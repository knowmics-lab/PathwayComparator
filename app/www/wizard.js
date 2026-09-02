/* ==========================================================================
   Search Wizard - stato e logica.

   Sostituisce la cascata di otto pickerInput di shinyWidgets con un unico
   componente che gestisce internamente tutto il flusso (Search by ->
   Output type -> Pathway/Node/Source/Destination) e comunica con Shiny
   SOLO in due momenti:
     1) quando serve la lista dinamica del secondo passo (richiesta),
     2) quando la selezione e' completa (un solo evento finale).
   Questo elimina l'intera classe di problemi incontrati con otto input
   Shiny separati che si osservavano a vicenda.
   ========================================================================== */

(function () {
  "use strict";

  var CONTAINER_ID = "searchWizard";
  var MAX_RENDERED_OPTIONS = 200;

  // ------------------------------------------------------------------
  // Stato del widget. "step" e' l'indice del passo attivo nella barra di
  // avanzamento: 1 = scelta modalita', 2 = primo passo di selezione,
  // 3 = secondo passo di selezione (solo se il piano ne prevede due).
  // ------------------------------------------------------------------
  var state = {
    step: 1,
    searchMode: "Pathway",
    outputType: "neighbors",
    maxHops: 1,
    limitPathLength: false,
    maxPathLength: 3,
    onlyPerturbedPaths: false,
    allNodes: [],
    allPathways: [],
    selections: {},       // es. {pathway: [...]},{gene:[...]},{sourceDest:{source:[...],dest:[...]}}
    dynamicChoices: null,  // scelte per il passo dinamico, una volta arrivate dal server
    loadingDynamic: false,
    requestSeq: 0,         // per ignorare risposte del server ormai superate
    searchText: {}         // testo di ricerca corrente per ciascun picker (per chiave)
  };

  // ------------------------------------------------------------------
  // Piano dei passi per la combinazione (searchMode, outputType)
  // corrente. Ogni voce descrive UNA selezione da fare.
  // ------------------------------------------------------------------
  function getStepPlan(searchMode, outputType) {
    if (searchMode === "Pathway") {
      if (outputType === "neighbors") {
        return [
          { key: "pathway", label: "Pathway", kind: "single", source: "staticPathways" },
          { key: "gene", label: "Node", kind: "single", source: "dynamicNodes" }
        ];
      }
      return [
        { key: "pathway", label: "Pathway", kind: "single", source: "staticPathways" },
        { key: "sourceDest", label: "Source & destination nodes", kind: "double", source: "dynamicNodes" }
      ];
    }
    if (outputType === "neighbors") {
      return [
        { key: "gene", label: "Node", kind: "single", source: "staticNodes" },
        { key: "pathway", label: "Pathway", kind: "single", source: "dynamicPathways" }
      ];
    }
    return [
      { key: "sourceDest", label: "Source & destination nodes", kind: "double", source: "staticNodes" },
      { key: "pathway", label: "Pathway", kind: "single", source: "dynamicPathways" }
    ];
  }

  function currentPlan() {
    return getStepPlan(state.searchMode, state.outputType);
  }

  // Una selezione (single o double) e' "completa" se ha almeno un
  // elemento (double: almeno un elemento in ENTRAMBE le liste).
  function isSelectionComplete(planEntry, value) {
    if (!value) return false;
    if (planEntry.kind === "double") {
      return !!(value.source && value.source.length > 0 && value.dest && value.dest.length > 0);
    }
    return Array.isArray(value) && value.length > 0;
  }

  // ------------------------------------------------------------------
  // Utility DOM
  // ------------------------------------------------------------------
  function el(tag, attrs, children) {
    var e = document.createElement(tag);
    attrs = attrs || {};
    Object.keys(attrs).forEach(function (k) {
      if (k === "class") e.className = attrs[k];
      else if (k === "text") e.textContent = attrs[k];
      else if (k.indexOf("data-") === 0) e.setAttribute(k, attrs[k]);
      else e[k] = attrs[k];
    });
    (children || []).forEach(function (c) { if (c) e.appendChild(c); });
    return e;
  }

  // ------------------------------------------------------------------
  // Riepilogo testuale di una selezione, per i chip dei passi completati.
  // ------------------------------------------------------------------
  function summaryLabel(planEntry, value, choiceList) {
    var labelFor = function (v) {
      var found = (choiceList || []).find(function (c) { return c.value === v; });
      return found ? found.label : v;
    };
    if (planEntry.kind === "double") {
      var srcTxt = (value.source || []).map(labelFor).join(", ");
      var dstTxt = (value.dest || []).map(labelFor).join(", ");
      return "Source: " + srcTxt + "  |  Destination: " + dstTxt;
    }
    return (value || []).map(labelFor).join(", ");
  }

  // ------------------------------------------------------------------
  // Sotto-componente: lista con ricerca e selezione multipla.
  // opts: {choices:[{value,label}], selected:[...], onChange:fn, placeholder}
  // ------------------------------------------------------------------
  function renderListPicker(key, choices, selected, onChange, placeholder) {
    var wrap = el("div", { class: "wiz-list-picker-inner" });
    var searchTxt = state.searchText[key] || "";

    var tagsRow = el("div", { class: "wiz-selected-tags" });
    selected.forEach(function (v) {
      var choice = choices.find(function (c) { return c.value === v; });
      var label = choice ? choice.label : v;
      var removeBtn = el("button", { class: "wiz-tag-remove", type: "button", text: "\u00d7" });
      removeBtn.addEventListener("click", function (ev) {
        ev.stopPropagation();
        onChange(selected.filter(function (x) { return x !== v; }));
      });
      tagsRow.appendChild(el("span", { class: "wiz-tag" }, [
        el("span", { class: "wiz-tag-text", text: label }),
        removeBtn
      ]));
    });
    if (selected.length > 0) wrap.appendChild(tagsRow);

    var searchInput = el("input", {
      class: "wiz-search-input", type: "text",
      placeholder: placeholder || "Search...", value: searchTxt
    });
    searchInput.addEventListener("input", function (ev) {
      state.searchText[key] = ev.target.value;
      render();
      // Ripristina il focus e la posizione del cursore: render() ridisegna
      // l'intero widget, quindi il campo va ri-selezionato esplicitamente.
      var refreshed = document.querySelector('[data-search-key="' + key + '"]');
      if (refreshed) {
        refreshed.focus();
        var pos = refreshed.value.length;
        refreshed.setSelectionRange(pos, pos);
      }
    });
    searchInput.setAttribute("data-search-key", key);
    wrap.appendChild(searchInput);

    var filtered = choices.filter(function (c) {
      return !searchTxt || c.label.toLowerCase().indexOf(searchTxt.toLowerCase()) !== -1;
    });
    var shown = filtered.slice(0, MAX_RENDERED_OPTIONS);

    var listEl = el("div", { class: "wiz-option-list" });
    if (shown.length === 0) {
      listEl.appendChild(el("div", { class: "wiz-option-empty", text: "No matching options" }));
    } else {
      shown.forEach(function (c) {
        var isSel = selected.indexOf(c.value) !== -1;
        var checkbox = el("input", { type: "checkbox" });
        checkbox.checked = isSel;
        var row = el("div", { class: "wiz-option-row" + (isSel ? " wiz-option-selected" : "") }, [
          checkbox,
          el("span", { class: "wiz-option-text", text: c.label })
        ]);
        row.addEventListener("click", function () {
          if (isSel) onChange(selected.filter(function (x) { return x !== c.value; }));
          else onChange(selected.concat([c.value]));
        });
        listEl.appendChild(row);
      });
      if (filtered.length > MAX_RENDERED_OPTIONS) {
        listEl.appendChild(el("div", {
          class: "wiz-option-more",
          text: "+" + (filtered.length - MAX_RENDERED_OPTIONS) + " more \u2014 keep typing to narrow down"
        }));
      }
    }
    wrap.appendChild(listEl);

    var actions = el("div", { class: "wiz-list-actions" }, [
      el("span", { class: "wiz-list-count", text: selected.length + " selected" })
    ]);
    if (selected.length > 0) {
      var clearBtn = el("button", { class: "wiz-clear-link", type: "button", text: "Clear" });
      clearBtn.addEventListener("click", function () { onChange([]); });
      actions.appendChild(clearBtn);
    }
    wrap.appendChild(actions);

    return wrap;
  }

  // ------------------------------------------------------------------
  // Corpo del passo 1: scelta di Search by / Output type.
  // ------------------------------------------------------------------
  function renderModeToggle() {
    function toggleGroup(label, options, currentValue, onSelect) {
      var btns = options.map(function (opt) {
        var btn = el("button", {
          class: "wiz-toggle-btn" + (opt.value === currentValue ? " selected" : ""),
          type: "button", text: opt.label
        });
        btn.addEventListener("click", function () { onSelect(opt.value); });
        return btn;
      });
      return el("div", { class: "wiz-toggle-group" }, [
        el("span", { class: "wiz-toggle-group-label", text: label }),
        el("div", { class: "wiz-toggle-btns" }, btns)
      ]);
    }

    var searchByGroup = toggleGroup("Search by", [
      { value: "Pathway", label: "Pathway" }, { value: "Node", label: "Node" }
    ], state.searchMode, function (v) {
      if (v === state.searchMode) return;
      state.searchMode = v;
      if (state.step === 1) { state.selections = {}; state.dynamicChoices = null; state.loadingDynamic = false; }
      else resetFromStep(2);
      render();
    });

    var outputTypeGroup = toggleGroup("Output type", [
      { value: "neighbors", label: "Ego-network" }, { value: "paths", label: "All paths" }
    ], state.outputType, function (v) {
      if (v === state.outputType) return;
      state.outputType = v;
      if (state.step === 1) { state.selections = {}; state.dynamicChoices = null; state.loadingDynamic = false; }
      else resetFromStep(2);
      render();
    });

    return el("div", { class: "wiz-toggle-row" }, [searchByGroup, outputTypeGroup]);
  }

  // ------------------------------------------------------------------
  // Riga "Max hops" / "Limit max path length to": SEPARATA dal flusso a
  // passi e SEMPRE visibile (in base al solo outputType corrente), cosi'
  // regolarle non tocca mai la selezione di pathway/nodi gia' fatta - se
  // la selezione e' gia' completa, il grafico si aggiorna dal vivo con il
  // nuovo valore.
  // ------------------------------------------------------------------
  function renderSettingsRow() {
    var row = el("div", { class: "wiz-toggle-row", style: "margin-top:6px;" });
    if (state.outputType === "neighbors") {
      var hopsInput = el("input", { type: "number", min: "1", max: "20", step: "1", value: state.maxHops });
      hopsInput.style.width = "60px";
      hopsInput.addEventListener("change", function (ev) {
        var v = parseInt(ev.target.value, 10);
        state.maxHops = (isNaN(v) || v < 1) ? 1 : v;
        maybeEmitSelection();
      });
      row.appendChild(el("div", { class: "wiz-toggle-group" }, [
        el("span", { class: "wiz-toggle-group-label", text: "Max hops" }),
        hopsInput
      ]));
    } else {
      var checkbox = el("input", { type: "checkbox" });
      checkbox.checked = state.limitPathLength;
      var lenInput = el("input", { type: "number", min: "1", max: "50", step: "1", value: state.maxPathLength });
      lenInput.style.width = "60px";
      lenInput.style.marginLeft = "6px";
      lenInput.disabled = !state.limitPathLength;
      checkbox.addEventListener("change", function (ev) {
        state.limitPathLength = ev.target.checked;
        render();
        maybeEmitSelection();
      });
      lenInput.addEventListener("change", function (ev) {
        var v = parseInt(ev.target.value, 10);
        state.maxPathLength = (isNaN(v) || v < 1) ? 1 : v;
        maybeEmitSelection();
      });
      var group = el("div", { class: "wiz-toggle-group" }, [
        el("label", { style: "display:flex; align-items:center; gap:6px; cursor:pointer;" }, [
          checkbox, el("span", { text: "Limit max path length to" })
        ]),
        lenInput
      ]);
      row.appendChild(group);

      // Opzione indipendente (si puo' combinare con il limite sulla
      // lunghezza sopra, non lo sostituisce): mostra solo i cammini che
      // passano per almeno un nodo con score diverso da zero, per
      // evidenziare solo i percorsi che coinvolgono nodi effettivamente
      // perturbati.
      var perturbedCheckbox = el("input", { type: "checkbox" });
      perturbedCheckbox.checked = state.onlyPerturbedPaths;
      perturbedCheckbox.addEventListener("change", function (ev) {
        state.onlyPerturbedPaths = ev.target.checked;
        maybeEmitSelection();
      });
      row.appendChild(el("div", { class: "wiz-toggle-group" }, [
        el("label", { style: "display:flex; align-items:center; gap:6px; cursor:pointer;" }, [
          perturbedCheckbox, el("span", { text: "Only show paths to perturbed nodes" })
        ])
      ]));
    }
    return row;
  }

  // ------------------------------------------------------------------
  // Richiede al server le scelte dinamiche per il secondo passo, dato il
  // valore scelto nel primo. Ignora risposte fuori sequenza (se l'utente
  // cambia rapidamente selezione prima che arrivi una risposta precedente).
  // ------------------------------------------------------------------
  function requestDynamicChoices(plan, firstStepValue) {
    state.loadingDynamic = true;
    state.dynamicChoices = null;
    state.requestSeq += 1;
    var seq = state.requestSeq;
    window.__wizardLastSeq = seq;
    Shiny.setInputValue("wizard_request_step2", {
      searchMode: state.searchMode,
      outputType: state.outputType,
      role: plan[1].source === "dynamicNodes" ? "nodesForPathway" : "pathwaysForNodes",
      value: firstStepValue,
      seq: seq
    }, { priority: "event" });
  }

  // ------------------------------------------------------------------
  // Azzera le selezioni a partire da un dato passo (incluso), quando
  // l'utente torna indietro o cambia modalita'. Se preserveOwnValue e'
  // true, mantiene il valore GIA' scelto per il passo a cui si torna
  // (cosi' l'utente lo ritrova ancora spuntato, invece di ripartire da
  // zero) - usato quando si clicca esplicitamente "change" su un passo
  // per modificarlo, non quando cambia la modalita' stessa (che rende
  // irrilevante qualunque selezione precedente).
  // ------------------------------------------------------------------
  function resetFromStep(fromStep, preserveOwnValue) {
    if (preserveOwnValue) {
      var plan = currentPlan();
      var planIdx = fromStep - 2;
      var key = (planIdx >= 0 && plan[planIdx]) ? plan[planIdx].key : null;
      var preserved = key ? state.selections[key] : null;
      state.selections = {};
      if (key) state.selections[key] = preserved;
    } else if (fromStep <= 2) {
      state.selections = {};
    }
    state.dynamicChoices = null;
    state.loadingDynamic = false;
    state.step = fromStep;
  }

  // ------------------------------------------------------------------
  // Chiamata quando una selezione (primo o secondo passo) cambia. NON
  // avanza mai automaticamente al passo successivo (lo farebbe gia' alla
  // prima voce selezionata, impedendo di sceglierne piu' di una - i menu
  // sono tutti a selezione multipla): l'avanzamento tra passi avviene solo
  // tramite il pulsante "Next" esplicito (vedi goToNextStep). Per l'ULTIMO
  // passo del piano, invece, non serve alcun pulsante: la selezione finale
  // viene inviata a Shiny a ogni modifica, cosi' il grafico si aggiorna
  // via via che l'utente affina la scelta (stesso comportamento della
  // versione precedente dell'app).
  // ------------------------------------------------------------------
  function onStepValueChange(stepIndex, planEntry, newValue) {
    state.selections[planEntry.key] = newValue;
    render();
    var plan = currentPlan();
    if (stepIndex === plan.length - 1) maybeEmitSelection();
  }

  // ------------------------------------------------------------------
  // Passa esplicitamente al passo successivo (pulsante "Next"), avviando
  // la richiesta delle scelte dinamiche per il passo che si apre, se
  // necessario.
  // ------------------------------------------------------------------
  function goToNextStep(stepIndex, planEntry) {
    var plan = currentPlan();
    var value = state.selections[planEntry.key];
    var nextEntry = plan[stepIndex + 1];
    if (nextEntry && (nextEntry.source === "dynamicNodes" || nextEntry.source === "dynamicPathways")) {
      requestDynamicChoices(plan, value);
    }
    state.step = stepIndex + 3; // stepIndex 0 -> passo 3 (il primo passo di selezione e' il 2)
    render();
  }

  // ------------------------------------------------------------------
  // Se tutti i passi del piano corrente sono completi, invia la
  // selezione finale a Shiny (un solo evento).
  // ------------------------------------------------------------------
  function maybeEmitSelection() {
    var plan = currentPlan();
    var allDone = plan.every(function (p) { return isSelectionComplete(p, state.selections[p.key]); });
    if (!allDone) {
      Shiny.setInputValue("wizardSelection", null);
      return;
    }
    var payload = {
      searchMode: state.searchMode,
      outputType: state.outputType,
      pathway: state.selections.pathway || null,
      gene: state.selections.gene || null,
      source: state.selections.sourceDest ? state.selections.sourceDest.source : null,
      dest: state.selections.sourceDest ? state.selections.sourceDest.dest : null,
      maxHops: state.maxHops,
      limitPathLength: state.limitPathLength,
      maxPathLength: state.maxPathLength,
      onlyPerturbedPaths: state.onlyPerturbedPaths
    };
    Shiny.setInputValue("wizardSelection", payload, { priority: "event" });
  }

  // ------------------------------------------------------------------
  // Render di un passo di tipo "single" (una sola lista).
  // ------------------------------------------------------------------
  function renderSingleStep(stepIndex, planEntry, choices) {
    var current = state.selections[planEntry.key] || [];
    return renderListPicker(planEntry.key, choices, current, function (newSel) {
      onStepValueChange(stepIndex, planEntry, newSel);
    }, "Search " + planEntry.label.toLowerCase() + "...");
  }

  // ------------------------------------------------------------------
  // Render di un passo di tipo "double" (Source + Destination affiancati).
  // ------------------------------------------------------------------
  function renderDoubleStep(stepIndex, planEntry, choices) {
    var current = state.selections[planEntry.key] || { source: [], dest: [] };
    var container = el("div", { class: "wiz-list-picker wiz-list-double" });

    var srcCol = el("div", { class: "wiz-list-picker-col" }, [
      el("div", { class: "wiz-list-col-label", text: "Source node" }),
      renderListPicker(planEntry.key + "_src", choices, current.source || [], function (newSel) {
        onStepValueChange(stepIndex, planEntry, { source: newSel, dest: current.dest || [] });
      }, "Search source node...")
    ]);
    var dstCol = el("div", { class: "wiz-list-picker-col" }, [
      el("div", { class: "wiz-list-col-label", text: "Destination node" }),
      renderListPicker(planEntry.key + "_dst", choices, current.dest || [], function (newSel) {
        onStepValueChange(stepIndex, planEntry, { source: current.source || [], dest: newSel });
      }, "Search destination node...")
    ]);
    container.appendChild(srcCol);
    container.appendChild(dstCol);
    return container;
  }

  function choicesForSource(source) {
    if (source === "staticNodes") return state.allNodes;
    if (source === "staticPathways") return state.allPathways;
    // dynamicNodes / dynamicPathways
    return state.dynamicChoices || [];
  }

  // ------------------------------------------------------------------
  // Render del corpo di un passo di selezione (single o double), con
  // stato di caricamento se le scelte sono dinamiche e non ancora arrivate.
  // Se NON e' l'ultimo passo del piano, aggiunge un pulsante "Next"
  // esplicito (abilitato solo a selezione completa) per avanzare - non si
  // avanza mai automaticamente alla prima voce scelta, per non impedire
  // la selezione multipla.
  // ------------------------------------------------------------------
  function renderSelectionStep(stepIndex, planEntry) {
    var isDynamic = planEntry.source === "dynamicNodes" || planEntry.source === "dynamicPathways";
    if (isDynamic && state.loadingDynamic) {
      return el("div", { class: "wiz-loading", text: "Loading " + planEntry.label.toLowerCase() + "..." });
    }
    var choices = choicesForSource(planEntry.source);
    var wrap = el("div", {});
    var body = planEntry.kind === "double"
      ? renderDoubleStep(stepIndex, planEntry, choices)
      : el("div", { class: "wiz-list-picker" }, [renderSingleStep(stepIndex, planEntry, choices)]);
    wrap.appendChild(body);

    var plan = currentPlan();
    if (stepIndex < plan.length - 1) {
      var value = state.selections[planEntry.key];
      var complete = isSelectionComplete(planEntry, value);
      var nextBtn = el("button", {
        class: "wiz-toggle-btn" + (complete ? " selected" : ""),
        type: "button", text: "Next", style: "margin-top:10px;"
      });
      if (!complete) { nextBtn.disabled = true; nextBtn.style.opacity = "0.5"; nextBtn.style.cursor = "default"; }
      else nextBtn.addEventListener("click", function () { goToNextStep(stepIndex, planEntry); });
      wrap.appendChild(nextBtn);
    }
    return wrap;
  }

  // ------------------------------------------------------------------
  // Riepilogo compatto (chip) per un passo gia' completato, cliccabile
  // per tornare indietro e modificarlo.
  // ------------------------------------------------------------------
  function renderCompletedSummary(stepIndex, planEntry) {
    var choices = choicesForSource(planEntry.source);
    var value = state.selections[planEntry.key];
    var text = summaryLabel(planEntry, value, choices);
    var row = el("div", { class: "wiz-completed-summary" }, [
      el("div", { class: "wiz-summary-chip" }, [
        el("span", { class: "wiz-chip-text", text: planEntry.label + ": " + text })
      ]),
      el("span", { class: "wiz-edit-hint", text: "change" })
    ]);
    row.addEventListener("click", function () {
      resetFromStep(stepIndex + 2, true);
      render();
    });
    return row;
  }

  // ------------------------------------------------------------------
  // Barra di avanzamento (cerchi numerati + linee di collegamento).
  // ------------------------------------------------------------------
  function renderProgress(plan) {
    var totalSteps = 1 + plan.length;
    var labels = ["Mode"].concat(plan.map(function (p) { return p.label; }));
    var bar = el("div", { class: "wiz-progress" });
    for (var i = 1; i <= totalSteps; i++) {
      var cls = "wiz-progress-step";
      if (i < state.step) cls += " completed";
      else if (i === state.step) cls += " active";
      var circleText = i < state.step ? "\u2713" : String(i);
      var stepDiv = el("div", { class: cls }, [
        el("div", { class: "wiz-progress-circle", text: circleText }),
        el("div", { class: "wiz-progress-label", text: labels[i - 1] })
      ]);
      if (i < state.step) {
        stepDiv.style.cursor = "pointer";
        stepDiv.addEventListener("click", function (clickedStep) {
          return function () { resetFromStep(clickedStep, clickedStep > 1); render(); };
        }(i));
      }
      bar.appendChild(stepDiv);
      if (i < totalSteps) bar.appendChild(el("div", { class: "wiz-progress-line" }));
    }
    return bar;
  }

  // ------------------------------------------------------------------
  // Render completo del widget.
  // ------------------------------------------------------------------
  function render() {
    var container = document.getElementById(CONTAINER_ID);
    if (!container) return;
    var plan = currentPlan();

    // Se il passo attivo supera il numero di passi disponibili (es. dopo
    // un cambio di modalita' che accorcia il piano), riportalo all'ultimo
    // passo valido.
    var totalSteps = 1 + plan.length;
    if (state.step > totalSteps) state.step = totalSteps;

    var root = el("div", { class: "search-wizard" });
    root.appendChild(renderProgress(plan));

    var body = el("div", { class: "wiz-step-body" });

    // Passo 1: Mode - riepilogo se gia' completato (sempre "completo" per
    // definizione, dato che ha dei default), altrimenti i toggle.
    if (state.step === 1) {
      body.appendChild(renderModeToggle());
      body.appendChild(renderSettingsRow());
      var nextBtn = el("button", {
        class: "wiz-toggle-btn selected", type: "button", text: "Next",
        style: "margin-top:12px;"
      });
      nextBtn.addEventListener("click", function () { state.step = 2; render(); });
      body.appendChild(nextBtn);
    } else {
      // Passo 1 completato: riepilogo cliccabile.
      var modeRow = el("div", { class: "wiz-completed-summary" }, [
        el("div", { class: "wiz-summary-chip" }, [
          el("span", {
            class: "wiz-chip-text",
            text: "Search by " + state.searchMode + " \u2014 " +
              (state.outputType === "neighbors" ? "Ego-network" : "All paths")
          })
        ]),
        el("span", { class: "wiz-edit-hint", text: "change" })
      ]);
      modeRow.addEventListener("click", function () { resetFromStep(1); render(); });
      body.appendChild(modeRow);
      body.appendChild(renderSettingsRow());

      // Passi di selezione: quelli PRIMA del passo attivo sono riepiloghi;
      // quello CORRENTE e' interattivo.
      for (var i = 0; i < plan.length; i++) {
        var stepNum = i + 2;
        if (stepNum < state.step) {
          body.appendChild(renderCompletedSummary(i, plan[i]));
        } else if (stepNum === state.step) {
          body.appendChild(el("div", { class: "wiz-step-title", text: plan[i].label }));
          body.appendChild(renderSelectionStep(i, plan[i]));
        }
      }
    }

    root.appendChild(body);

    var resetBtn = el("button", { class: "wiz-reset-link", type: "button", text: "Start over" });
    resetBtn.addEventListener("click", function () {
      state.step = 1;
      state.selections = {};
      state.dynamicChoices = null;
      state.loadingDynamic = false;
      state.searchText = {};
      render();
      Shiny.setInputValue("wizardSelection", null);
    });
    root.appendChild(resetBtn);

    container.innerHTML = "";
    container.appendChild(root);
  }

  // ------------------------------------------------------------------
  // Handler dei messaggi da Shiny.
  // ------------------------------------------------------------------
  document.addEventListener("DOMContentLoaded", function () {
    if (typeof Shiny === "undefined") return;

    Shiny.addCustomMessageHandler("wizard_init", function (msg) {
      state.allNodes = msg.allNodes || [];
      state.allPathways = msg.allPathways || [];
      state.step = 1;
      state.selections = {};
      state.dynamicChoices = null;
      state.loadingDynamic = false;
      state.searchText = {};
      render();
    });

    Shiny.addCustomMessageHandler("wizard_step2_choices", function (msg) {
      if (msg.seq !== window.__wizardLastSeq) return; // risposta superata
      state.dynamicChoices = msg.choices || [];
      state.loadingDynamic = false;
      render();
    });

    render();
  });
})();
