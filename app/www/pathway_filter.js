/* ==========================================================================
   Pathway Filter - selettore multi-scelta per la scheda "Node comparison",
   nello stesso stile visivo del wizard di ricerca (tag, ricerca testuale,
   lista con checkbox) - riusa DIRETTAMENTE le classi CSS gia' definite in
   wizard.css (.wiz-tag, .wiz-search-input, .wiz-option-list, ecc.), senza
   bisogno di un foglio di stile separato. Componente autonomo (non fa
   parte dello stato del wizard di ricerca): una sola lista persistente,
   senza passi, che invia la selezione a Shiny a ogni modifica.
   ========================================================================== */

(function () {
  "use strict";

  var CONTAINER_ID = "pathwayFilterWidget";
  var MAX_RENDERED_OPTIONS = 200;

  var state = {
    choices: [],
    selected: [],
    searchText: ""
  };

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

  function onChange(newSelected) {
    state.selected = newSelected;
    render();
    Shiny.setInputValue("nodeComparisonPathwayFilter", state.selected, { priority: "event" });
  }

  function render() {
    var container = document.getElementById(CONTAINER_ID);
    if (!container) return;

    var wrap = el("div", { class: "wiz-list-picker-inner" });

    var tagsRow = el("div", { class: "wiz-selected-tags" });
    state.selected.forEach(function (v) {
      var choice = state.choices.find(function (c) { return c.value === v; });
      var label = choice ? choice.label : v;
      var removeBtn = el("button", { class: "wiz-tag-remove", type: "button", text: "\u00d7" });
      removeBtn.addEventListener("click", function (ev) {
        ev.stopPropagation();
        onChange(state.selected.filter(function (x) { return x !== v; }));
      });
      tagsRow.appendChild(el("span", { class: "wiz-tag" }, [
        el("span", { class: "wiz-tag-text", text: label }),
        removeBtn
      ]));
    });
    if (state.selected.length > 0) wrap.appendChild(tagsRow);

    var searchInput = el("input", {
      class: "wiz-search-input", type: "text",
      placeholder: "Search pathway...", value: state.searchText
    });
    searchInput.addEventListener("input", function (ev) {
      state.searchText = ev.target.value;
      render();
      // render() ridisegna tutto il widget: il campo di ricerca va
      // ri-selezionato esplicitamente per non perdere il focus mentre si
      // digita (stesso accorgimento gia' usato nel wizard di ricerca).
      var refreshed = document.querySelector("#" + CONTAINER_ID + " .wiz-search-input");
      if (refreshed) {
        refreshed.focus();
        var pos = refreshed.value.length;
        refreshed.setSelectionRange(pos, pos);
      }
    });
    wrap.appendChild(searchInput);

    var filtered = state.choices.filter(function (c) {
      return !state.searchText || c.label.toLowerCase().indexOf(state.searchText.toLowerCase()) !== -1;
    });
    var shown = filtered.slice(0, MAX_RENDERED_OPTIONS);

    var listEl = el("div", { class: "wiz-option-list" });
    if (shown.length === 0) {
      listEl.appendChild(el("div", { class: "wiz-option-empty", text: "No matching pathways" }));
    } else {
      shown.forEach(function (c) {
        var isSel = state.selected.indexOf(c.value) !== -1;
        var checkbox = el("input", { type: "checkbox" });
        checkbox.checked = isSel;
        var row = el("div", { class: "wiz-option-row" + (isSel ? " wiz-option-selected" : "") }, [
          checkbox,
          el("span", { class: "wiz-option-text", text: c.label })
        ]);
        row.addEventListener("click", function () {
          if (isSel) onChange(state.selected.filter(function (x) { return x !== c.value; }));
          else onChange(state.selected.concat([c.value]));
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
      el("span", { class: "wiz-list-count", text: state.selected.length + " selected" })
    ]);
    if (state.selected.length > 0) {
      var clearBtn = el("button", { class: "wiz-clear-link", type: "button", text: "Clear" });
      clearBtn.addEventListener("click", function () { onChange([]); });
      actions.appendChild(clearBtn);
    }
    wrap.appendChild(actions);

    var outer = el("div", { class: "wiz-list-picker" }, [wrap]);
    container.innerHTML = "";
    container.appendChild(outer);
  }

  document.addEventListener("DOMContentLoaded", function () {
    if (typeof Shiny === "undefined") return;

    Shiny.addCustomMessageHandler("pathway_filter_init", function (msg) {
      state.choices = msg.choices || [];
      // Non azzera la selezione corrente quando arrivano scelte
      // aggiornate (es. cambio di networkSel non tocca le pathway
      // comuni) - ma rimuove eventuali valori selezionati che non sono
      // piu' tra le scelte disponibili (es. dopo l'eliminazione di un
      // file che cambia l'insieme di pathway comuni).
      var validValues = state.choices.map(function (c) { return c.value; });
      var stillValid = state.selected.filter(function (v) { return validValues.indexOf(v) !== -1; });
      if (stillValid.length !== state.selected.length) {
        state.selected = stillValid;
        Shiny.setInputValue("nodeComparisonPathwayFilter", state.selected, { priority: "event" });
      }
      render();
    });

    render();
  });
})();