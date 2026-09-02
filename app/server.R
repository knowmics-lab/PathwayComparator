function(input, output, session) {
  
  options(shiny.maxRequestSize=200*1024^2)
  
  rv <- reactiveValues(data.list = list(), common.pathways = NULL, list.all.nodes = NULL, compared = FALSE)
  output$hasFiles <- reactive({
    length(rv$data.list) > 0
  })
  outputOptions(output, "hasFiles", suspendWhenHidden = FALSE)
  
  # TRUE solo dopo che almeno un file e' stato caricato ed elaborato con
  # successo - controlla la visibilita' del wizard di ricerca.
  output$compared <- reactive({ rv$compared })
  outputOptions(output, "compared", suspendWhenHidden = FALSE)
  
  is.empty <- function(x) is.null(x) || length(x) == 0
  
  # Mostra/nasconde l'indicatore di attesa in ENTRAMBE le schede durante
  # un caricamento/eliminazione file. Usa shinyjs::show/hide (esecuzione
  # JS IMMEDIATA lato client) invece di un reactiveValues + conditionalPanel:
  # un valore reattivo che passa da TRUE a FALSE all'INTERNO dello stesso
  # ciclo reattivo (come accade qui, dato che l'intero
  # caricamento/eliminazione avviene in un solo observeEvent sincrono) non
  # produce MAI un cambiamento osservabile lato client - Shiny raggruppa
  # tutti i cambiamenti di un reactiveValues e invia al browser solo lo
  # stato FINALE del ciclo, mai gli stati intermedi. shinyjs::show/hide,
  # al contrario, invia il comando al client nel momento stesso in cui
  # viene chiamato (stesso meccanismo gia' usato con successo altrove in
  # questa app per "Updating visualization...").
  start.loading <- function() {
    shinyjs::hide("nodeComparisonContent")
    shinyjs::hide("networkViewContent")
    shinyjs::show("nodeComparisonLoading")
    shinyjs::show("networkViewLoading")
  }
  stop.loading <- function() {
    shinyjs::hide("nodeComparisonLoading")
    shinyjs::hide("networkViewLoading")
    if(isTRUE(rv$compared)) {
      shinyjs::show("nodeComparisonContent")
      shinyjs::show("networkViewContent")
    }
  }
  
  # Normalizza un campo arrivato da input$wizardSelection o
  # input$wizard_request_step2 (inviati dal widget JS) in un vero vettore
  # di caratteri, qualunque sia la forma con cui Shiny lo ha deserializzato
  # (un array JSON puo' arrivare come vettore o come lista a seconda del
  # percorso di codice interno di Shiny - un test isolato su
  # shiny:::safeFromJSON non si e' rivelato un indicatore affidabile del
  # comportamento REALE di input$xxx, come confermato da un errore
  # riscontrato in produzione: strsplit() falliva perche' il campo non
  # era di tipo character). Meglio normalizzare una volta sola qui, al
  # confine tra i dati del widget e il resto del codice R, che sperare che
  # ogni funzione a valle gestisca correttamente ogni possibile forma.
  to.char.vec <- function(x) {
    if(is.null(x)) return(NULL)
    x <- unlist(x, use.names = FALSE)
    if(length(x) == 0) return(NULL)
    as.character(x)
  }
  
  # TRUE quando il wizard ha inviato una selezione completa - controlla la
  # visibilita' dei parametri di visualizzazione (Layer, Show ...), che
  # devono comparire insieme al grafico, non subito dopo il caricamento.
  output$readyToPlot <- reactive({ !is.null(input$wizardSelection) })
  outputOptions(output, "readyToPlot", suspendWhenHidden = FALSE)
  
  uploaded.files <- reactive(input$uploadedFiles)
  
  shinyjs::hide("updatingMsg")
  
  shinyInput <- function(FUN, len, id, ...) {
    inputs <- character(len)
    for (i in seq_len(len)) {
      inputs[i] <- as.character(FUN(paste0(id, i), ...))
    }
    inputs
  }
  
  # Traduce le checkbox "Show compounds/drugs/miRNAs" (di default tutte NON
  # selezionate, quindi tutto nascosto) nel formato "Hide ..." atteso da
  # apply.hide.filters() in PathwaysFunctions.R.
  resolve.hide.elements <- function(show.elements) {
    hide <- c()
    if(!("Show compounds" %in% show.elements)) hide <- c(hide, "Hide chemical entities")
    if(!("Show drugs" %in% show.elements))     hide <- c(hide, "Hide drugs")
    if(!("Show miRNAs" %in% show.elements))    hide <- c(hide, "Hide miRNAs")
    hide
  }
  
  #---------------------------------------------------------------------
  # Conversione delle scelte R (vettori con nome, o liste per network) nel
  # formato JSON-friendly [{value,label}, ...] atteso dal widget JS
  # (www/wizard.js).
  #---------------------------------------------------------------------
  
  # Per le PATHWAY: value = nome completo, label = nome abbreviato se
  # troppo lungo. Ri-tronchiamo qui direttamente sul valore (non ci
  # affidiamo ai nomi eventualmente persi da funzioni come intersect()),
  # dato che per le pathway value e label derivano dalla STESSA stringa.
  pathways.to.js.choices <- function(pathway.names) {
    if(is.null(pathway.names) || length(pathway.names) == 0) return(list())
    pathway.names <- unname(pathway.names)
    labels <- truncate.label(pathway.names)
    lapply(seq_along(pathway.names), function(i) list(value = pathway.names[i], label = labels[i]))
  }
  
  # Per i NODI: get.list.selectable.nodes()/filter.selectable.nodes()
  # restituiscono gia' una lista (una per network) di vettori con nome
  # (value = "nodeName\nnetwork\nnode.type", name = etichetta gia'
  # abbreviata) - qui li appiattiamo in un unico elenco.
  # Un nodo con lo stesso nome puo' comparire in piu' network (tipicamente
  # quando piu' file caricati si riferiscono allo stesso organismo, come
  # piu' condizioni/linee cellulari dello stesso esperimento) - in tal
  # caso selezionare l'occorrenza da UNO QUALUNQUE di quei network produce
  # lo stesso risultato finale, dato che l'espansione per ortologhi in
  # build.pathway.net()/expand.to.ortho.nodes() lavora a livello di
  # ORGANISMO (cerca il nodeName nei dati dell'organismo di riferimento),
  # non di singolo file - quindi mostrare la stessa etichetta piu' volte,
  # una per network, era ridondante e confondeva l'utente (come nella
  # versione precedente dell'app, un nodo compare una sola volta, senza
  # specificare il network tra parentesi). La deduplicazione avviene sul
  # nodeName ESTRATTO dal valore (prima di un'eventuale troncatura per la
  # lunghezza), non sull'etichetta gia' troncata, per evitare falsi
  # duplicati se due nomi diversi si troncassero nello stesso identico
  # modo.
  nodes.to.js.choices <- function(list.options) {
    result <- list()
    seen.names <- character(0)
    for(net in names(list.options)) {
      opts <- list.options[[net]]
      if(length(opts) == 0) next
      labels <- names(opts)
      for(i in seq_along(opts)) {
        value.i <- unname(opts[i])
        node.name <- strsplit(value.i, "\n", fixed = TRUE)[[1]][1]
        if(node.name %in% seen.names) next
        seen.names <- c(seen.names, node.name)
        result[[length(result)+1]] <- list(value = value.i, label = unname(labels[i]))
      }
    }
    result
  }
  
  # Invia al widget le scelte "statiche" correnti (tutti i nodi per i
  # layer selezionati, tutte le pathway comuni) - richiamata al
  # caricamento/eliminazione file E ogni volta che cambia networkSel
  # (dato che l'elenco dei nodi dipende da quali layer sono confrontati).
  # Il parametro "networks" e' esplicito (invece di leggere sempre
  # input$networkSel) perche' subito dopo updatePickerInput("networkSel",
  # selected=...) il valore scelto non e' ancora riflesso in
  # input$networkSel - serve un giro di andata e ritorno col client prima
  # che lo sia - quindi refresh.search.setup() passa qui direttamente
  # initial.networks, invece di rischiare di inviare un primo wizard_init
  # con la lista nodi vuota.
  send.static.choices <- function(networks = input$networkSel) {
    if(!isTRUE(rv$compared) || is.null(rv$list.all.nodes)) {
      session$sendCustomMessage("wizard_init", list(allNodes = list(), allPathways = list()))
      return(invisible(NULL))
    }
    all.nodes.filtered <- filter.selectable.nodes(rv$list.all.nodes, networks)
    session$sendCustomMessage("wizard_init", list(
      allNodes = nodes.to.js.choices(all.nodes.filtered),
      allPathways = pathways.to.js.choices(sort(rv$common.pathways))
    ))
  }
  
  # Ricalcola pathway comuni/nodi selezionabili in base al contenuto
  # CORRENTE di rv$data.list - richiamata sia dopo un caricamento file
  # riuscito sia dopo un'eliminazione, cosi' che il wizard compaia/si
  # aggiorni automaticamente non appena i file cambiano.
  refresh.search.setup <- function() {
    if(length(rv$data.list) == 0) {
      rv$common.pathways <- NULL
      rv$list.all.nodes <- NULL
      rv$compared <- FALSE
      updatePickerInput(session,"networkSel",choices=NULL,selected=NULL)
      send.static.choices(NULL)
      return(invisible(NULL))
    }
    
    ordered.data.list <- rv$data.list[order(sapply(rv$data.list, function(el) el$organism))]
    list.organism <- as.character(sapply(ordered.data.list,function(el){el$organism}))
    
    if(length(list.organism)>1){
      org.pairs <- apply(t(combn(sort(list.organism),2)),1,function(row){paste0(row,collapse="-")})
      for(pair in org.pairs) {
        if(!pair %in% names(ortho.list)) {
          ortho.list[[pair]] <<- readRDS(paste0("Data/Orthologs/",pair,".rds"))
        }
      }
    } else {
      pair <- paste0(list.organism,"-",list.organism)
      if(!pair %in% names(ortho.list)) {
        ortho.list[[pair]] <<- readRDS(paste0("Data/Orthologs/",pair,".rds"))
      }
    }
    
    list.pathways <- lapply(list.organism,function(org){unique(pathway.list[[org]]$pathwayName)})
    rv$common.pathways <- unique(as.character(Reduce(intersect,list.pathways)))
    rv$list.all.nodes <- get.list.selectable.nodes(sort(rv$common.pathways),names(rv$data.list),rv$data.list,pathway.list)
    rv$compared <- TRUE
    
    max.opt <- 3
    if(length(rv$data.list)==2)
      max.opt <- 2
    initial.networks <- names(rv$data.list)[1:max.opt]
    updatePickerInput(session,"networkSel",choices=names(rv$data.list),selected = initial.networks)
    send.static.choices(initial.networks)
  }
  
  # Quando cambia networkSel, i nodi "statici" vanno rimandati al widget
  # (dipendono da quali layer sono confrontati) - le pathway comuni non
  # dipendono da networkSel, ma le rimandiamo comunque per semplicita'
  # (send.static.choices() le include sempre insieme).
  observeEvent(input$networkSel, {
    send.static.choices()
  }, ignoreInit = TRUE)
  
  #---------------------------------------------------------------------
  # Risponde alle richieste del widget per le scelte dinamiche del
  # secondo passo: nodi contenuti in una pathway scelta (ricerca by
  # Pathway), oppure pathway contenenti i nodi scelti (ricerca by Node -
  # con la stessa logica di intersezione, non unione, gia' in uso prima:
  # una pathway compare solo se contiene ALMENO UN nodo sorgente E ALMENO
  # UN nodo destinazione, quando la modalita' e' "paths").
  #---------------------------------------------------------------------
  observeEvent(input$wizard_request_step2, {
    req <- input$wizard_request_step2
    if(is.null(req)) return()
    
    if(req$role == "nodesForPathway") {
      opts <- get.list.selectable.nodes(to.char.vec(req$value), input$networkSel, rv$data.list, pathway.list)
      choices <- nodes.to.js.choices(opts)
    } else {
      if(req$outputType == "neighbors") {
        genes <- to.char.vec(req$value)
        pw <- sort(get.list.selectable.pathways(genes, rv$data.list, pathway.list))
      } else {
        source.nodes <- to.char.vec(req$value$source)
        dest.nodes <- to.char.vec(req$value$dest)
        pathways.with.source <- get.list.selectable.pathways(source.nodes, rv$data.list, pathway.list)
        pathways.with.dest   <- get.list.selectable.pathways(dest.nodes, rv$data.list, pathway.list)
        pw <- sort(intersect(pathways.with.source, pathways.with.dest))
      }
      choices <- pathways.to.js.choices(pw)
    }
    
    session$sendCustomMessage("wizard_step2_choices", list(seq = req$seq, choices = choices))
  })
  
  output$listFiles <- DT::renderDT({
    if(length(uploaded.files())==0) {
      tab.data <- data.frame("  "=character()," "=character(),check.names = F)
    } else {
      tab.data <- data.frame("  "=uploaded.files()," "=rep("Remove",length(uploaded.files())),check.names = F)
      tab.data[," "] <- shinyInput(actionButton, length(uploaded.files()),'delete_',label = "",icon=icon("trash"),
                                   class = "btn-delete-row",
                                   onclick = paste0('Shiny.onInputChange( \"delete_button\" , this.id, {priority: \"event\"})'))
    }
    datatable(tab.data, options = list(dom="t",ordering=F, language = list(
      zeroRecords = "No files yet"),
      initComplete = JS(
        "function(settings, json) {",
        "$(this.api().table().header()).css({'background-color': '#2f2f38','color': '#e4e4e8'});",
        "$(this.api().table().body()).css({'background-color': '#2f2f38','color': '#e4e4e8'});",
        "}")
    ),rownames = F,escape=F)
  })
  
  observeEvent(input$delete_button, {
    start.loading()
    on.exit(stop.loading(), add = TRUE)
    selectedRow <- as.numeric(strsplit(input$delete_button, "_")[[1]][2])
    rem.organism <- rv$data.list[[selectedRow]]$organism
    freq.organisms <- table(sapply(rv$data.list,function(el){el$organism}))
    if(freq.organisms[rem.organism]==1) {
      metapathway.list[[rem.organism]] <<- NULL
      pathway.list[[rem.organism]] <<- NULL
    }
    rv$data.list[[selectedRow]] <- NULL
    updateSelectInput(session,"uploadedFiles",choices=input$uploadedFiles[-selectedRow],selected=input$uploadedFiles[-selectedRow])
    
    refresh.search.setup()
    shinyjs::hide("updatingMsg")
  })
  
  observeEvent(input$hiddenUpload,{
    file.list <- input$hiddenUpload
    req(file.list)
    start.loading()
    on.exit(stop.loading(), add = TRUE)
    dataname.list <- c()
    for(i in 1:nrow(file.list)) {
      file <- file.list[i,]
      dataname <- strsplit(file$name,"\\.")[[1]][1]
      file.data <- tryCatch(read.phensim.file(file$datapath), error = function(e) NULL)
      if(is.null(file.data) || length(file.data$organism) == 0 || is.na(file.data$organism)) {
        showNotification(
          paste0("Impossibile leggere '", file$name, "': formato non riconosciuto o organismo non supportato."),
          type = "error", duration = 8
        )
        next
      }
      existing.names <- c(names(rv$data.list), dataname.list)
      if(dataname %in% existing.names) {
        suffix <- 2
        candidate <- paste0(dataname,"_",suffix)
        while(candidate %in% existing.names) {
          suffix <- suffix + 1
          candidate <- paste0(dataname,"_",suffix)
        }
        dataname <- candidate
      }
      rv$data.list[[dataname]] <- file.data
      organism <- as.character(file.data$organism)
      if(!organism %in% names(metapathway.list)) {
        metapathway.list[[organism]] <<- readRDS(paste0("Data/Metapathways/",organism,".rds"))
        organism.pathways <- readRDS(paste0("Data/Pathways/",organism,".rds"))
        organism.pathways$node.type <- classify.node.type(organism.pathways$node)
        pathway.list[[organism]] <<- organism.pathways
      }
      dataname.list <- c(dataname.list,dataname)
    }
    updateSelectInput(session,"uploadedFiles",choices=c(input$uploadedFiles,dataname.list),selected=c(input$uploadedFiles,dataname.list))
    refresh.search.setup()
  })
  
  #---------------------------------------------------------------------
  # Selezione finale, ricevuta come UN SOLO evento dal widget quando tutti
  # i passi del flusso sono stati completati (vedi www/wizard.js).
  #---------------------------------------------------------------------
  plot.trigger <- reactive({
    sel <- input$wizardSelection
    list(pathways = if(is.null(sel)) NULL else to.char.vec(sel$pathway),
         genes = if(is.null(sel)) NULL else to.char.vec(sel$gene),
         source.genes = if(is.null(sel)) NULL else to.char.vec(sel$source),
         dest.genes = if(is.null(sel)) NULL else to.char.vec(sel$dest),
         networks = input$networkSel,
         hide = resolve.hide.elements(input$showElements),
         view.mode = if(is.null(sel)) "neighbors" else sel$outputType,
         max.hops = if(is.null(sel)) 1 else sel$maxHops,
         max.length = if(!is.null(sel) && isTRUE(sel$limitPathLength)) sel$maxPathLength else Inf,
         only.perturbed.paths = !is.null(sel) && isTRUE(sel$onlyPerturbedPaths),
         data.list = rv$data.list)
  })
  plot.trigger.d <- debounce(plot.trigger, 400)
  
  observeEvent(plot.trigger(), {
    shinyjs::disable("networkSel")
    shinyjs::disable("showElements")
    shinyjs::runjs("$('#networkSel').selectpicker('refresh');")
    shinyjs::show("updatingMsg")
  }, ignoreInit = TRUE)
  
  legend.info <- reactiveVal(NULL)
  
  output$plotPathway <- renderVisNetwork({
    trig <- plot.trigger.d()
    on.exit({
      shinyjs::enable("networkSel")
      shinyjs::enable("showElements")
      shinyjs::runjs("$('#networkSel').selectpicker('refresh');")
      shinyjs::hide("updatingMsg")
    })
    resolved.view.mode <- if(is.null(trig$view.mode)) "neighbors" else trig$view.mode
    resolved.max.hops <- if(is.null(trig$max.hops) || is.na(trig$max.hops) || trig$max.hops < 1) 1 else trig$max.hops
    resolved.max.length <- if(is.null(trig$max.length) || is.na(trig$max.length)) Inf else trig$max.length
    ready <- if(resolved.view.mode=="paths") {
      !is.empty(trig$pathways) && !is.empty(trig$source.genes) && !is.empty(trig$dest.genes) && !is.empty(trig$networks)
    } else {
      !is.empty(trig$pathways) && !is.empty(trig$genes) && !is.empty(trig$networks)
    }
    if(ready){
      multilayer.net <- build.pathway.net(trig$data.list,metapathway.list,pathway.list,ortho.list,
                                          trig$networks,trig$pathways,trig$genes,trig$hide,
                                          view.mode = resolved.view.mode,
                                          source.genes = trig$source.genes, dest.genes = trig$dest.genes,
                                          max.hops = resolved.max.hops, max.length = resolved.max.length,
                                          only.perturbed.paths = isTRUE(trig$only.perturbed.paths))
      pathway.plot <- plot.pathway(multilayer.net$nodes,multilayer.net$edges,view.mode = resolved.view.mode)
      legend.info(if(nrow(multilayer.net$nodes)>0) c(get.legend.info(multilayer.net$nodes), list(view.mode=resolved.view.mode)) else NULL)
      pathway.plot
    } else {
      legend.info(NULL)
      NULL
    }
  })
  
  output$graphLegend <- renderUI({
    info <- legend.info()
    if(is.null(info)) return(NULL)
    tagList(
      lapply(seq_along(info$net.names), function(i) {
        tags$div(style = "display:flex; align-items:center; gap:8px; margin-bottom:8px;",
                 tags$div(style = paste0("width:16px; height:16px; border-radius:50%; flex-shrink:0; background:grey; border:3px solid ", info$colors[i], ";")),
                 tags$span(info$net.names[i], style = "font-size:13px; word-break:break-word;")
        )
      }),
      tags$div(style = "display:flex; align-items:center; gap:8px; margin-bottom:14px;",
               tags$div(style = "width:14px; height:14px; flex-shrink:0; background:white; border:3px solid black;"),
               tags$span("Endpoint", style = "font-size:13px;")
      ),
      if(identical(info$view.mode, "paths")) tagList(
        tags$div(style = "display:flex; align-items:center; gap:8px; margin-bottom:8px;",
                 tags$div(style = "width:14px; height:14px; border-radius:50%; flex-shrink:0; background:grey; border:4px solid #00A651;"),
                 tags$span("Source node", style = "font-size:13px;")
        ),
        tags$div(style = "display:flex; align-items:center; gap:8px; margin-bottom:14px;",
                 tags$div(style = "width:14px; height:14px; border-radius:50%; flex-shrink:0; background:grey; border:4px solid #ED1C24;"),
                 tags$span("Destination node", style = "font-size:13px;")
        )
      ),
      if(identical(info$view.mode, "neighbors")) tagList(
        tags$div(style = "display:flex; align-items:center; gap:8px; margin-bottom:14px;",
                 tags$div(style = "width:14px; height:14px; border-radius:50%; flex-shrink:0; background:grey; border:4px solid #8E44AD;"),
                 tags$span("Reference node", style = "font-size:13px;")
        )
      ),
      tags$div(
        tags$div("Score", style = "font-size:13px; text-align:center; margin-bottom:3px;"),
        tags$div(style = "width:100%; height:14px; border:1px solid #999; background:linear-gradient(to right, blue, grey, red);"),
        tags$div(style = "display:flex; justify-content:space-between; font-size:11px; margin-top:2px;",
                 tags$span("MIN"), tags$span("MAX")
        )
      )
    )
  })
  
  #---------------------------------------------------------------------
  # Scheda "Node comparison": tabella nodo x network con delta e colore -
  # vedi build.node.comparison.table()/compute.abs.max() in
  # PathwaysFunctions.R. Ricalcolata automaticamente ogni volta che
  # rv$data.list cambia (caricamento/eliminazione file) o che il filtro
  # pathway (input$nodeComparisonPathwayFilter) cambia - mostra tutti i
  # nodi di tutti i file correntemente caricati (ristretti alle pathway
  # selezionate, se presenti), non solo quelli perturbati.
  #---------------------------------------------------------------------
  observeEvent(rv$common.pathways, {
    choices <- if(is.null(rv$common.pathways)) list() else pathways.to.js.choices(sort(rv$common.pathways))
    session$sendCustomMessage("pathway_filter_init", list(choices = choices))
  }, ignoreNULL = FALSE)
  
  node.comparison.data <- reactive({
    build.node.comparison.table(rv$data.list, pathway.list, pathways = to.char.vec(input$nodeComparisonPathwayFilter))
  })
  
  output$nodeComparisonTable <- DT::renderDT({
    tbl <- node.comparison.data()
    if(nrow(tbl) == 0) {
      return(datatable(data.frame(Message = "No nodes to show (no files loaded, or no node matches the selected pathways)"), options = list(dom = "t"), rownames = FALSE))
    }
    
    abs.max <- compute.abs.max(tbl)
    score.cols <- setdiff(colnames(tbl), c("nodeName","delta","pathwaysList"))
    #Indici di colonna 0-based per DT (rownames=FALSE, quindi la colonna 0
    #e' "nodeName"): servono per applicare la colorazione SOLO alle
    #colonne dei punteggi, non a "Node" o "Delta". "pathwaysList" resta
    #nella tabella dati (necessaria per il tooltip sul nome del nodo) ma
    #viene nascosta dalla visualizzazione.
    score.col.indices <- match(score.cols, colnames(tbl)) - 1
    delta.col.index <- match("delta", colnames(tbl)) - 1
    node.col.index <- match("nodeName", colnames(tbl)) - 1
    pathways.col.index <- match("pathwaysList", colnames(tbl)) - 1
    
    display.tbl <- tbl
    display.tbl[,c(score.cols,"delta")] <- lapply(display.tbl[,c(score.cols,"delta"),drop=FALSE], function(x) round(x,2))
    
    #Colorazione delle celle: stessa interpolazione blue/grey/red gia'
    #usata per i nodi nel grafico della rete (vedi build.pathway.net,
    #"Set node colors"), ma con una trasformazione a radice quadrata
    #(segno preservato) prima della normalizzazione - riduce il peso
    #visivo dei rari valori estremi e aumenta la distinguibilita' dei
    #valori comuni piu' piccoli, restando pero' una scala ASSOLUTA (stesso
    #valore = sempre stesso colore, indipendentemente da pagina/ricerca
    #correnti - abs.max e' calcolato sull'INTERA tabella, non sulla
    #pagina mostrata). Il numero non viene mostrato direttamente nella
    #cella (solo il colore): resta leggibile al passaggio del mouse
    #tramite un tooltip personalizzato (non l'attributo HTML "title": il
    #tooltip nativo del browser e' un elemento di sistema, il cui font
    #NON e' controllabile via CSS - da qui la scelta di un tooltip fatto
    #a mano, per poterne controllare la dimensione del testo). Un bordo
    #bianco attorno a ciascuna cella separa visivamente le celle
    #adiacenti nella stessa riga, che altrimenti si fondono quando hanno
    #colori simili.
    color.js <- sprintf("
      function(td, cellData, rowData, row, col) {
        $(td).css('border', '3px solid #fff');
        $(td).css('position', 'relative');
        if (cellData === null || cellData === undefined) {
          $(td).css('background-color', '#fff');
          td.innerHTML = '<span style=\"color:#aaa; font-style:italic; font-size:12px;\">Not present</span>';
          return;
        }
        var v = parseFloat(cellData);
        var ABS_MAX = %s;
        var sign = v < 0 ? -1 : (v > 0 ? 1 : 0);
        var transformed = sign * Math.sqrt(Math.abs(v));
        var transformedMax = Math.sqrt(ABS_MAX);
        var t = Math.max(0, Math.min(1, (transformed + transformedMax) / (2*transformedMax)));
        var r,g,b;
        if (t <= 0.5) {
          var k = t / 0.5;
          r = Math.round(k*190); g = Math.round(k*190); b = Math.round(255 + k*(190-255));
        } else {
          var k = (t-0.5) / 0.5;
          r = Math.round(190 + k*(255-190)); g = Math.round(190 - k*190); b = Math.round(190 - k*190);
        }
        $(td).css('background-color', 'rgb('+r+','+g+','+b+')');
        $(td).addClass('node-cmp-score-cell');
        td.innerHTML = '<span class=\"node-cmp-tooltip\">' + v.toFixed(2) + '</span>';
      }", abs.max)
    
    #Elenco delle pathway a cui appartiene il nodo (colonna "Node"),
    #letto dalla colonna nascosta "pathwaysList" (indice %s) tramite
    #rowData - DataTables include SEMPRE tutte le colonne in rowData,
    #anche quelle nascoste dalla visualizzazione (visible=FALSE non le
    #rimuove dai dati, solo dal rendering).
    #
    #BUG FIX (segnalato in produzione): un tooltip al passaggio del mouse
    #con contenuto scorrevole e' problematico quando si trova dentro un
    #contenitore che ha GIA' un proprio scroll orizzontale (qui, la
    #tabella) - il gesto di scorrimento del mouse viene conteso tra i due
    #scroll annidati, rendendo difficile scorrere l'elenco. Sostituito
    #con un popover attivato al CLICK, in position:fixed (ancorato alla
    #finestra del browser, non alla tabella): il suo scroll verticale e'
    #cosi' completamente indipendente da qualunque contenitore scorrevole
    #della pagina. Un unico elemento popover viene creato una volta sola
    #e riusato per tutte le celle (evita di duplicare l'elemento per ogni
    #riga), con un pulsante di chiusura esplicito e chiusura al click
    #fuori dal popover.
    node.tooltip.js <- sprintf("
      function(td, cellData, rowData, row, col) {
        var pathwaysList = rowData[%s];
        if (!pathwaysList) { return; }
        $(td).addClass('node-cmp-name-cell');
        var pathwayNames = pathwaysList.split('|||');
        td.innerHTML = cellData + ' <span class=\"node-cmp-pathway-badge\">' + pathwayNames.length + '</span>';
        td.onclick = function(ev) {
          ev.stopPropagation();
          var pop = document.getElementById('nodeCmpPathwayPopover');
          if (!pop) {
            pop = document.createElement('div');
            pop.id = 'nodeCmpPathwayPopover';
            pop.className = 'node-cmp-pathway-popover';
            pop.innerHTML = '<div class=\"node-cmp-popover-header\">' +
              '<span class=\"node-cmp-popover-title\"></span>' +
              '<button type=\"button\" class=\"node-cmp-popover-close\">\\u00d7</button></div>' +
              '<div class=\"node-cmp-popover-body\"></div>';
            document.body.appendChild(pop);
            pop.querySelector('.node-cmp-popover-close').addEventListener('click', function(e) {
              e.stopPropagation();
              pop.style.display = 'none';
            });
            document.addEventListener('click', function(e) {
              if (pop.style.display === 'block' && !pop.contains(e.target)) {
                pop.style.display = 'none';
              }
            });
          }
          pop.querySelector('.node-cmp-popover-title').textContent = cellData + ' \\u2014 pathways';
          pop.querySelector('.node-cmp-popover-body').innerHTML = pathwayNames.join('<br>');
          pop.style.display = 'block';
          pop.style.visibility = 'hidden';
          var rect = td.getBoundingClientRect();
          var popRect = pop.getBoundingClientRect();
          var left = Math.min(rect.left, window.innerWidth - popRect.width - 10);
          left = Math.max(10, left);
          pop.style.left = left + 'px';
          if (rect.bottom + popRect.height + 8 > window.innerHeight) {
            pop.style.top = 'auto';
            pop.style.bottom = (window.innerHeight - rect.top + 4) + 'px';
          } else {
            pop.style.bottom = 'auto';
            pop.style.top = (rect.bottom + 4) + 'px';
          }
          pop.style.visibility = 'visible';
        };
      }", pathways.col.index)
    
    datatable(display.tbl,
              colnames = c("Node" = "nodeName", "Delta" = "delta"),
              rownames = FALSE,
              filter = "none",
              selection = "none",
              options = list(
                pageLength = 20,
                lengthMenu = list(c(20,50,100), c("20","50","100")),
                order = list(list(delta.col.index, "desc")),
                columnDefs = list(
                  list(targets = score.col.indices, createdCell = JS(color.js)),
                  list(targets = node.col.index, createdCell = JS(node.tooltip.js)),
                  list(targets = pathways.col.index, visible = FALSE)
                )
              ))
  }, server = TRUE)
  
}
