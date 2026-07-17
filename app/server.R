function(input, output, session) {
  
  options(shiny.maxRequestSize=200*1024^2)
  
  rv <- reactiveValues(data.list = list(), common.pathways = NULL, list.all.nodes = NULL, compared = FALSE)
  output$hasFiles <- reactive({
    length(rv$data.list) > 0
  })
  outputOptions(output, "hasFiles", suspendWhenHidden = FALSE)
  output$compared <- reactive({ rv$compared })
  outputOptions(output, "compared", suspendWhenHidden = FALSE)
  is.empty <- function(x) is.null(x) || length(x) == 0
  
  output$readyToPlot <- reactive({
    if(is.null(input$outputType) || is.null(input$searchOpt)) return(FALSE)
    if(input$outputType == "paths") {
      !is.empty(last.sel.pathways()) && !is.empty(last.sel.source()) && !is.empty(last.sel.dest())
    } else {
      !is.empty(last.sel.pathways()) && !is.empty(last.sel.genes())
    }
  })
  outputOptions(output, "readyToPlot", suspendWhenHidden = FALSE)
  
  uploaded.files <- reactive(input$uploadedFiles)
  
  shinyjs::hide("updatingMsg")
  
  last.sel.pathways <- reactiveVal(NULL)
  last.sel.genes <- reactiveVal(NULL)
  last.sel.source <- reactiveVal(NULL)
  last.sel.dest <- reactiveVal(NULL)
  
  observe({
    if(is.null(input$searchOpt) || input$searchOpt == "Pathway") last.sel.pathways(input$pathwaySel)
    else last.sel.pathways(input$pathwaySelNode)
  })
  observe({
    if(is.null(input$searchOpt) || input$searchOpt == "Node") last.sel.genes(input$geneSel)
    else last.sel.genes(input$geneSelPathway)
  })
  observe({
    if(is.null(input$searchOpt) || input$searchOpt == "Node") last.sel.source(input$sourceSel)
    else last.sel.source(input$sourceSelPathway)
  })
  observe({
    if(is.null(input$searchOpt) || input$searchOpt == "Node") last.sel.dest(input$destSel)
    else last.sel.dest(input$destSelPathway)
  })
  
  shinyInput <- function(FUN, len, id, ...) {
    inputs <- character(len)
    for (i in seq_len(len)) {
      inputs[i] <- as.character(FUN(paste0(id, i), ...))
    }
    inputs
  }
  
  resolve.hide.elements <- function(show.elements) {
    hide <- c()
    if(!("Show compounds" %in% show.elements)) hide <- c(hide, "Hide chemical entities")
    if(!("Show drugs" %in% show.elements))     hide <- c(hide, "Hide drugs")
    if(!("Show miRNAs" %in% show.elements))    hide <- c(hide, "Hide miRNAs")
    hide
  }
  
  refresh.search.setup <- function() {
    if(length(rv$data.list) == 0) {
      rv$common.pathways <- NULL
      rv$list.all.nodes <- NULL
      rv$compared <- FALSE
      last.sel.pathways(NULL)
      last.sel.genes(NULL)
      last.sel.source(NULL)
      last.sel.dest(NULL)
      updatePickerInput(session,"pathwaySel",choices=NULL,selected=NULL)
      updatePickerInput(session,"pathwaySelNode",choices=NULL,selected=NULL)
      updatePickerInput(session,"geneSel",choices=NULL,selected=NULL)
      updatePickerInput(session,"geneSelPathway",choices=NULL,selected=NULL)
      updatePickerInput(session,"sourceSel",choices=NULL,selected=NULL)
      updatePickerInput(session,"sourceSelPathway",choices=NULL,selected=NULL)
      updatePickerInput(session,"destSel",choices=NULL,selected=NULL)
      updatePickerInput(session,"destSelPathway",choices=NULL,selected=NULL)
      updatePickerInput(session,"networkSel",choices=NULL,selected=NULL)
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
    
    last.sel.pathways(NULL)
    last.sel.genes(NULL)
    last.sel.source(NULL)
    last.sel.dest(NULL)
    updatePickerInput(session,"pathwaySel",choices=NULL,selected=NULL)
    updatePickerInput(session,"pathwaySelNode",choices=NULL,selected=NULL)
    updatePickerInput(session,"geneSel",choices=NULL,selected=NULL)
    updatePickerInput(session,"geneSelPathway",choices=NULL,selected=NULL)
    updatePickerInput(session,"sourceSel",choices=NULL,selected=NULL)
    updatePickerInput(session,"sourceSelPathway",choices=NULL,selected=NULL)
    updatePickerInput(session,"destSel",choices=NULL,selected=NULL)
    updatePickerInput(session,"destSelPathway",choices=NULL,selected=NULL)
    updatePickerInput(session,"networkSel",choices=names(rv$data.list),selected = initial.networks)
  }
  
  refresh.pickers <- function() {
    shinyjs::runjs("$('#pathwaySel, #pathwaySelNode, #geneSel, #geneSelPathway, #sourceSel, #sourceSelPathway, #destSel, #destSelPathway, #networkSel').selectpicker('refresh');")
  }
  
  all.nodes.choices <- reactive({
    if(is.empty(input$networkSel) || is.null(rv$list.all.nodes)) list()
    else filter.selectable.nodes(rv$list.all.nodes, input$networkSel)
  })
  
  all.pathways.choices <- reactive({
    if(is.null(rv$common.pathways)) NULL else {
      sorted <- sort(rv$common.pathways)
      setNames(sorted, truncate.label(sorted))
    }
  })
  
  nodes.for.selected.pathway <- reactive({
    if(!is.empty(last.sel.pathways())) {
      get.list.selectable.nodes(last.sel.pathways(), input$networkSel, rv$data.list, pathway.list)
    } else {
      list()
    }
  })
  
  pathways.for.selected.nodes <- reactive({
    if(input$outputType == "neighbors") {
      if(!is.empty(last.sel.genes())) sort(get.list.selectable.pathways(last.sel.genes(), rv$data.list, pathway.list)) else NULL
    } else {
      if(!is.empty(last.sel.source()) && !is.empty(last.sel.dest())) {
        pathways.with.source <- get.list.selectable.pathways(last.sel.source(), rv$data.list, pathway.list)
        pathways.with.dest   <- get.list.selectable.pathways(last.sel.dest(), rv$data.list, pathway.list)
        sort(intersect(pathways.with.source, pathways.with.dest))
      } else {
        NULL
      }
    }
  })
  
  observeEvent(all.pathways.choices(), {
    updatePickerInput(session, "pathwaySel", choices=all.pathways.choices(), selected=NULL)
    if(is.null(input$searchOpt) || input$searchOpt == "Pathway") last.sel.pathways(NULL)
    refresh.pickers()
  }, ignoreNULL = FALSE)
  
  observeEvent(pathways.for.selected.nodes(), {
    updatePickerInput(session, "pathwaySelNode", choices=pathways.for.selected.nodes(), selected=NULL)
    if(!is.null(input$searchOpt) && input$searchOpt == "Node") last.sel.pathways(NULL)
    refresh.pickers()
  }, ignoreNULL = FALSE)
  
  observeEvent(all.nodes.choices(), {
    choices <- all.nodes.choices()
    updatePickerInput(session, "geneSel", choices=choices, selected=NULL)
    updatePickerInput(session, "sourceSel", choices=choices, selected=NULL)
    updatePickerInput(session, "destSel", choices=choices, selected=NULL)
    if(is.null(input$searchOpt) || input$searchOpt == "Node") {
      last.sel.genes(NULL)
      last.sel.source(NULL)
      last.sel.dest(NULL)
    }
    refresh.pickers()
  }, ignoreNULL = FALSE)
  
  observeEvent(nodes.for.selected.pathway(), {
    choices <- nodes.for.selected.pathway()
    updatePickerInput(session, "geneSelPathway", choices=choices, selected=NULL)
    updatePickerInput(session, "sourceSelPathway", choices=choices, selected=NULL)
    updatePickerInput(session, "destSelPathway", choices=choices, selected=NULL)
    if(!is.null(input$searchOpt) && input$searchOpt == "Pathway") {
      last.sel.genes(NULL)
      last.sel.source(NULL)
      last.sel.dest(NULL)
    }
    refresh.pickers()
  }, ignoreNULL = FALSE)
  
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
        "$(this.api().table().header()).css({'background-color': 'black','color': 'white'});",
        "$(this.api().table().body()).css({'background-color': 'black','color': 'white'});",
        "}")
    ),rownames = F,escape=F)
  })
  
  observeEvent(input$delete_button, {
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
  
  plot.trigger <- reactive({
    list(pathways = last.sel.pathways(),
         genes = last.sel.genes(),
         source.genes = last.sel.source(),
         dest.genes = last.sel.dest(),
         networks = input$networkSel,
         hide = resolve.hide.elements(input$showElements),
         view.mode = input$outputType,
         max.hops = input$maxHopsEgo,
         max.length = if(isTRUE(input$limitPathLength)) input$maxPathLength else Inf,
         data.list = rv$data.list)
  })
  plot.trigger.d <- debounce(plot.trigger, 400)
  
  observeEvent(plot.trigger(), {
    shinyjs::disable("searchOpt")
    shinyjs::disable("outputType")
    shinyjs::disable("pathwaySel")
    shinyjs::disable("pathwaySelNode")
    shinyjs::disable("geneSel")
    shinyjs::disable("geneSelPathway")
    shinyjs::disable("sourceSel")
    shinyjs::disable("sourceSelPathway")
    shinyjs::disable("destSel")
    shinyjs::disable("destSelPathway")
    shinyjs::disable("networkSel")
    shinyjs::disable("showElements")
    shinyjs::disable("maxHopsEgo")
    shinyjs::disable("limitPathLength")
    shinyjs::disable("maxPathLength")
    refresh.pickers()
    shinyjs::show("updatingMsg")
  }, ignoreInit = TRUE)
  
  legend.info <- reactiveVal(NULL)
  
  output$plotPathway <- renderVisNetwork({
    trig <- plot.trigger.d()
    on.exit({
      shinyjs::enable("searchOpt")
      shinyjs::enable("outputType")
      shinyjs::enable("pathwaySel")
      shinyjs::enable("pathwaySelNode")
      shinyjs::enable("geneSel")
      shinyjs::enable("geneSelPathway")
      shinyjs::enable("sourceSel")
      shinyjs::enable("sourceSelPathway")
      shinyjs::enable("destSel")
      shinyjs::enable("destSelPathway")
      shinyjs::enable("networkSel")
      shinyjs::enable("showElements")
      shinyjs::enable("maxHopsEgo")
      shinyjs::enable("limitPathLength")
      shinyjs::enable("maxPathLength")
      refresh.pickers()
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
                                          max.hops = resolved.max.hops, max.length = resolved.max.length)
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
  
}
