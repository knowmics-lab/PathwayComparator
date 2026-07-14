function(input, output, session) {
  
  options(shiny.maxRequestSize=200*1024^2)
  
  rv <- reactiveValues(data.list = list(), common.pathways = NULL, list.all.nodes = NULL)
  output$hasFiles <- reactive({
    length(rv$data.list) > 0
  })
  outputOptions(output, "hasFiles", suspendWhenHidden = FALSE)
  
  uploaded.files <- reactive(input$uploadedFiles)
  
  shinyjs::hide("searchOpt")
  shinyjs::hide("startSel")
  shinyjs::hide("networkSel")
  shinyjs::hide("hideElements")
  shinyjs::hide("viewMode")
  shinyjs::hide("updatingMsg")
  
  last.sel.pathways <- reactiveVal(NULL)
  last.sel.genes <- reactiveVal(NULL)
  
  shinyInput <- function(FUN, len, id, ...) {
    inputs <- character(len)
    for (i in seq_len(len)) {
      inputs[i] <- as.character(FUN(paste0(id, i), ...))
    }
    inputs
  }
  
  observe({
    last.sel.pathways(input$pathwaySel)
  })
  
  observe({
    last.sel.genes(input$geneSel)
  })
  
  output$startSel <- renderUI({
    if (input$searchOpt == "Pathway") {
      # Prima pathway, poi gene
      fluidRow(
        column(6, pickerInput("pathwaySel","Pathway",choices=NULL,multiple = T,width="100%",options=pickerOptions(actionsBox=T,liveSearch = T))),
        column(6, pickerInput("geneSel","Node",choices=NULL,multiple = T,width="100%",options=pickerOptions(actionsBox=T,liveSearch = T)))
      )
    } else {
      # Prima gene, poi pathway
      fluidRow(
        column(6, pickerInput("geneSel","Node",choices=NULL,multiple = T,width="100%",options=pickerOptions(actionsBox=T,liveSearch = T))),
        column(6, pickerInput("pathwaySel","Pathway",choices=NULL,multiple = T,width="100%",options=pickerOptions(actionsBox=T,liveSearch = T)))
      )
    }
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
    
    rv$common.pathways <- NULL
    rv$list.all.nodes <- NULL
    last.sel.pathways(NULL)
    last.sel.genes(NULL)
    updatePickerInput(session,"pathwaySel",choices=NULL,selected=NULL)
    updatePickerInput(session,"geneSel",choices=NULL,selected=NULL)
    updatePickerInput(session,"networkSel",choices=NULL,selected=NULL)
    shinyjs::hide("searchOpt")
    shinyjs::hide("startSel")
    shinyjs::hide("networkSel")
    shinyjs::hide("hideElements")
    shinyjs::hide("viewMode")
    shinyjs::hide("updatingMsg")
    
    if(length(rv$data.list)==1) {
      updateActionButton(session,"compare",label="Show",disabled=F)
    } else if(length(rv$data.list)>1) {
      updateActionButton(session,"compare",label="Compare",disabled=F)
    } else {
      updateActionButton(session,"compare",label="Show",disabled=T)
    }
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
    if(length(rv$data.list)==1) {
      updateActionButton(session,"compare",label="Show",disabled=F)
    } else if(length(rv$data.list)>1) {
      updateActionButton(session,"compare",label="Compare",disabled=F)
    } else {
      updateActionButton(session,"compare",label="Show",disabled=T)
    }
  })
  
  observeEvent(input$compare, {
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
    rv$list.all.nodes <- get.list.selectable.nodes(sort(rv$common.pathways),names(rv$data.list),character(0),rv$data.list,pathway.list)
    
    #Make elements visible
    shinyjs::show("searchOpt")
    shinyjs::show("startSel")
    shinyjs::show("networkSel")
    shinyjs::show("hideElements")
    shinyjs::show("viewMode")
    
    max.opt <- 3
    if(length(rv$data.list)==2)
      max.opt <- 2
    initial.networks <- names(rv$data.list)[1:max.opt]
    
    #Update picker lists
    if(input$searchOpt=="Pathway"){
      updatePickerInput(session,"pathwaySel",choices=sort(rv$common.pathways),selected=NULL)
      updatePickerInput(session,"geneSel",choices=list(),selected=NULL)
      last.sel.genes(NULL)
    } else {
      initial.gene.choices <- filter.selectable.nodes(rv$list.all.nodes, initial.networks, input$hideElements)
      updatePickerInput(session,"geneSel",choices=initial.gene.choices,selected=NULL)
      updatePickerInput(session,"pathwaySel",choices=NULL,selected=NULL)
      last.sel.pathways(NULL)
    }
    updatePickerInput(session,"networkSel",choices=names(rv$data.list),selected = initial.networks)
  })
  
  observeEvent(input$searchOpt,{
    if(!is.null(rv$common.pathways)){
      if(input$searchOpt=="Pathway"){
        updatePickerInput(session,"pathwaySel",choices=sort(rv$common.pathways),selected=NULL)
        last.sel.pathways(NULL)
        updatePickerInput(session,"geneSel",choices=list(),selected=NULL)
        last.sel.genes(NULL)
      } else {
        list.selectable.nodes <- filter.selectable.nodes(rv$list.all.nodes,input$networkSel,input$hideElements)
        updatePickerInput(session,"geneSel",choices=list.selectable.nodes,selected=NULL)
        last.sel.genes(NULL)
        updatePickerInput(session,"pathwaySel",choices=NULL,selected=NULL)
        last.sel.pathways(NULL)
      }
    }
  })
  
  observeEvent(last.sel.pathways(),{
    if (input$searchOpt == "Pathway"){
      if(!is.null(last.sel.pathways()) && !is.null(input$networkSel)) {
        updated.list.nodes <- get.list.selectable.nodes(last.sel.pathways(),input$networkSel,input$hideElements,rv$data.list,pathway.list)
        updatePickerInput(session,"geneSel",choices=updated.list.nodes)
      } else {
        updatePickerInput(session,"geneSel",choices=list(),selected=NULL)
        last.sel.genes(NULL)
      }
    }
  }, ignoreNULL = FALSE)
  
  observeEvent(last.sel.genes(),{
    if(input$searchOpt == "Node"){
      if(!is.null(last.sel.genes()) && !is.null(input$networkSel)) {
        updated.list.pathways <- get.list.selectable.pathways(last.sel.genes(),rv$data.list,pathway.list)
        updatePickerInput(session,"pathwaySel",choices=sort(updated.list.pathways),selected=NULL)
      } else {
        updatePickerInput(session,"pathwaySel",choices=NULL,selected=NULL)
        last.sel.pathways(NULL)
      }
    }
  }, ignoreNULL = FALSE)
  
  observeEvent({
    input$networkSel
    input$hideElements
  },{
    if(input$searchOpt=="Pathway"){
      if(!is.null(last.sel.pathways()) && !is.null(input$networkSel)) {
        updated.list <- get.list.selectable.nodes(last.sel.pathways(),input$networkSel,input$hideElements,rv$data.list,pathway.list)
        updated.gene.list <- unname(unlist(updated.list))
        if(any(last.sel.genes() %in% updated.gene.list)) {
          updatePickerInput(session, "geneSel", selected=last.sel.genes()[last.sel.genes() %in% updated.gene.list], choices=updated.list)
        } else{
          updatePickerInput(session, "geneSel", selected=NULL, choices=updated.list)
          last.sel.genes(NULL)
        }
      } else {
        updatePickerInput(session,"geneSel", selected=NULL, choices=list())
        last.sel.genes(NULL)
      }
    } else {
      if(!is.null(last.sel.genes()) && !is.null(input$networkSel)) {
        #Filter list of selectable genes
        updated.gene.list <- filter.selectable.nodes(rv$list.all.nodes,input$networkSel,input$hideElements)
        updated.gene.flat <- unname(unlist(updated.gene.list))
        #Update list of already selected genes (if needed)
        if(any(last.sel.genes() %in% updated.gene.flat)) {
          updatePickerInput(session, "geneSel", selected=last.sel.genes()[last.sel.genes() %in% updated.gene.flat], choices=updated.gene.list)
          updated.list.pathways <- get.list.selectable.pathways(last.sel.genes(),rv$data.list,pathway.list)
          pathway.selected <- if(!is.null(last.sel.pathways()) && any(last.sel.pathways() %in% updated.list.pathways)) {
            last.sel.pathways()[last.sel.pathways() %in% updated.list.pathways]
          } else {
            NULL
          }
          updatePickerInput(session,"pathwaySel",choices=sort(updated.list.pathways),selected=pathway.selected)
        } else{
          updatePickerInput(session, "geneSel", selected=NULL, choices=updated.gene.list)
          last.sel.genes(NULL)
          updatePickerInput(session,"pathwaySel",choices=NULL,selected=NULL)
          last.sel.pathways(NULL)
        }
      } else {
        updatePickerInput(session,"geneSel", selected=NULL, choices=list())
        last.sel.genes(NULL)
        updatePickerInput(session,"pathwaySel",choices=NULL,selected=NULL)
        last.sel.pathways(NULL)
      }
    }
  }, ignoreNULL = F)
  
  plot.trigger <- reactive({
    list(pathways = last.sel.pathways(), genes = last.sel.genes(),
         networks = input$networkSel, hide = input$hideElements,
         view.mode = input$viewMode, max.hops = input$maxHops,
         data.list = rv$data.list)
  })
  plot.trigger.d <- debounce(plot.trigger, 400)
  
  refresh.pickers <- function() {
    shinyjs::runjs("$('#pathwaySel, #geneSel, #networkSel').selectpicker('refresh');")
  }
  
  observeEvent(plot.trigger(), {
    shinyjs::disable("searchOpt")
    shinyjs::disable("pathwaySel")
    shinyjs::disable("geneSel")
    shinyjs::disable("networkSel")
    shinyjs::disable("hideElements")
    shinyjs::disable("viewMode")
    shinyjs::disable("maxHops")
    refresh.pickers()
    shinyjs::show("updatingMsg")
  }, ignoreInit = TRUE)
  
  legend.info <- reactiveVal(NULL)
  
  output$plotPathway <- renderVisNetwork({
    trig <- plot.trigger.d()
    on.exit({
      shinyjs::enable("searchOpt")
      shinyjs::enable("pathwaySel")
      shinyjs::enable("geneSel")
      shinyjs::enable("networkSel")
      shinyjs::enable("hideElements")
      shinyjs::enable("viewMode")
      shinyjs::enable("maxHops")
      refresh.pickers()
      shinyjs::hide("updatingMsg")
    })
    if(!is.null(trig$pathways) && !is.null(trig$genes) && !is.null(trig$networks)){
      resolved.view.mode <- if(is.null(trig$view.mode)) "neighbors" else trig$view.mode
      resolved.max.hops <- if(resolved.view.mode!="paths" || is.null(trig$max.hops) || is.na(trig$max.hops)) Inf else trig$max.hops
      multilayer.net <- build.pathway.net(trig$data.list,metapathway.list,pathway.list,ortho.list,
                                          trig$networks,trig$pathways,trig$genes,trig$hide,
                                          view.mode = resolved.view.mode, max.hops = resolved.max.hops)
      pathway.plot <- plot.pathway(multilayer.net$nodes,multilayer.net$edges,view.mode = resolved.view.mode)
      legend.info(if(nrow(multilayer.net$nodes)>0) get.legend.info(multilayer.net$nodes) else NULL)
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
