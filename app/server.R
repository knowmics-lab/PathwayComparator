function(input, output, session) {
  
  options(shiny.maxRequestSize=200*1024^2)
  
  #last.choice <- reactive(input$geneSel)
  uploaded.files <- reactive(input$uploadedFiles)
  
  shinyjs::hide("searchOpt")
  shinyjs::hide("startSel")
  shinyjs::hide("networkSel")
  shinyjs::hide("hideElements")
  
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
  
  # Renderizza dinamicamente i picker in base alla selezione
  output$startSel <- renderUI({
    if (input$searchOpt == "Pathway") {
      # Prima pathway, poi gene
      tagList(
        pickerInput("pathwaySel","Pathway",choices=NULL,multiple = T,width="100%",options=pickerOptions(actionsBox=T,liveSearch = T)),
        pickerInput("geneSel","Node",choices=NULL,multiple = T,options=pickerOptions(actionsBox=T,liveSearch = T))
      )
    } else {
      # Prima gene, poi pathway
      tagList(
        pickerInput("geneSel","Node",choices=NULL,multiple = T,options=pickerOptions(actionsBox=T,liveSearch = T)),
        pickerInput("pathwaySel","Pathway",choices=NULL,multiple = T,width="100%",options=pickerOptions(actionsBox=T,liveSearch = T))
      )
    }
  })
  
  output$listFiles <- DT::renderDT({
   #print(uploaded.files())
    if(length(uploaded.files())==0) {
      tab.data <- data.frame("  "=character()," "=character(),check.names = F)
    } else {
      tab.data <- data.frame("  "=uploaded.files()," "=rep("Remove",length(uploaded.files())),check.names = F)
      tab.data[," "] <- shinyInput(actionButton, length(uploaded.files()),'delete_',label = "",icon=icon("trash"),
                              style = "color: white; background-color: black",
                              onclick = paste0('Shiny.onInputChange( \"delete_button\" , this.id, {priority: \"event\"})'))
    }
    datatable(tab.data, options = list(dom="t",ordering=F, language = list(
        zeroRecords = "No pathway files"),
                                         initComplete = JS(
                                           "function(settings, json) {",
                                           "$(this.api().table().header()).css({'background-color': 'black','color': 'white'});",
                                           "$(this.api().table().body()).css({'background-color': 'black','color': 'white'});",
                                           "}")
                                         ),rownames = F,escape=F)
  })
  
  observeEvent(input$delete_button, {
    selectedRow <- as.numeric(strsplit(input$delete_button, "_")[[1]][2])
    rem.organism <- data.list[[selectedRow]]$organism
    freq.organisms <- table(sapply(data.list,function(el){el$organism}))
    if(freq.organisms[rem.organism]==1) {
      metapathway.list[[rem.organism]] <<- NULL
      pathway.list[[rem.organism]] <<- NULL
    }
    data.list[[selectedRow]] <<- NULL
    #if(length(data.list)>0) {}
    updateSelectInput(session,"uploadedFiles",choices=input$uploadedFiles[-selectedRow],selected=input$uploadedFiles[-selectedRow])
    if(length(data.list)==1) {
      updateActionButton(session,"compare",label="Show",disabled=F)
    } else if(length(data.list)>1) {
      updateActionButton(session,"compare",label="Compare",disabled=F)
    } else {
      updateActionButton(session,"compare",label="Show",disabled=T)
    }
  })
  
  observeEvent(input$hiddenUpload,{
    file.list <- input$hiddenUpload
    req(file.list)
    #print(file.list)
    dataname.list <- c()
    for(i in 1:nrow(file.list)) {
      file <- file.list[i,]
      dataname <- strsplit(file$name,"\\.")[[1]][1]
      file.data <- read.phensim.file(file$datapath)
      if(dataname %in% names(data.list)) {
        dataname <- paste0(dataname,"_",(length(input$uploadedFiles)+1))
      }
      data.list[[dataname]] <<- file.data
      organism <- as.character(file.data$organism)
      if(!organism %in% names(metapathway.list)) {
        metapathway.list[[organism]] <<- readRDS(paste0("Data/Metapathways/",organism,".rds"))
        pathway.list[[organism]] <<- readRDS(paste0("Data/Pathways/",organism,".rds"))
      }
      dataname.list <- c(dataname.list,dataname)
    }
    updateSelectInput(session,"uploadedFiles",choices=c(input$uploadedFiles,dataname.list),selected=c(input$uploadedFiles,dataname.list))
    if(length(data.list)==1) {
      updateActionButton(session,"compare",label="Show",disabled=F)
    } else if(length(data.list)>1) {
      updateActionButton(session,"compare",label="Compare",disabled=F)
    } else {
      updateActionButton(session,"compare",label="Show",disabled=T)
    }
  })
  
  observeEvent(input$compare, {
    data.list <- data.list[order(sapply(data.list, function(el) el$organism))]
    list.organism <- as.character(sapply(data.list,function(el){el$organism}))
    #print(list.organism)
    if(length(list.organism)>1){
      org.pairs <- apply(t(combn(sort(list.organism),2)),1,function(row){paste0(row,collapse="-")})
      for(pair in org.pairs) {
        ortho.list[[pair]] <<- readRDS(paste0("Data/Orthologs/",pair,".rds"))
      } 
    } else {
      pair <- paste0(list.organism,"-",list.organism)
      ortho.list[[pair]] <<- readRDS(paste0("Data/Orthologs/",pair,".rds"))
    }
    
    #Update list of selectable pathways and nodes
    list.pathways <- sapply(list.organism,function(org){unique(pathway.list[[org]]$pathwayName)})
    if(is.list(list.pathways)) {
      common.pathways <<- Reduce(intersect,list.pathways)
    } else {
      common.pathways <<- list.pathways
    }
    list.all.nodes <<- get.list.selectable.nodes(sort(common.pathways),names(data.list),input$hideElements)
    
    #Make elements visible
    shinyjs::show("searchOpt")
    shinyjs::show("startSel")
    shinyjs::show("networkSel")
    shinyjs::show("hideElements")
    
    #Update picker lists
    if(input$searchOpt=="Pathway"){
      updatePickerInput(session,"pathwaySel",choices=sort(common.pathways),selected=NULL)
      updatePickerInput(session,"geneSel",choices=list(),selected=NULL)
      last.sel.genes(NULL)
    } else {
      updatePickerInput(session,"geneSel",choices=list.all.nodes,selected=NULL)
      updatePickerInput(session,"pathwaySel",choices=NULL,selected=NULL)
      last.sel.pathways(NULL)
    }
    max.opt <- 3
    if(length(data.list)==2)
      max.opt <- 2
    updatePickerInput(session,"networkSel",choices=names(data.list),selected = names(data.list)[1:max.opt])
  })
  
  observeEvent(input$searchOpt,{
    if(exists("common.pathways")){
      if(input$searchOpt=="Pathway"){
        updatePickerInput(session,"pathwaySel",choices=sort(common.pathways),selected=NULL)
        last.sel.pathways(NULL)
        updatePickerInput(session,"geneSel",choices=list(),selected=NULL)
        last.sel.genes(NULL)
      } else {
        list.selectable.nodes <- filter.selectable.nodes(list.all.nodes,input$networkSel,input$hideElements)
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
        updated.list.nodes <- get.list.selectable.nodes(last.sel.pathways(),input$networkSel,input$hideElements)
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
        updated.list.pathways <- get.list.selectable.pathways(last.sel.genes())
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
        updated.list <- get.list.selectable.nodes(last.sel.pathways(),input$networkSel,input$hideElements)
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
        updated.gene.list <- filter.selectable.nodes(list.all.nodes,input$networkSel,input$hideElements)
        #Update list of already selected genes (if needed)
        if(any(last.sel.genes() %in% updated.gene.list)) {
          updatePickerInput(session, "geneSel", selected=last.sel.genes()[last.sel.genes() %in% updated.gene.list], choices=updated.gene.list)
        } else{
          updatePickerInput(session, "geneSel", selected=NULL, choices=updated.gene.list)
          last.sel.genes(NULL)
        }
      } else {
        updatePickerInput(session,"geneSel", selected=NULL, choices=list())
        last.sel.genes(NULL)
      }
    }
  }, ignoreNULL = F)
  
  output$plotPathway <- renderVisNetwork({
    print(paste0("PATHWAYS: ",last.sel.pathways()))
    print(paste0("GENES: ",last.sel.genes()))
    if(!is.null(last.sel.pathways()) && !is.null(last.sel.genes()) && !is.null(input$networkSel)){
      print("prova")
      multilayer.net <- build.pathway.net(data.list,metapathway.list,pathway.list,ortho.list,
                            input$networkSel,last.sel.pathways(),last.sel.genes(),input$hideElements)
      pathway.plot <- plot.pathway(multilayer.net$nodes,multilayer.net$edges)
      pathway.plot
    }
  })
  
}