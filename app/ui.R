library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(shinyjs)
library(data.table)
library(DT)
library(visNetwork)
library(igraph)
library(shinycssloaders)
library(sass)
library(fresh)

source("PathwaysFunctions.R")

app.theme <- create_theme(
  adminlte_sidebar(
    width = "250px",
    dark_bg = "black"
  )
)

header <- dashboardHeader(tags$li(a(icon("github"), "Github", target = "_blank", href = "https://github.com/GMicale/PathwayComparator",
                                    style = "font-size: 14px;"),class= 'dropdown'),
                          title = span("PAthway COmparator", style = "color: white; font-size: 17px; font-family: 'Segoe Print', 'Bradley Hand', cursive;"), titleWidth = 250)

sidebar <- dashboardSidebar(
  sidebarMenu(
    hidden(selectInput("uploadedFiles",label=NULL,choices=NULL,multiple=T)),
    p(HTML("<b>Pathway files</b> &nbsp;"),
      span(icon("circle-question"), id="listPath",
           title="Supported organisms: Human, Mouse, Rat, Worm, Fly, Zebrafish, Arabidopsis"),
      class="section-label"),
    DTOutput("listFiles"),
    fileInput("hiddenUpload",label=NULL,multiple=T,buttonLabel = "Add file")
  )
)

body <- dashboardBody(
  use_theme(app.theme),
  useShinyjs(),
  
  conditionalPanel(
    condition = "!output.hasFiles",
    box(width = 12, status = "primary", solidHeader = FALSE, title = "How to use?",
        tags$ol(
          tags$li("Upload one or more node perturbation score files (PHENSIM, MITHrIL o custom) from left panel."),
          tags$li("Choose search mode (by pathway o by node) and output type (Ego-network o All paths)."),
          tags$li("Select pathways/nodes progressively from the drop-down menus, then adjust layers and node types to show for the visualization.")
        )
    )
  ),
  
  conditionalPanel(
    condition = "output.compared",
    fluidRow(
      column(12,
             div(class = "search-mode-row",
                 radioButtons("searchOpt","Search by",choices=c("Pathway","Node"),inline=T,selected="Pathway"),
                 radioButtons("outputType","Output type",
                              choices = c("Ego-network" = "neighbors", "All paths" = "paths"),
                              selected = "neighbors", inline = TRUE),
                 conditionalPanel(
                   condition = "input.outputType=='neighbors'",
                   div(class = "inline-numeric-input",
                       numericInput("maxHopsEgo", "Max hops", value = 1, min = 1, max = 20, step = 1))
                 ),
                 conditionalPanel(
                   condition = "input.outputType=='paths'",
                   div(class = "limit-path-length-row",
                       checkboxInput("limitPathLength", "Limit max path length to", value = FALSE),
                       conditionalPanel(
                         condition = "input.limitPathLength",
                         numericInput("maxPathLength", NULL, value = 3, min = 1, max = 50, step = 1)
                       )
                   )
                 )
             )
      )
    ),
    fluidRow(
      column(12,
             div(class = "search-dropdowns-row",
                 conditionalPanel(
                   condition = "input.searchOpt=='Pathway'",
                   div(class = "search-dropdown-item",
                       pickerInput("pathwaySel","Pathway",choices=NULL,multiple=T,width="auto",options=pickerOptions(actionsBox=T,liveSearch=T)))
                 ),
                 conditionalPanel(
                   condition = "input.searchOpt=='Pathway' && input.outputType=='neighbors' && input.pathwaySel && input.pathwaySel.length>0",
                   div(class = "search-dropdown-item",
                       pickerInput("geneSelPathway","Node",choices=NULL,multiple=T,width="auto",options=pickerOptions(actionsBox=T,liveSearch=T)))
                 ),
                 conditionalPanel(
                   condition = "input.searchOpt=='Pathway' && input.outputType=='paths' && input.pathwaySel && input.pathwaySel.length>0",
                   div(class = "search-dropdown-item",
                       pickerInput("sourceSelPathway","Source node",choices=NULL,multiple=T,width="auto",options=pickerOptions(actionsBox=T,liveSearch=T)))
                 ),
                 conditionalPanel(
                   condition = "input.searchOpt=='Pathway' && input.outputType=='paths' && input.pathwaySel && input.pathwaySel.length>0",
                   div(class = "search-dropdown-item",
                       pickerInput("destSelPathway","Destination node",choices=NULL,multiple=T,width="auto",options=pickerOptions(actionsBox=T,liveSearch=T)))
                 ),
                 conditionalPanel(
                   condition = "input.searchOpt=='Node' && input.outputType=='neighbors'",
                   div(class = "search-dropdown-item",
                       pickerInput("geneSel","Node",choices=NULL,multiple=T,width="auto",options=pickerOptions(actionsBox=T,liveSearch=T)))
                 ),
                 conditionalPanel(
                   condition = "input.searchOpt=='Node' && input.outputType=='paths'",
                   div(class = "search-dropdown-item",
                       pickerInput("sourceSel","Source node",choices=NULL,multiple=T,width="auto",options=pickerOptions(actionsBox=T,liveSearch=T)))
                 ),
                 conditionalPanel(
                   condition = "input.searchOpt=='Node' && input.outputType=='paths'",
                   div(class = "search-dropdown-item",
                       pickerInput("destSel","Destination node",choices=NULL,multiple=T,width="auto",options=pickerOptions(actionsBox=T,liveSearch=T)))
                 ),
                 conditionalPanel(
                   condition = "input.searchOpt=='Node' &&
                         ((input.outputType=='neighbors' && input.geneSel && input.geneSel.length>0) ||
                          (input.outputType=='paths' && input.sourceSel && input.sourceSel.length>0 && input.destSel && input.destSel.length>0))",
                   div(class = "search-dropdown-item",
                       pickerInput("pathwaySelNode","Pathway",choices=NULL,multiple=T,width="auto",options=pickerOptions(actionsBox=T,liveSearch=T)))
                 )
             )
      )
    )
  ),
  conditionalPanel(
    condition = "output.readyToPlot",
    fluidRow(
      column(2,
             pickerInput("networkSel",
                         label = tagList("Layer ", icon("circle-question",
                                                        title="You can compare up to 3 layers at a time.")),
                         choices=NULL,multiple = T,width="100%",options=list(`max-options`=3))
      ),
      column(4,
             tags$p("", style="margin-left:20px; margin-top:26px; margin-bottom:4px; font-weight:bold;"),
             checkboxGroupInput("showElements", NULL, c("Show compounds","Show drugs","Show miRNAs"),
                                selected = character(0), inline=T)
      )
    )
  ),
  fluidRow(
    column(12, tags$p(icon("spinner", class="fa-spin"), " Updating visualization...",
                      id="updatingMsg", class="updating-msg"))
  ),
  fluidRow(
    column(10,visNetworkOutput("plotPathway",height="80vh") %>% withSpinner(color = "#8a7fd1")),
    column(2, tags$div(style = "margin-top:20px;", uiOutput("graphLegend")))
  ),
  
  tags$head(tags$style(sass(sass_file("www/bootswatch-cyborg.scss")))),
  
  tags$head(tags$style(HTML("
    .navbar-custom-menu { position: absolute; display: inline-block; }
    .content-wrapper { background-color: white; }
    .box {
      -webkit-box-shadow: none; -moz-box-shadow: none; box-shadow: none;
      font-size: 14px;
      border: 1px solid rgba(0,0,0,0.08);
    }
    .main-sidebar { font-size: 14px; }
    .sidebar { background-color: black; }
    .control-label { font-size: 14px; }
    .form-group { font-size: 14px; margin-bottom: 8px; }
    .section-label { margin-left: 20px; margin-top: 12px; margin-bottom: 4px; font-size: 15px; font-weight: bold; }
    .step-label { margin-left: 0; margin-top: 4px; margin-bottom: 6px; font-size: 15px; font-weight: bold; }
    .search-mode-row { display: flex; align-items: center; gap: 40px; flex-wrap: wrap; margin-bottom: 4px; }
    .search-mode-row #searchOpt, .search-mode-row #outputType { margin-bottom: 0; }
    #searchOpt, #outputType { display: flex; align-items: center; gap: 15px; flex-wrap: wrap; margin-bottom: 4px; }
    .inline-numeric-input .form-group { display: flex; align-items: center; gap: 8px; margin-bottom: 0; }
    .inline-numeric-input .form-group label { margin-bottom: 0; white-space: nowrap; }
    .inline-numeric-input .form-group input { width: 70px; }
    .limit-path-length-row::after { content: ''; display: table; clear: both; }
    .limit-path-length-row > div { float: left; margin: 0 !important; }
    .limit-path-length-row .form-group { margin: 0 !important; width: auto !important; }
    .limit-path-length-row input[type='number'] { width: 60px !important; margin-left: 6px; margin-top: -6px; }
    .shiny-label-null { display: none; }
    #searchOpt > label, #outputType > label { font-size: 15px; font-weight: bold; margin-bottom: 0; }
    .updating-msg { margin-left: 0; margin-top: 4px; margin-bottom: 4px; font-size: 13px; color: #8a7fd1; }
    input[type=checkbox] { transform: scale(1.15); }
    .search-dropdowns-row { display: flex; gap: 20px; flex-wrap: wrap; align-items: flex-start; margin-bottom: 8px; }
    .search-dropdown-item { flex: 1 1 240px; min-width: 200px; max-width: 340px; }
    .search-dropdown-item .bootstrap-select { max-width: 100%; }
    .search-dropdown-item .dropdown-toggle {
      max-width: 100%; overflow: hidden; text-overflow: ellipsis; white-space: nowrap;
    }
    .search-dropdown-item .filter-option,
    .search-dropdown-item .filter-option-inner-inner {
      max-width: 100%; overflow: hidden; text-overflow: ellipsis; white-space: nowrap;
    }
    .selectize-input { font-size: 14px; line-height: 16px; }
    .selectize-dropdown { font-size: 13px; line-height: 15px; }
    .dropdown-menu ul li:nth-child(n) a { color: black !important; font-size: 13px; }
    .picker { font-size: 13px; color: black; }
    .btn-file { font-size: 13px; white-space: nowrap; overflow: hidden; text-overflow: ellipsis; }
    .input-group, .form-group .input-group { max-width: 100%; }
    #listFiles thead { display: none; }
    #listFiles table.dataTable td { padding: 3px 8px; font-size: 13px; vertical-align: middle; }
    .btn-delete-row { color: #fff; background-color: #2a2a2a; border: 1px solid #444; padding: 1px 8px; font-size: 12px; }
    .btn-delete-row:hover { background-color: #c0392b; border-color: #c0392b; color: #fff; }
  ")))
)

dashboardPage(header,sidebar,body,title = "PathwayComparator",skin = "purple")
