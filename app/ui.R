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
    width = "300px",
    dark_bg = "black"
  )
)

header <- dashboardHeader(tags$li(a(icon("github"), "Github", target = "_blank", href = "https://github.com/GMicale/PathwayComparator",
                                    style = "font-size: 18px;"),class= 'dropdown'),
                          title = span("PAthway COmparator", style = "color: white; font-size: 23px; font-family: 'Segoe Print', 'Bradley Hand', cursive;"), titleWidth = 350)

sidebar <- dashboardSidebar(
  sidebarMenu(
    hidden(selectInput("uploadedFiles",label=NULL,choices=NULL,multiple=T)),
    p(HTML("<b>1. Pathway files</b> &nbsp;"),
      span(icon("circle-question"), id="listPath",
           title="Supported organisms: Human, Mouse, Rat, Worm, Fly, Zebrafish, Arabidopsis"),
      class="section-label"),
    DTOutput("listFiles"),
    fileInput("hiddenUpload",label=NULL,multiple=T,buttonLabel = "Add file"),
    actionButton("compare","Show",disabled = T)
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
          tags$li("Click \"Show\"/\"Compare\"."),
          tags$li("Choose search mode (by pathway o by node), select perturbation score files to compare (max 3) and filter node types to hide for the visualization."),
          tags$li("Select one or more pathways and nodes for the visualization from the drop-down menus.")
        )
    )
  ),
  
  fluidRow(
    column(12,radioButtons("searchOpt","2. Search by:",choices=c("Pathway","Node"),inline=T,selected="Pathway"))
  ),
  fluidRow(
    column(12, tags$p("3. Compare & filter", id="filterStepLabel", class="step-label"))
  ),
  fluidRow(
    column(6,uiOutput("startSel")),
    column(2,
           pickerInput("networkSel",
                       label = tagList("Layer ", icon("circle-question",
                                                      title="You can compare up to 3 layers at a time.")),
                       choices=NULL,multiple = T,width="100%",options=list(`max-options`=3))
    ),
    column(4,
           tags$p("", style="margin-left:20px; margin-top:30px; margin-bottom:5px; font-weight:bold;"),
           checkboxGroupInput("hideElements", NULL, c("Hide chemical entities","Hide drugs","Hide miRNAs"),
                              selected = c("Hide chemical entities","Hide drugs","Hide miRNAs"), inline=T)
    )
  ),
  fluidRow(
    column(12,
           div(class = "view-mode-row",
               radioButtons("viewMode", NULL,
                            choices = c("Show selected nodes in selected pathways" = "neighbors",
                                        "Show all paths to selected nodes" = "paths"),
                            selected = "neighbors", inline = TRUE),
               conditionalPanel(
                 condition = "input.viewMode == 'paths'",
                 numericInput("maxHops", "Max hops", value = 2, min = 1, max = 50, step = 1, width = "140px")
               )
           )
    )
  ),
  fluidRow(
    column(12, tags$p(icon("spinner", class="fa-spin"), " Updating visualization...",
                      id="updatingMsg", class="updating-msg"))
  ),
  fluidRow(
    column(12,visNetworkOutput("plotPathway",height="80vh") %>% withSpinner(color = "#8a7fd1"))
  ),
  
  tags$head(tags$style(sass(sass_file("www/bootswatch-cyborg.scss")))),
  
  tags$head(tags$style(HTML("
    .navbar-custom-menu { position: absolute; display: inline-block; }
    .content-wrapper { background-color: white; }
    .box {
      -webkit-box-shadow: none; -moz-box-shadow: none; box-shadow: none;
      font-size: 18px;
      border: 1px solid rgba(0,0,0,0.08);
    }
    .main-sidebar { font-size: 18px; }
    .sidebar { background-color: black; }
    .control-label { font-size: 18px; }
    .form-group { font-size: 18px; }
    .section-label { margin-left: 20px; margin-top: 20px; font-size: 18px; font-weight: bold; }
    .step-label { margin-left: 0; margin-top: 5px; margin-bottom: 10px; font-size: 18px; font-weight: bold; }
    .updating-msg { margin-left: 0; margin-top: 5px; margin-bottom: 5px; font-size: 16px; color: #8a7fd1; }
    input[type=checkbox] { transform: scale(1.4); }
    .view-mode-row { display: flex; align-items: center; gap: 25px; flex-wrap: wrap; }
    .view-mode-row .form-group { margin-bottom: 0; }
    .selectize-input { font-size: 18px; line-height: 20px; }
    .selectize-dropdown { font-size: 16px; line-height: 18px; }
    .dropdown-menu ul li:nth-child(n) a { color: black !important; font-size: 18px; }
    .picker { font-size: 16px; color: black; }
    .btn-file { font-size: 14px; white-space: nowrap; overflow: hidden; text-overflow: ellipsis; }
    .input-group, .form-group .input-group { max-width: 100%; }
    .btn-delete-row { color: #fff; background-color: #2a2a2a; border: 1px solid #444; }
    .btn-delete-row:hover { background-color: #c0392b; border-color: #c0392b; color: #fff; }
  ")))
)

dashboardPage(header,sidebar,body,title = "PathwayComparator",skin = "purple")
