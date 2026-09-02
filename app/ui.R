library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(shinyjs)
library(data.table)
library(DT)
library(visNetwork)
library(igraph)
library(shinycssloaders)
library(fresh)

source("PathwaysFunctions.R")

app.theme <- create_theme(
  adminlte_sidebar(
    width = "250px",
    dark_bg = "#24242b"
  )
)

header <- dashboardHeader(tags$li(a(icon("github"), "Github", target = "_blank", href = "https://github.com/GMicale/PathwayComparator",
                                    style = "font-size: 14px;"),class= 'dropdown'),
                          title = span(HTML("P<span style='opacity:0.85; font-weight:400;'>athway</span> C<span style='opacity:0.85; font-weight:400;'>omparator</span>"),
                                       style = "color: white; font-size: 17px; font-family: Georgia, 'Times New Roman', serif; font-weight: 700; letter-spacing: 0.01em;"), titleWidth = 250)

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
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "wizard.css"),
    tags$script(src = "wizard.js"),
    tags$script(src = "pathway_filter.js")
  ),
  
  conditionalPanel(
    condition = "!output.hasFiles",
    box(width = 12, status = "primary", solidHeader = FALSE, title = "How to use?",
        tags$ol(
          tags$li("Upload one or more node perturbation score files (PHENSIM, MITHrIL o custom) from left panel."),
          tags$li("Choose search mode (by pathway o by node) and output type (Ego-network o All paths)."),
          tags$li("Select pathways/nodes progressively, then adjust layers and node types to show for the visualization.")
        )
    )
  ),
  
  # ===================== TAB: Network view / Node comparison =====================
  # Le due modalita' di visualizzazione (grafico della rete vs tabella di
  # confronto punteggi tra i file caricati) condividono lo stesso
  # caricamento file, ma non la stessa selezione: il wizard di ricerca
  # (Search by Pathway/Node) resta specifico della scheda "Network view",
  # dato che "Node comparison" mostra SEMPRE tutti i nodi di tutti i file
  # caricati, senza bisogno di restringere a una pathway.
  #
  # Il pannello a tab compare non appena rv$data.list riceve il primo file
  # (output.hasFiles) - questo non soffre dello stesso problema di
  # tempistica dei pannelli di caricamento sotto, dato che deve solo
  # diventare vero UNA VOLTA e restare tale, non passare da vero a falso
  # all'interno dello stesso ciclo reattivo.
  #
  # Contenuto reale e indicatore di attesa in ciascuna scheda sono
  # avvolti in shinyjs::hidden() con ID propri, mostrati/nascosti
  # ESPLICITAMENTE lato server (vedi start.loading()/stop.loading() in
  # server.R) - un conditionalPanel basato su un reactiveValues NON
  # funzionerebbe qui, dato che l'intero caricamento/eliminazione avviene
  # in un solo observeEvent sincrono: un valore che passa da TRUE a FALSE
  # nello stesso ciclo reattivo non produce mai un aggiornamento
  # osservabile lato client.
  conditionalPanel(
    condition = "output.hasFiles",
    tabsetPanel(
      id = "mainTabs", type = "tabs",
      
      tabPanel("Node comparison",
               shinyjs::hidden(
                 tags$div(id = "nodeComparisonLoading", class = "tab-loading-panel",
                          icon("spinner", class="fa-spin"), " Loading files, please wait...")
               ),
               shinyjs::hidden(
                 tags$div(id = "nodeComparisonContent",
                          # Tabella nodo x network con delta e colore in scala assoluta
                          # (identica a quella usata nel grafico della rete, con
                          # trasformazione a radice quadrata per dare piu' risalto ai
                          # valori comuni piu' piccoli rispetto ai rari outlier estremi) -
                          # vedi build.node.comparison.table()/compute.abs.max() in
                          # PathwaysFunctions.R e output$nodeComparisonTable in server.R.
                          # Il filtro pathway (facoltativo, selezione multipla) restringe
                          # la tabella ai nodi che appartengono ad ALMENO UNA delle pathway
                          # scelte (unione, non intersezione); il tooltip sul nome del
                          # nodo (colonna "Node") continua pero' a mostrare SEMPRE tutte le
                          # pathway di appartenenza, non solo quella usata per filtrare.
                          fluidRow(
                            column(6, style="margin-top:16px;",
                                   tags$div(style="font-weight:bold; margin-bottom:6px;",
                                            "Filter by pathway ", icon("circle-question",
                                                                       title="Shows only nodes that belong to at least one of the selected pathways. Leave empty to show all nodes.")),
                                   tags$div(id = "pathwayFilterWidget")
                            )
                          ),
                          fluidRow(
                            column(12, tags$div(style="margin-top:8px;", DTOutput("nodeComparisonTable")))
                          )
                 )
               )
      ),
      
      tabPanel("Network view",
               shinyjs::hidden(
                 tags$div(id = "networkViewLoading", class = "tab-loading-panel",
                          icon("spinner", class="fa-spin"), " Loading files, please wait...")
               ),
               shinyjs::hidden(
                 tags$div(id = "networkViewContent",
                          fluidRow(
                            column(12, tags$div(id = "searchWizard"))
                          ),
                          # Layer da confrontare ed elementi della pathway da mostrare (di
                          # default nascosti - vedi selected=character(0) sotto): non
                          # influenzano cosa e' ricercabile/selezionabile nel wizard sopra,
                          # solo cosa viene disegnato nel grafico finale. Compaiono solo
                          # quando la selezione e' completa a sufficienza da produrre un
                          # grafico (output.readyToPlot).
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
                            column(10,visNetworkOutput("plotPathway",height="80vh") %>% withSpinner(color = "#6f63b8")),
                            column(2, tags$div(style = "margin-top:20px;", uiOutput("graphLegend")))
                          )
                 )
               )
      )
    )
  ),
  
  tags$head(tags$style(HTML("
    .navbar-custom-menu { position: absolute; display: inline-block; }
    .content-wrapper { background-color: #fafafc !important; }
    .skin-purple .main-header .logo, .skin-purple .main-header .navbar {
      background: linear-gradient(135deg, #8b7fd6, #453d7c) !important;
    }
    .skin-purple .main-header .logo:hover {
      background: linear-gradient(135deg, #9b8fe0, #4f4589) !important;
    }
    .box {
      -webkit-box-shadow: 0 1px 4px rgba(0,0,0,0.08); -moz-box-shadow: 0 1px 4px rgba(0,0,0,0.08); box-shadow: 0 1px 4px rgba(0,0,0,0.08);
      font-size: 14px;
      border: 1px solid rgba(0,0,0,0.06);
      border-radius: 6px;
    }
    .main-sidebar { font-size: 14px; }
    .sidebar { background-color: #24242b !important; }
    .control-label { font-size: 14px; }
    .form-group { font-size: 14px; margin-bottom: 8px; }
    .section-label { margin-left: 20px; margin-top: 12px; margin-bottom: 4px; font-size: 15px; font-weight: bold; }
    .updating-msg { margin-left: 0; margin-top: 4px; margin-bottom: 4px; font-size: 13px; color: #6f63b8; }
    .tab-loading-panel {
      padding: 60px 20px; text-align: center; font-size: 15px; color: #6f63b8; font-weight: 600;
    }
    input[type=checkbox] { transform: scale(1.15); }
    .selectize-input { font-size: 14px; line-height: 16px; }
    .selectize-dropdown { font-size: 13px; line-height: 15px; }
    .dropdown-menu ul li:nth-child(n) a { color: black !important; font-size: 13px; }
    .picker { font-size: 13px; color: black; }
    .btn-file { font-size: 13px; white-space: nowrap; overflow: hidden; text-overflow: ellipsis; }
    .input-group, .form-group .input-group { max-width: 100%; }
    #listFiles thead { display: none; }
    #listFiles table.dataTable td { padding: 3px 8px; font-size: 13px; vertical-align: middle; }
    .btn-delete-row { color: #fff; background-color: transparent; border: 1px solid #4a4a56; padding: 1px 8px; font-size: 12px; border-radius: 4px; transition: background-color 0.15s, border-color 0.15s; }
    .btn-delete-row:hover { background-color: #c0392b; border-color: #c0392b; color: #fff; }
    #mainTabs.nav-tabs {
      background: #f2f0fa; border: none; border-radius: 8px 8px 0 0;
      padding: 6px 6px 0; margin-bottom: 0; display: flex; gap: 4px;
    }
    #mainTabs.nav-tabs > li { float: none; }
    #mainTabs.nav-tabs > li > a {
      color: #6f63b8; font-weight: 600; font-size: 13.5px;
      border: none; border-radius: 6px 6px 0 0; background: transparent;
      padding: 11px 20px; margin: 0; transition: background-color 0.15s ease, color 0.15s ease;
    }
    #mainTabs.nav-tabs > li > a:hover {
      color: #453d7c; background: rgba(255,255,255,0.5); border-color: transparent;
    }
    #mainTabs.nav-tabs > li.active > a, #mainTabs.nav-tabs > li.active > a:hover, #mainTabs.nav-tabs > li.active > a:focus {
      color: #453d7c; background: #fff; border: none;
      box-shadow: 0 -1px 4px rgba(0,0,0,0.08);
    }
    .tab-content {
      background: #fff; border: 1px solid #e0dff0; border-top: none;
      border-radius: 0 0 8px 8px; padding: 18px 16px;
    }
    #nodeComparisonTable td { text-align: center; vertical-align: middle; }
    #nodeComparisonTable td:first-child { text-align: left; font-weight: 600; }
    #nodeComparisonTable td:last-child { font-size: 15px; font-weight: 600; }
    .node-cmp-tooltip {
      display: none; position: absolute; bottom: calc(100% + 6px); left: 50%; transform: translateX(-50%);
      background: #2c2c2c; color: #fff; font-size: 16px; padding: 6px 12px; border-radius: 5px;
      white-space: nowrap; z-index: 50; font-variant-numeric: tabular-nums;
    }
    .node-cmp-tooltip::after {
      content: ''; position: absolute; top: 100%; left: 50%; transform: translateX(-50%);
      border: 6px solid transparent; border-top-color: #2c2c2c;
    }
    .node-cmp-score-cell:hover .node-cmp-tooltip { display: block; }
    .node-cmp-name-cell { cursor: pointer; }
    .node-cmp-name-cell:hover { text-decoration: underline; text-decoration-color: #cfcbe8; }
    .node-cmp-pathway-badge {
      display: inline-block; background: #ece8f8; color: #6f63b8; font-weight: 700;
      font-size: 11px; border-radius: 9px; padding: 1px 7px; margin-left: 4px; vertical-align: middle;
    }
    .node-cmp-pathway-popover {
      display: none; position: fixed; z-index: 200;
      background: #fff; border-radius: 6px; box-shadow: 0 4px 16px rgba(0,0,0,0.25);
      min-width: 240px; max-width: 400px; border: 1px solid #e0dff0;
    }
    .node-cmp-popover-header {
      display: flex; align-items: center; justify-content: space-between; gap: 10px;
      padding: 8px 10px; border-bottom: 1px solid #f0eff8; background: #f8f7fc;
      border-radius: 6px 6px 0 0;
    }
    .node-cmp-popover-title { font-weight: 700; font-size: 13px; color: #2c2c2c; }
    .node-cmp-popover-close {
      background: none; border: none; font-size: 18px; line-height: 1; color: #767676;
      cursor: pointer; padding: 0 2px;
    }
    .node-cmp-popover-close:hover { color: #2c2c2c; }
    .node-cmp-popover-body {
      padding: 10px 12px; font-size: 13.5px; color: #2c2c2c; line-height: 1.6;
      max-height: 260px; overflow-y: auto;
    }
  ")))
)

dashboardPage(header,sidebar,body,title = "PathwayComparator",skin = "purple")
