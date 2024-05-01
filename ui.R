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
    #dark_hover_bg = "#81A1C1",
    #dark_color = "#2E3440"
  )
)

header <- dashboardHeader(title = span("PathwayComparator", style = "color: white; font-size: 25px; font-family: 'Segoe Print'"), titleWidth = 350)

sidebar <- dashboardSidebar(
  sidebarMenu(
    hidden(selectInput("uploadedFiles",label=NULL,choices=NULL,multiple=T)),
    DTOutput("listFiles"),
    fileInput("hiddenUpload",label=NULL,multiple=T,buttonLabel = "Add file"),
    actionButton("compare","Compare")
  )
)

body <- dashboardBody(
  use_theme(app.theme),
  useShinyjs(),
  fluidRow(
    column(4,selectizeInput("pathwaySel","Pathway",choices=NULL,width="100%")),
    column(2,selectizeInput("geneSel","Node",choices=NULL)),
    column(2,pickerInput("networkSel","Layer",choices=NULL,multiple = T,options=list(`max-options`=3))),
    column(4,checkboxGroupInput("hideElements", "", c("Hide chemical entities","Hide drugs","Hide miRNAs"),
                                selected = c("Hide chemical entities","Hide drugs","Hide miRNAs"), inline=T))
  ),
  fluidRow(
    column(12,visNetworkOutput("plotPathway",height="900px") %>% withSpinner())
  ),
  tags$head(tags$style(sass(sass_file("www/bootswatch-cyborg.scss")))),
  tags$head(tags$style(HTML('.box{-webkit-box-shadow: none; -moz-box-shadow: none;box-shadow: none; font-size: 18px;}'))),
  tags$head(tags$style(HTML(".main-sidebar { font-size: 18px;}"))),
  tags$head(tags$style(HTML(".sidebar { background-color: black;}"))),
  tags$head(tags$style(HTML(".control-label { font-size: 18px;}"))),
  tags$head(tags$style(HTML(".content-wrapper { background-color: white;}"))),
  tags$head(tags$style(HTML(".form-group { font-size: 18px;}"))),
  tags$head(tags$style(HTML("p { margin-left: 20px;}"))),
  tags$style("input[type=checkbox] { transform: scale(1.4); }"),
  tags$head(tags$style(HTML(".selectize-input { font-size: 18px; line-height: 20px;} .selectize-dropdown { font-size: 16px; line-height: 18px; }"))),
  tags$head(tags$style(HTML(".dropdown-menu ul li:nth-child(n) a { color: black !important; font-size: 18px;}"))),
  tags$head(tags$style(HTML(".picker {font-size: 16px; color: black;}")))
)

dashboardPage(header,sidebar,body,title = "PathwayComparator",skin = "purple")

