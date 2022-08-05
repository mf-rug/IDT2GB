# a shiny app to convert sequences with IDT-like (for the moment) modifications to .gb format
# load libraries
library(shiny)
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(DT))
library(shinyWidgets)
suppressPackageStartupMessages(library(shinyjs))
library(rclipboard)


ui <- fluidPage(
  useShinyjs(),
  titlePanel(title = "IDT2GB"),
  rclipboardSetup(),
  sidebarLayout(
    sidebarPanel(
      HTML('<strong><big>Input</strong></big>'), hr(),
      textAreaInput('seqs', label = NULL, placeholder = 'Enter names and sequences', height= '50vh'),
      HTML('<strong><big>Options</strong></big>'), hr(),
      fluidRow(column(12,
                      div(style="display: inline-block;vertical-align:top",HTML("<strong>name-sequence separator:&nbsp</strong>")),
                      div(style="display: inline-block;vertical-align:top",textInput('sep', NULL, value = '\\t', width = '100px'))
                      )),
      materialSwitch('remr', "remove 'r' from RNA sequences", value = TRUE, right = TRUE, status = 'primary'),
      materialSwitch('remp', "remove '+' from LNA sequences", value = TRUE, right = TRUE, status = 'primary'),
      materialSwitch('rems', "remove '*' from phosphorothioate sequences", right = TRUE, value = TRUE, status = 'primary'),
      materialSwitch('remm', "remove 'm' from 2'-O-methyl sequences", right = TRUE, value = TRUE, status = 'primary')
    ),
    mainPanel(
      HTML('<br><strong><big>Output</strong></big>'), hr(),
      DTOutput('outtable'), hr(),
      fluidRow(div(style="display: inline-block;vertical-align:top", HTML("&nbsp&nbsp&nbsp")),
               div(style="display: inline-block;vertical-align:top", hidden(uiOutput('clip'))), 
               div(style="display: inline-block;vertical-align:top", hidden(downloadButton('downloadgb', HTML('<small><font color="grey">Download .gb file</font></small>'), icon = icon('download'))))
      ), br(),
      htmlOutput('outtext')
    )
  )
)