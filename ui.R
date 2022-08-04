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
      textAreaInput('seqs', label = NULL, placeholder = 'Enter names and sequences', height= '50vh')
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