
library(shiny)
library(ggplot2)
library(ggridges)
library(gridExtra)
library(grid)
library(tidyverse)
library(ggpubr)
library(stringr)
library(scales)
library(egg)
library(zoo)
library(ggrepel)
library(reshape2)


source(file = '/home/sbp29/R/Projects/ybr/scripts/pull_down_results/analysis.R')
load(file = '/home/sbp29/R/Projects/ybr/scripts/pull_down_results/pddata.RData')

# Define UI for application
ui <- fluidPage(
  
  # Application title
  titlePanel("YBR Pull Down Results"),
  
  sidebarLayout(
    sidebarPanel(
      selectInput(inputId = "combine_results", 
                  label = "1. Combine replicate-wise results?",  # Give the input a label to be displayed in the app
                  choices = list("Yes" = "Y",
                                 "No" = "N"),
                  selected = "N"),
      selectInput(inputId = "test_group", 
                  label = "2. Select Test Group/s",  # Give the input a label to be displayed in the app
                  choices = unique(pddat$ID),
                  selected = unique(pddat$ID), multiple = T),
      selectInput(inputId = "control_group",
                  label = "3. Select Control Group/s",  # Give the input a label to be displayed in the app
                  choices = unique(pddat$ID),
                  selected = unique(pddat$ID), multiple = T),
      # radioButtons(inputId = "gc",
      #              label = "4. Growth Curves",
      #              choices = list("Yes" = "Y",
      #                             "No" = "N"),
      #              selected = "N"),
      # sliderInput(inputId = "n_col",
      #             label = "5. Modify Growth Curve Plot",
      #             min = 1, max = 5,
      #             value = 1),
      
      submitButton(text = "Analyze", icon = NULL, width = NULL),
      width = 3),
    
    # Show the output
    mainPanel(
      plotOutput("upset")
      # br(),
      # plotOutput("boxPlot")
    )
    # fluidRow(
    #   column(12,tableOutput('overlap_results'))
    # )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  overlap_results <- reactive(ms_overlap(input$combine_results, input$test_group, input$control_group))
  # overlap_results <- overlap_results()
  output$upset <- renderPlot({overlap_results()}, height = 900, width = 1200)
  

  # output$overlap_results <- renderTable(overlap_results[2])
}

# Run the application 
shinyApp(ui = ui, server = server)


