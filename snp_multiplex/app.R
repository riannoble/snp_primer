# Load required libraries
library(shiny)
library(dplyr)

# Source your functions script
source("C:/Users/riano/Documents/acornfinder/acornfinder/R/ria_fxns_draft.R")

# Define UI for application
ui <- fluidPage(

  # Application title
  titlePanel("SNP Analysis App"),

  # Sidebar layout with input and output definitions
  sidebarLayout(
    sidebarPanel(
      # Input for SNP ID (string)
      textInput("primer", "Enter Primer:", value = ""),

      # Input for shift (number)
      numericInput("shift", "Enter Shift Number:", value = 0, min = 0),

      # Action button to trigger calculation
      actionButton("calculate_button", "Calculate")
    ),

    # Main panel for displaying output table
    mainPanel(
      tableOutput("output_table")
    )
  )
)

# Define server logic
server <- function(input, output) {

  # Function to calculate and render output table
  observeEvent(input$calculate_button, {
    req(input$primer, input$shift)  # Require both inputs to be filled

    df <- get_primer_candidates(input$primer, input$shift)
    df <- get_self_filter(df)
    df <- get_cross_filter(df)
    output_data <- get_final_list(df)

    output$output_table <- DT::renderDT({
      datatable(output_data)  # Render output_data as a datatable
    })
  })
}

# Run the application
shinyApp(ui = ui, server = server)
