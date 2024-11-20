library(shiny)
library(dplyr)
library(purrr)
library(here)

#print(getwd())

#setwd(".../acornfinder.Rproj")

source("C:/Users/riano/Documents/acornfinder/acornfinder/R/fxns_draft.R")

# UI
ui <- fluidPage(
  titlePanel("Primer Analysis App"),
  sidebarLayout(
    sidebarPanel(
      textInput("primer", "Enter SNPs:", value = ""),
      actionButton("run_analysis", "Run Analysis"),
      downloadButton("download_data", "Download Results")
    ),
    mainPanel(
      tableOutput("results_table")
    )
  )
)

# Server
server <- function(input, output, session) {

  # Reactive expression to store the dataset
  results <- eventReactive(input$run_analysis, {
    primer_string <- input$primer
    shift = 100

    # Assuming the input string is used to generate a dataframe (e.g., `df`)
    # Replace this with actual data processing logic to generate `df`

    df <- get_primer_candidates(primer, shift)
    df <- get_self_filter(df)
    new_df <- get_cross_filter(df)
    final_combinations <- get_final_list(new_df)

    return(final_combinations)
  })

  # Render the results table
  output$results_table <- renderTable({
    req(results())
    head(results(), 10)  # Display the first 10 rows of the dataset
  })

  # Download handler
  output$download_data <- downloadHandler(
    filename = function() {
      paste("primer_analysis_results-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(results(), file, row.names = FALSE)
    }
  )
}

# Run the app
shinyApp(ui = ui, server = server)
