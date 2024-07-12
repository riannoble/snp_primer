library(shiny)

shinyServer(function(input, output, session) {
  # Copy the initial split_dfs from global.R
  split_dfs <- reactiveVal(split_dfs)

  observeEvent({
    input$shift
    input$desired_tm
    input$diff
    input$Heterodimer_tm
    input$Homodimer
    input$top
    input$hairpin
  }, {
    # Recompute the data frame based on new inputs
    df <- get_primer_candidates(primer, input$shift)
    df <- get_self_filter(df)
    df <- get_cross_filter(df)
    df <- get_complex_filter(df)
    split_dfs(split(df, df$groupID))
  })

  output$dynamic_tabs <- renderUI({
    current_split_dfs <- split_dfs()
    tabsetPanel(
      lapply(names(current_split_dfs), function(groupID) {
        tabPanel(
          title = paste("Group", groupID),
          DT::dataTableOutput(paste("table", groupID, sep = "_"))
        )
      })
    )
  })

  observe({
    current_split_dfs <- split_dfs()
    lapply(names(current_split_dfs), function(groupID) {
      output[[paste("table", groupID, sep = "_")]] <- DT::renderDataTable({
        DT::datatable(current_split_dfs[[groupID]])
      })
    })
  })
})
