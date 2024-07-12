library(shiny)

shinyUI(fluidPage(
  titlePanel("Primer Candidates Viewer"),
  sidebarLayout(
    sidebarPanel(
      # Add inputs for parameters if needed
      numericInput("shift", "Shift:", value = 100),
      numericInput("desired_tm", "Desired Tm:", value = 64),
      numericInput("diff", "Tm Difference:", value = 3),
      numericInput("Heterodimer_tm", "Heterodimer Tm:", value = 50),
      numericInput("Homodimer", "Homodimer:", value = 45),
      numericInput("top", "Top Candidates:", value = 2),
      numericInput("hairpin", "Hairpin:", value = 45)
    ),
    mainPanel(
      uiOutput("dynamic_tabs")
    )
  )
))
