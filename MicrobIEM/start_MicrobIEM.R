start_MicrobIEM <- function() {
  # load required packages
  if("shiny" %in% rownames(installed.packages()) == FALSE) 
    install.packages("shiny")
  library(shiny)
  
  # load ui and server files
  source("MicrobIEM/ui.R"); 
  source("MicrobIEM/server.R"); 
  # start the shiny app
  shinyApp(ui,server)
}