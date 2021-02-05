start_MicrobIEM <- function() {
  # load required packages
  if("shiny" %in% rownames(installed.packages()) == FALSE) 
    install.packages("shiny")
  library(shiny)
  if("shinyjs" %in% rownames(installed.packages()) == FALSE) 
    install.packages("shinyjs")
  library(shinyjs)
  if("ggplot2" %in% rownames(installed.packages()) == FALSE) 
    install.packages("ggplot2")
  library(ggplot2)
  
  # load ui and server files
  source("MicrobIEM/ui.R"); 
  source("MicrobIEM/server.R"); 
  # start the shiny app
  shinyApp(ui,server)
}
