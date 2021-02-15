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
  if("vegan" %in% rownames(installed.packages()) == FALSE) 
    install.packages("vegan")
  library(vegan)
  if("dplyr" %in% rownames(installed.packages()) == FALSE) 
    install.packages("dplyr")
  library(dplyr)
  if("shinyWidgets" %in% rownames(installed.packages()) == FALSE) 
    install.packages("shinyWidgets")
  library(shinyWidgets)
  if("reshape2" %in% rownames(installed.packages()) == FALSE) 
    install.packages("reshape2")
  library(reshape2)
  if("plotly" %in% rownames(installed.packages()) == FALSE) 
    install.packages("plotly")
  library(plotly)
  
  
  # load ui and server files
  source("MicrobIEM/ui.R"); 
  source("MicrobIEM/server.R"); 
  # start the shiny app
  shinyApp(ui,server)
}
