################################################################################
#
# This is the user-interface definition of MicrobIEM 
#
################################################################################

# ------------------------------------------------------------------------------
# Install and load required packages 
# ------------------------------------------------------------------------------
# Install packages
packages_ui <- c("shinyjs", "DT", "plotly")
install.packages(setdiff(packages_ui, rownames(installed.packages())))
# Load packages
library(shinyjs)
library(DT)
library(plotly)

# ------------------------------------------------------------------------------
# Define re-used parameters
# ------------------------------------------------------------------------------
# Possible filter criteria for contamination filter
neg_ratio_steps <- c("ignore" = Inf, 
                      "2" = 2,
                      "1.5" = 1.5,
                      "1" = 1,
                      "0.5" = 0.5,
                      "0.1" = 0.1)
#	Allow upload of files with max. size of 1000 Mb
options(shiny.maxRequestSize = 1000*1024^2)

# ------------------------------------------------------------------------------
# UI main function
# ------------------------------------------------------------------------------
ui <- fluidPage(
  # Application title
  titlePanel("MicrobIEM v0.7"),
  tags$div("Need help? Check out the ", 
           a("documentation", target = "_blank", 
             href = "https://github.com/LuiseRauer/MicrobIEM"), 
           " and the ", 
           a("example dataset", target = "_blank", 
             href = "https://github.com/LuiseRauer/MicrobIEM/tree/main/MicrobIEM/test-data"), 
           "on Github."), 
  br(),
  useShinyjs(), # Use shinyjs package 
  
  # Creating a sidebar
  sidebarLayout(
    sidebarPanel( # Sidebar
      tabsetPanel( # Tabset
        id = "tabset",
        tabPanel( # First tab
          title = "Filtering",
          br(), # One empty line for aesthetics
          
          # --------------------------------------------------------------------
          # All input fields for data and filtering options
          # --------------------------------------------------------------------
          # Uploading the meta file
          fileInput(inputId = "metafile", label = "Choose meta table",
                    accept = c("text/csv", 
                               "text/comma-separated-values,text/plain",
                               ".csv")),
          # Uploading the feature table
          fileInput(inputId = "featurefile", label = "Choose feature table",
                    accept = c("text/csv", 
                               "text/comma-separated-values,text/plain",
                               ".csv")),
          # Choose visualization type for filtering
          selectInput(inputId = "visualization_type",
                      label = "Visualization", 
                      choices = c("Correlation of reads and features",
                                  "Change in feature abundance",
                                  "Contamination removal - NEG1",
                                  "Contamination removal - NEG2",
                                  "Reduction of total reads"),
                      selected = "Correlation of reads and features"),
          # Choose parameters for sample filter
          textInput(inputId = "req_reads_per_sample",
                    label = "Minimum reads per sample",
                    value = "1"),
          # Choose parameters for feature abundance filter
          textInput(inputId = "req_reads_per_feature",
                    label = "Minimum reads per feature",
                    value = "1"),
          # Choose parameters for feature relative frequency filter
          numericInput(inputId = "req_ratio_per_feature",
                       label = "Minimum relative frequency per feature",
                       value = 0, min = 0, max = 1, step = 0.0005),
          # Choose parameters for contamination filter - NEG1
          h4(tags$b("Contaminant filter based on NEG1"), id = "header_neg1"), 
          selectInput(inputId = "req_ratio_neg1",
                      label = "Frequency mean ratio (NEG1/SAMPLE)",		
                      choices = neg_ratio_steps),
          tags$div(id = "placeholder_span_1"),
          # Choose parameters for contamination filter - NEG2
          h4(tags$b("Contaminant filter based on NEG2"), id = "header_neg2"), 
          selectInput(inputId = "req_ratio_neg2",
                      label = "Frequency mean ratio (NEG2/SAMPLE)",		
                      choices = neg_ratio_steps),
          tags$div(id = "placeholder_span_2"),

          # --------------------------------------------------------------------
          # All buttons for filtering data 
          # --------------------------------------------------------------------
          actionButton(inputId = "update_button", label = "Update plot"),
          actionButton(inputId = "back_button", label = "Back"),
          actionButton(inputId = "next_button", label = "Next")
        ) # Close first tab
      ) # Close tabset
    ), # Close sidebar

    # --------------------------------------------------------------------------
    # Output a plot, a text output, and a table in the main panel
    # --------------------------------------------------------------------------
    mainPanel(
      plotlyOutput("plot"),
      textOutput("text"),
      DT::dataTableOutput("table")
    )
  )
)
