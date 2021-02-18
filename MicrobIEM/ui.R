################################################################################
#
# This is the user-interface definition of MicrobIEM 
#
################################################################################


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
neg_span_steps <- c("ignore" = 0.0001,
                    "100 %" = 1, 
                    "80 %" = 0.8, 
                    "75 %" = 0.75, 
                    "60 %" = 0.6, 
                    "50 %" = 0.5, 
                    "40 %" = 0.4, 
                    "25 %" = 0.25, 
                    "20 %" = 0.2)

#	Allow upload of files with max. size of 50 Mb
options(shiny.maxRequestSize = 50*1024^2)

# ------------------------------------------------------------------------------
# UI main function
# ------------------------------------------------------------------------------

ui <- fluidPage(
  # Application title
  titlePanel("MicrobIEM"),
  useShinyjs(), # Use shinyjs package 
  
  # Creating a sidebar
  sidebarLayout(
    sidebarPanel( # Sidebar
      tabsetPanel( # Tabset
        id = "tabset",
        tabPanel( # First tab
          title = "Filtering",
          br(), # One empty line for aesthetics
          
          # -------------------------------------------------------------------
          # All input fields for data and filtering options
          # -------------------------------------------------------------------
          
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
          
          # Start button for analysis
          #actionButton(inputId = "start_button", label = "Start"),
          
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
          tags$div(id = "placeholder"), # NEWLINE
          selectInput(inputId = "req_ratio_neg1",
                      label = "Frequency mean ratio (NEG1/SAMPLE)",		
                      choices = neg_ratio_steps),
          selectInput(inputId = "req_span_neg1", 
                      label = "Minimum span threshold (NEG1)",
                      choices = neg_span_steps),
          
          # Choose parameters for contamination filter - NEG2
          h4(tags$b("Contaminant filter based on NEG2"), id = "header_neg2"), 
          selectInput(inputId = "req_ratio_neg2",
                      label = "Frequency mean ratio (NEG2/SAMPLE)",		
                      choices = neg_ratio_steps),
          selectInput(inputId = "req_span_neg2", 
                      label = "Minimum span threshold (NEG2)",
                      choices = neg_span_steps),
          
          # -------------------------------------------------------------------
          # All buttons for filtering data 
          # -------------------------------------------------------------------
          
          actionButton(inputId = "update_button", label = "Update plot"),
          actionButton(inputId = "back_button", label = "Back"),
          actionButton(inputId = "next_button", label = "Next")
          
        ) # Close first tab
      ) # Close tabset
    ), # Close sidebar
    
    # Show a plot and a table in the main panel
    mainPanel(
      plotlyOutput("plot"),
      textOutput("text"),
      DT::dataTableOutput("table")
    )
  )
)
