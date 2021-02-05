################################################################################
#
# This is the server logic of MicrobIEM
#
################################################################################

sample_types_allowed <- c("SAMPLE", "POS1", "POS2", "NEG1", "NEG2")

# ------------------------------------------------------------------------------
# Server logic main function
# ------------------------------------------------------------------------------

server <- function(input, output, session) {
  
  # Define reactive values (count steps etc.)
  reactives <- reactiveValues(step_var = 1, metadata = NA, featuredata = NA,
                              output_dir = NA)
  
  # Hide buttons and input fields in the beginning
  shinyjs::hide("visualization_type")
  shinyjs::hide("req_reads_per_sample") 
  shinyjs::hide("req_reads_per_feature")
  shinyjs::hide("req_ratio_per_feature")
  shinyjs::hide("header_neg1") 
  shinyjs::hide("req_ratio_neg1")
  shinyjs::hide("req_span_neg1")
  shinyjs::hide("header_neg2")
  shinyjs::hide("req_ratio_neg2")
  shinyjs::hide("req_span_neg2")
  shinyjs::hide("update_button")
  shinyjs::hide("back_button")
  shinyjs::hide("next_button")
  
  # Read and check the metafile
  observeEvent(input$metafile, {
    if (!(tools::file_ext(input$metafile$datapath) %in% c("csv", "txt"))) {
      showModal(modalDialog(
        title = "Error1", "Please choose a csv or txt file as meta table."))
    } else {
      reactives$metadata <- read.csv(input$metafile$datapath, sep = "\t", 
                                     header = TRUE, check.names = FALSE)
      #print(isolate(colnames(reactives$metadata)))
      colnames_md <- colnames(reactives$metadata)
      if (colnames_md[1] != "Sample_ID" || !("Sample_type" %in% colnames_md)) {
        showModal(modalDialog(
          title = "Error2", 
          "Please provide a meta information table with the first column 
          'Sample_ID' and a column 'Sample_type' to define samples and controls"))
        reactives$metadata <- NA
      } else {
        sample_types_observed <- unique(reactives$metadata[, "Sample_type"])
        if ((sum(sample_types_observed %in% sample_types_allowed) != 
             length(sample_types_observed))) {
          showModal(modalDialog(
            title = "Error3", 
            paste0("Please use only the following labels to define samples and 
                   controls: ", 
                   paste(sample_types_allowed, collapse = ", "))))
          reactives$metadata <- NA
        } else if (!is.na(reactives$featuredata) && !is.na(reactives$metadata)) {
          withProgress(message = "Data upload", file_open_success())
        }
      }
    }
  })
  
  # Read and check the feature table
  observeEvent(input$featurefile, {
    if (!(tools::file_ext(input$featurefile$datapath) %in% c("csv", "txt"))) {
      showModal(modalDialog(
        title = "Error4", "Please choose a csv or txt file as feature table."))
    } else {
      reactives$featuredata <- read.csv(input$featurefile$datapath, sep = "\t", 
                                        header = TRUE, check.names = FALSE)
      print(isolate(colnames(reactives$featuredata)))
      colnames_fd <- colnames(reactives$featuredata)
      if(length(colnames_fd) < 3 || colnames_fd[1] != "OTU_ID" || 
         colnames_fd[length(colnames_fd)] != "Taxonomy") {
        showModal(modalDialog(
          title = "Error5", "Please provide a feature table with the first 
          column 'OTU_ID', at least one sample, and the last column 'Taxonomy'"))
        reactives$featuredata <- NA
      } else if (!is.na(reactives$featuredata) && !is.na(reactives$metadata)) {
        withProgress(message = "Data upload", file_open_success())
      }
    }
  })
  
  # Define the file_open_success function
  file_open_success <- function() {
    sample_names_md <- reactives$metadata[, "Sample_ID"]
    sample_names_fd <- head(colnames(reactives$featuredata)[-1], -1)
    if (identical(sort(sample_names_md), sort(sample_names_fd))) {
      showModal(modalDialog(
        title = "Success2", "Meta table and feature table were uploaded. 
        Please note that filtering can take several minutes."))
      make_plot_filtering()
      shinyjs::show("update_button")
      shinyjs::show("back_button")
      shinyjs::show("next_button")
      shinyjs::show("visualization_type")
      # Create output directory
      if(is.na(reactives$output_dir)){
        reactives$output_dir=gsub(":","_",format(Sys.time(), "%Y_%m_%d_%a_%X"))
      }
      # Replace NA values in the metafile
      if (sum(is.na(reactives$metadata)) > 0) {
        reactives$metadata[is.na(reactives$metadata)] <- "n.a." 
        print(c("INFO|server|replaced ", sum(is.na(reactives$metadata)), 
                " in columns: ", 
                colnames(reactives$metadata)[colSums(is.na(reactives$metadata)) > 0]))
      }
    } else {
      showModal(modalDialog(
        title = "Error6", "Sample IDs do not match between meta table and 
        feature table."))
      #reactives$metadata <- NA
      #reactives$featuredata <- NA
    }
  }

  make_plot_filtering <- function() { # This will create the plot
    filter_feature_table(reactives$step_var, 
                         input$req_reads_per_sample,
                         input$req_reads_per_feature,
                         input$req_ratio_per_feature,
                         input$req_ratio_neg1,
                         input$req_span_neg1,
                         input$req_ratio_neg2,
                         input$req_span_neg2) # This will return the filtered feature table+params
    
    # Load data set depending on step_var
    # Design plot based on chosen plot
  }
  
  filter_feature_table <- function(a, b, c, d, e, f, g, h){
    showModal(modalDialog(title = "BLA", "filter_feature_table"))
  }
  
  observeEvent(input$update_button, {
    if(isolate(is.na(reactives$metadata))) {
      print("No file")
    } else {
      print("Success")
    }
    showModal(modalDialog(title = "BLA", "update_button"))
  })
  
  observeEvent(input$next_button, {
    print(paste0("INFO|server::nextstep",Sys.time()))
    reactives$step_var <- reactives$step_var +1
    showModal(modalDialog(title = "BLA", paste0("next_button", reactives$step_var)))
  })

  observeEvent(input$back_button, {
    print(paste0("INFO|server::backstep",Sys.time()))
    reactives$step_var <- reactives$step_var - 1
    showModal(modalDialog(title = "BLA", paste0("back_button", reactives$step_var)))
  })
  
  output$plot <- renderPlot({
    set.seed(1)  
    plot(rnorm(10, 0, 2), rnorm(10, as.numeric(input$req_reads_per_sample), 
                                isolate(reactives$step_var)))
  })
  
}
