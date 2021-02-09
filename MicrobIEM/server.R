################################################################################
#
# This is the server logic of MicrobIEM
#
################################################################################

# ------------------------------------------------------------------------------
# Define re-used parameters
# ------------------------------------------------------------------------------

sample_types_allowed <- c("SAMPLE", "POS1", "POS2", "NEG1", "NEG2")

# ------------------------------------------------------------------------------
# Server logic main function
# ------------------------------------------------------------------------------

server <- function(input, output, session) {
  
  # Define reactive values (count steps etc.)
  reactives <- reactiveValues(step_var = 1, output_dir = NA,
                              metadata = NA, featuredata = NA,
                              metadata_current = NA, featuredata_current = NA,
                              metadata_1 = NA, featuredata_1 = NA,
                              metadata_2 = NA, featuredata_2 = NA,
                              metadata_3 = NA, featuredata_3 = NA,
                              req_reads_per_sample_old = "1", 
                              req_reads_per_feature_old = "1",
                              req_ratio_per_feature_old = 0,
                              req_ratio_neg1_old = -1, req_span_neg1_old = -1,
                              req_ratio_neg2_old = -1, req_span_neg2_old = -1)
  
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

  # ----------------------------------------------------------------------------
  # Read and check the metafile
  # ----------------------------------------------------------------------------
  
  observeEvent(input$metafile, {
    # Check correct file extension
    if (!(tools::file_ext(input$metafile$datapath) %in% c("csv", "txt"))) {
      showModal(modalDialog(
        title = "Error1", "Please choose a csv or txt file as meta table."))
    } else {
      reactives$metadata <- read.csv(input$metafile$datapath, sep = "\t", 
                                     header = TRUE, check.names = FALSE)
      colnames_md <- colnames(reactives$metadata)
      # Check for columns Sample_ID and Sample_type
      if (colnames_md[1] != "Sample_ID" || !("Sample_type" %in% colnames_md)) {
        showModal(modalDialog(
          title = "Error2", 
          "Please provide a meta information table with the first column 
          'Sample_ID' and a column 'Sample_type' to define samples and controls"))
        reactives$metadata <- NA
      } else {
        sample_types_observed <- unique(reactives$metadata[, "Sample_type"])
        # Check for correct sample type definitions
        if ((sum(sample_types_observed %in% sample_types_allowed) != 
             length(sample_types_observed))) {
          showModal(modalDialog(
            title = "Error3", 
            paste0("Please use only the following labels to define samples and 
                   controls: ", 
                   paste(sample_types_allowed, collapse = ", "))))
          reactives$metadata <- NA
        } else if (!is.na(reactives$featuredata) && !is.na(reactives$metadata)) {
          # If everything is correct, start the file_open_success function
          withProgress(message = "Data upload", file_open_success())
        }
      }
    }
  })
  
  # ----------------------------------------------------------------------------
  # Read and check the featurefile
  # ----------------------------------------------------------------------------
  observeEvent(input$featurefile, {
    # Check correct file extension
    if (!(tools::file_ext(input$featurefile$datapath) %in% c("csv", "txt"))) {
      showModal(modalDialog(
        title = "Error4", "Please choose a csv or txt file as feature table."))
    } else {
      reactives$featuredata <- read.csv(input$featurefile$datapath, sep = "\t", 
                                        header = TRUE, check.names = FALSE)
      print(isolate(colnames(reactives$featuredata)))
      colnames_fd <- colnames(reactives$featuredata)
      # Check for columns OTU_ID and Taxonomy
      if(length(colnames_fd) < 3 || colnames_fd[1] != "OTU_ID" || 
         colnames_fd[length(colnames_fd)] != "Taxonomy") {
        showModal(modalDialog(
          title = "Error5", "Please provide a feature table with the first 
          column 'OTU_ID', at least one sample, and the last column 'Taxonomy'"))
        reactives$featuredata <- NA
      } else if (!is.na(reactives$featuredata) && !is.na(reactives$metadata)) {
        # If everything is correct, start the file_open_success function
        withProgress(message = "Data upload", file_open_success())
      }
    }
  })
  
  # ----------------------------------------------------------------------------
  # Define the file_open_success function
  # ----------------------------------------------------------------------------
  file_open_success <- function() {
    print("Start open_success_function")
    sample_names_md <- reactives$metadata[, "Sample_ID"]
    sample_names_fd <- head(colnames(reactives$featuredata)[-1], -1)
    # If sample names match in feature file and meta file, start filtering
    if (identical(sort(sample_names_md), sort(sample_names_fd))) {
      showModal(modalDialog(
        title = "Success2", "Meta table and feature table were uploaded. 
        Please note that filtering can take several minutes."))
      # Replace NA values in the metafile
      if (sum(is.na(reactives$metadata)) > 0) {
        reactives$metadata[is.na(reactives$metadata)] <- "n.a." 
        print(c("INFO|server|replaced ", sum(is.na(reactives$metadata)), 
                " in columns: ", 
                colnames(reactives$metadata)[colSums(is.na(reactives$metadata)) > 0]))
      }
      # Create output directory
      if (is.na(reactives$output_dir)) {
        reactives$output_dir <- gsub(":", "_", format(Sys.time(), 
                                                      "%Y_%m_%d_%a_%X"))
        if (!dir.exists(paste0(reactives$output_dir, "/output"))){
          dir.create(paste0(reactives$output_dir, "/output"), recursive = TRUE)
          print(paste0("INFO|server|output dir created", Sys.time()))
        } else {
          print(paste0("INFO|server|output dir not created", Sys.time()))
        }
      }
      # Save original files in output directory as ..._0?
      
      # Define rownames of featuredata and metadata
      rownames(reactives$metadata) <- reactives$metadata[, "Sample_ID"]
      rownames(reactives$featuredata) <- reactives$featuredata[, "OTU_ID"]
      # Remove controls from the data set
      # Define samples without controls
      real_sample_IDs <- 
        reactives$metadata[which(reactives$metadata[, "Sample_type"] == 
                                   sample_types_allowed[1]), "Sample_ID"]
      # Filter featuredata and metadata
      reactives$featuredata_1 <- 
        reactives$featuredata[, c("OTU_ID", real_sample_IDs, "Taxonomy")]
      reactives$metadata_1 <- 
        reactives$metadata[reactives$metadata$Sample_ID %in% real_sample_IDs, ]
      # Proceed one step - set step var to 2
      reactives$step_var <- 2
      # Start filtering the data set
      filter_feature_table(input$req_reads_per_sample,
                           input$req_reads_per_feature,
                           input$req_ratio_per_feature,
                           input$req_ratio_neg1,
                           input$req_span_neg1,
                           input$req_ratio_neg2,
                           input$req_span_neg2)
      # Make buttons appear
      shinyjs::show("update_button")
      shinyjs::show("back_button")
      shinyjs::show("next_button")
      shinyjs::show("visualization_type")
      shinyjs::show("req_reads_per_sample")
    } else {
      showModal(modalDialog(
        title = "Error6", "Sample IDs do not match between meta table and 
        feature table."))
    }
  }

  # ----------------------------------------------------------------------------
  # Define the filter_feature_table function
  # ----------------------------------------------------------------------------
  filter_feature_table <- function(req_reads_per_sample,
                                   req_reads_per_feature,
                                   req_ratio_per_feature,
                                   req_ratio_neg1,
                                   req_span_neg1,
                                   req_ratio_neg2,
                                   req_span_neg2){

    # --------------------------------------------------------------------------
    # Filter step 2: remove samples
    # --------------------------------------------------------------------------
    if(reactives$step_var == 2 || 
       reactives$req_reads_per_sample_old != input$req_reads_per_sample) {
      print("filtering_2 branch started")
      # Define samples to keep
      sample_read_sums <- colSums(
        reactives$featuredata_1[, 2:(ncol(reactives$featuredata_1)-1)])
      samples_to_keep <- names(
        sample_read_sums[sample_read_sums >= as.numeric(req_reads_per_sample)])
      # Filter featuredata and metadata
      reactives$featuredata_2 <- 
        reactives$featuredata_1[, c("OTU_ID", samples_to_keep, "Taxonomy")]
      reactives$metadata_2 <- 
        reactives$metadata_1[reactives$metadata_1$Sample_ID %in% samples_to_keep, ]
      reactives$featuredata_current <- reactives$featuredata_2
      reactives$metadata_current <- reactives$metadata_2
    }
    
    # --------------------------------------------------------------------------
    # Filter step 3: remove features by abundance
    # --------------------------------------------------------------------------
    if(reactives$step_var == 3 || reactives$req_reads_per_feature_old != input$req_reads_per_feature) {
      print("filtering_3 branch started")
      # Define feature to keep
      feature_read_sums <- rowSums(
        reactives$featuredata_2[, 2:(ncol(reactives$featuredata_2)-1)])
      feature_to_keep_abund <- names(
        feature_read_sums[feature_read_sums >= as.numeric(req_reads_per_feature)])
      # Filter featuredata and metadata
      reactives$featuredata_3 <- 
        reactives$featuredata_2[feature_to_keep_abund, ]
      # Remove empty samples:
      sample_read_sums <- colSums(
        reactives$featuredata_3[, 2:(ncol(reactives$featuredata_3)-1)])
      print(sample_read_sums)
      samples_to_keep <- names(sample_read_sums[sample_read_sums > 0])
      print(samples_to_keep)
      reactives$featuredata_3 <- 
        reactives$featuredata_3[, c("OTU_ID", samples_to_keep, "Taxonomy")]
      reactives$metadata_3 <- 
        reactives$metadata_2[reactives$metadata_2$Sample_ID %in% 
                               colnames(reactives$featuredata_3), ]
      reactives$featuredata_current <- reactives$featuredata_3
      reactives$metadata_current <- reactives$metadata_3
    }
    
    # --------------------------------------------------------------------------
    # Filter step 4: remove features by frequency
    # --------------------------------------------------------------------------
    if(reactives$step_var == 4 ||
       reactives$req_ratio_per_feature_old != input$req_ratio_per_feature) {
      print("filtering_4 branch started")
      # Convert feature table to frequencies
      featuredata_3_freq <- 
        t(decostand(t(reactives$featuredata_3[, 2:(ncol(reactives$featuredata_3)-1)]), 
                    method = "total"))
      # Define feature to keep
      feature_read_freqs <- apply(
        featuredata_3_freq[, 2:(ncol(featuredata_3_freq)-1)], 1, max)
      feature_to_keep_freq <- names(
        feature_read_freqs[feature_read_freqs >= as.numeric(req_ratio_per_feature)])
      # Filter featuredata and metadata
      reactives$featuredata_4 <- 
        reactives$featuredata_3[feature_to_keep_freq, ]
      # Remove empty samples:
      sample_read_sums <- colSums(
        reactives$featuredata_4[, 2:(ncol(reactives$featuredata_4)-1)])
      print(sample_read_sums)
      samples_to_keep <- names(sample_read_sums[sample_read_sums > 0])
      print(samples_to_keep)
      reactives$featuredata_4 <- 
        reactives$featuredata_4[, c("OTU_ID", samples_to_keep, "Taxonomy")]
      reactives$metadata_4 <- 
        reactives$metadata_3[reactives$metadata_3$Sample_ID %in% 
                               colnames(reactives$featuredata_4), ]
      reactives$featuredata_current <- reactives$featuredata_4
      reactives$metadata_current <- reactives$metadata_4
    }
    
    # --------------------------------------------------------------------------
    # Filter step 5: remove contaminants
    # --------------------------------------------------------------------------
    if(reactives$step_var == 5 ||
       reactives$req_ratio_neg1_old != input$req_ratio_neg1 ||
       reactives$req_span_neg1_old != input$req_span_neg1 ||
       reactives$req_ratio_neg2_old != input$req_ratio_neg2 ||
       reactives$req_span_neg2_old != input$req_span_neg2) {
      print("filtering_5 branch started")
      # Convert feature table to frequencies
      featuredata_4_freq <- 
        t(decostand(t(reactives$featuredata_4[, 2:(ncol(reactives$featuredata_4)-1)]), 
                    method = "total"))
      # Define types of controls (again :/)
      sample_IDs_neg1 <- 
        reactives$metadata[which(reactives$metadata[, "Sample_type"] == 
                                   sample_types_allowed[4]), "Sample_ID"]
      sample_IDs_neg2 <- 
        reactives$metadata[which(reactives$metadata[, "Sample_type"] == 
                                   sample_types_allowed[5]), "Sample_ID"]
      real_sample_IDs <- 
        reactives$metadata[which(reactives$metadata[, "Sample_type"] == 
                                   sample_types_allowed[1]), "Sample_ID"]
      # Define mean values
      sample_mean <- apply(
        featuredata_4_freq[, colnames(featuredata_4_freq) %in% real_sample_IDs], 
        1, mean)
      # NEG1:
      neg1_data <-
        reactives$featuredata[, colnames(reactives$featuredata) %in% sample_IDs_neg1]
      neg1_data_rel <- 
        t(decostand(t(neg1_data), method = "total"))
      neg1_mean <- apply(neg1_data_rel, 1, mean)
      neg1_span <- apply(neg1_data, 1, function(x) sum(x > 0)/length(x))
      # NEG2:
      neg2_data <-
        reactives$featuredata[, colnames(reactives$featuredata) %in% sample_IDs_neg2]
      neg2_data_rel <- 
        t(decostand(t(neg2_data), method = "total"))
      neg2_mean <- apply(neg2_data_rel, 1, mean)
      neg2_span <- apply(neg2_data, 1, function(x) sum(x > 0)/length(x))
      # Actual filtering
      # Save mean and span values in data frame
      # Calculate mean ratios
      # Apply filter criteria, return Sample IDs that should be kept
    }
    
    # provide information for visualisation
    if(input$visualization_type == "Correlation of reads and features") {
      data_to_plot <- data.frame(
        reads = colSums(reactives$featuredata_current[, 2:(
          ncol(reactives$featuredata_current)-1)]),
        features = apply(reactives$featuredata_current[, 2:(
          ncol(reactives$featuredata_current)-1)], 2, function(x) sum(x > 0)))
      print(str(data_to_plot))
      output$plot <- renderPlot({
        ggplot(data = data_to_plot, aes(x = reads, y = features)) +
          geom_point()
      }) 
    }
    # Save current input:
    reactives$req_reads_per_sample_old <- input$req_reads_per_sample
    reactives$req_reads_per_feature_old <- input$req_reads_per_feature
    reactives$req_ratio_per_feature_old <- input$req_ratio_per_feature
    reactives$req_ratio_neg1_old <- input$req_ratio_neg1
    reactives$req_span_neg1_old <- input$req_span_neg1
    reactives$req_ratio_neg2_old <- input$req_ratio_neg2
    reactives$req_span_neg2_old <- input$req_span_neg2
  }
  
  observeEvent(input$update_button, {
    filter_feature_table(input$req_reads_per_sample,
                         input$req_reads_per_feature,
                         input$req_ratio_per_feature,
                         input$req_ratio_neg1,
                         input$req_span_neg1,
                         input$req_ratio_neg2,
                         input$req_span_neg2)
  })
  
  observeEvent(input$next_button, {
    print(paste0("INFO|server::nextstep",Sys.time()))
    reactives$step_var <- reactives$step_var + 1
    filter_feature_table(input$req_reads_per_sample,
                         input$req_reads_per_feature,
                         input$req_ratio_per_feature,
                         input$req_ratio_neg1,
                         input$req_span_neg1,
                         input$req_ratio_neg2,
                         input$req_span_neg2)
    showModal(modalDialog(title = "BLA", paste0("next_button", reactives$step_var)))
    step_var_UIchange()
    print(paste0("OBJECT SIZE: ", object.size(reactives)))
  })

  observeEvent(input$back_button, {
    print(paste0("INFO|server::backstep",Sys.time()))
    reactives$step_var <- reactives$step_var - 1
    showModal(modalDialog(title = "BLA", paste0("back_button", reactives$step_var)))
    step_var_UIchange()
  })
  
  # Step_var_change function: build UI (and disable buttons depending on step)
  step_var_UIchange <- function(){
    if(reactives$step_var == 1){
      shinyjs::enable("metafile")
      shinyjs::enable("featurefile")
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
    }
    if(reactives$step_var == 2){
      shinyjs::disable("metafile")
      shinyjs::disable("featurefile")
      shinyjs::show("req_reads_per_sample") 
      shinyjs::enable("req_reads_per_sample") 
      shinyjs::hide("req_reads_per_feature")
      shinyjs::hide("req_ratio_per_feature")
      shinyjs::hide("header_neg1") 
      shinyjs::hide("req_ratio_neg1")
      shinyjs::hide("req_span_neg1")
      shinyjs::hide("header_neg2")
      shinyjs::hide("req_ratio_neg2")
      shinyjs::hide("req_span_neg2")
    }
    if(reactives$step_var == 3) {
      shinyjs::disable("metafile")
      shinyjs::disable("featurefile")
      shinyjs::disable("req_reads_per_sample") 
      shinyjs::show("req_reads_per_feature")
      shinyjs::enable("req_reads_per_feature")
      shinyjs::hide("req_ratio_per_feature")
      shinyjs::hide("header_neg1") 
      shinyjs::hide("req_ratio_neg1")
      shinyjs::hide("req_span_neg1")
      shinyjs::hide("header_neg2")
      shinyjs::hide("req_ratio_neg2")
      shinyjs::hide("req_span_neg2")
    }
    if(reactives$step_var == 4) {
      shinyjs::disable("metafile")
      shinyjs::disable("featurefile")
      shinyjs::disable("req_reads_per_sample") 
      shinyjs::disable("req_reads_per_feature")
      shinyjs::show("req_ratio_per_feature")
      shinyjs::enable("req_ratio_per_feature")
      shinyjs::hide("header_neg1") 
      shinyjs::hide("req_ratio_neg1")
      shinyjs::hide("req_span_neg1")
      shinyjs::hide("header_neg2")
      shinyjs::hide("req_ratio_neg2")
      shinyjs::hide("req_span_neg2")
    }
    if(reactives$step_var == 5) {
      shinyjs::disable("metafile")
      shinyjs::disable("featurefile")
      shinyjs::disable("req_reads_per_sample") 
      shinyjs::disable("req_reads_per_feature")
      shinyjs::disable("req_ratio_per_feature")
      shinyjs::show("header_neg1") 
      shinyjs::show("req_ratio_neg1")
      shinyjs::show("req_span_neg1")
      shinyjs::show("header_neg2")
      shinyjs::show("req_ratio_neg2")
      shinyjs::show("req_span_neg2")
    }
  }
  
  #output$plot <- renderPlot({
  #  set.seed(1)  
  #  plot(rnorm(10, 0, 2), rnorm(10, as.numeric(input$req_reads_per_sample), 
  #                              isolate(reactives$step_var)))
  #})
  
}
