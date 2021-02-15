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
                              featuredata_taxonomy = NA,
                              featuredata_neg1 = NA, featuredata_neg2 = NA,
                              metadata_current = NA, featuredata_current = NA,
                              metadata_1 = NA, featuredata_1 = NA,
                              metadata_2 = NA, featuredata_2 = NA,
                              metadata_3 = NA, featuredata_3 = NA,
                              metadata_4 = NA, featuredata_4 = NA,
                              metadata_5 = NA, featuredata_5 = NA,
                              req_reads_per_sample_old = "1", 
                              req_reads_per_feature_old = "1",
                              req_ratio_per_feature_old = 0,
                              req_ratio_neg1_old = Inf, req_span_neg1_old = 0.0001,
                              req_ratio_neg2_old = Inf, req_span_neg2_old = 0.0001)
  
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
      # Create output directory
      if (is.na(reactives$output_dir)) {
        reactives$output_dir <- gsub(":", "_", format(Sys.time(), 
                                                      "%Y_%m_%d_%a_%X"))
        if (!dir.exists(reactives$output_dir)){
          dir.create(paste0(reactives$output_dir, "/output"), recursive = TRUE)
          dir.create(paste0(reactives$output_dir, "/interim"), recursive = TRUE)
          print(paste0("INFO|server|output dir created", Sys.time()))
        } else {
          print(paste0("INFO|server|output dir not created", Sys.time()))
        }
      }
      
      # ------------------------------------------------------------------------
      # Prepare the data for filtering 
      # ------------------------------------------------------------------------
      # Replace NA values in the metafile
      if (sum(is.na(reactives$metadata)) > 0) {
        reactives$metadata[is.na(reactives$metadata)] <- "n.a." 
        print(c("INFO|server|replaced ", sum(is.na(reactives$metadata)), 
                " in columns: ", 
                colnames(reactives$metadata)[colSums(is.na(reactives$metadata)) > 0]))
      }
      # Define rownames of featuredata and metadata
      rownames(reactives$metadata) <- reactives$metadata[, "Sample_ID"]
      rownames(reactives$featuredata) <- reactives$featuredata[, "OTU_ID"]
      # Remove OTU_ID column and separate taxonomy columns in the feature file
      reactives$featuredata_taxonomy <- data.frame(
        Taxonomy = reactives$featuredata[, "Taxonomy", drop = FALSE])
      reactives$featuredata[, c("OTU_ID", "Taxonomy")] <- NULL
      # Subset real samples
      ID_SAMPLE <- 
        reactives$metadata[which(reactives$metadata[, "Sample_type"] == 
                                   sample_types_allowed[1]), "Sample_ID"]
      reactives$featuredata_1 <- reactives$featuredata[, ID_SAMPLE]
      reactives$metadata_1 <- 
        reactives$metadata[reactives$metadata$Sample_ID %in% ID_SAMPLE, ]
      # Subset NEG1 control samples
      # ADD CHECK IF THERE ARE NEG1 SAMPLES!
      ID_NEG1 <- 
        reactives$metadata[which(reactives$metadata[, "Sample_type"] == 
                                   sample_types_allowed[4]), "Sample_ID"]
      reactives$featuredata_neg1 <- reactives$featuredata[, ID_NEG1]
      # Subset NEG2 control samples 
      # ADD CHECK IF THERE ARE NEG2 SAMPLES!
      ID_NEG2 <- 
        reactives$metadata[which(reactives$metadata[, "Sample_type"] == 
                                   sample_types_allowed[5]), "Sample_ID"]
      reactives$featuredata_neg2 <- reactives$featuredata[, ID_NEG2]
      
      # ------------------------------------------------------------------------
      # Proceed one step and start the filtering 
      # ------------------------------------------------------------------------
      reactives$step_var <- 2
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
      print(colnames(reactives$featuredata_1))
      sample_read_sums <- colSums(reactives$featuredata_1)
      samples_to_keep <- names(
        sample_read_sums[sample_read_sums >= as.numeric(req_reads_per_sample)])
      # Filter featuredata and metadata
      reactives$featuredata_2 <- reactives$featuredata_1[, samples_to_keep]
      reactives$metadata_2 <- 
        reactives$metadata_1[reactives$metadata_1$Sample_ID %in% samples_to_keep, ]
      reactives$featuredata_current <- reactives$featuredata_2
      reactives$metadata_current <- reactives$metadata_2
      print(colnames(reactives$featuredata_2))
    }
    
    # --------------------------------------------------------------------------
    # Filter step 3: remove features by abundance
    # --------------------------------------------------------------------------
    if(reactives$step_var == 3 || reactives$req_reads_per_feature_old != input$req_reads_per_feature) {
      print("filtering_3 branch started")
      # Define feature to keep
      feature_read_sums <- rowSums(reactives$featuredata_2)
      feature_to_keep_abund <- names(
        feature_read_sums[feature_read_sums >= as.numeric(req_reads_per_feature)])
      # Filter featuredata and metadata
      reactives$featuredata_3 <- 
        reactives$featuredata_2[feature_to_keep_abund, ]
      # Remove empty samples:
      sample_read_sums <- colSums(reactives$featuredata_3)
      print(sample_read_sums) #
      samples_to_keep <- names(sample_read_sums[sample_read_sums > 0])
      print(samples_to_keep) #
      reactives$featuredata_3 <- reactives$featuredata_3[, samples_to_keep]
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
        t(decostand(t(reactives$featuredata_3), method = "total"))
      # Define feature to keep
      feature_read_freqs <- apply(featuredata_3_freq, 1, max)
      feature_to_keep_freq <- names(
        feature_read_freqs[feature_read_freqs >= as.numeric(req_ratio_per_feature)])
      # Filter featuredata and metadata
      reactives$featuredata_4 <- 
        reactives$featuredata_3[feature_to_keep_freq, ]
      # Remove empty samples:
      sample_read_sums <- colSums(reactives$featuredata_4)
      print(sample_read_sums) #
      samples_to_keep <- names(sample_read_sums[sample_read_sums > 0])
      print(samples_to_keep) #
      reactives$featuredata_4 <- reactives$featuredata_4[, samples_to_keep]
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
      # Convert current filtered feature table to frequencies
      featuredata_4_freq <- 
        t(decostand(t(reactives$featuredata_4), method = "total"))
      # Calculate sample mean per feature
      sample_mean <- data.frame(sample_mean = rowMeans(featuredata_4_freq))
      # Transform neg1-data to frequencies and calculate mean/span per feature
      # ADD CHECKS FOR NEG1/NEG2 IF THEY EXIST!
      neg1_mean <- 
        rowMeans(t(decostand(t(reactives$featuredata_neg1), method = "total")))
      neg1_span <- 
        apply(reactives$featuredata_neg1, 1, function(x) sum(x > 0)/length(x))
      # Transform neg2-data to frequencies and calculate mean/span per feature
      neg2_mean <-  
        rowMeans(t(decostand(t(reactives$featuredata_neg2), method = "total")))
      neg2_span <- 
        apply(reactives$featuredata_neg2, 1, function(x) sum(x > 0)/length(x))
      # Merge the different pieces of information for each feature
      reactives$filter_basis <- data.frame(neg1_mean = neg1_mean[names(neg1_mean)],
                                           neg1_span = neg1_span[names(neg1_span)],
                                           neg2_mean = neg2_mean[names(neg2_mean)],
                                           neg2_span = neg2_span[names(neg2_span)])
      reactives$filter_basis <- 
        merge(reactives$filter_basis, sample_mean, by = 0, all = TRUE)
      rownames(reactives$filter_basis) <- reactives$filter_basis[, "Row.names"]
      reactives$filter_basis["Row.names"] <- NULL
      # Calculate the ratios
      reactives$filter_basis["ratio_neg1"] <- 
        reactives$filter_basis$neg1_mean / reactives$filter_basis$sample_mean
      reactives$filter_basis["ratio_neg2"] <-
        reactives$filter_basis$neg2_mean / reactives$filter_basis$sample_mean
      # Apply filter criteria and return Sample IDs that should be removed
      if(as.numeric(input$req_ratio_neg1) == Inf && 
         as.numeric(input$req_span_neg1) != 0.0001) {
        feature_removed_neg1 <- reactives$filter_basis %>%
          filter(neg1_span >= as.numeric(input$req_span_neg1)) %>%
          rownames()
      } else {
        feature_removed_neg1 <- reactives$filter_basis %>%
          filter(neg1_span >= as.numeric(input$req_span_neg1)) %>%
          filter(ratio_neg1 > as.numeric(input$req_ratio_neg1)) %>%
          rownames()
      }
      if(as.numeric(input$req_ratio_neg2) == Inf && 
         as.numeric(input$req_span_neg2) != 0.0001) {
        feature_removed_neg2 <- reactives$filter_basis %>%
          filter(neg2_span >= as.numeric(input$req_span_neg2)) %>%
          rownames()
      } else {
        feature_removed_neg2 <- reactives$filter_basis %>%
          filter(neg2_span >= as.numeric(input$req_span_neg2)) %>%
          filter(ratio_neg2 > as.numeric(input$req_ratio_neg2)) %>%
          rownames()
      }
      # Filter featuredata and metadata
      reactives$featuredata_5 <- 
        reactives$featuredata_4[!rownames(reactives$featuredata_4) %in% unique(
          c(feature_removed_neg1, feature_removed_neg2)), ]
      # Remove empty samples
      sample_read_sums <- colSums(reactives$featuredata_5)
      print(sample_read_sums) #
      samples_to_keep <- names(sample_read_sums[sample_read_sums > 0])
      print(samples_to_keep) #
      reactives$featuredata_5 <- reactives$featuredata_5[, samples_to_keep]
      reactives$metadata_5 <- 
        reactives$metadata_4[reactives$metadata_4$Sample_ID %in% 
                               colnames(reactives$featuredata_5), ]
      reactives$featuredata_current <- reactives$featuredata_5
      reactives$metadata_current <- reactives$metadata_5
    }
    
    # --------------------------------------------------------------------------
    # Step 6: save final files
    # --------------------------------------------------------------------------
    if(reactives$step_var == 6) {
      # Save final featuretable with absolute counts
      featuredata_current <- merge(reactives$featuredata_current,
                                   reactives$featuredata_taxonomy, 
                                   by = 0, all.x = TRUE)
      colnames(featuredata_current)[1] <- "OTU_ID"
      write.table(
        featuredata_current, 
        file = paste0(reactives$output_dir, "/output/Featuretable_final.txt"),
        row.names = FALSE, sep = "\t", quote = FALSE)
      # Save final featuretable with relative abundance
      featuredata_current_rel <- 
        t(decostand(t(reactives$featuredata_current), method = "total"))
      featuredata_current_rel <- merge(featuredata_current_rel,
                                       reactives$featuredata_taxonomy, 
                                       by = 0, all.x = TRUE)
      colnames(featuredata_current_rel)[1] <- "OTU_ID"
      write.table(
        featuredata_current_rel, 
        file = paste0(reactives$output_dir, "/output/Featuretable_final_frequency.txt"),
        row.names = FALSE, sep = "\t", quote = FALSE)
      # Save final metatable
      write.table(
        data.frame(Sample_ID = rownames(reactives$metadata_current), 
                   reactives$metadata_current), 
        file = paste0(reactives$output_dir, "/output/Metatable_final.txt"),
        row.names = FALSE, sep = "\t", quote = FALSE)
      # Save contamination filter basis
      write.table(
        data.frame(OTU_ID = rownames(reactives$filter_basis), 
                   reactives$filter_basis), 
        file = paste0(reactives$output_dir, "/interim/Contamination_filter_basis.txt"),
        row.names = FALSE, sep = "\t", quote = FALSE)
      # Save reduction of reads per OTU by filter step
      Reduction_per_feature <- 
        Reduce(function(x, y) merge(x = x, y = y, by.x = "Row.names", by.y = 0, 
                                    all.x = TRUE), 
               list(data.frame(Row.names = names(rowSums(reactives$featuredata)),
                               Original = rowSums(reactives$featuredata)), 
                    data.frame(Step1 = rowSums(reactives$featuredata_1)), 
                    data.frame(Step2 = rowSums(reactives$featuredata_2)),
                    data.frame(Step3 = rowSums(reactives$featuredata_3)),
                    data.frame(Step4 = rowSums(reactives$featuredata_4)),
                    data.frame(Step5 = rowSums(reactives$featuredata_5))))
      colnames(Reduction_per_feature)[1] <- "OTU_ID"
      write.table(
        Reduction_per_feature, 
        file = paste0(reactives$output_dir, "/interim/Read_reduction_per_feature.txt"),
        row.names = FALSE, sep = "\t", quote = FALSE)
      # Save reduction of reads per species by filter step 
      Reduction_per_taxonomy <- 
        merge(Reduction_per_feature, reactives$featuredata_taxonomy,
              by.x = "OTU_ID", by.y = 0, all = TRUE)
      Reduction_per_taxonomy <- Reduction_per_taxonomy %>% 
        group_by(Taxonomy) %>%
        summarise_at(vars(Original, Step1, Step2, Step3, Step4, Step5),
                     sum, na.rm = TRUE)
      write.table(
        Reduction_per_taxonomy, 
        file = paste0(reactives$output_dir, "/interim/Read_reduction_per_taxonomy.txt"),
        row.names = FALSE, sep = "\t", quote = FALSE)
      # Save reduction of reads per metatable category by filter step
      
      # ------------------------------------------------------------------------
      # Create overview file
      # ------------------------------------------------------------------------
      writeLines(paste0(c(
        "MicrobIEM - quality control and analysis tool for microbiome data",
        "-----------------------------------------------------------------",
        "", "Input meta table:", input$metafile$name,
        "", "Input feature table:", input$featurefile$name,
        "", "Minimum reads per sample:", input$req_reads_per_sample,
        "", "Minimum reads per feature:", input$req_reads_per_feature,
        "", "Minimum relative frequency per feature:", input$req_ratio_per_feature,
        "", "Frequency mean ratio (NEG1 by SAMPLE):", 
        names(neg_ratio_steps[neg_ratio_steps == input$req_ratio_neg1]),
        "", "Span threshold (NEG1):", 
        names(neg_span_steps[neg_span_steps == input$req_span_neg1]),
        "", "Frequency mean ratio (NEG2 by SAMPLE):", 
        names(neg_ratio_steps[neg_ratio_steps == input$req_ratio_neg2]),
        "", "Span threshold (NEG2)", 
        names(neg_span_steps[neg_span_steps == input$req_span_neg2])), 
        sep = " "), paste0(reactives$output_dir, "/output/Filter_criteria.txt"))
    }
    # --------------------------------------------------------------------------
    # Pre-calculate values for analysis
    # --------------------------------------------------------------------------
    # Define alpha diversity indices
    Richness.function <- function(x) {sum(x > 0)}
    Shannon.function <- function(x) {
      -sum(scale(x, center = FALSE, scale = sum(x))*
             log(scale(x, center = FALSE, scale = sum(x))), na.rm = TRUE)
    }
    InvSimpson.function <- function(x) {
      1/(sum((scale(x, center = FALSE, scale = sum(x)))^2))
    }
    Simpson.function <- function(x) {
      sum((scale(x, center = FALSE, scale = sum(x)))^2)
    }
    Evenness.function <- function(x) {
      (-sum(scale(x, center = FALSE, scale = sum(x)) *
              log(scale(x, center = FALSE, scale = sum(x))), na.rm = TRUE)) /
        log(sum(x > 0))
    }
    # Calculate alpha diversity indices for current metadata
    reactives$alpha_diversity_values <- data.frame(
      Richness = apply(reactives$featuredata_current, 2, Richness.function),
      Shannon = apply(reactives$featuredata_current, 2, Shannon.function),
      Inv.Simpson = apply(reactives$featuredata_current, 2, InvSimpson.function),
      Simpson = apply(reactives$featuredata_current, 2, Simpson.function),
      Evenness = apply(reactives$featuredata_current, 2, Evenness.function))
    rownames(reactives$alpha_diversity_values) <- 
      colnames(reactives$featuredata_current)
    
    # ------------------------------------------------------------------------
    # Provide information for visualisation and build filtering plots
    # ------------------------------------------------------------------------
    if(input$visualization_type == "Correlation of reads and features") {
      data_to_plot <- data.frame(
        reads = colSums(reactives$featuredata_current),
        features = apply(reactives$featuredata_current, 2, function(x) sum(x > 0)))
      print(str(data_to_plot))
      output$plot <- renderPlot({
        ggplot(data = data_to_plot, aes(x = reads, y = features)) +
          geom_point()
      }) 
    }
    if(input$visualization_type == "Change in feature abundance") {
      data_to_plot <- data.frame(
        reads_before = rowSums(reactives$featuredata_1))
      data_to_plot <- merge(data_to_plot, data.frame(
        reads_now = rowSums(reactives$featuredata_current)), by = 0, all = TRUE)
      data_to_plot[is.na(data_to_plot)] <- 0
      print(str(data_to_plot))
      output$plot <- renderPlot({
        ggplot(data = data_to_plot, aes(x = reads_before, y = reads_now)) +
          geom_point()
      }) 
    }
    if(input$visualization_type == "Contamination removal - NEG1" ||
       input$visualization_type == "Contamination removal - NEG2") {
      if(reactives$step_var != 5) {
        showModal(modalDialog(title = "Error", "Please wait for contaminant
                              removal analysis."))
        updateSelectInput(session, inputId = "visualization_type",
                          selected = "Correlation of reads and features")
      } else {
        if(input$visualization_type == "Contamination removal - NEG1") {
          contamination_plot <- 
            ggplot(data = reactives$filter_basis, 
                   aes(x = ratio_neg1, y = neg1_span, size = sample_mean)) +
            geom_point() +
            scale_x_log10()
          if(as.numeric(input$req_span_neg1) != 0.0001) {
            contamination_plot <- contamination_plot +
              geom_hline(aes(alpha = "Span_threshold", 
                              yintercept = as.numeric(input$req_span_neg1))) +
              scale_alpha_manual(values = 1)
          }
          if(as.numeric(input$req_ratio_neg1) != Inf) {
            contamination_plot <- contamination_plot +
              geom_vline(aes(linetype = "Ratio_threshold",
                             xintercept = as.numeric(input$req_ratio_neg1)))
          }
          output$plot <- renderPlot({
            contamination_plot
          })
        }
        if(input$visualization_type == "Contamination removal - NEG2") {
          output$plot <- renderPlot({
            ggplot(data = reactives$filter_basis, 
                   aes(x = ratio_neg2, y = neg2_span, size = sample_mean)) +
              geom_point() +
              scale_x_log10()
          })
        }
      }
    }
    if(input$visualization_type == "Reduction of total reads") {
      sum_of_reads_function <- function(x) {
        if(is.na(x)) {return(NA)} else {return(sum(colSums(x)))}
      }
      data_to_plot <- data.frame(
        step = c("0 - original", "1 - without controls", "2 - sample filter",
                 "3 - feature abundance filter", "4 - feature frequency filter",
                 "5 - contamination filter"),
        sum_of_reads = c(sum_of_reads_function(reactives$featuredata),
                         sum_of_reads_function(reactives$featuredata_1),
                         sum_of_reads_function(reactives$featuredata_2),
                         sum_of_reads_function(reactives$featuredata_3),
                         sum_of_reads_function(reactives$featuredata_4),
                         sum_of_reads_function(reactives$featuredata_5)))
      print(str(data_to_plot))
      output$plot <- renderPlot({
        ggplot(data = data_to_plot, aes(x = step, y = sum_of_reads)) +
          geom_bar(stat = "identity")
      }) 
    }
    print("ROWNAMES")
    print(rownames(reactives$featuredata_current))
    print(rownames(reactives$metadata_current))
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
      shinyjs::hide("placeholder") # NEWLINE
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
      # NEWLINE :
      shinyjs::show("placeholder") # NEWLINE
      insertUI(selector = "#placeholder",
               ui = selectInput(inputId = "test_ID",
                                label = "Some value",		
                                choices = seq(0:ncol(reactives$featuredata_neg1))))
      
      shinyjs::show("header_neg1") 
      shinyjs::show("req_ratio_neg1")
      shinyjs::show("req_span_neg1")
      shinyjs::show("header_neg2")
      shinyjs::show("req_ratio_neg2")
      shinyjs::show("req_span_neg2")
    }
    if(reactives$step_var == 6) {
      appendTab(
        inputId = "tabset",
        tabPanel(
          title = "Alpha diversity",
          br(),
          selectInput(inputId = "metavar_alpha", 
                      label = "Metainformation variable",
                      choices = colnames(reactives$metadata_current)),
          selectInput(inputId = "subvar_alpha",
                      label = "Variable for subgroup analysis",
                      choices = c("ignore", colnames(reactives$metadata_current))),
          radioButtons(inputId = "scaling", 
                       label = "Y axis scaling",
                       choices = c("Z-normalized values", "Raw values")),
          actionButton(inputId = "alpha_analysis", 
                       label = "Update plot")
        )
      )
      appendTab(
        inputId = "tabset",
        tabPanel(
          title = "Beta diversity",
          br(),
          selectInput(inputId = "metavar_beta", 
                      label = "Metainformation variable",
                      choices = colnames(reactives$metadata_current)),
          radioButtons(inputId = "plot_beta", 
                       label = "Visualization type", 
                       choices = c(
                         "Principal coordinates analysis (PCoA)" = "beta_pcoa",
                         "Non-metric multidimensional scaling (nMDS)" = "beta_nmds")),
          actionButton(inputId = "beta_analysis", 
                       label = "Update plot")
        )
      )
      appendTab(
        inputId = "tabset",
        tabPanel(
          title = "Taxonomy",
          br(),
          selectInput(inputId = "metavar_taxonomy", 
                      label = "Metainformation variable",
                      choices = colnames(reactives$metadata_current)),
          radioButtons(inputId = "taxonomy_level", 
                       label = "Taxonomic level", 
                       choices = c("Domain", "Phylum", "Class", "Order", 
                                   "Family", "Genus", "Species")),
          actionButton(inputId = "tax_analysis", 
                       label = "Update plot")
        )
      )
      appendTab(
        inputId = "tabset",
        tabPanel(
          title = "Sample selection",
          br(),
          lapply(colnames(reactives$metadata_current), function(x) {
            pickerInput(inputId = paste0("sampleselector_", x),
                        label = x,
                        choices = unique(reactives$metadata_current[, x]),
                        selected = unique(reactives$metadata_current[, x]),
                        multiple = TRUE,
                        options = list(`actions-box` = TRUE,
                                       `selected-text-format` = "count",
                                       `count-selected-text` = "{0}/{1} selected"))
          }) # Close lapply function
        )
      )
      updateTabsetPanel(session, "tabset", selected = "Alpha diversity")
    }    
  }
  
  # ----------------------------------------------------------------------------
  # Prevent the user from messing up the files
  # ----------------------------------------------------------------------------
  observeEvent(input$tabset, {
    if(reactives$step_var == 6 && input$tabset == "Filtering") {
      showModal(modalDialog(
        title="Are you sure?",
        "Returning to filtering may overwrite your final filtered files",
        footer = tagList(actionButton("refilter_button", "Go back to filtering"),
                         modalButton("Stay on analysis site"))
      ))
      updateTabsetPanel(session, "tabset", selected = "Alpha diversity")
    }
    if(reactives$step_var == 5 && input$tabset != "Filtering") {
      showModal(modalDialog(
        title="Are you sure?",
        "Returning to analysis will not save the changes you made in filtering.
        If you want to save the changes and overwrite your final filtered files,
        press the 'Next' button at the bottom of the page.",
        footer = tagList(actionButton("reanalysis_button", 
                                      "Proceed to analysis without saving changes"),
                         modalButton("Stay on filtering site")
        )
      ))
      updateTabsetPanel(session, "tabset", selected = "Filtering")
    }
  })
  
  observeEvent(input$refilter_button, {
    reactives$step_var <- reactives$step_var - 1
    updateTabsetPanel(session, "tabset", selected = "Filtering")
    removeModal()
  })
  
  observeEvent(input$reanalysis_button, {
    reactives$step_var <- reactives$step_var + 1
    updateTabsetPanel(session, "tabset", selected = "Alpha diversity") 
    removeModal()
  })
  
  # ----------------------------------------------------------------------------
  # Alpha diversity analysis
  # ----------------------------------------------------------------------------
  observeEvent(input$alpha_analysis, {
    print(paste0("INFO|server::alpha_diversity",Sys.time()))
    # Filter feature and metadata according to sample selection
    samples_deselected <- Sample_selection_check()
    Metadata <- 
      reactives$metadata_current[!rownames(reactives$metadata_current) %in% 
                                   samples_deselected, ]
    # Select alpha diversity values for selected samples
    alpha_diversity_result <- 
      reactives$alpha_diversity_values[!rownames(reactives$metadata_current) %in% 
                                         samples_deselected ,]
    # Perform z-transformation if needed
    print(str(alpha_diversity_result))
    if (input$scaling == "Z-normalized values") {
      alpha_diversity_result <- 
        as.data.frame(apply(alpha_diversity_result, 2, function(x)
          scale(x, center = TRUE, scale = TRUE)), 
          row.names = rownames(alpha_diversity_result))
    }
    # Attach variables of interest for alpha diversity
    alpha_diversity_result["Group"] <- Metadata[paste0(input$metavar_alpha)]
    if(input$subvar_alpha != "ignore") {
      alpha_diversity_result["Subgroup"] <- Metadata[paste0(input$subvar_alpha)]
    }
    print(str(melt(alpha_diversity_result)))
    # Generate the alpha diversity plot
    alpha_diversity_result_plot <- melt(alpha_diversity_result)
    alpha_diversity_plot <- 
      ggplot(data = alpha_diversity_result_plot, 
             aes(x = Group, y = value, colour = Group)) +
      geom_boxplot()
    if (input$subvar_alpha != "ignore" && input$scaling == "Z-normalized values") {
      alpha_diversity_plot <- alpha_diversity_plot + 
        facet_wrap(Subgroup ~ variable, 
                   nrow = length(unique(alpha_diversity_result_plot$Subgroup)),
                   scales = "free_x")
    } else if (input$subvar_alpha != "ignore") {
      alpha_diversity_plot <- alpha_diversity_plot +
        facet_wrap(Subgroup ~ variable, 
                   ncol = length(levels(alpha_diversity_result_plot$variable)), 
                   scales = "free")      
    } else if (input$scaling == "Z-normalized values") {
      alpha_diversity_plot <- alpha_diversity_plot +
        facet_wrap(. ~ variable, 
                   ncol = length(levels(alpha_diversity_result_plot$variable)))
    } else {
      alpha_diversity_plot <- alpha_diversity_plot +
        facet_wrap(. ~ variable, 
                   ncol = length(levels(alpha_diversity_result_plot$variable)), 
                   scales = "free_y")
    }
    
    output$plot <- renderPlot({
      alpha_diversity_plot
    })
    # Output table and significance
  })
  
  observeEvent(input$subvar_alpha, {
    if(input$subvar_alpha == input$metavar_alpha) {
      showModal(modalDialog(
        title = "Error10", "Please choose a second variable for alpha diversity
        that is different from your main variable."))
      updateSelectInput(session, inputId = "subvar_alpha", selected = "ignore")
    }
  })

  Sample_selection_check <- function() {
    # return all sample names that should not appear in the analysis
    samples_deselected <- rownames(reactives$metadata_current)[sort(
      unique(unlist(sapply(colnames(reactives$metadata_current), function(x) {
        which(!reactives$metadata_current[, x] %in% 
                input[[paste0("sampleselector_", x)]])
        }
      )))
    )]
    print("SAMPLE SELECTION CHECK")
    print(samples_deselected)
    return(samples_deselected)
  }
  
  # ----------------------------------------------------------------------------
  # Beta diversity analysis
  # ----------------------------------------------------------------------------
  observeEvent(input$beta_analysis, {
    print(paste0("INFO|server::beta_diversity",Sys.time()))
    
  })
  
  # ----------------------------------------------------------------------------
  # Taxonomy analysis
  # ----------------------------------------------------------------------------
  observeEvent(input$tax_analysis, {
    print(paste0("INFO|server::taxnonomy",Sys.time()))
    
  })
  
  #output$plot <- renderPlot({
  #  set.seed(1)  
  #  plot(rnorm(10, 0, 2), rnorm(10, as.numeric(input$req_reads_per_sample), 
  #                              isolate(reactives$step_var)))
  #})
  
}
