################################################################################
#
# This is the server logic of MicrobIEM
#
################################################################################

# ------------------------------------------------------------------------------
# Install and load required packages 
# ------------------------------------------------------------------------------
# Install packages
packages_server <- c("shinyjs", "DT", "plotly", "shinyWidgets",
                     "dplyr", "reshape2", "ggplot2", "vegan")
install.packages(setdiff(packages_server, rownames(installed.packages())))
# Load packages
library(shinyjs)
library(ggplot2)
library(vegan)
library(dplyr)
library(shinyWidgets)
library(reshape2)
library(plotly)
library(DT)

# ------------------------------------------------------------------------------
# Define re-used parameters
# ------------------------------------------------------------------------------
sample_types_allowed <- c("SAMPLE", "POS1", "POS2", "NEG1", "NEG2")
neg_ratio_steps <- c("ignore" = Inf, 
                     "2" = 2,
                     "1.5" = 1.5,
                     "1" = 1,
                     "0.5" = 0.5,
                     "0.1" = 0.1)
plot_theme <- theme(
  strip.background = element_rect(fill = "#f5f5f5", colour = "grey50"),
  panel.background = element_rect(fill = NA, colour = "grey50"),
  panel.border = element_rect(fill = NA, colour = "grey50"),
  panel.grid = element_blank(),
  legend.key = element_blank()) 
set.seed(1)
theme_colours <- c("#4669A5", "#e65d1e", "#5dc9cf", "#ffda0a", "#ba002b", 
                   "#6ed957", "#ffbdc7", "#c97b00", "#fff79c", "#277527", 
                   "#ff7891", "#c6f279", "#3f1163", "#f0b400", "#729ce8", 
                   "#8F1870", "#e3cb6f", "#248c8c", "#e61e1e", "#98d8fe",
                   "#a84c1b", "#D68CCC", "#065e5e", "#ff9c75", "#9DBD32", 
                   "#B957A3", sample(rainbow(250)))
contamination_neg1_colours <- c("kept" = "#2fa4e7", # blue
                                "removed (NEG1)" = "#d05517", # orange
                                "removed (NEG2)" = "#757575", # dark grey
                                "removed (NEG1 & NEG2)" = "#f1aa87") # light orange
contamination_neg2_colours <- c("kept" = "#2fa4e7", # blue
                                "removed (NEG1)" = "#757575", # dark grey
                                "removed (NEG2)" = "#d05517", # orange
                                "removed (NEG1 & NEG2)" = "#f1aa87") # light orange

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
                              req_ratio_neg2_old = Inf, req_span_neg2_old = 0.0001,
                              neg1_span_steps = NA, req_span_neg1 = NA,
                              neg2_span_steps = NA, req_span_neg2 = NA,
                              analysis_started = FALSE,
                              beta_diversity_analysis = FALSE)
  
  # Hide buttons and input fields defined in UI function
  shinyjs::hide("visualization_type")
  shinyjs::hide("req_reads_per_sample")
  shinyjs::hide("req_reads_per_feature")
  shinyjs::hide("req_ratio_per_feature")
  shinyjs::hide("header_neg1")
  shinyjs::hide("req_ratio_neg1")
  shinyjs::hide("header_neg2")
  shinyjs::hide("req_ratio_neg2")
  shinyjs::hide("update_button")
  shinyjs::hide("back_button")
  shinyjs::hide("next_button")

  # ----------------------------------------------------------------------------
  # Read and check the metafile
  # ----------------------------------------------------------------------------
  observeEvent(input$metafile, {
    print(paste0("INFO - input in metafile - ", Sys.time()))
    # Check correct file extension
    if (!(tools::file_ext(input$metafile$datapath) == "txt")) {
      showModal(modalDialog(
        title = "Error 1a", "Please choose a tab-separated txt file as meta table."))
    } else {
      reactives$metadata <- read.csv(input$metafile$datapath, sep = "\t", 
                                     header = TRUE, check.names = FALSE)
      reactives$metadata[] <- lapply(reactives$metadata, as.character)
      colnames_md <- colnames(reactives$metadata)
      # Check for columns Sample_ID and Sample_type
      if (colnames_md[1] != "Sample_ID" || !("Sample_type" %in% colnames_md)) {
        showModal(modalDialog(
          title = "Error 2", 
          "Please provide a meta information table with the first column 
          'Sample_ID' and a column 'Sample_type' to define samples and controls."))
        reactives$metadata <- NA
      } else {
        sample_types_observed <- unique(reactives$metadata[, "Sample_type"])
        # Check for correct sample type definitions
        if ((sum(sample_types_observed %in% sample_types_allowed) != 
             length(sample_types_observed))) {
          showModal(modalDialog(
            title = "Error 3", 
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
    print(paste0("INFO - input in featurefile - ", Sys.time()))
    # Check correct file extension
    if (!(tools::file_ext(input$featurefile$datapath) == "txt")) {
      showModal(modalDialog(
        title = "Error 1b", "Please choose a tab-separated txt file as feature table."))
    } else {
      reactives$featuredata <- read.csv(input$featurefile$datapath, sep = "\t", 
                                        header = TRUE, check.names = FALSE)
      colnames_fd <- colnames(reactives$featuredata)
      # Check for columns OTU_ID and Taxonomy
      if(length(colnames_fd) < 3 || colnames_fd[1] != "OTU_ID" || 
         colnames_fd[length(colnames_fd)] != "Taxonomy") {
        showModal(modalDialog(
          title = "Error 4", "Please provide a feature table with the first 
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
    print(paste0("INFO - file_open_success started - ", Sys.time()))
    sample_names_md <- reactives$metadata[, "Sample_ID"]
    sample_names_fd <- head(colnames(reactives$featuredata)[-1], -1)
    # If sample names match in feature file and meta file, start filtering
    if (identical(sort(sample_names_md), sort(sample_names_fd))) {
      print(paste0("INFO - meta and featurefile loaded - ", Sys.time()))
      print(paste0("INFO - ", length(sample_names_md), 
                   " samples and controls found in this analysis"))
      
      # ------------------------------------------------------------------------
      # Prepare the data for filtering 
      # ------------------------------------------------------------------------
      # Replace NA values in the metafile
      if (sum(is.na(reactives$metadata)) > 0) {
        print(paste0("INFO - replacing NA values in metafile - ", Sys.time()))
        reactives$metadata[is.na(reactives$metadata)] <- "n.a." 
        print(c("INFO - Replaced ", sum(is.na(reactives$metadata)), 
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
      print(paste0("INFO - filtering step 1 - separate samples and controls - ", 
                   Sys.time()))
      ID_SAMPLE <- 
        reactives$metadata[which(reactives$metadata[, "Sample_type"] == 
                                   sample_types_allowed[1]), "Sample_ID"]
      reactives$featuredata_1 <- reactives$featuredata[, ID_SAMPLE]
      reactives$metadata_1 <- 
        reactives$metadata[reactives$metadata$Sample_ID %in% ID_SAMPLE, ]
      # Remove empty features
      feature_read_sums <- rowSums(reactives$featuredata_1)
      features_to_keep <- names(feature_read_sums[feature_read_sums > 0])
      reactives$featuredata_1 <- reactives$featuredata_1[features_to_keep, ]
      # Subset NEG1 control samples
      ID_NEG1 <- 
        reactives$metadata[which(reactives$metadata[, "Sample_type"] == 
                                   sample_types_allowed[4]), "Sample_ID"]
      reactives$featuredata_neg1 <- reactives$featuredata[, ID_NEG1, drop = FALSE]
      # Create span UI values and default value for NEG1
      neg1_span_values <- seq(1, ncol(reactives$featuredata_neg1))
      reactives$neg1_span_steps <- setNames(
        c(0.0001, neg1_span_values/length(neg1_span_values)), 
        nm = c("ignore", paste0(
          round(neg1_span_values/length(neg1_span_values), 2), " (",
          neg1_span_values, "/", length(neg1_span_values), ")")))
      reactives$req_span_neg1 <- reactives$neg1_span_steps[1]
      reactives$req_span_neg1_old <- reactives$req_span_neg1
      # Subset NEG2 control samples 
      ID_NEG2 <- 
        reactives$metadata[which(reactives$metadata[, "Sample_type"] == 
                                   sample_types_allowed[5]), "Sample_ID"]
      reactives$featuredata_neg2 <- reactives$featuredata[, ID_NEG2, drop = FALSE]
      print(paste0("INFO - detected ", length(ID_SAMPLE), " real samples, ",
                   length(ID_NEG1), " NEG1 controls, and ",
                   length(ID_NEG2), " NEG2 controls in the data."))
      # Create span UI values and default value for NEG2
      neg2_span_values <- seq(1, ncol(reactives$featuredata_neg2))
      reactives$neg2_span_steps <- setNames(
        c(0.0001, neg2_span_values/length(neg2_span_values)), 
        nm = c("ignore", paste0(
          round(neg2_span_values/length(neg2_span_values), 2), " (",
          neg2_span_values, "/", length(neg2_span_values), ")")))
      reactives$req_span_neg2 <- reactives$neg2_span_steps[1]
      reactives$req_span_neg2_old <- reactives$req_span_neg2
      
      # ------------------------------------------------------------------------
      # Pre-calculate values for taxonomy analysis - split by taxonomic levels
      # ------------------------------------------------------------------------
      print(paste0("INFO - split taxonomy by taxonomic levels - ", Sys.time()))
      # Split taxonomy by levels
      reactives$taxonomy_data <- 
        apply(reactives$featuredata_taxonomy, 1, function(x)
          unlist(strsplit(x, split = ";|,")))
      # Iterate over vector items, apply the unique colname to them, rbind all
      # Reactives$taxonomy_data is stored as a list if missing entries exist
      if(is.list(reactives$taxonomy_data)) { 
        reactives$taxonomy_data <- data.frame(
          do.call(rbind, lapply(reactives$taxonomy_data, "[", 
                                unique(unlist(c(sapply(reactives$taxonomy_data, 
                                                       names)))))))
      } else { # Or as matrix if only complete taxonomies exits (e.g. "NoGenus")
        reactives$taxonomy_data <- data.frame(t(reactives$taxonomy_data))
      }
      # Prevent a bug when there is only one taxonomic level in some annotations
      if(any(colnames(reactives$taxonomy_data) == "Taxonomy") && 
         any(colnames(reactives$taxonomy_data) == "NA.")) {
        reactives$taxonomy_data$Taxonomy[is.na(reactives$taxonomy_data$Taxonomy)] <- 
          reactives$taxonomy_data$NA.[!is.na(reactives$taxonomy_data$NA.)]
        reactives$taxonomy_data$NA. <- NULL
      }
      colnames(reactives$taxonomy_data) <- 
        c("Domain", "Phylum", "Class", "Order", "Family", "Genus", 
          "Species")[1:length(colnames(reactives$taxonomy_data))]
      # Replace NA and empty values
      reactives$taxonomy_data[is.na(reactives$taxonomy_data)] <- "No_annotation"
      reactives$taxonomy_data[reactives$taxonomy_data == ""] <- "No_annotation"
      
      # ------------------------------------------------------------------------
      # Proceed one step and start the filtering 
      # ------------------------------------------------------------------------
      reactives$step_var <- 2
      filter_feature_table()
      # Make buttons appear
      step_var_UIchange()
      shinyjs::show("update_button")
      shinyjs::show("back_button")
      shinyjs::show("next_button")
    } else {
      showModal(modalDialog(
        title = "Error 5", "Sample IDs do not match between meta table and 
        feature table."))
    }
  }

  # ----------------------------------------------------------------------------
  # Define the filter_feature_table function
  # ----------------------------------------------------------------------------
  filter_feature_table <- function() {

    # --------------------------------------------------------------------------
    # Filter step 2: remove samples
    # --------------------------------------------------------------------------
    if(reactives$step_var == 2 || 
       reactives$req_reads_per_sample_old != input$req_reads_per_sample) {
      print(paste0("INFO - filtering step 2 - remove samples - ", 
                   Sys.time()))
      # Define samples to keep
      sample_read_sums <- colSums(reactives$featuredata_1)
      samples_to_keep <- names(
        sample_read_sums[sample_read_sums >= as.numeric(input$req_reads_per_sample)])
      # Filter featuredata and metadata
      reactives$featuredata_2 <- reactives$featuredata_1[, samples_to_keep]
      reactives$metadata_2 <- 
        reactives$metadata_1[reactives$metadata_1$Sample_ID %in% samples_to_keep, ]
      reactives$featuredata_current <- reactives$featuredata_2
      reactives$metadata_current <- reactives$metadata_2
      print(paste0("INFO - ", nrow(reactives$metadata_current), " samples and ",
                   sum(rowSums(reactives$featuredata_current)), " reads in ", 
                   nrow(reactives$featuredata_current), " features remaining"))
    }
    
    # --------------------------------------------------------------------------
    # Filter step 3: remove features by abundance
    # --------------------------------------------------------------------------
    if(reactives$step_var == 3 || reactives$req_reads_per_feature_old != input$req_reads_per_feature) {
      print(paste0("INFO - filtering step 3 - remove features by abundance - ", 
                   Sys.time()))
      # Define feature to keep
      feature_read_sums <- rowSums(reactives$featuredata_2)
      feature_to_keep_abund <- names(
        feature_read_sums[feature_read_sums >= as.numeric(input$req_reads_per_feature)])
      # Filter featuredata and metadata
      reactives$featuredata_3 <- 
        reactives$featuredata_2[feature_to_keep_abund, ]
      # Remove empty samples:
      sample_read_sums <- colSums(reactives$featuredata_3)
      samples_to_keep <- names(sample_read_sums[sample_read_sums > 0])
      reactives$featuredata_3 <- reactives$featuredata_3[, samples_to_keep]
      reactives$metadata_3 <- 
        reactives$metadata_2[reactives$metadata_2$Sample_ID %in% 
                               colnames(reactives$featuredata_3), ]
      reactives$featuredata_current <- reactives$featuredata_3
      reactives$metadata_current <- reactives$metadata_3      
      print(paste0("INFO - ", nrow(reactives$metadata_current), " samples and ",
                   sum(rowSums(reactives$featuredata_current)), " reads in ", 
                   nrow(reactives$featuredata_current), " features remaining"))
    }
    
    # --------------------------------------------------------------------------
    # Filter step 4: remove features by frequency
    # --------------------------------------------------------------------------
    if(reactives$step_var == 4 ||
       reactives$req_ratio_per_feature_old != input$req_ratio_per_feature) {
      print(paste0("INFO - filtering step 4 - remove features by frequency - ", 
                   Sys.time()))
      # Convert feature table to frequencies
      featuredata_3_freq <- 
        t(decostand(t(reactives$featuredata_3), method = "total"))
      # Define feature to keep
      feature_read_freqs <- apply(featuredata_3_freq, 1, max)
      feature_to_keep_freq <- names(
        feature_read_freqs[feature_read_freqs >= as.numeric(input$req_ratio_per_feature)])
      # Filter featuredata and metadata
      reactives$featuredata_4 <- 
        reactives$featuredata_3[feature_to_keep_freq, ]
      # Remove empty samples:
      sample_read_sums <- colSums(reactives$featuredata_4)
      samples_to_keep <- names(sample_read_sums[sample_read_sums > 0])
      reactives$featuredata_4 <- reactives$featuredata_4[, samples_to_keep]
      reactives$metadata_4 <- 
        reactives$metadata_3[reactives$metadata_3$Sample_ID %in% 
                               colnames(reactives$featuredata_4), ]
      reactives$featuredata_current <- reactives$featuredata_4
      reactives$metadata_current <- reactives$metadata_4
      print(paste0("INFO - ", nrow(reactives$metadata_current), " samples and ",
                   sum(rowSums(reactives$featuredata_current)), " reads in ", 
                   nrow(reactives$featuredata_current), " features remaining"))
    }
    
    # --------------------------------------------------------------------------
    # Filter step 5: remove contaminants
    # --------------------------------------------------------------------------
    if(reactives$step_var == 5 ||
       reactives$req_ratio_neg1_old != input$req_ratio_neg1 ||
       reactives$req_span_neg1_old != reactives$req_span_neg1 ||
       reactives$req_ratio_neg2_old != input$req_ratio_neg2 ||
       reactives$req_span_neg2_old != reactives$req_span_neg2) {
      print(paste0("INFO - filtering step 5 - remove contaminants - ", 
                   Sys.time()))
      # Convert current filtered feature table to frequencies
      featuredata_4_freq <- 
        t(decostand(t(reactives$featuredata_4), method = "total"))
      # Calculate sample mean per feature
      sample_mean <- data.frame(sample_mean = rowMeans(featuredata_4_freq))
      # Transform neg1-data to frequencies and calculate mean/span per feature
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
      reactives$filter_basis <- data.frame(
        neg1_mean = neg1_mean[names(neg1_mean)],
        neg1_span = round(neg1_span[names(neg1_span)], 4),
        neg2_mean = neg2_mean[names(neg2_mean)],
        neg2_span = round(neg2_span[names(neg2_span)], 4))
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
         as.numeric(reactives$req_span_neg1) != 0.0001) {
        feature_removed_neg1 <- reactives$filter_basis %>%
          filter(neg1_span >= round(as.numeric(reactives$req_span_neg1), 4)) %>%
          rownames()
      } else {
        feature_removed_neg1 <- reactives$filter_basis %>%
          filter(neg1_span >= round(as.numeric(reactives$req_span_neg1), 4)) %>%
          filter(ratio_neg1 > as.numeric(input$req_ratio_neg1)) %>%
          rownames()
      }
      if(as.numeric(input$req_ratio_neg2) == Inf && 
         as.numeric(reactives$req_span_neg2) != 0.0001) {
        feature_removed_neg2 <- reactives$filter_basis %>%
          filter(neg2_span >= round(as.numeric(reactives$req_span_neg2), 4)) %>%
          rownames()
      } else {
        feature_removed_neg2 <- reactives$filter_basis %>%
          filter(neg2_span >= round(as.numeric(reactives$req_span_neg2), 4)) %>%
          filter(ratio_neg2 > as.numeric(input$req_ratio_neg2)) %>%
          rownames()
      }
      # Add column for information on removal 
      reactives$filter_basis <- 
        transform(reactives$filter_basis, status = ifelse(
          rownames(reactives$filter_basis) %in% feature_removed_neg1 &
            rownames(reactives$filter_basis) %in% feature_removed_neg2, 
          "removed (NEG1 & NEG2)", ifelse(
            rownames(reactives$filter_basis) %in% feature_removed_neg1,
            "removed (NEG1)", ifelse(
              rownames(reactives$filter_basis) %in% feature_removed_neg2,
              "removed (NEG2)", "kept"))))
      # Attach last taxonomy
      reactives$filter_basis["Taxonomy"] <- 
        apply(reactives$taxonomy_data, 1, function(x) tail(x[x != "No_annotation"], 1))
      # Filter featuredata and metadata
      reactives$featuredata_5 <- 
        reactives$featuredata_4[!rownames(reactives$featuredata_4) %in% unique(
          c(feature_removed_neg1, feature_removed_neg2)), ]
      # Remove empty samples
      sample_read_sums <- colSums(reactives$featuredata_5)
      samples_to_keep <- names(sample_read_sums[sample_read_sums > 0])
      reactives$featuredata_5 <- reactives$featuredata_5[, samples_to_keep]
      reactives$metadata_5 <- 
        reactives$metadata_4[reactives$metadata_4$Sample_ID %in% 
                               colnames(reactives$featuredata_5), ]
      reactives$featuredata_current <- reactives$featuredata_5
      reactives$metadata_current <- reactives$metadata_5
      print(paste0("INFO - ", nrow(reactives$metadata_current), " samples and ",
                   sum(rowSums(reactives$featuredata_current)), " reads in ", 
                   nrow(reactives$featuredata_current), " features remaining"))
    }
    
    # --------------------------------------------------------------------------
    # Step 6: Precalculate values for analysis and create files for download
    # --------------------------------------------------------------------------
    if(reactives$step_var == 6) {
      print(paste0("INFO - filtering step 6 - create final files for analysis - ", 
                   Sys.time()))
      reactives$featuredata_6 <- reactives$featuredata_current
      reactives$metadata_6 <- reactives$metadata_current
      print(paste0("INFO - create final files for download - ", Sys.time()))
      
      # ------------------------------------------------------------------------
      # Create final feature files
      # ------------------------------------------------------------------------
      # Create final featuretable with absolute counts
      reactives$download_featuredata_abs <- merge(
        reactives$featuredata_6, reactives$featuredata_taxonomy, 
        by = 0, all.x = TRUE)
      colnames(reactives$download_featuredata_abs)[1] <- "OTU_ID"

      # Create final featuretable with relative abundance
      reactives$download_featuredata_rel <- 
        t(decostand(t(reactives$featuredata_6), method = "total"))
      reactives$download_featuredata_rel <- merge(
        reactives$download_featuredata_rel, reactives$featuredata_taxonomy, 
        by = 0, all.x = TRUE)
      colnames(reactives$download_featuredata_rel)[1] <- "OTU_ID"

      # ------------------------------------------------------------------------
      # Create overview file for filter settings
      # ------------------------------------------------------------------------
      reactives$download_filtersettings <- paste0(c(
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
        names(reactives$neg1_span_steps[reactives$neg1_span_steps == 
                                          reactives$req_span_neg1]),
        "", "Frequency mean ratio (NEG2 by SAMPLE):", 
        names(neg_ratio_steps[neg_ratio_steps == input$req_ratio_neg2]),
        "", "Span threshold (NEG2):", 
        names(reactives$neg2_span_steps[reactives$neg2_span_steps == 
                                          reactives$req_span_neg2])), 
        sep = " ")
      
      # ------------------------------------------------------------------------
      # Create files for quality control
      # ------------------------------------------------------------------------
      # Create reduction of reads per feature by filter step
      reactives$download_reduction_per_feature <- 
        Reduce(function(x, y) merge(x = x, y = y, by.x = "Row.names", by.y = 0, 
                                    all.x = TRUE), 
               list(data.frame(Row.names = names(rowSums(reactives$featuredata)),
                               Original = rowSums(reactives$featuredata)), 
                    data.frame(Step1 = rowSums(reactives$featuredata_1)), 
                    data.frame(Step2 = rowSums(reactives$featuredata_2)),
                    data.frame(Step3 = rowSums(reactives$featuredata_3)),
                    data.frame(Step4 = rowSums(reactives$featuredata_4)),
                    data.frame(Step5 = rowSums(reactives$featuredata_5))))
      colnames(reactives$download_reduction_per_feature)[1] <- "OTU_ID"
      
      # Create reduction of reads per highest taxonomy by filter step 
      reactives$download_reduction_per_taxonomy <- 
        merge(reactives$download_reduction_per_feature, 
              reactives$featuredata_taxonomy,
              by.x = "OTU_ID", by.y = 0, all = TRUE)
      reactives$download_reduction_per_taxonomy <- 
        reactives$download_reduction_per_taxonomy %>% 
        group_by(Taxonomy) %>%
        summarise_at(vars(Original, Step1, Step2, Step3, Step4, Step5),
                     sum, na.rm = TRUE)
      
      # ------------------------------------------------------------------------
      # Pre-calculate values for alpha and beta diversity analysis
      # ------------------------------------------------------------------------
      print(paste0("INFO - pre-calculate alpha diversity values - ", Sys.time()))
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
      reactives$alpha_diversity_data <- data.frame(
        Richness = apply(reactives$featuredata_6, 2, Richness.function),
        Shannon = apply(reactives$featuredata_6, 2, Shannon.function),
        Inv.Simpson = apply(reactives$featuredata_6, 2, InvSimpson.function),
        Simpson = apply(reactives$featuredata_6, 2, Simpson.function),
        Evenness = apply(reactives$featuredata_6, 2, Evenness.function))
      rownames(reactives$alpha_diversity_data) <- 
        colnames(reactives$featuredata_6)
      
      # Create beta diversity distance matrix of all samples
      print(paste0("INFO - pre-calculate beta diversity values - ", Sys.time()))
      reactives$distance_matrix <- as.data.frame(
        as.matrix(vegdist(decostand(t(reactives$featuredata_6), method = "total"), 
                          method = "bray")))
      
      # ------------------------------------------------------------------------
      # Information for the user for analysis
      # ------------------------------------------------------------------------
      showModal(modalDialog(
        title = "Info", 
        HTML(paste0("There are ", nrow(reactives$metadata_6), " samples and ",
               sum(rowSums(reactives$featuredata_6)), " reads in ", 
               nrow(reactives$featuredata_6), 
               " features remaining for analysis.", 
               "<br><br> Please do not forget to ", 
               HTML("<b>download the final filtered 
               featurefile, metafile, and the filter settings in the 'Download'
               tab</b>."))),
        footer = tagList(modalButton("Ok"))
      ))
    }
    
    # --------------------------------------------------------------------------
    # Provide information for visualisation and build filtering plots
    # --------------------------------------------------------------------------
    if(input$visualization_type == "Correlation of reads and features") {
      correlation_plot <- data.frame(
        Reads = colSums(reactives$featuredata_current),
        Features = apply(reactives$featuredata_current, 2, function(x) sum(x > 0)))
      output$plot <- renderPlotly({
        ggplot(data = correlation_plot, aes(x = Reads, y = Features)) +
          geom_point(aes(text = paste("Sample:", rownames(correlation_plot))),
                     colour = "#2fa4e7") +
          scale_x_log10() +
          plot_theme
      })
    }
    if(input$visualization_type == "Change in feature abundance") {
      data_to_plot <- data.frame(
        reads_before = rowSums(reactives$featuredata_1))
      data_to_plot <- merge(data_to_plot, data.frame(
        reads_now = rowSums(reactives$featuredata_current)), by = 0, all = TRUE)
      data_to_plot[is.na(data_to_plot)] <- 0.1
      data_to_plot[data_to_plot == 0] <- 0.1
      data_to_plot <- transform(data_to_plot, Status = ifelse(reads_now == 0.1, 
                                                              "removed", "kept"))
      # Attach last taxonomy
      data_to_plot <- 
        merge(data_to_plot, 
              data.frame(Taxonomy = apply(reactives$taxonomy_data, 1, function(x) 
                tail(x[x != "No_annotation"], 1))),
              by.x = "Row.names", by.y = 0, all.x = TRUE)
      abundance_plot <- 
        ggplot(data = data_to_plot, aes(x = reads_before, y = reads_now, 
                                        colour = Status)) +
        geom_abline(slope = 1, intercept = 0,  colour = "#9ed3f3") +
        geom_point(aes(text = paste("Feature:", Row.names,
                                    "\nTaxonomy:", Taxonomy,
                                    "\nReads before:", reads_before,
                                    "\nReads now:", round(reads_now, 0),
                                    "\nStatus:", Status))) +
        scale_colour_manual(values = c("#2fa4e7", "#757575")) +
        scale_x_log10() + scale_y_log10() +
        xlab("Reads before") + ylab("Reads now") +
        plot_theme
      output$plot <- renderPlotly({
        ggplotly(abundance_plot, tooltip = "text")
      }) 
    }
    if(input$visualization_type == "Contamination removal - NEG1" ||
       input$visualization_type == "Contamination removal - NEG2") {
      if(reactives$step_var < 5) {
        showModal(modalDialog(title = "Warning 1", paste0(
        "The plot for contamination removal is only available at step 5. You are 
        currently at step ", reactives$step_var, ". Please dismiss this message 
        and proceed with your analysis by clicking the Next-button.")))
        updateSelectInput(session, inputId = "visualization_type",
                          selected = "Correlation of reads and features")
      } else {
        if(input$visualization_type == "Contamination removal - NEG1") {
          req_span_neg1 <- reactives$req_span_neg1
          req_ratio_neg1 <- input$req_ratio_neg1
          contamination_plot <- 
            ggplot(data = reactives$filter_basis, 
                   aes(x = ratio_neg1, y = neg1_span, size = sample_mean,
                       colour = status))
          if(as.numeric(req_span_neg1) != 0.0001) {
            contamination_plot <- contamination_plot +
              geom_hline(aes(yintercept = as.numeric(req_span_neg1)))
          }
          if(as.numeric(req_ratio_neg1) != Inf) {
            contamination_plot <- contamination_plot +
              geom_vline(aes(xintercept = as.numeric(req_ratio_neg1)))
          }
          contamination_plot <- contamination_plot +
            geom_point(aes(text = paste("Feature:", rownames(reactives$filter_basis),
                                        "\nTaxonomy:", reactives$filter_basis$Taxonomy,
                                        "\nRatio (NEG1):", round(ratio_neg1, 3),
                                        "\nSpan (NEG1):", round(neg1_span, 2),
                                        "\nMean rel. abundance (samples):", signif(sample_mean, 2),
                                        "\nStatus:", status))) +
            scale_x_log10() +
            scale_colour_manual("Status", values = contamination_neg1_colours) +
            scale_size_continuous("") +
            xlab("Ratio (NEG1)") + ylab("Span (NEG1)") +
            plot_theme
          output$plot <- renderPlotly({
            ggplotly(contamination_plot, tooltip = "text")
          })
        }
        if(input$visualization_type == "Contamination removal - NEG2") {
          req_span_neg2 <- reactives$req_span_neg2
          req_ratio_neg2 <- input$req_ratio_neg2
          contamination_plot <- 
            ggplot(data = reactives$filter_basis, 
                   aes(x = ratio_neg2, y = neg2_span, size = sample_mean,
                       colour = status))
          if(as.numeric(req_span_neg2) != 0.0001) {
            contamination_plot <- contamination_plot +
              geom_hline(aes(yintercept = as.numeric(req_span_neg2)))
          }
          if(as.numeric(req_ratio_neg2) != Inf) {
            contamination_plot <- contamination_plot +
              geom_vline(aes(xintercept = as.numeric(req_ratio_neg2)))
          }
          contamination_plot <- contamination_plot +
            geom_point(aes(text = paste("Feature:", rownames(reactives$filter_basis),
                                        "\nTaxonomy:", reactives$filter_basis$Taxonomy,
                                        "\nRatio (NEG2):", round(ratio_neg2, 3),
                                        "\nSpan (NEG2):", round(neg2_span, 2),
                                        "\nMean rel. abundance (samples):", signif(sample_mean, 2),
                                        "\nStatus:", status))) +
            scale_x_log10() +
            scale_colour_manual("Status", values = contamination_neg2_colours) +
            scale_size_continuous("") +
            xlab("Ratio (NEG2)") + ylab("Span (NEG2)") +
            plot_theme
          output$plot <- renderPlotly({
            ggplotly(contamination_plot, tooltip = "text")
          })          
        }
      }
    }
    if(input$visualization_type == "Reduction of total reads") {
      sum_of_reads_function <- function(x) {
        if(is.na(x)) {return(NA)} else {return(sum(colSums(x)))}
      }
      data_to_plot <- data.frame(
        Step = c("0 - original", "1 - without controls", "2 - sample filter",
                 "3 - feature abundance filter", "4 - feature frequency filter",
                 "5 - contamination filter"),
        Reads = c(sum_of_reads_function(reactives$featuredata),
                  sum_of_reads_function(reactives$featuredata_1),
                  sum_of_reads_function(reactives$featuredata_2),
                  sum_of_reads_function(reactives$featuredata_3),
                  sum_of_reads_function(reactives$featuredata_4),
                  sum_of_reads_function(reactives$featuredata_5)))
      reduction_of_reads <- 
        ggplot(data = data_to_plot, aes(x = Step, y = Reads, fill = Step)) +
        geom_bar(aes(text = paste("step:", Step,
                                  "\nreads:", Reads)), stat = "identity") +
        scale_fill_manual(values = c("#E3F0F7", "#C4E3F4", "#A4D6F2", "#85C8EF",
                                     "#65BBEC", "#2FA4E7")) +
        scale_x_discrete(labels = c("0", "1", "2", "3", "4", "5")) +
        plot_theme
      output$plot <- renderPlotly({
        ggplotly(reduction_of_reads, tooltip = "text")
      }) 
    }
    # Save current input:
    reactives$req_reads_per_sample_old <- input$req_reads_per_sample
    reactives$req_reads_per_feature_old <- input$req_reads_per_feature
    reactives$req_ratio_per_feature_old <- input$req_ratio_per_feature
    reactives$req_ratio_neg1_old <- input$req_ratio_neg1
    reactives$req_span_neg1_old <- reactives$req_span_neg1
    reactives$req_ratio_neg2_old <- input$req_ratio_neg2
    reactives$req_span_neg2_old <- reactives$req_span_neg2
  }
  
  # ----------------------------------------------------------------------------
  # Changed input in span filter
  # ----------------------------------------------------------------------------
  observeEvent(input$req_span_neg1, {
    reactives$req_span_neg1 <- input$req_span_neg1
  })
  observeEvent(input$req_span_neg2, {
    reactives$req_span_neg2 <- input$req_span_neg2
  })
  
  # ----------------------------------------------------------------------------
  # Update button
  # ----------------------------------------------------------------------------
  observeEvent(input$update_button, {
    print(paste0("INFO - update plot - ", Sys.time()))
    filter_feature_table()
    shinyjs::hide("text")
    shinyjs::hide("table")
    output$table <- NULL
  })
  
  # ----------------------------------------------------------------------------
  # Next button
  # ----------------------------------------------------------------------------
  observeEvent(input$next_button, {
    reactives$step_var <- reactives$step_var + 1
    print(paste0("INFO - next to step ", reactives$step_var, " - ", Sys.time()))
    filter_feature_table()
    step_var_UIchange()
    shinyjs::hide("text")
    shinyjs::hide("table")
    output$table <- NULL
  })

  # ----------------------------------------------------------------------------
  # Back button
  # ----------------------------------------------------------------------------
  observeEvent(input$back_button, {
    # Delete the previous feature and metadata
    eval(parse(text = paste0("reactives$featuredata_", reactives$step_var, 
                             " <- NA")))
    eval(parse(text = paste0("reactives$metadata_", reactives$step_var, 
                             " <- NA")))
    reactives$step_var <- reactives$step_var - 1
    print(paste0("INFO - back to step ", reactives$step_var, " - ", Sys.time()))
    step_var_UIchange()
  })
  
  # ----------------------------------------------------------------------------
  # Step_var_change function: build UI (and disable buttons depending on step)
  # ----------------------------------------------------------------------------
  step_var_UIchange <- function(){
    print(paste0("INFO - change UI for step ", reactives$step_var, " - ", 
                 Sys.time()))
    if(reactives$step_var == 1){
      shinyjs::enable("metafile")
      shinyjs::enable("featurefile")
      shinyjs::hide("visualization_type")
      shinyjs::hide("req_reads_per_sample")
    }
    if(reactives$step_var == 2){
      shinyjs::disable("metafile")
      shinyjs::disable("featurefile")
      shinyjs::show("visualization_type")
      shinyjs::show("req_reads_per_sample") 
      shinyjs::enable("req_reads_per_sample") 
      shinyjs::hide("req_reads_per_feature")
    }
    if(reactives$step_var == 3) {
      shinyjs::disable("req_reads_per_sample") 
      shinyjs::show("req_reads_per_feature")
      shinyjs::enable("req_reads_per_feature")
      shinyjs::hide("req_ratio_per_feature")
    }
    if(reactives$step_var == 4) {
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
      shinyjs::disable("req_ratio_per_feature")
      shinyjs::show("header_neg1") 
      shinyjs::show("req_ratio_neg1")
      insertUI(selector = "#placeholder_span_1",
               ui = selectInput(inputId = "req_span_neg1",
                                label = "Minimum span threshold (NEG1)",		
                                choices = reactives$neg1_span_steps))
      shinyjs::show("header_neg2")
      shinyjs::show("req_ratio_neg2")
      insertUI(selector = "#placeholder_span_2",
               ui = selectInput(inputId = "req_span_neg2",
                                label = "Minimum span threshold (NEG2)",		
                                choices = reactives$neg2_span_steps))
    }

    # --------------------------------------------------------------------------
    # Build UI for data analysis after filtering is completed
    # --------------------------------------------------------------------------
    if(reactives$step_var == 6) {
      if(reactives$analysis_started == FALSE) {
        reactives$analysis_started <- TRUE
        appendTab(
          inputId = "tabset",
          tabPanel(
            title = "Alpha diversity",
            br(),
            selectInput(inputId = "metavar_alpha", 
                        label = "Metainformation variable",
                        choices = colnames(reactives$metadata_6)),
            selectInput(inputId = "subvar_alpha",
                        label = "Variable for subgroup analysis",
                        choices = c("ignore", colnames(reactives$metadata_6))),
            radioButtons(inputId = "scaling", 
                         label = "Y axis scaling",
                         choices = c("Raw values", "Z-normalized values")),
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
                        choices = colnames(reactives$metadata_6)),
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
                        choices = colnames(reactives$metadata_6)),
            radioButtons(inputId = "taxonomy_level", 
                         label = "Taxonomic level", 
                         choices = c("Domain", "Phylum", "Class", "Order", 
                                     "Family", "Genus", "Species")),
            numericInput(inputId = "top_number_taxa",
                         label = "Top number of taxa displayed",
                         value = 10, min = 1, max = 20, step = 1),
            radioButtons(inputId = "per_group_overall", 
                         label = "Top number of taxa based on", 
                         choices = c("Overall", "Per group")),
            actionButton(inputId = "tax_analysis", 
                         label = "Update plot")
          )
        )
        appendTab(
          inputId = "tabset",
          tabPanel(
            title = "Sample selection",
            br(),
            lapply(colnames(reactives$metadata_6), function(x) {
              pickerInput(inputId = paste0("sampleselector_", x),
                          label = x,
                          choices = unique(reactives$metadata_6[, x]),
                          selected = unique(reactives$metadata_6[, x]),
                          multiple = TRUE,
                          options = list(`actions-box` = TRUE,
                                         `selected-text-format` = "count",
                                         `count-selected-text` = "{0}/{1} selected"))
            }) # Close lapply function
          )
        )
        appendTab(
          inputId = "tabset",
          tabPanel(
            title = "Download",
            br(),
            h5(tags$b("Final filtered data output:")),
            downloadButton("download_feature_abs", "Final feature file (counts)"),
            downloadButton("download_feature_rel", "Final feature file (rel. abundance)"),
            downloadButton("download_metafile", "Final metafile"),
            br(), br(),
            h5(tags$b("Final filter settings:")),
            downloadButton("download_filtersettings", "Filter settings"),
            br(), br(),
            h5(tags$b("Files for quality control:")),
            downloadButton("download_contbasis", "Basis for contamination filtering"),
            downloadButton("download_redfeature", "Read reduction per feature"),
            downloadButton("download_redtaxonomy", "Read reduction per taxonomy"),
            br(), br(),
            h5(tags$b("Basic diversity analyses (all samples):")),
            downloadButton("download_allalpha", "Alpha diversity values"),
            downloadButton("download_allbeta", "Beta diversity distance matrix"),
            br(), br(),
            h5(tags$b("Current analysis values:")),
            h6(tags$b("Analysis-specific values will be available for download 
                      once the respective analysis was performed."), 
               id = "analysis_first"),
            downloadButton(outputId = "download_current_alpha",
                                label = "Current alpha diversity values"),
            downloadButton(outputId = "download_current_beta",
                                label = "Current beta diversity values"),
            downloadButton(outputId = "download_current_taxonomy",
                                label = "Current taxonomy values")
          )
        )
      }
      updateTabsetPanel(session, "tabset", selected = "Download")
    }    
  }
  
  # ----------------------------------------------------------------------------
  # Prevent the user from messing up the files between filtering and analysis,
  #   show download buttons only for available current analysis values
  # ----------------------------------------------------------------------------
  observeEvent(input$tabset, {
    if(reactives$step_var == 6 && input$tabset == "Filtering") {
      showModal(modalDialog(
        title = "Are you sure?",
        "Returning to filtering may overwrite your final filtered files.",
        footer = tagList(actionButton("refilter_button", "Go back to filtering"),
                         modalButton("Stay on analysis site"))
      ))
      updateTabsetPanel(session, "tabset", selected = "Alpha diversity")
    }
    if(reactives$step_var <= 5 && input$tabset != "Filtering") {
      showModal(modalDialog(
        title = "Are you sure?",
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
    # Hide download buttons for current analysis when analysis was not yet done
    if(is.null(reactives$alpha_diversity_download)) {
      shinyjs::hide("download_current_alpha")
    } else {shinyjs::show("download_current_alpha")}
    if(is.null(reactives$beta_diversity_download)) {
      shinyjs::hide("download_current_beta")
    } else {shinyjs::show("download_current_beta")}
    if(is.null(reactives$taxonomy_download)) {
      shinyjs::hide("download_current_taxonomy")
    } else {shinyjs::show("download_current_taxonomy")}
    if(!is.null(reactives$alpha_diversity_download) && 
       !is.null(reactives$beta_diversity_download) &&
       !is.null(reactives$taxonomy_download)) {
      shinyjs::hide("analysis_first")
    }
  })
  
  observeEvent(input$refilter_button, {
    reactives$step_var <- reactives$step_var - 1
    updateTabsetPanel(session, "tabset", selected = "Filtering")
    removeModal()
  })
  
  observeEvent(input$reanalysis_button, {
    reactives$step_var <- 6
    updateTabsetPanel(session, "tabset", selected = "Alpha diversity") 
    removeModal()
  })
  
  # ----------------------------------------------------------------------------
  # Alpha diversity analysis
  # ----------------------------------------------------------------------------
  observeEvent(input$alpha_analysis, {
    print(paste0("INFO - alpha diversity analysis - ", Sys.time()))
    # Save input as reactive variables for download
    reactives$metavar_alpha <- input$metavar_alpha
    reactives$subvar_alpha <- input$subvar_alpha
    # Filter feature and metadata according to sample selection
    samples_deselected <- Sample_selection_check()
    Metadata <- 
      reactives$metadata_6[!rownames(reactives$metadata_6) %in% 
                                   samples_deselected, ]
    # Select alpha diversity values for selected samples
    alpha_diversity_result <- 
      reactives$alpha_diversity_data[!rownames(reactives$metadata_6) %in% 
                                         samples_deselected ,]
    # Perform z-transformation if needed
    if (input$scaling == "Z-normalized values") {
      alpha_diversity_result <- 
        as.data.frame(apply(alpha_diversity_result, 2, function(x)
          scale(x, center = TRUE, scale = TRUE)), 
          row.names = rownames(alpha_diversity_result))
    }
    # Attach variables of interest for alpha diversity
    alpha_diversity_result["Group"] <- Metadata[paste0(reactives$metavar_alpha)]
    if(reactives$subvar_alpha != "ignore") {
      alpha_diversity_result["Subgroup"] <- Metadata[paste0(reactives$subvar_alpha)]
    }
    # Generate format for download
    reactives$alpha_diversity_download <- alpha_diversity_result
    # Calculate pvalue table for alpha diversity
    if(reactives$subvar_alpha == "ignore" ||
       (reactives$subvar_alpha != "ignore" && 
        nrow(unique(Metadata[paste0(reactives$subvar_alpha)])) == 1)) {
      alpha_diversity_table <- as.data.frame(
        sapply(colnames(reactives$alpha_diversity_data), function(y) 
          if(length(unique(alpha_diversity_result$Group)) == 1) {
            1 # print 1 if there is only one category in "Group"
          } else {
            kruskal.test(formula(paste(y, "~ Group")),
                         data = alpha_diversity_result)$p.value}))
      colnames(alpha_diversity_table) <- "Group"
      alpha_diversity_table <- as.data.frame(t(alpha_diversity_table))
    } else {
      alpha_diversity_table <- as.data.frame(
        sapply(colnames(reactives$alpha_diversity_data), function(y) 
          sapply(unique(alpha_diversity_result$Subgroup), function(x)
            if(length(unique(subset(alpha_diversity_result, 
                                    Subgroup == x)[, "Group"])) == 1) {
              1 # print 1 if there is only one category in "Group" per subgroup
            } else {
              kruskal.test(formula(paste(y, "~ Group")), 
                           data = subset(alpha_diversity_result, 
                                         Subgroup == x))$p.value})))
    }
    # Round pvalues of the table
    alpha_diversity_table <- signif(alpha_diversity_table, 2)
    # Generate the alpha diversity plot
    if(reactives$subvar_alpha != "ignore") {
      alpha_diversity_result <- melt(alpha_diversity_result,
                                        id.vars = c("Group", "Subgroup"))
    } else {
      alpha_diversity_result <- melt(alpha_diversity_result,
                                          id.vars = "Group")
    }
    # Build alpha diversity plot
    alpha_diversity_plot <- 
      ggplot(data = alpha_diversity_result, aes(x = Group, y = value)) +
      geom_boxplot(aes(colour = Group, fill = Group)) +
      geom_point(aes(text = paste0("Group: ", Group,
                                  "\n", variable, " value: ", round(value, 3))), size = 1.5) +
      theme(strip.text.x = element_text(size = 8)) +
      plot_theme +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
            axis.title = element_blank()) +
      scale_colour_manual(reactives$metavar_alpha, values = theme_colours) +
      scale_fill_manual(reactives$metavar_alpha, values = alpha(theme_colours, 0.5))
    # Build facets based on selected number of variables and scaling
    if (reactives$subvar_alpha != "ignore" && input$scaling == "Z-normalized values") {
      alpha_diversity_plot <- alpha_diversity_plot + 
        facet_wrap(Subgroup ~ variable, 
                   nrow = length(unique(alpha_diversity_result$Subgroup)),
                   scales = "free_x")
    } else if (reactives$subvar_alpha != "ignore") {
      alpha_diversity_plot <- alpha_diversity_plot +
        facet_wrap(Subgroup ~ variable, 
                   ncol = length(levels(alpha_diversity_result$variable)), 
                   scales = "free")      
    } else if (input$scaling == "Z-normalized values") {
      alpha_diversity_plot <- alpha_diversity_plot +
        facet_wrap(. ~ variable, 
                   ncol = length(levels(alpha_diversity_result$variable)))
    } else {
      alpha_diversity_plot <- alpha_diversity_plot +
        facet_wrap(. ~ variable, 
                   ncol = length(levels(alpha_diversity_result$variable)), 
                   scales = "free_y")
    }
    # Output plot, empty line and p-value table
    output$plot <- renderPlotly({
      ggplotly(alpha_diversity_plot, tooltip = "text")
    })
    output$text <- renderText({""})
    output$table <- DT::renderDataTable({alpha_diversity_table})
    shinyjs::show("text")
    shinyjs::show("table")
  })
  
  # ----------------------------------------------------------------------------
  # Beta diversity analysis
  # ----------------------------------------------------------------------------
  observeEvent(input$beta_analysis, {
    print(paste0("INFO - beta diversity analysis - ", Sys.time()))
    # Save input as reactive variables for download
    reactives$metavar_beta <- input$metavar_beta
    reactives$plot_beta <- input$plot_beta
    # Filter feature and metadata according to sample selection
    # Add a variable to check that at least two samples are selected
    reactives$beta_diversity_analysis <- TRUE
    samples_deselected <- Sample_selection_check()
    reactives$beta_diversity_analysis <- FALSE
    Metadata <- 
      reactives$metadata_6[!rownames(reactives$metadata_6) %in% 
                                   samples_deselected, ]
    Featuredata <- 
      reactives$featuredata_6[, !colnames(reactives$featuredata_6) %in%
                                      samples_deselected]
    Featuredata_freq <- decostand(t(Featuredata), method = "total")
    # Calculate visualization for nMDS
    if (reactives$plot_beta == "beta_nmds") {
      set.seed(1)
      beta_diversity_data <- 
        metaMDS(Featuredata_freq, distance = "bray", k = 2, 
                autotransform = FALSE)
      beta_diversity_result <- data.frame(
        Axis1 = beta_diversity_data$points[, 1],
        Axis2 = beta_diversity_data$points[, 2])
      rownames(beta_diversity_result) <- 
        dimnames(beta_diversity_data$points)[[1]]
    # Calculate visualization for PCoA
    } else if (reactives$plot_beta == "beta_pcoa") {
      set.seed(1)
      beta_diversity_data <- 
        cmdscale(vegdist(Featuredata_freq, distance = "bray"), k = 2, 
                 eig = FALSE, add = FALSE, x.ret = FALSE)
      beta_diversity_result <- data.frame(
        Axis1 = beta_diversity_data[, 1],
        Axis2 = beta_diversity_data[, 2])
      rownames(beta_diversity_result) <- dimnames(beta_diversity_data)[[1]]
    }
    # Attach metadata variable and calculate pvalue for beta diversity
    beta_diversity_result[reactives$metavar_beta] <- 
      Metadata[paste0(reactives$metavar_beta)]
    if(length(unique(beta_diversity_result[[reactives$metavar_beta]])) <= 1) {
      beta_diversity_pvalue <- NA
    } else {
      set.seed(1)
      beta_diversity_pvalue <- 
        signif(adonis(vegdist(Featuredata_freq, method = "bray") ~ 
                 as.factor(beta_diversity_result[[reactives$metavar_beta]]))$aov.tab$'Pr(>F)'[1], 2)
    }
    # Generate format for download
    reactives$beta_diversity_download <- beta_diversity_result
    # Output plot and pvalue
    output$plot <- renderPlotly({
      ggplot(data = beta_diversity_result, 
             aes(x = Axis1, y = Axis2, colour = !!sym(reactives$metavar_beta))) +
        geom_point(aes(text = paste("Sample:", rownames(beta_diversity_result)))) +
        stat_ellipse() +
        plot_theme +
        scale_colour_manual(values = theme_colours)
    })
    output$text <- renderText({
      paste0("p-value = ", beta_diversity_pvalue)
    })
    shinyjs::show("text")
    shinyjs::hide("table")
  })
  
  # ----------------------------------------------------------------------------
  # Taxonomy analysis
  # ----------------------------------------------------------------------------
  observeEvent(input$tax_analysis, {
    print(paste0("INFO - taxonomy analysis - ", Sys.time()))
    # Save input as reactive variables for download
    reactives$metavar_taxonomy <- input$metavar_taxonomy
    reactives$taxonomy_level <- input$taxonomy_level
    # Filter feature and metadata according to sample selection
    samples_deselected <- Sample_selection_check()
    Metadata <- 
      reactives$metadata_6[!rownames(reactives$metadata_6) %in% 
                                   samples_deselected, ]
    Featuredata <- 
      reactives$featuredata_6[, !colnames(reactives$featuredata_6) %in%
                                      samples_deselected, drop = FALSE]
    Featuredata_freq <- t(decostand(t(Featuredata), method = "total"))
    # Merge taxonomic level and featuredata
    taxonomy_result <- 
      merge(Featuredata_freq, reactives$taxonomy_data[reactives$taxonomy_level], 
            by = 0, all.x = TRUE)
    # Sum up frequencies with identical taxonomy
    taxonomy_result <- taxonomy_result %>% 
      group_by(.data[[reactives$taxonomy_level]]) %>%
      summarise(across(.cols = is.numeric, .fns = sum, .names = "{col}")) %>%
      as.data.frame()
    # Check that the chosen number of taxa is available
    if(input$top_number_taxa > nrow(taxonomy_result)) {
      updateNumericInput(session, inputId = "top_number_taxa", 
                         value = nrow(taxonomy_result))
      reactives$top_number_taxa <- nrow(taxonomy_result)
      showModal(modalDialog(
        title = "Warning 3", 
        paste0("There are only ", nrow(taxonomy_result), " different taxa at ", 
               reactives$taxonomy_level, " level. The top number of taxa 
               displayed is set to ", nrow(taxonomy_result), ".")))
    } else {
      reactives$top_number_taxa <- input$top_number_taxa
    }
    rownames(taxonomy_result) <- taxonomy_result[[reactives$taxonomy_level]]
    taxonomy_result[[reactives$taxonomy_level]] <- NULL
    taxonomy_result <- as.data.frame(t(taxonomy_result))
    taxonomy_result <- 
      merge(taxonomy_result, 
            Metadata[, colnames(Metadata) == reactives$metavar_taxonomy, drop = FALSE], 
            by = 0, all = TRUE)
    taxonomy_result[[reactives$metavar_taxonomy]] <- 
      as.character(taxonomy_result[[reactives$metavar_taxonomy]])
    # Build mean frequency per selected group
    taxonomy_result <- taxonomy_result %>%
      group_by(.data[[reactives$metavar_taxonomy]]) %>%
      summarise(across(.cols = is.numeric, .fns = mean, .names = "{col}")) %>%
      as.data.frame()
    # Select top taxa per group or overall
    if(input$per_group_overall == "Overall") {
      top_taxa <- taxonomy_result %>%
        select_if(is.numeric) %>%
        colMeans() %>% 
        sort(decreasing = TRUE)
      top_taxa <- names(top_taxa[1:reactives$top_number_taxa])
    } else if (input$per_group_overall == "Per group") {
      top_taxa <- taxonomy_result %>%
        group_by(.data[[reactives$metavar_taxonomy]]) %>%
        select_if(is.numeric) %>%
        group_map(~ names(sort(.x, decreasing = TRUE)[1:reactives$top_number_taxa]))
      top_taxa <- unique(unlist(top_taxa))
    }
    # Summarise other taxa into "Others"
    top_taxa <- sort(top_taxa)
    taxonomy_others <- taxonomy_result %>% 
      select_if(is.numeric) %>%
      select(-all_of(top_taxa)) %>%
      rowSums() %>% 
      data.frame
    colnames(taxonomy_others) <- "Others"
    # Prevent showing an empty "Others" category in the plot
    if(colSums(taxonomy_others) != 0) {
      taxonomy_result <- merge(taxonomy_others,
                             taxonomy_result[, c(top_taxa, reactives$metavar_taxonomy)],
                             by = 0, all = TRUE)
    } else { # Prevent a bug in melt() when "Others" don't exist
      taxonomy_result["Row.names"] <- rownames(taxonomy_result)
    }
    # Generate format for download
    reactives$taxonomy_download <- taxonomy_result[, order(ncol(taxonomy_result):1)]
    # Plot the result
    taxonomy_result <- melt(
      taxonomy_result, id.vars = c("Row.names", reactives$metavar_taxonomy), 
      value.name = "Relative abundance") # Space needed to avoid a plotly bug
    colnames(taxonomy_result)[colnames(taxonomy_result) == "variable"] <- "Taxonomy"
    taxonomy_plot <- ggplot(data = taxonomy_result, 
                            aes(x = !!sym(reactives$metavar_taxonomy), 
                                y = !!sym("Relative abundance"),  
                                fill = Taxonomy)) +
      geom_bar(position = "fill", stat = "identity") + 
      plot_theme +
      ylab("Relative abundance") +
      scale_y_continuous(expand = c(0, 0)) +
      xlab(reactives$metavar_taxonomy)
    if(colSums(taxonomy_others) != 0) {
      taxonomy_plot <- taxonomy_plot +
        scale_fill_manual("Taxonomy", values = c("#d1d1d1", theme_colours))
    } else {
      taxonomy_plot <- taxonomy_plot +
        scale_fill_manual("Taxonomy", values = theme_colours)
    }
    output$plot <- renderPlotly({
      taxonomy_plot 
    })
    shinyjs::hide("text")
    shinyjs::hide("table")
  })

  # Check that the user does not select a taxonomy level that is not available
  observeEvent(input$taxonomy_level, {
    if(!(input$taxonomy_level %in% colnames(reactives$taxonomy_data))) {
      showModal(modalDialog(
        title = "Warning 4", paste0("It seems that your taxonomic annotation does 
        not provide information on ", input$taxonomy_level, " level. Please 
        select a lower taxonomic level for analysis.")))
      updateRadioButtons(session, inputId = "taxonomy_level", 
                         selected = "Domain")
    }
  })
  
  # ----------------------------------------------------------------------------
  # Sample_selection_check returns deselected sample names from Sample selection
  # ----------------------------------------------------------------------------
  Sample_selection_check <- function() {
    print(paste0("INFO - Sample selection check - ", Sys.time()))
    samples_deselected <- rownames(reactives$metadata_6)[sort(
      unique(unlist(sapply(colnames(reactives$metadata_6), function(x) {
        which(!reactives$metadata_6[, x] %in% 
                input[[paste0("sampleselector_", x)]])
      }
      )))
    )]
    # Check that the user does not deselect all samples
    if (length(samples_deselected) == nrow(reactives$metadata_6)) {
      showModal(modalDialog(
        title = "Warning 2", "It seems that you deselected all samples. Sample
        selection is set back, and all samples are selected for analysis."))
      lapply(colnames(reactives$metadata_6), function(x) {
        updatePickerInput(session, inputId = paste0("sampleselector_", x),
                          selected = unique(reactives$metadata_6[, x]))
      }) # Close lapply function
      samples_deselected <- c()
    }
    # Check that the user does select more than 2 samples for beta diversity
    if ((nrow(reactives$metadata_6)-length(samples_deselected)) <= 2 && 
        reactives$beta_diversity_analysis == TRUE) {
      showModal(modalDialog(
        title = "Warning 5", "It seems that you selected two or less samples
        for beta diversity analysis. Sample selection is set back, and all 
        samples are selected for analysis."))
      lapply(colnames(reactives$metadata_6), function(x) {
        updatePickerInput(session, inputId = paste0("sampleselector_", x),
                          selected = unique(reactives$metadata_6[, x]))
      }) # Close lapply function
      samples_deselected <- c()
    }
    print(paste0("INFO - ", 
                 nrow(reactives$metadata_6) - length(samples_deselected), 
                 " samples selected for analysis"))
    return(samples_deselected)
  }
  
  # ----------------------------------------------------------------------------
  # Data download - final filtered files
  # ----------------------------------------------------------------------------
  # Download final filtered feature table as counts
  output$download_feature_abs <- downloadHandler(
    filename = function() {
      paste0("MicrobIEM_Featuretable_abs_final_", 
             format(Sys.time(), "%Y_%m_%d"), ".txt")
    },
    content = function(file) {
      write.table(reactives$download_featuredata_abs, file = file, 
                  row.names = FALSE, sep = "\t", quote = FALSE)
    }
  )
  # Download final filtered feature table as relative abundances
  output$download_feature_rel <- downloadHandler(
    filename = function() {
      paste0("MicrobIEM_Featuretable_rel_final_", 
             format(Sys.time(), "%Y_%m_%d"), ".txt")
    },
    content = function(file) {
      write.table(reactives$download_featuredata_rel, file = file, 
                  row.names = FALSE, sep = "\t", quote = FALSE)
    }
  )
  # Download final metafile
  output$download_metafile <- downloadHandler(
    filename = function() {
      paste0("MicrobIEM_Metatable_final_", 
             format(Sys.time(), "%Y_%m_%d"), ".txt")
    },
    content = function(file) {
      write.table(reactives$metadata_6, file = file, row.names = FALSE, 
                  sep = "\t", quote = FALSE)
    }
  )
  # ----------------------------------------------------------------------------
  # Data download - Filter settings
  # ----------------------------------------------------------------------------
  output$download_filtersettings <- downloadHandler(
    filename = function() {
      paste0("MicrobIEM_Filter_settings_", 
             format(Sys.time(), "%Y_%m_%d"), ".txt")
    },
    content = function(file) {
      write.table(reactives$download_filtersettings, file = file, 
                  col.names = FALSE, row.names = FALSE, quote = FALSE)
    }
  )
  
  # ----------------------------------------------------------------------------
  # Data download - Quality control
  # ----------------------------------------------------------------------------
  # Download contamination filter basis
  output$download_contbasis <- downloadHandler(
    filename = function() {
      paste0("MicrobIEM_Contamination_filter_basis_", 
             format(Sys.time(), "%Y_%m_%d"), ".txt")
    },
    content = function(file) {
      write.table(data.frame(OTU_ID = rownames(reactives$filter_basis), 
                             reactives$filter_basis), file = file,
                  row.names = FALSE, sep = "\t", quote = FALSE)
    }
  )
  # Download reduction of reads per feature by filter step
  output$download_redfeature <- downloadHandler(
    filename = function() {
      paste0("MicrobIEM_Read_reduction_per_feature_", 
             format(Sys.time(), "%Y_%m_%d"), ".txt")
    },
    content = function(file) {
      write.table(reactives$download_reduction_per_feature, file = file,
                  row.names = FALSE, sep = "\t", quote = FALSE)
    }
  )
  # Download reduction of reads per highest taxonomy by filter step 
  output$download_redtaxonomy <- downloadHandler(
    filename = function() {
      paste0("MicrobIEM_Read_reduction_per_taxonomy_", 
             format(Sys.time(), "%Y_%m_%d"), ".txt")
    },
    content = function(file) {
      write.table(reactives$download_reduction_per_taxonomy, file = file, 
        row.names = FALSE, sep = "\t", quote = FALSE)
    }
  )
  
  # ----------------------------------------------------------------------------
  # Data download - Basic analysis values
  # ----------------------------------------------------------------------------
  # Save alpha diversity values
  output$download_allalpha <- downloadHandler(
    filename = function() {
      paste0("MicrobIEM_Alpha_diversity_values_", 
             format(Sys.time(), "%Y_%m_%d"), ".txt")
    },
    content = function(file) {
      write.table(data.frame(Sample_ID = rownames(reactives$alpha_diversity_data), 
                   reactives$alpha_diversity_data), file = file,
        row.names = FALSE, sep = "\t", quote = FALSE)
    }
  )
  # Save beta diversity values
  output$download_allbeta <- downloadHandler(
    filename = function() {
      paste0("MicrobIEM_Beta_diversity_distance_matrix_", 
             format(Sys.time(), "%Y_%m_%d"), ".txt")
    },
    content = function(file) {
      write.table(data.frame(Sample_ID = rownames(reactives$distance_matrix), 
                             reactives$distance_matrix, check.names = FALSE), 
                  file = file, row.names = FALSE, col.names = TRUE, 
                  sep = "\t", quote = FALSE)
    }
  )
  
  # ----------------------------------------------------------------------------
  # Current data download - alpha diversity
  # ----------------------------------------------------------------------------
  output$download_current_alpha <- downloadHandler(
    filename = function() {
      paste0(
        "MicrobIEM_Alpha_diversity_",
        gsub("[^0-9a-zA-Z_-]", "", reactives$metavar_alpha), "_", # Variable name
        if(reactives$subvar_alpha != "ignore") {
          paste0(gsub("[^0-9a-zA-Z_-]", "", reactives$subvar_alpha), "_")
        },
        format(Sys.time(), "%Y_%m_%d"), ".txt") # Date of download
    }, 
    content = function(file) {
      write.table(
        data.frame(Sample_ID = rownames(reactives$alpha_diversity_download),
                   reactives$alpha_diversity_download), 
        file = file, sep = "\t", dec = ".", quote = FALSE, row.names = FALSE)
    }
  )

  # ----------------------------------------------------------------------------
  # Current data download - beta diversity
  # ----------------------------------------------------------------------------
  output$download_current_beta <- downloadHandler(
    filename = function() {
      paste0(
        "MicrobIEM_Beta_diversity_",
        substr(reactives$plot_beta, 6, 9), "_", # nMDS or PCoA
        gsub("[^0-9a-zA-Z_-]", "", reactives$metavar_beta), "_", # Variable name
        format(Sys.time(), "%Y_%m_%d"), ".txt") # Date of download
    },
    content = function(file) {
      write.table(
        data.frame(Sample_ID = rownames(reactives$beta_diversity_download),
                   reactives$beta_diversity_download), 
        file = file, sep = "\t", dec = ".", quote = FALSE, row.names = FALSE)
    }
  )

  # ----------------------------------------------------------------------------
  # Current data download - taxonomy
  # ----------------------------------------------------------------------------
  output$download_current_taxonomy <- downloadHandler(
    filename = function() {
      paste0(
        "MicrobIEM_Taxonomy_",
        reactives$taxonomy_level, "_", # Taxonomic level
        gsub("[^0-9a-zA-Z_-]", "", reactives$metavar_taxonomy), "_", # Variable name
        format(Sys.time(), "%Y_%m_%d"), ".txt") # Date of download
    },
    content = function(file) {
      write.table(
        select(reactives$taxonomy_download, -Row.names), 
        file = file, sep = "\t", dec = ".", quote = FALSE, row.names = FALSE)
    }
  )
}
