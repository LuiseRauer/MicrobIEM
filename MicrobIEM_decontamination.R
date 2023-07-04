################################################################################
#
# MicrobIEM - an algorithm to reduce contaminants in microbiome 16S amplicon 
#   sequencing data
#
################################################################################

# MicrobIEM decontamination function
MicrobIEM_decontamination <- function(
  # data.frame, a feature table with count data, with ASVs/OTUs/taxa in rows, ASV IDs as rownames, and samples in columns
  feature_table, 
  # character vector, names of samples to be filtered 
  SAMPLE, 
  # character vector, names of NEG1 controls, can by empty
  NEG1 = character(0), 
  # character vector of names of NEG2 controls, can be empty
  NEG2 = character(0), 
  # Numeric, thresholds for filtering contaminants, default is to not use any filter (NA)
  ratio_NEG1_threshold = NA, # Sensible thresholds: > 0
  ratio_NEG2_threshold = NA, # Sensible thresholds: > 0
  span_NEG1_threshold = NA, # Sensible thresholds: between 0 and 1
  span_NEG2_threshold = NA # Sensible thresholds: between 0 and 1
  ){
  
  # Check input data types
  if(!(is(feature_table, "data.frame")))
    stop("feature_table must be provided as a data.frame.")
  if(!is(SAMPLE, "character") | !all(SAMPLE %in% colnames(feature_table)))
    stop("Please provide a vector of sample names to be filtered. Sample names
         must appear in the feature table column names")
  if(!is(NEG1, "character") | !all(NEG1 %in% colnames(feature_table)) |
     !is(NEG2, "character") | !all(NEG2 %in% colnames(feature_table)))
    stop("Please provide NEG1 and/or NEG2 samples as character. Control sample 
         names must appear in the feature table column names")
  if(length(NEG1) == 0 & length(NEG2) == 0) 
    stop("Please provide a vector of control samples for at least one type
         of negative controls (NEG1 and/or NEG2)")
  if(length(intersect(NEG1, NEG2)) != 0)
    stop("NEG1 and NEG2 samples must not overlap. Please provide unique 
         control sample names")
  if(length(intersect(NEG1, SAMPLE)) != 0 | length(intersect(NEG2, SAMPLE)) != 0)
    stop("Control samples and real samples must not overlap. Please provide
         unique names for samples and controls")
  thresholds <- c(ratio_NEG1_threshold, ratio_NEG2_threshold, 
                  span_NEG1_threshold, span_NEG2_threshold)
  if(all(is.na(thresholds)))
    message("No thresholds set, no contamination filtering will be performed")
  else if(!is(thresholds, "numeric"))
    stop("Ratio and span thresholds must be numeric")
  if(length(ratio_NEG1_threshold) > 1 | length(ratio_NEG2_threshold) > 1 | 
     length(span_NEG1_threshold) > 1 | length(span_NEG2_threshold) > 1) 
    stop("Ratio and span values must be a single value")

  # Put together the relevant samples for filtering
  feature_table_formatted <- data.frame(feature_table[, c(SAMPLE, NEG1, NEG2)],
                                        check.names = FALSE)
  
  # Check for empty samples or NA values in the feature table
  if(any(colSums(feature_table_formatted) == 0) | any(is.na(feature_table_formatted)))
    stop("Please remove empty samples and NA values before running MicrobIEM decontamination")
  # Check whether the feature table needs to be transformed to relative abundances
  if(!all(colSums(feature_table_formatted) == 1)) {
    message("Normalise feature table to relative abundances")
    feature_table_formatted <- sweep(
      feature_table_formatted, 2, colSums(feature_table_formatted), "/")
  }

  # Calculate sample mean per feature
  sample_mean <- rowMeans(feature_table_formatted[, SAMPLE, drop = FALSE])
  # Calculate mean and span per feature in NEG1 
  if(!length(NEG1) == 0) {
    neg1_mean <- rowMeans(feature_table_formatted[, NEG1, drop = FALSE])
    neg1_span <- apply(feature_table_formatted[, NEG1, drop = FALSE], 1, 
                       function(x) sum(x > 0)/length(x))
  }
  # Calculate mean and span per feature in NEG2
  if(!length(NEG2) == 0) {
    neg2_mean <- rowMeans(feature_table_formatted[, NEG2, drop = FALSE])
    neg2_span <- apply(feature_table_formatted[, NEG2, drop = FALSE], 1,
                       function(x) sum(x > 0)/length(x))
  }
  # Merge the different pieces of information for each feature
  filter_basis <- data.frame(
    Feature.ID = row.names(feature_table_formatted),
    sample_mean = sample_mean,
    neg1_mean = if(exists("neg1_mean")) neg1_mean else NA,
    neg1_span = if(exists("neg1_span")) round(neg1_span, 4) else NA,
    neg2_mean = if(exists("neg2_mean")) neg2_mean else NA,
    neg2_span = if(exists("neg2_span")) round(neg2_span, 4) else NA)
  rownames(filter_basis) <- rownames(feature_table_formatted)
  # Calculate the ratios
  filter_basis["neg1_ratio"] <- 
    filter_basis$neg1_mean / filter_basis$sample_mean
  filter_basis["neg2_ratio"] <-
    filter_basis$neg2_mean / filter_basis$sample_mean
  
  # Convert filter criteria to numeric values
  if(all(is.na(ratio_NEG1_threshold))) ratio_NEG1_threshold <- Inf
  if(all(is.na(ratio_NEG2_threshold))) ratio_NEG2_threshold <- Inf
  if(all(is.na(span_NEG1_threshold))) span_NEG1_threshold <- 0.0001
  if(all(is.na(span_NEG2_threshold))) span_NEG2_threshold <- 0.0001
  
  # Round span filter
  span_NEG1_threshold <- round(span_NEG1_threshold, 4)
  span_NEG2_threshold <- round(span_NEG2_threshold, 4)
  
  # Apply filter criteria and return Sample IDs that should be removed
  feature_removed_neg1 <- character(0)
  if(ratio_NEG1_threshold == Inf && span_NEG1_threshold != 0.0001) {
    # NEG1: Use only span filter - remove everything that appears in controls
    feature_removed_neg1 <- row.names(
      filter_basis[filter_basis$neg1_span >= span_NEG1_threshold, ])
  } else { 
    # NEG1: Use no filter, only ratio filter, or ratio and span filter together
    feature_removed_neg1 <- row.names(
      filter_basis[filter_basis$neg1_span >= span_NEG1_threshold & 
                     filter_basis$neg1_ratio > ratio_NEG1_threshold, ])
  }
  features_removed_neg2 <- character(0)
  if(ratio_NEG2_threshold == Inf && span_NEG2_threshold != 0.0001) {
    # NEG2: Use only span filter - remove everything that appears in controls
    feature_removed_neg2 <- row.names(
      filter_basis[filter_basis$neg2_span >= span_NEG2_threshold, ])
  } else {
    # NEG1: Use no filter, only ratio filter, or ratio and span filter together
    feature_removed_neg2 <- row.names(
      filter_basis[(filter_basis$neg2_span >= span_NEG2_threshold &
                     filter_basis$neg2_ratio > ratio_NEG2_threshold), ])
  }
  # Add columns for information on contaminant status and removal 
  filter_basis["status"] <- ifelse(
    rownames(filter_basis) %in% feature_removed_neg1 &
      rownames(filter_basis) %in% feature_removed_neg2, 
    "removed (NEG1 & NEG2)", ifelse(
      rownames(filter_basis) %in% feature_removed_neg1,
      "removed (NEG1)", ifelse(
        rownames(filter_basis) %in% feature_removed_neg2,
        "removed (NEG2)", "kept")))
  filter_basis["is_contaminant"] <- grepl("removed", filter_basis[, "status"])
  return(filter_basis)
} 

