################################################################################
#
# MicrobIEM - an algorithm to reduce contaminants in microbiome 16S amplicon 
#   sequencing data
#
################################################################################

# Download the MicrobIEM test data from Github
# Set the directory to where you stored the test data
file_directory <- c("C:/Users/.../test-data/") # CHANGE ME

# Read example feature file
MicrobIEM_featurefile <- read.delim(
  paste0(file_directory, "MicrobIEM_test-data_featurefile.txt"),
  row.names = 1)

# EXAMPLE 1 - Run MicrobIEM decontamination function only for NEG2 filtering
View(MicrobIEM_decontamination(
  MicrobIEM_featurefile, 
  SAMPLE = colnames(MicrobIEM_featurefile)[1:29], 
  NEG2 = c("Sample_34", "Sample_35"), 
  ratio_NEG2_threshold = 2, span_NEG2_threshold = 1, 
))

# EXAMPLE 2 - Run MicrobIEM decontamination function for NEG1 & NEG2 filtering
View(MicrobIEM_decontamination(
  MicrobIEM_featurefile, 
  SAMPLE = colnames(MicrobIEM_featurefile)[1:29], 
  NEG1 = c("Sample_31", "Sample_32", "Sample_33"), 
  NEG2 = c("Sample_34", "Sample_35"), 
  ratio_NEG1_threshold = 1, ratio_NEG2_threshold = 1,
  span_NEG1_threshold = NA, span_NEG2_threshold = NA
))
