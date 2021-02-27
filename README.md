# MicrobIEM - draft

### Content
[1. Overview](#Overview)  
[2. Prerequisites](#Prerequisites)  
[3. Start the tool](#Start)  
[4. Prepare your data](#Prepare)  
[5. Perform quality control](#Qualitycontrol)  
[6. Explore your data](#Explore)  
[7. Select samples of interest](#Select)  
[8. Save your figures](#Savefigures)  
[9. Save your results](#Saveresults)

### 1. Overview {#Overview}

<img src="MicrobIEM/man/01_Interface.png"/> 

MicrobIEM is a user-friendly tool for quality control and interactive analysis of microbiome data. A feature table and a metafile of a microbiome study are loaded via the graphical user interface. For each step in quality control, interactive visualisations allow users to explore their data and help defining thresholds for filtering the data. The final data set can then be further investigated with statistical analysis common in microbiome research. Raw data, figures and corresponding p-values can be downloaded in the end. Thus, MicrobIEM is a fast tool that allows the user to explore microbiome data without any knowledge in coding.

#### Available quality control steps:

- Sample filter: by number of total reads
- Feature filter: by maximum absolute abundance and maximum relative abundance
- Contaminant filter: by abundance in negative controls compared to samples

#### Available statistical analyses:
- Alpha diversity (Richness, Shannon, Simpson, Inverse Simpson, Evenness)
- Beta diversity (PCoA, nMDS, on Bray-Curtis dissimilarities)
- Taxonomy analysis

### 2. Prerequisites {#Prerequisites}
#### Download R
Download and install the software package R from the [R project website](https://cran.r-project.org/bin/windows/base/).

#### Download MicrobIEM
Download MicrobIEM and save the folder on your machine. Unzip the folder and copy its directory.  

### 3. Start the tool {#Start}
Open R and paste the following code chunk:
``` r
setwd("C:/user/.../MicrobIEM") # Change to the directory of your MicrobIEM folder
source("start_MicrobIEM.R")
start_MicrobIEM()
```
The last item in the directory (first line) should be the unzipped MicrobIEM folder.  
When you start MicrobIEM for the first time, this step may take some minutes because additional packages may need to be installed.

### 4. Prepare your data {#Prepare}
MicrobIEM requires the input data - a feature file and a meta file - to be in a specific format.

#### Featurefile
The feature file can be an OTU table or an ASV table, and contains sequenced read counts for each sample. It should be .txt or .csv file. The first column must be called "OTU_ID" and must contain unique names of features. The next columns start with the name of each sample. The last column must be called "Taxonomy" and contain information on taxonomic classification.
<img src="MicrobIEM/man/02_Featurefile.png"/> 


#### Metafile
The metafile contains additional information on each sample. It should be .txt or .csv file. The first column must be called "Sample_ID" and contains the same sample names that are found in the featurefile. One column in the metafile must be called "Sample_type" and contain the classification of samples into real samples and positive and negative controls. Please use the following terms to define samples and controls:

- "SAMPLE" for real samples
- "NEG1" or "NEG2" for up to 2 different types of negative controls
- "POS1" for positive controls

<img src="MicrobIEM/man/03_Metafile.png"/> 

### 5. Perform quality control {#Qualitycontrol} 
<img src="MicrobIEM/man/04_Contaminant-removal.png"/> 

### 6. Explore your data {#Explore}
<img src="MicrobIEM/man/05_Beta-diversity.png"/> 

### 7. Select samples of interest {#Select}
<img src="MicrobIEM/man/06_Sample-selection.png"/> 

### 8. Save your figures {#Savefigures}
<img src="MicrobIEM/man/07_Save-figures.png"/> 

### 9. Save your results {#Saveresults}
<img src="MicrobIEM/man/08_Save-results.png"/> 
