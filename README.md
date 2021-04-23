---
title: "MicrobIEM"
output:
  html_document:
    toc: true
---

## 1. Overview

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

## 2. Prerequisites
#### Download R
Download and install the software package R from the [R project website](https://cran.r-project.org/bin/windows/base/). If you already used R on your machine, please update it to at least version 4.0.

#### Download RStudio
Download and install RStudio, an integrated development environment (IDE) for R, from the [RStudio website](https://rstudio.com/products/rstudio/download/).

#### Install Shiny
Once R and RStudio are properly installed, open RStudio and type or copy the following command in the console:
``` r
install.packages("shiny")
```
Pressing enter starts the installation.

#### Download MicrobIEM
Download MicrobIEM and save and unzip the folder on your machine. 

You can download MicrobIEM by clicking the green "Code" button in the top right of this repository and select "Download ZIP".  

## 3. Start the tool
In your unzipped MicrobIEM folder, open the file 'server' by doubleclicking - it should automatically open in RStudio. Press the 'Run App' button in the upper middle to start MicrobIEM.  
When you start MicrobIEM for the first time, this step may take some minutes because additional packages may need to be installed.
<img src="MicrobIEM/man/09_Start-MicrobIEM.png"/> 

## 4. Prepare your data
MicrobIEM requires the input data - a feature file and a meta file - to be in a specific format. You can see an example of the required formats in the folder /MicrobIEM-main/MicrobIEM/Test-Data.

#### Featurefile
The feature file can be an OTU table or an ASV table, and contains sequenced read counts for each sample. It should be .txt or .csv file. The first column must be called "OTU_ID" and must contain unique names of features. The next columns start with the name of each sample. The last column must be called "Taxonomy" and contain information on taxonomic classification.
<img src="MicrobIEM/man/02_Featurefile.png"/> 

#### Metafile
The metafile contains additional information on each sample. It should be .txt or .csv file. The first column must be called "Sample_ID" and contains the same sample names that are found in the featurefile. One column in the metafile must be called "Sample_type" and contain the classification of samples into real samples and positive and negative controls. Please use the following terms to define samples and controls:

- "SAMPLE" for real samples
- "NEG1" or "NEG2" for up to 2 different types of negative controls
- "POS1" for positive controls

<img src="MicrobIEM/man/03_Metafile.png"/> 

## 5. Perform quality control 
<img src="MicrobIEM/man/04_Contaminant-removal.png"/> 

## 6. Explore your data
<img src="MicrobIEM/man/05_Beta-diversity.png"/> 

## 7. Select samples of interest
<img src="MicrobIEM/man/06_Sample-selection.png"/> 

## 8. Save your figures
<img src="MicrobIEM/man/07_Save-figures.png"/> 

## 9. Save your results
<img src="MicrobIEM/man/08_Save-results.png"/> 
