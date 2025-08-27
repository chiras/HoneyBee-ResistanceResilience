################################################
# Based on: 
# R analysis for metabarcoding data
# Author: Alexander Keller (LMU München) | keller@bio.lmu.de
# Marker: ITS2
# Bioinformatics: https://github.com/chiras/metabarcoding_pipeline
# Ref: https://doi.org/10.1098/rstb.2021.0171
################################################
# Project: ITS2_HoneyBee-ResiResi
# Created: Thu Nov 16 19:39:58 CET 2023
################################################

# ------------------------------------------------
# Workspace setup
# ------------------------------------------------
rm(list = ls())

# Set path if not "here":
#setwd('/PATH/TO/SCRIPT/ON/YOUR/SYSTEM') 

# Load libraries (deduplicated and ordered)
library(phyloseq); library(speedyseq); library(tidyr); library(dplyr)
library(ggplot2); library(ggrepel); library(ggnewscale); library(ggsci)
library(gghighlight); library(viridis); library(scales); library(patchwork)
library(rnaturalearth); library(sf); library(geosphere); library(bipartite); 
library(gdm); library(vegan); library(foreach); library(doParallel); 
library(doSNOW); library(progress)

# Create directories
dir.create("plots", showWarnings = FALSE)
dir.create("plots.supplement", showWarnings = FALSE)
dir.create("intermediate.data", showWarnings = FALSE)

# ------------------------------------------------
# Load custom functions & data
# ------------------------------------------------
marker <- "ITS2"
source('Rscript/metabarcoding_tools_0-1a.R')

# Unzip data if necessary
unzip("Data/total.asv_table.bytable.sed.txt.zip", overwrite = FALSE, exdir = "Data")
unzip("Data/total.taxonomy.sed.vsearch.zip", overwrite = FALSE, exdir = "Data")

# ------------------------------------------------
# Analysis pipeline
# ------------------------------------------------

### Results collection
test_results=list()

#### R-1: Preprocessing ###########################################
source("Rscript/R-1-Preprocessing.R")

#### R-2: Diversity & Crops (depends on R-1) ######################
source("Rscript/R-2-DiversityMap.R", print.eval = TRUE)

#### R-3: KDE Estimation (depends on R-1) #########################
source("Rscript/R-3-Estimate_KDE.R", print.eval = TRUE)

#### R-4: Scenario Abundance Estimation (depends on R-3) ##########
# only needed to be run once for creating temporary files

# First for Temperature
measure <- "Temperature"
source("Rscript/R-4-EstimateScenariosAbundance.R") 

# Second for Precipitation
measure <- "Precipitation"
source("Rscript/R-4-EstimateScenariosAbundance.R") 

#### R-5/6: Climate Risk Projections (depends on one-time execution of R-4) #############
# First for Temperature
measure <- "Temperature"
source("Rscript/R-5-EstimateScenariosAbundanceII.R", print.eval = TRUE)
source("Rscript/R-6-TemporalProgression.R", print.eval = TRUE)

# Second for Precipitation
measure <- "Precipitation"
source("Rscript/R-5-EstimateScenariosAbundanceII.R", print.eval = TRUE)
source("Rscript/R-6-TemporalProgression.R", print.eval = TRUE)

#### R-7: Resistance & Resilience (depends on R-4) ################
source("Rscript/R-7-ResistanceResilience.R", print.eval = TRUE)     ####RE-CHECK
source("Rscript/R-7.1-ResistanceResiliencePlot.R", print.eval = TRUE)

#### R-8: Historical Changes (self-contained) #####################
# First for Temperature
measure <- "Temperature"
source("Rscript/R-8-HistoricalChanges.R", print.eval = TRUE)

# Second for Precipitation
measure <- "Precipitation"
source("Rscript/R-8-HistoricalChanges.R", print.eval = TRUE)

# Statistics for both
source("Rscript/R-8.1-HistoricalChangesStats.R", print.eval = TRUE)

#### R-9: GDM & NMDS Ordination (depends on R-1) ##################
source("Rscript/R-9_GDMs.R", print.eval = TRUE)

# ------------------------------------------------
# Summary statistics (reported in manuscript)
# ------------------------------------------------
names(test_results)
