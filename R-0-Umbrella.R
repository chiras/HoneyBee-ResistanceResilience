################################################
# R analysis for metabarcoding data
# by Alexander Keller (LMU München)
# keller@bio.lmu.de
#
# Created: Thu Nov 16 19:39:58 CET 2023
# Project: ITS2_HoneyBee-ResiResi
# Marker: ITS2
################################################
# Bioinformatics follow:
# https://github.com/chiras/metabarcoding_pipeline
# https://doi.org/10.1098/rstb.2021.0171
################################################

# Preparation
## clear workspace
rm(list = ls())

## load libraries
library(phyloseq)
library(tidyr)
library(speedyseq)
library(dplyr)
library(viridis)
library(bipartite)
library(ggplot2)
library(ggsci)
library(rnaturalearth)
library(sf)
library(ggplot2)
library(ggrepel)
library(geosphere)
library(scales)
library(foreach)
library(doParallel)
library(doSNOW)
library(progress)

## setting directories, variables and custom functions
setwd('/Users/ra39huv/TMP/Data_processing/ITS2_Insignia.24-01-30/Test_2/HoneyBee-ResistanceResilience')
dir.create("plots")
dir.create("plots.supplement")

marker="ITS2"
source('Rscript/metabarcoding_tools_0-1a.R')

# unzip data 
unzip("Data/total.asv_table.bytable.sed.txt.zip", overwrite = F, exdir="Data")
unzip("Data/total.taxonomy.sed.vsearch.zip", overwrite = F, exdir="Data")

# run code for data preparations and graphs
test_results=list()

# Prepare intermediate data and estimate distributions
source("Rscript/R-1-Preprocessing.R")
source("Rscript/R-2-DiversityMap.R", print.eval=TRUE)
source("Rscript/R-3-EstimateDistributions.R", print.eval=TRUE)
source("Rscript/R-4-EstimateScenariosAbundance.R")

# Analyse risks
toAnalyze = "results.temp" 
source("Rscript/R-5-EstimateScenariosAbundanceII.R", print.eval=TRUE)
source("Rscript/R-6-TemporalProgression.R", print.eval=TRUE)

toAnalyze = "results.prec" 
source("Rscript/R-5-EstimateScenariosAbundanceII.R", print.eval=TRUE)
source("Rscript/R-6-TemporalProgression.R", print.eval=TRUE)

# Calculate Resistance and Resilience 
source("Rscript/R-7-ResistanceResilience.R", print.eval=TRUE)
source("Rscript/R-7.1-ResistanceResiliencePlot.R", print.eval=TRUE)

# Historical changes 
measure="Temperature"
source("Rscript/R-8-HistoricalChanges.R", print.eval=TRUE)

measure="Precipitation"
source("Rscript/R-8-HistoricalChanges.R", print.eval=TRUE)

measure="Precipitation"
source("Rscript/R-9_GDMs.R", print.eval=TRUE)

# Further statistics

names(test_results)

# Historical Temperature min and max
# southern
mean(test_results$Temperature_historical_change_since_1973_within[test_results$Temperature_historical_change_since_1973_within$CoordY.x  < 45,"mean"])
max(mean(test_results$Temperature_historical_change_since_1973_within[test_results$Temperature_historical_change_since_1973_within$CoordY.x  < 45,"max"]))

# central
mean(test_results$Temperature_historical_change_since_1973_within[test_results$Temperature_historical_change_since_1973_within$CoordY.x  > 45 & test_results$Temperature_historical_change_since_1973_within$CoordY.x < 54,"mean"])
max(mean(test_results$Temperature_historical_change_since_1973_within[test_results$Temperature_historical_change_since_1973_within$CoordY.x  > 45 & test_results$Temperature_historical_change_since_1973_within$CoordY.x < 54,"max"]))

# northern
mean(test_results$Temperature_historical_change_since_1973_within[test_results$Temperature_historical_change_since_1973_within$CoordY.x  > 54,"mean"])
max(mean(test_results$Temperature_historical_change_since_1973_within[test_results$Temperature_historical_change_since_1973_within$CoordY.x  > 54,"max"]))

# Precipitation 
# southern
diff(range(test_results$Precipitation_historical_change_since_1973_within[test_results$Precipitation_historical_change_since_1973_within$CoordY.x  < 45,c("min","max")]))
min(test_results$Precipitation_historical_change_since_1973_within[test_results$Precipitation_historical_change_since_1973_within$CoordY.x  < 45,"min"])

# central
diff(range(test_results$Precipitation_historical_change_since_1973_within[test_results$Precipitation_historical_change_since_1973_within$CoordY.x  > 45 & test_results$Precipitation_historical_change_since_1973_within$CoordY.x < 54,c("min","max")]))
min(test_results$Precipitation_historical_change_since_1973_within[test_results$Precipitation_historical_change_since_1973_within$CoordY.x  > 45 & test_results$Precipitation_historical_change_since_1973_within$CoordY.x < 54,"min"])

# northern
diff(range(test_results$Precipitation_historical_change_since_1973_within[test_results$Precipitation_historical_change_since_1973_within$CoordY.x  > 54,c("min","max")]))
min(test_results$Precipitation_historical_change_since_1973_within[test_results$Precipitation_historical_change_since_1973_within$CoordY.x  > 54,"min"])


