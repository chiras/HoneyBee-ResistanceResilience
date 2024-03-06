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
library(geosphere)
library(scales)
library(foreach)
library(doParallel)
library(doSNOW)
library(progress)

## setting directories, variables and custom functions
setwd('/Users/ra39huv/TMP/Data_processing/ITS2_Insignia.24-01-30/TestFullPipeline')
dir.create("plots")
dir.create("plots.supplement")

marker="ITS2"
source('Rscript/metabarcoding_tools_0-1a.R')

# unzip data 
unzip("Data/total.asv_table.bytable.sed.txt.zip", overwrite = F, exdir="Data")
unzip("Data/total.taxonomy.sed.vsearch.zip", overwrite = F, exdir="Data")

# run code for data preparations and graphs
source("Rscript/R-1-Preprocessing.R")
source("Rscript/R-2-DiversityMap.R", print.eval=TRUE)
source("Rscript/R-3-EstimateDistributions.R")
source("Rscript/R-4-EstimateScenariosAbundance.R")

toAnalyze = "results.temp" 
source("Rscript/R-5-EstimateScenariosAbundanceII.R", print.eval=TRUE)
source("Rscript/R-6-TemporalProgression.R", print.eval=TRUE)

toAnalyze = "results.prec" 
source("Rscript/R-5-EstimateScenariosAbundanceII.R", print.eval=TRUE)
source("Rscript/R-6-TemporalProgression.R", print.eval=TRUE)

source("Rscript/R-7-ResistanceResilience.R", print.eval=TRUE)

# stats

