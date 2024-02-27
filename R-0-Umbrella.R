################################################
# R data analysis for metabarcoding
# by Alexander Keller (LMU München)
# keller@bio.lmu.de
#
# Created: Thu Nov 16 19:39:58 CET 2023
# Project: ITS2_Insignia
# Marker: ITS2
# For: INSIGNIA
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

marker="ITS2"
source('Rscript/metabarcoding_tools_0-1a.R')

# run code
source("Rscript/R-1-Preprocessing.R")
source("Rscript/R-2-DiversityMap.R", print.eval=TRUE)
source("Rscript/R-3-EstimateDistributions.R")
source("Rscript/R-4-EstimateScenariosAbundance.R")

toAnalyze = "results.temp" 
source("Rscript/R-5-EstimateScenariosAbundanceII.R", print.eval=TRUE)

toAnalyze = "results.prec" 
source("Rscript/R-5-EstimateScenariosAbundanceII.R", print.eval=TRUE)
