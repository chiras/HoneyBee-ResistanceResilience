# R code for Quaresma et al. 2024

## Citation
**Honey bee food resources under threat with climate change**

Andreia Quaresma, Hans Baveco, Robert Brodschneider, Bas Buddendorf, Norman Carreck, Kristina Gratzer, Fani Hatjina, Ole Kilpinen, Ivo Roessink, Flemming Vejsnaes, Jozef van der Steen, M. Alice Pinto, Alexander Keller

(in submission)


## Calling the script
Main call of the R code is in the ```R-0-Umbrella.R``` file, where all other files are sourced. 

The data folder holds files that need unzipping before the analyses (done in the script though).


## Dependencies

R packages loaded in ```R-0-Umbrella.R``` are 

```R
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
```

## Computational requirements

Consider these analyses are quite intense and require at least 20GB of RAM, with > 24GB recommended. Most of the computational intensive tasks are parallelized, which means the more CPUs are available, the faster the processes. The function ``` makeCluster(cores[1]-2) ``` is used to estimate cores available, and will use all -2 of those.

Analyses for the article were performed on MacOSX 13.6.7 with a MacBook Pro Max M1 with 64 GB of RAM and R version 4.4.1 (2024-06-14).

## Metabarcoding preprocessing: 

Preprocessing of metabarcoding was done according to: https://github.com/chiras/metabarcoding_pipeline
Some parts of the R Script also origin from this GitHub repository. 