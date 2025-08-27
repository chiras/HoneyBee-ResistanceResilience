##############################################################
# R-7-ResistanceResilience.R
# Combined effects: estimating site-level resistance and 
# resilience to climate extremes using KDEs.
#
# Dependencies:
#   - Requires:
#       - samples, tmp.otu, tmp.smp, tmp.tax from R-1
#       - plant.temp, plant.kde.funs from R-3
#       - data.map, data.species, plant_taxa
# Outputs:
#   - In-memory object: prediction.resiresi (used in R-7.1)
##############################################################
# For each species, compute the probability of occurrence under changing climatic conditions
# Using the cumulative distribution function (CDF) from plant-specific parameters

# Determine resistance: Fraction of currently present plants likely to persist
# Determine resilience: Fraction of potential replacement species satisfying spatial, climatic, and taxonomic constraints

print("Calculating site-level Resistance and Resilience")


# -------------------------
# Setup and definitions
# -------------------------

cores <- detectCores()
cl <- makeCluster(cores[1] - 2)
registerDoSNOW(cl)

temp_levels <- seq(0, 5, 0.5)
prec_levels <- seq(0, 50, 5)

# Distance thresholds
min.dist.geo <- 0.35   # ~500 km
min.prob <- 0.15       # minimum KDE probability (both temp and prec)
min.tax <- 0.5         # max taxonomic distance (same family or closer)

# Extract plant KDE functions
kde_temp_funs <- lapply(plant.kde.funs, `[[`, "temp_fun")
kde_prec_funs <- lapply(plant.kde.funs, `[[`, "prec_fun")
names(kde_temp_funs) <- names(plant.kde.funs)
names(kde_prec_funs) <- names(plant.kde.funs)

# -------------------------
# Precompute taxonomic distances between species
# -------------------------

print("Precompute plant taxonomic distances")
plant.df <- merge(plant.temp, data.frame(tax_table(samples)), by="species")      
rownames(plant.df) <- plant.df$species
plant.df$endemism <- rescale(sqrt(1/plant.df$countries), c(0,1))


plant_taxa2 <- plant_taxa

pb <- txtProgressBar(0, length(plant_taxa2), style = 3)
distance.tax <-foreach (plant = plant_taxa2, .combine=cbind, .options.snow=opts, .packages=c("geosphere","phyloseq","scales")) %dopar% {
  print(plant)
  plant.qs = matrix(nrow = length(plant_taxa2), ncol=1)
  colnames(plant.qs)[1] <- plant
  rownames(plant.qs) <- plant_taxa2
  
  for (plant.replacement in plant_taxa2){
    tax.dist <- NA
    if (!is.na(plant.df[plant.replacement,"genus"]) && !is.na(plant.df[plant,"genus"])){ 
      if (plant == plant.replacement){
        tax.dist <- 0
      }else if (plant.df[plant,"genus"] == plant.df[plant.replacement,"genus"]){
        tax.dist <- 0.2
      }else if (plant.df[plant,"family"] == plant.df[plant.replacement,"family"]){
        tax.dist <- 0.4
      }else if (plant.df[plant,"order"] == plant.df[plant.replacement,"order"]){
        tax.dist <- 0.6
      }else if (plant.df[plant,"phylum"] == plant.df[plant.replacement,"phylum"]){
        tax.dist <- 0.8
      }else{
        tax.dist <- 1
      }
      plant.qs[plant.replacement,1] <- tax.dist
    }
  }  
  plant.qs
}

# -------------------------
# Precompute geographic distances between sites
# -------------------------

print("Precompute geographic distances between sampling sites and species association")

geodist.raw <- distm(sample_data(samples)[, c("CoordX", "CoordY")], fun = distHaversine)
geodist <- rescale(sqrt(geodist.raw), to = c(0, 1))  # rescale to 0–1
rownames(geodist) <- rownames(sample_data(samples))
colnames(geodist) <- rownames(sample_data(samples))

### Verify 500km radius
max_d <- max(geodist.raw)
threshold_m <- (min.dist.geo) * (sqrt(max_d) - sqrt(0)) + sqrt(0)
threshold_m <- threshold_m^2   # undo sqrt
threshold_km <- threshold_m / 1000

# Create site identifiers
data.map$site <- interaction(data.map$CoordY, data.map$CoordX, sep = "-")
sites <- unique(rownames(sample_data(samples)))

pb <- txtProgressBar(0, length(sites), style = 3)

distance.geo <- foreach(sample=sites, .options.snow=opts , .combine="cbind",  .packages=c("phyloseq","tidyr","speedyseq","dplyr")) %dopar% {
  plant.qs = matrix(nrow = length(plant_taxa), ncol=1)
  colnames(plant.qs)[1] <- sample
  rownames(plant.qs) <- plant_taxa
  
  for (plant in plant_taxa){
    sites.present1 <- colnames(otu_table(samples)[,otu_table(samples)[plant,]>0])
    min.distance <- min(geodist[sites.present1,sample])
    
    plant.qs[plant,1] <- min.distance
  }
  plant.qs
}

# -------------------------
# Compute site-level resistance and resilience
# -------------------------

print("Calculating Resistance and Resilience")

####### Using now full grid approach for less for loops

# Create full parameter grid
climate_grid <- expand.grid(
  temp_increase = temp_levels,
  prec_increase = prec_levels
)

pb <- txtProgressBar(0, nrow(climate_grid), style = 3)
prediction.resiresi <- foreach(i = 1:nrow(climate_grid), .combine = "rbind", 
                               .packages = c("phyloseq", "dplyr", "geosphere", "scales")) %dopar% {
                                 
    temp_increase <- climate_grid$temp_increase[i]
    prec_increase <- climate_grid$prec_increase[i]
  
    print(paste0(temp_increase,"--",prec_increase))

    .GlobalEnv$samples <- samples
    .GlobalEnv$sample_names <- sample_names

    site_predictions <- data.frame()
    
    for (site_id in sites) {
      site.otu <- tmp.otu[, site_id]
      site.phyloseq <- merge_phyloseq(site.otu, tmp.smp, tmp.tax)
      site.melt <- psmelt(site.phyloseq)
      site.melt <- site.melt[site.melt$Abundance > 0, ]
      plants.in.plot <- site.melt$species
      
      # Compute expected climate values for site under change
      temp_future <- site.melt$Temp[1] + temp_increase
      prec_future <- max(0, site.melt$Prec[1] - prec_increase) # lower bound at 0
      
      # Use already set KDE intervals
      thresholds  <- plant.temp[plant.temp$species %in% plants.in.plot, c("species","temp.p90.kde","prec.p10.kde")]
      thresholds$temp.fallout <- thresholds$temp.p90.kde < temp_future
      thresholds$prec.fallout <- thresholds$prec.p10.kde > prec_future
      
      # Resistance = plants that remain viable
      resistance.plants <- thresholds[(thresholds$temp.fallout + thresholds$prec.fallout) == 0,"species"]
      dropout.plants <- thresholds[(thresholds$temp.fallout + thresholds$prec.fallout) > 0,"species"]
      resistance <- length(resistance.plants) / length(plants.in.plot)
      
      if (resistance == 1) {
        resilience=0
        resilience_raw=0
      }else{
      
      # Candidates for replacement from nearby sites (excluding same site)
      rep.geo <- distance.geo[,site_id]
      candidate_species <- names(rep.geo[rep.geo < min.dist.geo & rep.geo > 0])
      
      # Filter viable candidates using KDEs
      candidate.thresholds  <- plant.temp[plant.temp$species %in% candidate_species, c("species","temp.p90.kde","prec.p10.kde")]
      candidate.thresholds$temp.fallout <- candidate.thresholds$temp.p90.kde < temp_future
      candidate.thresholds$prec.fallout <- candidate.thresholds$prec.p10.kde > prec_future
      candidate.thresholds <- candidate.thresholds %>%
        filter(!temp.fallout) %>%
        filter(!prec.fallout)
      
      if (length(candidate.thresholds[,1]==0)){
        resilience = 0
        resilience_raw = 0
      }
        
      # not already in plot
      candidate.thresholds <- candidate.thresholds %>%
        filter(!species %in% plants.in.plot)
            
      # Taxonomic resilience among dropouts
      replacements <- data.frame(original=dropout.plants, replacement=NA)

      for (dropout in replacements$original) {
        if (!dropout %in% rownames(distance.tax)) next
        tax_distances <- distance.tax[dropout, candidate.thresholds$species, drop = FALSE]
        related_species <- colnames(tax_distances)[tax_distances < min.tax]
        if (length(related_species) > 0) {
          replacements$replacement[replacements$original==dropout]<-related_species[1]
          candidate.thresholds <- candidate.thresholds %>%
            filter(species != related_species[1])
        }
      }
      
      # Resilience_raw = proportion of viable replacements (used)
      # Resilience = proportion of viable replacements rescaled by proportion dropped out (not used)
      
      resilience <- (sum(!is.na(replacements$replacement)) / length(replacements$replacement))  * (1 - resistance)
      resilience_raw <- (sum(!is.na(replacements$replacement)) / length(replacements$replacement))
      } # end else 100% resistance
      
      site_predictions <- rbind(site_predictions, data.frame(
        site = site_id,
        richness = length(plants.in.plot),
        prec_increase = prec_increase,
        temp_increase = temp_increase,
        resistance = resistance,
        resilience_raw = resilience_raw,
        resilience = resilience
      ))
    } # end site loop
    
    
    return(site_predictions)
  }

stopCluster(cl)

print("Resistance and Resilience calculations complete")

#### Metadata aggregation
metadata <- sample_data(data.species)
metadata$site <- rownames(metadata)
prediction.all <- merge(prediction.resiresi,data.frame(metadata), by="site")

#### Writing to file
write.table(
  prediction.all,
  file = "intermediate.data/sites_resi-resi.csv",
  sep = ",",
  row.names = FALSE
)
