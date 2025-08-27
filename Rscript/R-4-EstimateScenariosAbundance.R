##############################################################
# R-4-EstimateScenariosAbundance.R
# Estimate local plant responses to projected temperature or 
# precipitation change, using KDE thresholds (from R-3).
#
# Dependencies:
#   - Requires:
#       - samples, tmp.otu, tmp.smp, tmp.tax (from R-1)
#       - plant.temp, plant.kde.funs (from R-3)
#   - Requires: measure <- "Temperature" or "Precipitation"
# Outputs:
#   - prediction.temperature.csv and prediction.precipitation.csv
##############################################################

print(paste("Estimating future scenarios for:", measure))

# Set climate change levels
if (measure == "Temperature") {
  change_levels <- seq(0, 5, by = 0.5)  # Temperature increase (°C)
} else if (measure == "Precipitation") {
  change_levels <- seq(0, 50, by = 5)   # Precipitation decrease (mm)
} else {
  stop("measure must be 'Temperature' or 'Precipitation'")
}


# Prepare parallel backend
cores <- detectCores()
cl <- makeCluster(cores[1] - 2)
registerDoSNOW(cl)

# Get all site IDs (sample names)
site_names <- sample_names(samples)

# Set up progress bar
pb <- txtProgressBar(0, length(change_levels), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)


# Run predictions in parallel
predictions <- foreach(change = change_levels, .combine = "rbind", .options.snow = opts,
                       .packages = c("phyloseq", "tidyr", "speedyseq", "dplyr")) %dopar% {
local({                        
 scenario_df <- data.frame()  # collects predictions for this delta
 
 for (site_id in site_names) {
   # Prepare site-specific phyloseq object and melt
   site_otu <- tmp.otu[, site_id]
   site_phy <- merge_phyloseq(site_otu, tmp.smp, tmp.tax)
   site_melt <- psmelt(site_phy)
   site_melt <- site_melt[site_melt$Abundance > 0, ]
   
   for (plant_id in site_melt$species) {
     # Get metadata for plant
     meta <- plant.temp[plant.temp$species == plant_id, ]
     if (nrow(meta) == 0) next
     
     # Build result row
     row <- data.frame(
       species = plant_id,
       site = site_id,
       date = site_melt$Month[1],
       country = site_melt$Country[1],
       coordY = site_melt$CoordY[1],
       abundance = site_melt[site_melt$species == plant_id, "Abundance"],
       prediction = NA,
       delta = change,
       crop = meta$crop
     )
     
     # Only evaluate if species is sufficiently sampled. Otherwise prediction will remain NA
     if (meta$sites2 > 5) {
       
       # Extract current climate value
       current_value <- if (measure == "Temperature") {
         site_melt$Temp[1]
       } else {
         site_melt$Prec[1]
       }
       
       # Apply change (positive for temp increase, negative for precip loss)
       future_value <- if (measure == "Temperature") {
         current_value + change
       } else {
         current_value - change
       }
       
       # Get KDE threshold values
       if (measure == "Temperature") {
         # Fallback if missing
         if (is.na(meta$temp.p99.kde)) meta$temp.p99.kde <- meta$temp.max
         
         if (meta$temp.p99.kde < future_value) {
           row$prediction <- ">p99"
         } else if (meta$temp.p95.kde < future_value) {
           row$prediction <- ">p95"
         } else if (meta$temp.p90.kde < future_value) {
           row$prediction <- ">p90"
         } else {
           row$prediction <- "<=p90"
         }
         
       } else if (measure == "Precipitation") {
         # Fallback if missing
         if (is.na(meta$prec.p01.kde)) meta$prec.p01.kde <- meta$prec.min
         
         if (meta$prec.p01.kde > future_value) {
           row$prediction <- "<p01"
         } else if (meta$prec.p05.kde > future_value) {
           row$prediction <- "<p05"
         } else if (meta$prec.p10.kde > future_value) {
           row$prediction <- "<p10"
         } else {
           row$prediction <- ">=p10"
         }
       }
     }
     
     scenario_df <- rbind(scenario_df, row)
   }
 }
 
 return(scenario_df)
 })
}

# Write output as file
out_path <- paste0("intermediate.data/prediction.",measure,".csv")
write.table(predictions, file = out_path, sep = ",", row.names = FALSE)

print(paste0("Finished writing predictions for: ", measure))
