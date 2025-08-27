##############################################################
# R-3-Estimate_KDE.R
# Estimates kernel density distributions for all plant taxa
# Dependencies:
#   - Requires `samples` object from R-1-Preprocessing
# Outputs:
#   - intermediate.data/plant_metadata_complete.csv
#   - intermediate.data/plant_metadata_qualified.csv
#   - In-memory: `plant.kde.funs` used by R-4 through R-7.1
##############################################################
oldw <- getOption("warn")
options(warn = -1) # suppress warnings temporarily

print("Estimating plant temperature and precipitation distributions...")

# ------------------------------------------------------------
# Quick diagnostics: distribution of climate in sample sites
# ------------------------------------------------------------

pdf("plots.supplement/temperature_samples.pdf", width = 7, height = 6.5)
hist(data.map$Temp)
dev.off()

pdf("plots.supplement/precipitation_samples.pdf", width = 7, height = 6.5)
hist(data.map$Prec)
dev.off()


# ------------------------------------------------------------
# Prepare for parallel computation of KDEs
# ------------------------------------------------------------

cores <- detectCores()
cl <- makeCluster(cores[1] - 2)  # Leave 1-2 cores free
registerDoSNOW(cl)

# Extract pieces from phyloseq object, since it can't be passed across cores
tmp.otu <- otu_table(samples)
tmp.smp <- sample_data(samples)
tmp.tax <- tax_table(samples)
plant_taxa <- taxa_names(samples)

# Progress bar
pb <- txtProgressBar(0, length(plant_taxa), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

# ------------------------------------------------------------
# Parallel loop over plant species: KDE estimation
# ------------------------------------------------------------

plant.Temp <- foreach(
  plant_i = seq_along(plant_taxa), .combine = 'c', .options.snow = opts,
  .packages = c("phyloseq", "dplyr")
) %dopar% {
  
  plant <- plant_taxa[plant_i]
  target.otu <- tmp.otu[plant, ]
  target <- merge_phyloseq(target.otu, tmp.smp, tmp.tax)
  target.melt <- psmelt(target)
  abundance.plant <- sum(target.otu) / length(target.otu)
  
  # Skip species with no observations
  if (sum(target.melt$Abundance) == 0) return(NULL)
  
  # Filter only non-zero entries and define site IDs
  target.melt <- target.melt[target.melt$Abundance > 0, ]
  target.melt$location <- interaction(target.melt$CoordY, target.melt$CoordX)
  
  # Metadata for species
  
  crop <- plant %in% crops$V1 # Crop status
  countries <- length(unique(target.melt$Country)) # Number of countries where plant occurs
  sites <- length(unique(target.melt$location)) # Number of unique spatial sites (lat/lon combinations)
  sites2 <- length(unique(target.melt$Sample)) # Number of unique samples where this plant was recorded
  
  # Aggregate abundance by temperature and precipitation
  target.temp <- target.melt %>%
    group_by(Temp) %>%
    summarize(accumulated = sum(Abundance), .groups = "drop")
  
  target.prec <- target.melt %>%
    group_by(Prec) %>%
    summarize(accumulated = sum(Abundance), .groups = "drop")
  
  # Generate distribution samples
  dist.temp <- unlist(mapply(rep, target.temp$Temp, round(target.temp$accumulated, 3) * 1000))
  dist.prec <- unlist(mapply(rep, target.prec$Prec, round(target.prec$accumulated, 3) * 1000))
  
  # ----- KDE: Temperature -----
  density_temp <- density(dist.temp, from = 0, to = 60, bw = 1) #2
  dx_temp <- mean(diff(density_temp$x))
  
  # Area-normalized: for computing cumulative distribution
  density_temp_area_y <- density_temp$y / sum(density_temp$y * dx_temp)
  cdf_temp_y <- cumsum(density_temp_area_y) * dx_temp
  cdf_temp <- approxfun(density_temp$x, cdf_temp_y, rule = 2)
  density_temp_fun_cdf <- approxfun(density_temp$x, density_temp$y / sum(density_temp$y * dx_temp), rule = 2)
  
  # Peak-normalized: for visual plotting
  density_temp_fun <- approxfun(density_temp$x, density_temp$y / sum(density_temp$y), rule = 2)
  
  # Estimate percentiles
  percentiles_temp_kde <- approx(cdf_temp_y, density_temp$x, xout = c(0.90, 0.95, 0.99))$y
  
  # ----- KDE: Precipitation -----
  density_prec <- density(dist.prec, from = 0, to = 500, bw = 25)
  dx_prec <- mean(diff(density_prec$x))
  
  density_prec_area_y <- density_prec$y / sum(density_prec$y * dx_prec)
  cdf_prec_y <- cumsum(density_prec_area_y) * dx_prec
  cdf_prec <- approxfun(density_prec$x, cdf_prec_y, rule = 2)
  density_prec_fun_cdf <- approxfun(density_prec$x, density_prec$y / sum(density_prec$y * dx_prec), rule = 2)
  
  density_prec_fun <- approxfun(density_prec$x, density_prec$y / sum(density_prec$y), rule = 2)
  
  percentiles_prec_kde <- approx(cdf_prec_y, density_prec$x, xout = c(0.01, 0.05, 0.10))$y
  
  # Return species metadata + functions
  list(
    list(
      data = data.frame(
        species = plant,                     # Species name (taxon ID from the phyloseq object)
        n_records = nrow(target.melt),       # Total number of non-zero records (sample × abundance > 0)
        crop = crop,                         # Logical: TRUE if species is in the 'crops' list
        countries = countries,               # Number of unique countries in which the species was detected
        sites = sites,                       # Number of unique spatial sampling locations (based on CoordX × CoordY)
        sites2 = sites2,                     # Number of unique samples in which the species was present (Sample-level resolution)
        abundance = abundance.plant,         # Mean relative abundance across all samples (including zeros)
        
        # KDE-derived temperature percentiles (from cumulative density function)
        temp.p90.kde = percentiles_temp_kde[1],  # 90th percentile of temperature preference
        temp.p95.kde = percentiles_temp_kde[2],  # 95th percentile
        temp.p99.kde = percentiles_temp_kde[3],  # 99th percentile
        temp.max = max(dist.temp),               # Max observed temperature value (raw, untransformed)
        
        # KDE-derived precipitation percentiles (from cumulative density function)
        prec.p01.kde = percentiles_prec_kde[1],  # 1st percentile of precipitation preference
        prec.p05.kde = percentiles_prec_kde[2],  # 5th percentile
        prec.p10.kde = percentiles_prec_kde[3],  # 10th percentile
        prec.min = min(dist.prec)                # Minimum observed precipitation value
      ),
      
      functions = list(
        species = plant,
        temp_fun = density_temp_fun,
        temp_fun_cdf = density_temp_fun_cdf,
        prec_fun = density_prec_fun,
        prec_fun_cdf = density_prec_fun_cdf
      )
    )
  )
}

stopCluster(cl)

# ------------------------------------------------------------
# Collect outputs into usable objects after parallel KDE loop
# ------------------------------------------------------------

# Combine metadata rows from all species into one data frame
plant.temp <- do.call(rbind, lapply(plant.Temp, `[[`, "data"))

# Fill missing precipitation 1st percentiles with observed minimum
plant.temp$prec.p01.kde[is.na(plant.temp$prec.p01.kde)] <- plant.temp$prec.min[is.na(plant.temp$prec.p01.kde)]

# Extract the kernel functions into a named list
plant.kde.funs <- lapply(plant.Temp, `[[`, "functions")
names(plant.kde.funs) <- sapply(plant.kde.funs, `[[`, "species")


# ------------------------------------------------------------
# Save metadata tables to disk
# ------------------------------------------------------------

# Full dataset including all species
write.table(
  plant.temp,
  file = "intermediate.data/plant_metadata_complete.csv",
  sep = ",",
  row.names = FALSE
)

# Filter to only well-sampled species (occurring in >5 samples)
plant.temp.qualified <- plant.temp[plant.temp$sites2 > 5, ]
write.table(
  plant.temp.qualified,
  file = "intermediate.data/plant_metadata_qualified.csv",
  sep = ",",
  row.names = FALSE
)


# ------------------------------------------------------------
# Diagnostic plot: most abundant species (top 25)
# ------------------------------------------------------------

top_abundant <- plant.temp.qualified %>%
  arrange(desc(abundance)) %>%
  slice(1:25)

top_abundant$crop <- ifelse(top_abundant$crop, "crop", "wild")

# Clean long hybrid names for readability
top_abundant$species[top_abundant$species == "Crataegus monogyna x Crataegus punctata"] <- "Crataegus monogyna"

# Barplot of most abundant species
pdf("plots.supplement/most_abundant_plants.pdf", width = 8, height = 6)
ggplot(top_abundant, aes(x = reorder(species, -abundance), y = abundance, fill = crop)) +
  geom_bar(stat = "identity") +
  labs(x = "Species", y = "Overall relative abundance") +
  scale_fill_viridis(discrete = TRUE, option = "viridis", begin = 0.3, end = 1) +
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )
dev.off()


# ------------------------------------------------------------
# KDE - temperature response curves of top abundant species
# ------------------------------------------------------------

top_species <- top_abundant$species
x_vals <- seq(0, 50, length.out = 200)

kde_data <- do.call(rbind, lapply(plant.Temp, function(entry) {
  sp_name <- entry$functions$species
  if (sp_name %in% top_species) {
    y_vals <- entry$functions$temp_fun(x_vals)
    data.frame(x = x_vals, y = y_vals, species = sp_name)
  } else {
    NULL
  }
}))

kde_modes <- kde_data %>%
  group_by(species) %>%
  filter(y == max(y, na.rm = TRUE)) %>%
  slice(1) %>%  # Handle ties
  ungroup() %>%
  rename(label_x = x, label_y = y)

max_y_label <- max(kde_modes$label_y + 0.01, na.rm = TRUE)

pdf("plots.supplement/most_abundant_plants_temperature_kde.pdf", width = 12, height = 6)
ggplot(kde_data, aes(x = x, y = y, color = species)) +
  geom_line(linewidth = 1) +
  labs(x = "Temperature", y = "Density") +
  theme_minimal() +
  scale_color_viridis(option = "turbo", discrete = TRUE) +
  xlim(5, 30) +
  ylim(0, max_y_label) +
  coord_cartesian(clip = "off") +
  geom_text_repel(
    data = kde_modes,
    aes(x = label_x, y = max_y_label - 0.01, label = species, color = species),
    size = 3.5,
    angle = 90,
    hjust = 0,
    nudge_y = 0.002,
    max.overlaps = Inf,
    box.padding = 0.25,
    point.padding = 0.1,
    direction = "x",
    show.legend = FALSE
  ) +
  theme(legend.position = "none")
dev.off()

# ------------------------------------------------------------
# KDE - precipitation response curves of top abundant species
# ------------------------------------------------------------

x_vals <- seq(0, 300, length.out = 600)

kde_data <- do.call(rbind, lapply(plant.Temp, function(entry) {
  sp_name <- entry$functions$species
  if (sp_name %in% top_species) {
    y_vals <- entry$functions$prec_fun(x_vals)
    data.frame(x = x_vals, y = y_vals, species = sp_name)
  } else {
    NULL
  }
}))

kde_modes <- kde_data %>%
  group_by(species) %>%
  filter(y == max(y, na.rm = TRUE)) %>%
  slice(1) %>%
  ungroup() %>%
  rename(label_x = x, label_y = y)

pdf("plots.supplement/most_abundant_plants_precipitation_kde.pdf", width = 12, height = 6)
ggplot(kde_data, aes(x = x, y = y, color = species)) +
  geom_line(linewidth = 1) +
  labs(x = "Precipitation", y = "Density") +
  theme_minimal() +
  scale_color_viridis(option = "turbo", discrete = TRUE) +
  xlim(0, 250) +
  ylim(0, max_y_label) +
  coord_cartesian(clip = "off") +
  geom_text_repel(
    data = kde_modes,
    aes(x = label_x, y = max_y_label - 0.01, label = species, color = species),
    size = 3.5,
    angle = 90,
    hjust = 0,
    nudge_y = 0.002,
    max.overlaps = Inf,
    box.padding = 0.25,
    point.padding = 0.1,
    direction = "x",
    show.legend = FALSE
  ) +
  theme(legend.position = "none")
dev.off()

#restore original warning settings
options(warn = oldw)
