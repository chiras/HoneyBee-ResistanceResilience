oldw <- getOption("warn")
options(warn = -1) # suppress warnings temporarily

################################################
# Diversity and Latitude Mapping
################################################
# Dependencies:
#   - Requires:
#       - samples (from R-1)
################################################

print("Create diversity distribution plot")

# Calculate Shannon diversity and exponential (Effective Species Richness)
shannon <- diversity(t(otu_table(samples)))
exp_shannon <- exp(shannon)
richness <- rowSums(t(otu_table(samples))>0)

diversity_df <- data.frame(ExpShannon = exp_shannon,
                        Shannon = shannon,
                        Richness = richness,
                        SampleID = colnames(otu_table(samples)))

sample_data(samples)$SampleID <- rownames(sample_data(samples))

# Merge diversity metrics with sample metadata
diversity_df <- merge(diversity_df, data.frame(sample_data(samples)), by = "SampleID")

# Extract and format time variables
diversity_df$date2 <- as.Date(diversity_df$Month, format = "%d/%m/%Y")
diversity_df$month <- format(diversity_df$date2, "%m")
diversity_df <- diversity_df[!is.na(diversity_df$month), ]

# Month names for faceting
month_names <- c("05" = "May", "06" = "June", "07" = "July", "08" = "August")
diversity_df$month <- factor(month_names[diversity_df$month],
                          levels = c("May", "June", "July", "August"),
                          ordered = TRUE)

# Plot 1: Spatial diversity map by month
world <- map_data("world", regions = c("Czech Republic", unique(diversity_df$Country)))

pdf("plots/diversity_map.pdf", width = 12, height = 4)
ggplot(diversity_df, aes(x = CoordX, y = CoordY, col = ExpShannon)) +
  geom_polygon(data = world, aes(x = long, y = lat, group = group),
               fill = "lightgray", colour = "white") +
  geom_jitter(width = 0.2, height = 0.2, size = 1.5, alpha = 0.5) +
  facet_wrap(~ month, ncol = 4) +
  scale_color_viridis(option = "magma", direction = -1) +
  theme_bw() +
  labs(x = "Longitude", y = "Latitude",
       colour = bquote("Effective Species Richness (" ~ e^H ~ ")")) +
  theme(strip.background = element_rect(fill = "white"),
        panel.grid.minor = element_blank(),
        legend.position = "bottom")
dev.off()

# Plot 2: Latitudinal diversity gradients
pdf("plots/diversity_latitude_quadratic.pdf", width = 12, height = 4)
ggplot(diversity_df, aes(y = ExpShannon, x = CoordY, fill = CoordY)) +
  geom_boxplot(aes(group = as.factor(CoordY)), linewidth = 0.2, outlier.size = 0.1) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), col = "black") +
  facet_wrap(~ month, ncol = 4) +
  scale_fill_viridis(option = "viridis", direction = -1) +
  scale_y_continuous(limits = c(0, 15)) +
  theme_bw() +
  labs(x = "Latitude",
       y = bquote("Effective Species Richness (" ~ e^H ~ ")"),
       fill = "Latitude") +
  theme(strip.background = element_rect(fill = "white"),
        panel.grid.minor = element_blank(),
        legend.position = "bottom")
dev.off()

# Plot 2.b: Latitudinal diversity gradients, normal richness
pdf("plots.supplement/richness_latitude_quadratic.pdf", width = 12, height = 4)
ggplot(diversity_df, aes(y = Richness, x = CoordY, fill = CoordY)) +
  geom_boxplot(aes(group = as.factor(CoordY)), linewidth = 0.2, outlier.size = 0.1) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), col = "black") +
  facet_wrap(~ month, ncol = 4) +
  scale_fill_viridis(option = "viridis", direction = -1) +
  scale_y_continuous(limits = c(0, 25)) +
  theme_bw() +
  labs(x = "Latitude",
       y = bquote("Species Richness"),
       fill = "Latitude") +
  theme(strip.background = element_rect(fill = "white"),
        panel.grid.minor = element_blank(),
        legend.position = "bottom")
dev.off()

# Fit quadratic diversity-latitude model
final_model <- lm(log(ExpShannon) ~ CoordY + I(CoordY^2) + as.numeric(month), data = diversity_df)

# Save model summaries per month
test_results[["diversity_latitude_quadratic"]] <- summary(final_model)

for (m in c("May", "June", "July", "August")) {
  subset_data <- diversity_df[diversity_df$month == m, ]
  model_name <- paste0("diversity_latitude_quadratic_", m)
  test_results[[model_name]] <- summary(
    lm(log(ExpShannon) ~ CoordY + I(CoordY^2), data = subset_data)
  )
}


################################################
# Crop vs. Wild Proportions
################################################

# Calculate total crop and wild abundance per sample
crop_taxa <- crops$V1
crop_otu <- colSums(otu_table(samples)[taxa_names(samples) %in% crop_taxa, ])
wild_otu <- colSums(otu_table(samples)[!(taxa_names(samples) %in% crop_taxa), ])

cw_otu <- otu_table(t(data.frame(crops = crop_otu, wild = wild_otu)), taxa_are_rows = TRUE)
proportion <- otu_table(cw_otu[1, ] / colSums(cw_otu), taxa_are_rows = TRUE)
proportion <- merge_phyloseq(proportion, sample_data(samples))
proportion2 <- psmelt(proportion)

# Prepare time metadata
proportion2$date2 <- as.Date(proportion2$Month, format = "%d/%m/%Y")
proportion2$month <- format(proportion2$date2, "%m")
proportion2 <- proportion2[!is.na(proportion2$month), ]
proportion2$month <- factor(month_names[proportion2$month],
                            levels = c("May", "June", "July", "August"),
                            ordered = TRUE)

# Plot 3: Crop proportions by latitude
pdf("plots/proportion_crops.pdf", width = 12, height = 4)
ggplot(proportion2, aes(y = Abundance, x = CoordY, fill = CoordY)) +
  geom_boxplot(aes(group = as.factor(CoordY)), linewidth = 0.2, outlier.size = 0.1) +
  facet_wrap(~ month, ncol = 4) +
  scale_fill_viridis(option = "viridis", direction = -1) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_bw() +
  labs(x = "Latitude", y = "Proportion crops (%)", fill = "Latitude") +
  theme(strip.background = element_rect(fill = "white"),
        panel.grid.minor = element_blank(),
        legend.position = "bottom")
dev.off()

options(warn = oldw)
