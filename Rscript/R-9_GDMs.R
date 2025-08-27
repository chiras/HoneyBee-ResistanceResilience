################################################
# Generalized Dissimilarity Modeling (GDM) and NMDS ordination
# Outputs: GDM model summary, NMDS plots (Temp & Prec)
################################################

# ----- GDM MODELING -----

# Aggregate samples by site identifier (BK) to get mean relative abundances
data.species.rel.filter2 <- merge_samples(data.species.rel.filter, "BK", mean)
testData1 <- data.frame(otu_table(data.species.rel.filter2))
testData1 <- testData1[rownames(testData1) != "0", ]
testData1$site <- rownames(testData1)

# Extract environmental metadata per site
meta.data <- data.frame(
  site  = rownames(sample_data(data.species.rel.filter2)),
  lat   = sample_data(data.species.rel.filter2)$CoordY,
  long  = sample_data(data.species.rel.filter2)$CoordX,
  temp  = sample_data(data.species.rel.filter2)$Temp,
  prec  = sample_data(data.species.rel.filter2)$Prec
)
meta.data <- meta.data[meta.data$site != "0", ]

# Format data for GDM
exFormat <- formatsitepair(
  bioData = testData1,
  bioFormat = 1,
  siteColumn = "site",
  XColumn = "long",
  YColumn = "lat",
  predData = meta.data
)

# Run GDM model with geographic distance
model <- gdm(exFormat, geo = TRUE)
summary(model)

# Save GDM plot
pdf("plots.supplement/GDMs.pdf", width = 8, height = 12)
plot(model, plot.layout = c(3, 2))
dev.off()

# ----- EXPLAINED VARIANCE (and scaling of raw deviance explained) -----

# Raw total deviance explained
total <- 48.162
t_r <- 1.815  # temperature
t_g <- 0.932  # geography
t_p <- 0.448  # precipitation

# Convert partial deviance to relative contributions
scale_factor <- total / (t_r + t_g + t_p)
t_r <- t_r * scale_factor
t_g <- t_g * scale_factor
t_p <- t_p * scale_factor

# check total 
t_p + t_g + t_r


# ----- NMDS Ordination -----

# Ordination based on Bray-Curtis dissimilarity
ordination_obj <- ordinate(data.species.rel.filter, method = "NMDS", distance = "bray", k = 5, trymax = 100)
cat("NMDS stress:", ordination_obj$stress, "\n")
test_results[["NMDS_stress"]] <- ordination_obj$stress

# Ordination plots for TEMP

data <- data.species.rel.filter
sample_data(data)$MonthShort <- format(as.Date(sample_data(data)$Month, format = "%d/%m/%Y"), "%B")
sample_data(data)$MonthShort <- factor(sample_data(data)$MonthShort, levels = c("May", "June", "July", "August"), ordered = TRUE)

make_nmds_plot <- function(month_filter = NULL, var = "Temp", color_label = "ºC", title = "All months") {
  show_stress=F
  if (is.null(month_filter)){
    plot.show.legend = T
  }else{
    plot.show.legend = F
  }
  
  if (var == "Prec"){
    plot.col.rev = -1
  } else {
    plot.col.rev = 1
  }
  
  p <- plot_ordination(data, ordination_obj) +
    geom_point(aes(color = !!rlang::sym(var)), size = 2.5) +
    theme_bw() +
    scale_color_viridis_c(option = "inferno", direction = plot.col.rev, limits = global_range) +  # continuous scale for Temp or Prec
    labs(colour = color_label, title = title) + 
    xlim(c(-2.5,2.5))+
    ylim(c(-1.1,1.5))+
    theme(
      legend.position = if (plot.show.legend) "right" else "none",
      legend.text = element_text(size = 16),
      legend.title = element_text(size = 18),
      axis.text = element_text(size = 16),
      axis.title = element_text(size = 16),
      plot.title = element_text(size = 16, hjust = 0.5),
      strip.text = element_text(size = 16),
      strip.background = element_rect(color = "black")
    )
  
  # Apply ellipse with second (discrete) color scale
  if (!is.null(month_filter)) {
    p <- p + gghighlight(MonthShort == month_filter)
    p <- p + ggnewscale::new_scale_color()
    p <- p + suppressWarnings(
      stat_ellipse(data = function(df) dplyr::filter(df, !is.na(MonthShort)), aes(color = MonthShort, group = MonthShort), size = 1.2, show.legend = FALSE)
    )
    p <- p + scale_color_manual(values = c("May" = "#E64B35FF", "June" = "#4DBBD5FF", "July" = "#00A087FF", "August" = "#3C5488FF"))
  } else {
    p <- p + gghighlight()
    p <- p + ggnewscale::new_scale_color()
    p <- p + suppressWarnings(
      stat_ellipse(data = function(df) dplyr::filter(df, !is.na(MonthShort)), aes(color = MonthShort, group = MonthShort), size = 1.2, show.legend = T)
    )
    p <- p + scale_color_manual(values = c("May" = "#E64B35FF", "June" = "#4DBBD5FF", "July" = "#00A087FF", "August" = "#3C5488FF"))
    show_stress = T
  }
  
  if (show_stress && !is.null(test_results[["NMDS_stress"]])) {
    stress_val <- round(test_results[["NMDS_stress"]], 3)
    p <- p + annotate("text", label = paste("Stress =", stress_val), x = -Inf, y = Inf, hjust = -0.1, vjust = 2.5, size = 4, col="darkgray")
  }
  
  return(p)
}

# Combine TEMP plots
global_range <- range(sample_data(data.species.rel.filter)$Temp, na.rm = TRUE)

p1 <- make_nmds_plot(var = "Temp", title = "All months")
p2 <- make_nmds_plot("May", var = "Temp", title = "May")
p3 <- make_nmds_plot("June", var = "Temp", title = "June")
p4 <- make_nmds_plot("July", var = "Temp", title = "July")
p5 <- make_nmds_plot("August", var = "Temp", title = "August")


# Compose patchwork layout
final_plot <- (
  (p1 + patchwork::plot_spacer()) / (p2 + p3) / (p4 + p5)
) +
  plot_layout(
    widths = c(1, 1),
    heights = c(1, 1, 1),
    guides = "collect"
  ) &
  theme(
    legend.position = "right",  # colorbar for Temp (continuous)
    legend.box = "vertical",
    legend.justification = "center"
  )

# Export as square grid
pdf("plots.supplement/NMDS_Temperature.pdf", width = 10, height = 12)
print(final_plot)
dev.off()


# Combine PRECIPITATION plots
global_range <- range(sample_data(data.species.rel.filter)$Prec, na.rm = TRUE)

p6 <- make_nmds_plot(var = "Prec", color_label = "mm", title = "All months")
p7 <- make_nmds_plot("May", var = "Prec", color_label = "mm", title = "May")
p8 <- make_nmds_plot("June", var = "Prec", color_label = "mm", title = "June")
p9 <- make_nmds_plot("July", var = "Prec", color_label = "mm", title = "July")
p10 <- make_nmds_plot("August", var = "Prec", color_label = "mm", title = "August")


# Compose patchwork layout
final_plot <- (
  (p6 + patchwork::plot_spacer()) / (p7 + p8) / (p9 + p10)
) +
  plot_layout(
    widths = c(1, 1),
    heights = c(1, 1, 1),
    guides = "collect"
  ) &
  theme(
    legend.position = "right",  # colorbar for Temp (continuous)
    legend.box = "vertical",
    legend.justification = "center"
  )

# Export as square grid
pdf("plots.supplement/NMDS_Precipitation.pdf", width = 10, height = 12)
print(final_plot)
dev.off()

