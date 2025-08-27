###############################################################
# Estimates long-term changes in Temperature or Precipitation
# Input: measure = "Temperature" or "Precipitation"
###############################################################

# Load historical data file for selected measure
historical.data <- read.table(paste0("Data/", measure, "_Data.csv"), sep = ",", header = TRUE)

# Merge with sample metadata
meta <- data.frame(sample_data(samples))[, c("BK", "Month", "CoordY")]
historical.data.merge <- merge(meta, historical.data, by = "BK")

# Reshape historical time series data
hist.melt <- reshape2::melt(historical.data.merge, id.vars = c("BK", "Month", "CoordY"), value.name = "Value")
hist.melt$variable <- gsub("^X", "", hist.melt$variable)
hist.melt$year_record <- gsub("\\..*", "", hist.melt$variable)
hist.melt$month_record <- gsub(".*\\.", "", hist.melt$variable)

# Extract sampling month
hist.melt$month_sample <- gsub("^\\d{1,2}/", "", hist.melt$Month)
hist.melt$month_sample <- gsub("/\\d{4}$", "", hist.melt$month_sample)

# Filter to valid month comparisons (historical month = sampling month)
hist.melt2 <- hist.melt[hist.melt$month_record == hist.melt$month_sample, ]

# Aggregate duplicates (same BK, year, month)
hist.melt3 <- hist.melt2 %>%
  group_by(BK, variable, CoordY, year_record, month_record, month_sample) %>%
  summarize(value = mean(as.numeric(Value), na.rm = TRUE), .groups = "drop")

# Extract current (2023) values
hist.melt3.current <- hist.melt3[hist.melt3$year_record == "2023", ]

# Merge historical and current
hist.melt4 <- merge(hist.melt3, hist.melt3.current, by = c("BK", "month_sample"))
hist.melt4$year_record <- as.numeric(hist.melt4$year_record.x)

# === Plot 1: Full time series ===
pdf(paste0("plots.supplement/historical_", measure, "_change.pdf"), width = 5, height = 4)
ggplot(hist.melt4, aes(x = year_record, y = value.x, col = as.factor(CoordY.y))) +
  geom_smooth(method = "loess", n = 100, aes(fill = as.factor(CoordY.y)), alpha = 0.25) +
  geom_smooth(method = "loess", n = 100, se = FALSE) +
  scale_color_viridis(discrete = TRUE, option = "viridis", direction = -1) +
  scale_fill_viridis(discrete = TRUE, option = "viridis", direction = -1, alpha = 0.1) +
  xlim(1900, 2023) +
  theme_bw() +
  labs(x = "Time", y = paste("Monthly average", measure))
dev.off()

# === Difference between 1970s and 2023 ===
hist.melt5 <- hist.melt4[hist.melt4$year_record.x > 1969 & hist.melt4$year_record.x < 1976, ]

# Calculate change between past and present
hist.melt5$difference <- hist.melt5$value.x - hist.melt5$value.y
hist.melt5$difference <- -hist.melt5$difference  # Flip sign so increase = warming/wetting

# Axis config
if (measure == "Temperature") {
  breaks_plot <- seq(-10, 10, 1)
  ymin <- 5
  ymin2 <- 2.5
} else {
  breaks_plot <- seq(-200, 250, 25)
  ymin <- -50
  ymin2 <- 50
}

# === Statistical Tests: Latitude correlation ===
hist.mean.lat <- hist.melt5 %>%
  group_by(CoordY.x, month_record.x) %>%
  summarize(mean = mean(difference, na.rm = TRUE), .groups = "drop")

test_results[[paste0(measure, "_historical_change_since_1973_latitude")]] <-
  cor.test(hist.mean.lat$mean, as.numeric(hist.mean.lat$CoordY.x))

# === Statistical Tests: Within-latitude changes ===
wilcox.tests <- do.call(rbind, lapply(unique(hist.melt5$month_record.x), function(month) {
  do.call(rbind, lapply(unique(hist.melt5$CoordY.x), function(latitude) {
    subset_vals <- dplyr::filter(hist.melt5, CoordY.x == latitude & month_record.x == month)$difference
    if (length(subset_vals) > 1) {
      test <- wilcox.test(subset_vals, mu = 0)
      data.frame(
        CoordY.x = latitude,
        month_record.y = month,
        p.value = round(test$p.value, 3),
        est = test$statistic,
        mean = mean(subset_vals),
        max = max(subset_vals),
        min = min(subset_vals),
        ue = quantile(subset_vals, 1),
        le = quantile(subset_vals, 0),
        ue2 = median(subset_vals) + 1.58 * IQR(subset_vals) / sqrt(length(subset_vals)),
        le2 = median(subset_vals) - 1.58 * IQR(subset_vals) / sqrt(length(subset_vals))
      )
    } else {
      NULL
    }
  }))
}))

# Annotate significance
wilcox.tests$signif2 <- ifelse(substr(as.character(wilcox.tests$mean), 1, 1) != "-", "+", "-")
wilcox.tests$signif <- ""
wilcox.tests$signif[wilcox.tests$p.value < 0.1] <- "'"
wilcox.tests$signif[wilcox.tests$p.value < 0.05] <- "*"
wilcox.tests$signif[wilcox.tests$p.value < 0.01] <- "**"
wilcox.tests$signif[wilcox.tests$p.value < 0.001] <- "***"

test_results[[paste0(measure, "_historical_change_since_1973_within")]] <- wilcox.tests

# === Plot 3: Boxplots with significance annotations ===
hist.annot <- merge(hist.melt5, wilcox.tests, by = c("CoordY.x", "month_record.y"))

# Reformat months
month_map <- c("05" = "May", "06" = "June", "07" = "July", "08" = "August")
hist.annot$month_sample <- factor(month_map[hist.annot$month_sample],
                                  levels = c("May", "June", "July", "August"),
                                  ordered = TRUE)

pdf(paste0("plots.supplement/historical_", measure, "_1973.pdf"), width = 8, height = 6)

ggplot(hist.annot, aes(x = CoordY.x, y = difference, fill = as.factor(CoordY.x))) +
  annotate("rect", xmin = min(t.tests$CoordY.x) - 1, xmax = max(t.tests$CoordY.x) + 1,
           ymin = 0, ymax = ymin, alpha = 0.4, fill = "steelblue") +
  geom_hline(yintercept = 0, linewidth = 1, col = "red") +
  geom_boxplot(linewidth = 0.2, outlier.size = 0.1) +
  geom_smooth(aes(group = month_sample), method = "loess", se = FALSE, color = "black", size = 0.8)+
  facet_wrap(~ month_sample) +
  scale_y_continuous(breaks = breaks_plot) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_color_viridis(discrete = TRUE, option = "viridis", direction = -1) +
  scale_fill_viridis(discrete = TRUE, option = "viridis", direction = -1) +
  geom_text(aes(y = le - ymin2, label = signif), size = 3, lineheight = 0.5) +
  geom_text(aes(y = ue + ymin2, label = signif2), size = 3) +
  theme_bw() +
  labs(x = "Latitude", y = paste(measure, "difference vs. 1970–1975")) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "white"))
dev.off()

