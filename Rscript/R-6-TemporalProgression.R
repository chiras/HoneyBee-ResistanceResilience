##############################################################
# R-6-TemporalProgression.R
# Visualizes seasonal shifts in risk (across months) per latitude
# under different climate change scenarios.
#
# Dependencies:
#   - Requires: `prediction` from R-5
#   - Requires: `criterion` defining low-risk class (e.g. "<=p90" or ">=p10")
# Outputs:
#   - seasonal risk plots (per delta)
#   - seasonal linear model stored in `test_results`
##############################################################

# ------------------------------------------------------------
# Extract temporal information
# ------------------------------------------------------------
# Ensure prediction dates are parsed
prediction$date2 <- as.Date(prediction$date, format = "%d/%m/%Y")
prediction$month <- format(prediction$date2, "%m")
prediction$week <- strftime(prediction$date2, format = "%V")


# ------------------------------------------------------------
# Scenario-specific labeling and axis title
# ------------------------------------------------------------
if (measure == "Precipitation") {
  criterion<-">=p10"
  ylabel <- "% <p10"
  prediction$delta2 <- paste0("-", prediction$delta, "mm")
  prediction$delta2 <- factor(
    prediction$delta2,
    levels = unique(prediction2$delta2[order(as.numeric(prediction2$delta))]),
    ordered = TRUE
  )
} else {
  criterion<-"<=p90"
  ylabel <- "% >p90"
  prediction$delta2 <- paste0("+", prediction$delta, "°C")
  prediction$delta2 <- factor(
    prediction$delta2,
    levels = unique(prediction2$delta2[order(as.numeric(prediction2$delta))]),
    ordered = TRUE
  )
}

minfilter=30

# ------------------------------------------------------------
# Filter to relevant subset: remove NAs and non-risk entries
# ------------------------------------------------------------
prediction.season <- prediction[!is.na(prediction$prediction), ]
prediction.season <- prediction.season[prediction.season$prediction != criterion, ]
prediction.season <- prediction.season[!is.na(prediction.season$month), ]

# ------------------------------------------------------------
# Aggregate by latitude, scenario, month, and site
# ------------------------------------------------------------
prediction.season2 <- prediction.season %>%
  group_by(coordY, delta, delta2, month, site) %>%
  summarize(accumulated = sum(abundance, na.rm = TRUE), .groups = "drop")

# Round latitude into broader zones (every 2 degrees)
prediction.season2$coordY2 <- round(prediction.season2$coordY / 2, digits = 0) * 2

# ------------------------------------------------------------
# Plot 1: Full range of climate scenarios across months
# ------------------------------------------------------------
pdf(paste0("plots.supplement/season_q90_distribution",measure,"full.pdf"), width = 15, height = 6)

ggplot(prediction.season2, aes(
  x = as.numeric(month),
  y = accumulated,
  col = factor(coordY2),
  fill = factor(coordY2)
)) +
  geom_smooth(
    data = prediction.season2 %>% group_by(coordY2, delta) %>% filter(n() > minfilter),
    method = "loess",
    se = TRUE,
    alpha = 0.1
  ) +
  facet_wrap(~delta2, ncol = 6) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_color_viridis(discrete = TRUE, option = "viridis", direction = -1) +
  scale_fill_viridis(discrete = TRUE, option = "viridis", direction = -1, alpha = 0.01) +
  theme_bw() +
  labs(x = "Month", y = ylabel) +
  theme(
    panel.margin.y = unit(0, "lines"),
    strip.background = element_rect(fill = "white"),
    panel.grid.minor = element_blank()
  )

dev.off()

# ------------------------------------------------------------
# Plot 2: Subset of scenarios for main figure
# ------------------------------------------------------------
prediction.season2.1 <- prediction.season2[prediction.season2$delta %in% subset_levels, ]
prediction.season2.1$delta <- factor(prediction.season2.1$delta)

pdf(paste0("plots/season_q90_distribution",measure,"_sub.pdf"), width = 15, height = 3)

ggplot(prediction.season2.1, aes(
  x = as.numeric(month),
  y = accumulated,
  col = factor(coordY2),
  fill = factor(coordY2)
)) +
  geom_smooth(
    data = prediction.season2.1 %>% group_by(coordY2, delta) %>% filter(n() > minfilter),
    method = "loess",
    se = TRUE,
    alpha = 0.1
  ) +
  facet_wrap(~delta2, ncol = 6) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_color_viridis(discrete = TRUE, option = "viridis", direction = -1) +
  scale_fill_viridis(discrete = TRUE, option = "viridis", direction = -1, alpha = 0.01) +
  theme_bw() +
  labs(x = "Month", y = ylabel) +
  theme(
    panel.margin.y = unit(0, "lines"),
    strip.background = element_rect(fill = "white"),
    panel.grid.minor = element_blank()
  )

dev.off()

# ------------------------------------------------------------
# Linear model: predict risk over time and latitude
# ------------------------------------------------------------
prediction.at.risk <- prediction.season[prediction.season$prediction != criterion,] %>%
  group_by(coordY, delta2, month) %>%
  summarize(accumulated.rel = sum(abundance, na.rm = TRUE), .groups = "drop")

res.man <- lm(
  accumulated.rel ~
    as.numeric(month) +
    as.numeric(month):as.numeric(delta2) +
    as.numeric(coordY) +
    as.numeric(delta2):as.numeric(coordY) +
    as.numeric(delta2),
  data = prediction.at.risk
)

summary(res.man)

test_results[[paste0(measure, ".model.risk")]] <- summary(res.man)

# Restore warning settings
options(warn = oldw)
