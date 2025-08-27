###############################################################
# Estimates long-term changes in Temperature or Precipitation
# Statistics part
###############################################################

## Historical temperature change — regional summaries
# Define region boundaries
regions <- list(
  Southern = expression(CoordY.x < 45),
  Central  = expression(CoordY.x >= 45 & CoordY.x < 54),
  Northern = expression(CoordY.x >= 54)
)


# Helper function for stats
get_stats <- function(df, var_mean, var_max, min_var = NULL) {
  lapply(regions, function(cond) {
    idx <- eval(cond, envir = df)
    list(
      mean_change = mean(df[[var_mean]][idx], na.rm = TRUE),
      max_change  = max(df[[var_max]][idx], na.rm = TRUE),
      range_diff  = if (!is.null(min_var)) diff(range(df[idx, c("min", "max")], na.rm = TRUE)) else NA,
      min_change  = if (!is.null(min_var)) min(df[[min_var]][idx], na.rm = TRUE) else NA
    )
  })
}

# Temperature stats
temp_stats <- get_stats(
  df = test_results$Temperature_historical_change_since_1973_within,
  var_mean = "mean",
  min_var  = "min",
  var_max  = "max"
)

# Precipitation stats
prec_stats <- get_stats(
  df = test_results$Precipitation_historical_change_since_1973_within,
  var_mean = "mean",  
  min_var  = "min",
  var_max  = "max"
)

# Combine into data frame
overview <- data.frame(
  Region        = names(regions),
  Temp_mean     = sapply(temp_stats, `[[`, "mean_change"),
  Temp_max      = sapply(temp_stats, `[[`, "max_change"),
  Temp_range    = sapply(temp_stats, `[[`, "range_diff"),
  Prec_mean    = sapply(prec_stats, `[[`, "mean_change"),
  Prec_min      = sapply(prec_stats, `[[`, "min_change"),
  Prec_range    = sapply(prec_stats, `[[`, "range_diff")
)

print(overview)
test_results[["historical_change_since_1973_latitude_overview"]] <- overview
