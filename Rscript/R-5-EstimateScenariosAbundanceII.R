
# ----------------------------------------------------------------
# Temporarily suppress warnings (e.g. missing directories or NA coercion)
# ----------------------------------------------------------------
oldw <- getOption("warn")
options(warn = -1)

# ----------------------------------------------------------------
# Load predicted abundance loss data for the selected climate variable
# ----------------------------------------------------------------
if (measure == "Precipitation") {
  print("Loading prediction.precipitation from file...")
  assign("prediction",
         read.csv("intermediate.data/prediction.Precipitation.csv", stringsAsFactors = FALSE),
         envir = .GlobalEnv)
}

if (measure == "Temperature") {
  print("Loading prediction.temperature from file...")
  assign("prediction",
         read.csv("intermediate.data/prediction.Temperature.csv", stringsAsFactors = FALSE),
         envir = .GlobalEnv)
}

print(paste("Estimate future scenarios:", measure))

# ----------------------------------------------------------------
# Clean: remove predictions with NA risk classifications
# ----------------------------------------------------------------
prediction <- prediction[!is.na(prediction$prediction), ]

# ----------------------------------------------------------------
# Compute total pollen abundance (baseline delta = 0) per latitude
# ----------------------------------------------------------------
predictionSum <- prediction[prediction$delta == 0, ] %>%
  group_by(coordY) %>%
  summarize(countrysum = sum(abundance, na.rm = TRUE))

# ----------------------------------------------------------------
# Load ISO country codes and merge with prediction data
# ----------------------------------------------------------------
countrycodes <- read.table("Data/countrycodes.csv", sep = ",", header = TRUE)
rownames(countrycodes) <- countrycodes$name
prediction$countrycodes <- countrycodes[prediction$country, "alpha.2"]

# ----------------------------------------------------------------
# Summarize countries per latitude zone (for labels/facets)
# ----------------------------------------------------------------
predictionSum2 <- prediction[prediction$delta == 0, ] %>%
  group_by(coordY) %>%
  summarize(countriesX = paste(unique(country), collapse = ", "))

predictionSum3 <- prediction[prediction$delta == 0, ] %>%
  group_by(coordY) %>%
  summarize(countrycodesX = paste(unique(countrycodes), collapse = ","))

# ----------------------------------------------------------------
# Aggregate pollen loss per latitude × scenario × risk class × crop status
# ----------------------------------------------------------------
prediction2 <- prediction %>%
  group_by(coordY, prediction, crop, delta) %>%
  summarize(accumulated = sum(abundance, na.rm = TRUE), .groups = "drop")

# ----------------------------------------------------------------
# Merge auxiliary metadata (total abundance, country names, country codes)
# ----------------------------------------------------------------
prediction2 <- prediction2 %>%
  merge(predictionSum, by = "coordY") %>%
  merge(predictionSum2, by = "coordY") %>%
  merge(predictionSum3, by = "coordY")

# ----------------------------------------------------------------
# Compute relative abundance loss within each latitude zone
# ----------------------------------------------------------------
prediction2$accumulated.rel <- prediction2$accumulated / prediction2$countrysum

# ----------------------------------------------------------------
# Format categorical variables for plotting consistency
# ----------------------------------------------------------------

# Standardize crop labels (TRUE → "crop", FALSE → "wild")
prediction2$crop <- ifelse(prediction2$crop, "crop", "wild")

# Ensure factors are properly ordered
prediction2$prediction <- factor(prediction2$prediction)
prediction2$crop <- factor(prediction2$crop)
prediction2$delta <- factor(prediction2$delta)

# Latitude factor: decreasing to plot north (top) to south (bottom)
prediction2$coordY <- factor(prediction2$coordY,
                             levels = sort(unique(as.numeric(prediction2$coordY)), decreasing = TRUE),
                             ordered = TRUE)

# Create interaction variables for plot facets and stacked bar fills
prediction2$iPC <- interaction(prediction2$prediction, prediction2$crop)
prediction2$iPC <- factor(prediction2$iPC,
                          levels = sort(unique(as.character(prediction2$iPC))),
                          ordered = TRUE)

prediction2$iCC <- interaction(prediction2$coordY, prediction2$countriesX, sep = ": ")
prediction2$iCC <- factor(prediction2$iCC,
                          levels = sort(unique(as.character(prediction2$iCC)), decreasing = TRUE),
                          ordered = TRUE)

prediction2$iCC2 <- interaction(prediction2$coordY, prediction2$countrycodesX, sep = ": ")
prediction2$iCC2 <- factor(prediction2$iCC2,
                           levels = sort(unique(as.character(prediction2$iCC2)), decreasing = TRUE),
                           ordered = TRUE)

# Explicit reordering of prediction risk classes (percentile thresholds)
# This ensures plot colors are consistent and logical
risk_levels <- c(">p99", ">p95", ">p90", "<=p90", ">=p10", "<p10", "<p05", "<p01")
prediction2$prediction <- factor(prediction2$prediction,
                                 levels = risk_levels[risk_levels %in% levels(prediction2$prediction)],
                                 ordered = TRUE)

#### prepare risk plot, adjust data with respect to temperature or precipitation
if(measure == "Precipitation"){
  reverse = F
  prediction2$delta <-  factor(round(as.numeric(as.character(prediction2$delta)), digits=1))
  prediction$delta <-  factor(round(as.numeric(as.character(prediction$delta)), digits=1))
  
  prediction2$prediction <- factor(prediction2$prediction, levels=sort(as.character(unique( prediction2$prediction)), decreasing =F), ordered=T)
  prediction2$iPC <- factor(prediction2$iPC, levels=sort(as.character(unique(prediction2$iPC)), decreasing =T), ordered=TRUE)
  criterion=">=p10"
  prediction2$delta2 <- paste("-",prediction2$delta,"mm",sep="")
  prediction2$delta2 <- factor(prediction2$delta2,levels=unique(prediction2$delta2[order(prediction2$delta)]), ordered=T)
  subset_levels <- seq(0,50,10)

}else{
  reverse = T
  criterion="<=p90" 
  prediction2$delta2 <- paste("+",prediction2$delta,"°C",sep="")
  prediction2$delta2 <- factor(prediction2$delta2,levels=unique(prediction2$delta2[order(prediction2$delta)]), ordered=T)
  subset_levels <- seq(0,5,1)
}

cols <- c(
  "#1f78b4", "#1f78b4",   # blue for low risk
  "#f7e379", "#f7e379",   # yellow for moderate
  "#f2a134", "#f2a134",   # orange for high
  "#e51f1f", "#e51f1f"    # red for critical
)

#### full plot, all scenarios (supplement)

# revert the crops to negative values for mirror plotting
prediction2.r <- prediction2 %>%
  mutate(accumulated.rel = ifelse(crop == "crop", -accumulated.rel, accumulated.rel))

pdf(paste0("plots.supplement/risk_losses_abundance_latitude_",measure,"_full.pdf"), width=15, height=8)

ggplot(prediction2.r, aes(fill=iPC, shape=crop,y=accumulated.rel, x=1)) + 
  geom_bar(position = position_stack(reverse = T),stat="identity",) +
  facet_grid(iCC ~delta2, scales="free", drop=T)+ 
  geom_hline(yintercept = 0, color = "white", size = 0.5) +
  coord_flip()+
  scale_fill_manual(values =c(cols))+
  scale_color_viridis(option="viridis", discrete=T) + 
  theme_bw()+
  xlim(c(-0.5,1))+
  scale_alpha_discrete(range = c(0.7, 1))+
  theme_bw()+ 
  theme(panel.margin.y = unit(0, "lines"),
        strip.text.y = element_text(angle = 360, hjust = 0),
        strip.background =element_rect(fill="white"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  xlab("Latitude") + ylab("Proportion of pollen")+
  scale_y_continuous(breaks = c(0,0.5))+
  scale_x_continuous(breaks = c(0,2))

dev.off()


#### subset for main text
prediction2.1.r <- prediction2.r[prediction2.r$delta %in% subset_levels,] 

pdf(paste0("plots/risk_losses_abundance_latitude_",measure,"_sub.pdf"), width=8, height=6.5)

ggplot(prediction2.1.r, aes(fill=iPC, shape=crop,y=accumulated.rel, x=1)) + 
  #geom_bar(position = position_stack(reverse = T),stat="identity",fill="black") +
  geom_bar(position = position_stack(reverse = T),stat="identity",) +
  facet_grid(iCC ~delta2, scales="free", drop=T)+ 
  geom_hline(yintercept = 0, color = "white", size = 0.5) +
  coord_flip()+
  scale_fill_manual(values =c(cols))+
  scale_color_viridis(option="viridis", discrete=T) + theme_bw()+ 
  scale_alpha_discrete(range = c(0.7, 1))+
  theme_bw()+ 
  xlim(c(-0.5,1))+
  theme(panel.margin.y = unit(0, "lines"),
        strip.text.y = element_text(angle = 360),
        strip.background =element_rect(fill="white"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  xlab("Latitude") + ylab("Proportion of pollen")+
  scale_y_continuous(breaks = c(0,0.5))+
  scale_x_continuous(breaks = c(0,2)) + 
  theme(legend.position="bottom")

dev.off()



options(warn = oldw)

# ----------------------------------------------------------------
# Linear model: test whether latitude and scenario magnitude explain risk
# ----------------------------------------------------------------
# This fits a model including main effects and their interaction.
# Note: coordY and delta2 are factors, so we coerce to numeric for linear trend.

prediction.at.risk <- prediction2[prediction2$prediction != criterion,]
prediction.at.risk <- prediction.at.risk %>%
  group_by(coordY,delta2) %>%
  summarize(accumulated.rel = sum(accumulated.rel, na.rm=T))

res.man <- lm(accumulated.rel ~ as.numeric(coordY)+as.numeric(coordY):as.numeric(delta2)+as.numeric(delta2), data = prediction.at.risk)
summary(res.man)

test_results[[paste0("LinearModelSzenarios_",measure)]] <- summary(res.man)