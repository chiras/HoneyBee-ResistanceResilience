
################################################
# Data analysis 2
## Estimate plant temperature and precipitation distributions

print("Estimate plant temperature and precipitation distributions")
pdf("plots.supplement/temperature_distribution_samples.pdf", width=7, height=6.5)
  hist(data.map$Temp)
dev.off()

pdf("plots.supplement/precipitation_distribution_samples.pdf", width=7, height=6.5)
  hist(data.map$Prec)
dev.off()

### use multiple cores to speed up the process
cores=detectCores()
cl <- makeCluster(cores[1]-2) #not to overload your computer
registerDoSNOW(cl)

### prepare data for parallelization (phyloseq objects do not work to be passed)
tmp.otu <-otu_table(samples)
tmp.smp <-sample_data(samples)
tmp.tax <-tax_table(samples)
plant_taxa<- taxa_names(samples)

### create progress bar
pb <- txtProgressBar(0, length(plant_taxa), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)

### parallelisation over all plants in the dataset
plant.temp <- data.frame()
plant.temp2 <-foreach (plant_i = 1:length(plant_taxa), .combine=rbind, .options.snow=opts, .packages=c("phyloseq","dplyr")) %dopar% {
 
  # subset to target plant and prepare data
  plant <- plant_taxa[plant_i]
  target.otu <- tmp.otu[plant,]
  target <- merge_phyloseq(target.otu, tmp.smp, tmp.tax)
  target.melt <- psmelt(target)
  abundance.plant <- (sum(target.otu)/length(target.otu))
  
  if (!(sum(target.melt$Abundance) == 0) ){
  
  target.melt <- target.melt[target.melt$Abundance > 0,]
  target.melt$location = interaction(target.melt$CoordY, target.melt$CoordX)
  
  # Calculate endemism
  countries <- length(data.frame(table(target.melt$Country))[,1])
  sites <- length(data.frame(table(target.melt$location))[,1])
  sites2 <- length(unique(target.melt$location))

  # Check whether it is a crop
  crop <- plant %in% crops$V1 
  
  # Temperature & Precipitation
    records <- sum(target.melt$Abundance!=0) 
    
    # round to 0.5s
    target.melt$Temp2 <- round(target.melt$Temp/0.5)*0.5

    # square root for Precipitation due to scewed distribution
    target.melt$Prec2 <- sqrt(target.melt$Prec)
    
    # accumulate over all sites
    target.temp <- target.melt %>%
      group_by(Temp2) %>%
      summarize(accumulated = sum(Abundance))

    target.prec <- target.melt %>%
      group_by(Prec2) %>%
      summarize(accumulated = sum(Abundance))  
    
    # prepare for distribution calculation 
    dist.temp <- c()
    dist.prec <- c()
    
    for (temp in target.temp$Temp2){
      dist.tmp <- rep(temp,round(target.temp$accumulated[target.temp$Temp2==temp], digits=3)*1000)
      dist.temp=c(dist.temp,dist.tmp)
    }

    for (prec in target.prec$Prec2){
      dist.tmp <- rep(prec,round(target.prec$accumulated[target.prec$Prec2==prec], digits=3)*1000)
      dist.prec=c(dist.prec,dist.tmp)
    }

    # calculate mean and standard deviation    
    mean.temp = mean(dist.temp)
    sd.temp = sd(dist.temp)

    mean.prec = mean(dist.prec)
    sd.prec = sd(dist.prec)
    
    # set minimum SD to 1
    if (sd.prec < 1){
      sd.prec =1
    }

    if (sd.temp < 1){
      sd.temp =1
    }

    # calculate distribution metrics    
    quantiles.temp <- qnorm(c(0.05,0.25,0.75,0.90,0.95,0.99), mean = mean.temp, sd = sd.temp, lower.tail = TRUE, log.p = FALSE)
    quantiles.prec <- qnorm(c(0.01,0.05,0.10,0.25,0.50,0.75), mean = mean.prec, sd = sd.prec, lower.tail = TRUE, log.p = FALSE)
    
    # create output line for the plant
    tmp <- data.frame(species = plant,
                      n_records = records,
                      crop = crop,
                      countries = countries,
                      sites = sites,
                      sites2 = sites2,
                      abundance = abundance.plant,
                      temp.mean = mean.temp, 
                      temp.sd = sd.temp,
                      temp.median= median(dist.temp),
                      temp.q05 = quantiles.temp[1],
                      temp.q25 = quantiles.temp[2],
                      temp.q75 = quantiles.temp[3],
                      temp.q90 = quantiles.temp[4],
                      temp.q95 = quantiles.temp[5],
                      temp.q99 = quantiles.temp[6],
                      temp.max = max(dist.temp),
                      prec.mean = mean.prec, 
                      prec.sd = sd.prec,
                      prec.median= median(dist.prec),
                      prec.q01 = quantiles.prec[1],
                      prec.q05 = quantiles.prec[2],
                      prec.q10 = quantiles.prec[3],
                      prec.q25 = quantiles.prec[4],
                      prec.q50 = quantiles.prec[5],
                      prec.q75 = quantiles.prec[6],
                      prec.min = min(dist.prec))
 
    plant.temp <- rbind(plant.temp,tmp) 
    tmp
  }
}
stopCluster(cl)

plant.temp <- plant.temp2 # when using parallelization
write.table(plant.temp, file="tmp.plant.csv", sep=",")


plant.temp.qualified <- plant.temp[plant.temp$sites2 > 5,]
plant.temp.qualified$prec.mean=plant.temp.qualified$prec.mean^2
plant.temp.qualified$prec.sd=plant.temp.qualified$prec.sd^2

write.table(plant.temp.qualified, file="plant_niche_distributions.csv", sep=",")

# plot most abundant plants and their climatic distributions
top_abundant <- plant.temp.qualified %>%
  arrange(desc(abundance)) %>%
  slice(1:25)

top_abundant$crop[top_abundant$crop] <- "crop"
top_abundant$crop[top_abundant$crop==FALSE] <- "wild"

top_abundant$species[top_abundant$species=="Crataegus monogyna x Crataegus punctata"] <- "Crataegus monogyna"

# Create the barplot
pdf("plots.supplement/most_abundant_plants.pdf", width=8, height=6)

ggplot(top_abundant, aes(x = reorder(species, -abundance), y = abundance, fill = crop)) +
  geom_bar(stat = "identity") +
  labs(x = "Species", 
       y = "Overall relative abundance") +
  scale_fill_viridis(discrete=T, option="viridis",begin = 0.3,end = 1)+
  theme_minimal() +
  theme(legend.title = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

dev.off()

# temperature distributions
generate_normal_data <- function(mean, sd, species) {
  x_vals <- seq(0, 40, length.out = 100)  # Fixed range from 0 to 40
  y_vals <- dnorm(x_vals, mean = mean, sd = 1.2*sd)  # Generate corresponding y values using normal distribution
  data.frame(x = x_vals, y = y_vals, species = species)  # Return a data frame with the species name
}

normal_data <- do.call(rbind, lapply(1:nrow(top_abundant), function(i) {
  generate_normal_data(top_abundant$temp.mean[i], top_abundant$temp.sd[i], top_abundant$species[i])
}))

# Create the plot
pdf("plots.supplement/most_abundant_plants_temperature.pdf", width=12, height=6)

ggplot(normal_data, aes(x = x, y = y, color = species)) +
  geom_line(linewidth = 1) +
  labs(x = "Temperature", 
       y = "Density") +
  theme_minimal() +
  scale_color_viridis(option="turbo", discrete=T)+
  xlim(0, 35) +  
  ylim(0, 0.4) +  
  theme(legend.title = element_blank(),
        legend.position = "bottom")  +
  guides(color = guide_legend(ncol = 3))+
  geom_text_repel(data = top_abundant, 
                   aes(x = temp.mean, y = 0.2, 
                       label = species, color = species),
                   size = 4,        
                   angle = 90,hjust = 0,
                   nudge_y = 0.03,  
                   box.padding = 0.55,  
                   point.padding = 0.5, 
                   direction = "x",  
                   show.legend = FALSE) + 
  theme(legend.position = "none") 

dev.off()

# precipitation distributions
generate_normal_data <- function(mean, sd, species) {
  x_vals <- seq(0, 150, length.out = 100)  
  y_vals <- dnorm(x_vals, mean = mean, sd = 1.2*sd) 
  data.frame(x = x_vals, y = y_vals, species = species)  
}

normal_data <- do.call(rbind, lapply(1:nrow(top_abundant), function(i) {
  generate_normal_data(top_abundant$prec.mean[i], top_abundant$prec.sd[i], top_abundant$species[i])
}))

# Create the plot
pdf("plots.supplement/most_abundant_plants_precipitation.pdf", width=8, height=6)
ggplot(normal_data, aes(x = x, y = y, color = species)) +
  geom_line(linewidth = 1) +
  labs(x = "Precipitation", 
       y = "Density") +
  theme_minimal() +
  scale_color_viridis(option="turbo", discrete=T)+
  xlim(0, 150) +  
  ylim(0, 0.3) +  
  theme(legend.title = element_blank(),
        legend.position = "bottom")  +
  guides(color = guide_legend(ncol = 3))+
  geom_text_repel(data = top_abundant, 
                   aes(x = prec.mean, y = 0.12, 
                       label = species, color = species),
                   size = 4,        
                   angle = 90, hjust = 0,
                   nudge_y = 0.03,  
                   box.padding = 0.55,  
                   point.padding = 0.5, 
                   direction = "x",  
                   show.legend = FALSE) + 
  theme(legend.position = "none")  

dev.off()
