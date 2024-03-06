
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
  
  if (!(sum(target.melt$Abundance) == 0) ){
  
  target.melt <- target.melt[target.melt$Abundance > 0,]
  target.melt$location = interaction(target.melt$CoordY, target.melt$CoordX)
  
  # Calculate endemism
  countries <- length(data.frame(table(target.melt$Country))[,1])
  sites <- length(data.frame(table(target.melt$location))[,1])

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
