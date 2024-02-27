
################################################
# Data analysis 3
## Estimate future scenarios, per site
print("Estimate future scenarios: Preparation")

sample_names <- sample_names(samples)#[c(1,1999,2111,322,1199)]
predictionTemp <- data.frame()

### setting parallel computing
cl <- makeCluster(cores[1]-2) 
registerDoSNOW(cl)

predictionPrec <- data.frame()
predictionTemp <- data.frame()

prec_levels <- seq(0,25,2.5)
temp_levels <- seq(0,5,0.5)

## TEMPERATURE
pb <- txtProgressBar(0, length(temp_levels), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)

prediction.temperature <- foreach(prec_increase=temp_levels, .options.snow=opts , .combine="rbind",  .packages=c("phyloseq","tidyr","speedyseq","dplyr")) %dopar% {
  # remove and set global variables
  rm("temp_increase","siteY","site.phyloseq","site.melt","plant.p","plant","prediction.tmp","plant.species.meta")
  .GlobalEnv$samples <- samples
  .GlobalEnv$sample_names <- sample_names
  temp_increase = prec_increase
  
  # iterate over sites
  predictionPlants=data.frame()
  site.count=c()
  for (siteY in sample_names){
    site.otu <- tmp.otu[,siteY]
    site.phyloseq <- merge_phyloseq(site.otu,tmp.smp,tmp.tax)
    
    site.melt <- psmelt(site.phyloseq)
    site.melt <- site.melt[site.melt$Abundance > 0,]
    site.count=c(site.count,siteY)

    # iterate over plants
    for (plant.p in site.melt$species){
      plant.species.meta = plant.temp[plant.temp$species == plant.p,]
      
      prediction.tmp <- data.frame(species = plant.p,
                                   site = site.melt$Sample[1],
                                   date = site.melt$Month[1],
                                   country = site.melt$Country[1],
                                   coordY = site.melt$CoordY[1],
                                   abundance = site.melt[site.melt$species == plant.p,"Abundance"],
                                   prediction = NA,
                                   p.dist.value = pnorm(site.melt$Temp[1]+temp_increase, mean=plant.species.meta$temp.mean, sd=plant.species.meta$temp.sd),
                                   endemism = plant.species.meta$countries,
                                   temp_increase = temp_increase,
                                   crop = plant.species.meta$crop)
      
      if (plant.species.meta$sites > 10){
        if (plant.species.meta$temp.q99 < site.melt$Temp[1]+temp_increase){
            prediction.tmp$prediction <- ">q99"
          }else if (plant.species.meta$temp.q95 < site.melt$Temp[1]+temp_increase){
            prediction.tmp$prediction <- ">q95"
          }else if (plant.species.meta$temp.q90 < site.melt$Temp[1]+temp_increase){
            prediction.tmp$prediction <- ">q90"
          }else{
            prediction.tmp$prediction = "<=q90"
          }
        
      }
      # maybe real q value? comes later...
      rm("plant.p","plant")

      predictionPrec=rbind(predictionPrec,prediction.tmp)
      predictionPlants=rbind(predictionPlants,prediction.tmp)
    }
  }
  return(predictionPlants)
}
write.table(prediction.temperature, file="tmp.prediction2.temp.csv", sep=",")


## PRECIPITATION
pb <- txtProgressBar(0, length(prec_levels), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)

prediction.precipitation <- foreach(prec_increase=prec_levels, .options.snow=opts , .combine="rbind",  .packages=c("phyloseq","tidyr","speedyseq","dplyr")) %dopar% {
  # Global variables
  rm("temp_increase","siteY","site.phyloseq","site.melt","plant.p","plant","prediction.tmp","plant.species.meta")
  .GlobalEnv$samples <- samples
  .GlobalEnv$sample_names <- sample_names
  temp_increase = sqrt(prec_increase)

  predictionPlants=data.frame()
  site.count=c()
  # Iterate over sites
  for (siteY in sample_names){
    site.otu <- tmp.otu[,siteY]
    site.phyloseq <- merge_phyloseq(site.otu,tmp.smp,tmp.tax)
    
    site.melt <- psmelt(site.phyloseq)
    site.melt <- site.melt[site.melt$Abundance > 0,]
    site.count=c(site.count,siteY)

    # iterate over plants
    for (plant.p in site.melt$species){
      plant.species.meta = plant.temp[plant.temp$species == plant.p,]
      
      prediction.tmp <- data.frame(species = plant.p,
                                   site = site.melt$Sample[1],
                                   date = site.melt$Month[1],
                                   country = site.melt$Country[1],
                                   coordY = site.melt$CoordY[1],
                                   abundance = site.melt[site.melt$species == plant.p,"Abundance"],
                                   prediction = NA,
                                   p.dist.value = pnorm(sqrt(site.melt$Prec[1])-temp_increase, mean=plant.species.meta$prec.mean, sd=plant.species.meta$prec.sd),
                                   endemism = plant.species.meta$countries,
                                   temp_increase = temp_increase,
                                   crop = plant.species.meta$crop)
      
      if (plant.species.meta$sites > 10){
        if (plant.species.meta$prec.q01 > sqrt(site.melt$Prec[1])-temp_increase){
          prediction.tmp$prediction <- "<q01"
        }else if (plant.species.meta$prec.q05 > sqrt(site.melt$Prec[1])-temp_increase){
          prediction.tmp$prediction <- "<q05"
        }else if (plant.species.meta$prec.q10> sqrt(site.melt$Prec[1])-temp_increase){
          prediction.tmp$prediction <- "<q10"
        }else{
          prediction.tmp$prediction = ">=q10"
        }
        # maybe real q value? comes later...
        rm("plant.p","plant")
        
      }
      predictionPrec=rbind(predictionPrec,prediction.tmp)
      predictionPlants=rbind(predictionPlants,prediction.tmp)
    }
  }
  return(predictionPlants)
}

stopCluster(cl)
write.table(prediction.precipitation, file="tmp.prediction2.prec.csv", sep=",")
