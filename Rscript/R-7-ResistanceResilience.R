################################################
# Data analysis 4
## Combined effects: Resistance & Resilience
print("Combined effects: Resistance & Resilience")

# getting the range of potential distributions
range.prec <- 0:(max(round(data.map$Prec, digits=0))/5)*5
range.temp <- 0:(max(round(data.map$Temp+ max(temp_levels), digits=0)))
  
cores=detectCores()
cl <- makeCluster(cores[1]-2) #not to overload your computer
registerDoSNOW(cl)

print("Plant precipitation probability distribution precalculation")
pb <- txtProgressBar(0, length(range.prec), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)  

prediction.prec.exact <- foreach(prec_increase=range.prec, .options.snow=opts , .combine="cbind",  .packages=c("phyloseq","tidyr","speedyseq","dplyr")) %dopar% {
  plant.qs = matrix(nrow = length(plant_taxa), ncol=1)
  colnames(plant.qs)[1] <- prec_increase
  rownames(plant.qs) <- plant_taxa
  
  plant.qs <- as.matrix(plant.qs)
  for (plant in plant_taxa){
    plant.species.meta = plant.temp[plant.temp$species == plant,]
    
    if(length(plant.species.meta$prec.mean)>0){
      p.dist.value = pnorm(sqrt(prec_increase), mean=plant.species.meta$prec.mean, sd=plant.species.meta$prec.sd)
    }else(
      p.dist.value <- NA
    )
    plant.qs[plant,1] <- p.dist.value
  }
  plant.qs
}

print("Plant temperature probability distribution precalculation")
pb <- txtProgressBar(0, length(range.temp), style = 3)
prediction.temp.exact <- foreach(prec_increase=range.temp, .options.snow=opts , .combine="cbind",  .packages=c("phyloseq","tidyr","speedyseq","dplyr")) %dopar% {
    plant.qs = matrix(nrow = length(plant_taxa), ncol=1)
    colnames(plant.qs)[1] <- prec_increase
    rownames(plant.qs) <- plant_taxa
    
    plant.qs <- as.matrix(plant.qs)
    for (plant in plant_taxa){
      plant.species.meta = plant.temp[plant.temp$species == plant,]
      
      if(length(plant.species.meta$prec.mean)>0){
        p.dist.value = pnorm((prec_increase), mean=plant.species.meta$temp.mean, sd=plant.species.meta$temp.sd,lower.tail = F)
      }else(
        p.dist.value <- NA
      )
      plant.qs[plant,1] <- p.dist.value
    }
    plant.qs
}


print("Plant geographic distances")
geodist <- distm(sample_data(data.species.rel.filter)[,c("CoordX","CoordY")], fun = distHaversine)
geodist <- rescale(sqrt(geodist), c(0,1))
row.names(geodist) <-rownames(sample_data(data.species.rel.filter))
colnames(geodist) <- rownames(sample_data(data.species.rel.filter))     
      

data.map$site = interaction(data.map$CoordY, data.map$CoordX, sep = "-" )
pb <- txtProgressBar(0, length(data.map$site ), style = 3)

sites <- unique(rownames(data.map))
distance.geo <- foreach(sample=sites, .options.snow=opts , .combine="cbind",  .packages=c("phyloseq","tidyr","speedyseq","dplyr")) %dopar% {
  plant.qs = matrix(nrow = length(plant_taxa), ncol=1)
  colnames(plant.qs)[1] <- sample
  rownames(plant.qs) <- plant_taxa
  
  for (plant in plant_taxa){
    sites.present1 <- colnames(otu_table(data.species.rel.filter)[,otu_table(data.species.rel.filter)[plant,]>0])
    min.distance <- min(geodist[sites.present1,sample])

    plant.qs[plant,1] <- min.distance
  }
  plant.qs
}

print("Plant taxonomic distances")
plant.df <- merge(plant.temp, data.frame(tax_table(data.species.rel.filter)), by="species")      
rownames(plant.df) <- plant.df$species
plant.df$endemism <- rescale(sqrt(1/plant.df$countries), c(0,1))
      

plant_taxa2 <- plant_taxa

pb <- txtProgressBar(0, length(plant_taxa2), style = 3)
distance.tax <-foreach (plant = plant_taxa2, .combine=cbind, .options.snow=opts, .packages=c("geosphere","phyloseq","scales")) %dopar% {
 print(plant)
  plant.qs = matrix(nrow = length(plant_taxa2), ncol=1)
  colnames(plant.qs)[1] <- plant
  rownames(plant.qs) <- plant_taxa2
  
  for (plant.replacement in plant_taxa2){
    tax.dist <- NA
    if (!is.na(plant.df[plant.replacement,"genus"]) && !is.na(plant.df[plant,"genus"])){ 
      if (plant == plant.replacement){
        tax.dist <- 0
      }else if (plant.df[plant,"genus"] == plant.df[plant.replacement,"genus"]){
        tax.dist <- 0.2
      }else if (plant.df[plant,"family"] == plant.df[plant.replacement,"family"]){
        tax.dist <- 0.4
      }else if (plant.df[plant,"order"] == plant.df[plant.replacement,"order"]){
        tax.dist <- 0.6
      }else if (plant.df[plant,"phylum"] == plant.df[plant.replacement,"phylum"]){
        tax.dist <- 0.8
      }else{
        tax.dist <- 1
      }
    plant.qs[plant.replacement,1] <- tax.dist
    }
  }  
  plant.qs
}

print("Abundance estimation") # not used 

pb <- txtProgressBar(0, length(range.prec), style = 3)
prediction.abun.exact <-foreach (plant = plant_taxa, .combine=rbind, .options.snow=opts, .packages=c("geosphere","phyloseq","scales")) %dopar% {
  if (sum(tmp.otu[plant,]) > 0){
    
    taxon.otu <- tmp.otu[plant,tmp.otu[plant,]>0]
    taxon.mean.abundance <-  mean(as.numeric(taxon.otu))
    
  }else{
    taxon.otu <- NA
  }
  data.frame(row.names=plant, mean.abund = taxon.mean.abundance)
  }
prediction.abun.exact$mean.abund <- rescale(prediction.abun.exact$mean.abund, to=c(0,1))
prediction.abun.exact2 <- rescale(t((prediction.abun.exact)), to=c(0,1))

stopCluster(cl)


### Predict Resi Resi
temp_levels = seq(0,5,0.5)
prec_levels = seq(0,40,5)
cl <- makeCluster(cores[1]-2) #not to overload your computer
registerDoSNOW(cl)

print("Calculate site resistance and resilience")  

pb <- txtProgressBar(0, length(prec_levels)*length(temp_levels), style = 3)
prediction.resiresi <- foreach(prec_increase=prec_levels, .options.snow=opts , .combine="rbind",  .packages=c("phyloseq","tidyr","speedyseq","dplyr")) %:% 
  foreach(temp_increase=temp_levels, .options.snow=opts , .combine="rbind",  .packages=c("phyloseq","tidyr","speedyseq","dplyr")) %dopar% {
    .GlobalEnv$samples <- samples
    .GlobalEnv$sample_names <- sample_names
    
    predictionPlants=data.frame()
    site.count=c()
    site.prediction = data.frame()
    for (siteY in sample_names){
      site.otu <- tmp.otu[,siteY]
      site.phyloseq <- merge_phyloseq(site.otu,tmp.smp,tmp.tax)
      
      site.melt <- psmelt(site.phyloseq)
      site.melt <- site.melt[site.melt$Abundance > 0,]
      site.expected.temp <- as.character(round(site.melt$Temp[1]+as.numeric(temp_increase), digits=0))
      site.expected.prec <- as.character(round((site.melt$Prec[1]-as.numeric(prec_increase))/5, digits=0)*5)
      
      if (site.expected.prec < 0){site.expected.prec <- 0}
           
      plants.in.plot <- site.melt$species
      scores.in.plot <- c()
      scores.in.plot.rep <- c()
      
      rep.temp <- prediction.temp.exact[,as.character(site.expected.temp)]
      rep.prec <- prediction.prec.exact[,as.character(site.expected.prec)]
      rep.geo <- distance.geo[,siteY]
      
      sum(names(rep.temp) == names(rep.prec)) == length(names(rep.temp))
      sum(names(rep.temp) == names(rep.geo)) == length(names(rep.temp))
      
      min.dist.geo <- 0.20 # equals 300 km radius
      min.prob <- 0.15 # equals 15 % probability for prec and temp
      min.tax <- 0.5 # equals same family
      
      rep.prec[plants.in.plot]
      
      # resitance: plants/all dropping out
       resistance.plants <- rep.temp[plants.in.plot] > min.prob & rep.prec[plants.in.plot] > min.prob
       resistance.plot <- sum(resistance.plants)
       resistance <- resistance.plot/length(plants.in.plot)
       
      # resilience potential: new candidates that might serve as repacements
       # closeby
       closeby.candidates <- names(rep.geo[rep.geo != 0 & rep.geo < min.dist.geo])
       resilience.closeby <- sum(rep.temp[closeby.candidates] > min.prob & rep.prec[closeby.candidates] > min.prob)
       
       resilience <- resilience.closeby/length(closeby.candidates) # without other criteria
       
       rep.full.temp <- rep.temp[closeby.candidates]
       rep.full.prec <- rep.prec[closeby.candidates] 

       rep.plant.resilience.total = 0
       rep.plant.resilience.dropouts = 0
       rep.plant.dropouts = plants.in.plot[!(resistance.plants)]
       
       
      for (plant in site.melt$species){ # find replacement taxa that match all criteria
        rep.tax <- distance.tax[plant,]
        
        rep.full.temp <- rep.full.temp[!(names(rep.full.temp) %in% plants.in.plot)]
        rep.full.prec <- rep.full.prec[!(names(rep.full.prec) %in% plants.in.plot)]
        
        family.members <- names(rep.tax[rep.tax < min.tax])
        rep.full.temp2 <- rep.full.temp[names(rep.full.temp) %in% family.members]
        rep.full.prec2 <- rep.full.prec[names(rep.full.prec) %in% family.members]
        
        closeby.family.candidates <- names(rep.full.prec2)
        resilience.family.closeby <- sum(rep.full.temp2[closeby.family.candidates] > min.prob & rep.prec[closeby.family.candidates] > min.prob)
        
        if (resilience.family.closeby > 0 ){
          rep.plant.resilience.total = rep.plant.resilience.total + 1
          if(!(resistance.plants[plant])){
            rep.plant.resilience.dropouts = rep.plant.resilience.dropouts + 1
          }
        }
        
      } # for plant
      
      #resilience.plant <- rep.plant.resilience/length(plants.in.plot)
      resilience.dropouts <- (rep.plant.resilience.dropouts/length(rep.plant.dropouts))*(1-resistance)
      resilience.dropouts2 <- (rep.plant.resilience.dropouts/length(rep.plant.dropouts))
      resilience2 <- resilience *(1-resistance)
      
      if (length(rep.plant.dropouts)==0){ resilience.dropouts2 <- resilience.dropouts <- 0}
      

      site.prediction <- rbind(site.prediction, data.frame(site=siteY,richness=length(plants.in.plot),prec_increase=prec_increase, temp_increase=temp_increase, resistance = resistance, resilience = resilience, resilience2=resilience2, resilience.dropouts2=resilience.dropouts2,resilience.dropouts= resilience.dropouts))  
    } # for site
    return(site.prediction)
  }


#### aggregation
metadata <- sample_data(data.species)
metadata.xy <- data.frame(metadata[,c("CoordX","CoordY")])
metadata.env <- data.frame(metadata[,c("Temp","Vegetation")])
metadata$site <- rownames(metadata)

prediction.all <- merge(prediction.resiresi,data.frame(metadata), by="site")
prediction.all$site2 <- interaction(prediction.all$CoordY,prediction.all$CoordX)

prediction.all$month2 <- factor(gsub("^\\d\\d","",prediction.all$Month))

prediction.all.mean <- prediction.all %>%
  group_by(month2, CoordY, temp_increase, prec_increase) %>%
  summarize(mean.resistance = mean(resistance, na.rm=T),mean.resilience.dropouts = mean(resilience.dropouts2, na.rm=T))  

prediction.all.mean$coordY2 = round(prediction.all.mean$CoordY*2/10, digits=0)*5



### plot
prediction.all.mean$prec_increase2 <- paste("-",prediction.all.mean$prec_increase,"mm",sep="")
prediction.all.mean$prec_increase2 <- factor(prediction.all.mean$prec_increase2,levels=unique(prediction.all.mean$prec_increase2[order(prediction.all.mean$prec_increase,decreasing =F)]), ordered=T)
prediction.all.mean$temp_increase2 <- paste("+",prediction.all.mean$temp_increase,"°C",sep="")
prediction.all.mean$temp_increase2 <- factor(prediction.all.mean$temp_increase2,levels=unique(prediction.all.mean$temp_increase2[order(prediction.all.mean$temp_increase,decreasing =T)]), ordered=T)

#### full
pdf("resistance_resilience_full.pdf", width=15, height=15)

ggplot(prediction.all.mean)+
  geom_point(aes(y=mean.resilience.dropouts, x=mean.resistance, col=mean.resistance+mean.resilience.dropouts ), size=2, alpha=0.8)+
  facet_grid(temp_increase2~prec_increase2)+
  scale_x_continuous(limits=c(0,1), breaks=c(0,0.5))+
  scale_y_continuous(limits=c(0,1), breaks=c(0,0.5))+
  scale_color_viridis(option="inferno", direction=-1, begin=0, end=0.8)+theme_bw()+  
  theme(panel.margin.y = unit(0, "lines"),
        strip.text.y = element_text(angle = 360, hjust = 0),
        strip.background =element_rect(fill="white"))+
  theme(legend.position="bottom")+  
  xlab("Resistance Potential")+
  ylab("Resilience Potential")+
  coord_fixed()

dev.off()

#### subset
pdf("resistance_resilience_sub.pdf", width=10, height=9)

ggplot(prediction.all.mean[prediction.all.mean$temp_increase %in% 1:5 & prediction.all.mean$prec_increase %in% subset_levels[2:6],])+
  geom_point(aes(y=mean.resilience.dropouts, x=mean.resistance, col=mean.resistance+mean.resilience.dropouts ), size=2, alpha=0.8)+
  facet_grid(temp_increase2~prec_increase2)+
  scale_x_continuous(limits=c(0,1), breaks=c(0,0.5))+
  scale_y_continuous(limits=c(0,1), breaks=c(0,0.5))+
  scale_color_viridis(option="inferno", direction=-1, begin=0, end=0.8)+theme_bw()+  
  theme(panel.margin.y = unit(0, "lines"),
        strip.text.y = element_text(angle = 360, hjust = 0),
        strip.background =element_rect(fill="white"))+
  theme(legend.position="bottom")+  
  xlab("Resistance Potential")+
  ylab("Resilience Potential")+
  coord_fixed()

dev.off()

file.rename("resistance_resilience_sub.pdf", paste("plots","resistance_resilience_sub.pdf",sep="/"))
file.rename("resistance_resilience_full.pdf", paste("plots.supplement","resistance_resilience_full.pdf",sep="/"))

