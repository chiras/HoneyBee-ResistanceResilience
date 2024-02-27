
### Setting output directory, summarize results and adjust meta-data
    # run once for temperature, once for precipitation
    # return point 1

#toAnalyze = "results.temp" # results.prec or results.temp
dir.create(toAnalyze)
print(paste("Estimate future scenarios:", toAnalyze))

if(toAnalyze == "results.prec" ){
  prediction <- prediction.precipitation[!is.na(prediction.precipitation$prediction),]
}else{
  prediction <- prediction.temperature[!is.na(prediction.temperature$prediction),]
}

prediction$incidence <- 1*(prediction$abundance>0)


#### summarize results over latitudes

predictionSum <- prediction[prediction$temp_increase==0,] %>%
  group_by(coordY) %>%
  summarize(countrysum = sum(abundance, na.rm=T))

countrycodes <- read.table("Data/countrycodes.csv", sep=",", header=T)
rownames(countrycodes)<- countrycodes$name
prediction$countrycodes <- countrycodes[prediction$country,"alpha.2"]


predictionSum2 <- prediction[prediction$temp_increase==0,] %>%
  group_by(coordY) %>%
  summarize(countriesX= paste(unique(country), sep=",", collapse=", "))

predictionSum3 <- prediction[prediction$temp_increase==0,] %>%
  group_by(coordY) %>%
  summarize(countrycodesX= paste(unique(countrycodes), sep=",", collapse=","))


prediction2 <- prediction %>%
  group_by(coordY,prediction,crop,temp_increase) %>%
  summarize(accumulated = sum(abundance, na.rm=T))

countryY <- prediction %>%
  group_by(country) %>%
  summarize(coordY = mean(coordY, na.rm=T))


#### collate with metadata
prediction2.1 <- merge(prediction2,predictionSum, by="coordY")
prediction2.2 <- merge(prediction2.1,predictionSum2, by="coordY")
prediction2 <- merge(prediction2.2,predictionSum3, by="coordY")

#### recalibrate relative data with overall counts
prediction2$accumulated.rel <- prediction2$accumulated/prediction2$countrysum

#### adjust metadata
prediction2$crop[prediction2$crop] <- "crop"
prediction2$crop[prediction2$crop=="FALSE"] <- "wild"

prediction2$prediction <- factor(prediction2$prediction)
prediction2$crop <- factor(prediction2$crop)
#prediction2$country <- factor(prediction2$country)
#prediction2$country <- factor(prediction2$country, levels=unique(factor(countryY$country[order(countryY$coordY)])), ordered=TRUE)

prediction2$temp_increase <- factor(prediction2$temp_increase)
prediction2$coordY <- factor(prediction2$coordY, levels=sort(as.numeric(unique(prediction2$coordY)), decreasing = T), ordered=TRUE)
prediction2$iPC <- interaction(prediction2$prediction, prediction2$crop)
prediction2$iPC <- factor(prediction2$iPC, levels=sort(as.character(unique(prediction2$iPC))), ordered=TRUE)
prediction2$iCC <- factor(interaction(prediction2$coordY,prediction2$countriesX, sep=": "))
prediction2$iCC <- factor(prediction2$iCC, levels=sort(as.character(unique(prediction2$iCC)), decreasing =T), ordered=TRUE)
prediction2$iCC2 <- factor(interaction(prediction2$coordY,prediction2$countrycodesX, sep=": "))
prediction2$iCC2 <- factor(prediction2$iCC2, levels=sort(as.character(unique(prediction2$iCC2)), decreasing =T), ordered=TRUE)

prediction2$prediction <- factor(prediction2$prediction, levels=sort(as.character(unique(prediction2$prediction)), decreasing =T), ordered=TRUE)

#### prepare risk plot, adjust data with respect to temperature or precipitation
if(toAnalyze == "results.prec"){
  reverse = F
  prediction2$temp_increase <-  factor(round(as.numeric(as.character(prediction2$temp_increase))^2, digits=0))
  prediction$temp_increase <-  factor(round(as.numeric(as.character(prediction$temp_increase))^2, digits=0))
  
  prediction2$prediction <- factor(prediction2$prediction, levels=sort(as.character(unique( prediction2$prediction)), decreasing =F), ordered=T)
  prediction2$iPC <- factor(prediction2$iPC, levels=sort(as.character(unique(prediction2$iPC)), decreasing =T), ordered=TRUE)
  criterion=">=q10"
  prediction2$temp_increase2 <- paste("-",prediction2$temp_increase,"mm",sep="")
  prediction2$temp_increase2 <- factor(prediction2$temp_increase2,levels=unique(prediction2$temp_increase2[order(prediction2$temp_increase)]), ordered=T)
  subset_levels <- seq(4,20,4)

}else{
  reverse = T
  criterion="<=q90" 
  prediction2$temp_increase2 <- paste("+",prediction2$temp_increase,"°C",sep="")
  prediction2$temp_increase2 <- factor(prediction2$temp_increase2,levels=unique(prediction2$temp_increase2[order(prediction2$temp_increase)]), ordered=T)
  subset_levels <- seq(0,5,1)
}
cols <- c("#44ce1b","#44ce1b","#f7e379", "#f7e379", "#f2a134", "#f2a134","#e51f1f","#e51f1f")

#### full plot, all scenarios (supplement)

pdf(paste(toAnalyze,"taxa_loss_abundance_latitude_full.pdf",sep="/"), width=15, height=6)
ggplot(prediction2, aes(fill=iPC, shape=crop,y=accumulated.rel, alpha=crop, col="black", x=1)) + 
  geom_bar(position = position_stack(reverse = T),stat="identity",fill="black") +
  geom_bar(position = position_stack(reverse = T),stat="identity",) +
  facet_grid(iCC ~temp_increase2, scales="free", drop=T)+ 
  coord_flip()+
  scale_fill_manual(values =c(cols))+
  scale_color_viridis(option="viridis", discrete=T) + theme_bw()+ 
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
prediction2.1 <- prediction2[prediction2$temp_increase %in% subset_levels,] 

pdf(paste(toAnalyze,"taxa_loss_abundance_latitude_sub.pdf",sep="/"), width=8, height=6.5)
ggplot(prediction2.1, aes(fill=iPC, shape=crop,y=accumulated.rel, alpha=crop, col="black", x=1)) + 
  geom_bar(position = position_stack(reverse = T),stat="identity",fill="black") +
  geom_bar(position = position_stack(reverse = T),stat="identity",) +
  facet_grid(coordY+countrycodesX ~temp_increase2, scales="free", drop=T)+ 
  coord_flip()+
  scale_fill_manual(values =c(cols))+
  scale_color_viridis(option="viridis", discrete=T) + theme_bw()+ 
  scale_alpha_discrete(range = c(0.7, 1))+
  theme_bw()+ 
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

##  repeat from return point 1 for precipitation
