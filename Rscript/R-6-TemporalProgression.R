################################################
# Data analysis 4
## Plot risk per week per latitude
print("Evaluation of temporal risk changes")
oldw <- getOption("warn")
options(warn = -1) # temporarily remove warnings, in case dir exists and graphics warnings

# prepare date formats
prediction$date2 <- as.Date(prediction$date, format = "%d/%m/%Y")
prediction$month  <- format(prediction$date2, "%m")
prediction$week  <- strftime(prediction$date2, format = "%V")

prediction.season <- prediction[!is.na(prediction$prediction),]
prediction.season <- prediction[prediction$prediction!=criterion ,]
prediction.season = prediction.season[!is.na(prediction.season$month),]

# aggregate
prediction.season2 <- prediction.season %>%
  group_by(coordY,temp_increase, month,site) %>%
  summarize(accumulated = sum(abundance, na.rm=T))

prediction.season2$coordY2 = round(prediction.season2$coordY*2/10, digits=0)*5


#fit <- loess(accumulated~as.numeric(month)+coordY2,data=prediction.season2)
#prediction.season2$pred <- predict(fit)
#prediction.season2=prediction.season2[prediction.season2$temp_increase!=0,]

if(toAnalyze == "results.prec" ){
  ylabel <- "% <q10"
  prediction.season2$temp_increase2 <- paste("-",prediction.season2$temp_increase,"mm",sep="")
  prediction.season2$temp_increase2 <- factor(prediction.season2$temp_increase2,levels=unique(prediction.season2$temp_increase2[order(prediction2$temp_increase)]), ordered=T)

}else{
  ylabel <- "% >q90"
  prediction.season2$temp_increase2 <- paste("+",prediction.season2$temp_increase,"°C",sep="")
  prediction.season2$temp_increase2 <- factor(prediction.season2$temp_increase2,levels=unique(prediction.season2$temp_increase2[order(prediction2$temp_increase)]), ordered=T)

}

##### full
pdf(paste(toAnalyze,"season_q90_distribution_full.pdf",sep="/"), width=15, height=10)
ggplot(prediction.season2, aes(x=as.numeric(month),y=accumulated, col=factor(coordY2),fill=factor(coordY2)))+
  geom_smooth(data = prediction.season2 %>% group_by(coordY2, temp_increase) %>% filter(n() > 50), method = "loess", se=T) + 
  facet_wrap(temp_increase2~., ncol=6)+
  coord_cartesian(ylim=c(0, 1))+
  scale_color_viridis(discrete=T, option = "viridis", end=0.8)+
  scale_fill_viridis(discrete=T, option = "viridis", end=0.8, alpha=0.1)+
  theme_bw()+
  scale_x_continuous(breaks = 5:8, expand = c(0, 0.2))+
  scale_y_continuous(expand = c(0, 0))+
  ylab(ylabel) + xlab("Month")+
  theme(panel.margin.y = unit(0, "lines"),
        strip.background =element_rect(fill="white"),
        panel.grid.minor = element_blank())
dev.off()   
     
##### subset

prediction.season2.1 <- prediction.season2[prediction.season2$temp_increase %in% subset_levels,] 

pdf(paste(toAnalyze,"season_q90_distribution_sub.pdf",sep="/"), width=10, height=3)
ggplot(prediction.season2.1, aes(x=as.numeric(month),y=accumulated, col=factor(coordY2),fill=factor(coordY2)))+
  geom_smooth(data = prediction.season2.1 %>% group_by(coordY2, temp_increase) %>% filter(n() > 25), method = "loess", se=T) + 
  facet_wrap(temp_increase2~., ncol=6)+
  coord_cartesian(ylim=c(0, 1))+
  scale_color_viridis(discrete=T, option = "viridis", end=0.8)+
  scale_fill_viridis(discrete=T, option = "viridis", end=0.8, alpha=0.1)+
  theme_bw()+
  scale_x_continuous(breaks = 5:8, expand = c(0, 0.2))+
  scale_y_continuous(expand = c(0, 0))+
  theme(panel.margin.y = unit(0, "lines"),
        strip.background =element_rect(fill="white"),
        panel.grid.minor = element_blank())+ 
  theme(legend.position="bottom")+  
  ylab(ylabel) + xlab("Month")+
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))
  
dev.off()   

file.copy(paste(toAnalyze,"season_q90_distribution_sub.pdf",sep="/"), paste("plots",paste(toAnalyze,"season_q90_distribution_sub.pdf", sep="."),sep="/"))
file.copy(paste(toAnalyze,"season_q90_distribution_full.pdf",sep="/"), paste("plots.supplement",paste(toAnalyze,"season_q90_distribution_full.pdf", sep="."),sep="/"))


options(warn = oldw)

# ANOVA prediction.season2.1
