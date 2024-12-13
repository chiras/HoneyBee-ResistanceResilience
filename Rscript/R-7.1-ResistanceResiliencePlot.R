################################################
# Data analysis 4
## Combined effects: Resistance & Resilience. Plotting
print("Combined effects: Resistance & Resilience Plot")

prediction.all.mean$prec_increase2 <- paste("-",prediction.all.mean$prec_increase,"mm",sep="")
prediction.all.mean$prec_increase2 <- factor(prediction.all.mean$prec_increase2,levels=unique(prediction.all.mean$prec_increase2[order(prediction.all.mean$prec_increase,decreasing =F)]), ordered=T)
prediction.all.mean$temp_increase2 <- paste("+",prediction.all.mean$temp_increase,"°C",sep="")
prediction.all.mean$temp_increase2 <- factor(prediction.all.mean$temp_increase2,levels=unique(prediction.all.mean$temp_increase2[order(prediction.all.mean$temp_increase,decreasing =T)]), ordered=T)

prediction.all.mean2 <- prediction.all.mean[prediction.all.mean$temp_increase %in% 1:5 & prediction.all.mean$prec_increase %in% c(5,15,25,35,45),]
pdf("resistance_resilience_sub_density.pdf", width=7, height=8)

ggplot(prediction.all.mean2, aes(y=(mean.resilience.dropouts), x=(mean.resistance))) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = F, n = 300,
    h = c(1, 1)) +
  scale_fill_viridis_b() +  # Optional: Change color scale to viridis
  #geom_point(alpha = 0.3, size = 1) +  # Add scatter points for clarity
  theme_minimal() +
  facet_grid(temp_increase2~prec_increase2)+
  scale_x_continuous(limits=c(0,1), breaks=c(0,0.5))+
  scale_y_continuous(limits=c(0,1), breaks=c(0,0.5,1))+
  theme_bw()+  
  theme(panel.margin.y = unit(0, "lines"),
        strip.text.y = element_text(angle = 360, hjust = 0),
        strip.background =element_rect(fill="white"))+
  labs(fill = "Density")+  
  theme(legend.position="bottom")+  
  xlab("Resistance Potential")+
  ylab("Resilience Potential")+
  geom_vline(xintercept = 0.5, color = "red", linetype = "dashed") +
  geom_hline(yintercept = 0.5, color = "red", linetype = "dashed")
dev.off()

for (i in 1:4){
month_num <- c("/05/2023" ,"/06/2023" ,"/07/2023", "/08/2023")[i]
month_char <- c("May" ,"June" ,"July", "August")[i]
  
prediction.all.mean2 <- prediction.all.mean[prediction.all.mean$temp_increase %in% 1:5 & prediction.all.mean$prec_increase %in% c(5,15,25,35,45) & prediction.all.mean$month2==month_num,]

pdf(paste("plots.supplement/resistance_resilience_sub_density_",month_char,".pdf",sep=""), width=7, height=8)

ggplot(prediction.all.mean2, aes(y=(mean.resilience.dropouts), x=(mean.resistance))) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = F, n = 300,
                  h = c(1, 1)) +
  scale_fill_viridis_b() +  # Optional: Change color scale to viridis
  #geom_point(alpha = 0.3, size = 1) +  # Add scatter points for clarity
  theme_minimal() +
  facet_grid(temp_increase2~prec_increase2)+
  scale_x_continuous(limits=c(0,1), breaks=c(0,0.5))+
  scale_y_continuous(limits=c(0,1), breaks=c(0,0.5,1))+
  theme_bw()+  
  theme(panel.margin.y = unit(0, "lines"),
        strip.text.y = element_text(angle = 360, hjust = 0),
        strip.background =element_rect(fill="white"))+
  labs(fill = "Density")+  
  theme(legend.position="bottom")+  
  xlab("Resistance Potential")+
  ylab("Resilience Potential")+
  geom_vline(xintercept = 0.5, color = "red", linetype = "dashed") +
  geom_hline(yintercept = 0.5, color = "red", linetype = "dashed")
dev.off()
}
