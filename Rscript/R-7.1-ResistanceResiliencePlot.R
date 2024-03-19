################################################
# Data analysis 4
## Combined effects: Resistance & Resilience. Plotting
print("Combined effects: Resistance & Resilience Plot")

prediction.all.mean$prec_increase2 <- paste("-",prediction.all.mean$prec_increase,"mm",sep="")
prediction.all.mean$prec_increase2 <- factor(prediction.all.mean$prec_increase2,levels=unique(prediction.all.mean$prec_increase2[order(prediction.all.mean$prec_increase,decreasing =F)]), ordered=T)
prediction.all.mean$temp_increase2 <- paste("+",prediction.all.mean$temp_increase,"°C",sep="")
prediction.all.mean$temp_increase2 <- factor(prediction.all.mean$temp_increase2,levels=unique(prediction.all.mean$temp_increase2[order(prediction.all.mean$temp_increase,decreasing =T)]), ordered=T)

#### full
pdf("resistance_resilience_full.pdf", width=12, height=13)

ggplot(prediction.all.mean)+
  geom_point(aes(y=(mean.resilience.dropouts), x=(mean.resistance), col=(mean.resistance)+(mean.resilience.dropouts),fill=(mean.resistance)+(mean.resilience.dropouts) ), shape = 21,size = 2,colour = "black", stroke = 0.2)+
  facet_grid(temp_increase2~prec_increase2)+
  scale_x_continuous(limits=c(0,1), breaks=c(0,0.5))+
  scale_y_continuous(limits=c(0,1), breaks=c(0,0.5,1))+
  scale_color_viridis(option="inferno", direction=-1, begin=0, end=0.8)+ #, end=0.8
   scale_fill_viridis(option="inferno", direction=-1, begin=0, end=0.8)+ #, end=0.8
 theme_bw()+  #, end=0.8
  theme(panel.margin.y = unit(0, "lines"),
        strip.text.y = element_text(angle = 360, hjust = 0),
        strip.background =element_rect(fill="white"))+
  theme(legend.position="bottom")+  
  xlab("Resistance Potential")+
  ylab("Resilience Potential")+
  labs(color = "Resistance Potential + Resilience Potential")+  
  coord_fixed()
dev.off()

#### subset
pdf("resistance_resilience_sub.pdf", width=7, height=8)
prediction.all.mean2 <- prediction.all.mean[prediction.all.mean$temp_increase %in% 1:5 & prediction.all.mean$prec_increase %in% subset_levels[2:6],]
ggplot(prediction.all.mean2)+
  geom_point(aes(y=(mean.resilience.dropouts), x=(mean.resistance), col=(mean.resistance)+(mean.resilience.dropouts),fill=(mean.resistance)+(mean.resilience.dropouts) ), shape = 21,size = 2,colour = "black", stroke = 0.2)+
  facet_grid(temp_increase2~prec_increase2)+
  scale_x_continuous(limits=c(0,1), breaks=c(0,0.5))+
  scale_y_continuous(limits=c(0,1), breaks=c(0,0.5,1))+
  scale_color_viridis(option="inferno", direction=-1, begin=0, end=0.8)+ #, end=0.8
  scale_fill_viridis(option="inferno", direction=-1, begin=0, end=0.8)+ #, end=0.8
  theme_bw()+  
  theme(panel.margin.y = unit(0, "lines"),
        strip.text.y = element_text(angle = 360, hjust = 0),
        strip.background =element_rect(fill="white"))+
  labs(fill = "Resistance Potential + Resilience Potential")+  
  theme(legend.position="bottom")+  
  xlab("Resistance Potential")+
  ylab("Resilience Potential")+
  coord_fixed()

dev.off()

file.rename("resistance_resilience_sub.pdf", paste("plots","resistance_resilience_sub.pdf",sep="/"))
file.rename("resistance_resilience_full.pdf", paste("plots.supplement","resistance_resilience_full.pdf",sep="/"))

# plot legend

intervall=0.001

cline <- function(sumxy, intervall){
  x=seq(0,sumxy,intervall)
  y=sumxy-x
  return(data.frame(x=x,y=y))
}


max = max(prediction.all.mean2$mean.resistance)+max(prediction.all.mean2$mean.resilience.dropouts)

x=seq(0,1,intervall)
xy= expand.grid(x,x)
xy= xy[rowSums(xy)<max,]
plot <- ggplot(xy, aes(x=Var1, y=Var2, col=(Var1+Var2)))+  
  geom_point(size=1)+
  theme_bw()

pdf("resistance_resilience_legend.pdf", width=3.5, height=3)
plot +  geom_path(data=cline(0.25, intervall), aes(x=x,y=y), inherit.aes = F, col="black", linewidth=1 )+
  geom_path(data=cline(0.5, intervall), aes(x=x,y=y), inherit.aes = F, col="black", linewidth=1)+
  geom_path(data=cline(0.75, intervall), aes(x=x,y=y), inherit.aes = F, col="black", linewidth=1)+
  geom_path(data=cline(1, intervall), aes(x=x,y=y), inherit.aes = F, col="black", linewidth=1)+
  geom_path(data=cline(1.25, intervall), aes(x=x,y=y), inherit.aes = F, col="black", linewidth=1)+
  geom_path(data=cline(1.5, intervall), aes(x=x,y=y), inherit.aes = F, col="black", linewidth=1)+
  geom_path(data=cline(1.75, intervall), aes(x=x,y=y), inherit.aes = F, col="black", linewidth=1)+
  scale_x_continuous(limits=c(0,1), breaks=c(0,0.5,1))+
  scale_y_continuous(limits=c(0,1), breaks=c(0,0.5))+
  scale_color_viridis(option="inferno", direction=-1, begin=0, end=0.8)+ #, end=0.8
  theme_bw()+  
  labs(color = "Resistance Potential + Resilience Potential")+  
  theme(panel.margin.y = unit(0, "lines"),
        strip.text.y = element_text(angle = 360, hjust = 0),
        strip.background =element_rect(fill="white"))+
  theme(legend.position="bottom")+  
  xlab("Resistance Potential")+
  ylab("Resilience Potential")+
  coord_fixed()
dev.off()

## 

file.rename("resistance_resilience_legend.pdf", paste("plots.supplement","resistance_resilience_legend.pdf",sep="/"))

  

  