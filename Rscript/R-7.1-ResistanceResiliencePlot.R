##############################################################
# R-7.1-ResistanceResiliencePlot.R
# Combined effects: estimating site-level resistance and 
# resilience to climate extremes using KDEs.
#
# Dependencies:
#   - Requires:
#       - one time running of R-7 for intermediate file intermediate.data/sites_resi-resi.csv
# Outputs:
#   - Plots (sub and full)
##############################################################

print("Combined effects: Resistance & Resilience Plot")

# ------------------------------------------------------------
# Load resistance/resilience data for all sites
# ------------------------------------------------------------
prediction.all <- read.csv(
  "intermediate.data/sites_resi-resi.csv",
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# Factorization of variables
# ------------------------------------------------------------
prediction.all$prec_increase2 <- paste("-",prediction.all$prec_increase,"mm",sep="")
prediction.all$prec_increase2 <- factor(prediction.all$prec_increase2,levels=unique(prediction.all$prec_increase2[order(prediction.all$prec_increase,decreasing =F)]), ordered=T)
prediction.all$temp_increase2 <- paste("+",prediction.all$temp_increase,"°C",sep="")
prediction.all$temp_increase2 <- factor(prediction.all$temp_increase2,levels=unique(prediction.all$temp_increase2[order(prediction.all$temp_increase,decreasing =T)]), ordered=T)

# ------------------------------------------------------------
# Subsetting for partial plots
# ------------------------------------------------------------
prediction.all.sub <- prediction.all[prediction.all$temp_increase %in% 0:5 & prediction.all$prec_increase %in% c(0,10,20,30,40,50),]

pdf("plots/resistance_resilience_sub_density_all_raw.pdf", width=7, height=8)
ggplot(prediction.all.sub, aes(y=(resilience_raw), x=(resistance))) +
  stat_density_2d(aes(fill = after_stat(density)), geom = "raster", contour = F, n = 300, 
    h = c(1.5, 1.5) ) +  
  scale_fill_viridis_b(n.breaks = 6, option="mako") +  
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

# ------------------------------------------------------------
# Full supplement plots
# ------------------------------------------------------------

pdf("plots.supplement/resistance_resilience_full_density_all_raw.pdf", width=14, height=15)
ggplot(prediction.all, aes(y=(resilience_raw), x=(resistance))) +
  stat_density_2d(aes(fill = after_stat(density)), geom = "raster", contour = F, n = 300, 
                  h = c(1.25, 1.25) ) + 
  scale_fill_viridis_b(n.breaks = 6) +  
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

##### Curves Plot Suppl. 9

prediction.all.mean2 <- prediction.all[prediction.all$temp_increase %in% 0:5 & prediction.all$prec_increase %in% c(0,10,20,30,40,50),]
prediction.all.mean2$resiresi =  prediction.all.mean2$resilience_raw+prediction.all.mean2$resistance
prediction.all.mean2$date2 <- as.Date(prediction.all.mean2$Month, format = "%d/%m/%Y")
prediction.all.mean2$month2 <- format(prediction.all.mean2$date2, "%m")

prediction.all.mean3 <- prediction.all.mean2 %>%
  group_by(Lat, resiresi, month2,temp_increase2,prec_increase2) %>%
  summarize(resiresi = mean(resiresi, na.rm = TRUE), .groups = "drop")

pdf("plots.supplement/resistance_resilience_sub_curves.pdf", width=7, height=8)
ggplot(prediction.all.mean3)+
  geom_point(aes(y=(resiresi), x=(Lat), col=(month2),fill=(month2) ), shape = 21,size = 0.5, alpha=0.2, stroke = 0.2)+
  facet_grid(temp_increase2~prec_increase2)+
  scale_y_continuous(limits=c(0,1), breaks=c(0,0.5,1))+
  scale_color_viridis(option="viridis",  discrete=T)+ #, end=0.8
  scale_fill_viridis(option="viridis",  discrete=T)+ #, end=0.8
  theme_bw()+  
  theme(panel.margin.y = unit(0, "lines"),
        strip.text.y = element_text(angle = 360, hjust = 0),
        strip.background =element_rect(fill="white"))+
  #labs(fill = "Resistance Potential + Resilience Potential")+  
  theme(legend.position="bottom")+  
  xlab("Latitude")+
  #ylab("Resistance + Resilience")+
  geom_smooth(aes(y=(resiresi), x=(Lat), col=factor(month2),fill=factor(month2) ))

#coord_fixed()

dev.off()


