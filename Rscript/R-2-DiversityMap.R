
################################################
# Data analysis 1 
## Create diversity distribution plot
print("Create diversity distribution plot")

### calcualate diversity
diversity <- diversity(t(otu_table(samples)))
diversity <- data.frame((diversity))
colnames(diversity)[1] <- "Shannon"

diversity$SampleID <- rownames(diversity)
sample_data(samples)$SampleID <- rownames(sample_data(samples))

### combine diversity with metadata
diversity <- merge(diversity,data.frame(sample_data(samples)), by = "SampleID")

### get map
world <- ne_countries(scale = 50, returnclass = 'sf')
world <- map_data("world", regions=c("Czech Republic",names(table(diversity$Country))))

### plot diversity
pdf("plots.supplement/diversity_map.pdf", width=7, height=6.5)
ggplot() + geom_polygon(data = world, aes(x = long, y = lat, group = group), fill = "lightgray", colour = "white") + 
  theme_bw()+
  geom_jitter(data=diversity,aes(x=CoordX,y=CoordY, col=Shannon), width=0.2, height=0.2, size=2, alpha=0.5)+
  scale_color_viridis(option="magma") +   
  ylab("Latitude") + xlab("Longitude")

dev.off()

### plot diversity & latitiude 

pdf("plots.supplement/diversity_latitude_quadratic.pdf", width=5.5, height=5)
ggplot(diversity, aes(y=exp(Shannon), x=CoordY, fill=CoordY))+
 geom_boxplot(aes(group=as.factor(CoordY)), linewidth=0.2, outlier.size=0.1)+
  #geom_point()+
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), col="black")+
  scale_color_viridis(discrete=F, option="turbo", direction=-1)+
  scale_fill_viridis(discrete=F, option="turbo", direction=-1)+
  scale_y_continuous(limits = c(0, 15)) +
  theme_bw() +   
  ylab(bquote("Effective Diversity ("~e^H[2]~")"))+
  xlab("Latitude")
dev.off()

test_results[["diversity_latitude_quadratic"]] <- summary(step(lm(exp(Shannon)~CoordY+I(CoordY^2),data=diversity)))
