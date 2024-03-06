
################################################
# Data analysis 1 
## Create diversity distribution plot
print("Create diversity distribution plot")

samples <- subset_samples(data.species.rel.filter, CoordY != 0 & Temp != 0 & Project=="INSIGNIA-EU")

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
  geom_jitter(data=diversity,aes(x=CoordX,y=CoordY, col=Temp), width=0.2, height=0.2, size=2)+
  scale_color_viridis()
dev.off()

