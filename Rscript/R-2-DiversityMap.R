oldw <- getOption("warn")
options(warn = -1) # temporarily remove warnings, in case dir exists and graphics warnings

################################################
# Data analysis 1 
## Create diversity distribution plot
print("Create diversity distribution plot")

### calcualate diversity
diversity <- exp(diversity(t(otu_table(samples))))
diversity2 <- (diversity(t(otu_table(samples))))
diversity <- data.frame((diversity))
colnames(diversity)[1] <- "ExpShannon"
diversity$Shannon <- diversity2
diversity$SampleID <- rownames(diversity)

sample_data(samples)$SampleID <- rownames(sample_data(samples))

### combine diversity with metadata
diversity <- merge(diversity,data.frame(sample_data(samples)), by = "SampleID")

### get map
world <- ne_countries(scale = 50, returnclass = 'sf')
world <- map_data("world", regions=c("Czech Republic",names(table(diversity$Country))))

diversity$date2 <- as.Date(diversity$Month, format = "%d/%m/%Y")
diversity$month  <- format(diversity$date2, "%m")
diversity$week  <- strftime(diversity$date2, format = "%V")

diversity <- diversity[!is.na(diversity$month),]

diversity$month[diversity$month=="05"] <- "May"
diversity$month[diversity$month=="06"] <- "June"
diversity$month[diversity$month=="07"] <- "July"
diversity$month[diversity$month=="08"] <- "August"

diversity$month <- ordered(factor(diversity$month, levels = c("May", "June", "July","August")))

### plot diversity
pdf("plots.supplement/diversity_map.pdf", width=8, height=6)
ggplot(data=diversity,aes(x=CoordX,y=CoordY, col=ExpShannon)) + geom_polygon(data = world, aes(x = long, y = lat, group = group), fill = "lightgray", colour = "white") + 
  theme_bw()+
  geom_jitter(width=0.2, height=0.2, size=2, alpha=0.5)+
  scale_color_viridis(option="magma") +   
  ylab("Latitude") + xlab("Longitude")+
  facet_wrap(month~., ncol=2)+
    theme(panel.margin.y = unit(0, "lines"),
                                         strip.background =element_rect(fill="white"),
                                         panel.grid.minor = element_blank())+ 
  labs(colour = bquote("Effective Diversity ("~e^H[2]~")"))



dev.off()

### plot diversity & latitiude 

pdf("plots.supplement/diversity_latitude_quadratic.pdf", width=8, height=6)
ggplot(diversity, aes(y=ExpShannon, x=CoordY, fill=(CoordY)))+
 geom_boxplot(aes(group=as.factor(CoordY)), linewidth=0.2, outlier.size=0.1)+
  #geom_point()+
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), col="black")+
  scale_color_viridis(discrete=F, option="turbo", direction=-1)+
  scale_fill_viridis(discrete=F, option="turbo", direction=-1)+
 scale_y_continuous(limits = c(0, 15)) +
  theme_bw() +     facet_wrap(month~., ncol=2)+
    theme(panel.margin.y = unit(0, "lines"),
                                         strip.background =element_rect(fill="white"),
                                         panel.grid.minor = element_blank())+

  ylab(bquote("Effective Diversity ("~e^H[2]~")"))+
  xlab("Latitude")+  labs(fill = "Latitude")

dev.off()

test_results[["diversity_latitude_quadratic"]] <- summary(step(lm(exp(Shannon)~CoordY+I(CoordY^2)*Temp+Prec+as.numeric(month),data=diversity)))
test_results[["diversity_latitude_quadratic_5"]] <- summary(step(lm(exp(Shannon)~CoordY+I(CoordY^2),data=diversity[diversity$month=="May",])))
test_results[["diversity_latitude_quadratic_6"]] <- summary(step(lm(exp(Shannon)~CoordY+I(CoordY^2),data=diversity[diversity$month=="June",])))
test_results[["diversity_latitude_quadratic_7"]] <- summary(step(lm(exp(Shannon)~CoordY+I(CoordY^2),data=diversity[diversity$month=="July",])))
test_results[["diversity_latitude_quadratic_8"]] <- summary(step(lm(exp(Shannon)~CoordY+I(CoordY^2),data=diversity[diversity$month=="August",])))
options(warn = oldw)
