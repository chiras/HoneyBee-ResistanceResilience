#todo: matching issues? see differences in 1973 plots and overall

historical.data <- read.table(paste("Data/",measure,"_Data.csv",sep=""), sep=",", header=T)

historical.data.merge <- merge(data.frame(sample_data(samples))[,c("BK","Month","CoordY")],historical.data,by="BK")

hist.melt <- reshape2::melt((historical.data.merge), value.name="Value", id.vars=c("BK","Month","CoordY"))

# get year month format
hist.melt$variable <- gsub("^X","",hist.melt$variable)
hist.melt$year_record <- gsub("\\..*","",hist.melt$variable)
hist.melt$month_record <- gsub(".*\\.","",hist.melt$variable)
hist.melt$month_sample <- gsub("\\/[0-9]*$","",gsub("^[0-9]*\\/","",hist.melt$Month))

# get only valid comparisons of same months
hist.melt2 <- hist.melt[hist.melt$month_record==hist.melt$month_sample,]

# aggregate duplicate point from melt
hist.melt3 <- hist.melt2 %>%
  group_by(BK,variable,CoordY,year_record,month_record,month_sample) %>%
  summarize(value = mean(as.numeric(Value), na.rm=T))

# get current
hist.melt3.current <- hist.melt3[hist.melt3$year_record=="2023",]

hist.melt4 <- merge(hist.melt3, hist.melt3.current,by=c("BK","month_sample"))
hist.melt4$year_record <- as.numeric(hist.melt4$year_record.x)

# historical change over years 
pdf(paste("plots.supplement/",measure,"_historical_change.pdf",sep=""), width=8, height=6)
ggplot(hist.melt4,aes(y=value.x,x=year_record, col=as.factor(CoordY.y)))+
  geom_smooth(method="loess",n=100, aes(fill=as.factor(CoordY.y)), alpha=0.25)+
  geom_smooth(method="loess",n=100, se=F)+
  #geom_point(alpha=0.2,size=1)+
  #facet_wrap(month_record.x ~.)+
  theme_bw()+
  scale_color_viridis(discrete=T)+ #, end=0.8
  scale_fill_viridis(discrete=T)+ #, end=0.8
  xlim(1900,2023)
dev.off()

# difference between 2023 and 1973

hist.melt5 <- hist.melt4[hist.melt4$year_record.x=="1973",]
hist.melt5 <- hist.melt4[hist.melt4$year_record.x<"1976" & hist.melt4$year_record.x>"1969",]

# historical (1973) values are value.x, current values are value.y
if(measure=="Temperature"){
  # difference past - present 
  #-> positive numbers: decrease of Temp, negative numbers: increase of Temp over time
  hist.melt5$difference <- (hist.melt5$value.x) - (hist.melt5$value.y)
  breaks_plot=seq(-10,10,1)
  ymin = 5
  ymin2 = 2.5
}else{
  breaks_plot=seq(-200,+250,25)
  #-> positive numbers: decrease of Prec, negative numbers: increase of Prec over time
  hist.melt5$difference <- ((hist.melt5$value.x) - (hist.melt5$value.y))
  ymin = -25
  ymin2 = +50
}
hist.melt5$difference<- -hist.melt5$difference

pdf(paste("plots.supplement/",measure,"_historical_1973.pdf",sep=""), width=8, height=6)
ggplot(hist.melt5,aes(y=difference,x=CoordY.x, col=as.factor(CoordY.x)))+
  geom_boxplot()+
  geom_smooth(aes(y=difference,x=CoordY.x), inherit.aes = F, col="black", method="loess")+
  #geom_point(alpha=0.2,size=1)+
  facet_wrap(month_record.x ~.)+
  theme_bw()+
  xlab("Latitude") + ylab(paste(measure,"difference since 1973"))+
  scale_y_continuous(breaks=breaks_plot)+
  scale_color_viridis(discrete=T)#+
  #xlim(1900,2023)
dev.off()

hist.melt.6 <- hist.melt5 %>%
  group_by(CoordY.x,month_record.x ) %>%
  summarize(mean = mean(as.numeric(difference), na.rm=T))


test_results[[paste(measure,"historical_change_since_1973_latitude", sep="_")]] <- cor.test(hist.melt.6$mean,as.numeric(hist.melt.6$CoordY.x))

t.tests <- data.frame()
for (month in unique(hist.melt5$month_record.x)){
  for (latitude in unique(hist.melt5$CoordY.x)){
    sub_difference <- hist.melt5[hist.melt5$CoordY.x==latitude & hist.melt5$month_record.x==month,"difference"]
    test <- wilcox.test(sub_difference,mu=0,alternative="two.sided")
    t.tests <- rbind(t.tests,data.frame(CoordY.x=latitude, 
                                        month_record.y=month, 
                                        p.value = round(test$p.value, digits=3), 
                                        est=test$statistic, 
                                        mean=mean(sub_difference), 
                                        max=max(sub_difference), 
                                        min=min(sub_difference),
                                        ue= quantile(sub_difference, probs = 1),
                                        le= quantile(sub_difference, probs = 0),
                                        ue2= median(sub_difference)+ 1.58*IQR(sub_difference)/sqrt(length(sub_difference)),
                                        le2= median(sub_difference)- 1.58*IQR(sub_difference)/sqrt(length(sub_difference))
    ))
    # lower edge of notch = median - 1.58 * IQR / sqrt(n).
    # upper edge of notch = median + 1.58 * IQR / sqrt(n).
  
}}

t.tests$signif2 <- substr(as.character(t.tests$mean),1,1)
t.tests$signif2[t.tests$signif2!="-"]<-"+"


t.tests$signif <- ""
for (i in 1:length(t.tests$p.value)){
if(t.tests$p.value[i] < 0.001){
  t.tests$signif[i] <- paste("*","*","*",sep="\n")
}else if(t.tests$p.value[i] < 0.01){
  t.tests$signif[i] <- paste("*","*",sep="\n")
}else if(t.tests$p.value[i] < 0.05){
  t.tests$signif[i] <- paste("*",sep="\n")
}else if(t.tests$p.value[i] < 0.1){
  t.tests$signif[i] <- paste("'",sep="\n")
}
}

test_results[[paste(measure,"historical_change_since_1973_within", sep="_")]] <- t.tests
hist.melt6 <- merge(hist.melt5,t.tests, by=c("CoordY.x","month_record.y"))

pdf(paste("plots.supplement/",measure,"_historical_1973.pdf",sep=""), width=8, height=6)
plot <- ggplot(hist.melt6,aes(y=difference,x=CoordY.x, col=as.factor(CoordY.x)))+
  annotate("rect", xmin = min(t.tests$lat)-1, xmax =  max(t.tests$lat)+1, ymin = 0, ymax = ymin,
           alpha = .4,fill = "steelblue")+
  facet_wrap(month_record.y~.)+
  geom_hline(yintercept=0,linewidth=1, col="red")+
  geom_boxplot(outlier.shape = NA)+
  geom_smooth(aes(y=difference,x=CoordY.x), inherit.aes = F, col="black", method="loess")+
  #geom_point(alpha=0.2,size=1)+
  #facet_wrap(month_record.x ~.)+
  theme_bw()+
  xlab("Latitude") + ylab(paste(measure,"difference since 1973"))+
  scale_y_continuous(breaks=breaks_plot)+
  scale_x_continuous(expand=c(0,0))+
  scale_color_viridis(discrete=T) + 
  theme(panel.margin.y = unit(0, "lines"),
                                         strip.background =element_rect(fill="white"),
                                         panel.grid.minor = element_blank())
  
plot <- plot + geom_text(aes(x=CoordY.x, y=le-ymin2, label=signif),inherit.aes = T,lineheight = .25)#+
#  facet_wrap(month_record.y~.)
  
plot +  geom_text(aes(y = ue+ymin2, label=signif2),inherit.aes = T)#+
#  facet_wrap(month_record.y~.)
#xlim(1900,2023)
dev.off()

