#todo: matching issues in merge? 

historical.data <- read.table(paste("Data/",measure,"_Data.csv",sep=""), sep=",", header=T)

historical.data.merge <- merge(data.frame(sample_data(data.species.rel.filter))[,c("BK","Month","CoordY")],historical.data,by="BK")

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
  scale_color_viridis(discrete=T, end=0.8)+
  scale_fill_viridis(discrete=T, end=0.8)+
  xlim(1900,2023)
dev.off()

# difference between 2023 and 1973

hist.melt5 <- hist.melt4[hist.melt4$year_record.x=="1973",]

# historical (1973) values are value.x, current values are value.y
if(measure=="Temperature"){
  # difference past - present 
  #-> positive numbers: decrease of Temp, negative numbers: increase of Temp over time
  hist.melt5$difference <- (hist.melt5$value.x) - (hist.melt5$value.y)
  breaks_plot=seq(-10,10,1)
}else{
  breaks_plot=seq(-200,+200,20)
  #-> positive numbers: decrease of Prec, negative numbers: increase of Prec over time
  hist.melt5$difference <- ((hist.melt5$value.x) - (hist.melt5$value.y))
}
hist.melt5$difference<- -hist.melt5$difference

pdf(paste("plots.supplement/",measure,"_historical_1973.pdf",sep=""), width=8, height=6)
ggplot(hist.melt5,aes(y=difference,x=CoordY.x, col=as.factor(CoordY.x)))+
  geom_boxplot()+
  geom_smooth(aes(y=difference,x=CoordY.x), inherit.aes = F, col="black", method="loess")+
  #geom_point(alpha=0.2,size=1)+
  #facet_wrap(month_record.x ~.)+
  theme_bw()+
  xlab("Latitude") + ylab(paste(measure,"difference since 1973"))+
  scale_y_continuous(breaks=breaks_plot)+
  scale_color_viridis(discrete=T)#+
  #xlim(1900,2023)
dev.off()

hist.melt.6 <- hist.melt5 %>%
  group_by(CoordY.x) %>%
  summarize(mean = mean(as.numeric(difference), na.rm=T))


test_results[[paste(measure,"historical_change_since_1973_latitude", sep="_")]] <- cor.test(hist.melt.6$mean,as.numeric(hist.melt.6$CoordY.x))

t.tests <- data.frame()
for (latitude in unique(hist.melt.6$CoordY.x)){
  sub_difference <- hist.melt5[hist.melt5$CoordY.x==latitude,"difference"]
  test <- t.test(sub_difference,mu=0,alternative="two.sided")
  t.tests <- rbind(t.tests,data.frame(lat=latitude, p.value = round(test$p.value, digits=3), est=test$statistic))
}

t.tests$signif2 <- substr(as.character(t.tests$est),1,1)
t.tests$signif2[t.tests$signif2!="-"]<-"+"

t.tests$signif <- ""
for (i in 1:length(t.tests$p.value)){
if(t.tests$p.value[i] < 0.001){
  t.tests$signif[i] <- paste("***",t.tests$signif2[i],sep="\n")
}else if(t.tests$p.value[i] < 0.01){
  t.tests$signif[i] <- paste("**",t.tests$signif2[i],sep="\n")
}else if(t.tests$p.value[i] < 0.05){
  t.tests$signif[i] <- paste("*",t.tests$signif2[i],sep="\n")
}else if(t.tests$p.value[i] < 0.1){
  t.tests$signif[i] <- paste("'",t.tests$signif2[i],sep="\n")
}
}

test_results[[paste(measure,"historical_change_since_1973_within", sep="_")]] <- t.tests

pdf(paste("plots.supplement/",measure,"_historical_1973.pdf",sep=""), width=8, height=6)
ggplot(hist.melt5,aes(y=difference,x=CoordY.x, col=as.factor(CoordY.x)))+
  geom_boxplot()+
  geom_hline(yintercept=1)+
  geom_smooth(aes(y=difference,x=CoordY.x), inherit.aes = F, col="black", method="loess")+
  #geom_point(alpha=0.2,size=1)+
  #facet_wrap(month_record.x ~.)+
  theme_bw()+
  xlab("Latitude") + ylab(paste(measure,"difference since 1973"))+
  scale_y_continuous(breaks=breaks_plot)+
  scale_color_viridis(discrete=T)+
  geom_text(data=t.tests,aes(x=lat, y = min(hist.melt5$difference)-2, label=signif),inherit.aes = F)
#xlim(1900,2023)
dev.off()

