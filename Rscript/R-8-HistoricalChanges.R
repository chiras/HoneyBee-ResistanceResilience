#todo: matching issues in merge? 

measure = "Precipitation"
historical.data <- read.table(paste("Data/",measure,"_Data.csv",sep=""), sep=",", header=T)

head(historical.data)

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

hist.melt4$difference <- hist.melt4$value.x - hist.melt4$value.y
hist.melt4$year_record <- as.numeric(hist.melt4$year_record.x)

pdf(paste("plots.supplement/",measure,"_historical_change.pdf",sep=""), width=8, height=6)
ggplot(hist.melt4,aes(y=difference,x=year_record, col=as.factor(CoordY.y)))+
  geom_smooth(method="loess",n=100)+
  #  geom_point(alpha=0.2,size=1)+
  #facet_wrap(month_record.x ~.)+
  theme_bw()+
  scale_color_viridis(discrete=T)+
  xlim(1900,2023)
dev.off()
