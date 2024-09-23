## GDMs 
library(gdm)

data.species.rel.filter2 <- merge_samples(data.species.rel.filter, "BK", mean)
testData1<- data.frame((otu_table(data.species.rel.filter2)))

testData1 <- testData1[rownames(testData1)!="0",]
testData1a <- testData1
testData1a$site <- rownames(testData1)
#testData1a$Long <- meta.data$Long
#testData1a$Lat <- meta.data$Lat

meta.data <- data.frame(site=data.frame(row.names(sample_data(data.species.rel.filter2)))[,1],
                                 lat=data.frame(sample_data(data.species.rel.filter2)[,c("CoordY")])[,1],
                                 long=data.frame(sample_data(data.species.rel.filter2)[,c("CoordX")])[,1],
                                 temp=data.frame(sample_data(data.species.rel.filter2)[,c("Temp")])[,1],
                                 prec=data.frame(sample_data(data.species.rel.filter2)[,c("Prec")])[,1])

meta.data <- meta.data[meta.data$site!="0",]

meta.data$site == testData1a$site

exFormat1a <- formatsitepair(testData1a, 
                             bioFormat = 1, 
                             siteColumn="site", 
                             XColumn="long", 
                             YColumn="lat", 
                             predData= meta.data)

model <- gdm(exFormat1a, geo=T)
summary(model)
plot(model, plot.layout=c(3,2))

total = 48.162
t_r = 1.815
t_g = 0.932
t_p = 0.448

t_t = t_r + t_g + t_p

# partial variances 
t_t = 48.162/t_t
t_r = t_r * t_t
t_g = t_g * t_t
t_p = t_p * t_t

t_p # prec
t_g # geo
t_r # temp

# check total 
t_p + t_g + t_r

# diversity
testData1<- data.frame((otu_table(data.species.rel.filter2)))

div = (diversity(testData1))
div2 = data.frame(site=names(div) , esn=div)

div3 = merge(div2, meta.data, by="site")
div3 = div3[div3$lat>30,]
summary(step(lm(esn~temp+long+prec+lat+I(lat^2), data=div3)))

plot(div3$esn~div3$lat)
