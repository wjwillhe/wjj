  

library(ncdf4)
library(raster)
library(rgdal)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(sp)
library(stringr)
library(magrittr)
library(maps)
library(geosphere)
library(plyr)
library(abind)
library(fields)
library(factoextra)
library(lattice)
library(trend)
library(clValid)
library(cowplot)
library(lubridate)
library(VIM)

library(NbClust)

######处理地图shp文件
world<-readOGR('E:/R/worldmap/world.shp')
world<-fortify(world)
world_lon<-world[,1]


for( i in 1:586238)
{if(world_lon[i] < 0)
{world_lon[i] <- (360-abs(world_lon[i]))}
  
}

world[,1] <- world_lon

world1 <- world
world1[which(world1$long<80), ]<-NA
world1[which(world1$long>260), ]<-NA
world1[which(world1$lat<(-10)), ]<-NA
world1<-na.omit(world1)





world2<-world
world2[which(world2$long>0), ]<-NA
world2<-na.omit(world2)


#####
data = read.table("e:/data/NCEP_TRACK/mg1955-2014/R-kmeans/mg1955-2014.txt",header = F,sep="\t")#读取文本数据
data = as.matrix(data)
dd<-nc_open("E:/data/track_data/TRACK_ff_trs_pos_ERA5_1950_202207_NH_CYC.nc")

lat<-ncvar_get(dd,"traj_lat")
lon<-ncvar_get(dd,"traj_lon")
vor<-ncvar_get(dd,"traj_rel_vor_pos")
time<-ncvar_get(dd,"traj_time")
sid<-ncvar_get(dd,"parcel")
nc_close(dd)


 meg_lat <- matrix(nr=8,nc=1159)
 meg_lon <- matrix(nr=8,nc=1159)
 meg_vor <- matrix(nr=8,nc=1159)
 for(i in 1:1159){
   meg_lon[,i]<-lon[1:8,data[i]]
   meg_lat[,i]<-lat[1:8,data[i]]
   meg_vor[,i]<-vor[1:8,data[i]]
 }
 


######计算终点
######计算起始的位置

M1_lon<-meg_lon[1,]
M1_lat<-meg_lat[1,]

######计算终点879
n <- sapply(1:1159,function(i){
  M2_1 =sum(complete.cases(meg_lat[ ,i])) })#############最后一点

N1_lon <- sapply(1:1159,function(i){
  M2_1 =meg_lon[n[i],i]})

N1_lat <- sapply(1:1159,function(i){
  M2_1 =meg_lat[n[i],i]})


c_lat<- N1_lat- M1_lat
c_lon<- N1_lon- M1_lon

######Calculate average longitude and latitude
meanlat <- sapply(1:1159,function(i){
  M2_1 =mean(meg_lat[ ,i],na.rm = TRUE) })

meanlon <- sapply(1:1159,function(i){
  M2_1 =mean(meg_lon[ ,i],na.rm = TRUE) })


#####
data_clean <- as.matrix(data.frame(lon= M1_lon, lat = M1_lat,sid=data))%>%na.omit()

#####Calculate the three covariance without considering the weight
M2_lon <- sapply(1:1159,function(i){
  M2_1 = sum((meg_lon[ ,i] - meanlon[i])**2,na.rm = TRUE)/n[i]})


M2_lat <- sapply(1:1159,function(i){
  M2_1 = sum((meg_lat[ ,i] - meanlat[i])**2,na.rm = TRUE)/n[i]})

M2_r <- sapply(1:1159,function(i){
  M2_1 = sum((meg_lon[ ,i] - meanlon[i])*(meg_lat[ ,i] - meanlat[i]),na.rm = TRUE)/n[i]})
######################################


#####
result <- as.matrix(data.frame( c_lat = c_lat,M2_lat = M2_lat , M2_r= M2_r)) %>% na.omit()

####
result_scale <- scale(result) 

#head(result_scale)


#####Calculate distance coefficient
dunn.distance <- dist(result_scale,method = "euclidean")

cluster.result <- sapply(1:4,function(i){
  result <- (kmeans(x = result_scale, centers = i+1, iter.max = 9))$cluster})

#####
dunn.result <- sapply(1:4,function(i){
  result <- dunn(dunn.distance, cluster.result[,i])})

###Calculate the optimal number of clusters
#windows（）
set.seed(1234)
devAskNewPage(ask=TRUE)
nc <- NbClust(result_scale,min.nc=2, max.nc=9,method="kmeans")

###############kmeans
km_result <- kmeans(result_scale,3, nstart = 25)
####
data.result <- cbind(data_clean, cluster = km_result$cluster)


for(i in 1:1159)
{clus1 <- subset(data.result,as.factor(data.result[,4])==1)
clus2 <- subset(data.result,as.factor(data.result[,4])==2)
clus3 <- subset(data.result,as.factor(data.result[,4])==3)

}


write(clus1[ ,3] ,file="e:/data/NCEP_TRACK/mg1955-2014/R-kmeans/mgc_clus1.txt",nc=1)
write(clus2[ ,3] ,file="e:/data/NCEP_TRACK/mg1955-2014/R-kmeans/mgc_clus2.txt",nc=1)
write(clus3[ ,3] ,file="e:/data/NCEP_TRACK/mg1955-2014/R-kmeans/mgc_clus3.txt",nc=1)

