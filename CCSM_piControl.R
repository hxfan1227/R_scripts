######
#setup
#all of these must be done for the rest to work properly
cd /home/lmkh201/Documents/R/model_data/pr/piControl/CCSM4/r1i1p1
R

.libPaths(new="/home/lmkh201/Documents/R/R_packages")


library(maps)
library(mapdata)
library(mapproj)
library(fBasics) #masks norm
library(abind)

library(ncdf)
library(akima) 
library(pracma) #masks akimaInterp, inv, kron, pascal
library(clim.pact) #masks mod, num2str, map
library(R.methodsS3)
library(R.oo) #masks getDescription, getTitle, attach, getClasses, getMethods, attach, detach, gc, load, save
library(R.utils) #masks insert, timestamp, cat, commandArgs, getOption, inherits, isOpen, parse, warnings



###PRE-INDUSTRIAL CONTROL



##calculate global climatological mean ####################################
###800-1000 data
nc.1<- retrieve.nc("pr_Amon_CCSM4_piControl_r1i1p1_080001-130012.nc",v.nam="pr",t.rng=c("800","1000")) ##get data from 800-1000AD
mm<-nc.1$mm ##vector rep(1:12), length 2400
lat<-length(nc.1$lat) ##latitude res
lon<-length(nc.1$lon) ##longitude res
yy<-length(nc.1$yy)/12 ##number of years
latval<-nc.1$lat
lonval<-nc.1$lon


nc.1<-86400*nc.1$dat ##mm/day

mean1<-array(NA,dim=c(12,lat,lon)) #monthly mean array
for (i in 1:lat)
{
for (j in 1:lon)
{
mean1[,i,j]<-tapply(nc.1[,i,j],mm,mean,na.rm=T) ##monthly mean for each grid box 
}}

rm(nc.1)


###1000-1200
nc.2<- retrieve.nc("pr_Amon_CCSM4_piControl_r1i1p1_080001-130012.nc",v.nam="pr",t.rng=c("1000","1200")) ##get data from 1000-1200AD
mm<-nc.2$mm ##vector rep(1:12), length 2400
yy<-length(nc.2$yy)/12 ##number of years
nc.2<-86400*nc.2$dat ##mm/day

mean2<-array(NA,dim=c(12,lat,lon)) #monthly mean array
for (i in 1:lat)
{
for (j in 1:lon)
{
mean2[,i,j]<-tapply(nc.2[,i,j],mm,mean,na.rm=T) 
}}

rm(nc.2)


###1200-1300
nc.3<- retrieve.nc("pr_Amon_CCSM4_piControl_r1i1p1_080001-130012.nc",v.nam="pr",t.rng=c("1200","1300 12 16")) ##get data from 1200-1301AD
mm<-nc.3$mm ##vector rep(1:12), length 1212
yy<-length(nc.3$yy)/12 ##number of years
nc.3<-86400*nc.3$dat ##mm/day

mean3<-array(NA,dim=c(12,lat,lon)) #monthly mean array
for (i in 1:lat)
{
for (j in 1:lon)
{
mean3[,i,j]<-tapply(nc.3[,i,j],mm,mean,na.rm=T) 
}}

rm(nc.3)


##Define dimensions
x <- dim.def.ncdf( "Lat", "degrees_north", latval)
y <- dim.def.ncdf( "Lon", "degrees_east", lonval)
m <- dim.def.ncdf( "Month", "months", 1:12)

##define variables
mean.08.10<-var.def.ncdf("Mean0810","mm/day",list(m,x,y),NA,longname="Monthly mean precipitation 800-1000")

mean.10.12<-var.def.ncdf("Mean1012","mm/day",list(m,x,y),NA,longname="Monthly mean precipitation 1000-1200")

mean.12.13<-var.def.ncdf("Mean1213","mm/day",list(m,x,y),NA,longname="Monthly mean precipitation 1200-1300")

##create netCDF file
mmeannc<-create.ncdf("CCSM4_pr_piC_monthly_mean.nc", list(mean.08.10,mean.10.12,mean.12.13))

##add data to netCDF file
put.var.ncdf(mmeannc,mean.08.10,mean1)
put.var.ncdf(mmeannc,mean.10.12,mean2)
put.var.ncdf(mmeannc,mean.12.13,mean3)


close.ncdf(mmeannc)

rm(mean.08.10,mean.10.12,mean.12.13,mean1,mean2,mean3,mmeannc)

###250-500 data
nc.1<- retrieve.nc("pr_Amon_CCSM4_piControl_r1i1p1_025001-050012.nc",v.nam="pr") ##get data from 250-500
mm<-nc.1$mm ##vector rep(1:12), length 2400
lat<-length(nc.1$lat) ##latitude res
lon<-length(nc.1$lon) ##longitude res
yy<-length(nc.1$yy)/12 ##number of years
latval<-nc.1$lat
lonval<-nc.1$lon

nc.1<-86400*nc.1$dat ##mm/day

mean1<-array(NA,dim=c(12,lat,lon)) #monthly mean array
for (i in 1:lat)
{
for (j in 1:lon)
{
mean1[,i,j]<-tapply(nc.1[,i,j],mm,mean,na.rm=T) ##monthly mean for each grid box 
}}

rm(nc.1)

##define variable 
mean.025.050<-var.def.ncdf("Mean025050","mm/day",list(m,x,y),1e+30,longname="Monthly mean precipitation 250-500")
##add variable
mmeannc<-open.ncdf("CCSM4_pr_piC_monthly_mean.nc", write=T)
mmeannc<-var.add.ncdf(mmeannc,mean.025.050)
close.ncdf(mmeannc)
##add data
mmeannc<-open.ncdf("CCSM4_pr_piC_monthly_mean.nc", write=T)
put.var.ncdf(mmeannc,mean.025.050,mean1)

rm(mean.025.050,mean1)


###501-799 data
nc.1<- retrieve.nc("pr_Amon_CCSM4_piControl_r1i1p1_050101-079912.nc",v.nam="pr") ##get data from 501-799 
mm<-nc.1$mm ##vector rep(1:12), length 2400
lat<-length(nc.1$lat) ##latitude res
lon<-length(nc.1$lon) ##longitude res
yy<-length(nc.1$yy)/12 ##number of years
latval<-nc.1$lat
lonval<-nc.1$lon

nc.1<-86400*nc.1$dat ##mm/day

mean1<-array(NA,dim=c(12,lat,lon)) #monthly mean array
for (i in 1:lat)
{
for (j in 1:lon)
{
mean1[,i,j]<-tapply(nc.1[,i,j],mm,mean,na.rm=T) ##monthly mean for each grid box 
}}

rm(nc.1)


##define variable 
mean.501.799<-var.def.ncdf("Mean050080","mm/day",list(m,x,y),NA,longname="Monthly mean precipitation 501-799")
##add variable
mmeannc<-open.ncdf("CCSM4_pr_piC_monthly_mean.nc", write=T)
mmeannc<-var.add.ncdf(mmeannc,mean.501.799)
close.ncdf(mmeannc)
##add data
mmeannc<-open.ncdf("CCSM4_pr_piC_monthly_mean.nc", write=T)
put.var.ncdf(mmeannc,mean.501.799,mean1)
close.ncdf(mmeannc)
rm(mean.501.799,mean1)


##get climatological mean over whole time period
nc<-open.ncdf("CCSM4_pr_piC_monthly_mean.nc")
nc1<-get.var.ncdf(nc,nc$var[[1]])
nc2<-get.var.ncdf(nc,nc$var[[2]])
nc3<-get.var.ncdf(nc,nc$var[[3]])
nc4<-get.var.ncdf(nc,nc$var[[4]])
nc5<-get.var.ncdf(nc,nc$var[[5]])

means<-abind(nc1,nc2,nc3,nc4,nc5,along=1) ###bind all means
close.ncdf(nc)

months<-gl(12,1,12*nc$nvars) #generate levels for months (#months=levels,#each level, months*number of means)

dim(means)<-c(length(means[,1,1]),lat*lon)





clmmn<-array(NA, dim=c(12,lat*lon))
for (i in 1:(lat*lon))
{
clmmn[,i]<-vaggregate(means[,i],months,mean)
}

dim(clmmn)<-c(length(clmmn[,1]),lat,lon)

rm(i,j,nc1,nc2,nc3,nc4,nc5,means,mm,months)

















mean<-array(NA, dim=c(12,lat,lon))
for (i in 1:lat)
{
for (j in 1:lon)
{
mean[,i,j]<-tapply(means[,i,j],months,mean,na.rm=T)
}}

rm(i,j,nc1,nc2,nc3,nc4,nc5,means,mm,months)




##define variable 
monmean<-var.def.ncdf("Climatology","mm/day",list(nc$dim$Month,nc$dim$Lat,nc$dim$Lon),NA,longname="Monthly mean precipitation")
##create netCDF file
mmeannc<-create.ncdf("CCSM4_pr_piC.nc", monmean)
close.ncdf(mmeannc)
##add data
mmeannc<-open.ncdf("CCSM4_pr_piC.nc", write=T)
put.var.ncdf(mmeannc,monmean,mean)
close.ncdf(mmeannc)

rm(list = ls())

9.969210e+36


###########tropics only ##########################################
#get lat indexes 
nc<-open.ncdf("pr_Amon_CCSM4_piControl_r1i1p1_025001-050012.nc",readunlim=F)
s<-match(nc$dim$lat$vals[floor(nc$dim$lat$vals)==-40],nc$dim$lat$vals) ##get index 40S
n<-match(nc$dim$lat$vals[ceiling(nc$dim$lat$vals)==40],nc$dim$lat$vals) ##index 40N
close.ncdf(nc)

##retrieve climatology
nc<-open.ncdf("CCSM4_pr_piC.nc")
clim<-get.var.ncdf(nc,"Climatology")
trop<-clim[,s:n,]
close.ncdf(nc)

##define tropics lat dimension
nc<-open.ncdf("CCSM4_pr_piC.nc")
lat<-nc$dim$Lat$vals[s:n] ##lat 40S-40N
lat<-as.vector(lat)
lat<-dim.def.ncdf("Latitude (Tropics)","degrees_north",lat)

##define variable
climtrop<-var.def.ncdf("Climatology (Tropics)","mm/day",list(nc$dim$Month,lat,nc$dim$Lon),NA,longname="Monthly mean precipitation (Tropics)")

##add variable 
nc<-open.ncdf("CCSM4_pr_piC.nc",write=T)
nc<-var.add.ncdf(nc,climtrop)
close.ncdf(nc)

##add data to variable
nc<-open.ncdf("CCSM4_pr_piC.nc",write=T)
trop<-put.var.ncdf(nc,"Climatology (Tropics)",trop)
close.ncdf(nc)

rm(list = ls())





#########RETRIEVE DATA#######################################

######define dimensions TROPICS


###define Time dimension by taking timval from each data set
nc<-open.ncdf("pr_Amon_CCSM4_piControl_r1i1p1_025001-050012.nc")
s<-match(nc$dim$lat$vals[floor(nc$dim$lat$vals)==-40],nc$dim$lat$vals) ##get index 40S
n<-match(nc$dim$lat$vals[ceiling(nc$dim$lat$vals)==40],nc$dim$lat$vals) ##index 40N
lat<-nc$dim$lat$vals[s:n] ##lat 40S-40N
lat<-as.vector(lat)
lon<-nc$dim$lon$vals ##lon
lon<-as.vector(lon)
t250<-nc$dim$time$vals ##t250
t250<-as.vector(t250)
close.ncdf(nc)

nc<-open.ncdf("pr_Amon_CCSM4_piControl_r1i1p1_050101-079912.nc")
t501<-nc$dim$time$vals ##t501
t501<-as.vector(t501)
close.ncdf(nc)

nc<-open.ncdf("pr_Amon_CCSM4_piControl_r1i1p1_080001-130012.nc")
t800<-nc$dim$time$vals ##t800
t800<-as.vector(t800)
close.ncdf(nc)

t<-abind(t250,t501,t800) #250-1300
t<-as.vector(t)

tdim <- dim.def.ncdf( "Time", "days since 0001-01-01 00:00:00", t,unlim=TRUE)
latdim<-dim.def.ncdf("Latitude (Tropics)","degrees_north",lat)
londim<-dim.def.ncdf("Longitude","degrees_east",lon)

rm(lat,lon,nc,t,t250,t501,t800)


#### define variables
anom<-var.def.ncdf("ClimAnom","mm/day",list(latdim,londim,tdim),NA,longname="Precipitation Anomaly 250-1300")

rm(latdim,londim,tdim,n,s)



###250-500 data#################################
nc.1<- retrieve.nc("pr_Amon_CCSM4_piControl_r1i1p1_025001-050012.nc",v.nam="pr",y.rng=c(-40,40))
lat<-length(nc.1$lat) ##latitude res
lon<-length(nc.1$lon) ##longitude res
yy<-length(nc.1$tim) ##number of timesteps

nc.1<-86400*nc.1$dat ##mm/day [time,lat,lon]

###change dimension order to [lat,lon,time]
res<-lat*lon ##resolution, i.e. #gridboxes
dim(nc.1)<-c(yy,(res)) ##make matrix [time,loc]
nc.1<-t(nc.1) ##[loc,time]
dim(nc.1)<-c(lat,lon,yy) 

##retrieve climatology
nc<-open.ncdf("CCSM4_pr_piC.nc")
clim<-get.var.ncdf(nc,"Climatology (Tropics)")
close.ncdf(nc)

##make clim array for data set
res<-lat*lon ##resolution, i.e. #gridboxes
dim(clim)<-c(12,(res)) ##make matrix [month,loc]
clim<-t(clim) ##[loc,month]
clim<-replicate(yy/12,clim) ##array[55296,12,250] =[loc,month,yrs]
dim(clim)<-c(lat,lon,yy) #monthly mean at each grid box over entire period

##remove climatological mean
nc.1<-nc.1-clim ##climatological anomaly 
rm(clim,lat,lon,nc,res)


##create anomaly ncdf
nc<-create.ncdf("CCSM4_pr_piC_anom_tropics.nc",anom)
close.ncdf(nc)

##add data to variable
nc<-open.ncdf("CCSM4_pr_piC_anom_tropics.nc", write=T)
put.var.ncdf( nc, nc$var$ClimAnom, nc.1, start=c(1,1,1), count=c(-1,-1,yy) )
close.ncdf(nc)

start<-yy+1 ###time starting point for next data set

rm(anom)



###501-799 data###########################
nc.1<- retrieve.nc("pr_Amon_CCSM4_piControl_r1i1p1_050101-079912.nc",v.nam="pr",y.rng=c(-40,40))
lat<-length(nc.1$lat) ##latitude res
lon<-length(nc.1$lon) ##longitude res
yy<-length(nc.1$tim) ##number of timesteps

nc.1<-86400*nc.1$dat ##mm/day [time,lat,lon]

###change dimension order to [lat,lon,time]
res<-lat*lon ##resolution, i.e. #gridboxes
dim(nc.1)<-c(yy,(res)) ##make matrix [time,loc]
nc.1<-t(nc.1) ##[loc,time]
dim(nc.1)<-c(lat,lon,yy) 

##retrieve climatology
nc<-open.ncdf("CCSM4_pr_piC.nc")
clim<-get.var.ncdf(nc,"Climatology (Tropics)")
close.ncdf(nc)

##make clim array for data set
res<-lat*lon ##resolution, i.e. #gridboxes
dim(clim)<-c(12,(res)) ##make matrix [month,loc]
clim<-t(clim) ##[loc,month]
clim<-replicate(yy/12,clim) ##array[55296,12,250] =[loc,month,yrs]
dim(clim)<-c(lat,lon,yy) #monthly mean at each grid box over entire period

##remove climatological mean
nc.1<-nc.1-clim ##climatological anomaly 
rm(clim,lat,lon,nc,res)

##add data to variable
nc<-open.ncdf("CCSM4_pr_piC_anom_tropics.nc", write=T)
put.var.ncdf( nc, nc$var$ClimAnom, nc.1, start=c(1,1,start), count=c(-1,-1,yy) )
close.ncdf(nc)

start<-start+yy ###time starting point for next data set


###800-1300 data############################
nc.1<- retrieve.nc("pr_Amon_CCSM4_piControl_r1i1p1_080001-130012.nc",v.nam="pr",y.rng=c(-40,40))
lat<-length(nc.1$lat) ##latitude res
lon<-length(nc.1$lon) ##longitude res
yy<-length(nc.1$tim) ##number of timesteps

nc.1<-86400*nc.1$dat ##mm/day [time,lat,lon]

###change dimension order to [lat,lon,time]
res<-lat*lon ##resolution, i.e. #gridboxes
dim(nc.1)<-c(yy,(res)) ##make matrix [time,loc]
nc.1<-t(nc.1) ##[loc,time]
dim(nc.1)<-c(lat,lon,yy) 

##retrieve climatology
nc<-open.ncdf("CCSM4_pr_piC.nc")
clim<-get.var.ncdf(nc,"Climatology (Tropics)")
close.ncdf(nc)

##make clim array for data set
res<-lat*lon ##resolution, i.e. #gridboxes
dim(clim)<-c(12,(res)) ##make matrix [month,loc]
clim<-t(clim) ##[loc,month]
clim<-replicate(yy/12,clim) ##array[55296,12,250] =[loc,month,yrs]
dim(clim)<-c(lat,lon,yy) #monthly mean at each grid box over entire period

##remove climatological mean
nc.1<-nc.1-clim ##climatological anomaly 
rm(clim,lat,lon,nc,res)

##add data to variable
nc<-open.ncdf("CCSM4_pr_piC_anom_tropics.nc", write=T)
put.var.ncdf( nc, nc$var$ClimAnom, nc.1, start=c(1,1,start), count=c(-1,-1,yy) )
close.ncdf(nc)



#############################################################################
###################EOF ANALYSIS####################

nc<-open.ncdf("CCSM4_pr_piC_anom_tropics.nc", write=T)
nc.1<-get.var.ncdf(nc,"annual_mean",start=c(1,1,1),count=c(-1,-1,1051)) 














#switch lon
east<-round(lon/2)
west<-lon-east
nc.e<-nc.1[,,1:east]
nc.w<-nc.1[,,west:lon]
nc1<-abind(nc.w,nc.e,along=3)

rm(nc.1,nc.e,nc.w,east,west)


## EOF
nc<-wrap(nc,map=list(1,NA)) #make a matrix [t,spat]
nc<-detrend(t(nc)) #detrend [spat,t]
nc<-t(nc) #matrix [t,spat]

#for p>>n EOF and PC
L<-nc%*%t(nc) #small covariance matrix [t,t]
eig<-eigen(L) #eig$vectors and eig$values
lam<-diag(eig$values, nrow = length(eig$values)) #diagonal matrix Lambda (=values); equal to those of R=t(nc)%*%nc
B<-eig$vectors #eigenvectors B
D<-t(nc)%*%B #proportional EOFs
eof<-sweep(D, 2, sqrt(colSums(D^2)), FUN="/") #normalise EOFs

var.pic<-diag(lam)/tr(lam) #variance explained by each EOF
var.pic[1:10]*100 #variance % for first 10 EOFs
# [1] 25.0434527 11.1934515  5.5163690  2.4622298  1.7217062  1.4426120
# [7]  1.3081065  1.2675816  1.0269953  0.9798475

rm(L,eig,B,D,lam)











##EOF1
eof.pic1<-eof[,1] #1st EOF
pc.pic1<-nc%*%eof.pic1 #PC1
corr.pic1<-cor(pc.pic1,nc,method="spearman") #correlation1
dim(corr.pic1)<-c(192,288) 
corr.pic1<-t(corr.pic1) #get map dimensions [lon,lat]
varexp.pic1<-(corr.pic1^2)*100


##EOF2
eof.2k<-eof[,2] #1st EOF
pc.2k<-nc%*%eof.2k #PC1
corr.2k<-cor(pc.2k,nc,method="spearman") #correlation1
dim(corr.2k)<-c(192,288) 
corr.2k<-t(corr.2k) #get map dimensions [lon,lat]
varexp.2k<-(corr.2k^2)*100




############MAPS
##change EOF dimensions
dim(eof.1k)<-c(192,288) #1st EOF, in [lat,lon] last millennium
eof.1k<-t(eof.1k) #change to [lon,lat]

dim(eof.2k)<-c(192,288) #1st EOF, in [lat,lon] last millennium
eof.2k<-t(eof.2k) #change to [lon,lat]


#map settings
ccsm<-open.ncdf("pr_Amon_CCSM4_midHolocene_r1i1p1_100001-130012.nc") #extract lat, lon vectors
lon<-get.var.ncdf(ccsm,"lon") 
lat<-get.var.ncdf(ccsm,"lat")

data(world2HiresMapEnv) #map contours
colour<-colorRampPalette(c("red","white","blue"), space="rgb")
colour2<-colorRampPalette(c("white","yellow","orange","red"), space="rgb")



##past millennium
#global EOF1
png(filename="EOF_1_1k.png", width=930,height=560) ##EOF1
filled.contour(lon, lat, eof.1k*100, nlevels=40, zlim=c(-4,4), color.palette=colour, plot.title=title(main="1st EOF Last Millennium"),plot.axes={maps::map(database="world2Hires",interior=T,add=T)})
dev.off()

png(filename="Correlation_1_1k.png", width=930,height=560) ## correlation1
filled.contour(lon, lat, corr.1k, nlevels=40, zlim=c(-1,1), color.palette=colour, plot.title=title(main="1st EOF Correlation Last Millennium"),plot.axes={maps::map(database="world2Hires",interior=T,add=T)})
dev.off()

png(filename="Variance_1_1k.png", width=930,height=560) ## variance explained1
filled.contour(lon, lat, varexp.1k, nlevels=40, zlim=c(0,100), color.palette=colour2, plot.title=title(main="Variance Explained by EOF1 Last Millennium"),plot.axes={maps::map(database="world2Hires",interior=T,add=T)})
dev.off()

png("PC1_Var_1k.png", width=930,height=560) # PC1 and variance
plot.new()
par(mfrow=c(2,1)) #2 figures in 2 rows, 1 column
plot.ts(pc.1k,xlab="Year AD",ylab="",main="PC1 Last Millennium")
plot(var.1k[1:10]*100, main="Variances",xlab="EOF",ylab="Variance (%)",type="l")
dev.off()

#global EOF2
png(filename="EOF_2_1k.png", width=930,height=560) ##EOF1
filled.contour(lon, lat, eof.2k*100, nlevels=40, zlim=c(-4,4), color.palette=colour, plot.title=title(main="2nd EOF Last Millennium"),plot.axes={maps::map(database="world2Hires",interior=T,add=T)})
dev.off()

png(filename="Correlation_2_1k.png", width=930,height=560) ## correlation1
filled.contour(lon, lat, corr.2k, nlevels=40, zlim=c(-1,1), color.palette=colour, plot.title=title(main="2nd EOF Correlation Last Millennium"),plot.axes={maps::map(database="world2Hires",interior=T,add=T)})
dev.off()

png(filename="Variance_2_1k.png", width=930,height=560) ## variance explained1
filled.contour(lon, lat, varexp.2k, nlevels=40, zlim=c(0,100), color.palette=colour2, plot.title=title(main="Variance Explained by EOF2 Last Millennium"),plot.axes={maps::map(database="world2Hires",interior=T,add=T)})948
dev.off()

png("PC2_Var_1k.png", width=930,height=560) # PC1 and variance
plot.new()
par(mfrow=c(2,1)) #2 figures in 2 rows, 1 column
plot.ts(pc.2k,xlab="Year AD",ylab="",main="PC2 Last Millennium")
plot(var.1k[1:10]*100, main="Variances",xlab="EOF",ylab="Variance (%)",type="l")
dev.off()


png("PC1_MCA_LIA_1k.png", width=930,height=560) # PC1 and variance
plot.new()
par(mfrow=c(2,1)) #2 figures in 2 rows, 1 column
plot.ts(pc.1k[150:450],xlim=c(0,450),ylim=c(-150,150),xlab="Year since AD 1000",ylab="",main="PC1 MCA")
plot.ts(pc.1k[550:1000],xlim=c(0,450),ylim=c(-150,150),xlab="Year since AD 1400",ylab="",main="PC1 LIA")
dev.off()

























