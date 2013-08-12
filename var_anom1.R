


ncname<-"CCSM4_pr_piC_clim1.nc"	##name of netCDF file
units<- "mm/day"		##units of variable
varname<-"climatology1"		##name of variable created
n<-2				##number of timesets

##get climatological mean over whole time period
nc<-open.ncdf(ncname,write=T)
#create 4D array to hold clim means for each timeset
means<-array(NA,dim=c(nc$dim$lon1$len, nc$dim$lat$len,12,n))
#create array to hold final climatology
clim<-array(NA,dim=c(nc$dim$lon1$len, nc$dim$lat$len,12))
#generate levels for each month in each timeset
gg<-gl(12, 1, 12*n)
#loop to put data in means array, then take mean over all timesets (4th dimension)
for (i in 1:n)
{
means[,,,i]<-get.var.ncdf(nc,nc$var[[i]],count=c(-1,-1,12))
}

for (j in 1:nc$dim$lon1$len)
{
for (k in 1:nc$dim$lat$len)
{
for (l in 1:12)
clim[j,k,]<-vaggregate(means[j,k,,],gg,mean)
}}

close.ncdf(nc)


##define variable 
monmean<-var.def.ncdf(varname,units,list(nc$dim$lon1,nc$dim$lat,nc$dim$time),NA,longname="Climatology")

##add variable and data
nc<-open.ncdf(ncname, write=T)
nc<-var.add.ncdf(nc,monmean)
close.ncdf(nc)

nc<-open.ncdf(ncname, write=T)
put.var.ncdf(nc,varname,clim,count=c(-1,-1,12))
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

























