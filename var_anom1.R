#### Calculate anomalies ###################
## this script will calculate the climate anomalies. Repeat this for all time sets.

ncname<-"CCSM4_pr_piC_clim1.nc"	##name of netCDF file with clims
units<- "mm/day"		##units of variable
varclim<-"climatology1"		##name of climatology variable created
varanom<-"anomaly250_500"		##name of anomaly variable created
varanomlong<-"Precipitation Anomaly250_500"	##anom variable long name

data<-"pr_Amon_CCSM4_piControl_r1i1p1_025001-050012.nc"	##dataset name
vb<-"pr"			##variable name
conv<-86400		##unit conversion (1 if NA)
	#precipitation mm/day: 86400

#Break dataset into spatial (lon) sections
set<- 3			##number of sections
n<- 1			##section being used

x.min<- 0		##start lon
x.max<-	96		##end lon
	#for lon 288: 0-96,97-192, 193-288




#########1. ANOMALIES ###################
nc<-open.ncdf(data)
#define lat indexes for tropics (40S-40N)
lat.min<-match(nc$dim$lat$vals[round(nc$dim$lat$vals)==-40],nc$dim$lat$vals) 
lat.max<-match(nc$dim$lat$vals[round(nc$dim$lat$vals)==40],nc$dim$lat$vals) 

lon<-nc$dim$lon$len/set
lat<-lat.max-lat.min+1
time<-nc$dim$time$len
yy<-time/12

#get variable data and dimensions for the tropics
var<-get.var.ncdf(nc, vb, start=c(1+((n-1)*lon),lat.min,1), count=c(lon,lat,-1))
close.ncdf(nc)

#convert to mm/day
var<-var*conv

##retrieve climatology
nc<-open.ncdf(ncname)
clim<-get.var.ncdf(nc,varclim,count=c(-1,-1,12))
close.ncdf(nc)
#replicate for each year in dataset
clim<-replicate(yy,abind(clim,along=3))
#match variable dimensions
(clim)<-dim(var)

##remove climatological mean from data to get anomaly
var<-var-clim

#####2. Add to netCDF file #####################
#define netCDF variable
nc<-open.ncdf(ncname,write=T)
anom<-var.def.ncdf(varanom,units,list(nc$dim$lon1,nc$dim$lat,nc$dim$time),NA,longname=varanomlong)
close.ncdf(nc)
#add variable to netCDF
nc<-open.ncdf(ncname,write=T)
nc<-var.add.ncdf(nc,anom)
close.ncdf(nc)
#add data to variable
nc<-open.ncdf(ncname,write=T)
put.var.ncdf(nc,varanom,var, start=c(1,1,1), count=c(-1,-1,time))
close.ncdf(nc)

