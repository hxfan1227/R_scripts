##calculating tropical climatological mean SECTION 1 ####################################


###section 1
data<-	"pr_Amon_CCSM4_midHolocene_r1i1p1_100001-130012.nc" ##name of dataset to be used
###section 2
varname<- "climatology2"
longname<-"Climatology"


var<-	"pr"		##name of variable extracted
conv<-86400		##unit conversion (1 if NA)
	#precipitation mm/day: 86400

#Break dataset into spatial (lon) sections
set<- 3			##number of sections
n<- 2			##section being used

x.min<- 97		##start lon
x.max<-	192		##end lon
	#for lon 288: 0-96,97-192, 193-288
totaltime<-3612		##full dataset time length (all time sets)

ncname<-"CCSM4_pr_mh_clim2.nc"
units<- "mm/day"



###1. Get precipitation data (tropics) ##############

nc<-open.ncdf(data)
#define lat indexes for tropics (40S-40N)
lat.min<-match(nc$dim$lat$vals[round(nc$dim$lat$vals)==-40],nc$dim$lat$vals) 
lat.max<-match(nc$dim$lat$vals[round(nc$dim$lat$vals)==40],nc$dim$lat$vals) 

lon<-nc$dim$lon$len/set
lat<-lat.max-lat.min+1
time<-nc$dim$time$len
yy<-time/12


#get variable data and dimensions for the tropics
var<-get.var.ncdf(nc, "pr", start=c(1+((n-1)*lon),lat.min,1), count=c(lon,lat,-1))

close.ncdf(nc)

#convert to mm/day
var<-var*conv

#create levels for each month (jan=1,feb=2,...)
mm<-gl(12, 1, time)
#create an array to hold monthly means
clim<-array(NA,dim=c(lon,lat,12)) 
#change variable dims and loop to calculate means for each month

for (i in 1:lon)
{
for (j in 1:lat)
{
clim[i,j,]<-vaggregate(var[i,j,],mm,mean,na.rm=T) 
}}


###2. Create ncetCDF file ################################
##Define dimensions
x <- dim.def.ncdf( "lon2", "degrees_east", nc$dim$lon$vals[x.min:x.max])
y <- dim.def.ncdf( "lat", "degrees_north", nc$dim$lat$vals[lat.min:lat.max])
t <- dim.def.ncdf( "time", "timesteps", 1:totaltime,unlim=TRUE)

##define variables
climvar<-var.def.ncdf(varname,units,list(x,y,t),NA,longname=longname)

##create netCDF file
nc<-create.ncdf(ncname, climvar)

##add data to netCDF file
put.var.ncdf(nc,varname,clim,start=c(1,1,1),count=c(-1,-1,12))

close.ncdf(nc)

##SUBSEQUENT TIMES: define variable
nc<-open.ncdf(ncname,write=T)
climvar<-var.def.ncdf(varname,units,list(nc$dim$lon2,nc$dim$lat,nc$dim$time),NA,longname=longname)
nc<-var.add.ncdf(nc,climvar)

close.ncdf(nc)

##add data to netCDF file
nc<-open.ncdf(ncname,write=T)
put.var.ncdf(nc,varname,clim,start=c(1,1,1),count=c(-1,-1,12))

close.ncdf(nc)


