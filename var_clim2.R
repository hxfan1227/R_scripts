##calculating tropical climatological mean ####################################


###section 1
data<-	"pr_Amon_CCSM4_piControl_r1i1p1_025001-050012.nc" ##name of dataset to be used
var<-	"pr"						##name of variable extracted
conv<-86400						##unit conversion (1 if NA)
	#precipitation mm/day: 86400

#Break dataset into spatial (lon) sections
set<- 3							##number of sections
n<- 1							##section being used

x.min<- 97						##start lon
x.max<-	192						##end lon
	#for lon 288: 0-96,97-192, 193-288



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
#convert to mm/day
var<-var*conv

#create levels for each month (jan=1,feb=2,...)
mm<-gl(12, 1, time)
#create an array to hold monthly means
clim<-array(NA,dim=c(lon,lat,12)) 
#change variable dims and loop to calculate means for each month
dim(var)<-c(lon,lat,time)

for (i in 1:lon)
{
for (j in 1:lat)
{
clim[i,j,]<-vaggregate(var[i,j,],mm,mean,na.rm=T) 
}}





