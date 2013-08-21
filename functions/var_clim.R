##calculating tropical climatological mean SECTION 1 ####################################
#This calculates climatologies for the first longitudinal section. Repeat this for all time sections. 


function(vb,ncname,data,varname,units="mm/day"


vb<-	"pr"		##name of variable extracted
conv<-86400		##unit conversion (1 if NA)
	#precipitation mm/day: 86400

#Break dataset into spatial (lon) sections
set<- 3			##number of sections
n<- 1			##section being used

x.min<- 0		##start lon
x.max<-	96		##end lon
	#for lon 288: 0-96,97-192, 193-288
totaltime<-3612		##full dataset time length (all time sets)

###section 2
ncname<-"CCSM4_pr_mh_clim1.nc"

data<-	"pr_Amon_CCSM4_midHolocene_r1i1p1_100001-130012.nc" ##name of dataset to be used
varname<- "climatology1"
units<- "mm/day"
longname<-"Climatology"

### parameters ########
#vb		Character string; climatological variable to analyse (as in raw dataset filenames)
#model		Character string; model that data is taken from (as in raw dataset filenames)
#period		Character string; PMIP3 timeslice (pic,h,1k,mh) or empirical data range (yy-yy)
#dx		Numerical value; dataset LON slice
#dy		Numerical value; dataset LAT slice



varclim<-function(vb,model,period,dx,dy,conv=86400)
{
	##Get filenames of all data timesets within a model timeslice
	filename<-paste(vb,"Amon",model,"*",sep="_") #creates search string "[var]_Amon_[model]_*"
	files<-Sys.glob(filename) #vector with names of all raw datafiles
	data<-files[1]

	###1. Create new netCDF file ##############
	nc<-open.ncdf(data,write=T)
		lon<-nc$dim$lon$len/dx
		lat<-nc$dim$lat$len/dy
		time<-nc$dim$time$len
		yy<-time/12

		#Define dimensions
		#Take lon,lat values for slice [dx,dy]
		londim <- dim.def.ncdf( "lon", "degrees_east", nc$dim$lon$vals[(1+((x-1)*lon)):(x*lon)]) 
		latdim <- dim.def.ncdf( "lat", "degrees_north", nc$dim$lat$vals[(1+((y-1)*lat)):(y*lat)])
		#Time is unlim so doesn't matter if there are multiple timesets
		timedim <- dim.def.ncdf( "time", "timesteps", 1:time,unlim=TRUE)
	close.ncdf(nc)

	##Define new variables
		#create variable names
	climvar<-paste("climatology", x,y, sep = ".") #gives "climatology.x.y"
	anomvar<-paste("anomaly", x,y, sep = ".") #gives "anomaly.x.y"
	annual<-paste("annual_mean", x,y, sep = ".") #gives "annual_mean.x.y"
	winter<-paste("winter_mean", x,y, sep = ".") #gives "winter_mean.x.y"
	summer<-paste("summer_mean", x,y, sep = ".") #gives "summer_mean.x.y"
		#define vars
	climvar<-var.def.ncdf(climvar,units,list(londim,latdim,timedim),NA,longname=climvar)
	anomvar<-var.def.ncdf(anomvar,units,list(londim,latdim,timedim),NA,longname=anomvar)
	annual<-var.def.ncdf(annual,units,list(londim,latdim,timedim),NA,longname=annual)
	winter<-var.def.ncdf(winter,units,list(londim,latdim,timedim),NA,longname=winter)
	summer<-var.def.ncdf(summer,units,list(londim,latdim,timedim),NA,longname=summer)

	##create netCDF file
	ncname<-paste(model,vb,period,".nc",sep=".") #gives "model.variable.timeslice.nc"
	nc<-create.ncdf(ncname, climvar,anomvar,annual,winter,summer)
	close.ncdf(nc)

	#### 2. Get variable data ######################
	#get all timesets filenames 
	filename<-paste(vb,"Amon",model,"*",sep="_") #creates search string "[var]_Amon_[model]_*"
	files<-Sys.glob(filename) #vector with names of all raw datafiles

	if (length(files)==1)
	{
		nc<-open.ncdf(data,write=T)
			var<-get.var.ncdf(nc, vb, start=c(1+((x-1)*lon),1+((y-1)*lat),1), count=c(lon,lat,-1))
			if (vb=="pr")
			{
				#convert to mm/day
				var<-var*86400
			}
			else
			{
				var<-var #can be changed if it's necessary for SST to convert...?
			}
		close.ncdf(nc)

		####### 3. Calculate climatology ############
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
			}
		}
		##add data to new netCDF file
		nc<-open.ncdf(ncname,write=T)
			put.var.ncdf(nc,climvar,clim,start=c(1,1,1),count=c(-1,-1,12))
		close.ncdf(nc)
	}
	else
	{
		for (f in 1:length(files))
		{
			nc<-open.ncdf(files[f],write=T)
				var<-get.var.ncdf(nc, vb, start=c(1+((x-1)*lon),1+((y-1)*lat),1), count=c(lon,lat,-1))
				#convert to mm/day
				var<-var*conv
			close.ncdf(nc)

			####### 3. Calculate climatology ############
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
				}
			}
			##add data to new netCDF file
			nc<-open.ncdf(ncname,write=T)
				put.var.ncdf(nc,climvar,clim,start=c(1,1,(1+((f-1)*12),count=c(-1,-1,12))
			close.ncdf(nc)
		}
	}


####this script will calculate the final climatology over the entire timeseries. Repeat this for each lon section!!

ncname<-"CCSM4_pr_piC_clim3.nc"	##name of netCDF file with clims
units<- "mm/day"		##units of variable
varclim<-"climatology3"		##name of climatology variable created
nc<-open.ncdf(ncname)
long <-nc$dim$lon3		##lon dimension (change lon1/2/3)
longlen<- nc$dim$lon3$len	##lon dimension length (lon1/2/3)
close.ncdf(nc)
nt<-3				##number of timesets in clim netCDF


###1. CLIMATOLOGICAL MEAN #########
nc<-open.ncdf(ncname,write=T)
#create 4D array to hold clim means for each timeset
means<-array(NA,dim=c(longlen, nc$dim$lat$len,12,nt))
#create array to hold final climatology
clim<-array(NA,dim=c(longlen, nc$dim$lat$len,12))
#generate levels for each month in each timeset
gg<-gl(12, 1, 12*nt)
#loop to put data in means array, then take mean over all timesets (4th dimension)
for (i in 1:nt)
{
means[,,,i]<-get.var.ncdf(nc,nc$var[[i]],count=c(-1,-1,12))
}

for (j in 1:longlen)
{
for (k in 1:nc$dim$lat$len)
{
for (l in 1:12)
clim[j,k,]<-vaggregate(means[j,k,,],gg,mean)
}}

close.ncdf(nc)

#######2. Add to netCDF file ##################
##define variable 
monmean<-var.def.ncdf(varclim,units,list(long,nc$dim$lat,nc$dim$time),NA,longname="Climatology")

##add variable and data
nc<-open.ncdf(ncname, write=T)
nc<-var.add.ncdf(nc,monmean)

close.ncdf(nc)
nc<-open.ncdf(ncname, write=T)
put.var.ncdf(nc,varclim,clim,count=c(-1,-1,12))

close.ncdf(nc)


}


















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
var<-get.var.ncdf(nc, vb, start=c(1+((n-1)*lon),lat.min,1), count=c(lon,lat,-1))

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
##repeat for all timesets if applicable

##FIRST TIME: Define dimensions
nc<-open.ncdf(data)
x <- dim.def.ncdf( "lon1", "degrees_east", nc$dim$lon$vals[x.min:x.max])
y <- dim.def.ncdf( "lat", "degrees_north", nc$dim$lat$vals[lat.min:lat.max])
t <- dim.def.ncdf( "time", "timesteps", 1:totaltime,unlim=TRUE)
close.ncdf(nc)

##define variables
climvar<-var.def.ncdf(varname,units,list(x,y,t),NA,longname=longname)

##create netCDF file
nc<-create.ncdf(ncname, climvar)

##add data to netCDF file
put.var.ncdf(nc,varname,clim,start=c(1,1,1),count=c(-1,-1,12))

close.ncdf(nc)



##SUBSEQUENT TIMES: define variable
nc<-open.ncdf(ncname,write=T)
climvar<-var.def.ncdf(varname,units,list(nc$dim$lon1,nc$dim$lat,nc$dim$time),NA,longname=longname)
nc<-var.add.ncdf(nc,climvar)

close.ncdf(nc)

##add data to netCDF file
nc<-open.ncdf(ncname,write=T)
put.var.ncdf(nc,varname,clim,start=c(1,1,1),count=c(-1,-1,12))

close.ncdf(nc)


