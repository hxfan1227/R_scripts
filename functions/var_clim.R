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



varclim<-function(vb,model,period,dx,dy)
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

	if (length(files)==1) #if there is only 1 raw data timeset
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
	else #if there are more than 1 raw data timesets that need to be combined to make a final climatology
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
				#add each timeset's clim onto time dimension
				put.var.ncdf(nc,climvar,clim,start=c(1,1,(1+((f-1)*12),count=c(-1,-1,12)) 
			close.ncdf(nc)
		}
		##compute final climatology
		nc<-open.ncdf(ncname,write=T)
			#create 4D array to hold clim means for each timeset
			means<-array(NA,dim=c(lon,lat,12,length(files)))
			#create array to hold final climatology
			clim<-array(NA,dim=c(lon,lat,12))
			#generate levels for each month in each timeset
			gg<-gl(12, 1, 12*length(files))
			#loop to put data in means array
			for (f in 1:length(files))
			{
				means[,,,f]<-get.var.ncdf(nc,climvar,start=c(1,1,(1+((f-1)*12),count=c(-1,-1,12))
			}
			#take mean over all timesets (4th dimension)
			for (i in 1:lon)
			{
				for (j in 1:lat)
				{
					clim[i,j,]<-vaggregate(means[i,j,,],gg,mean)
				}
			}
			##add data to netCDF file. Overwrites old values for first timeset, but other timesets will still be there
			put.var.ncdf(nc,climvar,clim,count=c(-1,-1,12))
		close.ncdf(nc)



	}
}

