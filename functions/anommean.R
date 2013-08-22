##### parameters
# ncname	Character string; NetCDF filename extract data from/write data to
# anomvar	Anomaly variable used. Default "anomaly"
# winter	Boolean. If TRUE, extracted months are NOT consecutive within a year(like in winter) and those in between will be removed.
# meanvar	new variable name, suggested format "xxx_mean" where xxx are months/seasons 
# varlong	new variable long name, suggested format "xxx Mean Precipitation [date/period]"
# units		variable units. Default "mm/day"
# start		first month to remove or select
# end		last month to remove or select
		#annual = 1:12
		#winter NDJFM = 4:10 (remove)
		#summer MJJAS = 5:9 (select)


anommean<-function(ncname,anomvar="anomaly",winter=FALSE,meanvar,varlong,units="mm/day"  start,end)
{
	# 1. ########define new variable
	nc<-open.ncdf(ncname, write=T)
		##define variable
		newmean<-var.def.ncdf(meanvar, units, list(nc$dim$lon,nc$dim$lat,nc$dim$time),missval=NA, 
		longname=varlong)  
		##add to netcdf file
		nc<-var.add.ncdf(nc,newmean) 
	close.ncdf(nc)

	###### Calculating mean #########################
	# 2. ########get data
	nc<-open.ncdf(ncname, write=T) 
		lon<-nc$dim$lon$len ##number of lons
		lat<-nc$dim$lat$len ##number of lats
		time<-nc$dim$time$len ##number of timesteps
		yy<-time/12 ## number of years 

		##retrieve data
		anom<-get.var.ncdf(nc,anomvar,start=c(1,1,1),count=c(-1,-1,-1)) ##data anomaly1[lon,lat,time]
	close.ncdf(nc)

	##label the time dimension to be able to select required months
	yr<-gl(12, 1, time,labels=c("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec")) #month levels
	##make nc.1 into matrix[loc,time]
	dim(anom)<-c(lon*lat,time)
	##append column names (months) to dataset matrix
	colnames(anom)<-yr

	######## Creating new netCDF variable ###################
	nc<-open.ncdf(ncname, write=T)
		##define variable
		newmean<-var.def.ncdf(meanvar, units, list(nc$dim$lon,nc$dim$lat,nc$dim$time),missval=NA, 
		longname=varlong)  
		##add to netcdf file
		nc<-var.add.ncdf(nc,newmean) 
	close.ncdf(nc)

	# 3. ####### make months selection
	if (winter==TRUE) ####a) WINTER: selecting months to remove
	{	
		remove<-colnames(anom)[start:end] 
		#remove cols of unwanted months
		anom<-anom[,!colnames(anom) %in% remove]
		###if the months aren't consecutive within a year(e.g. jan, feb, dec) it is arguably better to remove the first and last records to ensure you calculate averages based on actual seasons (e.g. dec-jan) rather than within calendar years.

		#remove first and/or last months 
		anom<-anom[,start:(ncol(anom)-(12-end))] 
		#update number of years, as you will have removed a year essentially
		mm<-12-length(start:end)
		yy<-ncol(anom)/mm
	}
	else ####b) ANNUAL/SUMMER: directly selecting months 
	{
		extract<-colnames(anom)[start:end] 
		#extract cols of wanted months
		anom<-anom[,colnames(anom) %in% extract]
	}

	##restore dataset dimensions
	dim(anom)<-c(lon,lat,yy*mm)

	# 4. ############ take mean
	##create one level for each year
	y<-gl(yy,mm)

	##loop taking mean for each [lat,lon] grid 
	ncmean<-array(0,dim=c(lon,lat,yy))

	for (i in 1:lon)
	{
		for (j in 1:lat)
		{
			ncmean[i,j,]<-vaggregate(anom[i,j,],y,mean)
		}
	}

	# 5. ###########writing to NetCDF
	nc<-open.ncdf(ncname, write=T)
		##put data into variable
		put.var.ncdf(nc,meanvar,ncmean, start=c(1,1,1),count=c(-1,-1,yy)) 
	close.ncdf(nc)
}


