##### parameters
# winter	Boolean. If TRUE, extracted months are NOT consecutive within a year(like in winter) and those in between will be removed.
# start		first month to remove or select
# end		last month to remove or select
		#annual = 1:12
		#winter NDJFM = 4:10 (remove)
		#summer MJJAS = 5:9 (select)


anommean<-function(winter=FALSE, start,end)
{
	ncname<-paste(model,vb,period,x,y,".nc",sep=".") #gives "model.variable.timeslice.x.y.nc"
	## 1. Get anomalies and dimensions ############
	nc<-open.ncdf(ncname, write=T) 
		#get dimensions	for [x,y]
		lon<-nc$dim$lon$len ##number of lons
		lat<-nc$dim$lat$len ##number of lats
		time<-nc$dim$time$len ##number of timesteps
		yy<-time/12 ## number of years 
		#retrieve anomaly data
		anomvar<-paste("anomaly", x,y, sep = ".") #"anomaly.x.y"
		anom<-get.var.ncdf(nc,anomvar,count=c(-1,-1,-1)) 	close.ncdf(nc)
	## 2. label the time dimension with months #########
	yr<-gl(12, 1, time,labels=c("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec")) #month levels
	#make anom into matrix[loc,time]
	dim(anom)<-c(lon*lat,time)
	#append column names (months) to dataset matrix
	colnames(anom)<-yr

	# 3. Make months selection ###########

nc<- open.ncdf(ncname, write=T) 
	








	

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


