##### parameters
		


anommean<-function(x,y)
{
	ncname<-paste(model,vb,period,x,y,"nc",sep=".") #gives "model.variable.timeslice.x.y.nc"
	## 1. Get anomalies and dimensions ############
	nc<-open.ncdf(ncname, write=T) 
		#get dimensions	for [x,y]
		lon<-nc$dim$lon$len ##number of lons
		lat<-nc$dim$lat$len ##number of lats
		time<-nc$dim$time$len ##number of timesteps
		yy<-time/12 ## number of years 
		#retrieve anomaly data
		anomvar<-paste("anomaly", x,y, sep = ".") #"anomaly.x.y"
		anom<-get.var.ncdf(nc,anomvar,count=c(-1,-1,-1)) 	
	close.ncdf(nc)
	## 2. label the time dimension with months #########
	yr<-gl(12, 1, time,labels=c("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec")) #month levels
	#make anom into matrix[loc,time]
	dim(anom)<-c(lon*lat,time)
	#append column names (months) to dataset matrix
	colnames(anom)<-yr

	# 3. Make months selection ###########
	nc<- open.ncdf(ncname, write=T) 
		#get mean variable names
		grx<-glob2rx("*mean*")
		varnames<-attributes(nc$var)$names
		meantypes<-grep(grx,varnames,value=T) ##all means
		grx<-glob2rx("winter*")
		winter<-grep(grx,meantypes,value=T) #winter
		grx<-glob2rx("summer*")
		summer<-grep(grx,meantypes,value=T) #summer
		grx<-glob2rx("annual*")
		annual<-grep(grx,meantypes,value=T) #annual
	close.ncdf(nc)	
	for (mt in 1:length(meantypes))
	{
		if (meantypes[mt]==annual) { #ANNUAL: using all months
			ann<-anom
			
			##restore dataset dimensions
			mm<-12
			dim(ann)<-c(lon,lat,yy*mm)

			# 4. ############ take mean
			#create one level for each year
			y<-gl(yy,mm)
			#loop taking mean for each [lat,lon] grid 
			ncmean<-array(0,dim=c(lon,lat,yy))
			for (i in 1:lon)
			{
				for (j in 1:lat)
				{
					ncmean[i,j,]<-vaggregate(ann[i,j,],y,mean)
				}
			}

			# 5. ###########writing to NetCDF
			nc<-open.ncdf(ncname, write=T)
				meanvar<-meantypes[mt]
				##put data into variable
				put.var.ncdf(nc,meanvar,ncmean, start=c(1,1,1),count=c(-1,-1,yy)) 
			close.ncdf(nc)

		} else if (meantypes[mt]==winter) { #WINTER: selecting months to remove
			remv<-colnames(anom)[4:10] 
			#remove cols of unwanted months
			wint<-anom[,!colnames(anom) %in% remv]
			#remove first and last months 
			wint<-wint[,4:(ncol(wint)-2)] 
			#update number of years, as you will have removed a year essentially
			mm<-12-length(4:10)
			yys<-ncol(wint)/mm

			##restore dataset dimensions
			dim(wint)<-c(lon,lat,(yys*mm))

			# 4. ############ take mean
			#create one level for each year
			mm<-length(5:9)
			y<-gl(yy,mm)
			#loop taking mean for each [lat,lon] grid 
			ncmean<-array(0,dim=c(lon,lat,yy))
			for (i in 1:lon)
			{
				for (j in 1:lat)
				{
					ncmean[i,j,]<-vaggregate(wint[i,j,],y,mean)
				}
			}

			# 5. ###########writing to NetCDF
			nc<-open.ncdf(ncname, write=T)
				meanvar<-meantypes[mt]
				##put data into variable
				put.var.ncdf(nc,meanvar,ncmean, start=c(1,1,1),count=c(-1,-1,yy)) 
			close.ncdf(nc)

		} else if (meantypes[mt]==summer) { #SUMMER: directly selecting months 
			extr<-colnames(anom)[5:9] 
			#extract cols of wanted months
			summ<-anom[,colnames(anom) %in% extr]

			##restore dataset dimensions
			mm<-length(5:9)
			dim(summ)<-c(lon,lat,(yy*mm))

			# 4. ############ take mean
			#create one level for each year
			y<-gl(yy,mm)
			#loop taking mean for each [lat,lon] grid 
			ncmean<-array(0,dim=c(lon,lat,yy))
			for (i in 1:lon)
			{
				for (j in 1:lat)
				{
					ncmean[i,j,]<-vaggregate(summ[i,j,],y,mean)
				}
			}

		# 5. ###########writing to NetCDF
		nc<-open.ncdf(ncname, write=T)
			meanvar<-meantypes[mt]
			##put data into variable
			put.var.ncdf(nc,meanvar,ncmean, start=c(1,1,1),count=c(-1,-1,yy)) 
		close.ncdf(nc)
		}
	}

}


