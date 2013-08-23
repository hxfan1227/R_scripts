##### parameters
# vb		Character string; climatological variable to analyse (as in raw dataset filenames)
# model		Character string; model that data is taken from (as in raw dataset filenames)
# period		Character string; PMIP3 timeslice (pic,h,1k,mh) or other data range (yy-yy)
# x		Numerical value; dataset LON slice
# y		Numerical value; dataset LAT slice


varanom<-function(x,y)
{
	#get dimensions	for [x,y]
	ncname<-paste(model,vb,period,x,y,"nc",sep=".") #gives "model.variable.timeslice.x.y.nc"
	nc<-open.ncdf(ncname)
		lon<-nc$dim$lon$len
		lat<-nc$dim$lat$len
	close.ncdf(nc)

	##Get filenames of all data timesets within a model timeslice
	filename<-paste(vb,"Amon",model,"*",sep="_") #creates search string "[var]_Amon_[model]_*"
	files<-Sys.glob(filename) #vector with names of all raw datafiles
		nxt<-1

	for (f in 1:length(files))
	{
		data<-files[f] #raw datafile
		### 1. Get variable data ###########
		nc<-open.ncdf(data,write=T)
			time<-nc$dim$time$len
			yy<-time/12
			st<-nxt	#timeset beginning
			nxt<-nxt+time #next timeset's beginning
			var<-get.var.ncdf(nc, vb, start=c(1+((x-1)*lon),1+((y-1)*lat),1), count=c(lon,lat,-1))
		close.ncdf(nc)
		if (vb=="pr") {
			var<-var*86400 #convert to mm/day
		} else {
			var<-var #can be changed if necessary for SST
		}

		##retrieve climatology
		nc<-open.ncdf(ncname)
			climvar<-paste("climatology",x,y,sep=".")
			clim<-get.var.ncdf(nc,climvar,count=c(-1,-1,12))
		close.ncdf(nc)
		#replicate for each year in dataset
		clim<-replicate(yy,abind(clim,along=3))
		#match variable dimensions
		dim(clim)<-dim(var)
		##remove climatological mean from data to get anomaly
		anom<-var-clim

		#####2. Add to netCDF file #####################
		#add data to variable
		nc<-open.ncdf(ncname,write=T)
			anomvar<-paste("anomaly", x,y, sep = ".") #"anomaly.x.y"
			put.var.ncdf(nc,anomvar,anom, start=c(1,1,st), count=c(-1,-1,length(anom[1,1,])))
			put.var.ncdf(nc,"time",st:(nxt-1),start=st,count=time)
		close.ncdf(nc)
	}
}


