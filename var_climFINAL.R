#### Calculate final climatology ###################
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

