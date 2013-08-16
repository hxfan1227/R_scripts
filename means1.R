############################# CALCULATING MEANS	####################
# you can calculate means for selected months using the code below. Annual means obviously require all months, whereas seasonal means are derived from certain months only. 
#instructions for specific adaptations are given along the way. Generally tab spaces indicate where things need to be completed 
# Repeat for all lon sections!!


####parameters
set<-	3	##if dataset is too large, this breaks it into n sets along lon. 
n<-1		##number of the set being calculated. 

data<-	"CCSM4_pr_mh_clim1.nc"		##netCDF filename (in quotation marks!) to extract data from/write data to
anomvar<-"anomaly1"		##anomaly variable used

nc<-open.ncdf(data, write=T)
londim<-nc$dim$lon1		##lon dimension (change lon[n])
close.ncdf(nc)

meanvar<- "annual_mean"		##new variable name, suggested format "xxx_mean" where xxx are months/seasons 
varlong<-"Annual Mean Precipitation 1850-2005"		##new variable long name, suggested format "xxx Mean Precipitation 250-1300"
units<-"mm/day" 		##variable units

mm<-	12	##number of months retained
start<- 1	##first month to remove/select
end<-	12	##last month to remove/select
	#annual = 1:12
	#winter NDJFM = 4:10 (remove)
	#summer MJJAS = 5:9 (select)


######## Creating new netCDF variable ###################
# 1. ########define new variable
nc<-open.ncdf(data, write=T)
##define variable
newmean<-var.def.ncdf(meanvar, units, list(londim,nc$dim$lat,nc$dim$time),missval=NA, 
longname=varlong)  

##add to netcdf file
nc<-var.add.ncdf(nc,newmean) 


close.ncdf(nc)

###### Calculating mean #########################
# 2. ########get data
nc<-open.ncdf(data, write=T) 
lon<-londim$len 
time<-nc$dim$time$len ##number of timesteps
yy<-time/12 ## number of years 
lat<-nc$dim$lat$len ##number of lats

##retrieve data
anom<-get.var.ncdf(nc,anomvar,start=c(1,1,1),count=c(-1,-1,-1)) ##data anomaly1[lon,lat,time]

close.ncdf(nc)


##label the time dimension to be able to select required months
yr<-gl(12, 1, time,labels=c("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec")) #month levels
##make nc.1 into matrix[loc,time]
dim(anom)<-c(lon*lat,time)
##append column names (months) to dataset matrix
colnames(anom)<-yr






# 3. ####### make months selection (if any)

###you can either:
	#a) remove months (if the required ones are not consecutive within a year); or 
	#b) directly extract a number of months. The former is applicable for e.g. winter, the latter can be used for e.g. summer.




####a) WINTER: selecting months to remove
	remove<-colnames(anom)[start:end] 
	#remove cols of unwanted months
	anom<-anom[,!colnames(anom) %in% remove]
	###if the months aren't consecutive within a year(e.g. jan, feb, dec) it is arguably better to remove the first and last 	records to ensure you calculate averages based on actual seasons (e.g. dec-jan) rather than within calendar years.

	#remove first and/or last months 
	anom<-anom[,4:(ncol(anom)-(12-end))] 
	#update number of years, as you will have removed a year essentially
	yy<-ncol(anom)/mm	





####b) ANNUAL/SUMMER: directly selecting months 
	extract<-colnames(anom)[start:end] 
	#extract cols of wanted months
	anom<-anom[,colnames(anom) %in% extract]







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
}}



# 5. ###########writing to NetCDF
nc<-open.ncdf(data, write=T)
##put data into variable
put.var.ncdf(nc,meanvar,ncmean, start=c(1,1,1),count=c(-1,-1,yy)) 


close.ncdf(nc)






