############################# CALCULATING MEANS	####################
# you can calculate means for selected months using the code below. Annual means obviously require all months, whereas seasonal means are derived from certain months only. 
#instructions for specific adaptations are given along the way. Generally tab spaces indicate where things need to be completed 

#section 1
data<-	"CCSM4_pr_piC_anom_tropics.nc"		##netCDF filename (in quotation marks!) to extract data from/write data to
meanvar<- "winter_mean"		##new variable name, suggested format "xxx_mean" where xxx are months/seasons 
varlong<-"Nov-Mar Mean Precipitation 250-1300"		##new variable long name, suggested format "xxx Mean Precipitation 250-1300"


######## ONLY DONE ONCE PER VARIABLE ###################
# 1. ########define new variable
nc<-open.ncdf(data, write=T)
##define variable
newmean<-var.def.ncdf(meanvar, "mm/day", list(nc$dim$Lat,nc$dim$Longitude,nc$dim$Time),missval=NA, 
longname=varlong)  

##add to netcdf file
nc<-var.add.ncdf(nc,newmean) 
close.ncdf(nc)
#######################################################

#section 2
set<-	3	##if dataset is too large, this breaks it into n sets along lon. 
n<-3		##number of the set being calculated. 

#section 3
mm<-	5	##number of months retained
start<- 4	##first month to remove/select
end<-	10	##last month to remove/select
	#winter NDJFM = 4:10 (remove)
	#summer MJJAS = 5:9 (select)

#section 4
#--> check in section after completing step 3

#section 5
#--> check in section after completing step 4 for all n sets





##REPEAT AS OFTEN AS NECESSARY (n TIMES) TO GET FULL DATASET
# 2. ########get data
nc<-open.ncdf(data, write=T) 
lon<-nc$dim$Lon$len/set 
##retrieve data
nc.1<-get.var.ncdf(nc,"ClimAnom",start=c(1,((n-1)*lon)+1,1),count=c(-1,lon,-1)) ##data ClimAnom[lat,lon,time]
##define parameters
time<-nc$dim$Time$len ##number of timesteps
yy<-time/12 ## number of years 
lat<-nc$dim$Lat$len ##number of lats
close.ncdf(nc)

##label the time dimension to be able to select required months
yr<-gl(12, 1, time,labels=c("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec")) #month levels
##make nc.1 into matrix[loc,time]
dim(nc.1)<-c(lat*lon,time)
##append column names (months) to dataset matrix
colnames(nc.1)<-yr


# 3. ####### make months selection (if any)

###you can either:
	#a) remove months (if the required ones are not consecutive within a year); or 
	#b) directly extract a number of months. The former is applicable for e.g. winter, the latter can be used for e.g. summer.

####a) selecting months to remove
	remove<-colnames(nc.1)[start:end] 
	#remove cols of unwanted months
	nc.1<-nc.1[,!colnames(nc.1) %in% remove]
	###if the months aren't consecutive within a year(e.g. jan, feb, dec) it is arguably better to remove the first and last 	records to ensure you calculate averages based on actual seasons (e.g. dec-jan) rather than within calendar years.

	#remove first and/or last months 
	nc.1<-nc.1[,4:(ncol(nc.1)-(12-end))] 
	#update number of years, as you will have removed a year essentially
	yy<-ncol(nc.1)/mm	


####b) directly selecting months 
	extract<-colnames(nc.1)[start:end] 
	#extract cols of wanted months
	nc.1<-nc.1[,colnames(nc.1) %in% extract]


# 4. ############ take mean
##create one level for each year
y<-gl(yy,mm)

library(plyr) #needed for vaggregate
library(reshape) #possibly needed???

##loop taking mean for each lat*lon grid 
ncmean<-array(0,dim=c(lat*lon,yy))

for (i in 1:(lat*lon))
{
ncmean[i,]<-vaggregate(nc.1[i,],y,mean)
}


##restore dataset dimensions
dim(ncmean)<-c(lat,lon,yy)

##rename dataset to secure it
ncmean3		<-ncmean #label according to number of dataset section (i.e. n)



##REPEAT STEPS 2-4 AS OFTEN AS NECESSARY

# 5. ###########combining all sections
ncmean<-abind(ncmean1		#section 1
,ncmean2	#section 2
,ncmean3	#section 3... etc - add more as necessary
,along=	2)



##WHEN ALL DATASET HAS BEEN COMBINED
# 6. ###########writing to NetCDF
nc<-open.ncdf(data, write=T)
##put data into variable
put.var.ncdf(nc,meanvar,ncmean, start=c(1,1,1),count=c(-1,-1,yy)) 
close.ncdf(nc)



