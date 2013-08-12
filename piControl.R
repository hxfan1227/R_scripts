######
#setup
#all of these must be done for the rest to work properly
cd /home/lmkh201/Documents/R/model_data/pr/piControl/CCSM4/r1i1p1
R

.libPaths(new="/home/lmkh201/Documents/R/R_packages")


library(maps)
library(mapdata)
library(mapproj)
library(fBasics) #masks norm
library(abind)

library(ncdf)
library(akima) 
library(pracma) #masks akimaInterp, inv, kron, pascal
library(clim.pact) #masks mod, num2str, map
library(R.methodsS3)
library(R.oo) #masks getDescription, getTitle, attach, getClasses, getMethods, attach, detach, gc, load, save
library(R.utils) #masks insert, timestamp, cat, commandArgs, getOption, inherits, isOpen, parse, warnings
library(plyr) #needed for vaggregate
library(reshape) #possibly needed??? masks rename, round_any






###PRE-INDUSTRIAL CONTROL
##calculate climatological mean
data<-"pr_Amon_CCSM4_piControl_r1i1p1_080001-130012.nc"





###################create netCDF files
##Define dimensions
x <- dim.def.ncdf( "Lat", "degrees_north", latval)
y <- dim.def.ncdf( "Lon", "degrees_east", lonval)
m <- dim.def.ncdf( "Month", "months", 1:12)




### data
nc.1<- retrieve.nc(data,v.nam="pr",t.rng=c("800","1000")) ##get data, adjust t.rng to max 300 years if necessary
mm<-nc.1$mm ##vector rep(1:12), 
lat<-length(nc.1$lat) ##latitude res
lon<-length(nc.1$lon) ##longitude res
yy<-length(nc.1$yy)/12 ##number of years
latval<-nc.1$lat
lonval<-nc.1$lon

nc.1<-86400*nc.1$dat ##mm/day

mean1<-array(NA,dim=c(12,lat,lon)) #monthly mean array
for (i in 1:lat)
{
for (j in 1:lon)
{
mean1[,i,j]<-tapply(nc.1[,i,j],mm,mean,na.rm=T) ##monthly mean for each grid box 
}}

rm(nc.1)

