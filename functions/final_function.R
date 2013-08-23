function(dx=4,dy=2) ##dx=4,dy=2 should work for all model dataset (checked in piControl)
# X,Y is full spatial grid
#dy,dy is number of blocks to divide X,Y into
#x,y is individual block
{
	for (x in 1:dx)
	{
		for (y in 1:dy)
		{
			print(paste("Starting varclim",x,y,sep="."))			
			varclim[x,y]

			print(paste("Starting varanom",x,y,sep="."))	
			varanom[x,y]

			print(paste("Starting anommean",x,y,sep="."))	
			anommean[x,y]
		}
	}

eof[X,Y]
}




vb<-"pr"
model<-"MIROC-ESM"
period<-"piC"
dx<-4
dy<-2
units<-"mm/day"

source("/home/lmkh201/Documents/R/R_scripts/functions/varclim.R")
source("/home/lmkh201/Documents/R/R_scripts/functions/varanom.R")
source("/home/lmkh201/Documents/R/R_scripts/functions/anommean.R")

varclim(vb=vb,model=model,period=period,x=x,y=y,units=units)

varanom(vb=vb,model=model,period=period,x=x,y=y)

anommean(vb=vb,model=model,period=period,x=x,y=y)





# varclim[x,y]: varclim<-function(vb,model,period,x,y,units)
1. retrieves raw data filenames
2. creates netCDF file [x,y] with variables 
	climatology
	anomaly
	annual_mean
	winter_mean
	summer_mean
3. 	a) if there is only 1 file (i.e. timeset), retrieves data[x,y] and calculates climatology and adds to netCDF file variable "climatology.x.y".
	b) else loops the above for all timesets, adding each clim sequentially in time dimension, then calculates final climatology and overwrites this into netCDF file variable "climatology.x.y" at start of time dimension.

# varanom[x,y]: varanom<-function(vb,model,period,x,y)
1. retrieves raw data filenames 
2. gets raw variable data and climatology
3. calculates anomalies for each timeset and adds to netCDF file variable "anomaly.x.y".

# anommean[x,y]: anommean<-function(model,vb,period,x,y)
1. gets anomaly data
2. labels time dimension with month names
3. if/else loop to extract appropriate months for 
	a) annual 
	b) winter
	c) summer
4. calculate means and writes them to netCDF file variables "*_mean".
