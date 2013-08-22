function(data, newdata,dx=4,dy=2) ##dx=4,dy=2 should work for all model dataset (checked in piControl)
# X,Y is full spatial grid
#dy,dy is number of blocks to divide X,Y into
#x,y is individual block
{
	for (x in 1:dx)
	{
		for (y in 1:dy)
		{
			varclim[x,y]
			varanom[x,y]
			anommean[x,y]
		}
	}

eof[X,Y]
}



# varclim[x,y]: varclim<-function(vb,model,period,x,y,units)
1. retrieves raw data filenames
2. creates netCDF file [x,y] with variables 
	clim
	anom
	annual mean
	winter mean
	summer mean
3. 	a) if there is only 1 file (i.e. timeset), retrieves data[x,y] and calculates climatology and adds to netCDF file.
	b) else loops the above for all timesets, adding each clim sequentially in time dimension, then calculates final climatology and overwrites this into netCDF file start of time dimension.




# anommean[dx,dy]:
