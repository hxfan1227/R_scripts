MEANS<-function(vb,model,period,dx=4,dy=2,units) ##dx=4,dy=2 should work for all model dataset (checked in piControl)
# X,Y is full spatial grid
#dy,dy is number of blocks to divide X,Y into
#x,y is individual block
{
#source baby functions
source("/home/lmkh201/Documents/R/R_scripts/functions/varclim.R")
source("/home/lmkh201/Documents/R/R_scripts/functions/varanom.R")
source("/home/lmkh201/Documents/R/R_scripts/functions/anommean.R")

#loop all datachunks
	for (x in 1:dx)
	{
		for (y in 1:dy)
		{
			varclim(x,y)
			varanom(x,y)
			anommean(x,y)
		}
	}
}





