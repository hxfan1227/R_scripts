MEANS<-function(vb,model,period,dx=4,dy=2,units) {##dx=4,dy=2 should work for all model dataset (checked in piControl)
# X,Y is full spatial grid
#dy,dy is number of blocks to divide X,Y into
#x,y is individual block


print(is.character(vb))
model<-as.character(model)
period<-as.character(period)
units<-as.character(units)

#source baby functions
source("/home/lmkh201/Documents/R/R_scripts/functions/varclim.R")
source("/home/lmkh201/Documents/R/R_scripts/functions/varanom.R")
source("/home/lmkh201/Documents/R/R_scripts/functions/anommean.R")

#loop all datachunks
	for (x in 1:dx)
	{
		for (y in 1:dy)
		{
			print(paste("Starting varclim",x,y,sep="."))			
			varclim(vb,model,period,x,y,units)

			print(paste("Starting varanom",x,y,sep="."))	
			varanom(vb,model,period,x,y)

			print(paste("Starting anommean",x,y,sep="."))	
			anommean(vb,model,period,x,y)
		}
	}

}





