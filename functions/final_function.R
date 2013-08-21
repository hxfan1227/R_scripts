function(data, newdata,xslice=4,yslice=2) ##xslice=4,yslice=2 should work for all model dataset (checked in piControl)
{
	for (x in 1:xslice)
	{
		for (y in 1:yslice)
		{
		varclim[x,y]
		varanom[x,y]
		anommean[x,y]
		}
	}

eof
}
