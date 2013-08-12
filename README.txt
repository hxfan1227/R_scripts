**** Instructions for pre-processing model data and calculating EOFs ****

Execute scripts in this order:
1. var_clim[n].R
	There are n scripts, depending on how the dataset is divided along lon. Run all sections; each time-section datafile will need its own variable which are then combined into a final climatology before calculating the anomaly (in step 2).

2. var_anom.R
	
3. means.R

