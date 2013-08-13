**** Instructions for pre-processing model data and calculating EOFs ****

Execute scripts in this order:
1. var_clim[n].R
	There are n scripts, depending on how the dataset is divided along lon. Run all lon sections; within each, run every timeset. These will all be saved as a separate variable to be combined into a final climatology (in step 2).

2. var_climFINAL.R
	This combines all climatologies from the previous script to create a final climatology over the entire model timeseries. It can be adapted for all lon sections.

var_anom[n].R
	calculates climate anomalies using the final climatology from step 2. There are n scripts, depending on how the dataset is divided along lon. Again, run all timesets for all lon sections. This is all saved in the same netCDF file as the climatologies.
	
3. means.R

