**** Instructions for pre-processing model data and calculating EOFs ****

Execute scripts in this order:
1. var_clim[n].R
	There are n scripts, depending on how the dataset is divided along lon. Run all lon sections; within each, run every timeset. These will all be saved as a separate variable to be combined into a final climatology (in step 2).

2. var_climFINAL.R
	ONLY NECESSARY IF >1 TIME SETS
This combines all climatologies from the previous script to create a final climatology over the entire model timeseries for each lon set. It can be adapted for all lon sections. 

3. var_anom[n].R
	calculates climate anomalies using the final climatology from step 2. There are n scripts, depending on how the dataset is divided along lon. Again, run all timesets for all lon sections. It creates one variable, saved in the same netCDF file as the climatologies.
	
4. means.R
	calculates annual and seasonal means of the anomalies from step 3. The type of mean can be selected in the header; generally run this three times PER LON SET to get annual, winter and summer means. Variables are stored in the same netCDF files as climatologies and anomaly.

5. EOF.R
	This script concatenates means for all lons and calculates the EOFs. EOFs, variances, PCs and correlations for the first 3 EOFs get saved in .csv files. It also produces maps and figures which get saved as .png files. 


FUNCTIONS:

means: anommean
