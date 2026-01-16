import xarray as xr
import netCDF4 as nc
import numpy as np
import time as tt
import warnings
import os
import math
import random
import pandas as pd
from scipy import stats,interpolate
from scipy.signal import butter, lfilter, freqz
import scipy
import matplotlib
matplotlib.interactive(True)
import matplotlib.pyplot as plt

'''
The purpose of this script is to streamline the CASTS update process.
All files and sources should be able to be produced from one central script.

This package will include all NAFC-Oceanography QA/QC functions.
'''

#Test info
input_file = '/gpfs/fs7/dfo/dpnm/joc000/Data/CASTS/Data_Input/NAFC_Oceanography/3_duplicates/1950.nc'
ds = xr.open_dataset(input_file)
temp = ds.temperature.values
saln = ds.salinity.values
levels = ds.level.values
lat = ds.latitude.values
lon = ds.longitude.values
month = ds['time.month'].values

def global_impossible_parameters(
	temp,saln,
	):
	#UNESCO 1990 2.1
	#Define the temp and saln range
	temp_range = [-2.5,35]
	saln_range = [0,40]
	#Isolate the temperatures
	temp[temp <= temp_range[0]] = np.nan
	temp[temp >= temp_range[1]] = np.nan
	#Isolate the salinities
	saln[saln <= saln_range[0]] = np.nan
	saln[saln >= saln_range[1]] = np.nan
	#Return the data
	return(temp,saln)

#Note, UNESCO 2.2 and 2.3 are not applicable for CASTS

def global_profile_envelope(
	temp,saln,
	):
	#UNESCO 1990 2.4
	#Define the temp and saln ranges
	temp_range = {
	'0-50': [0,50,-2.5,35],
	'51-100': [51,100,-2.5,30],
	'101-400': [101,400,-2.5,28],
	'401-1100': [401,1100,-2,27],
	'1101-3000': [1101,3000,-1.5,18],
	'3001-5500': [3001,5500,-1.5,7],}
	saln_range = {
	'0-50': [0,50,0,40],
	'51-100': [51,100,1,40],
	'101-400': [101,400,3,40],
	'401-1100': [401,1100,10,40],
	'1101-3000': [1101,3000,22,38],
	'3001-5500': [3001,5500,33,37],}
	#Cycle through each depth range and remove data
	for i in temp_range:
		#Slice the data and remove suspect data
		temp_slice = temp[:,temp_range[i][0]:temp_range[i][1]]
		temp_slice[temp_slice <= temp_range[i][2]] = np.nan
		temp_slice[temp_slice >= temp_range[i][3]] = np.nan
		temp[:,temp_range[i][0]:temp_range[i][1]] = temp_slice
		saln_slice = saln[:,saln_range[i][0]:saln_range[i][1]]
		saln_slice[saln_slice <= saln_range[i][2]] = np.nan
		saln_slice[saln_slice >= saln_range[i][3]] = np.nan
		saln[:,saln_range[i][0]:saln_range[i][1]] = saln_slice
	#Return the data
	return(temp,saln)

def constant_profile(
	temp,saln,
	):
	#UNESCO 1990 2.5
	#Determine the number of measurements per cast
	temp_nom = np.sum(~np.isnan(temp),axis=1)
	saln_nom = np.sum(~np.isnan(saln),axis=1)
	#Determine the number of unique measurements (minus 1 for nan)
	temp_unique = np.array([np.unique(row).size-1 for row in temp])
	constant_temp = (temp_nom >= 3)*(temp_unique == 1)
	temp[constant_temp] = np.nan
	saln_unique = np.array([np.unique(row).size-1 for row in saln])
	constant_saln = (saln_nom >= 3)*(saln_unique == 1)
	saln[constant_temp] = np.nan
	#Return the data
	return(temp,saln)

def freezing_point(
	temp,saln,
	):
	#UNESCO 1990 2.6
	#Determine if the temperature is below freezing point
	#Citation, UNESCO 1983, Algorithms for Computation of Fundamental Properties of Seawater
	#Note, only completed when temperature and salinity measurements are taken together
	#Isolate salinity between 27 and 35
	saln_freezing = np.copy(saln)
	saln_freezing[saln_freezing <= 27] = np.nan
	saln_freezing[saln_freezing >= 35] = np.nan
	#Determine the Freezing Temperature
	T_freezing = -0.0575*saln_freezing+1.71052e-3*(saln_freezing**(3/2))-2.154996e-4*(saln_freezing**2)-(7/(53e-4*levels))
	#If temperature is greater than Freezing Temperature, do not remove
	temp[temp <= T_freezing] = np.nan
	#Return the data
	return(temp)

def spike_test(
	temp,saln,
	):
	#UNESCO 1990 2.7
	#Based upon WMO/IOC Manuals and Guides #3
	#Define the temperature and salinity thresholds
	temp_th = 2
	saln_th = 0.3
	#Determine the profiles with temperature spikes
	row_1 = temp[:,:-2]
	row_2 = temp[:,1:-1]
	row_3 = temp[:,2:]
	temp_spike = np.abs(row_2 - (row_3+row_1)/2) - np.abs(row_1-row_3)/2
	temp_spike = np.sum(temp_spike >= temp_th,axis=1).astype(bool)
	temp[temp_spike] = np.nan
	#Determine the profiles with salinity spikes
	row_1 = saln[:,:-2]
	row_2 = saln[:,1:-1]
	row_3 = saln[:,2:]
	saln_spike = np.abs(row_2 - (row_3+row_1)/2) - np.abs(row_1-row_3)/2
	saln_spike = np.sum(saln_spike >= saln_th,axis=1).astype(bool)
	saln[saln_spike] = np.nan
	#Return the data
	return(temp,saln)

def top_bottom_spike(
	temp,saln,
	):
	#UNESCO 1990 2.8
	#Determine if spikes are present in the shallowest/deepest measurements
	#Test the surface temperature
	temp_dn,temp_up = -10,10
	temp_spike = temp[:,0] - temp[:,1]
	temp_spike = ((temp_spike < temp_dn)+(temp_spike > temp_up)).astype(bool)
	temp[temp_spike,0] = np.nan
	#Test the bottom temperature
	#Determine the bottom measurements and one above bottom
	for i,row in enumerate(temp):
		if np.where(~np.isnan(row))[0].size != 0:
			bottom_idx = np.where(~np.isnan(row))[0][-1]
			bottom_temp = temp[i,bottom_idx]
			bottom_temp_m1 = temp[i,bottom_idx-1]
			temp_spike = bottom_temp - bottom_temp_m1
			if temp_spike < temp_dn or temp_spike > temp_up:
				temp[i,bottom_idx] = np.nan
	#Test the surface salinity
	saln_dn,saln_up = -5,5
	saln_spike = saln[:,0] - saln[:,1]
	saln_spike = ((saln_spike < saln_dn)+(saln_spike > saln_up)).astype(bool)
	saln[saln_spike,0] = np.nan
	#Test the bottom salinity
	#Determine the bottom measurements and one above bottom
	for i,row in enumerate(saln):
		if np.where(~np.isnan(row))[0].size != 0:
			bottom_idx = np.where(~np.isnan(row))[0][-1]
			bottom_saln = saln[i,bottom_idx]
			bottom_saln_m1 = saln[i,bottom_idx-1]
			saln_spike = bottom_saln - bottom_saln_m1
			if saln_spike < saln_dn or saln_spike > saln_up:
				saln[i,bottom_idx] = np.nan
	#Return the data
	return(temp,saln)

def gradient_test(
	temp,saln,
	):
	#UNESCO 1990 2.9
	#Perform the temperature gradient test
	temp_th = 10
	row_1 = temp[:,:-2]
	row_2 = temp[:,1:-1]
	row_3 = temp[:,2:]
	temp_spike = np.full(temp.shape,np.nan)
	temp_spike[:,:-2] = np.abs(row_2 - (row_1+row_3)/2)
	temp[temp_spike >= temp_th] = np.nan
	#Perform the salinity gradient test
	saln_th = 5
	row_1 = saln[:,:-2]
	row_2 = saln[:,1:-1]
	row_3 = saln[:,2:]
	saln_spike = np.full(saln.shape,np.nan)
	saln_spike[:,:-2] = np.abs(row_2 - (row_1+row_3)/2)
	saln[saln_spike >= saln_th] = np.nan
	#Return the data
	return(temp,saln)

#Note, UNESCO 1993 2.10 Density Inversion Test was not included
#Temperature and salinity data required for test, not always present for each profile

#Note, UNESCO 1993 2.11 Bottom Test
#Implement this later, this is better than our bathymetry test

'''
def temperature_inversion_test(
	temp,saln,
	):
	#UNESCO 1993 2.12
	#Perform a temperature inversion test
'''

def levitus_seasonal_statistic(
	temp,
	saln,
	lat,
	lon,
	month,
	level,
	salnclim_mean='/gpfs/fs7/dfo/dpnm/joc000/Data/CASTS/Data_Input/Other/climatology_files/salt.anal1deg.nc',
	salnclim_stdv='/gpfs/fs7/dfo/dpnm/joc000/Data/CASTS/Data_Input/Other/climatology_files/salt.sd1deg.nc',
	):
	#UNESCO 1993 3.1
	#Perform an outlier test using levitus seasonal statistics atlas
	#Import the salinity mean and stdv
	ds_salnmean = xr.open_dataset(salnclim_mean,decode_times=False)
	ds_salnstdv = xr.open_dataset(salnclim_stdv,decode_times=False)
	ds_salnmean['salt'] = ds_salnmean.salt.where(ds_salnmean.salt > 0.2)
	#Define the months by number, not days
	ds_salnmean['time'] = np.arange(1,12+1)
	ds_salnstdv['time'] = np.arange(1,12+1)
	#Determine the mean and standard deviation at each cast point
	x = xr.DataArray(lon+360, dims='z') #account for difference in lon formatting
	y = xr.DataArray(lat, dims='z')
	t = xr.DataArray(month, dims='z')
	warnings.simplefilter(action='ignore', category=FutureWarning)
	ds_salnmean = ds_salnmean.salt.interp(time=t,lat=y,lon=x,method='linear')
	ds_salnstdv = ds_salnstdv.salt.interp(time=t,lat=y,lon=x,method='linear')
	#Ensure that vertical gridding is consistent
	ds_salnmean = ds_salnmean.interp(level=levels)
	ds_salnstdv = ds_salnstdv.interp(level=levels)
	#Determine the distance to land of each point

	#Compare the CASTS salinity with the Levitus climatology
	saln_mean = ds_salnmean.values
	saln_stdv = ds_salnstdv.values
	#saln_stdv[saln_stdv < 0.01] = np.nan
	upper = saln >= (ds_salnmean.values + (ds_salnstdv.values*5))
	lower = saln <= (ds_salnmean.values - (ds_salnstdv.values*5))
	#Make into one index
	idx = upper.sum(axis=1).astype(bool) + lower.sum(axis=1).astype(bool)


#TESTING AREA, seems too strict as of now...
test_casts = saln[idx]
test_mean = ds_salnmean.values[idx]
test_stdv = ds_salnstdv.values[idx]

i = 1
plt.plot(levels,test_casts[i],marker='.',color='tab:red')
plt.plot(levels,test_mean[i]+(test_stdv[i]*5),color='k')
plt.plot(levels,test_mean[i]-(test_stdv[i]*5),color='k')











def NAFC_climatology(
	input_folder,
	output_folder,
	years,
	max_depth=5000,
	):
	'''
	input_folder: path location to yearly netcdf files0
	output_folder: path output location
	years: years that data is from (str)
	max_depth (default 5000): maximum depth bin averaging to in meters
	'''

	#Define lat and lon on a constant grid across years
	dc = 1
	x = np.arange(-100, -42, dc)
	y = np.arange(35, 80, dc)
	lon, lat = np.meshgrid(x,y)

	#Cycle through each month

	#Set up the empty arrays
	temp_3D_mean = np.full((years.size,y.shape[0]-1,x.shape[0]-1),np.nan)
	saln_3D_mean = np.full((years.size,y.shape[0]-1,x.shape[0]-1),np.nan)
	#temp_3D_stdv = np.full((years.size,y.shape[0]-1,x.shape[0]-1),np.nan)
	#saln_3D_stdv = np.full((years.size,y.shape[0]-1,x.shape[0]-1),np.nan)

	#Cycle through each of the years
	for i,year in enumerate(years):
		#Open the year of interest
		ds = xr.open_dataset(input_folder+year+'.nc')
		#Isolate the month of interest
		ds = ds.isel(time = ds['time.month'] == 8)
		#Determine the gridded mean for temp and saln
		if ds.time.size > 0:
			temp_3D_mean[i] = stats.binned_statistic_2d(
				ds.longitude.values,
				ds.latitude.values,
				ds.temperature[:,0].values,
				'mean',bins=[x,y]
				).statistic.T
		print(year)


