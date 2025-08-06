import glob
import os
import sys
import re
import warnings
from scipy import stats,interpolate
import numpy as np
import xarray as xr
import netCDF4 as nc
import time as tt
import pandas as pd
import matplotlib
matplotlib.interactive(True)
import matplotlib.pyplot as plt



'''

The purpose of this script is to streamline the CASTS update process.
All files and sources should be able to be produced from one central script.

This package will include all IML QA/QC functions.

'''


def raw_to_netcdf(
	file_input,
	file_output,
	max_depth=5000,
	):

	'''
	file_input: input location for climate_ae.nc file from Jean-Luc
	file_output: output location for yearly netcdf files
	max_depth (default 5000): maximum depth bin averaging to in meters
	'''

	#Import the raw data
	ds = xr.open_dataset(file_input+'climate_ae.nc')

	#Determine the year range covered
	years = np.unique(ds['time.year'].values)

	#Determine the cast start location for each 
	cast_start = np.concatenate((np.array([0]),np.cumsum(ds.cast_size.values)[:-1]))

	#Set up empty 1D variables
	original_file = ds.original_file.values.astype(str)
	cruise_id = ds.cruise_id.values.astype(str)
	station_id = ds.station_id.values.astype(str)
	platform_id = ds.platform_id.values.astype(str)

	#Set up empty temperature and salinity arrays
	temp_2D = np.full((ds.time.size,max_depth),np.nan)
	saln_2D = np.full((ds.time.size,max_depth),np.nan)

	#Load in the temperature, salinity, pressure
	ds_temp = ds.temperature.values
	ds_saln = ds.salinity.values
	ds_pres = ds.pressure.values

	#Cycle through each of the casts
	for i,value in enumerate(cast_start):

		#Slice the temperature, pressure, and salinity
		if i+1 != cast_start.size:
			temp_slice = ds_temp[value:cast_start[i+1]]
			saln_slice = ds_saln[value:cast_start[i+1]]
			pres_slice = ds_pres[value:cast_start[i+1]]
		else:
			ending_index = value + ds.temperature_row_size[-1].values
			temp_slice = ds_temp[value:ending_index]
			saln_slice = ds_saln[value:ending_index]
			pres_slice = ds_pres[value:ending_index]	

		#If the cast is empty, skip
		if pres_slice.size == 0:
			continue

		#Populate the 2D arrays by binning first
		temp_slice = stats.binned_statistic(
			pres_slice,
			temp_slice,
			'mean',
			bins = np.arange(max_depth+1)
			).statistic
		saln_slice = stats.binned_statistic(
			pres_slice,
			saln_slice,
			'mean',
			bins = np.arange(max_depth+1)
			).statistic

		#Populate the 2D arrays
		temp_2D[i,:] = temp_slice
		saln_2D[i,:] = saln_slice 

	#Remove the -99 values, cutoffs 
	temp_2D[temp_2D <= -2] = np.nan
	temp_2D[temp_2D > 35] = np.nan
	saln_2D[saln_2D <= 0] = np.nan
	saln_2D[saln_2D > 45] = np.nan

	#Record the latitude, longitude, time
	years = ds['time.year'].values
	ds_time = ds.time.values
	ds_latitude = ds.lat.values
	ds_longitude = ds.lon.values

	#Cycle through each of the years now
	for year in np.unique(years):

		#Create a slice of the data for a specific year
		time_slice = ds_time[years == year]
		lat_slice = ds_latitude[years == year]
		lon_slice = ds_longitude[years == year]
		temp_slice = temp_2D[years == year,:]
		saln_slice = saln_2D[years == year,:]
		original_file_slice = original_file[years == year]
		cruise_id_slice = cruise_id[years == year]
		station_id_slice = station_id[years == year]
		platform_id_slice = platform_id[years == year]

		#Create a filter
		filt = np.argsort(time_slice)

		#Save the data as a netcdf file 
		#Set up the .nc file
		nc_out = nc.Dataset(file_output+str(year)+'.nc','w')
		
		#File information
		nc_out.Conventions = 'CF-1.6' #Ask Fred what this means
		nc_out.title = 'MLI Yearly Netcdf File' #Temporary title for the .nc file
		nc_out.institution = 'Northwest Atlantic Fisheries Centre, Fisheries and Oceans Canada'
		nc_out.description = 'Output by jonathan.coyne@dfo-mpo.gc.ca'
		nc_out.history = 'Created ' + tt.ctime(tt.time())

		#Create dimensions
		time = nc_out.createDimension('time', None) #use date2 for this
		level = nc_out.createDimension('level', max_depth) 

		#Create coordinate variables
		times = nc_out.createVariable('time', np.float64, ('time',))
		levels = nc_out.createVariable('level', np.int32, ('level',))

		#Create 1D variables
		latitudes = nc_out.createVariable('latitude', np.float32, ('time'), zlib=True)
		longitudes = nc_out.createVariable('longitude', np.float32, ('time'), zlib=True)
		file_names = nc_out.createVariable('file_names',str,('time'),zlib=True)
		cruise_ids = nc_out.createVariable('cruise_id',str,('time'),zlib=True)
		station_ids = nc_out.createVariable('station_id',str,('time'),zlib=True)
		platform_ids = nc_out.createVariable('platform_id',str,('time'),zlib=True)

		#Create 2D variables
		temp = nc_out.createVariable('temperature', np.float32, ('time', 'level'), zlib=True, fill_value=-9999)
		saln = nc_out.createVariable('salinity', np.float32, ('time', 'level'), zlib=True, fill_value=-9999)

		#Variable Attributes
		latitudes.units = 'degree_north'
		longitudes.units = 'degree_east'
		times.units = 'seconds since 1900-01-01 00:00:00'
		times.calendar = 'gregorian'
		levels.units = 'dbar'
		levels.standard_name = "pressure"
		#levels.valid_range = np.array((0.0, 5000.0))
		levels.valid_min = 0
		temp.units = 'Celsius'
		temp.long_name = "Water Temperature" # (may be use to label plots)
		temp.standard_name = "sea_water_temperature"
		saln.long_name = "Practical Salinity"
		saln.standard_name = "sea_water_salinity"
		saln.units = "1"
		saln.valid_min = 0

		#Fill in the 1D variables 
		latitudes[:] = lat_slice[filt]
		longitudes[:] = lon_slice[filt]
		times[:] = time_slice[filt]
		file_names[:] = original_file_slice[filt]
		cruise_ids[:] = cruise_id_slice[filt]
		station_ids[:] = station_id_slice[filt]
		platform_ids[:] = platform_id_slice[filt]

		#Fill 2D structure
		temp[:,:] = temp_slice[filt]
		saln[:,:] = saln_slice[filt]

		#Convert to time stamps
		time_stamps = [pd.Timestamp(i).to_pydatetime() for i in time_slice[filt]]
		times[:] = nc.date2num(time_stamps, units=times.units, calendar=times.calendar)
		levels[:] = np.arange(max_depth)

		#Save and close the .nc file
		nc_out.close()

