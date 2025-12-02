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

This package will include all WOD QA/QC functions.

'''

def merge_netcdf(
	file_input,
	file_output,
	max_depth=5000,
	instr_id=''
	):
	'''
	file_input: input location for .nc files
	file_output: output location for yearly netcdf files
	max_depth (default 5000): maximum depth bin averaging to in meters
	'''

	#Temporary
	#file_input = '/gpfs/fs7/dfo/dpnm/joc000/Data/CASTS/Data_Input/WOD/OSD_MBT/'
	#file_output = '/gpfs/fs7/dfo/dpnm/joc000/Data/CASTS/Data_Input/WOD/data_processed/'
	#max_depth = 5000

	#Determine all file names
	file_name = os.listdir(file_input)

	#Define variables of interest
	temp = np.full([np.size(file_name),max_depth],np.nan)
	saln = np.full([np.size(file_name),max_depth],np.nan)
	lat = np.full(np.size(file_name),np.nan)
	lon = np.full(np.size(file_name),np.nan)
	time = np.empty(np.size(file_name),dtype='datetime64[m]')
	years = np.empty(np.size(file_name),dtype=int)
	instrument_id = np.empty(np.size(file_name),dtype='U70')
	source = np.full(np.size(file_name),'WOD')
	cruise_id = np.empty(np.size(file_name),dtype='U40')
	platform_id = np.empty(np.size(file_name),dtype='U80')

	#Cycle through each file
	for i,file in enumerate(file_name):
		#Open the file
		ds = xr.open_dataset(file_input+file)
		#Record the variables of interest
		if np.isin('Temperature',list(ds.variables)):
			temp[i,:] = stats.binned_statistic(
				ds.z.values,
				ds.Temperature.values,
				'mean',bins = np.arange(max_depth+1)
				).statistic
		if np.isin('Salinity',list(ds.variables)):
			saln[i,:] = stats.binned_statistic(
				ds.z.values,
				ds.Salinity.values,
				'mean',bins = np.arange(max_depth+1)
				).statistic
		lat[i] = float(ds.lat.values)
		lon[i] = float(ds.lon.values)
		time[i] = ds.time.values
		years[i] = int(ds['time.year'].values)
		if np.isin('WOD_cruise_identifier',list(ds.variables)):
			cruise_id[i] = str(ds.WOD_cruise_identifier.values.astype(str))
		if np.isin('Temperature_Instrument',list(ds.variables)):
			instrument_id[i] = str(ds.Temperature_Instrument.values.astype(str))
		else:
			instrument_id[i] = instr_id
		if np.isin('Platform',list(ds.variables)):
			platform_id[i] = str(ds.Platform.values.astype(str))

	#Sort according to time
	time_filt = np.argsort(time)
	temp = temp[time_filt]
	saln = saln[time_filt]
	lat = lat[time_filt]
	lon = lon[time_filt]
	file_name = np.array(file_name)[time_filt]
	time = time[time_filt]
	years = years[time_filt]
	cruise_id = cruise_id[time_filt]
	instrument_id = instrument_id[time_filt]
	platform_id = platform_id[time_filt]


	#Isolate and save into yearly files
	for year in np.unique(years):
		#Create a year filter
		year_filt = years == year

		#Save the data as a netcdf file 
		#Set up the .nc file
		nc_out = nc.Dataset(file_output+str(year)+'.nc','w')
		
		#File information
		nc_out.Conventions = 'CF-1.6'
		nc_out.title = 'WOD Yearly Netcdf File'
		nc_out.institution = 'Northwest Atlantic Fisheries Centre, Fisheries and Oceans Canada'
		nc_out.description = 'Output by jonathan.coyne@dfo-mpo.gc.ca'
		nc_out.history = 'Created ' + tt.ctime(tt.time())

		#Create dimensions
		time_nc = nc_out.createDimension('time', None) #use date2 for this
		level = nc_out.createDimension('level', max_depth) 

		#Create coordinate variables
		times = nc_out.createVariable('time', np.float64, ('time',))
		levels = nc_out.createVariable('level', np.int32, ('level',))

		#Create 1D variables
		latitudes = nc_out.createVariable('latitude', np.float32, ('time'), zlib=True)
		longitudes = nc_out.createVariable('longitude', np.float32, ('time'), zlib=True)
		file_names = nc_out.createVariable('file_names',str,('time'),zlib=True)
		cruise_ids = nc_out.createVariable('cruise_id',str,('time'),zlib=True)
		instrument_ids = nc_out.createVariable('instrument_id',str,('time'),zlib=True)
		platform_ids = nc_out.createVariable('platform_id',str,('time'),zlib=True)

		#Create 2D variables
		temps = nc_out.createVariable('temperature', np.float32, ('time', 'level'), zlib=True, fill_value=-9999)
		salns = nc_out.createVariable('salinity', np.float32, ('time', 'level'), zlib=True, fill_value=-9999)

		#Variable Attributes
		latitudes.units = 'degree_north'
		longitudes.units = 'degree_east'
		times.units = 'seconds since 1900-01-01 00:00:00'
		times.calendar = 'gregorian'
		levels.units = 'dbar'
		levels.standard_name = "pressure"
		#levels.valid_range = np.array((0.0, 5000.0))
		levels.valid_min = 0
		temps.units = 'Celsius'
		temps.long_name = "Water Temperature" # (may be use to label plots)
		temps.standard_name = "sea_water_temperature"
		salns.long_name = "Practical Salinity"
		salns.standard_name = "sea_water_salinity"
		salns.units = "1"
		salns.valid_min = 0

		#Fill in the 1D variables 
		latitudes[:] = lat[year_filt]
		longitudes[:] = lon[year_filt]
		times[:] = time[year_filt]
		file_names[:] = np.array(file_name)[year_filt]
		cruise_ids[:] = cruise_id[year_filt]
		instrument_ids[:] = instrument_id[year_filt]
		platform_ids[:] = platform_id[year_filt]

		#Fill 2D structure
		temps[:,:] = temp[year_filt]
		salns[:,:] = saln[year_filt]

		#Convert to time stamps
		time_stamps = [pd.Timestamp(i).to_pydatetime() for i in time[year_filt]]
		times[:] = nc.date2num(time_stamps, units=times.units, calendar=times.calendar)
		levels[:] = np.arange(max_depth)

		#Save and close the .nc file
		nc_out.close()
		print(year)




