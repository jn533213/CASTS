import glob
import os
import sys
import re
import warnings
import urllib
from scipy import stats,interpolate
import numpy as np
import xarray as xr
import netCDF4 as nc
import time as tt
import pyreadr
import pandas as pd
import matplotlib
matplotlib.interactive(True)
import matplotlib.pyplot as plt



'''

The purpose of this script is to streamline the CASTS update process.
All files and sources should be able to be produced from one central script.

This package will include all NEFSC QA/QC functions.

'''


def download_NEFSC(
	file_output,
	):
	#NEFSC
	address = 'https://comet.nefsc.noaa.gov/erddap/tabledap/ocdbs_v_erddap1.nc?UTC_DATETIME%2Clatitude%2Clongitude%2Cdepth%2Cpressure_dbars%2Csea_water_temperature%2Csea_water_salinity%2Cdissolved_oxygen%2Cfluorescence%2Cpar_sensor%2Ccast_number%2Ccruise_id%2Cpurpose_code%2Cbottom_depth%2CGEAR_TYPE'

	#Download the file
	urllib.request.urlretrieve(address,file_output+'NEFSC_ERDDAP.nc')






'''
TESTING AREA
file_input = '/gpfs/fs7/dfo/dpnm/joc000/Data/CASTS/Data_Input/NEFSC/ERDDAP/ocdbs_v_erddap1_c300_5bd0_e2e1.nc'
file_output = '/gpfs/fs7/dfo/dpnm/joc000/Data/CASTS/Data_Input/NEFSC/individual_netcdf_files/'
file_yearly_output = '/gpfs/fs7/dfo/dpnm/joc000/Data/CASTS/Data_Input/NEFSC/yearly_netcdf_files/'
max_level = 5000
'''

def raw_to_netcdf(
	file_input,
	file_output,
	file_yearly_output,
	max_level=5000
	):


	#Import the file
	ds = xr.open_dataset(file_input)

	#Determine the years present
	years = ds['UTC_DATETIME.year'].values

	#Cycle through each of the years
	for year in np.unique(years):
		ds_year_isel = ds.sel(row = years==year)

		#Isolate the data based upon cruise id
		cruise_id = ds_year_isel.cruise_id.values.astype(str)

		#Cycle through each
		for cruise in np.unique(cruise_id):
			ds_cruise_isel = ds_year_isel.sel(row = cruise_id == cruise)
			cast_number = ds_cruise_isel.cast_number.values

			#Set up the 1D variables
			latitude = np.zeros(np.unique(cast_number).size)
			longitude = np.zeros(np.unique(cast_number).size)
			time_1D = np.zeros(np.unique(cast_number).size).astype('datetime64[s]')
			cruise_ID = np.full(np.unique(cast_number).size, cruise)
			purpose_code = np.zeros(np.unique(cast_number).size).astype(str)
			bottom_depth = np.zeros(np.unique(cast_number).size).astype(str)
			gear_type = np.zeros(np.unique(cast_number).size).astype(str)

			#Set up the 2D variables
			temperature = np.zeros((np.unique(cast_number).size,max_level))
			salinity = np.zeros((np.unique(cast_number).size,max_level))

			#Cycle through each of the casts
			for i,cast in enumerate(np.unique(cast_number)):
				ds_cast_isel = ds_cruise_isel.sel(row = cast_number == cast)

				#Record the 1D variables
				latitude[i] = ds_cast_isel.latitude.values.mean()
				longitude[i] = ds_cast_isel.longitude.values.mean()
				purpose_code[i] = ds_cast_isel.purpose_code.values[0]
				bottom_depth[i] = ds_cast_isel.bottom_depth.values[0]
				gear_type[i] = ds_cast_isel.GEAR_TYPE.values[0]

				#Format the datetime
				date = ds_cast_isel.UTC_DATETIME.values[0].astype('datetime64[m]')
				#decimal_hour = ds_cast_isel.UTC_DECIMAL_HOUR.values[0]
				#if np.isnan(decimal_hour):
				#	hours = 0
				#	minutes = 0
				#else:
				#	hours = int(decimal_hour)
				#	minutes = np.round((decimal_hour*60) % 60,0).astype(int)
				#time_1D[i] = date+np.timedelta64(hours,'h')+np.timedelta64(minutes,'m')
				time_1D[i] = date

				#Record the 2D variables
				temperature[i,:] = stats.binned_statistic(
					ds_cast_isel.depth.values.astype(float),
					ds_cast_isel.sea_water_temperature.values.astype(float),
					'mean',
					bins = np.arange(max_level+1)
					).statistic
				salinity[i,:] = stats.binned_statistic(
					ds_cast_isel.depth.values.astype(float),
					ds_cast_isel.sea_water_temperature.values.astype(float),
					'mean',
					bins = np.arange(max_level+1)
					).statistic

			#Form into a cruise-specific netcdf
			#Set up the .nc file
			nc_out = nc.Dataset(file_output+str(year)+'_'+cruise+'.nc','w')
			
			#File information
			nc_out.Conventions = 'CF-1.6' #Ask Fred what this means
			nc_out.title = 'NEFSC Individual Netcdf File' #Temporary title for the .nc file
			nc_out.institution = 'Northwest Atlantic Fisheries Centre, Fisheries and Oceans Canada'
			nc_out.description = 'Output by jonathan.coyne@dfo-mpo.gc.ca'
			nc_out.history = 'Created ' + tt.ctime(tt.time())

			#Create dimensions
			time = nc_out.createDimension('time', None) #use date2 for this
			level = nc_out.createDimension('level', max_level) 

			#Create coordinate variables
			times = nc_out.createVariable('time', np.float64, ('time',))
			levels = nc_out.createVariable('level', np.int32, ('level',))

			#Create 1D variables
			latitudes = nc_out.createVariable('latitude', np.float32, ('time'), zlib=True)
			longitudes = nc_out.createVariable('longitude', np.float32, ('time'), zlib=True)
			purpose_codes = nc_out.createVariable('purpose_code', str, ('time'), zlib=True)
			bottom_depths = nc_out.createVariable('station_bottom_depth', str, ('time'), zlib=True)
			gear_types = nc_out.createVariable('gear_type', str, ('time'), zlib=True)
			cruise_ids = nc_out.createVariable('cruise_id', str, ('time'), zlib=True)

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

			#Create a time filter
			filt = np.argsort(time_1D)

			#Fill in the 1D variables 
			latitudes[:] = latitude[filt]
			longitudes[:] = longitude[filt]
			times[:] = time_1D[filt]
			purpose_codes[:] = purpose_code[filt]
			bottom_depths[:] = bottom_depth[filt]
			gear_types[:] = gear_type[filt]
			cruise_ids[:] = cruise_ID[filt]

			#Fill 2D structure
			temp[:,:] = temperature[filt]
			saln[:,:] = salinity[filt]

			#Convert to time stamps
			time_stamps = [pd.Timestamp(i).to_pydatetime() for i in time_1D]
			times[:] = nc.date2num(time_stamps, units=times.units, calendar=times.calendar)
			levels[:] = np.arange(max_level)

			#Save and close the .nc file
			nc_out.close()
			#print(str(year)+', '+cruise+' done.')

	#Cycle through each of the years and import
	for year in np.unique(years):

		#Determine which files are present
		nc_files = []
		for i in os.listdir(file_output):
			if i.startswith(str(year)):
				nc_files.append(i)

		#Import the nc files of interest
		ds_merged = {}
		for i in nc_files:
			ds_merged[i.split('.')[0]] = xr.open_dataset(file_output+i)

		#Merge all of the sources together
		ds_merged = xr.concat([ds_merged[i] for i in ds_merged],dim='time',combine_attrs='override')
		ds_merged = ds_merged.sortby('time')

		#Save and close the file
		ds_merged.to_netcdf(file_yearly_output + str(year) + '.nc',
			encoding={'time':{'units': "seconds since 1900-01-01 00:00:00"}})
		ds_merged.close()
		#print(str(year)+' done!')