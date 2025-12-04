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

This package will include all Argo Float QA/QC functions.

'''


def download_Argo(
	file_output,
	years,
	):
	'''
	2000 onwards
	'''
	#Cycle through each year and month provided
	for year in years:

		#Argo
		#Format the time of the address
		year_1 = year
		year_2 = str(int(year))
		#Cycle through each month, otherwise download may time out
		months = np.arange(1,12+1).astype(str)
		#months = np.array(['12'])
		for month in months:
			month_1 = month.zfill(2)
			if month_1 == '12':
				month_2 = '01'
				year_2 = str(int(year)+1)
			else:
				month_2 = str(int(month)+1).zfill(2)
			address = 'https://erddap.ifremer.fr/erddap/tabledap/ArgoFloats.nc?fileNumber%2Cdata_type%2Cformat_version%2Chandbook_version%2Creference_date_time%2Cdate_creation%2Cdate_update%2Cplatform_number%2Cproject_name%2Cpi_name%2Ccycle_number%2Cdirection%2Cdata_center%2Cdc_reference%2Cdata_state_indicator%2Cdata_mode%2Cplatform_type%2Cfloat_serial_no%2Cfirmware_version%2Cwmo_inst_type%2Ctime%2Ctime_qc%2Ctime_location%2Clatitude%2Clongitude%2Cposition_qc%2Cpositioning_system%2Cprofile_pres_qc%2Cprofile_temp_qc%2Cprofile_psal_qc%2Cvertical_sampling_scheme%2Cconfig_mission_number%2Cpres%2Cpres_qc%2Cpres_adjusted%2Cpres_adjusted_qc%2Cpres_adjusted_error%2Ctemp%2Ctemp_qc%2Ctemp_adjusted%2Ctemp_adjusted_qc%2Ctemp_adjusted_error%2Cpsal%2Cpsal_qc%2Cpsal_adjusted%2Cpsal_adjusted_qc%2Cpsal_adjusted_error&time%3E='+year_1+'-'+month_1+'-01T00%3A00%3A00Z&time%3C='+year_2+'-'+month_2+'-01T00%3A00%3A00Z&latitude%3E=35&latitude%3C=80&longitude%3E=-100&longitude%3C=-42'

			#Download the file
			try:
				urllib.request.urlretrieve(address,file_output+'Argo_'+year+'_'+month_1+'.nc')
				print(year+'-'+month_1)
			except urllib.error.HTTPError:
				print(year+'-'+month_1+', download failed')


'''
TESTING AREA
file_input = '/gpfs/fs7/dfo/dpnm/joc000/Data/CASTS/Data_Input/Argo_Float/data_raw/'
file_output = '/gpfs/fs7/dfo/dpnm/joc000/Data/CASTS/Data_Input/Argo_Float/data_processed/'
years = np.arange(2000,2024+1).astype(str)
max_depth = 5000
'''


def merge_netcdf(
	file_input,
	file_output,
	years,
	max_depth=5000,
	):

	#Cycle through each year
	for year in years:

		#Define the names of the relevant files
		files = np.array(glob.glob(file_input+'Argo_'+year+'_*.nc'))
		files = np.sort(files)

		#Check that files are present
		if files.size > 0:

			#Import all the files of interest
			ds = {}
			var_1D_dict = {}
			var_2D_dict = {}
			for file in files:
				title = file.split('/')[-1][:-3]
				ds[title] = xr.open_dataset(file)

				#Isolate all the casts
				time_lat_lon = np.array([ds[title].time.values.astype(str),
					ds[title].latitude.values.astype(str),
					ds[title].longitude.values.astype(str)]).T
				time_lat_lon = np.array([','.join(row) for row in time_lat_lon])
				n = np.unique(time_lat_lon)

				#Start with the 1D variables
				var_1D = [
				'latitude',
				'longitude',
				'time',
				'project_name',
				'platform_number',
				'dc_reference',
				'platform_type',
				'data_mode',
				]
				var_1D_dict[title] = {}
				for var in var_1D:
					#Check to see if var is present
					if np.isin(var,list(ds[title].variables)) == True:
						var_1D_dict[title][var] = [list(ds[title][var].values[time_lat_lon==i]) for i in n]
					else:
						var_1D_dict[title][var] = [list(np.tile('',np.sum(time_lat_lon==i))) for i in n]

					#Ensure that 1D variable output is constant throughout all depths
					for i,value in enumerate(var_1D_dict[title][var]):
						if np.unique(value).size > 1:
							print(var+', profile #'+str(i)+' has more than one unique value in one profile.')
					var_1D_dict[title][var] = np.array([i[0] for i in var_1D_dict[title][var]])
				var_1D_dict[title]['file_name'] = np.tile(title, n.shape)
				var_1D_dict[title]['platform_type'][var_1D_dict[title]['platform_type'] == ''] = 'ARGO: TYPE UNKNOWN'

				#Process the 2D variables
				var_2D = ['temp','psal']
				var_2D_dict[title] = {}
				for var in var_2D:
					var_2D_dict[title][var] = [list(ds[title][var].values[time_lat_lon==i]) for i in n]
					depth = [list(ds[title]['pres'].values[time_lat_lon==i]) for i in n]
					var_2D_dict[title][var] = np.array([
					stats.binned_statistic(
						value,
						var_2D_dict[title][var][i],
						'mean',
						bins=np.arange(max_depth+1)
						).statistic for i,value in enumerate(depth)])


			#Combine all the files together, 1D first
			file_titles = list(ds.keys())
			var_1D_tog = {}
			for var in var_1D_dict[title].keys():
				var_1D_tog[var] = np.concatenate([var_1D_dict[i][var] for i in file_titles])

			#Now 2D
			var_2D_tog = {}
			for var in var_2D:
				var_2D_tog[var] = np.concatenate([var_2D_dict[i][var] for i in file_titles])

			#Create a time sort
			time_filt = np.argsort(var_1D_tog['time'])

			#Set up the .nc file
			nc_out = nc.Dataset(file_output+year+'.nc','w')
			
			#File information
			nc_out.Conventions = 'CF-1.6' #Ask Fred what this means
			nc_out.title = 'Argo Floats Yearly Netcdf File' #Temporary title for the .nc file
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
			file_names = nc_out.createVariable('file_names', str, ('time'), zlib=True)
			project_names = nc_out.createVariable('project_name', str, ('time'), zlib=True)
			platform_numbers = nc_out.createVariable('platform_number', str, ('time'), zlib=True)
			dc_references = nc_out.createVariable('dc_reference', str, ('time'), zlib=True)
			platform_types = nc_out.createVariable('platform_type', str, ('time'), zlib=True)
			data_modes = nc_out.createVariable('data_mode', str, ('time'), zlib=True)

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
			latitudes[:] = var_1D_tog['latitude'][time_filt]
			longitudes[:] = var_1D_tog['longitude'][time_filt]
			file_names[:] = var_1D_tog['file_name'][time_filt]
			times[:] = var_1D_tog['time'][time_filt]
			project_names[:] = var_1D_tog['project_name'][time_filt]
			platform_numbers[:] = var_1D_tog['platform_number'][time_filt]
			dc_references[:] = var_1D_tog['dc_reference'][time_filt]
			platform_types[:] = var_1D_tog['platform_type'][time_filt]
			data_modes[:] = var_1D_tog['data_mode'][time_filt]

			#Fill 2D structure
			temp[:,:] = var_2D_tog['temp'][time_filt]
			saln[:,:] = var_2D_tog['psal'][time_filt]

			#Convert to time stamps
			time_stamps = [pd.Timestamp(i).to_pydatetime() for i in var_1D_tog['time'][time_filt].astype('datetime64[m]')]
			times[:] = nc.date2num(time_stamps, units=times.units, calendar=times.calendar)
			levels[:] = np.arange(max_depth)

			#Save and close the .nc file
			nc_out.close()
			print(year+' done.')


