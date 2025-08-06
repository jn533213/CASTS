import glob
import os
import sys
import re
import warnings
import urllib
from urllib.error import HTTPError
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

This package will include all CIOOS-NAFC QA/QC functions.

'''


def download_AZMP(
	file_output,
	years,
	):
	'''
	1999 onwards
	'''
	#Cycle through each year and month provided
	for year in years:

		#AZMP
		#Format the time of the address
		year_1 = year
		year_2 = str(int(year)+1)
		address = 'https://cioosatlantic.ca/erddap/tabledap/nafc_azmp_ctd_profiles.nc?id%2Cplatform_name%2Cwmo_platform_code%2Cplatform_call_sign%2Cdfo_nafc_platform_code%2Cdfo_nafc_platform_name%2Cstation%2Clatitude%2Clongitude%2Ctime%2Csounder_depth%2Cinstrument%2Cset_number%2Ccast_type%2Ccomment%2Ccast_direction%2Csub_interval%2Cfishing_strata%2Ccloud_coverage%2Cwind_dir%2Cwind_speed_knots%2Cair_pressure%2Cair_dry_temp_celsius%2Cair_wet_temp_celsius%2Cwaves_period%2Cwaves_height%2Cswell_dir%2Cswell_height%2Cswell_period%2Cice_conc%2Cice_stage%2Cice_bergs%2Cice_sandt%2Cscan%2Cdepth%2CPRESPR01%2CTEMPS901%2CCNDCST01%2CPSALST01%2CSIGTEQST%2CCPHLPR01%2CDOXYZZ01%2CDWIRRXUD%2CPHXXZZXX%2COPTCPS01%2CATTNZS01%2CCCOMD002&time%3E='+year_1+'-01-01T00%3A00%3A00Z&time%3C='+year_2+'-01-01T00%3A00%3A00Z'

		#Download the file
		urllib.request.urlretrieve(address,file_output+'CIOOS_NAFC_AZMP_'+year+'.nc')

def download_multispecies(
	file_output,
	years,
	):
	'''
	1995 onwards
	'''
	#Cycle through each year and month provided
	for year in years:

		#Multispecies
		#Format the time of the address
		year_1 = year
		year_2 = str(int(year)+1)
		address = 'https://cioosatlantic.ca/erddap/tabledap/nafc_multispecies_ctd_profiles.nc?id%2Cplatform_name%2Cwmo_platform_code%2Cplatform_call_sign%2Cdfo_nafc_platform_code%2Cdfo_nafc_platform_name%2Cstation%2Clatitude%2Clongitude%2Ctime%2Csounder_depth%2Cinstrument%2Cset_number%2Ccast_type%2Ccomment%2Ccast_direction%2Csub_interval%2Cfishing_strata%2Ccloud_coverage%2Cwind_dir%2Cwind_speed_knots%2Cair_pressure%2Cair_dry_temp_celsius%2Cair_wet_temp_celsius%2Cwaves_period%2Cwaves_height%2Cswell_dir%2Cswell_height%2Cswell_period%2Cice_conc%2Cice_stage%2Cice_bergs%2Cice_sandt%2Cscan%2Cdepth%2CPRESPR01%2CTEMPS901%2CCNDCST01%2CPSALST01%2CSIGTEQST%2CCPHLPR01%2CDOXYZZ01%2CDWIRRXUD%2CPHXXZZXX%2COPTCPS01%2CATTNZS01%2CCCOMD002&time%3E='+year_1+'-01-01T00%3A00%3A00Z&time%3C='+year_2+'-01-01T00%3A00%3A00Z'

		#Download the file - account for website timeout
		urllib.request.urlretrieve(address,file_output+'CIOOS_NAFC_multispecies_'+year+'.nc')

def download_stn27(
	file_output,
	years,
	):
	'''
	1983 onwards
	'''
	#Cycle through each year and month provided
	for year in years:

		#STN-27
		#Format the time of the address
		year_1 = year
		year_2 = str(int(year)+1)
		address = 'https://cioosatlantic.ca/erddap/tabledap/nafc_station_27_ctd_profiles.nc?id%2Cplatform_name%2Cwmo_platform_code%2Cplatform_call_sign%2Cdfo_nafc_platform_code%2Cdfo_nafc_platform_name%2Cstation%2Clatitude%2Clongitude%2Ctime%2Csounder_depth%2Cinstrument%2Cset_number%2Ccast_type%2Ccomment%2Ccast_direction%2Csub_interval%2Cfishing_strata%2Ccloud_coverage%2Cwind_dir%2Cwind_speed_knots%2Cair_pressure%2Cair_dry_temp_celsius%2Cair_wet_temp_celsius%2Cwaves_period%2Cwaves_height%2Cswell_dir%2Cswell_height%2Cswell_period%2Cice_conc%2Cice_stage%2Cice_bergs%2Cice_sandt%2Cscan%2Cdepth%2CPRESPR01%2CTEMPS901%2CCNDCST01%2CPSALST01%2CSIGTEQST%2CCPHLPR01%2CDOXYZZ01%2CDWIRRXUD%2CPHXXZZXX%2COPTCPS01%2CATTNZS01%2CCCOMD002&time%3E='+year_1+'-01-01T00%3A00%3A00Z&time%3C='+year_2+'-01-01T00%3A00%3A00Z'

		#Download the file
		urllib.request.urlretrieve(address,file_output+'CIOOS_NAFC_stn27_'+year+'.nc')

def download_unsorted(
	file_output,
	years,
	):
	'''
	1983 onwards
	'''
	#Cycle through each year and month provided
	for year in years:

		#Unsorted
		#Format the time of the address
		year_1 = year
		year_2 = str(int(year)+1)
		address = 'https://cioosatlantic.ca/erddap/tabledap/nafc_bulk_unsorted_ctd_profiles.nc?id%2Cplatform_name%2Cwmo_platform_code%2Cplatform_call_sign%2Cdfo_nafc_platform_code%2Cdfo_nafc_platform_name%2Cstation%2Clatitude%2Clongitude%2Ctime%2Csounder_depth%2Cinstrument%2Cset_number%2Ccast_type%2Ccomment%2Ccast_direction%2Csub_interval%2Cfishing_strata%2Ccloud_coverage%2Cwind_dir%2Cwind_speed_knots%2Cair_pressure%2Cair_dry_temp_celsius%2Cair_wet_temp_celsius%2Cwaves_period%2Cwaves_height%2Cswell_dir%2Cswell_height%2Cswell_period%2Cice_conc%2Cice_stage%2Cice_bergs%2Cice_sandt%2Cscan%2Cdepth%2CPRESPR01%2CTEMPS901%2CCNDCST01%2CPSALST01%2CSIGTEQST%2CCPHLPR01%2CDOXYZZ01%2CDWIRRXUD%2CPHXXZZXX%2COPTCPS01%2CATTNZS01%2CCCOMD002&time%3E='+year_1+'-01-01T00%3A00%3A00Z&time%3C='+year_2+'-01-01T00%3A00%3A00Z'

		#Download the file
		urllib.request.urlretrieve(address,file_output+'CIOOS_NAFC_unsorted_'+year+'.nc')

def download_NSRF(
	file_output,
	years,
	):
	'''
	2005 onwards
	'''
	#Cycle through each year and month provided
	for year in years:

		#Unsorted
		#Format the time of the address
		year_1 = year
		year_2 = str(int(year)+1)
		address = 'https://cioosatlantic.ca/erddap/tabledap/nafc_nsrf_ctd_profiles.nc?id%2Cplatform_name%2Cwmo_platform_code%2Cplatform_call_sign%2Cdfo_nafc_platform_code%2Cdfo_nafc_platform_name%2Cstation%2Clatitude%2Clongitude%2Ctime%2Csounder_depth%2Cinstrument%2Cset_number%2Ccast_type%2Ccomment%2Ccast_direction%2Csub_interval%2Cfishing_strata%2Ccloud_coverage%2Cwind_dir%2Cwind_speed_knots%2Cair_pressure%2Cair_dry_temp_celsius%2Cair_wet_temp_celsius%2Cwaves_period%2Cwaves_height%2Cswell_dir%2Cswell_height%2Cswell_period%2Cice_conc%2Cice_stage%2Cice_bergs%2Cice_sandt%2Cscan%2Cdepth%2CPRESPR01%2CTEMPS901%2CCNDCST01%2CPSALST01%2CSIGTEQST%2CCPHLPR01%2CDOXYZZ01%2CDWIRRXUD%2CPHXXZZXX%2COPTCPS01%2CATTNZS01%2CCCOMD002&time%3E='+year_1+'-01-01T00%3A00%3A00Z&time%3C='+year_2+'-01-01T00%3A00%3A00Z'

		#Download the file
		urllib.request.urlretrieve(address,file_output+'CIOOS_NAFC_NSRF_'+year+'.nc')



'''
TESTING AREA
file_input = '/gpfs/fs7/dfo/dpnm/joc000/Data/CASTS/Data_Input/CIOOS_NAFC/data_raw/'
file_output = '/gpfs/fs7/dfo/dpnm/joc000/Data/CASTS/Data_Input/CIOOS_NAFC/data_processed/'
years = np.arange(1983,2023+1).astype(str)
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
		files = np.array(glob.glob(file_input+'*_'+year+'.nc'))
		files = files[np.array([i.endswith('.nc') for i in files])]

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
				n = np.unique(ds[title]['id'].values)

				#Start with the 1D variables
				var_1D = [
				'latitude',
				'longitude',
				'time',
				'instrument',
				'platform_name',
				'set_number',
				'station',
				'comment',
				'id',
				]
				var_1D_dict[title] = {}
				for var in var_1D:
					var_1D_dict[title][var] = [list(ds[title][var].values[ds[title]['id'].values==i]) for i in n]
					
					#Ensure that 1D variable output is constant throughout all depths
					for i,value in enumerate(var_1D_dict[title][var]):
						if np.unique(value).size > 1:
							print(var+', profile #'+str(i)+' has more than one unique value in one profile.')
					var_1D_dict[title][var] = np.array([i[0] for i in var_1D_dict[title][var]])
				var_1D_dict[title]['file_name'] = np.tile(title, n.shape)

				#Process the 2D variables
				var_2D = ['TEMPS901','PSALST01']
				var_2D_dict[title] = {}
				for var in var_2D:
					var_2D_dict[title][var] = [list(ds[title][var].values[ds[title]['id'].values==i]) for i in n]
					depth = [list(ds[title]['PRESPR01'].values[ds[title]['id'].values==i]) for i in n]
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
			nc_out.title = 'CIOOS_NAFC Yearly Netcdf File' #Temporary title for the .nc file
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
			platform_names = nc_out.createVariable('platform_name', str, ('time'), zlib=True)
			instruments = nc_out.createVariable('instrument', str, ('time'), zlib=True)
			set_numbers = nc_out.createVariable('set_number', str, ('time'), zlib=True)
			stations = nc_out.createVariable('station', str, ('time'), zlib=True)
			comments = nc_out.createVariable('comment', str, ('time'), zlib=True)
			IDs = nc_out.createVariable('ID', str, ('time'), zlib=True)

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
			instruments[:] = var_1D_tog['instrument'][time_filt]
			platform_names[:] = var_1D_tog['platform_name'][time_filt]
			set_numbers[:] = var_1D_tog['set_number'][time_filt].astype(str)
			stations[:] = var_1D_tog['station'][time_filt].astype(str)
			comments[:] = var_1D_tog['comment'][time_filt]
			IDs[:] = var_1D_tog['id'][time_filt]

			#Fill 2D structure
			temp[:,:] = var_2D_tog['TEMPS901'][time_filt]
			saln[:,:] = var_2D_tog['PSALST01'][time_filt]

			#Convert to time stamps
			time_stamps = [pd.Timestamp(i).to_pydatetime() for i in var_1D_tog['time'][time_filt].astype('datetime64[m]')]
			times[:] = nc.date2num(time_stamps, units=times.units, calendar=times.calendar)
			levels[:] = np.arange(max_depth)

			#Save and close the .nc file
			nc_out.close()
			print(year+' done.')

