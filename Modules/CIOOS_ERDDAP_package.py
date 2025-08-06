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

This package will include all CIOOS-ERDDAP QA/QC functions.

'''

def download_AZMP(
	file_output,
	years,
	):
	'''
	1997 onwards
	'''
	#Cycle through each year and month provided
	for year in years:

		#AZMP
		#Format the time of the address
		year_1 = year
		year_2 = str(int(year)+1)
		address = 'https://cioosatlantic.ca/erddap/tabledap/bio_atlantic_zone_monitoring_program_ctd.nc?project%2Cplatform_name%2Cplatform_id%2Cchief_scientist%2Ccruise_name%2Ccruise_number%2Cgeographic_area%2Cstation%2Cid%2Cevent_number%2Cprofile_direction%2Clatitude%2Clongitude%2CpreciseLat%2CpreciseLon%2CQLATD_01%2CQLOND_01%2Ctime%2Cprofile_start_time%2Cprofile_end_time%2Cmeasurement_time%2Cdepth%2CAHSFZZ01%2CQALTB_01%2CPRESPR01%2CQPRES_01%2CPSLTZZ01%2CPSALST01%2CQPSAL_01%2CPSLTZZ02%2CPSALST02%2CQPSAL_02%2CTEMPP681%2CTEMPP901%2CTEMPS601%2CTEMPS901%2CTEMPPR01%2CQTEMP_01%2CTEMPP682%2CTEMPP902%2CTEMPS602%2CTEMPS902%2CTEMPPR02%2CQTEMP_02%2CCNDCST01%2CQCNDC_01%2CCNDCST02%2CQCNDC_02%2CDOXYZZ01%2CQDOXY_01%2CDOXYZZ02%2CQDOXY_02%2CCDOMZZ01%2CQCDOMZZ01%2CCDOMZZ02%2CQCDOMZZ02%2CCPHLPR01%2CQCPHLPR01%2CCPHLPR02%2CQCPHLPR02%2CPHXXZZ01%2CQPHPH_01%2CIRRDSV01%2CQSPAR_01%2CIRRDUV01%2CQPSAR_01%2CVSCTXX01%2CQVSCTXX01%2CATTNZS01%2CQATTNZS01%2COPTCPS01%2CQOPTCPS01%2CScanNumber%2CRecPerBin&time%3E='+year_1+'-01-01T00%3A00%3A00Z&time%3C='+year_2+'-01-01T00%3A00%3A00Z'

		#Download the file
		urllib.request.urlretrieve(address,file_output+'CIOOS_ERDDAP_AZMP_'+year+'.nc')


def download_AZOMP(
	file_output,
	years,
	):
	'''
	1993 onwards
	'''
	#Cycle through each year and month provided
	for year in years:

		#AZOMP
		#Format the time of the address
		year_1 = year
		year_2 = str(int(year)+1)
		address = 'https://cioosatlantic.ca/erddap/tabledap/bio_atlantic_zone_off_shelf_monitoring_program_ctd.nc?project%2Cplatform_name%2Cplatform_id%2Cchief_scientist%2Ccruise_name%2Ccruise_number%2Cgeographic_area%2Cstation%2Cid%2Cevent_number%2Cprofile_direction%2Clatitude%2Clongitude%2CpreciseLat%2CpreciseLon%2CQLATD_01%2CQLOND_01%2Ctime%2Cprofile_start_time%2Cprofile_end_time%2Cmeasurement_time%2Cdepth%2CAHSFZZ01%2CQALTB_01%2CPRESPR01%2CQPRES_01%2CPSLTZZ01%2CPSALST01%2CQPSAL_01%2CPSLTZZ02%2CPSALST02%2CQPSAL_02%2CTEMPP681%2CTEMPP901%2CTEMPS601%2CTEMPS901%2CTEMPPR01%2CQTEMP_01%2CTEMPP682%2CTEMPP902%2CTEMPS602%2CTEMPS902%2CTEMPPR02%2CQTEMP_02%2CCNDCST01%2CQCNDC_01%2CCNDCST02%2CQCNDC_02%2CDOXYZZ01%2CQDOXY_01%2CDOXYZZ02%2CQDOXY_02%2CCDOMZZ01%2CQCDOMZZ01%2CCPHLPR01%2CQCPHLPR01%2CCPHLPR02%2CQCPHLPR02%2CPHXXZZ01%2CQPHPH_01%2CIRRDSV01%2CQSPAR_01%2CIRRDUV01%2CQPSAR_01%2CVSCTXX01%2CQVSCTXX01%2CATTNZS01%2CQATTNZS01%2COPTCPS01%2CQOPTCPS01%2CScanNumber%2CRecPerBin&time%3E='+year_1+'-01-01T00%3A00%3A00Z&time%3C='+year_2+'-01-01T00%3A00%3A00Z'

		#Download the file
		urllib.request.urlretrieve(address,file_output+'CIOOS_ERDDAP_AZOMP_'+year+'.nc')


def download_ecosystems(
	file_output,
	years,
	):
	'''
	1996 onwards
	'''
	#Cycle through each year and month provided
	for year in years:

		#Ecosystems
		#Format the time of the address
		year_1 = year
		year_2 = str(int(year)+1)
		address = 'https://cioosatlantic.ca/erddap/tabledap/bio_maritimes_region_ecosystem_survey_ctd.nc?project%2Cplatform_name%2Cplatform_id%2Cchief_scientist%2Ccruise_name%2Ccruise_number%2Cgeographic_area%2Cstation%2Cid%2Cevent_number%2Cprofile_direction%2Clatitude%2Clongitude%2CpreciseLat%2CpreciseLon%2CQLATD_01%2CQLOND_01%2Ctime%2Cprofile_start_time%2Cprofile_end_time%2Cmeasurement_time%2Cdepth%2CAHSFZZ01%2CQALTB_01%2CPRESPR01%2CQPRES_01%2CPSLTZZ01%2CPSALST01%2CQPSAL_01%2CPSLTZZ02%2CPSALST02%2CQPSAL_02%2CTEMPP681%2CTEMPP901%2CTEMPS601%2CTEMPS901%2CTEMPPR01%2CQTEMP_01%2CTEMPP682%2CTEMPP902%2CTEMPS602%2CTEMPS902%2CTEMPPR02%2CQTEMP_02%2CCNDCST01%2CQCNDC_01%2CCNDCST02%2CQCNDC_02%2CDOXYZZ01%2CQDOXY_01%2CDOXYZZ02%2CQDOXY_02%2CCDOMZZ01%2CQCDOMZZ01%2CCDOMZZ02%2CQCDOMZZ02%2CCPHLPR01%2CQCPHLPR01%2CCPHLPR02%2CQCPHLPR02%2CPHXXZZ01%2CQPHPH_01%2CIRRDSV01%2CQSPAR_01%2CIRRDUV01%2CQPSAR_01%2CVSCTXX01%2CQVSCTXX01%2CATTNZS01%2CQATTNZS01%2COPTCPS01%2CQOPTCPS01%2CScanNumber%2CRecPerBin&time%3E='+year_1+'-01-01T00%3A00%3A00Z&time%3C='+year_2+'-01-01T00%3A00%3A00Z'

		#Download the file
		urllib.request.urlretrieve(address,file_output+'CIOOS_ERDDAP_ecosystems_'+year+'.nc')


def download_historical_offshore_international(
	file_output,
	years,
	):
	'''
	1969 to 2023
	'''
	#Cycle through each year and month provided
	for year in years:

		#Format the time of the address
		year_1 = year
		year_2 = str(int(year)+1)
		address = 'https://cioosatlantic.ca/erddap/tabledap/bio_historical_offshore_international_ctd.nc?project%2Cplatform_name%2Cplatform_id%2Cchief_scientist%2Ccruise_name%2Ccruise_number%2Cgeographic_area%2Cstation%2Cid%2Cevent_number%2Cprofile_direction%2Clatitude%2Clongitude%2CpreciseLat%2CpreciseLon%2Ctime%2Cprofile_start_time%2Cprofile_end_time%2Cmeasurement_time%2Cdepth%2CAHSFZZ01%2CQALTB_01%2CPRESPR01%2CQPRES_01%2CPSLTZZ01%2CPSALST01%2CQPSAL_01%2CTEMPP681%2CTEMPP901%2CTEMPS601%2CTEMPS901%2CTEMPPR01%2CQTEMP_01%2CCNDCST01%2CQCNDC_01%2CDOXYZZ01%2CQDOXY_01%2CCDOMZZ01%2CQCDOMZZ01%2CCPHLPR01%2CQCPHLPR01%2CCPHLPR02%2CQCPHLPR02%2CPHXXZZ01%2CQPHPH_01%2CIRRDUV01%2CQPSAR_01%2CScanNumber%2CRecPerBin&time%3E='+year_1+'-01-01T00%3A00%3A00Z&time%3C='+year_2+'-01-01T00%3A00%3A00Z'

		#Download the file
		urllib.request.urlretrieve(address,file_output+'CIOOS_ERDDAP_historical_offshore_international_'+year+'.nc')


def download_historical_west_greenland(
	file_output,
	years,
	):
	'''
	1989 to 1994
	'''
	#Cycle through each year and month provided
	for year in years:

		#Format the time of the address
		year_1 = year
		year_2 = str(int(year)+1)
		address = 'https://cioosatlantic.ca/erddap/tabledap/bio_historical_west_greenland_ctd.nc?platform_name%2Cchief_scientist%2Ccruise_name%2Ccruise_number%2Cgeographic_area%2Cstation%2Cid%2Cevent_number%2Cprofile_direction%2Clatitude%2Clongitude%2Ctime%2Cprofile_start_time%2Cdepth%2CPRESPR01%2CPSLTZZ01%2CTEMPPR01&time%3E='+year_1+'-01-01T00%3A00%3A00Z&time%3C='+year_2+'-01-01T00%3A00%3A00Z'

		#Download the file
		urllib.request.urlretrieve(address,file_output+'CIOOS_ERDDAP_historical_west_greenland_'+year+'.nc')


def download_historical_coastal_rosette(
	file_output,
	years,
	):
	'''
	1969 to 2011
	'''
	#Cycle through each year and month provided
	for year in years:

		#Format the time of the address
		year_1 = year
		year_2 = str(int(year)+1)
		address = 'https://cioosatlantic.ca/erddap/tabledap/bio_historical_coastal_ctd.nc?project%2Cplatform_name%2Cplatform_id%2Cchief_scientist%2Ccruise_name%2Ccruise_number%2Cgeographic_area%2Cstation%2Cid%2Cevent_number%2Cprofile_direction%2Clatitude%2Clongitude%2Ctime%2Cprofile_start_time%2Cprofile_end_time%2Cmeasurement_time%2Cdepth%2CPRESPR01%2CPSLTZZ01%2CPSALST01%2CTEMPPR01%2CCNDCST01%2CDOXYZZ01%2CDOXYZZ02%2CCPHLPR01%2COPTCPS01%2CScanNumber&time%3E='+year_1+'-01-01T00%3A00%3A00Z&time%3C='+year_2+'-01-01T00%3A00%3A00Z'

		#Download the file
		urllib.request.urlretrieve(address,file_output+'CIOOS_ERDDAP_historical_coastal_rosette_'+year+'.nc')
		'''
		#Format the time of the address
		year_1 = year
		year_2 = str(int(year))
		#Cycle through each month, otherwise download may time out
		months = np.arange(1,12+1).astype(str)
		for month in months:
			month_1 = month.zfill(2)
			if month_1 == '12':
				month_2 = '01'
				year_2 = str(int(year)+1)
			else:
				month_2 = str(int(month)+1).zfill(2)
			address = 'https://cioosatlantic.ca/erddap/tabledap/bio_historical_coastal_ctd.nc?project%2Cplatform_name%2Cplatform_id%2Cchief_scientist%2Ccruise_name%2Ccruise_number%2Cgeographic_area%2Cstation%2Cid%2Cevent_number%2Cprofile_direction%2Clatitude%2Clongitude%2Ctime%2Cprofile_start_time%2Cprofile_end_time%2Cmeasurement_time%2Cdepth%2CPRESPR01%2CPSLTZZ01%2CPSALST01%2CTEMPPR01%2CCNDCST01%2CDOXYZZ01%2CDOXYZZ02%2CCPHLPR01%2COPTCPS01%2CScanNumber&time%3E='+year_1+'-'+month_1+'-01T00%3A00%3A00Z&time%3C='+year_2+'-'+month_2+'-01T00%3A00%3A00Z'

			#Download the file
			try:
				urllib.request.urlretrieve(address,file_output+'CIOOS_ERDDAP_historical_coastal_rosette_'+year+'_'+month_1+'.nc')
				print(year+'-'+month_1)
			except urllib.error.HTTPError:
				print(year+'-'+month_1+', download failed')
		'''



def download_historical_arctic_rosette(
	file_output,
	years,
	):
	'''
	1973 to 2013
	'''
	#Cycle through each year and month provided
	for year in years:

		#Format the time of the address
		year_1 = year
		year_2 = str(int(year)+1)
		address = 'https://cioosatlantic.ca/erddap/tabledap/bio_historical_arctic_ctd.nc?project%2Cplatform_name%2Cplatform_id%2Cchief_scientist%2Ccruise_name%2Ccruise_number%2Cgeographic_area%2Cstation%2Cid%2Cevent_number%2Cprofile_direction%2Clatitude%2Clongitude%2Ctime%2Cprofile_start_time%2Cprofile_end_time%2Cmeasurement_time%2Cdepth%2CPRESPR01%2CQPRES_01%2CPSLTZZ01%2CPSALST01%2CQPSAL_01%2CTEMPP681%2CTEMPP901%2CTEMPS601%2CTEMPS901%2CTEMPPR01%2CQTEMP_01%2CCNDCST01%2CQCNDC_01%2CDOXYZZ01%2CScanNumber%2CRecPerBin&time%3E='+year_1+'-01-01T00%3A00%3A00Z&time%3C='+year_2+'-01-01T00%3A00%3A00Z'
		#Download the file
		urllib.request.urlretrieve(address,file_output+'CIOOS_ERDDAP_historical_arctic_rosette_'+year+'.nc')


def download_barrow_strait_flow_observational(
	file_output,
	years,
	):
	'''
	1998 to 2010
	'''
	#Cycle through each year and month provided
	for year in years:

		#Format the time of the address
		year_1 = year
		year_2 = str(int(year)+1)
		address = 'https://cioosatlantic.ca/erddap/tabledap/bio_barrow_strait_program_ctd.nc?project%2Cplatform_name%2Cplatform_id%2Cchief_scientist%2Ccruise_name%2Ccruise_number%2Cgeographic_area%2Cstation%2Cid%2Cevent_number%2Cprofile_direction%2Clatitude%2Clongitude%2Ctime%2Cprofile_start_time%2Cprofile_end_time%2Cdepth%2CPRESPR01%2CPSLTZZ01%2CPSALST01%2CTEMPPR01%2CCNDCST01%2CScanNumber&time%3E='+year_1+'-01-01T00%3A00%3A00Z&time%3C='+year_2+'-01-01T00%3A00%3A00Z'
		#Download the file
		urllib.request.urlretrieve(address,file_output+'CIOOS_ERDDAP_barrow_strait_flow_observational_'+year+'.nc')


def download_annual_ice_forecast_hydrographic_survey(
	file_output,
	years,
	):
	'''
	1972 to 1995
	'''
	#Cycle through each year and month provided
	for year in years:

		#Format the time of the address
		year_1 = year
		year_2 = str(int(year)+1)
		address = 'https://cioosatlantic.ca/erddap/tabledap/bio_annual_ice_forecast_hydrographic_survey_st_lawrence_ctd.nc?platform_name%2Cplatform_id%2Cchief_scientist%2Ccruise_name%2Ccruise_number%2Cgeographic_area%2Cstation%2Cid%2Cevent_number%2Cprofile_direction%2Clatitude%2Clongitude%2Ctime%2Cprofile_start_time%2Cprofile_end_time%2Cmeasurement_time%2Cdepth%2CPRESPR01%2CPSLTZZ01%2CPSALST01%2CTEMPPR01%2CCNDCST01%2CDOXYZZ01%2CScanNumber&time%3E='+year_1+'-01-01T00%3A00%3A00Z&time%3C='+year_2+'-01-01T00%3A00%3A00Z'
		#Download the file
		urllib.request.urlretrieve(address,file_output+'CIOOS_ERDDAP_annual_ice_forecast_hydrographic_survey_'+year+'.nc')



'''
TESTING AREA
file_input = '/gpfs/fs7/dfo/dpnm/joc000/Data/CASTS/Data_Input/CIOOS_ERDDAP/data_raw/Mar_2025/'
file_output = '/gpfs/fs7/dfo/dpnm/joc000/Data/CASTS/Data_Input/CIOOS_ERDDAP/data_processed/'
years = np.arange(1969,2024+1).astype(str)
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
				'project',
				'platform_name',
				'platform_id',
				'instrument_ID',
				'cruise_name',
				'cruise_number',
				'station',
				'id',
				]
				var_1D_dict[title] = {}
				for var in var_1D:
					#Check to see if var is present
					if np.isin(var,list(ds[title].variables)) == True:
						var_1D_dict[title][var] = [list(ds[title][var].values[ds[title]['id'].values==i]) for i in n]
					elif var == 'instrument_ID':
						var_1D_dict[title][var] = [list(np.tile(ds[title].attrs['odf_data_type'],np.sum(ds[title]['id'].values==i))) for i in n]
					else:
						var_1D_dict[title][var] = [list(np.tile('',np.sum(ds[title]['id'].values==i))) for i in n]

					#Ensure that 1D variable output is constant throughout all depths
					for i,value in enumerate(var_1D_dict[title][var]):
						if np.unique(value).size > 1:
							print(var+', profile #'+str(i)+' has more than one unique value in one profile.')
					var_1D_dict[title][var] = np.array([i[0] for i in var_1D_dict[title][var]])
				var_1D_dict[title]['file_name'] = np.tile(title, n.shape)

				#Process the 2D variables
				var_2D = ['TEMPPR01','PSALST01']
				var_2D_dict[title] = {}
				for var in var_2D:
					if np.isin(var,list(ds[title].variables)) == True:
						var_2D_dict[title][var] = [list(ds[title][var].values[ds[title]['id'].values==i]) for i in n]
					else:
						var_2D_dict[title][var] = [list(ds[title]['PSLTZZ01'].values[ds[title]['id'].values==i]) for i in n]
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
			nc_out.title = 'CIOOS_ERDDAP Yearly Netcdf File' #Temporary title for the .nc file
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
			projects = nc_out.createVariable('project', str, ('time'), zlib=True)
			platform_names = nc_out.createVariable('platform_name', str, ('time'), zlib=True)
			platform_ids = nc_out.createVariable('platform_ids', str, ('time'), zlib=True)
			instrument_ids = nc_out.createVariable('instrument_ids', str, ('time'), zlib=True)
			cruise_names = nc_out.createVariable('cruise_name', str, ('time'), zlib=True)
			cruise_numbers = nc_out.createVariable('cruise_number', str, ('time'), zlib=True)
			stations = nc_out.createVariable('station', str, ('time'), zlib=True)
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
			projects[:] = var_1D_tog['project'][time_filt]
			platform_names[:] = var_1D_tog['platform_name'][time_filt]
			platform_ids[:] = var_1D_tog['platform_id'][time_filt]
			instrument_ids[:] = var_1D_tog['instrument_ID'][time_filt]
			cruise_names[:] = var_1D_tog['cruise_name'][time_filt]
			cruise_numbers[:] = var_1D_tog['cruise_number'][time_filt]
			stations[:] = var_1D_tog['station'][time_filt]
			IDs[:] = var_1D_tog['id'][time_filt]

			#Fill 2D structure
			temp[:,:] = var_2D_tog['TEMPPR01'][time_filt]
			saln[:,:] = var_2D_tog['PSALST01'][time_filt]

			#Convert to time stamps
			time_stamps = [pd.Timestamp(i).to_pydatetime() for i in var_1D_tog['time'][time_filt].astype('datetime64[m]')]
			times[:] = nc.date2num(time_stamps, units=times.units, calendar=times.calendar)
			levels[:] = np.arange(max_depth)

			#Save and close the .nc file
			nc_out.close()
			print(year+' done.')








