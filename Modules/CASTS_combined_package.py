import glob
import os
import sys
import re
import csv
import warnings
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

This script will combine all other sources, also contains QA/QC functions.

'''


'''
#TESTING AREA
path = {}
path['path1'] = directory+'Data_Input/BIO_Climate/NetCDF/BIO_Climate_Databases/' #BIO, 1913-2010
path['path2'] = directory+'Data_Input/BIO_Climate/NetCDF/database_2008-2017/' #BIO, 2008-2017
path['path3'] = directory+'Data_Input/BIO_Climate/NetCDF/database_2018/' #BIO, 2018
path['path4'] = directory+'Data_Input/NAFC_Oceanography/4_outliers/' #NAFC, 1999-2022
path['path5'] = directory+'Data_Input/CIOOS_NAFC/data_processed/' #CIOOS_NAFC, 1983-2023
path['path6'] = directory+'Data_Input/CIOOS_ERDDAP/data_processed/' #CIOOS-ERRDAP, 1996-2020
path['path7'] = directory+'Data_Input/NEFSC/yearly_netcdf_files/' #NEFSC, 1981-2021
path['path8'] = directory+'Data_Input/NAFC_Aquaculture/data_processed/' #CTD Andry, 2009-2020
path['path9'] = directory+'Data_Input/IML/data_processed/' # MLI-Sourced, 1972-2022
path['path10'] = directory+'Data_Input/NCEI_GTSPP/data_processed/bottom_flagged/' # NCEI-GTSPP, 1990-2019
path['path11'] = directory+'Data_Input/Polar_Data_Catalogue/data_product/netcdf_yearly/' # Polar Data Catalogue, 2002-2020
path['path12'] = directory+'Data_Input/IEO_Spain/netcdf_yearly/' # IEO Spain, 2019-2021
path['path13'] = directory+'Data_Input/IEO_Spain/netcdf_yearly/EU_NAFO/' # IEO Spain, 1988-2020
path['path14'] = directory+'Data_Input/Other/Freds_Files/netcdf_yearly/' #Marine Institute, 2000-2008
path['path15'] = directory+'Data_Input/WOD/data_processed/' #World Ocean Database, 1873-1910
path_source = {}
path_source['path1'] = 'Climate'
path_source['path2'] = 'BIO-OMM'
path_source['path3'] = 'BIO-OMM'
path_source['path4'] = 'NAFC-Oceanography'
path_source['path5'] = 'CIOOS_NAFC'
path_source['path6'] = 'CIOOS_BIO'
path_source['path7'] = 'NEFSC'
path_source['path8'] = 'NAFC-Aquaculture'
path_source['path9'] = 'MLI'
path_source['path10'] = 'NCEI'
path_source['path11'] = 'Polar-Data-Catalogue'
path_source['path12'] = 'EU-NAFO'
path_source['path13'] = 'EU-NAFO'
path_source['path14'] = 'Marine-Institute-NL'
path_source['path15'] = 'WOD'
years=np.arange(1873,2024+1).astype(str)
file_output=directory+'Data_Products/1_combined_raw/'
'''


def combine_netcdf(
	path,
	path_source,
	years,
	file_output,
	):


	#Cycle through each year
	for year in years[:]:

		#Import each of the files, if one exists
		ds = {}

		#Cycle through each path
		for i in path.keys():
			if os.path.isfile(path[i]+str(year)+'.nc') == True:
				ds[i] = xr.open_dataset(path[i]+str(year)+'.nc')

		#Determine if data is available
		if len(ds) == 0:
			continue

		#Some variables will need to be renamed for consistency
		for i in ds.keys():
			if i == 'path1':
				ds[i] = ds[i].rename({'datatype': 'instrument_ID',})
				ds[i] = ds[i].rename({'maximumdepth': 'sounder_depth'})
			if i == 'path2' or i == 'path3':
				ds[i] = ds[i].rename({'maximumdepth': 'sounder_depth'})
			if i == 'path3':
				ds[i] = ds[i].rename({'datatype': 'instrument_ID'})
			if i == 'path5':
				ds[i] = ds[i].rename({'instrument': 'instrument_ID'})
				ds[i] = ds[i].rename({'platform_name': 'trip_ID'})
				ds[i] = ds[i].rename({'comment': 'comments'})
			if i == 'path7':
				ds[i] = ds[i].rename({'station_bottom_depth': 'sounder_depth'})
				ds[i] = ds[i].rename({'gear_type': 'instrument_type'})
				ds[i] = ds[i].rename({'cruise_id': 'trip_ID'})
			if i == 'path10':
				ds[i] = ds[i].rename({'file_name': 'file_names'})
				ds[i] = ds[i].rename({'source_ID': 'comments'})
				ds[i] = ds[i].rename({'data_type': 'instrument_ID',})
			if i == 'path12':
				ds[i] = ds[i].rename({'datatype': 'instrument_ID'})
				ds[i] = ds[i].rename({'instrument_id': 'instrument_type'})
				ds[i] = ds[i].rename({'station_id': 'comments'})
			if i == 'path13':
				ds[i] = ds[i].rename({'Cruise': 'trip_ID'})
				ds[i] = ds[i].rename({'Station': 'comment'})
				ds[i] = ds[i].rename({'bottom_depth': 'sounder_depth'})
				ds[i] = ds[i].rename({'instrument_id': 'instrument_type'})
			if i == 'path14':
				ds[i] = ds[i].rename({'mission_id': 'trip_ID'})

		#Isolate for the variables of interest
		variables = np.array([
			'level',
			'latitude',
			'longitude',
			'instrument_ID',
			'instrument_type',
			'comments',
			'trip_ID',
			'file_names',
			'sounder_depth',
			'temperature',
			'salinity',
			])
		for i in ds.keys():

			#Check to see if the variable is present
			variables_temp = variables[np.isin(variables,list(ds[i].variables))]
			ds[i] = ds[i][variables_temp]

			#If the variable is not available make a nan one
			if np.isin(variables,list(ds[i].variables)).all() == False:
				for ii in variables[~np.isin(variables,list(ds[i].variables))]:
					ds[i][ii] = (['time'], np.full(ds[i].time.size, np.nan).astype(str))

		#Merge all of the variables together
		ds_merged = xr.concat([ds[i] for i in ds.keys()],dim='time',combine_attrs='override')
		ds_merged = ds_merged.rename({'comments': 'station_ID'})

		#Determine where each path is located
		path_location = np.concatenate([np.full(ds[i].time.size,path_source[i]) for i in ds.keys()])
		ds_merged['source'] = (['time'],path_location.astype(object))

		#Remove all casts which have a standard deviation of 0 throughout, if they have more than 2 measurements
		#Determine for both temperature and salinity
		warnings.simplefilter("ignore", category=RuntimeWarning)
		temp_std = ds_merged.temperature.std(axis=1).values == 0
		temp_nom = ds_merged.level.size - np.isnan(ds_merged['temperature']).sum(axis=1)
		saln_std = ds_merged.salinity.std(axis=1).values == 0
		saln_nom = ds_merged.level.size - np.isnan(ds_merged['salinity']).sum(axis=1)

		#Remove the constants
		ds_merged['temperature'][temp_std*(temp_nom > 3),:] = np.nan
		ds_merged['salinity'][saln_std*(saln_nom > 3),:] = np.nan

		#Remove the new empties
		temp_empties = np.isnan(ds_merged['temperature']).sum(axis=1) == ds_merged.level.size
		saln_empties = np.isnan(ds_merged['salinity']).sum(axis=1) == ds_merged.level.size
		ds_merged = ds_merged.sel(time = ~(temp_empties.values*saln_empties.values))

		#Define the time and distance radius around the point 
		time_cutoff = np.timedelta64(10,'m')
		dist_cutoff = 0.0125

		#Import the time, latitude, and longitude
		time = ds_merged.time.values
		latitude = ds_merged.latitude.values
		longitude = ds_merged.longitude.values

		#Create a variable marking duplicates
		duplicates_flag = np.full(ds_merged.time.shape,False) 

		#Organize instrument ID into categories
		instrument_call = {}
		instrument_call['CT'] = ['19899','19P-7','19P45','19P53','GL','OLABS','STD12','NBmun','0CTD','0CTD1','1CTD','1CTD1','CD','CT','CTD','CU','FAPCD','FAPCU','GLCTD','GLctd','MEDCD','S0000','S0002','S0036','S0038','S0078','S0097','S0103','S0118','S0129','S0212','S0277','S0281','S0344','S0394','S0454','S0466','S0580','S0582','S0583','S0674','S0688','S0689','S0690','S0691','S0692','S0758','S0803','S0845','S0846','S0847','S0910','S0911','S0912','S0935','S1003','S1019','S1020','S1021','S1098','S1101','S1145','S1146','S1221','S1237','S1238','S1247','S1257','S1308','S1309','S1310','S1311','S1312','S1313','S1314','S1315','S1316','S1317','S1318','S1319','S1460','S2245','S2246','S2247','S2248','S4016','S4017','S4018','S4019','S4020','S4578','S4579','S4580','S4581','S4582','S4776','S4777','S5117','S7849','S7853','S7854','S7855','S7992','S8013','S8014','S9999','SB', 'SB002', 'SB036', 'SB038','SB12','SB279','SB280','SB394','SB450','SB453','SB454','SB455','SB580','SB581','SB583','SB688','SB689','SB690','SB691','SB692','SB73','SB74','SB81','SB82','SB822','SB83','SB845','SB880','SB890','SB900','SB910','SB920','SB935','SBMUN']
		instrument_call['BO'] = ['BO','FAPBO','FAPbo','FRCBO']
		instrument_call['BF'] = ['BF','FAPBF']
		instrument_call['TO'] = ['TO','FAPTO']
		instrument_call['DT'] = ['DT','BT']
		instrument_call['MB'] = ['MB']
		instrument_call['XB'] = ['XB','XBT05','XBT06','XBT07','XBT10']
		instrument_call['PF'] = ['PF','B3']
		instrument_call[''] = ['nan','21573.0']

		#Define the new instrument_ID data
		ds_instrument_ID = ds_merged.instrument_ID.values
		new_instrument_ID = np.full(time.size,'',dtype='<U32')

		#Organize the measurements into categories
		for i,value in enumerate(ds_instrument_ID):
			for ii in instrument_call:
				if np.isin(value, instrument_call[ii]):
					new_instrument_ID[i] = ii
					break

		#Cycle through each point 
		for i in np.arange(time.size):

			#Determine if the cast is already flagged as a duplicate, if so, skip
			if duplicates_flag[i] == False:

				#Determine the distance and time between that point and every other
				time_diff = time - time[i]
				lon_diff = longitude - longitude[i]
				lat_diff = latitude - latitude[i]

				#Determine which points make the cutoffs
				time_diff = np.abs(time_diff) < time_cutoff
				lon_diff = np.abs(lon_diff) < dist_cutoff
				lat_diff = np.abs(lat_diff) < dist_cutoff

				#Determine where all three criteras are met
				flagged_casts = time_diff*lat_diff*lon_diff

				#Determine if there are duplicates present, if so proceed
				if flagged_casts.sum() > 1:

					#Determine the corresponding sources, index
					flagged_source = ds_merged.source[flagged_casts].values
					flagged_index = np.where(flagged_casts == True)[0]

					#Determine the number of temperature and salinity measurements in each
					temp = ds_merged.temperature[flagged_index].values
					saln = ds_merged.salinity[flagged_index].values
					#Determine the number of recordings in each cast
					nom_temp = (~np.isnan(temp)).sum(axis=1)
					nom_saln = (~np.isnan(saln)).sum(axis=1)

					#Record the instrument ID
					dup_instrID = new_instrument_ID[flagged_casts]
					no_id = np.unique(dup_instrID)
					no_id = no_id[no_id != ''] #remove instances of no instrument ID specified

					#Determine which instrument nans should included with
					if np.sum(dup_instrID == '').sum() == 0:
						nans_include = {}
						for ii in no_id:
							nans_include[ii] = []
					elif no_id.size == 0:
						nans_include = {'':flagged_index}
					elif no_id.size == 1:
						nans_include = {no_id[0]:flagged_index[dup_instrID == '']}
					else:
						#Determine where all the nans are
						nan_index = flagged_index[dup_instrID == '']
						nans_include = {}
						for ii in no_id:
							nans_include[ii] = []
						#Cycle through each
						for nan_id in nan_index:
							#Determine the number of measurements where nans are
							diff = np.zeros(no_id.size)
							for ii,grp in enumerate(no_id):
								Tdiff = np.abs(np.mean(nom_temp[dup_instrID == grp]) - nom_temp[flagged_index == nan_id])
								Sdiff = np.abs(np.mean(nom_saln[dup_instrID == grp]) - nom_saln[flagged_index == nan_id])
								#Determine which averages closest to 0
								diff[ii] = np.mean([Tdiff,Sdiff])
							nans_include[no_id[np.argmin(diff)]].append(nan_id)

					#Using the nans_include variable, determine what we're dealing with
					if list(nans_include.keys()) == ['']:
						no_id = np.array([''])

					#Cycle through each of the instruments
					for ids in no_id:
						#Isolate the relevant profiles (plus nans)
						dup_instrID[np.isin(flagged_index,nans_include[ids])] = ids
						fi_id = flagged_index[dup_instrID == ids]
						fs_id = flagged_source[dup_instrID == ids]

						#Determine the number of temperature and salinity measurements in each
						temp = ds_merged.temperature[fi_id].values
						saln = ds_merged.salinity[fi_id].values
						#Determine the number of recordings in each cast
						nom_temp = (~np.isnan(temp)).sum(axis=1)
						nom_saln = (~np.isnan(saln)).sum(axis=1)

						#If climate and NAFC_oceanography are present, take climate
						if fs_id.size == 1:
							continue
						if fs_id.size == 2:
							if np.isin('Climate',fs_id) and np.isin('NAFC-Oceanography',fs_id):
								duplicates_flag[fi_id[fs_id == 'NAFC-Oceanography']] = True
								continue
							#If NAFC_oceanography and CIOOS_NAFC are present, take CIOOS_NAFC
							elif np.isin('CIOOS_NAFC',fs_id) and np.isin('NAFC-Oceanography',fs_id):
								duplicates_flag[fi_id[fs_id == 'NAFC-Oceanography']] = True
								continue
							#If anything and BIO-OMM are present, take anything
							elif np.isin(fs_id,'BIO-OMM').sum() == 1:
								duplicates_flag[fi_id[fs_id == 'BIO-OMM']] = True
								continue
							else:
								None

						#If Climate is present
						if np.isin('Climate',fs_id):

							#If Climate is present once
							if np.isin(fs_id,'Climate').sum() == 1:
								#If the Climate location has at least 90% of recordings, keep
								clim_place1 = {}
								place1 = {}

								#Determine the number of recordings in all spots
								for ii,value in enumerate(fs_id):
									if value != 'Climate':
										place1[ii] = int(nom_temp[ii] + nom_saln[ii])
									else:
										clim_place1[value] = int(nom_temp[ii] + nom_saln[ii])

								#Determine if Climate has enough to be kept (90% of others)
								place1_check = np.zeros(fs_id.size).astype(bool)
								x = 0
								for ii,value in enumerate(fs_id):
									if value != 'Climate':

										#Check number of measurements
										if place1[ii]*0.9 > clim_place1['Climate']:
											place1_check[x] = True
										x += 1

								#If all aren't, keep Climate, discard others
								if not any(place1_check):
									duplicates_flag[fi_id[~(fs_id == 'Climate')]] = True

								#Else take the cast with the highest number
								else:
									duplicates_flag[fi_id[fi_id != fi_id[np.argmax(nom_temp+nom_saln)]]] = True

							#If Climate is the only source present
							elif np.isin(fs_id,'Climate').all():
								duplicates_flag[fi_id[fi_id != fi_id[np.argmax(nom_temp+nom_saln)]]] = True

							#If Climate is present more than twice, but other variables are present
							else:
								#If the Climate location has at least 90% of recordings, keep
								clim_place1 = {}
								place1 = {}

								#Determine the number of recordings in all  spots
								for ii,value in enumerate(fs_id):
									if value != 'Climate':
										place1[ii] = int(nom_temp[ii] + nom_saln[ii])
									else:
										clim_place1[ii] = int(nom_temp[ii] + nom_saln[ii])

								#Determine which group has the highest number
								place1_winner = np.max([place1[ii]*0.9 for ii in place1])
								clim_place1_winner = np.max([clim_place1[ii] for ii in clim_place1])

								#If place1 group is higher
								if place1_winner > clim_place1_winner:
									duplicates_flag[fi_id[fi_id != fi_id[np.argmax(nom_temp+nom_saln)]]] = True

								#If clim_place1 group is higher
								else:
									place1 = np.argmax((nom_temp+nom_saln)[fs_id == 'Climate'])
									#Remove the non-Climate casts
									duplicates_flag[fi_id[np.isin(fi_id,fi_id[fs_id == 'Climate'],invert=True)]] = True
									#Remove the Climate casts that didn't win
									duplicates_flag[fi_id[fs_id == 'Climate'][fi_id[fs_id == 'Climate'] != fi_id[fs_id == 'Climate'][place1]]] = True

						#If Climate is not present
						else:
							#Determine which cast has the highest number of measurements and keep that one
							duplicates_flag[fi_id[fi_id != fi_id[np.argmax(nom_temp+nom_saln)]]] = True

		#Remove the duplicates from the ds
		ds_merged = ds_merged.sel(time=~duplicates_flag)
		ds_merged = ds_merged.sortby('time')

		#Manually set stations of interest
		stns_of_interest = {
		'S27-01':	[47.54667,-52.58667],
		'HL-02':	[44.2670,-63.3170],
		'Prince-5':	[44.9300,-66.8500],
		}
		#Set the size box around which to isolate 
		padding = 0.0125

		#Cycle through each station
		for i in stns_of_interest:

			#Define the box coordinates
			lon_min = stns_of_interest[i][1]-padding 
			lon_max = stns_of_interest[i][1]+padding
			lat_min = stns_of_interest[i][0]-padding
			lat_max = stns_of_interest[i][0]+padding

			#Isolate the latitudes and longitudes of interest
			lat_flag = (ds_merged.latitude.values > lat_min)*(ds_merged.latitude.values < lat_max)
			lon_flag = (ds_merged.longitude.values > lon_min)*(ds_merged.longitude.values < lon_max)

			#Change the comments value for flagged casts
			ds_merged['station_ID'][lat_flag*lon_flag] = i

		#Check to ensure all points are within the desired year
		ds_merged = ds_merged.sel(time = ds_merged['time.year'] == int(year))

		#Re-declare the file attributes
		ds_merged.attrs = {
		'Conventions': 'CF-1.6',
		'title': 'Canadian Atlantic Shelf Temperature-Salinity (CASTS) Dataset',
		'institution': 'Northwest Atlantic Fisheries Centre, Fisheries and Oceans Canada',
		'source': 'https://github.com/OceanAccessLab/CASH',
		'references': 'Cyr et al., 2022',
		'description': 'Temperature and salinity cast records from the CASTS Dataset.',
		'history': 'Created ' + tt.ctime(tt.time())
		}

		#Flatten the arrays to ensure they save properly
		for i in ['trip_ID','source','instrument_ID','instrument_type','file_names','station_ID','sounder_depth']:
			ds_merged[i] = ds_merged[i].astype(str)

		#Save the merged dataset
		ds_merged.to_netcdf(file_output+str(year)+'.nc',\
			encoding={'time':{'units': 'seconds since 1900-01-01 00:00:00', 'calendar': 'gregorian'}})
		ds_merged.close()



'''
TESTING AREA
file_input='/gpfs/fs7/dfo/dpnm/joc000/Data/CASTS/Data_Products/1_combined_raw/'
file_output='/gpfs/fs7/dfo/dpnm/joc000/Data/CASTS/Data_Products/2_restrict_region/'
years=np.arange(1912,2023+1).astype(str)
'''



def area_isolate(
	file_input,
	file_output,
	years,
		):

	#Cycle through each year
	for year in years:

		#Open the netcdf file
		if os.path.exists(file_input+str(year)+'.nc'):
			ds = xr.open_dataset(file_input+str(year)+'.nc')
		else:
			continue

		#Perform the cutoff for the temperature and salinity
		#Pre-define cut-offs for each variable that are non-realistic
		cutoff = {
		'temperature':	[-1.825,35],
		'salinity':		[0,45],
		}

		#Isolate the data
		place1 = ds.where(ds['temperature'] >= cutoff['temperature'][0], drop=False)  
		place1 = place1.where(place1['temperature'] <= cutoff['temperature'][1], drop=False)  
		place2 = ds.where(ds['salinity'] >= cutoff['salinity'][0], drop=False)  
		place2 = place2.where(place2['salinity'] <= cutoff['salinity'][1], drop=False)  
		ds['temperature'] = place1.temperature
		ds['salinity'] = place2.salinity

		#Isolate for the region of interest
		ds = ds.sel(time = (ds.latitude >= 35)*(ds.latitude <= 80))
		ds = ds.sel(time = (ds.longitude >= -100)*(ds.longitude <= -42))
		#Stop if file is empty
		if ds.time.size == 0:
			continue

		#TEMPORARY
		#Remove NAFC-Oceanography casts startiong with 2017_79 during 2017 or 2012_39 for 2012
		if np.isin(year, ['2017','2018']):

			#Write down the file name
			file_names = ds.file_names.values
			place1 = np.zeros(file_names.size).astype(bool)
			for i,value in enumerate(file_names):
				if str(value).startswith('2017_79'):
					place1[i] = True

			#Remove the flagged casts
			ds = ds.sel(time = (~place1))

		#Bin-average the values over specified depths
		bin1 = np.arange(0,1000+1,1)
		bin2 = np.arange(1010,2000+10,10)
		bin3 = np.arange(2100,5000+100,100)
		bins = np.concatenate((bin1,bin2,bin3))

		#Bin-average the depths
		temp_bins = ds.temperature.groupby_bins('level',bins=bins,include_lowest=True).mean()
		saln_bins = ds.salinity.groupby_bins('level',bins=bins,include_lowest=True).mean()
		ds_new = temp_bins.to_dataset()
		ds_new['salinity'] = saln_bins
		for i in [
			'longitude',
			'latitude',
			'trip_ID',
			'source',
			'instrument_ID',
			'instrument_type',
			'file_names',
			'station_ID',
			'sounder_depth']:
			ds_new[i] = ds[i]
		ds_new.attrs = ds.attrs
		ds_new = ds_new.assign_coords({'level_bins':('level_bins',bins[:-1])})
		ds_new = ds_new.rename({'level_bins':'level'})

		#Flatten the arrays to ensure they save properly
		for i in ['trip_ID','source','instrument_ID','instrument_type','file_names','station_ID','sounder_depth']:
			ds_new[i] = ds_new[i].astype(str)

		#Save the isolated netcdf files
		ds_new.to_netcdf(file_output+\
			str(year)+'.nc',
			encoding={'time':{'units': 'seconds since 1900-01-01 00:00:00', 'calendar': 'gregorian'}})
		ds_new.close()


'''
TESTING AREA
file_input='/gpfs/fs7/dfo/dpnm/joc000/Data/CASTS/Data_Products/2_restrict_region/'
file_output='/gpfs/fs7/dfo/dpnm/joc000/Data/CASTS/Data_Products/3_depth_check/'
years=np.arange(1912,2023+1).astype(str)
res=5
bath_path='/gpfs/fs7/dfo/dpnm/joc000/Data/GEBCO_2022/gebco_2022_n80.0_s35.0_w-100.0_e-42.0.nc'
'''


def depth_filter(
	file_input,
	file_output,
	years,
	res=5,
	bath_path='/gpfs/fs7/dfo/dpnm/joc000/Data/GEBCO_2022/gebco_2022_n80.0_s35.0_w-100.0_e-42.0.nc'
	):

	#Import the bathymetry data
	ds_bath = xr.open_dataset(bath_path)

	#Choose a resolution (by skipping) for the bathymetry
	bath_lon,bath_lat = np.meshgrid(ds_bath.lon[::res].values,ds_bath.lat[::res].values)
	bath_elv = ds_bath.elevation[::res,::res].values.astype(float)
	bath_elv[bath_elv >= 0] = np.nan

	#Define the indices for the bathymetry
	indx_lon,indx_lat = np.meshgrid(np.arange(bath_lon.shape[1]),np.arange(bath_lat.shape[0]))

	#Cycle through each year
	for year in years:

		#Open up the relevant data set
		if os.path.exists(file_input+str(year)+'.nc'):
			ds = xr.open_dataset(file_input+year+'.nc')
		else:
			continue

		#Record all the latitudes and longitude
		ds_lon = ds.longitude.values
		ds_lat = ds.latitude.values

		#Determine the maximum depth of the casts
		#Load in the temperature and salinity
		temp = ds.temperature.values
		saln = ds.salinity.values

		#Determine the depth of each cast (deepest measurement)
		depth_max_temp = np.zeros(ds.time.size)
		depth_max_saln = np.zeros(ds.time.size)
		level = ds.level.values

		#Cycle through each cast
		for i in np.arange(ds.time.size):

			#Determine where the last measurement is
			if np.isnan(temp[i]).all() == False:
				depth_max_temp[i] = level[np.where(np.isnan(temp[i,:]) == False)[0][-1]]
			if np.isnan(saln[i]).all() == False:
				depth_max_saln[i] = level[np.where(np.isnan(saln[i,:]) == False)[0][-1]]

		#Determine where the temperature or salinity is the maximum
		depth_max = np.zeros(ds.time.size)
		place1 = depth_max_temp > depth_max_saln
		depth_max[place1] = depth_max_temp[place1]
		depth_max[~place1] = depth_max_saln[~place1]

		#Determine the index in bathymetry for each of the casts
		ds_xcoord = interpolate.griddata(
		np.array((bath_lon.flatten(),bath_lat.flatten())).T,
		indx_lon.flatten().T,
		(ds_lon,ds_lat),
		method='nearest'
		)
		ds_ycoord = interpolate.griddata(
		np.array((bath_lon.flatten(),bath_lat.flatten())).T,
		indx_lat.flatten().T,
		(ds_lon,ds_lat),
		method='nearest'
		)

		#Determine the GEBCO bathymetry at each cast point
		depth_GEBCO = np.abs(bath_elv[ds_ycoord,ds_xcoord])

		#Determine whether casts are exceeding the depth threshold
		#If the cast recorded depth is greater than 200m
		#If the cast depth is greater than 1.5 times the GEBCO depth
		#Then remove cast
		depth_filter = np.zeros(depth_max.size).astype(bool)
		for i in np.arange(depth_filter.size):
			if depth_max[i] > 200:
				if (depth_max[i] - depth_GEBCO[i]) > 300:
					depth_filter[i] = True

		#Save the casts with the depth-filtered casts removed
		ds = ds.sel(time = ~depth_filter)
		if ds.time.size == 0:
			continue

		#Flatten the arrays to ensure they save properly
		for i in ['trip_ID','source','instrument_ID','instrument_type','file_names','station_ID','sounder_depth']:
			ds[i] = ds[i].astype(str)
		ds.to_netcdf(file_output+year+'.nc')
		ds.close()



'''
TESTING AREA
file_input='/gpfs/fs7/dfo/dpnm/joc000/Data/CASTS/Data_Products/3_depth_check/'
file_output='/gpfs/fs7/dfo/dpnm/joc000/Data/CASTS/Data_Products/4_outlier_check/'
climatology_path='/gpfs/fs7/dfo/dpnm/joc000/Data/CASTS/Data_Products/climatology/'
years=np.arange(1912,2023+1).astype(str)
'''



def outlier_check(
	file_input,
	file_output,
	climatology_path,
	years,
	):

	#Cycle through the years
	for year in years[:]:

		#Import the netcdf dataset, one year at a time
		if os.path.exists(file_input+str(year)+'.nc'):
			ds = xr.open_dataset(file_input+str(year)+'.nc')
		else:
			continue

		#Determine the coresponding coordinates for each cast in the climatology
		#Define the lat/lon for the NAFC
		NAFC_lat = ds.latitude.values
		NAFC_lon = ds.longitude.values

		#Import the lat/lon from the climatology
		ds_clim = xr.open_dataset(climatology_path+'bootstrap_climatology_01.nc')
		clim_lon,clim_lat = np.meshgrid(ds_clim.longitude.values, ds_clim.latitude.values)
		indx_lon,indx_lat = np.meshgrid(np.arange(clim_lon.shape[1]),np.arange(clim_lat.shape[0]))

		#Determine the indices for each of the NAFC casts
		NAFC_xcoord = interpolate.griddata(
		np.array((clim_lon.flatten(),clim_lat.flatten())).T,
		indx_lon.flatten().T,
		(NAFC_lon,NAFC_lat),
		method='nearest'
		)
		NAFC_ycoord = interpolate.griddata(
		np.array((clim_lon.flatten(),clim_lat.flatten())).T,
		indx_lat.flatten().T,
		(NAFC_lon,NAFC_lat),
		method='nearest'
		)

		#Check to see that index is correct
		for i,value in enumerate(NAFC_lon):
			if value < clim_lon[NAFC_ycoord[i],NAFC_xcoord[i]]:
				NAFC_xcoord[i] = NAFC_xcoord[i]-1
		for i,value in enumerate(NAFC_lat):
			if value < clim_lat[NAFC_ycoord[i],NAFC_xcoord[i]]:
				NAFC_ycoord[i] = NAFC_ycoord[i]-1

		#Determine which of the cast measurements are outliers
		outlier_mask = {}

		#Run through month by month
		for i in np.unique(ds['time.month'].values):

			#Import the climatology data
			ds_clim = xr.open_dataset(climatology_path+'bootstrap_climatology_'+"%.2d" % (i)+'.nc')

			#Determine a monthly filter for the NAFC data
			month_filt = ds['time.month'].values == i
			outlier_mask[i] = {}

			#Cycle through temperature and salinity
			for ii in ['temperature','salinity']:

				#Create a variable to mask data
				outlier_mask[i][ii] = np.zeros((month_filt.sum(),ds_clim.level.size)).astype(bool)

				#Determine the climatology values at each cast location
				data_mean = ds_clim[ii+'_mean'].values[:,NAFC_ycoord[month_filt],NAFC_xcoord[month_filt]]
				data_ster = ds_clim[ii+'_ster'].values[:,NAFC_ycoord[month_filt],NAFC_xcoord[month_filt]]
				
				#Isolate the data
				data_NAFC = ds[ii][month_filt].values
				
				#Blanket 50 ster window for climatology
				place2 = \
				(data_NAFC > (data_mean.T + (data_ster.T*40)))+\
				(data_NAFC < (data_mean.T - (data_ster.T*40)))

				#Determine if there 25% of measurements are outliers
				nom_NAFC = data_NAFC.shape[1] - np.isnan(data_NAFC).sum(axis=1)
				nom_outliers = place2.sum(axis=1)
				np.seterr(divide='ignore', invalid='ignore')
				place2[nom_outliers/nom_NAFC > 0.25,:] = True
				if int(year) >= 1970:
					#For years after 1970, don't remove individual outliers
					place2[nom_outliers/nom_NAFC <= 0.25,:] = False
				outlier_mask[i][ii][:,:] = place2[:,:]

		#Create one universal filter for temperature and salinity
		outlier_combined = {}
		empties = {}
		max_outliers = {}

		#Cycle through temperature and salinity
		for ii in ['temperature','salinity']:

			#Create the empty mask
			outlier_combined[ii] = np.zeros((month_filt.size,ds_clim.level.size)).astype(bool)

			#Cycle through each month
			for i in np.unique(ds['time.month'].values):

				#Create the month filter, populate the outlier combined
				month_filt = ds['time.month'].values == i
				outlier_combined[ii][month_filt] = outlier_mask[i][ii]

			#Add back to the original combined dataset
			ds[ii][:,:] = ds[ii].values*(~outlier_combined[ii].astype(bool))

			#Determine where all the empties are
			empties[ii] = np.isnan(ds[ii]).sum(axis=1).values
			empties[ii] = empties[ii] == ds.level.size

			#Determine where the flagged cast for removal is 
			max_outliers[ii] = outlier_combined[ii].sum(axis=1) == ds.level.size


		#Determine where there are true empties (missing temperature and salinity)
		empties = empties['temperature']*empties['salinity']

		#Remove the cast with one or both variables exceeding the max number of outliers
		max_outliers = max_outliers['temperature']*max_outliers['salinity']

		#Combine the two
		cast_removal = empties+max_outliers

		#Change the temperature and salinity outliers to nans
		for ii in ['temperature','salinity']:
			place1 = ds[ii].values
			place1[(outlier_combined[ii])] = np.nan
			ds[ii][:,:] = place1

		#Remove the casts 
		ds = ds.sel(time = ~cast_removal)
		if ds.time.size == 0:
			continue

		#Save the combined dataset
		#Flatten the arrays to ensure they save properly
		for i in ['instrument_ID','file_names','sounder_depth','instrument_type','station_ID','trip_ID','source']:
			ds[i] = ds[i].astype(str)

		#Save the isolated netcdf files
		ds.to_netcdf(file_output+str(year)+'.nc',
			encoding={'time':{'units': 'seconds since 1900-01-01 00:00:00', 'calendar': 'gregorian'}})
		ds.close()




'''
TESTING AREA
file_input='/gpfs/fs7/dfo/dpnm/joc000/Data/CASTS/Data_Products/4_outlier_check/'
file_output='/gpfs/fs7/dfo/dpnm/joc000/Data/CASTS/Data_Products/5_final_product/'
station_path='/gpfs/fs7/dfo/dpnm/joc000/Data/CASTS/Data_Input/Other/station_coordinates/'
years=np.arange(1912,2023+1).astype(str)
padding=0.0125
'''



def station_check(
	file_input,
	file_output,
	station_path,
	years,
	padding=0.0125
	):

	#Import the station ID data
	station_ds = pd.read_excel(station_path+'Station_ID.xlsx')

	#Isolate the data of interest
	program = station_ds['Program'].values
	station_ID = station_ds['Station_ID'].values
	station_lat = station_ds['Latitude_decimal'].values
	station_lon = station_ds['Longitude_decimal'].values*-1

	#Cycle through each year
	for year in years:

		#Import the dataset
		if os.path.exists(file_input+str(year)+'.nc'):
			ds = xr.open_dataset(file_input+year+'.nc')
		else:
			continue
		ds = ds.sel(time = (ds.instrument_ID != 'ML'))

		#Import the latitude, longitude, station ID
		ds_lat = ds.latitude.values 
		ds_lon = ds.longitude.values
		ds_ID = ds.station_ID.values
		ds_instrument_ID = ds.instrument_ID.values

		#Define the new station_ID data
		new_ID = np.full(ds_ID.size,'',dtype='<U32')

		#Cycle through each station ID from the xlsx file
		for i in np.arange(station_ID.size):

			#Isolate which casts fall into them 
			place1 = (ds_lon >= (station_lon[i] - padding))*(ds_lon <= (station_lon[i] + padding))
			place2 = (ds_lat >= (station_lat[i] - padding))*(ds_lat <= (station_lat[i] + padding))

			#Write the new station ID
			new_ID[place1*place2] = program[i]+'_'+station_ID[i]

		#Record the new data
		ds['station_ID_manual'] = (('time'), new_ID)

		#Organize instrument ID into categories
		instrument_call = {}
		instrument_call['CT'] = ['19899','19P-7','19P45','19P53','GL','OLABS','STD12','NBmun','0CTD','0CTD1','1CTD','1CTD1','CD','CT','CTD','CU','FAPCD','FAPCU','GLCTD','GLctd','MEDCD','S0000','S0002','S0036','S0038','S0078','S0097','S0103','S0118','S0129','S0212','S0277','S0281','S0344','S0394','S0454','S0466','S0580','S0582','S0583','S0674','S0688','S0689','S0690','S0691','S0692','S0758','S0803','S0845','S0846','S0847','S0910','S0911','S0912','S0935','S1003','S1019','S1020','S1021','S1098','S1101','S1145','S1146','S1221','S1237','S1238','S1247','S1257','S1308','S1309','S1310','S1311','S1312','S1313','S1314','S1315','S1316','S1317','S1318','S1319','S1460','S2245','S2246','S2247','S2248','S4016','S4017','S4018','S4019','S4020','S4578','S4579','S4580','S4581','S4582','S4776','S4777','S5117','S7849','S7853','S7854','S7855','S7992','S8013','S8014','S9999','SB', 'SB002', 'SB036', 'SB038','SB12','SB279','SB280','SB394','SB450','SB453','SB454','SB455','SB580','SB581','SB583','SB688','SB689','SB690','SB691','SB692','SB73','SB74','SB81','SB82','SB822','SB83','SB845','SB880','SB890','SB900','SB910','SB920','SB935','SBMUN']
		instrument_call['BO'] = ['BO','FAPBO','FAPbo','FRCBO']
		instrument_call['BF'] = ['BF','FAPBF']
		instrument_call['TO'] = ['TO','FAPTO']
		instrument_call['DT'] = ['DT','BT']
		instrument_call['MB'] = ['MB']
		instrument_call['XB'] = ['XB','XBT05','XBT06','XBT07','XBT10']
		instrument_call['PF'] = ['PF','B3']
		instrument_call[''] = ['nan','21573.0']

		#Define the new instrument_ID data
		new_instrument_ID = np.full(ds_ID.size,'',dtype='<U32')

		#Organize the measurements into categories
		for i,value in enumerate(ds_instrument_ID):
			for ii in instrument_call:
				if np.isin(value, instrument_call[ii]):
					new_instrument_ID[i] = ii
					break

		#Record the new data
		ds['instrument_ID_manual'] = (('time'), new_instrument_ID)

		#LAST MINUTE CHANGE - THIS SHOULD BE MOVED TO FIRST STEP
		#Rename IEO 
		source = ds.source.values
		source[source == 'IEO-Spain'] = 'EU_NAFO'
		source[source == 'BIO-OMO'] = 'BIO-OMM'
		ds['source'][:] = source

		#Edit the attributes to be constant for each year
		#Time
		ds['time'].attrs['standard_name'] = 'time'
		ds['time'].attrs['variable_name'] = 'time'
		ds['time'].encoding['units'] = 'seconds since 1900-01-01 00:00:00'
		#ds['time'].attrs['calendar'] = 'gregorian'

		#Level
		ds['level'].attrs['standard_name'] = 'sea_water_pressure'
		ds['level'].attrs['maximum_level'] = '5000'
		ds['level'].attrs['units'] = 'dbar'
		ds['level'].attrs['variable_name'] = 'level'
		ds['level'].attrs['example_boundary'] = 'level 0: (0, 1], level 4900: (4900,5000]'

		#Temperature
		ds['temperature'].attrs['standard_name'] = 'sea_water_potential_temperature'
		ds['temperature'].attrs['variable_name'] = 'temperature'
		ds['temperature'].attrs['units'] = 'degree_C'
		ds['temperature'].attrs['valid_min'] = '-1.825'
		ds['temperature'].attrs['valid_max'] = '35'
		ds['temperature'].attrs['missing_value'] = 'nan'

		#Salinity
		ds['salinity'].attrs['standard_name'] = 'sea_water_practical_salinity'
		ds['salinity'].attrs['variable_name'] = 'salinity'
		ds['salinity'].attrs['units'] = 'psu'
		ds['salinity'].attrs['valid_min'] = '0'
		ds['salinity'].attrs['valid_max'] = '45'
		ds['salinity'].attrs['missing_value'] = 'nan'

		#Longitude
		ds['longitude'].attrs['standard_name'] = 'longitude'
		ds['longitude'].attrs['variable_name'] = 'longitude'
		ds['longitude'].attrs['units'] = 'degree_east'

		#Latitude
		ds['latitude'].attrs['standard_name'] = 'latitude'
		ds['latitude'].attrs['variable_name'] = 'latitude'
		ds['latitude'].attrs['units'] = 'degree_north'

		#Trip_ID
		ds['trip_ID'].attrs['standard_name'] = 'trip_identification'
		ds['trip_ID'].attrs['variable_name'] = 'trip_ID'

		#Source
		ds['source'].attrs['standard_name'] = 'cast_source'
		ds['source'].attrs['variable_name'] = 'source'

		#Instrument_ID
		ds['instrument_ID'].attrs['standard_name'] = 'instrument_identification'
		ds['instrument_ID'].attrs['variable_name'] = 'instrument_ID'

		#Instrument_type
		ds['instrument_type'].attrs['standard_name'] = 'instrument_type'
		ds['instrument_type'].attrs['variable_name'] = 'instrument_type'

		#File_names
		ds['file_names'].attrs['standard_name'] = 'file_names'
		ds['file_names'].attrs['variable_name'] = 'file_names'

		#Station_ID
		ds['station_ID'].attrs['standard_name'] = 'station_identification'
		ds['station_ID'].attrs['variable_name'] = 'station_ID'

		#Sounder_depth
		ds['sounder_depth'].attrs['standard_name'] = 'sounder_depth'
		ds['sounder_depth'].attrs['variable_name'] = 'sounder_depth'
		ds['sounder_depth'].attrs['units'] = 'dbar'

		#Station_ID_manual
		ds['station_ID_manual'].attrs['standard_name'] = 'station_identification_manual'
		ds['station_ID_manual'].attrs['variable_name'] = 'station_ID_manual'

		#Instrument_ID_manual
		ds['instrument_ID_manual'].attrs['standard_name'] = 'instrument_identification_manual'
		ds['instrument_ID_manual'].attrs['variable_name'] = 'instrument_ID_manual'

		#Global Attributes
		ds.attrs['title'] = 'Canadian Atlantic Shelf Temperature-Salinity (CASTS) Data Product, 2023, V2.0'
		ds.attrs['institution'] = 'Northwest Atlantic Fisheries Centre, Fisheries and Oceans Canada (DFO)'
		ds.attrs['source'] = 'https://doi.org/10.20383/102.0739'
		ds.attrs['project_url'] = 'https://github.com/OceanAccessLab/CASTS'
		ds.attrs['references'] = 'Currently NA'
		ds.attrs['description'] = 'Temperature and salinity cast records from the CASTS Dataset'
		ds.attrs['comment'] = 'This work is licensed under the Creative Commons Attribution 4.0 International License. To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/.'
		ds.attrs['Conventions'] = 'CF-1.6'
		ds.attrs['doi'] = '10.20383/102.0739'
		ds.attrs['creator_names'] = 'Jonathan Coyne, Frederic Cyr'
		ds.attrs['creator_emails'] = 'jonathan.coyne@dfo-mpo.gc.ca, frederic.cyr@dfo-mpo.gc.ca'
		ds.attrs['geospatial_lon_max'] = '-42degE'
		ds.attrs['geospatial_lon_min'] = '-100degE'
		ds.attrs['geospatial_lat_min'] = '35degN'
		ds.attrs['geospatial_lat_max'] = '80degN'
		ds.attrs['file_year'] = year
		ds.attrs['keywords'] = 'Ocean Science, Ocean Temperature, Ocean Salinity, Historical Data'

		#Flatten the arrays to ensure they save properly
		for i in ['trip_ID','source','instrument_ID','instrument_ID_manual','instrument_type','file_names','station_ID','station_ID_manual','sounder_depth']:
			ds[i] = ds[i].astype(str)

		#Save the isolated netcdf files
		ds.to_netcdf(file_output+str(year)+'.nc',
			encoding={'time':{'units': 'seconds since 1900-01-01 00:00:00', 'calendar': 'gregorian'}})
		ds.close()
		exp = 'ncks -O --mk_rec_dmn time '+file_output+year+'.nc '+file_output+year+'.nc'
		os.system(exp)
		#Compress all files for memory purposes
		exp = 'nccopy -d 9 '+file_output+year+'.nc '+file_output+'compressed/'+year+'.nc'
		os.system(exp)








