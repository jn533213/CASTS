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
sys.path.append('/gpfs/fs7/dfo/dpnm/joc000/Data/CASTS/Data_Input/Other/Charlies_Scripts')
import cnv_tk
import matplotlib
matplotlib.interactive(True)
import matplotlib.pyplot as plt



'''

The purpose of this script is to streamline the CASTS update process.
All files and sources should be able to be produced from one central script.

This package will include all NAFC-Oceanography QA/QC functions.

'''

'''
TEST INFO
input_folder = '/gpfs/fs7/dfo/dpnm/joc000/Data/CASTS/Data_Input/NAFC_Oceanography/pfiles/2011/'
output_folder = '/gpfs/fs7/dfo/dpnm/joc000/Data/CASTS/Data_Input/NAFC_Oceanography/1_files_raw/2011.nc'
year = '2011'
max_depth = 5000
'''


def convertpfile(
	input_folder,
	output_folder,
	year,
	max_depth=5000,
	):
	'''
	input_folder: path location to yearly .pfiles
	output_folder: path location including name for netcdf file
	year: year that data is from (str)
	max_depth (default 5000): maximum depth bin averaging to in meters
	'''

	#Determine the files present
	year_files = glob.glob(input_folder+'/*.p'+year)

	#Create empty output variables for each potential file, 1D
	latitude = np.array([])
	longitude = np.array([])
	time_1D = np.array([])
	trip_id = np.array([])
	comment = np.array([])
	instrument_id = np.array([])
	instrument_type = np.array([])
	sounder = np.array([])
	file_name = np.array([])

	#Create empty output variables for 2D variables
	temperature = []
	salinity = []
	density = []

	#Cycle through each of the files
	for file in year_files:

		#Import the file
		my_file = open(file,encoding='utf8',errors='ignore') #Default is utf-8 encoding
		my_data = my_file.read().split('\n')
		my_file.close()

		#Isolate and record all header info
		str_start = [0,10,13,19,24,29,40,46,51,57,62,64]
		str_end = [8,12,18,23,29,40,46,51,57,61,63,78]
		header_info = []
		for i,value in enumerate(str_start):
			try: 
				header_info.append(my_data[1][value:str_end[i]].strip())
			except Exception:
				header_info.append('nan')

		#Replace 24:00 with 00:00 for time
		if header_info[6] == '24:00':
			header_info[6] = '00:00'

		#This section will cycle through to determine problem meta-data
		#Check to see that time is properly formatted
		try:
			pd.Timestamp(header_info[5]+' '+header_info[6])
		except Exception:
			print(file.split('/')[-1]+': problem time encountered, skipped.')
			continue
		#Check to see if data is GTS origin or BA origin, skip if True
		if header_info[8].upper()[-2:] == 'TE' or header_info[8].upper()[-2:] == 'BA':
			continue
		#Check to see that header info is correct
		if my_data[0] != 'NAFC_Y2K_HEADER':
			print(file.split('/')[-1]+': problem header encountered, skipped.')
			continue
		#Check to see that lat and lon aren't 0
		if header_info[1] == '0' and header_info[3] == '0':
			print(file.split('/')[-1]+': problem coordinates encountered, skipped.')
			continue
		#Check to see that lat isn't greater than 90 or less than -90
		if float(header_info[1]) > 90 or float(header_info[1]) < -90:
			print(file.split('/')[-1]+': problem latitude encountered, skipped.')
			continue
		#Check to see that lon isn't greater than 180 or less than -180
		if float(header_info[3]) > 180 or float(header_info[3]) < -180:
			print(file.split('/')[-1]+': problem longitude encountered, skipped.')
			continue

		#Record the meta-data
		latitude = np.append(latitude, float(header_info[1]) + float(header_info[2])/60.0)
		longitude = np.append(longitude,  np.sign(float(header_info[3]))*(np.abs(float(header_info[3])) + float(header_info[4])/60.0))
		time_1D = np.append(time_1D, pd.Timestamp(header_info[5]+' '+header_info[6]))
		trip_id = np.append(trip_id, header_info[0])
		sounder = np.append(sounder, int(header_info[7]))
		instrument_id = np.append(instrument_id, header_info[8])
		instrument_type = np.append(instrument_type, header_info[10])
		comment = np.append(comment, header_info[11])
		file_name = np.append(file_name, file.split('/')[-1])

		#Isolate the 2D data and determine column names
		for i,value in enumerate(my_data):
			if value == '-- DATA --':
				columns = my_data[i-1]
				my_data_2D = my_data[i+1:-1]
				my_data_2D = list(filter(None, my_data_2D))
				break

		#Check to see that data has at least one measurement
		if np.size(my_data_2D) == 0 or my_data_2D[0][0] == 'Empty' or my_data_2D[0] == 'Empty DataFrame':
			print(file.split('/')[-1]+': no data available, filled empty.')
			temperature.append(np.full(max_depth,np.nan))
			salinity.append(np.full(max_depth,np.nan))
			density.append(np.full(max_depth,np.nan))
			continue

		#Format into dictionary of arrays
		columns = re.split(r'\s{2,}', columns)
		columns = [i for i in columns if i not in ('')]
		columns = [i.strip() for i in columns]
		for i,value in enumerate(my_data_2D):
			place1 = re.split(r'\s{1,}', value)
			place1 = [ii for ii in place1 if ii not in ('')]
			my_data_2D[i] = place1
		'''
		if np.shape(my_data_2D)[0] == 1:
			print(file.split('/')[-1]+': no data available, filled empty.')
			temperature.append(np.full(max_depth,np.nan))
			salinity.append(np.full(max_depth,np.nan))
			density.append(np.full(max_depth,np.nan))
			continue
		'''
		#my_data_2D = np.array(my_data_2D).T.astype(float)
		try:
			my_data_2D = np.array(my_data_2D).T.astype(float)
		except Exception:
			print(file.split('/')[-1]+': problem data encountered, filled empty.')
			temperature.append(np.full(max_depth,np.nan))
			salinity.append(np.full(max_depth,np.nan))
			density.append(np.full(max_depth,np.nan))
			continue


		#Check to see that data and column number matches
		if np.size(columns) != my_data_2D.shape[0]:
			columns = columns[:my_data_2D.shape[0]]

		#Determine the pressure column place
		pres_place = np.where(np.isin(columns, ['pres','depth','pres-db','prdM','prSM']))[0]

		#If no pressure available, skip (will be later removed as empty)
		if pres_place.size == 0:
			print(file.split('/')[-1]+': no pressure available, filled empty.')
			temperature.append(np.full(max_depth,np.nan))
			salinity.append(np.full(max_depth,np.nan))
			density.append(np.full(max_depth,np.nan))
			continue

		#Format into dictionary
		my_data_2D_dict = {}
		for i,value in enumerate(columns):
			#Bin-average the 2D results
			my_data_2D_dict[value.strip()] = stats.binned_statistic(
				my_data_2D[pres_place[0]].flatten(),
				my_data_2D[i].flatten(),
				'mean',
				bins=np.arange(max_depth+1)
				).statistic

		#Record each of the variables
		temp_names = ['temp','temp90-C','tv290C','t090C']
		if np.isin(columns, temp_names).sum() == 1:
			var_name = np.array(columns)[np.isin(columns, temp_names)]
			temperature.append(my_data_2D_dict[var_name[0]])
		else:
			temperature.append(np.full(max_depth,np.nan))
		saln_names = ['sal','saln-PSU','sal00']
		if np.isin(columns, saln_names).sum() == 1:
			var_name = np.array(columns)[np.isin(columns, saln_names)]
			salinity.append(my_data_2D_dict[var_name[0]])
		else:
			salinity.append(np.full(max_depth,np.nan))
		dens_names = ['sigt','sigma-t','sigma-t00']
		if np.isin(columns, dens_names).sum() == 1:
			var_name = np.array(columns)[np.isin(columns, dens_names)]
			density.append(my_data_2D_dict[var_name[0]])
		else:
			density.append(np.full(max_depth,np.nan))

	#Bring the 2D variables together
	temperature = np.array(temperature)
	salinity = np.array(salinity)
	density = np.array(density)

	#Save the new variables as a netcdf array
	nc_out = nc.Dataset(output_folder,'w')

	#File information
	nc_out.Conventions = 'CF-1.6' #Ask Fred what this means
	nc_out.title = 'NAFC Yearly Netcdf File' #Temporary title for the .nc file
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
	trip_ids = nc_out.createVariable('trip_ID',str,('time'),zlib=True)
	comments = nc_out.createVariable('comments',str,('time'),zlib=True)
	instrument_types = nc_out.createVariable('instrument_type',str,('time'),zlib=True)
	instrument_ids = nc_out.createVariable('instrument_ID',str,('time'),zlib=True)
	file_names = nc_out.createVariable('file_names',str,('time'),zlib=True)

	#Create 2D variables
	temp = nc_out.createVariable('temperature', np.float32, ('time', 'level'), zlib=True, fill_value=-9999)
	saln = nc_out.createVariable('salinity', np.float32, ('time', 'level'), zlib=True, fill_value=-9999)
	dens = nc_out.createVariable('density', np.float32, ('time', 'level'), zlib=True, fill_value=-9999)

	#Variable Attributes
	latitudes.units = 'degree_north'
	longitudes.units = 'degree_east'
	times.units = 'seconds since 1900-01-01 00:00:00'
	times.calendar = 'gregorian'
	levels.units = 'dbar'
	levels.standard_name = "pressure"
	temp.units = 'Celsius'
	temp.long_name = "Water Temperature" # (may be use to label plots)
	temp.standard_name = "sea_water_temperature"
	saln.long_name = "Practical Salinity"
	saln.standard_name = "sea_water_salinity"
	saln.units = "1"
	dens.units = 'Kg m-3'
	dens.long_name = 'Sigma-t'
	dens.standard_name = 'sigma_t'

	#Create a time filter
	filt = np.argsort(time_1D)

	#Fill in the 1D variables 
	latitudes[:] = latitude[filt]
	longitudes[:] = longitude[filt]
	trip_ids[:] = trip_id.astype(str)[filt]
	comments[:] = comment.astype(str)[filt]
	instrument_types[:] = instrument_type.astype(str)[filt]
	instrument_ids[:] = instrument_id.astype(str)[filt]
	file_names[:] = file_name.astype(str)[filt]

	#Fill 2D structure
	temp[:,:] = temperature[filt]
	saln[:,:] = salinity[filt]
	dens[:,:] = density[filt]

	#Convert to time stamps
	time_stamps = np.array([pd.Timestamp(i).to_pydatetime() for i in time_1D[filt]])
	times[:] = nc.date2num(time_stamps, units=times.units, calendar=times.calendar)
	levels[:] = np.arange(max_depth)

	#Save and close the .nc file
	nc_out.close()


'''
TEST INFO
input_folder = '/gpfs/fs7/dfo/dpnm/joc000/Data/CASTS/Data_Input/NAFC_Oceanography/pfiles/2024/'
output_folder = '/gpfs/fs7/dfo/dpnm/joc000/Data/CASTS/Data_Input/NAFC_Oceanography/1_files_raw/2024.nc'
year = '2024'
max_depth = 5000
'''

def convertcnv(
	input_folder,
	output_folder,
	year,
	max_depth=5000,
	):
	'''
	input_folder: path location to yearly .cnv
	output_folder: path location including name for netcdf file
	year: year that data is from (str)
	max_depth (default 5000): maximum depth bin averaging to in meters
	'''

	#Determine the files present
	year_files = glob.glob(input_folder+'/*.pcnv')

	#Create empty output variables for each potential file, 1D
	latitude = np.array([])
	longitude = np.array([])
	time_1D = np.array([])
	trip_id = np.array([])
	comment = np.array([])
	instrument_id = np.array([])
	instrument_type = np.array([])
	sounder = np.array([])
	file_name = np.array([])

	#Create empty output variables for 2D variables
	temperature = []
	salinity = []
	density = []

	#Cycle through each of the files, save as netcdf
	for file in year_files:

		#Create the cast object
		cast = cnv_tk.Cast(file)
		cnv_tk.cnv_meta(cast,file)
		#Create a dataframe
		df = cnv_tk.cnv_to_dataframe(cast)
		df = cnv_tk.df_press_depth(cast)
		df = cnv_tk.StandardizedDF(cast, df)

		#Bin average the variables in 1dbar bins
		max_depth = 5000
		if year == '2024':
			variables = ['Temperature','Secondary_Salinity']
			variable_titles = {'Temperature': 'Temperature','Secondary_Salinity': 'Salinity'}
		else:
			variables = ['Temperature','Salinity']
			variable_titles = {'Temperature': 'Temperature','Salinity': 'Salinity'}
		variables_binned = {}
		for variable in variables:
			if variable in df.columns:
				#Remove this special case once we have an updated file
				if file != '/gpfs/fs7/dfo/dpnm/joc000/Data/CASTS/Data_Input/NAFC_Oceanography/pfiles/2024/cr24st133.pcnv':
					variables_binned[variable_titles[variable]] = stats.binned_statistic(
						df['Pressure'].values.astype(float),
						df[variable].values.astype(float),
						'mean',
						bins=np.arange(max_depth+1)
						).statistic
				else:
					if variable == 'Temperature':
						variables_binned[variable_titles[variable]] = stats.binned_statistic(
							df['Temperature'].values.astype(float),
							df['Secondary_Temperature'].values.astype(float),
							'mean',
							bins=np.arange(max_depth+1)
							).statistic
					if variable == 'Secondary_Salinity':
						variables_binned[variable_titles[variable]] = stats.binned_statistic(
							df['Temperature'].values.astype(float),
							df['Secondary_Salinity'].values.astype(float),
							'mean',
							bins=np.arange(max_depth+1)
							).statistic
			else:
				variables_binned[variable_titles[variable]] = np.full(max_depth,np.nan)
				if variable == 'Secondary_Salinity':
					if 'Salinity' in df.columns:
						variables_binned[variable_titles[variable]] = stats.binned_statistic(
							df['Pressure'].values.astype(float),
							df['Salinity'].values.astype(float),
							'mean',
							bins=np.arange(max_depth+1)
							).statistic

		#Record the meta-data
		latitude = np.append(latitude, cast.Latitude)
		longitude = np.append(longitude, cast.Longitude)
		time_1D = np.append(time_1D, pd.Timestamp(cast.CastDatetime))
		trip_id = np.append(trip_id, cast.triptag.lstrip().rstrip()+'_'+cast.trip.lstrip().rstrip())
		sounder = np.append(sounder, cast.SounderDepth.lstrip().rstrip())
		instrument_id = np.append(instrument_id, cast.InstrumentName.lstrip().rstrip())
		instrument_type = np.append(instrument_type, cast.castType.lstrip().rstrip())
		comment = np.append(comment, cast.comment.lstrip().rstrip())
		file_name = np.append(file_name, file.split('/')[-1])

		#Record each of the variables
		temperature.append(variables_binned['Temperature'])
		salinity.append(variables_binned['Salinity'])

	#Bring the 2D variables together
	temperature = np.array(temperature)
	salinity = np.array(salinity)

	#Save the new variables as a netcdf array
	nc_out = nc.Dataset(output_folder,'w')

	#File information
	nc_out.Conventions = 'CF-1.6' #Ask Fred what this means
	nc_out.title = 'NAFC Yearly Netcdf File' #Temporary title for the .nc file
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
	trip_ids = nc_out.createVariable('trip_ID',str,('time'),zlib=True)
	comments = nc_out.createVariable('comments',str,('time'),zlib=True)
	instrument_types = nc_out.createVariable('instrument_type',str,('time'),zlib=True)
	instrument_ids = nc_out.createVariable('instrument_ID',str,('time'),zlib=True)
	file_names = nc_out.createVariable('file_names',str,('time'),zlib=True)    

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
	#levels.valid_min = 0
	temp.units = 'Celsius'
	temp.long_name = "Water Temperature" # (may be use to label plots)
	temp.standard_name = "sea_water_temperature"
	saln.long_name = "Practical Salinity"
	saln.standard_name = "sea_water_salinity"
	saln.units = "1"
	#saln.valid_min = 0

	#Create a time filter
	filt = np.argsort(time_1D)

	#Fill in the 1D variables 
	latitudes[:] = latitude[filt]
	longitudes[:] = longitude[filt]
	trip_ids[:] = trip_id.astype(str)[filt]
	comments[:] = comment.astype(str)[filt]
	instrument_types[:] = instrument_type.astype(str)[filt]
	instrument_ids[:] = instrument_id.astype(str)[filt]
	file_names[:] = file_name.astype(str)[filt]

	#Fill 2D structure
	temp[:,:] = temperature[filt]
	saln[:,:] = salinity[filt]

	#Convert to time stamps
	time_stamps = np.array([pd.Timestamp(i).to_pydatetime() for i in time_1D[filt]])
	times[:] = nc.date2num(time_stamps, units=times.units, calendar=times.calendar)
	levels[:] = np.arange(max_depth)

	#Save and close the .nc file
	nc_out.close()



def remove_empties(
	input_file,
	output_file,
	flagged_files,
	year,
	):
	'''
	input_file: path location to yearly netcdf file
	output_file: path location including name for netcdf file output
	flagged_files: path location to flagged files (MBT BO pre 1980)
	year (str): year of interest
	'''

	#Import the file of interest
	ds = xr.open_dataset(input_file)

	#Run through each variable and determine the number of nans (missing) in each
	#Pre-define cut-offs for each variable that are non-realistic
	cutoff = {
	'temperature':	[-2,35],
	'salinity':		[0,45],
	}

	#Load in each of the variables one at a time
	missing = {}
	for ii in cutoff.keys():
		place1 = ds[ii].values

		#Perform the cutoffs listed above
		place1[place1 <= cutoff[ii][0]] = np.nan
		place1[place1 >= cutoff[ii][1]] = np.nan

		#Calculate the sum of missing values in each cast
		missing[ii] = np.isnan(place1).sum(axis=1)

	#Determine if a cast is completely empty
	missing = np.stack([missing[ii] for ii in cutoff.keys()])
	missing = missing.sum(axis=0) == missing.shape[0]*ds.level.size

	#Remove the empty casts
	ds = ds.sel(time=~missing)

	#Remove BA casts
	ds = ds.sel(time=~(ds.instrument_ID == 'FAPBA'))
	ds = ds.sel(time=~(ds.instrument_ID == 'MEDBA'))

	#Only run bad MB BO test pre 1980
	if int(year) < 1980:

		#Determine if a problem file is already present
		path_prob = flagged_files
		try:
			flagged_files = np.load(path_prob+year+'_flagged_files.npy')
			#print(year+' Flagged File exists - removing pre-flagged files')

			#Create a mask
			flagged_files = np.array([ii.split('.')[0][5:] for ii in flagged_files])
			file_names = ds.file_names.values.astype(str)
			file_names = np.array([ii.split('.')[0] for ii in file_names])
			prob_flag = np.isin(file_names,flagged_files)

		except Exception:

			#Remove the problem MBT casts
			for ii in np.arange(prob_flag.size):

				#Check to see if cast is MB or BO
				if instrument_ID[ii].endswith('MB') or instrument_ID[ii].endswith('BO'):

					#Define the instrument
					inst = instrument_ID[ii][-2:]

					#Determine the temperature
					temp_slice = ds.temperature[ii].values
					levels = ds.level.values

					#Filter out the nans
					levels = levels[~np.isnan(temp_slice)]
					temp_slice = temp_slice[~np.isnan(temp_slice)]

					#Check to see that there are at least 5 measurements present
					if temp_slice.size >= 5:

						#Determine if any measurement is between 0 and 2
						if np.min(temp_slice) > 0 and np.min(temp_slice) <= 1:

							#Do a manual check to see if it's a problem cast
							plt.figure()
							plt.plot(temp_slice,levels,marker='.')
							plt.vlines(0,ymin=levels[0],ymax=levels[-1],colors='k')
							plt.gca().invert_yaxis()
							plt.title(inst+', '+str(ii)+' out of '+str(prob_flag.size))
							plt.show()
							plt.pause(0.5)
							while True:
								place1 = input('Keep cast? ["" for y,n]: ')
								if place1 not in ('','n'):
									print('Not a valid answer. Please type "" for y or "n".')
								else:
									break
							if place1 == '':
								plt.close()
							elif place1 == 'n':
								prob_flag[ii] = True
								plt.close()

			#Save the flagged files
			prob_flag_files = ds.file_names[prob_flag].values.astype(str)
			np.save(flagged_files+year+'_flagged_files',prob_flag_files)

		#Remove the flagged MBs
		ds = ds.sel(time=~prob_flag)

	#Flatten the values inside the string variables
	for ii in ['trip_ID','comments','instrument_ID','instrument_type','file_names']:
		ds[ii] = ds[ii].astype(str)

	#Save the data now with the empty casts removed
	ds.to_netcdf(output_file)
	ds.close()


def remove_duplicates(
	input_file,
	output_file,
	year,
	):
	'''
	input_file: path location to yearly netcdf file
	output_file: path location including name for netcdf file output
	year (str): year of interest
	'''

	#Load the .nc file, now with no empties
	ds = xr.open_dataset(input_file)

	spatial_rounder = 2 #number of decimal points
	temporal_rounder = 'm' #minutes

	#Determine the times where there are time duplicates
	time = ds.time.values.astype('datetime64['+temporal_rounder+']')
	longitude = ds.longitude.values.round(spatial_rounder)
	latitude = ds.latitude.values.round(spatial_rounder)

	#Round the time to the nearest 10th minute
	time_10min = np.zeros(time.size).astype(str)
	for ii,value in enumerate(time):
		minutes = int(str(value)[-1:])
		if minutes > 5:
			add = np.timedelta64(10-minutes,'m')
			time_10min[ii] = value+add
		else:
			subtract = np.timedelta64(minutes,'m')
			time_10min[ii] = value-subtract
	time_10min = np.array([np.datetime64(ii) for ii in time_10min])


	#Start to bring all the variables together
	spatiotemporal_comp = np.array([
		time_10min.astype(str),
		latitude.astype(str),
		longitude.astype(str)
		])

	a = spatiotemporal_comp.T

	#Merge the time,lat,lon as one string, consider as a whole
	a = [row[0]+','+row[1]+','+row[2] for row in spatiotemporal_comp.T]
	unique_dates,unique_index,unique_count = np.unique(a,return_counts=True,return_index=True)

	#Mark the times as duplicates
	duplicates_flag = np.full(time.shape,False) 

	for ii in np.arange(unique_dates.size):

		#Determine if there is a count higher than 1 (duplicate)
		if unique_count[ii] > 1:

			#Isolate for the casts
			spots = np.where(np.array(a) == unique_dates[ii])[0]

			#SPECIAL CASES
			#Determine if one of the casts is XBT and one is BO, if so skip
			if unique_count[ii] == 2:
				place1 = ds.instrument_ID[spots].values
				if 'bo' in str(place1[0]).lower() or 'bo' in str(place1[1]).lower():
					if 'xb' in str(place1[0]).lower() or 'xb' in str(place1[1]).lower():
						continue

			#Determine if one of the casts is CTD and one is BO, if so skip
			if unique_count[ii] == 2:
				place1 = ds.instrument_ID[spots].values
				if 'bo' in str(place1[0]).lower() or 'bo' in str(place1[1]).lower():
					if 'cd' in str(place1[0]).lower() or 'cd' in str(place1[1]).lower():
						continue

			#Determine if one of the casts is CTD and one is XBT, if so skip
			if unique_count[ii] == 2:
				place1 = ds.instrument_ID[spots].values
				if 'xb' in str(place1[0]).lower() or 'xb' in str(place1[1]).lower():
					if 'cd' in str(place1[0]).lower() or 'cd' in str(place1[1]).lower():
						continue

			#Determine if one of the casts is S1 and the other S0, and they're both fishsets
			if unique_count[ii] == 2:
				place1 = ds.instrument_ID[spots].values
				place2 = ds.comments[spots].values
				if 'fishset' in str(place2[0]).lower() and 'fishset' in str(place2[1]).lower():
					if 's0' in str(place1[0]).lower() or 's0' in str(place1[1]).lower():
						if 's1' in str(place1[0]).lower() or 's1' in str(place1[1]).lower():
	
							#Import the temperature
							temp = ds.temperature[spots].values
							saln = ds.salinity[spots].values

							#Determine the number of recordings in each cast
							nom_temp = (~np.isnan(temp)).sum(axis=1)
							nom_saln = (~np.isnan(saln)).sum(axis=1)							
							
							#Keep the cast with the highest number
							duplicates_flag[spots[np.argmin(nom_temp+nom_saln)]] = True
							continue

			#If repeat is the comment of each cast for two casts
			if unique_count[ii] == 2:
				place2 = ds.comments[spots].values
				if 'repeat' in str(place2[0]).lower() and 'repeat' in str(place2[1]).lower():

					#Import the temperature
					temp = ds.temperature[spots].values
					saln = ds.salinity[spots].values

					#Determine the number of recordings in each cast
					nom_temp = (~np.isnan(temp)).sum(axis=1)
					nom_saln = (~np.isnan(saln)).sum(axis=1)							
					
					#Keep the cast with the highest number
					duplicates_flag[spots[np.argmin(nom_temp+nom_saln)]] = True
					continue

			#If repetitive is the comment of each cast for two casts
			if unique_count[ii] == 2:
				place2 = ds.comments[spots].values
				if 'repetitive' in str(place2[0]).lower() and 'repetitive' in str(place2[1]).lower():

					#Import the temperature
					temp = ds.temperature[spots].values
					saln = ds.salinity[spots].values

					#Determine the number of recordings in each cast
					nom_temp = (~np.isnan(temp)).sum(axis=1)
					nom_saln = (~np.isnan(saln)).sum(axis=1)							
					
					#Keep the cast with the highest number
					duplicates_flag[spots[np.argmin(nom_temp+nom_saln)]] = True
					continue

			#Import the temperature
			temp = ds.temperature[spots].values
			saln = ds.salinity[spots].values

			#Determine the number of recordings in each cast
			nom_temp = (~np.isnan(temp)).sum(axis=1)
			nom_saln = (~np.isnan(saln)).sum(axis=1)

			#Plot the temperature, salinity
			plt.subplot(1,2,1)
			for iii in np.arange(temp.shape[0]):
				plt.scatter(temp[iii],np.arange(ds.level.size),
					label='Temp'+str(iii)+': '+str(nom_temp[iii]))
			plt.legend()
			plt.gca().invert_yaxis()
			plt.xlim(xmin=-2,xmax=35)
			plt.grid()

			plt.subplot(1,2,2)
			for iii in np.arange(temp.shape[0]):
				plt.scatter(saln[iii],np.arange(ds.level.size),
					label='Saln'+str(iii)+': '+str(nom_saln[iii]))
			plt.legend()
			plt.gca().invert_yaxis()
			plt.xlim(xmin=30,xmax=37)
			plt.grid()

			#Plot the cast info
			title = ''
			for iii in np.arange(temp.shape[0]):
				title = title+\
				'CAST '+str(iii)+': '+\
				ds.instrument_ID[spots].values[iii]+', '+\
				ds.comments[spots].values[iii]+', '+\
				ds.file_names[spots].values[iii]+', '+\
				ds.time[spots].values[iii].astype(str)+'\n'
			plt.suptitle(title,fontsize=8)

			#Determine the cast you would like to keep
			for test in range(1,10):
				user = input("Select profile number to keep (type 'all' to keep all).")
				while np.isin(user, np.concatenate((np.arange(spots.size),['all']))) == False:
					user = input('Not a valid choice. Please re-select a profile number.')
				if np.isin(user, np.arange(spots.size).astype(str)):
					#Record the kept cast
					x = int(user)
					duplicates_flag[spots[~np.isin(spots,spots[x])]] = True
					break
				elif user == 'all':
					break

			plt.close()
			print(str(ii)+' out of '+str(unique_dates.size)+' done.')

	#Remove the duplicates from the ds
	ds = ds.sel(time=~duplicates_flag)

	#Flatten the values inside the string variables
	for ii in ['trip_ID','comments','instrument_ID','instrument_type','file_names']:
		ds[ii] = ds[ii].astype(str)

	#Save the new file
	ds.to_netcdf(output_file,mode='w',\
		encoding={'time':{'units': "seconds since 1900-01-01 00:00:00",'calendar':	'gregorian'}})
	ds.close()

def remove_outliers(
	input_file,
	output_file,
	climatology,
	year,
	):
	'''
	input_file: path location to yearly netcdf file
	output_file: path location including name for netcdf file output
	climatology: path location to climatology files
	year (str): year of interest
	'''

	#Import the netcdf dataset
	ds = xr.open_dataset(input_file)

	#Determine the coresponding coordinates for each cast in the climatology
	#Define the lat/lon for the NAFC
	NAFC_lat = ds.latitude.values
	NAFC_lon = ds.longitude.values

	#Import the lat/lon from the climatology
	ds_clim = xr.open_dataset(climatology+'/bootstrap_climatology_01.nc')
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
		ds_clim = xr.open_dataset(climatology+'bootstrap_climatology_'+"%.2d" % (i)+'.nc')

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

	#Save the combined dataset
	#Flatten the arrays to ensure they save properly
	for i in ['trip_ID','comments','instrument_type','instrument_ID','file_names']:
		ds[i] = ds[i].astype(str)

	#Save the isolated netcdf files
	ds.to_netcdf(output_file,encoding={'time':{'units': 'seconds since 1900-01-01 00:00:00', 'calendar': 'gregorian'}})
	ds.close()





