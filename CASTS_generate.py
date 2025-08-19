import glob
import os
import re
import warnings
from scipy import stats
import numpy as np
import xarray as xr
import netCDF4 as nc
import time as tt
import pandas as pd

#Provide the path to where custom packages are saved
import sys
sys.path.append('/fs/vnas_Hdfo/dpnm/joc000/Documents/Scripts/CASTS/Modules')
import NAFC_Oceanography_package as NAFC_oceanography
import CIOOS_NAFC_package as CIOOS_NAFC
import NAFC_Aquaculture_package as NAFC_aquaculture
import IML_package as IML
import BIO_Climate_package as BIO_climate
import CIOOS_ERDDAP_package as CIOOS_ERDDAP
import IEO_Spain_package as IEO_Spain
import NCEI_GTSPP_package as NCEI_GTSPP
import NEFSC_package as NEFSC
import Polar_Data_Catalogue_package as Polar_Data_Catalogue
import WOD_package as WOD
import CASTS_combined_package as CASTS_combined

'''
This is the central CASTS script.
All files will be processed here by calling to custom-made modules.
Jonathan Coyne, May 2024
'''

##########################################################################################################
#STANDARD INFO
max_depth = 5000
directory = '/gpfs/fs7/dfo/dpnm/joc000/Data/CASTS/'


##########################################################################################################
'''----------------------------------------------------------------------------------------------------'''
##########################################################################################################
#1. NAFC-OCEANOGRAPHY
#Update yearly

#Define the years of interest
years = np.arange(1912,2023+1).astype(str)

#Cycle through each year
for year in years[:]:

	#Convert to raw netcdf files
	input_folder = directory+'Data_Input/NAFC_Oceanography/pfiles/'
	output_folder = directory+'Data_Input/NAFC_Oceanography/1_files_raw/'
	if int(year) < 2023:
		NAFC_oceanography.convertpfile(input_folder+year+'/',output_folder+year+'.nc',year)
	else:
		NAFC_oceanography.convertcnv(input_folder,output_folder+year+'.nc',year)
	print(year+' -> raw file done!')

	#Remove empties, problem MBT BO profiles
	input_file = directory+'Data_Input/NAFC_Oceanography/1_files_raw/'
	output_file = directory+'Data_Input/NAFC_Oceanography/2_empties/'
	flagged_files = directory+'Data_Input/NAFC_Oceanography/2_empties/flagged_files/'
	NAFC_oceanography.remove_empties(input_file+year+'.nc',output_file+year+'.nc',flagged_files,year)
	print(year+' -> removed empties, problem MBT BO done!')

	#Remove duplicates
	input_file = directory+'Data_Input/NAFC_Oceanography/2_empties/'
	output_file = directory+'Data_Input/NAFC_Oceanography/3_duplicates/'
	NAFC_oceanography.remove_duplicates(input_file+year+'.nc',output_file+year+'.nc',year)
	print(year+' -> removed duplicates done!')

	#Remove outliers
	#Assumes that climatology has already been completed
	input_file = directory+'Data_Input/NAFC_Oceanography/3_duplicates/'
	output_file = directory+'Data_Input/NAFC_Oceanography/4_outliers/'
	climatology = directory+'Data_Input/NAFC_Oceanography/4_outliers/climatology/'
	NAFC_oceanography.remove_outliers(input_file+year+'.nc',output_file+year+'.nc',climatology,year)
	print(year+' -> removed outliers done!')



##########################################################################################################
#1. CIOOS-NAFC
#Update yearly

#Download all files
file_input = directory+'Data_Input/CIOOS_NAFC/data_raw/'
years = np.arange(1983,2024+1).astype(str)
years_azmp = np.arange(1999,2024+1).astype(str)
CIOOS_NAFC.download_AZMP(file_input,years_azmp)
years_multispecies = np.arange(1995,2024+1).astype(str)
CIOOS_NAFC.download_multispecies(file_input,years_multispecies)
years_stn27 = np.arange(1983,2024+1).astype(str)
CIOOS_NAFC.download_stn27(file_input,years_stn27)
years_unsorted = np.arange(1983,2024+1).astype(str)
CIOOS_NAFC.download_unsorted(file_input,years_unsorted)
years_NSRF = np.arange(2005,2024+1).astype(str)
CIOOS_NAFC.download_NSRF(file_input,years_NSRF)

#Mesh individual files together
file_output = directory+'Data_Input/CIOOS_NAFC/data_processed/'
CIOOS_NAFC.merge_netcdf(file_input,file_output,years)
print('CIOOS-NAFC -> yearly netcdf done!')



##########################################################################################################
#2. NAFC-AQUACULTURE
#Update yearly (with info from Andry Ratsimandresy)

#Define the years of interest
years = np.arange(2009,2020+1)

#Convert folders to individual netcdf files
input_path = directory+'Data_Input/NAFC_Aquaculture/data_raw/'
NAFC_aquaculture.raw_to_netcdf(input_path)
print('NAFC Aquaculture -> raw to individual netcdf done!')

#Convert to yearly files
file_output = directory+'Data_Input/NAFC_Aquaculture/data_processed/'
NAFC_aquaculture.netcdf_gather(input_path,file_output,years)
print('NAFC Aquaculture -> yearly netcdf done!')



##########################################################################################################
#3. IML
#Update yearly (with info from Jean-Luc)

#Convert to yearly files
file_input = directory+'Data_Input/IML/data_raw/2025_Aug/'
file_output = directory+'Data_Input/IML/data_processed/'
IML.raw_to_netcdf(file_input,file_output)
print('IML -> yearly netcdf done!')



##########################################################################################################
#4. BIO-Climate
#Historical, no update required

#Convert historical (pre-2010) Climate
file_input = directory+'Data_Input/BIO_Climate/R_Files/'
file_name = 'BIO_Climate_Databases'
file_output = directory+'Data_Input/BIO_Climate/NetCDF/'
BIO_climate.rdata_to_netcdf_pre2010(file_input+file_name+'.RData',file_name,file_output)
print('BIO Climate Historical -> yearly netcdf done!')

#Convert newer (post-2010) Climate
file_input = directory+'Data_Input/BIO_Climate/R_Files/'
file_name = 'database_2008-2017'
file_output = directory+'Data_Input/BIO_Climate/NetCDF/'
BIO_climate.rdata_to_netcdf_post2010(file_input+file_name+'.RData',file_name,file_output)
file_name = 'database_2018'
BIO_climate.rdata_to_netcdf_post2010(file_input+file_name+'.RData',file_name,file_output)
print('BIO Climate post-2010 -> yearly netcdf done!')



##########################################################################################################
#5. BIO CIOOS-ERDDAP
#Update yearly

#Download all files
file_input = directory+'Data_Input/CIOOS_ERDDAP/data_raw/Jan_2025/'
years = np.arange(1993,2024+1).astype(str)
years_azmp = np.arange(1997,2024+1).astype(str)
CIOOS_ERDDAP.download_AZMP(file_input,years_azmp)
years_azomp = np.concatenate([np.arange(1993,2016+1),np.arange(2018,2020+1),np.arange(2022,2024+1)]).astype(str)
CIOOS_ERDDAP.download_AZOMP(file_input,years_azomp)
years_ecosystems = np.arange(1996,2024+1).astype(str)
CIOOS_ERDDAP.download_ecosystems(file_input,years_ecosystems)

#Mesh individual files together
file_output = directory+'Data_Input/CIOOS_ERDDAP/data_processed/'
CIOOS_ERDDAP.merge_netcdf(file_input,file_output,years)
print('CIOOS-ERDDAP -> yearly netcdf done!')



##########################################################################################################
#6. IEO-Spain
#Update yearly (depending on available data)

#Convert from .dat files to netcdf
file_input = directory+'Data_Input/IEO_Spain/data_raw/'
file_output = directory+'Data_Input/IEO_Spain/data_processed/'
file_yearly = directory+'Data_Input/IEO_Spain/netcdf_yearly/'
IEO_Spain.raw_to_netcdf(file_input,file_output)
IEO_Spain.netcdf_yearly(file_output,file_yearly)
print('IEO Spain -> yearly netcdf done!')



##########################################################################################################
#7. NCEI-GTSPP
#Update yearly

#Download the data of interest
file_input = directory+'Data_Input/NCEI_GTSPP/data_raw/'
year = '2024'
NCEI_GTSPP.download_data(year,file_input)

#Merge individual files
file_output = directory+'Data_Input/NCEI_GTSPP/data_processed/'
NCEI_GTSPP.merge_netcdf(file_input,file_output,year)
#Depth check
file_depthcheck = directory+'Data_Input/NCEI_GTSPP/data_processed/bottom_flagged/'
NCEI_GTSPP.depth_check(file_output,file_depthcheck,year)
print('NCEI_GTSPP -> yearly netcdf done!')



##########################################################################################################
#8. NEFSC
#Update yearly (there's a CIOOS page)

#Download the data
file_input = directory+'Data_Input/NEFSC/ERDDAP/'
NEFSC.download_NEFSC(file_input)

#Re-format the downloaded data
file_output = directory+'Data_Input/NEFSC/individual_netcdf_files/'
file_yearly_output = '/gpfs/fs7/dfo/dpnm/joc000/Data/CASTS/Data_Input/NEFSC/yearly_netcdf_files/'
NEFSC.raw_to_netcdf(file_input+'NEFSC_ERDDAP.nc',file_output,file_yearly_output)
print('NEFSC -> yearly netcdf done!')



##########################################################################################################
#9. Polar Data Catalogue
#Update yearly (there's a website)

#Convert the folders to netcdf
file_input = directory+'Data_Input/Polar_Data_Catalogue/data_raw/'
file_output = directory+'Data_Input/Polar_Data_Catalogue/data_product/netcdf_group/'
file_yearly_output = directory+'Data_Input/Polar_Data_Catalogue/data_product/netcdf_yearly/'
Polar_Data_Catalogue.raw_to_netcdf(file_input,file_output)
Polar_Data_Catalogue.yearly_netcdf(file_input,file_output,file_yearly_output)
print('Polar Data Catalogue -> yearly netcdf done!')



##########################################################################################################
#10. World Ocean Database
#Update yearly (there's a website)

#Convert the folders to netcdf
file_input = directory+'Data_Input/Data_Input/WOD/OSD_MBT/'
file_output = directory+'Data_Input/WOD/data_processed/'
WOD.merge_netcdf(file_input,file_output)
print('World Ocean Database (WOD) -> yearly netcdf done!')


##########################################################################################################
'''----------------------------------------------------------------------------------------------------'''
##########################################################################################################
#Combine all sources together

#Bring all the sources together
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
path['path13'] = directory+'Data_Input/Other/Freds_Files/netcdf_yearly/' #Marine Institute, 2000-2008
path['path14'] = directory+'Data_Input/WOD/data_processed/' #World Ocean Database, 1873-1910
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
path_source['path13'] = 'Marine-Institute-NL'
path_source['path14'] = 'WOD'
years=np.arange(1873,2024+1).astype(str)
file_output=directory+'Data_Products/1_combined_raw/'
CASTS_combined.combine_netcdf(path,path_source,years,file_output)

#Isolate the area and finalize vertical resolution
file_input=directory+'Data_Products/1_combined_raw/'
file_output=directory+'Data_Products/2_restrict_region/'
CASTS_combined.area_isolate(file_input,file_output,years)

#Do a depth check
file_input=directory+'Data_Products/2_restrict_region/'
file_output=directory+'Data_Products/3_depth_check/'
CASTS_combined.depth_filter(file_input,file_output,years)

#Do an outlier check
file_input=directory+'Data_Products/3_depth_check/'
file_output=directory+'Data_Products/4_outlier_check/'
climatology_path=directory+'Data_Products/climatology/'
CASTS_combined.outlier_check(file_input,file_output,climatology_path,years)

#Do a station check
file_input=directory+'Data_Products/4_outlier_check/'
file_output=directory+'Data_Products/5_final_product/'
station_path=directory+'Data_Input/Other/station_coordinates/'
CASTS_combined.station_check(file_input,file_output,station_path,years)
print('CASTS Combined -> yearly netcdf done!')


