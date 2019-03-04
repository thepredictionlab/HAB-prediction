# make_data_csv.py

# exec(open('make_data_csv.py').read())

import datetime as dt
import numpy as np
import pandas as pd
from py_help import *
import sys
paramVector,dataNames,WVARS,priorMeans,priorStds,paramType,paramList,muList,stdList,dataNameList = defns('l')
# dataNames[0] = np.array(dataNames[0])[[True,True,True,True,True,True,False,False,False,True,False,False,False]].tolist()
# dataNames[1] = np.array(dataNames[1])[[True,True,True,True,True,True,False,False,False,True,False,False,False]].tolist()
# dataNames[3] = np.array(dataNames[3])[[False,False,True,True]].tolist()
#minTox = 0
minTox = 2
maxTox = 4
data_file = 'Data_latest_interp.npz'
data = np.load(data_path + data_file)

########################################################################
# Locations: LOCS = ['BB','BO','HA','HT','LB','LBP','LBS']
# TIME: decimal years (daily increments, at midday)
# NUT1: NO3+NO2 (mg/L)
# NUT2: O-Phos (mg/L)
# NUT3: TN (mg/L)
# NUT4: T-Phos (mg/L)
# TOX1: LCMSMS Cylindro (ppb)
# TOX2: LCMSMS Microcystin (ppb)
# TOX3: ELISA Cylindro (ppb)
# TOX4: ELISA Microcystin (ppb)
# DIV: bacterial/algal family name
# DEN: algal concentration
# TBV: total biovolume
# MTBV: temporal max of total biovolume
# FBV: fractional biovolume
# TEMP: temperature
# HUM: humidity
# PWI: peak wind speed
# WIS: wind speed
# RAIN: rain
# PRES: barometric pressure
########################################################################

dates = data['TIME']
WVARS = ['TEMP','MAXT','MINT','WDEG','HUM','PWI','PRES','RAIN','DVR','WIS','RAD']

# Assemble priors and types. Same length as inclusionVector, compile
# entries corresponding to nonzeros in incVec.
dataNameList = [];
for i1 in np.arange(len(dataNames)):
	for i2 in np.where(dataNames[i1])[0]:
		dataNameList.extend([dataNames[i1][i2]])

dataNameList = [k for sublist in dataNameList for item in [sublist if isinstance(sublist,list) else [sublist]] for k in item]

dataFull = np.array(dates)
dataNameList.insert(0,'date')

# This is for Total Biovolume only.  Adding other abundance measures will require a rewrite.
# All relevant location data is averaged to get TBV level at a given time.
for div in np.arange(len(dataNames[0])):
	divFull = np.zeros((data['TIME'].shape[0],1)) * np.nan
	#for loc in inputs['locations']:
	#	ind1 = data['LOCS'].tolist().index(loc)
	divFull = np.append(divFull, 
		np.reshape(data['TBV'][:,div],(divFull.shape[0],1)), 1)
		#np.reshape(data[dc][:,ind1,div],(divFull.shape[0],1)), 1)
	dataFull = np.append(dataFull, np.nanmean(divFull,1))
# MAX TBV data
for div in np.arange(len(dataNames[1])):
	divFull = np.zeros((data['TIME'].shape[0],1)) * np.nan
	#for loc in inputs['locations']:
	#	ind1 = data['LOCS'].tolist().index(loc)
	divFull = np.append(divFull, 
		np.reshape(data['MTBV'][:,div],(divFull.shape[0],1)), 1)
		#np.reshape(data[dc][:,ind1,div],(divFull.shape[0],1)), 1)
	dataFull = np.append(dataFull, np.nanmean(divFull,1))

# Calculate the MEAN nutrient level across relevant sampling sites
for nut in np.arange(len(dataNames[2])):
	nutFull = np.zeros((data['TIME'].shape[0],1)) * np.nan
	for loc in LOCS:
		nutFull = np.append(nutFull, data['NUT'+str(nut+1)][:,np.where(data['LOCS'] == loc)[0]], axis = 1)
	dataFull = np.append(dataFull, np.nanmean(nutFull,1))

# Collect MAX toxin level across sample sites
for tox in np.arange(minTox,maxTox):
	item = 'TOX'+str(tox+1)
	dataFull = np.append(dataFull, np.nanmax(data[item],1))

# Collect weather data from detroit lake station
for wea in np.arange(len(dataNames[4])):
	item = WVARS[wea]
	dataFull = np.append(dataFull, data[item])
dataFull = dataFull.reshape((data['TEMP'].shape[0],-1), order = 'F')

print('dataFull shape = ' + str(dataFull.shape) + '.')

# save location
fileName = data_path + 'data_latest_matrix.csv'
# save dataFull in csv
df = pd.DataFrame(dataFull)
df.to_csv(fileName,header=dataNameList)

sys.exit()

