# make_data_csv.py

import datetime as dt
import numpy as np
import pandas as pd

data_path = '../Data/Preprocessed/'
data_file = 'Data_hist_interp.npz'
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
# FBV: fractional biovolume
# TEMP: temperature
# HUM: humidity
# PWI: peak wind speed
# WIS: wind speed
# RAIN: rain
# PRES: barometric pressure
########################################################################

dates = data['TIME']
WVARS = ['TEMP','MAXT','MINT','WDEG','HUM','PWI','PRES','RAIN','DVR','WIS']


inclusionVector = [\
	[1,1,1,1,1,1,1,1,1,1,1,1,1], #TBV
	[1,1,1,1], # NUT
	[1,1,1,1], # TOX
	[1,1,1,1,1,1,1,1,1,1]] #WEA: TEMP, MAXT, MINT, WDEG, HUM, PWI, PRES, RAIN, DVR, WIS

# mask unwanted entries: 
inputs = {'abundances':{'TBV'}, #},#'TBV'}, # need loc and div
	'divisions':np.where(inclusionVector[0])[0], #[0,1,2,4], #[0,1,2],
	'locations':['BB','BO','HA','HT','LB','LBP','LBS'], #['LB','BO'],
	'nutrients':np.where(inclusionVector[1])[0], #{},#'NUT1','NUT2','NUT3','NUT4'},		#{}, #{'NUT2'},		# need loc
	'toxins':np.where(inclusionVector[2])[0], #{}, # 'TOX4'},
	'weather':np.where(inclusionVector[3])[0]} #{'TEMP','MAXT','RAIN'}}#,'WDEG'}} # 'TEMP','MINT','MAXT'}}

dataNames = [['baciTbv','chlorTbv','chrysTbv','cryptTbv','cyanoTbv',\
			'eugleTbv','haptTbv','hapt2Tbv','none','pyrrTbv',\
			'raphiTbv','rhodoTbv','xanthTbv'],\
			['NO','OP','TN','TP'],\
			['cyanoLC','microLC','cyanoEL','microEL'],\
			['tem','maxTem','minTem','degDays','hum','peakWind',\
			'pres','rain','shanDiv','windSpeed']]

# Assemble priors and types. Same length as inclusionVector, compile
# entries corresponding to nonzeros in incVec.
dataNameList = [];
for i1 in np.arange(len(inclusionVector)):
	for i2 in np.where(inclusionVector[i1])[0]:
		dataNameList.extend([dataNames[i1][i2]])

dataNameList = [k for sublist in dataNameList for item in [sublist if isinstance(sublist,list) else [sublist]] for k in item]

dataFull = np.array(dates)
dataNameList.insert(0,'date')
# This is for Total Biovolume only.  Adding other abundance measures will require a rewrite.
# All relevant location data is averaged to get TBV level at a given time.
for div in np.where(inclusionVector[0])[0]:
	divFull = np.zeros((data['TIME'].shape[0],1)) * np.nan
	#for loc in inputs['locations']:
	#	ind1 = data['LOCS'].tolist().index(loc)
	divFull = np.append(divFull, 
		np.reshape(data['TBV'][:,div],(divFull.shape[0],1)), 1)
		#np.reshape(data[dc][:,ind1,div],(divFull.shape[0],1)), 1)
	dataFull = np.append(dataFull, np.nanmean(divFull,1))

# Calculate the MEAN nutrient level across relevant sampling sites
for nut in np.where(inclusionVector[1])[0]:
	nutFull = np.zeros((data['TIME'].shape[0],1)) * np.nan
	for loc in inputs['locations']:
		nutFull = np.append(nutFull, data['NUT'+str(nut+1)][:,np.where(data['LOCS'] == loc)[0]], axis = 1)
	dataFull = np.append(dataFull, np.nanmean(nutFull,1))

# Collect MAX toxin level across sample sites
for tox in np.where(inclusionVector[2])[0]:
	item = 'TOX'+str(tox+1)
	dataFull = np.append(dataFull, np.nanmax(data[item],1))

# Collect weather data from detroit lake station
for wea in np.where(inclusionVector[3])[0]:
	item = WVARS[wea]
	dataFull = np.append(dataFull, data[item])
dataFull = dataFull.reshape((data['TEMP'].shape[0],-1), order = 'F')

print('dataFull shape = ' + str(dataFull.shape) + '.')

# save location
fileName = data_path + 'data_hist_matrix.csv'
# save dataFull in csv
df = pd.DataFrame(dataFull)
df.to_csv(fileName,header=dataNameList)



