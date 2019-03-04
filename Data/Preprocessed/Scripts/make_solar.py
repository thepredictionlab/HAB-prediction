# make_solar.py - Data/Scripts/
# 
# https://nsrdb.nrel.gov/current-version

import pandas as pd
from datetime import datetime
import numpy as np

def datetime2year(dt):
    year_part = dt - datetime(year=dt.year, month=1, day=1)
    year_length = datetime(year=dt.year+1, month=1, day=1) - datetime(year=dt.year, month=1, day=1)
    return dt.year + year_part/year_length

# check for top directory level - 'HAB-prediction' used on GitHub, 'GitHub' used on my laptop.
# check for other top directories (containing 'BMA', 'Data' folders, etc.) by adding other try statements.
cpath = os.path.abspath('.').split('/')
try:
    fatherLevel = cpath.index('HAB-prediction')
except(ValueError):
    try:
        fatherLevel = cpath.index('GitHub')
    except(ValueError):
        raise Exception('Top directory not found.')

fatherDir = '/'.join(cpath[:fatherLevel+1]) + '/'
folderName = fatherDir+'/Data/Raw_solar/Historical/'
SolarData = np.reshape( np.array([np.nan,np.nan]), (1,2))
cols = []
for Y in np.arange(1998,2018):
	fileName = '205854_44.69_-122.18_' + str(Y) + '.csv'
	## Data
	df = pd.read_csv(folderName+fileName, header=2)
	dataH = df['DHI']
	dataN = df['DNI']
	# sum solar radiation from sun and non-sun sources:
	data = dataH + dataN
	## Dates
	date = []
	for k in np.arange(data.shape[0]):
		date.append( dt.datetime( year=df['Year'][k], month=df['Month'][k], day=df['Day'][k], \
			hour=df['Hour'][k], minute=df['Minute'][k] ) )
	deci_date = []
	for t in np.arange(len(date)):
		deci_date.append(datetime2year(date[t]))
	deci_date = np.reshape(deci_date,len(date),1)
	newyear = np.hstack((deci_date,data))
	newyear = np.reshape(newyear,(-1,2),'F')
	SolarData = np.vstack((SolarData, newyear))

SolarData = SolarData[np.abs(SolarData[:,0] - 2008)<12,:]
SolarData = pd.DataFrame( SolarData )
cols = ['Date','Radiation']
SolarData.to_csv(fatherDir+'Data/Preprocessed/data_solar_hist.csv',header=cols)

