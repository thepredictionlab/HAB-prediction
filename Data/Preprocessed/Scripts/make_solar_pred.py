# make_solar_pred.py - Data/Scripts/
# 
# Uses historical measurements to estimate solar activity at a given day/time.

import pandas as pd
from datetime import datetime
import numpy as np

def yearLength(yr):
    yr = np.reshape(yr,(-1,1))
    year_length = [datetime(year=yr[k]+1, month=1, day=1) - datetime(year=yr[k], month=1, day=1) for k in np.arange(len(yr))]
    year_length = [year_length[k].days for k in np.arange(len(yr))]
    return year_length

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
dataDir = fatherDir+'Data/Preprocessed/'
solar = pd.read_csv(dataDir+'data_solar_hist.csv')

date = solar['Date']
rad = solar['Radiation']
days = np.mod(date,1)

newdate = []
newrad = []

for year in np.arange(2018, 2028):
	yl = yearLength(year)[0]
	daysize = 1.0/yl
	hhsize = daysize/48 # length of 30 min window during that year
	hhest = daysize/36 # should be right at 48, but give a little leeway, 40 min window
	for day in np.arange(yl):
		for hh in np.arange(48):
			# given a particular day of a future year and a particular half hour
			# collect nearest half hours for days within 3 days of present
			# translated back over past years.  Average these to get radiation estimate.
			repdays = np.linspace(-3*daysize,3*daysize,7) + day*daysize + hh*hhsize
			repdays = np.mod(repdays,1) # in case translation wraps over the new year
			reps = np.zeros(len(days))
			for r in np.arange(len(repdays)):
				rep_r = np.where(np.abs(days - repdays[r]) < hhest)[0]
				reps[rep_r] = 1
			newdate.append(year + day*daysize + hh*hhsize)
			newrad.append(np.mean(rad[np.where(reps == 1)[0]]))

newsolar = np.hstack((newdate,newrad))
newsolar = np.reshape(newsolar,(-1,2),'F')
NewSolarData = pd.DataFrame(newsolar)
cols = ['Date','Radiation']
NewSolarData.to_csv(fatherDir+'Data/Preprocessed/data_solar_pred.csv',header=cols)







