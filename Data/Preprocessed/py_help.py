# py_help.py
## helper file for the BMA code base
## includes the sigma, quadratic, and timestamp functions
## defines the parameter types and their prior distributions across scripts

import os
import pandas as pd
import numpy as np
import openpyxl as px
import matplotlib.pylab as plt
import datetime as dt
from scipy import interpolate

EXP = 'PCTBV'

###############################################################################
#### Function for calculating decimal years
from datetime import datetime
def yearLength(yr):
    yr = np.reshape(yr,(-1,1))
    year_length = [datetime(year=yr[k]+1, month=1, day=1) - datetime(year=yr[k], month=1, day=1) for k in np.arange(len(yr))]
    year_length = [year_length[k].days for k in np.arange(len(yr))]
    return year_length

def datetime2year(dt):
    year_part = dt - datetime(year=dt.year, month=1, day=1)
    year_length = datetime(year=dt.year+1, month=1, day=1) - datetime(year=dt.year, month=1, day=1)
    return dt.year + year_part/year_length

def year2datetime(date):
    year = np.array(np.uint(np.floor(date)), ndmin=2)
    year_part = date - year
    Seconds = year_part*np.uint(yearLength(year))*24*60*60
    Date = [ dt.datetime(year=year[0][k] ,month=1,day=1) + dt.timedelta(seconds = Seconds[0][k]) for k in np.arange(year.size) ]
    return Date

def timestamp():
	now = dt.datetime.now()
	s = str(now.year) + str(now.month) + str(now.day) + "_" \
	    + str(now.hour) + str(now.minute) + str(now.second)
	return s

# for nut1 and 4, C = 0.055, L = 0.008
# for nut2 and 3, C = 0.025, L = 0.002
def sigmoid(x,L,C):
	S = np.exp((x-C)/(2*L))/(np.exp((x-C)/(2*L))+1)
	return S

# for temp max growth at 29 C = 84.5 F
# so use v = 8.5, c < 0, h s.t. Q(40 F) = 0
# also using nonnegative range cause cold is just cold
def quadratic(x,v,c):
	h = -(c*(4-v)**2)
	Q = c*np.power((x-v),2) + h
#	Q = [np.max(q,0) for q in Q]
	return Q

def getDist(data,bins):
	h = np.histogram(data,bins)
	return h[0]/np.sum(h[0])
###############################################################################

interp_kind = 'zero' # 'previous'
LOCS = ['BB','BO','HA','HT','LB','LBP','LBS','MG']

mSamples = 8000
mCores = 2
mTuning = 3000

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
dataDir = fatherDir + 'Data/Raw_lake/Historical/'

data_path = fatherDir + 'Data/Preprocessed/'
matrix_file = data_path + 'data_hist_matrix.csv'
latest_matrix_file = data_path + 'data_latest_matrix.csv'
df = pd.read_csv(matrix_file,prefix='X')
df2 = pd.read_csv(latest_matrix_file,prefix='X')
dates = df['date']
dates2 = df2['date']

# the folder holding the models predicting a given response variable
experiment = EXP+'/'

# Define predictand.
if (EXP == 'PCTBV') | (EXP=='TBVcat'):
	dataOut = df['cyanoTbv']
	predOut = np.diff(dataOut)/np.array(dataOut[:-1])
	predOut = np.insert(predOut,len(predOut),np.nan)
	pBins = np.linspace(-2,2,20)
	pBins = np.hstack((-np.inf, pBins, np.inf))
elif EXP == 'PCTBVS':
	dataOut = df['cyanoTbv']
	predOut = np.diff(dataOut)/np.array(dataOut[:-1])
	predOut = np.insert(predOut,len(predOut),np.nan)
	predOut[np.where(predOut < 0)[0]] = -1
	predOut[np.where(predOut > 0)[0]] = 1
	predOut[np.where(predOut==0)[0]] = np.nan
	predOut = predOut[1:]
	predOut = np.insert(predOut,len(predOut),np.nan)
	pBins = np.zeros(1)
	pBins = np.hstack((-np.inf, pBins, np.inf))
elif EXP == 'TBVZ':
	dataOut = np.array(df['cyanoTbv'])
	dataOut = (dataOut - np.nanmean(dataOut))/np.nanstd(dataOut)
	predOut = np.zeros(dataOut.shape)
	predOut[dataOut < 1.25] = -1
	predOut[dataOut >=1.25] = +1
	predOut = predOut[1:]
	predOut = np.insert(predOut,len(predOut),np.nan)
	pBins = np.zeros(1)
	pBins = np.hstack((-np.inf, pBins, np.inf))
	presDataOut = np.array(df2['cyanoTbv'])
	presDataOut = (presDataOut - np.nanmean(presDataOut))/np.nanstd(presDataOut)
	presPredOut = np.zeros(presDataOut.shape)
	presPredOut[presDataOut < 1.25] = -1
	presPredOut[presDataOut >=1.25] = +1
	presPredOut = presPredOut[1:]
	presPredOut = np.insert(presPredOut,len(presPredOut),np.nan)
elif EXP == 'TOXLVL':
	# POSSIBLE DEFINITION OF TOXIN DATA: MAX OF CYLINDRO AND 2*MICRO.
	dataOut = np.array(df['microEL'])
	predOut = np.zeros(dataOut.shape)
	# safe level
	predOut[dataOut < np.log(0.3)] = 0
	# unsafe for vulnerable pop
	predOut[(dataOut >= np.log(0.3)) & (dataOut < np.log(1.6))] = 1
	# unsafe for adults
	predOut[(dataOut >= np.log(1.6)) & (dataOut < np.log(4))] = 2
	# unsafe for recreation
	predOut[dataOut >= np.log(4)] = 3
	predOut = predOut[1:]
	predOut = np.insert(predOut,len(predOut),np.nan)
	pBins = np.log(np.array([0.3, 1.6, 4]))
	pBins = np.hstack((-np.inf, pBins, np.inf))
elif EXP == 'TOXC':
	dataOut = np.array(df['microEL'])
	predOut = np.zeros(dataOut.shape)
	predOut[dataOut < np.log(0.3)] = -1
	predOut[dataOut >= np.log(0.3)] = +1
	predOut = predOut[1:]
	predOut = np.insert(predOut,len(predOut),np.nan)
	pBins = np.zeros(1)
	pBins = np.hstack((-np.inf, pBins, np.inf))
else:
	raise Exception('The distribution of the posterior must be defined for EXP = '+EXP)


#### TRANSFORM predOut TO FORM OF POSTERIOR DISTRIBUTION ####
# This way we train the model against the actual 'categories',
# rather than finer resolution than we will check for at the 
# end of the algorithm
#############################################################


# F # I # L # L # # I # N #


def defns(mode):
	if mode == 'h':
		inclVec = [np.ones(13).tolist(),np.ones(13).tolist(),np.ones(4).tolist(),np.ones(4).tolist(),np.ones(11).tolist()]
	elif mode == 'l':
		inclVec = [np.ones(13).tolist(),np.ones(13).tolist(),np.ones(4).tolist(),np.ones(4).tolist(),np.ones(11).tolist()]
		inclVec[0] = [1,1,1,1,1,1,0,0,0,1,0,0,0]
		inclVec[1] = [1,1,1,1,1,1,0,0,0,1,0,0,0]
		inclVec[3] = [0,0,1,1]
	else:
		raise Exception('Not a valid mode for defns() in py_help.')
	paramVector = [\
		['l','l','l','l','l','l','l','l','l','l','l','l','l'], # BIO
		['l','l','l','l','l','l','l','l','l','l','l','l','l'], # MBIO
		['s','s','s','s'], # NUT
		['l','l','l','l'], # TOX
		['q','q','q','l','l','l','l','l','l','l','l']] # WEA
	dataNames = [['baciTbv','chlorTbv','chrysTbv','cryptTbv','cyanoTbv',\
				'eugleTbv','haptTbv','hapt2Tbv','none','pyrrTbv',\
				'raphiTbv','rhodoTbv','xanthTbv'],\
				['baciMTbv','chlorMTbv','chrysMTbv','cryptMTbv','cyanoMTbv',\
				'eugleMTbv','haptMTbv','hapt2MTbv','noneM','pyrrMTbv',\
				'raphiMTbv','rhodoMTbv','xanthMTbv'],\
				['NO','OP','TN','TP'],\
				['cyanoLC','microLC','cyanoEL','microEL'],\
				['tem','maxTem','minTem','degDays','hum','peakWind',\
				'pres','rain','shanDiv','windSpeed','rad']]
	WVARS = ['TEMP','MAXT','MINT','WDEG','HUM','PWI','PRES','RAIN','DVR','WIS']
	# Define priors for all variables
	priorMeans = [\
		[-1,-1,-1,-1,+1,-1,-1,-1,-1,-1,-1,-1,-1], # BIO
		[-1,-1,-1,-1,+1,-1,-1,-1,-1,-1,-1,-1,-1], # MBIO
		[[0.008, 0.055],[0.002, 0.025],[0.002, 0.025],[0.008, 0.055]], # NUT
		[1,1,1,1], # TOX
		[[8.3,-0.1],[8.5,-0.3],[8.5,-0.3],1,0,-1,0,1,0,-1,1] # WEA
		]
	priorStds = [\
		[1,1,1,1,1,1,1,1,1,1,1,1,1], # BIO
		[1,1,1,1,1,1,1,1,1,1,1,1,1], # MBIO
		#[10,10,10,10,10,10,10,10,10,10,10,10,10], # BIO
		[[0.1,0.1],[0.1,0.1],[0.1,0.1],[0.1,0.1]], # NUT
		[2,2,2,2], # TOX
		[[1,0.5],[2,0.5],[2,0.5],1,1,1,1,1,1,1,1] # WEA
		]
	paramType = []; paramList = [];
	muList = []; stdList = []; dataNameList = [];
	for i1 in np.arange(len(inclVec)):
		paramType.append([]);
		for i2 in np.where(inclVec[i1]!=0)[0]:
			paramType[i1].append(paramVector[i1][i2])
			paramList.extend(paramVector[i1][i2])
			muList.extend([priorMeans[i1][i2]])
			stdList.extend([priorStds[i1][i2]])
			dataNameList.extend([dataNames[i1][i2]])
	paramList = [k for sublist in paramList for item in [sublist if isinstance(sublist,list) else [sublist]] for k in item]
	dataNameList = [k for sublist in dataNameList for item in [sublist if isinstance(sublist,list) else [sublist]] for k in item]
	muList = [k for sublist in muList for item in [sublist if isinstance(sublist,list) else [sublist]] for k in item]
	stdList = [k for sublist in stdList for item in [sublist if isinstance(sublist,list) else [sublist]] for k in item]
	return paramVector,dataNames,WVARS,priorMeans,priorStds,paramType,paramList,muList,stdList,dataNameList

#paramVector,dataNames,WVARS,priorMeans,priorStds,paramType,paramList,muList,stdList,dataNameList = defns('h')











################################################################ SAVE ########
# Locations: LOCS = ['BB','BO','HA','HT','LB','LBP','LBS']
# TIME: decimal years (daily increments, at midday)
# NUT1: NO3+NO2 (mg/L)                                  sparse
# NUT2: O-Phos (mg/L)                                   sparse
# NUT3: TN (mg/L)                                       sparse
# NUT4: T-Phos (mg/L)                                   sparse
# TOX1: LCMSMS Cylindro (ppb)                           sparse
# TOX2: LCMSMS Microcystin (ppb)
# TOX3: ELISA Cylindro (ppb)                            sparse
# TOX4: ELISA Microcystin (ppb)
# DIV: bacterial/algal family name
# DEN: algal concentration
# TBV: total biovolume
# FBV: fractional biovolume
# TEMP: temperature 
# T2: temp squared
# MINT: minimum temperature over the period
# MAXT: maximum temperature over the period
# HUM: humidity
# PWI: peak wind speed
# WIS: wind speed
# RAIN: rain
# PRES: barometric pressure
# WDEG: degree days integrated over each period
# DVR: shannon diversity of algal populations           sparse
# RAD: integrated solar radiation around the lake
################################################################ SAVE ########



