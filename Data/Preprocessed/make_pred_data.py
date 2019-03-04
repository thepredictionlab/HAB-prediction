# make_pred_data.py

# exec(open('make_pred_data.py').read())

## loads most recent data from GitHub/Data/Preprocessed/Latest
## and, given an experiment (choice of posterior) compiles into 
## input data 

import pandas as pd
import numpy as np
import openpyxl as px
import sys
import os
import matplotlib.pylab as plt
from scipy import interpolate

# Define start and interval of time windows used in analysis
import datetime as dt
date = dt.datetime(2018,1,1,12,0,0)
timeInt = dt.timedelta(days = 7)

# Supporting function
def datetime2year(dat):
    year_part = dat - dt.datetime(year=dat.year, month=1, day=1)
    year_length = dt.datetime(year=dat.year+1, month=1, day=1) - dt.datetime(year=dat.year, month=1, day=1)
    return dat.year + year_part/year_length

## Create daily timeseries in decimal years
TIME = datetime2year(date) # this is the daily timeseries
for i in range(100): # range(2200):
    date += timeInt
    TIME = np.hstack((TIME,datetime2year(date)))

## Make an interpolating function
def interpolator(times, data, mode, domain):
    f = interpolate.interp1d(times,data,kind=mode,bounds_error=False,fill_value=np.nan)
    return f(domain)

def integrator(times, data, mode, domain):
    dataInt = np.zeros(len(domain))
    for k in np.arange(1,len(domain)):
        mask = (times >= domain[k-1])&(times < domain[k])
        data = np.array(data)
        if np.any(mask):
            print(k)
            if mode == 'mean':
                dataMean = np.nanmean(data[mask])
                dataInt[k] = dataMean
            elif mode == 'max':
                dataMax = np.nanmax(data[mask])
                dataInt[k] = dataMax
            elif mode == 'min':
                dataMin = np.nanmin(data[mask])
                dataInt[k] = dataMin
            elif mode == 'sum':
                dataSum = np.nansum(data[mask])
                dataInt[k] = dataSum
            else:
                raise Exception('Not an appropriate input for integrator.')
                dataInt[k] = np.nan
        else:
            dataInt[k] = np.nan
    return dataInt

def seasonalCleaner(times, data, funct, mode, domain):
    t = np.array(times)
    N = funct(times, data, mode, domain)
    yr = 2013
    while yr <2019:
        print(str(yr)+': '+str(np.where((t>=yr)&(t<yr+1))[0].shape[0]) + ' records.')
        if (np.where(t < yr+1)[0].shape[0] > 0)&(np.where(t >= yr+1)[0].shape[0] > 0):
            maxSample = np.max(t[t < yr+1])
            nextSample = np.min(t[t >= yr+1])
            yl = yearLength(yr)
            yl = yl[0]
            maxInt = domain[np.min(np.where(domain >= maxSample)[0])]
            nextYr = int(np.floor(nextSample))
            yl = yearLength(nextYr)
            yl = yl[0]
            nextInt = domain[np.min(np.where(domain >= nextSample)[0])]
            N[(domain > maxInt)&(domain < nextInt)] = np.nan
            yr = nextYr
        else:
            yr += 1
    return N


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
saveDir = fatherDir+'Data/Preprocessed/Latest/'

LOCS = ['BB','BO','HA','HT','LB','LBP','LBS','MG']
# Load all npz files of latest data
cyanoData = np.load(saveDir+'hist_cyano.npz')
nutsData = np.load(saveDir+'hist_nuts.npz')
algaeData = np.load(saveDir+'hist_algae.npz')
usgsData = np.load(saveDir+'hist_usgs.npz')
ysiData = np.load(saveDir+'hist_ysi.npz')
weatherData = np.load(saveDir+'hist_weather.npz')

# nut1 = nitrates, nut2 = ophos, nut3 = tn, nut4 = tphos

Time = nutsData['TIME']
loc = nutsData['LOCS']
nut1 = nutsData['NO']
nut2 = nutsData['OP']
nut3 = nutsData['TN']
nut4 = nutsData['TP']

# Interpolate to daily time series
# and rearrange so its timeseries of each nut for each location
NUT1 = np.ones((len(TIME),len(LOCS))) * np.nan
for i in np.arange(0,len(LOCS)): #1): #0,len(LOCS)):
    ID = np.where(loc == LOCS[i])[0]
    #ID = np.isin(loc, LOCS)
    if len(ID) > 1:
        n = nut1[ID]
        t = Time[ID] 
        print('Nutrient 1 (nitrate/nitrite):')
        N = seasonalCleaner(t,n,integrator,'mean',TIME)
        NUT1[:,i] = N 
    
NUT2 = np.ones((len(TIME),len(LOCS))) * np.nan
for i in np.arange(1): #0,len(LOCS)):
    #ID = np.where(loc == LOCS[i])[0]
    ID = np.isin(loc, LOCS)
    if len(ID) > 3:
        n = nut2[ID]
        t = Time[ID] 
        print('Nutrient 2 (O-Phos):')
        N = seasonalCleaner(t,n,integrator,'mean',TIME)
        NUT2[:,i] = N 
    
NUT3 = np.ones((len(TIME),len(LOCS))) * np.nan
for i in np.arange(1): #0,len(LOCS)):
    #ID = np.where(loc == LOCS[i])[0]
    ID = np.isin(loc, LOCS)
    if len(ID) > 1:
        n = nut3[ID]
        t = Time[ID] 
        print('Nutrient 3 (TN):')
        N = seasonalCleaner(t,n,integrator,'mean',TIME)
        NUT3[:,i] = N 
    
NUT4 = np.ones((len(TIME),len(LOCS))) * np.nan
for i in np.arange(1): #0,len(LOCS)):
    #ID = np.where(loc == LOCS[i])[0]
    ID = np.isin(loc, LOCS)
    if len(ID) > 1:
        n = nut4[ID]
        t = Time[ID] 
        print('Nutrient 4 (T-Phos):')
        N = seasonalCleaner(t,n,integrator,'mean',TIME)
        NUT4[:,i] = N 
    
time = cyanoData['TIME']
loc = cyanoData['LOCS']
tox3 = cyanoData['CYLIN']
tox4 = cyanoData['MICRO']

# Interpolate to daily time series
# and rearrange so its timeseries of each nut for each location
TOX3 = np.ones((len(TIME),len(LOCS))) * np.nan
for i in np.arange(0,len(LOCS)):#1): #0,len(LOCS)):
    ID = np.where(loc == LOCS[i])[0]
    #ID = np.isin(loc, LOCS)
    if len(ID)>1:
        n = tox3[ID]
        t = time[ID]
        print('ELISA Cylindro:')
        N = seasonalCleaner(t,n,integrator,'max',TIME)
        TOX3[:,i] = N

TOX4 = np.ones((len(TIME),len(LOCS))) * np.nan
for i in np.arange(0,len(LOCS)):
    ID = np.where(loc == LOCS[i])[0]
    #ID = np.isin(loc, LOCS)
    if len(ID)>1:
        n = tox4[ID]
        t = time[ID]
        print('ELISA Microcystin:')
        N = seasonalCleaner(t,n,integrator,'max',TIME)
        TOX4[:,i] = N

time = algaeData['TIME']
loc = algaeData['LOCS']
den = algaeData['DENSITY']
div = algaeData['DIV']
tbv = algaeData['TBV']
fbv = algaeData['FBV']

udiv = np.unique(div)
DIV = udiv

# Interpolate time series
# and rearrange so its timeseries of each div
DEN = np.ones((len(TIME),len(udiv))) * np.nan # DEN = np.ones((len(TIME),len(LOCS),len(udiv))) * np.nan
for i in np.arange(1): #0,len(LOCS)):
    #ID = np.where(loc == LOCS[i])[0]
    ID = np.isin(loc, LOCS)
    if len(ID)>1:
        for j in np.arange(0,len(udiv)):
            KD = np.where((div == udiv[j]) & ID)[0] # (loc == LOCS[i]))[0]
            if len(KD) > 1:
                n = den[KD]
                t = time[KD]
                n1 = []
                t1 = []
                tvals = np.unique(t)
                for s in np.arange(len(tvals)):
                    conc_div_samples = np.where(t == tvals[s])[0]
                    n1.append(np.nansum(n[conc_div_samples]))
                    t1.append(tvals[s])
                if (len(n1) != len(tvals))|(len(t1) != len(tvals)):
                    raise Exception("Cell density data compilation failed.")
                N = seasonalCleaner(t1,n1,integrator,'mean',TIME)
                DEN[:,j] = N

# TBV = np.ones((len(TIME),len(LOCS),len(udiv))) * np.nan
## make sure to sum across concurrent samples (same div, different genus)
## then average across locations.
TBV = np.ones((len(TIME),len(udiv))) * np.nan
MTBV = np.ones((len(TIME),len(udiv))) * np.nan
for i in np.arange(1): #0,len(LOCS)):
    #ID = np.where(loc == LOCS[i])[0]
    ID = np.isin(loc, LOCS)
    if len(ID)>1:
        for j in np.arange(0,len(udiv)):
            KD = np.where((div == udiv[j]) & ID)[0]#(loc == LOCS[i]))[0]
            if len(KD) > 1:
                n = tbv[KD]
                t = time[KD]
                n1 = []
                t1 = []
                tvals = np.unique(t)
                for s in np.arange(len(tvals)):
                    conc_div_samples = np.where(t == tvals[s])[0]
                    n1.append(np.nansum(n[conc_div_samples]))
                    t1.append(tvals[s])
                if (len(n1) != len(tvals))|(len(t1) != len(tvals)):
                    raise Exception("Biovolume data compilation failed.")
                print('Total biovolume of div '+str(j))
                N = seasonalCleaner(t1,n1,integrator,'mean',TIME)
                MN = seasonalCleaner(t1,n1,integrator,'max',TIME)
                #TBV[:,i,j] = N
                TBV[:,j] = N
                MTBV[:,j] = MN

FBV = np.ones((len(TIME),len(LOCS),len(udiv))) * np.nan
for i in np.arange(1): #0,len(LOCS)):
    #ID = np.where(loc == LOCS[i])[0]
    ID = np.isin(loc, LOCS)
    if len(ID)>1:
        for j in np.arange(0,len(udiv)):
            KD = np.where((div == udiv[j]) & ID)[0] # (loc == LOCS[i]))[0]
            if len(KD) > 1:
                n = fbv[KD]
                t = time[KD]
                n1 = []
                t1 = []
                tvals = np.unique(t)
                for s in np.arange(len(tvals)):
                    conc_div_samples = np.where(t == tvals[s])[0]
                    n1.append(np.nansum(n[conc_div_samples]))
                    t1.append(tvals[s])
                if (len(n1) != len(tvals))|(len(t1) != len(tvals)):
                    raise Exception("Fractional biovolume data compilation failed.")
                print('Frac biovolume of div '+str(j))
                N = seasonalCleaner(t1,n1,interpolator,interp_kind,TIME)
                f = interpolate.interp1d(t,n,kind=interp_kind,bounds_error=False,fill_value=np.nan)
                FBV[:,i,j] = N

time = weatherData['TIME']
tem = weatherData['TEMP']
tem = np.array([float(k) for k in tem])
tem = np.reshape(tem,(-1,1))
hum = weatherData['HUM']
hum = np.array([float(k) for k in hum])
hum = np.reshape(hum,(-1,1))
pwi = weatherData['PEAKW']
pwi = np.array([float(k) for k in pwi])
pwi = np.reshape(pwi,(-1,1))
wis = weatherData['WSPD']
wis = np.array([float(k) for k in wis])
wis = np.reshape(wis,(-1,1))
rain = weatherData['PRECIP']
rain = np.array([float(k) for k in rain])
rain = np.reshape(rain,(-1,1))
pres = weatherData['PRES']
pres = np.array([float(k) for k in pres])
pres = np.reshape(pres,(-1,1))

N = seasonalCleaner(time, tem, integrator, 'mean', TIME)
TEMP = N / 10

#N = seasonalCleaner(time, np.power(tem/10,2), integrator, 'mean', TIME)
#T2 = N

N = seasonalCleaner(time, tem, integrator, 'min', TIME)
MINT = N / 10

N = seasonalCleaner(time, tem, integrator, 'max', TIME)
MAXT = N / 10

f = interpolate.interp1d(time.flatten(),hum.flatten(),kind=interp_kind,bounds_error=False,fill_value=np.nan)
N = f(TIME)
JD = np.where(np.isnan(N)==0)[0]
KD = np.where(np.diff(N[JD])==0)    
N[JD[KD]] = np.nan
HUM = N

f = interpolate.interp1d(time.flatten(),pwi.flatten(),kind=interp_kind,bounds_error=False,fill_value=np.nan)
N = f(TIME)
JD = np.where(np.isnan(N)==0)[0]
KD = np.where(np.diff(N[JD])==0)
N[JD[KD]] = np.nan
PWI = N

f = interpolate.interp1d(time.flatten(),wis.flatten(),kind=interp_kind,bounds_error=False,fill_value=np.nan)
N = f(TIME)
JD = np.where(np.isnan(N)==0)[0]
KD = np.where(np.diff(N[JD])==0)
N[JD[KD]] = np.nan
WIS = N

# Make a degree day count - # of daily max temps exceeding a given value
tempThreshold = 65

## Get daily time series, divided at midnight
date = dt.datetime(2013,1,1,0,0,0)
timeInt = dt.timedelta(days = 1)
TIMED = datetime2year(date) # this is the daily timeseries
for i in range(2200):
    date += timeInt
    TIMED = np.hstack((TIMED,datetime2year(date)))

## Get daily max temp over previous day
N = seasonalCleaner(time, tem, integrator, 'max', TIMED)
DEG = [np.nanmax([N[k]-tempThreshold,0]) for k in np.arange(len(N))]
DDEG = np.zeros(len(DEG))
for yr in np.arange(2013,2020):
    days = (TIMED >= yr) & (TIMED < yr+1)
    DDEG[days] = np.cumsum(np.array(DEG)[days])

## Get weekly cumulative degree days
N = seasonalCleaner(TIMED,DDEG,integrator, 'sum', TIME)
WDEG = N
## And rescale
WDEG = WDEG/10000

N = seasonalCleaner(time,rain,integrator,'max',TIME)
RAIN = N

N = seasonalCleaner(time,pres,integrator,'mean',TIME)
PRES = N

# Calculate species diversity
#DEN = np.nansum(DEN,1)
PRCT = DEN/np.repeat(np.reshape(np.sum(DEN,1),(-1,1)), DEN.shape[1], 1)
LPRCT = np.log(PRCT)
LPRCT[np.where(LPRCT==-np.inf)] = 0
DVR = np.sum(PRCT * LPRCT,1)

# Calculate solar radiation (integrated over chosen intervals)

sol = pd.read_csv(fatherDir+'Data/Preprocessed/data_solar_hist.csv')
times = sol['Date']
rad = sol['Radiation']
psol = pd.read_csv(fatherDir+'Data/Preprocessed/data_solar_pred.csv')
ptimes = psol['Date']
prad = psol['Radiation']

times = np.vstack((np.reshape(times,(-1,1)),np.reshape(ptimes,(-1,1))))
rad = np.vstack((np.reshape(rad,(-1,1)),np.reshape(prad,(-1,1))))
N = seasonalCleaner(times, rad, integrator, 'sum', TIME)
RAD = N/10000


np.savez(fatherDir+"Data/Preprocessed/Data_latest_interp.npz",LOCS=LOCS,TIME=TIME,NUT1=NUT1,WDEG=WDEG,
        NUT2=NUT2,NUT3=NUT3,NUT4=NUT4,TOX3=TOX3,TOX4=TOX4,DIV=DIV,DEN=DEN,DVR=DVR,
        TBV=TBV,MTBV=MTBV,FBV=FBV,TEMP=TEMP,MINT=MINT,MAXT=MAXT,
        HUM=HUM,PWI=PWI,WIS=WIS,RAIN=RAIN,PRES=PRES,RAD=RAD)





