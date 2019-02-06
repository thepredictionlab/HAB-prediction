# make_something.py

import matplotlib.pyplot as plt
#plt.style.use(['seaborn-darkgrid'])
import pymc3 as pm
import numpy as np

folderName = './Results/'
data_path = '../Data/Preprocessed/'
data_file = 'Data_historical.npz'

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

# mask unwanted entries: 
inputs = {'abundances':{'TBV'}, # need loc and div
	'divisions':[4],
	'locations':['LB','BO'],
	'nutrients':{}, #{'NUT2'},		# need loc
	'weather':{'TEMP'}}	

inputOutput = inputs.copy()
inputOutput['toxins'] = ['TOX3']

dataFull = np.array([])
dataList = []
for dc in inputs['abundances']:
	for loc in inputs['locations']:
		for div in inputs['divisions']:
			ind1 = inputs['locations'].index(loc)
			dataFull = np.append(dataFull, data[dc][:,ind1,div])

for nut in inputs['nutrients']:
	for loc in inputs['locations']:
		dataFull = np.append(dataFull, data[nut][:,np.where(data['LOCS'] == loc)[0]])

for item in inputs['weather']:
	dataFull = np.append(dataFull, data[item])

dataFull = dataFull.reshape((data['TEMP'].shape[0],-1), order = 'F')
print('dataFull shape = ' + str(dataFull.shape) + '.')

# Data matrix made. Filter for appropriate dates. Process.

# dates/times with NaNs
nanDays = np.any(np.isnan(dataFull),1)
domain = ~nanDays
print('Domain consists of ' + str(dataFull[domain,0].shape[0]) + ' days.')

# response variable
y1 = np.array( data['TOX3'][:,4] )

# initiate model
# assume linear response to the variables, normal priors for coefficients

mSamples = 2000
mCores = 4
mTuning = 500
with pm.Model() as model:
	norm100 = pm.Normal('n', mu=0, sd=100) # too broad?
	coeffs = pm.Normal('c', mu=0, sd=10, shape=dataFull.shape[1])
	mu = 0
	for k in np.arange(dataFull.shape[1]):
		mu += dataFull[domain,k]*coeffs[k]
		#print('coeff'+str(k)+' = '+dataList[k])
	eta = pm.Normal('noise', mu=norm100, sd=10)
	#eta2 = pm.AR('ar_noise',eta,sd = 1.0)
	obs = pm.Normal('obs', mu = mu+eta, sd = 1, observed = y1[domain])
	trace = pm.sample(5000, cores = 4, tuning = 500)
saveTraceName = './trace'
np.save(trace, saveTraceName)

pm.traceplot(trace)
savePlotName = folderName+str(mSamples)+'cores'+str(mCores)+'tuning'+str(mTuning)
plt.savefig(savePlotName)

plt.show()

ppc_w = pm.sample_posterior_predictive_w(traces, 1000, models,
                        weights=comp.weight.sort_index(ascending=True),
                        progressbar=False)