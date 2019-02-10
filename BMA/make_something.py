# make_something.py

import matplotlib.pyplot as plt
#plt.style.use(['seaborn-darkgrid'])
import pymc3 as pm
import numpy as np

folderName = './Results/'
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

# mask unwanted entries: 
inputs = {'abundances':{},#'TBV'}, # need loc and div
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
			ind1 = data['LOCS'].tolist().index(loc)
			dataFull = np.append(dataFull, data[dc][:,ind1,div])

for nut in inputs['nutrients']:
	for loc in inputs['locations']:
		dataFull = np.append(dataFull, data[nut][:,np.where(data['LOCS'] == loc)[0]])

for item in inputs['weather']:
	dataFull = np.append(dataFull, data[item])

dataFull = dataFull.reshape((data['TEMP'].shape[0],-1), order = 'F')
#dataFull[7:] = dataFull[:-7]
#dataFull[:7] = np.reshape(np.zeros(7),(7,1))
print('dataFull shape = ' + str(dataFull.shape) + '.')

# Filter for appropriate dates. Process.

## Validation via LOO:

# dates/times with NaNs
nanDays = np.any(np.isnan(dataFull),1) | np.isnan(dataOut)
domain = ~nanDays
print('Domain consists of ' + str(dataFull[domain,0].shape[0]) + ' days.')

# Data matrix made. Make output matrix.

dataOut = data['TOX3'][domain,4]
y1 = np.array( dataOut )

for t in np.arange(len(domain)):
	subData = dataFull[TIMES THAT AREN'T t,:]
	subOut = y1[SAME TIMES]
	linModel = pm.Model()
	with linModel as model:
		norm100 = pm.Normal('n', mu=0, sd=1) # too broad?
		coeffs = pm.Normal('c', mu=0, sd=1, shape=dataFull.shape[1])
		mu = 0
		for k in np.arange(dataFull.shape[1]):
			mu += dataFull[domain,k]*coeffs[k]
			#print('coeff'+str(k)+' = '+dataList[k])
		eta = pm.Normal('noise', mu=norm100, sd=10)
		#eta2 = pm.AR('ar_noise',eta,sd = 1.0)
		obs = pm.Normal('obs', mu = mu+eta, sd = 1, observed = y1)#[domain])
		trace = pm.sample(mSamples, cores = 4, tuning = mTuning)

	# use posterior to estimate response for left out time
	testData = dataFull[t,:]
	validation = y1[t]
	testResults = np.zeros(mSamples)
	for k in np.arange(mSamples):
		output = np.dot( trace['mu'][k], testData ) + trace['eta'][k] 
		# OR WHATEVER, inner product of coefficients mu with data
		testResults[k] = output
	plt.hist(testResults)
	plt.title(str(validation))
	plt.show()

# response variable
y1 = np.array( data['TOX3'][domain,4] )

# initiate model
# assume linear response to the variables, normal priors for coefficients

mSamples = 100
mCores = 2
mTuning = 50
linModel = pm.Model()
with pm.Model() as model:
	norm100 = pm.Normal('n', mu=0, sd=1) # too broad?
	coeffs = pm.Normal('c', mu=0, sd=1, shape=dataFull.shape[1])
	mu = 0
	for k in np.arange(dataFull.shape[1]):
		mu += dataFull[domain,k]*coeffs[k]
		#print('coeff'+str(k)+' = '+dataList[k])
	eta = pm.Normal('noise', mu=norm100, sd=10)
	#eta2 = pm.AR('ar_noise',eta,sd = 1.0)
	obs = pm.Normal('obs', mu = mu+eta, sd = 1, observed = y1)#[domain])
	trace = pm.sample(mSamples, cores = 4, tuning = mTuning)

saveTraceName = './trace.npz'
np.savez(saveTraceName,TRACE=trace)

pm.summary(trace)
pm.traceplot(trace)
plt.show()

savePlotName = folderName+str(mSamples)+'cores'+str(mCores)+'tuning'+str(mTuning)
plt.savefig(savePlotName)

plt.show()

ppc_w = pm.sample_posterior_predictive_w(traces, 1000, models,
                        weights=comp.weight.sort_index(ascending=True),
                        progressbar=False)
