# make_something.py

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
	'divisions':{1,2,3},
	'locations':{'LB','BO'},
	'nutrients':{'NUT1'},		# need loc
	'weather':{'TEMP'}}	

inputOutput = inputs.copy()
inputOutput['toxins'] = ['TOX3']

dataFull = np.array([])
dataList = []
for dc in inputs['abundances']:
	for loc in inputs['locations']:
		for div in inputs['divisions']:
			dataFull = np.append(dataFull, data[dc][], 0)

for 

for dc in inputs:
	for item in inputs[dc]:
		dataLen = len(data2[dc][item])
		dataAdd = data2[dc][item].reshape((dataLen,1))
		dataFull = np.append(dataFull.reshape((dataLen,-1)),dataAdd,axis=1)
		dataList.append(item)

mask = np.ndarray((dataFull.shape[0],1), dtype=bool)
mask = mask|~mask
for dc in inputOutput:
	for item in inputOutput[dc]:
		itemMask = ~np.isnan(data2[dc][item])
		print(item + ' has ' + str(itemMask.sum()) + ' entries.')
		mask = mask & itemMask.reshape((dataLen,1))
		print(mask.shape)
y1 = np.array( data2['toxins']['cylL'] ).reshape((dataLen,1))[mask]
dataFull = dataFull[mask.flatten(),:]

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
		mu += dataFull[:,k]*coeffs[k]
		print('coeff'+str(k)+' = '+dataList[k])
	eta = pm.Normal('noise', mu=norm100, sd=10)
	obs = pm.Normal('obs', mu = mu+eta, sd = 1, observed = y1)
	trace = pm.sample(5000, cores = 4, tuning = 500)
save_obj(trace, saveTraceName)

pm.traceplot(trace)
savePlotName = folderName+str(mSamples)+'cores'+str(mCores)+'tuning'+str(mTuning)
plt.savefig(savePlotName)

plt.show()

