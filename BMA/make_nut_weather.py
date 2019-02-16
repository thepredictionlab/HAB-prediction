# make_something.py

## for plotting: in bash execute
## export DISPLAY=localhost:0.0
## and turn on XMing

import matplotlib.pyplot as plt
#plt.style.use(['seaborn-darkgrid'])
import datetime as dt
import pymc3 as pm
import numpy as np
def timestamp():
	now = dt.datetime.now()
	s = str(now.year) + str(now.month) + str(now.day) + "_" \
	    + str(now.hour) + str(now.minute) + str(now.second)
	return s


folderName = './Results/TBV/4nut_tem_maxt_rain/'
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
inputs = {'abundances':{'TBV'}, #},#'TBV'}, # need loc and div
	'divisions':[4], #[0,1,2],
	'locations':['BB','BO','HA','HT','LB','LBP','LBS'], #['LB','BO'],
	'nutrients':{'NUT1','NUT2','NUT3','NUT4'},		#{}, #{'NUT2'},		# need loc
	'toxins':{}, # 'TOX4'},
	'weather':{'TEMP','MAXT','RAIN'}} # 'TEMP','MINT','MAXT'}}

dataList = ['cyano_tbv',\
			'nut1','nut2','nut3','nut4',\
			'tem','max_tem','rain']
paramType = ['l',\
			 's','s','s','s',\
			 'q','l','l']
# Note dataOut will already be included in dataFull
dataOut = data['TBV'][:,4]
paramOut = 'l'
#dataOut = np.nanmax(data['TOX4'][:,:],1)
predOut = dataOut[1:]
predOut = np.insert(predOut,len(predOut),np.nan)

# define some priors - coefficients relating input to dataOut var.
mus = [0.05, 0.1, 0.02, 0.04, 0.02, 0.04, 0.05, 0.1,\
        8.5, -0.3,\
        0.5, 0.5, 0.5]
# note the final linear coeff for dataout (current tbv level)

# simple variance estimate
stds = 0.2*np.ones(len(mus))

dataFull = np.array([])
for dc in inputs['abundances']:
	for div in inputs['divisions']:
		divFull = np.zeros((data['TIME'].shape[0],1)) * np.nan
		#for loc in inputs['locations']:
		#	ind1 = data['LOCS'].tolist().index(loc)
		divFull = np.append(divFull, 
			np.reshape(data[dc][:,div],(divFull.shape[0],1)), 1)
			#np.reshape(data[dc][:,ind1,div],(divFull.shape[0],1)), 1)
		dataFull = np.append(dataFull, np.nanmean(divFull,1))
# Calculate the mean nutrient level across relevant sampling sites
nutFull = np.zeros((data['TIME'].shape[0],1)) * np.nan
for nut in inputs['nutrients']:
	for loc in inputs['locations']:
		nutFull = np.append(nutFull, data[nut][:,np.where(data['LOCS'] == loc)[0]], axis = 1)
	dataFull = np.append(dataFull, np.nanmean(nutFull,1))
# Collect MAX toxin level across sample sites
for item in inputs['toxins']:
	dataFull = np.append(dataFull, np.nanmax(data[item],1))
# Collect weather data from detroit lake station
for item in inputs['weather']:
	dataFull = np.append(dataFull, data[item])
dataFull = dataFull.reshape((data['TEMP'].shape[0],-1), order = 'F')
print('dataFull shape = ' + str(dataFull.shape) + '.')

## Validation via LOO:
# dates/times with NaNs
nanDays = np.any(np.isnan(dataFull),1) | np.isnan(predOut)
domain = ~nanDays
print('Domain consists of ' + str(dataFull[domain,0].shape[0]) + ' days.')
sampledTimes = TIME[domain]
print(str(np.where(sampledTimes < 2018)[0].shape[0]) + ' days before 2018.')

# Data matrix made. Make output matrix.
y1 = np.array( predOut )

# for nut1 and 4, C = 0.05, L = 0.1
# for nut2 and 3, C = 0.02, L = 0.04
# for nut3, 
def sigmoid(x,L,C):
	S = L*np.exp((x-C)/(2*L))/(np.exp((x-C)/(2*L))+1)
	return S

# for temp max growth at 29 C = 84.5 F
# so use v = 8.5, c < 0, h s.t. Q(40 F) = 0
# also using nonnegative range cause cold is just cold
def quadratic(x,v,c):
	h = -(c*(4-v)**2)
	Q = c*np.power((x-v),2) + h
#	Q = [np.max(q,0) for q in Q]
	return Q

mSamples = 1000
mCores = 4
mTuning = 1000
numVars = len(mus) # dataFull.shape[1]
for t in np.arange(0,len(sampledTimes)):
	v = t #np.random.randint(len(sampledTimes))
	subDomain = domain
	subDomain[np.where(domain)[0][v]] = False
	subData =dataFull[subDomain,:]
	obsOut = y1[subDomain]
	linModel = pm.Model()
	with linModel as model:
		norm100 = pm.Normal('n', mu=0, sd=1)
		coeffs = pm.Normal('c', mu=mus, sd=stds, shape=numVars)
		mu = 0
		pInd = 0
		for k in np.arange(len(paramType)):
			if paramType[k] == 's':
				mu += sigmoid(subData[:,k],coeffs[pInd+1],coeffs[pInd]) # (x,L,C)
				pInd += 2
			elif paramType[k] == 'q':
				mu += quadratic(subData[:,k],coeffs[pInd],coeffs[pInd+1])
			elif paramType[k] == 'l':
				mu += subData[:,k]*coeffs[k]
				pInd += 1
			else:
				raise Exception('Not a valid variable type [paramType].')
		eta = pm.Normal('noise', mu=norm100, sd=1)
		#eta2 = pm.AR('ar_noise',eta,sd = 1.0)
		obs = pm.Normal('obs', mu = mu+eta, sd = 1, observed = obsOut)
		trace = pm.sample(mSamples, cores = mCores, tuning = mTuning)
	print(np.std(trace['c'],0))
	np.savez(folderName+'Trace_'+str(v)+'_'+timestamp(),\
			TRACE=trace,DATA=dataFull,RESPONSE=obsOut)
	plt.figure()
	pm.traceplot(trace)
	plt.legend(dataList)
	plt.title(paramType)
	plt.savefig(folderName+'Traceplot_'+str(v)+'_'+timestamp())
	#plt.show()
	plt.close()
	# use posterior to estimate response for left out time
	testData = np.reshape(dataFull[np.where(domain)[0][v],:], -1, 1)
	print(testData)
	validation = y1[np.where(domain)[0][v]]
	testResults = np.zeros(mSamples)
	for j in np.arange(mSamples):
		output = 0
		pValues = trace['c'][j]
		pInd = 0
		for k in np.arange(len(paramType)):
			if paramType[k] == 's':
				output += sigmoid(testData[k],pValues[pInd+1],pValues[pInd])
				pInd += 2
			elif paramType[k] == 'q':
				output += quadratic(testData[k],pValues[pInd],pValues[pInd+1])
				pInd += 2
			elif paramType[k] == 'l':
				output += testData[k]*pValues[pInd]
				pInd += 1
			else:
				raise Exception('This error should never have happened.')
		testResults[j] = output
	plt.figure()
	plt.hist(testResults, bins=50)
	plt.title('LOO: '+str(v)+', valid\'n: '+str(np.round(validation,3)))
	plt.savefig(folderName+'Validation_'+str(v)+'_'+timestamp())
	#plt.show()
	plt.close()

#ppc_w = pm.sample_posterior_predictive_w(traces, 1000, models,
 #                       weights=comp.weight.sort_index(ascending=True),
 #                       progressbar=False)

# # response variable
# y1 = np.array( data['TOX3'][domain,4] )

# # initiate model
# # assume linear response to the variables, normal priors for coefficients

# mSamples = 100
# mCores = 2
# mTuning = 50
# linModel = pm.Model()
# with pm.Model() as model:
# 	norm100 = pm.Normal('n', mu=0, sd=1) # too broad?
# 	coeffs = pm.Normal('c', mu=0, sd=1, shape=dataFull.shape[1])
# 	mu = 0
# 	for k in np.arange(dataFull.shape[1]):
# 		mu += dataFull[domain,k]*coeffs[k]
# 		#print('coeff'+str(k)+' = '+dataList[k])
# 	eta = pm.Normal('noise', mu=norm100, sd=10)
# 	#eta2 = pm.AR('ar_noise',eta,sd = 1.0)
# 	obs = pm.Normal('obs', mu = mu+eta, sd = 1, observed = y1)#[domain])
# 	trace = pm.sample(mSamples, cores = 4, tuning = mTuning)

# saveTraceName = './trace.npz'
# np.savez(saveTraceName,TRACE=trace)

# pm.summary(trace)
# pm.traceplot(trace)
# plt.show()

# savePlotName = folderName+str(mSamples)+'cores'+str(mCores)+'tuning'+str(mTuning)
# plt.savefig(savePlotName)

# plt.show()



#np.where(np.any(~np.isnan(data['TOX4']),1))[0].shape
#np.where(np.any(~np.isnan(dataFull[:,0:3]),1))[0].shape