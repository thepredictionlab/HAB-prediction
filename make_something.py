# make_something.py

# File1 = '../Data/Raw_lake/make_historical.py'
# File2 = './make_something.py'
# exec(open(File2).read())
## for plotting: in bash execute
## export DISPLAY=localhost:0.0
## and turn on XMing

import matplotlib.pyplot as plt
#plt.style.use(['seaborn-darkgrid'])
import subprocess
import datetime as dt
import pymc3 as pm
import numpy as np
def timestamp():
	now = dt.datetime.now()
	s = str(now.year) + str(now.month) + str(now.day) + "_" \
	    + str(now.hour) + str(now.minute) + str(now.second)
	return s

data_path = '../Data/Preprocessed/'
data_file = 'Data_hist_interp.npz'
data = np.load(data_path + data_file)

inclusionVector = [\
	[0,0,0,0,1,0,0,0,0,0,0,0,0], #TBV
	[0,0,0,0], # NUT
	[0,0,0,0], # TOX
	[0,0,0,1,0,0,0,1,0,0]] #WEA: TEMP, MAXT, MINT, WDEG, HUM, PWI, PRES, RAIN, DVR, WIS
print(inclusionVector)
# Note dataOut will already be included in dataFull
dataOut = data['TBV'][:,4]
# Future values of predictand
## next value
#predOut = dataOut[1:]; predOut = np.insert(predOut,len(predOut),np.nan)
## tbv "return"
predOut = np.diff(dataOut)/np.array(dataOut[:-1])
predOut = np.insert(predOut,len(predOut),np.nan)
predOut[np.where(predOut < 0)[0]] = -1
predOut[np.where(predOut > 0)[0]] = 1
#predOut[np.where(predOut==0)[0]] = np.nan

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
# Param types for each variable
paramVector = [\
	['l','l','l','l','l','l','l','l','l','l','l','l','l'], # BIO
	['s','s','s','s'], # NUT
	['l','l','l','l'], # TOX
	['q','q','q','l','l','l','l','l','l','l']] # WEA
# TEMP, MAXT, MINT, WDEG, HUM, PWI, PRES, RAIN, DVR, WIS

# Define priors for all variables
priorMeans = [\
	[-1,-1,-1,-1,+1,-1,-1,-1,-1,-1,-1,-1,-1], # BIO
	[[0.008, 0.055],[0.002, 0.025],[0.002, 0.025],[0.008, 0.055]], # NUT
	[1,1,1,1], # TOX
	[[8.3,-0.1],[8.5,-0.3],[8.5,-0.3],1,0,-1,0,1,0,-1] # WEA
	]
priorStds = [\
	[1,1,1,1,1,1,1,1,1,1,1,1,1], # BIO
	#[10,10,10,10,10,10,10,10,10,10,10,10,10], # BIO
	[[0.1,0.1],[0.1,0.1],[0.1,0.1],[0.1,0.1]], # NUT
	[2,2,2,2], # TOX
	[[1,0.5],[2,0.5],[2,0.5],1,1,1,1,1,1,1] # WEA
	]


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
paramType = []; mus = []; stds=[]; paramList = [];
muList = []; stdList = []; dataNameList = [];
numVars = 0; numParams = 0;
for i1 in np.arange(len(inclusionVector)):
	paramType.append([]); mus.append([]); stds.append([]);
	for i2 in np.where(inclusionVector[i1])[0]:
		paramType[i1].append(paramVector[i1][i2])
		mus[i1].append(priorMeans[i1][i2])
		stds[i1].append(priorStds[i1][i2])
		numParams += len(np.array(priorMeans[i1][i2]).flatten())
		paramList.extend(paramVector[i1][i2])
		muList.extend([priorMeans[i1][i2]])
		stdList.extend([priorStds[i1][i2]])
		dataNameList.extend([dataNames[i1][i2]])
	numVars += len(np.array(paramType[i1]))

paramList = [k for sublist in paramList for item in [sublist if isinstance(sublist,list) else [sublist]] for k in item]
muList = [k for sublist in muList for item in [sublist if isinstance(sublist,list) else [sublist]] for k in item]
stdList = [k for sublist in stdList for item in [sublist if isinstance(sublist,list) else [sublist]] for k in item]
dataNameList = [k for sublist in dataNameList for item in [sublist if isinstance(sublist,list) else [sublist]] for k in item]

# Define target folder, create if it doesn't already exist
folderName = './Results/TBVcat/' #+str(inclusionVector)+'/' #4tbv_tem_maxt_rain/'
for k in dataNameList:
	folderName += k
	folderName += '_'

folderName = folderName.rstrip('_') + '/'
subprocess.call(["./folderMaker.sh",folderName])

dataFull = np.array([])
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
		nutFull = np.append(nutFull, data['NUT'+str(nut)][:,np.where(data['LOCS'] == loc)[0]], axis = 1)
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

## Validation via LOO:
# dates/times with NaNs
nanDays = np.any(np.isnan(dataFull),1) | np.isnan(predOut)
domain = ~nanDays
print('Domain consists of ' + str(dataFull[domain,0].shape[0]) + ' observations.')
sampledTimes = data['TIME'][domain]
numPre2018 = np.where(sampledTimes < 2018)[0].shape[0]
print(str(numPre2018) + ' days before 2018.')

# Data matrix made. Make output matrix.
y1 = np.array( predOut )
np.savez(folderName+'data_time_priors.npz',data=dataFull,times=data['TIME'],domain=domain,response=y1,\
	priorMeans=priorMeans,priorStds=priorStds)

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

mSamples = 5000
mCores = 2
mTuning = 2000
# numVars = len(mus) # dataFull.shape[1]
for t in np.arange(13,numPre2018): #len(sampledTimes)-15,len(sampledTimes)):
	v = t #np.random.randint(len(sampledTimes))
	subDomain = domain.copy()
	subDomain[np.where(domain)[0][v]] = False
	subData =dataFull[subDomain,:]
	obsOut = y1[subDomain]
	linModel = pm.Model()
	with linModel as model:
		norm100 = pm.Normal('n', mu=0, sd=1)
		coeffs = pm.Normal('c', mu=np.array([float(k) for k in muList]),\
		 sd=np.array([float(k) for k in stdList]), shape=numParams)
		mu = 0
		pInd = 0
		for k in np.arange(numVars):
			if paramList[k] == 's':
				mu += sigmoid(subData[:,k],coeffs[pInd+1],coeffs[pInd]) # (x,L,C)
				pInd += 2
			elif paramList[k] == 'q':
				mu += quadratic(subData[:,k],coeffs[pInd],coeffs[pInd+1])
				pInd += 2
			elif paramList[k] == 'l':
				mu += subData[:,k]*coeffs[k]
				pInd += 1
			else:
				raise Exception('Not a valid variable type [paramType].')
		eta = pm.Normal('noise', mu=norm100, sd=1)
		#eta2 = pm.AR('ar_noise',eta,sd = 1.0)
		obs = pm.Normal('obs', mu = mu+eta, sd = 1, observed = obsOut)
		trace = pm.sample(mSamples, init='jitter+adapt_diag', cores = mCores, tune = mTuning, step_scale=0.2)
	print(np.std(trace['c'],0))
	np.savez(folderName+'Trace_'+str(v)+'_'+timestamp(),\
			TRACE=trace,DATA=dataFull,RESPONSE=obsOut)
	plt.figure()
	pm.traceplot(trace)
	plt.legend(dataNameList)
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
		for k in np.arange(len(paramList)):
			if paramList[k] == 's':
				output += sigmoid(testData[k],pValues[pInd+1],pValues[pInd])
				pInd += 2
			elif paramList[k] == 'q':
				output += quadratic(testData[k],pValues[pInd],pValues[pInd+1])
				pInd += 2
			elif paramList[k] == 'l':
				output += testData[k]*pValues[pInd]
				pInd += 1
			else:
				raise Exception('This error should never have happened.')
		testResults[j] = output+trace['noise'][j]
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