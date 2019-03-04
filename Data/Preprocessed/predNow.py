# predNow.py
## collect likelihoods, compile model posteriors for latest data.

## TODO: deal with the presence/absence of dates from individual models....

import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import datetime as dt
import pymc3 as pm
from py_help import *

paramVector,dataNames,WVARS,priorMeans,priorStds,paramType,paramList,\
muList,stdList,dataNameList = defns('l')

limitedDataNameList = dataNameList

paramVector,dataNames,WVARS,priorMeans,priorStds,paramType,paramList,\
muList,stdList,dataNameList = defns('h')

rootFolder = fatherDir+'BMA/Results/'

# the folder where the individual predictions of each model will be saved
# to later be averaged together using allLikelihoods.npz
saveFolder = '2019Pred/'

dirlist = np.array(os.listdir(rootFolder+experiment))

# for each model, for each time step in/after 2018, make trace based on past data
# use trace and 'present' data to get posterior on 'present observation'
# normalize each posterior, then take weighted average

# filter non-experiment-folder items from dirlist
dirlist = dirlist[[np.all(np.isin(j.split('_'),dataNameList)) for j in dirlist]]

trainingDataSize = {}

for bmodel in dirlist:
	# load data, parameter types, and priors
	varlist = bmodel.split('_')
	if np.all(np.isin(varlist,limitedDataNameList)):
		varInds = np.isin(dataNameList,varlist)
		params = np.array(paramList)[varInds]
		parInds = varInds.copy().tolist()
		j = 0
		for k in np.arange(len(varInds)):
			if paramList[k] != 'l':
				parInds.insert(j+1,parInds[j])
				j += 2
			else:
				j += 1
		mus = np.array(muList)[parInds]
		stds = np.array(stdList)[parInds]
		numParams = len(mus); numVars = len(varlist);

		data = np.array(df[varlist])
		domain = ~np.any(np.isnan(data),1) & ~np.isnan(predOut)
		
		numPre2018 = np.where(dates[domain] < 2018)[0].shape[0]
		trainingDataSize[bmodel] = numPre2018
		
		data2 = np.array(df2[varlist])
		domain2 = ~np.any(np.isnan(data2),1) & ~np.isnan(presPredOut)
		numPresent = np.where(dates2[domain2] >= 2019)[0].shape[0]
		numPre2019 = np.where(dates2[domain2] < 2019)[0].shape[0]
		if numPresent > 0:
			predictions = np.zeros((len(dates2[domain2])-numPre2019, len(pBins)))
			for t in np.arange(numPre2019,len(dates2)):
				if domain2[t]:
					# train the model (using pre-2019 data, fixed across weeks)
					subDomain = domain.copy() & (dates < dates2[t])
					subData =data[subDomain,:]
					obsOut = predOut[subDomain]
					linModel = pm.Model()
					with linModel as model:
						norm100 = pm.Normal('n', mu=0, sd=1)
						coeffs = pm.Normal('c', mu=np.array([float(k) for k in mus]),\
						 sd=np.array([float(k) for k in stds]), shape=numParams)
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

					# compute the posterior distribution of the predictand
					predictors = data2[t,:]
					observed = presPredOut[t]
					testResults = np.zeros(mSamples)
					for j in np.arange(mSamples):
						output = 0
						pValues = trace[j]['c']
						pInd = 0
						for k in np.arange(len(params)):
							if paramList[k] == 's':
								output += sigmoid(predictors[k],pValues[pInd+1],pValues[pInd])
								pInd += 2
							elif paramList[k] == 'q':
								output += quadratic(predictors[k],pValues[pInd],pValues[pInd+1])
								pInd += 2
							elif paramList[k] == 'l':
								output += predictors[k]*pValues[pInd]
								pInd += 1
							else:
								raise Exception('This error should never have happened.')
						testResults[j] = output+trace[j]['noise']
					d = getDist(testResults,pBins)
					o = getDist(observed,pBins)
					likeli = np.dot(d,o)
					np.savez(rootFolder+experiment+saveFolder+bmodel+str(t)+'.npz',\
						 LIKELIHOOD=likeli, PREDICTION=testResults, POSTERIOR=d, OBS=observed, TRACE=trace, BINS=pBins)

# trainingDataSize records how many training data points each model has before the 2018 season.
# Lower values should lower the predictive value of the model, and may be used to influence the averaging process.
np.savez(rootFolder+experiment+saveFolder+'trainingDataSize.npz', TDS=trainingDataSize)




















