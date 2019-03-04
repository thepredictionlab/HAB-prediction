# make_posterior.py

## gives a posterior distribution (subdividing [-2,2] for tbv percent change) for the given model folder
## 

import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import datetime as dt
from py_help import *

def makePost(folder):
	paramVector,dataNames,WVARS,priorMeans,priorStds,paramType,paramList,muList,stdList,dataNameList = defns('l')
	dirlist = np.array(os.listdir(folder))
	# list all the Trace_v_....npz files
	traceList = dirlist[np.where([k[:6]=='Trace_' for k in dirlist])[0]]
	print(traceList)
	Post = np.zeros((len(traceList),len(pBins)-1))
	Like = np.zeros(len(traceList))
	t = 0
	dtp = np.load(folder+'data_time_priors.npz')
	# get variables
	varList = folder.split('/')[-2].split('_')
	varInds = np.isin(dataNameList,varList)
	paramList = np.array(paramList)[varInds]
	for elmt in traceList:
		# use posterior to estimate response for left out time
		traceData = np.load(folder + elmt)
		looNum = int(elmt.split('_')[1])
		trace = traceData['TRACE']
		data = traceData['DATA']
		data = np.reshape(data,(data.shape[0],-1))
		data = data[dtp['domain'],:]
		response = predOut[dtp['domain']]

		predictors = data[looNum,:]
		observed = response[looNum]
		mSamples = len(trace)
		testResults = np.zeros(mSamples)
		for j in np.arange(mSamples):
			output = 0
			pValues = trace[j]['c']
			pInd = 0
			for k in np.arange(len(paramList)):
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
		Post[t,:] = d
		Like[t] = np.dot(d,o)
		plt.hist(testResults, bins=30)
		plt.title('LOO: '+str(looNum)+', valid\'n: '+str(np.round(observed,3))+', lk: '+str(Like[t]))
		plt.savefig(folder+'Validation_'+str(looNum)+'_'+timestamp())
		#plt.show()
		plt.clf()
		t += 1
		print(t)

	np.savez(folder+'Posterior', posterior=Post, likelihood=Like)





