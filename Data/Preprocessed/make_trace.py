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
import pandas as pd
from py_help import *
paramVector,dataNames,WVARS,priorMeans,priorStds,paramType,paramList,muList,stdList,dataNameList = defns('h')

# inclusionVector = [\
# 	[0,0,0,0,1,0,0,0,0,0,0,0,0], #TBV
# 	[0,0,0,0], # NUT
# 	[0,0,0,0], # TOX
# 	[0,0,0,1,0,0,0,1,0,0,1]] #WEA: TEMP, MAXT, MINT, WDEG, HUM, PWI, PRES, RAIN, DVR, WIS, RAD
###print(inclusionVector)


# # mask unwanted entries: 
# inputs = {'abundances':{'TBV'}, #},#'TBV'}, # need div
# 	'divisions':np.where(inclusionVector[0])[0], #[0,1,2,4], #[0,1,2],
# 	'max_divs':np.where(inclusionVector[1])[0],
# 	#'locations':LOCS,
# 	'nutrients':np.where(inclusionVector[2])[0], #{},#'NUT1','NUT2','NUT3','NUT4'},		#{}, #{'NUT2'},		# need loc
# 	'toxins':np.where(inclusionVector[3])[0], #{}, # 'TOX4'},
# 	'weather':np.where(inclusionVector[4])[0]} #{'TEMP','MAXT','RAIN'}}#,'WDEG'}} # 'TEMP','MINT','MAXT'}}

def makeTrace(inclusionVector):
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

	dataFull = np.array(df[dataNameList])

	# Define target folder, create if it doesn't already exist
	folderName = fatherDir+'BMA/Results/'+experiment #+str(inclusionVector)+'/' #4tbv_tem_maxt_rain/'
	for k in dataNameList:
		folderName += k
		folderName += '_'

	folderName = folderName.rstrip('_') + '/'
	subprocess.call(["./folderMaker.sh",folderName])

	## Validation via LOO:
	# dates/times with NaNs
	nanDays = np.any(np.isnan(dataFull),1) | np.isnan(predOut)
	domain = ~nanDays
	print('Domain consists of ' + str(dataFull[domain,0].shape[0]) + ' observations.')
	sampledTimes = dates[domain]
	numPre2018 = np.where(sampledTimes < 2018)[0].shape[0]
	print(str(numPre2018) + ' days before 2018.')

	np.savez(folderName+'data_time_priors.npz',data=dataFull,times=dates,domain=domain,response=predOut,\
		priorMeans=priorMeans,priorStds=priorStds)

	plt.figure()
	for t in np.arange(0,numPre2018):
		v = t #np.random.randint(len(sampledTimes))
		subDomain = domain.copy()
		subDomain[np.where(domain)[0][v]] = False
		subData =dataFull[subDomain,:]
		obsOut = predOut[subDomain]
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
		pm.traceplot(trace)
		plt.legend(dataNameList)
		plt.title(paramType)
		plt.savefig(folderName+'Traceplot_'+str(v)+'_'+timestamp())
		plt.clf()

####################### JUNK CODE ###############################################

#ppc_w = pm.sample_posterior_predictive_w(traces, 1000, models,
 #                       weights=comp.weight.sort_index(ascending=True),
 #                       progressbar=False)

# # response variable
# y1 = np.array( data['TOX3'][domain,4] )

# pm.summary(trace)
# pm.traceplot(trace)
# plt.show()

# savePlotName = folderName+str(mSamples)+'cores'+str(mCores)+'tuning'+str(mTuning)
# plt.savefig(savePlotName)

################################################################################

#np.where(np.any(~np.isnan(data['TOX4']),1))[0].shape
#np.where(np.any(~np.isnan(dataFull[:,0:3]),1))[0].shape