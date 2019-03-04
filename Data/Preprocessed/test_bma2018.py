# test2018b.py - combines the posteriors saved in test2018, 
## weighted according to allLikelihoods.npz to create 
## the true posterior estimate for the predictor at each time step.

### exec(open('test_bma2018.py').read())
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import datetime as dt
import pymc3 as pm
from py_help import *

maxDateIndex = 500

rootFolder = fatherDir+'BMA/Results/'
# the folder holding the models predicting a given response variable

##################################
experiment = 'TBVcat/'
##################################

# the folder where the individual predictions of each model will be saved
# to later by averaged together using allLikelihoods.npz
saveFolder = '2018Pred/'

# load likelihoods
lh = np.load(rootFolder+experiment+'allLikelihoods.npz')
likelihood = lh['LIKELIHOOD']
dirlist = np.array(os.listdir(rootFolder+experiment+saveFolder))

# filter for the dated .npz files and their associated dates
dirlist = dirlist[[k[-5].isnumeric() for k in dirlist]]
datelist = [''.join([j[n] for n in np.arange(-7,-4) if j[n].isnumeric()]) for j in dirlist]

preds = preds = np.load(rootFolder+experiment+saveFolder+dirlist[0])
finalPost = {}
plt.figure(figsize=(6.4,5.8))
for date in np.arange(maxDateIndex):
	t = str(date)
	if (np.isin(t,datelist))&(dates[date]>2018):
		sublist = dirlist[[k[-(4+len(t)):-4] == t for k in dirlist]]
		sublist = dirlist[[k == t for k in datelist]]
		modelList = [k.split(t)[0] for k in sublist]
		like = [0]
		post = np.zeros((1,len(pBins)-1))
		for bmodel in modelList:
			varlist = bmodel.split('_')
			preds = np.load(rootFolder+experiment+saveFolder+bmodel+str(t)+'.npz')
			# LIKELIHOOD=likeli, PREDICTION=testResults, OBS=observed, TRACE=trace, BINS=pBins
			like.append(preds['LIKELIHOOD'])
			modelPred = np.reshape(preds['PREDICTION'],(1,-1))
			modelPost = preds['POSTERIOR']
			post = np.vstack((post,modelPost))
		like = np.array(like)
		like = like/np.sum(like)
		like = np.reshape(like,(1,-1))
		finalPost[date] = np.matmul(like,post)
		binCenters = (preds['BINS'][1:] + preds['BINS'][:-1])/2
		expectation = np.dot(binCenters[1:-1], finalPost[date][0][1:-1])
		date
		dates[date]
		# plt.title('Observed: '+str(predOut[date])+'\nMean pred: '+str(expectation)+'\nDate: '+str(year2datetime(dates[date])))
		# plt.plot(binCenters, finalPost[date][0])
		# # plt.scatter([-1,1],finalPost[date][0])
		# #plt.show()
		# plt.savefig(rootFolder+experiment+saveFolder+'/Figures/prediction'+str(date)+'.png')
		# plt.clf()

np.savez(rootFolder+experiment+saveFolder+'2018results.npz',\
	POSTERIOR=finalPost,OBSERVATION=predOut,DATE=dates,BINCENTERS=binCenters)



#s = np.sort(list(likelihood.tolist().values()))
#np.where(np.array(list(likelihood.tolist().values())) > 0.2)

#list(likelihood.tolist().keys())[50]

# Max occurs at model OP
# Next occurs at NO

lik = 0
for k in np.arange(279,290):
	if predOut[k] < 0:
		lik += np.sum(finalPost[k][0][:11])
	else:
		lik += 1 - np.sum(finalPost[k][0][:11])

lik = lik/10
