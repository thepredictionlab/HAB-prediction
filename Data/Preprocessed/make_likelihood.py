# make_likelihood.py

# exec(open('make_likelihood.py').read())

## collects posterior results from models and calculates their mean likelihood 
## (categorical, + / -) from validation data.  This results in a dictionary of the form fileName : likelihood

## data_time_priors.npz holds an element 'response' which has all the response var values
	## and a 'data' element which holds the input variables over time.
## Posterior.npz holds the posterior performance of the model on the validation set (LOO)

import os
import numpy as np
from make_posterior import *
from py_help import *

likelihood = dict()
rootFolder = fatherDir+'BMA/Results/'
dirlist = os.listdir(rootFolder+experiment)
dirlist = np.array(dirlist)
dirlist = dirlist[np.where([k[0]!='.' for k in dirlist])[0]]
dirlist = dirlist[np.where([k[:3]!='all' for k in dirlist])[0]]

for modelFolder in dirlist:
	if ~np.isin(modelFolder,list(likelihood)):
		subdir = rootFolder+experiment+modelFolder+'/'
		makePost(subdir)
		subdirlist = os.listdir(subdir)
		subdirlist = np.array(subdirlist)
		subdirlist = subdirlist[np.where([k[:9]=='Posterior' for k in subdirlist])[0]]
		if len(subdirlist) > 0:
			results = np.load(subdir+subdirlist[0])
			modelLikelihood = np.nanmean(results['likelihood'])
			likelihood[modelFolder] = modelLikelihood
		else:
			print('No posterior file exists for model '+modelFolder)

np.savez(rootFolder+experiment+'allLikelihoods.npz', LIKELIHOOD=likelihood)




