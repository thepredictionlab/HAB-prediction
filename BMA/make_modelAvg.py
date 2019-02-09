# make_modelAvg.py
# Mathew Titus, The Prediction Lab, 2/7/19
# 
## iterate through files in /Data/Posteriors/[experiment]/
## named data_posterior[k].npz with k an index ranging from 0 to some maximum value;
## calculate the likelihood of each given a set of observations 
## ('data_validation.npz', found in /Data/Posteriors/[experiment] folder)
## Each data_posterior 
## file contains an array giving the posterior distribution output by 
## a given model, each row corresponding to a different input value.
## We record the averaged posterior distribution across models as well as
## the vector of likelihoods.

## NB: adding other 'data_*' files to the experiment's directory 
##     will cause errors.

import numpy as np
import os

filePath = '../Data/Posteriors/test/'
fileList = os.listdir(filePath)
fileNumbers = []
for file in fileList:
	if file[:5] == 'data_'
		file = file[14:-4]
		fileNumbers.append(file)

# get spatial resolution (# colms) and number of observations (# rows)
file1 = np.load(filePath + 'data_posterior' + fileNumbers[0] + '.npz')
[N,X] = file1['data'].shape
# get 'answers' to the predictions made for the N validations
observed = np.load(filePath + 'data_validation.npz')

for k in fileNumbers:
	post = np.load(filePath + 'data_posterior' + k + '.npz')
	post = post['data']
	if (post.shape[0] ~= N)|(post.shape[1] ~= X):
		raise Exception('Dataset ' + k + ' does not match the format of file1, a ('+str(N)+', '+str(X)+') matrix.')
		continue
	else:
# Process data
		