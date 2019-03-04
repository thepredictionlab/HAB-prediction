# fnc_mcmc.py

## exec(open('fnc_mcmc.py').read())

import numpy as np
from make_trace import *
from py_help import *
paramVector,dataNames,WVARS,priorMeans,priorStds,paramType,paramList,muList,stdList,dataNameList = defns('h')

##############################################################################################################

# Define the experiment to be run:
### PCTBV = percent change Total biovolume, resolved posterior
### PCTBVS = sign of the percent change of cyanobacteria TBV
### TBVZ = total biovolume is either high ( > 1.25sigma ) or low, 
		### corresponding to +1 or -1 for the posterior (classify positive output as high, negative as low)
### TOXC = toxin level, resolved into categories: low (-1), medium (0), high (+1)
### TOXLVL = microcystin ELISA toxin level, sorted into categories (0,1,2,3) as defined by the state

EXP = 'PCTBV' # go to py_help and edit EXP

##############################################################################################################

inclusionVector = [\
	[0,0,0,0,0,0,0,0,0,0,0,0,0], #TBV
	[0,0,0,0,0,0,0,0,0,0,0,0,0], #MTBV
	[0,0,0,0], # NUT
	[0,0,0,0], # TOX
	[0,0,0,0,0,0,0,0,0,0,0]] #WEA: TEMP, MAXT, MINT, WDEG, HUM, PWI, PRES, RAIN, DVR, WIS, RAD

# component lengths
j0 = 0; j1 = 13; j2 = 26; j3 = 13 + 17; j4 = 13+21; j5 = 13+32;
vecLen = 0
for k in np.arange(len(inclusionVector)):
	vecLen += len(inclusionVector[k])


while True:
	inclusionList = np.zeros(vecLen)
	for itrn in np.arange(6):
		## Modify one entry of the inclusion list
		inclusionList[np.random.randint(j4,j5)] += 1 # np.random.randint(vecLen)] += 1
		inclusionList = np.mod(inclusionList,2)
		inclusionVector[0] = inclusionList[j0:j1]
		inclusionVector[1] = inclusionList[j1:j2]
		inclusionVector[2] = inclusionList[j2:j3]
		inclusionVector[3] = inclusionList[j3:j4]
		inclusionVector[4] = inclusionList[j4:j5]
		## ?????????? Fix probability to keep # of vars low
		print('work to be done.')
		#print(inclusionVector)
		print(inclusionList)
		## Check to make sure a variable is included
		if np.any(inclusionList):
			makeTrace(inclusionVector)
			# file = './make_trace.py'
			# exec(open(file).read())
			print('Posteriors generated for model.')

