import numpy as np
import matplotlib.pylab as plt
from matplotlib.pyplot import figure, show, rc
import matplotlib.cm as cm



## Data
data = np.load("../Data/Preprocessed/Data_historical_7day.npz")
div  = data['DIV']  # column 4 is cyanobacteria
locs = data['LOCS'] # column 4 is the log boom

fbv  = data['FBV'][:,4,4]
temp = data['TEMP']
pwi  = data['PWI']
nut1 = data['NUT1'][:,4]
nut2 = data['NUT2'][:,4]
nut3 = data['NUT3'][:,4]
nut4 = data['NUT4'][:,4]

## Normalize data
fbv  = fbv  / np.nanmax(fbv)
temp = temp / np.nanmax(temp)
pwi  = pwi  / np.nanmax(pwi)
nut1 = nut1 / np.nanmax(nut1)
nut2 = nut2 / np.nanmax(nut2)
nut3 = nut3 / np.nanmax(nut3)
nut4 = nut4 / np.nanmax(nut4)
data = np.vstack((fbv,temp,pwi,nut1,nut2,nut3,nut4)).transpose()

## Choose a time to plot
ds = np.sum(data,1)
JD = np.where(np.isnan(ds)==0)[0] # when all data is present

ID = JD[2] # choose a day
dp = data[ID,:]
N = len(dp) # number of petals
radii = 10*dp # length of petals

## Plot
# force square figure and square axes looks better for polar, IMO
fig = figure(figsize=(8,8))
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], polar=True)

theta = np.arange(0.0, 2*np.pi, 2*np.pi/N) # angle of petals
width = np.pi/4*np.ones(N) # width of petals
bars = ax.bar(theta, radii, width=width, bottom=0.0)
for r,bar in zip(radii, bars):
    bar.set_facecolor( cm.jet(r/N))
    bar.set_alpha(0.5)

ax.set_xticklabels(['Cyano', 'Temp', 'P. Wind', 'NO3/2', 'O-Phos', 'TN', 'I-Phos', ''])


show()




