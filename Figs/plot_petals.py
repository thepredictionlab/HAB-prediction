import numpy as np
import matplotlib.pylab as plt
from matplotlib.pyplot import figure, show, rc
import matplotlib.cm as cm


## Data
data = np.load("../Data/Preprocessed/Data_historical_7day.npz")
time = data['TIME']
div  = data['DIV']  # column 4 is cyanobacteria
locs = data['LOCS'] # column 4 is the log boom

tox  = data['TOX4'][:,4]
fbv  = data['FBV'][:,4,4]
temp = data['TEMP']
pwi  = data['PWI']
hum  = data['HUM']
pres  = data['PRES']
nut1 = data['NUT1'][:,4] # NO3+NO2
nut2 = data['NUT2'][:,4] # O-Phos
nut3 = data['NUT3'][:,4] # TN
nut4 = data['NUT4'][:,4] # T_Phos
# Algal diversity (Shannon)
# Light level
# Degree days

## Nuts can be negative?
nut1[nut1<0]=0
nut2[nut2<0]=0
nut3[nut3<0]=0
nut4[nut4<0]=0

## Log nuts
nut1 = np.log10(nut1+1);
nut2 = np.log10(nut2+1);
nut3 = np.log10(nut3+1);
nut4 = np.log10(nut4+1);

## Normalize data
tox  = tox  / np.nanmax(tox)
fbv  = fbv  / np.nanmax(fbv)
temp = temp / np.nanmax(temp)
pwi  = pwi  / np.nanmax(pwi)
hum  = hum  / np.nanmax(hum)
pres = pres / np.nanmax(pres)
nut1 = nut1 / np.nanmax(nut1)
nut2 = nut2 / np.nanmax(nut2)
nut3 = nut3 / np.nanmax(nut3)
nut4 = nut4 / np.nanmax(nut4)
D = np.vstack((fbv,temp,pwi,hum,pres,nut1,nut2,nut3,nut4)).transpose()

## Choose a time to plot
ds = np.sum(D,1)
JD = np.where(np.isnan(ds)==0)[0] # when all data is present

## Choose a day to test
i=1 # Choose which day to plot
ID = JD[i] # extract a day
dp = D[ID,:]

## Plot
arrCnts = (10*dp) # length of petals
iN = len(dp) # number of petals
lObjectsALLlbls = ['Cyano','Temp','Wind','Humidity','Pressure','NO3+NO2','O-Phos','TN','T-Phos']
theta=np.arange(0,2*np.pi,2*np.pi/iN)
width = (2*np.pi)/iN *0.9
bottom = 3.0

#! Make figu
fig = plt.figure(figsize=(8, 8))
ax = fig.add_axes([0.1, 0.1, 0.75, 0.75], polar=True)
bars = ax.bar(theta, arrCnts, width=width, bottom=bottom)
for r,bar in zip(arrCnts, bars):
    bar.set_facecolor(cm.jet(r/iN))
    bar.set_alpha(0.5)

# Clear up
ax.set_xticks(theta)
plt.axis('off')

# Add petal labels
rotations = np.rad2deg(theta)
y0,y1 = ax.get_ylim()

for x, bar, rotation, label in zip(theta, bars, rotations, lObjectsALLlbls):
     offset = (bottom+bar.get_height())/(y1-y0)
     lab = ax.text(0, 0, label, transform=None,
             ha='center', va='center')
     renderer = ax.figure.canvas.get_renderer()
     bbox = lab.get_window_extent(renderer=renderer)
     invb = ax.transData.inverted().transform([[0,0],[bbox.width,0] ])
     lab.set_position((x,offset+(invb[1][0]-invb[0][0])/2.*2.7 ) )
     lab.set_transform(ax.get_xaxis_transform())
     lab.set_rotation(rotation)

# Provide score
ax.annotate('Clear',
    xy=(0, 0),  # theta, radius
    horizontalalignment='center',
    verticalalignment='center',
    color='k',zorder=21,fontsize=30)    


#! Add current date
ax.set_theta_zero_location("N")
plt.show()
fig.tight_layout()
plt.show()
plt.savefig("./PNG/Fig_petal.png",dpi=600)

#a
## center point data
#bv = fbv[np.where(np.isnan(fbv)==0)]
#cn = fbv[ID]
#
### Plot
## force square figure and square axes looks better for polar, IMO
#fig = figure(figsize=(8,8))
#ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], polar=True)
#ax.grid(False)
##ax.spines['polar'].set_visible(False)
#
### Add colored central disk
##    if cn <= np.percentile(bv,33):
##        c = ax.scatter(0, 0, c=[.3,.5,1], s=4500, alpha=1,zorder=20)
##    elif (cn > np.percentile(bv,33)) & (cn <= np.percentile(bv,66)):
##        c = ax.scatter(0, 0, c=[.3,1,.3], s=4500, alpha=1,zorder=20)
##    else:
##        c = ax.scatter(0, 0, c=[0,1,0], s=4500, alpha=1,zorder=20)
#
### Add petals
#theta = np.arange(0.0, 2*np.pi, 2*np.pi/N) # angle of petals
#width = np.pi/4*np.ones(N) # width of petals
#bars = ax.bar(theta, radii, width=width, bottom=0.0)
#for r,bar in zip(radii, bars):
#    bar.set_facecolor(cm.jet(r/N))
#    #bar.set_facecolor( cm.jet(r/N))
#    bar.set_alpha(0.5)
#
#ax.set_rorigin(-2.5)
#ax.set_theta_direction(-1)
#ax.set_yticklabels([])
#ax.set_xticklabels([])
#ax.set_theta_zero_location('W', offset=10)
#
## Make a list of shifted thetas to place the labels at.
#Values = radii
#MetricLabels = ['Temp','Wind','Humid','Pressure','NO3/2', 'O-Phos', 'TN', 'I-Phos']
#Theta = theta
#ThetaShifted = np.copy(Theta)
#for i in range(N-1):
#    ThetaShifted[i] = (Theta[i] + Theta[i+1])/2.0
#ThetaShifted[-1] = (Theta[-1] + 2.0*np.pi)/2.0
#
#ax.set_xticklabels(['Temp','Wind','Humid','Pressure','NO3/2', 'O-Phos', 'TN', 'I-Phos'])
#
#ax.annotate('A',
#    xy=(0, 0),  # theta, radius
#    xytext=(0.05, 1.05), # fraction, fraction
#    horizontalalignment='center',
#    verticalalignment='center',
#    color='k',zorder=21,fontsize=40)    
#
#show()
#
#
#
#
#
#
#
#
#lObjectsALLcnts = [1, 1, 1, 2, 2, 3, 5, 14, 15, 20, 32, 33, 51, 1, 1, 2, 2, 3, 3, 3, 3, 3, 4, 6, 7, 7, 10, 10, 14, 14, 14, 17, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 5, 5, 6, 14, 14, 27, 27, 1, 1, 2, 3, 4, 4, 5]
#
#lObjectsALLlbls = ['DuctPipe', 'Column', 'Protrusion', 'Tree', 'Pole', 'Bar', 'Undefined', 'EarthingConductor', 'Grooves', 'UtilityPipe', 'Cables', 'RainPipe', 'Moulding', 'Intrusion', 'PowerPlug', 'UtilityBox', 'Balcony', 'Lighting', 'Lock', 'Doorbell', 'Alarm', 'LetterBox', 'Grate', 'Undefined', 'CableBox', 'Canopy', 'Vent', 'PowerBox', 'UtilityHole', 'Recess', 'Protrusion', 'Shutter', 'Handrail', 'Lock', 'Mirror', 'SecuritySpike', 'Bench', 'Intrusion', 'Picture', 'Showcase', 'Camera', 'Undefined', 'Stair', 'Protrusion', 'Alarm', 'Graffiti', 'Lighting', 'Ornaments', 'SecurityBar', 'Grate', 'Vent', 'Lighting', 'UtilityHole', 'Intrusion', 'Undefined', 'Protrusion']
#
#iN = len(lObjectsALLcnts)
#arrCnts = np.array(lObjectsALLcnts)
#
#theta=np.arange(0,2*np.pi,2*np.pi/iN)
#width = (2*np.pi)/iN *0.9
#
#fig = plt.figure(figsize=(8, 8))
#ax = fig.add_axes([0.1, 0.1, 0.75, 0.75], polar=True)
#bars = ax.bar(theta, arrCnts, width=width, bottom=20.0)
#for r,bar in zip(arrCnts, bars):
#    bar.set_facecolor(cm.jet(r/iN))
#    bar.set_alpha(0.5)
#
#ax.set_xticks(theta)
#plt.axis('off')
#
#
#bottom = 50
#rotations = np.rad2deg(theta)
#y0,y1 = ax.get_ylim()
#
#for x, bar, rotation, label in zip(theta, bars, rotations, lObjectsALLlbls):
#     offset = (bottom+bar.get_height())/(y1-y0)
#     lab = ax.text(0, 0, label, transform=None,
#             ha='center', va='center')
#     renderer = ax.figure.canvas.get_renderer()
#     bbox = lab.get_window_extent(renderer=renderer)
#     invb = ax.transData.inverted().transform([[0,0],[bbox.width,0] ])
#     lab.set_position((x,offset+(invb[1][0]-invb[0][0])/2.*2.7 ) )
#     lab.set_transform(ax.get_xaxis_transform())
#     lab.set_rotation(rotation)
#
#
