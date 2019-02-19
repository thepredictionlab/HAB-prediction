import numpy as np
import matplotlib.pylab as plt
from matplotlib.pyplot import figure, show, rc
import matplotlib.cm as cm



######## Climatology for long-term forecast
## Data
data = np.load("../Data/Preprocessed/Data_historical_7day.npz")
time = data['TIME']
tox  = data['TOX4'][:,4]
fbv  = data['FBV'][:,4,4]

TI = time - np.floor(time)
ID = np.where(np.isnan(fbv)==0)[0]
x  = TI[ID]
y  = fbv[ID]
ID = np.where(y>20)
x  = x[ID]
y  = y[ID]
ID = np.where(x != x.min())[0]
x  = x[ID]
y  = y[ID]

# Bin to monthly timescale
from scipy import interpolate
X = np.arange((1/12)/2,1,1/12)
Y  = np.arange(0,len(X))
f = interpolate.interp1d(X, Y,'nearest')
yi = f(x)

MU = np.zeros(len(X))
for i in np.unique(yi):
    II = len(np.where(yi==i)[0])
    MU[int(i)] = II
MU += 1

######### Plot
import seaborn as sns
sns.set(style="whitegrid")

fig, ax = plt.subplots(figsize=(12,5))
N = 12
x = np.arange(N).astype(float)
x += .5
y = MU / MU.max()
y = (1-y)*10
mn = np.asarray(['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'])

#! Plot bars
pal = sns.color_palette("Greens_r", len(MU))
rank = MU.argsort().argsort()   # http://stackoverflow.com/a/6266510/1628638
sns.barplot(x=mn, y=MU, palette=np.array(pal[::-1])[rank])

# Highlight which month we are in
colors = ['k','r', 'k','k','k','k','k','k','k','k','k','k']
for xtick, color in zip(ax.get_xticklabels(), colors):
    xtick.set_color(color)

#set ticks
T =  np.asarray([1,5,9,13])
ax.set_yticks(T)

#set labels
labels=['0','4','8','12']
ax.set_yticklabels(labels)

#! Adorn
plt.ylabel('Number of blooms (since 2013)',fontsize=16)
plt.xlabel('Month of the year',fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)


#! Add current date
fig.tight_layout()
plt.show()
plt.savefig("./PNG/Fig_climatology.png",dpi=600)



