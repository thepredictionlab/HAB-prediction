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
D = np.vstack((temp,pwi,hum,pres,nut1,nut2,nut3,nut4)).transpose()

## Choose a time to plot
ds = np.sum(D,1)
JD = np.where(np.isnan(ds)==0)[0] # when all data is present

bgcolor = '#222222'


for i in np.arange(0,1):
    ID = JD[i] # choose a day
    dp = D[ID,:]
    radii = (10*dp) # length of petals
    N = len(dp) # number of petals

# center point data
    bv = fbv[np.where(np.isnan(fbv)==0)]
    cn = fbv[ID]

## Plot
# force square figure and square axes looks better for polar, IMO
    fig = figure(figsize=(8,8))
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], polar=True)
    ax.grid(False)
    #ax.spines['polar'].set_visible(False)

## Add central point
#    if cn <= np.percentile(bv,33):
#        c = ax.scatter(0, 0, c=[.3,.5,1], s=4500, alpha=1,zorder=20)
#    elif (cn > np.percentile(bv,33)) & (cn <= np.percentile(bv,66)):
#        c = ax.scatter(0, 0, c=[.3,1,.3], s=4500, alpha=1,zorder=20)
#    else:
#        c = ax.scatter(0, 0, c=[0,1,0], s=4500, alpha=1,zorder=20)

## Add petals
    theta = np.arange(0.0, 2*np.pi, 2*np.pi/N) # angle of petals
    width = np.pi/4*np.ones(N) # width of petals
    bars = ax.bar(theta, radii, width=width, bottom=0.0)
    for r,bar in zip(radii, bars):
        bar.set_facecolor( rvb(r/N))
        #bar.set_facecolor( cm.jet(r/N))
        bar.set_alpha(0.5)

    ax.set_rorigin(-2.5)
    ax.set_theta_direction(-1)
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    ax.set_theta_zero_location('W', offset=10)

# Make a list of shifted thetas to place the labels at.
    Values = radii
    MetricLabels = ['Temp','Wind','Humid','Pressure','NO3/2', 'O-Phos', 'TN', 'I-Phos']
    Theta = theta
    ThetaShifted = np.copy(Theta)
    for i in range(N-1):
        ThetaShifted[i] = (Theta[i] + Theta[i+1])/2.0
    ThetaShifted[-1] = (Theta[-1] + 2.0*np.pi)/2.0

    ax.set_xticklabels(['Temp','Wind','Humid','Pressure','NO3/2', 'O-Phos', 'TN', 'I-Phos'])

    ax.annotate('A',
        xy=(0, 0),  # theta, radius
        xytext=(0.05, 1.05), # fraction, fraction
        horizontalalignment='center',
        verticalalignment='center',
        color='k',zorder=21,fontsize=40)    

show()



####### PREDICTION PLOT
import matplotlib.colors as mcolors
def make_colormap(seq):
    """Return a LinearSegmentedColormap
    seq: a sequence of floats and RGB-tuples. The floats should be increasing
    and in the interval (0,1).
    """
    seq = [(None,) * 3, 0.0] + list(seq) + [1.0, (None,) * 3]
    cdict = {'red': [], 'green': [], 'blue': []}
    for i, item in enumerate(seq):
        if isinstance(item, float):
            r1, g1, b1 = seq[i - 1]
            r2, g2, b2 = seq[i + 1]
            cdict['red'].append([item, r1, r2])
            cdict['green'].append([item, g1, g2])
            cdict['blue'].append([item, b1, b2])

    return mcolors.LinearSegmentedColormap('CustomMap', cdict)

c = mcolors.ColorConverter().to_rgb
rvb = make_colormap(
[c('red'), 0.125, c('red'), c('orange'), 0.25, c('orange'),c('green'),0.5, c('green'),0.7, c('green'), c('blue'), 0.75, c('blue')])

objects = ('Worse alot', 'Worse a little', 'Same', 'Better a little', 'Better alot')
N = len(objects)
y_pos = len(objects) - np.arange(len(objects)) - 1
probability = [.7,.1,.05,.025,.025]

plt.figure()
plt.barh(y_pos, probability, align='center', alpha=0.5,color=rvb((np.abs(y_pos-N)-1)/N))
plt.yticks(y_pos, objects)
plt.xlabel('Probability')
plt.title('HAB Condition Forecast')



######## Climatology for long-term forecast
# Get data
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

# Make colorscale
def make_colormap(seq):
    """Return a LinearSegmentedColormap
    seq: a sequence of floats and RGB-tuples. The floats should be increasing
    and in the interval (0,1).
    """
    seq = [(None,) * 3, 0.0] + list(seq) + [1.0, (None,) * 3]
    cdict = {'red': [], 'green': [], 'blue': []}
    for i, item in enumerate(seq):
        if isinstance(item, float):
            r1, g1, b1 = seq[i - 1]
            r2, g2, b2 = seq[i + 1]
            cdict['red'].append([item, r1, r2])
            cdict['green'].append([item, g1, g2])
            cdict['blue'].append([item, b1, b2])

    return mcolors.LinearSegmentedColormap('CustomMap', cdict)


c = mcolors.ColorConverter().to_rgb
rvb = make_colormap(
[c('red'), 0.125, c('red'), c('orange'), 0.25, c('orange'),c('green'),0.5, c('green'),0.7, c('green'), c('blue'), 0.75, c('blue')])

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

#set ticks
T =  np.asarray([1,5,9,13])
ax.set_yticks(T)

#set labels
labels=['0','4','8','12']
ax.set_yticklabels(labels)

#! Adorn
plt.ylabel('Number of blooms',fontsize=14)
plt.xlabel('Month of the year',fontsize=14)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)


#! Add current date
fig.tight_layout()
plt.show()




