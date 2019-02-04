from netCDF4 import Dataset
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import glob
from scipy.interpolate import griddata
from datetime import datetime


### Load data
data = np.load("./Data/Data_LS8_timeseries.npz")
COL = data['COL']
TIME = data['TIME']
data = np.load("./Data/Data_lake_locations.npz")
LOCS = data['LOCS']


### Plot timeseries
import seaborn as sns
#sns.set_style("whitegrid")
sns.set_style("ticks")
sns.despine(offset=10, trim=True);
sns.set_context("talk")

plt.figure(figsize=[12,5])
plt.plot(TIM,COL[:,0,1],color=[.15,.9,.15])
plt.plot(TIM,COL[:,1,1])
plt.plot(TIM,COL[:,2,1])
plt.plot(TIM,COL[:,3,1])
plt.plot(TIM,COL[:,4,1])
plt.plot(TIM,COL[:,5,1])
plt.plot(TIM,COL[:,6,1])
plt.xlabel("Time")
plt.ylabel("Reflectance")
plt.tight_layout()
plt.savefig("./Figs/Fig_spectral_ts.png",dpi=600)

plt.figure(figsize=[12,5])
x = np.asarray([443,483,561,655,865,1609,2201])
y = np.asarray(COL[-18,:,1])
plt.plot(x,y,color=[.15,.9,.15])
plt.xlabel("Band nm")
plt.ylabel("Reflectance")
plt.title("May 2018")
plt.tight_layout()
plt.savefig("./Figs/Fig_spectrum.png",dpi=600)

### Plot maps
plt.figure(figsize=[9,7])
plt.pcolormesh(Lon,Lat,Col,cmap="rainbow")
plt.colorbar()
plt.clim([0,0.1])
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.title("L8 09/03/2018")
plt.tight_layout()
plt.savefig("./Figs/Fig_map_raw655.png",dpi=600)


