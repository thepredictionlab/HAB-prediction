from netCDF4 import Dataset
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import glob
from scipy.interpolate import griddata
from datetime import datetime


### Files
files = sorted(glob.glob("./Data/LS8/*.nc"))
TIM = np.empty(len(files),dtype="object")
data = np.load("./Data/Data_lake_locations.npz")
LOCS = data['LOCS']
COL = np.zeros((len(files),7,len(LOCS)))


### Convert time to decimal years
from datetime import datetime, timedelta
def datetime2year(dt):
    year_part = dt - datetime(year=dt.year, month=1, day=1)
    year_length = datetime(year=dt.year+1, month=1, day=1) - datetime(year=dt.year, month=1, day=1)
    return dt.year + year_part/year_length
TIME = np.zeros(len(TIM))
for i in np.arange(0,len(TIM)):
    dt = TIM[i] + timedelta(hours=12)
    TIME[i] = datetime2year(TIM[i])

### Save info about bands
BANDS = np.asarray(["443","483","561","655","865","1609","2201"])

### Loop through Data
t = 0
for file in files:

    # get time
    yr = int(file[18:22])
    mn = int(file[23:25])
    dy = int(file[26:28])
    TIM[t] = datetime(yr,mn,dy)

    # data
    data = Dataset(file)
    Lon = np.asarray(data['lon'])
    Lat = np.asarray(data['lat'])
    Col1 = np.asarray(data['rhorc_443'])
    Col2 = np.asarray(data['rhorc_483'])
    Col3 = np.asarray(data['rhorc_561'])
    Col4 = np.asarray(data['rhorc_655'])
    Col5 = np.asarray(data['rhorc_865'])
    Col6 = np.asarray(data['rhorc_1609'])
    Col7 = np.asarray(data['rhorc_2201'])

    # get specific time series
    ID = np.arange(0,len(Lon.flatten()))
    ID = np.reshape(ID,Lon.shape)

    # place to interp to
    #lat_i = 44.695236
    #lon_i =-122.207523
    #lat_i = 44.704206
    #lon_i = -122.237333
    lon_i = LOCS[:,0]
    lat_i = LOCS[:,1]

    # find nearest pixel
    JD = griddata((Lon.flatten(),Lat.flatten()),ID.flatten(),(lon_i,lat_i),'nearest')

    # Loop through all locations
    #Mat = np.ones(Col1.shape) # for testing
    for i in np.arange(0,len(JD)):
        KD = np.where(ID==JD[i])
        Mat[KD] = 0 # for testing

        # Make timeseris
        COL[t,0,i] = np.float(Col1[KD])
        COL[t,1,i] = np.float(Col2[KD])
        COL[t,2,i] = np.float(Col3[KD])
        COL[t,3,i] = np.float(Col4[KD])
        COL[t,4,i] = np.float(Col5[KD])
        COL[t,5,i] = np.float(Col6[KD])
        COL[t,6,i] = np.float(Col7[KD])

    t+=1
    print(t)


### Clear (remove neg numbers)
COL[COL<=0] = 1e-10


### Save
np.savez("./Data/Data_LS8_timeseries.npz",COL=COL,TIME=TIME,BANDS=BANDS)
