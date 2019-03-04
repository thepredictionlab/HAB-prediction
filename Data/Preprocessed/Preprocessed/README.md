### Preprocessed but cleaned and interpolated data

## In Data_historical
LOCS: locations in lake only ['BB','BO','HA','HT','LB','LBP','LBS']

TIME: decimal years (daily increments, at midday, from 2013)

NUT1: NO3+NO2 (mg/L)... time x locs

NUT2: O-Phos (mg/L)... time x locs

NUT3: TN (mg/L)... time x locs

NUT4: T-Phos (mg/L)... time x locs

TOX1: LCMSMS Cylindro (ppb)... time x locs

TOX2: LCMSMS Microcystin (ppb)... time x locs

TOX3: ELISA Cylindro (ppb)... time x locs

TOX4: ELISA Microcystin (ppb)... time x locs

DIV: bacterial/algal family name (there are 13)

DEN: algal concentration... time x locs x div

TBV: total biovolume... time x locs x div

FBV: fractional biovolume... time x locs x div

TEMP: temperature... time

HUM: humidity... time

PWI: peak wind speed... time

WIS: wind speed... time

RAIN: rain... time

PRES: barometric pressure... time


## In Data_sat_locations
The lat/lons of 53 locations defining the skeleton of the lake, created using 
ginput in python

## In Data_LS8_timeseries.npz
COL: the *color* of seven spectral bands:
BANDS: the spectral bands:
- 443
- 483
- 561
- 655
- 865
- 1609
- 2201
TIME: Time (decimal years)
