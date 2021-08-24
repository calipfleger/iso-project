#!/usr/bin/env python
"""
Created on Thu Jul 22 11:22:16 2021
@author: andrew
Iso2k Coral PRYSM data prep and model run
"""

import numpy as np
import pandas as pd
from psm.coral.sensor import pseudocoral
from psm.agemodels.banded import bam_simul_perturb
from psm.aux_functions.analytical_error import analytical_error
from psm.aux_functions.analytical_err_simple import analytical_err_simple

#====================================================================
# 1. LOAD CLIMATE FIELDS
#====================================================================

#Data directory
datadir='./data/iCESM_full_forcing_afedits/'

df = pd.read_csv(datadir + 'iBOFP_mean_afedit.csv', sep = ',')
sst = np.array(df["TS"])
# Load SST anomalies [K] (NOTE: THIS SHOULD BE A 1-D VECTOR OF DATA!)
# yearly
df_sst_a = pd.read_csv(datadir + 'ssta_BOFPiM_binned.csv', sep = ',')
ssta = np.array(df_sst_a["y"])

# monthly
ssta_m = np.array(df["TS"])
# Export monthly data
ssta_m = pd.DataFrame(ssta_m)
ssta_m.to_csv(datadir+"SSTA_BOFPiM_M.csv", index = False, header = False)

# Load SSS anomalies [psu] (NOTE: THIS SHOULD BE A 1-D VECTOR OF DATA!)
# yearly
df_sss_a = pd.read_csv(datadir + 'ssta_BOFPiM_binned.csv', sep = ',')
sssa = np.array(df_sss_a["y"])

# monthly
sssa_m = np.array(df["SALT_mean"])
# Export monthly data
sssa_m = pd.DataFrame(sssa_m) 
sssa_m.to_csv(datadir+"SSSA_BOFPiM_M.csv", index = False, header = False)

# set time axis
step = 1
start = 1851
stop = 1981
num = (stop-start)*step
time=np.linspace(start, stop, num = num, endpoint = False)

print('Preparing data...')
# 2.1 Convert Lats and Lons to standard format (if they aren't already).

# lon: convert (-180: +180) to (0: 360)
# lat: convert (0: 180) to (-90: +90)

# Enter your coordinates here:
lon=210.17
lat=-17.5

# make sure there are no negative-longitude coordinates: [0 to 360] only.
# longitude
if (lon<0.):
    lon = lon+360.

# latitude
if (lat>90.):
    lat = lat-90.

# 2.2 Ensure that Sea-Surface Temperature is in degrees C, not Kelvin.
sst=ssta #AAA
sss=sssa #AAA
temp_flag = any(sst>200)

for i in range(len(sst)):
    if (temp_flag):
        sst[i] = sst[i]-274.15

#====================================================================
# 3. CALL SENSOR MODEL: Call function 'pseudocoral' (see doc) to compute pseudocoral
#    timeseries and spatial field.
#====================================================================

# NOTE: THIS SHOULD BE A 1-D VECTOR OF DATA!
print('Running sensor model...')
coral = np.zeros(len(time))  # this will initialize a [Time x Lat x Lon] matrix of coral values.

# Fill coral array with data same size as input vectors.
# Important: if you have d18O of seawater, add a 5th argument 'd18O' and input data as vector array (default is -1).
# Inputs: lat, lon, SST, SSS OR d18O
# Optional d18O/T slope parameters:
# a= -0.22,b1=0.3007062,b2=0.2619054,b3=0.436509,b4=0.1552032 (see doctring)

for i in range(len(time)):
    coral[i] = pseudocoral(lat,lon,sst[i],sss[i])

#======================================================================
# 4. CALL OBSERVATION MODEL
#======================================================================

# 4.1 Specify and model rate of annual layer miscount: BAM (see doctring)
print('Running observation model...')
X = coral
X = X.reshape(len(X),1)
tp, Xp, tmc=bam_simul_perturb(X,time,param=[0.02,0.02],name='poisson',ns=1000,resize=0)
#======================================================================
# 4.2: Analytical Uncertainty Model:

#4.2.1 Simple Model: just add uncertainty bands based on measurement precision
sigma=0.1 # permil, measurement  precision
coral_upper, coral_lower = analytical_err_simple(X,sigma)

#4.2.2 Gaussian Noise Model for analytical error:
sigma=0.1
#nsamples = ## enter number of samples here
coral_Xn=analytical_error(X,sigma)
#====================================================================
# Save coral timeseries fields as BOFPy arrays in current directory.
outdir = './results/'
print('Saving time series...')
np.save(outdir + "BOFPiM_coral_d18O.npy", coral)
np.save(outdir + "BOFPiM_coral_age_perturbed.npy", Xp)
coral = pd.DataFrame(coral)
coral.to_csv(outdir+"BOFPiM_pseudocoral.csv", index = False, header = False)
#coral_error_bounds=np.save('coral_error.npy',coral_upper, coral_lower)
#====================================================================
