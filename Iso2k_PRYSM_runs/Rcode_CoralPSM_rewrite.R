library(geoChronR)
library(magrittr)
library(tidyverse)
library(lubridate)
library(rgdal)
library(ncdf4)
library(ncdf4.helpers)

#
#### Import iCESM data ####
nrec = 7
setwd("~/GitHub/iso-project/Iso2k_PRYSM_runs/data/iCESM_full_forcing_afedits")
# SST
sst = matrix(data = NA, nrow = 130, ncol = nrec)
sst[,1] = read.csv("ssta_KUBEiM_binned.csv", head = T)[,3]
sst[,2] = read.csv("ssta_SWBBiM_binned.csv", head = T)[,3]
sst[,3] = read.csv("ssta_KIPRiM_binned.csv", head = T)[,3]
sst[,4] = read.csv("ssta_LIRAiM_binned.csv", head = T)[,3]
sst[,5] = read.csv("ssta_NUMPiM_binned.csv", head = T)[,3]
sst[,6] = read.csv("ssta_COPMiM_binned.csv", head = T)[,3]
sst[,7] = read.csv("ssta_BOFPiM_binned.csv", head = T)[,3]
# SSS
sss = matrix(data = NA, nrow = 130, ncol = nrec)
sss[,1] = read.csv("sssa_KUBEiM_binned.csv", head = T)[,3]
sss[,2] = read.csv("sssa_SWBBiM_binned.csv", head = T)[,3]
sss[,3] = read.csv("sssa_KIPRiM_binned.csv", head = T)[,3]
sss[,4] = read.csv("sssa_LIRAiM_binned.csv", head = T)[,3]
sss[,5] = read.csv("sssa_NUMPiM_binned.csv", head = T)[,3]
sss[,6] = read.csv("sssa_COPMiM_binned.csv", head = T)[,3]
sss[,7] = read.csv("sssa_BOFPiM_binned.csv", head = T)[,3]
#
#### Extract and bin Iso2k records####
# Time dimension
start = 1851
end = 1980
time = seq(from = start, to = end, by = 1)
lat = matrix(data = NA, nrow = 1, ncol = nrec)
lon = matrix(data = NA, nrow = 1, ncol = nrec)
setwd("C:/Users/andre/Box Sync/Konecky Lab/User Storage/Andrew Flaim/Rcode/Iso2k") # Laptop
load("iso2k1_0_0.RData")
# Remove extraneous objects
rm(D, TS)
start = 1851
end = 1981
binvec =  seq(from = start, to = end, by = 1)
# Apply geoChronR function to extract primary TS
# SWBB
SWBB_data = sTS[which(pullTsVariable(sTS, 
                                     variable = "dataSetName") == "CO96SWBB")][[1]]
lat[2] = SWBB_data$geo_latitude
lon[2] = SWBB_data$geo_longitude
SWBB_rec = na.omit(data.frame(year = SWBB_data$year, 
                              val = SWBB_data$paleoData_values))
binned_SWBB = geoChronR::bin(SWBB_rec$year, SWBB_rec$val, binvec)
SWBB_anom = binned_SWBB$y-mean(na.omit(binned_SWBB$y))

# KUBE
KUBE_data = sTS[which(pullTsVariable(sTS, 
                                     variable = "dataSetName") == "CO05KUBE")][[1]]
lat[1] = KUBE_data$geo_latitude
lon[1] = KUBE_data$geo_longitude
KUBE_rec = na.omit(data.frame(year = KUBE_data$year, 
                              val = KUBE_data$paleoData_values))
binned_KUBE = geoChronR::bin(KUBE_rec$year, KUBE_rec$val, binvec)
KUBE_anom = binned_KUBE$y-mean(na.omit(binned_KUBE$y))

# LIRA
LIRA_data = sTS[which(pullTsVariable(sTS, 
                                     variable = "dataSetName") == "CO04LIRA")][[1]]
lat[4] = LIRA_data$geo_latitude
lon[4] = LIRA_data$geo_longitude
LIRA_rec = na.omit(data.frame(year = LIRA_data$year, 
                              val = LIRA_data$paleoData_values))
binned_LIRA = geoChronR::bin(LIRA_rec$year, LIRA_rec$val, binvec)
LIRA_anom = binned_LIRA$y-mean(na.omit(binned_LIRA$y))

# NUMP
NUMP_data = sTS[which(pullTsVariable(sTS, 
                                     variable = "dataSetName") == "CO11NUPM")][[5]]
lat[5] = NUMP_data$geo_latitude
lon[5] = NUMP_data$geo_longitude
NUMP_rec = na.omit(data.frame(year = NUMP_data$year, 
                              val = NUMP_data$paleoData_values))
binned_NUMP = geoChronR::bin(NUMP_rec$year, NUMP_rec$val, binvec)
NUMP_anom = binned_NUMP$y-mean(na.omit(binned_NUMP$y))

# COPM
COPM_data = sTS[which(pullTsVariable(sTS, 
                                     variable = "dataSetName") == "CO03COPM")][[1]]
lat[6] = COPM_data$geo_latitude
lon[6] = COPM_data$geo_longitude
COPM_rec = na.omit(data.frame(year = COPM_data$year, 
                              val = COPM_data$paleoData_values))
binned_COPM = geoChronR::bin(COPM_rec$year, COPM_rec$val, binvec)
COPM_anom = binned_COPM$y-mean(na.omit(binned_COPM$y))

# BOFP
BOFP_data = sTS[which(pullTsVariable(sTS, 
                                     variable = "dataSetName") == "CO98BOFP")][[1]]
lat[7] = BOFP_data$geo_latitude
lon[7] = BOFP_data$geo_longitude
BOFP_rec = na.omit(data.frame(year = BOFP_data$year, 
                              val = BOFP_data$paleoData_values))
binned_BOFP = geoChronR::bin(BOFP_rec$year, BOFP_rec$val, binvec)
BOFP_anom = binned_BOFP$y-mean(na.omit(binned_BOFP$y))

# KIPR
KIPR_data = sTS[which(pullTsVariable(sTS, 
                                     variable = "dataSetName") == "CO08KIPR")][[1]]
lat[3] = KIPR_data$geo_latitude
lon[3] = KIPR_data$geo_longitude
KIPR_rec = na.omit(data.frame(year = KIPR_data$year, 
                              val = KIPR_data$paleoData_values))
binned_KIPR = geoChronR::bin(KIPR_rec$year, KIPR_rec$val, binvec)
KIPR_anom = binned_KIPR$y-mean(na.omit(binned_KIPR$y))

#

#### Preparing data ####
# 2.1 Convert Lats and Lons to standard format (if they aren't already).
# Enter your coordinates here:

# make sure there are no negative-longitude coordinates: [0 to 360] only.
# longitude
for (i in 1:nrec){
  if (lon[i] < 0.){
    lon[i] = lon[i]+360.
  }
  # latitude
  if (lat[i]>90.){
    lat[i] = lat[i]-90.
  }
}


# 2.2 Ensure that Sea-Surface Temperature is in degrees C, not Kelvin.
temp_flag = any(sst>200)

for (i in 1:nrec){
  for (j in 1:length(time)){
    if (temp_flag){
      sst[i, j] = sst[i, j]-274.15
    }
  }
}


#### Assign beta values ####
b1=0.3007062
b2=0.2619054
b3=0.436509
b4=0.1552032
b5=0.15

V_corr = 0.97002

slope1 = V_corr*b1      # Red Sea
slope2 = V_corr*b2      # Tropical Pacific
slope3 = V_corr*b3      # South Pacific
slope4 = V_corr*b4      # Indian Ocean
slope5 = V_corr*b5      # Tropical Atlantic

# Given Lat/Lon pair, assign correct value to 'b'
b = matrix(data = NA, nrow = 1, ncol = nrec)
for (i in 1:nrec){
  if (lon[i]>=32.83 & lon[i]<=43.5 & lat[i]>=12.38 & lat[i]<=28.5){
    b[i] = slope1 #Red Sea
  } else if (lon[i]<=120.){
    b[i] = slope4 #Indian Ocean
  } else if (lon[i]>=270. & lat[i]>=10.){
    b[i] = slope5 #Tropical Atlantic ~300 E longitude
  } else if (lat[i]> -5. & lat[i]<=13.){
    b[i] = slope2 #Tropical Pacific
  } else if (lat[i]<= -5.){
    b[i] = slope3 #South Pacific
  } else{
    b[i] = slope2 #Default = Trop. Pacific.
  }
}

#### PSM ####
a = -0.22
coral = matrix(data = NA, nrow = 130, ncol = nrec)
for (i in 1:nrec){
  for (j in 1:length(time)){
    coral[j, i] = a*sst[j, i] + b[i]*sss[j, i]
  }
}

#### Plot Pseudocorals ####
rec_col = c("blue", "red", "dark green", "purple", "black", "orange", "brown")
for(i in 1:nrec){
  plot(time, coral[,i], type = "l", frame.plot = F, ylim = c(-0.3,0.3),
       xlab = "", ylab = "", axes = F, col = rec_col[i])
  par(new = T)
}
plot(NA, frame.plot = F, ylim = c(-0.2,0.2), xlim = c(min(time), max(time)),
     ylab = "Pseudo d18O anomalies", xlab = "Years (CE)", axes = T, main = "Pseudocorals")
#### Plot Atlantic Output and PSMs ####
# Plot SST
for(i in 1:3){
  plot(time, sst[,i], type = "l", frame.plot = F, ylim = c(-1,1),
       xlab = "", ylab = "", axes = F, col = rec_col[i])
  par(new = T)
}
plot(NA, frame.plot = F, ylim = c(-1,1), xlim = c(min(time), max(time)),
     ylab = "SST anomalies", xlab = "Years (CE)", axes = T, main = "Atlantic iCESM full forcing SSTs")
legend("topright", col = rec_col, 
       legend = c("KUBE", "SWBB", "KIPR"),
       lty = rep(1, 3))

# Plot SSS
for(i in 1:3){
  plot(time, sss[,i], type = "l", frame.plot = F, ylim = c(-0.5,0.5),
       xlab = "", ylab = "", axes = F, col = rec_col[i])
  par(new = T)
}
plot(NA, frame.plot = F, ylim = c(-0.5,0.5), xlim = c(min(time), max(time)),
     ylab = "SSS anomalies", xlab = "Years (CE)", axes = T, main = "Atlantic iCESM full forcing SSS")

# Plot pseudoproxies
for(i in 1:3){
  plot(time, coral[,i], type = "l", frame.plot = F, ylim = c(-0.3,0.3),
       xlab = "", ylab = "", axes = F, col = rec_col[i])
  par(new = T)
}
plot(NA, frame.plot = F, ylim = c(-0.3,0.3), xlim = c(min(time), max(time)),
     ylab = "Pseudo d18O anomalies", xlab = "Years (CE)", axes = T, main = "Atlantic Pseudocorals")

# Plot proxies
plot(time, KUBE_anom, type = "l", frame.plot = F, xlim = c(min(time), max(time)),
     ylab = "", xlab = "", axes = T, main = "", col = rec_col[1],
     ylim = c(-0.8, 0.8))
par(new = T)
plot(time, SWBB_anom, type = "l", frame.plot = F, xlim = c(min(time), max(time)),
     ylab = "", xlab = "", axes = F, main = "", col = rec_col[2],
     ylim = c(-0.8, 0.8))
par(new = T)
plot(time, KIPR_anom, type = "l", frame.plot = F, xlim = c(min(time), max(time)),
     ylab = "d18O anomalies", xlab = "Years (CE)", axes = T, main = "Atlantic Iso2k coral anomalies",
     col = rec_col[3], ylim = c(-0.8, 0.8))
#
#### Plot Pacific Output and PSMs ####
# Plot SST
for(i in 4:nrec){
  plot(time, sst[,i], type = "l", frame.plot = F, ylim = c(-1,1),
       xlab = "", ylab = "", axes = F, col = rec_col[i])
  par(new = T)
}
plot(NA, frame.plot = F, ylim = c(-1,1), xlim = c(min(time), max(time)),
     ylab = "SST anomalies", xlab = "Years (CE)", axes = T, main = "Pacific iCESM full forcing SSTs")
legend("topright", col = rec_col[4:7], 
       legend = c("LIRA", "NUMP", "COPM", "BOFP"),
       lty = rep(1, 4))

# Plot SSS
for(i in 4:nrec){
  plot(time, sss[,i], type = "l", frame.plot = F, ylim = c(-0.5,0.5),
       xlab = "", ylab = "", axes = F, col = rec_col[i])
  par(new = T)
}
plot(NA, frame.plot = F, ylim = c(-0.5,0.5), xlim = c(min(time), max(time)),
     ylab = "SSS anomalies", xlab = "Years (CE)", axes = T, main = "Pacific iCESM full forcing SSS")

# Plot pseudoproxies
for(i in 4:nrec){
  plot(time, coral[,i], type = "l", frame.plot = F, ylim = c(-0.3,0.3),
       xlab = "", ylab = "", axes = F, col = rec_col[i])
  par(new = T)
}
plot(NA, frame.plot = F, ylim = c(-0.3,0.3), xlim = c(min(time), max(time)),
     ylab = "Pseudo d18O anomalies", xlab = "Years (CE)", axes = T, main = "Pacific Pseudocorals")

# Plot proxies
plot(time, LIRA_anom, type = "l", frame.plot = F, xlim = c(min(time), max(time)),
     ylab = "", xlab = "", axes = T, main = "", col = rec_col[4],
     ylim = c(-0.8, 0.8))
par(new = T)
plot(time, NUMP_anom, type = "l", frame.plot = F, xlim = c(min(time), max(time)),
     ylab = "", xlab = "", axes = F, main = "", col = rec_col[5],
     ylim = c(-0.8, 0.8))
par(new = T)
plot(time, COPM_anom, type = "l", frame.plot = F, xlim = c(min(time), max(time)),
     ylab = "", xlab = "", axes = F, main = "", col = rec_col[6],
     ylim = c(-0.8, 0.8))
par(new = T)
plot(time, BOFP_anom, type = "l", frame.plot = F, xlim = c(min(time), max(time)),
     ylab = "d18O anomalies", xlab = "Years (CE)", axes = T, main = "Pacific Iso2k coral anomalies",
     col = rec_col[7], ylim = c(-0.8, 0.8))
#
#### Extract observational SSTs ####
setwd("C:/Users/andre/Box Sync/Konecky Lab/User Storage/Andrew Flaim/Rcode/Data") # Laptop

data = nc_open("sst.mon.mean.nc")
sst = ncvar_get(data, "sst")

lon = ncvar_get(data, varid = "lon")
lat = ncvar_get(data, varid = "lat")
sst_time = nc.get.time.series(data, v = "sst",
                              time.dim.name = "time")
# Extract grided data
time_index = which(format(sst_time, "%Y-%m-%d") == "1980-12-01")
sst.slice = sst[,,time_index]

bofp = sst[210, -17.5, ]

sst.location = bofp
trim = 1:1573
trim_time = sst_time[trim]

test = tibble(time = trim_time, 
              sst = as.vector(sst.location[trim])) %>% 
  mutate(time = as.Date(format(time, "%Y-%m-%d")))

# Bin timeseries
binvec =  seq(from = 1850, to = 1980, by = 1)
year = seq(from = 1850, to = 1981, by = 1/12)
thisrec = na.omit(data.frame(year = year,
                             val = test$sst))
binned_ann_recs = geoChronR::bin(time = thisrec$year, 
                                 val = thisrec$val, bin.vec = binvec)
sst_obs = binned_ann_recs$y-mean(na.omit(binned_ann_recs$y))

nc_close("sst.mon.mean.nc")
#
#### BOFP Moorea, French Polynesia - Fixed d18O with observational temp ####
# d18Osw = 0.6 from Conroy 2017
d18O = 0.6
a = -0.22
fixed = matrix(data = NA, nrow = length(time), ncol = 1)

for (i in 1:length(time)){
  fixed[i,] = a*sst_obs[i] + d18O
}

# Fixed with modeled temp
fixed_mod = matrix(data = NA, nrow = length(time), ncol = 1)
for (i in 1:length(time)){
  fixed_mod[i,] = a*sst[i,7] + d18O
}
#### Plot fixed obs. input and output ####
par(mfrow = c(4,1), mai = c(0, 0, 0, 0), omi = c(0.75, 0.75, 0.25, 0))
plot(time, sst[,7], type = "l", frame.plot = F, xaxt = "n", xlab = "", xpd = NA,
     ylab = "iCESM SST", cex.lab = 1.5, cex.axis = 1.5)
plot(time, fixed_mod, type = "l", frame.plot = F, xaxt = "n", xlab = "", xpd = NA,
     ylab = "iCESM Pseudo d18O", cex.lab = 1.5, cex.axis = 1.5)
plot(time, sst_obs, type = "l", frame.plot = F, xaxt = "n", xlab = "", xpd = NA,
     ylab = "COBE-SST", cex.lab = 1.5, cex.axis = 1.5)
plot(time, fixed, type = "l", frame.plot = F, xlab = "Years (CE)", xpd = NA,
     ylab = "Obs. Pseudo d18O", cex.lab = 1.5, cex.axis = 1.5)

# Fixed model, obs, and proxy
par(mfrow = c(3,1), mai = c(0, 0, 0, 0), omi = c(0.75, 0.75, 0.25, 0))
plot(time, BOFP_anom, type = "l", frame.plot = F, xaxt = "n", xlab = "", xpd = NA,
     ylab = "Iso2k d18O", cex.lab = 1.5, cex.axis = 1.5)
plot(time, fixed_mod, type = "l", frame.plot = F, xaxt = "n", xlab = "", xpd = NA,
     ylab = "iCESM Pseudo d18O", cex.lab = 1.5, cex.axis = 1.5)
plot(time, fixed, type = "l", frame.plot = F, xlab = "Years (CE)", xpd = NA,
     ylab = "Obs. Pseudo d18O", cex.lab = 1.5, cex.axis = 1.5)