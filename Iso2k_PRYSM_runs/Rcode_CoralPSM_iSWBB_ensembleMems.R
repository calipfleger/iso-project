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
setwd("~/GitHub/iso-project/full_forcing_files/icesm_fullforcing_indv")

iSWBB_1 = read.csv("iSWBB_1.csv", head = T)
iSWBB_2 = read.csv("iSWBB_2.csv", head = T)
iSWBB_3 = read.csv("iSWBB_3.csv", head = T)

iSWBB_1 = iSWBB_1[12:1572,]
iSWBB_2 = iSWBB_2[12:1572,]
iSWBB_3 = iSWBB_3[12:1572,]

iSWBB_1_sss = iSWBB_1[,2]
iSWBB_2_sss = iSWBB_2[,2]
iSWBB_3_sss = iSWBB_3[,2]

iSWBB_1_sst = iSWBB_1[,3]
iSWBB_2_sst = iSWBB_2[,3]
iSWBB_3_sst = iSWBB_3[,3]

#### Bin sst and sss ####
start_bin = 1851
end_bin = 1981
time_bin = seq(from = start_bin, to = end_bin, by = 1)
# Binned sst anomaly matrix
sst_binned = matrix(data = NA, nrow = 130, ncol = 3)
# Bin SWBB_1_sst
binvec = seq(from = start_bin, to = end_bin, by = 1/12)
iSWBB_1_sst_thisrec = na.omit(data.frame(year = binvec, val = iSWBB_1_sst))
iSWBB_1_sst_bin = geoChronR::bin(iSWBB_1_sst_thisrec$year, 
                                 iSWBB_1_sst_thisrec$val, time_bin)
sst_binned[,1] = iSWBB_1_sst_bin$y - mean(iSWBB_1_sst_bin$y)
# Bin SWBB_2_sst
iSWBB_2_sst_thisrec = na.omit(data.frame(year = binvec, val = iSWBB_2_sst))
iSWBB_2_sst_bin = geoChronR::bin(iSWBB_2_sst_thisrec$year, 
                                 iSWBB_2_sst_thisrec$val, time_bin)
sst_binned[,2] = iSWBB_2_sst_bin$y - mean(iSWBB_2_sst_bin$y)
# Bin SWBB_1_sst
iSWBB_3_sst_thisrec = na.omit(data.frame(year = binvec, val = iSWBB_3_sst))
iSWBB_3_sst_bin = geoChronR::bin(iSWBB_3_sst_thisrec$year, 
                                 iSWBB_3_sst_thisrec$val, time_bin)
sst_binned[,3] = iSWBB_3_sst_bin$y - mean(iSWBB_3_sst_bin$y)

# Binned sss anomaly matrix
sss_binned = matrix(data = NA, nrow = 130, ncol = 3)
# Bin SWBB_1_sss
binvec = seq(from = start_bin, to = end_bin, by = 1/12)
iSWBB_1_sss_thisrec = na.omit(data.frame(year = binvec, val = iSWBB_1_sss))
iSWBB_1_sss_bin = geoChronR::bin(iSWBB_1_sss_thisrec$year, 
                                 iSWBB_1_sss_thisrec$val, time_bin)
sss_binned[,1] = iSWBB_1_sss_bin$y - mean(iSWBB_1_sss_bin$y)
# Bin SWBB_2_sss
iSWBB_2_sss_thisrec = na.omit(data.frame(year = binvec, val = iSWBB_2_sss))
iSWBB_2_sss_bin = geoChronR::bin(iSWBB_2_sss_thisrec$year, 
                                 iSWBB_2_sss_thisrec$val, time_bin)
sss_binned[,2] = iSWBB_2_sss_bin$y - mean(iSWBB_2_sss_bin$y)
# Bin SWBB_3_sss
iSWBB_3_sss_thisrec = na.omit(data.frame(year = binvec, val = iSWBB_3_sss))
iSWBB_3_sss_bin = geoChronR::bin(iSWBB_3_sss_thisrec$year, 
                                 iSWBB_3_sss_thisrec$val, time_bin)
sss_binned[,3] = iSWBB_3_sss_bin$y - mean(iSWBB_3_sss_bin$y)

#### PSMs ####
# Time dimension
start = 1851
end = 1980
time = seq(from = start, to = end, by = 1)
a = -0.22
b = 0.145503
SWBB_ensemble = matrix(data = NA, nrow = 130, ncol = 3)
for (i in 1:3){
  for (j in 1:length(time)){
    SWBB_ensemble[j, i] = a*sst_binned[j, i] + b*sss_binned[j, i]
  }
}
SWBB_ensemble_mean = rowMeans(SWBB_ensemble)
#### Plot PSM output ####
for(i in 1:3){
  plot(time, SWBB_ensemble[,i], type = "l", frame.plot = F, ylim = c(-0.4,0.4),
       xlab = "", ylab = "", axes = F)
  par(new = T)
}
plot(time, SWBB_ensemble_mean, type = "l", frame.plot = F, ylim = c(-0.4,0.4),
     xlab = "", ylab = "", axes = F, lwd = 3)
par(new = T)
plot(NA, frame.plot = F, ylim = c(-0.4,0.4), xlim = c(min(time), max(time)),
     ylab = "Pseudo d18O anomalies", xlab = "Years (CE)", axes = T, 
     main = "SWBB Pseudocorals")
