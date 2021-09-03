#### Packages and dimension initiation ####
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

iCOPM_1 = read.csv("iCOPM_1.csv", head = T)
iCOPM_2 = read.csv("iCOPM_2.csv", head = T)
iCOPM_3 = read.csv("iCOPM_3.csv", head = T)

iCOPM_1 = iCOPM_1[12:1572,]
iCOPM_2 = iCOPM_2[12:1572,]
iCOPM_3 = iCOPM_3[12:1572,]

iCOPM_1_sss = iCOPM_1[,2]
iCOPM_2_sss = iCOPM_2[,2]
iCOPM_3_sss = iCOPM_3[,2]

iCOPM_1_sst = iCOPM_1[,3]
iCOPM_2_sst = iCOPM_2[,3]
iCOPM_3_sst = iCOPM_3[,3]

#### Bin sst and sss ####
# Binned sst anomaly matrix
# Time dimension
start_bin = 1851
end_bin = 1981
time_bin = seq(from = start_bin, to = end_bin, by = 1)

sst_binned = matrix(data = NA, nrow = 130, ncol = 3)
# Bin COPM_1_sst
binvec = seq(from = start_bin, to = end_bin, by = 1/12)
iCOPM_1_sst_thisrec = na.omit(data.frame(year = binvec, val = iCOPM_1_sst))
iCOPM_1_sst_bin = geoChronR::bin(iCOPM_1_sst_thisrec$year, 
                                 iCOPM_1_sst_thisrec$val, time_bin)
sst_binned[,1] = iCOPM_1_sst_bin$y - mean(iCOPM_1_sst_bin$y)
# Bin COPM_2_sst
iCOPM_2_sst_thisrec = na.omit(data.frame(year = binvec, val = iCOPM_2_sst))
iCOPM_2_sst_bin = geoChronR::bin(iCOPM_2_sst_thisrec$year, 
                                 iCOPM_2_sst_thisrec$val, time_bin)
sst_binned[,2] = iCOPM_2_sst_bin$y - mean(iCOPM_2_sst_bin$y)
# Bin COPM_1_sst
iCOPM_3_sst_thisrec = na.omit(data.frame(year = binvec, val = iCOPM_3_sst))
iCOPM_3_sst_bin = geoChronR::bin(iCOPM_3_sst_thisrec$year, 
                                 iCOPM_3_sst_thisrec$val, time_bin)
sst_binned[,3] = iCOPM_3_sst_bin$y - mean(iCOPM_3_sst_bin$y)

# Binned sss anomaly matrix
sss_binned = matrix(data = NA, nrow = 130, ncol = 3)
# Bin COPM_1_sss
binvec = seq(from = start_bin, to = end_bin, by = 1/12)
iCOPM_1_sss_thisrec = na.omit(data.frame(year = binvec, val = iCOPM_1_sss))
iCOPM_1_sss_bin = geoChronR::bin(iCOPM_1_sss_thisrec$year, 
                                 iCOPM_1_sss_thisrec$val, time_bin)
sss_binned[,1] = iCOPM_1_sss_bin$y - mean(iCOPM_1_sss_bin$y)
# Bin COPM_2_sss
iCOPM_2_sss_thisrec = na.omit(data.frame(year = binvec, val = iCOPM_2_sss))
iCOPM_2_sss_bin = geoChronR::bin(iCOPM_2_sss_thisrec$year, 
                                 iCOPM_2_sss_thisrec$val, time_bin)
sss_binned[,2] = iCOPM_2_sss_bin$y - mean(iCOPM_2_sss_bin$y)
# Bin COPM_3_sss
iCOPM_3_sss_thisrec = na.omit(data.frame(year = binvec, val = iCOPM_3_sss))
iCOPM_3_sss_bin = geoChronR::bin(iCOPM_3_sss_thisrec$year, 
                                 iCOPM_3_sss_thisrec$val, time_bin)
sss_binned[,3] = iCOPM_3_sss_bin$y - mean(iCOPM_3_sss_bin$y)

#### PSMs ####
a = -0.22
b = 0.2540535
start = 1851
end = 1980
time = seq(from = start, to = end, by = 1)
COPM_ensemble = matrix(data = NA, nrow = 130, ncol = 3)
for (i in 1:3){
  for (j in 1:length(time)){
    COPM_ensemble[j, i] = a*sst_binned[j, i] + b*sss_binned[j, i]
  }
}
COPM_ensemble_mean = rowMeans(COPM_ensemble)
#### Plot PSM output ####
for(i in 1:3){
  plot(time, COPM_ensemble[,i], type = "l", frame.plot = F, ylim = c(-0.4,0.4),
       xlab = "", ylab = "", axes = F, lwd = 0.5)
  par(new = T)
}
plot(time, COPM_ensemble_mean, type = "l", frame.plot = F, ylim = c(-0.4,0.4),
     xlab = "", ylab = "", axes = F, lwd = 3)
par(new = T)
plot(NA, frame.plot = F, ylim = c(-0.4,0.4), xlim = c(min(time), max(time)),
     ylab = "Pseudo d18O anomalies", xlab = "Years (CE)", axes = T, 
     main = "COPM Pseudocorals")
