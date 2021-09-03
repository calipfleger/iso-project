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

iKIPR_1 = read.csv("iKIPR_1.csv", head = T)
iKIPR_2 = read.csv("iKIPR_2.csv", head = T)
iKIPR_3 = read.csv("iKIPR_3.csv", head = T)

iKIPR_1 = iKIPR_1[12:1572,]
iKIPR_2 = iKIPR_2[12:1572,]
iKIPR_3 = iKIPR_3[12:1572,]

iKIPR_1_sss = iKIPR_1[,2]
iKIPR_2_sss = iKIPR_2[,2]
iKIPR_3_sss = iKIPR_3[,2]

iKIPR_1_sst = iKIPR_1[,3]
iKIPR_2_sst = iKIPR_2[,3]
iKIPR_3_sst = iKIPR_3[,3]

#### Bin sst and sss ####
start_bin = 1851
end_bin = 1981
time_bin = seq(from = start_bin, to = end_bin, by = 1)
# Binned sst anomaly matrix
sst_binned = matrix(data = NA, nrow = 130, ncol = 3)
# Bin KIPR_1_sst
binvec = seq(from = start_bin, to = end_bin, by = 1/12)
iKIPR_1_sst_thisrec = na.omit(data.frame(year = binvec, val = iKIPR_1_sst))
iKIPR_1_sst_bin = geoChronR::bin(iKIPR_1_sst_thisrec$year, 
                                 iKIPR_1_sst_thisrec$val, time_bin)
sst_binned[,1] = iKIPR_1_sst_bin$y - mean(iKIPR_1_sst_bin$y)
# Bin KIPR_2_sst
iKIPR_2_sst_thisrec = na.omit(data.frame(year = binvec, val = iKIPR_2_sst))
iKIPR_2_sst_bin = geoChronR::bin(iKIPR_2_sst_thisrec$year, 
                                 iKIPR_2_sst_thisrec$val, time_bin)
sst_binned[,2] = iKIPR_2_sst_bin$y - mean(iKIPR_2_sst_bin$y)
# Bin KIPR_1_sst
iKIPR_3_sst_thisrec = na.omit(data.frame(year = binvec, val = iKIPR_3_sst))
iKIPR_3_sst_bin = geoChronR::bin(iKIPR_3_sst_thisrec$year, 
                                 iKIPR_3_sst_thisrec$val, time_bin)
sst_binned[,3] = iKIPR_3_sst_bin$y - mean(iKIPR_3_sst_bin$y)

# Binned sss anomaly matrix
sss_binned = matrix(data = NA, nrow = 130, ncol = 3)
# Bin KIPR_1_sss
binvec = seq(from = start_bin, to = end_bin, by = 1/12)
iKIPR_1_sss_thisrec = na.omit(data.frame(year = binvec, val = iKIPR_1_sss))
iKIPR_1_sss_bin = geoChronR::bin(iKIPR_1_sss_thisrec$year, 
                                 iKIPR_1_sss_thisrec$val, time_bin)
sss_binned[,1] = iKIPR_1_sss_bin$y - mean(iKIPR_1_sss_bin$y)
# Bin KIPR_2_sss
iKIPR_2_sss_thisrec = na.omit(data.frame(year = binvec, val = iKIPR_2_sss))
iKIPR_2_sss_bin = geoChronR::bin(iKIPR_2_sss_thisrec$year, 
                                 iKIPR_2_sss_thisrec$val, time_bin)
sss_binned[,2] = iKIPR_2_sss_bin$y - mean(iKIPR_2_sss_bin$y)
# Bin KIPR_3_sss
iKIPR_3_sss_thisrec = na.omit(data.frame(year = binvec, val = iKIPR_3_sss))
iKIPR_3_sss_bin = geoChronR::bin(iKIPR_3_sss_thisrec$year, 
                                 iKIPR_3_sss_thisrec$val, time_bin)
sss_binned[,3] = iKIPR_3_sss_bin$y - mean(iKIPR_3_sss_bin$y)

#### PSMs ####
# Time dimension
start = 1851
end = 1980
time = seq(from = start, to = end, by = 1)
a = -0.22
b = 0.145503
KIPR_ensemble = matrix(data = NA, nrow = 130, ncol = 3)
for (i in 1:3){
  for (j in 1:length(time)){
    KIPR_ensemble[j, i] = a*sst_binned[j, i] + b*sss_binned[j, i]
  }
}
KIPR_ensemble_mean = rowMeans(KIPR_ensemble)
#### Plot PSM output ####
for(i in 1:3){
  plot(time, KIPR_ensemble[,i], type = "l", frame.plot = F, ylim = c(-0.4,0.4),
       xlab = "", ylab = "", axes = F)
  par(new = T)
}
plot(time, KIPR_ensemble_mean, type = "l", frame.plot = F, ylim = c(-0.4,0.4),
     xlab = "", ylab = "", axes = F, lwd = 3)
par(new = T)
plot(NA, frame.plot = F, ylim = c(-0.4,0.4), xlim = c(min(time), max(time)),
     ylab = "Pseudo d18O anomalies", xlab = "Years (CE)", axes = T, 
     main = "KIPR Pseudocorals")
