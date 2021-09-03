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

iLIRA_1 = read.csv("iLIRA_1.csv", head = T)
iLIRA_2 = read.csv("iLIRA_2.csv", head = T)
iLIRA_3 = read.csv("iLIRA_3.csv", head = T)

iLIRA_1 = iLIRA_1[12:1572,]
iLIRA_2 = iLIRA_2[12:1572,]
iLIRA_3 = iLIRA_3[12:1572,]

iLIRA_1_sss = iLIRA_1[,2]
iLIRA_2_sss = iLIRA_2[,2]
iLIRA_3_sss = iLIRA_3[,2]

iLIRA_1_sst = iLIRA_1[,3]
iLIRA_2_sst = iLIRA_2[,3]
iLIRA_3_sst = iLIRA_3[,3]

#### Bin sst and sss ####
# Binned sst anomaly matrix
# Time dimension
start_bin = 1851
end_bin = 1981
time_bin = seq(from = start_bin, to = end_bin, by = 1)

sst_binned = matrix(data = NA, nrow = 130, ncol = 3)
# Bin LIRA_1_sst
binvec = seq(from = start_bin, to = end_bin, by = 1/12)
iLIRA_1_sst_thisrec = na.omit(data.frame(year = binvec, val = iLIRA_1_sst))
iLIRA_1_sst_bin = geoChronR::bin(iLIRA_1_sst_thisrec$year, 
                                 iLIRA_1_sst_thisrec$val, time_bin)
sst_binned[,1] = iLIRA_1_sst_bin$y - mean(iLIRA_1_sst_bin$y)
# Bin LIRA_2_sst
iLIRA_2_sst_thisrec = na.omit(data.frame(year = binvec, val = iLIRA_2_sst))
iLIRA_2_sst_bin = geoChronR::bin(iLIRA_2_sst_thisrec$year, 
                                 iLIRA_2_sst_thisrec$val, time_bin)
sst_binned[,2] = iLIRA_2_sst_bin$y - mean(iLIRA_2_sst_bin$y)
# Bin LIRA_1_sst
iLIRA_3_sst_thisrec = na.omit(data.frame(year = binvec, val = iLIRA_3_sst))
iLIRA_3_sst_bin = geoChronR::bin(iLIRA_3_sst_thisrec$year, 
                                 iLIRA_3_sst_thisrec$val, time_bin)
sst_binned[,3] = iLIRA_3_sst_bin$y - mean(iLIRA_3_sst_bin$y)

# Binned sss anomaly matrix
sss_binned = matrix(data = NA, nrow = 130, ncol = 3)
# Bin LIRA_1_sss
binvec = seq(from = start_bin, to = end_bin, by = 1/12)
iLIRA_1_sss_thisrec = na.omit(data.frame(year = binvec, val = iLIRA_1_sss))
iLIRA_1_sss_bin = geoChronR::bin(iLIRA_1_sss_thisrec$year, 
                                 iLIRA_1_sss_thisrec$val, time_bin)
sss_binned[,1] = iLIRA_1_sss_bin$y - mean(iLIRA_1_sss_bin$y)
# Bin LIRA_2_sss
iLIRA_2_sss_thisrec = na.omit(data.frame(year = binvec, val = iLIRA_2_sss))
iLIRA_2_sss_bin = geoChronR::bin(iLIRA_2_sss_thisrec$year, 
                                 iLIRA_2_sss_thisrec$val, time_bin)
sss_binned[,2] = iLIRA_2_sss_bin$y - mean(iLIRA_2_sss_bin$y)
# Bin LIRA_3_sss
iLIRA_3_sss_thisrec = na.omit(data.frame(year = binvec, val = iLIRA_3_sss))
iLIRA_3_sss_bin = geoChronR::bin(iLIRA_3_sss_thisrec$year, 
                                 iLIRA_3_sss_thisrec$val, time_bin)
sss_binned[,3] = iLIRA_3_sss_bin$y - mean(iLIRA_3_sss_bin$y)

#### PSMs ####
a = -0.22
b = 0.4234225
start = 1851
end = 1980
time = seq(from = start, to = end, by = 1)
LIRA_ensemble = matrix(data = NA, nrow = 130, ncol = 3)
for (i in 1:3){
  for (j in 1:length(time)){
    LIRA_ensemble[j, i] = a*sst_binned[j, i] + b*sss_binned[j, i]
  }
}
LIRA_ensemble_mean = rowMeans(LIRA_ensemble)
#### Plot PSM output ####
for(i in 1:3){
  plot(time, LIRA_ensemble[,i], type = "l", frame.plot = F, ylim = c(-0.4,0.4),
       xlab = "", ylab = "", axes = F, lwd = 0.5)
  par(new = T)
}
plot(time, LIRA_ensemble_mean, type = "l", frame.plot = F, ylim = c(-0.4,0.4),
     xlab = "", ylab = "", axes = F, lwd = 3)
par(new = T)
plot(NA, frame.plot = F, ylim = c(-0.4,0.4), xlim = c(min(time), max(time)),
     ylab = "Pseudo d18O anomalies", xlab = "Years (CE)", axes = T, 
     main = "LIRA Pseudocorals")
