library(geoChronR)
library(magrittr)
library(tidyverse)
library(lubridate)
library(rgdal)
library(ncdf4)
library(ncdf4.helpers)

# Time dimension
start = 1851
end = 1981
time = seq(from = start, to = end, by = 1)
#
#### Import iCESM data ####
nrec = 7
setwd("~/GitHub/iso-project/full_forcing_files/icesm_fullforcing_indv")

iBOFP_1 = read.csv("iBOFP_1.csv", head = T)
iBOFP_2 = read.csv("iBOFP_2.csv", head = T)
iBOFP_3 = read.csv("iBOFP_3.csv", head = T)

iBOFP_1 = iBOFP_1[12:1572,]
iBOFP_2 = iBOFP_2[12:1572,]
iBOFP_3 = iBOFP_3[12:1572,]

iBOFP_1_sss = iBOFP_1[,2]
iBOFP_2_sss = iBOFP_2[,2]
iBOFP_3_sss = iBOFP_3[,2]

iBOFP_1_sst = iBOFP_1[,3]
iBOFP_2_sst = iBOFP_2[,3]
iBOFP_3_sst = iBOFP_3[,3]

#### Bin sst and sss ####
# Binned sst anomaly matrix
sst_binned = matrix(data = NA, nrow = 130, ncol = 3)
# Bin BOFP_1_sst
binvec = seq(from = start, to = end, by = 1/12)
iBOFP_1_sst_thisrec = na.omit(data.frame(year = binvec, val = iBOFP_1_sst))
iBOFP_1_sst_bin = geoChronR::bin(iBOFP_1_sst_thisrec$year, 
                                 iBOFP_1_sst_thisrec$val, time)
sst_binned[,1] = iBOFP_1_sst_bin$y - mean(iBOFP_1_sst_bin$y)
# Bin BOFP_2_sst
iBOFP_2_sst_thisrec = na.omit(data.frame(year = binvec, val = iBOFP_2_sst))
iBOFP_2_sst_bin = geoChronR::bin(iBOFP_2_sst_thisrec$year, 
                                 iBOFP_2_sst_thisrec$val, time)
sst_binned[,2] = iBOFP_2_sst_bin$y - mean(iBOFP_2_sst_bin$y)
# Bin BOFP_1_sst
iBOFP_3_sst_thisrec = na.omit(data.frame(year = binvec, val = iBOFP_3_sst))
iBOFP_3_sst_bin = geoChronR::bin(iBOFP_3_sst_thisrec$year, 
                                 iBOFP_3_sst_thisrec$val, time)
sst_binned[,3] = iBOFP_3_sst_bin$y - mean(iBOFP_3_sst_bin$y)

# Binned sss anomaly matrix
sss_binned = matrix(data = NA, nrow = 130, ncol = 3)
# Bin BOFP_1_sss
binvec = seq(from = start, to = end, by = 1/12)
iBOFP_1_sss_thisrec = na.omit(data.frame(year = binvec, val = iBOFP_1_sss))
iBOFP_1_sss_bin = geoChronR::bin(iBOFP_1_sss_thisrec$year, 
                                 iBOFP_1_sss_thisrec$val, time)
sss_binned[,1] = iBOFP_1_sss_bin$y - mean(iBOFP_1_sss_bin$y)
# Bin BOFP_2_sss
iBOFP_2_sss_thisrec = na.omit(data.frame(year = binvec, val = iBOFP_2_sss))
iBOFP_2_sss_bin = geoChronR::bin(iBOFP_2_sss_thisrec$year, 
                                 iBOFP_2_sss_thisrec$val, time)
sss_binned[,2] = iBOFP_2_sss_bin$y - mean(iBOFP_2_sss_bin$y)
# Bin BOFP_3_sss
iBOFP_3_sss_thisrec = na.omit(data.frame(year = binvec, val = iBOFP_3_sss))
iBOFP_3_sss_bin = geoChronR::bin(iBOFP_3_sss_thisrec$year, 
                                 iBOFP_3_sss_thisrec$val, time)
sss_binned[,3] = iBOFP_3_sss_bin$y - mean(iBOFP_3_sss_bin$y)

#### PSMs ####
a = -0.22
b = 0.4234225
BOFP_ensemble = matrix(data = NA, nrow = 130, ncol = 3)
for (i in 1:3){
  for (j in 1:length(time)){
    BOFP_ensemble[j, i] = a*sst[j, i] + b*sss[j, i]
  }
}

#### Plot PSM output ####
for(i in 1:3){
  plot(time, BOFP_ensemble[,i], type = "l", frame.plot = F, ylim = c(-0.2,0.2),
       xlab = "", ylab = "", axes = F)
  par(new = T)
}
plot(NA, frame.plot = F, ylim = c(-0.2,0.2), xlim = c(min(time), max(time)),
     ylab = "Pseudo d18O anomalies", xlab = "Years (CE)", axes = T, 
     main = "BOFP Pseudocorals")
