library(geoChronR)
library(magrittr)
library(tidyverse)
library(lubridate)
library(rgdal)
rec_col = c("blue", "red", "dark green", "purple", "black", "orange", "brown")

#
#### Import iCESM data ####
nrec = 7
setwd("~/GitHub/iso-project/Iso2k_PRYSM_runs/data/iCESM_full_forcing_afedits")
# SST
sst = matrix(data = NA, nrow = 130, ncol = nrec) %>%
        set_colnames(c("KUBE", "SWBB", "KIPR", "LIRA", "NUMP", "COPM", "BOFP"))
sst[,1] = read.csv("ssta_KUBEiM_binned.csv", head = T)[,3]
sst[,2] = read.csv("ssta_SWBBiM_binned.csv", head = T)[,3]
sst[,3] = read.csv("ssta_KIPRiM_binned.csv", head = T)[,3]
sst[,4] = read.csv("ssta_LIRAiM_binned.csv", head = T)[,3]
sst[,5] = read.csv("ssta_NUMPiM_binned.csv", head = T)[,3]
sst[,6] = read.csv("ssta_COPMiM_binned.csv", head = T)[,3]
sst[,7] = read.csv("ssta_BOFPiM_binned.csv", head = T)[,3]
# SSS
sss = matrix(data = NA, nrow = 130, ncol = nrec) %>%
        set_colnames(c("KUBE", "SWBB", "KIPR", "LIRA", "NUMP", "COPM", "BOFP"))
sss[,1] = read.csv("sssa_KUBEiM_binned.csv", head = T)[,3]
sss[,2] = read.csv("sssa_SWBBiM_binned.csv", head = T)[,3]
sss[,3] = read.csv("sssa_KIPRiM_binned.csv", head = T)[,3]
sss[,4] = read.csv("sssa_LIRAiM_binned.csv", head = T)[,3]
sss[,5] = read.csv("sssa_NUMPiM_binned.csv", head = T)[,3]
sss[,6] = read.csv("sssa_COPMiM_binned.csv", head = T)[,3]
sss[,7] = read.csv("sssa_BOFPiM_binned.csv", head = T)[,3]

#
#### Import PRYSM pseudocoral result ####
setwd("~/GitHub/iso-project/Iso2k_PRYSM_runs/results")
pseudo = matrix(data = NA, nrow = 130, ncol = nrec) %>%
        set_colnames(c("KUBE", "SWBB", "KIPR", "LIRA", "NUMP", "COPM", "BOFP"))
pseudo[,1] = as.vector(unlist(read.csv("KUBEiM_pseudocoral.csv", head = F)))
pseudo[,2] = as.vector(unlist(read.csv("SWBBiM_pseudocoral.csv", head = F)))
pseudo[,3] = as.vector(unlist(read.csv("KIPRiM_pseudocoral.csv", head = F)))
pseudo[,4] = as.vector(unlist(read.csv("LIRAiM_pseudocoral.csv", head = F)))
pseudo[,5] = as.vector(unlist(read.csv("NUMPiM_pseudocoral.csv", head = F)))
pseudo[,6] = as.vector(unlist(read.csv("COPMiM_pseudocoral.csv", head = F)))
pseudo[,7] = as.vector(unlist(read.csv("BOFPiM_pseudocoral.csv", head = F)))

#
#### Extract and bin Iso2k records####
# Time dimension
start = 1851
end = 1980
time = seq(from = start, to = end, by = 1)

# Convert model temps
#sst = sst - 273.15

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
SWBB_rec = na.omit(data.frame(year = SWBB_data$year, 
                              val = SWBB_data$paleoData_values))
binned_SWBB = geoChronR::bin(SWBB_rec$year, SWBB_rec$val, binvec)
SWBB_anom = binned_SWBB$y-mean(na.omit(binned_SWBB$y))

# KUBE
KUBE_data = sTS[which(pullTsVariable(sTS, 
                                     variable = "dataSetName") == "CO05KUBE")][[1]]
KUBE_rec = na.omit(data.frame(year = KUBE_data$year, 
                              val = KUBE_data$paleoData_values))
binned_KUBE = geoChronR::bin(KUBE_rec$year, KUBE_rec$val, binvec)
KUBE_anom = binned_KUBE$y-mean(na.omit(binned_KUBE$y))

# LIRA
LIRA_data = sTS[which(pullTsVariable(sTS, 
                                     variable = "dataSetName") == "CO04LIRA")][[1]]
LIRA_rec = na.omit(data.frame(year = LIRA_data$year, 
                              val = LIRA_data$paleoData_values))
binned_LIRA = geoChronR::bin(LIRA_rec$year, LIRA_rec$val, binvec)
LIRA_anom = binned_LIRA$y-mean(na.omit(binned_LIRA$y))

# NUMP
NUMP_data = sTS[which(pullTsVariable(sTS, 
                                     variable = "dataSetName") == "CO11NUPM")][[5]]
NUMP_rec = na.omit(data.frame(year = NUMP_data$year, 
                              val = NUMP_data$paleoData_values))
binned_NUMP = geoChronR::bin(NUMP_rec$year, NUMP_rec$val, binvec)
NUMP_anom = binned_NUMP$y-mean(na.omit(binned_NUMP$y))

# COPM
COPM_data = sTS[which(pullTsVariable(sTS, 
                                     variable = "dataSetName") == "CO03COPM")][[1]]
COPM_rec = na.omit(data.frame(year = COPM_data$year, 
                              val = COPM_data$paleoData_values))
binned_COPM = geoChronR::bin(COPM_rec$year, COPM_rec$val, binvec)
COPM_anom = binned_COPM$y-mean(na.omit(binned_COPM$y))

# BOFP
BOFP_data = sTS[which(pullTsVariable(sTS, 
                                     variable = "dataSetName") == "CO98BOFP")][[1]]
BOFP_rec = na.omit(data.frame(year = BOFP_data$year, 
                              val = BOFP_data$paleoData_values))
binned_BOFP = geoChronR::bin(BOFP_rec$year, BOFP_rec$val, binvec)
BOFP_anom = binned_BOFP$y-mean(na.omit(binned_BOFP$y))

# KIPR
KIPR_data = sTS[which(pullTsVariable(sTS, 
                                     variable = "dataSetName") == "CO08KIPR")][[1]]
KIPR_rec = na.omit(data.frame(year = KIPR_data$year, 
                              val = KIPR_data$paleoData_values))
binned_KIPR = geoChronR::bin(KIPR_rec$year, KIPR_rec$val, binvec)
KIPR_anom = binned_KIPR$y-mean(na.omit(binned_KIPR$y))

#
#### Plot iCESM output ####
# Plot SST
for(i in 1:nrec){
        plot(time, sst[,i], type = "l", frame.plot = F, ylim = c(-1,1),
             xlab = "", ylab = "", axes = F, col = rec_col[i])
        par(new = T)
}
plot(NA, frame.plot = F, ylim = c(-1,1), xlim = c(min(time), max(time)),
     ylab = "SST anomalies", xlab = "Years (CE)", axes = T, main = "iCESM full forcing SSTs")
legend("topright", col = rec_col, 
       legend = c("KUBE", "SWBB", "KIPR", "LIRA", "NUMP", "COPM", "BOFP"),
       lty = rep(1, nrec))

# Plot SSS
for(i in 1:nrec){
        plot(time, sss[,i], type = "l", frame.plot = F, ylim = c(-0.5,0.5),
             xlab = "", ylab = "", axes = F, col = rec_col[i])
        par(new = T)
}
plot(NA, frame.plot = F, ylim = c(-0.5,0.5), xlim = c(min(time), max(time)),
     ylab = "SSS anomalies", xlab = "Years (CE)", axes = T, main = "iCESM full forcing SSS")

#### Plot pseudoproxies ####
# Plot pseudoproxies
for(i in 1:nrec){
        plot(time, pseudo[,i], type = "l", frame.plot = F, ylim = c(-0.2,0.2),
             xlab = "", ylab = "", axes = F, col = rec_col[i])
        par(new = T)
}
plot(NA, frame.plot = F, ylim = c(-0.2,0.2), xlim = c(min(time), max(time)),
     ylab = "Pseudo d18O anomalies", xlab = "Years (CE)", axes = T, main = "Pseudocorals")
#
#### Plot Iso2k records ####
plot(binned_KUBE, type = "l", frame.plot = F, xlim = c(min(time), max(time)),
     ylab = "", xlab = "", axes = F, main = "", col = rec_col[1],
     ylim = c(-5.5, -3))
par(new = T)
plot(binned_LIRA, type = "l", frame.plot = F, xlim = c(min(time), max(time)),
     ylab = "", xlab = "", axes = F, main = "", col = rec_col[4],
     ylim = c(-5.5, -3))
par(new = T)
plot(binned_NUMP, type = "l", frame.plot = F, xlim = c(min(time), max(time)),
     ylab = "", xlab = "", axes = F, main = "", col = rec_col[5],
     ylim = c(-5.5, -3))
par(new = T)
plot(binned_SWBB, type = "l", frame.plot = F, xlim = c(min(time), max(time)),
     ylab = "", xlab = "", axes = F, main = "", col = rec_col[2],
     ylim = c(-5.5, -3))
par(new = T)
plot(binned_COPM, type = "l", frame.plot = F, xlim = c(min(time), max(time)),
     ylab = "", xlab = "", axes = F, main = "",col = rec_col[6],
     ylim = c(-5.5, -3))
par(new = T)
plot(binned_BOFP, type = "l", frame.plot = F, xlim = c(min(time), max(time)),
     ylab = "", xlab = "", axes = F, main = "",col = rec_col[7],
     ylim = c(-5.5, -3))
par(new = T)
plot(binned_KIPR, type = "l", frame.plot = F, xlim = c(min(time), max(time)),
     ylab = "d18O", xlab = "Years (CE)", axes = T, main = "Iso2k corals",
     col = rec_col[3], ylim = c(-5.5, -3))

#### Plot Iso2k anomalies ####
plot(time, KUBE_anom, type = "l", frame.plot = F, xlim = c(min(time), max(time)),
     ylab = "", xlab = "", axes = T, main = "", col = rec_col[1],
     ylim = c(-0.8, 0.8))
par(new = T)
plot(time, LIRA_anom, type = "l", frame.plot = F, xlim = c(min(time), max(time)),
     ylab = "", xlab = "", axes = F, main = "", col = rec_col[4],
     ylim = c(-0.8, 0.8))
par(new = T)
plot(time, NUMP_anom, type = "l", frame.plot = F, xlim = c(min(time), max(time)),
     ylab = "", xlab = "", axes = F, main = "", col = rec_col[5],
     ylim = c(-0.8, 0.8))
par(new = T)
plot(time, SWBB_anom, type = "l", frame.plot = F, xlim = c(min(time), max(time)),
     ylab = "", xlab = "", axes = F, main = "", col = rec_col[2],
     ylim = c(-0.8, 0.8))
par(new = T)
plot(time, COPM_anom, type = "l", frame.plot = F, xlim = c(min(time), max(time)),
     ylab = "", xlab = "", axes = F, main = "",col = rec_col[6],
     ylim = c(-0.8, 0.8))
par(new = T)
plot(time, BOFP_anom, type = "l", frame.plot = F, xlim = c(min(time), max(time)),
     ylab = "", xlab = "", axes = F, main = "",col = rec_col[7],
     ylim = c(-0.8, 0.8))
par(new = T)
plot(time, KIPR_anom, type = "l", frame.plot = F, xlim = c(min(time), max(time)),
     ylab = "d18O anomalies", xlab = "Years (CE)", axes = T, main = "Iso2k coral anomalies",
     col = rec_col[3], ylim = c(-0.8, 0.8))
#
#### Plot Atlantic Output and PSMS ####
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
        plot(time, pseudo[,i], type = "l", frame.plot = F, ylim = c(-0.05,0.05),
             xlab = "", ylab = "", axes = F, col = rec_col[i])
        par(new = T)
}
plot(NA, frame.plot = F, ylim = c(-0.05,0.05), xlim = c(min(time), max(time)),
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
#### Plot Pacific Output and PSMS ####
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
        plot(time, pseudo[,i], type = "l", frame.plot = F, ylim = c(-0.2,0.2),
             xlab = "", ylab = "", axes = F, col = rec_col[i])
        par(new = T)
}
plot(NA, frame.plot = F, ylim = c(-0.2,0.2), xlim = c(min(time), max(time)),
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
