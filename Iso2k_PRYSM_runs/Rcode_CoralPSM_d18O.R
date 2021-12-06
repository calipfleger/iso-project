library(geoChronR)
library(magrittr)
library(tidyverse)
library(lubridate)
library(rgdal)
library(ncdf4)
library(ncdf4.helpers)

#
#### Import iCESM data ####
#setwd("D:/GitHub/iso-project/TS_1850")
setwd("~/GitHub/iso-project/TS_1850")               # Laptop
sst = read.csv("TS_full1_clipped.csv", head = T)
#setwd("D:/GitHub/iso-project/d18osw_coral")
setwd("~/GitHub/iso-project/d18osw_coral")          # Laptop
d18Osw = read.csv("test_full_post1.csv", head = T)

# Trim to uniform timespan
d18Osw = d18Osw[12:1872,]
sst = sst[13:1873,]
# Reformat time axis
old.time = d18Osw[,1]
d18Osw = d18Osw[,-1]
sst = sst[,-1]
start = 1851
end = 2006
monthly = seq(from = start, to = end, by = 1/12)

# Remove incorrect records
sst_trim = sst
remove = vector()
for (i in 1:ncol(sst)){
  #print(i)
  name = which(names(sst)[i] == names(d18Osw))
  if (!length(name)){
    remove = cbind(remove, i)}
}
sst_trim = sst[,-remove]
sst = sst_trim

# Convert sst to degrees C, not Kelvin.
temp_flag = any(sst>200)
for (i in 1:dim(sst)[2]){
  for (j in 1:length(monthly)){
    if (temp_flag){
      sst[j, i] = sst[j, i]-274.15
    }
  }
}

# Convert to matching column order
sst_rs = sst * 0
for (i in 1:ncol(sst)){
  col = which(names(sst) == names(d18Osw)[i])
  sst_rs[,i] = sst[,col]
}
names(sst_rs) = names(d18Osw)
sst = sst_rs

#### PSM ####
a = -0.22
coral = matrix(data = NA, nrow = length(monthly), ncol = ncol(sst))
for (i in 1:ncol(sst)){
  for (j in 1:length(monthly)){
    coral[j, i] = a*sst[j, i] + d18Osw[j, i]
  }
}
coral = data.frame(coral)
names(coral) = names(sst)

coral.out = cbind(old.time, coral)
names(coral.out)[1] = "time"
#
#### Bin PSM input and output to annual ####
binvec =  seq(from = 1850, to = 2006, by = 1)
binYears = rowMeans(cbind(binvec[-1], binvec[-length(binvec)]))
binned_coral = matrix(NA, length(binYears), length(coral))
binned_sst = matrix(NA, length(binYears), length(coral))
binned_d18Osw = matrix(NA, length(binYears), length(coral))

for(i in 1:(ncol(binned_coral))) {
  thisrec = na.omit(data.frame(year = monthly, val = coral[,i]))
  thisrec_sst = na.omit(data.frame(year = monthly, val = sst[,i]))
  thisrec_d18Osw = na.omit(data.frame(year = monthly, val = d18Osw[,i]))
  binned_coral[,i] = geoChronR::bin(thisrec$year, thisrec$val, binvec)$y
  binned_sst[,i] = geoChronR::bin(thisrec_sst$year, thisrec_sst$val, binvec)$y
  binned_d18Osw[,i] = geoChronR::bin(thisrec_d18Osw$year, thisrec_d18Osw$val, binvec)$y
}

#### Export PSM output ####
setwd("D:/GitHub/iso-project")
write.csv(coral.out, "pseudocoral_1850-2005.csv", row.names = F)

#
#### Import and bin Iso2k records ####
data("wrld_simpl", package = "maptools")
# Load and extract Iso2k
if (Sys.info()['sysname'] == "Darwin"){
  setwd("~/Box Sync/Konecky Lab/User Storage/Andrew Flaim/Rcode/Iso2k") # MacOS
} else if (Sys.info()['login'] == "andre"){
  setwd("C:/Users/andre/Box Sync/Konecky Lab/User Storage/Andrew Flaim/Rcode/Iso2k") # Laptop
} else{
  setwd("D:/Box/Box_Sync/Konecky Lab/User Storage/Andrew Flaim/Rcode/Iso2k") # Desktop
}

load("iso2k1_0_0.RData")
# Remove extraneous objects
rm(D, TS)
# Apply geoChronR function to extrat primary TS
all_iso_ts = sTS[which(pullTsVariable(sTS, 
                                      variable = "paleoData_iso2kPrimaryTimeseries") == "TRUE")]
iso2k_coral = all_iso_ts[which(pullTsVariable(all_iso_ts,
                                        variable = "archiveType") == "Coral")]

CO06LIFI = na.omit(data.frame(year = iso2k_coral[[24]]$year,
                              val = iso2k_coral[[24]]$paleoData_values))
binvec =  seq(from = 1850, to = 2006, by = 1)
binYears = rowMeans(cbind(binvec[-1], binvec[-length(binvec)]))
binned_rec = geoChronR::bin(CO06LIFI$year, CO06LIFI$val, binvec)$y

#
#### Plot inputs and pseudocorals ####
cex = 2
par(mfrow = c(4, 1), mai = c(0.5, 1,0.25,0.5))
plot(binYears, binned_sst[,1], type = "l", xlab = "", xaxt = "n", ylab = "",
     frame.plot = F, cex.axis = cex, xlim = c(1850, 2000))
mtext(side = 2, line = 3.5, "SST (deg. C)", cex = cex)
plot(binYears, binned_d18Osw[,1], type = "l", xlab = "", xaxt = "n", ylab = "",
     frame.plot = F, cex.axis = cex, xlim = c(1850, 2000))
mtext(side = 2, line = 3.5, "d18Osw (permil)", cex = cex)
plot(binYears, binned_coral[,1], type = "l", xlab = "", ylab = "",
     frame.plot = F, xaxt = "n", cex.axis = cex, xlim = c(1850, 2000))
mtext(side = 2, line = 3.5, "Pseudocoral d18O", cex = cex)
plot(iso2k_coral[[24]][["year"]], iso2k_coral[[24]][["paleoData_values"]], type = "l",
     frame.plot = F, xlim = c(1850, 2000), xlab = "", ylab = "", cex.axis = cex)
mtext(side = 1, line = 2.5, "Year (CE)", cex = cex)
mtext(side = 2, line = 3.5, "CO06LIFI d18O", cex = cex)
#

#### Stack all pseudocorals ####
for (i in 1:ncol(coral)){
  plot(monthly, coral[,i], type = "l", frame.plot = F, axes = F, xlab = "", 
       ylab = "", ylim = range(coral))
  par(new = T)
}
axis(1)
mtext(side = 1, line = 2, "Year (CE)")
axis(2)
mtext(side = 2, line = 2.4, "Pseudocoral d18O")

#
#### Remove seasonal cycle ####
require(graphics) 
coralts = ts(coral[,1], start = 1850, end = 2005, deltat = 1/12) 
stlcoral = stl(coralts, s.window = "periodic")
plot(stlcoral)
print(stlcoral$time.series)
plot(stlcoral$time.series[,2]+stlcoral$time.series[,3])
summary(stlcoral)
#
###############################################################################
#### Plot Pseudocorals ####
rec_col = c("blue", "red", "dark green", "purple", "black", "orange", "brown")
for(i in 1:ncol(sst)){
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