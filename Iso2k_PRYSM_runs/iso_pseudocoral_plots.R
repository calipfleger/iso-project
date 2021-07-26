library(geoChronR)

# Import iCESM and observational data
setwd("~/GitHub/iso-project/Iso2k_PRYSM_runs/data")
sst = read.csv("ssta_binned.csv", head = T)
sss = read.csv("sssa_binned.csv", head = T)
sst_obs = read.csv("binned_BermudaSST.csv", head = T)

# Import PRYSM pseudocoral result
setwd("~/GitHub/iso-project/Iso2k_PRYSM_runs/results")
coral = read.csv("pseudocoral.csv", head = F)
coral = as.vector(unlist(coral))
coral_obs = read.csv("pseudocoral_obs.csv", head = F)
coral_obs = as.vector(unlist(coral_obs))

# Time dimension
start = 1851
end = 1980
time = seq(from = start, to = end, by = 1)

# Plot pseudocoral TS
plot(time, coral, type = "l", frame.plot = F, ylab = "pseudo-coral")

# Convert model temps
sst_c = sst$y - 273.15

#### Load and extract Iso2k ####
setwd("C:/Users/andre/Box Sync/Konecky Lab/User Storage/Andrew Flaim/Rcode/Iso2k") # Laptop
load("iso2k1_0_0.RData")
# Remove extraneous objects
rm(D, TS)
# Apply geoChronR function to extract primary TS
bermuda = sTS[which(pullTsVariable(sTS, 
                                      variable = "dataSetName") == "CO05KUBE")]
bermuda = bermuda[[1]] # This is the d18O subset

# Plot raw proxy data
plot(bermuda$year, bermuda$paleoData_values, type = "l",
     ylab = bermuda$paleoData_units, xlab = bermuda$yearUnits,
     main = bermuda$geo_siteName)

#### Bin proxy data ####
start = 1851
end = 1981
#time = seq(from = start, to = end, by = 1/12)
binvec =  seq(from = start, to = end, by = 1)
thisrec = na.omit(data.frame(year = bermuda$year, 
                             val = bermuda$paleoData_values))
binned_proxy = geoChronR::bin(thisrec$year, thisrec$val, binvec)

# Plot binned proxy data
plot(binned_proxy$x, binned_proxy$y, type = "l", frame.plot = F,
     ylab = bermuda$paleoData_units, xlab = bermuda$yearUnits,
     main = bermuda$geo_siteName)

#### Plot modeled SST and SSS data ####
par(mfrow = c(2, 1), omi = c(0.7, 0.4, 0, 0), mai = c(0.2, 0.5, 0.2, 0))
plot(sst$x, sst_c, type = "l", frame.plot = F, xaxt = "n", xlab = "", lwd = 2,
     ylab = "Temp (C)", xpd = NA, main = "iCESM Full Forcing Lat:32.5 Lon:-64.7")
plot(sss$x, sss$y, type = "l", frame.plot = F, xlab = "Year (CE)",
     ylab = "Salinity (PSU)", xpd = NA, lwd = 2)

#### Plot modeled and observed SST ####
par(mfrow = c(2, 1), omi = c(0.7, 0.4, 0, 0), mai = c(0.2, 0.5, 0.2, 0))
plot(sst$x, sst_c, type = "l", frame.plot = F, xaxt = "n", xlab = "", lwd = 2,
     ylab = "iCESM Temp (C)", xpd = NA, ylim = c(22, 24),
     main = "Modeled vs. Observed Temp Lat:32.5 Lon:-64.7")
plot(sst_obs$x, sst_obs$y, type = "l", frame.plot = F, xlab = "Year (CE)",
     ylab = "COBE-SST2 Temp (C)", xpd = NA, lwd = 2, ylim = c(22, 24))

#### Plot pseudo and binned proxy ####
par(mfrow = c(2, 1), omi = c(0.7, 0.4, 0, 0), mai = c(0.2, 0.5, 0.2, 0))
# Pseudo
plot(time, coral, type = "l", frame.plot = F, ylab = "Pseudo-coral",
     xaxt = "n", xlab = "", xpd = NA, lwd = 2,
     main = bermuda$geo_siteName)
# Proxy
plot(binned_proxy$x, binned_proxy$y, type = "l", frame.plot = F,
     ylab = "Proxy data (permil)", xlab = bermuda$yearUnits, xpd = NA,
     lwd = 2)

#### Plot correlations ####
cor = cor.test(coral, binned_proxy$y)
plot(coral, binned_proxy$y, pch = 16, xlab = "Pseudo-coral",
     ylab = "Proxy Data", main = c("r = 0.033", "p = 0.71"))

#### Plot modeled and obs psuedo ####
par(mfrow = c(3, 1), omi = c(0.7, 0.4, 0, 0), mai = c(0.2, 0.5, 0.2, 0))
# modeled
plot(time, coral, type = "l", frame.plot = F, ylab = "Pseudo-coral (modeled)",
     xaxt = "n", xlab = "", xpd = NA, lwd = 2,
     main = "Modeled vs. Observed SST Pseudo-corals")
# observed
plot(time, coral_obs, type = "l", frame.plot = F, xaxt = "n",
     ylab = "Pseudo-coral (obs)", xlab = "", xpd = NA,
     lwd = 2)
# Proxy
plot(binned_proxy$x, binned_proxy$y, type = "l", frame.plot = F,
     ylab = "Proxy data (permil)", xlab = "Year (CE)", xpd = NA,
     lwd = 2)

#### Plot modeled vs observed correlations ####
cor = cor.test(coral, coral_obs)
plot(coral, coral_obs, pch = 16, xlab = "Modeled SST pseudo-coral",
     ylab = "Observed SST pseudo-coral", main = c("r = 0.133", "p = 0.13"))
