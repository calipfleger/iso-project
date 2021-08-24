library(geoChronR)

#### Import iCESM and observational data ####
setwd("~/GitHub/iso-project/Iso2k_PRYSM_runs/data/iCESM_full_forcing_afedits")
sst = read.csv("ssta_COPMiM_binned.csv", head = T)
sss = read.csv("sssa_COPMiM_binned.csv", head = T)
#sst_obs = read.csv("binned_BermudaSST.csv", head = T)

# Import PRYSM pseudocoral result
setwd("~/GitHub/iso-project/Iso2k_PRYSM_runs/results")
coral = read.csv("COPMiM_pseudocoral.csv", head = F)
coral = as.vector(unlist(coral))
#coral_obs = read.csv("pseudocoral_obs.csv", head = F)
#coral_obs = as.vector(unlist(coral_obs))

# Time dimension
start = 1851
end = 1980
time = seq(from = start, to = end, by = 1)

# Plot pseudocoral TS
#plot(time, coral, type = "l", frame.plot = F, ylab = "pseudo-coral")

# Convert model temps
sst_c = sst$y - 273.15

#### Load and extract Iso2k ####
setwd("C:/Users/andre/Box Sync/Konecky Lab/User Storage/Andrew Flaim/Rcode/Iso2k") # Laptop
load("iso2k1_0_0.RData")
# Remove extraneous objects
rm(D, TS)
# Apply geoChronR function to extract primary TS
bermuda = sTS[which(pullTsVariable(sTS, 
                                      variable = "dataSetName") == "CO03COPM")]
bermuda = bermuda[[1]] # This is the d18O subset

# Plot raw proxy data
plot(bermuda$year, bermuda$paleoData_values, type = "l",
     ylab = bermuda$paleoData_units, xlab = bermuda$yearUnits,
     main = bermuda$geo_siteName, xlim = c(1880, 1950), frame.plot = F)

#
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

#### Full data stack ####
# Create location name variables
kube = "KUBEiM Lat:32.467, Lon:-64.7"
swbb = "SWBBiM Lat:25.3903, lon:-80.1715"
nump = "NUMPiM Lat:5.87 Lon:-162.13"
copm = "COPMiM Lat:5.87, Lon:-162.13"
bofp = "BOFPiM Lat:-17.5, Lon:-149.83"
lira = "LIRAiM Lat:-21.24, Lon:-159.83"

par(mfrow = c(4, 1), omi = c(0.5, 0.3, 0, 0), mai = c(0.2, 0.5, 0.2, 0))
# SST
plot(sst$x, sst_c, type = "l", frame.plot = F, xaxt = "n", xlab = "", lwd = 2,
     ylab = "iCESM SST (C)", xpd = NA, main = copm, 
     cex.lab = 1.5, cex.axis = 1.5)
# SSS
plot(sss$x, sss$y, type = "l", frame.plot = F, xlab = "", xaxt = "n",
     ylab = "iCESM SSS (PSU)", xpd = NA, lwd = 2,
     cex.lab = 1.5, cex.axis = 1.5, main = "")
# Modeled
plot(time, coral, type = "l", frame.plot = F, ylab = "Pseudo-coral (modeled)",
     xaxt = "n", xlab = "", xpd = NA, lwd = 2,
     main = "", cex.lab = 1.5, cex.axis = 1.5)
# Proxy
plot(binned_proxy$x, binned_proxy$y, type = "l", frame.plot = F,
     ylab = "Proxy data (permil)", xlab = "Year (CE)", xpd = NA,
     lwd = 2, cex.lab = 1.5, cex.axis = 1.5)

#### Plot pseudo vs. proxy correlation ####
#(cor = cor.test(coral, binned_proxy$y))
#plot(coral, binned_proxy$y, pch = 16, xlab = "Pseudo-coral",
#     ylab = "Proxy Data", main = c("r = -0.10", "p = 0.25"))

#