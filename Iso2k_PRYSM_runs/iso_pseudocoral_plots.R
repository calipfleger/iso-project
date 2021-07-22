# Import PRYSM pseudocoral result
setwd("~/GitHub/iso-project/Iso2k_PRYSM_runs/results")
coral = read.csv("pseudocoral.csv", head = F)
coral = as.vector(unlist(coral))

# Time dimension
start = 1851
end = 1980
time = seq(from = start, to = end, by = 1)

# Plot pseudocoral TS
plot(time, coral, type = "l", frame.plot = F, ylab = "pseudo-coral")

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

