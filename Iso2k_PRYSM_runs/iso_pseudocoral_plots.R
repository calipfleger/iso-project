library(geoChronR)

#### Import iCESM and observational data ####
setwd("~/GitHub/iso-project/Iso2k_PRYSM_runs/data/iCESM_full_forcing_afedits")
sst = read.csv("ssta_BOFPiM_binned.csv", head = T)
sss = read.csv("sssa_BOFPiM_binned.csv", head = T)
#sst_obs = read.csv("binned_BermudaSST.csv", head = T)

# Import PRYSM pseudocoral result
setwd("~/GitHub/iso-project/Iso2k_PRYSM_runs/results")
coral = read.csv("BOFPiM_pseudocoral.csv", head = F)
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
                                      variable = "dataSetName") == "CO98BOFP")]
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
     ylab = "iCESM SST (C)", xpd = NA, main = bofp, 
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
(cor = cor.test(coral, binned_proxy$y))
plot(coral, binned_proxy$y, pch = 16, xlab = "Pseudo-coral",
     ylab = "Proxy Data", main = c("r = -0.10", "p = 0.25"))

#
##############################################################################
#
#### Plot modeled SST and SSS data ####
par(mfrow = c(2, 1), omi = c(0.7, 0.4, 0, 0), mai = c(0.2, 0.5, 0.2, 0))
plot(sst$x, sst_c, type = "l", frame.plot = F, xaxt = "n", xlab = "", lwd = 2,
     ylab = "SST (C)", xpd = NA, main = "KUBE1 Lat:32.5 Lon:-64.7")
plot(sss$x, sss$y, type = "l", frame.plot = F, xlab = "Year (CE)",
     ylab = "SSS (PSU)", xpd = NA, lwd = 2)

#
#### Plot pseudo and binned proxy ####
par(mfrow = c(2, 1), omi = c(0.7, 0.4, 0, 0), mai = c(0.2, 0.5, 0.2, 0))
# Pseudo
plot(time, coral, type = "l", frame.plot = F, ylab = "KUBE2 Pseudo-coral",
     xaxt = "n", xlab = "", xpd = NA, lwd = 2,
     main = "KUBE1 Lat:32.5 Lon:-64.7")
# Proxy
plot(binned_proxy$x, binned_proxy$y, type = "l", frame.plot = F,
     ylab = "Proxy data (permil)", xlab = "CE", xpd = NA,
     lwd = 2)

#### Plot modeled and observed SST ####
par(mfrow = c(2, 1), omi = c(0.7, 0.4, 0, 0), mai = c(0.2, 0.5, 0.2, 0))
plot(sst$x, sst_c, type = "l", frame.plot = F, xaxt = "n", xlab = "", lwd = 2,
     ylab = "iCESM Temp (C)", xpd = NA, ylim = c(22, 24),
     main = "Modeled vs. Observed Temp Lat:32.5 Lon:-64.7")
plot(sst_obs$x, sst_obs$y, type = "l", frame.plot = F, xlab = "Year (CE)",
     ylab = "COBE-SST2 Temp (C)", xpd = NA, lwd = 2, ylim = c(22, 24))


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

cor = cor.test(coral_obs, binned_proxy$y)
plot(coral_obs, binned_proxy$y, pch = 16, xlab = "Observed SST pseudo-coral",
     ylab = "Proxy Data", main = c("r = 0.133", "p = 0.13"))

#### Map proxy locations ####
# Define archive point types
library(magrittr)
library(dplyr)
library(ggplot2)
library(ggtext)
data("wrld_simpl", package = "maptools")

shapes = c("Coral" = 15, "GlacierIce" = 20, "Sclerosponge" = 19, "Wood" = 18,
           "LakeSediment" = 17, "Speleothem" = 17, "MolluskShells" = 19,
           "MarineSediment" = 8)

#proxy = c("KUBE" = 15, "KUBE" = 20, "KIPR" = 18)
proxy = rep(15, 6)

setwd("~/GitHub/iso-project/Iso2k_PRYSM_runs/data")
coral_loc = read.csv("coral_proxy_locations.csv", head = T)

ggplot() +
        geom_polygon(data = wrld_simpl, aes(x = long, y = lat, group = group), 
                     fill = "grey", colour = "black", alpha = 0.2) +
        geom_point(data = coral_loc, mapping = aes(x = lon, y = lat, 
                                                    shape = archive), color = "blue", size = 7) +
        scale_shape_manual(values = proxy) + 
        geom_text(data = coral_loc, aes(x = lon, y = lat, label = archive),
                  hjust = -0.3, cex = 4) +
        # Removes Axes and labels
        scale_x_continuous(breaks = NULL) +
        scale_y_continuous(breaks = NULL) +
        xlab("") + 
        ylab("") +
        coord_equal(xlim = c(-180, 0), ylim = c(-90, 90)) +
        # Change theme to remove axes and ticks
        theme(panel.background = element_blank(),
              axis.ticks=element_blank(),
              panel.border = element_blank(),
              panel.grid.major = element_blank())
#