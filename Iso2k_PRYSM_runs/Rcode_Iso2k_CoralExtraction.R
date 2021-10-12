#### FULL CIRCUM-ATLANTIC SUBSET ####
library(magrittr)
library(tidyverse)
library(geoChronR)
library(lubridate)
library(rgdal)


#### Extracting Atlantic Annual Series ####
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
coral = all_iso_ts[which(pullTsVariable(all_iso_ts,
                                       variable = "archiveType") == "Coral")]

test = all_iso_ts[which(pullTsVariable(all_iso_ts,
                                               variable = "archiveType") == "Coral")]
allRecs = matrix(NA, length(coral), 14) %>%
  set_colnames(c("record", "resolution", "duration", 
                 "archive", "lat", "lon", "infMat", "var", "interp", 
                 "season", "start_year", "end_year", "site", "pub"))

for(i in 1:(length(coral))) {
  thisYearVec = na.omit(coral[[i]]$year)
    allRecs[i, 1] = coral[[i]]$dataSetName
    allRecs[i, 2] = (max(thisYearVec)-min(thisYearVec))/length(thisYearVec)
    allRecs[i, 3] = max(thisYearVec)-min(thisYearVec)
    allRecs[i, 4] = coral[[i]]$archiveType
    allRecs[i, 5] = coral[[i]]$geo_latitude
    allRecs[i, 6] = coral[[i]]$geo_longitude
    allRecs[i, 7] = coral[[i]]$paleoData_inferredMaterial
    allRecs[i, 8] = coral[[i]]$paleoData_variableName
    allRecs[i, 9] = coral[[i]]$isotopeInterpretation1_variableGroup
    allRecs[i, 10] = if(!is.null(coral[[i]]$isotopeInterpretation1_seasonality)){
      coral[[i]]$isotopeInterpretation1_seasonality} else {NA}
    allRecs[i, 11] = min(thisYearVec)
    allRecs[i, 12] = max(thisYearVec)
    allRecs[i, 13] = coral[[i]]$geo_siteName
    allRecs[i, 14] = if(!is.null(coral[[i]]$pub1_citation)){
      coral[[i]]$pub1_citation} else {NA}
}

allRecs = as.data.frame(allRecs) %>%
  mutate_at(c("resolution", "lat", "lon", "duration"), as.numeric)

#
#### Plot coral locations ####
# Import completed location list
setwd("C:/Users/andre/Box Sync/Konecky Lab/User Storage/Andrew Flaim")
completed = read.csv("old_iso2k_loc.csv", head = T)

shapes = c("Coral" = 20, "GlacierIce" = 20, "Sclerosponge" = 19, "Wood" = 18,
           "LakeSediment" = 17, "Speleothem" = 17, "MolluskShells" = 19,
           "MarineSediment" = 8)

ggplot() +
  geom_polygon(data = wrld_simpl, aes(x = long, y = lat, group = group), 
               fill = "grey", colour = "black", alpha = 0.2) +
  geom_point(data = allRecs, mapping = aes(x = lon, y = lat, 
                                               shape = archive), color = "black", size = 5) +
  geom_point(data = allRecs, mapping = aes(x = lon, y = lat, 
                                               shape = archive), color = "salmon", size = 4) +
  geom_point(data = completed, mapping = aes(x = lon, y = lat, 
                                           shape = archive), color = "forest green", size = 4) +
  scale_shape_manual(values = shapes) + 
  # Removes Axes and labels
  scale_x_continuous(breaks = NULL) +
  xlab("") + 
  ylab("") +
  # Change theme to remove axes and ticks
  theme(panel.background = element_blank(),
        axis.ticks=element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank())

#
#### Export coral dataframe ####
setwd("~/GitHub/iso-project/Iso2k_PRYSM_runs")
write.csv(allRecs, "iso2k_corals.csv", row.names = F)
