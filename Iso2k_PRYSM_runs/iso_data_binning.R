library(geoChronR)

#### Bin data ####
# Bin monthly SST
setwd("~/GitHub/iso-project/Iso2k_PRYSM_runs/data/KIPR/KIPR_afedits")
ssta_m = read.csv("KIPR_1_afedit.csv", head = T)[3]
start = 1851
end = 1981
time = seq(from = start, to = end, by = 1/12)
ssta_m = as.vector(unlist(ssta_m))
binvec =  seq(from = start, to = end, by = 1)

binned_ann_recs = matrix(NA, length(binvec), length(ssta_m)+1)

thisrec = na.omit(data.frame(year = time, val = ssta_m))

binned_ann_recs = geoChronR::bin(thisrec$year, thisrec$val, binvec)

write.csv(binned_ann_recs, "ssta_KIPR1_binned.csv")

# Bin monthly sss
sssa_m = read.csv("KIPR_1_afedit.csv", head = T)[2]
start = 1851
end = 1981
time = seq(from = start, to = end, by = 1/12)
sssa_m = as.vector(unlist(sssa_m))
binvec =  seq(from = start, to = end, by = 1)

binned_ann_recs = matrix(NA, length(binvec), length(sssa_m)+1)

thisrec = na.omit(data.frame(year = time, val = sssa_m))

binned_ann_recs = geoChronR::bin(thisrec$year, thisrec$val, binvec)

write.csv(binned_ann_recs, "sssa_KIPR1_binned.csv")
