# Import PRYSM pseudocoral result
setwd("~/GitHub/iso-project/Iso2k_PRYSM_runs/results")
coral = read.csv("pseudocoral.csv", head = F)
coral = as.vector(unlist(coral))

# Time dimension
start = 1851
end = 1950
time = seq(from = start, to = end, by = 1)

# Plot pseudocoral TS
plot(time, coral, type = "l", frame.plot = F)
