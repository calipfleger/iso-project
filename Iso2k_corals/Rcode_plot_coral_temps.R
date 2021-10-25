setwd("D:/GitHub/iso-project/Iso2k_corals")
coral_temps = read.csv("full_forcing1_coral_TS.csv", head = T)
coral_meta = read.csv("iso2k_corals.csv", head = T)
western = read.csv("iso2k_corals_western_pacific.csv", head = T)
north_atlantic = read.csv("iso2k_corals_north_atlantic.csv", head = T)
south_china = read.csv("iso2k_corals_south_china_sea.csv", head = T)
indian = read.csv("iso2k_corals_indian_ocean.csv", head = T)
north_pacific = read.csv("iso2k_corals_north_pacific.csv", head = T)
red = read.csv("iso2k_corals_red_sea.csv", head = T)
south_pacific = read.csv("iso2k_corals_south_pacfic.csv", head = T)

# Create new time variable and remove the model output time
time = seq(from = 1850+(1/12), to = 2006, by = 1/12)
coral_temps = coral_temps[,-1]
coral_temps = coral_temps-273.15

# Plot all temperature series
for (i in 1:ncol(coral_temps)){
  print(i)
  plot(time, coral_temps[,i], axes = F, xlab = "", ylab = "", frame.plot = F,
       type = "l", col = "grey", ylim = range(coral_temps))
  par(new = T)
}
axis(1)
axis(2)

#### Extract timeseries for each region ####
# Western Pacific
which_west = vector(length = nrow(western))
for (i in 1:nrow(western)){
  which_west[i] = which(coral_meta$record == western$record[i])
}
west_ts = coral_temps[,which_west]

# North Atlantic
which_NAtl = vector(length = nrow(north_atlantic))
for (i in 1:nrow(north_atlantic)){
  which_NAtl[i] = which(coral_meta$record == north_atlantic$record[i])
}
Natl_ts = coral_temps[,which_NAtl]

# South China Sea
which_china = vector(length = nrow(south_china))
for (i in 1:nrow(south_china)){
  which_china[i] = which(coral_meta$record == south_china$record[i])
}
china_ts = coral_temps[,which_china]

# Indian ocean
which_india = vector(length = nrow(indian))
for (i in 1:nrow(indian)){
  which_india[i] = which(coral_meta$record == indian$record[i])
}
india_ts = coral_temps[,which_india]

# North Pacific
which_NPac = vector(length = nrow(north_pacific))
for (i in 1:nrow(north_pacific)){
  which_NPac[i] = which(coral_meta$record == north_pacific$record[i])
}
NPac_ts = coral_temps[,which_NPac]

# Red Sea
which_red = vector(length = nrow(red))
for (i in 1:nrow(red)){
  which_red[i] = which(coral_meta$record == red$record[i])
}
red_ts = coral_temps[,which_red]

# South Pacific
which_SPac = vector(length = nrow(south_pacific))
for (i in 1:nrow(south_pacific)){
  which_SPac[i] = which(coral_meta$record == south_pacific$record[i])
}
SPac_ts = coral_temps[,which_SPac]

#### Plot each region ####
par(mfrow = c(4, 2), mai = c(0, 0.5, 0.5, 0.5), omi = c(0.5, 0.5, 0.5, 0))
# Western Pacific
recs = west_ts
for (i in 1:ncol(recs)){
  print(i)
  plot(time, recs[,i], axes = F, xlab = "", ylab = "", frame.plot = F,
       type = "l", col = "grey", ylim = range(coral_temps))
  par(new = T)
}
plot(time, rowMeans(recs),  axes = F, xlab = "", ylab = "", frame.plot = F,
     type = "l", col = "black", ylim = range(coral_temps), lwd = 2)
axis(2)
mtext(side = 3, line = 2, "Western Pacific")
mtext(side = 3, line = 1, text = c("n =           ", ncol(recs)))

# North Atlantic
recs = Natl_ts
for (i in 1:ncol(recs)){
  print(i)
  plot(time, recs[,i], axes = F, xlab = "", ylab = "", frame.plot = F,
       type = "l", col = "grey", ylim = range(coral_temps))
  par(new = T)
}
plot(time, rowMeans(recs),  axes = F, xlab = "", ylab = "", frame.plot = F,
     type = "l", col = "black", ylim = range(coral_temps), lwd = 2)
axis(2)
mtext(side = 3, line = 2, "North Atlantic")
mtext(side = 3, line = 1, text = c("n =           ", ncol(recs)))

# South China Sea
recs = china_ts
for (i in 1:ncol(recs)){
  print(i)
  plot(time, recs[,i], axes = F, xlab = "", ylab = "", frame.plot = F,
       type = "l", col = "grey", ylim = range(coral_temps))
  par(new = T)
}
plot(time, rowMeans(recs),  axes = F, xlab = "", ylab = "", frame.plot = F,
     type = "l", col = "black", ylim = range(coral_temps), lwd = 2)
axis(2)
mtext(side = 3, line = 2, "South China Sea")
mtext(side = 3, line = 1, text = c("n =           ", ncol(recs)))

# Indian Ocean
recs = india_ts
for (i in 1:ncol(recs)){
  print(i)
  plot(time, recs[,i], axes = F, xlab = "", ylab = "", frame.plot = F,
       type = "l", col = "grey", ylim = range(coral_temps))
  par(new = T)
}
plot(time, rowMeans(recs),  axes = F, xlab = "", ylab = "", frame.plot = F,
     type = "l", col = "black", ylim = range(coral_temps), lwd = 2)
axis(2)
mtext(side = 3, line = 2, "Indian Ocean")
mtext(side = 3, line = 1, text = c("n =           ", ncol(recs)))

# North Pacific
recs = NPac_ts
for (i in 1:ncol(recs)){
  print(i)
  plot(time, recs[,i], axes = F, xlab = "", ylab = "", frame.plot = F,
       type = "l", col = "grey", ylim = range(coral_temps))
  par(new = T)
}
plot(time, rowMeans(recs),  axes = F, xlab = "", ylab = "", frame.plot = F,
     type = "l", col = "black", ylim = range(coral_temps), lwd = 2)
axis(2)
mtext(side = 3, line = 2, "North Pacific")
mtext(side = 3, line = 1, text = c("n =           ", ncol(recs)))

# Red Sea
recs = red_ts
for (i in 1:ncol(recs)){
  print(i)
  plot(time, recs[,i], axes = F, xlab = "", ylab = "", frame.plot = F,
       type = "l", col = "grey", ylim = range(coral_temps))
  par(new = T)
}
plot(time, rowMeans(recs),  axes = F, xlab = "", ylab = "", frame.plot = F,
     type = "l", col = "black", ylim = range(coral_temps), lwd = 2)
axis(2)
mtext(side = 3, line = 2, "Red Sea")
mtext(side = 3, line = 1, text = c("n =           ", ncol(recs)))

# South Pacific
recs = SPac_ts
for (i in 1:ncol(recs)){
  print(i)
  plot(time, recs[,i], axes = F, xlab = "", ylab = "", frame.plot = F,
       type = "l", col = "grey", ylim = range(coral_temps))
  par(new = T)
}
plot(time, rowMeans(recs), axes = F, xlab = "", ylab = "", frame.plot = F,
     type = "l", col = "black", ylim = range(coral_temps), lwd = 2)
axis(1)
axis(2)
mtext(side = 3, line = 2, "South Pacific")
mtext(side = 3, line = 1, text = c("n =           ", ncol(recs)))

# Means
ylims = c(22, 29)
plot(time, rowMeans(NPac_ts), axes = F, xlab = "", ylab = "", frame.plot = F,
     type = "l", col = "#2596be75", ylim = ylims, lwd = 1)
par(new = T)
plot(time, rowMeans(china_ts), axes = F, xlab = "", ylab = "", frame.plot = F,
     type = "l", col = "#1bff0075", ylim = ylims, lwd = 1)
par(new = T)
plot(time, rowMeans(india_ts), axes = F, xlab = "", ylab = "", frame.plot = F,
     type = "l", col = "#ff000075", ylim = ylims, lwd = 1)
par(new = T)
plot(time, rowMeans(red_ts), axes = F, xlab = "", ylab = "", frame.plot = F,
     type = "l", col = "#0012ff75", ylim = ylims, lwd = 1)
par(new = T)
plot(time, rowMeans(SPac_ts), axes = F, xlab = "", ylab = "", frame.plot = F,
     type = "l", col = "#ff00e975", ylim = ylims, lwd = 1)
par(new = T)
plot(time, rowMeans(Natl_ts), axes = F, xlab = "", ylab = "", frame.plot = F,
     type = "l", col = "#ffce0075", ylim = ylims, lwd = 1)
par(new = T)
plot(time, rowMeans(west_ts), axes = F, xlab = "", ylab = "", frame.plot = F,
     type = "l", col = "#00000050", ylim = ylims, lwd = 1)
axis(1)
axis(2)
mtext(side = 3, line = 2, "Regional Means")
legend("bottomright", col = c("#2596be75", 
                              "#1bff0075",
                              "#ff000075",
                              "#0012ff75",
                              "#ff00e975",
                              "#ffce0075",
                              "#00000050"),
       lty = rep(1, 7), lwd = rep(3, 7),
       legend = c("North Pacific", "South China Sea",
                  "Indian Ocean", "Red Sea",
                  "South Pacific", "North Atlantic", "Western Pacific"))
