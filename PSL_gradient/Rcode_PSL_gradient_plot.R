# Import and plot PSL gradients from iCESM
setwd("~/GitHub/iso-project/PSL_gradient")

time = seq(from = 1850, to = 2006, by = 1/12)
time = time[-1]

#### Full forcing ####
f1 = read.csv("PSL_gradient_F1.csv", head = T)[,2]
f2 = read.csv("PSL_gradient_F2.csv", head = T)[,2]
f3 = read.csv("PSL_gradient_F3.csv", head = T)[,2]

full = cbind(f1, f2, f3)
mean = rowMeans(full)

plot(time, f1, type = "l", frame.plot = F, axes = F, xlab = "", ylab = "",
     ylim = c(-4, 4), col = "grey")
par(new = T)
plot(time, f2, type = "l", frame.plot = F, axes = F, xlab = "", ylab = "",
     ylim = c(-4, 4), col = "grey")
par(new = T)
plot(time, f3, type = "l", frame.plot = F, axes = F, xlab = "", ylab = "",
     ylim = c(-4, 4), col = "grey")
par(new = T)
plot(time, mean, type = "l", frame.plot = F, axes = F, xlab = "", ylab = "",
     ylim = c(-4, 4), col = "black", main = "iCESM full forcing ensemble MPSL index")
axis(1)
axis(2)

#### GHG forcing ####
ghg1 = read.csv("PSL_gradient_GHG.csv", head = T)[,2]
plot(time, ghg1, type = "l", frame.plot = F, axes = F, xlab = "", ylab = "",
     ylim = c(-4, 4), col = "black", main = "iCESM GHG MPSL index")
axis(1)
axis(2)


#### Full & GHG ####
plot(time, mean, type = "l", frame.plot = F, axes = F, xlab = "", ylab = "",
     ylim = c(-4, 4), col = "black", main = "iCESM full forcing & GHG MPSL index", lwd = 2)
par(new = T)
plot(time, ghg1, type = "l", frame.plot = F, axes = F, xlab = "", ylab = "",
     ylim = c(-4, 4), col = "forest green", main = "")
axis(1)
axis(2)
