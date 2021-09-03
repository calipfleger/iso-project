rec_col = c("blue", "red", "dark green", "purple", "black", "orange", "brown")
rec_col_temp = c("red", "dark green", "purple", "orange", "brown")

plot(time, SWBB_ensemble_mean, type = "l", frame.plot = F, ylim = c(-0.4,0.4),
     xlab = "", ylab = "", axes = F, lwd = 2, col = rec_col[2])
par(new = T)
plot(time, KIPR_ensemble_mean, type = "l", frame.plot = F, ylim = c(-0.4,0.4),
     xlab = "", ylab = "", axes = F, lwd = 2, col = rec_col[3])
par(new = T)
plot(time, LIRA_ensemble_mean, type = "l", frame.plot = F, ylim = c(-0.4,0.4),
     xlab = "", ylab = "", axes = F, lwd = 2, col = rec_col[4])
par(new = T)
plot(time, COPM_ensemble_mean, type = "l", frame.plot = F, ylim = c(-0.4,0.4),
     xlab = "", ylab = "", axes = F, lwd = 2, col = rec_col[6])
par(new = T)
plot(time, BOFP_ensemble_mean, type = "l", frame.plot = F, ylim = c(-0.4,0.4),
     xlab = "", ylab = "", axes = F, lwd = 2, col = rec_col[7])
par(new = T)
plot(NA, frame.plot = F, ylim = c(-0.4,0.4), xlim = c(min(time), max(time)),
     ylab = "Pseudo d18O anomalies", xlab = "Years (CE)", axes = T, 
     main = "Pseudocoral Ensemble Means")
legend("topright", col = rec_col_temp, 
       legend = c("SWBB", "KIPR", "LIRA", "COPM", "BOFP"),
       lty = rep(1, 3))
