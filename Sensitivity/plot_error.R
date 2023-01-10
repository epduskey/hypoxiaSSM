# Load and plot sensitivity analysis results

# Created: December 3, 2021
# Last modified: December 30, 2022

# Set working directory
setwd(paste(mypath, "hypoxiaSSM-main", sep = ""))

# Contents (ctrl-f):
#	I. Load results


########## I. Load results ##########

# Load calibration results
cal_one = read.table("Sensitivity/Output/cal_one.txt", header = T)
cal_two = read.table("Sensitivity/Output/cal_two.txt", header = T)
cal_three = read.table("Sensitivity/Output/cal_three.txt", header = T)

# Load projection results
proj_one = read.table("Sensitivity/Output/proj_one.txt", header = T)
proj_two = read.table("Sensitivity/Output/proj_two.txt", header = T)
proj_three = read.table("Sensitivity/Output/proj_three.txt", header = T)


########## VI. Plot calibration error ##########

# Plot colors
error_colors = palette.colors(palette = "Okabe-Ito")[1:4]

# Axis label colors
color = col2rgb(rep(c("#88CCEE","#DDCC77","#CC6677","#AA4499","#332288","#117733","#44AA99","#000000"),2))

# Model labels
model_labels = cal_one$model

# Initialize plot
jpeg("Sensitivity/sens_error_cal.jpg", width = 22, height = 9, units = 'cm', res = 600)
par(mfrow = c(1,3), oma = c(3,3,0,0), mar = c(5,2,2,2))

# Plot calibration for first repetition
plot(cal_one[order(cal_one$rank),]$rank, pch = 21, bg = error_colors[1], cex = 1.5, ylim = c(0,1), xlab = "", ylab = "", axes = F); axis(2); box()
for(i in seq(16)) axis(1, at = i, labels = model_labels[order(cal_one$rank)][i], col.axis = color[order(cal_one$rank)][i], las = 2)
points(cal_one[order(cal_one$rank),]$SSB, pch = 22, bg = error_colors[2], cex = 1.5)
points(cal_one[order(cal_one$rank),]$yield, pch = 23, bg = error_colors[3], cex = 1.5)
points(cal_one[order(cal_one$rank),]$growth, pch = 24, bg = error_colors[4], cex = 1.5)

# Y-axis label
mtext("Weighted error", 2, cex = 1.2, line = 3)

# Plot calibration for second repetition
plot(cal_two[order(cal_two$rank),]$rank, pch = 21, bg = error_colors[1], cex = 1.5, ylim = c(0,1), xlab = "", ylab = "", main = "Calibration", cex.main = 1.5, axes = F); axis(2); box()
for(i in seq(16)) axis(1, at = i, labels = model_labels[order(cal_two$rank)][i], col.axis = color[order(cal_two$rank)][i], las = 2)
points(cal_two[order(cal_two$rank),]$SSB, pch = 22, bg = error_colors[2], cex = 1.5)
points(cal_two[order(cal_two$rank),]$yield, pch = 23, bg = error_colors[3], cex = 1.5)
points(cal_two[order(cal_two$rank),]$growth, pch = 24, bg = error_colors[4], cex = 1.5)

# Plot calibration for third repetition
plot(cal_three[order(cal_three$rank),]$rank, pch = 21, bg = error_colors[1], cex = 1.5, ylim = c(0,1), xlab = "", ylab = "", axes = F); axis(2); box()
for(i in seq(16)) axis(1, at = i, labels = model_labels[order(cal_three$rank)][i], col.axis = color[order(cal_three$rank)][i], las = 2)
points(cal_three[order(cal_three$rank),]$SSB, pch = 22, bg = error_colors[2], cex = 1.5)
points(cal_three[order(cal_three$rank),]$yield, pch = 23, bg = error_colors[3], cex = 1.5)
points(cal_three[order(cal_three$rank),]$growth, pch = 24, bg = error_colors[4], cex = 1.5)

# Legend
par(fig = c(0,1,0,1), oma = c(0,1,0,0), mar = c(0,0,0,0), new = T)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend("bottom", bty = 'n', c("SSB", "Yield", "Growth", "Total"), pt.bg = error_colors[c(2,3,4,1)], pch = c(22,23,24,21), pt.cex = 1.5, xpd = TRUE, horiz = TRUE, cex = 1.2, seg.len=1)

# Finish plot
dev.off()


########## VII. Plot projection error ##########

# Initialize plot
jpeg("Sensitivity/sens_error_sim.jpg", width = 22, height = 9, units = 'cm', res = 600)
par(mfrow = c(1,3), oma = c(3,3,0,0), mar = c(5,2,2,2))

# Plot calibration for first repetition
plot(proj_one[order(proj_one$rank),]$rank, pch = 21, bg = error_colors[1], cex = 1.5, ylim = c(0,1), xlab = "", ylab = "", axes = F); axis(2); box()
for(i in seq(16)) axis(1, at = i, labels = model_labels[order(proj_one$rank)][i], col.axis = color[order(proj_one$rank)][i], las = 2)
points(proj_one[order(proj_one$rank),]$SSB, pch = 22, bg = error_colors[2], cex = 1.5)
points(proj_one[order(proj_one$rank),]$yield, pch = 23, bg = error_colors[3], cex = 1.5)
points(proj_one[order(proj_one$rank),]$growth, pch = 24, bg = error_colors[4], cex = 1.5)

# Y-axis label
mtext("Weighted error", 2, cex = 1.2, line = 3)

# Plot calibration for second repetition
plot(proj_two[order(proj_two$rank),]$rank, pch = 21, bg = error_colors[1], cex = 1.5, ylim = c(0,1), xlab = "", ylab = "", main = "Projection", cex.main = 1.5, axes = F); axis(2); box()
for(i in seq(16)) axis(1, at = i, labels = model_labels[order(proj_two$rank)][i], col.axis = color[order(proj_two$rank)][i], las = 2)
points(proj_two[order(proj_two$rank),]$SSB, pch = 22, bg = error_colors[2], cex = 1.5)
points(proj_two[order(proj_two$rank),]$yield, pch = 23, bg = error_colors[3], cex = 1.5)
points(proj_two[order(proj_two$rank),]$growth, pch = 24, bg = error_colors[4], cex = 1.5)

# Plot calibration for third repetition
plot(proj_three[order(proj_three$rank),]$rank, pch = 21, bg = error_colors[1], cex = 1.5, ylim = c(0,1), xlab = "", ylab = "", axes = F); axis(2); box()
for(i in seq(16)) axis(1, at = i, labels = model_labels[order(proj_three$rank)][i], col.axis = color[order(proj_three$rank)][i], las = 2)
points(proj_three[order(proj_three$rank),]$SSB, pch = 22, bg = error_colors[2], cex = 1.5)
points(proj_three[order(proj_three$rank),]$yield, pch = 23, bg = error_colors[3], cex = 1.5)
points(proj_three[order(proj_three$rank),]$growth, pch = 24, bg = error_colors[4], cex = 1.5)

# Legend
par(fig = c(0,1,0,1), oma = c(0,1,0,0), mar = c(0,0,0,0), new = T)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend("bottom", bty = 'n', c("SSB", "Yield", "Growth", "Total"), pt.bg = error_colors[c(2,3,4,1)], pch = c(22,23,24,21), pt.cex = 1.5, xpd = TRUE, horiz = TRUE, cex = 1.2, seg.len=1)

# Finish plot
dev.off()
