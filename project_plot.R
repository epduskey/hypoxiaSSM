# Load and plot results

# Created: May 18, 2021
# Last modified: January 4, 2023 by EPD

# Set working directory
setwd(paste(mypath, "hypoxiaSSM-main", sep = ""))

# Load packages
library(mizer)
library(RColorBrewer)
library(png)
library(xtable)

# Contents (ctrl-f):
#	I. Source scripts
#	II. Load results
#	III. Prepare and run calibrations
#	IV. Prepare and run simulations
#	V. Conceptual figure
#	VI. Plot oxygen
#	VII. Scaling plots
#	VIII. Plot growth during calibration period
#	IX. Plot SSB during projection period
#	X. Plot yield during projection period
#	XI. Plot growth during projection period
#	XII. Plot growth during projection period for cod only
#	XIII. Plot growth during projection period for first eight years of life
#	XIV. Plot cod occupancy
#	XV. Calculate error
#	XVI. Plot calibration error
#	XVII. Plot calibration error rank
#	XVIII. Plot projection error
#	XIX. Plot projection error rank
#	XX. Plot final weighted error rank for calibration
#	XXI. Plot final weighted error rank for projection


########## I. Source scripts ##########

# Source scripts
#	spmat: creates species_params data frame
#	bpsetup: organizes parameters for benthic and pelagic habitats
#	metabscale: scaling maintenance costs i.e. metabolic rate with oxygen and temperature
#	oxysensitivity: organized parameters for sensitivity to hypoxia
#	otscale: scaling occupancy and physiological rates with oxygen and temperature
#	rate_encounter: scaled encounter rate
#	rate_erag: scaled availability of energy for reproduction and growth
#	rate_feeding: scaled feeding level
#	rate_mortality: mortality including direct mortality due to hypoxia exposure
#	rate_predation: scaled predation rate
#	rate_RDI: scaled reproduction
#	rate_resourcemort: scaled resource mortality rate
#	resource_benthic: benthic resource function
#	resource_pelagic: pelagic resource function
#	plot_funcs: assorted functions to process and plot mizer sim output
#	plot_spectra: plot spectra with both benthic and pelagic resources
#	plot_growth: plot growth at a given time step
#	plot_yield: plot yield at a given time step
#	plot_diet: custom plot diet functions
#	plot_fmsy: plot terminal yield for each species at a given F
#	cal_occupancy: calibrate occupancy sensitivity parameters
#	cal_biomass: calibrate for biomass on average and dynamically
#	cal_growth: calibrate for growth on average and dynamically
#	cal_yield: calibrate for yield on average and dynamically
#	mizer_replace: altered get and set rate functions to operate with custom functions
#	fsetup: get annual effort arrays
source("Code/spmat.R")
source("Code/bpsetup.R")
source("Code/oxysensitivity.R")
source("Code/otmscale.R")
source("Code/rate_encounter.R")
source("Code/rate_erag.R")
source("Code/rate_feeding.R")
source("Code/rate_mortality.R")
source("Code/rate_predation.R")
source("Code/rate_predmort.R")
source("Code/rate_RDI.R")
source("Code/rate_resourcemort.R")
source("Code/resource_benthic.R")
source("Code/resource_pelagic.R")
source("Code/plot_funcs.R")
source("Code/plot_spectra.R")
source("Code/plot_growth.R")
source("Code/plot_yield.R")
source("Code/plot_diet.R")
source("Code/plot_fmsy.R")
source("Code/cal_occupancy.R")
source("Code/cal_biomass.R")
source("Code/cal_growth.R")
source("Code/cal_yield.R")
source("Code/mizer_replace.R")
source("Code/fsetup.R")


########## II. Load results ##########

# Benthos, occupancy, mortality, and physiological scaling
load("Calibration/params_full_bomp.rda")
params_bomp = params_full

# Benthos, occupancy, and mortality scaling
load("Calibration/params_full_bom.rda")
params_bom = params_full

# Benthos, occupancy, and physiological scaling
load("Calibration/params_full_bop.rda")
params_bop = params_full

# Benthos, mortality, and physiological scaling
load("Calibration/params_full_bmp.rda")
params_bmp = params_full

# Occupancy, mortality, and physiological scaling
load("Calibration/params_full_omp.rda")
params_omp = params_full

# Benthos and occupancy scaling
load("Calibration/params_full_bo.rda")
params_bo = params_full

# Benthos and mortality scaling
load("Calibration/params_full_bm.rda")
params_bm = params_full

# Benthos and physiological scaling
load("Calibration/params_full_bp.rda")
params_bp = params_full

# Occupancy and mortality scaling
load("Calibration/params_full_om.rda")
params_om = params_full

# Occupancy and physiological scaling
load("Calibration/params_full_op.rda")
params_op = params_full

# Occupancy and physiological scaling
load("Calibration/params_full_mp.rda")
params_mp = params_full

# Benthos scaling
load("Calibration/params_full_b.rda")
params_b = params_full

# Occupancy scaling
load("Calibration/params_full_o.rda")
params_o = params_full

# Mortality scaling
load("Calibration/params_full_m.rda")
params_m = params_full

# Physiological scaling
load("Calibration/params_full_p.rda")
params_p = params_full

# None scaling (with left fish)
load("Calibration/params_full_none.rda")
params_none = params_full


########## III. Prepare and run calibrations ##########

# Load oxygen data
load("Data/oxy_model.rda")

# Create an oxygen vector
t_tune = 100
times = 0:t_tune
benthic_oxygen_cal = vector(mode = "numeric", length = length(times))
pelagic_oxygen_cal = vector(mode = "numeric", length = length(times))

# Varying oxygen scenario
benthic_oxygen_cal[1:length(times)] = c(rep(mean(preds_all$df[preds_all$df$Year %in% seq(1991,2000), ]$Oxygen_benthic),t_tune-dim(preds_all$df[preds_all$df$Year %in% seq(1991,2000), ])[1]+1), preds_all$df[preds_all$df$Year %in% seq(1991,2000), ]$Oxygen_benthic)
pelagic_oxygen_cal[1:length(times)] = c(rep(mean(preds_all$df[preds_all$df$Year %in% seq(1991,2000), ]$Oxygen_pelagic),t_tune-dim(preds_all$df[preds_all$df$Year %in% seq(1991,2000), ])[1]+1), preds_all$df[preds_all$df$Year %in% seq(1991,2000), ]$Oxygen_pelagic)

# Store these in other_params
params_bomp_cal = params_bomp
params_bomp_cal@other_params$benthic_oxygen = benthic_oxygen_cal
params_bomp_cal@other_params$pelagic_oxygen = pelagic_oxygen_cal

params_bom_cal = params_bom
params_bom_cal@other_params$benthic_oxygen = benthic_oxygen_cal
params_bom_cal@other_params$pelagic_oxygen = pelagic_oxygen_cal

params_bop_cal = params_bop
params_bop_cal@other_params$benthic_oxygen = benthic_oxygen_cal
params_bop_cal@other_params$pelagic_oxygen = pelagic_oxygen_cal

params_bmp_cal = params_bmp
params_bmp_cal@other_params$benthic_oxygen = benthic_oxygen_cal
params_bmp_cal@other_params$pelagic_oxygen = pelagic_oxygen_cal

params_omp_cal = params_omp
params_omp_cal@other_params$benthic_oxygen = benthic_oxygen_cal
params_omp_cal@other_params$pelagic_oxygen = pelagic_oxygen_cal

params_bo_cal = params_bo
params_bo_cal@other_params$benthic_oxygen = benthic_oxygen_cal
params_bo_cal@other_params$pelagic_oxygen = pelagic_oxygen_cal

params_bm_cal = params_bm
params_bm_cal@other_params$benthic_oxygen = benthic_oxygen_cal
params_bm_cal@other_params$pelagic_oxygen = pelagic_oxygen_cal

params_bp_cal = params_bp
params_bp_cal@other_params$benthic_oxygen = benthic_oxygen_cal
params_bp_cal@other_params$pelagic_oxygen = pelagic_oxygen_cal

params_om_cal = params_om
params_om_cal@other_params$benthic_oxygen = benthic_oxygen_cal
params_om_cal@other_params$pelagic_oxygen = pelagic_oxygen_cal

params_op_cal = params_op
params_op_cal@other_params$benthic_oxygen = benthic_oxygen_cal
params_op_cal@other_params$pelagic_oxygen = pelagic_oxygen_cal

params_mp_cal = params_mp
params_mp_cal@other_params$benthic_oxygen = benthic_oxygen_cal
params_mp_cal@other_params$pelagic_oxygen = pelagic_oxygen_cal

params_b_cal = params_b
params_b_cal@other_params$benthic_oxygen = benthic_oxygen_cal
params_b_cal@other_params$pelagic_oxygen = pelagic_oxygen_cal

params_o_cal = params_bomp
params_o_cal@other_params$benthic_oxygen = benthic_oxygen_cal
params_o_cal@other_params$pelagic_oxygen = pelagic_oxygen_cal

params_m_cal = params_m
params_m_cal@other_params$benthic_oxygen = benthic_oxygen_cal
params_m_cal@other_params$pelagic_oxygen = pelagic_oxygen_cal

params_p_cal = params_p
params_p_cal@other_params$benthic_oxygen = benthic_oxygen_cal
params_p_cal@other_params$pelagic_oxygen = pelagic_oxygen_cal

params_none_cal = params_none
params_none_cal@other_params$benthic_oxygen = benthic_oxygen_cal
params_none_cal@other_params$pelagic_oxygen = pelagic_oxygen_cal

# Scale rates by oxygen and temperature: see otmscale.R
params_bomp_cal = occupancy(params_bomp_cal, t_tune)
params_bomp_cal = rate_scale(params_bomp_cal, t_tune)
params_bomp_cal = mscale(params_bomp_cal, t_tune)

params_bom_cal = occupancy(params_bom_cal, t_tune)

params_bop_cal = occupancy(params_bop_cal, t_tune)
params_bop_cal = rate_scale(params_bop_cal, t_tune)
params_bop_cal = mscale(params_bop_cal, t_tune)

params_bmp_cal = rate_scale(params_bmp_cal, t_tune)
params_bmp_cal = mscale(params_bmp_cal, t_tune)

params_omp_cal = occupancy(params_omp_cal, t_tune)
params_omp_cal = rate_scale(params_omp_cal, t_tune)
params_omp_cal = mscale(params_omp_cal, t_tune)

params_bo_cal = occupancy(params_bo_cal, t_tune)

params_bp_cal = rate_scale(params_bp_cal, t_tune)
params_bp_cal = mscale(params_bp_cal, t_tune)

params_om_cal = occupancy(params_om_cal, t_tune)

params_op_cal = occupancy(params_op_cal, t_tune)
params_op_cal = rate_scale(params_op_cal, t_tune)
params_op_cal = mscale(params_op_cal, t_tune)

params_mp_cal = rate_scale(params_mp_cal, t_tune)
params_mp_cal = mscale(params_mp_cal, t_tune)

params_o_cal = occupancy(params_o_cal, t_tune)

params_p_cal = rate_scale(params_p_cal, t_tune)
params_p_cal = mscale(params_p_cal, t_tune)

# Get effort array
effort_cal =  farray(params_bomp, fdat, seq(1991,2000), c(2,2,2,2), t_tune)

# Run calibration model
cal_bomp = project(params_bomp_cal, t_max = 100, effort = effort_cal)
cal_bom = project(params_bom_cal, t_max = 100, effort = effort_cal)
cal_bop = project(params_bop_cal, t_max = 100, effort = effort_cal)
cal_bmp = project(params_bmp_cal, t_max = 100, effort = effort_cal)
cal_omp = project(params_omp_cal, t_max = 100, effort = effort_cal)
cal_bo = project(params_bo_cal, t_max = 100, effort = effort_cal)
cal_bm = project(params_bm_cal, t_max = 100, effort = effort_cal)
cal_bp = project(params_bp_cal, t_max = 100, effort = effort_cal)
cal_om = project(params_om_cal, t_max = 100, effort = effort_cal)
cal_op = project(params_op_cal, t_max = 100, effort = effort_cal)
cal_mp = project(params_mp_cal, t_max = 100, effort = effort_cal)
cal_b = project(params_b_cal, t_max = 100, effort = effort_cal)
cal_o = project(params_o_cal, t_max = 100, effort = effort_cal)
cal_m = project(params_m_cal, t_max = 100, effort = effort_cal)
cal_p = project(params_p_cal, t_max = 100, effort = effort_cal)
cal_none = project(params_none_cal, t_max = 100, effort = effort_cal)


########## IV. Prepare and run simulations ##########

# Create an oxygen vector
t_tune = 100
times = 0:t_tune
benthic_oxygen_sim = vector(mode = "numeric", length = length(times))
pelagic_oxygen_sim = vector(mode = "numeric", length = length(times))

# Varying oxygen scenario
benthic_oxygen_sim[1:length(times)] = c(rep(mean(preds_all$df[preds_all$df$Year %in% seq(1991,2000), ]$Oxygen_benthic),t_tune-dim(preds_all$df[preds_all$df$Year %in% seq(1991,2019), ])[1]+1), preds_all$df[preds_all$df$Year %in% seq(1991,2019), ]$Oxygen_benthic)
pelagic_oxygen_sim[1:length(times)] = c(rep(mean(preds_all$df[preds_all$df$Year %in% seq(1991,2000), ]$Oxygen_pelagic),t_tune-dim(preds_all$df[preds_all$df$Year %in% seq(1991,2019), ])[1]+1), preds_all$df[preds_all$df$Year %in% seq(1991,2019), ]$Oxygen_pelagic)

# Store these in other_params
params_bomp_sim = params_bomp
params_bomp_sim@other_params$benthic_oxygen = benthic_oxygen_sim
params_bomp_sim@other_params$pelagic_oxygen = pelagic_oxygen_sim

params_bom_sim = params_bom
params_bom_sim@other_params$benthic_oxygen = benthic_oxygen_sim
params_bom_sim@other_params$pelagic_oxygen = pelagic_oxygen_sim

params_bop_sim = params_bop
params_bop_sim@other_params$benthic_oxygen = benthic_oxygen_sim
params_bop_sim@other_params$pelagic_oxygen = pelagic_oxygen_sim

params_bmp_sim = params_bmp
params_bmp_sim@other_params$benthic_oxygen = benthic_oxygen_sim
params_bmp_sim@other_params$pelagic_oxygen = pelagic_oxygen_sim

params_omp_sim = params_omp
params_omp_sim@other_params$benthic_oxygen = benthic_oxygen_sim
params_omp_sim@other_params$pelagic_oxygen = pelagic_oxygen_sim

params_bo_sim = params_bo
params_bo_sim@other_params$benthic_oxygen = benthic_oxygen_sim
params_bo_sim@other_params$pelagic_oxygen = pelagic_oxygen_sim

params_bm_sim = params_bm
params_bm_sim@other_params$benthic_oxygen = benthic_oxygen_sim
params_bm_sim@other_params$pelagic_oxygen = pelagic_oxygen_sim

params_bp_sim = params_bp
params_bp_sim@other_params$benthic_oxygen = benthic_oxygen_sim
params_bp_sim@other_params$pelagic_oxygen = pelagic_oxygen_sim

params_om_sim = params_om
params_om_sim@other_params$benthic_oxygen = benthic_oxygen_sim
params_om_sim@other_params$pelagic_oxygen = pelagic_oxygen_sim

params_op_sim = params_op
params_op_sim@other_params$benthic_oxygen = benthic_oxygen_sim
params_op_sim@other_params$pelagic_oxygen = pelagic_oxygen_sim

params_mp_sim = params_mp
params_mp_sim@other_params$benthic_oxygen = benthic_oxygen_sim
params_mp_sim@other_params$pelagic_oxygen = pelagic_oxygen_sim

params_b_sim = params_b
params_b_sim@other_params$benthic_oxygen = benthic_oxygen_sim
params_b_sim@other_params$pelagic_oxygen = pelagic_oxygen_sim

params_o_sim = params_o
params_o_sim@other_params$benthic_oxygen = benthic_oxygen_sim
params_o_sim@other_params$pelagic_oxygen = pelagic_oxygen_sim

params_m_sim = params_m
params_m_sim@other_params$benthic_oxygen = benthic_oxygen_sim
params_m_sim@other_params$pelagic_oxygen = pelagic_oxygen_sim

params_p_sim = params_p
params_p_sim@other_params$benthic_oxygen = benthic_oxygen_sim
params_p_sim@other_params$pelagic_oxygen = pelagic_oxygen_sim

params_none_sim = params_none
params_none_sim@other_params$benthic_oxygen = benthic_oxygen_sim
params_none_sim@other_params$pelagic_oxygen = pelagic_oxygen_sim

# Scale rates by oxygen and temperature: see otmscale.R
params_bomp_sim = occupancy(params_bomp_sim, t_tune)
params_bomp_sim = rate_scale(params_bomp_sim, t_tune)
params_bomp_sim = mscale(params_bomp_sim, t_tune)

params_bom_sim = occupancy(params_bom_sim, t_tune)

params_bop_sim = occupancy(params_bop_sim, t_tune)
params_bop_sim = rate_scale(params_bop_sim, t_tune)
params_bop_sim = mscale(params_bop_sim, t_tune)

params_bmp_sim = rate_scale(params_bmp_sim, t_tune)
params_bmp_sim = mscale(params_bmp_sim, t_tune)

params_omp_sim = occupancy(params_omp_sim, t_tune)
params_omp_sim = rate_scale(params_omp_sim, t_tune)
params_omp_sim = mscale(params_omp_sim, t_tune)

params_bo_sim = occupancy(params_bo_sim, t_tune)

params_bp_sim = rate_scale(params_bp_sim, t_tune)
params_bp_sim = mscale(params_bp_sim, t_tune)

params_om_sim = occupancy(params_om_sim, t_tune)

params_op_sim = occupancy(params_op_sim, t_tune)
params_op_sim = rate_scale(params_op_sim, t_tune)
params_op_sim = mscale(params_op_sim, t_tune)

params_mp_sim = rate_scale(params_mp_sim, t_tune)
params_mp_sim = mscale(params_mp_sim, t_tune)

params_o_sim = occupancy(params_o_sim, t_tune)

params_p_sim = rate_scale(params_p_sim, t_tune)
params_p_sim = mscale(params_p_sim, t_tune)

# Get effort array
effort_full = farray(params_bomp, fdat, seq(1991,2019), c(2,2,2,2), t_tune)

# Run full model
sim_bomp = project(params_bomp_sim, t_max = 100, effort = effort_full)
sim_bom = project(params_bom_sim, t_max = 100, effort = effort_full)
sim_bop = project(params_bop_sim, t_max = 100, effort = effort_full)
sim_bmp = project(params_bmp_sim, t_max = 100, effort = effort_full)
sim_omp = project(params_omp_sim, t_max = 100, effort = effort_full)
sim_bo = project(params_bo_sim, t_max = 100, effort = effort_full)
sim_bm = project(params_bm_sim, t_max = 100, effort = effort_full)
sim_bp = project(params_bp_sim, t_max = 100, effort = effort_full)
sim_om = project(params_om_sim, t_max = 100, effort = effort_full)
sim_op = project(params_op_sim, t_max = 100, effort = effort_full)
sim_mp = project(params_mp_sim, t_max = 100, effort = effort_full)
sim_b = project(params_b_sim, t_max = 100, effort = effort_full)
sim_o = project(params_o_sim, t_max = 100, effort = effort_full)
sim_m = project(params_m_sim, t_max = 100, effort = effort_full)
sim_p = project(params_p_sim, t_max = 100, effort = effort_full)
sim_none = project(params_none_sim, t_max = 100, effort = effort_full)


########## V. Conceptual figure ##########

# Initialize plot
jpeg("Plots/conceptual.jpg", width = 28, height = 19, units = 'cm', res = 600)

# Read drawings
sad_drawing = readPNG("Plots/bd_tp.png")
cop_drawing = readPNG("Plots/pd_tp.png")
cod_drawing = readPNG("Plots/cd_tp.png")
fln_drawing = readPNG("Plots/fd_tp.png")
spt_drawing = readPNG("Plots/sd_tp.png")
her_drawing = readPNG("Plots/hd_tp.png")

# Matrix layout
concept.layout = layout(matrix(c(1,7,7,7,8,2,7,7,7,8,3,4,5,6,8), 3, 5, byrow = T))
layout.show(concept.layout)

# Growth plot
a = 0.1
b = 3
kvb = 0.2
Linf = 1200
t0 = -1
age = seq(0, 20, length.out = 1000)
size = a * (Linf * (1 - exp(-kvb * (age - t0)))) ^ b
par(mar = c(3,3,3,1))
plot(age, size, type = 'l', lwd = 10, axes = F, xlab = "", ylab = "", main = "Growth", cex.main = 2)
lines(age, size, lwd = 5, col = "gray")
box()
mtext("Age", side = 1, cex = 1.5, line = 0.5)
mtext("Weight", side = 2, cex = 1.5, line = 0.5)
text(0.05*20, 0.90*a*1200^b, "a", cex = 2)

# Feeding plot
set.seed(487)
rnorm.test = rnorm(10000)
par(mar = c(3,3,3,1))
hist(rnorm.test, axes = F, main = "Feeding", xlim = c(-4,4), ylim = c(0,2000), cex.main = 2)
box()
mtext("log(Pred w/Prey w)", side = 1, cex = 1.5, line = 0.5)
mtext("Frequency", side = 2, cex = 1.5, line = 0.5)
text(0.05*8-4, 0.90*2000, "b", cex = 2)

# Reproduction plot
c = 10
d = 0.003
stock = seq(0,1000)
recruits = (c*stock) / (1 + d*stock)
par(mar = c(3,3,3,1))
plot(stock, recruits, type = 'l', lwd = 10, axes = F, xlab = "", ylab = "", main = "Reproduction", cex.main = 2)
lines(stock, recruits, lwd = 5, col = "gray")
box()
mtext("Stock", side = 1, cex = 1.5, line = 0.5)
mtext("Recruits", side = 2, cex = 1.5, line = 0.5)
text(0.05*1000, 0.90*2500, "c", cex = 2)

# Selectivity
L25 = 50
L50 = 100
Length = seq(1,1500)
selectivity = 1/(1 + exp(L50*log(3)/(L50-L25) - (log(3)/(L50-L25))*Length))
par(mar = c(3,3,3,1))
plot(a*Length^b, selectivity, log = 'x', type = 'l', lwd = 10, axes = F, xlab = "", ylab = "", main = "Fishing Selectivity", cex.main = 2)
lines(a*Length^b, selectivity, lwd = 5, col = "gray")
box()
mtext("log(Weight)", side = 1, cex = 1.5, line = 0.5)
mtext("Probability", side = 2, cex = 1.5, line = 0.5)
text(0.2, 0.90*1, "d", cex = 2)

# Metabolism scaling
V = 0.5
m = 3
oxygen = seq(0,8,length.out=100)
mscale = 1 + exp(-V*(oxygen - m))
par(mar = c(3,3,3,1))
plot(oxygen, mscale, type = 'l', lwd = 10, axes = F, xlab = "", ylab = "", main = "Metabolic Scaling", cex.main = 2, col.main = "blue")
lines(oxygen, mscale, lwd = 5, col = "blue")
box(col = "blue")
mtext("Oxygen", side = 1, cex = 1.5, line = 0.5, col = "blue")
mtext("Scaling", side = 2, cex = 1.5, line = 0.5, col = "blue")
text(0.95*8, 0.90*(1 + exp(-V*(0 - m))), "e", cex = 2)

# Other (occupancy, consumption, assimilation scaling)
U = 3
k = 3
oscale = 1 / (1 + exp(-U*(oxygen - k)))
par(mar = c(3,3,3,1))
plot(oxygen, oscale, type = 'l', lwd = 10, axes = F, xlab = "", ylab = "", main = "Other Scaling", cex.main = 2, col.main = "blue")
lines(oxygen, oscale, lwd = 5, col = "blue")
box(col = "blue")
mtext("Oxygen", side = 1, cex = 1.5, line = 0.5, col = "blue")
mtext("Scaling", side = 2, cex = 1.5, line = 0.5, col = "blue")
text(0.05*8, 0.90*(1 / (1 + exp(-U*(8 - k)))), "f", cex = 2)

# Spectrum plot
source("Code/plot_spectra.R")
concept_spectra(sim_bomp, 100)
text(1.7*1e5, 5*1e16, "g", cex = 2)

# Fish legend
plot(c(0,1), c(0,1), type = "n", axes = F, xlab = "", ylab = "")
text(0.52, 0.940, "Benthos", pos = 2, cex = 1.4); rasterImage(sad_drawing, 0.5, 0.90, 0.70, 0.98)
text(0.52, 0.845, "Plankton", pos = 2, cex = 1.4); rasterImage(cop_drawing, 0.5, 0.82, 0.70, 0.87)
text(0.52, 0.760, "Cod", pos = 2, cex = 1.4); rasterImage(cod_drawing, 0.5, 0.72, 0.99, 0.80)
text(0.52, 0.660, "Flounder", pos = 2, cex = 1.4); rasterImage(fln_drawing, 0.5, 0.62, 0.92, 0.70)
text(0.52, 0.565, "Sprat", pos = 2, cex = 1.4); rasterImage(spt_drawing, 0.5, 0.55, 0.75, 0.58)
text(0.52, 0.475, "Herring", pos = 2, cex = 1.4); rasterImage(her_drawing, 0.5, 0.45, 0.85, 0.50)

# Finish plot
dev.off()


########## VI. Plot oxygen ##########

# Plot benthic and pelagic
jpeg("Plots/oxygen.jpg", width = 23, height = 16, units = 'cm', res = 600)
par(mfrow = c(1,1), mar = c(5,5,3,1))
plot(Oxygen_benthic ~ Year, preds_all$df, type = 'l', lwd = 5, ylim = c(0,7), xlab = "Year", ylab = "Oxygen (mL/L)", cex.lab = 2, cex.axis = 1.5)
lines(preds_all$df$Year, preds_all$df$Oxygen_pelagic, type = 'l', lwd = 5, lty = 2, col = "gray")
legend("bottomleft", bty = 'n', c("pelagic","benthic"), lty = c(2,1), lwd = 5, cex = 1.5, col = c("gray","black"))
dev.off()


########## VII. Scaling plots ##########

# Get scaling parameters for cod and flounder
bomp_scaling = cbind(U_res=rep(resource_params(params_bomp)$U_oxy,2), k_res=rep(resource_params(params_bomp)$k_oxy,2), species_params(params_bomp)[c("cod","flounder"), c("U_hab","a_hab","U_crit","a_crit","U_met","a_met","z_mort","b_mort")])
bom_scaling = cbind(U_res=rep(resource_params(params_bom)$U_oxy,2), k_res=rep(resource_params(params_bom)$k_oxy,2), species_params(params_bom)[c("cod","flounder"), c("U_hab","a_hab","U_crit","a_crit","U_met","a_met","z_mort","b_mort")])
bop_scaling = cbind(U_res=rep(resource_params(params_bop)$U_oxy,2), k_res=rep(resource_params(params_bop)$k_oxy,2), species_params(params_bop)[c("cod","flounder"), c("U_hab","a_hab","U_crit","a_crit","U_met","a_met","z_mort","b_mort")])
bmp_scaling = cbind(U_res=rep(resource_params(params_bmp)$U_oxy,2), k_res=rep(resource_params(params_bmp)$k_oxy,2), species_params(params_bmp)[c("cod","flounder"), c("U_hab","a_hab","U_crit","a_crit","U_met","a_met","z_mort","b_mort")])
omp_scaling = cbind(U_res=rep(resource_params(params_omp)$U_oxy,2), k_res=rep(resource_params(params_omp)$k_oxy,2), species_params(params_omp)[c("cod","flounder"), c("U_hab","a_hab","U_crit","a_crit","U_met","a_met","z_mort","b_mort")])
bo_scaling = cbind(U_res=rep(resource_params(params_bo)$U_oxy,2), k_res=rep(resource_params(params_bo)$k_oxy,2), species_params(params_bo)[c("cod","flounder"), c("U_hab","a_hab","U_crit","a_crit","U_met","a_met","z_mort","b_mort")])
bm_scaling = cbind(U_res=rep(resource_params(params_bm)$U_oxy,2), k_res=rep(resource_params(params_bm)$k_oxy,2), species_params(params_bm)[c("cod","flounder"), c("U_hab","a_hab","U_crit","a_crit","U_met","a_met","z_mort","b_mort")])
bp_scaling = cbind(U_res=rep(resource_params(params_bp)$U_oxy,2), k_res=rep(resource_params(params_bp)$k_oxy,2), species_params(params_bp)[c("cod","flounder"), c("U_hab","a_hab","U_crit","a_crit","U_met","a_met","z_mort","b_mort")])
om_scaling = cbind(U_res=rep(resource_params(params_om)$U_oxy,2), k_res=rep(resource_params(params_om)$k_oxy,2), species_params(params_om)[c("cod","flounder"), c("U_hab","a_hab","U_crit","a_crit","U_met","a_met","z_mort","b_mort")])
op_scaling = cbind(U_res=rep(resource_params(params_op)$U_oxy,2), k_res=rep(resource_params(params_op)$k_oxy,2), species_params(params_op)[c("cod","flounder"), c("U_hab","a_hab","U_crit","a_crit","U_met","a_met","z_mort","b_mort")])
mp_scaling = cbind(U_res=rep(resource_params(params_mp)$U_oxy,2), k_res=rep(resource_params(params_mp)$k_oxy,2), species_params(params_mp)[c("cod","flounder"), c("U_hab","a_hab","U_crit","a_crit","U_met","a_met","z_mort","b_mort")])
b_scaling = cbind(U_res=rep(resource_params(params_b)$U_oxy,2), k_res=rep(resource_params(params_b)$k_oxy,2), species_params(params_b)[c("cod","flounder"), c("U_hab","a_hab","U_crit","a_crit","U_met","a_met","z_mort","b_mort")])
o_scaling = cbind(U_res=rep(resource_params(params_o)$U_oxy,2), k_res=rep(resource_params(params_o)$k_oxy,2), species_params(params_o)[c("cod","flounder"), c("U_hab","a_hab","U_crit","a_crit","U_met","a_met","z_mort","b_mort")])
m_scaling = cbind(U_res=rep(resource_params(params_m)$U_oxy,2), k_res=rep(resource_params(params_m)$k_oxy,2), species_params(params_m)[c("cod","flounder"), c("U_hab","a_hab","U_crit","a_crit","U_met","a_met","z_mort","b_mort")])
p_scaling = cbind(U_res=rep(resource_params(params_p)$U_oxy,2), k_res=rep(resource_params(params_p)$k_oxy,2), species_params(params_p)[c("cod","flounder"), c("U_hab","a_hab","U_crit","a_crit","U_met","a_met","z_mort","b_mort")])
none_scaling = cbind(U_res=rep(resource_params(params_none)$U_oxy,2), k_res=rep(resource_params(params_none)$k_oxy,2), species_params(params_none)[c("cod","flounder"), c("U_hab","a_hab","U_crit","a_crit","U_met","a_met","z_mort","b_mort")])

# Create data frame
scale_tab = rbind(
	bomp_scaling,
	bom_scaling,
	bop_scaling,
	bmp_scaling,
	omp_scaling,
	bo_scaling,
	bm_scaling,
	bp_scaling,
	om_scaling,
	op_scaling,
	mp_scaling,
	b_scaling,
	o_scaling,
	m_scaling,
	p_scaling,
	none_scaling)
scale_tab = cbind(Model = rep(c("BOMP","BOM","BOP","BMP","OMP","BO","BM","BP","OM","OP","MP","B","O","M","P","None"), each = 2), scale_tab)
scale_tab = cbind(Species = rep(c("Cod","Flounder"), 16), scale_tab)
scale_tab = scale_tab[c(seq(1,31,by=2), seq(2,32,by=2)), ]

# Create latex table using xtable
print(xtable(scale_tab), include.rownames = F)

# Initialize plot
jpeg("Plots/all_scaling.jpg", width = 32, height = 14, units = 'cm', res = 600)

# Choose weights to plot
pwt = c(1,20,40,60,80,100)

# Colors and transparency for plotting
color = col2rgb(rep(c("#D55E00","#0072B2","#F0E442","#CC79A7","#E69F00","#009E73","#000000","#D4D4D4"),2))
lty = rep(c(2,1),each=8)
alphas = seq(1.0*255, 0.2*255, length.out = length(params_bomp@w))

# Plot all scaling
oxygen = seq(0.001, 4, length.out = 100)
par(mfrow = c(1,4), mar = c(2,2,2,2), oma = c(5,5,3,10))
plot(oxygen, 1/(1+exp(-bomp_scaling["cod","U_res"]*(oxygen-bomp_scaling["cod","k_res"]))), ylim = c(0,1), xlab = "", ylab = "", type = 'l', lwd = 2, cex.axis = 1.6, cex.main = 2, lty = lty[1], col = rgb(color[1,1],color[2,1],color[3,1],alphas[1],maxColorValue=255), main = "Benthos")
lines(oxygen, 1/(1+exp(-bop_scaling["cod","U_res"]*(oxygen-bop_scaling["cod","k_res"]))), lwd = 2, lty = lty[2], col = rgb(color[1,2],color[2,2],color[3,2],alphas[1],maxColorValue=255))
lines(oxygen, 1/(1+exp(-bmp_scaling["cod","U_res"]*(oxygen-bmp_scaling["cod","k_res"]))), lwd = 2, lty = lty[3], col = rgb(color[1,3],color[2,3],color[3,3],alphas[1],maxColorValue=255))
lines(oxygen, 1/(1+exp(-omp_scaling["cod","U_res"]*(oxygen-omp_scaling["cod","k_res"]))), lwd = 2, lty = lty[4], col = rgb(color[1,4],color[2,4],color[3,4],alphas[1],maxColorValue=255))
lines(oxygen, 1/(1+exp(-bp_scaling["cod","U_res"]*(oxygen-bp_scaling["cod","k_res"]))), lwd = 2, lty = lty[5], col = rgb(color[1,5],color[2,5],color[3,5],alphas[1],maxColorValue=255))
lines(oxygen, 1/(1+exp(-op_scaling["cod","U_res"]*(oxygen-op_scaling["cod","k_res"]))), lwd = 2, lty = lty[6], col = rgb(color[1,6],color[2,6],color[3,6],alphas[1],maxColorValue=255))
lines(oxygen, 1/(1+exp(-mp_scaling["cod","U_res"]*(oxygen-mp_scaling["cod","k_res"]))), lwd = 2, lty = lty[7], col = rgb(color[1,7],color[2,7],color[3,7],alphas[1],maxColorValue=255))
lines(oxygen, 1/(1+exp(-p_scaling["cod","U_res"]*(oxygen-p_scaling["cod","k_res"]))), lwd = 2, lty = lty[8], col = rgb(color[1,8],color[2,8],color[3,8],alphas[1],maxColorValue=255))
lines(oxygen, 1/(1+exp(-b_scaling["cod","U_res"]*(oxygen-b_scaling["cod","k_res"]))), lwd = 2, lty = lty[9], col = rgb(color[1,9],color[2,9],color[3,9],alphas[1],maxColorValue=255))
lines(oxygen, 1/(1+exp(-o_scaling["cod","U_res"]*(oxygen-o_scaling["cod","k_res"]))), lwd = 2, lty = lty[10], col = rgb(color[1,10],color[2,10],color[3,10],alphas[1],maxColorValue=255))
lines(oxygen, 1/(1+exp(-m_scaling["cod","U_res"]*(oxygen-m_scaling["cod","k_res"]))), lwd = 2, lty = lty[11], col = rgb(color[1,11],color[2,11],color[3,11],alphas[1],maxColorValue=255))
lines(oxygen, 1/(1+exp(-bo_scaling["cod","U_res"]*(oxygen-bo_scaling["cod","k_res"]))), lwd = 2, lty = lty[12], col = rgb(color[1,12],color[2,12],color[3,12],alphas[1],maxColorValue=255))
lines(oxygen, 1/(1+exp(-bm_scaling["cod","U_res"]*(oxygen-bm_scaling["cod","k_res"]))), lwd = 2, lty = lty[13], col = rgb(color[1,13],color[2,13],color[3,13],alphas[1],maxColorValue=255))
lines(oxygen, 1/(1+exp(-om_scaling["cod","U_res"]*(oxygen-om_scaling["cod","k_res"]))), lwd = 2, lty = lty[14], col = rgb(color[1,14],color[2,14],color[3,14],alphas[1],maxColorValue=255))
lines(oxygen, 1/(1+exp(-bom_scaling["cod","U_res"]*(oxygen-bom_scaling["cod","k_res"]))), lwd = 2, lty = lty[15], col = rgb(color[1,15],color[2,15],color[3,15],alphas[1],maxColorValue=255))
lines(oxygen, 1/(1+exp(-none_scaling["cod","U_res"]*(oxygen-none_scaling["cod","k_res"]))), lwd = 2, lty = lty[16], col = rgb(color[1,16],color[2,16],color[3,16],alphas[1],maxColorValue=255))
plot(oxygen, 1/(1+exp(-bomp_scaling["cod","U_hab"]*(oxygen-bomp_scaling["cod","a_hab"]*params_bomp@other_params$P_crit["cod",1]))), ylim = c(0,1), xlab = "", ylab = "", type = 'l', lwd = 2, cex.axis = 1.6, cex.main = 2, lty = lty[1], col = rgb(color[1,1],color[2,1],color[3,1],alphas[1],maxColorValue=255), main = "Occupancy")
for(i in pwt) lines(oxygen, 1/(1+exp(-bop_scaling["cod","U_hab"]*(oxygen-bop_scaling["cod","a_hab"]*params_bop@other_params$P_crit["cod",i]))), lwd = 2, lty = lty[2], col = rgb(color[1,2],color[2,2],color[3,2],alphas[i],maxColorValue=255))
for(i in pwt) lines(oxygen, 1/(1+exp(-bmp_scaling["cod","U_hab"]*(oxygen-bmp_scaling["cod","a_hab"]*params_bmp@other_params$P_crit["cod",i]))), lwd = 2, lty = lty[3], col = rgb(color[1,3],color[2,3],color[3,3],alphas[i],maxColorValue=255))
for(i in pwt) lines(oxygen, 1/(1+exp(-omp_scaling["cod","U_hab"]*(oxygen-omp_scaling["cod","a_hab"]*params_omp@other_params$P_crit["cod",i]))), lwd = 2, lty = lty[4], col = rgb(color[1,4],color[2,4],color[3,4],alphas[i],maxColorValue=255))
for(i in pwt) lines(oxygen, 1/(1+exp(-bp_scaling["cod","U_hab"]*(oxygen-bp_scaling["cod","a_hab"]*params_bp@other_params$P_crit["cod",i]))), lwd = 2, lty = lty[5], col = rgb(color[1,5],color[2,5],color[3,5],alphas[i],maxColorValue=255))
for(i in pwt) lines(oxygen, 1/(1+exp(-op_scaling["cod","U_hab"]*(oxygen-op_scaling["cod","a_hab"]*params_op@other_params$P_crit["cod",i]))), lwd = 2, lty = lty[6], col = rgb(color[1,6],color[2,6],color[3,6],alphas[i],maxColorValue=255))
for(i in pwt) lines(oxygen, 1/(1+exp(-mp_scaling["cod","U_hab"]*(oxygen-mp_scaling["cod","a_hab"]*params_mp@other_params$P_crit["cod",i]))), lwd = 2, lty = lty[7], col = rgb(color[1,7],color[2,7],color[3,7],alphas[i],maxColorValue=255))
for(i in pwt) lines(oxygen, 1/(1+exp(-p_scaling["cod","U_hab"]*(oxygen-p_scaling["cod","a_hab"]*params_p@other_params$P_crit["cod",i]))), lwd = 2, lty = lty[8], col = rgb(color[1,8],color[2,8],color[3,8],alphas[i],maxColorValue=255))
for(i in pwt) lines(oxygen, 1/(1+exp(-b_scaling["cod","U_hab"]*(oxygen-b_scaling["cod","a_hab"]*params_b@other_params$P_crit["cod",i]))), lwd = 2, lty = lty[9], col = rgb(color[1,9],color[2,9],color[3,9],alphas[i],maxColorValue=255))
for(i in pwt) lines(oxygen, 1/(1+exp(-o_scaling["cod","U_hab"]*(oxygen-o_scaling["cod","a_hab"]*params_o@other_params$P_crit["cod",i]))), lwd = 2, lty = lty[10], col = rgb(color[1,10],color[2,10],color[3,10],alphas[i],maxColorValue=255))
for(i in pwt) lines(oxygen, 1/(1+exp(-m_scaling["cod","U_hab"]*(oxygen-m_scaling["cod","a_hab"]*params_m@other_params$P_crit["cod",i]))), lwd = 2, lty = lty[11], col = rgb(color[1,11],color[2,11],color[3,11],alphas[i],maxColorValue=255))
for(i in pwt) lines(oxygen, 1/(1+exp(-bo_scaling["cod","U_hab"]*(oxygen-bo_scaling["cod","a_hab"]*params_bo@other_params$P_crit["cod",i]))), lwd = 2, lty = lty[12], col = rgb(color[1,12],color[2,12],color[3,12],alphas[i],maxColorValue=255))
for(i in pwt) lines(oxygen, 1/(1+exp(-bm_scaling["cod","U_hab"]*(oxygen-bm_scaling["cod","a_hab"]*params_bm@other_params$P_crit["cod",i]))), lwd = 2, lty = lty[13], col = rgb(color[1,13],color[2,13],color[3,13],alphas[i],maxColorValue=255))
for(i in pwt) lines(oxygen, 1/(1+exp(-om_scaling["cod","U_hab"]*(oxygen-om_scaling["cod","a_hab"]*params_om@other_params$P_crit["cod",i]))), lwd = 2, lty = lty[14], col = rgb(color[1,14],color[2,14],color[3,14],alphas[i],maxColorValue=255))
for(i in pwt) lines(oxygen, 1/(1+exp(-bom_scaling["cod","U_hab"]*(oxygen-bom_scaling["cod","a_hab"]*params_bom@other_params$P_crit["cod",i]))), lwd = 2, lty = lty[15], col = rgb(color[1,15],color[2,15],color[3,15],alphas[i],maxColorValue=255))
for(i in pwt) lines(oxygen, 1/(1+exp(-none_scaling["cod","U_hab"]*(oxygen-none_scaling["cod","a_hab"]*params_none@other_params$P_crit["cod",i]))), lwd = 2, lty = lty[16], col = rgb(color[1,16],color[2,16],color[3,16],alphas[i],maxColorValue=255))
plot(oxygen, 1/(1+exp(-bomp_scaling["cod","U_crit"]*(oxygen-bomp_scaling["cod","a_crit"]*params_bomp@other_params$P_crit["cod",1]))), ylim = c(0.0001,100), axes = F, xlab = "", ylab = "", type = 'l', log = 'y', lwd = 2, cex.main = 2, lty = lty[1], col = rgb(color[1,1],color[2,1],color[3,1],alphas[1],maxColorValue=255), main = "Physiology - cod")
axis(1, cex.axis = 1.6); axis(2, at = c(1e-04,1e-03,1e-02,1e-01,1e+00,1e+01,1e+02), labels = c(0.0001,0.001,0.01,0.1,1,10,100), cex.axis = 1.6); box()
for(i in pwt) lines(oxygen, 1/(1+exp(-bomp_scaling["cod","U_crit"]*(oxygen-bomp_scaling["cod","a_crit"]*params_bomp@other_params$P_crit["cod",i]))), lwd = 2, lty = lty[1], col = rgb(color[1,1],color[2,1],color[3,1],alphas[i],maxColorValue=255))
for(i in pwt) lines(oxygen, 1/(1+exp(-bop_scaling["cod","U_crit"]*(oxygen-bop_scaling["cod","a_crit"]*params_bop@other_params$P_crit["cod",i]))), lwd = 2, lty = lty[2], col = rgb(color[1,2],color[2,2],color[3,2],alphas[i],maxColorValue=255))
for(i in pwt) lines(oxygen, 1/(1+exp(-bmp_scaling["cod","U_crit"]*(oxygen-bmp_scaling["cod","a_crit"]*params_bmp@other_params$P_crit["cod",i]))), lwd = 2, lty = lty[3], col = rgb(color[1,3],color[2,3],color[3,3],alphas[i],maxColorValue=255))
for(i in pwt) lines(oxygen, 1/(1+exp(-omp_scaling["cod","U_crit"]*(oxygen-omp_scaling["cod","a_crit"]*params_omp@other_params$P_crit["cod",i]))), lwd = 2, lty = lty[4], col = rgb(color[1,4],color[2,4],color[3,4],alphas[i],maxColorValue=255))
for(i in pwt) lines(oxygen, 1/(1+exp(-bp_scaling["cod","U_crit"]*(oxygen-bp_scaling["cod","a_crit"]*params_bm@other_params$P_crit["cod",i]))), lwd = 2, lty = lty[5], col = rgb(color[1,5],color[2,5],color[3,5],alphas[i],maxColorValue=255))
for(i in pwt) lines(oxygen, 1/(1+exp(-op_scaling["cod","U_crit"]*(oxygen-op_scaling["cod","a_crit"]*params_om@other_params$P_crit["cod",i]))), lwd = 2, lty = lty[6], col = rgb(color[1,6],color[2,6],color[3,6],alphas[i],maxColorValue=255))
for(i in pwt) lines(oxygen, 1/(1+exp(-mp_scaling["cod","U_crit"]*(oxygen-mp_scaling["cod","a_crit"]*params_op@other_params$P_crit["cod",i]))), lwd = 2, lty = lty[7], col = rgb(color[1,7],color[2,7],color[3,7],alphas[i],maxColorValue=255))
for(i in pwt) lines(oxygen, 1/(1+exp(-p_scaling["cod","U_crit"]*(oxygen-p_scaling["cod","a_crit"]*params_b@other_params$P_crit["cod",i]))), lwd = 2, lty = lty[8], col = rgb(color[1,8],color[2,8],color[3,8],alphas[i],maxColorValue=255))
for(i in pwt) lines(oxygen, 1/(1+exp(-b_scaling["cod","U_crit"]*(oxygen-b_scaling["cod","a_crit"]*params_o@other_params$P_crit["cod",i]))), lwd = 2, lty = lty[9], col = rgb(color[1,9],color[2,9],color[3,9],alphas[i],maxColorValue=255))
for(i in pwt) lines(oxygen, 1/(1+exp(-o_scaling["cod","U_crit"]*(oxygen-o_scaling["cod","a_crit"]*params_m@other_params$P_crit["cod",i]))), lwd = 2, lty = lty[10], col = rgb(color[1,10],color[2,10],color[3,10],alphas[i],maxColorValue=255))
for(i in pwt) lines(oxygen, 1/(1+exp(-m_scaling["cod","U_crit"]*(oxygen-m_scaling["cod","a_crit"]*params_mp@other_params$P_crit["cod",i]))), lwd = 2, lty = lty[11], col = rgb(color[1,11],color[2,11],color[3,11],alphas[i],maxColorValue=255))
for(i in pwt) lines(oxygen, 1/(1+exp(-bo_scaling["cod","U_crit"]*(oxygen-bo_scaling["cod","a_crit"]*params_p@other_params$P_crit["cod",i]))), lwd = 2, lty = lty[12], col = rgb(color[1,12],color[2,12],color[3,12],alphas[i],maxColorValue=255))
for(i in pwt) lines(oxygen, 1/(1+exp(-bm_scaling["cod","U_crit"]*(oxygen-bm_scaling["cod","a_crit"]*params_bo@other_params$P_crit["cod",i]))), lwd = 2, lty = lty[13], col = rgb(color[1,13],color[2,13],color[3,13],alphas[i],maxColorValue=255))
for(i in pwt) lines(oxygen, 1/(1+exp(-om_scaling["cod","U_crit"]*(oxygen-om_scaling["cod","a_crit"]*params_bp@other_params$P_crit["cod",i]))), lwd = 2, lty = lty[14], col = rgb(color[1,14],color[2,14],color[3,14],alphas[i],maxColorValue=255))
for(i in pwt) lines(oxygen, 1/(1+exp(-bom_scaling["cod","U_crit"]*(oxygen-bom_scaling["cod","a_crit"]*params_bom@other_params$P_crit["cod",i]))), lwd = 2, lty = lty[15], col = rgb(color[1,15],color[2,15],color[3,15],alphas[i],maxColorValue=255))
for(i in pwt) lines(oxygen, 1/(1+exp(-none_scaling["cod","U_crit"]*(oxygen-none_scaling["cod","a_crit"]*params_none@other_params$P_crit["cod",i]))), lty = lty[16], lwd = 2, col = rgb(color[1,16],color[2,16],color[3,16],alphas[i],maxColorValue=255))
for(i in pwt) lines(oxygen, 1+exp(-bomp_scaling["cod","U_met"]*(oxygen-bomp_scaling["cod","a_met"]*params_bomp@other_params$P_crit["cod",i])), lwd = 2, lty = lty[1], col = rgb(color[1,1],color[2,1],color[3,1],alphas[i],maxColorValue=255))
for(i in pwt) lines(oxygen, 1+exp(-bop_scaling["cod","U_met"]*(oxygen-bop_scaling["cod","a_met"]*params_bop@other_params$P_crit["cod",i])), lwd = 2, lty = lty[2], col = rgb(color[1,2],color[2,2],color[3,2],alphas[i],maxColorValue=255))
for(i in pwt) lines(oxygen, 1+exp(-bmp_scaling["cod","U_met"]*(oxygen-bmp_scaling["cod","a_met"]*params_bmp@other_params$P_crit["cod",i])), lwd = 2, lty = lty[3], col = rgb(color[1,3],color[2,3],color[3,3],alphas[i],maxColorValue=255))
for(i in pwt) lines(oxygen, 1+exp(-omp_scaling["cod","U_met"]*(oxygen-omp_scaling["cod","a_met"]*params_omp@other_params$P_crit["cod",i])), lwd = 2, lty = lty[4], col = rgb(color[1,4],color[2,4],color[3,4],alphas[i],maxColorValue=255))
for(i in pwt) lines(oxygen, 1+exp(-bp_scaling["cod","U_met"]*(oxygen-bp_scaling["cod","a_met"]*params_bm@other_params$P_crit["cod",i])), lwd = 2, lty = lty[5], col = rgb(color[1,5],color[2,5],color[3,5],alphas[i],maxColorValue=255))
for(i in pwt) lines(oxygen, 1+exp(-op_scaling["cod","U_met"]*(oxygen-op_scaling["cod","a_met"]*params_om@other_params$P_crit["cod",i])), lwd = 2, lty = lty[6], col = rgb(color[1,6],color[2,6],color[3,6],alphas[i],maxColorValue=255))
for(i in pwt) lines(oxygen, 1+exp(-mp_scaling["cod","U_met"]*(oxygen-mp_scaling["cod","a_met"]*params_op@other_params$P_crit["cod",i])), lwd = 2, lty = lty[7], col = rgb(color[1,7],color[2,7],color[3,7],alphas[i],maxColorValue=255))
for(i in pwt) lines(oxygen, 1+exp(-p_scaling["cod","U_met"]*(oxygen-p_scaling["cod","a_met"]*params_b@other_params$P_crit["cod",i])), lwd = 2, lty = lty[8], col = rgb(color[1,8],color[2,8],color[3,8],alphas[i],maxColorValue=255))
for(i in pwt) lines(oxygen, 1+exp(-b_scaling["cod","U_met"]*(oxygen-b_scaling["cod","a_met"]*params_o@other_params$P_crit["cod",i])), lwd = 2, lty = lty[9], col = rgb(color[1,9],color[2,9],color[3,9],alphas[i],maxColorValue=255))
for(i in pwt) lines(oxygen, 1+exp(-o_scaling["cod","U_met"]*(oxygen-o_scaling["cod","a_met"]*params_m@other_params$P_crit["cod",i])), lwd = 2, lty = lty[10], col = rgb(color[1,10],color[2,10],color[3,10],alphas[i],maxColorValue=255))
for(i in pwt) lines(oxygen, 1+exp(-m_scaling["cod","U_met"]*(oxygen-m_scaling["cod","a_met"]*params_mp@other_params$P_crit["cod",i])), lwd = 2, lty = lty[11], col = rgb(color[1,11],color[2,11],color[3,11],alphas[i],maxColorValue=255))
for(i in pwt) lines(oxygen, 1+exp(-bo_scaling["cod","U_met"]*(oxygen-bo_scaling["cod","a_met"]*params_p@other_params$P_crit["cod",i])), lwd = 2, lty = lty[12], col = rgb(color[1,12],color[2,12],color[3,12],alphas[i],maxColorValue=255))
for(i in pwt) lines(oxygen, 1+exp(-bm_scaling["cod","U_met"]*(oxygen-bm_scaling["cod","a_met"]*params_bo@other_params$P_crit["cod",i])), lwd = 2, lty = lty[13], col = rgb(color[1,13],color[2,13],color[3,13],alphas[i],maxColorValue=255))
for(i in pwt) lines(oxygen, 1+exp(-om_scaling["cod","U_met"]*(oxygen-om_scaling["cod","a_met"]*params_bp@other_params$P_crit["cod",i])), lwd = 2, lty = lty[14], col = rgb(color[1,14],color[2,14],color[3,14],alphas[i],maxColorValue=255))
for(i in pwt) lines(oxygen, 1+exp(-bom_scaling["cod","U_met"]*(oxygen-bom_scaling["cod","a_met"]*params_bom@other_params$P_crit["cod",i])), lwd = 2, lty = lty[15], col = rgb(color[1,15],color[2,15],color[3,15],alphas[i],maxColorValue=255))
for(i in pwt) lines(oxygen, 1+exp(-none_scaling["cod","U_met"]*(oxygen-none_scaling["cod","a_met"]*params_none@other_params$P_crit["cod",i])), lwd = 2, lty = lty[16], col = rgb(color[1,16],color[2,16],color[3,16],alphas[i],maxColorValue=255))
plot(oxygen, 1/(1+exp(-bomp_scaling["flounder","U_crit"]*(oxygen-bomp_scaling["flounder","a_crit"]*params_bomp@other_params$P_crit["flounder",1]))), ylim = c(0.5,1.5), xlab = "", ylab = "", type = 'l', lwd = 2, cex.axis = 1.6, cex.main = 2, col = rgb(color[1,1],color[2,1],color[3,1],alphas[1],maxColorValue=255), main = "Physiology - flounder")
for(i in pwt) lines(oxygen, 1/(1+exp(-bomp_scaling["flounder","U_crit"]*(oxygen-bomp_scaling["flounder","a_crit"]*params_bomp@other_params$P_crit["flounder",i]))), lwd = 2, lty = lty[1], col = rgb(color[1,1],color[2,1],color[3,1],alphas[i],maxColorValue=255))
for(i in pwt) lines(oxygen, 1/(1+exp(-bop_scaling["flounder","U_crit"]*(oxygen-bop_scaling["flounder","a_crit"]*params_bop@other_params$P_crit["flounder",i]))), lwd = 2, lty = lty[2], col = rgb(color[1,2],color[2,2],color[3,2],alphas[i],maxColorValue=255))
for(i in pwt) lines(oxygen, 1/(1+exp(-bmp_scaling["flounder","U_crit"]*(oxygen-bmp_scaling["flounder","a_crit"]*params_bmp@other_params$P_crit["flounder",i]))), lwd = 2, lty = lty[3], col = rgb(color[1,3],color[2,3],color[3,3],alphas[i],maxColorValue=255))
for(i in pwt) lines(oxygen, 1/(1+exp(-omp_scaling["flounder","U_crit"]*(oxygen-omp_scaling["flounder","a_crit"]*params_omp@other_params$P_crit["flounder",i]))), lwd = 2, lty = lty[4], col = rgb(color[1,4],color[2,4],color[3,4],alphas[i],maxColorValue=255))
for(i in pwt) lines(oxygen, 1/(1+exp(-bp_scaling["flounder","U_crit"]*(oxygen-bp_scaling["flounder","a_crit"]*params_bm@other_params$P_crit["flounder",i]))), lwd = 2, lty = lty[5], col = rgb(color[1,5],color[2,5],color[3,5],alphas[i],maxColorValue=255))
for(i in pwt) lines(oxygen, 1/(1+exp(-op_scaling["flounder","U_crit"]*(oxygen-op_scaling["flounder","a_crit"]*params_om@other_params$P_crit["flounder",i]))), lwd = 2, lty = lty[6], col = rgb(color[1,6],color[2,6],color[3,6],alphas[i],maxColorValue=255))
for(i in pwt) lines(oxygen, 1/(1+exp(-mp_scaling["flounder","U_crit"]*(oxygen-mp_scaling["flounder","a_crit"]*params_op@other_params$P_crit["flounder",i]))), lwd = 2, lty = lty[7], col = rgb(color[1,7],color[2,7],color[3,7],alphas[i],maxColorValue=255))
for(i in pwt) lines(oxygen, 1/(1+exp(-p_scaling["flounder","U_crit"]*(oxygen-p_scaling["flounder","a_crit"]*params_b@other_params$P_crit["flounder",i]))), lwd = 2, lty = lty[8], col = rgb(color[1,8],color[2,8],color[3,8],alphas[i],maxColorValue=255))
for(i in pwt) lines(oxygen, 1/(1+exp(-b_scaling["flounder","U_crit"]*(oxygen-b_scaling["flounder","a_crit"]*params_o@other_params$P_crit["flounder",i]))), lwd = 2, lty = lty[9], col = rgb(color[1,9],color[2,9],color[3,9],alphas[i],maxColorValue=255))
for(i in pwt) lines(oxygen, 1/(1+exp(-o_scaling["flounder","U_crit"]*(oxygen-o_scaling["flounder","a_crit"]*params_m@other_params$P_crit["flounder",i]))), lwd = 2, lty = lty[10], col = rgb(color[1,10],color[2,10],color[3,10],alphas[i],maxColorValue=255))
for(i in pwt) lines(oxygen, 1/(1+exp(-m_scaling["flounder","U_crit"]*(oxygen-m_scaling["flounder","a_crit"]*params_mp@other_params$P_crit["flounder",i]))), lwd = 2, lty = lty[11], col = rgb(color[1,11],color[2,11],color[3,11],alphas[i],maxColorValue=255))
for(i in pwt) lines(oxygen, 1/(1+exp(-bo_scaling["flounder","U_crit"]*(oxygen-bo_scaling["flounder","a_crit"]*params_p@other_params$P_crit["flounder",i]))), lwd = 2, lty = lty[12], col = rgb(color[1,12],color[2,12],color[3,12],alphas[i],maxColorValue=255))
for(i in pwt) lines(oxygen, 1/(1+exp(-bm_scaling["flounder","U_crit"]*(oxygen-bm_scaling["flounder","a_crit"]*params_bo@other_params$P_crit["flounder",i]))), lwd = 2, lty = lty[13], col = rgb(color[1,13],color[2,13],color[3,13],alphas[i],maxColorValue=255))
for(i in pwt) lines(oxygen, 1/(1+exp(-om_scaling["flounder","U_crit"]*(oxygen-om_scaling["flounder","a_crit"]*params_bp@other_params$P_crit["flounder",i]))), lwd = 2, lty = lty[14], col = rgb(color[1,14],color[2,14],color[3,14],alphas[i],maxColorValue=255))
for(i in pwt) lines(oxygen, 1/(1+exp(-bom_scaling["flounder","U_crit"]*(oxygen-bom_scaling["flounder","a_crit"]*params_bom@other_params$P_crit["flounder",i]))), lwd = 2, lty = lty[15], col = rgb(color[1,15],color[2,15],color[3,15],alphas[i],maxColorValue=255))
for(i in pwt) lines(oxygen, 1/(1+exp(-none_scaling["flounder","U_crit"]*(oxygen-none_scaling["flounder","a_crit"]*params_none@other_params$P_crit["flounder",i]))), lwd = 2, lty = lty[16], col = rgb(color[1,16],color[2,16],color[3,16],alphas[i],maxColorValue=255))
for(i in pwt) lines(oxygen, 1+exp(-bomp_scaling["flounder","U_met"]*(oxygen-bomp_scaling["flounder","a_met"]*params_bomp@other_params$P_crit["flounder",i])), lwd = 2, lty = lty[1], col = rgb(color[1,1],color[2,1],color[3,1],alphas[i],maxColorValue=255))
for(i in pwt) lines(oxygen, 1+exp(-bop_scaling["flounder","U_met"]*(oxygen-bop_scaling["flounder","a_met"]*params_bop@other_params$P_crit["flounder",i])), lwd = 2, lty = lty[2], col = rgb(color[1,2],color[2,2],color[3,2],alphas[i],maxColorValue=255))
for(i in pwt) lines(oxygen, 1+exp(-bmp_scaling["flounder","U_met"]*(oxygen-bmp_scaling["flounder","a_met"]*params_bmp@other_params$P_crit["flounder",i])), lwd = 2, lty = lty[3], col = rgb(color[1,3],color[2,3],color[3,3],alphas[i],maxColorValue=255))
for(i in pwt) lines(oxygen, 1+exp(-omp_scaling["flounder","U_met"]*(oxygen-omp_scaling["flounder","a_met"]*params_omp@other_params$P_crit["flounder",i])), lwd = 2, lty = lty[4], col = rgb(color[1,4],color[2,4],color[3,4],alphas[i],maxColorValue=255))
for(i in pwt) lines(oxygen, 1+exp(-bp_scaling["flounder","U_met"]*(oxygen-bp_scaling["flounder","a_met"]*params_bm@other_params$P_crit["flounder",i])), lwd = 2, lty = lty[5], col = rgb(color[1,5],color[2,5],color[3,5],alphas[i],maxColorValue=255))
for(i in pwt) lines(oxygen, 1+exp(-op_scaling["flounder","U_met"]*(oxygen-op_scaling["flounder","a_met"]*params_om@other_params$P_crit["flounder",i])), lwd = 2, lty = lty[6], col = rgb(color[1,6],color[2,6],color[3,6],alphas[i],maxColorValue=255))
for(i in pwt) lines(oxygen, 1+exp(-mp_scaling["flounder","U_met"]*(oxygen-mp_scaling["flounder","a_met"]*params_op@other_params$P_crit["flounder",i])), lwd = 2, lty = lty[7], col = rgb(color[1,7],color[2,7],color[3,7],alphas[i],maxColorValue=255))
for(i in pwt) lines(oxygen, 1+exp(-p_scaling["flounder","U_met"]*(oxygen-p_scaling["flounder","a_met"]*params_b@other_params$P_crit["flounder",i])), lwd = 2, lty = lty[8], col = rgb(color[1,8],color[2,8],color[3,8],alphas[i],maxColorValue=255))
for(i in pwt) lines(oxygen, 1+exp(-b_scaling["flounder","U_met"]*(oxygen-b_scaling["flounder","a_met"]*params_o@other_params$P_crit["flounder",i])), lwd = 2, lty = lty[9], col = rgb(color[1,9],color[2,9],color[3,9],alphas[i],maxColorValue=255))
for(i in pwt) lines(oxygen, 1+exp(-o_scaling["flounder","U_met"]*(oxygen-o_scaling["flounder","a_met"]*params_m@other_params$P_crit["flounder",i])), lwd = 2, lty = lty[10], col = rgb(color[1,10],color[2,10],color[3,10],alphas[i],maxColorValue=255))
for(i in pwt) lines(oxygen, 1+exp(-m_scaling["flounder","U_met"]*(oxygen-m_scaling["flounder","a_met"]*params_mp@other_params$P_crit["flounder",i])), lwd = 2, lty = lty[11], col = rgb(color[1,11],color[2,11],color[3,11],alphas[i],maxColorValue=255))
for(i in pwt) lines(oxygen, 1+exp(-bo_scaling["flounder","U_met"]*(oxygen-bo_scaling["flounder","a_met"]*params_p@other_params$P_crit["flounder",i])), lwd = 2, lty = lty[12], col = rgb(color[1,12],color[2,12],color[3,12],alphas[i],maxColorValue=255))
for(i in pwt) lines(oxygen, 1+exp(-bm_scaling["flounder","U_met"]*(oxygen-bm_scaling["flounder","a_met"]*params_bo@other_params$P_crit["flounder",i])), lwd = 2, lty = lty[13], col = rgb(color[1,13],color[2,13],color[3,13],alphas[i],maxColorValue=255))
for(i in pwt) lines(oxygen, 1+exp(-om_scaling["flounder","U_met"]*(oxygen-om_scaling["flounder","a_met"]*params_bp@other_params$P_crit["flounder",i])), lwd = 2, lty = lty[14], col = rgb(color[1,14],color[2,14],color[3,14],alphas[i],maxColorValue=255))
for(i in pwt) lines(oxygen, 1+exp(-bom_scaling["flounder","U_met"]*(oxygen-bom_scaling["flounder","a_met"]*params_bom@other_params$P_crit["flounder",i])), lwd = 2, lty = lty[15], col = rgb(color[1,15],color[2,15],color[3,15],alphas[i],maxColorValue=255))
for(i in pwt) lines(oxygen, 1+exp(-none_scaling["flounder","U_met"]*(oxygen-none_scaling["flounder","a_met"]*params_none@other_params$P_crit["flounder",i])), lwd = 2, lty = lty[16], col = rgb(color[1,16],color[2,16],color[3,16],alphas[i],maxColorValue=255))
mtext("Oxygen", side = 1, line = 2, outer = T, cex = 2)
mtext("Scaling", side = 2, line = 2, outer = T, cex = 2)

# Add legend
par(fig = c(0,1,0,1), oma = c(0,0,0,0), mar = c(0,0,0,0), new = T)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend("right", bty = 'n', c("BOMP","BOP","BMP","OMP","BP","OP","MP","P","B","O","M","BO","BM","OM","BOM","None"), col = apply(color, 2, function(x) rgb(x[1],x[2],x[3],maxColorValue=255)), lty = lty, lwd = 2, xpd = TRUE, cex = 1.5, seg.len = 4)
par(mfrow = c(1,1))

# Finish plot
dev.off()


########## VIII. Plot growth during calibration period ##########

# Calibration period
cal_ts = seq(1991,2000)

# Colors for plotting
color = rep(c("#D55E00","#0072B2","#F0E442","#CC79A7","#E69F00","#009E73","#000000","#D4D4D4"),2)
lty = rep(c(2,1),each=8)

# Plot growth time series
jpeg("Plots/cal_growth.jpg", width = 44, height = 20, units = 'cm', res = 600)
par(mfrow = c(4,10), mar = c(0,0,5,0), oma = c(5,5,3,12))
for(i in 1:length(cal_ts)) {
	temp_bomp = myGrowthCurves(cal_bomp, (t_tune-length(cal_ts))+i, "cod", max_age = 15)
	temp_bop = myGrowthCurves(cal_bop, (t_tune-length(cal_ts))+i, "cod", max_age = 15)
	temp_bmp = myGrowthCurves(cal_bmp, (t_tune-length(cal_ts))+i, "cod", max_age = 15)
	temp_omp = myGrowthCurves(cal_omp, (t_tune-length(cal_ts))+i, "cod", max_age = 15)
	temp_bp = myGrowthCurves(cal_bp, (t_tune-length(cal_ts))+i, "cod", max_age = 15)
	temp_op = myGrowthCurves(cal_op, (t_tune-length(cal_ts))+i, "cod", max_age = 15)
	temp_mp = myGrowthCurves(cal_mp, (t_tune-length(cal_ts))+i, "cod", max_age = 15)
	temp_p = myGrowthCurves(cal_p, (t_tune-length(cal_ts))+i, "cod", max_age = 15)
	temp_b = myGrowthCurves(cal_b, (t_tune-length(cal_ts))+i, "cod", max_age = 15)
	temp_o = myGrowthCurves(cal_o, (t_tune-length(cal_ts))+i, "cod", max_age = 15)
	temp_m = myGrowthCurves(cal_m, (t_tune-length(cal_ts))+i, "cod", max_age = 15)
	temp_bo = myGrowthCurves(cal_bo, (t_tune-length(cal_ts))+i, "cod", max_age = 15)
	temp_bm = myGrowthCurves(cal_bm, (t_tune-length(cal_ts))+i, "cod", max_age = 15)
	temp_om = myGrowthCurves(cal_om, (t_tune-length(cal_ts))+i, "cod", max_age = 15)
	temp_bom = myGrowthCurves(cal_bom, (t_tune-length(cal_ts))+i, "cod", max_age = 15)
	temp_none = myGrowthCurves(cal_none, (t_tune-length(cal_ts))+i, "cod", max_age = 15)
	tempdat = cod_ivb[cod_ivb$Year == cal_ts[i], ]
	tempdat$Prop = tempdat$IndWgt/species_params(params_bomp_cal)["cod",]$w_inf
	tempdat$Age = tempdat$Age + 0.5
	plot(Prop ~ Age, tempdat, pch = 24, cex = 0.5, axes = F, xlim = c(-2,17), ylim = c(0,1.2)); axis(1,at=c(0,5,10,15)); box(); mtext(cal_ts[i], cex = 0.75)
	if(i == 1) {axis(2,at=c(0,1),labels=c(0,round(species_params(params_bomp_cal)["cod",]$w_inf)))}
	temp_bomp = temp_bomp/species_params(params_bomp_cal)["cod",]$w_inf
	lines(dimnames(temp_bomp)$Age, temp_bomp, lwd = 2, lty = lty[1], col = color[1])
	temp_bop = temp_bop/species_params(params_bop_cal)["cod",]$w_inf
	lines(dimnames(temp_bop)$Age, temp_bop, lwd = 2, lty = lty[2], col = color[2])
	temp_bmp = temp_bmp/species_params(params_bmp_cal)["cod",]$w_inf
	lines(dimnames(temp_bmp)$Age, temp_bmp, lwd = 2, lty = lty[3], col = color[3])
	temp_omp = temp_omp/species_params(params_omp_cal)["cod",]$w_inf
	lines(dimnames(temp_omp)$Age, temp_omp, lwd = 2, lty = lty[4], col = color[4])
	temp_bp = temp_bp/species_params(params_bp_cal)["cod",]$w_inf
	lines(dimnames(temp_bp)$Age, temp_bp, lwd = 2, lty = lty[5], col = color[5])
	temp_op = temp_op/species_params(params_op_cal)["cod",]$w_inf
	lines(dimnames(temp_op)$Age, temp_op, lwd = 2, lty = lty[6], col = color[6])
	temp_mp = temp_mp/species_params(params_mp_cal)["cod",]$w_inf
	lines(dimnames(temp_mp)$Age, temp_mp, lwd = 2, lty = lty[7], col = color[7])
	temp_p = temp_p/species_params(params_p_cal)["cod",]$w_inf
	lines(dimnames(temp_p)$Age, temp_p, lwd = 2, lty = lty[8], col = color[8])
	temp_b = temp_b/species_params(params_b_cal)["cod",]$w_inf
	lines(dimnames(temp_b)$Age, temp_b, lwd = 2, lty = lty[9], col = color[9])
	temp_o = temp_o/species_params(params_o_cal)["cod",]$w_inf
	lines(dimnames(temp_o)$Age, temp_o, lwd = 2, lty = lty[10], col = color[10])
	temp_m = temp_m/species_params(params_m_cal)["cod",]$w_inf
	lines(dimnames(temp_m)$Age, temp_m, lwd = 2, lty = lty[11], col = color[11])
	temp_bo = temp_bo/species_params(params_bo_cal)["cod",]$w_inf
	lines(dimnames(temp_bo)$Age, temp_bo, lwd = 2, lty = lty[12], col = color[12])
	temp_bm = temp_bm/species_params(params_bm_cal)["cod",]$w_inf
	lines(dimnames(temp_bm)$Age, temp_bm, lwd = 2, lty = lty[13], col = color[13])
	temp_om = temp_om/species_params(params_om_cal)["cod",]$w_inf
	lines(dimnames(temp_om)$Age, temp_om, lwd = 2, lty = lty[14], col = color[14])
	temp_bom = temp_bom/species_params(params_bom_cal)["cod",]$w_inf
	lines(dimnames(temp_bom)$Age, temp_bom, lwd = 2, lty = lty[15], col = color[15])
	temp_none = temp_none/species_params(params_none_cal)["cod",]$w_inf
	lines(dimnames(temp_none)$Age, temp_none, lwd = 2, lty = lty[16], col = color[16])	
	tempobs = cod_ivb_cal[i,]$a*(cod_ivb_cal[i,]$L_inf*(1-exp(-cod_ivb_cal[i,]$k*(seq(0,15,length.out=50)))))^cod_ivb_cal[i,]$b
	tempobs = tempobs/species_params(params_bomp_cal)["cod",]$w_inf
	lines(seq(0,15,length.out=50), tempobs, lwd = 2, lty = 3, col = "dodgerblue")
}
for(i in 1:length(cal_ts)) {
	temp_bomp = myGrowthCurves(cal_bomp, (t_tune-length(cal_ts))+i, "flounder", max_age = 26)
	temp_bop = myGrowthCurves(cal_bop, (t_tune-length(cal_ts))+i, "flounder", max_age = 26)
	temp_bmp = myGrowthCurves(cal_bmp, (t_tune-length(cal_ts))+i, "flounder", max_age = 26)
	temp_omp = myGrowthCurves(cal_omp, (t_tune-length(cal_ts))+i, "flounder", max_age = 26)
	temp_bp = myGrowthCurves(cal_bp, (t_tune-length(cal_ts))+i, "flounder", max_age = 26)
	temp_op = myGrowthCurves(cal_op, (t_tune-length(cal_ts))+i, "flounder", max_age = 26)
	temp_mp = myGrowthCurves(cal_mp, (t_tune-length(cal_ts))+i, "flounder", max_age = 26)
	temp_p = myGrowthCurves(cal_p, (t_tune-length(cal_ts))+i, "flounder", max_age = 26)
	temp_b = myGrowthCurves(cal_b, (t_tune-length(cal_ts))+i, "flounder", max_age = 26)
	temp_o = myGrowthCurves(cal_o, (t_tune-length(cal_ts))+i, "flounder", max_age = 26)
	temp_m = myGrowthCurves(cal_m, (t_tune-length(cal_ts))+i, "flounder", max_age = 26)
	temp_bo = myGrowthCurves(cal_bo, (t_tune-length(cal_ts))+i, "flounder", max_age = 26)
	temp_bm = myGrowthCurves(cal_bm, (t_tune-length(cal_ts))+i, "flounder", max_age = 26)
	temp_om = myGrowthCurves(cal_om, (t_tune-length(cal_ts))+i, "flounder", max_age = 26)
	temp_bom = myGrowthCurves(cal_bom, (t_tune-length(cal_ts))+i, "flounder", max_age = 26)
	temp_none = myGrowthCurves(cal_none, (t_tune-length(cal_ts))+i, "flounder", max_age = 26)
	tempdat = flounder_ivb
	tempdat$Prop = tempdat$IndWgt/species_params(params_bomp_cal)["flounder",]$w_inf
	tempdat$Age = tempdat$Age + 0.5
	plot(Prop ~ Age, tempdat, pch = 24, cex = 0.5, axes = F, xlim = c(-2,28), ylim = c(0,1)); axis(1,at=c(0,10,20)); box()
	if(i == 1) {axis(2,at=c(0,1),labels=c(0,round(species_params(params_bomp_cal)["flounder",]$w_inf)))}
	temp_bomp = temp_bomp/species_params(params_bomp_cal)["flounder",]$w_inf
	lines(dimnames(temp_bomp)$Age, temp_bomp, lwd = 2, lty = lty[1], col = color[1])
	temp_bop = temp_bop/species_params(params_bop_cal)["flounder",]$w_inf
	lines(dimnames(temp_bop)$Age, temp_bop, lwd = 2, lty = lty[2], col = color[2])
	temp_bmp = temp_bmp/species_params(params_bmp_cal)["flounder",]$w_inf
	lines(dimnames(temp_bmp)$Age, temp_bmp, lwd = 2, lty = lty[3], col = color[3])
	temp_omp = temp_omp/species_params(params_omp_cal)["flounder",]$w_inf
	lines(dimnames(temp_omp)$Age, temp_omp, lwd = 2, lty = lty[4], col = color[4])
	temp_bp = temp_bp/species_params(params_bp_cal)["flounder",]$w_inf
	lines(dimnames(temp_bp)$Age, temp_bp, lwd = 2, lty = lty[5], col = color[5])
	temp_op = temp_op/species_params(params_op_cal)["flounder",]$w_inf
	lines(dimnames(temp_op)$Age, temp_op, lwd = 2, lty = lty[6], col = color[6])
	temp_mp = temp_mp/species_params(params_mp_cal)["flounder",]$w_inf
	lines(dimnames(temp_mp)$Age, temp_mp, lwd = 2, lty = lty[7], col = color[7])
	temp_p = temp_p/species_params(params_p_cal)["flounder",]$w_inf
	lines(dimnames(temp_p)$Age, temp_p, lwd = 2, lty = lty[8], col = color[8])
	temp_b = temp_b/species_params(params_b_cal)["flounder",]$w_inf
	lines(dimnames(temp_b)$Age, temp_b, lwd = 2, lty = lty[9], col = color[9])
	temp_o = temp_o/species_params(params_o_cal)["flounder",]$w_inf
	lines(dimnames(temp_o)$Age, temp_o, lwd = 2, lty = lty[10], col = color[10])
	temp_m = temp_m/species_params(params_m_cal)["flounder",]$w_inf
	lines(dimnames(temp_m)$Age, temp_m, lwd = 2, lty = lty[11], col = color[11])
	temp_bo = temp_bo/species_params(params_bo_cal)["flounder",]$w_inf
	lines(dimnames(temp_bo)$Age, temp_bo, lwd = 2, lty = lty[12], col = color[12])
	temp_bm = temp_bm/species_params(params_bm_cal)["flounder",]$w_inf
	lines(dimnames(temp_bm)$Age, temp_bm, lwd = 2, lty = lty[13], col = color[13])
	temp_om = temp_om/species_params(params_om_cal)["flounder",]$w_inf
	lines(dimnames(temp_om)$Age, temp_om, lwd = 2, lty = lty[14], col = color[14])
	temp_bom = temp_bom/species_params(params_bom_cal)["flounder",]$w_inf
	lines(dimnames(temp_bom)$Age, temp_bom, lwd = 2, lty = lty[15], col = color[15])
	temp_none = temp_none/species_params(params_none_cal)["flounder",]$w_inf
	lines(dimnames(temp_none)$Age, temp_none, lwd = 2, lty = lty[16], col = color[16])
	tempobs = flounder_lh["a","stan"]*(flounder_lh["L_inf","stan"]*(1-exp(-flounder_lh["k","stan"]*(seq(0,26,length.out=50)))))^flounder_lh["b","stan"]
	tempobs = tempobs/species_params(params_bomp_cal)["flounder",]$w_inf
	lines(seq(0,26,length.out=50), tempobs, lwd = 2, lty = 3, col = "dodgerblue")
}
for(i in 1:length(cal_ts)) {
	temp_bomp = myGrowthCurves(cal_bomp, (t_tune-length(cal_ts))+i, "sprat", max_age = 16)
	temp_bomp = temp_bomp/species_params(params_bomp_cal)["sprat",]$w_inf
	temp_bop = myGrowthCurves(cal_bop, (t_tune-length(cal_ts))+i, "sprat", max_age = 16)
	temp_bop = temp_bop/species_params(params_bop_cal)["sprat",]$w_inf
	temp_bmp = myGrowthCurves(cal_bmp, (t_tune-length(cal_ts))+i, "sprat", max_age = 16)
	temp_bmp = temp_bmp/species_params(params_bmp_cal)["sprat",]$w_inf
	temp_omp = myGrowthCurves(cal_omp, (t_tune-length(cal_ts))+i, "sprat", max_age = 16)
	temp_omp = temp_omp/species_params(params_omp_cal)["sprat",]$w_inf
	temp_bp = myGrowthCurves(cal_bp, (t_tune-length(cal_ts))+i, "sprat", max_age = 16)
	temp_bp = temp_bp/species_params(params_bp_cal)["sprat",]$w_inf
	temp_op = myGrowthCurves(cal_op, (t_tune-length(cal_ts))+i, "sprat", max_age = 16)
	temp_op = temp_op/species_params(params_op_cal)["sprat",]$w_inf
	temp_mp = myGrowthCurves(cal_mp, (t_tune-length(cal_ts))+i, "sprat", max_age = 16)
	temp_mp = temp_mp/species_params(params_mp_cal)["sprat",]$w_inf
	temp_p = myGrowthCurves(cal_p, (t_tune-length(cal_ts))+i, "sprat", max_age = 16)
	temp_p = temp_p/species_params(params_p_cal)["sprat",]$w_inf
	temp_b = myGrowthCurves(cal_b, (t_tune-length(cal_ts))+i, "sprat", max_age = 16)
	temp_b = temp_b/species_params(params_b_cal)["sprat",]$w_inf
	temp_o = myGrowthCurves(cal_o, (t_tune-length(cal_ts))+i, "sprat", max_age = 16)
	temp_o = temp_o/species_params(params_o_cal)["sprat",]$w_inf
	temp_m = myGrowthCurves(cal_m, (t_tune-length(cal_ts))+i, "sprat", max_age = 16)
	temp_m = temp_m/species_params(params_m_cal)["sprat",]$w_inf
	temp_bo = myGrowthCurves(cal_bo, (t_tune-length(cal_ts))+i, "sprat", max_age = 16)
	temp_bo = temp_bo/species_params(params_bo_cal)["sprat",]$w_inf
	temp_bm = myGrowthCurves(cal_bm, (t_tune-length(cal_ts))+i, "sprat", max_age = 16)
	temp_bm = temp_bm/species_params(params_bm_cal)["sprat",]$w_inf
	temp_om = myGrowthCurves(cal_om, (t_tune-length(cal_ts))+i, "sprat", max_age = 16)
	temp_om = temp_om/species_params(params_om_cal)["sprat",]$w_inf
	temp_bom = myGrowthCurves(cal_bom, (t_tune-length(cal_ts))+i, "sprat", max_age = 16)
	temp_bom = temp_bom/species_params(params_bom_cal)["sprat",]$w_inf
	temp_none = myGrowthCurves(cal_none, (t_tune-length(cal_ts))+i, "sprat", max_age = 16)
	temp_none = temp_none/species_params(params_none_cal)["sprat",]$w_inf
	tempobs = sprat_ivb_cal[i,]$a*(sprat_ivb_cal[i,]$L_inf*(1-exp(-sprat_ivb_cal[i,]$k*(seq(0,16,length.out=50)))))^sprat_ivb_cal[i,]$b
	tempobs = tempobs/species_params(params_bomp_cal)["sprat",]$w_inf
	tempdat = sprat_ivb[sprat_ivb$Year == cal_ts[i], ]
	tempdat$Age = tempdat$Age + 0.5
	if(nrow(tempdat)>0) {
		tempdat$Prop = tempdat$IndWgt/species_params(params_bomp_cal)["sprat",]$w_inf
		plot(Prop ~ Age, tempdat, pch = 24, cex = 0.5, axes = F, xlim = c(-2,18), ylim = c(0,1.3)); axis(1,at=c(0,5,10,15)); box()
		lines(dimnames(temp_bomp)$Age, temp_bomp, lwd = 2, lty = lty[1], col = color[1])
		lines(dimnames(temp_bop)$Age, temp_bop, lwd = 2, lty = lty[2], col = color[2])
		lines(dimnames(temp_bmp)$Age, temp_bmp, lwd = 2, lty = lty[3], col = color[3])
		lines(dimnames(temp_omp)$Age, temp_omp, lwd = 2, lty = lty[4], col = color[4])
		lines(dimnames(temp_bp)$Age, temp_bp, lwd = 2, lty = lty[5], col = color[5])
		lines(dimnames(temp_op)$Age, temp_op, lwd = 2, lty = lty[6], col = color[6])
		lines(dimnames(temp_mp)$Age, temp_mp, lwd = 2, lty = lty[7], col = color[7])
		lines(dimnames(temp_p)$Age, temp_p, lwd = 2, lty = lty[8], col = color[8])
		lines(dimnames(temp_b)$Age, temp_b, lwd = 2, lty = lty[9], col = color[9])
		lines(dimnames(temp_o)$Age, temp_o, lwd = 2, lty = lty[10], col = color[10])
		lines(dimnames(temp_m)$Age, temp_m, lwd = 2, lty = lty[11], col = color[11])
		lines(dimnames(temp_bo)$Age, temp_bo, lwd = 2, lty = lty[12], col = color[12])
		lines(dimnames(temp_bm)$Age, temp_bm, lwd = 2, lty = lty[13], col = color[13])
		lines(dimnames(temp_om)$Age, temp_om, lwd = 2, lty = lty[14], col = color[14])
		lines(dimnames(temp_bom)$Age, temp_bom, lwd = 2, lty = lty[15], col = color[15])
		lines(dimnames(temp_none)$Age, temp_none, lwd = 2, lty = lty[16], col = color[16])
		lines(seq(0,16,length.out=50), tempobs, lwd = 2, lty = 3, col = "dodgerblue")
	} else {
		plot(seq(0,16,length.out=50), tempobs, type = 'l', lwd = 2, lty = 1, col = "red", axes = F, xlim = c(-2,18), ylim = c(0,1.3)); axis(1,at=c(0,5,10,15)); box()
		lines(dimnames(temp_bomp)$Age, temp_bomp, lwd = 2, lty = lty[1], col = color[1])
		lines(dimnames(temp_bop)$Age, temp_bop, lwd = 2, lty = lty[2], col = color[2])
		lines(dimnames(temp_bmp)$Age, temp_bmp, lwd = 2, lty = lty[3], col = color[3])
		lines(dimnames(temp_omp)$Age, temp_omp, lwd = 2, lty = lty[4], col = color[4])
		lines(dimnames(temp_bp)$Age, temp_bp, lwd = 2, lty = lty[5], col = color[5])
		lines(dimnames(temp_op)$Age, temp_op, lwd = 2, lty = lty[6], col = color[6])
		lines(dimnames(temp_mp)$Age, temp_mp, lwd = 2, lty = lty[7], col = color[7])
		lines(dimnames(temp_p)$Age, temp_p, lwd = 2, lty = lty[8], col = color[8])
		lines(dimnames(temp_b)$Age, temp_b, lwd = 2, lty = lty[9], col = color[9])
		lines(dimnames(temp_o)$Age, temp_o, lwd = 2, lty = lty[10], col = color[10])
		lines(dimnames(temp_m)$Age, temp_m, lwd = 2, lty = lty[11], col = color[11])
		lines(dimnames(temp_bo)$Age, temp_bo, lwd = 2, lty = lty[12], col = color[12])
		lines(dimnames(temp_bm)$Age, temp_bm, lwd = 2, lty = lty[13], col = color[13])
		lines(dimnames(temp_om)$Age, temp_om, lwd = 2, lty = lty[14], col = color[14])
		lines(dimnames(temp_bom)$Age, temp_bom, lwd = 2, lty = lty[15], col = color[15])
		lines(dimnames(temp_none)$Age, temp_none, lwd = 2, lty = lty[16], col = color[16])
		lines(seq(0,16,length.out=50), tempobs, lwd = 2, lty = 3, col = "dodgerblue")
	}
	if(i == 1) {axis(2,at=c(0,1),labels=c(0,round(species_params(params_bomp_cal)["sprat",]$w_inf)))}
}
for(i in 1:length(cal_ts)) {
	temp_bomp = myGrowthCurves(cal_bomp, (t_tune-length(cal_ts))+i, "herring", max_age = 13)
	temp_bomp = temp_bomp/species_params(params_bomp_cal)["herring",]$w_inf
	temp_bop = myGrowthCurves(cal_bop, (t_tune-length(cal_ts))+i, "herring", max_age = 13)
	temp_bop = temp_bop/species_params(params_bop_cal)["herring",]$w_inf
	temp_bmp = myGrowthCurves(cal_bmp, (t_tune-length(cal_ts))+i, "herring", max_age = 13)
	temp_bmp = temp_bmp/species_params(params_bmp_cal)["herring",]$w_inf
	temp_omp = myGrowthCurves(cal_omp, (t_tune-length(cal_ts))+i, "herring", max_age = 13)
	temp_omp = temp_omp/species_params(params_omp_cal)["herring",]$w_inf
	temp_bp = myGrowthCurves(cal_bp, (t_tune-length(cal_ts))+i, "herring", max_age = 13)
	temp_bp = temp_bp/species_params(params_bp_cal)["herring",]$w_inf
	temp_op = myGrowthCurves(cal_op, (t_tune-length(cal_ts))+i, "herring", max_age = 13)
	temp_op = temp_op/species_params(params_op_cal)["herring",]$w_inf
	temp_mp = myGrowthCurves(cal_mp, (t_tune-length(cal_ts))+i, "herring", max_age = 13)
	temp_mp = temp_mp/species_params(params_mp_cal)["herring",]$w_inf
	temp_p = myGrowthCurves(cal_p, (t_tune-length(cal_ts))+i, "herring", max_age = 13)
	temp_p = temp_p/species_params(params_p_cal)["herring",]$w_inf
	temp_b = myGrowthCurves(cal_b, (t_tune-length(cal_ts))+i, "herring", max_age = 13)
	temp_b = temp_b/species_params(params_b_cal)["herring",]$w_inf
	temp_o = myGrowthCurves(cal_o, (t_tune-length(cal_ts))+i, "herring", max_age = 13)
	temp_o = temp_o/species_params(params_o_cal)["herring",]$w_inf
	temp_m = myGrowthCurves(cal_m, (t_tune-length(cal_ts))+i, "herring", max_age = 13)
	temp_m = temp_m/species_params(params_m_cal)["herring",]$w_inf
	temp_bo = myGrowthCurves(cal_bo, (t_tune-length(cal_ts))+i, "herring", max_age = 13)
	temp_bo = temp_bo/species_params(params_bo_cal)["herring",]$w_inf
	temp_bm = myGrowthCurves(cal_bm, (t_tune-length(cal_ts))+i, "herring", max_age = 13)
	temp_bm = temp_bm/species_params(params_bm_cal)["herring",]$w_inf
	temp_om = myGrowthCurves(cal_om, (t_tune-length(cal_ts))+i, "herring", max_age = 13)
	temp_om = temp_om/species_params(params_om_cal)["herring",]$w_inf
	temp_bom = myGrowthCurves(cal_bom, (t_tune-length(cal_ts))+i, "herring", max_age = 13)
	temp_bom = temp_bom/species_params(params_bom_cal)["herring",]$w_inf
	temp_none = myGrowthCurves(cal_none, (t_tune-length(cal_ts))+i, "herring", max_age = 13)
	temp_none = temp_none/species_params(params_none_cal)["herring",]$w_inf
	tempobs = herring_ivb_cal[i,]$a*(herring_ivb_cal[i,]$L_inf*(1-exp(-herring_ivb_cal[i,]$k*(seq(0,13,length.out=50)))))^herring_ivb_cal[i,]$b
	tempobs = tempobs/species_params(params_bomp_cal)["herring",]$w_inf
	tempdat = herring_ivb[herring_ivb$Year == cal_ts[i], ]
	tempdat$Age = tempdat$Age + 0.5
	if(nrow(tempdat)>0) {
		tempdat$Prop = tempdat$IndWgt/species_params(params_bomp_cal)["herring",]$w_inf
		plot(Prop ~ Age, tempdat, pch = 24, cex = 0.5, axes = F, xlim = c(-2,18), ylim = c(0,1)); axis(1,at=c(0,5,10,15)); box()
		lines(dimnames(temp_bomp)$Age, temp_bomp, lwd = 2, lty = lty[1], col = color[1])
		lines(dimnames(temp_bop)$Age, temp_bop, lwd = 2, lty = lty[2], col = color[2])
		lines(dimnames(temp_bmp)$Age, temp_bmp, lwd = 2, lty = lty[3], col = color[3])
		lines(dimnames(temp_omp)$Age, temp_omp, lwd = 2, lty = lty[4], col = color[4])
		lines(dimnames(temp_bp)$Age, temp_bp, lwd = 2, lty = lty[5], col = color[5])
		lines(dimnames(temp_op)$Age, temp_op, lwd = 2, lty = lty[6], col = color[6])
		lines(dimnames(temp_mp)$Age, temp_mp, lwd = 2, lty = lty[7], col = color[7])
		lines(dimnames(temp_p)$Age, temp_p, lwd = 2, lty = lty[8], col = color[8])
		lines(dimnames(temp_b)$Age, temp_b, lwd = 2, lty = lty[9], col = color[9])
		lines(dimnames(temp_o)$Age, temp_o, lwd = 2, lty = lty[10], col = color[10])
		lines(dimnames(temp_m)$Age, temp_m, lwd = 2, lty = lty[11], col = color[11])
		lines(dimnames(temp_bo)$Age, temp_bo, lwd = 2, lty = lty[12], col = color[12])
		lines(dimnames(temp_bm)$Age, temp_bm, lwd = 2, lty = lty[13], col = color[13])
		lines(dimnames(temp_om)$Age, temp_om, lwd = 2, lty = lty[14], col = color[14])
		lines(dimnames(temp_bom)$Age, temp_bom, lwd = 2, lty = lty[15], col = color[15])
		lines(dimnames(temp_none)$Age, temp_none, lwd = 2, lty = lty[16], col = color[16])
		lines(seq(0,13,length.out=50), tempobs, lwd = 2, lty = 3, col = "dodgerblue")
	} else {
		plot(seq(0,13,length.out=50), tempobs, type = 'l', lwd = 2, lty = 1, col = "red", axes = F, xlim = c(-2,18), ylim = c(0,1)); axis(1,at=c(0,5,10,15)); box()
		lines(dimnames(temp_bomp)$Age, temp_bomp, lwd = 2, lty = lty[1], col = color[1])
		lines(dimnames(temp_bop)$Age, temp_bop, lwd = 2, lty = lty[2], col = color[2])
		lines(dimnames(temp_bmp)$Age, temp_bmp, lwd = 2, lty = lty[3], col = color[3])
		lines(dimnames(temp_omp)$Age, temp_omp, lwd = 2, lty = lty[4], col = color[4])
		lines(dimnames(temp_bp)$Age, temp_bp, lwd = 2, lty = lty[5], col = color[5])
		lines(dimnames(temp_op)$Age, temp_op, lwd = 2, lty = lty[6], col = color[6])
		lines(dimnames(temp_mp)$Age, temp_mp, lwd = 2, lty = lty[7], col = color[7])
		lines(dimnames(temp_p)$Age, temp_p, lwd = 2, lty = lty[8], col = color[8])
		lines(dimnames(temp_b)$Age, temp_b, lwd = 2, lty = lty[9], col = color[9])
		lines(dimnames(temp_o)$Age, temp_o, lwd = 2, lty = lty[10], col = color[10])
		lines(dimnames(temp_m)$Age, temp_m, lwd = 2, lty = lty[11], col = color[11])
		lines(dimnames(temp_bo)$Age, temp_bo, lwd = 2, lty = lty[12], col = color[12])
		lines(dimnames(temp_bm)$Age, temp_bm, lwd = 2, lty = lty[13], col = color[13])
		lines(dimnames(temp_om)$Age, temp_om, lwd = 2, lty = lty[14], col = color[14])
		lines(dimnames(temp_bom)$Age, temp_bom, lwd = 2, lty = lty[15], col = color[15])
		lines(dimnames(temp_none)$Age, temp_none, lwd = 2, lty = lty[16], col = color[16])
		lines(seq(0,13,length.out=50), tempobs, lwd = 2, lty = 3, col = "dodgerblue")
	}
	if(i == 1) {axis(2,at=c(0,1),labels=c(0,round(species_params(params_bomp_cal)["herring",]$w_inf)))}
}
mtext("Age", side = 1, line = 3, outer = T, cex = 1.5)
mtext("Weight (g)", side = 2, line = 3, outer = T, cex = 1.5)
par(fig = c(0,1,0,1), oma = c(0,0,0,0), mar = c(0,0,0,0), new = T)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend("right", bty = 'n', c("BOMP","BOP","BMP","OMP","BP","OP","MP","P","B","O","M","BO","BM","OM","BOM","None","Data","von B"), col = c(color, "black", "dodgerblue"), lty = c(lty,NA,3), lwd = c(rep(2,16),NA,2), pch = c(rep(NA,16),24,NA), pt.cex = c(rep(NA,16),0.5,NA), xpd = TRUE, cex = 1.5, seg.len = 4)
par(fig = c(0,1,0,1), oma = c(0,0,0,10), mar = c(5,5,3,1), new = T)
plot(0, 0, type = 'l', bty = 'n', xlab = "", ylab = "", xaxt = 'n', yaxt = 'n')
text(0, 0.95, "Cod", adj = 0.5, cex = 1.5)
text(0, seq(0.95,-0.68,length.out=4)[2], "Flounder", adj = 0.5, cex = 1.5)
text(0, seq(0.95,-0.68,length.out=4)[3], "Sprat", adj = 0.5, cex = 1.5)
text(0, -0.68, "Herring", adj = 0.5, cex = 1.5)
par(mfrow = c(1,1))
dev.off()


########## IX. Plot SSB during projection period ##########

# Analysis period
sim_ts = seq(1991,2019)

# Colors and line types for plotting
color = rep(c("#D55E00","#0072B2","#F0E442","#CC79A7","#E69F00","#009E73","#000000","#D4D4D4"),2)
lty = rep(c(2,1),each=8)

# Plot biomass time series
jpeg("Plots/sim_biomass.jpg", width = 24, height = 22, units = 'cm', res = 600)
par(mfrow = c(2,2), oma = c(3,3,3,12), mar = c(5,3,3,1))
plot(sim_ts, rowSums((sweep(sim_bomp@n[(t_tune-length(sim_ts)+1):t_tune,1,], 2, params_bomp_sim@w*params_bomp_sim@dw, "*"))[,params_bomp_sim@w >= species_params(params_bomp_sim)$w_mat[1]])/1000/1000/1000, ylim = c(0,200), cex.main = 1.5, cex.axis = 1.2, col = color[1], lty = lty[1], type = 'l', lwd = 2, xlab = "", ylab = "", main = "Cod")
polygon(x = c(1991,2000,2000,1991), y = c(0,0,200,200), border = NA, col = rgb(0.5,0.5,0.5,0.2))
points(cod_ssb$Year, cod_ssb$SSB/1000/1000/1000, pch = 17, cex = 1.5, col = "dodgerblue")
lines(sim_ts, rowSums((sweep(sim_bop@n[(t_tune-length(sim_ts)+1):t_tune,1,], 2, params_bop_sim@w*params_bop_sim@dw, "*"))[,params_bop_sim@w >= species_params(params_bop_sim)$w_mat[1]])/1000/1000/1000, lwd = 2, col = color[2], lty = lty[2])
lines(sim_ts, rowSums((sweep(sim_bmp@n[(t_tune-length(sim_ts)+1):t_tune,1,], 2, params_bmp_sim@w*params_bmp_sim@dw, "*"))[,params_bmp_sim@w >= species_params(params_bmp_sim)$w_mat[1]])/1000/1000/1000, lwd = 2, col = color[3], lty = lty[3])
lines(sim_ts, rowSums((sweep(sim_omp@n[(t_tune-length(sim_ts)+1):t_tune,1,], 2, params_omp_sim@w*params_omp_sim@dw, "*"))[,params_omp_sim@w >= species_params(params_omp_sim)$w_mat[1]])/1000/1000/1000, lwd = 2, col = color[4], lty = lty[4])
lines(sim_ts, rowSums((sweep(sim_bp@n[(t_tune-length(sim_ts)+1):t_tune,1,], 2, params_bp_sim@w*params_bp_sim@dw, "*"))[,params_bp_sim@w >= species_params(params_bp_sim)$w_mat[1]])/1000/1000/1000, lwd = 2, col = color[5], lty = lty[5])
lines(sim_ts, rowSums((sweep(sim_op@n[(t_tune-length(sim_ts)+1):t_tune,1,], 2, params_op_sim@w*params_op_sim@dw, "*"))[,params_op_sim@w >= species_params(params_op_sim)$w_mat[1]])/1000/1000/1000, lwd = 2, col = color[6], lty = lty[6])
lines(sim_ts, rowSums((sweep(sim_mp@n[(t_tune-length(sim_ts)+1):t_tune,1,], 2, params_mp_sim@w*params_mp_sim@dw, "*"))[,params_mp_sim@w >= species_params(params_mp_sim)$w_mat[1]])/1000/1000/1000, lwd = 2, col = color[7], lty = lty[7])
lines(sim_ts, rowSums((sweep(sim_p@n[(t_tune-length(sim_ts)+1):t_tune,1,], 2, params_p_sim@w*params_p_sim@dw, "*"))[,params_p_sim@w >= species_params(params_p_sim)$w_mat[1]])/1000/1000/1000, lwd = 2, col = color[8], lty = lty[8])
lines(sim_ts, rowSums((sweep(sim_b@n[(t_tune-length(sim_ts)+1):t_tune,1,], 2, params_b_sim@w*params_b_sim@dw, "*"))[,params_b_sim@w >= species_params(params_b_sim)$w_mat[1]])/1000/1000/1000, lwd = 2, col = color[9], lty = lty[9])
lines(sim_ts, rowSums((sweep(sim_o@n[(t_tune-length(sim_ts)+1):t_tune,1,], 2, params_o_sim@w*params_o_sim@dw, "*"))[,params_o_sim@w >= species_params(params_o_sim)$w_mat[1]])/1000/1000/1000, lwd = 2, col = color[10], lty = lty[10])
lines(sim_ts, rowSums((sweep(sim_m@n[(t_tune-length(sim_ts)+1):t_tune,1,], 2, params_m_sim@w*params_m_sim@dw, "*"))[,params_m_sim@w >= species_params(params_m_sim)$w_mat[1]])/1000/1000/1000, lwd = 2, col = color[11], lty = lty[11])
lines(sim_ts, rowSums((sweep(sim_bo@n[(t_tune-length(sim_ts)+1):t_tune,1,], 2, params_bo_sim@w*params_bo_sim@dw, "*"))[,params_bo_sim@w >= species_params(params_bo_sim)$w_mat[1]])/1000/1000/1000, lwd = 2, col = color[12], lty = lty[12])
lines(sim_ts, rowSums((sweep(sim_bm@n[(t_tune-length(sim_ts)+1):t_tune,1,], 2, params_bm_sim@w*params_bm_sim@dw, "*"))[,params_bm_sim@w >= species_params(params_bm_sim)$w_mat[1]])/1000/1000/1000, lwd = 2, col = color[13], lty = lty[13])
lines(sim_ts, rowSums((sweep(sim_om@n[(t_tune-length(sim_ts)+1):t_tune,1,], 2, params_om_sim@w*params_om_sim@dw, "*"))[,params_om_sim@w >= species_params(params_om_sim)$w_mat[1]])/1000/1000/1000, lwd = 2, col = color[14], lty = lty[14])
lines(sim_ts, rowSums((sweep(sim_bom@n[(t_tune-length(sim_ts)+1):t_tune,1,], 2, params_bom_sim@w*params_bom_sim@dw, "*"))[,params_bom_sim@w >= species_params(params_bom_sim)$w_mat[1]])/1000/1000/1000, lwd = 2, col = color[15], lty = lty[15])
lines(sim_ts, rowSums((sweep(sim_none@n[(t_tune-length(sim_ts)+1):t_tune,1,], 2, params_none_sim@w*params_none_sim@dw, "*"))[,params_none_sim@w >= species_params(params_none_sim)$w_mat[1]])/1000/1000/1000, lwd = 2, col = color[16], lty = lty[16])
plot(sim_ts, rowSums((sweep(sim_bomp@n[(t_tune-length(sim_ts)+1):t_tune,2,], 2, params_bomp_sim@w*params_bomp_sim@dw, "*"))[,params_bomp_sim@w >= species_params(params_bomp_sim)$w_mat[2]])/1000/1000/1000, ylim = c(0,30), cex.main = 1.5, axes = F, col = color[1], type = 'l', lty = lty[1], lwd = 2, xlab = "", ylab = "", main = "Flounder"); axis(1, cex.axis = 1.2); axis(2, at = seq(0,30,10), cex.axis = 1.2); box()
polygon(x = c(1991,2000,2000,1991), y = c(0,0,30,30), border = NA, col = rgb(0.5,0.5,0.5,0.2))
points(flounder_ssb$Year, flounder_ssb$SSB/1000/1000/1000, pch = 17, cex = 1.5, col = "dodgerblue")
lines(sim_ts, rowSums((sweep(sim_bop@n[(t_tune-length(sim_ts)+1):t_tune,2,], 2, params_bop_sim@w*params_bop_sim@dw, "*"))[,params_bop_sim@w >= species_params(params_bop_sim)$w_mat[2]])/1000/1000/1000, lwd = 2, col = color[2], lty = lty[2])
lines(sim_ts, rowSums((sweep(sim_bmp@n[(t_tune-length(sim_ts)+1):t_tune,2,], 2, params_bmp_sim@w*params_bmp_sim@dw, "*"))[,params_bmp_sim@w >= species_params(params_bmp_sim)$w_mat[2]])/1000/1000/1000, lwd = 2, col = color[3], lty = lty[3])
lines(sim_ts, rowSums((sweep(sim_omp@n[(t_tune-length(sim_ts)+1):t_tune,2,], 2, params_omp_sim@w*params_omp_sim@dw, "*"))[,params_omp_sim@w >= species_params(params_omp_sim)$w_mat[2]])/1000/1000/1000, lwd = 2, col = color[4], lty = lty[4])
lines(sim_ts, rowSums((sweep(sim_bp@n[(t_tune-length(sim_ts)+1):t_tune,2,], 2, params_bp_sim@w*params_bp_sim@dw, "*"))[,params_bp_sim@w >= species_params(params_bp_sim)$w_mat[2]])/1000/1000/1000, lwd = 2, col = color[5], lty = lty[5])
lines(sim_ts, rowSums((sweep(sim_op@n[(t_tune-length(sim_ts)+1):t_tune,2,], 2, params_op_sim@w*params_op_sim@dw, "*"))[,params_op_sim@w >= species_params(params_op_sim)$w_mat[2]])/1000/1000/1000, lwd = 2, col = color[6], lty = lty[6])
lines(sim_ts, rowSums((sweep(sim_mp@n[(t_tune-length(sim_ts)+1):t_tune,2,], 2, params_mp_sim@w*params_mp_sim@dw, "*"))[,params_mp_sim@w >= species_params(params_mp_sim)$w_mat[2]])/1000/1000/1000, lwd = 2, col = color[7], lty = lty[7])
lines(sim_ts, rowSums((sweep(sim_p@n[(t_tune-length(sim_ts)+1):t_tune,2,], 2, params_p_sim@w*params_p_sim@dw, "*"))[,params_p_sim@w >= species_params(params_p_sim)$w_mat[2]])/1000/1000/1000, lwd = 2, col = color[8], lty = lty[8])
lines(sim_ts, rowSums((sweep(sim_b@n[(t_tune-length(sim_ts)+1):t_tune,2,], 2, params_b_sim@w*params_b_sim@dw, "*"))[,params_b_sim@w >= species_params(params_b_sim)$w_mat[2]])/1000/1000/1000, lwd = 2, col = color[9], lty = lty[9])
lines(sim_ts, rowSums((sweep(sim_o@n[(t_tune-length(sim_ts)+1):t_tune,2,], 2, params_o_sim@w*params_o_sim@dw, "*"))[,params_o_sim@w >= species_params(params_o_sim)$w_mat[2]])/1000/1000/1000, lwd = 2, col = color[10], lty = lty[10])
lines(sim_ts, rowSums((sweep(sim_m@n[(t_tune-length(sim_ts)+1):t_tune,2,], 2, params_m_sim@w*params_m_sim@dw, "*"))[,params_m_sim@w >= species_params(params_m_sim)$w_mat[2]])/1000/1000/1000, lwd = 2, col = color[11], lty = lty[11])
lines(sim_ts, rowSums((sweep(sim_bo@n[(t_tune-length(sim_ts)+1):t_tune,2,], 2, params_bo_sim@w*params_bo_sim@dw, "*"))[,params_bo_sim@w >= species_params(params_bo_sim)$w_mat[2]])/1000/1000/1000, lwd = 2, col = color[12], lty = lty[12])
lines(sim_ts, rowSums((sweep(sim_bm@n[(t_tune-length(sim_ts)+1):t_tune,2,], 2, params_bm_sim@w*params_bm_sim@dw, "*"))[,params_bm_sim@w >= species_params(params_bm_sim)$w_mat[2]])/1000/1000/1000, lwd = 2, col = color[13], lty = lty[13])
lines(sim_ts, rowSums((sweep(sim_om@n[(t_tune-length(sim_ts)+1):t_tune,2,], 2, params_om_sim@w*params_om_sim@dw, "*"))[,params_om_sim@w >= species_params(params_om_sim)$w_mat[2]])/1000/1000/1000, lwd = 2, col = color[14], lty = lty[14])
lines(sim_ts, rowSums((sweep(sim_bom@n[(t_tune-length(sim_ts)+1):t_tune,2,], 2, params_bom_sim@w*params_bom_sim@dw, "*"))[,params_bom_sim@w >= species_params(params_bom_sim)$w_mat[2]])/1000/1000/1000, lwd = 2, col = color[15], lty = lty[15])
lines(sim_ts, rowSums((sweep(sim_none@n[(t_tune-length(sim_ts)+1):t_tune,2,], 2, params_none_sim@w*params_none_sim@dw, "*"))[,params_none_sim@w >= species_params(params_none_sim)$w_mat[2]])/1000/1000/1000, lwd = 2, col = color[16], lty = lty[16])
plot(sim_ts, rowSums((sweep(sim_bomp@n[(t_tune-length(sim_ts)+1):t_tune,3,], 2, params_bomp_sim@w*params_bomp_sim@dw, "*"))[,params_bomp_sim@w >= species_params(params_bomp_sim)$w_mat[3]])/1000/1000/1000, ylim = c(0,3000), cex.main = 1.5, axes = F, col = color[1], type = 'l', lty = lty[1], lwd = 2, xlab = "", ylab = "", main = "Sprat"); axis(1, cex.axis = 1.2); axis(2, at = seq(0,3000,1000), cex.axis = 1.2); box()
polygon(x = c(1991,2000,2000,1991), y = c(0,0,3000,3000), border = NA, col = rgb(0.5,0.5,0.5,0.2))
points(sprat_ssb$Year, sprat_ssb$SSB/1000/1000/1000, pch = 17, cex = 1.5, col = "dodgerblue")
lines(sim_ts, rowSums((sweep(sim_bop@n[(t_tune-length(sim_ts)+1):t_tune,3,], 2, params_bop_sim@w*params_bop_sim@dw, "*"))[,params_bop_sim@w >= species_params(params_bop_sim)$w_mat[3]])/1000/1000/1000, lwd = 2, col = color[2], lty = lty[2])
lines(sim_ts, rowSums((sweep(sim_bmp@n[(t_tune-length(sim_ts)+1):t_tune,3,], 2, params_bmp_sim@w*params_bmp_sim@dw, "*"))[,params_bmp_sim@w >= species_params(params_bmp_sim)$w_mat[3]])/1000/1000/1000, lwd = 2, col = color[3], lty = lty[3])
lines(sim_ts, rowSums((sweep(sim_omp@n[(t_tune-length(sim_ts)+1):t_tune,3,], 2, params_omp_sim@w*params_omp_sim@dw, "*"))[,params_omp_sim@w >= species_params(params_omp_sim)$w_mat[3]])/1000/1000/1000, lwd = 2, col = color[4], lty = lty[4])
lines(sim_ts, rowSums((sweep(sim_bp@n[(t_tune-length(sim_ts)+1):t_tune,3,], 2, params_bp_sim@w*params_bp_sim@dw, "*"))[,params_bp_sim@w >= species_params(params_bp_sim)$w_mat[3]])/1000/1000/1000, lwd = 2, col = color[5], lty = lty[5])
lines(sim_ts, rowSums((sweep(sim_op@n[(t_tune-length(sim_ts)+1):t_tune,3,], 2, params_op_sim@w*params_op_sim@dw, "*"))[,params_op_sim@w >= species_params(params_op_sim)$w_mat[3]])/1000/1000/1000, lwd = 2, col = color[6], lty = lty[6])
lines(sim_ts, rowSums((sweep(sim_mp@n[(t_tune-length(sim_ts)+1):t_tune,3,], 2, params_mp_sim@w*params_mp_sim@dw, "*"))[,params_mp_sim@w >= species_params(params_mp_sim)$w_mat[3]])/1000/1000/1000, lwd = 2, col = color[7], lty = lty[7])
lines(sim_ts, rowSums((sweep(sim_p@n[(t_tune-length(sim_ts)+1):t_tune,3,], 2, params_p_sim@w*params_p_sim@dw, "*"))[,params_p_sim@w >= species_params(params_p_sim)$w_mat[3]])/1000/1000/1000, lwd = 2, col = color[8], lty = lty[8])
lines(sim_ts, rowSums((sweep(sim_b@n[(t_tune-length(sim_ts)+1):t_tune,3,], 2, params_b_sim@w*params_b_sim@dw, "*"))[,params_b_sim@w >= species_params(params_b_sim)$w_mat[3]])/1000/1000/1000, lwd = 2, col = color[9], lty = lty[9])
lines(sim_ts, rowSums((sweep(sim_o@n[(t_tune-length(sim_ts)+1):t_tune,3,], 2, params_o_sim@w*params_o_sim@dw, "*"))[,params_o_sim@w >= species_params(params_o_sim)$w_mat[3]])/1000/1000/1000, lwd = 2, col = color[10], lty = lty[10])
lines(sim_ts, rowSums((sweep(sim_m@n[(t_tune-length(sim_ts)+1):t_tune,3,], 2, params_m_sim@w*params_m_sim@dw, "*"))[,params_m_sim@w >= species_params(params_m_sim)$w_mat[3]])/1000/1000/1000, lwd = 2, col = color[11], lty = lty[11])
lines(sim_ts, rowSums((sweep(sim_bo@n[(t_tune-length(sim_ts)+1):t_tune,3,], 2, params_bo_sim@w*params_bo_sim@dw, "*"))[,params_bo_sim@w >= species_params(params_bo_sim)$w_mat[3]])/1000/1000/1000, lwd = 2, col = color[12], lty = lty[12])
lines(sim_ts, rowSums((sweep(sim_bm@n[(t_tune-length(sim_ts)+1):t_tune,3,], 2, params_bm_sim@w*params_bm_sim@dw, "*"))[,params_bm_sim@w >= species_params(params_bm_sim)$w_mat[3]])/1000/1000/1000, lwd = 2, col = color[13], lty = lty[13])
lines(sim_ts, rowSums((sweep(sim_om@n[(t_tune-length(sim_ts)+1):t_tune,3,], 2, params_om_sim@w*params_om_sim@dw, "*"))[,params_om_sim@w >= species_params(params_om_sim)$w_mat[3]])/1000/1000/1000, lwd = 2, col = color[14], lty = lty[14])
lines(sim_ts, rowSums((sweep(sim_bom@n[(t_tune-length(sim_ts)+1):t_tune,3,], 2, params_bom_sim@w*params_bom_sim@dw, "*"))[,params_bom_sim@w >= species_params(params_bom_sim)$w_mat[3]])/1000/1000/1000, lwd = 2, col = color[15], lty = lty[15])
lines(sim_ts, rowSums((sweep(sim_none@n[(t_tune-length(sim_ts)+1):t_tune,3,], 2, params_none_sim@w*params_none_sim@dw, "*"))[,params_none_sim@w >= species_params(params_none_sim)$w_mat[3]])/1000/1000/1000, lwd = 2, col = color[16], lty = lty[16])
plot(sim_ts, rowSums((sweep(sim_bomp@n[(t_tune-length(sim_ts)+1):t_tune,4,], 2, params_bomp_sim@w*params_bomp_sim@dw, "*"))[,params_bomp_sim@w >= species_params(params_bomp_sim)$w_mat[4]])/1000/1000/1000, ylim = c(0,1200), cex.main = 1.5, axes = F, col = color[1], type = 'l', lty = lty[1], lwd = 2, xlab = "", ylab = "", main = "Herring"); axis(1, cex.axis = 1.2); axis(2, at = seq(0,1200,400), cex.axis = 1.2); box()
polygon(x = c(1991,2000,2000,1991), y = c(0,0,1200,1200), border = NA, col = rgb(0.5,0.5,0.5,0.2))
points(herring_ssb$Year, herring_ssb$SSB/1000/1000/1000, pch = 17, cex = 1.5, col = "dodgerblue")
lines(sim_ts, rowSums((sweep(sim_bop@n[(t_tune-length(sim_ts)+1):t_tune,4,], 2, params_bop_sim@w*params_bop_sim@dw, "*"))[,params_bop_sim@w >= species_params(params_bop_sim)$w_mat[4]])/1000/1000/1000, lwd = 2, col = color[2], lty = lty[2])
lines(sim_ts, rowSums((sweep(sim_bmp@n[(t_tune-length(sim_ts)+1):t_tune,4,], 2, params_bmp_sim@w*params_bmp_sim@dw, "*"))[,params_bmp_sim@w >= species_params(params_bmp_sim)$w_mat[4]])/1000/1000/1000, lwd = 2, col = color[3], lty = lty[3])
lines(sim_ts, rowSums((sweep(sim_omp@n[(t_tune-length(sim_ts)+1):t_tune,4,], 2, params_omp_sim@w*params_omp_sim@dw, "*"))[,params_omp_sim@w >= species_params(params_omp_sim)$w_mat[4]])/1000/1000/1000, lwd = 2, col = color[4], lty = lty[4])
lines(sim_ts, rowSums((sweep(sim_bp@n[(t_tune-length(sim_ts)+1):t_tune,4,], 2, params_bp_sim@w*params_bp_sim@dw, "*"))[,params_bp_sim@w >= species_params(params_bp_sim)$w_mat[4]])/1000/1000/1000, lwd = 2, col = color[5], lty = lty[5])
lines(sim_ts, rowSums((sweep(sim_op@n[(t_tune-length(sim_ts)+1):t_tune,4,], 2, params_op_sim@w*params_op_sim@dw, "*"))[,params_op_sim@w >= species_params(params_op_sim)$w_mat[4]])/1000/1000/1000, lwd = 2, col = color[6], lty = lty[6])
lines(sim_ts, rowSums((sweep(sim_mp@n[(t_tune-length(sim_ts)+1):t_tune,4,], 2, params_mp_sim@w*params_mp_sim@dw, "*"))[,params_mp_sim@w >= species_params(params_mp_sim)$w_mat[4]])/1000/1000/1000, lwd = 2, col = color[7], lty = lty[7])
lines(sim_ts, rowSums((sweep(sim_p@n[(t_tune-length(sim_ts)+1):t_tune,4,], 2, params_p_sim@w*params_p_sim@dw, "*"))[,params_p_sim@w >= species_params(params_p_sim)$w_mat[4]])/1000/1000/1000, lwd = 2, col = color[8], lty = lty[8])
lines(sim_ts, rowSums((sweep(sim_b@n[(t_tune-length(sim_ts)+1):t_tune,4,], 2, params_b_sim@w*params_b_sim@dw, "*"))[,params_b_sim@w >= species_params(params_b_sim)$w_mat[4]])/1000/1000/1000, lwd = 2, col = color[9], lty = lty[9])
lines(sim_ts, rowSums((sweep(sim_o@n[(t_tune-length(sim_ts)+1):t_tune,4,], 2, params_o_sim@w*params_o_sim@dw, "*"))[,params_o_sim@w >= species_params(params_o_sim)$w_mat[4]])/1000/1000/1000, lwd = 2, col = color[10], lty = lty[10])
lines(sim_ts, rowSums((sweep(sim_m@n[(t_tune-length(sim_ts)+1):t_tune,4,], 2, params_m_sim@w*params_m_sim@dw, "*"))[,params_m_sim@w >= species_params(params_m_sim)$w_mat[4]])/1000/1000/1000, lwd = 2, col = color[11], lty = lty[11])
lines(sim_ts, rowSums((sweep(sim_bo@n[(t_tune-length(sim_ts)+1):t_tune,4,], 2, params_bo_sim@w*params_bo_sim@dw, "*"))[,params_bo_sim@w >= species_params(params_bo_sim)$w_mat[4]])/1000/1000/1000, lwd = 2, col = color[12], lty = lty[12])
lines(sim_ts, rowSums((sweep(sim_bm@n[(t_tune-length(sim_ts)+1):t_tune,4,], 2, params_bm_sim@w*params_bm_sim@dw, "*"))[,params_bm_sim@w >= species_params(params_bm_sim)$w_mat[4]])/1000/1000/1000, lwd = 2, col = color[13], lty = lty[13])
lines(sim_ts, rowSums((sweep(sim_om@n[(t_tune-length(sim_ts)+1):t_tune,4,], 2, params_om_sim@w*params_om_sim@dw, "*"))[,params_om_sim@w >= species_params(params_om_sim)$w_mat[4]])/1000/1000/1000, lwd = 2, col = color[14], lty = lty[14])
lines(sim_ts, rowSums((sweep(sim_bom@n[(t_tune-length(sim_ts)+1):t_tune,4,], 2, params_bom_sim@w*params_bom_sim@dw, "*"))[,params_bom_sim@w >= species_params(params_bom_sim)$w_mat[4]])/1000/1000/1000, lwd = 2, col = color[15], lty = lty[15])
lines(sim_ts, rowSums((sweep(sim_none@n[(t_tune-length(sim_ts)+1):t_tune,4,], 2, params_none_sim@w*params_none_sim@dw, "*"))[,params_none_sim@w >= species_params(params_none_sim)$w_mat[4]])/1000/1000/1000, lwd = 2, col = color[16], lty = lty[16])
mtext("Year", 1, outer = T, cex = 1.5)
mtext("Spawning stock biomass, SSB (kt)", 2, outer = T, cex = 1.5, line = 1)
par(fig = c(0,1,0,1), oma = c(0,0,0,0), mar = c(0,0,0,0), new = T)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend("right", bty = 'n', c("BOMP","BOP","BMP","OMP","BP","OP","MP","P","B","O","M","BO","BM","OM","BOM","None","Data"), col = c(color, "dodgerblue"), lty = c(lty,NA), lwd = c(rep(2,16),NA), pch = c(rep(NA,16),17), pt.cex = c(rep(NA,16),1.5), xpd = TRUE, cex = 1.5, seg.len = 4)
par(mfrow = c(1,1))
dev.off()

# Plot biomass time series
jpeg("Plots/sim_biomass_best.jpg", width = 24, height = 22, units = 'cm', res = 600)
par(mfrow = c(2,2), oma = c(3,3,3,12), mar = c(5,3,3,1))
lty = c("dashed","dotted","longdash","dotdash","twodash","9212","solid","solid")
plot(sim_ts, rowSums((sweep(sim_bomp@n[(t_tune-length(sim_ts)+1):t_tune,1,], 2, params_bomp_sim@w*params_bomp_sim@dw, "*"))[,params_bomp_sim@w >= species_params(params_bomp_sim)$w_mat[1]])/1000/1000/1000, ylim = c(0,200), cex.main = 1.5, cex.axis = 1.2, col = color[1], lty = lty[1], type = 'n', lwd = 4, xlab = "", ylab = "", main = "Cod")
polygon(x = c(1991,2000,2000,1991), y = c(0,0,200,200), border = NA, col = rgb(0.5,0.5,0.5,0.2))
points(cod_ssb$Year, cod_ssb$SSB/1000/1000/1000, pch = 17, cex = 1.5, col = "dodgerblue")
lines(sim_ts, rowSums((sweep(sim_none@n[(t_tune-length(sim_ts)+1):t_tune,1,], 2, params_none_sim@w*params_none_sim@dw, "*"))[,params_none_sim@w >= species_params(params_none_sim)$w_mat[1]])/1000/1000/1000, lwd = 4, col = color[16], lty = lty[8])
lines(sim_ts, rowSums((sweep(sim_b@n[(t_tune-length(sim_ts)+1):t_tune,1,], 2, params_b_sim@w*params_b_sim@dw, "*"))[,params_b_sim@w >= species_params(params_b_sim)$w_mat[1]])/1000/1000/1000, lwd = 4, col = color[9], lty = lty[1])
lines(sim_ts, rowSums((sweep(sim_o@n[(t_tune-length(sim_ts)+1):t_tune,1,], 2, params_o_sim@w*params_o_sim@dw, "*"))[,params_o_sim@w >= species_params(params_o_sim)$w_mat[1]])/1000/1000/1000, lwd = 4, col = color[10], lty = lty[2])
lines(sim_ts, rowSums((sweep(sim_m@n[(t_tune-length(sim_ts)+1):t_tune,1,], 2, params_m_sim@w*params_m_sim@dw, "*"))[,params_m_sim@w >= species_params(params_m_sim)$w_mat[1]])/1000/1000/1000, lwd = 4, col = color[11], lty = lty[3])
lines(sim_ts, rowSums((sweep(sim_bo@n[(t_tune-length(sim_ts)+1):t_tune,1,], 2, params_bo_sim@w*params_bo_sim@dw, "*"))[,params_bo_sim@w >= species_params(params_bo_sim)$w_mat[1]])/1000/1000/1000, lwd = 4, col = color[12], lty = lty[4])
lines(sim_ts, rowSums((sweep(sim_bm@n[(t_tune-length(sim_ts)+1):t_tune,1,], 2, params_bm_sim@w*params_bm_sim@dw, "*"))[,params_bm_sim@w >= species_params(params_bm_sim)$w_mat[1]])/1000/1000/1000, lwd = 4, col = color[13], lty = lty[5])
lines(sim_ts, rowSums((sweep(sim_om@n[(t_tune-length(sim_ts)+1):t_tune,1,], 2, params_om_sim@w*params_om_sim@dw, "*"))[,params_om_sim@w >= species_params(params_om_sim)$w_mat[1]])/1000/1000/1000, lwd = 4, col = color[14], lty = lty[6])
lines(sim_ts, rowSums((sweep(sim_bom@n[(t_tune-length(sim_ts)+1):t_tune,1,], 2, params_bom_sim@w*params_bom_sim@dw, "*"))[,params_bom_sim@w >= species_params(params_bom_sim)$w_mat[1]])/1000/1000/1000, lwd = 4, col = color[15], lty = lty[7])
plot(sim_ts, rowSums((sweep(sim_bomp@n[(t_tune-length(sim_ts)+1):t_tune,2,], 2, params_bomp_sim@w*params_bomp_sim@dw, "*"))[,params_bomp_sim@w >= species_params(params_bomp_sim)$w_mat[2]])/1000/1000/1000, ylim = c(0,30), cex.main = 1.5, axes = F, col = color[1], type = 'n', lty = lty[1], lwd = 4, xlab = "", ylab = "", main = "Flounder"); axis(1, cex.axis = 1.2); axis(2, at = seq(0,30,10), cex.axis = 1.2); box()
polygon(x = c(1991,2000,2000,1991), y = c(0,0,30,30), border = NA, col = rgb(0.5,0.5,0.5,0.2))
points(flounder_ssb$Year, flounder_ssb$SSB/1000/1000/1000, pch = 17, cex = 1.5, col = "dodgerblue")
lines(sim_ts, rowSums((sweep(sim_none@n[(t_tune-length(sim_ts)+1):t_tune,2,], 2, params_none_sim@w*params_none_sim@dw, "*"))[,params_none_sim@w >= species_params(params_none_sim)$w_mat[2]])/1000/1000/1000, lwd = 4, col = color[16], lty = lty[8])
lines(sim_ts, rowSums((sweep(sim_b@n[(t_tune-length(sim_ts)+1):t_tune,2,], 2, params_b_sim@w*params_b_sim@dw, "*"))[,params_b_sim@w >= species_params(params_b_sim)$w_mat[2]])/1000/1000/1000, lwd = 4, col = color[9], lty = lty[1])
lines(sim_ts, rowSums((sweep(sim_o@n[(t_tune-length(sim_ts)+1):t_tune,2,], 2, params_o_sim@w*params_o_sim@dw, "*"))[,params_o_sim@w >= species_params(params_o_sim)$w_mat[2]])/1000/1000/1000, lwd = 4, col = color[10], lty = lty[2])
lines(sim_ts, rowSums((sweep(sim_m@n[(t_tune-length(sim_ts)+1):t_tune,2,], 2, params_m_sim@w*params_m_sim@dw, "*"))[,params_m_sim@w >= species_params(params_m_sim)$w_mat[2]])/1000/1000/1000, lwd = 4, col = color[11], lty = lty[3])
lines(sim_ts, rowSums((sweep(sim_bo@n[(t_tune-length(sim_ts)+1):t_tune,2,], 2, params_bo_sim@w*params_bo_sim@dw, "*"))[,params_bo_sim@w >= species_params(params_bo_sim)$w_mat[2]])/1000/1000/1000, lwd = 4, col = color[12], lty = lty[4])
lines(sim_ts, rowSums((sweep(sim_bm@n[(t_tune-length(sim_ts)+1):t_tune,2,], 2, params_bm_sim@w*params_bm_sim@dw, "*"))[,params_bm_sim@w >= species_params(params_bm_sim)$w_mat[2]])/1000/1000/1000, lwd = 4, col = color[13], lty = lty[5])
lines(sim_ts, rowSums((sweep(sim_om@n[(t_tune-length(sim_ts)+1):t_tune,2,], 2, params_om_sim@w*params_om_sim@dw, "*"))[,params_om_sim@w >= species_params(params_om_sim)$w_mat[2]])/1000/1000/1000, lwd = 4, col = color[14], lty = lty[6])
lines(sim_ts, rowSums((sweep(sim_bom@n[(t_tune-length(sim_ts)+1):t_tune,2,], 2, params_bom_sim@w*params_bom_sim@dw, "*"))[,params_bom_sim@w >= species_params(params_bom_sim)$w_mat[2]])/1000/1000/1000, lwd = 4, col = color[15], lty = lty[7])
plot(sim_ts, rowSums((sweep(sim_bomp@n[(t_tune-length(sim_ts)+1):t_tune,3,], 2, params_bomp_sim@w*params_bomp_sim@dw, "*"))[,params_bomp_sim@w >= species_params(params_bomp_sim)$w_mat[3]])/1000/1000/1000, ylim = c(0,3000), cex.main = 1.5, axes = F, col = color[1], type = 'n', lty = lty[1], lwd = 4, xlab = "", ylab = "", main = "Sprat"); axis(1, cex.axis = 1.2); axis(2, at = seq(0,3000,1000), cex.axis = 1.2); box()
polygon(x = c(1991,2000,2000,1991), y = c(0,0,3000,3000), border = NA, col = rgb(0.5,0.5,0.5,0.2))
points(sprat_ssb$Year, sprat_ssb$SSB/1000/1000/1000, pch = 17, cex = 1.5, col = "dodgerblue")
lines(sim_ts, rowSums((sweep(sim_none@n[(t_tune-length(sim_ts)+1):t_tune,3,], 2, params_none_sim@w*params_none_sim@dw, "*"))[,params_none_sim@w >= species_params(params_none_sim)$w_mat[3]])/1000/1000/1000, lwd = 4, col = color[16], lty = lty[8])
lines(sim_ts, rowSums((sweep(sim_b@n[(t_tune-length(sim_ts)+1):t_tune,3,], 2, params_b_sim@w*params_b_sim@dw, "*"))[,params_b_sim@w >= species_params(params_b_sim)$w_mat[3]])/1000/1000/1000, lwd = 4, col = color[9], lty = lty[1])
lines(sim_ts, rowSums((sweep(sim_o@n[(t_tune-length(sim_ts)+1):t_tune,3,], 2, params_o_sim@w*params_o_sim@dw, "*"))[,params_o_sim@w >= species_params(params_o_sim)$w_mat[3]])/1000/1000/1000, lwd = 4, col = color[10], lty = lty[2])
lines(sim_ts, rowSums((sweep(sim_m@n[(t_tune-length(sim_ts)+1):t_tune,3,], 2, params_m_sim@w*params_m_sim@dw, "*"))[,params_m_sim@w >= species_params(params_m_sim)$w_mat[3]])/1000/1000/1000, lwd = 4, col = color[11], lty = lty[3])
lines(sim_ts, rowSums((sweep(sim_bo@n[(t_tune-length(sim_ts)+1):t_tune,3,], 2, params_bo_sim@w*params_bo_sim@dw, "*"))[,params_bo_sim@w >= species_params(params_bo_sim)$w_mat[3]])/1000/1000/1000, lwd = 4, col = color[12], lty = lty[4])
lines(sim_ts, rowSums((sweep(sim_bm@n[(t_tune-length(sim_ts)+1):t_tune,3,], 2, params_bm_sim@w*params_bm_sim@dw, "*"))[,params_bm_sim@w >= species_params(params_bm_sim)$w_mat[3]])/1000/1000/1000, lwd = 4, col = color[13], lty = lty[5])
lines(sim_ts, rowSums((sweep(sim_om@n[(t_tune-length(sim_ts)+1):t_tune,3,], 2, params_om_sim@w*params_om_sim@dw, "*"))[,params_om_sim@w >= species_params(params_om_sim)$w_mat[3]])/1000/1000/1000, lwd = 4, col = color[14], lty = lty[6])
lines(sim_ts, rowSums((sweep(sim_bom@n[(t_tune-length(sim_ts)+1):t_tune,3,], 2, params_bom_sim@w*params_bom_sim@dw, "*"))[,params_bom_sim@w >= species_params(params_bom_sim)$w_mat[3]])/1000/1000/1000, lwd = 4, col = color[15], lty = lty[7])
plot(sim_ts, rowSums((sweep(sim_bomp@n[(t_tune-length(sim_ts)+1):t_tune,4,], 2, params_bomp_sim@w*params_bomp_sim@dw, "*"))[,params_bomp_sim@w >= species_params(params_bomp_sim)$w_mat[4]])/1000/1000/1000, ylim = c(0,1200), cex.main = 1.5, axes = F, col = color[1], type = 'n', lty = lty[1], lwd = 4, xlab = "", ylab = "", main = "Herring"); axis(1, cex.axis = 1.2); axis(2, at = seq(0,1200,400), cex.axis = 1.2); box()
polygon(x = c(1991,2000,2000,1991), y = c(0,0,1200,1200), border = NA, col = rgb(0.5,0.5,0.5,0.2))
points(herring_ssb$Year, herring_ssb$SSB/1000/1000/1000, pch = 17, cex = 1.5, col = "dodgerblue")
lines(sim_ts, rowSums((sweep(sim_none@n[(t_tune-length(sim_ts)+1):t_tune,4,], 2, params_none_sim@w*params_none_sim@dw, "*"))[,params_none_sim@w >= species_params(params_none_sim)$w_mat[4]])/1000/1000/1000, lwd = 4, col = color[16], lty = lty[8])
lines(sim_ts, rowSums((sweep(sim_b@n[(t_tune-length(sim_ts)+1):t_tune,4,], 2, params_b_sim@w*params_b_sim@dw, "*"))[,params_b_sim@w >= species_params(params_b_sim)$w_mat[4]])/1000/1000/1000, lwd = 4, col = color[9], lty = lty[1])
lines(sim_ts, rowSums((sweep(sim_o@n[(t_tune-length(sim_ts)+1):t_tune,4,], 2, params_o_sim@w*params_o_sim@dw, "*"))[,params_o_sim@w >= species_params(params_o_sim)$w_mat[4]])/1000/1000/1000, lwd = 4, col = color[10], lty = lty[2])
lines(sim_ts, rowSums((sweep(sim_m@n[(t_tune-length(sim_ts)+1):t_tune,4,], 2, params_m_sim@w*params_m_sim@dw, "*"))[,params_m_sim@w >= species_params(params_m_sim)$w_mat[4]])/1000/1000/1000, lwd = 4, col = color[11], lty = lty[3])
lines(sim_ts, rowSums((sweep(sim_bo@n[(t_tune-length(sim_ts)+1):t_tune,4,], 2, params_bo_sim@w*params_bo_sim@dw, "*"))[,params_bo_sim@w >= species_params(params_bo_sim)$w_mat[4]])/1000/1000/1000, lwd = 4, col = color[12], lty = lty[4])
lines(sim_ts, rowSums((sweep(sim_bm@n[(t_tune-length(sim_ts)+1):t_tune,4,], 2, params_bm_sim@w*params_bm_sim@dw, "*"))[,params_bm_sim@w >= species_params(params_bm_sim)$w_mat[4]])/1000/1000/1000, lwd = 4, col = color[13], lty = lty[5])
lines(sim_ts, rowSums((sweep(sim_om@n[(t_tune-length(sim_ts)+1):t_tune,4,], 2, params_om_sim@w*params_om_sim@dw, "*"))[,params_om_sim@w >= species_params(params_om_sim)$w_mat[4]])/1000/1000/1000, lwd = 4, col = color[14], lty = lty[6])
lines(sim_ts, rowSums((sweep(sim_bom@n[(t_tune-length(sim_ts)+1):t_tune,4,], 2, params_bom_sim@w*params_bom_sim@dw, "*"))[,params_bom_sim@w >= species_params(params_bom_sim)$w_mat[4]])/1000/1000/1000, lwd = 4, col = color[15], lty = lty[7])
mtext("Year", 1, outer = T, cex = 1.5)
mtext("Spawning stock biomass, SSB (kt)", 2, outer = T, cex = 1.5, line = 1)
par(fig = c(0,1,0,1), oma = c(0,0,0,0), mar = c(0,0,0,0), new = T)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend("right", bty = 'n', c("B","O","M","BO","BM","OM","BOM","None","Data"), col = c(color[9:16],"dodgerblue"), lty = c(lty,NA), lwd = c(rep(4,8),NA), pch = c(rep(NA,8),17), pt.cex = c(rep(NA,8),1.5), xpd = TRUE, cex = 1.5, seg.len = 4)
par(mfrow = c(1,1))
dev.off()


########## X. Plot yield during projection period ##########

# Colors and line types for plotting
color = rep(c("#D55E00","#0072B2","#F0E442","#CC79A7","#E69F00","#009E73","#000000","#D4D4D4"),2)
lty = rep(c(2,1),each=8)

# Plot yield
yield_bomp = getYield(sim_bomp)
yield_bom = getYield(sim_bom)
yield_bop = getYield(sim_bop)
yield_bmp = getYield(sim_bmp)
yield_omp = getYield(sim_omp)
yield_bo = getYield(sim_bo)
yield_bm = getYield(sim_bm)
yield_bp = getYield(sim_bp)
yield_om = getYield(sim_om)
yield_op = getYield(sim_op)
yield_mp = getYield(sim_mp)
yield_b = getYield(sim_b)
yield_o = getYield(sim_o)
yield_m = getYield(sim_m)
yield_p = getYield(sim_p)
yield_none = getYield(sim_none)
jpeg("Plots/sim_yield.jpg", width = 24, height = 22, units = 'cm', res = 600)
par(mfrow = c(2,2), oma = c(3,3,3,12), mar = c(5,3,3,1))
plot(sim_ts, yield_bomp[(t_tune-length(sim_ts)+1):t_tune, "cod"]/1000/1000/1000, ylim = c(0,150), cex.main = 1.5, cex.axis = 1.2, col = color[1], type = 'l', lty = lty[1], lwd = 2, xlab = "", ylab = "", main = "Cod")
polygon(x = c(1991,2000,2000,1991), y = c(0,0,150,150), border = NA, col = rgb(0.5,0.5,0.5,0.2))
lines(sim_ts, yield_bop[(t_tune-length(sim_ts)+1):t_tune, "cod"]/1000/1000/1000, lwd = 2, col = color[2], lty = lty[2])
lines(sim_ts, yield_bmp[(t_tune-length(sim_ts)+1):t_tune, "cod"]/1000/1000/1000, lwd = 2, col = color[3], lty = lty[3])
lines(sim_ts, yield_omp[(t_tune-length(sim_ts)+1):t_tune, "cod"]/1000/1000/1000, lwd = 2, col = color[4], lty = lty[4])
lines(sim_ts, yield_bp[(t_tune-length(sim_ts)+1):t_tune, "cod"]/1000/1000/1000, lwd = 2, col = color[5], lty = lty[5])
lines(sim_ts, yield_op[(t_tune-length(sim_ts)+1):t_tune, "cod"]/1000/1000/1000, lwd = 2, col = color[6], lty = lty[6])
lines(sim_ts, yield_mp[(t_tune-length(sim_ts)+1):t_tune, "cod"]/1000/1000/1000, lwd = 2, col = color[7], lty = lty[7])
lines(sim_ts, yield_p[(t_tune-length(sim_ts)+1):t_tune, "cod"]/1000/1000/1000, lwd = 2, col = color[8], lty = lty[8])
lines(sim_ts, yield_b[(t_tune-length(sim_ts)+1):t_tune, "cod"]/1000/1000/1000, lwd = 2, col = color[9], lty = lty[9])
lines(sim_ts, yield_o[(t_tune-length(sim_ts)+1):t_tune, "cod"]/1000/1000/1000, lwd = 2, col = color[10], lty = lty[10])
lines(sim_ts, yield_m[(t_tune-length(sim_ts)+1):t_tune, "cod"]/1000/1000/1000, lwd = 2, col = color[11], lty = lty[11])
lines(sim_ts, yield_bo[(t_tune-length(sim_ts)+1):t_tune, "cod"]/1000/1000/1000, lwd = 2, col = color[12], lty = lty[12])
lines(sim_ts, yield_bm[(t_tune-length(sim_ts)+1):t_tune, "cod"]/1000/1000/1000, lwd = 2, col = color[13], lty = lty[13])
lines(sim_ts, yield_om[(t_tune-length(sim_ts)+1):t_tune, "cod"]/1000/1000/1000, lwd = 2, col = color[14], lty = lty[14])
lines(sim_ts, yield_bom[(t_tune-length(sim_ts)+1):t_tune, "cod"]/1000/1000/1000, lwd = 2, col = color[15], lty = lty[15])
lines(sim_ts, yield_none[(t_tune-length(sim_ts)+1):t_tune, "cod"]/1000/1000/1000, lwd = 2, col = color[16], lty = lty[16])
points(sim_ts, cod_catch[cod_catch$Year %in% sim_ts, ]$Catch/1000, pch = 17, cex = 1.5, col = "dodgerblue")
plot(sim_ts, yield_bomp[(t_tune-length(sim_ts)+1):t_tune, "flounder"]/1000/1000/1000, ylim = c(0,30), cex.main = 1.5, axes = F, col = color[1], type = 'l', lty = lty[1], lwd = 2, xlab = "", ylab = "", main = "Flounder"); axis(1, cex.axis = 1.2); axis(2, at = seq(0,30,10), cex.axis = 1.2); box()
polygon(x = c(1991,2000,2000,1991), y = c(0,0,30,30), border = NA, col = rgb(0.5,0.5,0.5,0.2))
lines(sim_ts, yield_bop[(t_tune-length(sim_ts)+1):t_tune, "flounder"]/1000/1000/1000, lwd = 2, col = color[2], lty = lty[2])
lines(sim_ts, yield_bmp[(t_tune-length(sim_ts)+1):t_tune, "flounder"]/1000/1000/1000, lwd = 2, col = color[3], lty = lty[3])
lines(sim_ts, yield_omp[(t_tune-length(sim_ts)+1):t_tune, "flounder"]/1000/1000/1000, lwd = 2, col = color[4], lty = lty[4])
lines(sim_ts, yield_bp[(t_tune-length(sim_ts)+1):t_tune, "flounder"]/1000/1000/1000, lwd = 2, col = color[5], lty = lty[5])
lines(sim_ts, yield_om[(t_tune-length(sim_ts)+1):t_tune, "flounder"]/1000/1000/1000, lwd = 2, col = color[6], lty = lty[6])
lines(sim_ts, yield_op[(t_tune-length(sim_ts)+1):t_tune, "flounder"]/1000/1000/1000, lwd = 2, col = color[7], lty = lty[7])
lines(sim_ts, yield_mp[(t_tune-length(sim_ts)+1):t_tune, "flounder"]/1000/1000/1000, lwd = 2, col = color[8], lty = lty[8])
lines(sim_ts, yield_p[(t_tune-length(sim_ts)+1):t_tune, "flounder"]/1000/1000/1000, lwd = 2, col = color[9], lty = lty[9])
lines(sim_ts, yield_b[(t_tune-length(sim_ts)+1):t_tune, "flounder"]/1000/1000/1000, lwd = 2, col = color[10], lty = lty[10])
lines(sim_ts, yield_o[(t_tune-length(sim_ts)+1):t_tune, "flounder"]/1000/1000/1000, lwd = 2, col = color[11], lty = lty[11])
lines(sim_ts, yield_m[(t_tune-length(sim_ts)+1):t_tune, "flounder"]/1000/1000/1000, lwd = 2, col = color[12], lty = lty[12])
lines(sim_ts, yield_bo[(t_tune-length(sim_ts)+1):t_tune, "flounder"]/1000/1000/1000, lwd = 2, col = color[13], lty = lty[13])
lines(sim_ts, yield_bm[(t_tune-length(sim_ts)+1):t_tune, "flounder"]/1000/1000/1000, lwd = 2, col = color[14], lty = lty[14])
lines(sim_ts, yield_bom[(t_tune-length(sim_ts)+1):t_tune, "flounder"]/1000/1000/1000, lwd = 2, col = color[15], lty = lty[15])
lines(sim_ts, yield_none[(t_tune-length(sim_ts)+1):t_tune, "flounder"]/1000/1000/1000, lwd = 2, col = color[16], lty = lty[16])
points(sim_ts, flounder_catch[flounder_catch$Year %in% sim_ts, ]$Catch/1000, pch = 17, cex = 1.5, col = "dodgerblue")
plot(sim_ts, yield_bomp[(t_tune-length(sim_ts)+1):t_tune, "sprat"]/1000/1000/1000, ylim = c(0,900), cex.main = 1.5, axes = F, col = color[1], type = 'l', lty = lty[1], lwd = 2, xlab = "", ylab = "", main = "Sprat"); axis(1, cex.axis = 1.2); axis(2, at = seq(0,900,300), cex.axis = 1.2); box()
polygon(x = c(1991,2000,2000,1991), y = c(0,0,900,900), border = NA, col = rgb(0.5,0.5,0.5,0.2))
lines(sim_ts, yield_bop[(t_tune-length(sim_ts)+1):t_tune, "sprat"]/1000/1000/1000, lwd = 2, col = color[2], lty = lty[2])
lines(sim_ts, yield_bmp[(t_tune-length(sim_ts)+1):t_tune, "sprat"]/1000/1000/1000, lwd = 2, col = color[3], lty = lty[3])
lines(sim_ts, yield_omp[(t_tune-length(sim_ts)+1):t_tune, "sprat"]/1000/1000/1000, lwd = 2, col = color[4], lty = lty[4])
lines(sim_ts, yield_bp[(t_tune-length(sim_ts)+1):t_tune, "sprat"]/1000/1000/1000, lwd = 2, col = color[5], lty = lty[5])
lines(sim_ts, yield_op[(t_tune-length(sim_ts)+1):t_tune, "sprat"]/1000/1000/1000, lwd = 2, col = color[6], lty = lty[6])
lines(sim_ts, yield_mp[(t_tune-length(sim_ts)+1):t_tune, "sprat"]/1000/1000/1000, lwd = 2, col = color[7], lty = lty[7])
lines(sim_ts, yield_p[(t_tune-length(sim_ts)+1):t_tune, "sprat"]/1000/1000/1000, lwd = 2, col = color[8], lty = lty[8])
lines(sim_ts, yield_b[(t_tune-length(sim_ts)+1):t_tune, "sprat"]/1000/1000/1000, lwd = 2, col = color[9], lty = lty[9])
lines(sim_ts, yield_o[(t_tune-length(sim_ts)+1):t_tune, "sprat"]/1000/1000/1000, lwd = 2, col = color[10], lty = lty[10])
lines(sim_ts, yield_m[(t_tune-length(sim_ts)+1):t_tune, "sprat"]/1000/1000/1000, lwd = 2, col = color[11], lty = lty[11])
lines(sim_ts, yield_bo[(t_tune-length(sim_ts)+1):t_tune, "sprat"]/1000/1000/1000, lwd = 2, col = color[12], lty = lty[12])
lines(sim_ts, yield_bm[(t_tune-length(sim_ts)+1):t_tune, "sprat"]/1000/1000/1000, lwd = 2, col = color[13], lty = lty[13])
lines(sim_ts, yield_om[(t_tune-length(sim_ts)+1):t_tune, "sprat"]/1000/1000/1000, lwd = 2, col = color[14], lty = lty[14])
lines(sim_ts, yield_bom[(t_tune-length(sim_ts)+1):t_tune, "sprat"]/1000/1000/1000, lwd = 2, col = color[15], lty = lty[15])
lines(sim_ts, yield_none[(t_tune-length(sim_ts)+1):t_tune, "sprat"]/1000/1000/1000, lwd = 2, col = color[16], lty = lty[16])
points(sim_ts, sprat_catch[sprat_catch$Year %in% sim_ts, ]$Catch/1000, pch = 17, cex = 1.5, col = "dodgerblue")
plot(sim_ts, yield_bomp[(t_tune-length(sim_ts)+1):t_tune, "herring"]/1000/1000/1000, ylim = c(0,600), cex.main = 1.5, axes = F, col = color[1], type = 'l', lty = lty[1], lwd = 2, xlab = "", ylab = "", main = "Herring"); axis(1, cex.axis = 1.2); axis(2, at = seq(0,600,200), cex.axis = 1.2); box()
polygon(x = c(1991,2000,2000,1991), y = c(0,0,600,600), border = NA, col = rgb(0.5,0.5,0.5,0.2))
lines(sim_ts, yield_bop[(t_tune-length(sim_ts)+1):t_tune, "herring"]/1000/1000/1000, lwd = 2, col = color[2], lty = lty[2])
lines(sim_ts, yield_bmp[(t_tune-length(sim_ts)+1):t_tune, "herring"]/1000/1000/1000, lwd = 2, col = color[3], lty = lty[3])
lines(sim_ts, yield_omp[(t_tune-length(sim_ts)+1):t_tune, "herring"]/1000/1000/1000, lwd = 2, col = color[4], lty = lty[4])
lines(sim_ts, yield_bp[(t_tune-length(sim_ts)+1):t_tune, "herring"]/1000/1000/1000, lwd = 2, col = color[5], lty = lty[5])
lines(sim_ts, yield_op[(t_tune-length(sim_ts)+1):t_tune, "herring"]/1000/1000/1000, lwd = 2, col = color[6], lty = lty[6])
lines(sim_ts, yield_mp[(t_tune-length(sim_ts)+1):t_tune, "herring"]/1000/1000/1000, lwd = 2, col = color[7], lty = lty[7])
lines(sim_ts, yield_p[(t_tune-length(sim_ts)+1):t_tune, "herring"]/1000/1000/1000, lwd = 2, col = color[8], lty = lty[8])
lines(sim_ts, yield_b[(t_tune-length(sim_ts)+1):t_tune, "herring"]/1000/1000/1000, lwd = 2, col = color[9], lty = lty[9])
lines(sim_ts, yield_o[(t_tune-length(sim_ts)+1):t_tune, "herring"]/1000/1000/1000, lwd = 2, col = color[10], lty = lty[10])
lines(sim_ts, yield_m[(t_tune-length(sim_ts)+1):t_tune, "herring"]/1000/1000/1000, lwd = 2, col = color[11], lty = lty[11])
lines(sim_ts, yield_bo[(t_tune-length(sim_ts)+1):t_tune, "herring"]/1000/1000/1000, lwd = 2, col = color[12], lty = lty[12])
lines(sim_ts, yield_bm[(t_tune-length(sim_ts)+1):t_tune, "herring"]/1000/1000/1000, lwd = 2, col = color[13], lty = lty[13])
lines(sim_ts, yield_om[(t_tune-length(sim_ts)+1):t_tune, "herring"]/1000/1000/1000, lwd = 2, col = color[14], lty = lty[14])
lines(sim_ts, yield_bom[(t_tune-length(sim_ts)+1):t_tune, "herring"]/1000/1000/1000, lwd = 2, col = color[15], lty = lty[15])
lines(sim_ts, yield_none[(t_tune-length(sim_ts)+1):t_tune, "herring"]/1000/1000/1000, lwd = 2, col = color[16], lty = lty[16])
points(sim_ts, herring_catch[herring_catch$Year %in% sim_ts, ]$Catch/1000, pch = 17, cex = 1.5, col = "dodgerblue")
mtext("Year", 1, outer = T, cex = 1.5)
mtext("Yield (kt)", 2, outer = T, cex = 1.5, line = 1)
par(fig = c(0,1,0,1), oma = c(0,0,0,0), mar = c(0,0,0,0), new = T)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend("right", bty = 'n', c("BOMP","BOP","BMP","OMP","BP","OP","MP","P","B","O","M","BO","BM","OM","BOM","None","Data"), col = c(color, "dodgerblue"), lty = c(lty,NA), lwd = c(rep(2,16),NA), pch = c(rep(NA,16),17), pt.cex = c(rep(NA,16),1.5), xpd = TRUE, cex = 1.5, seg.len = 4)
par(mfrow = c(1,1))
dev.off()

jpeg("Plots/sim_yield_best.jpg", width = 24, height = 22, units = 'cm', res = 600)
par(mfrow = c(2,2), oma = c(3,3,3,12), mar = c(5,3,3,1))
lty = c("dashed","dotted","longdash","dotdash","twodash","9212","solid","solid")
plot(sim_ts, yield_bomp[(t_tune-length(sim_ts)+1):t_tune, "cod"]/1000/1000/1000, ylim = c(0,150), cex.main = 1.5, cex.axis = 1.2, col = color[1], type = 'n', lty = lty[1], lwd = 4, xlab = "", ylab = "", main = "Cod")
polygon(x = c(1991,2000,2000,1991), y = c(0,0,150,150), border = NA, col = rgb(0.5,0.5,0.5,0.2))
points(sim_ts, cod_catch[cod_catch$Year %in% sim_ts, ]$Catch/1000, pch = 17, cex = 1.5, col = "dodgerblue")
lines(sim_ts, yield_none[(t_tune-length(sim_ts)+1):t_tune, "cod"]/1000/1000/1000, lwd = 4, col = color[16], lty = lty[8])
lines(sim_ts, yield_b[(t_tune-length(sim_ts)+1):t_tune, "cod"]/1000/1000/1000, lwd = 4, col = color[9], lty = lty[1])
lines(sim_ts, yield_o[(t_tune-length(sim_ts)+1):t_tune, "cod"]/1000/1000/1000, lwd = 4, col = color[10], lty = lty[2])
lines(sim_ts, yield_m[(t_tune-length(sim_ts)+1):t_tune, "cod"]/1000/1000/1000, lwd = 4, col = color[11], lty = lty[3])
lines(sim_ts, yield_bo[(t_tune-length(sim_ts)+1):t_tune, "cod"]/1000/1000/1000, lwd = 4, col = color[12], lty = lty[4])
lines(sim_ts, yield_bm[(t_tune-length(sim_ts)+1):t_tune, "cod"]/1000/1000/1000, lwd = 4, col = color[13], lty = lty[5])
lines(sim_ts, yield_om[(t_tune-length(sim_ts)+1):t_tune, "cod"]/1000/1000/1000, lwd = 4, col = color[14], lty = lty[6])
lines(sim_ts, yield_bom[(t_tune-length(sim_ts)+1):t_tune, "cod"]/1000/1000/1000, lwd = 4, col = color[15], lty = lty[7])
plot(sim_ts, yield_bomp[(t_tune-length(sim_ts)+1):t_tune, "flounder"]/1000/1000/1000, ylim = c(0,30), cex.main = 1.5, axes = F, col = color[1], type = 'n', lty = lty[1], lwd = 4, xlab = "", ylab = "", main = "Flounder"); axis(1, cex.axis = 1.2); axis(2, at = seq(0,30,10), cex.axis = 1.2); box()
polygon(x = c(1991,2000,2000,1991), y = c(0,0,30,30), border = NA, col = rgb(0.5,0.5,0.5,0.2))
points(sim_ts, flounder_catch[flounder_catch$Year %in% sim_ts, ]$Catch/1000, pch = 17, cex = 1.5, col = "dodgerblue")
lines(sim_ts, yield_none[(t_tune-length(sim_ts)+1):t_tune, "flounder"]/1000/1000/1000, lwd = 4, col = color[16], lty = lty[8])
lines(sim_ts, yield_b[(t_tune-length(sim_ts)+1):t_tune, "flounder"]/1000/1000/1000, lwd = 4, col = color[9], lty = lty[1])
lines(sim_ts, yield_o[(t_tune-length(sim_ts)+1):t_tune, "flounder"]/1000/1000/1000, lwd = 4, col = color[10], lty = lty[2])
lines(sim_ts, yield_m[(t_tune-length(sim_ts)+1):t_tune, "flounder"]/1000/1000/1000, lwd = 4, col = color[11], lty = lty[3])
lines(sim_ts, yield_bo[(t_tune-length(sim_ts)+1):t_tune, "flounder"]/1000/1000/1000, lwd = 4, col = color[12], lty = lty[4])
lines(sim_ts, yield_bm[(t_tune-length(sim_ts)+1):t_tune, "flounder"]/1000/1000/1000, lwd = 4, col = color[13], lty = lty[5])
lines(sim_ts, yield_om[(t_tune-length(sim_ts)+1):t_tune, "flounder"]/1000/1000/1000, lwd = 4, col = color[14], lty = lty[6])
lines(sim_ts, yield_bom[(t_tune-length(sim_ts)+1):t_tune, "flounder"]/1000/1000/1000, lwd = 4, col = color[15], lty = lty[7])
plot(sim_ts, yield_bomp[(t_tune-length(sim_ts)+1):t_tune, "sprat"]/1000/1000/1000, ylim = c(0,900), cex.main = 1.5, axes = F, col = color[1], type = 'n', lty = lty[1], lwd = 4, xlab = "", ylab = "", main = "Sprat"); axis(1, cex.axis = 1.2); axis(2, at = seq(0,900,300), cex.axis = 1.2); box()
polygon(x = c(1991,2000,2000,1991), y = c(0,0,900,900), border = NA, col = rgb(0.5,0.5,0.5,0.2))
points(sim_ts, sprat_catch[sprat_catch$Year %in% sim_ts, ]$Catch/1000, pch = 17, cex = 1.5, col = "dodgerblue")
lines(sim_ts, yield_none[(t_tune-length(sim_ts)+1):t_tune, "sprat"]/1000/1000/1000, lwd = 4, col = color[16], lty = lty[8])
lines(sim_ts, yield_b[(t_tune-length(sim_ts)+1):t_tune, "sprat"]/1000/1000/1000, lwd = 4, col = color[9], lty = lty[1])
lines(sim_ts, yield_o[(t_tune-length(sim_ts)+1):t_tune, "sprat"]/1000/1000/1000, lwd = 4, col = color[10], lty = lty[2])
lines(sim_ts, yield_m[(t_tune-length(sim_ts)+1):t_tune, "sprat"]/1000/1000/1000, lwd = 4, col = color[11], lty = lty[3])
lines(sim_ts, yield_bo[(t_tune-length(sim_ts)+1):t_tune, "sprat"]/1000/1000/1000, lwd = 4, col = color[12], lty = lty[4])
lines(sim_ts, yield_bm[(t_tune-length(sim_ts)+1):t_tune, "sprat"]/1000/1000/1000, lwd = 4, col = color[13], lty = lty[5])
lines(sim_ts, yield_om[(t_tune-length(sim_ts)+1):t_tune, "sprat"]/1000/1000/1000, lwd = 4, col = color[14], lty = lty[6])
lines(sim_ts, yield_bom[(t_tune-length(sim_ts)+1):t_tune, "sprat"]/1000/1000/1000, lwd = 4, col = color[15], lty = lty[7])
plot(sim_ts, yield_bomp[(t_tune-length(sim_ts)+1):t_tune, "herring"]/1000/1000/1000, ylim = c(0,600), cex.main = 1.5, axes = F, col = color[1], type = 'n', lty = lty[1], lwd = 4, xlab = "", ylab = "", main = "Herring"); axis(1, cex.axis = 1.2); axis(2, at = seq(0,600,200), cex.axis = 1.2); box()
polygon(x = c(1991,2000,2000,1991), y = c(0,0,600,600), border = NA, col = rgb(0.5,0.5,0.5,0.2))
points(sim_ts, herring_catch[herring_catch$Year %in% sim_ts, ]$Catch/1000, pch = 17, cex = 1.5, col = "dodgerblue")
lines(sim_ts, yield_none[(t_tune-length(sim_ts)+1):t_tune, "herring"]/1000/1000/1000, lwd = 4, col = color[16], lty = lty[8])
lines(sim_ts, yield_b[(t_tune-length(sim_ts)+1):t_tune, "herring"]/1000/1000/1000, lwd = 4, col = color[9], lty = lty[1])
lines(sim_ts, yield_o[(t_tune-length(sim_ts)+1):t_tune, "herring"]/1000/1000/1000, lwd = 4, col = color[10], lty = lty[2])
lines(sim_ts, yield_m[(t_tune-length(sim_ts)+1):t_tune, "herring"]/1000/1000/1000, lwd = 4, col = color[11], lty = lty[3])
lines(sim_ts, yield_bo[(t_tune-length(sim_ts)+1):t_tune, "herring"]/1000/1000/1000, lwd = 4, col = color[12], lty = lty[4])
lines(sim_ts, yield_bm[(t_tune-length(sim_ts)+1):t_tune, "herring"]/1000/1000/1000, lwd = 4, col = color[13], lty = lty[5])
lines(sim_ts, yield_om[(t_tune-length(sim_ts)+1):t_tune, "herring"]/1000/1000/1000, lwd = 4, col = color[14], lty = lty[6])
lines(sim_ts, yield_bom[(t_tune-length(sim_ts)+1):t_tune, "herring"]/1000/1000/1000, lwd = 4, col = color[15], lty = lty[7])
mtext("Year", 1, outer = T, cex = 1.5)
mtext("Yield (kt)", 2, outer = T, cex = 1.5, line = 1)
par(fig = c(0,1,0,1), oma = c(0,0,0,0), mar = c(0,0,0,0), new = T)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend("right", bty = 'n', c("B","O","M","BO","BM","OM","BOM","None","Data"), col = c(color[9:16], "dodgerblue"), lty = c(lty,NA), lwd = c(rep(4,8),NA), pch = c(rep(NA,8),17), pt.cex = c(rep(NA,8),1.5), xpd = TRUE, cex = 1.5, seg.len = 4)
par(mfrow = c(1,1))
dev.off()


########## XI. Plot growth during projection period ##########

# Colors and line types for plotting
color = rep(c("#D55E00","#0072B2","#F0E442","#CC79A7","#E69F00","#009E73","#000000","#D4D4D4"),2)
lty = rep(c(2,1),each=8)

# Plot growth time series
prj_ts = seq(2001,2019)
jpeg("Plots/sim_growth.jpg", width = 44, height = 20, units = 'cm', res = 600)
par(mfrow = c(4,19), mar = c(0,0,5,0), oma = c(5,5,3,12))
for(i in 1:length(prj_ts)) {
	temp_bomp = myGrowthCurves(sim_bomp, (t_tune-length(prj_ts))+i, "cod", max_age = 15)
	temp_bop = myGrowthCurves(sim_bop, (t_tune-length(prj_ts))+i, "cod", max_age = 15)
	temp_bmp = myGrowthCurves(sim_bmp, (t_tune-length(prj_ts))+i, "cod", max_age = 15)
	temp_omp = myGrowthCurves(sim_omp, (t_tune-length(prj_ts))+i, "cod", max_age = 15)
	temp_bp = myGrowthCurves(sim_bp, (t_tune-length(prj_ts))+i, "cod", max_age = 15)
	temp_op = myGrowthCurves(sim_op, (t_tune-length(prj_ts))+i, "cod", max_age = 15)
	temp_mp = myGrowthCurves(sim_mp, (t_tune-length(prj_ts))+i, "cod", max_age = 15)
	temp_p = myGrowthCurves(sim_p, (t_tune-length(prj_ts))+i, "cod", max_age = 15)
	temp_b = myGrowthCurves(sim_b, (t_tune-length(prj_ts))+i, "cod", max_age = 15)
	temp_o = myGrowthCurves(sim_o, (t_tune-length(prj_ts))+i, "cod", max_age = 15)
	temp_m = myGrowthCurves(sim_m, (t_tune-length(prj_ts))+i, "cod", max_age = 15)
	temp_bo = myGrowthCurves(sim_bo, (t_tune-length(prj_ts))+i, "cod", max_age = 15)
	temp_bm = myGrowthCurves(sim_bm, (t_tune-length(prj_ts))+i, "cod", max_age = 15)
	temp_om = myGrowthCurves(sim_om, (t_tune-length(prj_ts))+i, "cod", max_age = 15)
	temp_bom = myGrowthCurves(sim_bom, (t_tune-length(prj_ts))+i, "cod", max_age = 15)
	temp_none = myGrowthCurves(sim_none, (t_tune-length(prj_ts))+i, "cod", max_age = 15)
	tempdat = cod_ivb_sim[cod_ivb_sim$Year == prj_ts[i], ]
	tempdat$Prop = tempdat$IndWgt/species_params(params_bomp_sim)["cod",]$w_inf
	tempdat$Age = tempdat$Age + 0.5
	plot(Prop ~ Age, tempdat, axes = F, pch = 24, cex = 0.5, xlim = c(-2,17), ylim = c(0,1.2)); axis(1,at=c(0,5,10,15)); box(); mtext(prj_ts[i], cex = 0.75)
	if(i == 1) {axis(2,at=c(0,1),labels=c(0,round(species_params(params_bomp_sim)["cod",]$w_inf)))}
	temp_bomp = temp_bomp/species_params(params_bomp_sim)["cod",]$w_inf
	lines(dimnames(temp_bomp)$Age, temp_bomp, lwd = 2, lty = lty[1], col = color[1])
	temp_bop = temp_bop/species_params(params_bop_sim)["cod",]$w_inf
	lines(dimnames(temp_bop)$Age, temp_bop, lwd = 2, lty = lty[2], col = color[2])
	temp_bmp = temp_bmp/species_params(params_bmp_sim)["cod",]$w_inf
	lines(dimnames(temp_bmp)$Age, temp_bmp, lwd = 2, lty = lty[3], col = color[3])
	temp_omp = temp_omp/species_params(params_omp_sim)["cod",]$w_inf
	lines(dimnames(temp_omp)$Age, temp_omp, lwd = 2, lty = lty[4], col = color[4])
	temp_bp = temp_bp/species_params(params_bp_sim)["cod",]$w_inf
	lines(dimnames(temp_bp)$Age, temp_bp, lwd = 2, lty = lty[5], col = color[5])
	temp_op = temp_op/species_params(params_op_sim)["cod",]$w_inf
	lines(dimnames(temp_op)$Age, temp_op, lwd = 2, lty = lty[6], col = color[6])
	temp_mp = temp_mp/species_params(params_mp_sim)["cod",]$w_inf
	lines(dimnames(temp_mp)$Age, temp_mp, lwd = 2, lty = lty[7], col = color[7])
	temp_p = temp_p/species_params(params_p_sim)["cod",]$w_inf
	lines(dimnames(temp_p)$Age, temp_p, lwd = 2, lty = lty[8], col = color[8])
	temp_b = temp_b/species_params(params_b_sim)["cod",]$w_inf
	lines(dimnames(temp_b)$Age, temp_b, lwd = 2, lty = lty[9], col = color[9])
	temp_o = temp_o/species_params(params_o_sim)["cod",]$w_inf
	lines(dimnames(temp_o)$Age, temp_o, lwd = 2, lty = lty[10], col = color[10])
	temp_m = temp_m/species_params(params_m_sim)["cod",]$w_inf
	lines(dimnames(temp_m)$Age, temp_m, lwd = 2, lty = lty[11], col = color[11])
	temp_bo = temp_bo/species_params(params_bo_sim)["cod",]$w_inf
	lines(dimnames(temp_bo)$Age, temp_bo, lwd = 2, lty = lty[12], col = color[12])
	temp_bm = temp_bm/species_params(params_bm_sim)["cod",]$w_inf
	lines(dimnames(temp_bm)$Age, temp_bm, lwd = 2, lty = lty[13], col = color[13])
	temp_om = temp_om/species_params(params_om_sim)["cod",]$w_inf
	lines(dimnames(temp_om)$Age, temp_om, lwd = 2, lty = lty[14], col = color[14])
	temp_bom = temp_bom/species_params(params_bom_sim)["cod",]$w_inf
	lines(dimnames(temp_bom)$Age, temp_bom, lwd = 2, lty = lty[15], col = color[15])
	temp_none = temp_none/species_params(params_none_sim)["cod",]$w_inf
	lines(dimnames(temp_none)$Age, temp_none, lwd = 2, lty = lty[16], col = color[16])	
	tempobs = cod_ivb_prj[i,]$a*(cod_ivb_prj[i,]$L_inf*(1-exp(-cod_ivb_prj[i,]$k*(seq(0,15,length.out=50)))))^cod_ivb_prj[i,]$b
	tempobs = tempobs/species_params(params_bomp_sim)["cod",]$w_inf
	lines(seq(0,15,length.out=50), tempobs, lwd = 2, lty = 3, col = "dodgerblue")
}
for(i in 1:length(prj_ts)) {
	temp_bomp = myGrowthCurves(sim_bomp, (t_tune-length(prj_ts))+i, "flounder", max_age = 26)
	temp_bomp = temp_bomp/species_params(params_bomp_sim)["flounder",]$w_inf
	temp_bop = myGrowthCurves(sim_bop, (t_tune-length(prj_ts))+i, "flounder", max_age = 26)
	temp_bop = temp_bop/species_params(params_bop_sim)["flounder",]$w_inf
	temp_bmp = myGrowthCurves(sim_bmp, (t_tune-length(prj_ts))+i, "flounder", max_age = 26)
	temp_bmp = temp_bmp/species_params(params_bmp_sim)["flounder",]$w_inf
	temp_omp = myGrowthCurves(sim_omp, (t_tune-length(prj_ts))+i, "flounder", max_age = 26)
	temp_omp = temp_omp/species_params(params_omp_sim)["flounder",]$w_inf
	temp_bp = myGrowthCurves(sim_bp, (t_tune-length(prj_ts))+i, "flounder", max_age = 26)
	temp_bp = temp_bp/species_params(params_bp_sim)["flounder",]$w_inf
	temp_op = myGrowthCurves(sim_op, (t_tune-length(prj_ts))+i, "flounder", max_age = 26)
	temp_op = temp_op/species_params(params_op_sim)["flounder",]$w_inf
	temp_mp = myGrowthCurves(sim_mp, (t_tune-length(prj_ts))+i, "flounder", max_age = 26)
	temp_mp = temp_mp/species_params(params_mp_sim)["flounder",]$w_inf
	temp_p = myGrowthCurves(sim_p, (t_tune-length(prj_ts))+i, "flounder", max_age = 26)
	temp_p = temp_p/species_params(params_p_sim)["flounder",]$w_inf
	temp_b = myGrowthCurves(sim_b, (t_tune-length(prj_ts))+i, "flounder", max_age = 26)
	temp_b = temp_b/species_params(params_b_sim)["flounder",]$w_inf
	temp_o = myGrowthCurves(sim_o, (t_tune-length(prj_ts))+i, "flounder", max_age = 26)
	temp_o = temp_o/species_params(params_o_sim)["flounder",]$w_inf
	temp_m = myGrowthCurves(sim_m, (t_tune-length(prj_ts))+i, "flounder", max_age = 26)
	temp_m = temp_m/species_params(params_m_sim)["flounder",]$w_inf
	temp_bo = myGrowthCurves(sim_bo, (t_tune-length(prj_ts))+i, "flounder", max_age = 26)
	temp_bo = temp_bo/species_params(params_bo_sim)["flounder",]$w_inf
	temp_bm = myGrowthCurves(sim_bm, (t_tune-length(prj_ts))+i, "flounder", max_age = 26)
	temp_bm = temp_bm/species_params(params_bm_sim)["flounder",]$w_inf
	temp_om = myGrowthCurves(sim_om, (t_tune-length(prj_ts))+i, "flounder", max_age = 26)
	temp_om = temp_om/species_params(params_om_sim)["flounder",]$w_inf
	temp_bom = myGrowthCurves(sim_bom, (t_tune-length(prj_ts))+i, "flounder", max_age = 26)
	temp_bom = temp_bom/species_params(params_bom_sim)["flounder",]$w_inf
	temp_none = myGrowthCurves(sim_none, (t_tune-length(prj_ts))+i, "flounder", max_age = 26)
	temp_none = temp_none/species_params(params_none_sim)["flounder",]$w_inf
	tempobs = flounder_ivb_prj[i,]$a*(flounder_ivb_prj[i,]$L_inf*(1-exp(-flounder_ivb_prj[i,]$k*(seq(0,26,length.out=50)))))^flounder_ivb_prj[i,]$b
	tempobs = tempobs/species_params(params_bomp_sim)["flounder",]$w_inf
	tempdat = flounder_ivb_sim[flounder_ivb_sim$Year == prj_ts[i], ]
	tempdat$Age = tempdat$Age + 0.5
	if(nrow(tempdat)>0) {
		tempdat$Prop = tempdat$IndWgt/species_params(params_bomp_sim)["flounder",]$w_inf
		plot(Prop ~ Age, tempdat, pch = 24, cex = 0.5, axes = F, xlim = c(-2,28), ylim = c(0,1)); axis(1,at=c(0,5,10,15)); box()
		lines(dimnames(temp_bomp)$Age, temp_bomp, lwd = 2, lty = lty[1], col = color[1])
		lines(dimnames(temp_bop)$Age, temp_bop, lwd = 2, lty = lty[2], col = color[2])
		lines(dimnames(temp_bmp)$Age, temp_bmp, lwd = 2, lty = lty[3], col = color[3])
		lines(dimnames(temp_omp)$Age, temp_omp, lwd = 2, lty = lty[4], col = color[4])
		lines(dimnames(temp_bp)$Age, temp_bp, lwd = 2, lty = lty[5], col = color[5])
		lines(dimnames(temp_op)$Age, temp_op, lwd = 2, lty = lty[6], col = color[6])
		lines(dimnames(temp_mp)$Age, temp_mp, lwd = 2, lty = lty[7], col = color[7])
		lines(dimnames(temp_p)$Age, temp_p, lwd = 2, lty = lty[8], col = color[8])
		lines(dimnames(temp_b)$Age, temp_b, lwd = 2, lty = lty[9], col = color[9])
		lines(dimnames(temp_o)$Age, temp_o, lwd = 2, lty = lty[10], col = color[10])
		lines(dimnames(temp_m)$Age, temp_m, lwd = 2, lty = lty[11], col = color[11])
		lines(dimnames(temp_bo)$Age, temp_bo, lwd = 2, lty = lty[12], col = color[12])
		lines(dimnames(temp_bm)$Age, temp_bm, lwd = 2, lty = lty[13], col = color[13])
		lines(dimnames(temp_om)$Age, temp_om, lwd = 2, lty = lty[14], col = color[14])
		lines(dimnames(temp_bom)$Age, temp_bom, lwd = 2, lty = lty[15], col = color[15])
		lines(dimnames(temp_none)$Age, temp_none, lwd = 2, lty = lty[16], col = color[16])
		lines(seq(0,26,length.out=50), tempobs, lwd = 2, lty = 3, col = "dodgerblue")
	} else {
		plot(seq(0,26,length.out=50), tempobs, type = 'n', lwd = 2, axes = F, xlim = c(-2,28), ylim = c(0,1)); axis(1,at=c(0,5,10,15)); box()
		lines(dimnames(temp_bomp)$Age, temp_bomp, lwd = 2, lty = lty[1], col = color[1])
		lines(dimnames(temp_bop)$Age, temp_bop, lwd = 2, lty = lty[2], col = color[2])
		lines(dimnames(temp_bmp)$Age, temp_bmp, lwd = 2, lty = lty[3], col = color[3])
		lines(dimnames(temp_omp)$Age, temp_omp, lwd = 2, lty = lty[4], col = color[4])
		lines(dimnames(temp_bp)$Age, temp_bp, lwd = 2, lty = lty[5], col = color[5])
		lines(dimnames(temp_op)$Age, temp_op, lwd = 2, lty = lty[6], col = color[6])
		lines(dimnames(temp_mp)$Age, temp_mp, lwd = 2, lty = lty[7], col = color[7])
		lines(dimnames(temp_p)$Age, temp_p, lwd = 2, lty = lty[8], col = color[8])
		lines(dimnames(temp_b)$Age, temp_b, lwd = 2, lty = lty[9], col = color[9])
		lines(dimnames(temp_o)$Age, temp_o, lwd = 2, lty = lty[10], col = color[10])
		lines(dimnames(temp_m)$Age, temp_m, lwd = 2, lty = lty[11], col = color[11])
		lines(dimnames(temp_bo)$Age, temp_bo, lwd = 2, lty = lty[12], col = color[12])
		lines(dimnames(temp_bm)$Age, temp_bm, lwd = 2, lty = lty[13], col = color[13])
		lines(dimnames(temp_om)$Age, temp_om, lwd = 2, lty = lty[14], col = color[14])
		lines(dimnames(temp_bom)$Age, temp_bom, lwd = 2, lty = lty[15], col = color[15])
		lines(dimnames(temp_none)$Age, temp_none, lwd = 2, lty = lty[16], col = color[16])
		lines(seq(0,26,length.out=50), tempobs, lwd = 2, lty = 3, col = "dodgerblue")
	}
	if(i == 1) {axis(2,at=c(0,1),labels=c(0,round(species_params(params_bomp_sim)["flounder",]$w_inf)))}
}
for(i in 1:length(prj_ts)) {
	temp_bomp = myGrowthCurves(sim_bomp, (t_tune-length(prj_ts))+i, "sprat", max_age = 16)
	temp_bomp = temp_bomp/species_params(params_bomp_sim)["sprat",]$w_inf
	temp_bop = myGrowthCurves(sim_bop, (t_tune-length(prj_ts))+i, "sprat", max_age = 16)
	temp_bop = temp_bop/species_params(params_bop_sim)["sprat",]$w_inf
	temp_bmp = myGrowthCurves(sim_bmp, (t_tune-length(prj_ts))+i, "sprat", max_age = 16)
	temp_bmp = temp_bmp/species_params(params_bmp_sim)["sprat",]$w_inf
	temp_omp = myGrowthCurves(sim_omp, (t_tune-length(prj_ts))+i, "sprat", max_age = 16)
	temp_omp = temp_omp/species_params(params_omp_sim)["sprat",]$w_inf
	temp_bp = myGrowthCurves(sim_bp, (t_tune-length(prj_ts))+i, "sprat", max_age = 16)
	temp_bp = temp_bp/species_params(params_bp_sim)["sprat",]$w_inf
	temp_op = myGrowthCurves(sim_op, (t_tune-length(prj_ts))+i, "sprat", max_age = 16)
	temp_op = temp_op/species_params(params_op_sim)["sprat",]$w_inf
	temp_mp = myGrowthCurves(sim_mp, (t_tune-length(prj_ts))+i, "sprat", max_age = 16)
	temp_mp = temp_mp/species_params(params_mp_sim)["sprat",]$w_inf
	temp_p = myGrowthCurves(sim_p, (t_tune-length(prj_ts))+i, "sprat", max_age = 16)
	temp_p = temp_p/species_params(params_p_sim)["sprat",]$w_inf
	temp_b = myGrowthCurves(sim_b, (t_tune-length(prj_ts))+i, "sprat", max_age = 16)
	temp_b = temp_b/species_params(params_b_sim)["sprat",]$w_inf
	temp_o = myGrowthCurves(sim_o, (t_tune-length(prj_ts))+i, "sprat", max_age = 16)
	temp_o = temp_o/species_params(params_o_sim)["sprat",]$w_inf
	temp_m = myGrowthCurves(sim_m, (t_tune-length(prj_ts))+i, "sprat", max_age = 16)
	temp_m = temp_m/species_params(params_m_sim)["sprat",]$w_inf
	temp_bo = myGrowthCurves(sim_bo, (t_tune-length(prj_ts))+i, "sprat", max_age = 16)
	temp_bo = temp_bo/species_params(params_bo_sim)["sprat",]$w_inf
	temp_bm = myGrowthCurves(sim_bm, (t_tune-length(prj_ts))+i, "sprat", max_age = 16)
	temp_bm = temp_bm/species_params(params_bm_sim)["sprat",]$w_inf
	temp_om = myGrowthCurves(sim_om, (t_tune-length(prj_ts))+i, "sprat", max_age = 16)
	temp_om = temp_om/species_params(params_om_sim)["sprat",]$w_inf
	temp_bom = myGrowthCurves(sim_bom, (t_tune-length(prj_ts))+i, "sprat", max_age = 16)
	temp_bom = temp_bom/species_params(params_bom_sim)["sprat",]$w_inf
	temp_none = myGrowthCurves(sim_none, (t_tune-length(prj_ts))+i, "sprat", max_age = 16)
	temp_none = temp_none/species_params(params_none_sim)["sprat",]$w_inf
	tempobs = sprat_ivb_prj[i,]$a*(sprat_ivb_prj[i,]$L_inf*(1-exp(-sprat_ivb_prj[i,]$k*(seq(0,16,length.out=50)))))^sprat_ivb_prj[i,]$b
	tempobs = tempobs/species_params(params_bomp_sim)["sprat",]$w_inf
	tempdat = sprat_ivb_sim[sprat_ivb_sim$Year == prj_ts[i], ]
	tempdat$Age = tempdat$Age + 0.5
	if(nrow(tempdat)>0) {
		tempdat$Prop = tempdat$IndWgt/species_params(params_bomp_sim)["sprat",]$w_inf
		plot(Prop ~ Age, tempdat, pch = 24, cex = 0.5, axes = F, xlim = c(-2,18), ylim = c(0,1.3)); axis(1,at=c(0,5,10,15)); box()
		lines(dimnames(temp_bomp)$Age, temp_bomp, lwd = 2, lty = lty[1], col = color[1])
		lines(dimnames(temp_bop)$Age, temp_bop, lwd = 2, lty = lty[2], col = color[2])
		lines(dimnames(temp_bmp)$Age, temp_bmp, lwd = 2, lty = lty[3], col = color[3])
		lines(dimnames(temp_omp)$Age, temp_omp, lwd = 2, lty = lty[4], col = color[4])
		lines(dimnames(temp_bp)$Age, temp_bp, lwd = 2, lty = lty[5], col = color[5])
		lines(dimnames(temp_op)$Age, temp_op, lwd = 2, lty = lty[6], col = color[6])
		lines(dimnames(temp_mp)$Age, temp_mp, lwd = 2, lty = lty[7], col = color[7])
		lines(dimnames(temp_p)$Age, temp_p, lwd = 2, lty = lty[8], col = color[8])
		lines(dimnames(temp_b)$Age, temp_b, lwd = 2, lty = lty[9], col = color[9])
		lines(dimnames(temp_o)$Age, temp_o, lwd = 2, lty = lty[10], col = color[10])
		lines(dimnames(temp_m)$Age, temp_m, lwd = 2, lty = lty[11], col = color[11])
		lines(dimnames(temp_bo)$Age, temp_bo, lwd = 2, lty = lty[12], col = color[12])
		lines(dimnames(temp_bm)$Age, temp_bm, lwd = 2, lty = lty[13], col = color[13])
		lines(dimnames(temp_om)$Age, temp_om, lwd = 2, lty = lty[14], col = color[14])
		lines(dimnames(temp_bom)$Age, temp_bom, lwd = 2, lty = lty[15], col = color[15])
		lines(dimnames(temp_none)$Age, temp_none, lwd = 2, lty = lty[16], col = color[16])
		lines(seq(0,16,length.out=50), tempobs, lwd = 2, lty = 3, col = "dodgerblue")
	} else {
		plot(seq(0,16,length.out=50), tempobs, type = 'n', axes = F, xlim = c(-2,18), ylim = c(0,1.3)); axis(1,at=c(0,5,10,15)); box()
		lines(dimnames(temp_bomp)$Age, temp_bomp, lwd = 2, lty = lty[1], col = color[1])
		lines(dimnames(temp_bop)$Age, temp_bop, lwd = 2, lty = lty[2], col = color[2])
		lines(dimnames(temp_bmp)$Age, temp_bmp, lwd = 2, lty = lty[3], col = color[3])
		lines(dimnames(temp_omp)$Age, temp_omp, lwd = 2, lty = lty[4], col = color[4])
		lines(dimnames(temp_bp)$Age, temp_bp, lwd = 2, lty = lty[5], col = color[5])
		lines(dimnames(temp_op)$Age, temp_op, lwd = 2, lty = lty[6], col = color[6])
		lines(dimnames(temp_mp)$Age, temp_mp, lwd = 2, lty = lty[7], col = color[7])
		lines(dimnames(temp_p)$Age, temp_p, lwd = 2, lty = lty[8], col = color[8])
		lines(dimnames(temp_b)$Age, temp_b, lwd = 2, lty = lty[9], col = color[9])
		lines(dimnames(temp_o)$Age, temp_o, lwd = 2, lty = lty[10], col = color[10])
		lines(dimnames(temp_m)$Age, temp_m, lwd = 2, lty = lty[11], col = color[11])
		lines(dimnames(temp_bo)$Age, temp_bo, lwd = 2, lty = lty[12], col = color[12])
		lines(dimnames(temp_bm)$Age, temp_bm, lwd = 2, lty = lty[13], col = color[13])
		lines(dimnames(temp_om)$Age, temp_om, lwd = 2, lty = lty[14], col = color[14])
		lines(dimnames(temp_bom)$Age, temp_bom, lwd = 2, lty = lty[15], col = color[15])
		lines(dimnames(temp_none)$Age, temp_none, lwd = 2, lty = lty[16], col = color[16])
		lines(seq(0,16,length.out=50), tempobs, lwd = 2, lty = 3, col = "dodgerblue")
	}
	if(i == 1) {axis(2,at=c(0,1),labels=c(0,round(species_params(params_bomp_sim)["sprat",]$w_inf)))}
}
for(i in 1:length(prj_ts)) {
	temp_bomp = myGrowthCurves(sim_bomp, (t_tune-length(prj_ts))+i, "herring", max_age = 13)
	temp_bomp = temp_bomp/species_params(params_bomp_sim)["herring",]$w_inf
	temp_bop = myGrowthCurves(sim_bop, (t_tune-length(prj_ts))+i, "herring", max_age = 13)
	temp_bop = temp_bop/species_params(params_bop_sim)["herring",]$w_inf
	temp_bmp = myGrowthCurves(sim_bmp, (t_tune-length(prj_ts))+i, "herring", max_age = 13)
	temp_bmp = temp_bmp/species_params(params_bmp_sim)["herring",]$w_inf
	temp_omp = myGrowthCurves(sim_omp, (t_tune-length(prj_ts))+i, "herring", max_age = 13)
	temp_omp = temp_omp/species_params(params_omp_sim)["herring",]$w_inf
	temp_bp = myGrowthCurves(sim_bp, (t_tune-length(prj_ts))+i, "herring", max_age = 13)
	temp_bp = temp_bp/species_params(params_bp_sim)["herring",]$w_inf
	temp_op = myGrowthCurves(sim_op, (t_tune-length(prj_ts))+i, "herring", max_age = 13)
	temp_op = temp_op/species_params(params_op_sim)["herring",]$w_inf
	temp_mp = myGrowthCurves(sim_mp, (t_tune-length(prj_ts))+i, "herring", max_age = 13)
	temp_mp = temp_mp/species_params(params_mp_sim)["herring",]$w_inf
	temp_p = myGrowthCurves(sim_p, (t_tune-length(prj_ts))+i, "herring", max_age = 13)
	temp_p = temp_p/species_params(params_p_sim)["herring",]$w_inf
	temp_b = myGrowthCurves(sim_b, (t_tune-length(prj_ts))+i, "herring", max_age = 13)
	temp_b = temp_b/species_params(params_b_sim)["herring",]$w_inf
	temp_o = myGrowthCurves(sim_o, (t_tune-length(prj_ts))+i, "herring", max_age = 13)
	temp_o = temp_o/species_params(params_o_sim)["herring",]$w_inf
	temp_m = myGrowthCurves(sim_m, (t_tune-length(prj_ts))+i, "herring", max_age = 13)
	temp_m = temp_m/species_params(params_m_sim)["herring",]$w_inf
	temp_bo = myGrowthCurves(sim_bo, (t_tune-length(prj_ts))+i, "herring", max_age = 13)
	temp_bo = temp_bo/species_params(params_bo_sim)["herring",]$w_inf
	temp_bm = myGrowthCurves(sim_bm, (t_tune-length(prj_ts))+i, "herring", max_age = 13)
	temp_bm = temp_bm/species_params(params_bm_sim)["herring",]$w_inf
	temp_om = myGrowthCurves(sim_om, (t_tune-length(prj_ts))+i, "herring", max_age = 13)
	temp_om = temp_om/species_params(params_om_sim)["herring",]$w_inf
	temp_bom = myGrowthCurves(sim_bom, (t_tune-length(prj_ts))+i, "herring", max_age = 13)
	temp_bom = temp_bom/species_params(params_bom_sim)["herring",]$w_inf
	temp_none = myGrowthCurves(sim_none, (t_tune-length(prj_ts))+i, "herring", max_age = 13)
	temp_none = temp_none/species_params(params_none_sim)["herring",]$w_inf
	tempobs = herring_ivb_prj[i,]$a*(herring_ivb_prj[i,]$L_inf*(1-exp(-herring_ivb_prj[i,]$k*(seq(0,13,length.out=50)))))^herring_ivb_prj[i,]$b
	tempobs = tempobs/species_params(params_bomp_sim)["herring",]$w_inf
	tempdat = herring_ivb_sim[herring_ivb_sim$Year == prj_ts[i], ]
	tempdat$Age = tempdat$Age + 0.5
	if(nrow(tempdat)>0) {
		tempdat$Prop = tempdat$IndWgt/species_params(params_bomp_sim)["herring",]$w_inf
		plot(Prop ~ Age, tempdat, pch = 24, cex = 0.5, axes = F, xlim = c(-2,18), ylim = c(0,1)); axis(1,at=c(0,5,10,15)); box()
		lines(dimnames(temp_bomp)$Age, temp_bomp, lwd = 2, lty = lty[1], col = color[1])
		lines(dimnames(temp_bop)$Age, temp_bop, lwd = 2, lty = lty[2], col = color[2])
		lines(dimnames(temp_bmp)$Age, temp_bmp, lwd = 2, lty = lty[3], col = color[3])
		lines(dimnames(temp_omp)$Age, temp_omp, lwd = 2, lty = lty[4], col = color[4])
		lines(dimnames(temp_bp)$Age, temp_bp, lwd = 2, lty = lty[5], col = color[5])
		lines(dimnames(temp_op)$Age, temp_op, lwd = 2, lty = lty[6], col = color[6])
		lines(dimnames(temp_mp)$Age, temp_mp, lwd = 2, lty = lty[7], col = color[7])
		lines(dimnames(temp_p)$Age, temp_p, lwd = 2, lty = lty[8], col = color[8])
		lines(dimnames(temp_b)$Age, temp_b, lwd = 2, lty = lty[9], col = color[9])
		lines(dimnames(temp_o)$Age, temp_o, lwd = 2, lty = lty[10], col = color[10])
		lines(dimnames(temp_m)$Age, temp_m, lwd = 2, lty = lty[11], col = color[11])
		lines(dimnames(temp_bo)$Age, temp_bo, lwd = 2, lty = lty[12], col = color[12])
		lines(dimnames(temp_bm)$Age, temp_bm, lwd = 2, lty = lty[13], col = color[13])
		lines(dimnames(temp_om)$Age, temp_om, lwd = 2, lty = lty[14], col = color[14])
		lines(dimnames(temp_bom)$Age, temp_bom, lwd = 2, lty = lty[15], col = color[15])
		lines(dimnames(temp_none)$Age, temp_none, lwd = 2, lty = lty[16], col = color[16])
		lines(seq(0,13,length.out=50), tempobs, lwd = 2, lty = 3, col = "dodgerblue")
	} else {
		plot(seq(0,13,length.out=50), tempobs, type = 'n', axes = F, xlim = c(-2,18), ylim = c(0,1)); axis(1,at=c(0,5,10,15)); box()
		lines(dimnames(temp_bomp)$Age, temp_bomp, lwd = 2, lty = lty[1], col = color[1])
		lines(dimnames(temp_bop)$Age, temp_bop, lwd = 2, lty = lty[2], col = color[2])
		lines(dimnames(temp_bmp)$Age, temp_bmp, lwd = 2, lty = lty[3], col = color[3])
		lines(dimnames(temp_omp)$Age, temp_omp, lwd = 2, lty = lty[4], col = color[4])
		lines(dimnames(temp_bp)$Age, temp_bp, lwd = 2, lty = lty[5], col = color[5])
		lines(dimnames(temp_op)$Age, temp_op, lwd = 2, lty = lty[6], col = color[6])
		lines(dimnames(temp_mp)$Age, temp_mp, lwd = 2, lty = lty[7], col = color[7])
		lines(dimnames(temp_p)$Age, temp_p, lwd = 2, lty = lty[8], col = color[8])
		lines(dimnames(temp_b)$Age, temp_b, lwd = 2, lty = lty[9], col = color[9])
		lines(dimnames(temp_o)$Age, temp_o, lwd = 2, lty = lty[10], col = color[10])
		lines(dimnames(temp_m)$Age, temp_m, lwd = 2, lty = lty[11], col = color[11])
		lines(dimnames(temp_bo)$Age, temp_bo, lwd = 2, lty = lty[12], col = color[12])
		lines(dimnames(temp_bm)$Age, temp_bm, lwd = 2, lty = lty[13], col = color[13])
		lines(dimnames(temp_om)$Age, temp_om, lwd = 2, lty = lty[14], col = color[14])
		lines(dimnames(temp_bom)$Age, temp_bom, lwd = 2, lty = lty[15], col = color[15])
		lines(dimnames(temp_none)$Age, temp_none, lwd = 2, lty = lty[16], col = color[16])
		lines(seq(0,13,length.out=50), tempobs, lwd = 2, lty = 3, col = "dodgerblue")
	}
	if(i == 1) {axis(2,at=c(0,1),labels=c(0,round(species_params(params_bomp_sim)["herring",]$w_inf)))}
}
mtext("Age", side = 1, line = 3, outer = T, cex = 1.5)
mtext("Weight (g)", side = 2, line = 3, outer = T, cex = 1.5)
par(fig = c(0,1,0,1), oma = c(0,0,0,0), mar = c(0,0,0,0), new = T)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend("right", bty = 'n', c("BOMP","BOP","BMP","OMP","BP","OP","MP","P","B","O","M","BO","BM","OM","BOM","None","Data","von B"), col = c(color, "black", "dodgerblue"), lty = c(lty,NA,3), lwd = c(rep(2,16),NA,2), pch = c(rep(NA,16),24,NA), pt.cex = c(rep(NA,16),0.5,NA), xpd = TRUE, cex = 1.5, seg.len = 4)
par(fig = c(0,1,0,1), oma = c(0,0,0,10), mar = c(5,5,3,1), new = T)
plot(0, 0, type = 'l', bty = 'n', xlab = "", ylab = "", xaxt = 'n', yaxt = 'n')
text(0, 0.95, "Cod", adj = 0.5, cex = 1.5)
text(0, seq(0.95,-0.68,length.out=4)[2], "Flounder", adj = 0.5, cex = 1.5)
text(0, seq(0.95,-0.68,length.out=4)[3], "Sprat", adj = 0.5, cex = 1.5)
text(0, -0.68, "Herring", adj = 0.5, cex = 1.5)
par(mfrow = c(1,1))
dev.off()


########## XII. Plot growth during projection period for cod only ##########

# Colors and line types for plotting
color = rep(c("#D55E00","#0072B2","#F0E442","#CC79A7","#E69F00","#009E73","#000000","#D4D4D4"),2)
lty = rep(c(2,1),each=8)

# Plot growth time series for cod only
prj_ts = seq(2001,2019)
jpeg("Plots/cod_growth_eight.jpg", width = 24, height = 22, units = 'cm', res = 600)
par(mfrow = c(4,5), mar = c(0,0,4,0), oma = c(5,5,3,12))
for(i in 1:length(prj_ts)) {
	temp_bomp = myGrowthCurves(sim_bomp, (t_tune-length(prj_ts))+i, "cod", max_age = 8)
	temp_bop = myGrowthCurves(sim_bop, (t_tune-length(prj_ts))+i, "cod", max_age = 8)
	temp_bmp = myGrowthCurves(sim_bmp, (t_tune-length(prj_ts))+i, "cod", max_age = 8)
	temp_omp = myGrowthCurves(sim_omp, (t_tune-length(prj_ts))+i, "cod", max_age = 8)
	temp_bp = myGrowthCurves(sim_bp, (t_tune-length(prj_ts))+i, "cod", max_age = 8)
	temp_op = myGrowthCurves(sim_op, (t_tune-length(prj_ts))+i, "cod", max_age = 8)
	temp_mp = myGrowthCurves(sim_mp, (t_tune-length(prj_ts))+i, "cod", max_age = 8)
	temp_p = myGrowthCurves(sim_p, (t_tune-length(prj_ts))+i, "cod", max_age = 8)
	temp_b = myGrowthCurves(sim_b, (t_tune-length(prj_ts))+i, "cod", max_age = 8)
	temp_o = myGrowthCurves(sim_o, (t_tune-length(prj_ts))+i, "cod", max_age = 8)
	temp_m = myGrowthCurves(sim_m, (t_tune-length(prj_ts))+i, "cod", max_age = 8)
	temp_bo = myGrowthCurves(sim_bo, (t_tune-length(prj_ts))+i, "cod", max_age = 8)
	temp_bm = myGrowthCurves(sim_bm, (t_tune-length(prj_ts))+i, "cod", max_age = 8)
	temp_om = myGrowthCurves(sim_om, (t_tune-length(prj_ts))+i, "cod", max_age = 8)
	temp_bom = myGrowthCurves(sim_bom, (t_tune-length(prj_ts))+i, "cod", max_age = 8)
	temp_none = myGrowthCurves(sim_none, (t_tune-length(prj_ts))+i, "cod", max_age = 8)
	tempdat = cod_ivb_sim[cod_ivb_sim$Year == prj_ts[i], ]
	tempdat = tempdat[tempdat$Age <= 7, ]
	tempdat$Prop = tempdat$IndWgt/species_params(params_bomp_sim)["cod",]$w_inf
	tempdat$Age = tempdat$Age + 0.5
	plot(Prop ~ Age, tempdat, pch = 24, cex = 0.5, axes = F, xlim = c(-1,9), ylim = c(0,0.7)); axis(1,at=seq(0,8,by=2),cex.axis=1.5); box(); mtext(prj_ts[i], cex = 1.2)
	if(i %in% seq(1,19,by=5)) {axis(2,at=c(0,10000/species_params(params_bomp_sim)["cod",]$w_inf),labels=c(0,10000),cex.axis=1.5)}
	temp_bomp = temp_bomp/species_params(params_bomp_sim)["cod",]$w_inf
	lines(dimnames(temp_bomp)$Age, temp_bomp, lwd = 2, lty = lty[1], col = color[1])
	temp_bop = temp_bop/species_params(params_bop_sim)["cod",]$w_inf
	lines(dimnames(temp_bop)$Age, temp_bop, lwd = 2, lty = lty[2], col = color[2])
	temp_bmp = temp_bmp/species_params(params_bmp_sim)["cod",]$w_inf
	lines(dimnames(temp_bmp)$Age, temp_bmp, lwd = 2, lty = lty[3], col = color[3])
	temp_omp = temp_omp/species_params(params_omp_sim)["cod",]$w_inf
	lines(dimnames(temp_omp)$Age, temp_omp, lwd = 2, lty = lty[4], col = color[4])
	temp_bp = temp_bp/species_params(params_bp_sim)["cod",]$w_inf
	lines(dimnames(temp_bp)$Age, temp_bp, lwd = 2, lty = lty[5], col = color[5])
	temp_op = temp_op/species_params(params_op_sim)["cod",]$w_inf
	lines(dimnames(temp_op)$Age, temp_op, lwd = 2, lty = lty[6], col = color[6])
	temp_mp = temp_mp/species_params(params_mp_sim)["cod",]$w_inf
	lines(dimnames(temp_mp)$Age, temp_mp, lwd = 2, lty = lty[7], col = color[7])
	temp_p = temp_p/species_params(params_p_sim)["cod",]$w_inf
	lines(dimnames(temp_p)$Age, temp_p, lwd = 2, lty = lty[8], col = color[8])
	temp_b = temp_b/species_params(params_b_sim)["cod",]$w_inf
	lines(dimnames(temp_b)$Age, temp_b, lwd = 2, lty = lty[9], col = color[9])
	temp_o = temp_o/species_params(params_o_sim)["cod",]$w_inf
	lines(dimnames(temp_o)$Age, temp_o, lwd = 2, lty = lty[10], col = color[10])
	temp_m = temp_m/species_params(params_m_sim)["cod",]$w_inf
	lines(dimnames(temp_m)$Age, temp_m, lwd = 2, lty = lty[11], col = color[11])
	temp_bo = temp_bo/species_params(params_bo_sim)["cod",]$w_inf
	lines(dimnames(temp_bo)$Age, temp_bo, lwd = 2, lty = lty[12], col = color[12])
	temp_bm = temp_bm/species_params(params_bm_sim)["cod",]$w_inf
	lines(dimnames(temp_bm)$Age, temp_bm, lwd = 2, lty = lty[13], col = color[13])
	temp_om = temp_om/species_params(params_om_sim)["cod",]$w_inf
	lines(dimnames(temp_om)$Age, temp_om, lwd = 2, lty = lty[14], col = color[14])
	temp_bom = temp_bom/species_params(params_bom_sim)["cod",]$w_inf
	lines(dimnames(temp_bom)$Age, temp_bom, lwd = 2, lty = lty[15], col = color[15])
	temp_none = temp_none/species_params(params_none_sim)["cod",]$w_inf
	lines(dimnames(temp_none)$Age, temp_none, lwd = 2, lty = lty[16], col = color[16])	
	tempobs = cod_ivb_prj[i,]$a*(cod_ivb_prj[i,]$L_inf*(1-exp(-cod_ivb_prj[i,]$k*(seq(0,15,length.out=50)))))^cod_ivb_prj[i,]$b
	tempobs = tempobs/species_params(params_bomp_sim)["cod",]$w_inf
	lines(seq(0,8,length.out=50), tempobs, lwd = 2, lty = 3, col = "dodgerblue")
}
mtext("Age", side = 1, line = 3, outer = T, cex = 1.5)
mtext("Weight (g)", side = 2, line = 3, outer = T, cex = 1.5)
par(fig = c(0,1,0,1), oma = c(0,0,0,0), mar = c(0,0,0,0), new = T)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend("right", bty = 'n', c("BOMP","BOP","BMP","OMP","BP","OP","MP","P","B","O","M","BO","BM","OM","BOM","None","Data","von B"), col = c(color,"black","dodgerblue"), lty = c(lty,NA,3), lwd = c(rep(2,16),NA,2), pch = c(rep(NA,16),24,NA), pt.cex = c(rep(NA,16),0.5,NA), xpd = TRUE, cex = 1.5, seg.len = 4)
par(fig = c(0,1,0,1), oma = c(0,0,0,0), mar = c(5,5,3,1), new = T)
plot(0, 0, type = 'l', bty = 'n', xlab = "", ylab = "", xaxt = 'n', yaxt = 'n')
par(mfrow = c(1,1))
dev.off()

# Plot growth time series for cod only
prj_ts = seq(2001,2019)
jpeg("Plots/cod_growth_eight_best.jpg", width = 24, height = 22, units = 'cm', res = 600)
par(mfrow = c(4,5), mar = c(0,0,4,0), oma = c(5,5,3,12))
lty = c("dashed","dotted","longdash","dotdash","twodash","9212","solid","solid")

for(i in 1:length(prj_ts)) {
	temp_none = myGrowthCurves(sim_none, (t_tune-length(prj_ts))+i, "cod", max_age = 8)
	temp_b = myGrowthCurves(sim_b, (t_tune-length(prj_ts))+i, "cod", max_age = 8)
	temp_o = myGrowthCurves(sim_o, (t_tune-length(prj_ts))+i, "cod", max_age = 8)
	temp_m = myGrowthCurves(sim_m, (t_tune-length(prj_ts))+i, "cod", max_age = 8)
	temp_bo = myGrowthCurves(sim_bo, (t_tune-length(prj_ts))+i, "cod", max_age = 8)
	temp_bm = myGrowthCurves(sim_bm, (t_tune-length(prj_ts))+i, "cod", max_age = 8)
	temp_om = myGrowthCurves(sim_om, (t_tune-length(prj_ts))+i, "cod", max_age = 8)
	temp_bom = myGrowthCurves(sim_bom, (t_tune-length(prj_ts))+i, "cod", max_age = 8)
	tempdat = cod_ivb_sim[cod_ivb_sim$Year == prj_ts[i], ]
	tempdat = tempdat[tempdat$Age <= 7, ]
	tempdat$Prop = tempdat$IndWgt/species_params(params_bomp_sim)["cod",]$w_inf
	tempdat$Age = tempdat$Age + 0.5
	plot(Prop ~ Age, tempdat, pch = 24, cex = 0.5, axes = F, xlim = c(-1,9), ylim = c(0,0.7)); axis(1,at=seq(0,8,by=2),cex.axis=1.5); box(); mtext(prj_ts[i], cex = 1.2)
	if(i %in% seq(1,19,by=5)) {axis(2,at=c(0,10000/species_params(params_bomp_sim)["cod",]$w_inf),labels=c(0,10000),cex.axis=1.5)}
	temp_none = temp_none/species_params(params_none_sim)["cod",]$w_inf
	lines(dimnames(temp_none)$Age, temp_none, lwd = 4, lty = lty[8], col = color[16])	
	temp_b = temp_b/species_params(params_b_sim)["cod",]$w_inf
	lines(dimnames(temp_b)$Age, temp_b, lwd = 4, lty = lty[1], col = color[9])
	temp_o = temp_o/species_params(params_o_sim)["cod",]$w_inf
	lines(dimnames(temp_o)$Age, temp_o, lwd = 4, lty = lty[2], col = color[10])
	temp_m = temp_m/species_params(params_m_sim)["cod",]$w_inf
	lines(dimnames(temp_m)$Age, temp_m, lwd = 4, lty = lty[3], col = color[11])
	temp_bo = temp_bo/species_params(params_bo_sim)["cod",]$w_inf
	lines(dimnames(temp_bo)$Age, temp_bo, lwd = 4, lty = lty[4], col = color[12])
	temp_bm = temp_bm/species_params(params_bm_sim)["cod",]$w_inf
	lines(dimnames(temp_bm)$Age, temp_bm, lwd = 4, lty = lty[5], col = color[13])
	temp_om = temp_om/species_params(params_om_sim)["cod",]$w_inf
	lines(dimnames(temp_om)$Age, temp_om, lwd = 4, lty = lty[6], col = color[14])
	temp_bom = temp_bom/species_params(params_bom_sim)["cod",]$w_inf
	lines(dimnames(temp_bom)$Age, temp_bom, lwd = 4, lty = lty[7], col = color[15])
	tempobs = cod_ivb_prj[i,]$a*(cod_ivb_prj[i,]$L_inf*(1-exp(-cod_ivb_prj[i,]$k*(seq(0,15,length.out=50)))))^cod_ivb_prj[i,]$b
	tempobs = tempobs/species_params(params_bomp_sim)["cod",]$w_inf
	lines(seq(0,8,length.out=50), tempobs, lwd = 4, lty = 3, col = "dodgerblue")
}
mtext("Age", side = 1, line = 3, outer = T, cex = 1.5)
mtext("Weight (g)", side = 2, line = 3, outer = T, cex = 1.5)
par(fig = c(0,1,0,1), oma = c(0,0,0,0), mar = c(0,0,0,0), new = T)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend("right", bty = 'n', c("B","O","M","BO","BM","OM","BOM","None","Data","von B"), col = c(color[9:16],"black","dodgerblue"), lty = c(lty,NA,"dotted"), lwd = c(rep(4,8),NA,4), pch = c(rep(NA,8),24,NA), pt.cex = c(rep(NA,8),0.5,NA), xpd = TRUE, cex = 1.5, seg.len = 4)
par(fig = c(0,1,0,1), oma = c(0,0,0,0), mar = c(5,5,3,1), new = T)
plot(0, 0, type = 'l', bty = 'n', xlab = "", ylab = "", xaxt = 'n', yaxt = 'n')
par(mfrow = c(1,1))
dev.off()


########## XIII. Plot growth during projection period for first eight years of life ##########

# Colors and line types for plotting
color = rep(c("#D55E00","#0072B2","#F0E442","#CC79A7","#E69F00","#009E73","#000000","#D4D4D4"),2)
lty = rep(c(2,1),each=8)

# Plot growth time series of the first eight years of life
prj_ts = seq(2001,2019)
jpeg("Plots/sim_growth_eight.jpg", width = 44, height = 20, units = 'cm', res = 600)
par(mfrow = c(4,19), mar = c(0,0,5,0), oma = c(5,5,3,12))
for(i in 1:length(prj_ts)) {
	temp_bomp = myGrowthCurves(sim_bomp, (t_tune-length(prj_ts))+i, "cod", max_age = 8)
	temp_bop = myGrowthCurves(sim_bop, (t_tune-length(prj_ts))+i, "cod", max_age = 8)
	temp_bmp = myGrowthCurves(sim_bmp, (t_tune-length(prj_ts))+i, "cod", max_age = 8)
	temp_omp = myGrowthCurves(sim_omp, (t_tune-length(prj_ts))+i, "cod", max_age = 8)
	temp_bp = myGrowthCurves(sim_bp, (t_tune-length(prj_ts))+i, "cod", max_age = 8)
	temp_op = myGrowthCurves(sim_op, (t_tune-length(prj_ts))+i, "cod", max_age = 8)
	temp_mp = myGrowthCurves(sim_mp, (t_tune-length(prj_ts))+i, "cod", max_age = 8)
	temp_p = myGrowthCurves(sim_p, (t_tune-length(prj_ts))+i, "cod", max_age = 8)
	temp_b = myGrowthCurves(sim_b, (t_tune-length(prj_ts))+i, "cod", max_age = 8)
	temp_o = myGrowthCurves(sim_o, (t_tune-length(prj_ts))+i, "cod", max_age = 8)
	temp_m = myGrowthCurves(sim_m, (t_tune-length(prj_ts))+i, "cod", max_age = 8)
	temp_bo = myGrowthCurves(sim_bo, (t_tune-length(prj_ts))+i, "cod", max_age = 8)
	temp_bm = myGrowthCurves(sim_bm, (t_tune-length(prj_ts))+i, "cod", max_age = 8)
	temp_om = myGrowthCurves(sim_om, (t_tune-length(prj_ts))+i, "cod", max_age = 8)
	temp_bom = myGrowthCurves(sim_bom, (t_tune-length(prj_ts))+i, "cod", max_age = 8)
	temp_none = myGrowthCurves(sim_none, (t_tune-length(prj_ts))+i, "cod", max_age = 8)
	tempdat = cod_ivb_sim[cod_ivb_sim$Year == prj_ts[i], ]
	tempdat$Prop = tempdat$IndWgt/species_params(params_bomp_sim)["cod",]$w_inf
	tempdat$Age = tempdat$Age + 0.5
	plot(Prop ~ Age, tempdat[tempdat$Age <= 8, ], pch = 24, cex = 0.5, axes = F, xlim = c(-2,10), ylim = c(0,0.5)); axis(1,at=c(0,4,8)); box(); mtext(prj_ts[i], cex = 0.75)
	if(i == 1) {axis(2,at=c(0,0.5),labels=c(0,round(species_params(params_bomp_sim)["cod",]$w_inf/2)))}
	temp_bomp = temp_bomp/species_params(params_bomp_sim)["cod",]$w_inf
	lines(dimnames(temp_bomp)$Age, temp_bomp, lwd = 2, lty = lty[1], col = color[1])
	temp_bop = temp_bop/species_params(params_bop_sim)["cod",]$w_inf
	lines(dimnames(temp_bop)$Age, temp_bop, lwd = 2, lty = lty[2], col = color[2])
	temp_bmp = temp_bmp/species_params(params_bmp_sim)["cod",]$w_inf
	lines(dimnames(temp_bmp)$Age, temp_bmp, lwd = 2, lty = lty[3], col = color[3])
	temp_omp = temp_omp/species_params(params_omp_sim)["cod",]$w_inf
	lines(dimnames(temp_omp)$Age, temp_omp, lwd = 2, lty = lty[4], col = color[4])
	temp_bp = temp_bp/species_params(params_bp_sim)["cod",]$w_inf
	lines(dimnames(temp_bp)$Age, temp_bp, lwd = 2, lty = lty[5], col = color[5])
	temp_op = temp_op/species_params(params_op_sim)["cod",]$w_inf
	lines(dimnames(temp_op)$Age, temp_op, lwd = 2, lty = lty[6], col = color[6])
	temp_mp = temp_mp/species_params(params_mp_sim)["cod",]$w_inf
	lines(dimnames(temp_mp)$Age, temp_mp, lwd = 2, lty = lty[7], col = color[7])
	temp_p = temp_p/species_params(params_p_sim)["cod",]$w_inf
	lines(dimnames(temp_p)$Age, temp_p, lwd = 2, lty = lty[8], col = color[8])
	temp_b = temp_b/species_params(params_b_sim)["cod",]$w_inf
	lines(dimnames(temp_b)$Age, temp_b, lwd = 2, lty = lty[9], col = color[9])
	temp_o = temp_o/species_params(params_o_sim)["cod",]$w_inf
	lines(dimnames(temp_o)$Age, temp_o, lwd = 2, lty = lty[10], col = color[10])
	temp_m = temp_m/species_params(params_m_sim)["cod",]$w_inf
	lines(dimnames(temp_m)$Age, temp_m, lwd = 2, lty = lty[11], col = color[11])
	temp_bo = temp_bo/species_params(params_bo_sim)["cod",]$w_inf
	lines(dimnames(temp_bo)$Age, temp_bo, lwd = 2, lty = lty[12], col = color[12])
	temp_bm = temp_bm/species_params(params_bm_sim)["cod",]$w_inf
	lines(dimnames(temp_bm)$Age, temp_bm, lwd = 2, lty = lty[13], col = color[13])
	temp_om = temp_om/species_params(params_om_sim)["cod",]$w_inf
	lines(dimnames(temp_om)$Age, temp_om, lwd = 2, lty = lty[14], col = color[14])
	temp_bom = temp_bom/species_params(params_bom_sim)["cod",]$w_inf
	lines(dimnames(temp_bom)$Age, temp_bom, lwd = 2, lty = lty[15], col = color[15])
	temp_none = temp_none/species_params(params_none_sim)["cod",]$w_inf
	lines(dimnames(temp_none)$Age, temp_none, lwd = 2, lty = lty[16], col = color[16])	
	tempobs = cod_ivb_prj[i,]$a*(cod_ivb_prj[i,]$L_inf*(1-exp(-cod_ivb_prj[i,]$k*(seq(0,8,length.out=50)))))^cod_ivb_prj[i,]$b
	tempobs = tempobs/species_params(params_bomp_sim)["cod",]$w_inf
	lines(seq(0,8,length.out=50), tempobs, lwd = 2, lty = 3, col = "dodgerblue")
}
for(i in 1:length(prj_ts)) {
	temp_bomp = myGrowthCurves(sim_bomp, (t_tune-length(prj_ts))+i, "flounder", max_age = 8)
	temp_bomp = temp_bomp/species_params(params_bomp_sim)["flounder",]$w_inf
	temp_bop = myGrowthCurves(sim_bop, (t_tune-length(prj_ts))+i, "flounder", max_age = 8)
	temp_bop = temp_bop/species_params(params_bop_sim)["flounder",]$w_inf
	temp_bmp = myGrowthCurves(sim_bmp, (t_tune-length(prj_ts))+i, "flounder", max_age = 8)
	temp_bmp = temp_bmp/species_params(params_bmp_sim)["flounder",]$w_inf
	temp_omp = myGrowthCurves(sim_omp, (t_tune-length(prj_ts))+i, "flounder", max_age = 8)
	temp_omp = temp_omp/species_params(params_omp_sim)["flounder",]$w_inf
	temp_bp = myGrowthCurves(sim_bp, (t_tune-length(prj_ts))+i, "flounder", max_age = 8)
	temp_bp = temp_bp/species_params(params_bp_sim)["flounder",]$w_inf
	temp_op = myGrowthCurves(sim_op, (t_tune-length(prj_ts))+i, "flounder", max_age = 8)
	temp_op = temp_op/species_params(params_op_sim)["flounder",]$w_inf
	temp_mp = myGrowthCurves(sim_mp, (t_tune-length(prj_ts))+i, "flounder", max_age = 8)
	temp_mp = temp_mp/species_params(params_mp_sim)["flounder",]$w_inf
	temp_p = myGrowthCurves(sim_p, (t_tune-length(prj_ts))+i, "flounder", max_age = 8)
	temp_p = temp_p/species_params(params_p_sim)["flounder",]$w_inf
	temp_b = myGrowthCurves(sim_b, (t_tune-length(prj_ts))+i, "flounder", max_age = 8)
	temp_b = temp_b/species_params(params_b_sim)["flounder",]$w_inf
	temp_o = myGrowthCurves(sim_o, (t_tune-length(prj_ts))+i, "flounder", max_age = 8)
	temp_o = temp_o/species_params(params_o_sim)["flounder",]$w_inf
	temp_m = myGrowthCurves(sim_m, (t_tune-length(prj_ts))+i, "flounder", max_age = 8)
	temp_m = temp_m/species_params(params_m_sim)["flounder",]$w_inf
	temp_bo = myGrowthCurves(sim_bo, (t_tune-length(prj_ts))+i, "flounder", max_age = 8)
	temp_bo = temp_bo/species_params(params_bo_sim)["flounder",]$w_inf
	temp_bm = myGrowthCurves(sim_bm, (t_tune-length(prj_ts))+i, "flounder", max_age = 8)
	temp_bm = temp_bm/species_params(params_bm_sim)["flounder",]$w_inf
	temp_om = myGrowthCurves(sim_om, (t_tune-length(prj_ts))+i, "flounder", max_age = 8)
	temp_om = temp_om/species_params(params_om_sim)["flounder",]$w_inf
	temp_bom = myGrowthCurves(sim_bom, (t_tune-length(prj_ts))+i, "flounder", max_age = 8)
	temp_bom = temp_bom/species_params(params_bom_sim)["flounder",]$w_inf
	temp_none = myGrowthCurves(sim_none, (t_tune-length(prj_ts))+i, "flounder", max_age = 8)
	temp_none = temp_none/species_params(params_none_sim)["flounder",]$w_inf
	tempobs = flounder_ivb_prj[i,]$a*(flounder_ivb_prj[i,]$L_inf*(1-exp(-flounder_ivb_prj[i,]$k*(seq(0,8,length.out=50)))))^flounder_ivb_prj[i,]$b
	tempobs = tempobs/species_params(params_bomp_sim)["flounder",]$w_inf
	tempdat = flounder_ivb_sim[flounder_ivb_sim$Year == prj_ts[i], ]
	tempdat$Age = tempdat$Age + 0.5
	if(nrow(tempdat)>0) {
		tempdat$Prop = tempdat$IndWgt/species_params(params_bomp_sim)["flounder",]$w_inf
		plot(Prop ~ Age, tempdat[tempdat$Age <= 8, ], pch = 24, cex = 0.5, axes = F, xlim = c(-2,10), ylim = c(0,1)); axis(1,at=c(0,4,8)); box()
		lines(dimnames(temp_bomp)$Age, temp_bomp, lwd = 2, lty = lty[1], col = color[1])
		lines(dimnames(temp_bop)$Age, temp_bop, lwd = 2, lty = lty[2], col = color[2])
		lines(dimnames(temp_bmp)$Age, temp_bmp, lwd = 2, lty = lty[3], col = color[3])
		lines(dimnames(temp_omp)$Age, temp_omp, lwd = 2, lty = lty[4], col = color[4])
		lines(dimnames(temp_bp)$Age, temp_bp, lwd = 2, lty = lty[5], col = color[5])
		lines(dimnames(temp_op)$Age, temp_op, lwd = 2, lty = lty[6], col = color[6])
		lines(dimnames(temp_mp)$Age, temp_mp, lwd = 2, lty = lty[7], col = color[7])
		lines(dimnames(temp_p)$Age, temp_p, lwd = 2, lty = lty[8], col = color[8])
		lines(dimnames(temp_b)$Age, temp_b, lwd = 2, lty = lty[9], col = color[9])
		lines(dimnames(temp_o)$Age, temp_o, lwd = 2, lty = lty[10], col = color[10])
		lines(dimnames(temp_m)$Age, temp_m, lwd = 2, lty = lty[11], col = color[11])
		lines(dimnames(temp_bo)$Age, temp_bo, lwd = 2, lty = lty[12], col = color[12])
		lines(dimnames(temp_bm)$Age, temp_bm, lwd = 2, lty = lty[13], col = color[13])
		lines(dimnames(temp_om)$Age, temp_om, lwd = 2, lty = lty[14], col = color[14])
		lines(dimnames(temp_bom)$Age, temp_bom, lwd = 2, lty = lty[15], col = color[15])
		lines(dimnames(temp_none)$Age, temp_none, lwd = 2, lty = lty[16], col = color[16])
		lines(seq(0,8,length.out=50), tempobs, lwd = 2, lty = 3, col = "dodgerblue")
	} else {
		plot(seq(0,8,length.out=50), tempobs, type = 'n', axes = F, xlim = c(-2,10), ylim = c(0,1)); axis(1,at=c(0,4,8)); box()
		lines(dimnames(temp_bomp)$Age, temp_bomp, lwd = 2, lty = lty[1], col = color[1])
		lines(dimnames(temp_bop)$Age, temp_bop, lwd = 2, lty = lty[2], col = color[2])
		lines(dimnames(temp_bmp)$Age, temp_bmp, lwd = 2, lty = lty[3], col = color[3])
		lines(dimnames(temp_omp)$Age, temp_omp, lwd = 2, lty = lty[4], col = color[4])
		lines(dimnames(temp_bp)$Age, temp_bp, lwd = 2, lty = lty[5], col = color[5])
		lines(dimnames(temp_op)$Age, temp_op, lwd = 2, lty = lty[6], col = color[6])
		lines(dimnames(temp_mp)$Age, temp_mp, lwd = 2, lty = lty[7], col = color[7])
		lines(dimnames(temp_p)$Age, temp_p, lwd = 2, lty = lty[8], col = color[8])
		lines(dimnames(temp_b)$Age, temp_b, lwd = 2, lty = lty[9], col = color[9])
		lines(dimnames(temp_o)$Age, temp_o, lwd = 2, lty = lty[10], col = color[10])
		lines(dimnames(temp_m)$Age, temp_m, lwd = 2, lty = lty[11], col = color[11])
		lines(dimnames(temp_bo)$Age, temp_bo, lwd = 2, lty = lty[12], col = color[12])
		lines(dimnames(temp_bm)$Age, temp_bm, lwd = 2, lty = lty[13], col = color[13])
		lines(dimnames(temp_om)$Age, temp_om, lwd = 2, lty = lty[14], col = color[14])
		lines(dimnames(temp_bom)$Age, temp_bom, lwd = 2, lty = lty[15], col = color[15])
		lines(dimnames(temp_none)$Age, temp_none, lwd = 2, lty = lty[16], col = color[16])
		lines(seq(0,8,length.out=50), tempobs, lwd = 2, lty = 3, col = "dodgerblue")
	}
	if(i == 1) {axis(2,at=c(0,1),labels=c(0,round(species_params(params_bomp_sim)["flounder",]$w_inf)))}
}
for(i in 1:length(prj_ts)) {
	temp_bomp = myGrowthCurves(sim_bomp, (t_tune-length(prj_ts))+i, "sprat", max_age = 8)
	temp_bomp = temp_bomp/species_params(params_bomp_sim)["sprat",]$w_inf
	temp_bop = myGrowthCurves(sim_bop, (t_tune-length(prj_ts))+i, "sprat", max_age = 8)
	temp_bop = temp_bop/species_params(params_bop_sim)["sprat",]$w_inf
	temp_bmp = myGrowthCurves(sim_bmp, (t_tune-length(prj_ts))+i, "sprat", max_age = 8)
	temp_bmp = temp_bmp/species_params(params_bmp_sim)["sprat",]$w_inf
	temp_omp = myGrowthCurves(sim_omp, (t_tune-length(prj_ts))+i, "sprat", max_age = 8)
	temp_omp = temp_omp/species_params(params_omp_sim)["sprat",]$w_inf
	temp_bp = myGrowthCurves(sim_bp, (t_tune-length(prj_ts))+i, "sprat", max_age = 8)
	temp_bp = temp_bp/species_params(params_bp_sim)["sprat",]$w_inf
	temp_op = myGrowthCurves(sim_op, (t_tune-length(prj_ts))+i, "sprat", max_age = 8)
	temp_op = temp_op/species_params(params_op_sim)["sprat",]$w_inf
	temp_mp = myGrowthCurves(sim_mp, (t_tune-length(prj_ts))+i, "sprat", max_age = 8)
	temp_mp = temp_mp/species_params(params_mp_sim)["sprat",]$w_inf
	temp_p = myGrowthCurves(sim_p, (t_tune-length(prj_ts))+i, "sprat", max_age = 8)
	temp_p = temp_p/species_params(params_p_sim)["sprat",]$w_inf	
	temp_b = myGrowthCurves(sim_b, (t_tune-length(prj_ts))+i, "sprat", max_age = 8)
	temp_b = temp_b/species_params(params_b_sim)["sprat",]$w_inf
	temp_o = myGrowthCurves(sim_o, (t_tune-length(prj_ts))+i, "sprat", max_age = 8)
	temp_o = temp_o/species_params(params_o_sim)["sprat",]$w_inf
	temp_m = myGrowthCurves(sim_m, (t_tune-length(prj_ts))+i, "sprat", max_age = 8)
	temp_m = temp_m/species_params(params_m_sim)["sprat",]$w_inf
	temp_bo = myGrowthCurves(sim_bo, (t_tune-length(prj_ts))+i, "sprat", max_age = 8)
	temp_bo = temp_bo/species_params(params_bo_sim)["sprat",]$w_inf
	temp_bm = myGrowthCurves(sim_bm, (t_tune-length(prj_ts))+i, "sprat", max_age = 8)
	temp_bm = temp_bm/species_params(params_bm_sim)["sprat",]$w_inf
	temp_om = myGrowthCurves(sim_om, (t_tune-length(prj_ts))+i, "sprat", max_age = 8)
	temp_om = temp_om/species_params(params_om_sim)["sprat",]$w_inf
	temp_bom = myGrowthCurves(sim_bom, (t_tune-length(prj_ts))+i, "sprat", max_age = 8)
	temp_bom = temp_bom/species_params(params_bom_sim)["sprat",]$w_inf
	temp_none = myGrowthCurves(sim_none, (t_tune-length(prj_ts))+i, "sprat", max_age = 8)
	temp_none = temp_none/species_params(params_none_sim)["sprat",]$w_inf
	tempobs = sprat_ivb_prj[i,]$a*(sprat_ivb_prj[i,]$L_inf*(1-exp(-sprat_ivb_prj[i,]$k*(seq(0,8,length.out=50)))))^sprat_ivb_prj[i,]$b
	tempobs = tempobs/species_params(params_bomp_sim)["sprat",]$w_inf
	tempdat = sprat_ivb_sim[sprat_ivb_sim$Year == prj_ts[i], ]
	tempdat$Age = tempdat$Age + 0.5
	if(nrow(tempdat)>0) {
		tempdat$Prop = tempdat$IndWgt/species_params(params_bomp_sim)["sprat",]$w_inf
		plot(Prop ~ Age, tempdat[tempdat$Age <= 8, ], pch = 24, cex = 0.5, axes = F, xlim = c(-2,10), ylim = c(0,1)); axis(1,at=c(0,4,8)); box()
		lines(dimnames(temp_bomp)$Age, temp_bomp, lwd = 2, lty = lty[1], col = color[1])
		lines(dimnames(temp_bop)$Age, temp_bop, lwd = 2, lty = lty[2], col = color[2])
		lines(dimnames(temp_bmp)$Age, temp_bmp, lwd = 2, lty = lty[3], col = color[3])
		lines(dimnames(temp_omp)$Age, temp_omp, lwd = 2, lty = lty[4], col = color[4])
		lines(dimnames(temp_bp)$Age, temp_bp, lwd = 2, lty = lty[5], col = color[5])
		lines(dimnames(temp_op)$Age, temp_op, lwd = 2, lty = lty[6], col = color[6])
		lines(dimnames(temp_mp)$Age, temp_mp, lwd = 2, lty = lty[7], col = color[7])
		lines(dimnames(temp_p)$Age, temp_p, lwd = 2, lty = lty[8], col = color[8])
		lines(dimnames(temp_b)$Age, temp_b, lwd = 2, lty = lty[9], col = color[9])
		lines(dimnames(temp_o)$Age, temp_o, lwd = 2, lty = lty[10], col = color[10])
		lines(dimnames(temp_m)$Age, temp_m, lwd = 2, lty = lty[11], col = color[11])
		lines(dimnames(temp_bo)$Age, temp_bo, lwd = 2, lty = lty[12], col = color[12])
		lines(dimnames(temp_bm)$Age, temp_bm, lwd = 2, lty = lty[13], col = color[13])
		lines(dimnames(temp_om)$Age, temp_om, lwd = 2, lty = lty[14], col = color[14])
		lines(dimnames(temp_bom)$Age, temp_bom, lwd = 2, lty = lty[15], col = color[15])
		lines(dimnames(temp_none)$Age, temp_none, lwd = 2, lty = lty[16], col = color[16])
		lines(seq(0,8,length.out=50), tempobs, lwd = 2, lty = 3, col = "dodgerblue")
	} else {
		plot(seq(0,8,length.out=50), tempobs, type = 'l', lwd = 2, lty = lty[1], col = "red", axes = F, xlim = c(-2,10), ylim = c(0,1)); axis(1,at=c(0,4,8)); box()
		lines(dimnames(temp_bomp)$Age, temp_bomp, lwd = 2, lty = lty[1], col = color[1])
		lines(dimnames(temp_bop)$Age, temp_bop, lwd = 2, lty = lty[2], col = color[2])
		lines(dimnames(temp_bmp)$Age, temp_bmp, lwd = 2, lty = lty[3], col = color[3])
		lines(dimnames(temp_omp)$Age, temp_omp, lwd = 2, lty = lty[4], col = color[4])
		lines(dimnames(temp_bp)$Age, temp_bp, lwd = 2, lty = lty[5], col = color[5])
		lines(dimnames(temp_op)$Age, temp_op, lwd = 2, lty = lty[6], col = color[6])
		lines(dimnames(temp_mp)$Age, temp_mp, lwd = 2, lty = lty[7], col = color[7])
		lines(dimnames(temp_p)$Age, temp_p, lwd = 2, lty = lty[8], col = color[8])
		lines(dimnames(temp_b)$Age, temp_b, lwd = 2, lty = lty[9], col = color[9])
		lines(dimnames(temp_o)$Age, temp_o, lwd = 2, lty = lty[10], col = color[10])
		lines(dimnames(temp_m)$Age, temp_m, lwd = 2, lty = lty[11], col = color[11])
		lines(dimnames(temp_bo)$Age, temp_bo, lwd = 2, lty = lty[12], col = color[12])
		lines(dimnames(temp_bm)$Age, temp_bm, lwd = 2, lty = lty[13], col = color[13])
		lines(dimnames(temp_om)$Age, temp_om, lwd = 2, lty = lty[14], col = color[14])
		lines(dimnames(temp_bom)$Age, temp_bom, lwd = 2, lty = lty[15], col = color[15])
		lines(dimnames(temp_none)$Age, temp_none, lwd = 2, lty = lty[16], col = color[16])
		lines(seq(0,8,length.out=50), tempobs, lwd = 2, lty = 3, col = "dodgerblue")
	}
	if(i == 1) {axis(2,at=c(0,1),labels=c(0,round(species_params(params_bomp_sim)["sprat",]$w_inf)))}
}
for(i in 1:length(prj_ts)) {
	temp_bomp = myGrowthCurves(sim_bomp, (t_tune-length(prj_ts))+i, "herring", max_age = 8)
	temp_bomp = temp_bomp/species_params(params_bomp_sim)["herring",]$w_inf
	temp_bop = myGrowthCurves(sim_bop, (t_tune-length(prj_ts))+i, "herring", max_age = 8)
	temp_bop = temp_bop/species_params(params_bop_sim)["herring",]$w_inf
	temp_bmp = myGrowthCurves(sim_bmp, (t_tune-length(prj_ts))+i, "herring", max_age = 8)
	temp_bmp = temp_bmp/species_params(params_bmp_sim)["herring",]$w_inf
	temp_omp = myGrowthCurves(sim_omp, (t_tune-length(prj_ts))+i, "herring", max_age = 8)
	temp_omp = temp_omp/species_params(params_omp_sim)["herring",]$w_inf
	temp_bp = myGrowthCurves(sim_bp, (t_tune-length(prj_ts))+i, "herring", max_age = 8)
	temp_bp = temp_bp/species_params(params_bp_sim)["herring",]$w_inf
	temp_op = myGrowthCurves(sim_op, (t_tune-length(prj_ts))+i, "herring", max_age = 8)
	temp_op = temp_op/species_params(params_op_sim)["herring",]$w_inf
	temp_mp = myGrowthCurves(sim_mp, (t_tune-length(prj_ts))+i, "herring", max_age = 8)
	temp_mp = temp_mp/species_params(params_mp_sim)["herring",]$w_inf
	temp_p = myGrowthCurves(sim_p, (t_tune-length(prj_ts))+i, "herring", max_age = 8)
	temp_p = temp_p/species_params(params_p_sim)["herring",]$w_inf
	temp_b = myGrowthCurves(sim_b, (t_tune-length(prj_ts))+i, "herring", max_age = 8)
	temp_b = temp_b/species_params(params_b_sim)["herring",]$w_inf
	temp_o = myGrowthCurves(sim_o, (t_tune-length(prj_ts))+i, "herring", max_age = 8)
	temp_o = temp_o/species_params(params_o_sim)["herring",]$w_inf
	temp_m = myGrowthCurves(sim_m, (t_tune-length(prj_ts))+i, "herring", max_age = 8)
	temp_m = temp_m/species_params(params_m_sim)["herring",]$w_inf
	temp_bo = myGrowthCurves(sim_bo, (t_tune-length(prj_ts))+i, "herring", max_age = 8)
	temp_bo = temp_bo/species_params(params_bo_sim)["herring",]$w_inf
	temp_bm = myGrowthCurves(sim_bm, (t_tune-length(prj_ts))+i, "herring", max_age = 8)
	temp_bm = temp_bm/species_params(params_bm_sim)["herring",]$w_inf
	temp_om = myGrowthCurves(sim_om, (t_tune-length(prj_ts))+i, "herring", max_age = 8)
	temp_om = temp_om/species_params(params_om_sim)["herring",]$w_inf
	temp_bom = myGrowthCurves(sim_bom, (t_tune-length(prj_ts))+i, "herring", max_age = 8)
	temp_bom = temp_bom/species_params(params_bom_sim)["herring",]$w_inf
	temp_none = myGrowthCurves(sim_none, (t_tune-length(prj_ts))+i, "herring", max_age = 8)
	temp_none = temp_none/species_params(params_none_sim)["herring",]$w_inf
	tempobs = herring_ivb_prj[i,]$a*(herring_ivb_prj[i,]$L_inf*(1-exp(-herring_ivb_prj[i,]$k*(seq(0,8,length.out=50)))))^herring_ivb_prj[i,]$b
	tempobs = tempobs/species_params(params_bomp_sim)["herring",]$w_inf
	tempdat = herring_ivb_sim[herring_ivb_sim$Year == prj_ts[i], ]
	tempdat$Age = tempdat$Age + 0.5
	if(nrow(tempdat)>0) {
		tempdat$Prop = tempdat$IndWgt/species_params(params_bomp_sim)["herring",]$w_inf
		plot(Prop ~ Age, tempdat[tempdat$Age <= 8, ], pch = 24, cex = 0.5, axes = F, xlim = c(-2,10), ylim = c(0,1)); axis(1,at=c(0,4,8)); box()
		lines(dimnames(temp_bomp)$Age, temp_bomp, lwd = 2, lty = lty[1], col = color[1])
		lines(dimnames(temp_bop)$Age, temp_bop, lwd = 2, lty = lty[2], col = color[2])
		lines(dimnames(temp_bmp)$Age, temp_bmp, lwd = 2, lty = lty[3], col = color[3])
		lines(dimnames(temp_omp)$Age, temp_omp, lwd = 2, lty = lty[4], col = color[4])
		lines(dimnames(temp_bp)$Age, temp_bp, lwd = 2, lty = lty[5], col = color[5])
		lines(dimnames(temp_op)$Age, temp_op, lwd = 2, lty = lty[6], col = color[6])
		lines(dimnames(temp_mp)$Age, temp_mp, lwd = 2, lty = lty[7], col = color[7])
		lines(dimnames(temp_p)$Age, temp_p, lwd = 2, lty = lty[8], col = color[8])
		lines(dimnames(temp_b)$Age, temp_b, lwd = 2, lty = lty[9], col = color[9])
		lines(dimnames(temp_o)$Age, temp_o, lwd = 2, lty = lty[10], col = color[10])
		lines(dimnames(temp_m)$Age, temp_m, lwd = 2, lty = lty[11], col = color[11])
		lines(dimnames(temp_bo)$Age, temp_bo, lwd = 2, lty = lty[12], col = color[12])
		lines(dimnames(temp_bm)$Age, temp_bm, lwd = 2, lty = lty[13], col = color[13])
		lines(dimnames(temp_om)$Age, temp_om, lwd = 2, lty = lty[14], col = color[14])
		lines(dimnames(temp_bom)$Age, temp_bom, lwd = 2, lty = lty[15], col = color[15])
		lines(dimnames(temp_none)$Age, temp_none, lwd = 2, lty = lty[16], col = color[16])
		lines(seq(0,8,length.out=50), tempobs, lwd = 2, lty = 3, col = "dodgerblue")
	} else {
		plot(seq(0,8,length.out=50), tempobs, type = 'l', lwd = 2, lty = lty[1], col = "red", axes = F, xlim = c(-2,10), ylim = c(0,1)); axis(1,at=c(0,4,8)); box()
		lines(dimnames(temp_bomp)$Age, temp_bomp, lwd = 2, lty = lty[1], col = color[1])
		lines(dimnames(temp_bop)$Age, temp_bop, lwd = 2, lty = lty[2], col = color[2])
		lines(dimnames(temp_bmp)$Age, temp_bmp, lwd = 2, lty = lty[3], col = color[3])
		lines(dimnames(temp_omp)$Age, temp_omp, lwd = 2, lty = lty[4], col = color[4])
		lines(dimnames(temp_bp)$Age, temp_bp, lwd = 2, lty = lty[5], col = color[5])
		lines(dimnames(temp_op)$Age, temp_op, lwd = 2, lty = lty[6], col = color[6])
		lines(dimnames(temp_mp)$Age, temp_mp, lwd = 2, lty = lty[7], col = color[7])
		lines(dimnames(temp_p)$Age, temp_p, lwd = 2, lty = lty[8], col = color[8])
		lines(dimnames(temp_b)$Age, temp_b, lwd = 2, lty = lty[9], col = color[9])
		lines(dimnames(temp_o)$Age, temp_o, lwd = 2, lty = lty[10], col = color[10])
		lines(dimnames(temp_m)$Age, temp_m, lwd = 2, lty = lty[11], col = color[11])
		lines(dimnames(temp_bo)$Age, temp_bo, lwd = 2, lty = lty[12], col = color[12])
		lines(dimnames(temp_bm)$Age, temp_bm, lwd = 2, lty = lty[13], col = color[13])
		lines(dimnames(temp_om)$Age, temp_om, lwd = 2, lty = lty[14], col = color[14])
		lines(dimnames(temp_bom)$Age, temp_bom, lwd = 2, lty = lty[15], col = color[15])
		lines(dimnames(temp_none)$Age, temp_none, lwd = 2, lty = lty[16], col = color[16])
		lines(seq(0,8,length.out=50), tempobs, lwd = 2, lty = 3, col = "dodgerblue")
	}
	if(i == 1) {axis(2,at=c(0,1),labels=c(0,round(species_params(params_bomp_sim)["herring",]$w_inf)))}
}
mtext("Age", side = 1, line = 3, outer = T, cex = 1.5)
mtext("Weight (g)", side = 2, line = 3, outer = T, cex = 1.5)
par(fig = c(0,1,0,1), oma = c(0,0,0,0), mar = c(0,0,0,0), new = T)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend("right", bty = 'n', c("BOMP","BOP","BMP","OMP","BP","OP","MP","P","B","O","M","BO","BM","OM","BOM","None","Data","von B"), col = c(color,"black","dodgerblue"), lty = c(lty,NA,3), lwd = c(rep(2,12),4,4,4,2,NA,2), pch = c(rep(NA,16),24,NA), pt.cex = c(rep(NA,16),0.5,NA), xpd = TRUE, cex = 1.5, seg.len = 4)
par(fig = c(0,1,0,1), oma = c(0,0,0,10), mar = c(5,5,3,1), new = T)
plot(0, 0, type = 'l', bty = 'n', xlab = "", ylab = "", xaxt = 'n', yaxt = 'n')
text(0, 0.95, "Cod", adj = 0.5, cex = 1.5)
text(0, seq(0.95,-0.68,length.out=4)[2], "Flounder", adj = 0.5, cex = 1.5)
text(0, seq(0.95,-0.68,length.out=4)[3], "Sprat", adj = 0.5, cex = 1.5)
text(0, -0.68, "Herring", adj = 0.5, cex = 1.5)
par(mfrow = c(1,1))
dev.off()


########## XIV. Plot cod occupancy ##########

# A function to get occupancy values 
#	params: mizer params object
#	t: time steps to run model
#	effort: a vector of length four containing effort for each species
#	lty: line type
#	color: line color
#	lwd: line width
#	returns nothing, plots a line of occupancy values during 1991-2019
getocc = function(params, t = 100, effort, lty, color, lwd) {
	# Set up params object
	ret = params
	
	# Run to steady state
	sim = project(ret, t_max = t, effort = effort)
	
	# Set up params object
	ret = sim@params
	ret@initial_n[] <- sim@n[dim(sim@n)[1],, ]
    	ret@initial_n_pp[] <- sim@n_pp[dim(sim@n)[1], ]
    	ret@initial_n_other[] <- sim@n_other[dim(sim@n)[1], ]
    	ret@initial_effort[] <- sim@effort[dim(sim@n)[1], ]
    	ret@time_modified <- lubridate::now()
	
	# Calculate occupancy
	bocc = ret@other_params$benthic_occupancy[,"cod",params@w >= cod_lh["a","stan"]*species_params(ret)$l50[1]^cod_lh["b","stan"]]
	bocc_wm = rowSums(sweep(bocc, 2, initialN(ret)["cod",params@w >= cod_lh["a","stan"]*species_params(ret)$l50[1]^cod_lh["b","stan"]], "*")) / sum(initialN(ret)["cod",params@w >= cod_lh["a","stan"]*species_params(ret)$l50[1]^cod_lh["b","stan"]])

	# Add lines
	lines(sim_ts, bocc_wm[73:101], lty = lty, lwd = lwd, col = color)
}

# Save plot of occupancy over the years
jpeg("Plots/occupancy.jpg", width = 24, height = 16, units = 'cm', res = 600)
par(oma = c(0,0,0,10))
plot(sim_ts, type = 'n', lwd = 3, xlim = c(1991,2019), ylim = c(0.60,1.00), col = "red", pch=16, cex.lab = 1.5, cex.axis = 1.2, xlab = "Year", ylab = "Benthic occupancy")
getocc(params_bomp, effort = c(1,1,1,1), lty = lty[1], color = color[1], lwd = 2)
getocc(params_bop, effort = c(1,1,1,1), lty = lty[2], color = color[2], lwd = 2)
getocc(params_bmp, effort = c(1,1,1,1), lty = lty[3], color = color[3], lwd = 2)
getocc(params_omp, effort = c(1,1,1,1), lty = lty[4], color = color[4], lwd = 2)
getocc(params_bp, effort = c(1,1,1,1), lty = lty[5], color = color[5], lwd = 2)
getocc(params_op, effort = c(1,1,1,1), lty = lty[6], color = color[6], lwd = 2)
getocc(params_mp, effort = c(1,1,1,1), lty = lty[7], color = color[7], lwd = 2)
getocc(params_p, effort = c(1,1,1,1), lty = lty[8], color = color[8], lwd = 2)
getocc(params_b, effort = c(1,1,1,1), lty = lty[9], color = color[9], lwd = 2)
getocc(params_o, effort = c(1,1,1,1), lty = lty[10], color = color[10], lwd = 2)
getocc(params_m, effort = c(1,1,1,1), lty = lty[11], color = color[11], lwd = 2)
getocc(params_bo, effort = c(1,1,1,1), lty = lty[12], color = color[12], lwd = 2)
getocc(params_bm, effort = c(1,1,1,1), lty = lty[13], color = color[13], lwd = 2)
getocc(params_om, effort = c(1,1,1,1), lty = lty[14], color = color[14], lwd = 2)
getocc(params_bom, effort = c(1,1,1,1), lty = lty[15], color = color[15], lwd = 2)
getocc(params_none, effort = c(1,1,1,1), lty = lty[16], color = color[16], lwd = 2)
polygon(x = c(1991,2000,2000,1991), y = c(0.6,0.6,1,1), border = NA, col = rgb(0.5,0.5,0.5,0.2))
points(sim_ts[sim_ts <= 2015], cod_benthic[cod_benthic$Year %in% sim_ts, ]$Benthic_occ, pch = 17, col = "dodgerblue", cex = 1.5)	
par(fig = c(0,1,0,1), oma = c(0,0,0,0), mar = c(0,0,0,0), new = T)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend("right", bty = 'n', c("BOMP","BOP","BMP","OMP","BP","OP","MP","P","B","O","M","BO","BM","OM","BOM","None","Data"), col = c(color,"dodgerblue"), lty = c(lty,NA), lwd = c(rep(2,12),4,4,4,2,NA), pch = c(rep(NA,16),17), pt.cex = c(rep(NA,16),1.5), xpd = TRUE, cex = 1.2, seg.len = 4)
dev.off()


########## XV. Calculate error ##########

# Calibration years
cal_ts = seq(1991,2000)

# Projection years
prj_ts = seq(2001,2019)

# Model labels
model_labels = c(
	"BOMP (Full)",
	"BOM",
	"BOP",
	"BMP",
	"OMP",
	"BO",
	"BM",
	"BP",
	"OM",
	"OP",
	"MP",
	"B",
	"O",
	"M",
	"P",
	"None"
)

# SSB calibration error
cbe_cod_bomp = sum((log(cod_ssb[cod_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_bomp@n[(t_tune-length(cal_ts)+1):t_tune,1,], 2, params_bomp@w*params_bomp@dw, "*"))[,params_bomp@w >= species_params(params_bomp)$w_mat[1]])))^2)
cbe_cod_bom = sum((log(cod_ssb[cod_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_bom@n[(t_tune-length(cal_ts)+1):t_tune,1,], 2, params_bom@w*params_bom@dw, "*"))[,params_bom@w >= species_params(params_bom)$w_mat[1]])))^2)
cbe_cod_bop = sum((log(cod_ssb[cod_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_bop@n[(t_tune-length(cal_ts)+1):t_tune,1,], 2, params_bop@w*params_bop@dw, "*"))[,params_bop@w >= species_params(params_bop)$w_mat[1]])))^2)
cbe_cod_bmp = sum((log(cod_ssb[cod_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_bmp@n[(t_tune-length(cal_ts)+1):t_tune,1,], 2, params_bmp@w*params_bmp@dw, "*"))[,params_bmp@w >= species_params(params_bmp)$w_mat[1]])))^2)
cbe_cod_omp = sum((log(cod_ssb[cod_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_omp@n[(t_tune-length(cal_ts)+1):t_tune,1,], 2, params_omp@w*params_omp@dw, "*"))[,params_omp@w >= species_params(params_omp)$w_mat[1]])))^2)
cbe_cod_bo = sum((log(cod_ssb[cod_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_bo@n[(t_tune-length(cal_ts)+1):t_tune,1,], 2, params_bo@w*params_bo@dw, "*"))[,params_bo@w >= species_params(params_bo)$w_mat[1]])))^2)
cbe_cod_bm = sum((log(cod_ssb[cod_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_bm@n[(t_tune-length(cal_ts)+1):t_tune,1,], 2, params_bm@w*params_bm@dw, "*"))[,params_bm@w >= species_params(params_bm)$w_mat[1]])))^2)
cbe_cod_bp = sum((log(cod_ssb[cod_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_bp@n[(t_tune-length(cal_ts)+1):t_tune,1,], 2, params_bp@w*params_bp@dw, "*"))[,params_bp@w >= species_params(params_bp)$w_mat[1]])))^2)
cbe_cod_om = sum((log(cod_ssb[cod_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_om@n[(t_tune-length(cal_ts)+1):t_tune,1,], 2, params_om@w*params_om@dw, "*"))[,params_om@w >= species_params(params_om)$w_mat[1]])))^2)
cbe_cod_op = sum((log(cod_ssb[cod_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_op@n[(t_tune-length(cal_ts)+1):t_tune,1,], 2, params_op@w*params_op@dw, "*"))[,params_op@w >= species_params(params_op)$w_mat[1]])))^2)
cbe_cod_mp = sum((log(cod_ssb[cod_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_mp@n[(t_tune-length(cal_ts)+1):t_tune,1,], 2, params_mp@w*params_mp@dw, "*"))[,params_mp@w >= species_params(params_mp)$w_mat[1]])))^2)
cbe_cod_b = sum((log(cod_ssb[cod_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_b@n[(t_tune-length(cal_ts)+1):t_tune,1,], 2, params_b@w*params_b@dw, "*"))[,params_b@w >= species_params(params_b)$w_mat[1]])))^2)
cbe_cod_o = sum((log(cod_ssb[cod_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_o@n[(t_tune-length(cal_ts)+1):t_tune,1,], 2, params_o@w*params_o@dw, "*"))[,params_o@w >= species_params(params_o)$w_mat[1]])))^2)
cbe_cod_m = sum((log(cod_ssb[cod_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_m@n[(t_tune-length(cal_ts)+1):t_tune,1,], 2, params_m@w*params_m@dw, "*"))[,params_m@w >= species_params(params_m)$w_mat[1]])))^2)
cbe_cod_p = sum((log(cod_ssb[cod_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_p@n[(t_tune-length(cal_ts)+1):t_tune,1,], 2, params_p@w*params_p@dw, "*"))[,params_p@w >= species_params(params_p)$w_mat[1]])))^2)
cbe_cod_none = sum((log(cod_ssb[cod_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_none@n[(t_tune-length(cal_ts)+1):t_tune,1,], 2, params_none@w*params_none@dw, "*"))[,params_none@w >= species_params(params_none)$w_mat[1]])))^2)

cbe_flounder_bomp = sum((log(flounder_ssb[flounder_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_bomp@n[(t_tune-length(cal_ts)+1):t_tune,2,], 2, params_bomp@w*params_bomp@dw, "*"))[,params_bomp@w >= species_params(params_bomp)$w_mat[2]])))^2)
cbe_flounder_bom = sum((log(flounder_ssb[flounder_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_bom@n[(t_tune-length(cal_ts)+1):t_tune,2,], 2, params_bom@w*params_bom@dw, "*"))[,params_bom@w >= species_params(params_bom)$w_mat[2]])))^2)
cbe_flounder_bop = sum((log(flounder_ssb[flounder_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_bop@n[(t_tune-length(cal_ts)+1):t_tune,2,], 2, params_bop@w*params_bop@dw, "*"))[,params_bop@w >= species_params(params_bop)$w_mat[2]])))^2)
cbe_flounder_bmp = sum((log(flounder_ssb[flounder_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_bmp@n[(t_tune-length(cal_ts)+1):t_tune,2,], 2, params_bmp@w*params_bmp@dw, "*"))[,params_bmp@w >= species_params(params_bmp)$w_mat[2]])))^2)
cbe_flounder_omp = sum((log(flounder_ssb[flounder_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_omp@n[(t_tune-length(cal_ts)+1):t_tune,2,], 2, params_omp@w*params_omp@dw, "*"))[,params_omp@w >= species_params(params_omp)$w_mat[2]])))^2)
cbe_flounder_bo = sum((log(flounder_ssb[flounder_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_bo@n[(t_tune-length(cal_ts)+1):t_tune,2,], 2, params_bo@w*params_bo@dw, "*"))[,params_bo@w >= species_params(params_bo)$w_mat[2]])))^2)
cbe_flounder_bm = sum((log(flounder_ssb[flounder_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_bm@n[(t_tune-length(cal_ts)+1):t_tune,2,], 2, params_bm@w*params_bm@dw, "*"))[,params_bm@w >= species_params(params_bm)$w_mat[2]])))^2)
cbe_flounder_bp = sum((log(flounder_ssb[flounder_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_bp@n[(t_tune-length(cal_ts)+1):t_tune,2,], 2, params_bp@w*params_bp@dw, "*"))[,params_bp@w >= species_params(params_bp)$w_mat[2]])))^2)
cbe_flounder_om = sum((log(flounder_ssb[flounder_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_om@n[(t_tune-length(cal_ts)+1):t_tune,2,], 2, params_om@w*params_om@dw, "*"))[,params_om@w >= species_params(params_om)$w_mat[2]])))^2)
cbe_flounder_op = sum((log(flounder_ssb[flounder_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_op@n[(t_tune-length(cal_ts)+1):t_tune,2,], 2, params_op@w*params_op@dw, "*"))[,params_op@w >= species_params(params_op)$w_mat[2]])))^2)
cbe_flounder_mp = sum((log(flounder_ssb[flounder_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_mp@n[(t_tune-length(cal_ts)+1):t_tune,2,], 2, params_mp@w*params_mp@dw, "*"))[,params_mp@w >= species_params(params_mp)$w_mat[2]])))^2)
cbe_flounder_b = sum((log(flounder_ssb[flounder_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_b@n[(t_tune-length(cal_ts)+1):t_tune,2,], 2, params_b@w*params_b@dw, "*"))[,params_b@w >= species_params(params_b)$w_mat[2]])))^2)
cbe_flounder_o = sum((log(flounder_ssb[flounder_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_o@n[(t_tune-length(cal_ts)+1):t_tune,2,], 2, params_o@w*params_o@dw, "*"))[,params_o@w >= species_params(params_o)$w_mat[2]])))^2)
cbe_flounder_m = sum((log(flounder_ssb[flounder_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_m@n[(t_tune-length(cal_ts)+1):t_tune,2,], 2, params_m@w*params_m@dw, "*"))[,params_m@w >= species_params(params_m)$w_mat[2]])))^2)
cbe_flounder_p = sum((log(flounder_ssb[flounder_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_p@n[(t_tune-length(cal_ts)+1):t_tune,2,], 2, params_p@w*params_p@dw, "*"))[,params_p@w >= species_params(params_p)$w_mat[2]])))^2)
cbe_flounder_none = sum((log(flounder_ssb[flounder_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_none@n[(t_tune-length(cal_ts)+1):t_tune,2,], 2, params_none@w*params_none@dw, "*"))[,params_none@w >= species_params(params_none)$w_mat[2]])))^2)

cbe_sprat_bomp = sum((log(sprat_ssb[sprat_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_bomp@n[(t_tune-length(cal_ts)+1):t_tune,3,], 2, params_bomp@w*params_bomp@dw, "*"))[,params_bomp@w >= species_params(params_bomp)$w_mat[3]])))^2)
cbe_sprat_bom = sum((log(sprat_ssb[sprat_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_bom@n[(t_tune-length(cal_ts)+1):t_tune,3,], 2, params_bom@w*params_bom@dw, "*"))[,params_bom@w >= species_params(params_bom)$w_mat[3]])))^2)
cbe_sprat_bop = sum((log(sprat_ssb[sprat_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_bop@n[(t_tune-length(cal_ts)+1):t_tune,3,], 2, params_bop@w*params_bop@dw, "*"))[,params_bop@w >= species_params(params_bop)$w_mat[3]])))^2)
cbe_sprat_bmp = sum((log(sprat_ssb[sprat_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_bmp@n[(t_tune-length(cal_ts)+1):t_tune,3,], 2, params_bmp@w*params_bmp@dw, "*"))[,params_bmp@w >= species_params(params_bmp)$w_mat[3]])))^2)
cbe_sprat_omp = sum((log(sprat_ssb[sprat_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_omp@n[(t_tune-length(cal_ts)+1):t_tune,3,], 2, params_omp@w*params_omp@dw, "*"))[,params_omp@w >= species_params(params_omp)$w_mat[3]])))^2)
cbe_sprat_bo = sum((log(sprat_ssb[sprat_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_bo@n[(t_tune-length(cal_ts)+1):t_tune,3,], 2, params_bo@w*params_bo@dw, "*"))[,params_bo@w >= species_params(params_bo)$w_mat[3]])))^2)
cbe_sprat_bm = sum((log(sprat_ssb[sprat_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_bm@n[(t_tune-length(cal_ts)+1):t_tune,3,], 2, params_bm@w*params_bm@dw, "*"))[,params_bm@w >= species_params(params_bm)$w_mat[3]])))^2)
cbe_sprat_bp = sum((log(sprat_ssb[sprat_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_bp@n[(t_tune-length(cal_ts)+1):t_tune,3,], 2, params_bp@w*params_bp@dw, "*"))[,params_bp@w >= species_params(params_bp)$w_mat[3]])))^2)
cbe_sprat_om = sum((log(sprat_ssb[sprat_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_om@n[(t_tune-length(cal_ts)+1):t_tune,3,], 2, params_om@w*params_om@dw, "*"))[,params_om@w >= species_params(params_om)$w_mat[3]])))^2)
cbe_sprat_op = sum((log(sprat_ssb[sprat_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_op@n[(t_tune-length(cal_ts)+1):t_tune,3,], 2, params_op@w*params_op@dw, "*"))[,params_op@w >= species_params(params_op)$w_mat[3]])))^2)
cbe_sprat_mp = sum((log(sprat_ssb[sprat_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_mp@n[(t_tune-length(cal_ts)+1):t_tune,3,], 2, params_mp@w*params_mp@dw, "*"))[,params_mp@w >= species_params(params_mp)$w_mat[3]])))^2)
cbe_sprat_b = sum((log(sprat_ssb[sprat_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_b@n[(t_tune-length(cal_ts)+1):t_tune,3,], 2, params_b@w*params_b@dw, "*"))[,params_b@w >= species_params(params_b)$w_mat[3]])))^2)
cbe_sprat_o = sum((log(sprat_ssb[sprat_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_o@n[(t_tune-length(cal_ts)+1):t_tune,3,], 2, params_o@w*params_o@dw, "*"))[,params_o@w >= species_params(params_o)$w_mat[3]])))^2)
cbe_sprat_m = sum((log(sprat_ssb[sprat_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_m@n[(t_tune-length(cal_ts)+1):t_tune,3,], 2, params_m@w*params_m@dw, "*"))[,params_m@w >= species_params(params_m)$w_mat[3]])))^2)
cbe_sprat_p = sum((log(sprat_ssb[sprat_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_p@n[(t_tune-length(cal_ts)+1):t_tune,3,], 2, params_p@w*params_p@dw, "*"))[,params_p@w >= species_params(params_p)$w_mat[3]])))^2)
cbe_sprat_none = sum((log(sprat_ssb[sprat_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_none@n[(t_tune-length(cal_ts)+1):t_tune,3,], 2, params_none@w*params_none@dw, "*"))[,params_none@w >= species_params(params_none)$w_mat[3]])))^2)

cbe_herring_bomp = sum((log(herring_ssb[herring_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_bomp@n[(t_tune-length(cal_ts)+1):t_tune,4,], 2, params_bomp@w*params_bomp@dw, "*"))[,params_bomp@w >= species_params(params_bomp)$w_mat[4]])))^2)
cbe_herring_bom = sum((log(herring_ssb[herring_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_bom@n[(t_tune-length(cal_ts)+1):t_tune,4,], 2, params_bom@w*params_bom@dw, "*"))[,params_bom@w >= species_params(params_bom)$w_mat[4]])))^2)
cbe_herring_bop = sum((log(herring_ssb[herring_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_bop@n[(t_tune-length(cal_ts)+1):t_tune,4,], 2, params_bop@w*params_bop@dw, "*"))[,params_bop@w >= species_params(params_bop)$w_mat[4]])))^2)
cbe_herring_bmp = sum((log(herring_ssb[herring_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_bmp@n[(t_tune-length(cal_ts)+1):t_tune,4,], 2, params_bmp@w*params_bmp@dw, "*"))[,params_bmp@w >= species_params(params_bmp)$w_mat[4]])))^2)
cbe_herring_omp = sum((log(herring_ssb[herring_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_omp@n[(t_tune-length(cal_ts)+1):t_tune,4,], 2, params_omp@w*params_omp@dw, "*"))[,params_omp@w >= species_params(params_omp)$w_mat[4]])))^2)
cbe_herring_bo = sum((log(herring_ssb[herring_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_bo@n[(t_tune-length(cal_ts)+1):t_tune,4,], 2, params_bo@w*params_bo@dw, "*"))[,params_bo@w >= species_params(params_bo)$w_mat[4]])))^2)
cbe_herring_bm = sum((log(herring_ssb[herring_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_bm@n[(t_tune-length(cal_ts)+1):t_tune,4,], 2, params_bm@w*params_bm@dw, "*"))[,params_bm@w >= species_params(params_bm)$w_mat[4]])))^2)
cbe_herring_bp = sum((log(herring_ssb[herring_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_bp@n[(t_tune-length(cal_ts)+1):t_tune,4,], 2, params_bp@w*params_bp@dw, "*"))[,params_bp@w >= species_params(params_bp)$w_mat[4]])))^2)
cbe_herring_om = sum((log(herring_ssb[herring_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_om@n[(t_tune-length(cal_ts)+1):t_tune,4,], 2, params_om@w*params_om@dw, "*"))[,params_om@w >= species_params(params_om)$w_mat[4]])))^2)
cbe_herring_op = sum((log(herring_ssb[herring_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_op@n[(t_tune-length(cal_ts)+1):t_tune,4,], 2, params_op@w*params_op@dw, "*"))[,params_op@w >= species_params(params_op)$w_mat[4]])))^2)
cbe_herring_mp = sum((log(herring_ssb[herring_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_mp@n[(t_tune-length(cal_ts)+1):t_tune,4,], 2, params_mp@w*params_mp@dw, "*"))[,params_mp@w >= species_params(params_mp)$w_mat[4]])))^2)
cbe_herring_b = sum((log(herring_ssb[herring_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_b@n[(t_tune-length(cal_ts)+1):t_tune,4,], 2, params_b@w*params_b@dw, "*"))[,params_b@w >= species_params(params_b)$w_mat[4]])))^2)
cbe_herring_o = sum((log(herring_ssb[herring_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_o@n[(t_tune-length(cal_ts)+1):t_tune,4,], 2, params_o@w*params_o@dw, "*"))[,params_o@w >= species_params(params_o)$w_mat[4]])))^2)
cbe_herring_m = sum((log(herring_ssb[herring_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_m@n[(t_tune-length(cal_ts)+1):t_tune,4,], 2, params_m@w*params_m@dw, "*"))[,params_m@w >= species_params(params_m)$w_mat[4]])))^2)
cbe_herring_p = sum((log(herring_ssb[herring_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_p@n[(t_tune-length(cal_ts)+1):t_tune,4,], 2, params_p@w*params_p@dw, "*"))[,params_p@w >= species_params(params_p)$w_mat[4]])))^2)
cbe_herring_none = sum((log(herring_ssb[herring_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_none@n[(t_tune-length(cal_ts)+1):t_tune,4,], 2, params_none@w*params_none@dw, "*"))[,params_none@w >= species_params(params_none)$w_mat[4]])))^2)

# Sum across species
cbe_cod = c(cbe_cod_bomp, cbe_cod_bom, cbe_cod_bop, cbe_cod_bmp, cbe_cod_omp, cbe_cod_bo, cbe_cod_bm, cbe_cod_bp, cbe_cod_om, cbe_cod_op, cbe_cod_mp, cbe_cod_b, cbe_cod_o, cbe_cod_m, cbe_cod_p, cbe_cod_none)
cbe_flounder = c(cbe_flounder_bomp, cbe_flounder_bom, cbe_flounder_bop, cbe_flounder_bmp, cbe_flounder_omp, cbe_flounder_bo, cbe_flounder_bm, cbe_flounder_bp, cbe_flounder_om, cbe_flounder_op, cbe_flounder_mp, cbe_flounder_b, cbe_flounder_o, cbe_flounder_m, cbe_flounder_p, cbe_flounder_none)
cbe_sprat = c(cbe_sprat_bomp, cbe_sprat_bom, cbe_sprat_bop, cbe_sprat_bmp, cbe_sprat_omp, cbe_sprat_bo, cbe_sprat_bm, cbe_sprat_bp, cbe_sprat_om, cbe_sprat_op, cbe_sprat_mp, cbe_sprat_b, cbe_sprat_o, cbe_sprat_m, cbe_sprat_p, cbe_sprat_none)
cbe_herring = c(cbe_herring_bomp, cbe_herring_bom, cbe_herring_bop, cbe_herring_bmp, cbe_herring_omp, cbe_herring_bo, cbe_herring_bm, cbe_herring_bp, cbe_herring_om, cbe_herring_op, cbe_herring_mp, cbe_herring_b, cbe_herring_o, cbe_herring_m, cbe_herring_p, cbe_herring_none)

# Total and ranked error
cbe_cod_weight = 0.3
cbe_flounder_weight = 0.1 
cbe_sprat_weight = 0.3
cbe_herring_weight = 0.3
cbe_error = cbe_cod_weight*cbe_cod + cbe_flounder_weight*cbe_flounder + cbe_sprat_weight*cbe_sprat + cbe_herring_weight*cbe_herring
cbe_rank = (cbe_error - min(cbe_error)) / (max(cbe_error) - min(cbe_error))

# SSB projection error
pbe_cod_bomp = sum((log(cod_ssb[cod_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_bomp@n[(t_tune-length(prj_ts)+1):t_tune,1,], 2, params_bomp@w*params_bomp@dw, "*"))[,params_bomp@w >= species_params(params_bomp)$w_mat[1]])))^2)
pbe_cod_bom = sum((log(cod_ssb[cod_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_bom@n[(t_tune-length(prj_ts)+1):t_tune,1,], 2, params_bom@w*params_bom@dw, "*"))[,params_bom@w >= species_params(params_bom)$w_mat[1]])))^2)
pbe_cod_bop = sum((log(cod_ssb[cod_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_bop@n[(t_tune-length(prj_ts)+1):t_tune,1,], 2, params_bop@w*params_bop@dw, "*"))[,params_bop@w >= species_params(params_bop)$w_mat[1]])))^2)
pbe_cod_bmp = sum((log(cod_ssb[cod_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_bmp@n[(t_tune-length(prj_ts)+1):t_tune,1,], 2, params_bmp@w*params_bmp@dw, "*"))[,params_bmp@w >= species_params(params_bmp)$w_mat[1]])))^2)
pbe_cod_omp = sum((log(cod_ssb[cod_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_omp@n[(t_tune-length(prj_ts)+1):t_tune,1,], 2, params_omp@w*params_omp@dw, "*"))[,params_omp@w >= species_params(params_omp)$w_mat[1]])))^2)
pbe_cod_bo = sum((log(cod_ssb[cod_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_bo@n[(t_tune-length(prj_ts)+1):t_tune,1,], 2, params_bo@w*params_bo@dw, "*"))[,params_bo@w >= species_params(params_bo)$w_mat[1]])))^2)
pbe_cod_bm = sum((log(cod_ssb[cod_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_bm@n[(t_tune-length(prj_ts)+1):t_tune,1,], 2, params_bm@w*params_bm@dw, "*"))[,params_bm@w >= species_params(params_bm)$w_mat[1]])))^2)
pbe_cod_bp = sum((log(cod_ssb[cod_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_bp@n[(t_tune-length(prj_ts)+1):t_tune,1,], 2, params_bp@w*params_bp@dw, "*"))[,params_bp@w >= species_params(params_bp)$w_mat[1]])))^2)
pbe_cod_om = sum((log(cod_ssb[cod_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_om@n[(t_tune-length(prj_ts)+1):t_tune,1,], 2, params_om@w*params_om@dw, "*"))[,params_om@w >= species_params(params_om)$w_mat[1]])))^2)
pbe_cod_op = sum((log(cod_ssb[cod_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_op@n[(t_tune-length(prj_ts)+1):t_tune,1,], 2, params_op@w*params_op@dw, "*"))[,params_op@w >= species_params(params_op)$w_mat[1]])))^2)
pbe_cod_mp = sum((log(cod_ssb[cod_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_mp@n[(t_tune-length(prj_ts)+1):t_tune,1,], 2, params_mp@w*params_mp@dw, "*"))[,params_mp@w >= species_params(params_mp)$w_mat[1]])))^2)
pbe_cod_b = sum((log(cod_ssb[cod_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_b@n[(t_tune-length(prj_ts)+1):t_tune,1,], 2, params_b@w*params_b@dw, "*"))[,params_b@w >= species_params(params_b)$w_mat[1]])))^2)
pbe_cod_o = sum((log(cod_ssb[cod_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_o@n[(t_tune-length(prj_ts)+1):t_tune,1,], 2, params_o@w*params_o@dw, "*"))[,params_o@w >= species_params(params_o)$w_mat[1]])))^2)
pbe_cod_m = sum((log(cod_ssb[cod_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_m@n[(t_tune-length(prj_ts)+1):t_tune,1,], 2, params_m@w*params_m@dw, "*"))[,params_m@w >= species_params(params_m)$w_mat[1]])))^2)
pbe_cod_p = sum((log(cod_ssb[cod_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_p@n[(t_tune-length(prj_ts)+1):t_tune,1,], 2, params_p@w*params_p@dw, "*"))[,params_p@w >= species_params(params_p)$w_mat[1]])))^2)
pbe_cod_none = sum((log(cod_ssb[cod_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_none@n[(t_tune-length(prj_ts)+1):t_tune,1,], 2, params_none@w*params_none@dw, "*"))[,params_none@w >= species_params(params_none)$w_mat[1]])))^2)

pbe_flounder_bomp = sum((log(flounder_ssb[flounder_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_bomp@n[(t_tune-length(prj_ts)+1):t_tune,2,], 2, params_bomp@w*params_bomp@dw, "*"))[,params_bomp@w >= species_params(params_bomp)$w_mat[2]])))^2)
pbe_flounder_bom = sum((log(flounder_ssb[flounder_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_bom@n[(t_tune-length(prj_ts)+1):t_tune,2,], 2, params_bom@w*params_bom@dw, "*"))[,params_bom@w >= species_params(params_bom)$w_mat[2]])))^2)
pbe_flounder_bop = sum((log(flounder_ssb[flounder_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_bop@n[(t_tune-length(prj_ts)+1):t_tune,2,], 2, params_bop@w*params_bop@dw, "*"))[,params_bop@w >= species_params(params_bop)$w_mat[2]])))^2)
pbe_flounder_bmp = sum((log(flounder_ssb[flounder_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_bmp@n[(t_tune-length(prj_ts)+1):t_tune,2,], 2, params_bmp@w*params_bmp@dw, "*"))[,params_bmp@w >= species_params(params_bmp)$w_mat[2]])))^2)
pbe_flounder_omp = sum((log(flounder_ssb[flounder_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_omp@n[(t_tune-length(prj_ts)+1):t_tune,2,], 2, params_omp@w*params_omp@dw, "*"))[,params_omp@w >= species_params(params_omp)$w_mat[2]])))^2)
pbe_flounder_bo = sum((log(flounder_ssb[flounder_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_bo@n[(t_tune-length(prj_ts)+1):t_tune,2,], 2, params_bo@w*params_bo@dw, "*"))[,params_bo@w >= species_params(params_bo)$w_mat[2]])))^2)
pbe_flounder_bm = sum((log(flounder_ssb[flounder_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_bm@n[(t_tune-length(prj_ts)+1):t_tune,2,], 2, params_bm@w*params_bm@dw, "*"))[,params_bm@w >= species_params(params_bm)$w_mat[2]])))^2)
pbe_flounder_bp = sum((log(flounder_ssb[flounder_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_bp@n[(t_tune-length(prj_ts)+1):t_tune,2,], 2, params_bp@w*params_bp@dw, "*"))[,params_bp@w >= species_params(params_bp)$w_mat[2]])))^2)
pbe_flounder_om = sum((log(flounder_ssb[flounder_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_om@n[(t_tune-length(prj_ts)+1):t_tune,2,], 2, params_om@w*params_om@dw, "*"))[,params_om@w >= species_params(params_om)$w_mat[2]])))^2)
pbe_flounder_op = sum((log(flounder_ssb[flounder_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_op@n[(t_tune-length(prj_ts)+1):t_tune,2,], 2, params_op@w*params_op@dw, "*"))[,params_op@w >= species_params(params_op)$w_mat[2]])))^2)
pbe_flounder_mp = sum((log(flounder_ssb[flounder_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_mp@n[(t_tune-length(prj_ts)+1):t_tune,2,], 2, params_mp@w*params_mp@dw, "*"))[,params_mp@w >= species_params(params_mp)$w_mat[2]])))^2)
pbe_flounder_b = sum((log(flounder_ssb[flounder_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_b@n[(t_tune-length(prj_ts)+1):t_tune,2,], 2, params_b@w*params_b@dw, "*"))[,params_b@w >= species_params(params_b)$w_mat[2]])))^2)
pbe_flounder_o = sum((log(flounder_ssb[flounder_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_o@n[(t_tune-length(prj_ts)+1):t_tune,2,], 2, params_o@w*params_o@dw, "*"))[,params_o@w >= species_params(params_o)$w_mat[2]])))^2)
pbe_flounder_m = sum((log(flounder_ssb[flounder_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_m@n[(t_tune-length(prj_ts)+1):t_tune,2,], 2, params_m@w*params_m@dw, "*"))[,params_m@w >= species_params(params_m)$w_mat[2]])))^2)
pbe_flounder_p = sum((log(flounder_ssb[flounder_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_p@n[(t_tune-length(prj_ts)+1):t_tune,2,], 2, params_p@w*params_p@dw, "*"))[,params_p@w >= species_params(params_p)$w_mat[2]])))^2)
pbe_flounder_none = sum((log(flounder_ssb[flounder_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_none@n[(t_tune-length(prj_ts)+1):t_tune,2,], 2, params_none@w*params_none@dw, "*"))[,params_none@w >= species_params(params_none)$w_mat[2]])))^2)

pbe_sprat_bomp = sum((log(sprat_ssb[sprat_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_bomp@n[(t_tune-length(prj_ts)+1):t_tune,3,], 2, params_bomp@w*params_bomp@dw, "*"))[,params_bomp@w >= species_params(params_bomp)$w_mat[3]])))^2)
pbe_sprat_bom = sum((log(sprat_ssb[sprat_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_bom@n[(t_tune-length(prj_ts)+1):t_tune,3,], 2, params_bom@w*params_bom@dw, "*"))[,params_bom@w >= species_params(params_bom)$w_mat[3]])))^2)
pbe_sprat_bop = sum((log(sprat_ssb[sprat_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_bop@n[(t_tune-length(prj_ts)+1):t_tune,3,], 2, params_bop@w*params_bop@dw, "*"))[,params_bop@w >= species_params(params_bop)$w_mat[3]])))^2)
pbe_sprat_bmp = sum((log(sprat_ssb[sprat_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_bmp@n[(t_tune-length(prj_ts)+1):t_tune,3,], 2, params_bmp@w*params_bmp@dw, "*"))[,params_bmp@w >= species_params(params_bmp)$w_mat[3]])))^2)
pbe_sprat_omp = sum((log(sprat_ssb[sprat_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_omp@n[(t_tune-length(prj_ts)+1):t_tune,3,], 2, params_omp@w*params_omp@dw, "*"))[,params_omp@w >= species_params(params_omp)$w_mat[3]])))^2)
pbe_sprat_bo = sum((log(sprat_ssb[sprat_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_bo@n[(t_tune-length(prj_ts)+1):t_tune,3,], 2, params_bo@w*params_bo@dw, "*"))[,params_bo@w >= species_params(params_bo)$w_mat[3]])))^2)
pbe_sprat_bm = sum((log(sprat_ssb[sprat_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_bm@n[(t_tune-length(prj_ts)+1):t_tune,3,], 2, params_bm@w*params_bm@dw, "*"))[,params_bm@w >= species_params(params_bm)$w_mat[3]])))^2)
pbe_sprat_bp = sum((log(sprat_ssb[sprat_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_bp@n[(t_tune-length(prj_ts)+1):t_tune,3,], 2, params_bp@w*params_bp@dw, "*"))[,params_bp@w >= species_params(params_bp)$w_mat[3]])))^2)
pbe_sprat_om = sum((log(sprat_ssb[sprat_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_om@n[(t_tune-length(prj_ts)+1):t_tune,3,], 2, params_om@w*params_om@dw, "*"))[,params_om@w >= species_params(params_om)$w_mat[3]])))^2)
pbe_sprat_op = sum((log(sprat_ssb[sprat_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_op@n[(t_tune-length(prj_ts)+1):t_tune,3,], 2, params_op@w*params_op@dw, "*"))[,params_op@w >= species_params(params_op)$w_mat[3]])))^2)
pbe_sprat_mp = sum((log(sprat_ssb[sprat_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_mp@n[(t_tune-length(prj_ts)+1):t_tune,3,], 2, params_mp@w*params_mp@dw, "*"))[,params_mp@w >= species_params(params_mp)$w_mat[3]])))^2)
pbe_sprat_b = sum((log(sprat_ssb[sprat_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_b@n[(t_tune-length(prj_ts)+1):t_tune,3,], 2, params_b@w*params_b@dw, "*"))[,params_b@w >= species_params(params_b)$w_mat[3]])))^2)
pbe_sprat_o = sum((log(sprat_ssb[sprat_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_o@n[(t_tune-length(prj_ts)+1):t_tune,3,], 2, params_o@w*params_o@dw, "*"))[,params_o@w >= species_params(params_o)$w_mat[3]])))^2)
pbe_sprat_m = sum((log(sprat_ssb[sprat_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_m@n[(t_tune-length(prj_ts)+1):t_tune,3,], 2, params_m@w*params_m@dw, "*"))[,params_m@w >= species_params(params_m)$w_mat[3]])))^2)
pbe_sprat_p = sum((log(sprat_ssb[sprat_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_p@n[(t_tune-length(prj_ts)+1):t_tune,3,], 2, params_p@w*params_p@dw, "*"))[,params_p@w >= species_params(params_p)$w_mat[3]])))^2)
pbe_sprat_none = sum((log(sprat_ssb[sprat_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_none@n[(t_tune-length(prj_ts)+1):t_tune,3,], 2, params_none@w*params_none@dw, "*"))[,params_none@w >= species_params(params_none)$w_mat[3]])))^2)

pbe_herring_bomp = sum((log(herring_ssb[herring_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_bomp@n[(t_tune-length(prj_ts)+1):t_tune,4,], 2, params_bomp@w*params_bomp@dw, "*"))[,params_bomp@w >= species_params(params_bomp)$w_mat[4]])))^2)
pbe_herring_bom = sum((log(herring_ssb[herring_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_bom@n[(t_tune-length(prj_ts)+1):t_tune,4,], 2, params_bom@w*params_bom@dw, "*"))[,params_bom@w >= species_params(params_bom)$w_mat[4]])))^2)
pbe_herring_bop = sum((log(herring_ssb[herring_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_bop@n[(t_tune-length(prj_ts)+1):t_tune,4,], 2, params_bop@w*params_bop@dw, "*"))[,params_bop@w >= species_params(params_bop)$w_mat[4]])))^2)
pbe_herring_bmp = sum((log(herring_ssb[herring_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_bmp@n[(t_tune-length(prj_ts)+1):t_tune,4,], 2, params_bmp@w*params_bmp@dw, "*"))[,params_bmp@w >= species_params(params_bmp)$w_mat[4]])))^2)
pbe_herring_omp = sum((log(herring_ssb[herring_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_omp@n[(t_tune-length(prj_ts)+1):t_tune,4,], 2, params_omp@w*params_omp@dw, "*"))[,params_omp@w >= species_params(params_omp)$w_mat[4]])))^2)
pbe_herring_bo = sum((log(herring_ssb[herring_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_bo@n[(t_tune-length(prj_ts)+1):t_tune,4,], 2, params_bo@w*params_bo@dw, "*"))[,params_bo@w >= species_params(params_bo)$w_mat[4]])))^2)
pbe_herring_bm = sum((log(herring_ssb[herring_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_bm@n[(t_tune-length(prj_ts)+1):t_tune,4,], 2, params_bm@w*params_bm@dw, "*"))[,params_bm@w >= species_params(params_bm)$w_mat[4]])))^2)
pbe_herring_bp = sum((log(herring_ssb[herring_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_bp@n[(t_tune-length(prj_ts)+1):t_tune,4,], 2, params_bp@w*params_bp@dw, "*"))[,params_bp@w >= species_params(params_bp)$w_mat[4]])))^2)
pbe_herring_om = sum((log(herring_ssb[herring_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_om@n[(t_tune-length(prj_ts)+1):t_tune,4,], 2, params_om@w*params_om@dw, "*"))[,params_om@w >= species_params(params_om)$w_mat[4]])))^2)
pbe_herring_op = sum((log(herring_ssb[herring_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_op@n[(t_tune-length(prj_ts)+1):t_tune,4,], 2, params_op@w*params_op@dw, "*"))[,params_op@w >= species_params(params_op)$w_mat[4]])))^2)
pbe_herring_mp = sum((log(herring_ssb[herring_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_mp@n[(t_tune-length(prj_ts)+1):t_tune,4,], 2, params_mp@w*params_mp@dw, "*"))[,params_mp@w >= species_params(params_mp)$w_mat[4]])))^2)
pbe_herring_b = sum((log(herring_ssb[herring_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_b@n[(t_tune-length(prj_ts)+1):t_tune,4,], 2, params_b@w*params_b@dw, "*"))[,params_b@w >= species_params(params_b)$w_mat[4]])))^2)
pbe_herring_o = sum((log(herring_ssb[herring_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_o@n[(t_tune-length(prj_ts)+1):t_tune,4,], 2, params_o@w*params_o@dw, "*"))[,params_o@w >= species_params(params_o)$w_mat[4]])))^2)
pbe_herring_m = sum((log(herring_ssb[herring_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_m@n[(t_tune-length(prj_ts)+1):t_tune,4,], 2, params_m@w*params_m@dw, "*"))[,params_m@w >= species_params(params_m)$w_mat[4]])))^2)
pbe_herring_p = sum((log(herring_ssb[herring_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_p@n[(t_tune-length(prj_ts)+1):t_tune,4,], 2, params_p@w*params_p@dw, "*"))[,params_p@w >= species_params(params_p)$w_mat[4]])))^2)
pbe_herring_none = sum((log(herring_ssb[herring_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_none@n[(t_tune-length(prj_ts)+1):t_tune,4,], 2, params_none@w*params_none@dw, "*"))[,params_none@w >= species_params(params_none)$w_mat[4]])))^2)

# Sum across species
pbe_cod = c(pbe_cod_bomp, pbe_cod_bom, pbe_cod_bop, pbe_cod_bmp, pbe_cod_omp, pbe_cod_bo, pbe_cod_bm, pbe_cod_bp, pbe_cod_om, pbe_cod_op, pbe_cod_mp, pbe_cod_b, pbe_cod_o, pbe_cod_m, pbe_cod_p, pbe_cod_none)
pbe_flounder = c(pbe_flounder_bomp, pbe_flounder_bom, pbe_flounder_bop, pbe_flounder_bmp, pbe_flounder_omp, pbe_flounder_bo, pbe_flounder_bm, pbe_flounder_bp, pbe_flounder_om, pbe_flounder_op, pbe_flounder_mp, pbe_flounder_b, pbe_flounder_o, pbe_flounder_m, pbe_flounder_p, pbe_flounder_none)
pbe_sprat = c(pbe_sprat_bomp, pbe_sprat_bom, pbe_sprat_bop, pbe_sprat_bmp, pbe_sprat_omp, pbe_sprat_bo, pbe_sprat_bm, pbe_sprat_bp, pbe_sprat_om, pbe_sprat_op, pbe_sprat_mp, pbe_sprat_b, pbe_sprat_o, pbe_sprat_m, pbe_sprat_p, pbe_sprat_none)
pbe_herring = c(pbe_herring_bomp, pbe_herring_bom, pbe_herring_bop, pbe_herring_bmp, pbe_herring_omp, pbe_herring_bo, pbe_herring_bm, pbe_herring_bp, pbe_herring_om, pbe_herring_op, pbe_herring_mp, pbe_herring_b, pbe_herring_o, pbe_herring_m, pbe_herring_p, pbe_herring_none)

# Total and ranked error
pbe_cod_weight = 0.3
pbe_flounder_weight = 0.1 
pbe_sprat_weight = 0.3
pbe_herring_weight = 0.3
pbe_error = pbe_cod_weight*pbe_cod + pbe_flounder_weight*pbe_flounder + pbe_sprat_weight*pbe_sprat + pbe_herring_weight*pbe_herring
pbe_rank = (pbe_error - min(pbe_error)) / (max(pbe_error) - min(pbe_error))

# Yield projection error
yield_bomp = getYield(cal_bomp)
yield_bom = getYield(cal_bom)
yield_bop = getYield(cal_bop)
yield_bmp = getYield(cal_bmp)
yield_omp = getYield(cal_omp)
yield_bo = getYield(cal_bo)
yield_bm = getYield(cal_bm)
yield_bp = getYield(cal_bp)
yield_om = getYield(cal_om)
yield_op = getYield(cal_op)
yield_mp = getYield(cal_mp)
yield_b = getYield(cal_b)
yield_o = getYield(cal_o)
yield_m = getYield(cal_m)
yield_p = getYield(cal_p)
yield_none = getYield(cal_none)

cye_cod_bomp = sum((log(cod_catch[cod_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_bomp[(t_tune-length(cal_ts)+1):t_tune, "cod"]))^2)
cye_cod_bom = sum((log(cod_catch[cod_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_bom[(t_tune-length(cal_ts)+1):t_tune, "cod"]))^2)
cye_cod_bop = sum((log(cod_catch[cod_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_bop[(t_tune-length(cal_ts)+1):t_tune, "cod"]))^2)
cye_cod_bmp = sum((log(cod_catch[cod_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_bmp[(t_tune-length(cal_ts)+1):t_tune, "cod"]))^2)
cye_cod_omp = sum((log(cod_catch[cod_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_omp[(t_tune-length(cal_ts)+1):t_tune, "cod"]))^2)
cye_cod_bo = sum((log(cod_catch[cod_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_bo[(t_tune-length(cal_ts)+1):t_tune, "cod"]))^2)
cye_cod_bm = sum((log(cod_catch[cod_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_bm[(t_tune-length(cal_ts)+1):t_tune, "cod"]))^2)
cye_cod_bp = sum((log(cod_catch[cod_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_bp[(t_tune-length(cal_ts)+1):t_tune, "cod"]))^2)
cye_cod_om = sum((log(cod_catch[cod_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_om[(t_tune-length(cal_ts)+1):t_tune, "cod"]))^2)
cye_cod_op = sum((log(cod_catch[cod_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_op[(t_tune-length(cal_ts)+1):t_tune, "cod"]))^2)
cye_cod_mp = sum((log(cod_catch[cod_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_mp[(t_tune-length(cal_ts)+1):t_tune, "cod"]))^2)
cye_cod_b = sum((log(cod_catch[cod_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_b[(t_tune-length(cal_ts)+1):t_tune, "cod"]))^2)
cye_cod_o = sum((log(cod_catch[cod_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_o[(t_tune-length(cal_ts)+1):t_tune, "cod"]))^2)
cye_cod_m = sum((log(cod_catch[cod_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_m[(t_tune-length(cal_ts)+1):t_tune, "cod"]))^2)
cye_cod_p = sum((log(cod_catch[cod_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_p[(t_tune-length(cal_ts)+1):t_tune, "cod"]))^2)
cye_cod_none = sum((log(cod_catch[cod_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_none[(t_tune-length(cal_ts)+1):t_tune, "cod"]))^2)

cye_flounder_bomp = sum((log(flounder_catch[flounder_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_bomp[(t_tune-length(cal_ts)+1):t_tune, "flounder"]))^2)
cye_flounder_bom = sum((log(flounder_catch[flounder_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_bom[(t_tune-length(cal_ts)+1):t_tune, "flounder"]))^2)
cye_flounder_bop = sum((log(flounder_catch[flounder_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_bop[(t_tune-length(cal_ts)+1):t_tune, "flounder"]))^2)
cye_flounder_bmp = sum((log(flounder_catch[flounder_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_bmp[(t_tune-length(cal_ts)+1):t_tune, "flounder"]))^2)
cye_flounder_omp = sum((log(flounder_catch[flounder_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_omp[(t_tune-length(cal_ts)+1):t_tune, "flounder"]))^2)
cye_flounder_bo = sum((log(flounder_catch[flounder_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_bo[(t_tune-length(cal_ts)+1):t_tune, "flounder"]))^2)
cye_flounder_bm = sum((log(flounder_catch[flounder_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_bm[(t_tune-length(cal_ts)+1):t_tune, "flounder"]))^2)
cye_flounder_bp = sum((log(flounder_catch[flounder_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_bp[(t_tune-length(cal_ts)+1):t_tune, "flounder"]))^2)
cye_flounder_om = sum((log(flounder_catch[flounder_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_om[(t_tune-length(cal_ts)+1):t_tune, "flounder"]))^2)
cye_flounder_op = sum((log(flounder_catch[flounder_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_op[(t_tune-length(cal_ts)+1):t_tune, "flounder"]))^2)
cye_flounder_mp = sum((log(flounder_catch[flounder_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_mp[(t_tune-length(cal_ts)+1):t_tune, "flounder"]))^2)
cye_flounder_b = sum((log(flounder_catch[flounder_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_b[(t_tune-length(cal_ts)+1):t_tune, "flounder"]))^2)
cye_flounder_o = sum((log(flounder_catch[flounder_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_o[(t_tune-length(cal_ts)+1):t_tune, "flounder"]))^2)
cye_flounder_m = sum((log(flounder_catch[flounder_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_m[(t_tune-length(cal_ts)+1):t_tune, "flounder"]))^2)
cye_flounder_p = sum((log(flounder_catch[flounder_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_p[(t_tune-length(cal_ts)+1):t_tune, "flounder"]))^2)
cye_flounder_none = sum((log(flounder_catch[flounder_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_none[(t_tune-length(cal_ts)+1):t_tune, "flounder"]))^2)

cye_sprat_bomp = sum((log(sprat_catch[sprat_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_bomp[(t_tune-length(cal_ts)+1):t_tune, "sprat"]))^2)
cye_sprat_bom = sum((log(sprat_catch[sprat_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_bom[(t_tune-length(cal_ts)+1):t_tune, "sprat"]))^2)
cye_sprat_bop = sum((log(sprat_catch[sprat_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_bop[(t_tune-length(cal_ts)+1):t_tune, "sprat"]))^2)
cye_sprat_bmp = sum((log(sprat_catch[sprat_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_bmp[(t_tune-length(cal_ts)+1):t_tune, "sprat"]))^2)
cye_sprat_omp = sum((log(sprat_catch[sprat_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_omp[(t_tune-length(cal_ts)+1):t_tune, "sprat"]))^2)
cye_sprat_bo = sum((log(sprat_catch[sprat_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_bo[(t_tune-length(cal_ts)+1):t_tune, "sprat"]))^2)
cye_sprat_bm = sum((log(sprat_catch[sprat_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_bm[(t_tune-length(cal_ts)+1):t_tune, "sprat"]))^2)
cye_sprat_bp = sum((log(sprat_catch[sprat_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_bp[(t_tune-length(cal_ts)+1):t_tune, "sprat"]))^2)
cye_sprat_om = sum((log(sprat_catch[sprat_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_om[(t_tune-length(cal_ts)+1):t_tune, "sprat"]))^2)
cye_sprat_op = sum((log(sprat_catch[sprat_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_op[(t_tune-length(cal_ts)+1):t_tune, "sprat"]))^2)
cye_sprat_mp = sum((log(sprat_catch[sprat_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_mp[(t_tune-length(cal_ts)+1):t_tune, "sprat"]))^2)
cye_sprat_b = sum((log(sprat_catch[sprat_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_b[(t_tune-length(cal_ts)+1):t_tune, "sprat"]))^2)
cye_sprat_o = sum((log(sprat_catch[sprat_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_o[(t_tune-length(cal_ts)+1):t_tune, "sprat"]))^2)
cye_sprat_m = sum((log(sprat_catch[sprat_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_m[(t_tune-length(cal_ts)+1):t_tune, "sprat"]))^2)
cye_sprat_p = sum((log(sprat_catch[sprat_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_p[(t_tune-length(cal_ts)+1):t_tune, "sprat"]))^2)
cye_sprat_none = sum((log(sprat_catch[sprat_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_none[(t_tune-length(cal_ts)+1):t_tune, "sprat"]))^2)

cye_herring_bomp = sum((log(herring_catch[herring_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_bomp[(t_tune-length(cal_ts)+1):t_tune, "herring"]))^2)
cye_herring_bom = sum((log(herring_catch[herring_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_bom[(t_tune-length(cal_ts)+1):t_tune, "herring"]))^2)
cye_herring_bop = sum((log(herring_catch[herring_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_bop[(t_tune-length(cal_ts)+1):t_tune, "herring"]))^2)
cye_herring_bmp = sum((log(herring_catch[herring_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_bmp[(t_tune-length(cal_ts)+1):t_tune, "herring"]))^2)
cye_herring_omp = sum((log(herring_catch[herring_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_omp[(t_tune-length(cal_ts)+1):t_tune, "herring"]))^2)
cye_herring_bo = sum((log(herring_catch[herring_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_bo[(t_tune-length(cal_ts)+1):t_tune, "herring"]))^2)
cye_herring_bm = sum((log(herring_catch[herring_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_bm[(t_tune-length(cal_ts)+1):t_tune, "herring"]))^2)
cye_herring_bp = sum((log(herring_catch[herring_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_bp[(t_tune-length(cal_ts)+1):t_tune, "herring"]))^2)
cye_herring_om = sum((log(herring_catch[herring_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_om[(t_tune-length(cal_ts)+1):t_tune, "herring"]))^2)
cye_herring_op = sum((log(herring_catch[herring_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_op[(t_tune-length(cal_ts)+1):t_tune, "herring"]))^2)
cye_herring_mp = sum((log(herring_catch[herring_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_mp[(t_tune-length(cal_ts)+1):t_tune, "herring"]))^2)
cye_herring_b = sum((log(herring_catch[herring_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_b[(t_tune-length(cal_ts)+1):t_tune, "herring"]))^2)
cye_herring_o = sum((log(herring_catch[herring_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_o[(t_tune-length(cal_ts)+1):t_tune, "herring"]))^2)
cye_herring_m = sum((log(herring_catch[herring_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_m[(t_tune-length(cal_ts)+1):t_tune, "herring"]))^2)
cye_herring_p = sum((log(herring_catch[herring_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_p[(t_tune-length(cal_ts)+1):t_tune, "herring"]))^2)
cye_herring_none = sum((log(herring_catch[herring_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_none[(t_tune-length(cal_ts)+1):t_tune, "herring"]))^2)

# Sum across species
cye_cod = c(cye_cod_bomp, cye_cod_bom, cye_cod_bop, cye_cod_bmp, cye_cod_omp, cye_cod_bo, cye_cod_bm, cye_cod_bp, cye_cod_om, cye_cod_op, cye_cod_mp, cye_cod_b, cye_cod_o, cye_cod_m, cye_cod_p, cye_cod_none)
cye_flounder = c(cye_flounder_bomp, cye_flounder_bom, cye_flounder_bop, cye_flounder_bmp, cye_flounder_omp, cye_flounder_bo, cye_flounder_bm, cye_flounder_bp, cye_flounder_om, cye_flounder_op, cye_flounder_mp, cye_flounder_b, cye_flounder_o, cye_flounder_m, cye_flounder_p, cye_flounder_none)
cye_sprat = c(cye_sprat_bomp, cye_sprat_bom, cye_sprat_bop, cye_sprat_bmp, cye_sprat_omp, cye_sprat_bo, cye_sprat_bm, cye_sprat_bp, cye_sprat_om, cye_sprat_op, cye_sprat_mp, cye_sprat_b, cye_sprat_o, cye_sprat_m, cye_sprat_p, cye_sprat_none)
cye_herring = c(cye_herring_bomp, cye_herring_bom, cye_herring_bop, cye_herring_bmp, cye_herring_omp, cye_herring_bo, cye_herring_bm, cye_herring_bp, cye_herring_om, cye_herring_op, cye_herring_mp, cye_herring_b, cye_herring_o, cye_herring_m, cye_herring_p, cye_herring_none)

# Total and ranked error
cye_cod_weight = 0.25
cye_flounder_weight = 0.25 
cye_sprat_weight = 0.25
cye_herring_weight = 0.25
cye_error = cye_cod_weight*cye_cod + cye_flounder_weight*cye_flounder + cye_sprat_weight*cye_sprat + cye_herring_weight*cye_herring
cye_rank = (cye_error - min(cye_error)) / (max(cye_error) - min(cye_error))

# Yield projection error
yield_bomp = getYield(sim_bomp)
yield_bom = getYield(sim_bom)
yield_bop = getYield(sim_bop)
yield_bmp = getYield(sim_bmp)
yield_omp = getYield(sim_omp)
yield_bo = getYield(sim_bo)
yield_bm = getYield(sim_bm)
yield_bp = getYield(sim_bp)
yield_om = getYield(sim_om)
yield_op = getYield(sim_op)
yield_mp = getYield(sim_mp)
yield_b = getYield(sim_b)
yield_o = getYield(sim_o)
yield_m = getYield(sim_m)
yield_p = getYield(sim_p)
yield_none = getYield(sim_none)

pye_cod_bomp = sum((log(cod_catch[cod_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_bomp[(t_tune-length(prj_ts)+1):t_tune, "cod"]))^2)
pye_cod_bom = sum((log(cod_catch[cod_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_bom[(t_tune-length(prj_ts)+1):t_tune, "cod"]))^2)
pye_cod_bop = sum((log(cod_catch[cod_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_bop[(t_tune-length(prj_ts)+1):t_tune, "cod"]))^2)
pye_cod_bmp = sum((log(cod_catch[cod_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_bmp[(t_tune-length(prj_ts)+1):t_tune, "cod"]))^2)
pye_cod_omp = sum((log(cod_catch[cod_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_omp[(t_tune-length(prj_ts)+1):t_tune, "cod"]))^2)
pye_cod_bo = sum((log(cod_catch[cod_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_bo[(t_tune-length(prj_ts)+1):t_tune, "cod"]))^2)
pye_cod_bm = sum((log(cod_catch[cod_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_bm[(t_tune-length(prj_ts)+1):t_tune, "cod"]))^2)
pye_cod_bp = sum((log(cod_catch[cod_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_bp[(t_tune-length(prj_ts)+1):t_tune, "cod"]))^2)
pye_cod_om = sum((log(cod_catch[cod_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_om[(t_tune-length(prj_ts)+1):t_tune, "cod"]))^2)
pye_cod_op = sum((log(cod_catch[cod_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_op[(t_tune-length(prj_ts)+1):t_tune, "cod"]))^2)
pye_cod_mp = sum((log(cod_catch[cod_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_mp[(t_tune-length(prj_ts)+1):t_tune, "cod"]))^2)
pye_cod_b = sum((log(cod_catch[cod_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_b[(t_tune-length(prj_ts)+1):t_tune, "cod"]))^2)
pye_cod_o = sum((log(cod_catch[cod_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_o[(t_tune-length(prj_ts)+1):t_tune, "cod"]))^2)
pye_cod_m = sum((log(cod_catch[cod_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_m[(t_tune-length(prj_ts)+1):t_tune, "cod"]))^2)
pye_cod_p = sum((log(cod_catch[cod_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_p[(t_tune-length(prj_ts)+1):t_tune, "cod"]))^2)
pye_cod_none = sum((log(cod_catch[cod_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_none[(t_tune-length(prj_ts)+1):t_tune, "cod"]))^2)

pye_flounder_bomp = sum((log(flounder_catch[flounder_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_bomp[(t_tune-length(prj_ts)+1):t_tune, "flounder"]))^2)
pye_flounder_bom = sum((log(flounder_catch[flounder_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_bom[(t_tune-length(prj_ts)+1):t_tune, "flounder"]))^2)
pye_flounder_bop = sum((log(flounder_catch[flounder_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_bop[(t_tune-length(prj_ts)+1):t_tune, "flounder"]))^2)
pye_flounder_bmp = sum((log(flounder_catch[flounder_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_bmp[(t_tune-length(prj_ts)+1):t_tune, "flounder"]))^2)
pye_flounder_omp = sum((log(flounder_catch[flounder_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_omp[(t_tune-length(prj_ts)+1):t_tune, "flounder"]))^2)
pye_flounder_bo = sum((log(flounder_catch[flounder_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_bo[(t_tune-length(prj_ts)+1):t_tune, "flounder"]))^2)
pye_flounder_bm = sum((log(flounder_catch[flounder_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_bm[(t_tune-length(prj_ts)+1):t_tune, "flounder"]))^2)
pye_flounder_bp = sum((log(flounder_catch[flounder_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_bp[(t_tune-length(prj_ts)+1):t_tune, "flounder"]))^2)
pye_flounder_om = sum((log(flounder_catch[flounder_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_om[(t_tune-length(prj_ts)+1):t_tune, "flounder"]))^2)
pye_flounder_op = sum((log(flounder_catch[flounder_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_op[(t_tune-length(prj_ts)+1):t_tune, "flounder"]))^2)
pye_flounder_mp = sum((log(flounder_catch[flounder_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_mp[(t_tune-length(prj_ts)+1):t_tune, "flounder"]))^2)
pye_flounder_b = sum((log(flounder_catch[flounder_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_b[(t_tune-length(prj_ts)+1):t_tune, "flounder"]))^2)
pye_flounder_o = sum((log(flounder_catch[flounder_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_o[(t_tune-length(prj_ts)+1):t_tune, "flounder"]))^2)
pye_flounder_m = sum((log(flounder_catch[flounder_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_m[(t_tune-length(prj_ts)+1):t_tune, "flounder"]))^2)
pye_flounder_p = sum((log(flounder_catch[flounder_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_p[(t_tune-length(prj_ts)+1):t_tune, "flounder"]))^2)
pye_flounder_none = sum((log(flounder_catch[flounder_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_none[(t_tune-length(prj_ts)+1):t_tune, "flounder"]))^2)

pye_sprat_bomp = sum((log(sprat_catch[sprat_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_bomp[(t_tune-length(prj_ts)+1):t_tune, "sprat"]))^2)
pye_sprat_bom = sum((log(sprat_catch[sprat_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_bom[(t_tune-length(prj_ts)+1):t_tune, "sprat"]))^2)
pye_sprat_bop = sum((log(sprat_catch[sprat_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_bop[(t_tune-length(prj_ts)+1):t_tune, "sprat"]))^2)
pye_sprat_bmp = sum((log(sprat_catch[sprat_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_bmp[(t_tune-length(prj_ts)+1):t_tune, "sprat"]))^2)
pye_sprat_omp = sum((log(sprat_catch[sprat_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_omp[(t_tune-length(prj_ts)+1):t_tune, "sprat"]))^2)
pye_sprat_bo = sum((log(sprat_catch[sprat_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_bo[(t_tune-length(prj_ts)+1):t_tune, "sprat"]))^2)
pye_sprat_bm = sum((log(sprat_catch[sprat_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_bm[(t_tune-length(prj_ts)+1):t_tune, "sprat"]))^2)
pye_sprat_bp = sum((log(sprat_catch[sprat_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_bp[(t_tune-length(prj_ts)+1):t_tune, "sprat"]))^2)
pye_sprat_om = sum((log(sprat_catch[sprat_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_om[(t_tune-length(prj_ts)+1):t_tune, "sprat"]))^2)
pye_sprat_op = sum((log(sprat_catch[sprat_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_op[(t_tune-length(prj_ts)+1):t_tune, "sprat"]))^2)
pye_sprat_mp = sum((log(sprat_catch[sprat_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_mp[(t_tune-length(prj_ts)+1):t_tune, "sprat"]))^2)
pye_sprat_b = sum((log(sprat_catch[sprat_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_b[(t_tune-length(prj_ts)+1):t_tune, "sprat"]))^2)
pye_sprat_o = sum((log(sprat_catch[sprat_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_o[(t_tune-length(prj_ts)+1):t_tune, "sprat"]))^2)
pye_sprat_m = sum((log(sprat_catch[sprat_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_m[(t_tune-length(prj_ts)+1):t_tune, "sprat"]))^2)
pye_sprat_p = sum((log(sprat_catch[sprat_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_p[(t_tune-length(prj_ts)+1):t_tune, "sprat"]))^2)
pye_sprat_none = sum((log(sprat_catch[sprat_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_none[(t_tune-length(prj_ts)+1):t_tune, "sprat"]))^2)

pye_herring_bomp = sum((log(herring_catch[herring_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_bomp[(t_tune-length(prj_ts)+1):t_tune, "herring"]))^2)
pye_herring_bom = sum((log(herring_catch[herring_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_bom[(t_tune-length(prj_ts)+1):t_tune, "herring"]))^2)
pye_herring_bop = sum((log(herring_catch[herring_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_bop[(t_tune-length(prj_ts)+1):t_tune, "herring"]))^2)
pye_herring_bmp = sum((log(herring_catch[herring_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_bmp[(t_tune-length(prj_ts)+1):t_tune, "herring"]))^2)
pye_herring_omp = sum((log(herring_catch[herring_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_omp[(t_tune-length(prj_ts)+1):t_tune, "herring"]))^2)
pye_herring_bo = sum((log(herring_catch[herring_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_bo[(t_tune-length(prj_ts)+1):t_tune, "herring"]))^2)
pye_herring_bm = sum((log(herring_catch[herring_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_bm[(t_tune-length(prj_ts)+1):t_tune, "herring"]))^2)
pye_herring_bp = sum((log(herring_catch[herring_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_bp[(t_tune-length(prj_ts)+1):t_tune, "herring"]))^2)
pye_herring_om = sum((log(herring_catch[herring_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_om[(t_tune-length(prj_ts)+1):t_tune, "herring"]))^2)
pye_herring_op = sum((log(herring_catch[herring_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_op[(t_tune-length(prj_ts)+1):t_tune, "herring"]))^2)
pye_herring_mp = sum((log(herring_catch[herring_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_mp[(t_tune-length(prj_ts)+1):t_tune, "herring"]))^2)
pye_herring_b = sum((log(herring_catch[herring_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_b[(t_tune-length(prj_ts)+1):t_tune, "herring"]))^2)
pye_herring_o = sum((log(herring_catch[herring_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_o[(t_tune-length(prj_ts)+1):t_tune, "herring"]))^2)
pye_herring_m = sum((log(herring_catch[herring_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_m[(t_tune-length(prj_ts)+1):t_tune, "herring"]))^2)
pye_herring_p = sum((log(herring_catch[herring_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_p[(t_tune-length(prj_ts)+1):t_tune, "herring"]))^2)
pye_herring_none = sum((log(herring_catch[herring_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_none[(t_tune-length(prj_ts)+1):t_tune, "herring"]))^2)

# Sum across species
pye_cod = c(pye_cod_bomp, pye_cod_bom, pye_cod_bop, pye_cod_bmp, pye_cod_omp, pye_cod_bo, pye_cod_bm, pye_cod_bp, pye_cod_om, pye_cod_op, pye_cod_mp, pye_cod_b, pye_cod_o, pye_cod_m, pye_cod_p, pye_cod_none)
pye_flounder = c(pye_flounder_bomp, pye_flounder_bom, pye_flounder_bop, pye_flounder_bmp, pye_flounder_omp, pye_flounder_bo, pye_flounder_bm, pye_flounder_bp, pye_flounder_om, pye_flounder_op, pye_flounder_mp, pye_flounder_b, pye_flounder_o, pye_flounder_m, pye_flounder_p, pye_flounder_none)
pye_sprat = c(pye_sprat_bomp, pye_sprat_bom, pye_sprat_bop, pye_sprat_bmp, pye_sprat_omp, pye_sprat_bo, pye_sprat_bm, pye_sprat_bp, pye_sprat_om, pye_sprat_op, pye_sprat_mp, pye_sprat_b, pye_sprat_o, pye_sprat_m, pye_sprat_p, pye_sprat_none)
pye_herring = c(pye_herring_bomp, pye_herring_bom, pye_herring_bop, pye_herring_bmp, pye_herring_omp, pye_herring_bo, pye_herring_bm, pye_herring_bp, pye_herring_om, pye_herring_op, pye_herring_mp, pye_herring_b, pye_herring_o, pye_herring_m, pye_herring_p, pye_herring_none)

# Total and ranked error
pye_cod_weight = 0.25
pye_flounder_weight = 0.25 
pye_sprat_weight = 0.25
pye_herring_weight = 0.25
pye_error = pye_cod_weight*pye_cod + pye_flounder_weight*pye_flounder + pye_sprat_weight*pye_sprat + pye_herring_weight*pye_herring
pye_rank = (pye_error - min(pye_error)) / (max(pye_error) - min(pye_error))

# Growth calibration error
cge_cod_bomp = proj_growth(params_bomp, cal_ts, 100, "cod")
cge_cod_bom = proj_growth(params_bom, cal_ts, 100, "cod")
cge_cod_bop = proj_growth(params_bop, cal_ts, 100, "cod")
cge_cod_bmp = proj_growth(params_bmp, cal_ts, 100, "cod")
cge_cod_omp = proj_growth(params_omp, cal_ts, 100, "cod")
cge_cod_bo = proj_growth(params_bo, cal_ts, 100, "cod")
cge_cod_bm = proj_growth(params_bm, cal_ts, 100, "cod")
cge_cod_bp = proj_growth(params_bp, cal_ts, 100, "cod")
cge_cod_om = proj_growth(params_om, cal_ts, 100, "cod")
cge_cod_op = proj_growth(params_op, cal_ts, 100, "cod")
cge_cod_mp = proj_growth(params_mp, cal_ts, 100, "cod")
cge_cod_b = proj_growth(params_b, cal_ts, 100, "cod")
cge_cod_o = proj_growth(params_o, cal_ts, 100, "cod")
cge_cod_m = proj_growth(params_m, cal_ts, 100, "cod")
cge_cod_p = proj_growth(params_p, cal_ts, 100, "cod")
cge_cod_none = proj_growth(params_none, cal_ts, 100, "cod")

cge_flounder_bomp = proj_growth(params_bomp, cal_ts, 100, "flounder")
cge_flounder_bom = proj_growth(params_bom, cal_ts, 100, "flounder")
cge_flounder_bop = proj_growth(params_bop, cal_ts, 100, "flounder")
cge_flounder_bmp = proj_growth(params_bmp, cal_ts, 100, "flounder")
cge_flounder_omp = proj_growth(params_omp, cal_ts, 100, "flounder")
cge_flounder_bo = proj_growth(params_bo, cal_ts, 100, "flounder")
cge_flounder_bm = proj_growth(params_bm, cal_ts, 100, "flounder")
cge_flounder_bp = proj_growth(params_bp, cal_ts, 100, "flounder")
cge_flounder_om = proj_growth(params_om, cal_ts, 100, "flounder")
cge_flounder_op = proj_growth(params_op, cal_ts, 100, "flounder")
cge_flounder_mp = proj_growth(params_mp, cal_ts, 100, "flounder")
cge_flounder_b = proj_growth(params_b, cal_ts, 100, "flounder")
cge_flounder_o = proj_growth(params_o, cal_ts, 100, "flounder")
cge_flounder_m = proj_growth(params_m, cal_ts, 100, "flounder")
cge_flounder_p = proj_growth(params_p, cal_ts, 100, "flounder")
cge_flounder_none = proj_growth(params_none, cal_ts, 100, "flounder")

cge_sprat_bomp = proj_growth(params_bomp, cal_ts, 100, "sprat")
cge_sprat_bom = proj_growth(params_bom, cal_ts, 100, "sprat")
cge_sprat_bop = proj_growth(params_bop, cal_ts, 100, "sprat")
cge_sprat_bmp = proj_growth(params_bmp, cal_ts, 100, "sprat")
cge_sprat_omp = proj_growth(params_omp, cal_ts, 100, "sprat")
cge_sprat_bo = proj_growth(params_bo, cal_ts, 100, "sprat")
cge_sprat_bm = proj_growth(params_bm, cal_ts, 100, "sprat")
cge_sprat_bp = proj_growth(params_bp, cal_ts, 100, "sprat")
cge_sprat_om = proj_growth(params_om, cal_ts, 100, "sprat")
cge_sprat_op = proj_growth(params_op, cal_ts, 100, "sprat")
cge_sprat_mp = proj_growth(params_mp, cal_ts, 100, "sprat")
cge_sprat_b = proj_growth(params_b, cal_ts, 100, "sprat")
cge_sprat_o = proj_growth(params_o, cal_ts, 100, "sprat")
cge_sprat_m = proj_growth(params_m, cal_ts, 100, "sprat")
cge_sprat_p = proj_growth(params_p, cal_ts, 100, "sprat")
cge_sprat_none = proj_growth(params_none, cal_ts, 100, "sprat")

cge_herring_bomp = proj_growth(params_bomp, cal_ts, 100, "herring")
cge_herring_bom = proj_growth(params_bom, cal_ts, 100, "herring")
cge_herring_bop = proj_growth(params_bop, cal_ts, 100, "herring")
cge_herring_bmp = proj_growth(params_bmp, cal_ts, 100, "herring")
cge_herring_omp = proj_growth(params_omp, cal_ts, 100, "herring")
cge_herring_bo = proj_growth(params_bo, cal_ts, 100, "herring")
cge_herring_bm = proj_growth(params_bm, cal_ts, 100, "herring")
cge_herring_bp = proj_growth(params_bp, cal_ts, 100, "herring")
cge_herring_om = proj_growth(params_om, cal_ts, 100, "herring")
cge_herring_op = proj_growth(params_op, cal_ts, 100, "herring")
cge_herring_mp = proj_growth(params_mp, cal_ts, 100, "herring")
cge_herring_b = proj_growth(params_b, cal_ts, 100, "herring")
cge_herring_o = proj_growth(params_o, cal_ts, 100, "herring")
cge_herring_m = proj_growth(params_m, cal_ts, 100, "herring")
cge_herring_p = proj_growth(params_p, cal_ts, 100, "herring")
cge_herring_none = proj_growth(params_none, cal_ts, 100, "herring")

# Sum across species
cge_cod = c(cge_cod_bomp, cge_cod_bom, cge_cod_bop, cge_cod_bmp, cge_cod_omp, cge_cod_bo, cge_cod_bm, cge_cod_bp, cge_cod_om, cge_cod_op, cge_cod_mp, cge_cod_b, cge_cod_o, cge_cod_m, cge_cod_p, cge_cod_none)
cge_flounder = c(cge_flounder_bomp, cge_flounder_bom, cge_flounder_bop, cge_flounder_bmp, cge_flounder_omp, cge_flounder_bo, cge_flounder_bm, cge_flounder_bp, cge_flounder_om, cge_flounder_op, cge_flounder_mp, cge_flounder_b, cge_flounder_o, cge_flounder_m, cge_flounder_p, cge_flounder_none)
cge_sprat = c(cge_sprat_bomp, cge_sprat_bom, cge_sprat_bop, cge_sprat_bmp, cge_sprat_omp, cge_sprat_bo, cge_sprat_bm, cge_sprat_bp, cge_sprat_om, cge_sprat_op, cge_sprat_mp, cge_sprat_b, cge_sprat_o, cge_sprat_m, cge_sprat_p, cge_sprat_none)
cge_herring = c(cge_herring_bomp, cge_herring_bom, cge_herring_bop, cge_herring_bmp, cge_herring_omp, cge_herring_bo, cge_herring_bm, cge_herring_bp, cge_herring_om, cge_herring_op, cge_herring_mp, cge_herring_b, cge_herring_o, cge_herring_m, cge_herring_p, cge_herring_none)

# Total and ranked error
cge_cod_weight = 0.25
cge_flounder_weight = 0.25
cge_sprat_weight = 0.25
cge_herring_weight = 0.25
cge_error = cge_cod_weight*cge_cod + cge_flounder_weight*cge_flounder + cge_sprat_weight*cge_sprat + cge_herring_weight*cge_herring
cge_rank = (cge_error - min(cge_error)) / (max(cge_error) - min(cge_error))

# Growth projection error
pge_cod_bomp = proj_growth(params_bomp, prj_ts, 100, "cod")
pge_cod_bom = proj_growth(params_bom, prj_ts, 100, "cod")
pge_cod_bop = proj_growth(params_bop, prj_ts, 100, "cod")
pge_cod_bmp = proj_growth(params_bmp, prj_ts, 100, "cod")
pge_cod_omp = proj_growth(params_omp, prj_ts, 100, "cod")
pge_cod_bo = proj_growth(params_bo, prj_ts, 100, "cod")
pge_cod_bm = proj_growth(params_bm, prj_ts, 100, "cod")
pge_cod_bp = proj_growth(params_bp, prj_ts, 100, "cod")
pge_cod_om = proj_growth(params_om, prj_ts, 100, "cod")
pge_cod_op = proj_growth(params_op, prj_ts, 100, "cod")
pge_cod_mp = proj_growth(params_mp, prj_ts, 100, "cod")
pge_cod_b = proj_growth(params_b, prj_ts, 100, "cod")
pge_cod_o = proj_growth(params_o, prj_ts, 100, "cod")
pge_cod_m = proj_growth(params_m, prj_ts, 100, "cod")
pge_cod_p = proj_growth(params_p, prj_ts, 100, "cod")
pge_cod_none = proj_growth(params_none, prj_ts, 100, "cod")

pge_flounder_bomp = proj_growth(params_bomp, prj_ts, 100, "flounder")
pge_flounder_bom = proj_growth(params_bom, prj_ts, 100, "flounder")
pge_flounder_bop = proj_growth(params_bop, prj_ts, 100, "flounder")
pge_flounder_bmp = proj_growth(params_bmp, prj_ts, 100, "flounder")
pge_flounder_omp = proj_growth(params_omp, prj_ts, 100, "flounder")
pge_flounder_bo = proj_growth(params_bo, prj_ts, 100, "flounder")
pge_flounder_bm = proj_growth(params_bm, prj_ts, 100, "flounder")
pge_flounder_bp = proj_growth(params_bp, prj_ts, 100, "flounder")
pge_flounder_om = proj_growth(params_om, prj_ts, 100, "flounder")
pge_flounder_op = proj_growth(params_op, prj_ts, 100, "flounder")
pge_flounder_mp = proj_growth(params_mp, prj_ts, 100, "flounder")
pge_flounder_b = proj_growth(params_b, prj_ts, 100, "flounder")
pge_flounder_o = proj_growth(params_o, prj_ts, 100, "flounder")
pge_flounder_m = proj_growth(params_m, prj_ts, 100, "flounder")
pge_flounder_p = proj_growth(params_p, prj_ts, 100, "flounder")
pge_flounder_none = proj_growth(params_none, prj_ts, 100, "flounder")

pge_sprat_bomp = proj_growth(params_bomp, prj_ts, 100, "sprat")
pge_sprat_bom = proj_growth(params_bom, prj_ts, 100, "sprat")
pge_sprat_bop = proj_growth(params_bop, prj_ts, 100, "sprat")
pge_sprat_bmp = proj_growth(params_bmp, prj_ts, 100, "sprat")
pge_sprat_omp = proj_growth(params_omp, prj_ts, 100, "sprat")
pge_sprat_bo = proj_growth(params_bo, prj_ts, 100, "sprat")
pge_sprat_bm = proj_growth(params_bm, prj_ts, 100, "sprat")
pge_sprat_bp = proj_growth(params_bp, prj_ts, 100, "sprat")
pge_sprat_om = proj_growth(params_om, prj_ts, 100, "sprat")
pge_sprat_op = proj_growth(params_op, prj_ts, 100, "sprat")
pge_sprat_mp = proj_growth(params_mp, prj_ts, 100, "sprat")
pge_sprat_b = proj_growth(params_b, prj_ts, 100, "sprat")
pge_sprat_o = proj_growth(params_o, prj_ts, 100, "sprat")
pge_sprat_m = proj_growth(params_m, prj_ts, 100, "sprat")
pge_sprat_p = proj_growth(params_p, prj_ts, 100, "sprat")
pge_sprat_none = proj_growth(params_none, prj_ts, 100, "sprat")

pge_herring_bomp = proj_growth(params_bomp, prj_ts, 100, "herring")
pge_herring_bom = proj_growth(params_bom, prj_ts, 100, "herring")
pge_herring_bop = proj_growth(params_bop, prj_ts, 100, "herring")
pge_herring_bmp = proj_growth(params_bmp, prj_ts, 100, "herring")
pge_herring_omp = proj_growth(params_omp, prj_ts, 100, "herring")
pge_herring_bo = proj_growth(params_bo, prj_ts, 100, "herring")
pge_herring_bm = proj_growth(params_bm, prj_ts, 100, "herring")
pge_herring_bp = proj_growth(params_bp, prj_ts, 100, "herring")
pge_herring_om = proj_growth(params_om, prj_ts, 100, "herring")
pge_herring_op = proj_growth(params_op, prj_ts, 100, "herring")
pge_herring_mp = proj_growth(params_mp, prj_ts, 100, "herring")
pge_herring_b = proj_growth(params_b, prj_ts, 100, "herring")
pge_herring_o = proj_growth(params_o, prj_ts, 100, "herring")
pge_herring_m = proj_growth(params_m, prj_ts, 100, "herring")
pge_herring_p = proj_growth(params_p, prj_ts, 100, "herring")
pge_herring_none = proj_growth(params_none, prj_ts, 100, "herring")

# Sum across species
pge_cod = c(pge_cod_bomp, pge_cod_bom, pge_cod_bop, pge_cod_bmp, pge_cod_omp, pge_cod_bo, pge_cod_bm, pge_cod_bp, pge_cod_om, pge_cod_op, pge_cod_mp, pge_cod_b, pge_cod_o, pge_cod_m, pge_cod_p, pge_cod_none)
pge_flounder = c(pge_flounder_bomp, pge_flounder_bom, pge_flounder_bop, pge_flounder_bmp, pge_flounder_omp, pge_flounder_bo, pge_flounder_bm, pge_flounder_bp, pge_flounder_om, pge_flounder_op, pge_flounder_mp, pge_flounder_b, pge_flounder_o, pge_flounder_m, pge_flounder_p, pge_flounder_none)
pge_sprat = c(pge_sprat_bomp, pge_sprat_bom, pge_sprat_bop, pge_sprat_bmp, pge_sprat_omp, pge_sprat_bo, pge_sprat_bm, pge_sprat_bp, pge_sprat_om, pge_sprat_op, pge_sprat_mp, pge_sprat_b, pge_sprat_o, pge_sprat_m, pge_sprat_p, pge_sprat_none)
pge_herring = c(pge_herring_bomp, pge_herring_bom, pge_herring_bop, pge_herring_bmp, pge_herring_omp, pge_herring_bo, pge_herring_bm, pge_herring_bp, pge_herring_om, pge_herring_op, pge_herring_mp, pge_herring_b, pge_herring_o, pge_herring_m, pge_herring_p, pge_herring_none)

# Total and ranked error
pge_cod_weight = 0.25
pge_flounder_weight = 0.25 
pge_sprat_weight = 0.25
pge_herring_weight = 0.25
pge_error = pge_cod_weight*pge_cod + pge_flounder_weight*pge_flounder + pge_sprat_weight*pge_sprat + pge_herring_weight*pge_herring
pge_rank = (pge_error - min(pge_error)) / (max(pge_error) - min(pge_error))

# Calculate final weighted error
ssb_weight = 0.3
yield_weight = 0.2
growth_weight = 0.4
final_error_cal = ssb_weight*cbe_rank + yield_weight*cye_rank + growth_weight*cge_rank
final_error_sim = ssb_weight*pbe_rank + yield_weight*pye_rank + growth_weight*pge_rank
model_labels[order(final_error_cal)]
model_labels[order(final_error_sim)]


########## XVI. Plot calibration error ##########

# Plot colors
color_cod = rgb(0, 166, 174, maxColorValue = 255)
color_flounder = rgb(252, 255, 221, maxColorValue = 255)
color_sprat = rgb(205, 240, 203, maxColorValue = 255)
color_herring = rgb(119, 209, 181, maxColorValue = 255)

# Axis label colors
color = rep(c("#D55E00","#0072B2","#F0E442","#CC79A7","#E69F00","#009E73","#000000","#D4D4D4"),2)

# Initialize plot
jpeg("Plots/calibration_error.jpg", width = 30, height = 13, units = 'cm', res = 600)

# Set up plot
par(mfrow = c(1,3), mar = c(2,2,2,2), oma = c(9,5,3,1))

# SSB error
plot(cbe_error[order(cbe_error)], cex = 2, pch = 25, bg = "black", ylim = c(0,1), xlab = "", ylab = "", cex.main = 2, main = "ln(SSB)", axes = F); axis(2, cex.axis = 1.5, at = seq(0,1,0.2)); box()
for(i in seq(16)) axis(1, cex.axis = 1.5, at = i, labels = model_labels[order(cbe_error)][i], col.axis = color[order(cbe_error)][i], las = 2)
points(seq(16), 0.3*cbe_cod[order(cbe_cod)], cex = 2, pch = 21, bg = color_cod)
points(seq(16), 0.1*cbe_flounder[order(cbe_flounder)], cex = 2, pch = 22, bg = color_flounder)
points(seq(16), 0.3*cbe_sprat[order(cbe_sprat)], cex = 2, pch = 23, bg = color_sprat)
points(seq(16), 0.3*cbe_herring[order(cbe_herring)], cex = 2, pch = 24, bg = color_herring)

# Catch error
plot(cye_error[order(cye_error)], cex = 2, pch = 25, bg = "black", ylim = c(0,4), xlab = "", ylab = "", cex.main = 2, main = "ln(Catch)", axes = F); axis(2, cex.axis = 1.5, at = seq(0,4,1)); box()
for(i in seq(16)) axis(1, cex.axis = 1.5, at = i, labels = model_labels[order(cye_error)][i], col.axis = color[order(cye_error)][i], las = 2)
points(seq(16), 0.25*cye_cod[order(cye_cod)], cex = 2, pch = 21, bg = color_cod)
points(seq(16), 0.25*cye_flounder[order(cye_flounder)], cex = 2, pch = 22, bg = color_flounder)
points(seq(16), 0.25*cye_sprat[order(cye_sprat)], cex = 2, pch = 23, bg = color_sprat)
points(seq(16), 0.25*cye_herring[order(cye_herring)], cex = 2, pch = 24, bg = color_herring)

# Growth error
plot(cge_error[order(cge_error)], cex = 2, pch = 25, bg = "black", ylim = c(0,80000), xlab = "", ylab = "", cex.main = 2, main = "Growth (% max)", axes = F); axis(2, cex.axis = 1.5, at = seq(0,80000,20000)); box()
for(i in seq(16)) axis(1, cex.axis = 1.5, at = i, labels = model_labels[order(cge_error)][i], col.axis = color[order(cge_error)][i], las = 2)
points(seq(16), 0.25*cge_cod[order(cge_cod)], cex = 2, pch = 21, bg = color_cod)
points(seq(16), 0.25*cge_flounder[order(cge_flounder)], cex = 2, pch = 22, bg = color_flounder)
points(seq(16), 0.25*cge_sprat[order(cge_sprat)], cex = 2, pch = 23, bg = color_sprat)
points(seq(16), 0.25*cge_herring[order(cge_herring)], cex = 2, pch = 24, bg = color_herring)

# Error label
mtext("Weighted error", 2, line = 2, outer = T, cex = 1.5)

# Legend
par(fig = c(0,1,0,1), oma = c(0,0,0,0), mar = c(0,0,0,0), new = T)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend("bottom", bty = 'n', c("Cod","Flounder","Sprat","Herring","Total"), pt.bg = c(color_cod, color_flounder, color_sprat, color_herring, "black"), pch = c(21,22,23,24,25), pt.cex = 2, xpd = TRUE, horiz = TRUE, cex = 1.5, seg.len=1)

# Finish plot
dev.off()


########## XVII. Plot calibration error rank ##########

# Initialize plot
jpeg("Plots/calibration_error_rank.jpg", width = 30, height = 13, units = 'cm', res = 600)

# Set up plot
par(mfrow = c(1,3), mar = c(2,2,2,2), oma = c(9,5,3,1))

# SSB error
plot(cbe_rank[order(cbe_rank)], cex = 2, pch = 25, bg = "black", ylim = c(0,1), xlab = "", ylab = "", cex.main = 2, main = "ln(SSB)", axes = F); axis(2, cex.axis = 1.5, at = seq(0,1,0.2)); box()
for(i in seq(16)) axis(1, cex.axis = 1.5, at = i, labels = model_labels[order(cbe_rank)][i], col.axis = color[order(cbe_rank)][i], las = 2)

# Catch error
plot(cye_rank[order(cye_rank)], cex = 2, pch = 25, bg = "black", ylim = c(0,1), xlab = "", ylab = "", cex.main = 2, main = "ln(Catch)", axes = F); axis(2, cex.axis = 1.5, at = seq(0,1,0.2)); box()
for(i in seq(16)) axis(1, cex.axis = 1.5, at = i, labels = model_labels[order(cye_rank)][i], col.axis = color[order(cye_rank)][i], las = 2)

# Growth error
plot(cge_rank[order(cge_rank)], cex = 2, pch = 25, bg = "black", ylim = c(0,1), xlab = "", ylab = "", cex.main = 2, main = "Growth (% max)", axes = F); axis(2, cex.axis = 1.5, at = seq(0,1,0.2)); box()
for(i in seq(16)) axis(1, cex.axis = 1.5, at = i, labels = model_labels[order(cge_rank)][i], col.axis = color[order(cge_rank)][i], las = 2)

# Error label
mtext("Weighted error (normalized)", 2, line = 2, outer = T, cex = 1.5)

# Finish plot
dev.off()


########## XVIII. Plot projection error ##########

# Plot colors
color_cod = rgb(0, 166, 174, maxColorValue = 255)
color_flounder = rgb(252, 255, 221, maxColorValue = 255)
color_sprat = rgb(205, 240, 203, maxColorValue = 255)
color_herring = rgb(119, 209, 181, maxColorValue = 255)

# Axis label colors
color = rep(c("#D55E00","#0072B2","#F0E442","#CC79A7","#E69F00","#009E73","#000000","#D4D4D4"),2)

# Initialize plot
jpeg("Plots/simulation_error.jpg", width = 30, height = 13, units = 'cm', res = 600)

# Set up plot
par(mfrow = c(1,3), mar = c(2,2,2,2), oma = c(9,5,3,1))

# SSB error
plot(pbe_error[order(pbe_error)], cex = 2, pch = 25, bg = "black", ylim = c(0,100), xlab = "", ylab = "", cex.main = 2, main = "ln(SSB)", axes = F); axis(2, cex.axis = 1.5, at = seq(0,100,25)); box()
for(i in seq(16)) axis(1, cex.axis = 1.5, at = i, labels = model_labels[order(pbe_error)][i], col.axis = color[order(pbe_error)][i], las = 2)
points(seq(16), 0.3*pbe_cod[order(pbe_cod)], cex = 2, pch = 21, bg = color_cod)
points(seq(16), 0.1*pbe_flounder[order(pbe_flounder)], cex = 2, pch = 22, bg = color_flounder)
points(seq(16), 0.3*pbe_sprat[order(pbe_sprat)], cex = 2, pch = 23, bg = color_sprat)
points(seq(16), 0.3*pbe_herring[order(pbe_herring)], cex = 2, pch = 24, bg = color_herring)

# Catch error
plot(pye_error[order(pye_error)], cex = 2, pch = 25, bg = "black", ylim = c(0,80), xlab = "", ylab = "", cex.main = 2, main = "ln(Catch)", axes = F); axis(2, cex.axis = 1.5, at = seq(0,80,20)); box()
for(i in seq(16)) axis(1, cex.axis = 1.5, at = i, labels = model_labels[order(pye_error)][i], col.axis = color[order(pye_error)][i], las = 2)
points(seq(16), 0.25*pye_cod[order(pye_cod)], cex = 2, pch = 21, bg = color_cod)
points(seq(16), 0.25*pye_flounder[order(pye_flounder)], cex = 2, pch = 22, bg = color_flounder)
points(seq(16), 0.25*pye_sprat[order(pye_sprat)], cex = 2, pch = 23, bg = color_sprat)
points(seq(16), 0.25*pye_herring[order(pye_herring)], cex = 2, pch = 24, bg = color_herring)

# Growth error
plot(pge_error[order(pge_error)], cex = 2, pch = 25, bg = "black", ylim = c(0,200000), xlab = "", ylab = "", cex.main = 2, main = "Growth (% max; 1000's)", axes = F); axis(2, cex.axis = 1.5, at = seq(0,200000,50000), labels = seq(0,200,50)); box()
for(i in seq(16)) axis(1, cex.axis = 1.5, at = i, labels = model_labels[order(pge_error)][i], col.axis = color[order(pge_error)][i], las = 2)
points(seq(16), 0.25*pge_cod[order(pge_cod)], cex = 2, pch = 21, bg = color_cod)
points(seq(16), 0.25*pge_flounder[order(pge_flounder)], cex = 2, pch = 22, bg = color_flounder)
points(seq(16), 0.25*pge_sprat[order(pge_sprat)], cex = 2, pch = 23, bg = color_sprat)
points(seq(16), 0.25*pge_herring[order(pge_herring)], cex = 2, pch = 24, bg = color_herring)

# Error label
mtext("Weighted error", 2, line = 2, outer = T, cex = 1.5)

# Legend
par(fig = c(0,1,0,1), oma = c(0,0,0,0), mar = c(0,0,0,0), new = T)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend("bottom", bty = 'n', c("Cod", "Flounder", "Sprat", "Herring", "Total"), pt.bg = c(color_cod, color_flounder, color_sprat, color_herring, "black"), pch = c(21,22,23,24,25), pt.cex = 2, xpd = TRUE, horiz = TRUE, cex = 1.5, seg.len=1)

# Finish plot
dev.off()


########## XIX. Plot projection error rank ##########

# Initialize plot
jpeg("Plots/simulation_error_rank.jpg", width = 30, height = 13, units = 'cm', res = 600)

# Set up plot
par(mfrow = c(1,3), mar = c(2,2,2,2), oma = c(9,5,3,1))

# SSB error
plot(pbe_rank[order(pbe_rank)], cex = 2, pch = 25, bg = "black", ylim = c(0,1), xlab = "", ylab = "", cex.main = 2, main = "ln(SSB)", axes = F); axis(2, cex.axis = 1.5, at = seq(0,1,0.2)); box()
for(i in seq(16)) axis(1, cex.axis = 1.5, at = i, labels = model_labels[order(pbe_rank)][i], col.axis = color[order(pbe_rank)][i], las = 2)

# Catch error
plot(pye_rank[order(pye_rank)], cex = 2, pch = 25, bg = "black", ylim = c(0,1), xlab = "", ylab = "", cex.main = 2, main = "ln(Catch)", axes = F); axis(2, cex.axis = 1.5, at = seq(0,1,0.2)); box()
for(i in seq(16)) axis(1, cex.axis = 1.5, at = i, labels = model_labels[order(pye_rank)][i], col.axis = color[order(pye_rank)][i], las = 2)

# Growth error
plot(pge_rank[order(pge_rank)], cex = 2, pch = 25, bg = "black", ylim = c(0,1), xlab = "", ylab = "", cex.main = 2, main = "Growth (% max)", axes = F); axis(2, cex.axis = 1.5, at = seq(0,1,0.2)); box()
for(i in seq(16)) axis(1, cex.axis = 1.5, at = i, labels = model_labels[order(pge_rank)][i], col.axis = color[order(pge_rank)][i], las = 2)

# Error label
mtext("Weighted error (normalized)", 2, line = 2, outer = T, cex = 1.5)

# Finish plot
dev.off()


########## XX. Plot final weighted error rank for calibration ##########

# Plot colors
error_colors = palette.colors(palette = "Okabe-Ito")[1:4]

# Initialize plot
jpeg("Plots/final_error_cal.jpg", width = 13, height = 15, units = 'cm', res = 600)
par(oma = c(3,1,0,0))

# Plot calibration
plot(final_error_cal[order(final_error_cal)], pch = 21, bg = error_colors[1], cex = 1.5, ylim = c(0,1), xlab = "", ylab = "", main = "Calibration", cex.main = 1.5, axes = F); axis(2); box()
for(i in seq(16)) axis(1, at = i, labels = model_labels[order(final_error_cal)][i], col.axis = color[order(final_error_cal)][i], las = 2)
points(ssb_weight*cbe_rank[order(final_error_cal)], pch = 22, bg = error_colors[2], cex = 1.5)
points(yield_weight*cye_rank[order(final_error_cal)], pch = 23, bg = error_colors[3], cex = 1.5)
points(growth_weight*cge_rank[order(final_error_cal)], pch = 24, bg = error_colors[4], cex = 1.5)

# Y-axis label
mtext("Weighted error", 2, cex = 1.5, line = -1, outer = T)

# Legend
par(fig = c(0,1,0,1), oma = c(0,1,0,0), mar = c(0,0,0,0), new = T)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend("bottom", bty = 'n', c("SSB", "Yield", "Growth", "Total"), pt.bg = error_colors[c(2,3,4,1)], pch = c(22,23,24,21), pt.cex = 1.5, xpd = TRUE, horiz = TRUE, cex = 1.2, seg.len=1)

# Finish plot
dev.off()


########## XXI. Plot final weighted error rank for projection ##########

# Initialize plot
jpeg("Plots/final_error_sim.jpg", width = 13, height = 15, units = 'cm', res = 600)
par(oma = c(3,1,0,0))

# Plot projection
plot(final_error_sim[order(final_error_sim)], pch = 21, bg = error_colors[1], cex = 1.5, ylim = c(0,1), xlab = "", ylab = "", main = "Projection", cex.main = 1.5, axes = F); axis(2); box()
for(i in seq(16)) axis(1, at = i, labels = model_labels[order(final_error_sim)][i], col.axis = color[order(final_error_sim)][i], las = 2)
points(ssb_weight*pbe_rank[order(final_error_sim)], pch = 22, bg = error_colors[2], cex = 1.5)
points(yield_weight*pye_rank[order(final_error_sim)], pch = 23, bg = error_colors[3], cex = 1.5)
points(growth_weight*pge_rank[order(final_error_sim)], pch = 24, bg = error_colors[4], cex = 1.5)

# Y-axis label
mtext("Weighted error", 2, cex = 1.5, line = -1, outer = T)

# Legend
par(fig = c(0,1,0,1), oma = c(0,1,0,0), mar = c(0,0,0,0), new = T)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend("bottom", bty = 'n', c("SSB", "Yield", "Growth", "Total"), pt.bg = error_colors[c(2,3,4,1)], pch = c(22,23,24,21), pt.cex = 1.5, xpd = TRUE, horiz = TRUE, cex = 1.2, seg.len=1)

# Finish plot
dev.off()