# Load and plot calibration results

# Created: February 24, 2022
# Last modified: January 4, 2023

# Set working directory
setwd(paste(mypath, "hypoxiaSSM-main", sep = ""))

# Load packages
library(mizer)
library(jpeg)
library(png)
library(xtable)
library(diagram)
library(shape)
library(gridExtra)
library(grid)
library(ggpubr)
library(colourvalues)
library(fields)

# Contents (ctrl-f):
#	I. Source scripts
#	II. Load results
#	III. Run simulations
#	IV. Results table
#	IV. Plot results table
#	V. Get occupancy
#	VI. Get diet
#	VII. Plot food web diagrams
#	VIII. Plot diet changes


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
#	scenarios: run ideal and poor oxygen scenarios for a given set of parameters
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
source("Code/scenarios.R")


########## II. Load results ##########

# Full model
load("Calibration/params_full_bomp.rda")
params_bomp = params_full

# Full model
load("Calibration/params_full_bom.rda")
params_bom = params_full

# Full model
load("Calibration/params_full_bop.rda")
params_bop = params_full

# Full model
load("Calibration/params_full_bmp.rda")
params_bmp = params_full

# Full model
load("Calibration/params_full_omp.rda")
params_omp = params_full

# Occupancy and physiological scaling
load("Calibration/params_full_bo.rda")
params_bo = params_full

# Occupancy and physiological scaling
load("Calibration/params_full_bm.rda")
params_bm = params_full

# Occupancy and physiological scaling
load("Calibration/params_full_bp.rda")
params_bp = params_full

# Occupancy and physiological scaling
load("Calibration/params_full_om.rda")
params_om = params_full

# Occupancy and physiological scaling
load("Calibration/params_full_op.rda")
params_op = params_full

# Occupancy and physiological scaling
load("Calibration/params_full_mp.rda")
params_mp = params_full

# Occupancy scaling
load("Calibration/params_full_b.rda")
params_b = params_full

# Physiological scaling
load("Calibration/params_full_o.rda")
params_o = params_full

# Physiological scaling
load("Calibration/params_full_m.rda")
params_m = params_full

# Physiological scaling
load("Calibration/params_full_p.rda")
params_p = params_full

# None scaling (with left fish)
load("Calibration/params_full_none.rda")
params_none = params_full


########## III. Run simulations ##########

# Choose time
t_max = 500

# # Use scenarios function to scale objects, and run simulations for best model and components
# sim_bom = oxy_scenario(params_bom, t = t_max, p = F)
# sim_bo = oxy_scenario(params_bo, t = t_max, p = F)
# sim_bm = oxy_scenario(params_bm, t = t_max, p = F)
# sim_om = oxy_scenario(params_om, t = t_max, p = F)
# sim_b = oxy_scenario(params_b, t = t_max, p = F)
# sim_o = oxy_scenario(params_o, t = t_max, p = F)
# sim_m = oxy_scenario(params_m, t = t_max, p = F)
# sim_none = oxy_scenario(params_none, t = t_max, p = F)

# # Use scenarios function to scale objects, and run simulations for all other models
# sim_bomp = oxy_scenario(params_bomp, t = t_max, p = T)
# sim_bop = oxy_scenario(params_bop, t = t_max, p = T)
# sim_bmp = oxy_scenario(params_bmp, t = t_max, p = T)
# sim_omp = oxy_scenario(params_omp, t = t_max, p = T)
# sim_bp = oxy_scenario(params_bp, t = t_max, p = T)
# sim_op = oxy_scenario(params_op, t = t_max, p = T)
# sim_mp = oxy_scenario(params_mp, t = t_max, p = T)
# sim_p = oxy_scenario(params_p, t = t_max, p = T)

# # Save all scenario output
# save(object = sim_bom, file = "Scenarios/sim_bom.rda")
# save(object = sim_bo, file = "Scenarios/sim_bo.rda")
# save(object = sim_bm, file = "Scenarios/sim_bm.rda")
# save(object = sim_om, file = "Scenarios/sim_om.rda")
# save(object = sim_b, file = "Scenarios/sim_b.rda")
# save(object = sim_o, file = "Scenarios/sim_o.rda")
# save(object = sim_m, file = "Scenarios/sim_m.rda")
# save(object = sim_none, file = "Scenarios/sim_none.rda")
# save(object = sim_bomp, file = "Scenarios/sim_bomp.rda")
# save(object = sim_bop, file = "Scenarios/sim_bop.rda")
# save(object = sim_bmp, file = "Scenarios/sim_bmp.rda")
# save(object = sim_omp, file = "Scenarios/sim_omp.rda")
# save(object = sim_bp, file = "Scenarios/sim_bp.rda")
# save(object = sim_op, file = "Scenarios/sim_op.rda")
# save(object = sim_mp, file = "Scenarios/sim_mp.rda")
# save(object = sim_p, file = "Scenarios/sim_p.rda")

# Load scenarios for best model, components of best model, and none model
load("Scenarios/sim_bom.rda")
load("Scenarios/sim_bo.rda")
load("Scenarios/sim_bm.rda")
load("Scenarios/sim_om.rda")
load("Scenarios/sim_b.rda")
load("Scenarios/sim_o.rda")
load("Scenarios/sim_m.rda")
load("Scenarios/sim_none.rda")

# Load scenarios for all other models
load("Scenarios/sim_bomp.rda")
load("Scenarios/sim_bop.rda")
load("Scenarios/sim_bmp.rda")
load("Scenarios/sim_omp.rda")
load("Scenarios/sim_bp.rda")
load("Scenarios/sim_op.rda")
load("Scenarios/sim_mp.rda")
load("Scenarios/sim_p.rda")


########## IV. Results table ##########

# Create data frame
scen = data.frame(Model = rep(c("BOM","BO","BM","OM","B","O","M","None"), each = 3))
scen$O2 = round(c(sim_bom$oxy, sim_bo$oxy, sim_bm$oxy, sim_om$oxy, sim_b$oxy, sim_o$oxy, sim_m$oxy, sim_none$oxy), 1)
scen$Cod_SSB = rep(NA, nrow(scen))
scen$Cod_Yield = rep(NA, nrow(scen))
scen$Cod_Size = rep(NA, nrow(scen))
scen$Flounder_SSB = rep(NA, nrow(scen))
scen$Flounder_Yield = rep(NA, nrow(scen))
scen$Flounder_Size = rep(NA, nrow(scen))
scen$Sprat_SSB = rep(NA, nrow(scen))
scen$Sprat_Yield = rep(NA, nrow(scen))
scen$Sprat_Size = rep(NA, nrow(scen))
scen$Herring_SSB = rep(NA, nrow(scen))
scen$Herring_Yield = rep(NA, nrow(scen))
scen$Herring_Size = rep(NA, nrow(scen))

# Fill SSB
scen[1, seq(3,12,by=3)] = round(get_ssb(sim_bom$yay)/1000/1000/1000, 1)
scen[2, seq(3,12,by=3)] = round(get_ssb(sim_bom$okay)/1000/1000/1000, 1)
scen[3, seq(3,12,by=3)] = round(get_ssb(sim_bom$death)/1000/1000/1000, 1)
scen[4, seq(3,12,by=3)] = round(get_ssb(sim_bo$yay)/1000/1000/1000, 1)
scen[5, seq(3,12,by=3)] = round(get_ssb(sim_bo$okay)/1000/1000/1000, 1)
scen[6, seq(3,12,by=3)] = round(get_ssb(sim_bo$death)/1000/1000/1000, 1)
scen[7, seq(3,12,by=3)] = round(get_ssb(sim_bm$yay)/1000/1000/1000, 1)
scen[8, seq(3,12,by=3)] = round(get_ssb(sim_bm$okay)/1000/1000/1000, 1)
scen[9, seq(3,12,by=3)] = round(get_ssb(sim_bm$death)/1000/1000/1000, 1)
scen[10, seq(3,12,by=3)] = round(get_ssb(sim_om$yay)/1000/1000/1000, 1)
scen[11, seq(3,12,by=3)] = round(get_ssb(sim_om$okay)/1000/1000/1000, 1)
scen[12, seq(3,12,by=3)] = round(get_ssb(sim_om$death)/1000/1000/1000, 1)
scen[13, seq(3,12,by=3)] = round(get_ssb(sim_b$yay)/1000/1000/1000, 1)
scen[14, seq(3,12,by=3)] = round(get_ssb(sim_b$okay)/1000/1000/1000, 1)
scen[15, seq(3,12,by=3)] = round(get_ssb(sim_b$death)/1000/1000/1000, 1)
scen[16, seq(3,12,by=3)] = round(get_ssb(sim_o$yay)/1000/1000/1000, 1)
scen[17, seq(3,12,by=3)] = round(get_ssb(sim_o$okay)/1000/1000/1000, 1)
scen[18, seq(3,12,by=3)] = round(get_ssb(sim_o$death)/1000/1000/1000, 1)
scen[19, seq(3,12,by=3)] = round(get_ssb(sim_m$yay)/1000/1000/1000, 1)
scen[20, seq(3,12,by=3)] = round(get_ssb(sim_m$okay)/1000/1000/1000, 1)
scen[21, seq(3,12,by=3)] = round(get_ssb(sim_m$death)/1000/1000/1000, 1)
scen[22, seq(3,12,by=3)] = round(get_ssb(sim_none$yay)/1000/1000/1000, 1)
scen[23, seq(3,12,by=3)] = round(get_ssb(sim_none$yay)/1000/1000/1000, 1)
scen[24, seq(3,12,by=3)] = round(get_ssb(sim_none$yay)/1000/1000/1000, 1)

# Fill Yield
scen[1, seq(4,13,by=3)] = round(getYield(sim_bom$yay)[t_max+1, ]/1000/1000/1000, 1)
scen[2, seq(4,13,by=3)] = round(getYield(sim_bom$okay)[t_max+1, ]/1000/1000/1000, 1)
scen[3, seq(4,13,by=3)] = round(getYield(sim_bom$death)[t_max+1, ]/1000/1000/1000, 1)
scen[4, seq(4,13,by=3)] = round(getYield(sim_bo$yay)[t_max+1, ]/1000/1000/1000, 1)
scen[5, seq(4,13,by=3)] = round(getYield(sim_bo$okay)[t_max+1, ]/1000/1000/1000, 1)
scen[6, seq(4,13,by=3)] = round(getYield(sim_bo$death)[t_max+1, ]/1000/1000/1000, 1)
scen[7, seq(4,13,by=3)] = round(getYield(sim_bm$yay)[t_max+1, ]/1000/1000/1000, 1)
scen[8, seq(4,13,by=3)] = round(getYield(sim_bm$okay)[t_max+1, ]/1000/1000/1000, 1)
scen[9, seq(4,13,by=3)] = round(getYield(sim_bm$death)[t_max+1, ]/1000/1000/1000, 1)
scen[10, seq(4,13,by=3)] = round(getYield(sim_om$yay)[t_max+1, ]/1000/1000/1000, 1)
scen[11, seq(4,13,by=3)] = round(getYield(sim_om$okay)[t_max+1, ]/1000/1000/1000, 1)
scen[12, seq(4,13,by=3)] = round(getYield(sim_om$death)[t_max+1, ]/1000/1000/1000, 1)
scen[13, seq(4,13,by=3)] = round(getYield(sim_b$yay)[t_max+1, ]/1000/1000/1000, 1)
scen[14, seq(4,13,by=3)] = round(getYield(sim_b$okay)[t_max+1, ]/1000/1000/1000, 1)
scen[15, seq(4,13,by=3)] = round(getYield(sim_b$death)[t_max+1, ]/1000/1000/1000, 1)
scen[16, seq(4,13,by=3)] = round(getYield(sim_o$yay)[t_max+1, ]/1000/1000/1000, 1)
scen[17, seq(4,13,by=3)] = round(getYield(sim_o$okay)[t_max+1, ]/1000/1000/1000, 1)
scen[18, seq(4,13,by=3)] = round(getYield(sim_o$death)[t_max+1, ]/1000/1000/1000, 1)
scen[19, seq(4,13,by=3)] = round(getYield(sim_m$yay)[t_max+1, ]/1000/1000/1000, 1)
scen[20, seq(4,13,by=3)] = round(getYield(sim_m$okay)[t_max+1, ]/1000/1000/1000, 1)
scen[21, seq(4,13,by=3)] = round(getYield(sim_m$death)[t_max+1, ]/1000/1000/1000, 1)
scen[22, seq(4,13,by=3)] = round(getYield(sim_none$yay)[t_max+1, ]/1000/1000/1000, 1)
scen[23, seq(4,13,by=3)] = round(getYield(sim_none$okay)[t_max+1, ]/1000/1000/1000, 1)
scen[24, seq(4,13,by=3)] = round(getYield(sim_none$death)[t_max+1, ]/1000/1000/1000, 1)

# Fill size
scen[1, seq(5,14,by=3)] = round(c(myGrowthCurves(sim_bom$yay, t = t_max, species = "cod", max_age = 15)[1,50], myGrowthCurves(sim_bom$yay, t = t_max, species = "flounder", max_age = 26)[1,50], myGrowthCurves(sim_bom$yay, t = t_max, species = "sprat", max_age = 16)[1,50], myGrowthCurves(sim_bom$yay, t = t_max, species = "herring", max_age = 13)[1,50]))
scen[2, seq(5,14,by=3)] = round(c(myGrowthCurves(sim_bom$okay, t = t_max, species = "cod", max_age = 15)[1,50], myGrowthCurves(sim_bom$okay, t = t_max, species = "flounder", max_age = 26)[1,50], myGrowthCurves(sim_bom$okay, t = t_max, species = "sprat", max_age = 16)[1,50], myGrowthCurves(sim_bom$okay, t = t_max, species = "herring", max_age = 13)[1,50]))
scen[3, seq(5,14,by=3)] = round(c(myGrowthCurves(sim_bom$death, t = t_max, species = "cod", max_age = 15)[1,50], myGrowthCurves(sim_bom$death, t = t_max, species = "flounder", max_age = 26)[1,50], myGrowthCurves(sim_bom$death, t = t_max, species = "sprat", max_age = 16)[1,50], myGrowthCurves(sim_bom$death, t = t_max, species = "herring", max_age = 13)[1,50]))
scen[4, seq(5,14,by=3)] = round(c(myGrowthCurves(sim_bo$yay, t = t_max, species = "cod", max_age = 15)[1,50], myGrowthCurves(sim_bo$yay, t = t_max, species = "flounder", max_age = 26)[1,50], myGrowthCurves(sim_bo$yay, t = t_max, species = "sprat", max_age = 16)[1,50], myGrowthCurves(sim_bo$yay, t = t_max, species = "herring", max_age = 13)[1,50]))
scen[5, seq(5,14,by=3)] = round(c(myGrowthCurves(sim_bo$okay, t = t_max, species = "cod", max_age = 15)[1,50], myGrowthCurves(sim_bo$okay, t = t_max, species = "flounder", max_age = 26)[1,50], myGrowthCurves(sim_bo$okay, t = t_max, species = "sprat", max_age = 16)[1,50], myGrowthCurves(sim_bo$okay, t = t_max, species = "herring", max_age = 13)[1,50]))
scen[6, seq(5,14,by=3)] = round(c(myGrowthCurves(sim_bo$death, t = t_max, species = "cod", max_age = 15)[1,50], myGrowthCurves(sim_bo$death, t = t_max, species = "flounder", max_age = 26)[1,50], myGrowthCurves(sim_bo$death, t = t_max, species = "sprat", max_age = 16)[1,50], myGrowthCurves(sim_bo$death, t = t_max, species = "herring", max_age = 13)[1,50]))
scen[7, seq(5,14,by=3)] = round(c(myGrowthCurves(sim_bm$yay, t = t_max, species = "cod", max_age = 15)[1,50], myGrowthCurves(sim_bm$yay, t = t_max, species = "flounder", max_age = 26)[1,50], myGrowthCurves(sim_bm$yay, t = t_max, species = "sprat", max_age = 16)[1,50], myGrowthCurves(sim_bm$yay, t = t_max, species = "herring", max_age = 13)[1,50]))
scen[8, seq(5,14,by=3)] = round(c(myGrowthCurves(sim_bm$okay, t = t_max, species = "cod", max_age = 15)[1,50], myGrowthCurves(sim_bm$okay, t = t_max, species = "flounder", max_age = 26)[1,50], myGrowthCurves(sim_bm$okay, t = t_max, species = "sprat", max_age = 16)[1,50], myGrowthCurves(sim_bm$okay, t = t_max, species = "herring", max_age = 13)[1,50]))
scen[9, seq(5,14,by=3)] = round(c(myGrowthCurves(sim_bm$death, t = t_max, species = "cod", max_age = 15)[1,50], myGrowthCurves(sim_bm$death, t = t_max, species = "flounder", max_age = 26)[1,50], myGrowthCurves(sim_bm$death, t = t_max, species = "sprat", max_age = 16)[1,50], myGrowthCurves(sim_bm$death, t = t_max, species = "herring", max_age = 13)[1,50]))
scen[10, seq(5,14,by=3)] = round(c(myGrowthCurves(sim_om$yay, t = t_max, species = "cod", max_age = 15)[1,50], myGrowthCurves(sim_om$yay, t = t_max, species = "flounder", max_age = 26)[1,50], myGrowthCurves(sim_om$yay, t = t_max, species = "sprat", max_age = 16)[1,50], myGrowthCurves(sim_om$yay, t = t_max, species = "herring", max_age = 13)[1,50]))
scen[11, seq(5,14,by=3)] = round(c(myGrowthCurves(sim_om$okay, t = t_max, species = "cod", max_age = 15)[1,50], myGrowthCurves(sim_om$okay, t = t_max, species = "flounder", max_age = 26)[1,50], myGrowthCurves(sim_om$okay, t = t_max, species = "sprat", max_age = 16)[1,50], myGrowthCurves(sim_om$okay, t = t_max, species = "herring", max_age = 13)[1,50]))
scen[12, seq(5,14,by=3)] = round(c(myGrowthCurves(sim_om$death, t = t_max, species = "cod", max_age = 15)[1,50], myGrowthCurves(sim_om$death, t = t_max, species = "flounder", max_age = 26)[1,50], myGrowthCurves(sim_om$death, t = t_max, species = "sprat", max_age = 16)[1,50], myGrowthCurves(sim_om$death, t = t_max, species = "herring", max_age = 13)[1,50]))
scen[13, seq(5,14,by=3)] = round(c(myGrowthCurves(sim_b$yay, t = t_max, species = "cod", max_age = 15)[1,50], myGrowthCurves(sim_b$yay, t = t_max, species = "flounder", max_age = 26)[1,50], myGrowthCurves(sim_b$yay, t = t_max, species = "sprat", max_age = 16)[1,50], myGrowthCurves(sim_b$yay, t = t_max, species = "herring", max_age = 13)[1,50]))
scen[14, seq(5,14,by=3)] = round(c(myGrowthCurves(sim_b$okay, t = t_max, species = "cod", max_age = 15)[1,50], myGrowthCurves(sim_b$okay, t = t_max, species = "flounder", max_age = 26)[1,50], myGrowthCurves(sim_b$okay, t = t_max, species = "sprat", max_age = 16)[1,50], myGrowthCurves(sim_b$okay, t = t_max, species = "herring", max_age = 13)[1,50]))
scen[15, seq(5,14,by=3)] = round(c(myGrowthCurves(sim_b$death, t = t_max, species = "cod", max_age = 15)[1,50], myGrowthCurves(sim_b$death, t = t_max, species = "flounder", max_age = 26)[1,50], myGrowthCurves(sim_b$death, t = t_max, species = "sprat", max_age = 16)[1,50], myGrowthCurves(sim_b$death, t = t_max, species = "herring", max_age = 13)[1,50]))
scen[16, seq(5,14,by=3)] = round(c(myGrowthCurves(sim_o$yay, t = t_max, species = "cod", max_age = 15)[1,50], myGrowthCurves(sim_o$yay, t = t_max, species = "flounder", max_age = 26)[1,50], myGrowthCurves(sim_o$yay, t = t_max, species = "sprat", max_age = 16)[1,50], myGrowthCurves(sim_o$yay, t = t_max, species = "herring", max_age = 13)[1,50]))
scen[17, seq(5,14,by=3)] = round(c(myGrowthCurves(sim_o$okay, t = t_max, species = "cod", max_age = 15)[1,50], myGrowthCurves(sim_o$okay, t = t_max, species = "flounder", max_age = 26)[1,50], myGrowthCurves(sim_o$okay, t = t_max, species = "sprat", max_age = 16)[1,50], myGrowthCurves(sim_o$okay, t = t_max, species = "herring", max_age = 13)[1,50]))
scen[18, seq(5,14,by=3)] = round(c(myGrowthCurves(sim_o$death, t = t_max, species = "cod", max_age = 15)[1,50], myGrowthCurves(sim_o$death, t = t_max, species = "flounder", max_age = 26)[1,50], myGrowthCurves(sim_o$death, t = t_max, species = "sprat", max_age = 16)[1,50], myGrowthCurves(sim_o$death, t = t_max, species = "herring", max_age = 13)[1,50]))
scen[19, seq(5,14,by=3)] = round(c(myGrowthCurves(sim_m$yay, t = t_max, species = "cod", max_age = 15)[1,50], myGrowthCurves(sim_m$yay, t = t_max, species = "flounder", max_age = 26)[1,50], myGrowthCurves(sim_m$yay, t = t_max, species = "sprat", max_age = 16)[1,50], myGrowthCurves(sim_m$yay, t = t_max, species = "herring", max_age = 13)[1,50]))
scen[20, seq(5,14,by=3)] = round(c(myGrowthCurves(sim_m$okay, t = t_max, species = "cod", max_age = 15)[1,50], myGrowthCurves(sim_m$okay, t = t_max, species = "flounder", max_age = 26)[1,50], myGrowthCurves(sim_m$okay, t = t_max, species = "sprat", max_age = 16)[1,50], myGrowthCurves(sim_m$okay, t = t_max, species = "herring", max_age = 13)[1,50]))
scen[21, seq(5,14,by=3)] = round(c(myGrowthCurves(sim_m$death, t = t_max, species = "cod", max_age = 15)[1,50], myGrowthCurves(sim_m$death, t = t_max, species = "flounder", max_age = 26)[1,50], myGrowthCurves(sim_m$death, t = t_max, species = "sprat", max_age = 16)[1,50], myGrowthCurves(sim_m$death, t = t_max, species = "herring", max_age = 13)[1,50]))
scen[22, seq(5,14,by=3)] = round(c(myGrowthCurves(sim_none$yay, t = t_max, species = "cod", max_age = 15)[1,50], myGrowthCurves(sim_none$yay, t = t_max, species = "flounder", max_age = 26)[1,50], myGrowthCurves(sim_none$yay, t = t_max, species = "sprat", max_age = 16)[1,50], myGrowthCurves(sim_none$yay, t = t_max, species = "herring", max_age = 13)[1,50]))
scen[23, seq(5,14,by=3)] = round(c(myGrowthCurves(sim_none$okay, t = t_max, species = "cod", max_age = 15)[1,50], myGrowthCurves(sim_none$okay, t = t_max, species = "flounder", max_age = 26)[1,50], myGrowthCurves(sim_none$okay, t = t_max, species = "sprat", max_age = 16)[1,50], myGrowthCurves(sim_none$okay, t = t_max, species = "herring", max_age = 13)[1,50]))
scen[24, seq(5,14,by=3)] = round(c(myGrowthCurves(sim_none$death, t = t_max, species = "cod", max_age = 15)[1,50], myGrowthCurves(sim_none$death, t = t_max, species = "flounder", max_age = 26)[1,50], myGrowthCurves(sim_none$death, t = t_max, species = "sprat", max_age = 16)[1,50], myGrowthCurves(sim_none$death, t = t_max, species = "herring", max_age = 13)[1,50]))

# Remove last two rows
scen_table = scen[-c(23,24), ]

# Create table
print(xtable(scen_table), include.rownames = F)


########## IV. Plot results table ##########

color_b = "#D55E00"
color_o = "#0072B2"
color_m = "#F0E442"
color_bo = "#CC79A7"
color_bm = "#E69F00"
color_om = "#009E73"
color_bom = "#000000"
color_none = "#D4D4D4"

# Initialize plot
jpeg("Plots/table_scenarios.jpg", width = 16, height = 18, units = 'cm', res = 600)

# Layout plot
par(mfrow = c(4,3), mar = c(1.7,1.7,1.7,1.7), oma = c(4,3,0,8))

# Cod SSB
plot(seq(3), ylim = c(40,100), type = 'n', main = "SSB (kt)", cex.main = 1.5, axes = F, ylab = "Cod", cex.lab = 2)
axis(1, seq(3), rev(seq(3)), cex.axis = 1.2)
axis(2, seq(40,100,20), cex.axis = 1.2)
box()
mtext("Cod", side = 2, line = 3, cex = 1.2)
lines(seq(3), rep(scen_table[scen_table$Model == "None", ]$Cod_SSB,3), lwd = 4, lty = 1, col = color_none)
lines(seq(3), scen_table[scen_table$Model == "B", ]$Cod_SSB, lwd = 4, lty = 2, col = color_b)
lines(seq(3), scen_table[scen_table$Model == "O", ]$Cod_SSB, lwd = 4, lty = 3, col = color_o)
lines(seq(3), scen_table[scen_table$Model == "M", ]$Cod_SSB, lwd = 4, lty = 5, col = color_m)
lines(seq(3), scen_table[scen_table$Model == "BO", ]$Cod_SSB, lwd = 4, lty = 4, col = color_bo)
lines(seq(3), scen_table[scen_table$Model == "BM", ]$Cod_SSB, lwd = 4, lty = 6, col = color_bm)
lines(seq(3), scen_table[scen_table$Model == "OM", ]$Cod_SSB, lwd = 4, lty = "9212", col = color_om)
lines(seq(3), scen_table[scen_table$Model == "BOM", ]$Cod_SSB, lwd = 4, lty = 1, col = color_bom)

# Cod Yield
plot(seq(3), ylim = c(20,80), type = 'n', main = "Yield (kt)", cex.main = 1.5, axes = F)
axis(1, seq(3), rev(seq(3)), cex.axis = 1.2)
axis(2, seq(20,80,20), cex.axis = 1.2)
box()
lines(seq(3), rep(scen_table[scen_table$Model == "None", ]$Cod_Yield,3), lwd = 4, lty = 1, col = color_none)
lines(seq(3), scen_table[scen_table$Model == "B", ]$Cod_Yield, lwd = 4, lty = 2, col = color_b)
lines(seq(3), scen_table[scen_table$Model == "O", ]$Cod_Yield, lwd = 4, lty = 3, col = color_o)
lines(seq(3), scen_table[scen_table$Model == "M", ]$Cod_Yield, lwd = 4, lty = 5, col = color_m)
lines(seq(3), scen_table[scen_table$Model == "BO", ]$Cod_Yield, lwd = 4, lty = 4, col = color_bo)
lines(seq(3), scen_table[scen_table$Model == "BM", ]$Cod_Yield, lwd = 4, lty = 6, col = color_bm)
lines(seq(3), scen_table[scen_table$Model == "OM", ]$Cod_Yield, lwd = 4, lty = "9212", col = color_om)
lines(seq(3), scen_table[scen_table$Model == "BOM", ]$Cod_Yield, lwd = 4, lty = 1, col = color_bom)

# Cod Size
plot(seq(3), ylim = c(5000,15000), type = 'n', main = "Size (g)", cex.main = 1.5, axes = F)
axis(1, seq(3), rev(seq(3)), cex.axis = 1.2)
axis(2, c(5000,10000,15000), cex.axis = 1.2)
box()
lines(seq(3), rep(scen_table[scen_table$Model == "None", ]$Cod_Size,3), lwd = 4, lty = 1, col = color_none)
lines(seq(3), scen_table[scen_table$Model == "B", ]$Cod_Size, lwd = 4, lty = 2, col = color_b)
lines(seq(3), scen_table[scen_table$Model == "O", ]$Cod_Size, lwd = 4, lty = 3, col = color_o)
lines(seq(3), scen_table[scen_table$Model == "M", ]$Cod_Size, lwd = 4, lty = 5, col = color_m)
lines(seq(3), scen_table[scen_table$Model == "BO", ]$Cod_Size, lwd = 4, lty = 4, col = color_bo)
lines(seq(3), scen_table[scen_table$Model == "BM", ]$Cod_Size, lwd = 4, lty = 6, col = color_bm)
lines(seq(3), scen_table[scen_table$Model == "OM", ]$Cod_Size, lwd = 4, lty = "9212", col = color_om)
lines(seq(3), scen_table[scen_table$Model == "BOM", ]$Cod_Size, lwd = 4, lty = 1, col = color_bom)

# Flounder SSB
plot(seq(3), ylim = c(0,12), type = 'n', cex.main = 1.5, axes = F)
axis(1, seq(3), rev(seq(3)), cex.axis = 1.2)
axis(2, seq(0,12,3), cex.axis = 1.2)
box()
mtext("Flounder", side = 2, line = 3, cex = 1.2)
lines(seq(3), rep(scen_table[scen_table$Model == "None", ]$Flounder_SSB,3), lwd = 4, lty = 1, col = color_none)
lines(seq(3), scen_table[scen_table$Model == "B", ]$Flounder_SSB, lwd = 4, lty = 2, col = color_b)
lines(seq(3), scen_table[scen_table$Model == "O", ]$Flounder_SSB, lwd = 4, lty = 3, col = color_o)
lines(seq(3), scen_table[scen_table$Model == "M", ]$Flounder_SSB, lwd = 4, lty = 5, col = color_m)
lines(seq(3), scen_table[scen_table$Model == "BO", ]$Flounder_SSB, lwd = 4, lty = 4, col = color_bo)
lines(seq(3), scen_table[scen_table$Model == "BM", ]$Flounder_SSB, lwd = 4, lty = 6, col = color_bm)
lines(seq(3), scen_table[scen_table$Model == "OM", ]$Flounder_SSB, lwd = 4, lty = "9212", col = color_om)
lines(seq(3), scen_table[scen_table$Model == "BOM", ]$Flounder_SSB, lwd = 4, lty = 1, col = color_bom)

# Flounder Yield
plot(seq(3), ylim = c(0,12), type = 'n', cex.main = 1.5, axes = F)
axis(1, seq(3), rev(seq(3)), cex.axis = 1.2)
axis(2, seq(0,12,3), cex.axis = 1.2)
box()
lines(seq(3), rep(scen_table[scen_table$Model == "None", ]$Flounder_Yield,3), lwd = 4, lty = 1, col = color_none)
lines(seq(3), scen_table[scen_table$Model == "B", ]$Flounder_Yield, lwd = 4, lty = 2, col = color_b)
lines(seq(3), scen_table[scen_table$Model == "O", ]$Flounder_Yield, lwd = 4, lty = 3, col = color_o)
lines(seq(3), scen_table[scen_table$Model == "M", ]$Flounder_Yield, lwd = 4, lty = 5, col = color_m)
lines(seq(3), scen_table[scen_table$Model == "BO", ]$Flounder_Yield, lwd = 4, lty = 4, col = color_bo)
lines(seq(3), scen_table[scen_table$Model == "BM", ]$Flounder_Yield, lwd = 4, lty = 6, col = color_bm)
lines(seq(3), scen_table[scen_table$Model == "OM", ]$Flounder_Yield, lwd = 4, lty = "9212", col = color_om)
lines(seq(3), scen_table[scen_table$Model == "BOM", ]$Flounder_Yield, lwd = 4, lty = 1, col = color_bom)

# Flounder Size
plot(seq(3), ylim = c(800,1200), type = 'n', cex.main = 1.5, axes = F)
axis(1, seq(3), rev(seq(3)), cex.axis = 1.2)
axis(2, seq(800,1200,200), cex.axis = 1.2)
box()
lines(seq(3), rep(scen_table[scen_table$Model == "None", ]$Flounder_Size,3), lwd = 4, lty = 1, col = color_none)
lines(seq(3), scen_table[scen_table$Model == "B", ]$Flounder_Size, lwd = 4, lty = 2, col = color_b)
lines(seq(3), scen_table[scen_table$Model == "O", ]$Flounder_Size, lwd = 4, lty = 3, col = color_o)
lines(seq(3), scen_table[scen_table$Model == "M", ]$Flounder_Size, lwd = 4, lty = 5, col = color_m)
lines(seq(3), scen_table[scen_table$Model == "BO", ]$Flounder_Size, lwd = 4, lty = 4, col = color_bo)
lines(seq(3), scen_table[scen_table$Model == "BM", ]$Flounder_Size, lwd = 4, lty = 6, col = color_bm)
lines(seq(3), scen_table[scen_table$Model == "OM", ]$Flounder_Size, lwd = 4, lty = "9212", col = color_om)
lines(seq(3), scen_table[scen_table$Model == "BOM", ]$Flounder_Size, lwd = 4, lty = 1, col = color_bom)

# Sprat SSB
plot(seq(3), ylim = c(800,2000), type = 'n', cex.main = 1.5, axes = F)
axis(1, seq(3), rev(seq(3)), cex.axis = 1.2)
axis(2, seq(800,2000,600), cex.axis = 1.2)
box()
mtext("Sprat", side = 2, line = 3, cex = 1.2)
lines(seq(3), rep(scen_table[scen_table$Model == "None", ]$Sprat_SSB,3), lwd = 4, lty = 1, col = color_none)
lines(seq(3), scen_table[scen_table$Model == "B", ]$Sprat_SSB, lwd = 4, lty = 2, col = color_b)
lines(seq(3), scen_table[scen_table$Model == "O", ]$Sprat_SSB, lwd = 4, lty = 3, col = color_o)
lines(seq(3), scen_table[scen_table$Model == "M", ]$Sprat_SSB, lwd = 4, lty = 5, col = color_m)
lines(seq(3), scen_table[scen_table$Model == "BO", ]$Sprat_SSB, lwd = 4, lty = 4, col = color_bo)
lines(seq(3), scen_table[scen_table$Model == "BM", ]$Sprat_SSB, lwd = 4, lty = 6, col = color_bm)
lines(seq(3), scen_table[scen_table$Model == "OM", ]$Sprat_SSB, lwd = 4, lty = "9212", col = color_om)
lines(seq(3), scen_table[scen_table$Model == "BOM", ]$Sprat_SSB, lwd = 4, lty = 1, col = color_bom)

# Sprat Yield
plot(seq(3), ylim = c(200,600), type = 'n', cex.main = 1.5, axes = F)
axis(1, seq(3), rev(seq(3)), cex.axis = 1.2)
axis(2, seq(200,600,200), cex.axis = 1.2)
box()
lines(seq(3), rep(scen_table[scen_table$Model == "None", ]$Sprat_Yield,3), lwd = 4, lty = 1, col = color_none)
lines(seq(3), scen_table[scen_table$Model == "B", ]$Sprat_Yield, lwd = 4, lty = 2, col = color_b)
lines(seq(3), scen_table[scen_table$Model == "O", ]$Sprat_Yield, lwd = 4, lty = 3, col = color_o)
lines(seq(3), scen_table[scen_table$Model == "M", ]$Sprat_Yield, lwd = 4, lty = 5, col = color_m)
lines(seq(3), scen_table[scen_table$Model == "BO", ]$Sprat_Yield, lwd = 4, lty = 4, col = color_bo)
lines(seq(3), scen_table[scen_table$Model == "BM", ]$Sprat_Yield, lwd = 4, lty = 6, col = color_bm)
lines(seq(3), scen_table[scen_table$Model == "OM", ]$Sprat_Yield, lwd = 4, lty = "9212", col = color_om)
lines(seq(3), scen_table[scen_table$Model == "BOM", ]$Sprat_Yield, lwd = 4, lty = 1, col = color_bom)

# Sprat Size
plot(seq(3), ylim = c(15,25), type = 'n', cex.main = 1.5, axes = F)
axis(1, seq(3), rev(seq(3)), cex.axis = 1.2)
axis(2, seq(15,25,5), cex.axis = 1.2)
box()
lines(seq(3), rep(scen_table[scen_table$Model == "None", ]$Sprat_Size,3), lwd = 4, lty = 1, col = color_none)
lines(seq(3), scen_table[scen_table$Model == "B", ]$Sprat_Size, lwd = 4, lty = 2, col = color_b)
lines(seq(3), scen_table[scen_table$Model == "O", ]$Sprat_Size, lwd = 4, lty = 3, col = color_o)
lines(seq(3), scen_table[scen_table$Model == "M", ]$Sprat_Size, lwd = 4, lty = 5, col = color_m)
lines(seq(3), scen_table[scen_table$Model == "BO", ]$Sprat_Size, lwd = 4, lty = 4, col = color_bo)
lines(seq(3), scen_table[scen_table$Model == "BM", ]$Sprat_Size, lwd = 4, lty = 6, col = color_bm)
lines(seq(3), scen_table[scen_table$Model == "OM", ]$Sprat_Size, lwd = 4, lty = "9212", col = color_om)
lines(seq(3), scen_table[scen_table$Model == "BOM", ]$Sprat_Size, lwd = 4, lty = 1, col = color_bom)

# Herring SSB
plot(seq(3), ylim = c(400,600), type = 'n', cex.main = 1.5, axes = F)
axis(1, seq(3), rev(seq(3)), cex.axis = 1.2)
axis(2, seq(400,600,100), cex.axis = 1.2)
box()
mtext("Herring", side = 2, line = 3, cex = 1.2)
lines(seq(3), rep(scen_table[scen_table$Model == "None", ]$Herring_SSB,3), lwd = 4, lty = 1, col = color_none)
lines(seq(3), scen_table[scen_table$Model == "B", ]$Herring_SSB, lwd = 4, lty = 2, col = color_b)
lines(seq(3), scen_table[scen_table$Model == "O", ]$Herring_SSB, lwd = 4, lty = 3, col = color_o)
lines(seq(3), scen_table[scen_table$Model == "M", ]$Herring_SSB, lwd = 4, lty = 5, col = color_m)
lines(seq(3), scen_table[scen_table$Model == "BO", ]$Herring_SSB, lwd = 4, lty = 4, col = color_bo)
lines(seq(3), scen_table[scen_table$Model == "BM", ]$Herring_SSB, lwd = 4, lty = 6, col = color_bm)
lines(seq(3), scen_table[scen_table$Model == "OM", ]$Herring_SSB, lwd = 4, lty = "9212", col = color_om)
lines(seq(3), scen_table[scen_table$Model == "BOM", ]$Herring_SSB, lwd = 4, lty = 1, col = color_bom)

# Herring Yield
plot(seq(3), ylim = c(150,300), type = 'n', cex.main = 1.5, axes = F)
axis(1, seq(3), rev(seq(3)), cex.axis = 1.2)
axis(2, seq(150,300,50), cex.axis = 1.2)
box()
lines(seq(3), rep(scen_table[scen_table$Model == "None", ]$Herring_Yield,3), lwd = 4, lty = 1, col = color_none)
lines(seq(3), scen_table[scen_table$Model == "B", ]$Herring_Yield, lwd = 4, lty = 2, col = color_b)
lines(seq(3), scen_table[scen_table$Model == "O", ]$Herring_Yield, lwd = 4, lty = 3, col = color_o)
lines(seq(3), scen_table[scen_table$Model == "M", ]$Herring_Yield, lwd = 4, lty = 5, col = color_m)
lines(seq(3), scen_table[scen_table$Model == "BO", ]$Herring_Yield, lwd = 4, lty = 4, col = color_bo)
lines(seq(3), scen_table[scen_table$Model == "BM", ]$Herring_Yield, lwd = 4, lty = 6, col = color_bm)
lines(seq(3), scen_table[scen_table$Model == "OM", ]$Herring_Yield, lwd = 4, lty = "9212", col = color_om)
lines(seq(3), scen_table[scen_table$Model == "BOM", ]$Herring_Yield, lwd = 4, lty = 1, col = color_bom)

# Herring Size
plot(seq(3), ylim = c(170,180), type = 'n', cex.main = 1.5, axes = F)
axis(1, seq(3), rev(seq(3)), cex.axis = 1.2)
axis(2, seq(170,180,5), cex.axis = 1.2)
box()
lines(seq(3), rep(scen_table[scen_table$Model == "None", ]$Herring_Size,3), lwd = 4, lty = 1, col = color_none)
lines(seq(3), scen_table[scen_table$Model == "B", ]$Herring_Size, lwd = 4, lty = 2, col = color_b)
lines(seq(3), scen_table[scen_table$Model == "O", ]$Herring_Size, lwd = 4, lty = 3, col = color_o)
lines(seq(3), scen_table[scen_table$Model == "M", ]$Herring_Size, lwd = 4, lty = 5, col = color_m)
lines(seq(3), scen_table[scen_table$Model == "BO", ]$Herring_Size, lwd = 4, lty = 4, col = color_bo)
lines(seq(3), scen_table[scen_table$Model == "BM", ]$Herring_Size, lwd = 4, lty = 6, col = color_bm)
lines(seq(3), scen_table[scen_table$Model == "OM", ]$Herring_Size, lwd = 4, lty = "9212", col = color_om)
lines(seq(3), scen_table[scen_table$Model == "BOM", ]$Herring_Size, lwd = 4, lty = 1, col = color_bom)

# Axis label
mtext("Oxygen (mL/L)", side = 1, line = 2, outer = T, cex = 1.8)

# Add legend
par(fig = c(0,1,0,1), oma = c(0,0,0,0), mar = c(0,0,0,0), new = T)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend("right", bty = 'n', c("B","O","M","BO","BM","OM","BOM","None"), col = c(color_b,color_o,color_m,color_bo,color_bm,color_om,color_bom,color_none), lwd = 4, lty = c("dashed","dotted","longdash","dotdash","twodash","9212","solid","solid"), xpd = TRUE, cex = 1.2, seg.len=4)

# Finish plot
dev.off()


########## V. Results table for supplement ##########

# Create data frame
supp = data.frame(Model = rep(c("BOMP","BOM","BOP","BMP","OMP","BO","BM","BP","OM","OP","MP","B","O","M","P","None"), each = 3))
supp$O2 = round(c(sim_bomp$oxy, sim_bom$oxy, sim_bop$oxy, sim_bmp$oxy, sim_omp$oxy, sim_bo$oxy, sim_bm$oxy, sim_bp$oxy, sim_om$oxy, sim_op$oxy, sim_mp$oxy, sim_b$oxy, sim_o$oxy, sim_m$oxy, sim_p$oxy, sim_none$oxy), 1)
supp$Cod_SSB = rep(NA, nrow(supp))
supp$Cod_Yield = rep(NA, nrow(supp))
supp$Cod_Size = rep(NA, nrow(supp))
supp$Flounder_SSB = rep(NA, nrow(supp))
supp$Flounder_Yield = rep(NA, nrow(supp))
supp$Flounder_Size = rep(NA, nrow(supp))
supp$Sprat_SSB = rep(NA, nrow(supp))
supp$Sprat_Yield = rep(NA, nrow(supp))
supp$Sprat_Size = rep(NA, nrow(supp))
supp$Herring_SSB = rep(NA, nrow(supp))
supp$Herring_Yield = rep(NA, nrow(supp))
supp$Herring_Size = rep(NA, nrow(supp))

# Fill SSB
supp[1, seq(3,12,by=3)] = round(get_ssb(sim_bomp$yay)/1000/1000/1000, 1)
supp[2, seq(3,12,by=3)] = round(get_ssb(sim_bomp$okay)/1000/1000/1000, 1)
supp[3, seq(3,12,by=3)] = round(get_ssb(sim_bomp$death)/1000/1000/1000, 1)
supp[4, seq(3,12,by=3)] = round(get_ssb(sim_bom$yay)/1000/1000/1000, 1)
supp[5, seq(3,12,by=3)] = round(get_ssb(sim_bom$okay)/1000/1000/1000, 1)
supp[6, seq(3,12,by=3)] = round(get_ssb(sim_bom$death)/1000/1000/1000, 1)
supp[7, seq(3,12,by=3)] = round(get_ssb(sim_bop$yay)/1000/1000/1000, 1)
supp[8, seq(3,12,by=3)] = round(get_ssb(sim_bop$okay)/1000/1000/1000, 1)
supp[9, seq(3,12,by=3)] = round(get_ssb(sim_bop$death)/1000/1000/1000, 1)
supp[10, seq(3,12,by=3)] = round(get_ssb(sim_bmp$yay)/1000/1000/1000, 1)
supp[11, seq(3,12,by=3)] = round(get_ssb(sim_bmp$okay)/1000/1000/1000, 1)
supp[12, seq(3,12,by=3)] = round(get_ssb(sim_bmp$death)/1000/1000/1000, 1)
supp[13, seq(3,12,by=3)] = round(get_ssb(sim_omp$yay)/1000/1000/1000, 1)
supp[14, seq(3,12,by=3)] = round(get_ssb(sim_omp$okay)/1000/1000/1000, 1)
supp[15, seq(3,12,by=3)] = round(get_ssb(sim_omp$death)/1000/1000/1000, 1)
supp[16, seq(3,12,by=3)] = round(get_ssb(sim_bo$yay)/1000/1000/1000, 1)
supp[17, seq(3,12,by=3)] = round(get_ssb(sim_bo$okay)/1000/1000/1000, 1)
supp[18, seq(3,12,by=3)] = round(get_ssb(sim_bo$death)/1000/1000/1000, 1)
supp[19, seq(3,12,by=3)] = round(get_ssb(sim_bm$yay)/1000/1000/1000, 1)
supp[20, seq(3,12,by=3)] = round(get_ssb(sim_bm$okay)/1000/1000/1000, 1)
supp[21, seq(3,12,by=3)] = round(get_ssb(sim_bm$death)/1000/1000/1000, 1)
supp[22, seq(3,12,by=3)] = round(get_ssb(sim_bp$yay)/1000/1000/1000, 1)
supp[23, seq(3,12,by=3)] = round(get_ssb(sim_bp$okay)/1000/1000/1000, 1)
supp[24, seq(3,12,by=3)] = round(get_ssb(sim_bp$death)/1000/1000/1000, 1)
supp[25, seq(3,12,by=3)] = round(get_ssb(sim_om$yay)/1000/1000/1000, 1)
supp[26, seq(3,12,by=3)] = round(get_ssb(sim_om$okay)/1000/1000/1000, 1)
supp[27, seq(3,12,by=3)] = round(get_ssb(sim_om$death)/1000/1000/1000, 1)
supp[28, seq(3,12,by=3)] = round(get_ssb(sim_op$yay)/1000/1000/1000, 1)
supp[29, seq(3,12,by=3)] = round(get_ssb(sim_op$okay)/1000/1000/1000, 1)
supp[30, seq(3,12,by=3)] = round(get_ssb(sim_op$death)/1000/1000/1000, 1)
supp[31, seq(3,12,by=3)] = round(get_ssb(sim_mp$yay)/1000/1000/1000, 1)
supp[32, seq(3,12,by=3)] = round(get_ssb(sim_mp$okay)/1000/1000/1000, 1)
supp[33, seq(3,12,by=3)] = round(get_ssb(sim_mp$death)/1000/1000/1000, 1)
supp[34, seq(3,12,by=3)] = round(get_ssb(sim_b$yay)/1000/1000/1000, 1)
supp[35, seq(3,12,by=3)] = round(get_ssb(sim_b$okay)/1000/1000/1000, 1)
supp[36, seq(3,12,by=3)] = round(get_ssb(sim_b$death)/1000/1000/1000, 1)
supp[37, seq(3,12,by=3)] = round(get_ssb(sim_o$yay)/1000/1000/1000, 1)
supp[38, seq(3,12,by=3)] = round(get_ssb(sim_o$okay)/1000/1000/1000, 1)
supp[39, seq(3,12,by=3)] = round(get_ssb(sim_o$death)/1000/1000/1000, 1)
supp[40, seq(3,12,by=3)] = round(get_ssb(sim_m$yay)/1000/1000/1000, 1)
supp[41, seq(3,12,by=3)] = round(get_ssb(sim_m$okay)/1000/1000/1000, 1)
supp[42, seq(3,12,by=3)] = round(get_ssb(sim_m$death)/1000/1000/1000, 1)
supp[43, seq(3,12,by=3)] = round(get_ssb(sim_p$yay)/1000/1000/1000, 1)
supp[44, seq(3,12,by=3)] = round(get_ssb(sim_p$okay)/1000/1000/1000, 1)
supp[45, seq(3,12,by=3)] = round(get_ssb(sim_p$death)/1000/1000/1000, 1)
supp[46, seq(3,12,by=3)] = round(get_ssb(sim_none$yay)/1000/1000/1000, 1)
supp[47, seq(3,12,by=3)] = round(get_ssb(sim_none$okay)/1000/1000/1000, 1)
supp[48, seq(3,12,by=3)] = round(get_ssb(sim_none$death)/1000/1000/1000, 1)

# Fill Yield
supp[1, seq(4,13,by=3)] = round(getYield(sim_bomp$yay)[t_max+1, ]/1000/1000/1000, 1)
supp[2, seq(4,13,by=3)] = round(getYield(sim_bomp$okay)[t_max+1, ]/1000/1000/1000, 1)
supp[3, seq(4,13,by=3)] = round(getYield(sim_bomp$death)[t_max+1, ]/1000/1000/1000, 1)
supp[4, seq(4,13,by=3)] = round(getYield(sim_bom$yay)[t_max+1, ]/1000/1000/1000, 1)
supp[5, seq(4,13,by=3)] = round(getYield(sim_bom$okay)[t_max+1, ]/1000/1000/1000, 1)
supp[6, seq(4,13,by=3)] = round(getYield(sim_bom$death)[t_max+1, ]/1000/1000/1000, 1)
supp[7, seq(4,13,by=3)] = round(getYield(sim_bop$yay)[t_max+1, ]/1000/1000/1000, 1)
supp[8, seq(4,13,by=3)] = round(getYield(sim_bop$okay)[t_max+1, ]/1000/1000/1000, 1)
supp[9, seq(4,13,by=3)] = round(getYield(sim_bop$death)[t_max+1, ]/1000/1000/1000, 1)
supp[10, seq(4,13,by=3)] = round(getYield(sim_bmp$yay)[t_max+1, ]/1000/1000/1000, 1)
supp[11, seq(4,13,by=3)] = round(getYield(sim_bmp$okay)[t_max+1, ]/1000/1000/1000, 1)
supp[12, seq(4,13,by=3)] = round(getYield(sim_bmp$death)[t_max+1, ]/1000/1000/1000, 1)
supp[13, seq(4,13,by=3)] = round(getYield(sim_omp$yay)[t_max+1, ]/1000/1000/1000, 1)
supp[14, seq(4,13,by=3)] = round(getYield(sim_omp$okay)[t_max+1, ]/1000/1000/1000, 1)
supp[15, seq(4,13,by=3)] = round(getYield(sim_omp$death)[t_max+1, ]/1000/1000/1000, 1)
supp[16, seq(4,13,by=3)] = round(getYield(sim_bo$yay)[t_max+1, ]/1000/1000/1000, 1)
supp[17, seq(4,13,by=3)] = round(getYield(sim_bo$okay)[t_max+1, ]/1000/1000/1000, 1)
supp[18, seq(4,13,by=3)] = round(getYield(sim_bo$death)[t_max+1, ]/1000/1000/1000, 1)
supp[19, seq(4,13,by=3)] = round(getYield(sim_bm$yay)[t_max+1, ]/1000/1000/1000, 1)
supp[20, seq(4,13,by=3)] = round(getYield(sim_bm$okay)[t_max+1, ]/1000/1000/1000, 1)
supp[21, seq(4,13,by=3)] = round(getYield(sim_bm$death)[t_max+1, ]/1000/1000/1000, 1)
supp[22, seq(4,13,by=3)] = round(getYield(sim_bp$yay)[t_max+1, ]/1000/1000/1000, 1)
supp[23, seq(4,13,by=3)] = round(getYield(sim_bp$okay)[t_max+1, ]/1000/1000/1000, 1)
supp[24, seq(4,13,by=3)] = round(getYield(sim_bp$death)[t_max+1, ]/1000/1000/1000, 1)
supp[25, seq(4,13,by=3)] = round(getYield(sim_om$yay)[t_max+1, ]/1000/1000/1000, 1)
supp[26, seq(4,13,by=3)] = round(getYield(sim_om$okay)[t_max+1, ]/1000/1000/1000, 1)
supp[27, seq(4,13,by=3)] = round(getYield(sim_om$death)[t_max+1, ]/1000/1000/1000, 1)
supp[28, seq(4,13,by=3)] = round(getYield(sim_op$yay)[t_max+1, ]/1000/1000/1000, 1)
supp[29, seq(4,13,by=3)] = round(getYield(sim_op$okay)[t_max+1, ]/1000/1000/1000, 1)
supp[30, seq(4,13,by=3)] = round(getYield(sim_op$death)[t_max+1, ]/1000/1000/1000, 1)
supp[31, seq(4,13,by=3)] = round(getYield(sim_mp$yay)[t_max+1, ]/1000/1000/1000, 1)
supp[32, seq(4,13,by=3)] = round(getYield(sim_mp$okay)[t_max+1, ]/1000/1000/1000, 1)
supp[33, seq(4,13,by=3)] = round(getYield(sim_mp$death)[t_max+1, ]/1000/1000/1000, 1)
supp[34, seq(4,13,by=3)] = round(getYield(sim_b$yay)[t_max+1, ]/1000/1000/1000, 1)
supp[35, seq(4,13,by=3)] = round(getYield(sim_b$okay)[t_max+1, ]/1000/1000/1000, 1)
supp[36, seq(4,13,by=3)] = round(getYield(sim_b$death)[t_max+1, ]/1000/1000/1000, 1)
supp[37, seq(4,13,by=3)] = round(getYield(sim_o$yay)[t_max+1, ]/1000/1000/1000, 1)
supp[38, seq(4,13,by=3)] = round(getYield(sim_o$okay)[t_max+1, ]/1000/1000/1000, 1)
supp[39, seq(4,13,by=3)] = round(getYield(sim_o$death)[t_max+1, ]/1000/1000/1000, 1)
supp[40, seq(4,13,by=3)] = round(getYield(sim_m$yay)[t_max+1, ]/1000/1000/1000, 1)
supp[41, seq(4,13,by=3)] = round(getYield(sim_m$okay)[t_max+1, ]/1000/1000/1000, 1)
supp[42, seq(4,13,by=3)] = round(getYield(sim_m$death)[t_max+1, ]/1000/1000/1000, 1)
supp[43, seq(4,13,by=3)] = round(getYield(sim_p$yay)[t_max+1, ]/1000/1000/1000, 1)
supp[44, seq(4,13,by=3)] = round(getYield(sim_p$okay)[t_max+1, ]/1000/1000/1000, 1)
supp[45, seq(4,13,by=3)] = round(getYield(sim_p$death)[t_max+1, ]/1000/1000/1000, 1)
supp[46, seq(4,13,by=3)] = round(getYield(sim_none$yay)[t_max+1, ]/1000/1000/1000, 1)
supp[47, seq(4,13,by=3)] = round(getYield(sim_none$okay)[t_max+1, ]/1000/1000/1000, 1)
supp[48, seq(4,13,by=3)] = round(getYield(sim_none$death)[t_max+1, ]/1000/1000/1000, 1)

# Fill size
supp[1, seq(5,14,by=3)] = round(c(myGrowthCurves(sim_bomp$yay, t = t_max, species = "cod", max_age = 15)[1,50], myGrowthCurves(sim_bomp$yay, t = t_max, species = "flounder", max_age = 26)[1,50], myGrowthCurves(sim_bomp$yay, t = t_max, species = "sprat", max_age = 16)[1,50], myGrowthCurves(sim_bomp$yay, t = t_max, species = "herring", max_age = 13)[1,50]))
supp[2, seq(5,14,by=3)] = round(c(myGrowthCurves(sim_bomp$okay, t = t_max, species = "cod", max_age = 15)[1,50], myGrowthCurves(sim_bomp$okay, t = t_max, species = "flounder", max_age = 26)[1,50], myGrowthCurves(sim_bomp$okay, t = t_max, species = "sprat", max_age = 16)[1,50], myGrowthCurves(sim_bomp$okay, t = t_max, species = "herring", max_age = 13)[1,50]))
supp[3, seq(5,14,by=3)] = round(c(myGrowthCurves(sim_bomp$death, t = t_max, species = "cod", max_age = 15)[1,50], myGrowthCurves(sim_bomp$death, t = t_max, species = "flounder", max_age = 26)[1,50], myGrowthCurves(sim_bomp$death, t = t_max, species = "sprat", max_age = 16)[1,50], myGrowthCurves(sim_bomp$death, t = t_max, species = "herring", max_age = 13)[1,50]))
supp[4, seq(5,14,by=3)] = round(c(myGrowthCurves(sim_bom$yay, t = t_max, species = "cod", max_age = 15)[1,50], myGrowthCurves(sim_bom$yay, t = t_max, species = "flounder", max_age = 26)[1,50], myGrowthCurves(sim_bom$yay, t = t_max, species = "sprat", max_age = 16)[1,50], myGrowthCurves(sim_bom$yay, t = t_max, species = "herring", max_age = 13)[1,50]))
supp[5, seq(5,14,by=3)] = round(c(myGrowthCurves(sim_bom$okay, t = t_max, species = "cod", max_age = 15)[1,50], myGrowthCurves(sim_bom$okay, t = t_max, species = "flounder", max_age = 26)[1,50], myGrowthCurves(sim_bom$okay, t = t_max, species = "sprat", max_age = 16)[1,50], myGrowthCurves(sim_bom$okay, t = t_max, species = "herring", max_age = 13)[1,50]))
supp[6, seq(5,14,by=3)] = round(c(myGrowthCurves(sim_bom$death, t = t_max, species = "cod", max_age = 15)[1,50], myGrowthCurves(sim_bom$death, t = t_max, species = "flounder", max_age = 26)[1,50], myGrowthCurves(sim_bom$death, t = t_max, species = "sprat", max_age = 16)[1,50], myGrowthCurves(sim_bom$death, t = t_max, species = "herring", max_age = 13)[1,50]))
supp[7, seq(5,14,by=3)] = round(c(myGrowthCurves(sim_bop$yay, t = t_max, species = "cod", max_age = 15)[1,50], myGrowthCurves(sim_bop$yay, t = t_max, species = "flounder", max_age = 26)[1,50], myGrowthCurves(sim_bop$yay, t = t_max, species = "sprat", max_age = 16)[1,50], myGrowthCurves(sim_bop$yay, t = t_max, species = "herring", max_age = 13)[1,50]))
supp[8, seq(5,14,by=3)] = round(c(myGrowthCurves(sim_bop$okay, t = t_max, species = "cod", max_age = 15)[1,50], myGrowthCurves(sim_bop$okay, t = t_max, species = "flounder", max_age = 26)[1,50], myGrowthCurves(sim_bop$okay, t = t_max, species = "sprat", max_age = 16)[1,50], myGrowthCurves(sim_bop$okay, t = t_max, species = "herring", max_age = 13)[1,50]))
supp[9, seq(5,14,by=3)] = round(c(myGrowthCurves(sim_bop$death, t = t_max, species = "cod", max_age = 15)[1,50], myGrowthCurves(sim_bop$death, t = t_max, species = "flounder", max_age = 26)[1,50], myGrowthCurves(sim_bop$death, t = t_max, species = "sprat", max_age = 16)[1,50], myGrowthCurves(sim_bop$death, t = t_max, species = "herring", max_age = 13)[1,50]))
supp[10, seq(5,14,by=3)] = round(c(myGrowthCurves(sim_bmp$yay, t = t_max, species = "cod", max_age = 15)[1,50], myGrowthCurves(sim_bmp$yay, t = t_max, species = "flounder", max_age = 26)[1,50], myGrowthCurves(sim_bmp$yay, t = t_max, species = "sprat", max_age = 16)[1,50], myGrowthCurves(sim_bmp$yay, t = t_max, species = "herring", max_age = 13)[1,50]))
supp[11, seq(5,14,by=3)] = round(c(myGrowthCurves(sim_bmp$okay, t = t_max, species = "cod", max_age = 15)[1,50], myGrowthCurves(sim_bmp$okay, t = t_max, species = "flounder", max_age = 26)[1,50], myGrowthCurves(sim_bmp$okay, t = t_max, species = "sprat", max_age = 16)[1,50], myGrowthCurves(sim_bmp$okay, t = t_max, species = "herring", max_age = 13)[1,50]))
supp[12, seq(5,14,by=3)] = round(c(myGrowthCurves(sim_bmp$death, t = t_max, species = "cod", max_age = 15)[1,50], myGrowthCurves(sim_bmp$death, t = t_max, species = "flounder", max_age = 26)[1,50], myGrowthCurves(sim_bmp$death, t = t_max, species = "sprat", max_age = 16)[1,50], myGrowthCurves(sim_bmp$death, t = t_max, species = "herring", max_age = 13)[1,50]))
supp[13, seq(5,14,by=3)] = round(c(myGrowthCurves(sim_omp$yay, t = t_max, species = "cod", max_age = 15)[1,50], myGrowthCurves(sim_omp$yay, t = t_max, species = "flounder", max_age = 26)[1,50], myGrowthCurves(sim_omp$yay, t = t_max, species = "sprat", max_age = 16)[1,50], myGrowthCurves(sim_omp$yay, t = t_max, species = "herring", max_age = 13)[1,50]))
supp[14, seq(5,14,by=3)] = round(c(myGrowthCurves(sim_omp$okay, t = t_max, species = "cod", max_age = 15)[1,50], myGrowthCurves(sim_omp$okay, t = t_max, species = "flounder", max_age = 26)[1,50], myGrowthCurves(sim_omp$okay, t = t_max, species = "sprat", max_age = 16)[1,50], myGrowthCurves(sim_omp$okay, t = t_max, species = "herring", max_age = 13)[1,50]))
supp[15, seq(5,14,by=3)] = round(c(myGrowthCurves(sim_omp$death, t = t_max, species = "cod", max_age = 15)[1,50], myGrowthCurves(sim_omp$death, t = t_max, species = "flounder", max_age = 26)[1,50], myGrowthCurves(sim_omp$death, t = t_max, species = "sprat", max_age = 16)[1,50], myGrowthCurves(sim_omp$death, t = t_max, species = "herring", max_age = 13)[1,50]))
supp[16, seq(5,14,by=3)] = round(c(myGrowthCurves(sim_bo$yay, t = t_max, species = "cod", max_age = 15)[1,50], myGrowthCurves(sim_bo$yay, t = t_max, species = "flounder", max_age = 26)[1,50], myGrowthCurves(sim_bo$yay, t = t_max, species = "sprat", max_age = 16)[1,50], myGrowthCurves(sim_bo$yay, t = t_max, species = "herring", max_age = 13)[1,50]))
supp[17, seq(5,14,by=3)] = round(c(myGrowthCurves(sim_bo$okay, t = t_max, species = "cod", max_age = 15)[1,50], myGrowthCurves(sim_bo$okay, t = t_max, species = "flounder", max_age = 26)[1,50], myGrowthCurves(sim_bo$okay, t = t_max, species = "sprat", max_age = 16)[1,50], myGrowthCurves(sim_bo$okay, t = t_max, species = "herring", max_age = 13)[1,50]))
supp[18, seq(5,14,by=3)] = round(c(myGrowthCurves(sim_bo$death, t = t_max, species = "cod", max_age = 15)[1,50], myGrowthCurves(sim_bo$death, t = t_max, species = "flounder", max_age = 26)[1,50], myGrowthCurves(sim_bo$death, t = t_max, species = "sprat", max_age = 16)[1,50], myGrowthCurves(sim_bo$death, t = t_max, species = "herring", max_age = 13)[1,50]))
supp[19, seq(5,14,by=3)] = round(c(myGrowthCurves(sim_bm$yay, t = t_max, species = "cod", max_age = 15)[1,50], myGrowthCurves(sim_bm$yay, t = t_max, species = "flounder", max_age = 26)[1,50], myGrowthCurves(sim_bm$yay, t = t_max, species = "sprat", max_age = 16)[1,50], myGrowthCurves(sim_bm$yay, t = t_max, species = "herring", max_age = 13)[1,50]))
supp[20, seq(5,14,by=3)] = round(c(myGrowthCurves(sim_bm$okay, t = t_max, species = "cod", max_age = 15)[1,50], myGrowthCurves(sim_bm$okay, t = t_max, species = "flounder", max_age = 26)[1,50], myGrowthCurves(sim_bm$okay, t = t_max, species = "sprat", max_age = 16)[1,50], myGrowthCurves(sim_bm$okay, t = t_max, species = "herring", max_age = 13)[1,50]))
supp[21, seq(5,14,by=3)] = round(c(myGrowthCurves(sim_bm$death, t = t_max, species = "cod", max_age = 15)[1,50], myGrowthCurves(sim_bm$death, t = t_max, species = "flounder", max_age = 26)[1,50], myGrowthCurves(sim_bm$death, t = t_max, species = "sprat", max_age = 16)[1,50], myGrowthCurves(sim_bm$death, t = t_max, species = "herring", max_age = 13)[1,50]))
supp[22, seq(5,14,by=3)] = round(c(myGrowthCurves(sim_bp$yay, t = t_max, species = "cod", max_age = 15)[1,50], myGrowthCurves(sim_bp$yay, t = t_max, species = "flounder", max_age = 26)[1,50], myGrowthCurves(sim_bp$yay, t = t_max, species = "sprat", max_age = 16)[1,50], myGrowthCurves(sim_bp$yay, t = t_max, species = "herring", max_age = 13)[1,50]))
supp[23, seq(5,14,by=3)] = round(c(myGrowthCurves(sim_bp$okay, t = t_max, species = "cod", max_age = 15)[1,50], myGrowthCurves(sim_bp$okay, t = t_max, species = "flounder", max_age = 26)[1,50], myGrowthCurves(sim_bp$okay, t = t_max, species = "sprat", max_age = 16)[1,50], myGrowthCurves(sim_bp$okay, t = t_max, species = "herring", max_age = 13)[1,50]))
supp[24, seq(5,14,by=3)] = round(c(myGrowthCurves(sim_bp$death, t = t_max, species = "cod", max_age = 15)[1,50], myGrowthCurves(sim_bp$death, t = t_max, species = "flounder", max_age = 26)[1,50], myGrowthCurves(sim_bp$death, t = t_max, species = "sprat", max_age = 16)[1,50], myGrowthCurves(sim_bp$death, t = t_max, species = "herring", max_age = 13)[1,50]))
supp[25, seq(5,14,by=3)] = round(c(myGrowthCurves(sim_om$yay, t = t_max, species = "cod", max_age = 15)[1,50], myGrowthCurves(sim_om$yay, t = t_max, species = "flounder", max_age = 26)[1,50], myGrowthCurves(sim_om$yay, t = t_max, species = "sprat", max_age = 16)[1,50], myGrowthCurves(sim_om$yay, t = t_max, species = "herring", max_age = 13)[1,50]))
supp[26, seq(5,14,by=3)] = round(c(myGrowthCurves(sim_om$okay, t = t_max, species = "cod", max_age = 15)[1,50], myGrowthCurves(sim_om$okay, t = t_max, species = "flounder", max_age = 26)[1,50], myGrowthCurves(sim_om$okay, t = t_max, species = "sprat", max_age = 16)[1,50], myGrowthCurves(sim_om$okay, t = t_max, species = "herring", max_age = 13)[1,50]))
supp[27, seq(5,14,by=3)] = round(c(myGrowthCurves(sim_om$death, t = t_max, species = "cod", max_age = 15)[1,50], myGrowthCurves(sim_om$death, t = t_max, species = "flounder", max_age = 26)[1,50], myGrowthCurves(sim_om$death, t = t_max, species = "sprat", max_age = 16)[1,50], myGrowthCurves(sim_om$death, t = t_max, species = "herring", max_age = 13)[1,50]))
supp[28, seq(5,14,by=3)] = round(c(myGrowthCurves(sim_op$yay, t = t_max, species = "cod", max_age = 15)[1,50], myGrowthCurves(sim_op$yay, t = t_max, species = "flounder", max_age = 26)[1,50], myGrowthCurves(sim_op$yay, t = t_max, species = "sprat", max_age = 16)[1,50], myGrowthCurves(sim_op$yay, t = t_max, species = "herring", max_age = 13)[1,50]))
supp[29, seq(5,14,by=3)] = round(c(myGrowthCurves(sim_op$okay, t = t_max, species = "cod", max_age = 15)[1,50], myGrowthCurves(sim_op$okay, t = t_max, species = "flounder", max_age = 26)[1,50], myGrowthCurves(sim_op$okay, t = t_max, species = "sprat", max_age = 16)[1,50], myGrowthCurves(sim_op$okay, t = t_max, species = "herring", max_age = 13)[1,50]))
supp[30, seq(5,14,by=3)] = round(c(myGrowthCurves(sim_op$death, t = t_max, species = "cod", max_age = 15)[1,50], myGrowthCurves(sim_op$death, t = t_max, species = "flounder", max_age = 26)[1,50], myGrowthCurves(sim_op$death, t = t_max, species = "sprat", max_age = 16)[1,50], myGrowthCurves(sim_op$death, t = t_max, species = "herring", max_age = 13)[1,50]))
supp[31, seq(5,14,by=3)] = round(c(myGrowthCurves(sim_mp$yay, t = t_max, species = "cod", max_age = 15)[1,50], myGrowthCurves(sim_mp$yay, t = t_max, species = "flounder", max_age = 26)[1,50], myGrowthCurves(sim_mp$yay, t = t_max, species = "sprat", max_age = 16)[1,50], myGrowthCurves(sim_mp$yay, t = t_max, species = "herring", max_age = 13)[1,50]))
supp[32, seq(5,14,by=3)] = round(c(myGrowthCurves(sim_mp$okay, t = t_max, species = "cod", max_age = 15)[1,50], myGrowthCurves(sim_mp$okay, t = t_max, species = "flounder", max_age = 26)[1,50], myGrowthCurves(sim_mp$okay, t = t_max, species = "sprat", max_age = 16)[1,50], myGrowthCurves(sim_mp$okay, t = t_max, species = "herring", max_age = 13)[1,50]))
supp[33, seq(5,14,by=3)] = round(c(myGrowthCurves(sim_mp$death, t = t_max, species = "cod", max_age = 15)[1,50], myGrowthCurves(sim_mp$death, t = t_max, species = "flounder", max_age = 26)[1,50], myGrowthCurves(sim_mp$death, t = t_max, species = "sprat", max_age = 16)[1,50], myGrowthCurves(sim_mp$death, t = t_max, species = "herring", max_age = 13)[1,50]))
supp[34, seq(5,14,by=3)] = round(c(myGrowthCurves(sim_b$yay, t = t_max, species = "cod", max_age = 15)[1,50], myGrowthCurves(sim_b$yay, t = t_max, species = "flounder", max_age = 26)[1,50], myGrowthCurves(sim_b$yay, t = t_max, species = "sprat", max_age = 16)[1,50], myGrowthCurves(sim_b$yay, t = t_max, species = "herring", max_age = 13)[1,50]))
supp[35, seq(5,14,by=3)] = round(c(myGrowthCurves(sim_b$okay, t = t_max, species = "cod", max_age = 15)[1,50], myGrowthCurves(sim_b$okay, t = t_max, species = "flounder", max_age = 26)[1,50], myGrowthCurves(sim_b$okay, t = t_max, species = "sprat", max_age = 16)[1,50], myGrowthCurves(sim_b$okay, t = t_max, species = "herring", max_age = 13)[1,50]))
supp[36, seq(5,14,by=3)] = round(c(myGrowthCurves(sim_b$death, t = t_max, species = "cod", max_age = 15)[1,50], myGrowthCurves(sim_b$death, t = t_max, species = "flounder", max_age = 26)[1,50], myGrowthCurves(sim_b$death, t = t_max, species = "sprat", max_age = 16)[1,50], myGrowthCurves(sim_b$death, t = t_max, species = "herring", max_age = 13)[1,50]))
supp[37, seq(5,14,by=3)] = round(c(myGrowthCurves(sim_o$yay, t = t_max, species = "cod", max_age = 15)[1,50], myGrowthCurves(sim_o$yay, t = t_max, species = "flounder", max_age = 26)[1,50], myGrowthCurves(sim_o$yay, t = t_max, species = "sprat", max_age = 16)[1,50], myGrowthCurves(sim_o$yay, t = t_max, species = "herring", max_age = 13)[1,50]))
supp[38, seq(5,14,by=3)] = round(c(myGrowthCurves(sim_o$okay, t = t_max, species = "cod", max_age = 15)[1,50], myGrowthCurves(sim_o$okay, t = t_max, species = "flounder", max_age = 26)[1,50], myGrowthCurves(sim_o$okay, t = t_max, species = "sprat", max_age = 16)[1,50], myGrowthCurves(sim_o$okay, t = t_max, species = "herring", max_age = 13)[1,50]))
supp[39, seq(5,14,by=3)] = round(c(myGrowthCurves(sim_o$death, t = t_max, species = "cod", max_age = 15)[1,50], myGrowthCurves(sim_o$death, t = t_max, species = "flounder", max_age = 26)[1,50], myGrowthCurves(sim_o$death, t = t_max, species = "sprat", max_age = 16)[1,50], myGrowthCurves(sim_o$death, t = t_max, species = "herring", max_age = 13)[1,50]))
supp[40, seq(5,14,by=3)] = round(c(myGrowthCurves(sim_m$yay, t = t_max, species = "cod", max_age = 15)[1,50], myGrowthCurves(sim_m$yay, t = t_max, species = "flounder", max_age = 26)[1,50], myGrowthCurves(sim_m$yay, t = t_max, species = "sprat", max_age = 16)[1,50], myGrowthCurves(sim_m$yay, t = t_max, species = "herring", max_age = 13)[1,50]))
supp[41, seq(5,14,by=3)] = round(c(myGrowthCurves(sim_m$okay, t = t_max, species = "cod", max_age = 15)[1,50], myGrowthCurves(sim_m$okay, t = t_max, species = "flounder", max_age = 26)[1,50], myGrowthCurves(sim_m$okay, t = t_max, species = "sprat", max_age = 16)[1,50], myGrowthCurves(sim_m$okay, t = t_max, species = "herring", max_age = 13)[1,50]))
supp[42, seq(5,14,by=3)] = round(c(myGrowthCurves(sim_m$death, t = t_max, species = "cod", max_age = 15)[1,50], myGrowthCurves(sim_m$death, t = t_max, species = "flounder", max_age = 26)[1,50], myGrowthCurves(sim_m$death, t = t_max, species = "sprat", max_age = 16)[1,50], myGrowthCurves(sim_m$death, t = t_max, species = "herring", max_age = 13)[1,50]))
supp[43, seq(5,14,by=3)] = round(c(myGrowthCurves(sim_p$yay, t = t_max, species = "cod", max_age = 15)[1,50], myGrowthCurves(sim_p$yay, t = t_max, species = "flounder", max_age = 26)[1,50], myGrowthCurves(sim_p$yay, t = t_max, species = "sprat", max_age = 16)[1,50], myGrowthCurves(sim_p$yay, t = t_max, species = "herring", max_age = 13)[1,50]))
supp[44, seq(5,14,by=3)] = round(c(myGrowthCurves(sim_p$okay, t = t_max, species = "cod", max_age = 15)[1,50], myGrowthCurves(sim_p$okay, t = t_max, species = "flounder", max_age = 26)[1,50], myGrowthCurves(sim_p$okay, t = t_max, species = "sprat", max_age = 16)[1,50], myGrowthCurves(sim_p$okay, t = t_max, species = "herring", max_age = 13)[1,50]))
supp[45, seq(5,14,by=3)] = round(c(myGrowthCurves(sim_p$death, t = t_max, species = "cod", max_age = 15)[1,50], myGrowthCurves(sim_p$death, t = t_max, species = "flounder", max_age = 26)[1,50], myGrowthCurves(sim_p$death, t = t_max, species = "sprat", max_age = 16)[1,50], myGrowthCurves(sim_p$death, t = t_max, species = "herring", max_age = 13)[1,50]))
supp[46, seq(5,14,by=3)] = round(c(myGrowthCurves(sim_none$yay, t = t_max, species = "cod", max_age = 15)[1,50], myGrowthCurves(sim_none$yay, t = t_max, species = "flounder", max_age = 26)[1,50], myGrowthCurves(sim_none$yay, t = t_max, species = "sprat", max_age = 16)[1,50], myGrowthCurves(sim_none$yay, t = t_max, species = "herring", max_age = 13)[1,50]))
supp[47, seq(5,14,by=3)] = round(c(myGrowthCurves(sim_none$okay, t = t_max, species = "cod", max_age = 15)[1,50], myGrowthCurves(sim_none$okay, t = t_max, species = "flounder", max_age = 26)[1,50], myGrowthCurves(sim_none$okay, t = t_max, species = "sprat", max_age = 16)[1,50], myGrowthCurves(sim_none$okay, t = t_max, species = "herring", max_age = 13)[1,50]))
supp[48, seq(5,14,by=3)] = round(c(myGrowthCurves(sim_none$death, t = t_max, species = "cod", max_age = 15)[1,50], myGrowthCurves(sim_none$death, t = t_max, species = "flounder", max_age = 26)[1,50], myGrowthCurves(sim_none$death, t = t_max, species = "sprat", max_age = 16)[1,50], myGrowthCurves(sim_none$death, t = t_max, species = "herring", max_age = 13)[1,50]))

# Create table
print(xtable(supp), include.rownames = F, digits = 0)


########## VI. Get occupancy ##########

# Get all occupancy values
yay_occ = c(sim_bom$yay@params@other_params$benthic_occupancy[,"cod",])
okay_occ = c(sim_bom$okay@params@other_params$benthic_occupancy[,"cod",])
death_occ = c(sim_bom$death@params@other_params$benthic_occupancy[,"cod",])

# Check that all values are equal
var(yay_occ) == 0
var(okay_occ) == 0
var(death_occ) == 0

# If var is zero for all, just grab one value
yay_occ = yay_occ[1]
okay_occ = okay_occ[1]
death_occ = death_occ[1]


########## VII. Get diet ##########

# Get weight categories and weight category widths
w_all = sim_bom$yay@params@w
dw_all = sim_bom$yay@params@dw

# Use mature fish for diet
idx_cod = w_all >= species_params(params_bom)["cod", ]$w_mat
idx_flounder = w_all >= species_params(params_bom)["flounder", ]$w_mat
idx_sprat = w_all >= species_params(params_bom)["sprat", ]$w_mat
idx_herring = w_all >= species_params(params_bom)["herring", ]$w_mat

# Get maintenance costs as a proportion for all species
cod_metab = params_bom@metab["cod",idx_cod] / rowSums(getHypoxiaDiet(sim_bom$yay,t_max+1,proportion=F)["cod",idx_cod,])
flounder_metab = params_bom@metab["flounder",idx_flounder] / rowSums(getHypoxiaDiet(sim_bom$yay,t_max+1,proportion=F)["flounder",idx_flounder,])
sprat_metab = params_bom@metab["sprat",idx_sprat] / rowSums(getHypoxiaDiet(sim_bom$yay,t_max+1,proportion=F)["sprat",idx_sprat,])
herring_metab = params_bom@metab["herring",idx_herring] / rowSums(getHypoxiaDiet(sim_bom$yay,t_max+1,proportion=F)["herring",idx_herring,])

# Get 3 mL/L weighted diet (g/year)
yay_cod = colSums(sweep(getHypoxiaDiet(sim_bom$yay, t_max+1, proportion = F)["cod",idx_cod,],1,cod_metab,"*") * (sim_bom$yay@n[dim(sim_bom$yay@n)[1],,] * w_all * dw_all)["cod",idx_cod]) / sum((sim_bom$yay@n[dim(sim_bom$yay@n)[1],,] * w_all * dw_all)["cod",idx_cod])
yay_flounder = colSums(sweep(getHypoxiaDiet(sim_bom$yay, t_max+1, proportion = F)["flounder",idx_flounder,],1,flounder_metab,"*") * (sim_bom$yay@n[dim(sim_bom$yay@n)[1],,] * w_all * dw_all)["flounder",idx_flounder]) / sum((sim_bom$yay@n[dim(sim_bom$yay@n)[1],,] * w_all * dw_all)["flounder",idx_flounder])
yay_sprat = colSums(sweep(getHypoxiaDiet(sim_bom$yay, t_max+1, proportion = F)["sprat",idx_sprat,],1,sprat_metab,"*") * (sim_bom$yay@n[dim(sim_bom$yay@n)[1],,] * w_all * dw_all)["sprat",idx_sprat]) / sum((sim_bom$yay@n[dim(sim_bom$yay@n)[1],,] * w_all * dw_all)["sprat",idx_sprat])
yay_herring = colSums(sweep(getHypoxiaDiet(sim_bom$yay, t_max+1, proportion = F)["herring",idx_herring,],1,herring_metab,"*") * (sim_bom$yay@n[dim(sim_bom$yay@n)[1],,] * w_all * dw_all)["herring",idx_herring]) / sum((sim_bom$yay@n[dim(sim_bom$yay@n)[1],,] * w_all * dw_all)["herring",idx_herring])

# Relative values for plotting arrows
yayarw_cod = 0.5 * log(yay_cod + 1)
yayarw_flounder = 0.5 * log(yay_flounder + 1)
yayarw_sprat = 0.5 * log(yay_sprat + 1)
yayarw_herring = 0.5 * log(yay_herring + 1)

# Get 2 mL/L weighted diet (g/year)
okay_cod = colSums(sweep(getHypoxiaDiet(sim_bom$okay, t_max+1, proportion = F)["cod",idx_cod,],1,cod_metab,"*") * (sim_bom$okay@n[dim(sim_bom$okay@n)[1],,] * w_all * dw_all)["cod",idx_cod]) / sum((sim_bom$okay@n[dim(sim_bom$okay@n)[1],,] * w_all * dw_all)["cod",idx_cod])
okay_flounder = colSums(sweep(getHypoxiaDiet(sim_bom$okay, t_max+1, proportion = F)["flounder",idx_flounder,],1,flounder_metab,"*") * (sim_bom$okay@n[dim(sim_bom$okay@n)[1],,] * w_all * dw_all)["flounder",idx_flounder]) / sum((sim_bom$okay@n[dim(sim_bom$okay@n)[1],,] * w_all * dw_all)["flounder",idx_flounder])
okay_sprat = colSums(sweep(getHypoxiaDiet(sim_bom$okay, t_max+1, proportion = F)["sprat",idx_sprat,],1,sprat_metab,"*") * (sim_bom$okay@n[dim(sim_bom$okay@n)[1],,] * w_all * dw_all)["sprat",idx_sprat]) / sum((sim_bom$okay@n[dim(sim_bom$okay@n)[1],,] * w_all * dw_all)["sprat",idx_sprat])
okay_herring = colSums(sweep(getHypoxiaDiet(sim_bom$okay, t_max+1, proportion = F)["herring",idx_herring,],1,herring_metab,"*") * (sim_bom$okay@n[dim(sim_bom$okay@n)[1],,] * w_all * dw_all)["herring",idx_herring]) / sum((sim_bom$okay@n[dim(sim_bom$okay@n)[1],,] * w_all * dw_all)["herring",idx_herring])

# Relative values for plotting arrows
okayarw_cod = 0.5 * log(okay_cod + 1)
okayarw_flounder = 0.5 * log(okay_flounder + 1)
okayarw_sprat = 0.5 * log(okay_sprat + 1)
okayarw_herring = 0.5 * log(okay_herring + 1)

# Get 1 mL/L weighted diet (g/year)
death_cod = colSums(sweep(getHypoxiaDiet(sim_bom$death, t_max+1, proportion = F)["cod",idx_cod,],1,cod_metab,"*") * (sim_bom$death@n[dim(sim_bom$death@n)[1],,] * w_all * dw_all)["cod",idx_cod]) / sum((sim_bom$death@n[dim(sim_bom$death@n)[1],,] * w_all * dw_all)["cod",idx_cod])
death_flounder = colSums(sweep(getHypoxiaDiet(sim_bom$death, t_max+1, proportion = F)["flounder",idx_flounder,],1,flounder_metab,"*") * (sim_bom$death@n[dim(sim_bom$death@n)[1],,] * w_all * dw_all)["flounder",idx_flounder]) / sum((sim_bom$death@n[dim(sim_bom$death@n)[1],,] * w_all * dw_all)["flounder",idx_flounder])
death_sprat = colSums(sweep(getHypoxiaDiet(sim_bom$death, t_max+1, proportion = F)["sprat",idx_sprat,],1,sprat_metab,"*") * (sim_bom$death@n[dim(sim_bom$death@n)[1],,] * w_all * dw_all)["sprat",idx_sprat]) / sum((sim_bom$death@n[dim(sim_bom$death@n)[1],,] * w_all * dw_all)["sprat",idx_sprat])
death_herring = colSums(sweep(getHypoxiaDiet(sim_bom$death, t_max+1, proportion = F)["herring",idx_herring,],1,herring_metab,"*") * (sim_bom$death@n[dim(sim_bom$death@n)[1],,] * w_all * dw_all)["herring",idx_herring]) / sum((sim_bom$death@n[dim(sim_bom$death@n)[1],,] * w_all * dw_all)["herring",idx_herring])

# Relative values for plotting arrows
deatharw_cod = 0.5 * log(death_cod + 1)
deatharw_flounder = 0.5 * log(death_flounder + 1)
deatharw_sprat = 0.5 * log(death_sprat + 1)
deatharw_herring = 0.5 * log(death_herring + 1)

# Colors
okay_all = c((okay_cod/yay_cod)[1:5], (okay_flounder/yay_flounder)[c(2,5)], (okay_sprat/yay_sprat)[6], (okay_herring/yay_herring)[6])
death_all = c((death_cod/yay_cod)[1:5], (death_flounder/yay_flounder)[c(2,5)], (death_sprat/yay_sprat)[6], (death_herring/yay_herring)[6])
rc = color_values(c(okay_all, death_all, 3.5, -1.5), palette = "rdbu")

# Relative fish sizes
wt_cod = scen_table[scen_table$Model == "BOM", "Cod_Size"] / max(scen_table[scen_table$Model == "BOM", "Cod_Size"])
l_cod = (scen_table[scen_table$Model == "BOM", "Cod_Size"] / cod_lh["a","stan"]) ^ (1/cod_lh["b","stan"]); l_cod = l_cod / max(l_cod)
wt_sprat = scen_table[scen_table$Model == "BOM", "Sprat_Size"] / max(scen_table[scen_table$Model == "BOM", "Sprat_Size"])
l_sprat = (scen_table[scen_table$Model == "BOM", "Sprat_Size"] / sprat_lh["a","stan"]) ^ (1/sprat_lh["b","stan"]); l_sprat = l_sprat / max(l_sprat)
wt_herring = scen_table[scen_table$Model == "BOM", "Herring_Size"] / max(scen_table[scen_table$Model == "BOM", "Herring_Size"])
l_herring = (scen_table[scen_table$Model == "BOM", "Herring_Size"] / herring_lh["a","stan"]) ^ (1/herring_lh["b","stan"]); l_herring = l_herring / max(l_herring)
wt_flounder = scen_table[scen_table$Model == "BOM", "Flounder_Size"] / max(scen_table[scen_table$Model == "BOM", "Flounder_Size"])
l_flounder = (scen_table[scen_table$Model == "BOM", "Flounder_Size"] / flounder_lh["a","stan"]) ^ (1/flounder_lh["b","stan"]); l_flounder = l_flounder / max(l_flounder)


########## VIII. Plot food web diagrams ##########

# Initialize plot
jpeg("Plots/foodweb_scenarios.jpg", width = 20, height = 14, units = 'cm', res = 600)

# Read drawings
sad_drawing = readPNG("Plots/bd_tp.png")
cop_drawing = readPNG("Plots/pd_tp.png")
cod_drawing = readPNG("Plots/cd_tp.png")
fln_drawing = readPNG("Plots/fd_tp.png")
spt_drawing = readPNG("Plots/sd_tp.png")
her_drawing = readPNG("Plots/hd_tp.png")

# Layout plot
par(mfrow = c(1,3), mar = c(2,2,2,2), oma = c(6,0,3,0))

# Yay scenario
plot(1, type = 'n', xlim = c(0,1), ylim = c(0,1), axes = F, xlab = "", ylab = "", main = "3 mL/L", cex.main = 2); box(lwd = 3)
waterplot(100, "Blues 2")
abline(0.5, 0, lwd = 3)
rasterImage(sad_drawing, 0.08, 0.00, 0.23, 0.15)
rasterImage(cop_drawing, 0.42, 0.90, 0.58, 1.00)
rasterImage(cod_drawing, 0.20+((0.80-0.20)-l_cod[1]*(0.80-0.20))/2, 0.50-(0.20*wt_cod[1])*yay_occ, 0.80-((0.80-0.20)-l_cod[1]*(0.80-0.20))/2, 0.50+(0.20*wt_cod[1])*(1-yay_occ))
rasterImage(fln_drawing, 0.62+((0.92-0.62)-l_flounder[1]*(0.92-0.62))/2, 0.00+((0.15-0.00)-wt_flounder[1]*(0.15-0.00))/2, 0.92-((0.92-0.62)-l_flounder[1]*(0.92-0.62))/2, 0.15-((0.15-0.00)-wt_flounder[1]*(0.15-0.00))/2)
rasterImage(spt_drawing, 0.08+((0.23-0.08)-l_sprat[1]*(0.23-0.08))/2, 0.75+((0.80-0.75)-wt_sprat[1]*(0.80-0.75))/2, 0.23-((0.23-0.08)-l_sprat[1]*(0.23-0.08))/2, 0.80-((0.80-0.75)-wt_sprat[1]*(0.80-0.75))/2)
rasterImage(her_drawing, 0.67+((0.92-0.67)-l_herring[1]*(0.92-0.67))/2, 0.74+((0.81-0.74)-wt_herring[1]*(0.81-0.74))/2, 0.92-((0.92-0.67)-l_sprat[1]*(0.92-0.67))/2, 0.81-((0.81-0.74)-wt_herring[1]*(0.81-0.74))/2)
selfarrow(c(0.83,0.42), path = "R", curve = c(0.03,0.03), lwd = yayarw_cod["cod"]*1.5, arr.lwd = yayarw_cod["cod"]*1.5, lcol = rgb(0,0,0), arr.pos = 0.5, arr.length = 0.2, arr.type = "triangle")
selfarrow(c(0.83,0.42), path = "R", curve = c(0.03,0.03), lwd = yayarw_cod["cod"], arr.lwd = yayarw_cod["cod"], lcol = rgb(1,1,1), arr.pos = 0.5, arr.length = 0.2, arr.type = "triangle")
polygon(c(0.83,0.83-yayarw_cod["cod"]*0.0075,0.83+yayarw_cod["cod"]*0.0075), c(0.395,0.437,0.437), col = rgb(1,1,1), lwd = 1.5)
arrows(0.75, 0.16, 0.65, 0.30, length = 0.1, lwd = yayarw_cod["flounder"]*1.2, col = rgb(0,0,0))
arrows(0.75, 0.16, 0.65, 0.30, length = 0.1, lwd = yayarw_cod["flounder"], col = rgb(1,1,1))
arrows(0.18, 0.72, 0.30, 0.52, length = 0.1, lwd = yayarw_cod["sprat"]*1.2, col = rgb(0,0,0))
arrows(0.18, 0.72, 0.30, 0.52, length = 0.1, lwd = yayarw_cod["sprat"], col = rgb(1,1,1))
arrows(0.82, 0.72, 0.70, 0.52, length = 0.1, lwd = yayarw_cod["herring"]*1.2, col = rgb(0,0,0))
arrows(0.82, 0.72, 0.70, 0.52, length = 0.1, lwd = yayarw_cod["herring"], col = rgb(1,1,1))
arrows(0.25, 0.16, 0.35, 0.30, length = 0.1, lwd = yayarw_cod["Benthic"]*1.2, col = rgb(0,0,0))
arrows(0.25, 0.16, 0.35, 0.30, length = 0.1, lwd = yayarw_cod["Benthic"], col = rgb(1,1,1))
selfarrow(c(0.93,0.07), path = "R", curve = c(0.03,0.03), lwd = yayarw_flounder["flounder"]*1.5, arr.lwd = yayarw_flounder["flounder"], lcol = rgb(0,0,0), arr.pos = 0.5, arr.length = 0.2, arr.type = "triangle")
selfarrow(c(0.93,0.07), path = "R", curve = c(0.03,0.03), lwd = yayarw_flounder["flounder"], arr.lwd = yayarw_flounder["flounder"], lcol = rgb(1,1,1), arr.pos = 0.5, arr.length = 0.2, arr.type = "triangle")
polygon(c(0.93,0.93-yayarw_flounder["flounder"]*0.05,0.93+yayarw_flounder["flounder"]*0.05), c(0.060,0.082,0.082), col = rgb(1,1,1), lwd = 1.5)
arrows(0.25, 0.09, 0.55, 0.09, length = 0.1, lwd = yayarw_flounder["Benthic"]*1.2, col = rgb(0,0,0))
arrows(0.25, 0.09, 0.55, 0.09, length = 0.1, lwd = yayarw_flounder["Benthic"], col = rgb(1,1,1))
arrows(0.40, 0.92, 0.20, 0.82, length = 0.1, lwd = yayarw_sprat["Pelagic"]*1.2, col = rgb(0,0,0))
arrows(0.40, 0.92, 0.20, 0.82, length = 0.1, lwd = yayarw_sprat["Pelagic"], col = rgb(1,1,1))
arrows(0.60, 0.92, 0.80, 0.82, length = 0.1, lwd = yayarw_herring["Pelagic"]*1.2, col = rgb(0,0,0))
arrows(0.60, 0.92, 0.80, 0.82, length = 0.1, lwd = yayarw_herring["Pelagic"], col = rgb(1,1,1))

# Okay scenario
plot(1, type = 'n', xlim = c(0,1), ylim = c(0,1), axes = F, xlab = "", ylab = "", main = "2 mL/L", cex.main = 2); box(lwd = 3)
waterplot(100, "Blues 2")
abline(0.5, 0, lwd = 3)
rasterImage(sad_drawing, 0.08, 0.00, 0.23, 0.15)
rasterImage(cop_drawing, 0.42, 0.90, 0.58, 1.00)
rasterImage(cod_drawing, 0.20+((0.80-0.20)-l_cod[2]*(0.80-0.20))/2, 0.50-(0.20*wt_cod[2])*okay_occ, 0.80-((0.80-0.20)-l_cod[2]*(0.80-0.20))/2, 0.50+(0.20*wt_cod[2])*(1-okay_occ))
rasterImage(fln_drawing, 0.62+((0.92-0.62)-l_flounder[2]*(0.92-0.62))/2, 0.00+((0.15-0.00)-wt_flounder[2]*(0.15-0.00))/2, 0.92-((0.92-0.62)-l_flounder[2]*(0.92-0.62))/2, 0.15-((0.15-0.00)-wt_flounder[2]*(0.15-0.00))/2)
rasterImage(spt_drawing, 0.08+((0.23-0.08)-l_sprat[2]*(0.23-0.08))/2, 0.75+((0.80-0.75)-wt_sprat[2]*(0.80-0.75))/2, 0.23-((0.23-0.08)-l_sprat[2]*(0.23-0.08))/2, 0.80-((0.80-0.75)-wt_sprat[2]*(0.80-0.75))/2)
rasterImage(her_drawing, 0.67+((0.92-0.67)-l_herring[2]*(0.92-0.67))/2, 0.74+((0.81-0.74)-wt_herring[2]*(0.81-0.74))/2, 0.92-((0.92-0.67)-l_herring[2]*(0.92-0.67))/2, 0.81-((0.81-0.74)-wt_herring[2]*(0.81-0.74))/2)
selfarrow(c(0.83,0.44), path = "R", curve = c(0.03,0.03), lwd = okayarw_cod["cod"]*1.5, arr.lwd = okayarw_cod["cod"]*1.5, lcol = rgb(0,0,0), arr.pos = 0.5, arr.length = 0.2, arr.type = "triangle")
selfarrow(c(0.83,0.44), path = "R", curve = c(0.03,0.03), lwd = okayarw_cod["cod"], arr.lwd = okayarw_cod["cod"], lcol = rc[1], arr.pos = 0.5, arr.length = 0.2, arr.type = "triangle")
polygon(c(0.83,0.83-okayarw_cod["cod"]*0.0075,0.83+okayarw_cod["cod"]*0.0075), c(0.415,0.457,0.457), col = rc[1], lwd = 1.5)
arrows(0.75, 0.16, 0.65, 0.32, length = 0.1, lwd = okayarw_cod["flounder"]*1.2, col = rgb(0,0,0))
arrows(0.75, 0.16, 0.65, 0.32, length = 0.1, lwd = okayarw_cod["flounder"], col = rc[2])
arrows(0.18, 0.72, 0.30, 0.54, length = 0.1, lwd = okayarw_cod["sprat"]*1.2, col = rgb(0,0,0))
arrows(0.18, 0.72, 0.30, 0.54, length = 0.1, lwd = okayarw_cod["sprat"], col = rc[3])
arrows(0.82, 0.72, 0.70, 0.54, length = 0.1, lwd = okayarw_cod["herring"]*1.2, col = rgb(0,0,0))
arrows(0.82, 0.72, 0.70, 0.54, length = 0.1, lwd = okayarw_cod["herring"], col = rc[4])
arrows(0.25, 0.16, 0.35, 0.32, length = 0.1, lwd = okayarw_cod["Benthic"]*1.2, col = rgb(0,0,0))
arrows(0.25, 0.16, 0.35, 0.32, length = 0.1, lwd = okayarw_cod["Benthic"], col = rc[5])
selfarrow(c(0.93,0.07), path = "R", curve = c(0.03,0.03), lwd = okayarw_flounder["flounder"]*1.5, arr.lwd = okayarw_flounder["flounder"]*1.5, lcol = rgb(0,0,0), arr.pos = 0.5, arr.length = 0.2, arr.type = "triangle")
selfarrow(c(0.93,0.07), path = "R", curve = c(0.03,0.03), lwd = okayarw_flounder["flounder"], arr.lwd = okayarw_flounder["flounder"], lcol = rc[6], arr.pos = 0.5, arr.length = 0.2, arr.type = "triangle")
polygon(c(0.93,0.93-okayarw_flounder["flounder"]*0.05,0.93+okayarw_flounder["flounder"]*0.05), c(0.060,0.082,0.082), col = rc[6], lwd = 1.5)
arrows(0.25, 0.09, 0.55, 0.09, length = 0.1, lwd = okayarw_flounder["Benthic"]*1.2, col = rgb(0,0,0))
arrows(0.25, 0.09, 0.55, 0.09, length = 0.1, lwd = okayarw_flounder["Benthic"], col = rc[7])
arrows(0.40, 0.92, 0.20, 0.82, length = 0.1, lwd = okayarw_sprat["Pelagic"]*1.2, col = rgb(0,0,0))
arrows(0.40, 0.92, 0.20, 0.82, length = 0.1, lwd = okayarw_sprat["Pelagic"], col = rc[8])
arrows(0.60, 0.92, 0.80, 0.82, length = 0.1, lwd = okayarw_herring["Pelagic"]*1.2, col = rgb(0,0,0))
arrows(0.60, 0.92, 0.80, 0.82, length = 0.1, lwd = okayarw_herring["Pelagic"], col = rc[9])

# Death scenario
plot(1, type = 'n', xlim = c(0,1), ylim = c(0,1), axes = F, xlab = "", ylab = "", main = "1 mL/L", cex.main = 2); box(lwd = 3)
waterplot(100, "Blues 2")
abline(0.5, 0, lwd = 3)
rasterImage(sad_drawing, 0.08, 0.00, 0.23, 0.15)
rasterImage(cop_drawing, 0.42, 0.90, 0.58, 1.00)
rasterImage(cod_drawing, 0.20+((0.80-0.20)-l_cod[3]*(0.80-0.20))/2, 0.50-(0.20*wt_cod[3])*death_occ, 0.80-((0.80-0.20)-l_cod[3]*(0.80-0.20))/2, 0.50+(0.20*wt_cod[3])*(1-death_occ))
rasterImage(fln_drawing, 0.62+((0.92-0.62)-l_flounder[3]*(0.92-0.62))/2, 0.00+((0.15-0.00)-wt_flounder[3]*(0.15-0.00))/2, 0.92-((0.92-0.62)-l_flounder[3]*(0.92-0.62))/2, 0.15-((0.15-0.00)-wt_flounder[3]*(0.15-0.00))/2)
rasterImage(spt_drawing, 0.08+((0.23-0.08)-l_sprat[3]*(0.23-0.08))/2, 0.75+((0.80-0.75)-wt_sprat[3]*(0.80-0.75))/2, 0.23-((0.23-0.08)-l_sprat[3]*(0.23-0.08))/2, 0.80-((0.80-0.75)-wt_sprat[3]*(0.80-0.75))/2)
rasterImage(her_drawing, 0.67+((0.92-0.67)-l_herring[3]*(0.92-0.67))/2, 0.74+((0.81-0.74)-wt_herring[3]*(0.81-0.74))/2, 0.92-((0.92-0.67)-l_sprat[3]*(0.92-0.67))/2, 0.81-((0.81-0.74)-wt_herring[3]*(0.81-0.74))/2)
selfarrow(c(0.83,0.46), path = "R", curve = c(0.03,0.03), lwd = deatharw_cod["cod"]*1.5, arr.lwd = deatharw_cod["cod"]*1.5, lcol = rgb(0,0,0), arr.pos = 0.5, arr.length = 0.2, arr.type = "triangle")
selfarrow(c(0.83,0.46), path = "R", curve = c(0.03,0.03), lwd = deatharw_cod["cod"], arr.lwd = deatharw_cod["cod"], lcol = rc[10], arr.pos = 0.5, arr.length = 0.2, arr.type = "triangle")
polygon(c(0.83,0.83-deatharw_cod["cod"]*0.0075,0.83+deatharw_cod["cod"]*0.0075), c(0.435,0.477,0.477), col = rc[10], lwd = 1.5)
arrows(0.75, 0.16, 0.65, 0.39, length = 0.1, lwd = deatharw_cod["flounder"]*1.2, col = rgb(0,0,0))
arrows(0.75, 0.16, 0.65, 0.39, length = 0.1, lwd = deatharw_cod["flounder"], col = rc[11])
arrows(0.18, 0.72, 0.30, 0.55, length = 0.1, lwd = deatharw_cod["sprat"]*1.2, col = rgb(0,0,0))
arrows(0.18, 0.72, 0.30, 0.55, length = 0.1, lwd = deatharw_cod["sprat"], col = rc[12])
arrows(0.82, 0.72, 0.70, 0.55, length = 0.1, lwd = deatharw_cod["herring"]*1.2, col = rgb(0,0,0))
arrows(0.82, 0.72, 0.70, 0.55, length = 0.1, lwd = deatharw_cod["herring"], col = rc[13])
arrows(0.25, 0.16, 0.35, 0.39, length = 0.1, lwd = deatharw_cod["Benthic"]*1.2, col = rgb(0,0,0))
arrows(0.25, 0.16, 0.35, 0.39, length = 0.1, lwd = deatharw_cod["Benthic"], col = rc[14])
selfarrow(c(0.93,0.07), path = "R", curve = c(0.03,0.03), lwd = deatharw_flounder["flounder"]*1.5, arr.lwd = deatharw_flounder["flounder"]*1.5, lcol = rgb(0,0,0), arr.pos = 0.5, arr.length = 0.2, arr.type = "triangle")
selfarrow(c(0.93,0.07), path = "R", curve = c(0.03,0.03), lwd = deatharw_flounder["flounder"], arr.lwd = deatharw_flounder["flounder"], lcol = rc[15], arr.pos = 0.5, arr.length = 0.2, arr.type = "triangle")
polygon(c(0.93,0.93-deatharw_flounder["flounder"]*0.05,0.93+deatharw_flounder["flounder"]*0.05), c(0.060,0.082,0.082), col = rc[15], lwd = 1.5)
arrows(0.25, 0.09, 0.55, 0.09, length = 0.1, lwd = deatharw_flounder["Benthic"]*1.2, col = rgb(0,0,0))
arrows(0.25, 0.09, 0.55, 0.09, length = 0.1, lwd = deatharw_flounder["Benthic"], col = rc[16])
arrows(0.40, 0.92, 0.20, 0.82, length = 0.1, lwd = deatharw_sprat["Pelagic"]*1.2, col = rgb(0,0,0))
arrows(0.40, 0.92, 0.20, 0.82, length = 0.1, lwd = deatharw_sprat["Pelagic"], col = rc[17])
arrows(0.60, 0.92, 0.80, 0.82, length = 0.1, lwd = deatharw_herring["Pelagic"]*1.2, col = rgb(0,0,0))
arrows(0.60, 0.92, 0.80, 0.82, length = 0.1, lwd = deatharw_herring["Pelagic"], col = rc[18])

# Add title
mtext("Oxygen", 3, line = 1, outer = T, cex = 1.5)

# Add title
par(fig = c(0,1,0,1), oma = c(0,0,0,0), mar = c(0,0,0,0), new = T)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
par(mfrow = c(1,1))

# Add color scale
image.plot(zlim = c(0,3.5), legend.only = T, horizontal = T, smallplot = c(0.10,0.90,0.08,0.11), cex = 0.1, col = hcl.colors(1001, "rdbu")[301:1000], border = NULL)
text(0, -0.80, "Proportion change")

# Finish plot
dev.off()


########## IX. Plot diet changes ##########

# Get diets
diet_yay = getHypoxiaDiet(sim_bom$yay, t_max)["cod",,]
diet_okay = getHypoxiaDiet(sim_bom$okay, t_max)["cod",,]
diet_death = getHypoxiaDiet(sim_bom$death, t_max)["cod",,]
colnames(diet_yay) = c("cod", "flounder", "sprat", "herring", "benthos", "plankton")
colnames(diet_okay) = c("cod", "flounder", "sprat", "herring", "benthos", "plankton")
colnames(diet_death) = c("cod", "flounder", "sprat", "herring", "benthos", "plankton")
        
# Construct diet data frame
prey = dimnames(diet_yay)$prey
prey = factor(prey, levels = rev(prey))
plot_yay = data.frame(w = params_bom@w, Proportion = c(diet_yay), Prey = rep(prey, each = length(params_bom@w)))
plot_okay = data.frame(w = params_bom@w, Proportion = c(diet_okay), Prey = rep(prey, each = length(params_bom@w)))
plot_death = data.frame(w = params_bom@w, Proportion = c(diet_death), Prey = rep(prey, each = length(params_bom@w)))
plot_yay = plot_yay[plot_yay$Proportion > 0.001, ]
plot_okay = plot_okay[plot_okay$Proportion > 0.001, ]
plot_death = plot_death[plot_death$Proportion > 0.001, ]
        
# Plot colors
color_benthos = rgb(38, 24, 95, maxColorValue=255)
color_pelagic = rgb(0, 106, 168, maxColorValue = 255)
color_cod = rgb(0, 166, 174, maxColorValue = 255)
color_flounder = rgb(252, 255, 221, maxColorValue = 255)
color_sprat = rgb(205, 240, 203, maxColorValue = 255)
color_herring = rgb(119, 209, 181, maxColorValue = 255)
        
# Plot diet
jpeg("Plots/diet_scenarios.jpg", width = 20, height = 8, units = 'cm', res = 600)
legend_colors = c(benthos=color_benthos, plankton=color_pelagic, cod=color_cod, flounder=color_flounder, sprat=color_sprat, herring=color_herring)
legend_which = intersect(plot_yay$Prey, names(legend_colors))
g1 = ggplot(plot_yay) + 
	geom_area(aes(x = w, y = Proportion, fill = Prey)) + 
	scale_x_log10() + labs(x = "", title = "3 mL/L") + 
	theme(plot.title = element_text(hjust = 0.5)) + 
	scale_fill_manual(values = legend_colors[legend_which])
g2 = ggplot(plot_okay) + 
	geom_area(aes(x = w, y = Proportion, fill = Prey)) + 
	scale_x_log10() + 
	labs(x = "Size [g]", y = "", title = "2 mL/L") + scale_fill_manual(values = legend_colors[legend_which]) +
	theme(plot.title = element_text(hjust = 0.5))
g3 = ggplot(plot_death) + 
	geom_area(aes(x = w, y = Proportion, fill = Prey)) + 
	scale_x_log10() + 
	labs(x = "", y = "", title = "1 mL/L") + scale_fill_manual(values = legend_colors[legend_which]) +
	theme(plot.title = element_text(hjust = 0.5))
fig = ggarrange(g1, g2, g3, ncol = 3, common.legend = T, legend = "right") +
	theme(plot.margin = margin(0.4,0,0,0, "cm"))
fig = annotate_figure(fig, top = text_grob("Oxygen", size = 14)) #+
	geom_segment(aes(x = 0.2, y = 0.9, xend = 0.75, yend = 0.9), size = 1, lineend = "round", linejoin = "round", arrow = arrow(length = unit(0.2, "cm")))
plot(fig)
dev.off()
