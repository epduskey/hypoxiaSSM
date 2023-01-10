# Fitting growth

# A function to fit growth with h or gamma
#	par: input parameters c(U_crit, a_crit, U_met, a_met) in that order
#	params: mizer params object
#	species: which species to fit
#	yr: years to which to fit
#	t: time steps of model to run
#	effort: a vector containing effort values for cod, flounder, sprat, and herring in that order
#	sim: mizer simulation object, if already run
#	returns sum of relative errors in all species' SSB 
time_growth = function(par, params, species, yr, t, effort, sim = NULL) {
	
	# New params object
	ret = params
	
	# Other variables
	uc_cod = species_params(ret)["cod",]$U_crit
	uc_flounder = species_params(ret)["flounder",]$U_crit
	ac_cod = species_params(ret)["cod",]$a_crit
	ac_flounder = species_params(ret)["flounder",]$a_crit
	um_cod = species_params(ret)["cod",]$U_met
	um_flounder = species_params(ret)["flounder",]$U_met
	am_cod = species_params(ret)["cod",]$a_met
	am_flounder = species_params(ret)["flounder",]$a_met
	zm_cod = species_params(ret)["cod",]$z_mort
	zm_flounder = species_params(ret)["flounder",]$z_mort
	bm_cod = species_params(ret)["cod",]$b_mort
	bm_flounder = species_params(ret)["flounder",]$z_mort
	
	# Set up occupancy optim object if absent
	if(!("occ.optim" %in% ls())) { 
		occ.optim = list()
		occ.optim$par = c(10^5, 0)
	}
	
	# Run simulation if not already provided
	if(is.null(sim)) {
		# Assign parameters
		if(species == "cod") {
			ret = oxy_sensitivity(ret, 
				U_hab = c(occ.optim$par[1], 10^5, -10^5, -10^5),
				a_hab = c(occ.optim$par[2], 0, 0, 0),
				U_crit = c(par[1], uc_flounder, 10^5, 10^5),
				a_crit = c(par[2], ac_flounder, 0, 0),
				U_met = c(par[3], um_flounder, 10^5, 10^5),
				a_met = c(par[4], am_flounder, -10^5, -10^5),
				z_mort = c(zm_cod, zm_flounder, 0, 0),
				b_mort = c(bm_cod, bm_flounder, 0, 0))
		} else if(species == "flounder") {
			ret = oxy_sensitivity(ret, 
				U_hab = c(occ.optim$par[1], 10^5, -10^5, -10^5),
				a_hab = c(occ.optim$par[2], 0, 0, 0),
				U_crit = c(uc_cod, par[1], 10^5, 10^5),
				a_crit = c(ac_cod, par[2], 0, 0),
				U_met = c(um_cod, par[3], 10^5, 10^5),
				a_met = c(am_cod, par[4], -10^5, -10^5),
				z_mort = c(zm_cod, zm_flounder, 0, 0),
				b_mort = c(bm_cod, bm_flounder, 0, 0))
		}
		
		# Create an oxygen vector
		times = 0:t
		benthic_oxygen = vector(mode = "numeric", length = length(times))
		pelagic_oxygen = vector(mode = "numeric", length = length(times))

		# Time series oxygen
		benthic_oxygen[1:length(times)] = c(rep(mean(preds_all$df[preds_all$df$Year %in% yr, ]$Oxygen_benthic),t-length(yr)+1), preds_all$df[preds_all$df$Year %in% yr, ]$Oxygen_benthic)
		pelagic_oxygen[1:length(times)] = c(rep(mean(preds_all$df[preds_all$df$Year %in% yr, ]$Oxygen_pelagic),t-length(yr)+1), preds_all$df[preds_all$df$Year %in% yr, ]$Oxygen_pelagic)

		# Store these in other_params
		ret@other_params$benthic_oxygen = benthic_oxygen
		ret@other_params$pelagic_oxygen = pelagic_oxygen

		# Choose reference temperature
		Tref = 19
		ret@other_params$Tref = Tref

		# Create a temperature vector
		benthic_temp = vector(mode = "numeric", length = length(times))
		pelagic_temp = vector(mode = "numeric", length = length(times))

		# Constant temperature scenario
		benthic_temp[1:length(times)] = 19
		pelagic_temp[1:length(times)] = 19

		# Store these in other_params
		ret@other_params$benthic_temp = benthic_temp
		ret@other_params$pelagic_temp = pelagic_temp

		# Scale rates by oxygen and temperature: see otmscale.R
		ret = rate_scale(ret, t)
		ret = mscale(ret, t)
	
		# Run to steady state
		sim = project(ret, t_max = t, effort = effort)
	}
		
	# Growth matrices
	growth_mod = array(dim = c(length(yr), length(species), 50))
	species_all = c("cod","flounder","sprat","herring")
	dimnames(growth_mod) = list(Year = yr, Species = species_all[species_all %in% species], Age = seq(0,1,length.out=50))
	sp_maxage = c(15,26,16,13)
	growth_obs = growth_mod
	
	# Growth data
	sp_obs = list(cod = cod_ivb_cal, flounder = flounder_ivb_cal, sprat = sprat_ivb_cal, herring = herring_ivb_cal)[species]
	
	for(i in 1:dim(growth_mod)[1]) {
		
		# Model growth
	
		# von Bertalanffy growth
		sp = species_params(ret)[rownames(species_params(ret)) %in% species, ]
		L_inf = (sp$w_inf/sp$a)^(1/sp$b)
		for(j in 1:dim(growth_obs)[2]) {
			growth_mod[i,j,] = myGrowthCurves(sim, (t-length(yr))+i, species = rownames(sp)[j], max_age = sp_maxage[j], percentage = T)
			growth_obs[i,j,] = sp_obs[[j]][i,]$a * (sp_obs[[j]][i,]$L_inf * (1 - exp(-sp_obs[[j]][i,]$k * (as.numeric(dimnames(growth_mod)$Age)*sp_maxage[j])))) ^ sp_obs[[j]][i,]$b
			growth_obs[i,j,] = (growth_obs[i,j,]/sp$w_inf[j])*100
		}
	}
		
	error = sum((growth_mod - growth_obs)^2)
	return(error)
}

# A function to return growth error
#	params: mizer params object
#	species: which species to fit
#	yr: years to which to fit
#	t: time steps of model to run
#	effort: a vector containing effort values for cod, flounder, sprat, and herring in that order
#	sim: mizer simulation object, if already run
#	returns sum of relative errors in all species' SSB 
proj_growth = function(params, yr, t, species, sim = NULL) {
	
	# Run simulation if not already provided
	if(is.null(sim)) {
		# Create effort array
		effort = farray(params, fdat, yr, c(2,2,2,2), t)

		# Run to steady state
		sim = project(params, t_max = t, effort = effort)
	}
		
	# Growth matrices
	growth_mod = array(dim = c(length(yr), length(species), 50))
	dimnames(growth_mod) = list(Year = yr, Species = species, Age = seq(0,1,length.out=50))
	sp_maxage = c(15,26,16,13)
	growth_obs = growth_mod
	
	# Growth data
	sp_obs = list(cod = rbind(cod_ivb_cal,cod_ivb_prj), flounder = rbind(flounder_ivb_cal,flounder_ivb_prj), sprat = rbind(sprat_ivb_cal,sprat_ivb_prj), herring = rbind(herring_ivb_cal,herring_ivb_prj))[species]
	
	for(i in 1:dim(growth_mod)[1]) {
		
		# Model growth
	
		# von Bertalanffy growth
		sp = species_params(params)[rownames(species_params(params)) %in% species, ]
		L_inf = (sp$w_inf/sp$a)^(1/sp$b)
		for(j in 1:dim(growth_obs)[2]) {
			sp_temp = sp_obs[[j]][sp_obs[[j]]$Year %in% yr[i], 2:6]
			growth_mod[i,j,] = myGrowthCurves(sim, (t-length(yr))+i, species = rownames(sp)[j], max_age = sp_maxage[j], percentage = T)
			if(sum(is.na(sp_temp)) == 0) {
				growth_obs[i,j,] = sp_temp$a * (sp_temp$L_inf * (1 - exp(-sp_temp$k * (as.numeric(dimnames(growth_mod)$Age)*sp_maxage[j])))) ^ sp_temp$b
			} else {
				growth_obs[i,j,] = NA
				growth_mod[i,j,] = NA
			}
			growth_obs[i,j,] = (growth_obs[i,j,]/sp$w_inf[j])*100
		}
	}
		
	error = sum((growth_mod - growth_obs)^2, na.rm = T)
	return(error)
}

# A function to fit growth with h or gamma
#	par: benthic oxygen sensitivity input parameters c(U_oxy, k_oxy) in that order
#	params: mizer params object
#	yr: years to which to fit
#	t: time steps of model to run
#	effort: a vector containing effort values for cod, flounder, sprat, and herring in that order
#	sim: mizer simulation object, if already run
#	returns sum of relative errors in all species' growth over time 
benthic_growth = function(par, params, yr, t, effort, sim = NULL) {
	
	# New params object
	ret = params
	
	# Run simulation if not already provided
	if(is.null(sim)) {	
		# Assign parameters
		resource_params(ret)$U_oxy = par[1]
		resource_params(ret)$k_oxy = par[2]
		
		# Create an oxygen vector
		times = 0:t
		benthic_oxygen = vector(mode = "numeric", length = length(times))
		pelagic_oxygen = vector(mode = "numeric", length = length(times))

		# Time series oxygen
		benthic_oxygen[1:length(times)] = c(rep(mean(preds_all$df[preds_all$df$Year %in% yr, ]$Oxygen_benthic),t-length(yr)+1), preds_all$df[preds_all$df$Year %in% yr, ]$Oxygen_benthic)
		pelagic_oxygen[1:length(times)] = c(rep(mean(preds_all$df[preds_all$df$Year %in% yr, ]$Oxygen_pelagic),t-length(yr)+1), preds_all$df[preds_all$df$Year %in% yr, ]$Oxygen_pelagic)

		# Store these in other_params
		ret@other_params$benthic_oxygen = benthic_oxygen
		ret@other_params$pelagic_oxygen = pelagic_oxygen

		# Choose reference temperature
		Tref = 19
		ret@other_params$Tref = Tref

		# Create a temperature vector
		benthic_temp = vector(mode = "numeric", length = length(times))
		pelagic_temp = vector(mode = "numeric", length = length(times))

		# Constant temperature scenario
		benthic_temp[1:length(times)] = 19
		pelagic_temp[1:length(times)] = 19

		# Store these in other_params
		ret@other_params$benthic_temp = benthic_temp
		ret@other_params$pelagic_temp = pelagic_temp

		# Scale rates by oxygen and temperature: see otmscale.R
		ret = rate_scale(ret, t)
		ret = mscale(ret, t)
	
		# Run to steady state
		sim = project(ret, t_max = t, effort = effort)
	}
		
	# Growth matrices
	growth_mod = array(dim = c(length(yr), nrow(params@species_params), 50))
	species_all = c("cod","flounder","sprat","herring")
	dimnames(growth_mod) = list(Year = yr, Species = species_all, Age = seq(0,1,length.out=50))
	sp_maxage = c(15,26,16,13)
	growth_obs = growth_mod
	
	# Growth data
	sp_obs = list(cod = cod_ivb_cal, flounder = flounder_ivb_cal, sprat = sprat_ivb_cal, herring = herring_ivb_cal)
	
	for(i in 1:dim(growth_mod)[1]) {
		
		# Model growth
	
		# von Bertalanffy growth
		sp = species_params(ret)
		L_inf = (sp$w_inf/sp$a)^(1/sp$b)
		for(j in 1:dim(growth_obs)[2]) {
			growth_mod[i,j,] = myGrowthCurves(sim, (t-length(yr))+i, species = rownames(sp)[j], max_age = sp_maxage[j], percentage = T)
			growth_obs[i,j,] = sp_obs[[j]][i,]$a * (sp_obs[[j]][i,]$L_inf * (1 - exp(-sp_obs[[j]][i,]$k * (as.numeric(dimnames(growth_mod)$Age)*sp_maxage[j])))) ^ sp_obs[[j]][i,]$b
			growth_obs[i,j,] = (growth_obs[i,j,]/sp$w_inf[j])*100
		}
	}
		
	error = sum((growth_mod - growth_obs)^2, na.rm = T)
	return(error)
}