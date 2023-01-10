# Density independent reproduction scaled by hypoxia exposure
#	What is scaled?
#		1. Reproduction efficiency is scaled by physiological hypoxia response
hypoxiaRDI = function(params, n, n_pp, n_other, t, e_growth, mort, e_repro, ...) {

	# Calculate total energy from per capita energy
	e_repro_pop <- drop((e_repro * n * params@other_params$rate_scale[t+1,,]) %*% params@dw)

	# Assume sex_ratio = 0.5
	rdi <- 0.5 * (e_repro_pop * params@species_params$erepro) / params@w[params@w_min_idx]
	return(rdi)
}