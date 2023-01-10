# Load oxygen data
load("Data/oxy_model.rda")

# A function to find ideal vs poor oxygen and run respective simulations
#	params: a Mizer params object
#	t: time steps to run model
#	b: is there benthic resource carrying capacity scaling
#	o: is there occupancy scaling
#	p: is there physiological scaling
#	returns a list with ideal (yay) and poor (death) simulations
oxy_scenario = function(params, t, p = F) {
	
	# Create params objects
	ret_yay = params
	ret_okay = params
	ret_death = params
	
	# Create an oxygen vector
	times = 0:t
	bo_yay = vector(mode = "numeric", length = length(times))
	po_yay = vector(mode = "numeric", length = length(times))
	bo_okay = bo_yay
	po_okay = po_yay
	bo_death = bo_yay
	po_death = po_yay
	
	# Extend temperature vectors
	ret_yay@other_params$benthic_temp = rep(ret_yay@other_params$benthic_temp[1], t+1)
	ret_yay@other_params$pelagic_temp = rep(ret_yay@other_params$pelagic_temp[1], t+1)
	ret_okay@other_params$benthic_temp = rep(ret_yay@other_params$benthic_temp[1], t+1)
	ret_okay@other_params$pelagic_temp = rep(ret_yay@other_params$pelagic_temp[1], t+1)
	ret_death@other_params$benthic_temp = rep(ret_yay@other_params$benthic_temp[1], t+1)
	ret_death@other_params$pelagic_temp = rep(ret_yay@other_params$pelagic_temp[1], t+1)
	
	# Use 3 mL/L as ideal oxygen
	oxy_yay = 3
	
	# Use 2 mL/L as okay oxygen
	oxy_okay = 2
	
	# Use 1 mL/L as low oxygen
	oxy_death = 1

	# Use 6 mL/l as pelagic oxygen
	oxy_pelagic = 6
	
	# Store oxygen values in a vector
	oxy = c(yay = oxy_yay, okay = oxy_okay, death = oxy_death)
	
	# Fill oxygen vectors
	bo_yay[1:length(times)] = oxy_yay
	po_yay[1:length(times)] = oxy_pelagic
	bo_okay[1:length(times)] = oxy_okay
	po_okay[1:length(times)] = oxy_pelagic
	bo_death[1:length(times)] = oxy_death
	po_death[1:length(times)] = oxy_pelagic
	
	# Store these in other params
	ret_yay@other_params$benthic_oxygen = bo_yay
	ret_yay@other_params$pelagic_oxygen = po_yay
	ret_okay@other_params$benthic_oxygen = bo_okay
	ret_okay@other_params$pelagic_oxygen = po_okay
	ret_death@other_params$benthic_oxygen = bo_death
	ret_death@other_params$pelagic_oxygen = po_death
	
	# Shorten temperature object if necessary
	ret_yay@other_params$benthic_temp = ret_yay@other_params$benthic_temp[1:length(bo_yay)]
	ret_yay@other_params$pelagic_temp = ret_yay@other_params$pelagic_temp[1:length(po_yay)]
	ret_okay@other_params$benthic_temp = ret_okay@other_params$benthic_temp[1:length(bo_okay)]
	ret_okay@other_params$pelagic_temp = ret_okay@other_params$pelagic_temp[1:length(po_okay)]
	ret_death@other_params$benthic_temp = ret_death@other_params$benthic_temp[1:length(bo_death)]
	ret_death@other_params$pelagic_temp = ret_death@other_params$pelagic_temp[1:length(po_death)]
	
	# Scale occupancy
	ret_yay = occupancy(ret_yay, t)
	ret_okay = occupancy(ret_okay, t)
	ret_death = occupancy(ret_death, t)
	
	# Scale objects
	if(p) {
		# Scale physiological rates
		ret_yay = rate_scale(ret_yay, t)
		ret_yay = mscale(ret_yay, t)
		ret_okay = rate_scale(ret_okay, t)
		ret_okay = mscale(ret_okay, t)
		ret_death = rate_scale(ret_death, t)
		ret_death = mscale(ret_death, t)
	}
	if(!p) {
		# Scale physiological rates
		ret_yay = rate_scale(ret_yay, t)
		ret_yay = mscale(ret_yay, t)
		ret_okay = rate_scale(ret_okay, t)
		ret_okay = mscale(ret_okay, t)
		ret_death = rate_scale(ret_death, t)
		ret_death = mscale(ret_death, t)
		
		# Make sure scaling is set to zero
		ret_yay@other_params$rate_scale[,,] = 1
		ret_yay@other_params$metab_scale[,,] = 1
		ret_okay@other_params$rate_scale[,,] = 1
		ret_okay@other_params$metab_scale[,,] = 1
		ret_death@other_params$rate_scale[,,] = 1
		ret_death@other_params$metab_scale[,,] = 1
	}
	
	# Run scenarios
	sim_yay = project(ret_yay, t_max = t, effort = c(1,1,1,1))
	sim_okay = project(ret_okay, t_max = t, effort = c(1,1,1,1))
	sim_death = project(ret_death, t_max = t, effort = c(1,1,1,1))
	
	# Return scenarios
	return(list(yay = sim_yay, okay = sim_okay, death = sim_death, oxy = oxy))
}

# A function to get SSB
#	sim: a Mizer sim object
#	returns ssb for each species at the last time step of the model
get_ssb = function(sim) {

	# Set up SSB output matrix
	ssb = matrix(NA, nrow = 1, ncol = 4)
	colnames(ssb) = c("cod", "flounder", "sprat", "herring")
	n_all = sim@n[dim(sim@n)[1],,]
	
	# Get SSB for each species
	for(i in 1:ncol(ssb)) {
		
		# Use weights greater than w_mat
		idx_w = sim@params@w >= species_params(sim@params)$w_mat[i]
		
		# Model
		n = n_all[i,idx_w]
		mod = n * sim@params@w[idx_w] * sim@params@dw[idx_w]
		ssb[,i] = sum(mod)
	}
	
	return(ssb)
}

# A function to plot a gradient of colors on the food web diagrams
#	n: number of rectangles
#	pal: palette
#	returns nothing, plots a vertical gradient on the food web plots for ***aesthetics***
waterplot = function(n, pal) {
	
	# Divide the vertical areas of the plot
	y = seq(0, 1, length.out = n+1)
	
	# Get colors
	colors = col2rgb(hcl.colors(n, pal, rev = F))
	
	# Plot the polygons
	for(i in 1:n) {
		rd = colors[1,i]
		gr = colors[2,i]
		bl = colors[3,i]
		polygon(c(0,0,1,1), c(y[i],y[i+1],y[i+1],y[i]), col = rgb(rd,gr,bl,alpha=0.25*255,maxColorValue=255), border = NA)
	}
}