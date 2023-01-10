# Universal feeding level scaled by hypoxia exposure
#	What is scaled?
#		1. Maximum intake level scaled by physiological hypoxia response
hypoxiaFeedingLevel <- function(params, n, n_pp, n_other, t, encounter, ...) {
	return(encounter / (encounter + (params@intake_max * params@other_params$rate_scale[t+1,,])))
}