# Functions for calculating composite log-likelihoods under all models 
#	of convergent selection: neutral model, independent mutations,
#	standing variant model with source, standing variant model without source,
#	migration model, combinations of convergent selection models (mixed)
#
# Args:
#	positions: vector of genomic positions for region


##neutral model
calcLikelihood_bin_neutral = function(site, det_FOmegas, inv_FOmegas, selSiteLoc) {
	# Calculates log-likelihood of data at a given position for 
	#	neutral model
	#
	# Args:
	#	site: element of vector "positions" of position for log-likelihood to be calculated  
	#	det_FOmegas: list of length numBins of determinants of variance/covariance matrices for
	#		given model
	#	inv_FOmegas: list of length numBins of inverses of variance/covariance matrices for
	#		given model
	#	selSiteLoc: element of vector "selSite" of proposed location of selected site
	#
	# Returns:
	#	log-likelihood of data at a given position
	bin = distBins[site, j]
	my.x = as.matrix(freqs_MC[ , site])
	my.e = epsilons[site]*(1 - epsilons[site])

	likelihood = 1 / (sqrt((2 * pi)^k * (det_FOmegas * my.e^rank))) * exp(-1 / 2 * t(my.x - mu) %*%
		(inv_FOmegas / my.e) %*% (my.x - mu))
	return(log(likelihood))
}

calcCompLikelihood_neutral = function(selSiteLoc, det_FOmegas, inv_FOmegas) {
	# Calculates composite log-likelihood of all data for neutral model
	#
	# Args: 
	#	selSiteLoc: element of vector "selSite" of proposed location of selected site
	#	det_FOmegas: list of length numBins of determinants of variance/covariance matrices for
	#		given model
	#	inv_FOmegas: list of length numBins of inverses of variance/covariance matrices for
	#		given model
	#
	# Returns:
	#	composite log-likelihood of data under neutral model
	all = sapply(1 : length(positions), {
			function(i) calcLikelihood_bin_neutral(i, det_FOmegas, inv_FOmegas, j)
		})
	return(sum(all))
}

calcLikelihood_bin.1par = function(site, selSiteLoc, det_FOmegas, inv_FOmegas, par1) {
	# Calculates log-likelihood of data at a given position for models
	#	with a single parameter (independent mutations model) 
	#
	# Args:
	#	site: element of vector "positions" of position for log-likelihood to be calculated  
	#	selSiteLoc: element of vector "selSite" of proposed location of selected site
	#	det_FOmegas: list of length numBins of determinants of variance/covariance matrices for
	#		given model of convergent adaptation
	#	inv_FOmegas: list of length numBins of inverses of variance/covariance matrices for
	#		given model of convergent adaptation
	#	par1: element of vector of parameter used to specify given model of convergent adaptation
	#		("sels" for independent mutations model)
	#
	# Returns:
	#	log-likelihood of data at a given position
	bin = distBins[site, selSiteLoc]
	my.x = as.matrix(freqs_MC[ , site])
	my.e = epsilons[site] * (1 - epsilons[site])
	likelihood = 1 / (sqrt((2 * pi)^k * (det_FOmegas[[par1]][[bin]] * my.e^rank))) * exp(-1 / 2 * t(my.x -
		mu) %*% (inv_FOmegas[[par1]][[bin]] / my.e) %*% (my.x - mu))
	return(log(likelihood))
}

calcCompLikelihood.1par = function(selSiteLoc, det_FOmegas, inv_FOmegas, par1) {
	# Calculates composite log-likelihood of all data for models with a 
	#	single parameter (independent mutations model) 
	#
	# Args: 
	#	selSiteLoc: element of vector "selSite" of proposed location of selected site
	#	det_FOmegas: list of length numBins of determinants of variance/covariance matrices for
	#		given model of convergent adaptation
	#	inv_FOmegas: list of length numBins of inverses of variance/covariance matrices for
	#		given model of convergent adaptation
	#	par1: element of vector of parameter used to specify given model of convergent adaptation
	#		("sels" for independent mutations model)
	#
	# Returns:
	#	composite log-likelihood of data under model
	all = sapply(1 : length(positions), {
			function(i) calcLikelihood_bin.1par(i, selSiteLoc, det_FOmegas, inv_FOmegas, par1)
	})
	return(sum(all))
}

calcLikelihood_bin.3par = function(site, selSiteLoc, det_FOmegas, inv_FOmegas, par1, par2, par3) {
	# Calculates log-likelihood of data at a given position for models
	#	with three parameters (standing model without source, migration model) 
	#
	# Args:
	#	site: element of vector "positions" of position for log-likelihood to be calculated  
	#	selSiteLoc: element of vector "selSite" of proposed location of selected site
	#	det_FOmegas: list of length numBins of determinants of variance/covariance matrices for
	#		given model of convergent adaptation
	#	inv_FOmegas: list of length numBins of inverses of variance/covariance matrices for
	#		given model of convergent adaptation
	#	par1: element of vector of parameter used to specify given model of convergent adaptation
	#		("sels" for standing model without source and migration model)
	#	par2: element of vector of parameter used to specify given model of convergent adaptation
	#		("gs" for standing model without source, "migs" for migration model)
	#	par3: element of vector of parameter used to specify given model of convergent adaptation
	#		("times" for standing model without source, "source" for migration model)
	#
	# Returns:
	#	log-likelihood of data at a given position
	bin = distBins[site, selSiteLoc]
	my.x = as.matrix(freqs_MC[ , site])
	my.e = epsilons[site] * (1-epsilons[site])

	likelihood = 1 / (sqrt((2 * pi)^k * (det_FOmegas[[par1]][[par2]][[par3]][[bin]] * my.e^rank))) * exp(-1 / 2
		* t(my.x - mu) %*% (inv_FOmegas[[par1]][[par2]][[par3]][[bin]] / my.e) %*% (my.x - mu))
	
	return(log(likelihood))
}

calcCompLikelihood_bin.3par = function(selSiteLoc, det_FOmegas, inv_FOmegas, par1, par2, par3) {
	# Calculates composite log-likelihood of all data for models with a 
	#	three parameters (standing model without source, migration model) 
	#
	# Args: 
	#	selSiteLoc: element of vector "selSite" of proposed location of selected site
	#	det_FOmegas: list of length numBins of determinants of variance/covariance matrices for
	#		given model of convergent adaptation
	#	inv_FOmegas: list of length numBins of inverses of variance/covariance matrices for
	#		given model of convergent adaptation
	#	par1: element of vector of parameter used to specify given model of convergent adaptation
	#		("sels" for standing model without source and migration model)
	#	par2: element of vector of parameter used to specify given model of convergent adaptation
	#		("gs" for standing model without source, "migs" for migration model)
	#	par3: element of vector of parameter used to specify given model of convergent adaptation
	#		("times" for standing model without source, "source" for migration model)
	#
	# Returns:
	#	composite log-likelihood of data under model
	all = sapply(1 : nrow(distances), {
			function(i) calcLikelihood_bin.3par(i, selSiteLoc, det_FOmegas, inv_FOmegas, par1, par2, par3)
	})
	return(sum(all))
}

calcLikelihood_bin.4par = function(site, selSiteLoc, det_FOmegas, inv_FOmegas, par1, par2, par3, par4) {
	# Calculates log-likelihood of data at a given position for models
	#	with four parameters (standing model with source) 
	#
	# Args:
	#	site: element of vector "positions" of position for log-likelihood to be calculated  
	#	selSiteLoc: element of vector "selSite" of proposed location of selected site
	#	det_FOmegas: list of length numBins of determinants of variance/covariance matrices for
	#		given model of convergent adaptation
	#	inv_FOmegas: list of length numBins of inverses of variance/covariance matrices for
	#		given model of convergent adaptation
	#	par1: element of vector of parameter used to specify given model of convergent adaptation
	#		("sels" for standing model with source)
	#	par2: element of vector of parameter used to specify given model of convergent adaptation
	#		("gs" for standing model with source)
	#	par3: element of vector of parameter used to specify given model of convergent adaptation
	#		("times" for standing model with source)
	#	par4: element of vector of parameter used to specify given model of convergent adaptation
	#		("sources" for standing model with source)
	#
	# Returns:
	#	log-likelihood of data at a given position
	bin = distBins[site, selSiteLoc]
	my.x = as.matrix(freqs_MC[ , site])
	my.e = epsilons[site]*(1-epsilons[site])

	likelihood = 1/(sqrt((2 * pi)^k * (det_FOmegas[[par1]][[par2]][[par3]][[par4]][[bin]] * my.e^rank))) *
		exp(-1 / 2 * t(my.x-mu) %*% (inv_FOmegas[[par1]][[par2]][[par3]][[par4]][[bin]] / my.e) %*% (my.x - mu))
	
	return(log(likelihood))
}

calcCompLikelihood_bin.4par = function(selSiteLoc, det_FOmegas, inv_FOmegas, par1, par2, par3, par4) {
	# Calculates composite log-likelihood of all data for models with a 
	#	three parameters (standing model with source, migration model) 
	#
	# Args: 
	#	selSiteLoc: element of vector "selSite" of proposed location of selected site
	#	det_FOmegas: list of length numBins of determinants of variance/covariance matrices for
	#		given model of convergent adaptation
	#	inv_FOmegas: list of length numBins of inverses of variance/covariance matrices for
	#		given model of convergent adaptation
	#	par1: element of vector of parameter used to specify given model of convergent adaptation
	#		("sels" for standing model with source)
	#	par2: element of vector of parameter used to specify given model of convergent adaptation
	#		("gs" for standing model with source)
	#	par3: element of vector of parameter used to specify given model of convergent adaptation
	#		("times" for standing model with source)
	#	par4: element of vector of parameter used to specify given model of convergent adaptation
	#		("sources" for standing model with source)
	#
	# Returns:
	#	composite log-likelihood of data under model
	all = sapply(1 : length(positions), {
			function(i) calcLikelihood_bin.4par(i, selSiteLoc, det_FOmegas, inv_FOmegas, par1, par2, par3, par4)
		})
	return(sum(all))
}

calcLikelihood_bin.5par = function(site, selSiteLoc, det_FOmegas, inv_FOmegas, par1, par2, par3, par4, par5) {
	# Calculates log-likelihood of data at a given position for models
	#	with four parameters (model with multiple modes) 
	#
	# Args:
	#	site: element of vector "positions" of position for log-likelihood to be calculated  
	#	selSiteLoc: element of vector "selSite" of proposed location of selected site
	#	det_FOmegas: list of length numBins of determinants of variance/covariance matrices for
	#		given model of convergent adaptation
	#	inv_FOmegas: list of length numBins of inverses of variance/covariance matrices for
	#		given model of convergent adaptation
	#	par1: element of vector of parameter used to specify given model of convergent adaptation
	#		("sels" for model with multiple modes)
	#	par2: element of vector of parameter used to specify given model of convergent adaptation
	#		("gs" for model with multiple modes)
	#	par3: element of vector of parameter used to specify given model of convergent adaptation
	#		("times" for model with multiple modes)
	#	par4: element of vector of parameter used to specify given model of convergent adaptation
	#		("migs" for model with multiple modes)
	#	par5: element of vector of parameter used to specify given model of convergent adaptation
	#		("sources" for model with multiple modes)
	#
	# Returns:
	#	log-likelihood of data at a given position
	bin = distBins[site, selSiteLoc]
	my.x = as.matrix(freqs_MC[ ,site])
	my.e = epsilons[site] * (1-epsilons[site])

	likelihood = 1/(sqrt((2 * pi)^k * (det_FOmegas[[par1]][[par2]][[par3]][[par4]][[par5]][[bin]] * my.e^rank)))
		* exp(-1 / 2 * t(my.x - mu) %*% (inv_FOmegas[[par1]][[par2]][[par3]][[par4]][[par5]][[bin]] / my.e) %*%
		(my.x - mu))
	
	return(log(likelihood))
}

calcCompLikelihood_bin.5par = function(selSiteLoc, det_FOmegas, inv_FOmegas, par1, par2, par3, par4, par5) {
	# Calculates composite log-likelihood of all data for models with a 
	#	three parameters (model with multiple modes) 
	#
	# Args: 
	#	selSiteLoc: element of vector "selSite" of proposed location of selected site
	#	det_FOmegas: list of length numBins of determinants of variance/covariance matrices for
	#		given model of convergent adaptation
	#	inv_FOmegas: list of length numBins of inverses of variance/covariance matrices for
	#		given model of convergent adaptation
	#	par1: element of vector of parameter used to specify given model of convergent adaptation
	#		("sels" for model with multiple modes)
	#	par2: element of vector of parameter used to specify given model of convergent adaptation
	#		("gs" for model with multiple modes)
	#	par3: element of vector of parameter used to specify given model of convergent adaptation
	#		("times" for model with multiple modes)
	#	par4: element of vector of parameter used to specify given model of convergent adaptation
	#		("migs" for model with multiple modes)
	#	par5: element of vector of parameter used to specify given model of convergent adaptation
	#		("sources" for model with multiple modes)
	#
	# Returns:
	#	composite log-likelihood of data under model
	all = sapply(1 : length(positions), {
			function(i) calcLikelihood_bin.5par(i, selSiteLoc, det_FOmegas, inv_FOmegas, par1, par2, par3, par4, par5)
	})
	return(sum(all))
}