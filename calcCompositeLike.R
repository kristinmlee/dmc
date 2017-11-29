# Functions for calculating composite log-likelihoods under all models 
#	of convergent selection: neutral model, independent mutations,
#	standing variant model with source, standing variant model without source,
#	migration model, combinations of convergent selection models (mixed)
#
# Args:
#	numPops: number of populations sampled (both selected and non-selected)
#
#	positions: vector of genomic positions for region
#
#freqs: matrix of allele frequencies for window of interest with	
#		dimension numberOfPopulations x numberOfSites
#
# numBins: the number of bins in which to bin alleles a given distance from the proposed
#		selected sites
#   NOTE: must be same as number specified in generating covariance matrices
#
#	*Specify parameter spaces for likelihood calculations*
# *NOTE: must be same as number specified in generating covariance matrices*
#	selSite: vector of positions of proposed selected sites
#	sels: vector of proposed selection coefficients
#	times: vector of proposed time in generations the variant is standing 
#		in populations before selection occurs and prior to migration from
#		source population
#	gs: vector of proposed frequencies of the standing variant
#	migs: migration rate (proportion of individuals from source each generation)
#		*Note: cannot be 0
#	sources: vector of proposed source population of the beneficial allele
#		for both migration and standing variant with source models
#		*Note: the source must be a selected population in selPops
#
#
# Inverses and determinants for all relevant models of convergent selection



#calculate distances from proposed selected sites and bin
distances = sapply(1:length(selSite), function(i) abs(positions - selSite[i]))
numBins = 1000
my.seq = seq(min(distances) - 0.001, max(distances) + 0.001, length.out = (numBins + 1))
distBins = apply(distances, 2, function(i) as.numeric(cut(i, my.seq)))

#mean centering
M = numPops
Tmatrix = matrix(data = rep(-1 / M, (M - 1) * M), nrow = M - 1, ncol = M)
diag(Tmatrix) = (M - 1) / M 

#get site-specific mean allele frequencies across populations and mean-centered population allele frequencies
freqs = t(freqs)
epsilons = rowMeans(freqs)
freqs_MC = sapply(1 : nrow(freqs), function(i) Tmatrix %*% freqs[i,])

#MVN parameters
k = numPops - 1
mu = as.matrix(rep(0, k))
rank = numPops - 1

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
	bin = distBins[site, selSiteLoc]
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
			function(i) calcLikelihood_bin_neutral(i, det_FOmegas, inv_FOmegas, selSiteLoc)
		})
	return(sum(all))
}

calcLikelihood_bin_1par = function(site, selSiteLoc, det_FOmegas, inv_FOmegas, par1) {
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

calcCompLikelihood_1par = function(selSiteLoc, det_FOmegas, inv_FOmegas, par1) {
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
			function(i) calcLikelihood_bin_1par(i, selSiteLoc, det_FOmegas, inv_FOmegas, par1)
	})
	return(sum(all))
}

calcLikelihood_bin_3par = function(site, selSiteLoc, det_FOmegas, inv_FOmegas, par1, par2, par3) {
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

calcCompLikelihood_3par = function(selSiteLoc, det_FOmegas, inv_FOmegas, par1, par2, par3) {
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
			function(i) calcLikelihood_bin_3par(i, selSiteLoc, det_FOmegas, inv_FOmegas, par1, par2, par3)
	})
	return(sum(all))
}

calcLikelihood_bin_4par = function(site, selSiteLoc, det_FOmegas, inv_FOmegas, par1, par2, par3, par4) {
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

calcCompLikelihood_4par = function(selSiteLoc, det_FOmegas, inv_FOmegas, par1, par2, par3, par4) {
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
			function(i) calcLikelihood_bin_4par(i, selSiteLoc, det_FOmegas, inv_FOmegas, par1, par2, par3, par4)
		})
	return(sum(all))
}

calcLikelihood_bin_5par = function(site, selSiteLoc, det_FOmegas, inv_FOmegas, par1, par2, par3, par4, par5) {
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

	likelihood = 1/(sqrt((2 * pi)^k * (det_FOmegas[[par1]][[par2]][[par3]][[par4]][[par5]][[bin]] * my.e^rank))) *
	  exp(-1 / 2 * t(my.x - mu) %*% (inv_FOmegas[[par1]][[par2]][[par3]][[par4]][[par5]][[bin]] / my.e) %*%
		(my.x - mu))
	
	return(log(likelihood))
}

calcCompLikelihood_5par = function(selSiteLoc, det_FOmegas, inv_FOmegas, par1, par2, par3, par4, par5) {
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
			function(i) calcLikelihood_bin_5par(i, selSiteLoc, det_FOmegas, inv_FOmegas, par1, par2, par3, par4, par5)
	})
	return(sum(all))
}