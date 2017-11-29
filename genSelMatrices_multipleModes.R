# Functions for generating variance/covariance matrices (F^(S)) with
#	multiple modes of convergent adaptation with examples at bottom
#
# Args:
#	rec: per base pair recombination rate estimate for the region
#	Ne: effective population size estimate
#	numPops: number of populations sampled (both selected and non-selected)
#	sampleSizes: vector of sample sizes of length numPops (i.e. twice the number
#		of diploid individuals sampled in each population)
#
#	F_estimate: estimate of neutral variance/covariance matrix
#		generated with "calcNeutralF.R"
#
#	positions: vector of genomic positions for region
#
#	sets: list of length number of different modes of convergence to be specified
#		vector "modes" where each element in list contains vector of populations with
#   a given single mode of convergence
#     i.e. if populations 2 and 6 share a mode and populations 3 has another,
#     sets = list(c(2,6), 3)
#
# numBins: the number of bins in which to bin alleles a given distance from the proposed
#		selected sites
#
#	*Specify parameter spaces for likelihood calculations*
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
#		*Note: the source must be a selected population in "set"


selPops = unlist(sets)
nonSelPops = seq(1, numPops)[- selPops]
distances = sapply(1:length(selSite), function(i) abs(positions - selSite[i]))

sampleErrorMatrix = diag(1/sampleSizes, nrow = numPops, ncol = numPops)

##bin distances
numBins = 1000
my.seq = seq(min(distances) - 0.001, max(distances) + 0.001, length.out = (numBins + 1))
midDistances = sapply(1:numBins, function(i) mean(c(my.seq[i], my.seq[i+1])))

##MVN parameters
k = numPops - 1
mu = as.matrix(rep(0, k))
rank = numPops - 1

##mean centering
M = numPops
Tmatrix = matrix(data = rep(-1 / M, (M - 1) * M), nrow = M - 1, ncol = M)
diag(Tmatrix) = (M - 1) / M 

calctotAddF_ind.set = function(y, set){
	# Computes additional variance/covariance induced
	#	by convergent adaptation due to independent
	# 	mutations in a set of populations
	#
	# Args:
	#	y: the probability of recombining off the 
	# 		beneficial background of the sweep (a
	#		function of recombination distance)
	#	set: the set of selected populations to
	#		apply this operation to
	#
	# Returns:
	#	additional variance/covariance to be added as numPops x numPops matrix

	selMatrix = F_estimate
	temp.selPops = sets[[set]]
	for(i  in temp.selPops) {
		selMatrix[i,i] = F_estimate[i,i] + y^2 * (1 - F_estimate[i,i])
	}
	return(selMatrix - F_estimate)
}

calctotAddF_sv_source.set = function(y, Rf, rt, p_no, p_one, my.source, set) {
	# Computes additional variance/covariance induced
	#	by convergent adaptation due to selection on
	# 	shared ancestral standing variation in a set
	#	of populations *with a source of the standing 
	#	variant specified*
	#
	# Args:
	#	y: the probability of recombining off the 
	# 		beneficial background of the sweep (a
	#		function of recombination distance)
	#	Rf: a function of the recombination distance
	#		that is used to specify whether two neutral lineages 
	# 		coalesce or recombine off the standing phase
	#	rt: the probability that a single lineage does not
	#		recombine off onto the non-beneficial background
	#		during the standing phase for t generations (a
	#		function of recombination distance)
	#	p_no: the probability no lineages recombine off or coalesce
	#		during time t (a function of recombination distance)
	#	p_one: the probability one lineage recombines off
	#		during time t (a function of recombination distance)
	#	my.source: the proposed source of the beneficial allele
	#		*Note: the source must be a selected population in "set"
	#	set: the set of selected populations to
	#		apply this operation to
	#
	# Returns:
	#	additional variance/covariance to be added as numPops x numPops matrix
	selMatrix = F_estimate
	temp.selPops = sets[[set]]
	temp.nonSelPops = seq(1, numPops)[-temp.selPops]
	
	selMatrix[my.source, my.source] = (1 - y^2) * (F_estimate[my.source, my.source]) + y^2 * (1 / (1 + Rf) + 
		Rf / (1 + Rf) * (F_estimate[my.source, my.source]))
	
	for(i in temp.selPops[temp.selPops != my.source]) {
		selMatrix[i,i] = (1 - y)^2 * (F_estimate[i, i]) + y^2 * (p_no*(1 / (1+Rf) + Rf / (1 + Rf) *
			F_estimate[my.source, my.source]) + (1 - p_no) * (1 / (1 + Rf)) + ((1 - p_no) * Rf / (1 + Rf) -
			Rf / (1 + Rf / 2) * (1 - p_one) * rt) * F_estimate[i,i] + (1 - p_one) * Rf / (1 + Rf / 2) *
			rt * F_estimate[i, my.source]) + 2 * y * (1 - y) * (rt * F_estimate[i, my.source] + (1 - rt) *
			F_estimate[i, i])
			
		selMatrix[i, my.source] = (1 - y)^2 * F_estimate[i, my.source] + y^2 * (rt^2 * (1 / (1 + Rf) + 
			Rf / (1 + Rf) * F_estimate[my.source, my.source]) + rt * (1 - rt) * 
			(F_estimate[my.source, my.source]) + rt * (1 - rt) * F_estimate[i, my.source] + (1 - rt)^2 *
			F_estimate[i, my.source]) + y * (1 - y) * (rt * F_estimate[my.source, my.source] + (1 - rt) *
			F_estimate[i, my.source]) + y * (1 - y) * F_estimate[i, my.source]
			
		selMatrix[my.source, i] = selMatrix[i, my.source]			
		
		for(k in temp.nonSelPops) {
			selMatrix[k, i] = y * rt * F_estimate[k, my.source] + (1 - y) * F_estimate[i, k] + 
				y * (1 - rt) * F_estimate[i, k]
			selMatrix[i, k] = y * rt * F_estimate[k, my.source] + (1 - y) * F_estimate[i, k] + 
				y * (1 - rt) * F_estimate[i, k]
		}
			

		for(j  in temp.selPops[temp.selPops != my.source]) {
			if(i != j)
				selMatrix[i,j] = (1 - y)^2 * F_estimate[i, j] + y^2 * ((rt^2 * (1 / (1 + Rf) + 
					Rf / (1 + Rf) * F_estimate[my.source, my.source]) + (1 - rt)^2 * F_estimate[i, j]) + 
					rt * (1 - rt) * (F_estimate[i, my.source] + F_estimate[j, my.source])) + (1 - y) * y *
					(2 * (1 - rt) * F_estimate[i, j] + rt * (F_estimate[i, my.source] + 
					F_estimate[j, my.source]))
		}

	}
	return(selMatrix - F_estimate)
}

calctotAddF_mig.set = function(y, e_delta, my.Q, my.source, set){
	# Computes additional variance/covariance induced
	#	by convergent adaptation due to migration
	# 	in a set of populations
	#
	# Args:
	#	y: the probability of recombining off the 
	# 		beneficial background of the sweep (a
	#		function of recombination distance)
	#	e_delta: the probability of recombining off the 
	# 		beneficial background of the sweep in the
	#		source population for time delta (a
	#		function of recombination distance)
	#	my.Q: the probability of coalescing before recombination
	#		at the selected site
	#	my.source: the proposed migration source of the
	#		beneficial allele
	#		*Note: the source must be a selected population in "set"
	#	set: the set of selected populations to
	#		apply this operation to
	#
	# Returns:
	#	additional variance/covariance to be added as numPops x numPops matrix
	selMatrix = F_estimate
	temp.selPops = sets[[set]]
	temp.nonSelPops = seq(1, numPops)[-temp.selPops]
		
	if(is.element(my.source, temp.selPops)) {
		selMatrix[my.source, my.source] = (F_estimate[my.source, my.source]) +
			y^2 * (1 - (F_estimate[my.source, my.source]))
				
		for(i in temp.selPops[temp.selPops != my.source]) {
		  selMatrix[i, i] = my.Q * (y^2 + (1 - y^2) * (F_estimate[i, i]) + 2 * y * (1 - y) * (F_estimate[i, my.source])) +
		    (1 - my.Q) * (y^2 * e_delta^2 + (1 - y)^2 * (F_estimate[i, i]) + 2 * y * (1 - y) * 
		    F_estimate[my.source, i] + y^2 * (1 - e_delta^2) * (F_estimate[my.source, my.source]))
			selMatrix[i, my.source] = y^2 * e_delta + (1 - y) * F_estimate[my.source, i] + 
				y * (1 - y * e_delta) * (F_estimate[my.source, my.source])
			selMatrix[my.source, i] = y^2 * e_delta + (1 - y) * F_estimate[my.source, i] +
				y * (1 - y * e_delta) * (F_estimate[my.source, my.source])
			
			for(k in temp.nonSelPops) {
				selMatrix[k, i] = (1 - y) * F_estimate[i, k] + y * F_estimate[my.source, k]
				selMatrix[i, k] = (1 - y) * F_estimate[i, k] + y * F_estimate[my.source, k]
			}
			
			for(j  in temp.selPops[temp.selPops != my.source]) {
				if(i != j)
				selMatrix[i,j] = y^2 * e_delta^2 + y^2 * (1 - e_delta^2) * (F_estimate[my.source, my.source]) +
					(1 - y)^2 * F_estimate[i,j]	+ (1 - y) * y * (F_estimate[i, my.source] +
					F_estimate[j, my.source])
			}
		}
	}
	return(selMatrix - F_estimate)
}

calcFOmegas_mixed = function(sel, g, time, mig, my.source, modes) {
	# Computes mean-centered and sample-size corrected 
	#	variance/covariance matrices with convergent
	#	selection (F^(S)) for all bins for a model with multiple
	#	modes of convergence (independent mutations, standing variant 
	#	with source, and/or migration) for a given set of parameters
	#
	# Args:
	#	sel: strength of selection
	#	g: frequency of the standing variant
	#	t: time in generations the variant is standing in populations
	#		before selection occurs and prior to migration from source population
	#	mig: migration rate (proportion of individuals from source each generation)
	#	my.source: source population of the beneficial allele for both migration
	#		and standing variant with source models
	#		*Note: the source must be a selected population in "set"
	#	modes: character vector of length sets defining mode for each set
	#		of selected populations ("ind", "sv", and/or "mig")
	#
	# Returns:
	#	list of variance/covariance matrices with convergent
	#		selection (F^(S)) of length numBins
	y = exp(-rec * midDistances / sel * log(4 * Ne * sel))

	FOmegas = vector("list", length(midDistances))
	FOmegas = lapply(FOmegas, function(i) matrix(0, nrow = numPops, ncol = numPops))
	
	
	for(set in 1 : length(sets)) {
		if(modes[set] == "sv") {
			y = exp(-rec * midDistances / sel * log(1 / g))
			Rf = 4 * Ne * rec * midDistances * g
			rt = exp(-rec * midDistances * time)
			p_no = exp(-time * (2 * rec * midDistances + 1 / (2 * Ne * g)))
			p_one = exp(-time * (rec * midDistances + 1 / (2 * Ne * g)))
			FOmegas = lapply(1 : length(midDistances), function(i) calctotAddF_sv_source.set(y[[i]], Rf[[i]],
				rt[[i]], p_no[[i]], p_one[[i]], my.source, set) + FOmegas[[i]])	
		}
		else if(modes[set] == "ind") {
			FOmegas = lapply(1 : length(midDistances), function(i) calctotAddF_ind.set(y[[i]], set) +
				FOmegas[[i]])
		}
		else if(modes[set] == "mig") {
			delta = 1 / sel * log(1 + sel/(mig))
			e_delta = exp(-rec * midDistances * delta)
			my.Q = 1 / (1 + 4 * Ne * mig)
			FOmegas = lapply(1 : length(midDistances), function(i) calctotAddF_mig.set(y[[i]], e_delta[[i]],
				my.Q, my.source, set) + FOmegas[[i]])
		}
	}
	
	FOmegas = lapply(FOmegas, function(i) Tmatrix %*% (i + F_estimate + sampleErrorMatrix) %*% t(Tmatrix))
	return(FOmegas)
}