# Functions for calculating composite log-likelihoods under all models 
#	of convergent selection. Function used depends on how many parameters
# there are in the model you are calculating composite likelihood for

##neutral model
calcLikelihood_bin_neutral = function(site, selSiteLoc, det_FOmegas, inv_FOmegas) {
  # Calculates log-likelihood of data at a given position for 
  #	neutral model
  #
  # Args:
  #	site: element of vector "positions" of position for log-likelihood to be calculated 
  #	selSiteLoc: element of vector "selSite" of proposed location of selected site
  #	det_FOmegas: list of length numBins of determinants of variance/covariance matrices for
  #		given model
  #	inv_FOmegas: list of length numBins of inverses of variance/covariance matrices for
  #		given model
  #	selSiteLoc: element of vector "selSite" of proposed location of selected site
  #
  # Returns:
  #	log-likelihood of data at a given position
  bin = distBins[site,selSiteLoc]
  my.x = as.matrix(freqs_MC[, site ])
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
    function(i) calcLikelihood_bin_neutral(i, selSiteLoc, det_FOmegas, inv_FOmegas)
  })
  return(sum(all))
}

saveRDS(calcCompLikelihood_neutral, file="calcCompLikelihood_neutral.RDS")
saveRDS(calcLikelihood_bin_neutral, file="calcLikelihood_bin_neutral.RDS")

##Any model with one parameter
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
  #
  # Returns:
  #	likelihood of data at a given position
  bin = distBins[site, selSiteLoc]
  my.x = as.matrix(freqs_MC[, site ])
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
  #
  # Returns:
  #	composite log-likelihood of data under model
  all = sapply(1 : length(positions), {
    function(i) calcLikelihood_bin.1par(i, selSiteLoc, det_FOmegas, inv_FOmegas, par1)
  })
  return(sum(all))
}

saveRDS(calcLikelihood_bin.1par, file="calcLikelihood_bin_1par.RDS")
saveRDS(calcCompLikelihood.1par, file="calcCompLikelihood_1par.RDS")

##Any model with two parameters
calcLikelihood_bin.2par = function(site, selSiteLoc, det_FOmegas, inv_FOmegas, par1, par2) {
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
  #	par2: element of vector of parameter used to specify given model of convergent adaptation
  #
  # Returns:
  #	likelihood of data at a given position
  bin = distBins[site, selSiteLoc]
  my.x = as.matrix(freqs_MC[ , site])
  my.e = epsilons[site] * (1-epsilons[site])
  
  likelihood = 1 / (sqrt((2 * pi)^k * (det_FOmegas[[par1]][[par2]][[bin]] * my.e^rank))) * exp(-1 / 2
                                                                                               * t(my.x - mu) %*% (inv_FOmegas[[par1]][[par2]][[bin]] / my.e) %*% (my.x - mu))
  
  return(log(likelihood))
}

calcCompLikelihood.2par = function(selSiteLoc, det_FOmegas, inv_FOmegas, par1, par2) {
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
  #	par2: element of vector of parameter used to specify given model of convergent adaptation
  #
  # Returns:
  #	composite log-likelihood of data under model
  all = sapply(1 : length(positions), {
    function(i) calcLikelihood_bin.2par(i, selSiteLoc, det_FOmegas, inv_FOmegas, par1, par2)
  })
  return(sum(all))
}

saveRDS(calcLikelihood_bin.2par, file="calcLikelihood_bin_2par.RDS")
saveRDS(calcCompLikelihood.2par ,file="calcCompLikelihood_2par.RDS")


##Any model with three parameters
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
  #	par2: element of vector of parameter used to specify given model of convergent adaptation
  #	par3: element of vector of parameter used to specify given model of convergent adaptation
  #
  # Returns:
  #	likelihood of data at a given position
  bin = distBins[site, selSiteLoc]
  my.x = as.matrix(freqs_MC[ , site])
  my.e = epsilons[site] * (1-epsilons[site])
  
  likelihood = 1 / (sqrt((2 * pi)^k * (det_FOmegas[[par1]][[par2]][[par3]][[bin]] * my.e^rank))) * exp(-1 / 2
                                                                                                       * t(my.x - mu) %*% (inv_FOmegas[[par1]][[par2]][[par3]][[bin]] / my.e) %*% (my.x - mu))
  
  return(log(likelihood))
}

calcCompLikelihood.3par = function(selSiteLoc, det_FOmegas, inv_FOmegas, par1, par2, par3) {
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
  #	par2: element of vector of parameter used to specify given model of convergent adaptation
  #	par3: element of vector of parameter used to specify given model of convergent adaptation
  #
  # Returns:
  #	composite log-likelihood of data under model
  all = sapply(1 : length(positions), {
    function(i) calcLikelihood_bin.3par(i, selSiteLoc, det_FOmegas, inv_FOmegas, par1, par2, par3)
  })
  return(sum(all))
}

saveRDS(calcLikelihood_bin.3par ,file="calcLikelihood_bin_3par.RDS")
saveRDS(calcCompLikelihood.3par,file="calcCompLikelihood_3par.RDS")

