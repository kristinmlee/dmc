# Functions for get maximum composite log-likelihood estimates
# for parameters under all models 



getMCLEind = function(compLike, selSite, sels) {
  # Gets MCLE for location of selected site and strength of selection (s)
  #	  for independent mutations model 
  #
  # Args: 
  #	compLike: nested lists of composite likelihoods (list of length length(selSite) where each
  #   element is a list of length length(sels) containing composite likelihood values)
  #	selSite: selected sites used to calculate compLike
  #	sels: selection coefficients used to calculate compLike
  #
  # Returns:
  #	MCLE for location of selected site and s
  maxLocIndex = ceiling(which.max(unlist(compLike)) / length(sels))
  maxSelIndex = which.max(unlist(compLike[[maxLocIndex]]))
                          
  maxLoc = selSite[maxLocIndex]
  maxSel = sels[maxSelIndex]
  return(cbind(maxLoc, maxSel))
}

getMCLEmig = function(compLike, selSite, sels, migs, sources) {
  # Gets MCLE for location of selected site strength of selection (s), migration rate (m),
  #   and source population for migration model
  #
  # Args: 
  #	compLike: nested lists of composite likelihoods (list of length length(selSite) where each
  #   element is a list of length length(sels) ... of length length(migs) ... of length length(sources)
  #   containing composite likelihood values)
  #	selSite: selected sites used to calculate compLike
  #	sels: selection coefficients used to calculate compLike
  # migs: migration rates used to calculate compLike
  # sources: source populations used to calculate compLike
  #
  # Returns:
  #	MCLE for location of selected site, s, m, and source population
  maxLocIndex = ceiling(which.max(unlist(compLike)) / (length(sels) * length(migs) * length(sources)))
  maxSelIndex = ceiling(which.max(unlist(compLike[[maxLocIndex]])) / (length(migs) * length(sources)))
  maxMigIndex = ceiling(which.max(unlist(compLike[[maxLocIndex]][[maxSelIndex]])) / length(sources))
  maxSourceIndex = which.max(unlist(compLike[[maxLocIndex]][[maxSelIndex]][[maxMigIndex]]))
  
  maxLoc = selSite[maxLocIndex]
  maxSel = sels[maxSelIndex]
  maxMig = migs[maxMigIndex]
  maxSource = sources[maxSourceIndex]
  
  return(cbind(maxLoc, maxSel, maxMig, maxSource))
}

getMCLEsv_source = function(compLike, selSite, sels, gs, times, sources) {
  # Gets MCLE for location of selected site strength of selection (s), frequency of standing variant (g),
  #   standing time (t), and source population for standing variant model
  #
  # Args: 
  #	compLike: nested lists of composite likelihoods (list of length length(selSite) where each
  #   element is a list of length length(sels) ... of length length(gs) ... of length length(times) ...
  #   of length length(sources) containing composite likelihood values)
  #	selSite: selected sites used to calculate compLike
  #	sels: selection coefficients used to calculate compLike
  # gs: frequencies of standing variant used to calculate compLike
  # times: standing times used to calculate compLike
  # sources: source populations used to calculate compLike
  #
  # Returns:
  #	MCLE for location of selected site, s, g, t, and source population
  maxLocIndex = ceiling(which.max(unlist(compLike)) / (length(sels) * length(gs) * length(times) * length(sources)))
  maxSelIndex = ceiling(which.max(unlist(compLike[[maxLocIndex]])) / (length(gs) * length(times) * length(sources)))
  maxGIndex = ceiling(which.max(unlist(compLike[[maxLocIndex]][[maxSelIndex]])) / (length(times) * length(sources)))
  maxTimeIndex = ceiling(which.max(unlist(compLike[[maxLocIndex]][[maxSelIndex]][[maxGIndex]])) / length(sources))
  maxSourceIndex = which.max(unlist(compLike[[maxLocIndex]][[maxSelIndex]][[maxGIndex]][[maxTimeIndex]]))
  
  maxLoc = selSite[maxLocIndex]
  maxSel = sels[maxSelIndex]
  maxG = gs[maxGIndex]
  maxTime = times[maxTimeIndex]
  maxSource = sources[maxSourceIndex]
  
  return(cbind(maxLoc, maxSel, maxG, maxTime, maxSource))
}

getMCLEmixed = function(compLike, selSite, sels, gs, times, migs, sources) {
  # Gets MCLE for location of selected site strength of selection (s), frequency of standing variant (g),
  #   standing time (t), migration rate (m), and source population for models with mixed modes
  #
  # Args: 
  #	compLike: nested lists of composite likelihoods (list of length length(selSite) where each
  #   element is a list of length length(sels) ... of length length(gs) ... of length length(times) ...
  #   of length length(migs) ... of length length(sources) containing composite likelihood values)
  #	selSite: selected sites used to calculate compLike
  #	sels: selection coefficients used to calculate compLike
  # gs: frequencies of standing variant used to calculate compLike
  # times: standing times used to calculate compLike
  # migs: migration rates used to calculate compLike
  # sources: source populations used to calculate compLike
  #
  # Returns:
  #	MCLE for location of selected site, s, g, t, m, and source population
  # Note: some of these values will have no meaning depending on which modes are specified
  maxLocIndex = ceiling(which.max(unlist(compLike)) / (length(sels) * length(gs) * length(times) * length(migs) * length(sources)))
  maxSelIndex = ceiling(which.max(unlist(compLike[[maxLocIndex]])) / (length(gs) * length(times) * length(migs) * length(sources)))
  maxGIndex = ceiling(which.max(unlist(compLike[[maxLocIndex]][[maxSelIndex]])) / (length(times) * length(migs) * length(sources)))
  maxTimeIndex = ceiling(which.max(unlist(compLike[[maxLocIndex]][[maxSelIndex]][[maxGIndex]])) / (length(migs) * length(sources)))
  maxMigIndex = ceiling(which.max(unlist(compLike[[maxLocIndex]][[maxSelIndex]][[maxGIndex]][[maxTimeIndex]])) / length(sources))
  maxSourceIndex = which.max(unlist(compLike[[maxLocIndex]][[maxSelIndex]][[maxMigIndex]][[maxTimeIndex]][[maxMigIndex]]))
  
  maxLoc = selSite[maxLocIndex]
  maxSel = sels[maxSelIndex]
  maxG = gs[maxGIndex]
  maxTime = times[maxTimeIndex]
  maxMig = migs[maxMigIndex]
  maxSource = sources[maxSourceIndex]
  
  return(cbind(maxLoc, maxSel, maxG, maxTime, maxSource))
}
  