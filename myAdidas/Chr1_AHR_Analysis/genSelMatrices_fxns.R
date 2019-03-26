## genSelMatrices_fxns.R
##
## These are a set of functions that are saved and used in genSelMatrices_exec.R
## They generate expected selection matrices for each model according to the
## distance bin away from selected site and given set of parameters 

# Function to calculate probability of not recombining out during the sweep
prob_no_rec_out=function(x,s,r) {
  # x is starting frequency
  # s is strength of selection against homozygote for intact haplotypes
  # r is recombination rate between site of interest and selected site
  freq=x # starting frequency of beneficial variant
  t=1 #starting generation, used as a tracker for index to store frequency
  while(freq[t]<0.85) {
    p=freq[t]
    w.bar=p^2 + (1-0.5*s)*2*p*(1-p) + (1-s)*(1-p)^2
    freq[t+1]=(p^2 + (1-0.5*s)*p*(1-p))/w.bar
    t=t+1
  }
  no_rec_out=1-r*(1-freq) # vector of probability not recombining out in a generation
  return(prod(no_rec_out))
}

saveRDS(prob_no_rec_out,"prob_no_rec_out.RDS")

### Standing Variant Model (ILS) ### 
calctotAddF_stdVar = function(y, r, rt, g){
  # Computes mean-centered and sample-size corrected 
  #	coancestry matrices with convergent adaptation (F^(S))
  #	due to selection on shared ancestral standing variation
  #	*WITHOUT a source of the standing variant specified* for a 
  #	single distance bin
  #
  # Args:
  #	y: vector or length numPops for the probability of recombining off the 
  # 		beneficial background of the sweep (calculated by function prob_no_rec_out)
  # r: recombation rate between given distance bin and selected site
  #	rt: the probability that a single lineage does not
  #		recombine off onto the non-beneficial background
  #		during the standing phase for t generations (a
  #		function of recombination distance)
  #	g: frequency of the standing variant
  #
  # Global variables used (need to assign outside of function):
  #   selPops: vector of names of populations under selection
  #
  # Calculated within function:
  #	Rf: a function of the recombination distance & Ne
  #		that is used to specify whether two neutral lineages 
  # 		coalesce or recombine off the standing phase
  #
  # Returns:
  #	mean-centered and sample-size corrected coancestry matrix
  
  # Rf varies by population, depending on Ne
  Rf = sapply(Ne, function(N) { 4 * N * r * g })
  
  selMatrix = F_estimate
  for(i  in selPops) {
    selMatrix[i,i] = (1 - y^2) * (F_estimate[i, i]) + y^2 * (1 / (1 + Rf[i]) + Rf[i] / (1 + Rf[i]) * (F_estimate[i, i]))
    
    for(j  in selPops) {
      if(i != j) {
        RF = mean(Rf[i],Rf[j])
        selMatrix[i, j] = (1 - y^2) * F_estimate[i, j] + y^2 * ((rt^2 *(1 / (1 + RF) + 
                                                                          RF / (1 + RF) * F_estimate[i, j])) + (1 - rt^2) * F_estimate[i, j])
      }
    }
  }
  return(Tmatrix %*% (selMatrix + sampleErrorMatrix) %*% t(Tmatrix))
}
saveRDS(calctotAddF_stdVar, file="calctotAddF_stdVar.RDS")


calcFOmegas_stdVar = function(sel, G, time) {
  # Generates mean-centered and sample-size corrected 
  #	coancestry matrices with convergent
  #	selection (F^(S)) for all bins under standing variant model
  #	*WITHOUT source* for a given set of parameters
  #
  # Args:
  #	sel: strength of selection
  #	G: frequency of the standing variant
  #	time: time in generations the variant is standing in populations
  #		between when the selected populations split and the onset of selection
  #
  # Returns:
  #	list of length numBins of coancestry matrices 
  #  with convergent selection (F^(S)) 
  
  s.index=which(sapply(1:length(sels), function(i) {isTRUE(all.equal(sel,sels[i]))}))
  g.index=which(sapply(1:length(gs), function(i) {isTRUE(all.equal(G,gs[i]))}))
  y=ys.stand[[s.index]][[g.index]]
  rt = exp(-rec*midDistances * time)
  R = rec*midDistances
  FOmegas = lapply(1 : length(midDistances), function(i) calctotAddF_stdVar(y=y[[i]], r=R[[i]], rt=rt[[i]], g=G))
  return(FOmegas)
}

saveRDS(calcFOmegas_stdVar,file="calcFOmegas_stdVar.RDS")

### Standing Variant Source Model ### 
calctotAddF_stdVar.source = function(y, rt, g, r, p_no, p_one ) {
  # Computes mean-centered and sample-size corrected 
  #	coancestry matrices with convergent adaptation (F^(S))
  #	due to selection on shared ancestral standing variation
  #	*with a source of the standing variant specified* for a 
  #	single distance bin
  #
  # Args:
  #	y: the probability of recombining off the beneficial background of the sweep 
  #    (calculated by function prob_no_rec_out). (Single number -- not a vector)
  #	Rf: a function of the recombination distance
  #		that is used to specify whether two neutral lineages 
  # 		coalesce or recombine off the standing phase; 
  #   changes depending on which pop(s) are considered (because of differences in Ne)
  #	rt: the probability that a single lineage does not
  #		recombine off onto the non-beneficial background
  #		during the standing phase for t generations (a
  #		function of recombination distance)
  # g: frequency of the standing variant
  # r: recombation rate between given distance bin and selected site
  #	p_no: the probability no lineages recombine off or coalesce
  #		during time t (a function of recombination distance) *within the recipient selected pop
  #	p_one: the probability one lineage recombines off
  #		during time t (a function of recombination distance) *within the recipient selected pop
  #
  # Global variables used (need to assign outside of function):
  #	  source: source population of the beneficial allele for migration model
  #   recipient: grandis pop that sweeps deletion
  #
  # Calculated within function:
  #	Rf: a function of the recombination distance & Ne
  #		that is used to specify whether two neutral lineages 
  # 		coalesce or recombine off the standing phase
  #
  # Returns:
  #	mean-centered and sample-size corrected coancestry matrix
  
  #	my.source: the proposed source of the beneficial allele
  #		*Note: the source must be a selected population in "selPops"
  selMatrix = F_estimate
  
  # within source pop (that is selected)
  Rf_s = 4 * Ne[source] * r * g # Rf source
  selMatrix[source, source] = (1 - y^2) * (F_estimate[source, source]) + y^2 * (1 / (1 + Rf_s) + 
                                                                                              Rf_s / (1 + Rf_s) * (F_estimate[source, source]))
  for(i in recipient) {
    Rf_r = 4 * Ne[i] * r * g # Rf recipient
    
    selMatrix[i,i] = (1 - y)^2 * (F_estimate[i, i]) + y^2 * (p_no*(1 / (1+Rf_s) + Rf_s / (1 + Rf_s) *
                                                                     F_estimate[source, source]) + (1 - p_no) * (1 / (1 + Rf_r)) + ((1 - p_no) * Rf_r / (1 + Rf_r) -
                                                                                                                                            Rf_r / (1 + Rf_r / 2) * (1 - p_one) * rt) * F_estimate[i,i] + (1 - p_one) * Rf_r / (1 + Rf_r / 2) *
                                                               rt * F_estimate[i, source]) + 2 * y * (1 - y) * (rt * F_estimate[i, source] + (1 - rt) *
                                                                                                                     F_estimate[i, i])
    
    selMatrix[i, source] = (1 - y)^2 * F_estimate[i, source] + y^2 * (rt^2 * (1 / (1 + Rf_s) + 
                                                                                      Rf_s / (1 + Rf_s) * F_estimate[source, source]) + rt * (1 - rt) * 
                                                                              (F_estimate[source, source]) + rt * (1 - rt) * F_estimate[i, source] + (1 - rt)^2 *
                                                                              F_estimate[i, source]) + y * (1 - y) * (rt * F_estimate[source, source] + (1 - rt) *
                                                                                                                           F_estimate[i, source]) + y * (1 - y) * F_estimate[i, source]
    
    selMatrix[source, i] = selMatrix[i, source]			
    
    for(k in nonSelPops) {
      selMatrix[k, i] = y * rt * F_estimate[k, source] + (1 - y) * F_estimate[i, k] + 
        y * (1 - rt) * F_estimate[i, k]
      selMatrix[i, k] = y * rt * F_estimate[k, source] + (1 - y) * F_estimate[i, k] + 
        y * (1 - rt) * F_estimate[i, k]
    }
  }
  return(Tmatrix %*% (selMatrix + sampleErrorMatrix) %*% t(Tmatrix))
}

saveRDS(calctotAddF_stdVar.source,file="calctotAddF_stdVar_source.RDS")

calcFOmegas_stdVar.source = function(sel, G, time) {
  # Generates mean-centered and sample-size corrected 
  #	coancestry matrices with convergent
  #	selection (F^(S)) for all bins under standing variant model
  #	*with source* for a given set of parameters
  #
  # Args:
  #	sel: strength of selection
  #	G: frequency of the standing variant
  #	time: time in generations the variant is standing in populations
  #		between migration from source population and onset of selection
  #
  # Returns:
  #	list of length numBins of coancestry matrices 
  #		with convergent selection (F^(S)) 
  
  s.index=which(sapply(1:length(sels), function(i) {isTRUE(all.equal(sel,sels[i]))}))
  g.index=which(sapply(1:length(gs), function(i) {isTRUE(all.equal(G,gs[i]))}))
  y=ys.stand[[s.index]][[g.index]]
  
  rt = exp(-rec*midDistances * time) # prob of not recombining out in time t within recipient pop
  R = rec*midDistances
  
  p_No = exp(-time * (2 * rec*midDistances + 1/(2 * Ne[recipient] * G)))
  p_One = exp(-time * (rec*midDistances + 1/(2 * Ne[recipient] * G)))
  
  FOmegas = lapply(1 : length(midDistances), function(i) calctotAddF_stdVar.source(y[[i]], rt[[i]], g=G, r=R[[i]], p_no=p_No[[i]],p_one=p_One[[i]]))
  return(FOmegas)
}

saveRDS(calcFOmegas_stdVar.source,file="calcFOmegas_stdVar_source.RDS")

### Staggered Sweeps Model ###
calctotAddF_mig.stagSweeps = function (y, g, t_m, r) {
  # Computes mean-centered and sample-size corrected 
  #	coancestry matrices with convergent adaptation (F^(S))
  #	due to a sweep on a new variant in the source pop, followed by
  # migration to the recipient pop and sweep, for a single distance bin
  #
  # Args:
  # y: vector or length numPops for the probability of recombining off the 
  # 		beneficial background of the sweep (calculated by function prob_no_rec_out)
  # g: proportion of alleles in recipient pop that came from source pop
  # t_m: number of generations between end of sweep and migration (looking pastwards)
  # r: recombination rate for that bin (between bin and selected site)
  #
  # Global variables used (need to assign outside of function):
  #	  source: source population of the beneficial allele for migration model
  #   recipient: grandis pop that sweeps deletion
  #
  # Returns:
  #	mean-centered and sample-size corrected coancestry matrix
  e_m = exp(-r*t_m)
  
  selMatrix = F_estimate
  
  # within selected heteroclitus pop
  selMatrix[source,source] = F_estimate[source,source] + y[source]^2 * (1 - F_estimate[source, source])
  
  # probabilities that a single lineage sampled in the recipient pop does or doesn't migrate back to source
  p_recip_mig=y[recipient]*(e_m + (1-e_m)*g) +(1-y[recipient])*g # prob lineage in recipient migrates to source
  p_recip_no_mig=y[recipient]*(1-e_m)*(1-g) + (1-y[recipient])*(1-g) # prob lineage in source doesn't migrate to source
  
  # between source pop and recipient pop ** make sure you modify sel matrix for within source pop first
  selMatrix[source,recipient] =selMatrix[source,source]*p_recip_mig + F_estimate[recipient,source]*p_recip_no_mig
  selMatrix[recipient,source] = selMatrix[source, recipient]
  
  # within recipient pop ** make sure you modify sel matrix for within source pop first
  z1 = g^2 *selMatrix[source,source] + 2*g*(1-g)*F_estimate[source,recipient] + (1-g)^2*F_estimate[recipient,recipient]
  z2 = g*selMatrix[source,source] + (1-g)*F_estimate[source,recipient]
  w = 1 - exp(-t_m*(2*r + 1/(2*Ne[recipient]*g)))
  term1a= 1/(1+4*Ne[recipient]*g*r) + ((4*Ne[recipient]*g*r)/(1+4*Ne[recipient]*g*r))*(e_m*z2 + (1-e_m)*z1)
  term1=w*term1a + (1-w)*selMatrix[source,source]
  term2=e_m*z2 + (1-e_m)*z1
  selMatrix[recipient,recipient] = y[recipient]^2*term1 + (1-y[recipient])^2*z1 + 2*y[recipient]*(1-y[recipient])*term2 
  
  #between recipient pop and non-source pops:
  for (pop in setdiff(1:numPops,c(source,recipient))){
    selMatrix[recipient,pop] = F_estimate[recipient,pop]*p_recip_no_mig + F_estimate[source,pop]*p_recip_mig
    selMatrix[pop,recipient]=selMatrix[recipient,pop]
  }
  
  return(Tmatrix %*% (selMatrix + sampleErrorMatrix) %*% t(Tmatrix))
  
}

saveRDS(calctotAddF_mig.stagSweeps,file="calctotAddF_mig_stagSweeps.RDS")

calcFOmegas_mig.stagSweeps = function(sel, G, standing.time) {
  # Generates mean-centered and sample-size corrected 
  #	coancestry matrices with convergent adaptation (F^(S)) for all bins
  #	under model: migration from a source pop before the sweep in 
  # either pop, followed by same deletion + sweep 
  #
  # Args:
  #	sel: strength of selection
  # G: proportion of alleles in recipient pop that came from source pop
  # standing.time: probability of not recombining out during standing time if linked to deletion
  #      (function of recombination distance)
  #
  #
  # Returns:
  #	list of length numBins of coancestry matrices 
  #		with convergent selection (F^(S)) 
  s.index=which(sapply(1:length(sels), function(i) {isTRUE(all.equal(sel,sels[i]))}))
  g.index=which(sapply(1:length(gs), function(i) {isTRUE(all.equal(G,gs[i]))}))
  y.ind=ys.ind[[s.index]] # returns list of y corresponding to each rec distance for the selection coefficient
  y.stand=ys.stand[[s.index]][[g.index]]
  y = lapply(1:length(midDistances), function(i) {
    sapply(1:length(pops),function(p) {
      # if Ne = pop size of heteroclitus: start from freq 1/2N, otherwise start from freq G
      # this is because we know that the heteroclitus pop is the source and thus starts from 1/2N
      ifelse(Ne[p]==8.3e6,y.ind[[i]][p],y.stand[[i]])
    })
  })
  
  t_M = standing.time
  R = rec*midDistances
  
  FOmegas = lapply(1 : length(midDistances), function(i) calctotAddF_mig.stagSweeps(y=y[[i]], g=G, t_m=t_M, r=R[[i]]))
  return(FOmegas)
}

saveRDS(calcFOmegas_mig.stagSweeps,file="calcFOmegas_mig_stagSweeps.RDS")

### Concurrent Sweeps Model ###
# DMC originally calls this migration
calctotAddF_mig.concSweeps = function(y, e_delta, my.Q){
  # Computes mean-centered and sample-size corrected 
  #	coancestry matrices with convergent adaptation (F^(S))
  #	due to migration for a single distance bin
  #
  # Args:
  #	y: vector or length numPops for the probability of recombining off the 
  # 		beneficial background of the sweep (calculated by function prob_no_rec_out)
  #	e_delta: the probability of recombining off the 
  # 		beneficial background of the sweep in the
  #		source population for time delta (a
  #		function of recombination distance)
  #	my.Q: the probability of coalescing before recombination
  #		at the selected site
  #
  # Global variables used (need to assign outside of function):
  #	  source: source population of the beneficial allele for migration model
  #   recipient: grandis pop that sweeps deletion
  #   selPops: vector of names of populations under selection
  #
  # Returns:
  #	mean-centered and sample-size corrected coancestry matrix
  
  selMatrix = F_estimate
  if(is.element(source, selPops)) {
    selMatrix[source, source] = (F_estimate[source, source]) +
      y[source]^2 * (1 - (F_estimate[source, source]))
    
    for(i in selPops[selPops != source]) {
      selMatrix[i, i] = my.Q * (y[i]^2 + (1 - y[i]^2) * (F_estimate[i, i])) + (1 - my.Q) * 
        (y[i]^2 * e_delta^2 + (1 - y[i])^2 * (F_estimate[i, i]) + 2 * y[i] * (1 - y[i]) * 
           F_estimate[source, i] + y[i]^2 * (1 - e_delta^2) * (F_estimate[source, source]))
      selMatrix[i, source] = y[i]*y[source] * e_delta + (1 - y[i]) * F_estimate[source, i] + 
        y[i] * (1 - y[source] * e_delta) * (F_estimate[source, source])
      selMatrix[source, i] = selMatrix[i, source]
      
      for(k in nonSelPops) {
        selMatrix[k, i] = (1 - y[i]) * F_estimate[i, k] + y[i] * F_estimate[source, k]
        selMatrix[i, k] = (1 - y[i]) * F_estimate[i, k] + y[i] * F_estimate[source, k]
      }
      
      for(j  in selPops[selPops != source]) {
        if(i != j)
          selMatrix[i, j] = y[i]*y[j] * e_delta^2 + y[i]*y[j] * (1 - e_delta^2) * (F_estimate[source, source]) +
            (1 - y[i])*(1-y[j]) * F_estimate[i,j]	+ (1 - y[i]) * y[j] * F_estimate[i, source] +
            (1 - y[j]) * y[i] * F_estimate[j, source]
      }
    }
  }
  return(Tmatrix %*% (selMatrix + sampleErrorMatrix) %*% t(Tmatrix))
}

saveRDS(calctotAddF_mig.concSweeps,file="calctotAddF_mig_concSweeps.RDS")

calcFOmegas_mig.concSweeps = function(sel, mig) {
  # Generates mean-centered and sample-size corrected 
  #	coancestry matrices with convergent
  #	selection (F^(S)) for all bins under concurrent sweeps model 
  #	for a given set of parameters
  #
  # Args:
  #	sel: strength of selection
  #	mig: migration rate (proportion of individuals from source each generation)
  #
  # Returns:
  #	list of coancestry matrices with convergent
  #		selection (F^(S)) of length numBins
  s.index=which(sapply(1:length(sels), function(i) {isTRUE(all.equal(sel,sels[i]))}))
  y=ys.ind[[s.index]] # returns list of y corresponding to each rec distance for the selection coefficient, where each element is vector of length num pops
  
  delta = 1 / sel * log(1 + sel/(mig))
  e_delta = exp(-rec * midDistances * delta)
  my.Q = 1 / (1 + 4 * Ne[recipient] * mig)
  FOmegas = lapply(1 : length(midDistances), function(i) calctotAddF_mig.concSweeps(y[[i]], e_delta[[i]], my.Q))
  return(FOmegas)
}

saveRDS(calcFOmegas_mig.concSweeps,file="calcFOmegas_mig_concSweeps.RDS")
