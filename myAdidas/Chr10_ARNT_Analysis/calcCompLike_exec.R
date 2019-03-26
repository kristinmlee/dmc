# Implementation of calculating composite log-likelihoods under all models of convergent selection
#
# This script uses functions saved in calcCompLike_fxns.R
#
# Args:
#	F_estimate: estimate of neutral variance/covariance matrix
#		generated with "calcNeutralF.R"
#	sampleSizes: vector of sample sizes (# chromosomes sampled) of length numPops + same order of populations as F_estimate 
#	numPops: number of populations sampled (both selected and non-selected)
#	positions: vector of genomic positions for region
# freqs: matrix of allele frequencies for window of interest with	
#		dimension numberOfPopulations x numberOfSites
# numBins: the number of bins in which to bin alleles a given distance from the proposed
#		selected sites
#   NOTE: must be same as number specified in generating covariance matrices
# Inverses and determinants for all relevant models of convergent selection
#	Parameter spaces for likelihood calculations *must be same as those in genSelMatrices.R

F_estimate=readRDS("neutralF.RDS")
sampleSizes=readRDS("sampleSizes.RDS")
numPops=length(sampleSizes)
pops=colnames(F_estimate)

#loading in parameters used in genSelMatrices_exec.R
selSite=readRDS("selSite.RDS")
sels=readRDS("sels.RDS")
times=readRDS("times.RDS")
stdVar_times=readRDS("stdVar_times.RDS")
gs=readRDS("gs.RDS")
migs=readRDS("migs.RDS")

sampleErrorMatrix = diag(1/sampleSizes, nrow = numPops, ncol = numPops)

##Multivariate Normal parameters
k = numPops - 1
mu = as.matrix(rep(0, k))
rank = numPops - 1

##mean centering
M = numPops
Tmatrix = matrix(data = rep(-1 / M, (M - 1) * M), nrow = M - 1, ncol = M)
diag(Tmatrix) = (M - 1) / M 

# positions and frequencies correspond to each other
positions=readRDS("target_freqs.RDS")$pos

freqs=readRDS("target_freqs.RDS")[,pops]
epsilons=rowMeans(freqs)
freqs_MC = sapply(1:nrow(freqs), function(i) Tmatrix %*% t(freqs[i,])) # MC = mean centered
mu=rep(0,k) # because we mean_centered

## Categorize neutral sites into bins of distances to each selected site.
#load distances and my.seq from genSelMatrices
distances = readRDS("distances.RDS") # #distance of each position to selected site (row), depending on which site is the selected one (column)
my.seq = readRDS("myseq.RDS")
# distBins: matrix of dimensions numSites x numSelSites (proposed)
# for each sel site (column), get distance bin that each position 
# you are calculating likelihood of belongs to
distBins = apply(distances, 2, function(i) as.numeric(cut(i, my.seq)))

# Load generic functions to calculate composite likelihoods
calcCompLikelihood_neutral=readRDS(file="calcCompLikelihood_neutral.RDS")
calcLikelihood_bin_neutral=readRDS(file="calcLikelihood_bin_neutral.RDS")
calcLikelihood_bin.2par=readRDS(file="calcLikelihood_bin_2par.RDS")
calcCompLikelihood.2par=readRDS(file="calcCompLikelihood_2par.RDS")
calcLikelihood_bin.3par=readRDS(file="calcLikelihood_bin_3par.RDS")
calcCompLikelihood.3par=readRDS(file="calcCompLikelihood_3par.RDS")

# save image so that you can actually make new files to run these separately
save.image("calcCompLike_ARNT.Rdata")

#################### IMPLEMENTATION #################### (make separate scripts if you want to run in parallel)

### Neutral model 
# returns same value in vector of length selSite
det_FOmegas_neutral=readRDS(file="det_FOmegas_neutral_ARNT.RDS")
inv_FOmegas_neutral=readRDS(file="inv_FOmegas_neutral_ARNT.RDS")
compLikelihood_neutral = sapply(1 : ncol(distances), function(sel_site) calcCompLikelihood_neutral(sel_site, det_FOmegas_neutral,inv_FOmegas_neutral))
saveRDS(compLikelihood_neutral,"compLikelihood_neutral_ARNT.RDS")

## Standing variant model
# three parameters: sels, gs, times (+ location of selected site)
det_FOmegas_stdVar=readRDS(file="det_FOmegas_stdVar_ARNT.RDS")
inv_FOmegas_stdVar=readRDS(file="inv_FOmegas_stdVar_ARNT.RDS")
compLikelihood_stdVar= lapply(1 : ncol(distances), function(sel_site) {
  lapply(1:length(sels), function(s){
    lapply(1:length(gs), function(g){
      lapply(1:length(stdVar_times), function(t){
        calcCompLikelihood.3par(sel_site, det_FOmegas_stdVar, inv_FOmegas_stdVar, s, g, t)
      })
    })
  })
})
saveRDS(compLikelihood_stdVar,"compLikelihood_stdVar_45.RDS")

## Standing variant source model
# three parameters: sels, gs, times (+ location of selected site)
det_FOmegas_stdVar.source=readRDS(file="det_FOmegas_stdVar_source_ARNT.RDS")
inv_FOmegas_stdVar.source=readRDS(file="inv_FOmegas_stdVar_source_ARNT.RDS")
compLikelihood_stdVar.source= lapply(1 : ncol(distances), function(sel_site) {
  lapply(1:length(sels), function(s){
    lapply(1:length(gs), function(g){
      lapply(1:length(times), function(t){
        calcCompLikelihood.3par(sel_site, det_FOmegas_stdVar.source, inv_FOmegas_stdVar.source, s, g, t)
      })
    })
  })
})
saveRDS(compLikelihood_stdVar.source,"compLikelihood_stdVar_source_45.RDS")

#### Migration with staggered sweeps ####
# three parameters: sels, gs, times (+ location of selected site)
det_FOmegas_mig.stagSweeps=readRDS(file="det_FOmegas_mig_stagSweeps_ARNT.RDS")
inv_FOmegas_mig.stagSweeps=readRDS(file="inv_FOmegas_mig_stagSweeps_ARNT.RDS")
compLikelihood_mig.stagSweeps=lapply(1 : ncol(distances), function(sel_site) {
  lapply(1:length(sels), function(s){
    lapply(1:length(gs), function(g){
      lapply(1:length(times), function(t){
        calcCompLikelihood.3par(sel_site, det_FOmegas_mig.stagSweeps, inv_FOmegas_mig.stagSweeps, s, g, t)
      })
    })
  })
})
saveRDS(compLikelihood_mig.stagSweeps,"compLikelihood_mig_stagSweeps_ARNT.RDS")

### Migration with concurrent sweeps ####
# two parameters: sels, migs (+ location of selected site)
det_FOmegas_mig.concSweeps=readRDS(file="det_FOmegas_mig_concSweeps_ARNT.RDS")
inv_FOmegas_mig.concSweeps=readRDS(file="inv_FOmegas_mig_concSweeps_ARNT.RDS")
compLikelihood_mig.concSweeps=lapply(1 : ncol(distances), function(sel_site) {
  lapply(1:length(sels), function(s){
    lapply(1:length(migs), function(m){
      calcCompLikelihood.2par(sel_site,det_FOmegas_mig.concSweeps, inv_FOmegas_mig.concSweeps, s, m)
    })
  })
})
saveRDS(compLikelihood_mig.concSweeps,"compLikelihood_mig_concSweeps_ARNT.RDS")


