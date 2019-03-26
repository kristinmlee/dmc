## Script to generate and save estimate of matrix of neutral probabilities of coalescing (F)
##
## Uncommented lines show how we did this for the AHR analysis
## There are comments speicifying when the code was different for the ARNT analysis
## The neutral F matrix is different for each region's analysis because we included
##   different individuals in the populations, based on who contained selected haplotypes
##   (the sweeps are not complete and the model assumes that alleles sampled in the selected
##    population are linked to the selected allele)
##
## Args:
##	neutral_freqs: matrix of allele frequencies at putatively neutral sites with	
##		dimension numberOfPopulations x numberOfSites
##    -- allele frequencies should be randomized wrt the reference allele
##	sampleSizes: vector of sample sizes of length numberOfPopulations
##		(i.e. number of chromosomes sampled in each population)
##    Note: should be in same order as populations listed in neutral_freqs

# matrix of allele frequencies at putatively neutral sites that includes monomorphic sites
# the reference allele is already randomized
# column names are the population labels
allFreqs=readRDS("~/Documents/grandis/popAlleleRandFreqs_grandis_neutral.RDS") # allele freq matrix for individuals analyzed for AHR analysis
# allFreqs=readRDS("~/Documents/grandis/popAlleleRandFreqs_grandis_neutral_chr10indiv.RDS")  # allele freq matrix for individuals analyzed for ARNT analysis

# get dimensions numPops x numSites
allRunFreq = t(allFreqs)

# these are the populations we are including in the analysis
# AHR name format (below): [population initials]_[species][number of copies of deletion]
keep_pops=c("VB_g2","SP_g0","ER_h2","KC_h0","BP_h0","F_h0")
# ARNT allele freq matrix only includes populations of interest; name format: [population initials]; c("VB","SP","ER","KC","BP","F")

# subset allele frequency matrix to include only the populations used in analysis
allRunFreq=allRunFreq[keep_pops,]

# Identify and exclude monomorphic sites from allRunFreq:
colsWithAllFixed = which(apply(allRunFreq, 2, function(j) sum(j) == nrow(allRunFreq)))
colsWithAllLost = which(apply(allRunFreq, 2, function(j) sum(j) == 0))
allRunFreq = allRunFreq[, -c(colsWithAllFixed,colsWithAllLost)]

# Calculate Neutral F Matrix
numLoci = ncol(allRunFreq)
my.means.rand = (allRunFreq %*% t(allRunFreq)) / numLoci
sampleSizes = readRDS("sampleSizes.RDS")

diag(my.means.rand) = diag(my.means.rand) * sampleSizes / (sampleSizes - 1) - rowMeans(allRunFreq) /
  (sampleSizes - 1)

dist.ij = which(my.means.rand == min(my.means.rand), arr.ind = TRUE)[1, ]

A.rand = mean(allRunFreq[dist.ij[1], ] * allRunFreq[dist.ij[2], ])
C.rand = mean(allRunFreq[dist.ij[1], ] * (1 - allRunFreq[dist.ij[2], ]))
F_estimate = (my.means.rand - A.rand) / C.rand

# populations on opposite sides of the root have negative covariance (but practically zero), 
# even though my.means.rand-A.rand should equal zero so we just set those values to be zero
F_estimate[dist.ij[1],dist.ij[2]]=0
F_estimate[dist.ij[2],dist.ij[1]]=0

# save F_estimate 
saveRDS(F_estimate, "neutralF.RDS")