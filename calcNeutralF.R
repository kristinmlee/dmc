# Script to generate and save estimate of neutral variance/covariance matrix (F)
#
# Args:
#	allFreqs: matrix of allele frequencies at putatively neutral sites with	
#		dimension numberOfPopulations x numberOfSites
#	sampleSizes: vector of sample sizes of length numberOfPopulations
#		(i.e. twice the number of diploid individuals sampled in each population)
#	neutralF_filename: name/path as string for neutral variance/covariance matrices
#		to be saved as an R object

#randomize reference allele
allRunFreq = apply(allFreqs, 2, function(my.freqs) {
	if(runif(1) < 0.5) {
		my.freqs = 1 - my.freqs
	}
	my.freqs
})

numLoci = ncol(allRunFreq)
my.means.rand = (allRunFreq %*% t(allRunFreq)) / numLoci

diag(my.means.rand) = diag(my.means.rand) * sampleSizes / (sampleSizes - 1) - rowMeans(allRunFreq) /
	(sampleSizes - 1)

dist.ij = which(my.means.rand == min(my.means.rand), arr.ind = TRUE)[1, ]

A.rand = mean(allRunFreq[dist.ij[1], ] * allRunFreq[dist.ij[2], ])
C.rand = mean(allRunFreq[dist.ij[1], ] * (1 - allRunFreq[dist.ij[2], ]))

F_estimate = (my.means.rand - A.rand) / C.rand
saveRDS(F_estimate, paste(neutralF_filename, ".RDS", sep = "", collapse = ""))