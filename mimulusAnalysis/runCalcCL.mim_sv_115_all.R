setwd("/Users/kristinlee/Desktop/mimulusAnalysis")

det_FOmegas_stdVar_e = readRDS("det_FOmegas_stdVar_e.mim.RDS")

inv_FOmegas_stdVar_e = readRDS("inv_FOmegas_stdVar_e.mim.RDS")

selSiteInterest = 1:15

rec = 4.72e-8
Ne = 750000

numPops = 4
sets = c(1,3)
sisterPops = c(2, 4)
selPops = unlist(sets)

#mean centering
M = numPops
Tmatrix = matrix(data = rep(-1/M, (M-1)*M), nrow = M-1, ncol = M)
diag(Tmatrix) = (M-1)/M

freqs = readRDS("rand_freqs_all.mim.RDS")

epsilons = rowMeans(freqs)
freqs_MC = sapply(1:nrow(freqs), function(i) Tmatrix %*% freqs[i,])


positions = readRDS("mimPositions_all.RDS")
selSite = sort(c(seq(min(positions), max(positions), length.out = 30), 309000))[selSiteInterest]
distances = sapply(1:length(selSite), function(i) abs(positions - selSite[i]))

sels = c(seq(0.001,0.01,length=10),seq(0.011,0.06,length=20), 0.08, seq(0.1, 0.6, by = 0.05))
times = c(5, seq(10,1000,length=15),seq(1500,3000,length=15))
gs = 10^(-(2:10))
migs = c(0.00001, 0.0001, 0.0005, 0.001, 0.005, 0.01, seq(0.1, 1, by = 0.1))
sources = c(1,3)

sampleSizes = c(31, 21, 20, 25)*2
sampleErrorMatrix = diag(1/sampleSizes, nrow = numPops, ncol = numPops)

##bin distances
numBins = 1000
my.seq = seq(min(distances) - 0.001, max(distances) + 0.001, length.out = (numBins + 1))
midDistances = sapply(1:numBins, function(i) mean(c(my.seq[i], my.seq[i+1])))

distBins = apply(distances, 2, function(i) as.numeric(cut(i, my.seq)))

#MVN parameters
k = numPops-1
mu = as.matrix(rep(0, k))
rank = numPops - 1

calcLikelihood_bin.3par = function(site, j, det_FOmegas, inv_FOmegas, par1, par2, par3) {
	bin = distBins[site, j]
	my.x = as.matrix(freqs_MC[ ,site])
	my.e = epsilons[site]*(1-epsilons[site])

	likelihood = 1/(sqrt((2*pi)^k*(det_FOmegas[[par1]][[par2]][[par3]][[bin]]*my.e^rank))) * exp(-1/2*t(my.x-mu) %*% (inv_FOmegas[[par1]][[par2]][[par3]][[bin]]/my.e) %*% (my.x-mu))
	
	return(log(likelihood))
}

calcCompLikelihood_bin.3par = function(j, det_FOmegas, inv_FOmegas, par1, par2, par3) {
	all = sapply(1:nrow(distances), {
			function(i) calcLikelihood_bin.3par(i, j, det_FOmegas, inv_FOmegas, par1, par2, par3)
		})
	return(sum(all))
}

compLikelihood_stdVar = lapply(1 : ncol(distances), function(j) {
	lapply(1 : length(sels), function(sel) {
		lapply(1 : length(gs), function(g) {
			lapply(1 : length(times), function(t) {
				calcCompLikelihood_bin.3par(j, det_FOmegas_stdVar_e, inv_FOmegas_stdVar_e, sel, g, t)
			})
		})
	})
})

saveRDS(compLikelihood_stdVar, paste("compLikelihood_stdVar_selSites_", paste(selSiteInterest, collapse = "_"), ".mim.RDS", sep = ""))