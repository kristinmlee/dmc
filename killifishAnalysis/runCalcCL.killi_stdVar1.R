setwd("/Users/kristinlee/Desktop/killifishAnalysis")

det_FOmegas_stdVar_e = readRDS("det_FOmegas_stdVar_e.killifish_1.RDS")

inv_FOmegas_stdVar_e = readRDS("inv_FOmegas_stdVar_e.killifish_1.RDS")

rec = 2.17e-08

Ne = 8300000

numPops = 8
sets = list(c(2, 4, 6), 8)

selSiteInterest = 1:30

selPops = unlist(sets)
nonSelPops = seq(1, numPops)[- selPops]

positions = readRDS("positions.killifish.RDS")

selSite = seq(min(positions), max(positions), length.out = 30)[selSiteInterest]
distances = sapply(1:length(selSite), function(i) abs(positions - selSite[i]))

#mean centering
M = numPops
Tmatrix = matrix(data = rep(-1/M, (M-1)*M), nrow = M-1, ncol = M)
diag(Tmatrix) = (M-1)/M

freqs = readRDS("rand_freqs.killi.RDS")
epsilons = rowMeans(freqs)
freqs_MC = sapply(1:nrow(freqs), function(i) Tmatrix %*% freqs[i,])

sels = c(0.001, 0.005, seq(0.01, 0.05, by = 0.01), seq(0.06, 0.2, by = 0.02), 0.3, 0.4, 0.5, 0.6)
times = c(0, 5, 50, 100, 500, 1000, 5000, 10e6)
gs = 10^(-(2:10))
migs = c(0.00001, 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.1, 0.3, 0.5, 0.9, 1)
sources = c(2,4,6,8)[1]

sampleSizes = c(48, 48, 49, 50, 50, 43, 47, 49)*2
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

calcLikelihood_bin.4par = function(site, j, det_FOmegas, inv_FOmegas, par1, par2, par3, par4) {
	bin = distBins[site, j]
	my.x = as.matrix(freqs_MC[ ,site])
	my.e = epsilons[site]*(1-epsilons[site])

	likelihood = 1/(sqrt((2*pi)^k*(det_FOmegas[[par1]][[par2]][[par3]][[par4]][[bin]]*my.e^rank))) * exp(-1/2*t(my.x-mu) %*% (inv_FOmegas[[par1]][[par2]][[par3]][[par4]][[bin]]/my.e) %*% (my.x-mu))
	
	return(log(likelihood))
}

calcCompLikelihood_bin.4par = function(j, det_FOmegas, inv_FOmegas, par1, par2, par3, par4) {
	all = sapply(1:nrow(distances), {
			function(i) calcLikelihood_bin.4par(i, j, det_FOmegas, inv_FOmegas, par1, par2, par3, par4)
		})
	return(sum(all))
}

compLikelihood_stdVar = lapply(1 : ncol(distances), function(j) {
	lapply(1 : length(sels), function(sel) {
		lapply(1 : length(gs), function(g) {
			lapply(1 : length(times), function(t) {
				lapply(1: length(sources), function(my.source) {
				calcCompLikelihood_bin.4par(j, det_FOmegas_stdVar_e, inv_FOmegas_stdVar_e, sel, g, t, my.source)
				})
			})
		})
	})
})

saveRDS(compLikelihood_stdVar, paste("compLikelihood_stdVar_selSites_source", sources, ".killi.RDS", sep = ""))