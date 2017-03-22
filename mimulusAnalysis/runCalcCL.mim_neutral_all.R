setwd("/Users/kristinlee/Desktop/mimulusAnalysis")
library("MASS")

numPops = 4
M = numPops
Tmatrix = matrix(data = rep(-1/M, (M-1)*M), nrow = M-1, ncol = M)
diag(Tmatrix) = (M-1)/M

sampleSizes = c(31, 21, 20, 25)*2
sampleErrorMatrix = diag(1/sampleSizes, nrow = numPops, ncol = numPops)

F_estimate = readRDS("F_estimate.mim.RDS")

det_FOmegas_neutral = det(Tmatrix %*% (F_estimate + sampleErrorMatrix) %*% t(Tmatrix))
inv_FOmegas_neutral = ginv(Tmatrix %*% (F_estimate + sampleErrorMatrix) %*% t(Tmatrix))
saveRDS(det_FOmegas_neutral, "det_FOmegas_neutral.mim_all.RDS")
saveRDS(inv_FOmegas_neutral, "inv_FOmegas_neutral.mim_all.RDS")

selSiteInterest = 1:31

rec = 4.72e-8
Ne = 750000

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


##bin distances
numBins = 1000
my.seq = seq(min(distances) - 0.001, max(distances) + 0.001, length.out = (numBins + 1))
midDistances = sapply(1:numBins, function(i) mean(c(my.seq[i], my.seq[i+1])))

distBins = apply(distances, 2, function(i) as.numeric(cut(i, my.seq)))

#MVN parameters
k = numPops-1
mu = as.matrix(rep(0, k))
rank = numPops - 1

calcLikelihood_bin_neutral = function(site, det_FOmegas, inv_FOmegas, j) {
	bin = distBins[site, j]
	my.x = as.matrix(freqs_MC[ ,site])
	my.e = epsilons[site]*(1-epsilons[site])

	likelihood = 1/(sqrt((2*pi)^k*(det_FOmegas*my.e^rank))) * exp(-1/2*t(my.x-mu) %*% (inv_FOmegas/my.e) %*% (my.x-mu))
	return(log(likelihood))
}

calcCompLikelihood_bin_neutral = function(j, det_FOmegas, inv_FOmegas) {
	all = sapply(1:nrow(distances), {
			function(i) calcLikelihood_bin_neutral(i, det_FOmegas, inv_FOmegas, j)
		})
	return(sum(all))
}

compLikelihood_neutral = sapply(1 : ncol(distances), function(j) calcCompLikelihood_bin_neutral(j, det_FOmegas_neutral, inv_FOmegas_neutral))

saveRDS(compLikelihood_neutral, "compLikelihood_neutral.mim.RDS")
