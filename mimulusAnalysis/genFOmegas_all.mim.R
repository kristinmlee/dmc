library("MASS")

F_estimate = readRDS("F_estimate.mim.RDS")

rec = 4.72e-8
Ne = 750000

numPops = 4
sets = c(1,3)
sisterPops = c(2, 4)
selPops = unlist(sets)

positions = readRDS("mimPositions_all.RDS")
selSite = seq(min(positions), max(positions), length.out = 30)

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

#MVN parameters
k = numPops-1
mu = as.matrix(rep(0, k))
rank = numPops - 1

#mean centering
M = numPops
Tmatrix = matrix(data = rep(-1/M, (M-1)*M), nrow = M-1, ncol = M)
diag(Tmatrix) = (M-1)/M

det_FOmegas_neutral = det(Tmatrix %*% (F_estimate + sampleErrorMatrix) %*% t(Tmatrix))
inv_FOmegas_neutral = ginv(Tmatrix %*% (F_estimate + sampleErrorMatrix) %*% t(Tmatrix))

calcFOmegas_indSweeps_e = function(sel) {
	y = exp(-rec*midDistances/sel*log(4*Ne*sel))
	FOmegas = lapply(1 : length(y), function(i) (calctotAddF_indSweeps_e(y[[i]])))
	return(FOmegas)
}

calctotAddF_indSweeps_e = function(y){
	selMatrix = F_estimate
	for(i  in selPops) {
		selMatrix[i,i] = F_estimate[i,i] + y^2 * (1 - F_estimate[i,i])
	}
	return(Tmatrix %*% (selMatrix + sampleErrorMatrix) %*% t(Tmatrix))
}

calcFOmegas_stdVar_e = function(sel, g, time) {
	y = exp(-rec*midDistances/sel*log(1/g))
	Rf = 4*Ne*rec*midDistances*g
	rt = exp(-rec*midDistances*time)
	FOmegas = lapply(1 : length(midDistances), function(i) calctotAddF_stdVar_e(y[[i]], Rf[[i]], rt[[i]]))
	return(FOmegas)
}

calctotAddF_stdVar_e = function(y, Rf, rt){
	selMatrix = F_estimate
		
	for(i  in selPops) {
		selMatrix[i,i] = (1-y^2) * (F_estimate[i,i]) + y^2 * (1/(1+Rf) +Rf/(1+Rf) * (F_estimate[i,i]))
		for(j  in selPops) {
			if(i != j)
			selMatrix[i,j] = (1-y^2)*F_estimate[i,j] + y^2 * ((rt^2 *(1/(1+Rf) +Rf/(1+Rf) * F_estimate[i,j])) + (1-rt^2)*F_estimate[i,j])
		}
	}
	#return(selMatrix)
	return(Tmatrix %*% (selMatrix + sampleErrorMatrix) %*% t(Tmatrix))
}

calcFOmegas_mig_e = function(sel, mig, source) {
	y = exp(-rec*midDistances/sel*log(4*Ne*sel))
	delta = (1/sel*log(1+sel/mig))/2
	e_delta = exp(-rec*midDistances*delta)
	my.Q = 1 / (1+4*Ne*mig)
	FOmegas = lapply(1 : length(midDistances), function(i) calctotAddF_mig_e(y[[i]], e_delta[[i]], my.Q, source))
	return(FOmegas)
}

calctotAddF_mig_e = function(y, e_delta, my.Q, source){
	selMatrix = F_estimate
	
	if(is.element(source, selPops)) {
		selMatrix[source, source] = (F_estimate[source, source]) + y^2 * (1 - (F_estimate[source, source]))
		
		sisterSource = sisterPops[which(selPops == source)]
		
		for(i in selPops[selPops != source]) {
			selMatrix[i,i] = my.Q * (y^2 + (1-y^2) * (F_estimate[i,i])) + (1-my.Q) * (y^2*e_delta^2 + (1-y)^2 * (F_estimate[i, i])+ 2*y*(1-y)*F_estimate[source, i] + y^2*(1-e_delta^2) * (F_estimate[source, source]))
			selMatrix[i,source] = y^2 * e_delta + (1-y)* F_estimate[source, i] + y*(1-y*e_delta)*(F_estimate[source, source])
			selMatrix[source, i] = y^2 * e_delta + (1-y)* F_estimate[source, i] + y*(1-y*e_delta)*(F_estimate[source, source])

			for(j  in selPops[selPops != source]) {
				if(i != j)
				selMatrix[i,j] = y^2 * e_delta^2 + y^2(1-e_delta)*(F_estimate[source, source]) + (1-y^2)*F_estimate[i,j]	
			}
			selMatrix[i, sisterSource] = (1 - y) * F_estimate[i, sisterSource] + y * F_estimate[source, sisterSource]
			selMatrix[sisterSource, i] = selMatrix[i, sisterSource]
			
			selMatrix[i, sisterPops[which(selPops == i)]] = y * F_estimate[source, sisterPops[which(selPops == i)]] + (1-y) * F_estimate[i, sisterPops[which(selPops == i)]]
			selMatrix[sisterPops[which(selPops == i)], i] = selMatrix[i, sisterPops[which(selPops == i)]]	

		}
	}
	else {
		for(i in selPops) {
			selMatrix[i,i] = my.Q * (y^2 + (1-y^2) * (F_estimate[i,i])) + (1-my.Q) * (y^2*e_delta^2 + (1-y)^2 * (F_estimate[i, i])+ 2*y(1-y)*F_estimate[source, i] + y^2*(1-e_delta^2) * (F_estimate[source, source]))

			for(j  in selPops) {
				if(i != j)
				selMatrix[i,j] = y^2 * e_delta^2 + y^2(1-e_delta)*(F_estimate[source, source]) + (1-y^2)*F_estimate[i,j]	
				#should be f_uu instead of f_ss if s is proxy for unsampled
			}
			selMatrix[i, source] = y * (F_estimate[source, source])
			selMatrix[source, i] = selMatrix[i, source]
			#should be f_us instead of f_ss if s is proxy for unsampled
			
			selMatrix[i, sisterPops[which(selPops == i)]] = (1-y) * F_estimate[i, sisterPops[which(selPops == i)]]
			selMatrix[sisterPops[which(selPops == i)], i] = selMatrix[i, sisterPops[which(selPops == i)]]
		}	
	}
	return(Tmatrix %*% (selMatrix + sampleErrorMatrix) %*% t(Tmatrix))
}


FOmegas_indSweeps_e = lapply(sels, function(sel) calcFOmegas_indSweeps_e(sel))
inv_FOmegas_indSweeps_e = lapply(FOmegas_indSweeps_e, function(i) lapply(i, function(j) ginv(j)))
det_FOmegas_indSweeps_e = lapply(FOmegas_indSweeps_e, function(i) lapply(i, function(j) det(j)))

saveRDS(FOmegas_indSweeps_e, "FOmegas_indSweeps_e.mim.RDS")
saveRDS(inv_FOmegas_indSweeps_e, "inv_FOmegas_indSweeps_e.mim.RDS")
saveRDS(det_FOmegas_indSweeps_e, "det_FOmegas_indSweeps_e.mim.RDS")



FOmegas_stdVar_e = lapply(sels ,function(sel) {
	lapply(gs, function(g) {
		lapply(times, function(time) {
			calcFOmegas_stdVar_e(sel, g, time)
		})
	})
})

inv_FOmegas_stdVar_e = lapply(FOmegas_stdVar_e, function(sel) {
	lapply(sel, function(g) {
		lapply(g, function(time) {
			lapply(time, function(dist) {
				ginv(dist)
			})
		})
	})
})

det_FOmegas_stdVar_e = lapply(FOmegas_stdVar_e, function(sel) {
	lapply(sel, function(g) {
		lapply(g, function(time) {
			lapply(time, function(dist) {
				det(dist)
			})
		})
	})
})

saveRDS(FOmegas_stdVar_e, "FOmegas_stdVar_e.mim.RDS")
saveRDS(inv_FOmegas_stdVar_e, "inv_FOmegas_stdVar_e.mim.RDS")
saveRDS(det_FOmegas_stdVar_e, "det_FOmegas_stdVar_e.mim.RDS")

FOmegas_mig_e = lapply(sels ,function(sel) {
	lapply(migs, function(mig) {
		lapply(sources, function(source) {
			calcFOmegas_mig_e(sel, mig, source)
		})
	})
})

inv_FOmegas_mig_e = lapply(FOmegas_mig_e, function(sel) {
	lapply(sel, function(mig) {
		lapply(mig, function(source) {
			lapply(source, function(dist) {
				ginv(dist)
			})
		})
	})
})

det_FOmegas_mig_e = lapply(FOmegas_mig_e, function(sel) {
	lapply(sel, function(mig) {
		lapply(mig, function(source) {
			lapply(source, function(dist) {
				det(dist)
			})
		})
	})
})

saveRDS(FOmegas_mig_e, "FOmegas_mig_e.mim.RDS")
saveRDS(inv_FOmegas_mig_e, "inv_FOmegas_mig_e.mim.RDS")
saveRDS(det_FOmegas_mig_e, "det_FOmegas_mig_e.mim.RDS")

