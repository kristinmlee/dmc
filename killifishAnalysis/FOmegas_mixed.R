sets = list(c(2,4,6), 8)
my.modes = c("sv", "ind")

calctotAddF_indSweeps.set_e = function(y, set){
	selMatrix = F_estimate
	temp.selPops = sets[[set]]
	for(i  in temp.selPops) {
		selMatrix[i,i] = F_estimate[i,i] + y^2 * (1 - F_estimate[i,i])
	}
	return(selMatrix - F_estimate)
}

calctotAddF_mig.set_e = function(y, e_delta, my.Q, my.source, set){
	selMatrix = F_estimate
	temp.selPops = sets[[set]]
	temp.nonSelPops = seq(1, numPops)[- temp.selPops]
		
	if(is.element(my.source, temp.selPops)) {
		selMatrix[my.source, my.source] = (F_estimate[my.source, my.source]) + y^2 * (1 - (F_estimate[my.source, my.source]))
				
		for(i in selPops[temp.selPops != my.source]) {
			selMatrix[i,i] = my.Q * (y^2 + (1-y^2) * (F_estimate[i,i])) + (1-my.Q) * (y^2*e_delta^2 + (1-y)^2 * (F_estimate[i, i])+ 2*y*(1-y)*F_estimate[my.source, i] + y^2*(1-e_delta^2) * (F_estimate[my.source, my.source]))
			selMatrix[i, my.source] = y^2 * e_delta + (1-y)* F_estimate[my.source, i] + y*(1-y*e_delta)*(F_estimate[my.source, my.source])
			selMatrix[my.source, i] = y^2 * e_delta + (1-y)* F_estimate[my.source, i] + y*(1-y*e_delta)*(F_estimate[my.source, my.source])
			
			for(k in temp.nonSelPops) {
				selMatrix[k, i] = (1 - y) * F_estimate[i, k] + y * F_estimate[my.source, k]
				selMatrix[i, k] = (1 - y) * F_estimate[i, k] + y * F_estimate[my.source, k]
			}
			

			for(j  in temp.selPops[temp.selPops != my.source]) {
				if(i != j)
				selMatrix[i,j] = y^2 * e_delta^2 + y^2*(1-e_delta^2)*(F_estimate[my.source, my.source]) + (1-y)^2*F_estimate[i,j]	+ (1-y)*y*(F_estimate[i, my.source] + F_estimate[j, my.source])
			}

		}
	}

	return(selMatrix - F_estimate)
}

calctotAddF_stdVar.set_e = function(y, Rf, rt, p_no, p_one, my.source, set) {
	selMatrix = F_estimate
	temp.selPops = sets[[set]]
	temp.nonSelPops = seq(1, numPops)[- temp.selPops]
	
	selMatrix[my.source, my.source] = (1-y^2) * (F_estimate[my.source, my.source]) + y^2 * (1/(1+Rf) +Rf/(1+Rf) * (F_estimate[my.source, my.source]))
	
	for(i in temp.selPops[temp.selPops != my.source]) {
		selMatrix[i,i] = (1-y)^2 * (F_estimate[i, i]) + y^2 * (p_no*(1/(1+Rf) +Rf/(1+Rf)*F_estimate[my.source, my.source]) + (1 - p_no) * (1/(1+Rf)) + ((1 - p_no) * Rf/(1+Rf) - Rf/(1+Rf/2) * (1-p_one) * rt) * F_estimate[i,i] + (1-p_one)*Rf/(1+Rf/2)*rt* F_estimate[i, my.source]) + 2*y*(1-y)*(rt*F_estimate[i, my.source] + (1-rt)*F_estimate[i, i])
			
		selMatrix[i, my.source] = (1-y)^2*F_estimate[i, my.source] + y^2 * (rt^2 *(1/(1+Rf) +Rf/(1+Rf) * F_estimate[my.source, my.source]) + rt*(1-rt)*(F_estimate[my.source, my.source]) + rt*(1-rt)*F_estimate[i, my.source] + (1-rt)^2*F_estimate[i, my.source]) + y*(1-y)*(rt*F_estimate[my.source, my.source] + (1-rt)*F_estimate[i, my.source]) + y*(1-y)*F_estimate[i, my.source]
		selMatrix[my.source, i] = selMatrix[i, my.source]			
		for(k in temp.nonSelPops) {
				selMatrix[k, i] = y * rt * F_estimate[k, my.source] + (1-y) * F_estimate[i, k] + y*(1-rt)*F_estimate[i, k]
				selMatrix[i, k] = y * rt * F_estimate[k, my.source] + (1-y) * F_estimate[i, k] + y*(1-rt)*F_estimate[i, k]
			}
			

			for(j  in temp.selPops[temp.selPops != my.source]) {
				if(i != j)
				selMatrix[i,j] = (1-y)^2*F_estimate[i, j] + y^2 * ((rt^2 *(1/(1+Rf) + Rf/(1+Rf) * F_estimate[my.source, my.source]) + (1-rt)^2 * F_estimate[i, j]) + rt*(1-rt)*(F_estimate[i, my.source] + F_estimate[j, my.source])) + (1-y)*y*(2*(1-rt)*F_estimate[i, j] + rt*(F_estimate[i, my.source] + F_estimate[j, my.source]))
			}

		}

	return(selMatrix - F_estimate)
}

calcFOmegas_mixed_e = function(sel, g, time, mig, my.source, modes) {
	y = exp(-rec*midDistances/sel*log(4*Ne*sel))

	FOmegas = vector("list", length(midDistances))
	FOmegas = lapply(FOmegas, function(i) matrix(0, nrow = numPops, ncol = numPops))
	
	
	for(set in 1 : length(sets)) {
		if(modes[set] == "sv") {
			y = exp(-rec*midDistances/sel*log(1/g))
			Rf = 4*Ne*rec*midDistances*g
			rt = exp(-rec*midDistances*time)
			p_no = exp(-time*(2*rec*midDistances + 1/(2*Ne*g)))
			p_one = exp(-time*(rec*midDistances + 1/(2*Ne*g)))
			FOmegas = lapply(1 : length(midDistances), function(i) calctotAddF_stdVar.set_e(y[[i]], Rf[[i]], rt[[i]], p_no[[i]], p_one[[i]], my.source, set) + FOmegas[[i]])	
		}
		else if(modes[set] == "ind") {
			FOmegas = lapply(1 : length(midDistances), function(i) calctotAddF_indSweeps.set_e(y[[i]], set) + FOmegas[[i]])
		}
		else if(modes[set] == "mig") {
			delta = 1/sel*log(1+sel/(mig))
			e_delta = exp(-rec*midDistances*delta)
			my.Q = 1 / (1+4*Ne*mig)
			FOmegas = lapply(1 : length(midDistances), function(i) calctotAddF_mig.set_e(y[[i]], e_delta[[i]], my.Q, my.source, set) + FOmegas[[i]])
		}
	}
	
	FOmegas = lapply(FOmegas, function(i) Tmatrix %*% (i + F_estimate + sampleMatrix) %*% t(Tmatrix))


	return(FOmegas)
}

calcFOmegas_mixed_e = function(sel, g, time, mig, source, modes) {
	y = exp(-rec*midDistances/sel*log(4*Ne*sel))

	FOmegas = vector("list", length(midDistances))
	FOmegas = lapply(FOmegas, function(i) matrix(0, nrow = numPops, ncol = numPops))
	
	
	for(set in 1 : length(sets)) {
		if(modes[set] == "sv") {
			y = exp(-rec*midDistances/sel*log(1/g*sel))
			Rf = 4*Ne*rec*midDistances*g
			rt = exp(-rec*midDistances*time)
			FOmegas = lapply(1 : length(midDistances), function(i) calctotAddF_stdVar.set_e(y[[i]], Rf[[i]], rt[[i]], set) + FOmegas[[i]])	
		}
		else if(modes[set] == "ind") {
			FOmegas = lapply(1 : length(midDistances), function(i) calctotAddF_indSweeps.set_e(y[[i]], set) + FOmegas[[i]])
		}
		else if(modes[set] == "mig") {
			delta = 1/sel*log(1+sel/(mig))
			e_delta = exp(-rec*midDistances*delta)
			my.Q = 1 / (1+4*Ne*mig)
			FOmegas = lapply(1 : length(midDistances), function(i) calctotAddF_mig.set_e(y[[i]], e_delta[[i]], my.Q, source, set) + FOmegas[[i]])
		}
	}
	
	FOmegas = lapply(FOmegas, function(i) Tmatrix %*% (i + F_estimate + sampleMatrix) %*% t(Tmatrix))


	return(FOmegas)
}

my.modes_svInd = c("sv", "ind")
FOmegas_mixed.svInd = lapply(sels ,function(sel) {
	lapply(gs, function(g) {
		lapply(times, function(time) {
			lapply(migs[1], function(mig) {
				lapply(sources, function(my.source) {
					calcFOmegas_mixed_e(sel, g, time, mig, my.source, my.modes_svInd)
				})
			})
		})
	})
})


my.modes_migInd = c("mig", "ind")
FOmegas_mixed.svInd = lapply(sels ,function(sel) {
	lapply(gs[1], function(g) {
		lapply(times[1], function(time) {
			lapply(migs, function(mig) {
				lapply(sources, function(my.source) {
					calcFOmegas_mixed_e(sel, g, time, mig, my.source, my.modes_migInd)
				})
			})
		})
	})
})
