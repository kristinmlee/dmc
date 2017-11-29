# Scripts to generate inverses and determinants of variance/covariance matrices
#	with any model of convergent selection (F^(S)): neutral model (F), independent mutations,
#	standing variant model with source, standing variant model without source, migration
# 	model, combinations of convergent selection models (mixed)
#
# Args:
# numPops: number of populations sampled (both selected and non-selected)
#	sampleSizes: vector of sample sizes of length numberOfPopulations
#		(i.e. twice the number of diploid individuals sampled in each population)
#		Note: needed only for neutral model since other models are already sample-size
#		corrected
#	*variance/covariance matrices to apply operations to*
#	F_estimate: estimate of neutral variance/covariance matrix
#		generated with "calcNeutralF.R" (for neutral model)
#	FOmegas_indSweeps: list of length numBins of variance/covariance matrices with convergent
#		adapatation due to independent mutations model
#	FOmegas_stdVar.source: list of length numBins of variance/covariance matrices with convergent
#		adapatation due to standing variant with source model
#	FOmegas_stdVar: list of length numBins of variance/covariance matrices with convergent
#		adapatation due to standing variant without source model
#	FOmegas_mig: list of length numBins of variance/covariance matrices with convergent
#		adapatation due to migration model
#	FOmegas_mixed: 	list of length numBins of variance/covariance matrices with convergent
#		adapatation due to multiple specified modes of convergent adaptation


##Neutral model
sampleErrorMatrix = diag(1/sampleSizes, nrow = numPops, ncol = numPops)
det_FOmegas_neutral = det(Tmatrix %*% (F_estimate + sampleErrorMatrix) %*% t(Tmatrix))
inv_FOmegas_neutral = ginv(Tmatrix %*% (F_estimate + sampleErrorMatrix) %*% t(Tmatrix))

##Independent mutations model
det_FOmegas_indSweeps = lapply(FOmegas_indSweeps, function(sel) {
	lapply(sel, function(dist) {
		det(dist)
	})
})
inv_FOmegas_indSweeps = lapply(FOmegas_indSweeps, function(sel) {
	lapply(sel, function(dist) {
		ginv(dist)
	})
})

##Standing variant with source model
det_FOmegas_stdVar.source = lapply(FOmegas_stdVar.source, function(sel) {
	lapply(sel, function(g) {
		lapply(g, function(time) {
			lapply(time, function(my.source) {
				lapply(my.source, function(dist) {
					det(dist)
				})
			})
		})
	})
})
inv_FOmegas_stdVar.source = lapply(FOmegas_stdVar.source, function(sel) {
	lapply(sel, function(g) {
		lapply(g, function(time) {
			lapply(time, function(my.source) {
				lapply(my.source, function(dist) {
					ginv(dist)
				})
			})
		})
	})
})

##Standing variant without source model
det_FOmegas_stdVar = lapply(FOmegas_stdVar, function(sel) {
	lapply(sel, function(g) {
		lapply(g, function(time) {
			lapply(time, function(dist) {
				det(dist)
			})
		})
	})
})
inv_FOmegas_stdVar = lapply(FOmegas_stdVar, function(sel) {
	lapply(sel, function(g) {
		lapply(g, function(time) {
			lapply(time, function(dist) {
				ginv(dist)
			})
		})
	})
})

##Migration model
det_FOmegas_mig = lapply(FOmegas_mig, function(sel) {
	lapply(sel, function(mig) {
		lapply(mig, function(source) {
			lapply(source, function(dist) {
				det(dist)
			})
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

##Multiple convergent modes specified model
detFOmegas_mixed = lapply(FOmegas_mixed, function(sel) {
	lapply(sel, function(g) {
		lapply(g, function(time) {
			lapply(time, function(mig) {
				lapply(mig, function(source) {
					lapply(source, function(dist) {
						det(dist)
					})
				})	
			})
		})
	})
})
invFOmegas_mixed = lapply(FOmegas_mixed, function(sel) {
	lapply(sel, function(g) {
		lapply(g, function(time) {
			lapply(time, function(mig) {
				lapply(mig, function(source) {
					lapply(source, function(dist) {
						ginv(dist)
					})
				})	
			})
		})
	})
})