# Calculate inverses and determinants of selection matrices generated with genSelMatrices_exec.R
# Remember that there is a selection matrix for each parameter + distance bin combination for each model
# These results are used in calcCompLike_exec.R

library(MASS) #to get ginv() function

## Load FOmegas saved in genSelMatrices_exec.R
FOmegas_stdVar= readRDS("FOmegas_stdVar.RDS") #pars: sels, gs, std var times
FOmegas_stdVar.source= readRDS("FOmegas_stdVar_source.RDS") #pars: sels, gs, times
FOmegas_mig.stagSweeps=readRDS("FOmegas_mig_stagSweeps.RDS") # pars: sels, gs, times
FOmegas_mig.concSweeps=readRDS("FOmegas_mig_concSweeps.RDS") # pars: sels, migs

## Neutral model
sampleSizes=readRDS("sampleSizes.RDS")
F_estimate=readRDS("neutralF.RDS")
numPops=ncol(F_estimate)
M = numPops
Tmatrix = matrix(data = rep(-1 / M, (M - 1) * M), nrow = M - 1, ncol = M)
diag(Tmatrix) = (M - 1) / M 
sampleErrorMatrix = diag(1/sampleSizes, nrow = numPops, ncol = numPops)
det_FOmegas_neutral = det(Tmatrix %*% (F_estimate + sampleErrorMatrix) %*% t(Tmatrix))
inv_FOmegas_neutral = ginv(Tmatrix %*% (F_estimate + sampleErrorMatrix) %*% t(Tmatrix))
saveRDS(det_FOmegas_neutral, file="det_FOmegas_neutral_AHR.RDS")
saveRDS(inv_FOmegas_neutral, file="inv_FOmegas_neutral_AHR.RDS")

## Standing Variant Model (ILS)
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
saveRDS(det_FOmegas_stdVar, file="det_FOmegas_stdVar_AHR.RDS")
saveRDS(inv_FOmegas_stdVar, file="inv_FOmegas_stdVar_AHR.RDS")

## Standing Variant Source Model
det_FOmegas_stdVar.source = lapply(FOmegas_stdVar.source, function(sel) {
  lapply(sel, function(g) {
    lapply(g, function(time) {
      lapply(time, function(dist) {
        det(dist)
      })
    })
  })
})
inv_FOmegas_stdVar.source = lapply(FOmegas_stdVar.source, function(sel) {
  lapply(sel, function(g) {
    lapply(g, function(time) {
      lapply(time, function(dist) {
        ginv(dist)
      })
    })
  })
})
saveRDS(det_FOmegas_stdVar.source, file="det_FOmegas_stdVar_source_AHR.RDS")
saveRDS(inv_FOmegas_stdVar.source, file="inv_FOmegas_stdVar_source_AHR.RDS")


## Staggered Sweeps Model
det_FOmegas_mig.stagSweeps = lapply(FOmegas_mig.stagSweeps, function(sel){
  lapply(sel, function(g){
    lapply(g, function(time){
      lapply(time, function(dist){
        det(dist)
      })
    })
  })
})
saveRDS(det_FOmegas_mig.stagSweeps, file="det_FOmegas_mig_stagSweeps_AHR.RDS")

inv_FOmegas_mig.stagSweeps = lapply(FOmegas_mig.stagSweeps, function(sel){
  lapply(sel, function(g){
    lapply(g, function(time){
      lapply(time, function(dist){
        ginv(dist)
      })
    })
  })
})
saveRDS(inv_FOmegas_mig.stagSweeps, file="inv_FOmegas_mig_stagSweeps_AHR.RDS")

## Concurrent Sweeps Model
det_FOmegas_mig.concSweeps = lapply(FOmegas_mig.concSweeps, function(sel) {
  lapply(sel, function(mig) {
    lapply(mig, function(dist) {
      det(dist)
    })
  })
})
inv_FOmegas_mig.concSweeps = lapply(FOmegas_mig.concSweeps, function(sel) {
  lapply(sel, function(mig) {
    lapply(mig, function(dist) {
      ginv(dist)
    })
  })
})
saveRDS(det_FOmegas_mig.concSweeps, file="det_FOmegas_mig_concSweeps_AHR.RDS")
saveRDS(inv_FOmegas_mig.concSweeps, file="inv_FOmegas_mig_concSweeps_AHR.RDS")
