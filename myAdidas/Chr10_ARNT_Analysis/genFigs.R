#### Load results from each model ####

#neutral is a list of length selsite
cl_neutral = readRDS("compLikelihood_neutral_ARNT.RDS")[1]

#standing var without source
# three parameters: sels, gs, std var times
cl_sv = readRDS("compLikelihood_stdVar_ARNT.RDS")

#standing var with source
# three parameters: sels, gs, times
cl_svSource = readRDS("compLikelihood_stdVar_source_ARNT.RDS")

#mig with staggered sweeps (because of standing time)
# three parameters: sels, gs, times
cl_mig_stagSweeps = readRDS("compLikelihood_mig_stagSweeps_ARNT.RDS")

#mig with concurrent / immediate sweeps
# two parameters: sels, migs
cl_mig_concSweeps = readRDS("compLikelihood_mig_concSweeps_ARNT.RDS")

# load parameters over which likelihoods were calculated
selSite=readRDS("selSite.RDS")
gs = readRDS("gs.RDS")
migs = readRDS("migs.RDS")
sels = readRDS("sels.RDS")
stdVar_times = readRDS("stdVar_times.RDS")
times = readRDS("times.RDS")

#### Identify Profile Composite Log-Likelihood Surface ####

## Composite log-likelihood surface for selected site
maxCL_sv.selSite =  sapply(1:length(cl_sv), function(i) max(unlist(cl_sv[[i]]))  )
maxCL_svSource.selSite =  sapply(1:length(cl_svSource), function(i) max(unlist(cl_svSource[[i]]),na.rm=T)  )
maxCL_mig_stagSweeps.selSite =  sapply(1:length(cl_mig_stagSweeps), function(i) max(unlist(cl_mig_stagSweeps[[i]]))  )
maxCL_mig_concSweeps.selSite =  sapply(1:length(cl_mig_concSweeps), function(i) max(unlist(cl_mig_concSweeps[[i]]))  )

## Composite log-likelihood surface for selection coefficient
cl_sv_bySel =  lapply(1:length(sels), function(sel)  lapply(1:length(selSite), function(i) lapply(1:length(gs) ,function(j)  sapply(1:length(stdVar_times), function(k) max(unlist(cl_sv[[i]][[sel]][[j]][[k]]))) ) ))
cl_svSource_bySel =  lapply(1:length(sels), function(sel)  lapply(1:length(selSite), function(i) lapply(1:length(gs) ,function(j)  sapply(1:length(times), function(k) max(unlist(cl_svSource[[i]][[sel]][[j]][[k]]))) ) ))
cl_mig_stagSweeps_bySel =  lapply(1:length(sels), function(sel)  lapply(1:length(selSite), function(i) lapply(1:length(gs) ,function(j)  sapply(1:length(times), function(k) max(unlist(cl_mig_stagSweeps[[i]][[sel]][[j]][[k]]))) ) ))
cl_mig_concSweeps_bySel =  lapply(1:length(sels), function(sel)  lapply(1:length(selSite), function(i) lapply(1:length(migs) ,function(j) max(unlist(cl_mig_concSweeps[[i]][[sel]][[j]])) ) ))

maxCL_sv.sels = sapply(1:length(sels), function(s) max(unlist(cl_sv_bySel[[s]])))
maxCL_svSource.sels = sapply(1:length(sels), function(s) max(unlist(cl_svSource_bySel[[s]])))
maxCL_mig_stagSweeps.sels = sapply(1:length(sels), function(s) max(unlist(cl_mig_stagSweeps_bySel[[s]])))
maxCL_mig_concSweeps.sels = sapply(1:length(sels), function(s) max(unlist(cl_mig_concSweeps_bySel[[s]])))

## Composite log-likelihood surface for amount of time the variant was standing after migration and before selection begins.
cl_mig_stagSweeps_byTime = lapply(1 : length(times), function(time) lapply(1:length(selSite), function(site) lapply(1 : length(sels), function(sel) lapply(1 : length(gs), function(g) cl_mig_stagSweeps[[site]][[sel]][[g]][[time]]))) )
cl_svSource_byTime = lapply(1 : length(times), function(time) lapply(1:length(selSite), function(site)  lapply(1 : length(sels), function(sel) lapply(1 : length(gs), function(g) cl_svSource[[site]][[sel]][[g]][[time]]))))

maxCL_mig_stagSweeps.times = sapply(1 : length(times), function(i) max(unlist(cl_mig_stagSweeps_byTime[[i]])))
maxCL_svSource.times = sapply(1 : length(times), function(i) max(unlist(cl_svSource_byTime[[i]])))

## Composite log-likelihood surface for amount of time the variant was standing after selected populations split and before selection begins.
cl_sv_byTime = lapply(1 : length(stdVar_times), function(time) lapply(1:length(selSite), function(site)  lapply(1 : length(sels), function(sel) lapply(1 : length(gs), function(g) cl_sv[[site]][[sel]][[g]][[time]]))))
maxCL_sv.times = sapply(1 : length(stdVar_times), function(i) max(unlist(cl_svSource_byTime[[i]])))

## Composite log-likelihood surface for the frequency of the selected variant when standing
cl_mig_stagSweeps_byG = lapply(1 : length(gs), function(g) lapply(1:length(selSite), function(site) lapply(1 : length(sels), function(sel) lapply(1 : length(times), function(time) cl_mig_stagSweeps[[site]][[sel]][[g]][[time]]))))
cl_svSource_byG = lapply(1 : length(gs), function(g) lapply(1:length(selSite), function(site) lapply(1 : length(sels), function(sel) lapply(1 : length(times), function(time) cl_svSource[[site]][[sel]][[g]][[time]]))))
cl_sv_byG = lapply(1 : length(gs), function(g) lapply(1:length(selSite), function(site) lapply(1 : length(sels), function(sel) lapply(1 : length(stdVar_times), function(time) cl_sv[[site]][[sel]][[g]][[time]]))))

maxCL_mig_stagSweeps.gs = sapply(1 : length(gs), function(i) max(unlist(cl_mig_stagSweeps_byG[[i]])))
maxCL_svSource.gs = sapply(1 : length(gs), function(i) max(unlist(cl_svSource_byG[[i]])))
maxCL_sv.gs = sapply(1 : length(gs), function(i) max(unlist(cl_sv_byG[[i]])))

## Composite log-likelihood surface for the migration rate
cl_mig_concSweeps_byMig =  lapply(1:length(migs), function(m)  lapply(1:length(selSite), function(i) lapply(1:length(sels) ,function(j) max(unlist(cl_mig_concSweeps[[i]][[j]][[m]])) ) ))
maxCL_mig_concSweeps.migs = sapply(1:length(migs), function(m) max(unlist(cl_mig_concSweeps_byMig[[m]])))

#### Generate Plots ####

## plotting colors for each model
sv = "midnightblue"
src = "aquamarine4"
stag = "chocolate2"
conc = "violetred4"

pdf("compLike_grandis_selSite.pdf", width = 7, height = 5)
par(oma = c(0.8,1,0.5,0.3))
yMin_site = min(c(maxCL_mig_stagSweeps.selSite, maxCL_svSource.selSite, maxCL_sv.selSite, maxCL_mig_concSweeps.selSite) - cl_neutral)
yMax_site = max(c(maxCL_mig_stagSweeps.selSite, maxCL_svSource.selSite, maxCL_sv.selSite, maxCL_mig_concSweeps.selSite) - cl_neutral)
plot(selSite/1e6, maxCL_svSource.selSite - cl_neutral, ylab = "Composite log-likelihood (model - neutral)", ylim = c(yMin_site, yMax_site), type = "l", col = src, xlim = c(min(selSite)/1e6, max(selSite)/1e6), xlab = "position of selected site (Mbp)", lwd = 3, cex.lab = 1.1, xaxt="n")
lines(selSite/1e6, maxCL_sv.selSite - cl_neutral, col = sv, lwd = 3, lty = 2)
lines(selSite/1e6, maxCL_mig_concSweeps.selSite - cl_neutral, col = conc, lwd = 3, lty = 1)
lines(selSite/1e6, maxCL_mig_stagSweeps.selSite - cl_neutral, col = stag, lwd = 3, lty = 2)
axis(side=1,at=seq(18,34,by=2))
legend("topright", lty = c(1,2, 1, 2), legend = c("Introgression (sv source)", "Introgression (staggered sweeps)", "Introgression (concurrent sweeps)", "ILS"), lwd = 3, col = c(src,stag,conc,sv), cex=0.73)
dev.off() 

pdf("compLike_grandis_sels_all.pdf", width = 6, height = 5)
par(oma = c(0.8,1,0.5,0.3))
yMin_all = min(c(maxCL_mig_stagSweeps.sels, maxCL_svSource.sels, maxCL_sv.sels, maxCL_mig_concSweeps.sels) - cl_neutral)
yMax_all = max(c(maxCL_mig_stagSweeps.sels, maxCL_svSource.sels, maxCL_sv.sels, maxCL_mig_concSweeps.sels) - cl_neutral)
plot(sels, maxCL_svSource.sels - cl_neutral, ylab = "Composite log-likelihood (model - neutral)", ylim = c(yMin_all, yMax_all), type = "l", col = src, xlim = c(0, 1), xlab = "selection strength", lwd = 3, cex.lab = 1.1)
lines(sels, maxCL_sv.sels - cl_neutral, col = sv, lwd = 3, lty = 2)
lines(sels, maxCL_mig_concSweeps.sels - cl_neutral, col = conc, lwd = 3, lty = 1)
lines(sels, maxCL_mig_stagSweeps.sels - cl_neutral, col = stag, lwd = 3, lty = 2)
legend("bottomleft", lty = c(1,2, 1, 2), legend = c("Introgression (sv source)", "Introgression (staggered sweeps)", "Introgression (concurrent sweeps)", "ILS"), lwd = 3, col = c(src,stag,conc,sv))
dev.off() 

pdf("compLike_grandis_sels_topTwo.pdf", width = 6, height = 5)
par(oma = c(0.8,1,0.5,0.3))
yMin_two = min(c(maxCL_mig_stagSweeps.sels[7:13], maxCL_svSource.sels[7:13]) - cl_neutral)
yMax_two = max(c(maxCL_mig_stagSweeps.sels[7:13], maxCL_svSource.sels[7:13]) - cl_neutral)
plot(sels[7:13], maxCL_svSource.sels[7:13] - cl_neutral, ylab = "Composite log-likelihood (model - neutral)", ylim = c(yMin_two, yMax_two), type = "l", col = src, xlab = "selection strength", lwd = 3, cex.lab = 1.1)
lines(sels[7:13], maxCL_mig_stagSweeps.sels[7:13] - cl_neutral, col = stag, lwd = 3, lty = 2)
legend("bottomright", lty = c(1,2), legend = c("Introgression (sv source)", "Introgression (staggered sweeps)"), lwd = 3, col = c(src, stag))
dev.off() 

pdf("compLike_grandis_times.pdf", width = 6, height = 5)
par(oma = c(0.8,1,0.5,0.3))
plot(times[1:20], maxCL_svSource.times[1:20] - cl_neutral, ylab = "Composite log-likelihood (model - neutral)", type = "l", col = src, xlab = "time between introgression and selection (generations)", lwd = 3, cex.lab = 1.1)
lines(times[1:20], maxCL_mig_stagSweeps.times[1:20] - cl_neutral, col = stag, lwd = 3, lty = 2)
legend("topright", lty = c(1,2), legend = c("Introgression (sv source)", "Introgression (staggered sweeps)"), lwd = 3, col = c(src,stag))
dev.off() 

pdf("compLike_grandis_timeSV.pdf", width = 6, height = 5)
par(oma = c(0.8,1,0.5,0.3))
plot(log10(stdVar_times), maxCL_sv.times - cl_neutral, ylab = "Composite log-likelihood (model - neutral)", type = "l", col = sv, lty = 2, xlab = expression(paste(log[10], "(standing time (generations))")), lwd = 3, cex.lab = 1.1)
legend(x="topright",lty=2,col=sv,legend="ILS",lwd=3)
dev.off() 

pdf("compLike_grandis_g.pdf", width = 6, height = 5)
yMin_3 = min(c(maxCL_mig_stagSweeps.sels, maxCL_svSource.sels, maxCL_sv.sels) - cl_neutral)
yMax_3 = max(c(maxCL_mig_stagSweeps.sels, maxCL_svSource.sels, maxCL_sv.sels) - cl_neutral)
par(oma = c(0.8,1,0.5,0.3))
plot(log10(gs), maxCL_sv.gs - cl_neutral, ylim = c(yMin_3-25000, yMax_3), ylab = "Composite log-likelihood (model - neutral)", type = "l", col = sv, lty = 2, xlab = expression(paste(log[10], "(frequency of variant prior to selection)")), lwd = 3, cex.lab = 1.1)
lines(log10(gs), maxCL_svSource.gs - cl_neutral, col = src, lwd = 3)
lines(log10(gs), maxCL_mig_stagSweeps.gs - cl_neutral, col = stag, lty = 2, lwd = 3)
legend("bottomleft", lty = c(1,2, 2), legend = c("Introgression (sv source)", "Introgression (staggered sweeps)", "ILS"), lwd = 3, col = c(src,stag, sv))
dev.off() 



