noah_0<-read.table("~/Documents/test_data_for_Kristin/Noah_data/Scaffold0.freq.txt",sep="\t",head=TRUE,as.is=TRUE)
noah_1<-read.table("~/Documents/test_data_for_Kristin/Noah_data/Scaffold1.freq.txt",sep="\t",head=TRUE,as.is=TRUE)
noah_2<-read.table("~/Documents/test_data_for_Kristin/Noah_data/Scaffold2.freq.txt",sep="\t",head=TRUE,as.is=TRUE)
noah_3<-read.table("~/Documents/test_data_for_Kristin/Noah_data/Scaffold3.freq.txt",sep="\t",head=TRUE,as.is=TRUE)

allFreqs = rbind(noah_0[,-1], noah_1[,-1], noah_2[,-1], noah_3[,-1])

samplesPerPop = c(48, 48, 49, 50, 50, 43, 47, 49)*2

dFreq = t(allFreqs)

allRunFreq = apply(dFreq, 2, function(my.freqs) {
	if(runif(1)<0.5){my.freqs<-1-my.freqs};my.freqs
})


numLoci = ncol(allRunFreq)
my.means.rand = (allRunFreq %*% t(allRunFreq)) / numLoci

diag(my.means.rand) = diag(my.means.rand) * samplesPerPop / (samplesPerPop - 1) - rowMeans(allRunFreq)/ (samplesPerPop - 1)

dist.ij = which(my.means.rand == min(my.means.rand), arr.ind = TRUE)[1, ]

A.rand = mean(allRunFreq[dist.ij[1],]*allRunFreq[dist.ij[2],])
C.rand = mean(allRunFreq[dist.ij[1],]*(1-allRunFreq[dist.ij[2],]))

f_estimate.rand = (my.means.rand - A.rand) / C.rand
saveRDS(f_estimate.rand, "F_estimate.killi.RDS")

