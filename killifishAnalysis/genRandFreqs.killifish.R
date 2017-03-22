freqs = read.table("~/Documents/test_data_for_Kristin/Noah_data/Scaffold9893.freq.txt",sep="\t",head=TRUE,as.is=TRUE)

freqs = t(freqs[ ,-1])

allRunFreq = apply(t(freqs), 2, function(my.freqs) {
	if(runif(1)<0.5){my.freqs<-1-my.freqs};my.freqs
})

freqs = t(allRunFreq)

saveRDS(freqs, "rand_freqs.killi.RDS")

