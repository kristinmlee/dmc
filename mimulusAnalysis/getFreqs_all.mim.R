copper.mine<-read.table("~/Documents/test_data_for_Kristin/COMS_top10_targetRegions.SNPcount",as.is=TRUE)
copper.mine_adj <-read.table("~/Documents/test_data_for_Kristin/COMS_top10_adjacentRegions.SNPcount",as.is=TRUE)

names(copper.mine)<-c("scaf",	"position",	"RefPop1",	"AltPop1",	"RefPop2",	"AltPop2",	"RefPop3",	"AltPop3", "RefPop4",	"AltPop4")

copper.mine$FreqPop1 <- copper.mine$RefPop1 / (copper.mine$RefPop1 + copper.mine$AltPop1)
copper.mine$FreqPop2 <- copper.mine$RefPop2 / (copper.mine$RefPop2 + copper.mine$AltPop2)
copper.mine$FreqPop3 <- copper.mine$RefPop3 / (copper.mine$RefPop3 + copper.mine$AltPop3)
copper.mine$FreqPop4 <- copper.mine$RefPop4 / (copper.mine$RefPop4 + copper.mine$AltPop4)

names(copper.mine_adj)<-c("scaf",	"position",	"RefPop1",	"AltPop1",	"RefPop2",	"AltPop2",	"RefPop3",	"AltPop3", "RefPop4",	"AltPop4")

copper.mine_adj$FreqPop1 <- copper.mine_adj$RefPop1 / (copper.mine_adj$RefPop1 + copper.mine_adj$AltPop1)
copper.mine_adj$FreqPop2 <- copper.mine_adj$RefPop2 / (copper.mine_adj$RefPop2 + copper.mine_adj$AltPop2)
copper.mine_adj$FreqPop3 <- copper.mine_adj$RefPop3 / (copper.mine_adj$RefPop3 + copper.mine_adj$AltPop3)
copper.mine_adj$FreqPop4 <- copper.mine_adj$RefPop4 / (copper.mine_adj$RefPop4 + copper.mine_adj$AltPop4)

scaf8 = rbind(copper.mine[copper.mine$scaf == "scaffold_8", ], copper.mine_adj[copper.mine_adj$scaf == "scaffold_8", ])
copper.mine_freqs = as.matrix(scaf8[ ,grep("FreqPop", names(scaf8))])

saveRDS(copper.mine_freqs, "mimFreqs_all.RDS")

allPositions = as.matrix(scaf8[ ,2])
saveRDS(allPositions, "mimPositions_all.RDS")