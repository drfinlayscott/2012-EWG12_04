library(FLAdvice)
library(FLBioDym)
library(plyr)
source("../funs.R")
scn <- rbind(
	expand.grid(Btrig=0.5, Ftar=0.75, maxHR=0.35, aLag=c(1,3,5), srvBias=c(0.5, 1, 1.5), cthBias=c(0.5,1), IEM="linear", clag=c(5,10,20))
	)
load("RData.MSEsims01cth")
MSEcth.summ <- mseSumm(MSEsims01cth, scn)
save(MSEcth.summ, file="../RData.MSEcth.summ")


