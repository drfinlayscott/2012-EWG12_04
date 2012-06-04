library(FLAdvice)
library(FLBioDym)
library(plyr)
source("../funs.R")
scn <- rbind(
	expand.grid(Btrig=0.5, Ftar=0.75, maxHR=0.35, aLag=c(1,3,5), srvBias=c(0.5, 1, 1.5), cthBias=c(0.5,1), IEM="linear", slag=c(5,10,20), clag=c(5,10,20))
	)
load("RData.MSEsims01srv")
MSEsrv.summ <- mseSumm(MSEsims01srv, scn)
save(MSEsrv.summ, scn, file="../RData.MSEsrv.summ")


