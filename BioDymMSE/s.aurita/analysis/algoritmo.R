#####################################################################
# first runs for testing
# pelat: B[t+1] = B[t] + r/p*B[t]*(1-(B[t]/k)^p)) - C[t]
#--------------------------------------------------------------------
# year definitions: 
#	aLag - assessment lag in years (user set)
#	lastPyr - last year for projections (user set)
#	iniPyr - first year for projections (user set)
#	nPyr - number of year to project (computed)
#	aYr - assessment year (computed)
#	iYr - i year in the loop, it coincides with the intermediate year
#			in the usual ICES settings (loop) 
#	dtYr - last year with data which is the last year for which there 
#			are estimates (loop)
#	advYr - advice year, the year for which advice is being given (loop)     
#--------------------------------------------------------------------
#####################################################################

library(FLAdvice)
library(FLBioDym)
library(plyr)
#library(multicore)
#library(foreach)
#library(doMC)
#registerDoMC(3)
source("funs.R")

#====================================================================
# Simulation settings
#====================================================================

nits <- 5				# number of iterations
iniyr <- 2011 			# first year in projections
lastyr <- 2061 			# last year in projections
npyr <- lastyr-iniyr+1 	# number of years to project
srsd <- 0.2 			# sd for S/R

#====================================================================
# Rebuild the stock history with simulated data
#====================================================================

#--------------------------------------------------------------------
# Sardinella data
#--------------------------------------------------------------------
nc <- read.table('../data/nc_sar.dat', header=TRUE)
ia <- read.table('../data/cpue_sar.dat', header=TRUE)
ca <- FLQuant(nc$catch, dimnames=list(age='all', year=nc$year))
cp <- FLQuant(c(rep(NA, 7), ia$cpue), dimnames=list(age='all', year=nc$year))

#--------------------------------------------------------------------
# gislasim and projection
#--------------------------------------------------------------------
# Use CECAF SA B0 as starting point
par <- as(data.frame(linf=28.5, k=0.8, t0=-0.1, s=0.8, v=1750, a50=1), 'FLPar')
sar <- lh(gislasim(par), range=c(min=1, max=8, minfbar=1, maxfbar=6))

# stk with initial F closest to estimated HR
saa <- as(sar, 'FLStock')
saa <- saa[,7]
dimnames(saa) <- list(year=1989)

# prepare for projection
saa <- window(saa, FLBRP=sar, end=iniyr-1)

# projection control
trg <- fwdControl(data.frame(year=1990:2010, val=c(ca), quantity="catch"))

# catch projection
saa <- fwd(saa, trg, sr=list(model="bevholt", params=params(sar)))

name(saa) <- "SAA"
desc(saa) <- "Simulated Sardinella aurita, partly conditioned on CECAF SA 2010"

# CECAF SA results
sa <- read.table('../data/sa.dat', sep='\t', header=T)

par(mfrow=c(2,1))
plot(sa$hr~sa$year, type="l", ylim=c(0,1), main="Harvest Rate")
lines(c(catch(saa)/stock(saa))[-1]~sa$year, col=2)
plot(sa$ssb~sa$year, type="l", ylim=c(0,2000), main="SSB")
lines(c(ssb(saa))[-1]~sa$year, col=2)


#====================================================================
# Conditioning
#====================================================================

om <- saa
srBH <- as.FLSR(om,model="bevholt")
params(srBH) <- params(sar)

# Residuals - simulate residuals as lognormal with sd=srsd
set.seed(123)
srRsdl <- FLQuant(rlnorm(npyr*nits, 0, srsd), dimnames=list(year=srBH@range["minyear"]:lastyr, iter=1:nits)) 

#--------------------------------------------------------------------
# create OM object
# Note: this object is projected at Fsq and the values for the first
#	intermediate year are used in the projections
#--------------------------------------------------------------------

# window with FLBRP expands the object including weights, etc, the 
# brp doesn't seem to do anything except dispatching. it replaces "stf". 
OM <- window(om, FLBRP=sar, end=lastyr)

# trick to get iterations, start with M and fwd will add to other slots
m(OM) <- propagate(m(OM),nits)

# project to the end of projections at last year F level
ctrl <- fwdControl(data.frame(year=iniyr:lastyr, quantity="f", val=c(fbar(om)[,ac(om@range["maxyear"])])*rep(1,npyr)))
OM <- fwd(OM, ctrl=ctrl, sr=srBH, sr.residuals=srRsdl)
plot(OM)

#====================================================================
# MP settings
#====================================================================

# these bounds are related with ADMB adjusted to CECAF usage
bounds <- bounds(FLBioDym())
bounds["sigma","start"]=0.50
bounds["q",    "start"]=1.0
bounds["q",    1]      =1.0
bounds["b0",   c("phase","start")]=c(-1,0.2)
bounds[,"lower"]=bounds[,"start"]*0.1
bounds[,"upper"]=bounds[,"start"]*10.0
bounds["p",    c("phase","lower","upper","start")] = c(1,1,1.1,1)
bounds["r",    c("phase","lower","upper","start")] = c(1, 0.32, 1.28, 0.64)
bounds["K",    c("phase","lower","upper","start")] = c(2, 1750, 7000, 2500)

#====================================================================
# Scenarios and simulations
#====================================================================

# scenarios #1
scn01 <- expand.grid(Btrig=0.5, Ftar=1, aLag=c(1), srvBias=c(1), cthBias=c(1), IEM=c("linear", "cst")

MSEsims00 <- mseBD(OM=OM, start=iniyr, sr=srBH, srRsdl=srRsdl, bounds=bounds, CV=0.2, Btrig=0.5, Ftar=1, aLag=1, srvBias=1, cthBias=1, IEM="linear")

MSEsims00fixed <- mseBD(OM=OM, start=iniyr, sr=srBH, srRsdl=srRsdl, bounds=bounds, CV=0.2, Btrig=0.5, Ftar=1, aLag=1, srvBias=1.5, cthBias=1, IEM="linear")


save(MSEsims01, file="RData.MSEsims01")

# scenarios #2
scn02 <- expand.grid(Btrig=0.5, CV=0.2, Ftar=1, aLag=c(1,3,5), srvBias=c(0.5, 1, 1.5), cthBias=c(0.5, 1), IEM="cst")

MSEsims02 <- split(scn02, scn02)
for(i in 1:nrow(scn02)){
	MSEsims02[[i]] <- try(mseBD(OM=OM, start=iniyr, sr=srBH, srRsdl=srRsdl, bounds=bounds, CV=scn[i,"CV"], Btrig=scn[i,"Btrig"], Ftar=scn[i,"Ftar"], aLag=scn[i, "aLag"], srvBias=scn[i,"srvBias"], cthBias=scn[i,"cthBias"], IEM=scn[i,"IEM"], seed=22))
}

save(MSEsims02, file="RData.MSEsims02")

# scenarios #3
scn03 <- expand.grid(Btrig=0.5, CV=0.2, Ftar=1.5, aLag=c(1,3,5), srvBias=c(0.5, 1, 1.5), cthBias=c(0.5, 1), IEM="linear")

MSEsims03 <- split(scn03, scn03)
for(i in 1:nrow(scn03)){
	MSEsims03[[i]] <- try(mseBD(OM=OM, start=iniyr, sr=srBH, srRsdl=srRsdl, bounds=bounds, CV=scn[i,"CV"], Btrig=scn[i,"Btrig"], Ftar=scn[i,"Ftar"], aLag=scn[i, "aLag"], srvBias=scn[i,"srvBias"], cthBias=scn[i,"cthBias"], IEM=scn[i,"IEM"], seed=22))
}

save(MSEsims03, file="RData.MSEsims03")

# scenarios #4
scn04 <- expand.grid(Btrig=0.5, CV=0.2, Ftar=1.5, aLag=c(1,3,5), srvBias=c(0.5, 1, 1.5), cthBias=c(0.5, 1), IEM="cst")

MSEsims04 <- split(scn04, scn04)
for(i in 1:nrow(scn04)){
	MSEsims04[[i]] <- try(mseBD(OM=OM, start=iniyr, sr=srBH, srRsdl=srRsdl, bounds=bounds, CV=scn[i,"CV"], Btrig=scn[i,"Btrig"], Ftar=scn[i,"Ftar"], aLag=scn[i, "aLag"], srvBias=scn[i,"srvBias"], cthBias=scn[i,"cthBias"], IEM=scn[i,"IEM"], seed=22))
}

save(MSEsims04, file="RData.MSEsims04")

# scenarios #5 b0=0.5
scn <- expand.grid(Btrig=0.5, CV=0.2, Ftar=1, aLag=c(1,3,5), srvBias=c(0.5, 1, 1.5), cthBias=c(0.5, 1), IEM="linear")

MSEsims05 <- split(scn, scn)
for(i in 1:nrow(scn)){
	MSEsims05[[i]] <- try(mseBD(OM=OM, start=iniyr, sr=srBH, srRsdl=srRsdl, bounds=bounds, CV=scn[i,"CV"], Btrig=scn[i,"Btrig"], Ftar=scn[i,"Ftar"], aLag=scn[i, "aLag"], srvBias=scn[i,"srvBias"], cthBias=scn[i,"cthBias"], IEM=scn[i,"IEM"], seed=22))
}

# scenarios #6 p=1
scn <- expand.grid(Btrig=0.5, CV=0.2, Ftar=1, aLag=c(1,3,5), srvBias=c(0.5, 1, 1.5), cthBias=c(0.5, 1), IEM="linear")

MSEsims06 <- split(scn, scn)
for(i in 1:nrow(scn)){
	MSEsims06[[i]] <- try(mseBD(OM=OM, start=iniyr, sr=srBH, srRsdl=srRsdl, bounds=bounds, CV=scn[i,"CV"], Btrig=scn[i,"Btrig"], Ftar=scn[i,"Ftar"], aLag=scn[i, "aLag"], srvBias=scn[i,"srvBias"], cthBias=scn[i,"cthBias"], IEM=scn[i,"IEM"], seed=22))
}

# paralel
res0 <- foreach(i = 1:nrow(scn), .packages=c("FLCore", "FLash", "FLEDA", "MASS", "FLBioDym", "FLAdvice", "plyr"), .export="an") %dopar% mseBD(OM=OM, start=iniyr, sr=srBHAR1, srRsdl=srRsdlAR, bounds=bounds, CV=scn[i,"CV"], Btrig=scn[i,"Btrig"], Ftar=scn[i,"Ftar"], aLag=scn[i, "aLag"], srvBias=scn[i,"srvBias"], cthBias=scn[i,"cthBias"], IEM=scn[i,"IEM"])

#====================================================================
# Analysis of results
#====================================================================

#--------------------------------------------------------------------
# Read results and process
#--------------------------------------------------------------------

load("out01/RData.MSEsims01")
MSEsumm01 <- mseSumm(MSEsims01, scn01)


xyplot(data~factor(year)|par, groups=qtl, data=subset(MSEsumm01, runid==15), scales=list(y=list(relation="free")), par.settings=list(superpose.line=list(lty=c(2,1,2), col=c(2,1,2))), type="l")

xyplot(data~factor(year), groups=par, data=subset(MSEsumm01, runid==13 & par %in% c("TAC", "catch") & qtl==0.5), type="l")

xyplot(data~factor(year)|scn, groups=par, data=subset(MSEsumm01, runid<19 & par %in% c("TAC", "catch") & qtl==0.5 & year > 2000), type="l", par.strip.text=list(cex=0.5), layout=c(6,3))




> scn01 <- expand.grid(Btrig=0.5, CV=0.2, Ftar=1, aLag=c(1), srvBias=c(1), cthBias=c(1), IEM="linear")
> mseSumm(FLStocks(MSEsims00fixed), scn01)->df0
> xyplot(data~factor(year)|par, groups=qtl, data=subset(MSEsumm01, runid==15), scales=list(y=list(relation="free")), par.settings=list(superpose.line=list(lty=c(2,1,2), col=c(2,1,2))), type="l")
Error in xyplot(data ~ factor(year) | par, groups = qtl, data = subset(MSEsumm01,  : 
  error in evaluating the argument 'data' in selecting a method for function 'xyplot': Error in subset(MSEsumm01, runid == 15) : object 'MSEsumm01' not found
> xyplot(data~factor(year)|par, groups=qtl, data=df0, scales=list(y=list(relation="free")), par.settings=list(superpose.line=list(lty=c(2,1,2), col=c(2,1,2))), type="l")


