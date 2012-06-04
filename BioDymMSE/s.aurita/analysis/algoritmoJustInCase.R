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
source("funs.R")

#====================================================================
# Simulation settings
#====================================================================

nits <- 150				# number of iterations
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
par <- as(data.frame(linf=30, k=0.4, t0=-1, s=0.8, v=1750, a50=1), 'FLPar')
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

#====================================================================
# MP settings
#====================================================================

# these bounds are related with ADMB adjusted to CECAF usage
bounds <- bounds(FLBioDym())
bounds["sigma","start"]=0.50
bounds["q",    "start"]=1.0
bounds["q",    1]      =1.0
bounds["b0",   c("phase","start")]=c(-1,0.2)
bounds["p",   c("phase","start")]=c(-1,1)
bounds[,"lower"]=bounds[,"start"]*0.1
bounds[,"upper"]=bounds[,"start"]*10.0
bounds["r",    c("phase","lower","upper","start")] = c(1, 0.32, 1.28, 0.64)
bounds["K",    c("phase","lower","upper","start")] = c(2, 1750, 7000, 2500)

#====================================================================
# Scenarios and simulations
#====================================================================

scn <- mkScn(ref="base") 

base <- mseBD(OM=OM, start=iniyr, sr=srBH, srRsdl=srRsdl, bounds=bounds, Btrig=0.5, Ftar=0.75, maxHR=0.35, aLag=1, srvBias=1, cthBias=1, IEM="linear")

#--------------------------------------------------------------------
# assessment lag effect
#--------------------------------------------------------------------

scn[2,] <- mkScn(ref="aLag3", aLag=3, runid=2) 
scn[3,] <- mkScn(ref="aLag5", aLag=5, runid=3) 

base.aLag3 <- mseBD(OM=OM, start=iniyr, sr=srBH, srRsdl=srRsdl, bounds=bounds, Btrig=0.5, Ftar=0.75, maxHR=0.35, aLag=3, srvBias=1, cthBias=1, IEM="linear")

base.aLag5 <- mseBD(OM=OM, start=iniyr, sr=srBH, srRsdl=srRsdl, bounds=bounds, Btrig=0.5, Ftar=0.75, maxHR=0.35, aLag=5, srvBias=1, cthBias=1, IEM="linear")

#--------------------------------------------------------------------
# underreporting effect
#--------------------------------------------------------------------

scn[4,] <- mkScn(ref="cthBias0.5", runid=4, cthBias=0.5) 

base.cthBias0.5 <- mseBD(OM=OM, start=iniyr, sr=srBH, srRsdl=srRsdl, bounds=bounds, Btrig=0.5, Ftar=0.75, maxHR=0.35, aLag=1, srvBias=1, cthBias=0.5, IEM="linear")

#--------------------------------------------------------------------
# survey coverage
#--------------------------------------------------------------------
scn[5,] <- mkScn(ref="srvBias0.5", runid=5, srvBias=0.5) 
scn[6,] <- mkScn(ref="srvBias1.5", runid=6, srvBias=0.5) 

base.srvBias0.5 <- mseBD(OM=OM, start=iniyr, sr=srBH, srRsdl=srRsdl, bounds=bounds, Btrig=0.5, Ftar=0.75, maxHR=0.35, aLag=1, srvBias=0.5, cthBias=1, IEM="linear")

base.srvBias1.5 <- mseBD(OM=OM, start=iniyr, sr=srBH, srRsdl=srRsdl, bounds=bounds, Btrig=0.5, Ftar=0.75, maxHR=0.35, aLag=1, srvBias=1.5, cthBias=1, IEM="linear")

#--------------------------------------------------------------------
# mp effect
#--------------------------------------------------------------------

scn[7,] <- mkScn(ref="srv", runid=7, Btrig=NA, Ftar=NA, maxHR=NA, slag=5, clag=5, am="srv", b0=NA) 
scn[8,] <- mkScn(ref="cth", runid=8, Btrig=NA, Ftar=NA, maxHR=NA, clag=5, am="cth", b0=NA) 

base.srv <- mseBD(OM=OM, start=iniyr, sr=srBH, srRsdl=srRsdl, bounds=bounds, Btrig=0.5, Ftar=0.75, maxHR=0.35, aLag=1, srvBias=1, cthBias=1, IEM="linear", slag=5, clag=5, am="srv")

base.cth <- mseBD(OM=OM, start=iniyr, sr=srBH, srRsdl=srRsdl, bounds=bounds, Btrig=0.5, Ftar=0.75, maxHR=0.35, aLag=1, srvBias=1, cthBias=1, IEM="linear", clag=5, am="cth")

#--------------------------------------------------------------------
# HCR effect
#--------------------------------------------------------------------
scn[9,] <- mkScn(ref="maxHR1", runid=9, maxHR=1) 
scn[10,] <- mkScn(ref="Btrig1", runid=10, Btrig=1) 
scn[11,] <- mkScn(ref="Ftar1", runid=11, Ftar=1) 

base.maxHR1 <- mseBD(OM=OM, start=iniyr, sr=srBH, srRsdl=srRsdl, bounds=bounds, Btrig=0.5, Ftar=0.75, maxHR=1, aLag=1, srvBias=1, cthBias=1, IEM="linear")

base.Btrig1 <- mseBD(OM=OM, start=iniyr, sr=srBH, srRsdl=srRsdl, bounds=bounds, Btrig=1, Ftar=0.75, maxHR=0.35, aLag=1, srvBias=1, cthBias=1, IEM="linear")

base.Ftar1 <- mseBD(OM=OM, start=iniyr, sr=srBH, srRsdl=srRsdl, bounds=bounds, Btrig=0.5, Ftar=1, maxHR=0.35, aLag=1, srvBias=1, cthBias=1, IEM="linear")

#--------------------------------------------------------------------
# b0 effect
#--------------------------------------------------------------------

scn[12,] <- mkScn(ref="b0.5", runid=12, b0=0.5) 

bounds["b0",   c("phase","start")]=c(-1,0.5)

base.b0.5 <- mseBD(OM=OM, start=iniyr, sr=srBH, srRsdl=srRsdl, bounds=bounds, Btrig=0.5, Ftar=0.75, maxHR=0.35, aLag=1, srvBias=1, cthBias=1, IEM="linear")

bounds["b0",   c("phase","start")]=c(-1,0.2)

#--------------------------------------------------------------------
# survey slag=1
#--------------------------------------------------------------------

scn[13,] <- mkScn(ref="srv.slag1", runid=13, Btrig=NA, Ftar=NA, maxHR=NA, slag=1, clag=5, am="srv", b0=NA) 
scn[14,] <- mkScn(ref="srv.alag5", runid=14, Btrig=NA, Ftar=NA, maxHR=NA, aLag=5, slag=5, clag=5, am="srv", b0=NA) 
scn[15,] <- mkScn(ref="srv.alag5bias0.5", runid=15, Btrig=NA, Ftar=NA, maxHR=NA, srvBias=0.5, aLag=5, slag=5, clag=5, am="srv", b0=NA) 

base.srvslag1 <- mseBD(OM=OM, start=iniyr, sr=srBH, srRsdl=srRsdl, bounds=bounds, Btrig=0.5, Ftar=0.75, maxHR=0.35, aLag=1, srvBias=1, cthBias=1, IEM="linear", slag=1, clag=5, am="srv")

#--------------------------------------------------------------------
# survey alag=5
#--------------------------------------------------------------------
base.srvslag5alag5 <- mseBD(OM=OM, start=iniyr, sr=srBH, srRsdl=srRsdl, bounds=bounds, Btrig=0.5, Ftar=0.75, maxHR=0.35, aLag=5, srvBias=1, cthBias=1, IEM="linear", slag=5, clag=5, am="srv")

#--------------------------------------------------------------------
# survey slag=5 srvBias=0.5
#--------------------------------------------------------------------
base.srvslag5alag5srvBias0.5 <- mseBD(OM=OM, start=iniyr, sr=srBH, srRsdl=srRsdl, bounds=bounds, Btrig=0.5, Ftar=0.75, maxHR=0.35, aLag=5, srvBias=0.5, cthBias=1, IEM="linear", slag=5, clag=5, am="srv")

#--------------------------------------------------------------------
# survey simetric multiplier NOTE: requires changing the code 
# multipliers are hardcoded
#--------------------------------------------------------------------

base.srvSim <- mseBD(OM=OM, start=iniyr, sr=srBH, srRsdl=srRsdl, bounds=bounds, Btrig=0.5, Ftar=0.75, maxHR=0.35, aLag=1, srvBias=1, cthBias=1, IEM="linear", slag=5, clag=5, am="srv")

scn[16,] <- mkScn(ref="srvSimetric", runid=16) 

#--------------------------------------------------------------------
# IEM="cst"
#--------------------------------------------------------------------

base.iemCts <- mseBD(OM=OM, start=iniyr, sr=srBH, srRsdl=srRsdl, bounds=bounds, Btrig=0.5, Ftar=0.75, maxHR=0.35, aLag=1, srvBias=1, cthBias=0.5, IEM="cst")

scn[17,] <- mkScn(ref="iemCts", runid=17, cthBias=0.5, IEM="cts") 

#====================================================================
# Analysis
#====================================================================

lst <- list(base, base.aLag3, base.aLag5, base.cthBias0.5, base.srvBias0.5, base.srvBias1.5, base.srv, base.cth, base.maxHR1, base.Btrig1, base.Ftar1, base.b0.5, base.srvslag1, base.srvslag5alag5, base.srvslag5alag5srvBias0.5, base.srvSim, base.iemCts)

base.summ <- mseSumm(lst, scn)

save(base, base.aLag3, base.aLag5, base.b0.5, base.Btrig1, base.cth, base.cthBias0.5, base.Ftar1, base.iemCts, base.maxHR1, base.srv, base.srvBias0.5, base.srvBias1.5, base.srvSim, base.srvslag1, base.srvslag5alag5, base.srvslag5alag5srvBias0.5, base.iemCts, base.summ, scn, file="../report/RData.mse")

xyplot(data~factor(year)|par, groups=qtl, data=subset(base.summ, ref=="base"), scales=list(y=list(relation="free")), par.settings=list(superpose.line=list(lty=c(2,1,2), col=c(2,1,2))), type="l", main="base case")

xyplot(data~factor(year), groups=par, data=subset(base.summ, ref=="base" & par %in% c("HCR:TAC", "catch") & qtl==0.5), type="l")

q("yes")

