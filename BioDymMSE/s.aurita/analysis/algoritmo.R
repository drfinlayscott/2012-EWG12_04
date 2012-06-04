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
#registerDoMC(2)
source("funs.R")

#====================================================================
# Simulation settings
#====================================================================

nits <- 500				# number of iterations
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
bounds["p",   c("phase","start")]=c(-1,1)
bounds[,"lower"]=bounds[,"start"]*0.1
bounds[,"upper"]=bounds[,"start"]*10.0
#bounds["p",    c("phase","lower","upper","start")] = c(1,1,1.1,1)
bounds["r",    c("phase","lower","upper","start")] = c(1, 0.32, 1.28, 0.64)
bounds["K",    c("phase","lower","upper","start")] = c(2, 1750, 7000, 2500)

#====================================================================
# Scenarios and simulations
#====================================================================


#====================================================================
# Analysis of results
#====================================================================

#--------------------------------------------------------------------
# Read results and process
#--------------------------------------------------------------------

load("RData.MSEcth.summ")
load("RData.MSEsrv.summ")
load("RData.MSEbd.b0.2.unif.summ")

# all scenarios
scn <- rbind(
	expand.grid(Btrig=0.5, Ftar=0.75, maxHR=0.35, aLag=c(1,3,5), srvBias=c(0.5, 1, 1.5), cthBias=c(0.5,1), IEM="linear", slag=NA, clag=NA, am="bd"),
	expand.grid(Btrig=0.3, Ftar=1, maxHR=c(0.35, 1), aLag=c(1,3,5), srvBias=c(0.5, 1, 1.5), cthBias=c(0.5,1), IEM="linear", slag=NA, clag=NA, am="bd"),
	expand.grid(Btrig=NA, Ftar=NA, maxHR=NA, aLag=c(1,3,5), srvBias=c(0.5, 1, 1.5), cthBias=c(0.5,1), IEM="linear", slag=NA, clag=c(5,10,20), am="cth"),
	expand.grid(Btrig=NA, Ftar=NA, maxHR=NA, aLag=c(1,3,5), srvBias=c(0.5, 1, 1.5), cthBias=c(0.5,1), IEM="linear", slag=c(5,10,20), clag=c(5,10,20), am="srv")
	)

# puting all together

MSEsumm <- as.data.frame(matrix(NA, nrow=nrow(MSEbd.b0.2.unif.summ)+nrow(MSEsrv.summ)+nrow(MSEcth.summ), ncol=ncol(MSEsrv.summ)))
names(MSEsumm) <- names(MSEsrv.summ)
MSEsumm$am <- "bd"
MSEsumm$b0 <- NA

n1 <- nrow(MSEbd.b0.2.unif.summ)
n2 <- nrow(MSEsrv.summ) 
n3 <- nrow(MSEcth.summ)

MSEsumm[1:n1,names(MSEbd.b0.2.unif.summ)] <- MSEbd.b0.2.unif.summ
MSEsumm[1:n1,"b0"] <- 0.2

MSEsumm[(n1+1):(n1+n2),names(MSEsrv.summ)] <- MSEsrv.summ
MSEsumm[(n1+1):(n1+n2),"am"] <- "srv"
MSEsumm[(n1+1):(n1+n2),"Btrig"] <- NA
MSEsumm[(n1+1):(n1+n2),"Ftar"] <- NA
MSEsumm[(n1+1):(n1+n2),"maxHR"] <- NA

MSEsumm[(n1+n2+1):(n1+n2+n3),names(MSEcth.summ)] <- MSEcth.summ
MSEsumm[(n1+n2+1):(n1+n2+n3),"am"] <- "cth"
MSEsumm[(n1+1):(n1+n2),"Btrig"] <- NA
MSEsumm[(n1+1):(n1+n2),"Ftar"] <- NA
MSEsumm[(n1+1):(n1+n2),"maxHR"] <- NA
MSEsumm[(n1+1):(n1+n2),"slag"] <- NA

MSEsumm[,"season"] <- NULL
MSEsumm[,"area"] <- NULL
MSEsumm[,"unit"] <- NULL
MSEsumm[,"iter"] <- NULL

#--------------------------------------------------------------------
# Example plots
#--------------------------------------------------------------------

xyplot(data~factor(year)|par, groups=qtl, data=subset(MSEsumm, Ftar==0.75 & Btrig==0.5 & maxHR==0.35 & aLag==1 & srvBias==1 & cthBias==1), scales=list(y=list(relation="free")), par.settings=list(superpose.line=list(lty=c(2,1,2), col=c(2,1,2))), type="l", main="base case")

xyplot(data~factor(year), groups=par, data=subset(MSEsumm, Ftar==0.75 & Btrig==0.5 & maxHR==0.35 & aLag==1 & srvBias==1 & cthBias==1 & par %in% c("HCR:TAC", "catch") & qtl==0.5), type="l")

xyplot(data~factor(year)|scn, groups=par, data=subset(MSEsummBaseb0.5, par %in% c("HCR:TAC", "catch") & qtl==0.5 & year > 2000), type="l", par.strip.text=list(cex=0.5), layout=c(3,2))

subset(MSEsumm, runid==55)

xyplot(data~factor(year)|par, groups=qtl, data=subset(MSEsumm, runid==55), scales=list(y=list(relation="free")), par.settings=list(superpose.line=list(lty=c(2,1,2), col=c(2,1,2))), type="l", main=unique(subset(MSEsumm, runid==55)$scn))

