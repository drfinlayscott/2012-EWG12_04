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
library(multicore)
library(foreach)
library(doMC)
registerDoMC(3)
source("funs.R")

#====================================================================
# General settings
#====================================================================

nits <- 250				# number of iterations
iniyr <- 2009 			# first year in projections
lastyr <- 2039 			# last year in projections
npyr <- lastyr-iniyr+1 	# number of years to project
srsd <- 0.3 			# sd for S/R

#====================================================================
# Conditioning
#====================================================================

#--------------------------------------------------------------------
# VPA Results for OM
#--------------------------------------------------------------------

data(ple4)
om <- ple4

#--------------------------------------------------------------------
# SRR 
#--------------------------------------------------------------------
# models
srBH <- fmle(as.FLSR(om,model="bevholt"))
srBHAR1 <- fmle(as.FLSR(om,model="bevholtAR1"))
srSG <- fmle(as.FLSR(om,model="segreg"))
srGM <- fmle(as.FLSR(om,model="geomean"))
    
# Residuals - simulate residuals as lognormal with sd=srsd
# I use a hacked simAR method to control residuals cv
set.seed(123)
srRsdl <- FLQuant(rlnorm(npyr*nits, 0, srsd), dimnames=list(year=srBH@range["minyear"]:lastyr, iter=1:nits)) 
set.seed(123)
srRsdlAR <- exp(simAR(residuals(srBH), nits, end=lastyr, cv=srsd))

#--------------------------------------------------------------------
# BRPs
#--------------------------------------------------------------------

omBrp <- FLBRPs("Mean" = brp(FLBRP(om, nyears=dims(om)$year, sr=srGM)), "Hockey Stick"=brp(FLBRP(om,nyears=dims(om)$year,sr=srSG)))

omBrp <- FLBRPs(llply(omBrp, function(x) {
	fbar(x) <- seq(0, 0.5, length.out=101) # fixed max fbar as 0.5 ?
	brp(x)
}))

#--------------------------------------------------------------------
# create OM object
# Note: this object is projected at Fsq and the values for the first
#	intermediate year are used in the projections
#--------------------------------------------------------------------

# window with FLBRP expands the object including weights, etc, the 
# brp doesn't seem to do anything except dispatching. it replaces "stf". 
OM <- window(om, FLBRP=omBrp[[1]], end=lastyr)

# trick to get iterations, start with M and fwd will add to other slots
m(OM) <- propagate(m(OM),nits)

# project to the end of projections at last year F level
ctrl <- fwdControl(data.frame(year=iniyr:lastyr, quantity="f", val=c(fbar(om)[,ac(om@range["maxyear"])])*rep(1,npyr)))
OM <- fwd(OM, ctrl=ctrl, sr=srBHAR1, sr.residuals=srRsdl)
plot(OM)

#====================================================================
# MP settings
#====================================================================

# these bounds are related with ADMB 
bounds <- bounds(FLBioDym())
bounds["r",    "start"]=0.35
bounds["K",    "start"]=1.6e+8
bounds["sigma","start"]=0.5
bounds["q",    "start"]=1.0
bounds["q",    1]      =1.0
bounds["b0",   c("phase","start")]=c(-1,0.2)
bounds[,"lower"]=bounds[,"start"]*0.1
bounds[,"upper"]=bounds[,"start"]*10.0
bounds["p",    c("phase","lower","upper","start")] = c(1,1,1.1,1)

#====================================================================
# Scenarios and simulations
#====================================================================

# scenarios
scn <- expand.grid(Btrig=0.5, CV=0.2, Ftar=1, aLag=c(1,3,5), srvBias=c(0.5, 1, 1.5), cthBias=c(0.5, 1), IEM=c("linear", "cst"))

MSEsims01 <- split(scn, scn)
for(i in 1:nrow(scn)){
	MSEsims01[[i]] <- mseBD(OM=OM, start=iniyr, sr=srBHAR1, srRsdl=srRsdlAR, bounds=bounds, CV=scn[i,"CV"], Btrig=scn[i,"Btrig"], Ftar=scn[i,"Ftar"], aLag=scn[i, "aLag"], srvBias=scn[i,"srvBias"], cthBias=scn[i,"cthBias"], IEM=scn[i,"IEM"])

}


# paralel
res0 <- foreach(i = 1:nrow(scn), .packages=c("FLCore", "FLash", "FLEDA", "MASS", "FLBioDym", "FLAdvice", "plyr"), .export="an") %dopar% mseBD(OM=OM, start=iniyr, sr=srBHAR1, srRsdl=srRsdlAR, bounds=bounds, CV=scn[i,"CV"], Btrig=scn[i,"Btrig"], Ftar=scn[i,"Ftar"], aLag=scn[i, "aLag"], srvBias=scn[i,"srvBias"], cthBias=scn[i,"cthBias"], IEM=scn[i,"IEM"])

#====================================================================
# Analysis of results
#====================================================================

MSEsumm01 <- mseSumm(MSEsims01, scn)

xyplot(data~factor(year)|par, groups=qtl, data=subset(MSEsumm01, runid==13), scales=list(y=list(relation="free")), par.settings=list(superpose.line=list(lty=c(2,1,2), col=c(2,1,2))), type="l")

xyplot(data~factor(year), groups=par, data=subset(MSEsumm01, runid==13 & par %in% c("TAC", "catch") & qtl==0.5), type="l")

xyplot(data~factor(year)|scn, groups=par, data=subset(MSEsumm01, runid<19 & par %in% c("TAC", "catch") & qtl==0.5 & year > 2000), type="l", par.strip.text=list(cex=0.3))



