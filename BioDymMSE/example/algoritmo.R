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
source("funs.R")

#====================================================================
# General settings
#====================================================================

nits <- 250				# number of iterations
iniyr <- 2009 			# first year in projections
lastyr <- 2059 			# last year in projections
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
scn <- expand.grid(Btrig=0.5, CV=0.2, Ftar=1)

# MSE
set.seed(123)
res00 <- mlply(scn, function(Btrig, CV, Ftar) 
	mseBD(OM, iniyr, sr=srSG, srRsdl, CV=CV, Btrig=Btrig, Ftar=Ftar, bounds=bounds, aLag=1))

set.seed(123)
res1 <- mlply(scn, function(Btrig, CV, Ftar) 
	mseBD(OM, iniyr, sr=srSG, srRsdl, CV=CV, Btrig=Btrig, Ftar=Ftar, bounds=bounds, aLag=3))

set.seed(123)
res2 <- mlply(scn, function(Btrig, CV, Ftar) 
	mseBD(OM, iniyr, sr=srSG, srRsdl, CV=CV, Btrig=Btrig, Ftar=Ftar, bounds=bounds, aLag=5))

set.seed(123)
res3 <- mlply(scn, function(Btrig, CV, Ftar) 
	mseBD(OM, iniyr, sr=srBH, srRsdl, CV=CV, Btrig=Btrig, Ftar=Ftar, bounds=bounds, aLag=5))

set.seed(123)
res40 <- mlply(scn, function(Btrig, CV, Ftar) 
	mseBD(OM, iniyr, sr=srBHAR1, srRsdlAR, CV=CV, Btrig=Btrig, Ftar=Ftar, bounds=bounds, aLag=5))

set.seed(123)
res5 <- mlply(scn, function(Btrig, CV, Ftar) 
	mseBD(OM, iniyr, sr=srBHAR1, srRsdlAR, CV=CV, Btrig=Btrig, Ftar=Ftar, bounds=bounds, aLag=5, cthBias=0.5))

set.seed(123)
res6 <- mlply(scn, function(Btrig, CV, Ftar) 
	mseBD(OM, iniyr, sr=srBHAR1, srRsdlAR, CV=CV, Btrig=Btrig, Ftar=Ftar, bounds=bounds, aLag=5, srvBias=0.5))

# new code with control on HR and limit on F
set.seed(123)
res7 <- mlply(scn, function(Btrig, CV, Ftar) 
	mseBD(OM, iniyr, sr=srBHAR1, srRsdlAR, CV=CV, Btrig=Btrig, Ftar=Ftar, bounds=bounds, aLag=5))





plot(FLStocks(lag1=res0[[1]], lag3=res1[[1]], lag5=res2[[1]]))
plot(FLStocks(sg=res2[[1]], bh=res3[[1]], bhAR=res4[[1]]))
plot(FLStocks(nobias=res4[[1]], cthBias0.5=res5[[1]], srvBias0.5=res6[[1]]))




#====================================================================
# proj without management loop            
#====================================================================

OMs <- FLStocks()
OMs[[1]] <- fwd(OM, ctrl=ctrl, sr=srs[["Mean"]], sr.residuals=srRsdl)
OMs[[2]] <- fwd(OM, ctrl=ctrl, sr=srs[["Mean"]], sr.residuals=srRsdlAR)
OMs[[3]] <- fwd(OM, ctrl=ctrl, sr=srs[["Hockey Stick"]], sr.residuals=srRsdl)
OMs[[4]] <- fwd(OM, ctrl=ctrl, sr=srs[["Hockey Stick"]], sr.residuals=srRsdlAR)
names(OMs)=1:4

q("yes")

