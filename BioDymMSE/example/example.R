#####################################################################################################################
#  Problem for fixing
# 1) still need to source ggplotFL plots; problem is that it tries to plot FLPBRP object
# 2) need to add flag to FLSR to indicate fit status
# 3) Be able to create FLSRs object etc
# *) pass names from FLlst to FLlst objects when using llply, also should create FLSRs as required
# 5) FLSR plot should be able to plot FLSRs
# *) FLBRP does not work with iters, e.g. see ssb(omBrp[[1]]) and parans(omBrp[[1]])
# *) Broken recycling in fwdSetup
# *) fwd Not tested yet
# *) No qmax() for FLPar
# 10) need to be able to flag iters to re-run, i.e. those where b-ve
# 11) Arithmetic refpts * FLQuant not working
# 12) Replace hard coded years by variables
# 13) rec sd is fixed on 0.3, should be user setup
#####################################################################################################################

library(FLAdvice)
library(FLBioDym)
source("funs.R")

#### OM Conditioning ##############################################################
nits <- 5

#### VPA Results for OM############################################################
data(ple4)
om <- ple4
plot(om)

#### SRR ########################################################################
## 1st fit
sr <- fmle(as.FLSR(om,model="bevholt"))
plot(sr)

## AR(1)
srAR1 <- fmle(as.FLSR(om,model="bevholtAR1"))
plot(srAR1)
vcov(srAR1)

#sv(srAR1,spr0=spr0(FLBRP(om)))
sv(params(srAR1)[1:2], model=SRModelName(model(srAR1)), spr0=spr0(FLBRP(om)))["s"]
         
srSG <- fmle(as.FLSR(om,model="segreg"))
srGM <- fmle(as.FLSR(om,model="geomean"))

srs <- FLSRs("Mean" = srGM, "Hockey Stick" = srSG, "Beverton & Holt" = sr, "Beverton & Holt AR(1)" = srAR1)


#srDiags <- ldply(srs, diags) # not working yet

## Plot Functional Forms
getFF <- function(sr){
	ssb <- FLQuant(seq(0, signif(max(ssb(sr)*1.1), 2), length.out=101))
	return(model.frame(FLQuants("Recruits"=predict(sr,ssb=ssb), "SSB"=ssb), drop=TRUE))
}

fig1 <- ggplot(transform(ldply(srs, getFF), Recruits = Recruits/1000, SSB = SSB/1000, SRR = factor(.id, levels=c("Mean","Hockey Stick","Beverton & Holt","Beverton & Holt AR(1)")))) + geom_line(aes(SSB, Recruits, group=SRR, colour=SRR)) + geom_point(aes(ssb,rec), data=transform(model.frame(srs[[1]]), rec=rec/1000, ssb=ssb/1000), colour="blue") 


# not working yet, still to decide on diagnostics
#fig2 <- ggplot(subset(diags(srs[["Beverton & Holt"]]), name=="AR(1) Residuals")) + geom_point(aes(x,y)) + geom_smooth(aes(x,y),method = "lm") + xlab(expression(Residual[t])) + ylab(expression(Residual[t+1]))
    
#### Residuals
# meth 1 - bootstrap residuals simulated has lognormal with sd=historical sd
srRsdlAR <- exp(simAR(sr, nits, 2004, 2020, 15))
srRsdl <- srRsdlAR
n <- prod(dim(srRsdlAR[,ac(2004:2020)]))
srRsdl[,ac(2004:2020)] <- c(srRsdlAR[,ac(2004:2020)])[sample(1:n)] # why rearranging residuals here ? 

# meth 2 - simulate residuals has lognormal with sd=0.3
srRsdl <- FLQuant(rlnorm(17*nits, 0, .3), dimnames=list(year=2004:2020, iter=1:nits)) 

## BRPs ###################################################################################
omBrp <- FLBRPs("Mean" = brp(FLBRP(om, nyears=dims(om)$year, sr=srGM)), "Hockey Stick"=brp(FLBRP(om,nyears=dims(om)$year,sr=srSG)))

omBrp <- FLBRPs(llply(omBrp, function(x) {
	fbar(x) <- seq(0, 0.5, length.out=101) # fixed max fbar as 0.5
	brp(x)
}))

fig3a <- plot(omBrp[[1]], obs=TRUE)
fig3b <- plot(omBrp)
    
fig4 <- ggplot(as.data.frame(omBrp[[1]][[c("catch.sel","catch.wt","mat","m")]])) + geom_line(aes(age,data)) + facet_wrap(~qname,scales="free")

save(om,srs,srRsdlAR,srRsdl,omBrp,file="om.RData")

#################################################################################################################

##### MP ########################################################################################################
bd <- as(om,"FLBioDym")
bd <- admbBD(bd)
plot(bd)
plot(residuals(bd), type="l")
hist(residuals(bd))
qqnorm(residuals(bd))

plotSP(bd)
kobe(model.frame(mcf(bd[[c("stock","harvest")]]))) + geom_path(aes(stock/c(bmsy(bd)), harvest/c(fmsy(bd)))) + geom_point(aes(stock/c(bmsy(bd)),harvest/c(fmsy(bd)),colour=year),size=1.25)

cf <- transform(rbind(cbind(MSE="MP", model.frame(bd[[c("stock","catch")]], drop=TRUE)), cbind(MSE="OM", model.frame(om[[c("stock","catch")]], drop=TRUE))), harvest=catch/stock)

rp <- rbind(data.frame(MSE="MP", "bmsy"=c(bmsy(bd)), "msy"=c(msy(bd))), data.frame(MSE="OM", "bmsy"=c(refpts(omBrp[[1]])["msy","biomass"]), "msy"=c(refpts(omBrp[[1]])["msy","yield"][,"yield"])))

rel <- merge(cf,rp)
rel <- transform(rel,stock=stock/bmsy,catch=catch/msy,harvest=harvest*bmsy/msy)[,1:5]

fig5 <- ggplot(melt(rel,id=c("MSE","year"))) + geom_line(aes(year,value,group=MSE,colour=MSE)) + facet_wrap(~variable,scale="free",nrow=3) 

fig6 <- ggplot(melt(cf,id=c("MSE","year"))) + geom_line(aes(year,value,group=MSE,colour=MSE)) + facet_wrap(~variable,scale="free",nrow=3) 

fig7 <- kobe(rel) + geom_path(aes(stock,harvest,group=MSE,colour=year,size=MSE))

scl <- ddply(rel,.(MSE),function(x) c(sMn=mean(x$stock),yMn=mean(x$catch),fMn=mean(x$harvest)))

fig8 <- kobe(transform(merge(rel,scl),stock=stock/sMn,harvest=harvest/fMn)) + geom_path(aes(stock,harvest,group=MSE,colour=MSE),size=2)

save(fig1,fig3a,fig3b,fig4,fig5,fig6,fig7,fig8, file="figs-OM.RData")

### Sims ###################################################################################
OM <- window(om, FLBRP=omBrp[[1]], end=2020)
m(OM) <- propagate(m(OM),nits)
ctrl <- fwdControl(data.frame(year=2005:2020,quantity="f",val=c(fbar(om)[,"2004"])*seq(0.75,0.75,length.out=16)))
OM <- fwd(OM, ctrl=ctrl, sr=srs[["Mean"]], sr.residuals=srRsdl)
plot(OM)

# these bounds are related with ADMB 
bounds <- bounds(FLBioDym())
bounds["r",    "start"]=0.35
bounds["K",    "start"]=1.6e+8
bounds["sigma","start"]=0.5
bounds["q",    "start"]=1.0
bounds["q",    1]      =1.0
bounds["p",    c("phase","start")]=c(-1,1)
bounds["b0",   c("phase","start")]=c(-1,0.5)
bounds[,"lower"]=bounds[,"start"]*0.1
bounds[,"upper"]=bounds[,"start"]*10.0

# scenarios
scn <- expand.grid(Btrig=c(0.4,0.75), CV=c(0.3,0.5), Ftar=0.75)

# MSE
set.seed(123)
res1 <- mlply(scn, function(Btrig, CV, Ftar) 
	mseFLBioDym(OM, 2004, sr=srs[["Hockey Stick"]], srRsdl, CV=CV, Btrig=Btrig, Ftar=Ftar, bounds=bounds))
save(res1,file="res1.RData")

# EJ tests to introduce distinct lag in advice
set.seed(123)
res10 <- mlply(scn, function(Btrig, CV, Ftar) 
	mseBD(OM, 2004, sr=srs[["Hockey Stick"]], srRsdl, CV=CV, Btrig=Btrig, Ftar=Ftar, bounds=bounds))

set.seed(123)
res11 <- mlply(scn, function(Btrig, CV, Ftar) 
	mseBD(OM, 2004, sr=srs[["Hockey Stick"]], srRsdl, CV=CV, Btrig=Btrig, Ftar=Ftar, bounds=bounds, intYr=2))

set.seed(123)
res12 <- mlply(scn, function(Btrig, CV, Ftar) 
	mseBD(OM, 2004, sr=srs[["Hockey Stick"]], srRsdl, CV=CV, Btrig=Btrig, Ftar=Ftar, bounds=bounds, intYr=3))


res2 <- mlply(scn, function(Btrig,CV,Ftar) 
   mseFLBioDym(OM, 2004, sr=srs[["Hockey Stick"]], srRsdlAR, CV=CV, Btrig=Btrig, Ftar=Ftar, bounds=bounds))
save(res2,file="res2.RData")

res3=mlply(scn, function(Btrig,CV,Ftar) 
   mseFLBioDym(OM, 2004, sr=srs[["Mean"]], srRsdl, CV=CV, Btrig=Btrig, Ftar=Ftar, bounds=bounds))
save(res3,file="res3.RData")

res4=mlply(scn, function(Btrig,CV,Ftar) 
   mseFLBioDym(OM, 2004, sr=srs[["Mean"]], srRsdlAR, CV=CV, Btrig=Btrig, Ftar=Ftar, bounds=bounds))
save(res4,file="res4.RData")

# proj without management loop ??            
OMs <- FLStocks()
OMs[[1]] <- fwd(OM, ctrl=ctrl, sr=srs[["Mean"]], sr.residuals=srRsdl)
OMs[[2]] <- fwd(OM, ctrl=ctrl, sr=srs[["Mean"]], sr.residuals=srRsdlAR)
OMs[[3]] <- fwd(OM, ctrl=ctrl, sr=srs[["Hockey Stick"]], sr.residuals=srRsdl)
OMs[[4]] <- fwd(OM, ctrl=ctrl, sr=srs[["Hockey Stick"]], sr.residuals=srRsdlAR)
save(OMs,file="proj.RData")
names(OMs)=1:4

# plots
plot(FLStocks(llply(OMs, window,start=2004)))
plot(FLStocks(llply(res1,window,start=2004)))
plot(FLStocks(llply(res2,window,start=2004)))
plot(FLStocks(llply(res3,window,start=2004)))
plot(FLStocks(llply(res4,window,start=2004)))

