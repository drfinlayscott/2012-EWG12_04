# Base scenarios
scn <- rbind(
	expand.grid(Btrig=0.5, Ftar=0.75, maxHR=0.35, aLag=5, srvBias=1, cthBias=1, IEM="linear"),
	expand.grid(Btrig=0.3, Ftar=1, maxHR=0.35, aLag=5, srvBias=1, cthBias=1, IEM="linear"),
	expand.grid(Btrig=0.5, Ftar=0.75, maxHR=0.35, aLag=3, srvBias=1, cthBias=1, IEM="linear"),
	expand.grid(Btrig=0.3, Ftar=1, maxHR=0.35, aLag=3, srvBias=1, cthBias=1, IEM="linear")#,
	#expand.grid(Btrig=0.5, Ftar=0.75, maxHR=0.35, aLag=1, srvBias=1, cthBias=1, IEM="linear"),
	#expand.grid(Btrig=0.3, Ftar=1, maxHR=0.35, aLag=1, srvBias=1, cthBias=1, IEM="linear")
	)

# b0=0.2
MSEsimsBaseb0.2 <- foreach(i = 1:nrow(scn), .packages=c("FLCore", "FLash", "MASS", "FLBioDym", "FLAdvice", "plyr"), .export="an") %dopar% mseBD(OM=OM, start=iniyr, sr=srBH, srRsdl=srRsdl, bounds=bounds, Btrig=scn[i,"Btrig"], Ftar=scn[i,"Ftar"], aLag=scn[i, "aLag"], srvBias=scn[i,"srvBias"], cthBias=scn[i,"cthBias"], IEM=scn[i,"IEM"])

MSEsimBaseb0.2lag1Ftar1 <- mseBD(OM=OM, start=iniyr, sr=srBH, srRsdl=srRsdl, bounds=bounds, Btrig=0.3, Ftar=1, aLag=1, srvBias=1, cthBias=1, IEM="linear")

objmse <- mseBD(OM=OM, start=iniyr, sr=srBH, srRsdl=srRsdl, am="srv")
objmse2 <- mseBD(OM=OM, start=iniyr, sr=srBH, srRsdl=srRsdl, am="cth", aLag=5)
objmse3 <- mseBD(OM=OM, start=iniyr, sr=srBH, srRsdl=srRsdl, am="srv", aLag=5)
objmse4 <- mseBD(OM=OM, start=iniyr, sr=srBH, srRsdl=srRsdl, bounds=bounds, Btrig=0.3, Ftar=1, aLag=1, srvBias=1, cthBias=1, IEM="linear")


# b0=0.5

bounds["b0",   c("phase","start")]=c(-1,0.5)
bounds[,"lower"]=bounds[,"start"]*0.1
bounds[,"upper"]=bounds[,"start"]*10.0

MSEsimsBaseb0.5 <- foreach(i = 1:nrow(scn), .packages=c("FLCore", "FLash", "MASS", "FLBioDym", "FLAdvice", "plyr"), .export="an") %dopar% mseBD(OM=OM, start=iniyr, sr=srBH, srRsdl=srRsdl, bounds=bounds, Btrig=scn[i,"Btrig"], Ftar=scn[i,"Ftar"], aLag=scn[i, "aLag"], srvBias=scn[i,"srvBias"], cthBias=scn[i,"cthBias"], IEM=scn[i,"IEM"])

