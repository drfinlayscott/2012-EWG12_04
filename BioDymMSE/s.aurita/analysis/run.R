scn <- mkScn(ref="base") 
scn[2,] <- mkScn(ref="aLag3", aLag=3, runid=2) 
scn[3,] <- mkScn(ref="aLag5", aLag=5, runid=3) 
scn[4,] <- mkScn(ref="cthBias0.5", runid=4, cthBias=0.5) 
scn[5,] <- mkScn(ref="srvBias0.5", runid=5, srvBias=0.5) 
scn[6,] <- mkScn(ref="srvBias1.5", runid=6, srvBias=0.5) 
scn[7,] <- mkScn(ref="srv", runid=7, Btrig=NA, Ftar=NA, maxHR=NA, slag=5, clag=5, am="srv", b0=NA) 
scn[8,] <- mkScn(ref="cth", runid=8, Btrig=NA, Ftar=NA, maxHR=NA, clag=5, am="cth", b0=NA) 
scn[9,] <- mkScn(ref="maxHR1", runid=9, maxHR=1) 
scn[10,] <- mkScn(ref="Btrig1", runid=10, Btrig=1) 
scn[11,] <- mkScn(ref="Ftar1", runid=11, Ftar=1) 
scn[12,] <- mkScn(ref="b0.5", runid=12, b0=0.5) 
scn[13,] <- mkScn(ref="srv.slag1", runid=13, Btrig=NA, Ftar=NA, maxHR=NA, slag=1, clag=5, am="srv", b0=NA) 
scn[14,] <- mkScn(ref="srv.alag5", runid=14, Btrig=NA, Ftar=NA, maxHR=NA, aLag=5, slag=5, clag=5, am="srv", b0=NA) 
scn[15,] <- mkScn(ref="srv.alag5bias0.5", runid=15, Btrig=NA, Ftar=NA, maxHR=NA, srvBias=0.5, aLag=5, slag=5, clag=5, am="srv", b0=NA) 
scn[16,] <- mkScn(ref="srvSimetric", runid=16) 
scn[17,] <- mkScn(ref="iemCts", runid=17, cthBias=0.5, IEM="cts") 


base <- mseBD(OM=OM, start=iniyr, sr=srBH, srRsdl=srRsdl, bounds=bounds, Btrig=0.5, Ftar=0.75, maxHR=0.35, aLag=1, srvBias=1, cthBias=1, IEM="linear")

base.aLag3 <- mseBD(OM=OM, start=iniyr, sr=srBH, srRsdl=srRsdl, bounds=bounds, Btrig=0.5, Ftar=0.75, maxHR=0.35, aLag=3, srvBias=1, cthBias=1, IEM="linear")

base.aLag5 <- mseBD(OM=OM, start=iniyr, sr=srBH, srRsdl=srRsdl, bounds=bounds, Btrig=0.5, Ftar=0.75, maxHR=0.35, aLag=5, srvBias=1, cthBias=1, IEM="linear")

base.cthBias0.5 <- mseBD(OM=OM, start=iniyr, sr=srBH, srRsdl=srRsdl, bounds=bounds, Btrig=0.5, Ftar=0.75, maxHR=0.35, aLag=1, srvBias=1, cthBias=0.5, IEM="linear")

base.srvBias0.5 <- mseBD(OM=OM, start=iniyr, sr=srBH, srRsdl=srRsdl, bounds=bounds, Btrig=0.5, Ftar=0.75, maxHR=0.35, aLag=1, srvBias=0.5, cthBias=1, IEM="linear")

base.srvBias1.5 <- mseBD(OM=OM, start=iniyr, sr=srBH, srRsdl=srRsdl, bounds=bounds, Btrig=0.5, Ftar=0.75, maxHR=0.35, aLag=1, srvBias=1.5, cthBias=1, IEM="linear")

base.srv <- mseBD(OM=OM, start=iniyr, sr=srBH, srRsdl=srRsdl, bounds=bounds, Btrig=0.5, Ftar=0.75, maxHR=0.35, aLag=1, srvBias=1, cthBias=1, IEM="linear", slag=5, clag=5, am="srv")

base.cth <- mseBD(OM=OM, start=iniyr, sr=srBH, srRsdl=srRsdl, bounds=bounds, Btrig=0.5, Ftar=0.75, maxHR=0.35, aLag=1, srvBias=1, cthBias=1, IEM="linear", clag=5, am="cth")

base.maxHR1 <- mseBD(OM=OM, start=iniyr, sr=srBH, srRsdl=srRsdl, bounds=bounds, Btrig=0.5, Ftar=0.75, maxHR=1, aLag=1, srvBias=1, cthBias=1, IEM="linear")

base.Btrig1 <- mseBD(OM=OM, start=iniyr, sr=srBH, srRsdl=srRsdl, bounds=bounds, Btrig=1, Ftar=0.75, maxHR=0.35, aLag=1, srvBias=1, cthBias=1, IEM="linear")

base.Ftar1 <- mseBD(OM=OM, start=iniyr, sr=srBH, srRsdl=srRsdl, bounds=bounds, Btrig=0.5, Ftar=1, maxHR=0.35, aLag=1, srvBias=1, cthBias=1, IEM="linear")

bounds["b0",   c("phase","start")]=c(-1,0.5)

base.b0.5 <- mseBD(OM=OM, start=iniyr, sr=srBH, srRsdl=srRsdl, bounds=bounds, Btrig=0.5, Ftar=0.75, maxHR=0.35, aLag=1, srvBias=1, cthBias=1, IEM="linear")

bounds["b0",   c("phase","start")]=c(-1,0.2)

base.srvslag1 <- mseBD(OM=OM, start=iniyr, sr=srBH, srRsdl=srRsdl, bounds=bounds, Btrig=0.5, Ftar=0.75, maxHR=0.35, aLag=1, srvBias=1, cthBias=1, IEM="linear", slag=1, clag=5, am="srv")

base.srvslag5alag5 <- mseBD(OM=OM, start=iniyr, sr=srBH, srRsdl=srRsdl, bounds=bounds, Btrig=0.5, Ftar=0.75, maxHR=0.35, aLag=5, srvBias=1, cthBias=1, IEM="linear", slag=5, clag=5, am="srv")

base.srvslag5alag5srvBias0.5 <- mseBD(OM=OM, start=iniyr, sr=srBH, srRsdl=srRsdl, bounds=bounds, Btrig=0.5, Ftar=0.75, maxHR=0.35, aLag=5, srvBias=0.5, cthBias=1, IEM="linear", slag=5, clag=5, am="srv")

base.iemCts <- mseBD(OM=OM, start=iniyr, sr=srBH, srRsdl=srRsdl, bounds=bounds, Btrig=0.5, Ftar=0.75, maxHR=0.35, aLag=1, srvBias=1, cthBias=0.5, IEM="cst")

lst <- list(base, base.aLag3, base.aLag5, base.cthBias0.5, base.srvBias0.5, base.srvBias1.5, base.srv, base.cth, base.maxHR1, base.Btrig1, base.Ftar1, base.b0.5, base.srvslag1, base.srvslag5alag5, base.srvslag5alag5srvBias0.5, base.srvSim, base.iemCts)

base.summ <- mseSumm(lst, scn)

save(base, base.aLag3, base.aLag5, base.b0.5, base.Btrig1, base.cth, base.cthBias0.5, base.Ftar1, base.iemCts, base.maxHR1, base.srv, base.srvBias0.5, base.srvBias1.5, base.srvSim, base.srvslag1, base.srvslag5alag5, base.srvslag5alag5srvBias0.5, base.iemCts, base.summ, scn, file="../report/RData.mse2")

q("no")



