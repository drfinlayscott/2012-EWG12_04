library(FLBioDym)
data(ple4)
bd <- as(ple4, "FLBioDym")
bd <- admbBD(bd)
bd <- fwd(bd, catch=FLQuant(300000, dimnames=list(year=2009)))
plot(bd) # no 2009 information ...

