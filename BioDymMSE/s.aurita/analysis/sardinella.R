# explore.R - DESC
# explore.R

# Copyright 2003-2012 FLR Team. Distributed under the GPL 2 or later
# Maintainer: Iago Mosqueira, JRC
# $Id: $

# http://www.fao.org/fishery/rfb/cecaf/en


library(FLAdvice)
library(FLBioDym)
library(FLAssess)
library(plyr)

# Sardinella data
nc <- read.table('data/nc_sar.dat', header=TRUE)
ia <- read.table('data/cpue_sar.dat', header=TRUE)

ca <- FLQuant(nc$catch, dimnames=list(age='all', year=nc$year))
cp <- FLQuant(c(rep(NA, 7), ia$cpue), dimnames=list(age='all', year=nc$year))

bd <- FLBioDym(catch=ca, index=cp)

# bounds & start
# r
bd@bounds[1,] <- c(1, 0.32, 1.28, 0.64)
# K
bd@bounds[2,] <- c(2, 1750, 7000, 2500)
# p
bd@bounds[3,] <- c(-1, 1, 1, 1)
# b0
bd@bounds[4,] <- c(-1, 0.25, 0.95, 0.5)
# q
bd@bounds[5,] <- c(1, 0.1, 10, 1)
# sigma
bd@bounds[6,] <- c(1, 0.01, 1, 0.1)

# admbBD assessment
res <- admbBD(bd)

# LH parameters from FishBase & FLBioDym
par <- as(data.frame(linf=28.5, k=0.4, t0=-0.1, s=0.8, v=4500, a50=1), 'FLPar')

# Use CECAF SA B0 as starting point
par <- as(data.frame(linf=28.5, k=0.8, t0=-0.1, s=0.8, v=1750, a50=1), 'FLPar')

# gislasim
sar <- lh(gislasim(par), range=c(min=1, max=8, minfbar=1, maxfbar=6))

# stk with initial F closest to estimated HR
stk <- as(sar, 'FLStock')
stk <- stk[,12]
dimnames(stk) <- list(year=1989)

# prepare for projection
saa <- stf(stk, 21, 1)

# projection control
trg <- fwdControl(data.frame(year=1990:2010, val=c(ca), quantity="catch"))

# catch projection
saa <- fwd(saa, trg, sr=list(model="bevholt", params=params(sar)))

name(saa) <- "SAA"
desc(saa) <- "Simulated Sardinella aurita, partly conditioned on CECAF SA 2010"

save(saa, file=paste("SAA_", format(Sys.time(), "%Y%m%d_%H%M"), ".RData", sep=""))

# CECAF SA results
sa <- read.table('data/sa.dat', sep='\t', header=T)
