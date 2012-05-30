mseBD <- function (OM, start, sr, srRsdl = FLQuant(1, dimnames = dimnames(window(rec(OM),
    start = start))), CV = 0.2, Ftar = 0.75, Btrig = 0.75, Fmin = Ftar *
    0.1, Blim = Btrig * 1e-04, Bpct = 0.5, Fpct = 0.5, jk = FALSE,
    bounds = NULL, aLag=1, maxHR=0.35, maxF=2, IEM="linear", cthBias=1, srvBias=1, seed=1234){

	#--------------------------------------------------------------------
	# set year's info  
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
	#	maxF=max(fbar(OM))*1.1
	#--------------------------------------------------------------------


	#--------------------------------------------------------------------
	# general settings
	#--------------------------------------------------------------------
	# define tempdir for admb to use for paralelization
	bdir <- paste(tempdir(), runif(1,1,2), sep="/")
	system(paste("mkdir", bdir, sep=" "))
	# other variables required
	set.seed(seed)
	iniPyr <- start
	lastPyr <- OM@range["maxyear"]
	nPyr <- lastPyr-iniPyr+1
	# assessment years
	aYrs <- seq(iniPyr, lastPyr, aLag)

	#--------------------------------------------------------------------
	# introduce OEM on historical data
	#--------------------------------------------------------------------
	bd <- as(OM, "FLBioDym") # <<< this bit is required
      # abundance index variability
      index(bd) <- index(bd)*rlnorm(prod(dim(index(bd))), 0, CV)
      # abundance index bias
	index(bd) <- index(bd)*runif(length(index(bd)), srvBias*0.95, srvBias*1.05)
	# catch bias
	catch(bd) <- catch(bd)*runif(length(catch(bd)), cthBias*0.95, cthBias*1.05)
	# add bounds for BioDym	
    	if (!is.null(bounds)) bd@bounds <- bounds

	#--------------------------------------------------------------------
	# IEM setup
	#--------------------------------------------------------------------
	cmin <- min(catch(OM))
	cmax <- max(catch(OM))
	#iem <- lm(iem~catch, data=data.frame(catch=c(cmin,cmax),iem=c(cthBias,1)))

	#--------------------------------------------------------------------
	# object to register TACs, start year catch taken from OM 
	#--------------------------------------------------------------------
	PARrec <- catch(OM)
	PARrec <- FLCore::expand(PARrec, age=c("all", "HCR:hr", paste("BDM", dimnames(params(bd))$params), sep=":")) # hack to overcome exported method by reshape
	PARrec["hr"] <- catch(OM)/stock(OM) 
	dimnames(PARrec)[[1]][1] <- "HCR:TAC"	

	#--------------------------------------------------------------------
	# Go !!
	#--------------------------------------------------------------------
    for (iYr in iniPyr:(lastPyr-aLag)) {
        cat("===================", iYr, "===================\n")
		# Is it assessment year ? if yes do your magic !
		if(iYr %in% aYrs){
			dtaYr <- ac(iYr-1)
			dtaYrs <- ac((iYr-aLag):(iYr-1))
			advYrs <- ac((iYr+1):(iYr+aLag)) 

			# OEM	
	        bd <- window(bd, end = an(dtaYr))
	        # abundance index variability
	        index(bd)[,dtaYrs] <- stock(OM)[,dtaYrs]*rlnorm(prod(dim(index(bd)[,dtaYrs])), 0, CV)
	        # abundance index bias
       	  index(bd)[,dtaYrs] <- index(bd)[,dtaYrs]*runif(length(index(bd)[,dtaYrs]), srvBias*0.95, srvBias*1.05)
		  # catch bias
	        catch(bd)[,dtaYrs] <- computeCatch(OM)[, dtaYrs]
	        catch(bd)[,dtaYrs] <- catch(bd)[,dtaYrs]*runif(length(catch(bd)[,dtaYrs]), cthBias*0.95, cthBias*1.05)

            # MP
            # assessment
            bd <- admbBD(bd, dir=bdir)
            PARrec[-c(1,2), ac(iYr)] <- c(params(bd)) 

			# projections considering TAC was caught in iYr
			# to deal with a missing feature in fwd the last year of data must also be included
			ct <- PARrec["TAC",c(dtaYr, iYr)]
			ct[,dtaYr] <- catch(bd)[,dtaYr]
			bd <- fwd(bd, catch=ct)
			# HCR
		    if (jk) {
		        #hv <- hcrJK(bd, Ftar, Btrig, Fmin, Blim, Fpct, Bpct)
		        warning("Using jk hcr requires checking.", immediate.=T)
		    } else {
		    	# the lag on the hcr method is between the advice year and the data year, so it's 2
		        hv <- hcr(bd, FLPar(Ftar = Ftar, Btrig = Btrig, Fmin = Fmin, Blim = Blim), lag=2)
				# check F is not above maxF and replace if it is
				hv[hv>maxHR] <- maxHR
				PARrec["hr", ac(iYr)] <- hv
				tac <- TAC(bd, hv) # this is just hv*b
		    }
			# update TAC record, all advYrs get the same TAC, may need more options
			PARrec["TAC",advYrs] <- tac[,rep(1, length(advYrs))]
			
		    # IEM
		    ctrl <- fwdControl(data.frame(year = an(c(dtaYr, iYr, iYr, rep(advYrs, rep(2, aLag)))), max = NA, quantity = c("catch", rep(c("catch","f"), aLag+1))))
		    dms <- dimnames(ctrl@trgtArray)
		    dms$iter <- 1:nits
		    ctrl@trgtArray <- array(NA, lapply(dms, length), dms)
		    ctrl@trgtArray[1, "val", ] <- catch(OM)[,ac(dtaYr)]
		    ctrl@trgtArray[2*(1:(aLag+1)), "val", ] <- PARrec["TAC",c(iYr, advYrs)]
			# catch will be corrected by the bias level introduced on the OEM ?!
			# There may be better options like linking through effort

			ie <- ctrl@trgtArray[2*(1:(aLag+1)), "val", ]
			if(IEM=="linear" & cthBias!=1){
				ie[ie<cmin & !is.na(ie)] <- cthBias
				ie[ie>cmax & !is.na(ie)] <- 1
				#ie[ie>=cmin & ie<=cmax & !is.na(ie)] <- predict(iem, newdata=data.frame(catch=ie[ie>=cmin & ie<=cmax & !is.na(ie)]))
				ie[ie>=cmin & ie<=cmax & !is.na(ie)] <- cthBias + (ie[ie>=cmin & ie<=cmax & !is.na(ie)] - cmin)*(1-cthBias)/(cmax-cmin)
			} else {
				ie[] <- cthBias			
			}
			ctrl@trgtArray[2*(1:(aLag+1)), "val", ] <- ctrl@trgtArray[2*(1:(aLag+1)), "val", ]/runif(length(ie), ie*0.95, ie*1.05)

			# limit on F to stabilize the simulations
		    ctrl@trgtArray[2*(1:(aLag+1))+1, "max", ] <- maxF

			# OM
		    OM <- fwd(OM, ctrl = ctrl, sr = sr, sr.residuals = srRsdl) 

		} else {
			#tac <- TACrec[,c(iYr+1)]
	        cat("=============== No assessment ==============\n")
		}
    }
	attr(OM, "PARs") <- PARrec
    return(OM)
}


# OTHER FUNCTIONS

setGeneric("an", function(object) standardGeneric("an"))
setMethod("an", "character", function(object) as.numeric(object))
setMethod("an", "factor", function(object) as.numeric(as.character(object)))

setMethod("simAR", "FLQuant", function (object, nits, start = dims(object)$minyear, end = dims(object)$maxyear, n = 1, cv=sd(object)) 
    {
        arobject = acf(c(object), lag.max = n)
        plot(arobject)
        cvobject = cv
        objectRsdl = window(object, end = end)
        objectRsdl[] = 0
        objectRsdl = log(rlnorm(nits, objectRsdl, cvobject))
        for (i in end:(dims(objectRsdl)$minyear + n)) objectRsdl[, 
            ac(i)] = apply(sweep(objectRsdl[, rev(ac(i - 1:n))], 
            2, arobject$acf[n:1], "*"), 6, sum)
        return(objectRsdl)
    })

setGeneric("mseSumm", function(object, ...) standardGeneric("mseSumm"))
setMethod("mseSumm", c("list"), function(object, scn, statSumm=TRUE, probs=c(0.2,0.5,0.8)){
	lst <- split(scn, scn)
	for(i in 1:nrow(scn)){
		x <- object[[i]]
		df0 <- as.data.frame(fbar(x))
		df0[,1] <- "fbar"
		df1 <- as.data.frame(ssb(x))
		df1[,1] <- "ssb"
		df2 <- as.data.frame(rec(x))
		df2[,1] <- "rec"
		df3 <- as.data.frame(catch(x))
		df3[,1] <- "catch"
		df4 <- as.data.frame(attr(x, "PARs"))
		out <- rbind(df0, df1, df2, df3, df4)
		out <- cbind(out, scn[i,][rep(1,nrow(out)),])
		out[,"runid"] <- i
		out[,"scn"] <- paste(paste(names(scn), scn[i,], sep="="), collapse="; ")
		names(out)[1] <- "par"
		if(statSumm==TRUE){
			out <- lapply(split(out, out[,c("par","year")]),function(x){
				xx <- x[1:3,]
				xx$qtl <- probs
				xx$data <- quantile(x$data, probs=probs, na.rm=TRUE)
				xx
			})
			out <- do.call("rbind", out)
		}
		lst[[i]] <- out
	}
	lst <- do.call("rbind", lst)
	row.names(lst) <- NULL
	lst
})



