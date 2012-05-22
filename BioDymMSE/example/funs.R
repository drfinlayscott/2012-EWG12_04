mseBD <- function (OM, start, sr, srRsdl = FLQuant(1, dimnames = dimnames(window(rec(OM),
    start = start))), CV = 0.3, Ftar = 0.75, Btrig = 0.75, Fmin = Ftar *
    0.1, Blim = Btrig * 1e-04, Bpct = 0.5, Fpct = 0.5, jk = FALSE,
    bounds = NULL, aLag=1, maxF=0.35){

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
	#--------------------------------------------------------------------

	iniPyr <- start
	lastPyr <- OM@range["maxyear"]
	nPyr <- lastPyr-iniPyr+1

	#--------------------------------------------------------------------
	# not so sure this block of code is needed ... to be checked

	# deal with iterations
    nits <- c(OM = dims(OM)$iter, sr = dims(params(sr))$iter, rsdl = dims(srRsdl)$iter)
    if (length(unique(nits)) >= 2 & !(1 %in% nits)) ("Stop, iters not '1 or n' in OM")
    nits <- max(nits)

	#stock(OM) = propagate(stock(OM), nits)
 	if(dims(OM)$iter!=nits) stock(OM) <- propagate(stock(OM), nits)

    bd <- as(OM, "FLBioDym") # <<< this bit is required
    if (!is.null(bounds)) bd@bounds <- bounds
    #bd = propagate(bd, nits)
    index(bd) <- index(bd) * rlnorm(prod(dim(index(bd))), 0, CV) # ME NOT LIKE !!
	#--------------------------------------------------------------------

	# object to register TACs, start year catch taken from OM 
	TACrec <- catch(bd)

	# assessment years
	aYrs <- seq(iniPyr, lastPyr, aLag)

	# Go !!
    for (iYr in iniPyr:(lastPyr-aLag)) {
        cat("===================", iYr, "===================\n")
		# Is it assessment year ? if yes do your magic !
		if(iYr %in% aYrs){
			dtaYr <- ac(iYr-1)
			dtaYrs <- ac((iYr-aLag):(iYr-1))
			advYrs <- ac((iYr+1):(iYr+aLag)) 
			# OEM
	        bd <- window(bd, end = an(dtaYr))
	        index(bd)[,dtaYrs] <- stock(OM)[,dtaYrs] * rlnorm(prod(dim(index(bd)[,dtaYrs])), 0, CV)
	        catch(bd)[,dtaYrs] <- computeCatch(OM)[, dtaYrs]
            # assessment
            bd <- admbBD(bd)
			# projections considering TAC was caught in iYr
			# to deal with a missing feature in fwd the last year of data must also be included
			ct <- TACrec[,c(dtaYr, iYr)]
			ct[,dtaYr] <- catch(bd)[,dtaYr]
			bd <- fwd(bd, catch=ct)
			# HCR
		    if (jk) {
		        #hv <- hcrJK(bd, Ftar, Btrig, Fmin, Blim, Fpct, Bpct)
		        warning("Using jk hcr requires checking.", immediate.=T)
		    } else {
		    	# the lag on the hcr method is between the advice year and the data year,
		    	# so in the traditional "assessments evry year" it's 2
		        hv <- hcr(bd, FLPar(Ftar = Ftar, Btrig = Btrig, Fmin = Fmin, Blim = Blim), lag=2)
				tac <- TAC(bd, hv) # this is just hv*b
		    }
			# update TAC record, all advYrs get the same TAC, may need more options 
			TACrec[,advYrs] <- c(tac)

		    # Project OM
		    
		    ctrl <- fwdControl(data.frame(year = an(c(dtaYr, iYr, iYr, rep(advYrs, rep(2, aLag)))), max = NA, quantity = c("catch", rep(c("catch","f"), aLag+1))))
		    dms <- dimnames(ctrl@trgtArray)
		    dms$iter <- 1:nits
		    ctrl@trgtArray <- array(NA, lapply(dms, length), dms)
			# good place to introduce implementation error ??!!
		    ctrl@trgtArray[1, "val", ] <- catch(OM)[,ac(dtaYr)]
		    ctrl@trgtArray[2*(1:(aLag+1)), "val", ] <- TACrec[,c(iYr, advYrs)]
		    ctrl@trgtArray[2*(1:(aLag+1))+1, "max", ] <- maxF
			# OM projection, need to check if the residuals are being allocated correcly
		    OM <- fwd(OM, ctrl = ctrl, sr = sr, sr.residuals = srRsdl) 

		} else {
			#tac <- TACrec[,c(iYr+1)]
	        cat("=============== No assessment ==============\n")
		}
    }
	attr(OM, "TAC") <- TACrec
    return(OM)
}


# OTHER FUNCTIONS

setGeneric("an", function(object) standardGeneric("an"))
setMethod("an", "character", function(object) as.numeric(object))
setMethod("an", "factor", function(object) as.numeric(as.character(object)))






