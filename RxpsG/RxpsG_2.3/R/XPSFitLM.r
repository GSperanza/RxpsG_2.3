## ===========================================================
## 		XPSfit
## call th XPSdofit AND plot the result
## ===========================================================




#---showListParam  prints a list of parameters in a formatted style
#'showListParam prints a list of parameters in a formatted style
#'
#'@param lista list of parameter to be plotted
#'@param decimals numeric value number of decimals to be printed
#'@param ParamNames list of parameter names
#'@return returns the formatted print of parameters
#'
#'@examples
#'
#'\dontrun{
#' XPSfit(SampData[["C1s"]])
#'}
#'
#'@export
#'

showListParam <- function(lista, decimals, ParamNames){
     maxLen <- 15  #max length reserved for list elements
     ParamNames <- unlist(ParamNames)
     ParamNames <- sapply(ParamNames, function(x) encodeString(x, width=15, justify="right")) #encodestring formats each of the parameters of the table
     cat("\n", ParamNames)
     lista <- lapply(lista, function(x) round(x, digits=decimals)) #rounds numbers to 2 decimals
     lista <- lapply(lista, function(x) as.character(x)) #trasform each element of the table in string
     lista <- lapply(lista, function(x) encodeString(x, width=17, justify="right"))
     NN <- length(lista)  #now this is a list composed by a single column of 1 strimg containing the frrmatted parameters
     for (ii in 1:NN){
         lista[[ii]] <- substr(lista[[ii]], 4, nchar(lista[[ii]]))  #skip the first element of the list (component identifier)
         xx <- paste(lista[[ii]], sep="", collapse="")
         xx <- paste("C", ii, ": ", xx, sep="")
         cat("\n", xx)
     }
}





## ===========================================================
## 		XPSFitLM
## Esegue il FIT della curva originale-BGND
## Salva le componenti del Fit
## ===========================================================
#'Perform the fit of a Core Line
#'
#'Calculate the fit of a XPSCoreLine object. The \code{XPSFitLM} performs the
#'fit and stop. This function is the main function for performing the fit. The
#'\code{XPSfit} is the main function which calls the \code{XPSFitLM} and
#'displays the result.
#'
#'
#'@param Object XPSCoreLine object
#'@param plt logical if TRUE residual plot is drawn
#'@param verbose logical by default TRUE to print fit interation information
#'@param \dots  further parameters to the fitting function
#'@return The Object slot \code{FitLM} will be filled with the result of the
#'calculation. All the values of the components will be update and the result
#'will be displyed.
#'@seealso \link{nlsLM}
#'@examples
#'
#'\dontrun{
#' showListParam(c(1271, 285, 1.4), 3, c("h", "mu", "sigma" ) )
#'}
#'
#'@export
#'

XPSFitLM <- function(Object, plt=TRUE, verbose=TRUE, ...) {

    import::here(nls.lm, nls.lm.control, .from=minpack.lm )

 	  # data to fit: curve - baseline
 	  datafit <- data.frame(x = Object@RegionToFit$x,
 						          y = Object@RegionToFit$y - Object@Baseline$y)
	   # for any FunctName sapply ceates the correct list of parameters
	   # next these parameters are modified if link are present
	   xnam <- sapply(names(Object@Components), function(x) {  #sapply runs on the fit components
		      number <- unlist(strsplit(x, "C"))[2]  #index of component
		      fnct <- slot(Object@Components[[x]],"funcName")
		      fnctArgs <- formalArgs(fnct)
		      # the first parameter is x, wont be changed
		      LL <- length(fnctArgs)-1  # LL number of remaining parameters
		      # for each of the remaining parameters add an index for component C1 index=1 ... for component C5 index=5 ecc.
		      ParamIdx <- c("", rep.int(number, LL)) # x unmodified + sequence of indexes ==
		      fparm <- paste(fnctArgs,ParamIdx, sep="")
		      funct <- paste(paste(fnct, "(", sep=""), paste(fparm, collapse=","),")", sep="")

		      #now modify parameters if links are present
		      if ( length(Object@Components[[x]]@link) ) {
   	     		for ( idx in seq_along(slot(Object@Components[[x]],"link"))) {
   	      		  	funct <- sub(Object@Components[[x]]@link[[idx]]$variable, Object@Components[[x]]@link[[idx]]$expr, funct)
   	     		}
   		   }
		      return(funct)
		  })

	   fmla <- as.formula(paste("y ~ ", paste(xnam, collapse= "+")))
    cat("\n----------------- FIT EXPRESSION ----------------\n")
    print(fmla)
#--- print parameters values
	   startpar <- lapply(Object@Components, function(x) getParam(x, parameter="start"))
	   ubound <- lapply(Object@Components, function(x) getParam(x, parameter="max"))
	   lbound <- lapply(Object@Components, function(x) getParam(x, parameter="min"))
	   ParamNames <- lapply(seq_along(Object@Components), function(x,y) paste(rownames(y[[x]]@param),x,sep=""), y = Object@Components )

    cat("\n-------------------------------------------------")
    cat("\n\n START:")
    showListParam(startpar, 3, ParamNames[[1]] )
    cat("\n\n UpLim:")
    showListParam(ubound, 3, ParamNames[[1]] )
    cat("\n\n LwLim:")
    showListParam(lbound, 3, ParamNames[[1]] )
    cat("\n-------------------------------------------------")
    startpar <- unlist(startpar)
    ubound <- unlist(ubound)
	   lbound <- unlist(lbound)
	   ParamNames <- unlist(ParamNames)

	   # parameter names
	   ParamNames <- unlist(lapply(seq_along(Object@Components), function(x,y) paste(rownames(y[[x]]@param),x,sep=""), y = Object@Components ))
	   names(startpar) <- ParamNames
	   names(ubound) <- ParamNames
	   names(lbound) <- ParamNames
	   # backup of start values: this vector will be set with new fit values and considering link conditions
	   tmpfit <- startpar

#--- link: index indicates the param to drop since linked
	   index <- NULL
    # drop the variable from fit paramter list
	   for (number in seq_along(Object@Components)) {
		     if (length(Object@Components[[number]]@link)) {
	   	   		for (idx in seq_along(slot(Object@Components[[number]],"link"))) {
					        index <- c(index,Object@Components[[number]]@link[[idx]]$position)
	   			   }
		     }
	   }
	   if (! is.null(index) ) {
		      startpar <- startpar[-c(index)]
		      ubound <- ubound[-c(index)]
		      lbound <- lbound[-c(index)]
	   }


#--- Fit minpack.LM
    NPrint <- 0
    if (verbose) {NPrint <- 1}
       fit <- NULL
       cat("\n\n >>> Start Fitting")
       fit <- try(XPSnlsLM(     #modified nlsLM: same algorithm but avoids fit processing to block upon NON-convergence
                      formula = fmla,
                      data = datafit,
                      start = startpar,
                      upper = ubound,
                      lower = lbound,
                      algorithm = "LM", #Levenberg-Marquardt
                      trace = FALSE, #prints residual sum of errors at each iteration
                      control = nls.lm.control(ftol = 1e-10, factor = 0.1, maxiter = 500, nprint = NPrint),
                      ... ), silent=TRUE)

      # if not-error update parameters

       if ( ! identical(class(fit), "try-error") )  {
#	        if (show_summary) print(summary(fit)) # summary
#--- fit coefficients
	         cfit <- coef(fit) # this is a list
#--- Update values of tmpfit with values cfit from XPSnlsLM

	         for(nomi in names(cfit)){
             num <- grep(nomi, ParamNames)
             tmpfit[num] <- cfit[nomi]
	         }
	         # if there are links modify the correspondent linked parameters
	         for (idx in seq_along(Object@Components) ) {
             for ( ii in seq_along(Object@Components[[idx]]@link) ) {
                 FUN <- Object@Components[[idx]]@link[[ii]]$FUN
                 if ( ! is.na(FUN) ) {
                    origname <- Object@Components[[idx]]@link[[ii]]$newvar # mu1
                    ## find the value correspondent to the varable in cfit
                    num <- grep(origname, names(cfit))
                    origvalue <- cfit[num]
                    newvalue <- Object@Components[[idx]]@link[[ii]]$value
                    value <- sapply(origvalue,FUN,newvalue)
                    names(value) <- Object@Components[[idx]]@link[[ii]]$variable
                 } else {
                    origname <- Object@Components[[idx]]@link[[ii]]$expr # sigma1
                    ## find the value correspondent to the varable in cfit
                    num <- grep(origname, names(cfit))
                    value <- cfit[num]
                    names(value) <- Object@Components[[idx]]@link[[ii]]$variable
                }
                ## change the value in the position indidcated by $position
                index <- Object@Components[[idx]]@link[[ii]]$position
                tmpfit[index] <- value
             }
	         }

#--- update start values
          IndexNamesOftmpfit <- as.numeric(gsub("[^0-9]", "", names(tmpfit)))
          for (idx in seq_along(Object@Components) ) {
#--- modified 7 Oct 2015 cgange grtep with wich to find indexes and avoid duplicates in the case Ncomp > 9
			           num <- which(idx == IndexNamesOftmpfit) #

              Object@Components[[idx]] <- setParam(Object@Components[[idx]],
                                                    parameter="start",
                                                    value=tmpfit[num]) # start
#--- y values
              Object@Components[[idx]] <- Ycomponent(Object@Components[[idx]],
                                                    x=Object@RegionToFit$x,
                                                    y=Object@Baseline$y)
          }

#--- update slot Fit
   	      Object@Fit$y <- fitted(fit)  #fitted(fit) computes the fitting function using the best fit parameters
	         Object@Fit$fit <- fit
	         if (plt==TRUE){
              XPSresidualPlot(Object)  #see XPSClass.r
          }
       }
       cat("\n\n >>> Fit done!")
       return(Object)
   }
