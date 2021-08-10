#-----------------------------------------------------------
#XPSConstrain: make constraints among fitting parameters and 
#stores informatin into the XPSSample slots
#-----------------------------------------------------------

# Eg:
# constraint on mu: mu2=mu1+1.5
# prova[[2]] <- XPSConstrain(prova[[2]], 2, "link", "mu2", expr="mu1+1.5")
# constraint on sigma: sigma2=sigma1
# prova[[2]] <- XPSConstrain(prova[[2]], 3, "link", "sigma2", expr="sigma1")

#'XPSConstrain function to set fit constraints
#'
#'@param Object XPSSample
#'@param ncomponent the fit component index
#'@param action  link, fix, set, remove the constraint to apply or remove
#'@param variable h, mu, sigma, etc.   the fit parameter to constrain
#'@param parameter one of the values start, min, max
#'@param value the value to set
#'@param expr the constraint expression
#'
#'@examples
#'\dontrun{
#' XPSConstrain(XPSCoreline[[3]], ncomponent=2, action="link", variable="mu", expr="mu1+1.5")
#' XPSConstrain(XPSCoreline[[2]], ncomponent=3, action="link", variable="sigma", expr="sigma1")
#'}
#'
#'@export
#'

XPSConstrain <- function(Object,         # XPSSample[[CoreLine]]
                         ncomponent=NA,  # the fit component number
                         action=c("show", "link", "fix", "set", "remove"), # the kind of contrait to apply
                         variable=NULL,  # one of the name of parameters eg: "h" "mu" "sigma" ...
                         parameter=NULL, # one of c("start","min", "max")
                         value=NA,       # value to set
                         expr )          # like "mu1 + 1.5" or "h1 / 2", "h1 * 0.5", "sigma1"
   {
      action <- match.arg(action)
      switch(action,
          "show" = {
             if ( is.na(ncomponent)) ncomponent <- seq_along(Object@Components)
             for (n in ncomponent ) {
                  cat("\n ----------------------------------------------------\n")
                  cat(" Component : ", slot(Object@Components[[n]],"label")," ", slot(Object@Components[[n]],"funcName"),"\n" )
                  cat(sprintf("%10s _ %25s \n"," ","parameter" ))
                  cat(sprintf("%10s | %12s  %12s  %12s \n","variable","start", "min", "max"))
                  for (par in seq_len(nrow(Object@Components[[n]]@param))) {
                       cat(sprintf("%10s : %12.5f  %12.5f  %12.5f \n", rownames(Object@Components[[n]]@param)[par], Object@Components[[n]]@param[par,"start"], Object@Components[[n]]@param[par,"min"], Object@Components[[n]]@param[par,"max"]) )
                  }
             }
          },
          "link" = {
		           ## requires: variable, expr., FUN, value
		           # increase the number of components link
	            if ( is.na(ncomponent) || is.null(variable) ) {
		               stop(" the component number 'ncomponent' as well as 'variable' are required.\n")
        	       }

		           num <- length(Object@Components[[ncomponent]]@link) + 1  #if other links are present, add the new one
		           var <- paste(variable,ncomponent, sep="")
		           expr <- gsub("\ ","",expr) ## remove blank spaces

		           ## get the position of the variable in the list of parameters
		           paramnames <- unlist(lapply(seq_along(Object@Components), function(x,y) paste(rownames(y[[x]]@param),x,sep=""), y = Object@Components ))
		           position <- charmatch(var, paramnames)
             ## set list for default case
		           Object@Components[[ncomponent]]@link[[num]] <- list(action="link",
		                                                               variable=var,
                                                                 expr=expr,
		                                                               FUN=NA,
		                                                               value=NA,
		                                                               newvar=NA,
		                                                               position=position)
		           funpos <- regexpr("[[:punct:]]",expr)
		           ## detect if an operation is present between linked parameters es: "mu1+1.5"
		           if (attr(funpos,"match.length") > 0) {
		              FUN <- regmatches(expr,funpos) ## [1] "+"
		              exprsplit <- unlist(strsplit(expr,FUN, fixed=TRUE)) ## [1] "mu1" "1.5"
		              Object@Components[[ncomponent]]@link[[num]]$FUN <- FUN
		              Object@Components[[ncomponent]]@link[[num]]$newvar <- exprsplit[1] # mu1
		              Object@Components[[ncomponent]]@link[[num]]$value <- as.numeric(exprsplit[2]) # 1.5
		           }
		        },
		        "fix" = {
		  	        if ( is.na(ncomponent) ) {
                stop(" the component number 'ncomponent' is required.\n")
             }
			          if ( is.null(variable) ) {
                stop(" 'variable' is required.\n")
             }
             # variable <- match.arg(variable, rownames(slot(Object@Components[[ncomponent]],"param")))
			          if ( is.na(value) ) {
			             value <- getParam(Object@Components[[ncomponent]],"start", variable)
             }
   			       ## now set c("start","min", "max") at the same value
   	         value <- rep(value,3)
   			       Object@Components[[ncomponent]] <- setParam(Object@Components[[ncomponent]],
                                                         parameter=NULL,
                                                         variable=variable,
                                                         value=value)
		           },
		        "set" = {
		           if ( length(ncomponent)==0 ) stop(" the component number 'ncomponent' is required.\n")
		           if ( is.null(variable) || is.na(value) ) {
                     stop(" 'variable to set' and 'value' are required.\n")
                 }
			          Object@Components[[ncomponent]] <- setParam(Object@Components[[ncomponent]],
											                                              parameter=parameter,
												                                             variable=variable,
												                                             value=value)
		        },
		        "remove" = {
		           if ( length(ncomponent)==0 ) {
                    cat("\n Please specify the Fit component!\n")
                    return(Object)
                 } else {
      	             Object@Components[[ncomponent]]@link <- list()
                    value <- Object@Components[[ncomponent]]@param[1,1] # h value
	                   Object@Components[[ncomponent]]@param[1,] <- c(value, 0, 5*value)
      	             value <- Object@Components[[ncomponent]]@param[2,1] # mu value
			                 Object@Components[[ncomponent]]@param[2,] <- c(value, value-1.5, value+1.5)
      	             value <- Object@Components[[ncomponent]]@param[3,1] # sigma value
		         	        Object@Components[[ncomponent]]@param[3,] <- c(value, 0.01, 10)
                 }
             }
	   )
	   return(invisible(Object))
}
