### $Id: readRPL.R,v 0.0.1
###
###             Read Genplot RPL data
###

#'read.vamas
#'
#'Read Genplot RPL data
#'
#'Read Genplot analyzed spectra in RPL format
#'a RPL folder containing analyzed data must be present
#'
#'@param Object an XPSSample object where to put the read data
#'@param file File
#'@param debug debug switch
#'
#' 
#'@export
#'


readGenplot <- function(Object, file=NULL, debug=FALSE) {
   get_short  <- function(con, POS=0, origin="start") {
            seek(con, where=POS, origin=origin)
            return(readBin(con,"int",1,2))
   }
        
   get_record_lenght  <- function(con, POS=0, origin="start") {
            seek(con, where=POS, origin=origin)
            return(readBin(con,"int",1,4))
  	}

	if ( is.null(file) ) {
      if ( "tcltk" %in% search() ) { 
         file <- tk_choose.files() 
      } else { 
        file <- file.choose() 
      }
	}
	fp <- file(file, open="rb")
	on.exit(close(fp))

#--- HEADER
	reclen <- get_record_lenght(fp,0) ## Version
	reclen <- get_record_lenght(fp, reclen-4, origin="current")# text
	reclen <- get_record_lenght(fp, reclen-4, origin="current")# text
	NPT <- get_short(fp, origin="current")
	Ncol <- get_short(fp, origin="current")
	reclen <- get_record_lenght(fp, origin="current")# text

#--- DATA
	## loop on columns
	for (ii in 1:Ncol) {
       reclen <- get_record_lenght(fp, origin="current") # Column descriptor Data column 1
       reclen <- get_record_lenght(fp, reclen-4, origin="current")# text
       Nbyte <- get_record_lenght(fp, reclen-4, origin="current")

       Y <- readBin(fp, "double", n=NPT, size=4) # read NPT data
       reclen <- get_record_lenght(fp, origin="current") ## end of data

       ## Transform data into list
       Parameters <- list()
       NMaxC <- Y[3]
       Parameters$Boundaries <- list(x=c(Y[1], Y[2]), y=c(Y[5], Y[6]))
		     Parameters$NMaxC <- Y[3]
		     Parameters$Shift <- Y[4]
		     Parameters$BGND <- Y[7]
		     Parameters$ncomp <- Y[8]
		     ncomp <- Y[8]
		     Parameters$Parm <-  Y[seq_len(NMaxC+1)]

		     if (ncomp > 0 ) {
            sidx <- NMaxC + 2
            idx <- seq(from=sidx,length.out=NMaxC)
            Parameters$Xcenter  <- Y[idx][1:ncomp]
            sidx <- sidx + NMaxC
            idx <- seq(from=sidx,length.out=NMaxC)
            Parameters$XcenterRL <- Y[idx][1:ncomp]
            sidx <- sidx + NMaxC
            idx <- seq(from=sidx,length.out=NMaxC)
            Parameters$XcenterRH <- Y[idx][1:ncomp]
            sidx <- sidx + NMaxC
            idx <- seq(from=sidx,length.out=NMaxC)
            Parameters$XcenterMD <- Y[idx][1:ncomp]
            sidx <- sidx + NMaxC
            idx <- seq(from=sidx,length.out=NMaxC)
            Parameters$XcenterLK <- Y[idx][1:ncomp]
            sidx <- sidx + NMaxC
            idx <- seq(from=sidx,length.out=NMaxC)
            Parameters$XcenterRT <- Y[idx][1:ncomp]

            sidx <- sidx + NMaxC
            idx <- seq(from=sidx,length.out=NMaxC)
            Parameters$Sigma <- Y[idx][1:ncomp]
            sidx <- sidx + NMaxC
            idx <- seq(from=sidx,length.out=NMaxC)
            Parameters$SigmaRL <- Y[idx][1:ncomp]
            sidx <- sidx + NMaxC
            idx <- seq(from=sidx,length.out=NMaxC)
            Parameters$SigmaRH <- Y[idx][1:ncomp]
            sidx <- sidx + NMaxC
            idx <- seq(from=sidx,length.out=NMaxC)
            Parameters$SigmaMD <- Y[idx][1:ncomp]
            sidx <- sidx + NMaxC
            idx <- seq(from=sidx,length.out=NMaxC)
            Parameters$SigmaLK <- Y[idx][1:ncomp]
            sidx <- sidx + NMaxC
            idx <- seq(from=sidx,length.out=NMaxC)
            Parameters$SigmaRT <- Y[idx][1:ncomp]

            sidx <- sidx + NMaxC
            idx <- seq(from=sidx,length.out=NMaxC)
            Parameters$Amplitude <- Y[idx][1:ncomp]
            sidx <- sidx + NMaxC
            idx <- seq(from=sidx,length.out=NMaxC)
            Parameters$AmplitudeRL <- Y[idx][1:ncomp]
            sidx <- sidx + NMaxC
            idx <- seq(from=sidx,length.out=NMaxC)
            Parameters$AmplitudeRH <- Y[idx][1:ncomp]
            sidx <- sidx + NMaxC
            idx <- seq(from=sidx,length.out=NMaxC)
            Parameters$AmplitudeMD <- Y[idx][1:ncomp]
            sidx <- sidx + NMaxC
            idx <- seq(from=sidx,length.out=NMaxC)
            Parameters$AmplitudeLK <- Y[idx][1:ncomp]
            sidx <- sidx + NMaxC
            idx <- seq(from=sidx,length.out=NMaxC)
            Parameters$AmplitudeRT <- Y[idx][1:ncomp]

            sidx <- sidx + NMaxC
            idx <- seq(from=sidx,length.out=NMaxC)
            Parameters$MixGL <- Y[idx][1:ncomp]
            sidx <- sidx + NMaxC
            idx <- seq(from=sidx,length.out=NMaxC)
            Parameters$MixGLRL <- Y[idx][1:ncomp]
            sidx <- sidx + NMaxC
            idx <- seq(from=sidx,length.out=NMaxC)
            Parameters$MixGLRH <- Y[idx][1:ncomp]
            sidx <- sidx + NMaxC
            idx <- seq(from=sidx,length.out=NMaxC)
            Parameters$MixGLMD <- Y[idx][1:ncomp]

            sidx <- sidx + NMaxC
            idx <- seq(from=sidx,length.out=NMaxC)
            Parameters$Asym <- Y[idx][1:ncomp]
            sidx <- sidx + NMaxC
            idx <- seq(from=sidx,length.out=NMaxC)
            Parameters$AsymRL <- Y[idx][1:ncomp]
            sidx <- sidx + NMaxC
            idx <- seq(from=sidx,length.out=NMaxC)
            Parameters$AsymRH <- Y[idx][1:ncomp]
            sidx <- sidx + NMaxC
            idx <- seq(from=sidx,length.out=NMaxC)
            Parameters$AsymMD <- Y[idx][1:ncomp]

            sidx <- sidx + NMaxC
            idx <- seq(from=sidx,length.out=NMaxC+4)
            Parameters$IntegralComponents <- Y[idx][1:ncomp]
            Parameters$Integral_RTF_BGND_SUM <- Y[idx][(NMaxC+1):(NMaxC+3)]

            sidx <- sidx + NMaxC + 4
            idx <- seq(from=sidx,length.out=NMaxC)
            Parameters$RSF <- Y[idx][1:ncomp]
		     }

		     if ( debug ) {
            cat("\n NPT = ",NPT)
            cat("\n N. column = ", ii)
            cat("\n N. data byte = ", Nbyte)
            cat("\n Parm : ",Parameters$Parm)
            cat("\n RSF  : ",Parameters$RSF)
            cat("\n ncomp : ",ncomp)
#            if ( "IntegralComponents" %in% names(Parameters) ) {
#               cat("\n ", Parameters$IntegralComponents )
#               cat("\n ", Parameters$Integral_RTF_BGND_SUM )
#            }
       }

#--- If coreline was analyzed ncomp==Y[8] > 1

       object <- Object[[ii]]
       if (object@Symbol %in% c("survey", "Survey")){
        	  next
       }
       if ( ncomp > -1 ) {
        	  if (Parameters$Shift != 0) {
        	      Parameters$Boundaries$x <- Parameters$Boundaries$x - Parameters$Shift
           }
           slot(object,"Boundaries") <- Parameters$Boundaries
           # set RegionToFit
           object <- XPSsetRegionToFit(object)
           # set RSF
           if ( length(Parameters$RSF) != 0 ) {
              slot(object,"RSF") <- Parameters$RSF[1] }
        	  # set baseline
        	  if ( Parameters$BGND == 1 ) {
              object <- XPSbaseline(object, bgtype="linear", BgDeg=NULL, Wgt=NULL, splinePoints=NULL)
           } else {
              object <- XPSbaseline(object, bgtype="Shirley", BgDeg=NULL, Wgt=NULL, splinePoints=NULL)
           }
        	  # set SHIFT: in old Genplot analysis Shift is applied also to the Boundaries
       	  if ( Parameters$Shift != 0 ) {
              # shift for original data: mod 7 Sept2015
              object[[1]] <- object[[1]] + Parameters$Shift
              # RegionToFit & Baseline
              if ( length(object@RegionToFit) ) object@RegionToFit$x <- object@RegionToFit$x + Parameters$Shift
              if ( length(object@Baseline) ) object@Baseline$x <- object@Baseline$x + Parameters$Shift
              slot(object,"Shift") <- Parameters$Shift
           }

#--- Fit components
        	  if (ncomp > 0) {
       	  	    ## add components
       		      for(cindex in 1:ncomp ) {
        			        ## function depends on asym
        			        if ((Parameters$Asym[cindex] != 0) || (Parameters$AsymMD[cindex] != 0) ) {
        				           funct <- "AsymmGaussLorentz" }
        			        else {
                       funct <- "GaussLorentzSum"
                   }
            					  ## add the new component
            					  object@Components[[cindex]] <- fitAlgorithms[[funct]]
            					  slot(object@Components[[cindex]],"rsf") <- Parameters$RSF[1]
            					  ## start h, mu, sigma, mix, asym
            					  startpar <- c(Parameters$Amplitude[cindex],
         								                Parameters$Xcenter[cindex],
         				              				  Parameters$Sigma[cindex],
         			              					  Parameters$MixGL[cindex]
  								                      )
  					            lowbond <- c(Parameters$AmplitudeRL[cindex],
  				                    				  Parameters$XcenterRL[cindex],
  			                    					  Parameters$SigmaRL[cindex],
  		                    						  Parameters$MixGLRL[cindex]
  	                    						  )
            					  highbond <- c(Parameters$AmplitudeRH[cindex],
                       								  Parameters$XcenterRH[cindex],
                       								  Parameters$SigmaRH[cindex],
                       								  Parameters$MixGLRH[cindex]
                      								  )

            					  if (funct == "AsymmGaussLorentz" && Parameters$Asym[cindex] != 0) {
	               					  startpar <- c(startpar, Parameters$Asym[cindex] )
						                 lowbond  <- c(lowbond, Parameters$AsymRL[cindex])
					               	  highbond <- c(highbond, Parameters$AsymRH[cindex])
             			   }
					              ## Xcenter == Start
	             	    slot(object,"Components")[[cindex]] <- setParam(object@Components[[cindex]],
          							                                                		parameter="start",
					         		                                                			value = startpar
										                                                        )

					              slot(object,"Components")[[cindex]] <- setParam(object@Components[[cindex]],
							                                                         			parameter="min",
										                                                         value = lowbond
										                                                        )

					              slot(object,"Components")[[cindex]] <- setParam(object@Components[[cindex]],
								                                                         		parameter="max",
										                                                         value = highbond
										                                                        )

					              # label (used in the plot) # #1, #2, ...
					              slot(object@Components[[cindex]],"label") <- paste("#",as.character(cindex),sep="")
					              # name of this component
				            	  names(object@Components)[cindex] <- paste("C",as.character(cindex),sep="") # C1, C2, ...
			            		  ## add the y values for this component
			            		  object@Components[[cindex]] <- Ycomponent(object@Components[[cindex]], x = object@RegionToFit$x, y = object@Baseline$y)
		             		  ## set RSF
					              if (length(Parameters$RSF) > 1 && Parameters$RSF[cindex] != 0 ) {
					            	     slot(object@Components[[cindex]],"rsf") <- Parameters$RSF[cindex]
                   } else {
     	            			  slot(object@Components[[cindex]],"rsf") <- slot(object,"RSF")
   	            	  }
        		     }# end loop on ncomp
        		     yFit <- do.call("rbind", lapply(object@Components, function(x) x@ycoor))
        		     object@Fit$y <- colSums(yFit) - (ncomp * object@Baseline$y)
        	   }# end if ncomp > 0
        }# end if ncomp > -1
        Object[[ii]] <- object
	 } # end for(ii in 1:Ncol)
  return(Object)
}

