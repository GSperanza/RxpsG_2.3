#Computes standard deviation of an XPS-Coreline. 
#Initially thought to be used for the MaxEnt procedure, but not used in the final MaxEnt() routine.


XPS.STDV <- function(){

        STDV <- 0
        Object <- NULL
        Noise <- NULL  #array containing the noise
        NTrials <- 20 #N synthesized spectra to compute the StDev
        Isynth <- NULL #spectral intensity to compute the StDev
        SyntSpect <- NULL #array to store the Synthesize Spectrum
        if (exists("activeFName", envir=.GlobalEnv)==FALSE){
           gmessage("Please load an XPSSAMPLE", title="XPS-SAMPLE LACKING", icon="error")
           return()
        }
        FName <- get(activeFName, envir=.GlobalEnv)
#W        if (exists("activeSpectName", envir=.GlobalEnv)==FALSE){
           SpectList <- XPSSpectList(activeFName)
           CLmainwin <- gwindow("CORELINE SELECTION", visible=FALSE)
           size(CLmainwin) <- c(200, 100)
           CLgroup <- ggroup(label="", container=CLmainwin)
           CLframe <- gframe(text="SELECT THE CORELINES TO ANALYZE", horizontal=FALSE, spacing=5, container=CLgroup)
           CLobj <- gcheckboxgroup(SpectList, selected=1, horizontal=FALSE, handler=function(h,...){
                           SpectList <<- svalue(CLobj)
                        }, container=CLframe)
           gbutton("Save & Exit", handler=function(h, ...){
                           dispose(CLmainwin)
                        }, container=CLframe)
           visible(CLmainwin)<-TRUE
           CLmainwin$set_modal(TRUE)
#        }
#        for(jj in 1:N.XS){ #this for runs on the N. XPSSample == N. tilt angles
#            Object <- FName[[CLname]]
#            RSF <- FName[[CLname]]@RSF
#        }
        N.CL <- length(SpectList)
        for(ii in 1:N.CL){     #this for runs on the N. CoreLines
                XPS.CL <- unlist(strsplit(SpectList[ii], "\\."))   #drop "NUMber." in component name
                Indx <- as.integer(XPS.CL[1])
                CLname <- XPS.CL[2]
                Object <- FName[[CLname]]
                RSF <- FName[[CLname]]@RSF
                Estep <- abs(FName[[CLname]]@.Data[[1]][2] - FName[[CLname]]@.Data[[1]][1])
                Object@Components <- list()
                Object@Fit <- list()
                NN <- length(Object@RegionToFit$y)
                Rx <- range(Object@RegionToFit$x)   #range returns always Rx[1} < Rx[2]
                NComp <- floor((Rx[2]-Rx[1])/1.5)   #each 2eV a Gaussian component
                for(kk in 1:NComp){
                    xx <- Rx[1]+kk*1.5
                    idx <- findXIndex(Object@RegionToFit$x, xx)
                    yy <- Object@RegionToFit$y[idx]
                    Object <- XPSaddComponent(Object, type = "Gauss", peakPosition = list(x = xx, y = yy) )
                }
                Object <- XPSFitLM(Object, verbose=FALSE)
                plot(Object)
                answ <- gconfirm("Is the peak fitting OK?", title="CONTROL THE FIT QUALITY", icon="info")
                if(answ == FALSE){
                   gmessage("Increase the Number of Fitting Components")
                   NComp <- floor((Rx[2]-Rx[1])/1.1)   #each 1.1eV a Gaussian component
                   for(kk in 1:NComp){
                       xx <- Rx[1]+kk*1.1
                       idx <- findXIndex(Object@RegionToFit$x, xx)
                       yy <- Object@RegionToFit$y[idx]
                       Object <- XPSaddComponent(Object, type = "Gauss", peakPosition = list(x = xx, y = yy) )
                   }
                   Object <- XPSFitLM(Object, verbose=FALSE)
                }
                MaxNoise <- max(Object@Fit$y - Object@RegionToFit$y + Object@Baseline$y) #max of the noise superposed to the coreline
                MinNoise <- min(Object@Fit$y - Object@RegionToFit$y + Object@Baseline$y) #max of the noise superposed to the coreline
                NoiseAmpli <- MaxNoise-MinNoise
                LL <- length(Object@RegionToFit$y)
                I0 <- sum(Object@RegionToFit$y - Object@Baseline$y)*Estep #/RSF
                for (kk in 1:NTrials){
                    Noise <- runif(LL, min=-NoiseAmpli/2, max=NoiseAmpli/2)
                    SyntSpect <- Object@Fit$y + Noise
                    Isynth <- sum(SyntSpect)*Estep #/RSF
                    STDV <- STDV + (Isynth-I0)^2
                }
                STDV <- sqrt( 1/20*STDV)
cat("\n***** Coreline", CLname, "   Area", I0)
cat("\n***** Noise", NoiseAmpli, "   STDV", STDV)
        }
}

