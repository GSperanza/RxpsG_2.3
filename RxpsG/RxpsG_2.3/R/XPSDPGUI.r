#Depth Profile GUI: ARXPS function   8-5-2020
#Variance of concentrations is made following article:
#K. Harrison, L. B. Hazel: "The Determination of Uncertainties in Quantitative
#  XPS/AES and its Impact on Data Acquisition Strategy", SIA, (1992),18, 368-376

#'Estimate the thickness of a layer by analyzing the trend of core line intensities
#'as a function of the tilt angle.
#'The Signal decay with depth as described by the Drude equation....
#'See Briggs...
#'
#'
#'@examples
#'
#'\dontrun{
#'	DepthPro()
#'}
#'
#'@export
#'

XPSDP <- function(){


   CheckCL <- function(){
       CLtestList <- CommonCL        #duplicate the CommonCL list
       N.XS <- length(CLlist)        #N. XPSSpectra selected
       for (ii in 1:N.XS){
           N.CL <- length(CommonCL)      #N. Coreline names
           for (jj in N.CL:1){
              CL <- unlist(strsplit(CommonCL[jj], "\\.")) #extract pattern to compare skipping the CoreLine index
              CL <- CL[2]
              xx <- grep(CL, CLlist[[ii]])#is CL string present in CLlist?
              if (length(xx)==0) {      #pattern CL not present in CLlist
                 CommonCL <<- CommonCL[-jj]  #drop the elements of CommonCL not present in CommonCL
              }
           }
       }
   }

   MakeBaseLine <- function(SampName, Object) {
        BasLinType <- NULL
        BgDeg <- NULL
        Symbol <- Object@Symbol

        Selectwin <- gwindow("DEFINE BASELINE", visible=FALSE, parent=window)
        SelectGrp <- ggroup(horizontal=FALSE, container=Selectwin)
        txt <- paste("NO Baseline found on ", Symbol, " of ", SampName, sep="")
        message <- glabel(txt, container=SelectGrp)
        font(message) <- list(weight="bold")
        gseparator(horizontal=TRUE, container=SelectGrp)
        plot(Object)
		      bg <- gradio(items=c("linear","Shirley", "polynomial", "spline"),  horizontal=TRUE, handler=function(h, ...){
                              BasLinType <- svalue(bg)
                              if (BasLinType=="polynomial"){
                                  enabled(PolyDeg) <- TRUE
                              } else {
                                  enabled(PolyDeg) <- FALSE
                              }
                      }, container=SelectGrp)
        PolyDeg <- gedit(initial.msg = "Polynom degree:", container=SelectGrp)
        enabled(PolyDeg) <- FALSE
		      gbutton("OK", handler=function(...){
		                        BasLinType <<- svalue(bg)
		                        BgDeg <<- as.numeric(svalue(PolyDeg))
                              dispose(Selectwin)
                      }, container=SelectGrp)
        visible(Selectwin) <- TRUE
        Selectwin$set_modal(TRUE)
        plot(Object)

        if (BasLinType == "linear" || BasLinType == "Shirley") {
           gmessage(msg="==> Set the Baseline Limits", title="HELP INFO", icon="info")
           pos<-locator(n=2, type="p", col="red", lwd=2)
           Object@Boundaries$x<-pos$x
           Object@Boundaries$y<-pos$y
           Object <- XPSsetRegionToFit(Object)
           Object <- XPSbaseline(Object, BasLinType, BgDeg, Wgt, splinePoints )
           Object <- XPSsetRSF(Object)
        } else if (BasLinType == "polynomial") {
           gmessage(msg="==> Set the Baseline Limits", title="HELP INFO", icon="info")
           pos<-locator(n=2, type="p", col="red", lwd=2)
           Object@Boundaries$x<-pos$x
           Object@Boundaries$y<-pos$y
           Object <- XPSbaseline(Object, BasLinType, BgDeg, Wgt, splinePoints )
        } else if (BasLinType == "spline") {
            splinePoints<-list(x=NULL, y=NULL)
            txt<-"==> LEFT click to set spline points; RIGHT to exit"
            gmessage(msg=txt, title="HELP INFO", icon="info")
            pos<-c(1,1) # only to enter in  the loop
            while (length(pos) > 0) {  #pos != NULL => mouse right button not pressed
                  pos<-locator(n=1, type="p", col=3, cex=1.5, lwd=2, pch=1)
                  if (length(pos) > 0) {
                      splinePoints$x <- c(splinePoints$x, pos$x)  # $x and $y must be separate to add new coord to splinePoints
                      splinePoints$y <- c(splinePoints$y, pos$y)
                  }
            }
            # Now make BaseLine
            decr<-FALSE #Kinetic energy set
            if (Object@Flags[1] == TRUE) {
               idx<-order(splinePoints$x, decreasing=decr)  #idx is the vector of permutations for ascending order
               splinePoints$x<-splinePoints$x[idx] #splinePoints$x[idx] = SplinePoints selected in ascending order
               splinePoints$y<-splinePoints$y[idx] #splinePoints$x[idx] = SplinPoints selected in ascending order
            }
            LL<-length(splinePoints$x)

            Object@Boundaries$x<-c(splinePoints$x[1],splinePoints$x[LL]) #set the boundaries of the baseline
            Object@Boundaries$y<-c(splinePoints$y[1],splinePoints$y[LL])
            Object <- XPSsetRegionToFit(Object)
            Object <- XPSbaseline(Object, BasLinType, BgDeg, Wgt, splinePoints )
            Object <- XPSsetRSF(Object)
        }
        plot(Object)
        return(Object)
   }


   Calc.Var <- function(){
   #function that calculate the variance associated to the spectrum SpectDataY[[jj]][[ii]]
   #the standard deviation and the variance of the spctral intensity is computed following the article
   #of K. Harrison L. B. Hazel1 "The Determination of Uncertainties in Quantitative XPS/AES and its
   #Impact on Data Acquisition Strategy",SIA (1002), 18, 368-376. In modern instruments the number
   #signal intensity is high and generally the SNR low. Then the effect of noise on the StDev is low.
   #Higher is the effect of baseline selection and ascillation of X-ray intensity which cannot ne accounted.
   #StDev(I) is computed using eq. 6 Harrison. Concerning the concentrations we DO NOT use eq. 2e Harrison.
   #VARIANCE PROPERTIES:
   #SUM uncorrelated variables:  Var[SUM_i(Xi)] = SUM_i[ Var(Xi) ]
   #SUM correlated variables:  Var[SUM_i(Xi)] = SUM_i SUM_j[ Cov(Xi,Xj) ] = SUM_i[Var(Xi)] +2*SUM_i,j[Cov(Xi,Xj)] where 1 <=i < j <=n
   #    where cov(X,Y)=the covariance is defined as the expected value (or mean) of the product of their deviations
   #    from their individual expected values:  E[(X-E(X))*(Y-E(Y))]
   #    for uncorrelated variables cov(X,Y)=0
   #FRACTION: the variance of a fraction cannot computed be only approximated:
   #    var(X/Y) ~ [E(X)/E(Y)]^2 * [var(X)/E(X)^2  + var(Y)/E(Y)^2 + cov(X,Y)/E(X)E(Y)]
   #    if X and Y are uncorrelated:
   #    var(X/Y) ~ [E(X)/E(Y)]^2 * [var(X)/E(X)^2  + var(Y)/E(Y)^2]
   #
   #Let us indicate the expectation value E of the intensity of element Xj with E(Xj)
   #For the total intensity Y the expectation value is E(Y) = E(X1+X2+...+Xn) = E(X1)+E(X2)+...+E(Xn)
   #The variance associated to the concentration of X will be:
   #
   #   var(Cj) = var(Xj/Y) ~ [E(Xj)/E(Y)]^2 * [var(Xj)/E(Xj)^2  + var(Y)/E(Y)^2 + cov(Xj,Y)/E(Xj)E(Y)]
   #
   #for independent variables Xj, Y the covariance cov(Xj,Y) = 0. Xj and Y are not because Y = SUM.i(Xi)
   #includes Xj. We can solve the problem considering: Z=Y-Xj and then the concentration Cj will be:
   #
   #   Cj = Xj/(Xj+Z) but now X and Z are independent and then
   #
   #   var(Cj) =  [E(Xj)/E(Z)]^2 * [var(Xj)/E(Xj)^2  + var(Z)/E(Z)^2 ]
   #
   #More details in  Error estimates in CASAXPS on line PDF version.
   #
   #Var.I = matrix[Nelmnts, Nangles] marix of variances associated to peak intensities
   #Var.C = matrix[Nelmnts, Nangles] marix of variances associated to element concentrations

        StDev.I <- matrix(data=0, nrow=N.CL, ncol=N.XS)#Standard deviation of coreline spectral intensities
        StDev.C <- matrix(data=0, nrow=N.CL, ncol=N.XS)#Standard deviation of coreline concentrations
        Var.I <- matrix(data=0, nrow=N.CL, ncol=N.XS)  #variance associated to the spectral intensities
        Var.C <- matrix(data=0, nrow=N.CL, ncol=N.XS)  #variance associated to the element concetrations
        RSF <- array(dim=N.CL)
        Estep <- array(dim=N.CL)
        normFact <- array(dim=N.CL)
        Itot <- NULL
        Imed <- NULL
        Object <- NULL
        CL.Fit <- list()
        NTrials <- 50 #N synthesized spectra to compute the StDev on the spectrum integral intensity
        Noise <- NULL  #array containing the noise
        SyntSpect <- NULL #array to store the Synthesize Spectrum
        Isynth <- matrix(data=NA, nrow=N.CL, ncol=NTrials)       #synthesized spectral intensity to compute the StDev
        IsynthMax <- matrix(data=NA, nrow=N.CL, ncol=NTrials)    #synthesized spectral intensity to compute the StDev
        IsynthMin <- matrix(data=NA, nrow=N.CL, ncol=NTrials)    #synthesized spectral intensity to compute the StDev
        SyntCLconc <- matrix(data=0, nrow=N.CL, ncol=NTrials)    #Concentration of synthesized corelines
        SyntCLconcMax <- matrix(data=0, nrow=N.CL, ncol=NTrials) #Max Concentration of synthesized corelines
        SyntCLconcMin <- matrix(data=0, nrow=N.CL, ncol=NTrials) #Min Concentration of synthesized corelines
        SyntCLconcMed <- array(dim=N.CL) #Mean Concentration of synthesized corelines

        rho <<- as.numeric(svalue(FilmDens))        
        for(tt in 1:N.XS){ #this for runs on the N. XPSSample == N. tilt angles
#For fixed tilt angle now all the corelines are fitted to synthezize a noise-free curve
            Itot <- 0
            for(ii in 1:N.CL){     #this for runs on the N. CoreLines
                CLname <- SelectedCL[ii]  #generally occurs not all the corelines of an XPSSample are selected also the order of selection si unknown
cat("\n XPSSample name: ", SelectedFName[tt], CLname)
                FName <- get(SelectedFName[tt], envir=.GlobalEnv)
                Object <- FName[[CLname]]
                RSF[ii] <- Object@RSF
                Estep[ii] <- abs(Object@.Data[[1]][2] - Object@.Data[[1]][1])
                Object@Components <- list()   #reset existing fit
                Object@Fit <- list()
                CL.Fit[[ii]] <- list(x=NULL, y=NULL)
                Rx <- range(Object@RegionToFit$x)   #range returns always Rx[1} < Rx[2]
                NComp <- floor((Rx[2]-Rx[1])/1.5)   #each 1.5eV a Gaussian component
#Noise affects the signal. However as pointed out by Harrison, this is not the main source of error
#when performing the quantification because generally the spectra have a good SNR. Errors come from background
#substraction, imprecision of sensitivity factors and fluctuation of X intensity. Then following eq. 6 harrison:
#Observe that Estep/RSF are outside sqrt() This prevent to directly correct A for these parameters

                for(kk in 1:NComp){
                    xx <- Rx[1]+kk*1.5
                    idx <- findXIndex(Object@RegionToFit$x, xx)
                    yy <- Object@RegionToFit$y[idx]
                    Object <- XPSaddComponent(Object, type = "Gauss", peakPosition = list(x = xx, y = yy) )
                }
                Object <- XPSFitLM(Object, verbose=FALSE)
                plot(Object)
#                answ <- gconfirm("Is the peak fitting OK?", title="CONTROL THE FIT QUALITY", icon="info")
answ <- TRUE
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
                CL.Fit[[ii]]$x <- Object@RegionToFit$x
                CL.Fit[[ii]]$y <- Object@Fit$y  #the fit does notg contain the baseline contribution
                MM <- max(CL.Fit[[ii]]$y)
                MM <- findYIndex(CL.Fit[[ii]]$y, MM, 0.001)[1]
                PosMax <- CL.Fit[[ii]]$x[MM]
                if(Object@Flags) {  #Binding energy scale
                   XEnergy <- get("XPSSettings", envir=.GlobalEnv)$General[5] #X radiation energy
                   XEnergy <- as.numeric(XEnergy)
                   PosMax <- XEnergy-Posmax  #transforms Binding energy in Kinetic energy
                }
                attLength <- (1000/rho)*(49/PosMax^2 + 0.11*PosMax^0.5)
                normFact[ii] <- (Estep[ii]/RSF[ii])
            }

#We use the core-line fit to generate NTrials independent spectra affected by white noise. It is then possible to
#compute the concentration using these spectra and the concentration mean values for each element for each tilt
            for (kk in 1:NTrials){
                Itot <- 0
                for(ii in 1:N.CL){
                    LL <- length(CL.Fit[[ii]]$y)
                    CLname <- SelectedCL[ii]   #generally occurs not all the corelines of an XPSSample are selected also the order of selection si unknown
                    RSF[ii] <- FName[[CLname]]@RSF
                    Estep[ii] <- abs(FName[[CLname]]@.Data[[1]][2] - FName[[CLname]]@.Data[[1]][1])
                    MaxNoise <- max(CL.Fit[[ii]]$y - (FName[[CLname]]@RegionToFit$y - FName[[CLname]]@Baseline$y)) #max of the noise superposed to the coreline
                    MinNoise <- min(CL.Fit[[ii]]$y - (FName[[CLname]]@RegionToFit$y - FName[[CLname]]@Baseline$y)) #max of the noise superposed to the coreline
                    NoiseAmpli <- MaxNoise-MinNoise
                    Noise <- runif(LL, min=-NoiseAmpli/2, max=NoiseAmpli/2) #white noise generation
                    Isynth[ii,kk] <- sum(CL.Fit[[ii]]$y + Noise) #synthesized spectrum = noise-free spectrum + noise
                    Noise <- runif(LL, min=0, max=NoiseAmpli) #white noise generation
                    IsynthMax[ii,kk] <- sum(CL.Fit[[ii]]$y + Noise) #synthesized spectrum = noise-free spectrum + all positive noise
                    IsynthMin[ii,kk] <- sum(CL.Fit[[ii]]$y - Noise) #synthesized spectrum = noise-free spectrum - all positive noise
                }
            }

            for (kk in 1:NTrials){
                 for(ii in 1:N.CL){
                     Itot <- sum(Isynth[ ,kk]*normFact[])  #each Isynth is normalized for its Estep/RSF and then the integral is computed
                     SyntCLconc[ii,kk] <- (Isynth[ii,kk]*normFact[ii])/Itot
                     ItotMin <- sum(IsynthMin[ ,kk]*normFact[])+(IsynthMax[ii,kk]-IsynthMin[ii,kk])*normFact[ii]
                     SyntCLconcMax[ii,kk] <- (IsynthMax[ii,kk]*normFact[ii])/ItotMin
                     ItotMax <- sum(IsynthMax[ ,kk]*normFact[])+(IsynthMin[ii,kk]-IsynthMax[ii,kk])*normFact[ii]
                     SyntCLconcMin[ii,kk] <- (IsynthMin[ii,kk]*normFact[ii])/ItotMax
                 }
            }

            for(ii in 1:N.CL){
                SyntCLconcMed[ii] <- sum(SyntCLconcMax[ii, ] + SyntCLconcMin[ii, ])/(2*NTrials)
                Var.C[ii,tt] <- (sum(SyntCLconcMax[ii,]-SyntCLconcMed[ii])^2 + sum(SyntCLconcMin[ii,]-SyntCLconcMed[ii])^2)/NTrials
                StDev.C[ii,tt] <- sqrt(Var.C[ii,tt])
                Imed <- sum(CL.Fit[[ii]]$y)
                StDev.I[ii,tt] <- normFact[ii]*sqrt(sum((Imed-Isynth[ii, ])^2)/NTrials) #CasaXPS: "Error estimates in CasaXPS" pp. 16
                Var.I[ii,tt] <- StDev.I[ii,tt]^2
            }

#print(StDev.C)
#scan(n=1)
        }
cat("\n \n")
cat("\n ************  Intensity StDev ************\n")
print(StDev.I)
cat("\n *********** Intensity Variance ***********\n")
print(Var.I)
cat("\n ************  Concentration StDev ************\n")
print(StDev.C)
cat("\n *********** Concentration Variance ***********\n")
print(Var.C)
cat("\n **********************************************\n")
        return(Var.C)
   }


   FitData <- function(){

            import::here(modFit, .from=FME)
            import::here(gradient, .from=rootSolve)
            model <- as.numeric(svalue(DP.Model, index=TRUE))

# FitResiduals() at beginning are evaluated using the Start parameters. Then modFit() generates new parameters
# and calls the FitResiduals() with this new Parms. FitResiduals(Parms) calls the FitFunction() which now is
# evaluated using the new parameters FitParms. Then residuals estimated using the new Parms are evaluated.
# Following the Start format, also the new parameters are generated using the same names of variables.

            FitFunction <- function(Start, FitExpr){
                       if (model == 1) {
                          d <- Start["d"]
                          bkg <- Start["bkg"]
                       } else if (model == 2) {
                          d <- Start["d"]
                          A <- Start["A"]
                          Rgh <- Start["Rgh"]
                          blr <- Start["blr"]
                       }
                       FitCurve <- eval(parse(text=FitExpr))
                       return(FitCurve)
            }

            FitResiduals <- function(Parms) {
                       residuals <- DataToFit$y-FitFunction(Parms, FitExpr)
                       return(residuals)
            }

            FitMtd <- svalue(FitMethod, index=TRUE)
            ptol <- maxiter <- nprint <- NULL
            if (FitMtd ==  1){
               FitMtd <-"Marq"
               ctrl<-list(ptol= 1e-6, maxiter=1000, nprint=1)   #ctrl for Marquartd optimization
            } else if (FitMtd == 2) {
               FitMtd <-"Newton"
               ctrl<-list(ptol= 1e-6, maxiter=1000, nprint=1)   #ctrl for Newton
            } else if (FitMtd == 3){
               FitMtd <-"Port"
               ctrl<-list(rel.tol=1e-6, eval.max=200, iter.max=1000, trace=1)  # ctrl for Port
            } else if (FitMtd == 4) {
               FitMtd <-"CG"
               ctrl<-list(reltol=1e-6, iter.max=1000, trace=1) # ctrl for Nelder-Mead CG BFGS L-BFGS-B SANN
            } else if (FitMtd == 5) {
               FitMtd <-"SANN"
               ctrl<-list(reltol=1e-6, iter.max=1000, trace=1) # ctrl for Nelder-Mead CG BFGS L-BFGS-B SANN
            } else if (FitMtd == 6) {
               FitMtd <-"Pseudo"
               ctrl<-list(varleft=1e-6, numiter=1000, verbose=TRUE) # ctrl per Pseudo
            }
            FitExpr <- NULL

#--- Due to limited experimental data ModFit is applied for fitting
#    To avoid defining an XPSSample and a Coreline with ARXPS data stored,
#    Modfit is integrated in this GUI and called directly

            if (model == 1) { ### Thickness estimated using Classic decay model: I1/I2 = R.Rsf * {exp[-d/L1 sin(t) ]} / {1 - exp[d/L2 sin(t) ]}
               FitCurve <<- NULL
               bkg <- NULL
               Start <- c(d=0.1, bkg=0)  # c collects the set of parameters with their names
               Lower <- c(0, 0)
               Upper <- c(10, DataToFit$y[1])
               IniParms <<- data.frame(
                           row.names = c("d", "bkg"),
                           start = Start,
                           min = Lower,
                           max = Upper)
#               DataToFit <- list(x=c(90,74.14,68.6,60.2,49,40.1,31.5,25.2,20,18.7,17.4,17,16.8)*pi/180, #example of ARXPSData
#                                 y=c(2.65,2.62,2.62,2.6,2.6,2.55,2.5,2.4,2.28,2.2,2,1.58,1.1))
               Tilt <- DataToFit$x

               FitExpr <- "bkg + R.Rsf * exp(-d/(La*sin(Tilt))) / (1-exp(-d/(Lb*sin(Tilt)) ) )"
               FitEstimation <<- modFit(f = FitResiduals, p = Start, lower=Lower, upper=Upper, method=FitMtd, control=ctrl)
               FitParms <<- FitEstimation$par
               d <- unlist(FitParms["d"])
               bkg <- unlist(FitParms["bkg"])

               cat("\n ----Best Fit Param----\n")
               print(FitParms)
               cat("\n ----------------------\n")

               FitCurve <<- bkg + R.Rsf * exp(-d/(La*sin(Tilt))) / (1-exp(-d/(Lb*sin(Tilt)) ))
               Tilt <- Tilt*180/pi
               x <- cbind(Tilt, Tilt)
               y <- cbind(DataToFit$y, FitCurve)
               matplot(x, y, type="b", pch=16, cex=2, col=c("blue", "red"), xlab="Tilt Angle (deg)", ylab="Fit")


            } else if (model == 2) { ### Rougness modified Classic Model
               FitCurve <<- NULL
               Start <- c(d=1, A=0.01, Rgh=0.0001, blr=0.1)  # c collects the set of parameters with their names
               Lower <- c(0, 1e-7, 0, 1e-2)
               Upper <- c(10, 1,   5, 0.1 )
               IniParms <<- data.frame(
                           row.names = c("d", "A", "Rgh", "blr"),
                           start = Start,
                           min = Lower,
                           max = Upper)
#               DataToFit <- list(x=c(90,74.14,68.6,60.2,49,40.1,31.5,25.2,20,18.7,17.4,17,16.8)*pi/180,   #example of ARXPSData
#                                 y=c(2.65,2.62,2.62,2.6,2.6,2.55,2.5,2.4,2.28,2.2,2,1.58,1.1))
               Tilt <- DataToFit$x
#                           +-----roughness factor-----+  +--------------classic function with blurred dependance on tilt--------------+
               FitExpr <- "(1-A*exp(1/(sin(Tilt)-Rgh) )) * R.Rsf * exp(-d/(La*(1-blr+blr*sin(Tilt))) ) / (1-exp(-d/(Lb*(1-blr+blr*sin(Tilt)) )))"
#                                                                           +-------Blurr------+                      +-------Blurr------+
               FitEstimation <<- modFit(f = FitResiduals, p = Start, lower=Lower, upper=Upper, method=FitMtd, control=ctrl)
               FitParms <<- FitEstimation$par

               cat("\n ----Best Fit Param----\n")
               print(FitParms)
               cat("\n ----------------------\n")

               d <- FitParms["d"]
               A <- FitParms["A"]
               Rgh <- FitParms["Rgh"]
               blr <- FitParms["blr"]
               FitCurve <<- (1-A*exp(1/(sin(Tilt)-Rgh) )) * R.Rsf * exp(-d/(La*(1-blr+blr*sin(Tilt))) ) / (1-exp(-d/(Lb*(1-blr+blr*sin(Tilt)) )))
               Tilt <- Tilt*180/pi
               x <- cbind(Tilt, Tilt)
               y <- cbind(DataToFit$y, FitCurve)
               matplot(x, y, type="b", pch=16, cex=2, col=c("blue", "red"), xlab="Tilt Angle (deg)", ylab="Fit")

            } else if (model == 3) { #Max Entropy Model
               ####### Max Entropy method
            }

            return(FitCurve)
   }

   SaveResults <- function(h,...){
         SaveWin<-gwindow("SAVE DEPTH PROFILE RESULTS", visible=FALSE)
         SaveGroup <- ggroup(label="", horizontal=FALSE, container=SaveWin)

         SaveFrame <-gframe(text=" OUTPUT DATAFILE ", spacing=5, container=SaveGroup)
         SaveObj1 <- glabel(SelectedFName, container=SaveFrame)
         SaveObj2 <- gedit("", initial.msg="Output File Name: ", container=SaveFrame)
         SaveObj3 <- gbutton(" OK " , handler=function(h, ...){
                           DestFile <- svalue(SaveObj2)
                           DestFile <- unlist(strsplit(DestFile, "\\."))
                           DestFile <- paste(DestFile[1], ".RData", sep="")  #Force the extension to ".RData"
                           dispose(SaveWin)

                           FNameOut <- new("XPSSample")
                           FNameList <- ""
                           FileNames <- NULL
                           CLnames <- NULL
                           kk <- 1
                           for (ii in 1:N.XS){  #Exstract the Selected Coreline from the list of selected XPSSamples
                               for (jj in 1:N.CL){
                                   FName <- get(SelectedFName[ii], envir=.GlobalEnv)
                                   FNameOut[[kk]] <- FName[[SelectedCL[jj]]] #load the selected coreline in a temporary XPSSample
                                   CLnames <- c(CLnames, SelectedCL[jj])
                                   kk <- kk+1
                               }
                               FileNames <- c(FileNames, FName@Filename) #for each XPSSample obtain the data FileName
                               FNameList <- paste(FNameList, SelectedFName[ii], sep=", ")
                           }
                           FNameOut[[kk]] <- new("XPSCoreLine") #create a new coreline structure to save ARXPS data (data to fit, best fit, thickness...

                           CLnames <- c(CLnames, "I1/I2")
                           if (N.CL == 1){      #initialize CLcomment when only 1 coreline selected: DepthProfile performed on CL fit components
                              CLcomment <- paste("coreline ", SelectedCL[1], "fit components: ", svalue(ComponentCK[[1]]), svalue(ComponentCK[[2]]), sep="")
                           } else if (N.CL == 2) {
                              CLcomment <- paste("corelines: ", SelectedCL[1], ", ", SelectedCL[2], sep="")  #initialize CLname when only 2 corelines selected
                           }
                           DPModel <- svalue(DP.Model)
                           FNameOut@Sample <- FName@Sample
                           FNameOut@Comments <- c(paste("==>ARXPS: depth profile analysis on files: ", FNameList, sep=""),
                                                     paste("==>Analysis performed on ", CLcomment, sep=""),
                                                     paste("==>",DPModel, " model used to fit the ratio of their spectral intensities", sep=""))
                           FNameOut@User <- FName@User
                           FNameOut@Filename <- FileNames
                           FNameOut@names <- CLnames
                           #Define Data Structure for Depth Profile fitting
                           DPModel <- "Classic Depth Profile Model"  #by default Classic method set as Model name
                           if (DPModel[1] == "Roughness") { newDPModel <- "Classic Depth Profile Model" }   #Roughness modified method set name DpthProfileRoughness
                           if (DPModel[1] == "Max") { DPModel <- "Maximum Entropy Depth Profile Model" } #Max Entropy method  set name DpthProfileEntropy
                           Tilt <- Tilt*180/pi
                           FNameOut[[kk]]@.Data[[1]] <- Tilt      #abscissa are the tilt angles
                           FNameOut[[kk]]@.Data[[2]] <- Ratio     #ordinate are the ratios
                           FNameOut[[kk]]@.Data[[3]] <- rep(1, length(Tilt))  #analyzer transfer function here unitary for all spectral data
                           FNameOut[[kk]]@RegionToFit$x <- Tilt      #abscissa are the tilt angles
                           FNameOut[[kk]]@RegionToFit$y <- Ratio     #ordinate are the ratios
                           FNameOut[[kk]]@Baseline$x <- Tilt
                           FNameOut[[kk]]@Baseline$y <- Bkg <- rep(0, N.XS) 
                           FNameOut[[kk]]@Baseline$baseline <- new("baseline",
                                          baseline = matrix(data=Tilt, nrow=1),   #matrix() required by class("baseline") see baseline() function of R
                                          corrected = matrix(data=Ratio-min(Ratio),nrow=1),
                                          spectra = matrix(data=Ratio, nrow=1),
                                          call = match.call()
                                        )
                           FNameOut[[kk]]@Baseline$type = "Linear"
                           FNameOut[[kk]]@Fit$y <- FitCurve
                           FNameOut[[kk]]@Fit$fit <- FitEstimation
                           FNameOut[[kk]]@Boundaries <- list(x=range(Tilt*180/pi), y= range(Ratio))
                           FNameOut[[kk]]@units <- c("Tilt Angle (deg)", paste(SelectedSpect[2], SelectedSpect[1], sep="/")) #adopt the units and flags of the last loaded FName[[ coreline1 ]]
                           FNameOut[[kk]]@Flags <- c(FALSE, FALSE, FALSE, FALSE) #FName[[SelectedCL[1]]]@Flags
                           FNameOut[[kk]]@Info <- paste("ARXPS analysis using ", DPModel, " on ", paste(FileNames, collapse=" "), " XPSSamples", sep="")
                           FNameOut[[kk]]@Symbol <- "I1/I2"

                           #Add Fit Function to experim. data
                           FNameOut[[kk]] <- XPSaddComponent(FNameOut[[kk]], type = "Generic")
                           names(FNameOut[[kk]]@Components) <- "DP"
                           FNameOut[[kk]]@Components[[1]]@funcName <- "Depth_Profile"
                           FNameOut[[kk]]@Components[[1]]@description <- DPModel
                           FNameOut[[kk]]@Components[[1]]@label = "DP"

                           #Setting Fit Parameters
                           FNameOut[[kk]]@Components[[1]]@param <- IniParms
                           FNameOut[[kk]]@Components[[1]]@param$Fit <- FitParms
                           FNameOut[[kk]]@Components[[1]]@ycoor <- FitCurve
                           assign(DestFile, FNameOut, envir=.GlobalEnv)
#print(str(FNameOut[[kk]]))
                           plot(FNameOut)
            }, container=SaveFrame)
            visible(SaveWin) <- TRUE
      XPSSaveRetrieveBkp("save")
   }


#-------variables---
#---load list of file ID and correspondent FileNames
      FNameList<-XPSFNameList()
      if (length(FNameList) == 0) { return() }
      SampID <- ""
      SpectList <- ""
      SourceFileList <- NULL
      SelectedFName <- NULL
      CLlist <- list()
      ComponentCK <- list() #list used to make widgets for selection of coreLine fit components
      CommonCL <- NULL #the list of Corelines of the first XPSSample used as reference
      SelectedCL <- NULL
      SelectedSpect <- NULL
      N.XS <- NULL     #number of selected XPSSamples
      N.CL <- NULL     #number of selected common lines


      Tilt <- NULL
      angles <- list()
      infoLab <- list()
      rho <- NULL
      XrayE <- NULL
      N.Layers <- NULL     #N. Layers used to model the depth profile - MaxEnt
      L.Thickness <- NULL  #Thickness of the individual layer - MaxEnt
      AreaComp1 <- NULL
      AreaComp2 <- NULL
      Ratio <- NULL
      R.Rsf <- NULL
      PosMax1 <- NULL
      PosMax2 <- NULL
      BE <- NULL
      La <- NULL    #attenuation length element1    cannot use the name L1 problems in XPSFitLM
      Lb <- NULL    #attenuation length element2    cannot use the name L2 problems in XPSFitLM
      FitEstimation <- NULL
      FitCurve <- NULL
      IniParms <- data.frame()
      FitParms <- NULL
      IniData <- FALSE

      SpectDataX <- list()
      SpectDataY <- list()
      SpectDataRSF <- list()
      Exp.Data <- NULL #matrix[Nelements, NTilt] = integral intensities of all elements for all the tilt angles - MaxEnt
      Att.L <- NULL    #vector[N.Elements] = Element attenuation lengths  - MaxEnt
      Variance <- NULL #matrix[Nelements, NTilt] = Variance of element concentration due to the spectral noise - MaxEnt
      Intensity <- NULL #matrix[Nelements, NTilt]= integral intensity of corelines for the various tilt angles
      Conc <- NULL #matrix[Nelements, NTilt]= element concentrations foe each tilt angle



#####----main---
      AddWin <- gwindow("ARXPS: DEPTH PROFILE ANALYSIS", visible=FALSE)
      size(AddWin) <- c(600, 400)
      Addgroup <- ggroup(horizontal=FALSE, container=AddWin)
      NoteBK<-gnotebook(expand=TRUE, container = Addgroup)

      SelGroup <- ggroup(label="CORELINE SELECTION", horizontal=TRUE, spacing=5, container=NoteBK)

      layoutS1 <- glayout(homogeneous=FALSE, spacing=5, container=SelGroup)

      layoutS1[1,1] <-  AddFrame1 <- gframe("SELECT THE SOURCE XPS-SAMPLES", spacing=5, horizontal=TRUE, container=layoutS1)
      SourceFiles <- gcheckboxgroup(FNameList, selected=-1, horizontal=FALSE, handler=function(h,...){
                                   CLlist <<- NULL  #reset the CoreLine list to be checked
                                   txt <- NULL
                                   SelectedFName <<- svalue(SourceFiles)
                                   N.XS <<- length(SelectedFName)
                                   for (ii in 1:N.XS){
                                       CLlist[[ii]] <<- XPSSpectList(SelectedFName[ii])
                                       names(CLlist)[[ii]] <<- as.character(ii)
                                       txt <- paste(txt, paste(CLlist[[ii]], collapse="  "), "\n")  #transform vector of strings CLlist[] in one string and add carriage return
                                       svalue(CoreLineList) <- txt
                                   }
                         }, container = AddFrame1)

      layoutS1[1,2] <-  AddFrame2 <- gframe("CORE LINES LIST", spacing=5, horizontal=FALSE, container=layoutS1)
      CoreLineList <- glabel(text="", container=AddFrame2)  #just to not collapse AddFrame2

      layoutS1[2,1] <- Compare <-gbutton("COMPARE CORE LINE NAMES", handler=function(h, ...){
                                   CommonCL <<- XPSSpectList(SelectedFName[1])
                                   CheckCL()   #check for corelines common to selected XPSSamples
                                   delete(AddFrame3,CommonCoreLines)
                                   CommonCoreLines<<-gcheckboxgroup(CommonCL, selected=-1, horizontal=TRUE, handler=function(h,...){
                                                                          enabled(Compare) <- FALSE
                                                                          enabled(SelectCL) <- TRUE
                                                                   }, container=AddFrame3)
                         }, container=layoutS1)


      layoutS1[3,1] <-  AddFrame3 <- gframe("COMMON CORE LINES", spacing=5,  container=layoutS1)
      CommonCoreLines <- glabel("  ", container=AddFrame3)  #just to not collapse AddFrame3

      layoutS1[3,2] <- SelectCL <-gbutton("SELECT CORELINES", handler=function(h,...){
                               SlctdCL <- svalue(CommonCoreLines)
                               N.CL <<- length(SlctdCL)
                               if (N.CL > 2) {
                                  gmessage(msg=" More than two core lines selected! \n Max Entropy method will be selected.", title="CORE LINE SELECTION", icon="warning")
                               } 
                               
# ===> N.CL = 1 analysis of fit components
                               if (N.CL == 1) {
                                  CLname <- unlist(strsplit(SlctdCL, "\\."))   #drop "NUMBER." at beginning of coreLine name
                                  CLname <- CLname[2]
                                  SelectedCL <<- CLname                  #save selected core lines
                                  txt <- paste("Only ", CLname, " selected! Do you want to analyze fitting components?", sep="")
                                  answ <- gconfirm(msg=txt, title="ANALYSIS ON FIT COMPONENTS", icon="warning")
                                  if (answ == "FALSE") {
                                      gmessage(msg="Please select the second core line and proceed", title="CORE LINE SELECTION", icon="warning")
                                      return()
                                  }
                                  enabled(SelectCL) <- FALSE
                                  tmp <- new("XPSSample",
                                             Project = " ",
                                             Comments = " ",
                                             User=Sys.getenv('USER'),
                                             Filename=" " )

                                  for(ii in 1:N.XS){ #Now load the selected coreline from XPSSamples at different tilt angle
                                     FName <- get(SelectedFName[ii], envir=.GlobalEnv)
                                     txt <- paste(SelectedFName[ii], " Fit Comp.: ", sep="")
                                     layoutFC[ii,1] <- glabel(txt, spacing=5, container=layoutFC)
                                     tmp[[ii]] <- FName[[CLname]] #load the selected coreline in a temporary XPSSample
                                     CompNames <- names(tmp[[ii]]@Components)  #List of names of the fitting components
                                     layoutFC[ii,2] <- ComponentCK[[ii]] <<- gcheckboxgroup(CompNames, horizontal = TRUE, checked=FALSE, container = layoutFC) #CoreLineComp is an array
                                  }
                                  E.Units <- "Kinetic Energy [eV]"
                                  BE <<- FALSE #Kinetic energy scale for the XPSSamples
                                  if (FName[[CLname]]@Flags[1] == TRUE) {
                                      BE <<- TRUE        #BE found for the last XPSSample holds also for the other loaded XPSSamples
                                      E.Units <- "Binding [eV]"
                                  }
                                  plot(tmp)
                                  enabled(SelectFC) <- TRUE
                                  for (ii in 1:N.XS){
                                      Tiltlayout[ii, 1] <<- glabel(SelectedFName[ii], container=Tiltlayout)
                                      Tiltlayout[ii, 2] <<- angles[[ii]] <<- gedit("", initial.msg="Tilt ang.", container=Tiltlayout)
#                                     Tiltlayout[ii, 3] <<- infoLab[[ii]] <<- glabel(" " , container=Tiltlayout)
                                  }

# ===> N.CL = 2 compare intensity of 2 corelines
                               } else if (N.CL == 2) {
                                  RSF <- NULL
                                  #Store spectral data and info in the SpectData list
                                  for (ii in 1:N.CL){   #N.CL must be == 2
                                      CLname <- SlctdCL[ii]
                                      CLname <- unlist(strsplit(CLname, "\\."))#skip the CoreLine index
                                      CLname <- CLname[2]
                                      SelectedCL[ii] <<- CLname                  #save selected core lines name for labeling axis and messages
                                      RSF[ii] <- paste("RSF_", CLname, sep="")
                                      BE <<- NULL
                                      for (jj in 1:N.XS){                          #N.XS = number of selected XPS Samples to analyze
                                          FName<-get(SelectedFName[jj], envir=.GlobalEnv)
                                          SpectDataX[[CLname]][[jj]] <<- FName[[CLname]]@RegionToFit$x  #abscissa of selected spectrum
                                          if ( ! hasBaseline(FName[[CLname]]) ){ #Has a Baseline the selected coreline in the list of XPSSample?
                                             FName[[CLname]] <- MakeBaseLine(SelectedFName[jj], FName[[CLname]])         #if not, define the Baseline
                                             assign(SelectedFName[jj], FName, envir=.GlobalEnv) #save defined baseline in the .GlobalEnv
                                          }
                                          SpectDataY[[CLname]][[jj]] <<- FName[[CLname]]@RegionToFit$y-FName[[CLname]]@Baseline$y  #spectrum without baseline
                                          SpectDataRSF[[CLname]][1] <<- FName[[SelectedCL[ii]]]@RSF  #Relative Sensitivity Factor of the selected Coreline
                                          if (is.null(BE) ){
                                             BE <<- FName[[CLname]]@Flags[1]  #get energy units from the forst XPSSample: BE=TRUE => binding energy scale used during acquisition
                                          }  else if ( FName[[CLname]]@Flags[1] != BE) { #compare energy units and exit if different
                                             gmessage("ATTENTION: different energy scales used for XPS corelines!", title="ENERGY UNITS ERROR", icon="error")
                                             return()
                                          }
                                      }
                                  }
                                  # Structure of SpectData is
                                  # SpectDataX[[CL1]][[jj]] = list of abscissa Corelines type 1 from XPSSamples at differet tilt angles jj
                                  # SpectDataX[[CL2]][[jj]] = list of abscissa Corelines type 2 from XPSSamples at differet tilt angles jj
                                  # SpectDataY[[CL1]][[jj]] = list of spectra Corelines type 1 from XPSSamples at differet tilt angles jj
                                  # SpectDataY[[CL2]][[jj]] = list of spectra Corelines type 2 from XPSSamples at differet tilt angles jj
                                  # SpectDataRSF[[CL1]][1] = list of RSF corresponding to Corelines type 1
                                  # SpectDataRSF[[CL2]][1] = list of RSF corresponding to Corelines type 2

                                  #Now is possible to construct the widget to associate the tilt angles
                                  for (ii in 1:N.XS){
                                      Tiltlayout[ii, 1] <<- glabel(SelectedFName[ii], container=Tiltlayout)
                                      Tiltlayout[ii, 2] <<- angles[[ii]] <<- gedit("", initial.msg="Tilt ang.", container=Tiltlayout)
#                                      Tiltlayout[ii, 3] <<- infoLab[[ii]] <<- glabel(" " , container=Tiltlayout)
                                  }
                                  enabled(SelectCL) <- FALSE
                                  enabled(FilmDens) <- TRUE
                                  enabled(Xenergy) <- TRUE
                                  enabled(FitMethod) <- TRUE
                                  enabled(StartCalc) <- TRUE
                                  enabled(Layers) <- FALSE
                                  enabled(LThckns) <- FALSE
                                  svalue(NoteBK) <- 2
 
# ===> N.CL > 2 apply Max-Entropy method
                               } else if (N.CL > 2) { #MAXENT: more than two coreline selected => elemental depth profile with MaxEnt
                                  RSF <- NULL
                                  #Store spectral data and info in the SpectData list
                                  for (ii in 1:N.CL){
                                      CLname <- SlctdCL[ii]
                                      CLname <- unlist(strsplit(CLname, "\\."))#skip the CoreLine index
                                      CLname <- CLname[2]
                                      SelectedCL[ii] <<- CLname                  #save selected core lines name for labeling axis and messages
                                      RSF[ii] <- paste("RSF_", CLname, sep="")
                                      BE <<- NULL
                                      for (jj in 1:N.XS){                          #N.XS = number of selected XPS Samples to analyze
                                          FName<-get(SelectedFName[jj], envir=.GlobalEnv)
                                          SpectDataX[[CLname]][[jj]] <<- FName[[CLname]]@RegionToFit$x  #abscissa of selected spectrum
                                          if ( ! hasBaseline(FName[[CLname]]) ){ #Has a Baseline the selected coreline in the list of XPSSample?
                                             FName[[CLname]] <- MakeBaseLine(SelectedFName[jj], FName[[CLname]])         #if not, define the Baseline
                                             assign(SelectedFName[jj], FName, envir=.GlobalEnv) #save defined baseline in the .GlobalEnv
                                          }
                                          SpectDataY[[CLname]][[jj]] <<- FName[[CLname]]@RegionToFit$y-FName[[CLname]]@Baseline$y  #spectrum without baseline
                                          SpectDataRSF[[CLname]][1] <<- FName[[SelectedCL[ii]]]@RSF  #Relative Sensitivity Factor of the selected Coreline
                                          if (is.null(BE) ){
                                             BE <<- FName[[CLname]]@Flags[1]  #get energy units from the forst XPSSample
                                          }  else if ( FName[[CLname]]@Flags[1] != BE) { #compare energy units and exit if different
                                             gmessage("ATTENTION: different energy scales used for XPS corelines!", title="ENERGY UNITS ERROR", icon="error")
                                             return()
                                          }
                                      }
                                  }
                                  # SpectDataY[[1]][[jj]] = Coreline1 from XPSSamples at differet tilt angles jj
                                  # SpectDataY[[2]][[jj]] = Coreline2 from XPSSamples at differet tilt angles jj
                                  # ...
                                  # SpectDataY[[N.CL]][[jj]] = Coreline N from XPSSamples at differet tilt angles jj

                                  # SpectDataRSF[[1]][[1]] = RSF1 corresponding to Coreline1 from XPSSamples at differet tilt angles jj
                                  # SpectDataRSF[[2]][[1]] = RSF2 corresponding to Coreline2 from XPSSamples at differet tilt angles jj
                                  # ...
                                  # SpectDataRSF[[N.CL]][[1]] = RSF.N corresponding to CorelineN from XPSSamples at differet tilt angles jj

                                  #Now is possible to construct the widget to associate the tilt angles
                                  for (ii in 1:N.XS){
                                      Tiltlayout[ii, 1] <<- glabel(SelectedFName[ii], container=Tiltlayout)
                                      Tiltlayout[ii, 2] <<- angles[[ii]] <<- gedit("", initial.msg="Tilt ang.", container=Tiltlayout)
#                                      Tiltlayout[ii, 3] <<- infoLab[[ii]] <<- glabel(" " , container=Tiltlayout)
                                  }
                                  enabled(SelectCL) <- FALSE
                                  enabled(FilmDens) <- TRUE
                                  enabled(Xenergy) <- TRUE
                                  enabled(FitMethod) <- FALSE
                                  enabled(StartCalc) <- TRUE
                                  enabled(Layers) <- TRUE
                                  enabled(LThckns) <- TRUE
                                  svalue(NoteBK) <- 2
                                  svalue(DP.Model, index=TRUE) <- 3

                               }
                         }, container=layoutS1)
                         enabled(SelectCL) <- FALSE


      layoutS1[4,1] <-  AddFrame4 <- gframe("FIT COMPONENTS", horizontal = FALSE, spacing=5,  container=layoutS1)
      layoutFC <- glayout(homogeneous=FALSE, spacing=5, container=AddFrame4)
      glabel("   ", container= AddFrame4)  #white spaces just to define the frame
      layoutS1[4,2] <- SelectFC <-gbutton("SELECT FIT COMPONENTS", handler=function(h,...){
                               #Store spectral data and info in the SpectData list
                               for(ii in 1:N.XS){ #Now load the selected coreline from XPSSamples at different tilt angle
                                  CmpFit <- svalue(ComponentCK[[ii]], index=TRUE)
                                  if (length(CmpFit) != 2) {
                                      gmessage(msg = "Two component for each XPSSample has to be selected. Please control!", title="FIT COMPONENT SELECTION", icon = "warning")
                                      return()
                                  }
                                  FName <- get(SelectedFName[ii], envir=.GlobalEnv)
                                  SpectDataX[["C1"]][[ii]] <<- FName[[SelectedCL]]@RegionToFit$x  #abscissa of selected spectrum
                                  SpectDataY[["C1"]][[ii]] <<- FName[[SelectedCL]]@Components[[CmpFit[1]]]@ycoor  #spectrum without baseline
                                  SpectDataX[["C2"]][[ii]] <<- FName[[SelectedCL]]@RegionToFit$x  #abscissa of selected spectrum
                                  SpectDataY[["C2"]][[ii]] <<- FName[[SelectedCL]]@Components[[CmpFit[2]]]@ycoor  #spectrum without baseline
                                  MM <- max(SpectDataY[["C1"]][[ii]])
                                  MM <- findYIndex(SpectDataY[["C1"]][[ii]], MM, 0.001) #index corresponding to the max
                                  MM <- MM[1]
                                  SpectDataX[["PosMaxC1"]][[ii]] <<- SpectDataX[["C1"]][[ii]][MM] #position of the component1 max
                                  MM <- max(SpectDataY[["C2"]][[ii]])
                                  MM <- findYIndex(SpectDataY[["C2"]][[ii]], MM, 0.001) #index corresponding to the max
                                  MM <- MM[1]
                                  SpectDataX[["PosMaxC2"]][[ii]] <<- SpectDataX[["C1"]][[ii]][MM] #position of the component1 max
                               }
                               SpectDataRSF[["C1"]][1]<<- FName[[SelectedCL]]@RSF  #Relative Sensitivity Factor of the selected Coreline
                               SpectDataRSF[["C2"]][1]<<- FName[[SelectedCL]]@RSF  #Relative Sensitivity Factor of the selected Coreline
                               if (SpectDataX[["PosMaxC1"]][[1]] < SpectDataX[["PosMaxC2"]][[1]]) {
                                  SelectedCL[3] <<- "Low Energy Fit Comp."      #SelectedCL[3], SelectedCL[4] labels for axis and messages
                                  SelectedCL[4] <<- "High Energy Fit Comp."     #SelectedCL[1], SelectedCL[2] reserved for CoreLine Names
                               } else {
                                  SelectedCL[3] <<- "High Energy Fit Comp."
                                  SelectedCL[4] <<- "Low Energy Fit Comp."
                               }

                               for (ii in 1:N.XS){
                                   Tiltlayout[ii, 1] <<- glabel(SelectedFName[ii], container=Tiltlayout)
                                   Tiltlayout[ii, 2] <<- angles[[ii]] <<- gedit("", initial.msg="Tilt ang.", container=Tiltlayout)
                               }
                               enabled(SelectFC) <- FALSE
                               enabled(FilmDens) <- TRUE
                               enabled(Xenergy) <- TRUE
                               enabled(FitMethod) <- TRUE
                               enabled(StartCalc) <- TRUE
                               enabled(Layers) <- FALSE
                               enabled(LThckns) <- FALSE
                               svalue(NoteBK) <- 2
                         }, container=layoutS1)
      enabled(SelectFC) <- FALSE

#---Thickness Estimation

      EstimGroup <- ggroup(label="THICKNESS ESTIMATION", horizontal=TRUE, spacing=5, container=NoteBK)

      TiltGroup <- ggroup(horizontal=FALSE, spacing=5, container=EstimGroup)
      message <- glabel("Set tilt angles of XPSSamples (degrees)", container=TiltGroup)
#      font(message) <- list(weight="bold")
      gseparator(horizontal=TRUE, container=TiltGroup)
      Tiltlayout <- glayout(homogeneous=FALSE, spacing=3, container=TiltGroup)
      #Tiltlayout defined at row 500 and at row 589 following N.Coreline selection
      gseparator(horizontal=FALSE, container=TiltGroup)
      Infolayout <- glayout(homogeneous=FALSE, spacing=3, container=TiltGroup)
      Infolayout[1, 1] <- message1 <- glabel("Film density [Kg/m^3]", spacing=8, container=Infolayout)
      Infolayout[1, 2] <- message2 <- glabel("X-ray Excitation Energy (eV)", spacing=8, container=Infolayout)
      Infolayout[2, 1] <- FilmDens <- gedit("", initial.msg="Coating Density", spacing=8, container=Infolayout)
      Infolayout[2, 2] <- Xenergy <- gedit("1486.6", initial.msg="", spacing=8, container=Infolayout)
      Infolayout[3, 1] <- message3 <- glabel("Model Profile on N. Layers", spacing=8, container=Infolayout)
      Infolayout[3, 2] <- message4 <- glabel("Layer Thickness [nm]", spacing=8, container=Infolayout)
      Infolayout[4, 1] <- Layers <- gedit("", initial.msg="50", spacing=8, container=Infolayout)
      Infolayout[4, 2] <- LThckns <- gedit("", initial.msg="0.2", spacing=8, container=Infolayout)
      enabled(FilmDens) <- FALSE
      enabled(Xenergy) <- FALSE
      enabled(Layers) <- FALSE
      enabled(LThckns) <- FALSE

      gseparator(horizontal=FALSE, container=TiltGroup)
      glabel("Thickness estimation methods", spacig=3, container=TiltGroup)
      DP.Model <- gradio(c("Classic", "Roughness Modified", "Max Entropy"), spacig=3, selected=1, horizontal=TRUE, handler = function(h,...){
                               enabled(StartCalc)<-TRUE
                         }, container=TiltGroup)
      glabel("Residual minimization methods", spacig=3, container=TiltGroup)
      FitMethod <- gradio(c("Marquardt", "Newton", "Port", "Conj.Gradient", "SANN", "Pseudo"), selected=1, spacig=3, horizontal=TRUE, container=TiltGroup)
      enabled(FitMethod) <- FALSE

      StartCalc <- gbutton("START ESTIMATION", handler=function(...){
                               Model <- svalue(DP.Model, index=TRUE)
                               if ( Model == 3 ){  #N.CL>3, only MaxEnt method available
                                  cat("\n ==> Preparing DataSets for MaxEnt Analysis.")
                                  rho <<- as.numeric(svalue(FilmDens))
                                  XrayE <<- as.numeric(svalue(Xenergy))
                                  N.Layers <<- as.numeric(svalue(Layers))
                                  L.Thickness <<- as.numeric(svalue(LThckns))
                                  Intensity <<- matrix(nrow=N.CL, ncol=N.XS) #dimensions[Elements, Tilt] are as in Livesey
                                  for (ii in 1:N.XS){
                                        Tilt[ii] <<- as.numeric(svalue(angles[[ii]]))
                                  }
                                  Tilt <<- Tilt*pi/180  #degree to radian conversion
                                  if (is.na(rho) || is.na(XrayE) || is.na(sum(Tilt)) || is.na(N.Layers) || is.na(L.Thickness)){  #if tilt contains NA (i.e. lacking angle) sum(tilt)==NA
                                        gmessage(msg="INFORMATION LACKING! Please control and fill in all the items", title = "WARNING", icon="warning")
                                        return()
                                  }
                                  idx <- order(Tilt, decreasing=TRUE)
                                  Tilt <<- Tilt[idx] #tilt angles ordered in decreasing order
                                  SelectedFName <<- SelectedFName[idx] #order the selected FNames as the Tilt order
                                  SpectDataX <- lapply(SpectDataX, function(z) z<-z[idx])  #apply tilt ordering to SpectDataX components
                                  SpectDataY <- lapply(SpectDataY, function(z) z<-z[idx])  #apply tilt ordering to all  SpectDataY components
                                  #SpectDataRSF: RSF is the same for all tilt angles, just 1 value for each element

                                  N.Tilt <- length(Tilt)
                                  if (N.Tilt != N.XS){
                                     gmessage(msg="PROBLEM! N. XPS-Samples does not correspond to N. tilt angles. Analysis stopped.", title="N. TILT ANGLES INCONGUENCE", icon="error")
                                     return()
                                  }
                                  N.Layers <- as.integer(svalue(Layers))
                                  Exp.Data <<- matrix(data=NA, nrow=N.CL, ncol=N.XS)
                                  cat("\n ==> Estimation of Element Attenuation Length")
                                  for(jj in 1:N.CL){
                                      PosMax <- NULL
                                      for(ii in 1:N.XS){ #this for runs on the N. XPSSample == N. tilt angles
                                          Intensity[jj, ii] <<- sum(SpectDataY[[jj]][[ii]])

                                          Exp.Data[jj, ii] <<- sum(SpectDataY[[jj]][[ii]]) #SpectDataY[[jj]][[ii]]: jj runs on N.Corelines = 2, ii runs on N XPSSample = N.XS == Ntilt angles
                                          MM <- max(SpectDataY[[jj]][[ii]])         #max of selected corelines1 or fit component1
                                          MM <- findYIndex(SpectDataY[[jj]][[ii]], MM, 0.001)
                                          MM <- MM[1]    #the index relative to the max CL can be a vector depending spectral noise and if precision 0.001 is enough
                                          PosMax[ii] <- SpectDataX[[jj]][[ii]][MM]  #energy correspondent to the max of the jjth CL or position of the jjth fit component 1
                                      }                                             #see pp. 36 in Briggs "Surface analysis of polymers and static SIMS" Cambridge Uni. Press 1998 Cambridge
                                      PosMax <- mean(PosMax) #this is the mean value of the energy position of coreline jj calculated for the different tilt angles
                                      if (BE == TRUE) {PosMax <- XrayE - PosMax} # If BE scale transform in KE!
                                      Att.L[jj] <<- (1000/rho)*( 49/PosMax^2 + 0.11*PosMax^0.5) #attenuation length at energy correspondent to jj coreline position.
                                  }
                                  cat("\n ==> Estimation of Core-Line noise.")
                                  Variance <<- Calc.Var()


                                  tmp <- NULL
                                  cat("\n ==> Estimation of Concentrations.")
                                  Conc <<- matrix(nrow=N.CL, ncol=N.XS)
                                  for(ii in 1:N.XS){
                                      FName <- get(SelectedFName[ii], envir=.GlobalEnv)
                                      tmp <- XPSquantify(FName, verbose=TRUE)
                                     	TotQuant <- sum(unlist(sapply(tmp, function(x) x$quant)))
        	                             Conc[,ii] <<- sapply(tmp, function(x) {sum(unlist(x$quant))/TotQuant*100})
                                  }
                                  cat("\n ==> Call Max-Entropy procedure")

MaxEntData <- list(Variance=Variance, Concent=Conc, Intensity=Intensity, RSF=unlist(SpectDataRSF), Att.L=Att.L, Tilt=Tilt, N.Layers=N.Layers, L.Thickness=L.Thickness, SelectedFName=SelectedFName, SelectedCL=SelectedCL)
cat("\n *************** MaxEnt List **************\n")
print(MaxEntData)
assign("MaxEntData", MaxEntData, envir=.GlobalEnv)  #save the xxx.Rdata XPSSample in the .GlobalEnv
save(list="MaxEntData", file="Z:/X/LAVORI/R/Analysis/test/DepthProf/MaxEntData.Rdata", compress=TRUE)
                                  XPSMaxEnt(Variance, Intensity, Conc, Att.L, Tilt, N.Layers, L.Thickness)
                               }

                               if (IniData == FALSE){ #initializes data for thickness estimation
                                  if (length(rho)==0 || length(XrayE)==0 || is.na(sum(Tilt))){  #if Tilt, rho, XrayE not already initialized
                                     rho <<- as.numeric(svalue(FilmDens))
                                     XrayE <<- as.numeric(svalue(Xenergy))
                                     for (ii in 1:N.XS){
                                        Tilt[ii] <<- as.numeric(svalue(angles[[ii]]))
                                     }
                                     Tilt <<- Tilt*pi/180  #degree to radian conversion
                                     if (is.na(rho) || is.na(XrayE) || is.na(sum(Tilt)) ){  #if tilt contains NA (i.e. lacking angle) sum(tilt)==NA
                                        gmessage(msg="INFORMATION LACKING! Please control and fill in all the items", title = "WARNING", icon="warning")
                                        return()
                                     }
                                  }
                                  idx <- order(Tilt, decreasing=TRUE)
                                  Tilt <<- Tilt[idx] #tilt angles ordered in decreasing order
                                  #order spectral data correspondently
                                  SpectDataX[[1]] <<- SpectDataX[[1]][idx] #correspondent SpectralData applying the permutations done for Tilt vector
                                  SpectDataX[[2]] <<- SpectDataX[[2]][idx]
                                  SpectDataY[[1]] <<- SpectDataY[[1]][idx]
                                  SpectDataY[[2]] <<- SpectDataY[[2]][idx]
                                  SelectedFName <<- SelectedFName[idx] #also the list of XPSSamples is ordered with correspondence to Tilt angles
                                  for (ii in 1:N.XS){
                                      AreaComp1[ii] <- sum(SpectDataY[[1]][[ii]]) # SpectData[jj, ii]: jj runs on N.Corelines = 2, ii runs on N XPSSample = N.XS
                                      AreaComp2[ii] <- sum(SpectDataY[[2]][[ii]]) # This is the integral of selected Coreline 2 of the XPSSample ii
                                      Ratio[ii] <<- AreaComp1[ii]/AreaComp2[ii]
                                      MM <- max(SpectDataY[[1]][[ii]])       #max of selected corelines1 or fit component1
                                      MM <- findYIndex(SpectDataY[[1]][[ii]], MM, 0.001)
                                      MM <- MM[1] #index relative to the max CL1 can be a vector depending spectral noise and if precision 0.001 is enough

                                      PosMax1[ii] <<- SpectDataX[[1]][[ii]][MM]  #energy correspondent to the max CL1 or position fit component 1
                                      MM <- max(SpectDataY[[2]][[ii]])       #max of selected corelines2 or fit component 2
                                      MM <- findYIndex(SpectDataY[[2]][[ii]], MM, 0.001)
                                      MM <- MM[1]

                                      PosMax2[ii] <<- SpectDataX[[2]][[ii]][MM]  #energy correspondent to the max CL2 or position fit component 2
                                  }
                                  PosMax1 <<- round(mean(PosMax1), 2)    # PosMax1 is the average position of the Coreline1 (element1) acquired at different tilt angles
                                  PosMax2 <<- round(mean(PosMax2), 2)    # PosMax2 is the average position of the Coreline2 (element2) acquired at different tilt angles
                                  RSF1 <- SpectDataRSF[[1]]         # Sensitivity Factor of element 1
                                  RSF2 <- SpectDataRSF[[2]]         # Sensitivity Factor of element 2

                                  # the estimation of the film thickness requires the intensities from surface and bulk elements
                                  # measured under identical conditions i.e. the ratio of intensities from samples composed by only
                                  # by element 1 or by element 2= R.Rsf  is represended by the ratio of the element RSF
                                  R.Rsf <<- RSF1/RSF2

                                  # the estimation of the film thickness is performed considering that the intensities
                                  # of element1 and element2 decay with an exponential law with the depth.
                                  # The ratio of the intensities I1/I2 with I1=intensity of element1 placed in the bulk
                                  # and I2=intensity of element2 placed in the film on the surface. Following
                                  # D. Briggs in "Surface analysis of polymers by XPS and static SIMS", Cambridge University
                                  # press (1998):
                                  #
                                  #   I1/I2 = R.Rsf * {exp[-d/L1 sin(t) ]} / {1 - exp[d/L2 sin(t) ]}   (1)
                                  #
                                  # where t=tilt angle. This equation has an easy solution only if element1 and element2 are
                                  # two components of the same chemical element for example CHx and CF in a fluorurated polymer.
                                  # In this case R.Rsf = 1  and  L1 = L2 and d may be estimaded.
                                  # If element1 differs from element2,  R.Rsf != 1,  L1 != L2
                                  # the equation can be solved numerically by fitting the intensity ratio I1/I2 by varying the
                                  # parameter d = film thickness.
                                  #Now control if the Ratio increases/decreases with the tilt angle
                                  #This tells us which element (Coreline) is on the surface anwhich below
                                  SelectedSpect <<- SelectedCL
                                  if (N.CL ==1) {
                                      SelectedSpect[1] <<- paste(SelectedCL[1], SelectedCL[3], sep=".")
                                      SelectedSpect[2] <<- paste(SelectedCL[1], SelectedCL[4], sep=".")
                                  }
                                  tmp <- sort(Ratio, decreasing=TRUE)
                                  tmp <- (sum(Ratio - tmp)) # sum == 0  Ratio is equal to tmp which decreases with decreasing the tilt angle 90 -> 10 deg
                                  if (tmp == 0 ) {          # tmp = 0 means  original CL2 is on the surface and original CL1 in the bulk
                                     txt <- paste(SelectedSpect[2], " is on the surface and ", SelectedSpect[1], " is below the surface", sep="")
                                  } else {
                                     tmp <- sort(Ratio, decreasing=FALSE)
                                     tmp <- (sum(Ratio - tmp)) # sum == 0 if Ratio increases with decreasing the tilt angle 90 -> 10 deg
                                     if (tmp == 0 ) {          # tmp = 0 means now original CL1 is on the surface and original CL2 in the bulk
                                        SelectedSpect[c(1,2)] <<- SelectedSpect[c(2,1)] #swaps the elements of SelectedSpect to modify the warning message
                                        txt <- paste(SelectedSpect[2], " is on the surface and ", SelectedSpect[1], " is below the surface", sep="")
                                        Ratio <<- 1/Ratio   #reciprocal Ratio => CL1 on the surface, CL2 in the bulk => can use the same fitting function as the previous case
                                     } else {
                                        txt <- paste("Bad Ratio ", SelectedSpect[1], "/", SelectedSpect[2], "! Not possile to identify which element is on the surface \n  RESTART THE ANALYSIS SELECTING DIFFERENT ELEMENTS!", sep="")
                                        gmessage(msg=txt, title="BAD CORE LINE SELECTION", icon="error")
                                        SelectedCL <<- NULL
                                        Tilt <<- NULL
                                        PosMax1 <<- NULL
                                        PosMax2 <<- NULL
                                        Ratio <<- NULL
                                        return()
                                     }
                                  }
                                  par(mfrow=c(1,1)) #set single panel plot
                                  plot(Tilt*180/pi, Ratio, type="b", pch=16, cex=2, col="blue", xlab="Tilt Angle (deg)", ylab=paste(SelectedSpect[1], "/", SelectedSpect[2]))
                                  txt <- paste("Look at the plot and confirm if: \n", txt, sep="")
                                  answ <- gconfirm(msg=txt, title="CONFIRM ELEMENT POSITION", icon="warning")
                                  if (answ == FALSE) {
                                     gmessage(msg="Oops... Wrong analysis. Unable to continue. Sorry!", title="WRONG ANALYSIS", icon="warning")
                                     return()
                                  } else {
                                     E.Units <- "Kinetic Energy [eV]"
                                     if (BE == TRUE) {              #binding energy set
                                         PosMax1 <- XrayE-PosMax1   #now posmax  is in kinetic energy
                                         PosMax2 <- XrayE-PosMax2   #now posmax  is in kinetic energy
                                         E.Units <- "Binding [eV]"
                                         enabled(SaveAndExit) <- TRUE # Saving data blocked: ctrls on Dest file needed
                                     }
                                     La <<- (1000/rho)*( 49/PosMax1^2 + 0.11*PosMax1^0.5) #attenuation length at energy correspondent to coreline1 position. Cannot use L1 name
                                     Lb <<- (1000/rho)*( 49/PosMax2^2 + 0.11*PosMax2^0.5) #attenuation length at energy correspondent to coreline1 position. Cannot use L2 name

                                     #call fitting routine
                                     DataToFit <<- data.frame(x=Tilt, y=Ratio)
                                     Fit <<- FitData()
                                     IniData <<- TRUE
                                  }
                               } else {
                                  Fit <<- FitData()   # This allows repeating the fit selecting other fitting models
                               }  # end of  if(IniData==FALSE)
                         }, container=TiltGroup)
      enabled(StartCalc) <- FALSE

#--- COMMON BUTTONS


      gbutton("RESET",  handler=function(h, ...) {
                               svalue(SourceFiles) <- NULL
                               svalue(CoreLineList) <- NULL
                               svalue(CommonCoreLines) <- NULL
                               svalue(FilmDens) <- NULL
                               svalue(Xenergy) <- NULL
                               for (ii in 1:N.XS){
                                   Tiltlayout[ii, 1] <<- glabel(" ", container=Tiltlayout)
                                   Tiltlayout[ii, 2] <<- angles[[ii]] <<- gedit("", initial.msg="Tilt ang.", container=Tiltlayout)
                               }

                               FNameList<-XPSFNameList()
                               if (length(FNameList) == 0) { return() }
                               SampID <<- ""
                               SpectList <<- ""
                               SourceFileList <<- NULL
                               SelectedFName <<- NULL
                               SpectDataX <<- list()
                               SpectDataY <<- list()
                               SpectDataRSF <<- list()
                               Exp.Data <<- NULL
                               Att.L <<- NULL
                               Variance <<- NULL

                               CLlist <<- list()
                               ComponentCK <<- list()
                               CommonCL <<- XPSSpectList(FNameList[1]) #the list of Corelines of the first XPSSample used as reference
                               SelectedCL <<- NULL
                               SelectedSpect <<- NULL
                               N.XS <<- NULL     #number of selected XPSSamples
                               N.CL <<- NULL     #number of selected common lines

                               Tilt <<- NULL
                               angles <<- list()
                               infoLab <<- list()
                               AreaComp1 <<- NULL
                               AreaComp2 <<- NULL
                               Ratio <<- NULL
                               PosMax1 <<- NULL
                               PosMax2 <<- NULL
                               rho <<- NULL
                               XrayE <<- NULL
                               N.Layers <<- NULL
                               BE <<- NULL
                               La <<- NULL    #attenuation length element1    cannot use the name L1 problems in XPSFitLM
                               Lb <<- NULL    #attenuation length element2    cannot use the name L2 problems in XPSFitLM
                               R.Rsf <<- NULL
                               FitEstimation <<- NULL
                               FitCurve <<- NULL
                               IniParms <<- data.frame()
                               FitParms <<- NULL
                               IniData <<- FALSE
                               svalue(NoteBK) <- 1 #go to the first page
                         },container=Addgroup)


      SaveAndExit<-gbutton("SAVE DEPTH PROFILE RESULTS", handler=function(h, ...){
                                  SaveResults()
                                  XPSSaveRetrieveBkp("save")
                         }, container=Addgroup)
      enabled(SaveAndExit) <- FALSE # Saving data blocked: ctrls on Dest file needed

      gbutton("EXIT",  handler=function(h, ...) {
                                  dispose(AddWin)
                                  XPSSaveRetrieveBkp("save")
                         },container=Addgroup)


#DestFName XPSSample modified with baselines if not present


      svalue(NoteBK) <- 1
      visible(AddWin) <- TRUE

}


