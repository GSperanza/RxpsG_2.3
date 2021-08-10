#function to remove noise from XPS-Sample spectra

#'Smoothing functions to remove noise from spectra.
#'
#'XPSFilter contains a list of different filters to remove noise from spectra.
#'
#'In Sawitzy Golay filter at point i is a weighted average of data i-n... i ... i+n
#'In Autoregressive filters the output variable depends linearly on its own previous values. 
#'In the Moving Average filters the output at point i is the average of data i-n ... i ... i+n .
#'The FFT filter applys the FFT transform to perform filtering.
#'The Wavelets fileter uses the wavelets to perform filtering.
#'The FIRfilter is a Finite Impulse Response with a zero distortion.
#'The Butterworth is an Infinite Impulse Response filter.
#'No parameters are passed to this function
#'@examples
#'
#'\dontrun{
#'	XPSFilter()
#'}
#'
#'@export
#'


XPSFilter <- function() {

   import::here(filter, sgolay, hamming, fir1, filtfilt, butter, .from=signal)
   import::here(mra, .from=wavelets)


   PlotData <- function() {
      XXX <- cbind(unlist(Object@.Data[1]))
      YYY <- cbind(unlist(Object@.Data[2]),Filtrato)
      if (BkgSubtr) {
         YYY <- cbind(YYY, BackGnd)
      }
      Xlim <- range(XXX)
      if (Object@Flags[1]==TRUE) {
         Xlim <- rev(Xlim)  ## reverse x-axis
      }
      Ylim <- range(YYY)

      matplot(x=XXX, y=YYY, type="l", lty="solid", col=c("black", "red", "green"),
              xlim=Xlim, ylim=Ylim, xlab=Object@units[1], ylab=Object@units[2])
      return()
   }



   BkgSubtraction <- function(data){
      BackGnd <<- NULL
      LL <- length(data)
      rng <- floor(LL/15)
      if (rng<3) {rng==3}
      if (rng>30) {rng==30}
      bkg1 <- mean(data[1:rng])
      bkg2 <- mean(data[(LL-rng):LL])
      stp <- (bkg1-bkg2)/LL
      Dta_Bkg <- NULL
      for (ii in 1:LL){
          Dta_Bkg[ii] <- data[ii]-(bkg1-ii*stp)
          BackGnd[ii] <<- bkg1-ii*stp
      }
      return(Dta_Bkg)
   }

   normalize <- function(data){
      coefficenti <- sgolay(p=1, n=(2*7+1), ts=1) #7 coefficient for Sawitzky-Golay filtering
      LL <- length(data)
      rng <- as.integer(LL/20)
      if (rng < 7) rng <- 7 #minimum interval length for filtering with 7 coeff.

      YYY <- unlist(data)
      Min <- min(YYY)       #compute averaged noise free value around the data minimum

      FFF <- filter(filt=coefficenti, x=YYY) #Savitzky-Golay filtering
      Min <- min(FFF)
      Max <- max(FFF)
      if (SaveAmpli) AmpliRng <<- c(Min, Max)
      YYY <- (YYY-Min)/(Max-Min)
      return(YYY)
   }

   resize <- function(data){
      SaveAmpli <<- FALSE         #preserve ampliture range of original data
      data <- normalize(data)
      Min <- AmpliRng[1]
      Max <- AmpliRng[2]
      data <- (data*(Max-Min)) + Min
      return(data)
   }


   SetAmplitude <- function(){
#amplitude filtered data does not correspond to that of original data: matching is needed
      LL <- length(Object@.Data[[2]])
      BkgSubtraction(Object@.Data[[2]])  #To Define just the background
      coefficenti <- sgolay(p=1, n=(2*7+1), ts=1) #7 coefficient for Sawitzky-Golay filtering
#Padding data at edges: adds a number of data corresponding to the N. coeff Sawitzki-Golay filter.
#to avoid negative indexing when min(data( corresponds to 1 or LL
      BackGndO <- BackGnd
      XXX <- Object@.Data[[2]]-BackGndO
      pad <- rep(mean(XXX[1:5], 7))
      XXX <- c(pad, XXX)
      pad <- rep(mean(XXX[(LL-5):LL], 7))
      XXX <- c(XXX, pad)

      Min <- min(XXX)      #minimum of BKG subtracted original data Object@.Data[[2]]
      indx <- which(XXX == Min)
      rng <- 10
      while ( (indx+rng)>LL || (indx-rng)<0) { rng <- rng-1 }
      FFF <- filter(filt=coefficenti, x=XXX[(indx-rng):(indx+rng)]) #Savitzky-Golay filtering  BKG subtr original data
      Min <- min(FFF)

      Max <- max(XXX)      #maximum of BKG subtracted original data Object@.Data[[2]]
      indx <- which(XXX == Max)
      rng <- 10
      while ( (indx+rng)>LL || (indx-rng)<0) { rng <- rng-1 }
      FFF <- filter(filt=coefficenti, x=XXX[(indx-rng):(indx+rng)]) #Savitzky-Golay filtering
      Max <- max(FFF)
      AmpliFact <- Max-Min              #amplification factor

#Do the same for filtered data
      XXX <- BkgSubtraction(Filtrato)  #To Define just the background
      BackGndF <- BackGnd
      minF <- min(XXX)      #minimum of BKG subtracted filtered data
      indx <- which(XXX == minF)
      rng <- 5
      while ( (indx+rng)>LL || (indx-rng)<0) { rng <- rng-1 }
      minF <- mean(XXX[(indx-rng):(indx+rng)])  #average of BKG subtracted original data around minimum value

      maxF <- max(XXX)      #maximum of BKG subtracted original data
      indx <- which(XXX == maxF)
      rng <- 5
      while ( (indx+rng)>LL || (indx-rng)<0) { rng <- rng-1 }
      maxF <- mean(XXX[(indx-rng):(indx+rng)])  #average of BKG subtracted original data around maximum value
      DampFact <-  maxF-minF                      #damping factor

#Filtered data may contain noise: is needed a mean value around the maximum
      if (! BkgSubtr){
         Dta_Amp <- XXX*AmpliFact/DampFact+BackGndO #match amplitude of filtered data on BKG-UNSUBTRACTED original data
      } else {
         Dta_Amp <- XXX*AmpliFact/DampFact #match amplitude of filtered data on UNSUBTRACTED original data
      }

      return(Dta_Amp)
   }

   SavGolay <- function() {
      BkgSubtr <<- svalue(BkgSubtr2)
      FiltOrder <- as.numeric(svalue(F2FiltDeg))
      DaFiltrare <<- unlist(Object@.Data[[2]])
      Filtrato <<- NULL
      coefficenti <- NULL
      coefficenti <- sgolay(p=1, n=(2*FiltOrder+1), ts=1)
      LL <- length(DaFiltrare)
      if (! BkgSubtr) { # background substraction NOT selected
         Filtrato <<- filter(filt=coefficenti, x=DaFiltrare) #Savitzky-Golay filtering
      }
      if (BkgSubtr) {
         DaFiltrare <<- BkgSubtraction(DaFiltrare)
         Filtrato <<- filter(filt=coefficenti, x=DaFiltrare) #Savitzky-Golay filtering
      }
      PlotData()
      return()
   }


   AutoRegMedMob <- function() {
      BkgSubtr <<- svalue(BkgSubtr3)
      FilterType <- svalue(F3FiltType)
      Ro <- as.numeric(svalue(F3FiltDeg))
      coefficenti <- NULL
      Filtrato <<- NULL
      #preprocessing
      DaFiltrare <<- unlist(Object@.Data[[2]])
      LL <- length(DaFiltrare)
      if (BkgSubtr) { DaFiltrare <<- BkgSubtraction(DaFiltrare) }
      SaveAmpli <<- TRUE
      DaFiltrare <<- normalize(DaFiltrare)
      #now filtering
      tmp1 <- NULL
      tmp2 <- NULL
      if (FilterType=="AutoRegressive"){
#--- Autoregressive filter
         DaFiltrare <<- c(rep(DaFiltrare[1],2*Ro), DaFiltrare, rep(DaFiltrare[LL],2*Ro)) #2*Ro Padding at edges
         Filtrato <<- stats::filter(x=DaFiltrare, filter=rep(1,2*Ro),  sides=2, method="convolution", circular=TRUE)  #voglio usare FILTER pacchetto STATS non pacchetto SIGNAL
         Filtrato <<- Filtrato[(2*Ro+1) : (LL+2*Ro)]
      } else if (FilterType=="MovingAverage") {
#--- Moving average filter
         RRo <- Ro/(15+1)  #input degree of moving average filter [0.1 : 0.99]
         LL <- length(DaFiltrare)
         tmp1[1] <- DaFiltrare[1]
         tmp2[LL] <- DaFiltrare[LL]
         for(ii in 2:LL){
            tmp1[ii]=RRo*tmp1[ii-1] + (1-RRo)*DaFiltrare[ii];              #forward filtering
            tmp2[LL-ii+1]=RRo*tmp2[LL-ii+2]+(1-RRo)*DaFiltrare[LL-ii+1];   #backward filtering
         }
         Filtrato <<- (tmp1+tmp2)/2
      }
      #now set correct amplitude of filtered data
      Filtrato <<- resize(Filtrato)   #match the amplitude of Filtered with the original data
      PlotData()
      return()
   }

   FFTfilt <- function() {
      BkgSubtr <<- svalue(BkgSubtr4)
      FiltOrder <- as.numeric(svalue(F4FiltDeg))
      DaFiltrare <<- NULL
      Filtrato <<- NULL
      #preprocessing
      DaFiltrare <<- unlist(Object@.Data[[2]])
      LL <- length(DaFiltrare)
      if (BkgSubtr) DaFiltrare <- BkgSubtraction(DaFiltrare) #background subtraction to avoid FFT spourious oscillations
      SaveAmpli <<- TRUE
      DaFiltrare <<- normalize(DaFiltrare)
      #now filtering
      stopF <- floor(LL*0.015)
      FiltOrder <- as.integer(LL-LL*((15-FiltOrder)/15)^2 - stopF ) #non linear progrssion: for low filt ord FFTfilt does not work
#      FiltOrder <- FiltOrder*floor((LL-5)/10)
      fftTransf <- fft(DaFiltrare)
      fftTransf[(LL-FiltOrder):LL] <- 0+0i #FFT high frequencies components are forced to zero
      Filtrato <<- InvFFT <- Re(fft(fftTransf, inverse = TRUE)/length(fftTransf))
      #now set correct amplitude of filtered data
      Filtrato <<- resize(Filtrato)   #match the amplitude of Filtered with the original data
      PlotData()
      return()
   }


msSmoothMRA <- function(){
      DaFiltrare <<- NULL
      Filtrato <<- NULL
      WTdeg <- 2*as.numeric(svalue(F5WTDeg))  #This defines the N. Daubechies wavelets of the transform
      Nlevels <- as.numeric(svalue(F5WTForce))
      #preprocessing
      tmp <- (Object@.Data[[2]])
      if (BkgSubtr) tmp <- BkgSubtraction(tmp)
      SaveAmpli <<- TRUE
      tmp <- normalize(tmp)
      LL <- length(tmp)
      DaFiltrare[1:WTdeg] <<- rep(tmp[1], WTdeg) #head padding with WTdeg values
      DaFiltrare[(WTdeg+1):(WTdeg+LL)] <<- tmp
      DaFiltrare[(WTdeg+LL+1):(LL+2*WTdeg)] <<- rep(tmp[LL], WTdeg) #end padding with WTdeg values
      #now filtering
      WTfilt <- paste("d", WTdeg, sep="")  #d2 means Daubechies length=2,  d8 means Daubechies length=8
      ## maximum level according to the length(x) and selected filter degree WTdeg
#      max.levels <- floor( log((LL-1)/(WTdeg-1) + 1, base=2) )  #following wavelets::mra suggestions
#      max.levels <- floor(log(LL,2))  # max.levels following modwt
      max.levels <- 10 #fixed n.levels of the decomposition. This allows increasing the Nlevels (filter force) up to 10!

      if (Nlevels > max.levels){
         gmessage("Attention: Nlevel must be < 1/2 * Filter Order", title="BAD DEGREE OF NOISE REJECTION", icon="error")
      }

      mra.modwt <- mra(DaFiltrare
                        , filter = WTfilt
                        , n.levels = max.levels
                        , method = "modwt"
                        , boundary="periodic")

      mra.mat <- sapply(c(mra.modwt@D[Nlevels], mra.modwt@S[Nlevels]), cbind)

      Filtrato <<- rowSums(mra.mat)
      Filtrato <<- Filtrato[(WTdeg+1):(WTdeg+LL)] #eliminates head and tail padded values
      Filtrato <<- resize(Filtrato)
      PlotData()
      return()
   }


   WTfilt <- function(){   #this filter does not work properly
      BkgSubtr <<- svalue(BkgSubtr5)
      DaFiltrare <<- NULL
      Filtrato <<- NULL
      Force <- as.numeric(svalue(F5WTForce))
      Threshold <- 1/(as.numeric(svalue(F5WTThreshold)))  #if the degree of filtering increases the treshold must decrease
      WaveletType <- "d6"  # 6 Daubechies functions are initially selected to model the original data
#The difference between 6 and 20 Daubechies in modeling the XPS spectra is negligible
#Also the difference between the different types of Waveletts is negligible => use of Daubechies functions
#For this reason: selection on the number and type of Waveletts was removed
#BY DEFAULT: WaveletType=Daubechies, WavletNum=6
      DaFiltrare <<- (Object@.Data[[2]])
      #preprocessing
      DaFiltrare <<- BkgSubtraction(DaFiltrare)
      if (BkgSubtr) { DaFiltrare <<- BkgSubtraction(DaFiltrare) }
      minD <- min(DaFiltrare)
      maxD <- max(DaFiltrare)
      DaFiltrare <<- (DaFiltrare-minD)/(maxD-minD)
      LL <- length(DaFiltrare)
      NLevel <- floor(log(LL,2))
#      if (Force=="Light") {Force <- 1} #Wavelet amplitude is decreased from level N.Level=1
#      if (Force=="Mild") {Force <- as.numeric(NLevel/2)} #Wavelet amplitude is decreased from level N.Level/2
#      if (Force=="Strong") {Force <- NLevel} #Wavelet amplitude is decreased from beginning level N.Level
      Force <- Force*NLevel/10
      DataWT <- modwt(DaFiltrare,filter=WaveletType,n.levels=NLevel,boundary="reflection", fast=TRUE)
#here wavelets are modeling the noise then
#setting the filter strength corresponds to set the level where to start decreasing the wavelets amplitude
#the filter strength varies from 1 to Nlevel
#Treshold describes how much the wavelets amplitude has to be reduced. Treshold varies from 0.1 to 1
#here set level in the range [1:Nlevel] and threshold in [0.1:1]
      for (ii in 1:NLevel) {
          if (ii >= Force){
             DataWT@W[[ii]] <- DataWT@W[[ii]]/Threshold
          }
      }
      Filtrato <<- imodwt(DataWT, fast=TRUE)  #After wavelet re-modulation perform the inverse transform

      minF <- min(Filtrato)
      maxF <- max(Filtrato)
      Filtrato <<- (Filtrato-minF)*(maxD-minD)/(maxF-minF)+minD#adjust the amplitude of the filtered data to that of original-BKG
      if (! BkgSubtr) {                 # if background subtraction option NOT selected
         Filtrato <<- Filtrato + BackGnd #add BKG if option BKG substraction not selected
      }
      PlotData()
      return()
   }


   FirFilt <- function() {
      BkgSubtr <<- svalue(BkgSubtr6)
      FiltOrder <- as.numeric(svalue(F6FiltOrder))
      CutOFF <- as.numeric(svalue(F6CutOFF))
      CutOFF <- (16-CutOFF)/15    #here cut-off were freq/Nyquist for example 40/500Hz (sampling freq = 1kHz) where 40Hz is adapted to the typical XPS noise
      coefficenti <- NULL
      DaFiltrare <<- NULL
      Filtrato <<- NULL
      #preprocessing
      tmp <- unlist(Object@.Data[[2]])
      if (BkgSubtr) tmp <- BkgSubtraction(tmp)
      SaveAmpli <<- TRUE
      tmp <- normalize(tmp)
      LL <- length(tmp)
      DaFiltrare[1:FiltOrder] <<- rep(tmp[1], FiltOrder) #initial padding with FiltOrder values
      DaFiltrare[(FiltOrder+1):(FiltOrder+LL)] <<- tmp
      DaFiltrare[(FiltOrder+LL+1):(LL+2*FiltOrder)] <<- rep(tmp[LL], FiltOrder) #end padding with FiltOrder values
      #now filtering
      coefficenti <- fir1(n=FiltOrder, w=CutOFF, type = "low", window = hamming(FiltOrder+1), scale = TRUE)
      Filtrato <<- filtfilt(filt=coefficenti, x=DaFiltrare)
      Filtrato <<- Filtrato[(FiltOrder+1):(FiltOrder+LL)]    #eliminates initial and end padded values
      #now set correct amplitude of filtered data
      Filtrato <<- resize(Filtrato)
      PlotData()
      return()
   }


   ButtFilt <- function() {
      BkgSubtr <<- svalue(BkgSubtr7)
      FiltOrder <- as.numeric(svalue(F7FiltOrder))
      CutOFF <- as.numeric(svalue(F7CutOFF))
      CutOFF <- (16-CutOFF)/15    #here cut-off were freq/Nyquist for example 40/500Hz (sampling freq = 1kHz) where 40Hz is adapted to the typical XPS noise
      coefficenti <- NULL
      DaFiltrare <<- NULL
      Filtrato <<- NULL
      #preprocessing
      tmp <- unlist(Object@.Data[[2]])
      if (BkgSubtr) tmp <- BkgSubtraction(tmp)
      LL <- length(tmp)
      ZP <- 2*FiltOrder
      DaFiltrare[1:ZP] <<- rep(tmp[1], ZP) #initial padding with FiltOrder values
      DaFiltrare[(ZP+1):(ZP+LL)] <<- tmp
      DaFiltrare[(ZP+LL+1):(LL+2*ZP)] <<- rep(tmp[LL], ZP) #end padding with FiltOrder values

      #now filtering
      coefficenti <- butter(n=FiltOrder, W=CutOFF, type = "low", window = hamming(FiltOrder+1), scale = TRUE)
      Filtrato <<- filtfilt(filt=coefficenti, x=DaFiltrare)
      Filtrato <<- Filtrato[(ZP+1):(ZP+LL)]    #eliminates initial and end padded values
      #now set correct amplitude of filtered data
      PlotData()
      return()
   }

   makeCombo <- function(){
      ActiveFName <<- svalue(F1XpsSpect)
      FName <<- get(ActiveFName,envir=.GlobalEnv)  #carico in SampID il relativo XPSSAmple
      SpectList <<- XPSSpectList(ActiveFName)
      SpectIndx <<- 1
      delete(F1frame2, F1CoreLine)
      F1CoreLine <<- gcombobox(SpectList, selected=-1, editable=FALSE, handler=function(h,...){
                        XPSComponent <- svalue(F1CoreLine)
                        XPSComponent <- unlist(strsplit(XPSComponent, "\\."))   #skip the ".NUMBER" at beginning CoreLine name
                        SpectIndx <- as.integer(XPSComponent[1])
                        SpectName <- XPSComponent[2]
                        assign("activeFName", ActiveFName,.GlobalEnv) #set the activeSpectral name to the actual CoreLine name
                        assign("activeSpectName", SpectName,.GlobalEnv) #set the activeSpectral name to the actual CoreLine name
                        assign("activeSpectIndx", SpectIndx,.GlobalEnv) #set the activeSpectIndex to the actual CoreLine
                        Object <<- FName[[SpectIndx]]
                        plot(FName[[SpectIndx]])
                        enabled(F2group1) <- TRUE
                        enabled(F3group1) <- TRUE
                        enabled(F4group1) <- TRUE
                        enabled(F5group1) <- TRUE
                        enabled(F6group1) <- TRUE
                        enabled(F7group1) <- TRUE
                        enabled(F2filter) <- FALSE
                        enabled(F3filter) <- FALSE
                        enabled(F4filter) <- FALSE
                        enabled(F5filter) <- FALSE
                        enabled(F6filter) <- FALSE
                        enabled(F7filter) <- FALSE
                     }, container=F1frame2)
      plot(FName)
   }



#======== Variabili ========
   if (is.na(activeFName)){
       gmessage("No data present: please load and XPS Sample", title="XPS SAMPLES MISSING", icon="error")
       return()
   }
   ActiveFName <- get("activeFName", envir=.GlobalEnv)  #cload the XPSSample name (string)
   FName <- get(ActiveFName, envir=.GlobalEnv)   #load the active XPSSample (data)
   FNameList <- XPSFNameList()                   #list of all XPSSamples loaded in .GlobalEnv
   FNameIdx <- grep(ActiveFName,FNameList)
   SpectIndx <- get("activeSpectIndx", envir=.GlobalEnv)#load the index of the active CoreLine
   Object <- FName[[SpectIndx]]
   SpectList <- XPSSpectList(ActiveFName)

   DaFiltrare <- NULL
   Filtrato <- NULL
   BackGnd <- NULL
   BkgSubtr <- FALSE
   AmpliRng <- NULL
   SaveAmpli <- FALSE
   filterType <- c("AutoRegressive", "MovingAverage")
   waveletType <- list("Daubechies", "LeastAsymmetric", "BestLocalized", "Coiflet")
   waveletType[["Daubechies"]] <- c(2,4,6,8,10,12,14,16,18,20)
   waveletType[["LeastAsymmetric"]] <- c(8,10,12,14,16,18,20)
   waveletType[["BestLocalized"]] <- c(14,18,20)
   waveletType[["Coiflet"]] <- c(6,12,18,24,30)
   waveNumber <- "  "
   FiltInfo <- NULL


   mainFwin <- gwindow("DATA FILTERING", visible=FALSE)
   mainGroup <- ggroup(spacing=1, horizontal=FALSE, container=mainFwin)

   mainFrame1 <- gframe(spacing=1, horizontal=FALSE, container=mainGroup)

#   F1frame1 <- gframe(label=" CORELINE SELECTION ", spacing=1, horizontal=FALSE, container=mainGroup)
   F1group2 <- ggroup(label="  ",horizontal=TRUE, container=mainFrame1)

   F1frame1 <- gframe(" Select XPS Sample ", spacing=5, container=F1group2)
   F1XpsSpect <- gcombobox(FNameList, selected=FNameIdx, editable=FALSE, handler=function(h,...){
                            makeCombo()
                 }, container=F1frame1)

   F1frame2 <- gframe(" Select Coreline ", spacing=5, container=F1group2)
   F1CoreLine <- gcombobox(SpectList, selected=-1, editable=FALSE, handler=function(h,...){
                           XPSComponent <- svalue(F1CoreLine)
                           XPSComponent <- unlist(strsplit(XPSComponent, "\\."))   #skip the ".NUMBER" at beginning CoreLine name
                           SpectIndx <- as.integer(XPSComponent[1])
                           SpectName <- XPSComponent[2]
                           assign("activeSpectName", SpectName,.GlobalEnv) #set the activeSpectral name to the actual CoreLine name
                           assign("activeSpectIndx", SpectIndx,.GlobalEnv) #set the activeSpectIndex to the actual CoreLine
                           Object <<- FName[[SpectIndx]]
                           plot(FName[[SpectIndx]])
                           enabled(F2group1) <- TRUE
                           enabled(F3group1) <- TRUE
                           enabled(F4group1) <- TRUE
                           enabled(F5group1) <- TRUE
                           enabled(F6group1) <- TRUE
                           enabled(F7group1) <- TRUE
                           enabled(F2filter) <- FALSE
                           enabled(F3filter) <- FALSE
                           enabled(F4filter) <- FALSE
                           enabled(F5filter) <- FALSE
                           enabled(F6filter) <- FALSE
                           enabled(F7filter) <- FALSE
                 }, container = F1frame2)

   F1frame3 <- gframe(   "Update XPS Sample List", spacing=5, container=F1group2)
   UpdateButt <- gbutton("            UPDATE            ", spacing=30, handler=function(h,...){
                          FNameList <<- XPSFNameList()                   #list of all XPSSamples loaded in .GlobalEnv
                          delete(F1frame1, F1XpsSpect)
                          F1XpsSpect <<- gcombobox(FNameList, selected=FNameIdx, editable=FALSE, handler=function(h,...){
                                                   makeCombo()
                                         }, container=F1frame1)
                          makeCombo()
                 }, container=F1frame3)
   addSpring(F1frame3)
   glabel(" ", container=mainFrame1)   #add line space


#===== NoteBook =====
   nb1 <- gnotebook(expand=FALSE, spacing=10, container = mainFrame1)

# --- Tab1 ---
   F2group1 <- ggroup(label=" SAVITZKY GOLAY ", spacing=1, horizontal=FALSE, container=nb1)
   F2group2 <- ggroup(label="  ", horizontal=TRUE, container=F2group1)
   F2frame1 <- gframe(text="Degree of Noise Rejection", spacing=5, container=F2group2)
   F2FiltDeg <- gcombobox(c(1:15), selected=-1, editable=FALSE, handler=function(h,...){
                    enabled(F2filter) <- TRUE
                 }, container = F2frame1)

   BkgSubtr2 <- gcheckbox("BKG Subtraction", checked=FALSE, container=F2group2)

   F2group3 <- ggroup(label="  ",horizontal=TRUE, container=F2group1)
   F2filter <- gbutton(" FILTER ", handler=function(h,...){
                     SavGolay()
                     FiltInfo <<- paste("Savitzky Golay, ", svalue(F2FiltDeg), sep="")
                     enabled(SaveNewButt) <- TRUE
                     enabled(SaveButt) <- TRUE
                     enabled(SaveTestButt) <- TRUE
                 }, container = F2group3)

# --- Tab2 ---
   F3group1 <- ggroup(label=" AUTOREGRESSIVE : MOVING AVERAGE ", spacing=1, horizontal=FALSE, container=nb1)
   enabled(F3group1) <- FALSE
   F3group2 <- ggroup(label="  ", horizontal=TRUE, container=F3group1)
   F3frame1 <- gframe(text=" Select Filter Type ", spacing=5, container=F3group2)
   F3FiltType <- gcombobox(filterType, selected=-1, editable=FALSE, handler=function(h,...){
                    enabled(F3FiltDeg) <- TRUE
                 }, container = F3frame1)

   F3frame2 <- gframe(text="Degree of Noise Rejection", spacing=5, container=F3group2)
   F3FiltDeg <- gcombobox(c(1:15), selected=-1, editable=FALSE, handler=function(h,...){
                     enabled(F3filter) <- TRUE
                 }, container = F3frame2)

   BkgSubtr3 <- gcheckbox("BKG Subtraction", checked=FALSE, container=F3group2)

   F3group3 <- ggroup(label="  ",horizontal=TRUE, container=F3group1)
   F3filter <- gbutton(" FILTER ", handler=function(h,...){
                     AutoRegMedMob()
                     FiltInfo <<- paste(svalue(F3FiltType), ", ", svalue(F3FiltDeg), sep="")
                     enabled(SaveNewButt) <- TRUE
                     enabled(SaveButt) <- TRUE
                     enabled(SaveTestButt) <- TRUE
                 }, container = F3group3)

# --- Tab3 ---
   F4group1 <- ggroup(label=" FFT FILTER ", spacing=1, horizontal=FALSE, container=nb1)
   F4group2 <- ggroup(label="  ", horizontal=TRUE, container=F4group1)
   F4frame1 <- gframe(text="Degree of Noise Rejection", spacing=5, container=F4group2)
   F4FiltDeg <- gcombobox(c(1:15), selected=-1, editable=FALSE, handler=function(h,...){
                    enabled(F4filter) <- TRUE
                 }, container = F4frame1)

   BkgSubtr4 <- gcheckbox("BKG Subtraction", checked=FALSE, container=F4group2)

   F4group3 <- ggroup(label="  ",horizontal=TRUE, container=F4group1)
   F4filter <- gbutton(" FILTER ", handler=function(h,...){
                     FFTfilt()
                     FiltInfo <<- paste("FFT filter, ", svalue(F4FiltDeg), sep="")
                     enabled(SaveNewButt) <- TRUE
                     enabled(SaveButt) <- TRUE
                     enabled(SaveTestButt) <- TRUE
                 }, container = F4group3)

# --- Tab4 ---
   F5group1 <- ggroup(label=" WAVELETS FILTERING ", spacing=1, horizontal=FALSE, container=nb1)
   F5group2 <- ggroup(label="  ",horizontal=TRUE, container=F5group1)
   F5frame1 <- gframe(text=" Select Filter Order ", spacing=5, container=F5group2)
   F5WTDeg <- gcombobox(c(1:10), selected=5, editable=FALSE, container = F5frame1)

   F5frame2 <- gframe(text=" Degree of Noise Rejection ", spacing=5, container=F5group2)
   F5WTForce <- gcombobox(c(1:10), selected=-1, editable=FALSE, handler=function(h,...){
                     enabled(F5filter) <- TRUE
                 }, container = F5frame2)

   BkgSubtr5 <- gcheckbox("BKG Subtraction", checked=FALSE, container=F5group2)

   F5group2 <- ggroup(label="  ",horizontal=TRUE, container=F5group1)
   F5filter <- gbutton(" FILTER ", handler=function(h,...){
#                     WTfilt()
                     msSmoothMRA()
                     FiltInfo <<- paste("Wavelets filter, ", svalue(F5WTDeg), ", ", svalue(F5WTForce), sep="")
                     enabled(SaveNewButt) <- TRUE
                     enabled(SaveButt) <- TRUE
                     enabled(SaveTestButt) <- TRUE
                 }, container = F5group2)

# --- Tab5 ---
   F6group1 <- ggroup(label=" FIR FILTER ", spacing=1, horizontal=FALSE, container=nb1)
   F6group2 <- ggroup(label="  ",horizontal=TRUE, container=F6group1)
   F6frame1 <- gframe(text=" Select Filter Order ", spacing=5, container=F6group2)
   F6FiltOrder <- gcombobox(c(20,40,60), selected=1, editable=FALSE, handler=function(h,...){
                     enabled(F6CutOFF) <- TRUE
                 }, container = F6frame1)
   F6frame2 <- gframe(text="Degree of Noise Rejection", spacing=5, container=F6group2)
   F6CutOFF <- gcombobox(c(1:15), selected=-1, editable=FALSE, handler=function(h,...){
                     enabled(F6filter) <- TRUE
                 }, container = F6frame2)

   BkgSubtr6 <- gcheckbox("BKG Subtraction", checked=FALSE, container=F6group2)

   F6group3 <- ggroup(label="  ",horizontal=TRUE, container=F6group1)
   F6filter <- gbutton(" FILTER ", handler=function(h,...){
                     FirFilt()
                     FiltInfo <<- paste("FIR filter, ", svalue(F6FiltOrder), ", ", svalue(F6CutOFF), sep="")
                     enabled(SaveNewButt) <- TRUE
                     enabled(SaveButt) <- TRUE
                     enabled(SaveTestButt) <- TRUE
                 }, container = F6group3)

# --- Tab6 ---
   F7group1 <- ggroup(label=" BUTTERWORTH FILTER ", spacing=1, horizontal=FALSE, container=nb1)
   F7group2 <- ggroup(label="  ",horizontal=TRUE, container=F7group1)
   F7frame1 <- gframe(text=" Select Filter Order ", spacing=5, container=F7group2)
   F7FiltOrder <- gcombobox(c(4,6,8,10), selected=1, editable=FALSE, handler=function(h,...){
                    enabled(F7CutOFF) <- TRUE
                 }, container = F7frame1)
   F7frame2 <- gframe(text="Degree of Noise Rejection", spacing=5, container=F7group2)
   F7CutOFF <- gcombobox(c(1:15), selected=-1, editable=FALSE, handler=function(h,...){
                     enabled(F7filter) <- TRUE
                 }, container = F7frame2)

   BkgSubtr7 <- gcheckbox("BKG Subtraction", checked=FALSE, container=F7group2)

   F7group3 <- ggroup(label="  ",horizontal=TRUE, container=F7group1)
   F7filter <- gbutton(" FILTER ", handler=function(h,...){
                     ButtFilt()
                     FiltInfo <<- paste("Butterworth filter, ", svalue(F7FiltOrder), ", ", svalue(F7CutOFF), sep="")
                     enabled(SaveNewButt) <- TRUE
                     enabled(SaveButt) <- TRUE
                     enabled(SaveTestButt) <- TRUE
                 }, container = F7group3)


#--- Common buttons

   FButtgroup1 <- ggroup(horizontal=TRUE, spacing=1, container= mainGroup)
   FButtgroup2 <- ggroup(horizontal=TRUE, spacing=1, container= mainGroup)


   SaveNewButt <- gbutton(" SAVE AS A NEW CORE LINE ", spacing=1, handler=function(h,...){ #Filtered data saved in a new coreline
                    SpectIndx <- get("activeSpectIndx",.GlobalEnv)
                    Symbol <- get("activeSpectName",.GlobalEnv)
                    NCL <- length(FName)  #number of XPSSample corelines
                    CLNames <- names(FName)
                    FName[[NCL+1]] <<- FName[[SpectIndx]]
                    FName[[NCL+1]]@Symbol <<- Symbol
                    Info <- FName[[NCL+1]]@Info
                    nI <- length(Info)
                    Info[nI+1] <- paste("   Smoothing: ", FiltInfo, sep="")
                    FName[[NCL+1]]@Info <<- Info
                    FName@names <<- c(CLNames,Symbol)
                    FName[[NCL+1]]@.Data[[2]] <<- Filtrato
                    if (length(FName[[SpectIndx]]@RegionToFit$y>0)){ #RegionToFit defined
                       rng <- range(FName[[SpectIndx]]@RegionToFit$x)
                       if (Object@Flags[1]==TRUE) {rng <- rev(rng)}    #Binding energy set
                       idx1 <- which(FName[[TestIdx]]@.Data[[1]] == rng[1])
                       idx2 <- which(FName[[TestIdx]]@.Data[[1]] == rng[2])
                       FName[[TestIdx]]@RegionToFit$y <<- Filtrato[idx1:idx2]
                    }
                    assign(ActiveFName, FName,envir=.GlobalEnv)   #save XPSSample with additional smoothed coreline
#                    assign("activeSpectIndx", (NCL+1), envir=.GlobalEnv)      #set the activeSpectIndx be the smoothed core line
                    DaFiltrare <<- NULL
                    Filtrato <<- NULL
                    BackGnd <<- NULL
                    plot(FName)
                    XPSSaveRetrieveBkp("save")
                 }, container = FButtgroup1)

   SaveButt <- gbutton(" SAVE IN THE ORIGINAL CORE LINE", spacing=1, handler=function(h,...){  #Original data replaced by Filterd data
                    SpectIndx <- get("activeSpectIndx",.GlobalEnv)
                    FName[[SpectIndx]]@.Data[[2]] <<- Filtrato
                    if (length(FName[[SpectIndx]]@RegionToFit$y>0)){ #RegionToFit defined
                       rng <- range(FName[[SpectIndx]]@RegionToFit$x)
                       if (Object@Flags[1]==TRUE) {rng <- rev(rng)}    #Binding energy set
                       idx1 <- which(FName[[TestIdx]]@.Data[[1]] == rng[1])
                       idx2 <- which(FName[[TestIdx]]@.Data[[1]] == rng[2])
                       FName[[TestIdx]]@RegionToFit$y <<- Filtrato[idx1:idx2]
                    }
                    Info <- FName[[SpectIndx]]@Info
                    xxx <- grep("Smoothing:", Info)     #is the Test Filtering already present?
                    if (length(xxx)>0) {              #Filtering Comments are present?
                         nI <- length(Info)
                    } else {
                         nI <- length(Info)+1
                    }
                    Info[nI] <- paste("   Smoothing: ", FiltInfo, sep="")  #Add or Change filter information
                    FName[[SpectIndx]]@Info <<- Info
                    assign(ActiveFName, FName,envir=.GlobalEnv) #setto lo spettro attivo eguale ell'ultimo file caricato
                    Object <<- FName[[SpectIndx]]       #update Object used for filtering
                    plot(FName)
                    XPSSaveRetrieveBkp("save")
                 }, container = FButtgroup1)

   SaveTestButt <- gbutton(" SAVE AS A SMOOTHING-TEST ", spacing=1, handler=function(h,...){  #in Test Coreline new filtering is overwitten
                    SpectIndx <- get("activeSpectIndx",.GlobalEnv)
                    Symbol <- get("activeSpectName",.GlobalEnv)
                    CLNames <- names(FName)
                    LL <- length(FName)
                    chrPos <- regexpr("ST.", CLNames[SpectIndx])    #Are we working on a Smoothing-Test core line?
                    if (chrPos>0) {                               #Smoothing-Test Coreline IS PRESENT
                       TestIdx <- SpectIndx
                       Info <- FName[[TestIdx]]@Info             #update filter information
                       LL <- length(Info)
                       Info[LL] <- paste("   ::: Smoothing Test: ", FiltInfo, sep="")
                       FName[[TestIdx]]@Info <<- Info
                    } else {
                       chrPos <- which(regexpr("ST.", CLNames) >0)     #Find the index if there is a "ST." identifying the Smoothing-Test coreline
                       if (length(chrPos)==0) {                      #Smoothing-Test coreline is NOT present
                          TestIdx <- LL+1                           #Smoothing-Test coreline is added to FName as a new coreline
                          FName[[TestIdx]] <<- FName[[SpectIndx]] #We are testing a filter on a coreline and save results in the Smoothing-Test coreline
                          FName@names <<- c(CLNames,paste("ST.", Symbol, sep=""))  #modify the names and info of the new Smoothing-Test coreline
                          Info <- FName[[TestIdx]]@Info
                          LL <- length(Info)
                          if (regexpr(":::", Info[LL]) < 0){ LL <- LL+1 } #Smoothing-Test Info still not present
                          Info[LL] <- paste("   ::: Smoothing Test: ", FiltInfo, sep="")
                          FName[[TestIdx]]@Info <<- Info
                          FName[[TestIdx]]@Symbol <<- paste("ST.", Symbol, sep="")
                       } else {
                          TestIdx <- chrPos                            #Smoothing-Test coreline is present
                          Info <- FName[[TestIdx]]@Info
                          LL <- length(Info)
                          Info[LL] <- paste("   ::: Smoothing Test: ", FiltInfo, sep="")
                          FName[[TestIdx]]@Info <<- Info
                       }
                    }
                    FName[[TestIdx]]@.Data[[2]] <<- Filtrato         #change the original data with the filtered data
                    if (length(FName[[SpectIndx]]@RegionToFit$y>0)){ #RegionToFit defined
                       rng <- range(FName[[SpectIndx]]@RegionToFit$x)
                       if (Object@Flags[1]==TRUE) {rng <- rev(rng)}    #Binding energy set
                       idx1 <- which(FName[[TestIdx]]@.Data[[1]] == rng[1])
                       idx2 <- which(FName[[TestIdx]]@.Data[[1]] == rng[2])
                       FName[[TestIdx]]@RegionToFit$y <<- Filtrato[idx1:idx2]
                    }
                    assign(ActiveFName, FName,envir=.GlobalEnv)      #Save the modified XPSSample in the globalEnv
#                    assign("activeSpectIndx", TestIdx, envir=.GlobalEnv)  #set the activeSpectIndx be the smoothed core line
                    DaFiltrare <<- NULL
                    Filtrato <<- NULL
                    BackGnd <<- NULL
                    plot(FName)
                    XPSSaveRetrieveBkp("save")
                 }, container = FButtgroup1)

   gbutton(" RESET ", spacing=1, handler=function(h,...){
                    DaFiltrare <<- NULL
                    Filtrato <<- NULL
                    BackGnd <<- NULL
                    plot(FName)
                 }, container = FButtgroup1)

   gbutton(" EXIT ", spacing=1, handler=function(h,...){
                    dispose(mainFwin)
                    XPSSaveRetrieveBkp("save")
                 }, container = FButtgroup1)


   enabled(F2filter) <- FALSE
   enabled(F3filter) <- FALSE
   enabled(F4filter) <- FALSE
   enabled(F5filter) <- FALSE
   enabled(F6filter) <- FALSE
   enabled(F7filter) <- FALSE
   enabled(SaveNewButt) <- FALSE
   enabled(SaveButt) <- FALSE
   enabled(SaveTestButt) <- FALSE

   for (ii in 6:1){
      svalue(nb1) <- ii #update notebook pages
   }
   visible(mainFwin) <- TRUE
}
