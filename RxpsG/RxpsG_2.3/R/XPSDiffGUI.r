#function to differentiate XPS-Sample spectra

#'Differentiate is a function to compute the derivative of the XPSSpectra.
#'
#'@examples
#'
#'\dontrun{
#'	XPSDiff()
#'}
#' 
#'@export
#'


XPSDiff <- function(){

   Differ <- function(Object){
        LL <- length(Object)
        tmp <- NULL
        for(ii in 2:LL){
           tmp[ii] <- Object[ii] - Object[ii-1]
        }
        tmp[1] <- tmp[2]
        return(tmp)
   }

   BkgSubtraction <- function(data){   #linear BKG subtraction
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

   plotData <- function() {
      SpectIndx <- get("activeSpectIndx", envir=.GlobalEnv) #set the activeSpectIndex to the actual CoreLine
      XXX <- cbind(unlist(FName[[SpectIndx]]@.Data[1]))
#      YYY <- cbind(unlist(FName[[SpectIndx]]@.Data[2]),BackGnd+unlist(Differentiated[2]))
#      YYY <- cbind(unlist(FName[[SpectIndx]]@.Data[2]),BackGnd, unlist(Differentiated[2]))
      YYY <- cbind(unlist(FName[[SpectIndx]]@.Data[2]), unlist(Differentiated[2]))
      Xlim <- range(XXX)
      if (FName[[SpectIndx]]@Flags[1]==TRUE) {
         Xlim <- rev(Xlim)  ## reverse x-axis
      }
      Ylim <- range(YYY)
      matplot(x=XXX, y=YYY, type="l", lty="solid", col=c("black", "blue", "red"),
              xlim=Xlim, ylim=Ylim, xlab=FName[[SpectIndx]]@units[1], ylab=FName[[SpectIndx]]@units[2])
      return()
   }

   makeCombo <- function(){
      SelectedFName <- svalue(D1XpsSpect)
      FName <<- get(SelectedFName,envir=.GlobalEnv)  #carico in SampID il relativo XPSSAmple
      SpectList <<- XPSSpectList(SelectedFName)
      SpectIndx <<- 1
      delete(D1Frame2, D1CoreLine)
      D1CoreLine <<- gcombobox(SpectList, selected=-1, editable=FALSE, handler=function(h,...){
                        SpectName <- svalue(D1CoreLine)
                        SpectName <- unlist(strsplit(SpectName, "\\."))
                        SpectIndx <<- as.numeric(SpectName[1])
                        SpectName <<- SpectName[2]
                        assign("activeSpectName",SpectName,envir=.GlobalEnv)
                        assign("activeSpectIndx",SpectIndx,envir=.GlobalEnv)
                        plot(FName[[SpectIndx]])
                        enabled(DiffDeg) <- TRUE
                     }, container=D1Frame2)
      plot(FName)
   }
   
   SaveFName <- function(saveMode){
      SpectIndx <- get("activeSpectIndx",.GlobalEnv)
      Symbol <- FName[[SpectIndx]]@Symbol
      LL <- length(FName)
      DDeg <- as.numeric(svalue(DiffDeg))
      Nchr <- nchar(Symbol)
      chrPos <- gregexpr(pattern ="D.",Symbol)  #the result of gregexpr is a listcontaining the character position of D.
      chrPos <- chrPos[[1]][1]
#save data in NEW CORELINE mode
      if (chrPos > 0 && saveMode=="NewCL") {  #"D." present, operanting on previously differentiated data
         idx <- length(FName)+1                 #number of corelines +1
         FName[[idx]] <<- FName[[SpectIndx]]  #Add the new Coreline.
         FName[[idx]]@.Data <<- Differentiated
         if (length(FName[[idx]]@RegionToFit$y>0)){ #RegionToFit defined: remove all in the new differentiateed coreline
            FName[[idx]] <<- XPSremove(FName[[idx]],"all")
         }
         DDeg <- DDeg+as.integer(substr(Symbol,chrPos+2,chrPos+2))  #sum of the previous differentiation degree with the actual one
         substr(Symbol,chrPos+2,chrPos+2) <- as.character(DDeg)     #update the differentiation degree with the actual one
         DDeg <- as.integer(DDeg)
         Info <- FName[[idx]]@Info
         nI <- length(Info)+1                   #add new Info line
         SpectIndx <- idx
      }
      if (chrPos < 0 && saveMode=="NewCL") {  #"D." NOT present, operanting on NON differentiated data
         idx <- length(FName) +1                #number of corelines +1
         FName[[idx]] <<- FName[[SpectIndx]]  #Add the new Coreline.
         FName[[idx]]@.Data <<- Differentiated
         if (length(FName[[idx]]@RegionToFit$y>0)){ #RegionToFit defined: remove all in the new differentiateed coreline
            FName[[idx]] <<- XPSremove(FName[[idx]],"all")
         }
         Symbol <- paste("D.",DDeg,".",Symbol, sep="")  #compose the new coreline Symbol
         Info <- FName[[idx]]@Info
         nI <- length(Info)+1                   #add new Info line
         SpectIndx <- idx
      }
#save data in OVERWRITE mode
      if (chrPos > 0 && saveMode=="Overwrite") { #"D." is present in Symbol: data were previously differentiated
         FName[[SpectIndx]]@.Data <<- Differentiated  #change the orinal data with the filtered data
         DDeg <- DDeg+as.integer(substr(Symbol,chrPos+2,chrPos+2))  #sum of the previous differentiation degree with the actual one
         substr(Symbol,chrPos+2,chrPos+2) <- as.character(DDeg)     #update the differentiation degree with the actual one
         DDeg <- as.integer(DDeg)
         Info <- FName[[SpectIndx]]@Info
         nI <- length(Info)
      }
      if (chrPos < 0 && saveMode=="Overwrite") { #"D." NOT present, operanting on NON differentiated data
         idx <- length(FName)+1                 #number of corelines +1
         FName[[idx]] <<- FName[[SpectIndx]]  #Add the new Coreline.
         FName[[idx]]@.Data <<- Differentiated  #change the original data with the differentiated data
         Symbol <- paste("D.",DDeg,".",Symbol, sep="")  #compose the new coreline Symbol
         Info <- FName[[SpectIndx]]@Info
         nI <- length(Info)+1                   #add new Info line
         SpectIndx <- idx
      }
      FName[[SpectIndx]]@Symbol <<- Symbol    #update core line Symbol
      FName@names[SpectIndx] <<- Symbol              #update/add core line name to list of names
      Info[nI] <- paste("   :::Differentiated ", Symbol, ": Diff. degree: ", DDeg, sep="") #Update Differentiation information
      FName[[SpectIndx]]@Info <<- Info
      assign(ActiveFName, FName,envir=.GlobalEnv) #setto lo spettro attivo eguale ell'ultimo file caricato
      assign("activeSpectName", Symbol,.GlobalEnv) #set the activeSpectral name to the actual CoreLine name
      assign("activeSpectIndx", SpectIndx,.GlobalEnv) #set the activeSpectral name to the actual CoreLine name

      return(SpectIndx)
   }

#--- Variables
    plot.new()
    if (is.na(activeFName)){
       gmessage("No data present: please load and XPS Sample", title="XPS SAMPLES MISSING", icon="error")
       return()
    }
    ActiveFName <- get("activeFName", envir=.GlobalEnv)  #cload the XPSSample name (string)
    FName <- get(ActiveFName, envir=.GlobalEnv)  #load the active XPSSample (data)
    FNameList <- XPSFNameList()                  #list of loaded XPSSamples
    FNameIdx <- grep(ActiveFName,FNameList)
    SpectIndx <- get("activeSpectIndx", envir=.GlobalEnv) #load the index of the active CoreLine
    SpectName <- get("activeSpectName", envir=.GlobalEnv) #namer of the active CoreLine
    SpectList <- XPSSpectList(ActiveFName)       #list of CoreLine spectra belonging to FName XPSSample
    Differentiated <- NULL
    DiffBkp <- NULL
    BackGnd <- NULL


#===== Main Panel
    MainWin <- gwindow("DATA DIFFERENTIATION", visible=FALSE)
    size(MainWin) <- c(320, 200)
    mainFrame <- gframe(text="Differentiate", horizontal=FALSE, spacing=5, container=MainWin)

# --- Line 1 ---
    D1group1 <- ggroup(label="     CORELINE SELECTION     ", spacing=5, horizontal=FALSE, container=mainFrame)
    D1group2 <- ggroup(label="  ",horizontal=TRUE, container=D1group1)

    D1frame1 <- gframe(" Select XPS Sample ", spacing=5, container=D1group2)
    D1XpsSpect <- gcombobox(FNameList, selected=FNameIdx, editable=FALSE, handler=function(h,...){
                            makeCombo()
                 }, container=D1frame1)



    D1Frame2 <- gframe(" Select Coreline ", spacing=5, container=D1group2)
    D1CoreLine <- gcombobox(SpectList, selected=-1, editable=FALSE, handler=function(h,...){
                           XPSComponent <- svalue(D1CoreLine)
                           XPSComponent <- unlist(strsplit(XPSComponent, "\\."))   #skip the ".NUMBER" at beginning CoreLine name
                           SpectIndx <<- as.integer(XPSComponent[1])
                           SpectName <<- XPSComponent[2]
                           assign("activeSpectName", SpectName,.GlobalEnv) #set the activeSpectral name to the actual CoreLine name
                           assign("activeSpectIndx", SpectIndx,.GlobalEnv) #set the activeSpectIndex to the actual CoreLine
                           plot(FName[[SpectIndx]])
                           enabled(DiffDeg) <- TRUE
                 }, container = D1Frame2)

# --- Line 2 ---
    D2group1 <- ggroup(label="     DIFFERENTIATE     ", spacing=5, horizontal=FALSE, container=mainFrame)
    D2group2 <- ggroup(label="  ",horizontal=TRUE, container=D2group1)
    D2frame1 <- gframe(" Select Differentiation Degree ", spacing=5, container=D2group2)
    DiffDeg <- gcombobox(c(1,2,3,4,5), selected=-1, editable=FALSE, handler=function(h,...){
                           enabled(DiffButt) <- TRUE
                           enabled(ResButt) <- TRUE
                 }, container = D2frame1)
    DiffButt <- gbutton("            Differentiate            ", handler=function(h,...){
                           Ndiff <- as.numeric(svalue(DiffDeg))
                           SpectIndx <- get("activeSpectIndx", envir=.GlobalEnv) #set the activeSpectIndex to the actual CoreLine
                           Differentiated <<- FName[[SpectIndx]]@.Data
                           Object <- FName[[SpectIndx]]@.Data[[2]]
                           for(ii in 1:Ndiff){
                             Object <- Differ(Object)
                           }
                           Differentiated[[2]] <<- Object
                           DiffBkp <<- Differentiated
                           BkgSubtraction(FName[[SpectIndx]]@.Data[[2]])  #plotData needs the backgrnd to be defined
                           plotData()
                           enabled(SaveNewButt) <- TRUE
                           enabled(OverWButt) <- TRUE
                 }, container = D2frame1)

    ResButt <- gbutton(" RESET ", handler=function(h,...){
                           Differentiated <<- NULL
                           DiffBkp <<- NULL
                           BackGnd <<- NULL
                           Differentiated <<- NULL
                           svalue(D1CoreLine) <- ""
                           svalkue(DiffDeg) <- ""
                           makeCombo()
#                           plot(FName[[SpectIndx]])
                           plot.new()
                           enabled(DiffButt) <- FALSE
                           enabled(ResButt) <- FALSE
                           enabled(SaveNewButt) <- FALSE
                           enabled(OverWButt) <- FALSE
                 }, container = D2frame1)

    D2frame2 <- gframe(" Amplify Differentiated Data ", spacing=5, container=D2group1)
    AmpliDeg <- gcombobox(c(1,5,10,30,50,100), selected=-1, editable=FALSE, handler=function(h,...){
                           AA <- as.integer(svalue(AmpliDeg))
                           Differentiated[[2]] <<- AA*DiffBkp[[2]]
                           plotData()
                 }, container = D2frame2)


#--- Common buttons 

    DButtgroup <- ggroup(horizontal=TRUE, container=mainFrame)

    SaveNewButt <- gbutton("SAVE AS A NEW CORE LINE", handler=function(h,...){ #Filtered data saved in a new coreline
                    SaveFName("NewCL")
                    Differentiated <<- NULL
                    BackGnd <<- NULL
                    makeCombo()
                 }, container = DButtgroup)

    OverWButt <- gbutton("OVERWRITE PREVIOUS DIFFERENTIATION", handler=function(h,...){ #Filtered data saved in a new coreline
                    SaveFName("Overwrite")
                    Differentiated <<- NULL
                    BackGnd <<- NULL
                    makeCombo()
                    XPSSaveRetrieveBkp("save")
                 }, container = DButtgroup)



    gbutton("EXIT", handler=function(h,...){
                    dispose(MainWin)
                    XPSSaveRetrieveBkp("save")
                 }, container = DButtgroup)


    enabled(DiffDeg) <- FALSE
    enabled(DiffButt) <- FALSE
    enabled(SaveNewButt) <- FALSE
    enabled(OverWButt) <- FALSE
    visible(MainWin) <- TRUE
}
