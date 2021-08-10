#CoreLine processing: adding/deleting corelines or Baseline and Fits to exiting Corelines

#'Add a CoreLine (with fit if existing)to an XPS-Sample, or just BaseLine and Fit
#'or remove a Coreline from a XPS-Sample
#'
#'Add a Coreline from a source XPS-Sample into a destination XPS-Sample or
#'Add a just a BaseLine and Fit from a source XPS-Sample into a destination XPS-Sample or
#'remove a CoreLine from a XPS-Sample. No parameters are passed to this function
#'
#'
#'@examples
#'
#'\dontrun{
#'	ProcessCoreLine()
#'}
#'
#'@export
#'

XPSProcessCoreLine <- function(){

#---

   CtrlRepCL <- function(destIndx, SpectName, DestSpectList){    #CTRL for repeated CL: search for core-lines with same name
      winCL <- gwindow("SELECT CORELINE", visible=FALSE)
      size(winCL) <- c(100, 200)
      groupCL <- ggroup(horizontal=FALSE, container=winCL)
      DestFile <- svalue(DestFileName)
      N.CL <- length(destIndx)
      msg <- paste("Found", N.CL, SpectName, "spectra in", DestFile,"\nPlease select the coreline \nto add Baseline and Fit", sep=" ")
      txt <- glabel(text=msg, container=groupCL)
      font(txt) <- list(family="sans",size=12)
      gseparator(horizontal=TRUE, container=groupCL)
      selectCL <- gradio(DestSpectList[destIndx], selected=1, horizontal=TRUE, container=groupCL)
      gbutton("    OK     ", handler=function(h, ...){
                     CoreLine <- svalue(selectCL)
                     CoreLine <- unlist(strsplit(CoreLine, "\\."))   #drop "NUMBER." at beginning of coreLine name
                     RepCLidx <<- as.numeric(CoreLine[1])
#                     svalue(infoWin) <- ""
                     dispose(winCL)
                  }, container=groupCL)
      visible(winCL) <- TRUE
      winCL$set_modal(TRUE)  #nothing can be done while running this macro
   }                         #modal mode takes the control of the 'return' value which CANNOT be used to return variable values

   updateObj <- function(SFpointer, FRMpointer, CLpointer,...){
      SelectedFName <- svalue(SFpointer)
      SpectList <<- XPSSpectList(SelectedFName)
      SourceFName <- get(SelectedFName,envir=.GlobalEnv)  #load the source XPSSample file
      delete(FRMpointer,CLpointer)
      CLpointer <- gcombobox(SpectList, selected=-1, handler=function(h, ...){
                                  SourceCoreline <- svalue(CLpointer)
                                  SourceCoreline <- unlist(strsplit(SourceCoreline, "\\."))   #slipt strin at character "."
                                  indx <- as.integer(SourceCoreline[1])
                                  plot(SourceFName[[indx]])
                                  enabled(DestFileName) <- TRUE # disable selection of the destination file
                             }, editable=FALSE, container=FRMpointer)
      add(FRMpointer,CLpointer)
      plot(SourceFName)
      enabled(CLpointer) <- TRUE #enable core line selection
      return(CLpointer)

   }


   MakeBaseLine <- function(SourceFile, indx, DestFile, destIndx) {
        BLinfo <- SourceFile[[indx]]@Baseline$type
        BasLinType <- BLinfo[1]
        splinePoints <- NULL
        deg <- NULL
        Wgt <- NULL

        if (BasLinType == "shirley"){BasLinType <- "Shirley"}       #different names for old/new RXPSG packages
        if (BasLinType == "2p.shirley"){BasLinType <- "2P.Shirley"} #transform to new BaseLineNames.
        if (BasLinType == "3p.Shirley"){BasLinType <- "3P.Shirley"} #Exact Baseline Names required to generate the Baseline see XPSClass
        if (BasLinType == "lp.shirley"){BasLinType <- "LP.Shirley"}
        if (BasLinType == "2p.tougaard"){BasLinType <- "2P.Tougaard"}
        if (BasLinType == "3p.tougaard"){BasLinType <- "3P.Tougaard"}
        if (BasLinType == "4p.tougaard"){BasLinType <- "4P.Tougaard"}

        if (BasLinType == "linear" || BasLinType == "Shirley" || BasLinType == "2P.Shirley" || BasLinType == "2P.Tougaard" || BasLinType == "3P.Tougaard") {
           txt <- paste(BLinfo[1], " background found!\n  ==> Set the Baseline Limits")
           gmessage(msg=txt, title="HELP INFO", icon="info")
           plot(DestFile[[destIndx]])
           pos <- locator(n=2, type="p", col="red", lwd=2)
           DestFile[[destIndx]]@Boundaries$x <- pos$x
           DestFile[[destIndx]]@Boundaries$y <- pos$y
           DestFile[[destIndx]] <- XPSsetRegionToFit(DestFile[[destIndx]])
           DestFile[[destIndx]] <- XPSbaseline(DestFile[[destIndx]], BasLinType, deg, Wgt, splinePoints )
           DestFile[[destIndx]]@RSF <- SourceFile[[indx]]@RSF
        } else if (BasLinType == "polynomial") {
           deg <- as.numeric(BLinfo[2])
           gmessage(msg="Polynomial backgound found!\n ==> Set the Baseline Limits", title="HELP INFO", icon="info")
           plot(DestFile[[destIndx]])
           pos <- locator(n=2, type="p", col="red", lwd=2)
           DestFile[[destIndx]]@Boundaries$x <- pos$x
           DestFile[[destIndx]]@Boundaries$y <- pos$y
           DestFile[[destIndx]] <- XPSsetRegionToFit(DestFile[[destIndx]])
           DestFile[[destIndx]] <- XPSbaseline(DestFile[[destIndx]], BasLinType, deg, Wgt, splinePoints )
           DestFile[[destIndx]]@RSF <- SourceFile[[indx]]@RSF
        } else if (BasLinType == "spline") {
            splinePoints <- list(x=NULL, y=NULL)
            txt <- "Spline background found! \n ==> LEFT click to set spline points; RIGHT to exit"
            gmessage(msg=txt, title="HELP INFO", icon="info")
            plot(DestFile[[destIndx]])
            pos <- c(1,1) # only to enter in  the loop
            while (length(pos) > 0) {  #pos != NULL => mouse right button not pressed
                  pos <- locator(n=1, type="p", col=3, cex=1.5, lwd=2, pch=1)
                  if (length(pos) > 0) {
                      splinePoints$x <- c(splinePoints$x, pos$x)  # $x and $y must be separate to add new coord to splinePoints
                      splinePoints$y <- c(splinePoints$y, pos$y)
                  }
            }
            # Now make BaseLine
            decr <- FALSE #Kinetic energy set
            if (DestFile[[destIndx]]@Flags[1] == TRUE) { decr <- TRUE }
            idx <- order(splinePoints$x, decreasing=decr)
            splinePoints$x <- splinePoints$x[idx] #splinePoints$x in ascending order
            splinePoints$y <- splinePoints$y[idx] #following same order select the correspondent splinePoints$y
            LL <- length(splinePoints$x)

            DestFile[[destIndx]]@Boundaries$x <- c(splinePoints$x[1],splinePoints$x[LL]) #set the boundaries of the baseline
            DestFile[[destIndx]]@Boundaries$y <- c(splinePoints$y[1],splinePoints$y[LL])
            DestFile[[destIndx]] <- XPSsetRegionToFit(DestFile[[destIndx]])
            DestFile[[destIndx]] <- XPSbaseline(DestFile[[destIndx]], BasLinType, deg, Wgt, splinePoints )
            DestFile[[destIndx]]@RSF <- SourceFile[[indx]]@RSF
        } else if (BasLinType == "3P.Shirley" || BasLinType == "LP.Shirley" || BasLinType == "4P.Tougaard") {
            Wgt <- as.numeric(BLinfo[2])
            txt <- paste(BLinfo[1], " background found!\n  ==> Set the Baseline Limits")
            gmessage(msg=txt, title="HELP INFO", icon="info")
            plot(DestFile[[destIndx]])
            pos <- locator(n=2, type="p", col="red", lwd=2)
            DestFile[[destIndx]]@Boundaries$x <- pos$x
            DestFile[[destIndx]]@Boundaries$y <- pos$y
            DestFile[[destIndx]] <- XPSsetRegionToFit(DestFile[[destIndx]])
            DestFile[[destIndx]] <- XPSbaseline(DestFile[[destIndx]], BasLinType, deg, Wgt, splinePoints )
            DestFile[[destIndx]]@RSF <- SourceFile[[indx]]@RSF
        }

        plot(DestFile[[destIndx]])
        return(DestFile[[destIndx]])
   }



#---

   SaveSpectrum <- function(){
      activeFName <- svalue(DestFileName)
      assign(activeFName, DestFName, envir=.GlobalEnv) #save changes in the destinationfile
      assign("activeFName", activeFName, envir=.GlobalEnv) #set the Active XPSSample == DestinationFile
      assign("activeSpectName", activeSpectName,envir=.GlobalEnv) #Set the activeSpect == added coreline
      assign("activeSpectIndx", activeSpectIndx,envir=.GlobalEnv) #Set the active Index == index of added coreline
      plot(DestFName)
      XPSSaveRetrieveBkp("save")
      enabled(SaveSpect) <- FALSE    #SaveSpectrum file disabled
      enabled(SaveNewSpect) <- FALSE #Save Destination file disabled
      enabled(SaveAndExit) <- FALSE  #Saving data blocked: ctrls on Dest file needed
   }

   SaveNewSpectrum <- function(){
      activeFName <- svalue(SourceFile11) #name of the manipulated XPSSample (one of the math operations was performed)
      LL <- length(activeFName)
      if ( length(activeFName)>0 ){     #Math operations performed on SourceFile11
         SourceFName <- DestFName
         DestFName <- get(activeFName, envir=.GlobalEnv) #retrieve original source XPSSample which will be the destination file
         SpectName <- activeSpectName
         SpectIndx <- activeSpectIndx
         Symbol <- paste(prefix,SpectName, sep="")
         CoreLineList <- names(DestFName)
         destIndx <- length(CoreLineList)+1
         DestFName[[destIndx]] <- SourceFName[[SpectIndx]]  #this is the manipulated core line
         DestFName[[destIndx]]@Symbol <- Symbol
         DestFName@names <- c(CoreLineList,Symbol)
         assign(activeFName, DestFName, envir=.GlobalEnv) #Save the changes in a new coreline in the destination file
         assign("activeFName", activeFName, envir=.GlobalEnv)
         assign("activeSpectName", SpectName,envir=.GlobalEnv)
         assign("activeSpectIndx", SpectIndx,envir=.GlobalEnv)
         msg <- paste(" Coreline: ",  Symbol, " saved in the XPS Sample ", DestFName@Filename, sep="")
         cat("\n ==> ", msg)
#         svalue(infoWin) <- msg
         plot(DestFName)
         enabled(SaveNewSpect) <- FALSE #Save Destination file disabled if SourceFile OK
      } else {
         activeFName <- svalue(SourceFile1)
         SourceFName <- DestFName
         DestFName <- get(activeFName, envir=.GlobalEnv) #retrieve original source XPSSample which will be the destination file
         SpectName <- activeSpectName
         SpectIndx <- activeSpectIndx
         CoreLineList <- names(DestFName)
         destIndx <- length(CoreLineList)+1
         DestFName[[destIndx]] <- SourceFName[[SpectIndx]]  #this is the manipulated core line
         DestFName[[destIndx]]@Symbol <- SpectName
         DestFName@names <- c(CoreLineList,SpectName)
         assign(activeFName, DestFName, envir=.GlobalEnv) #Save the changes in a new coreline in the destination file
         assign("activeFName", activeFName, envir=.GlobalEnv)
         assign("activeSpectName", SpectName,envir=.GlobalEnv)
         assign("activeSpectIndx", SpectIndx,envir=.GlobalEnv)
         msg <- paste(" Coreline: ",  SpectName, " saved in the XPS Sample ", DestFName@Filename, sep="")
         cat("\n ==> ", msg)
#         svalue(infoWin) <- msg
         plot(DestFName)
      }
      enabled(SaveSpect) <- FALSE    #SaveSpectrum file disabled
      enabled(SaveNewSpect) <- FALSE #Save Destination file disabled
      enabled(SaveAndExit) <- FALSE  #Saving data blocked: ctrls on Dest file needed
   }


#----- variables -----

#---load list of file ID and correspondent FileNames
      FNameList <- XPSFNameList()
      DestFName <- NULL
      RepCLidx <- NULL
      activeSpectIndx <- NULL
      activeSpectName <- NULL
      SampID <- ""
      SpectList <- ""
      CullData <- NULL  #rangeX of the region to cull
      prefix <- ""

      gmessage(msg=" Remember to save data after each operation \n otherwise you will loss the results", title="SAVE RESULTS", icon="warning")

#####----main---
      AddWin <- gwindow("CORELINE PROCESSING", parent=c(100, 0), visible=FALSE)

      Addgroup <- ggroup(horizontal=FALSE, container=AddWin)
      NoteBK <- gnotebook(expand=TRUE, container = Addgroup)

#--- TAB1
      T1group <- ggroup(label="CORELINE PROCESSING", horizontal=TRUE, spacing=5, container=NoteBK)

      layout1T1 <- glayout(homogeneous=FALSE, spacing=3, container=T1group)


      layout1T1[1,1] <-  AddFrame1 <- gframe("SELECT THE SOURCE XPS-SAMPLE", spacing=3, container=layout1T1)
      SourceFile1 <- gcombobox(FNameList, selected=-1, editable=FALSE, expand=FALSE, handler=function(h,...){
#                                   svalue(infoWin) <- ""
                                   SourceCoreline1 <<- updateObj(SourceFile1,AddFrame2,SourceCoreline1)
                             }, container = AddFrame1)

      layout1T1[2,1] <- AddFrame2 <- gframe(text=" SELECT CORELINE TO PROCESS ", spacing=3, container=layout1T1)
      SourceCoreline1 <- gcombobox(SpectList, selected=-1, editable=FALSE, spacing=5, handler=function(h,...) {
                                   enabled(DestFileName) <- TRUE # selection destination file enabled
                             }, container=AddFrame2)

      enabled(SourceCoreline1) <- FALSE  #CoreLine selection blocked

      layout1T1[3,1] <- AddFrame3 <- gframe(text="SELECT THE DESTINATION FILE NAME", spacing=3, container=layout1T1)
      DestFileName <- gcombobox(FNameList, selected=-1, editable=FALSE, handler=function(h,...) {
                               SourceFile <- svalue(SourceFile1)
                               DestFile <- svalue(DestFileName)
#                               svalue(infoWin) <- ""
                               SpectName <- svalue(SourceCoreline1)
                               SpectName <- unlist(strsplit(SpectName, "\\."))   #split the string at character "."
                               SpectIndx <- as.integer(SpectName[1])
                               SpectName <- SpectName[2]
                               if (SourceFile == DestFile) {
                                  gmessage(msg="Warning: Destination File Name == Source File Name" , title = "WARNING: BAD DESTINATION FILE NAME!",  icon = "warning")
                               } else {
                                  SourceFile <- get(SourceFile, envir=.GlobalEnv)
                                  DestFile <- get(DestFile, envir=.GlobalEnv)
                                  destIndx <- grep(SpectName, names(DestFile))
                                  if (length(destIndx)==0){
                                     txt <- paste("No ",SpectName, " core line present in the destination file. Operation stopped!")
                                     gmessage(msg=txt , title = "WARNING: BAD DESTINATION FILE!",  icon = "warning")
                                     return()
                                  }
                                  LL <- length(destIndx)
                                  for(ii in 1:LL){
                                     if (DestFile[[destIndx[ii] ]]@Flags[1] != SourceFile[[SpectIndx]]@Flags[1]) { #acquisitions made using different energy scale
                                         gmessage(msg="XPS-Samples have different Energy units. Operation stopped!" , title = "WARNING: BAD DESTINATION FILE!",  icon = "warning")
                                         return()
                                     }
                                  }
                                  plot(DestFile)
                                  enabled(SaveSpect) <- TRUE # enable saving data if Dest File OK
                                  enabled(SaveAndExit) <- TRUE
                               }
                             }, container = AddFrame3)

      enabled(DestFileName) <- FALSE # selection of destination file blocked

      layout1T1[4,1] <- AddFrame4 <- gframe(text=" ADD WHOLE CORELINE ", horizontal=FALSE, spacing=3, container=layout1T1)
      Addbutton1 <- gbutton("Add New Coreline and Fit",  handler=function(h, ...){
                               SourceFile <- svalue(SourceFile1)
                               SourceCoreline <- svalue(SourceCoreline1)
#                               svalue(infoWin) <- ""
                               DestFile <- svalue(DestFileName)
                               if (length(SourceCoreline)==0) {   #No coreline selected
                                   gmessage(mag="Please select the source coreline please!" , title = "SOURCE CORELINE SELECTION",  icon = "warning")
                               } else if (length(DestFile)==0) {   #No destination file selected
                                   gmessage(mag="Please select the destination file please!" , title = "DESTINATION FILE SELECTION",  icon = "warning")
                               } else {
                                  SourceFile <- get(SourceFile, envir=.GlobalEnv)
                                  DestFName <<- get(DestFile, envir=.GlobalEnv)
                                  SourceCoreline <- unlist(strsplit(SourceCoreline, "\\."))   #split string at character "."
                                  SpectName <- SourceCoreline[2]
                                  SpectIndx <- as.integer(SourceCoreline[1])
                                  CoreLineList <- names(DestFName)
                                  destIndx <- length(DestFName)+1
                                  DestFName[[destIndx]] <<- SourceFile[[SpectIndx]]
                                  DestFName@names <<- c(CoreLineList,SourceFile[[SpectIndx]]@Symbol)
                                  activeSpectIndx <<- destIndx
                                  activeSpectName <<- SpectName
                                  msg <- paste(" Coreline: ",  SourceCoreline[2], "added to XPS Sample. PLEASE SAVE DATA!", DestFName@Filename)
                                  cat("\n ==> ", msg)
#                                  svalue(infoWin) <- msg
                                  plot(DestFName)
                                  enabled(SaveSpect) <- TRUE    #Enable SaveSpectrum
                                  enabled(SaveAndExit) <- TRUE
                              }
                           }, container=AddFrame4)

      layout1T1[6,1] <- AddFrame5 <- gframe(text=" ADD BASELINE AND FIT ", horizontal=FALSE, spacing=3, container=layout1T1)
      Addbutton2 <- gbutton("Add Baseline and Fit to Original Coreline",  handler=function(h, ...){
                               SourceFile <- svalue(SourceFile1)
                               SourceCoreline <- svalue(SourceCoreline1)
                               DestFile <- svalue(DestFileName)
#                               svalue(infoWin) <- ""
                               if (length(SourceCoreline)==0) {   #no coreline selected
                                  gmessage(msg="Please select the source coreline please!" , title = "SOURCE CORELINE SELECTION",  icon = "warning")
                               } else if (length(DestFile)==0) {  #no destination file selected
                                  gmessage(msg="Please select the destination file please!" , title = "DESTINATION FILE SELECTION",  icon = "warning")
                               } else {
                                  SourceFile <- get(SourceFile, envir=.GlobalEnv)
                                  DestSpectList <- XPSSpectList(DestFile)
                                  DestFName <<- get(DestFile, envir=.GlobalEnv)
                                  SourceCoreline <- unlist(strsplit(SourceCoreline, "\\."))
                                  SpectIndx <- as.integer(SourceCoreline[1])
                                  SpectName <- SourceCoreline[2]
                                  CoreLineList <- names(DestFName)
                                  destIndx <- grep(SpectName, CoreLineList)   #The index of the coreline in the destinationFile can be different from that od sourceFile => grep()

                                  if (length(destIndx) > 1){                #The same coreline can be present more than one time
                                     CtrlRepCL(destIndx, SpectName, DestSpectList)
                                     destIndx <- RepCLidx
                                  }

                                  if (length(destIndx) == 0) {
                                      text <- paste(SpectName, "not present in the Destination File: Fit copy stopped")
                                      gmessage(msg=text , title = "FIT COPY TO DESTINATION",  icon = "warning")
                                  } else {
                                      if (length(SourceFile[[SpectIndx]]@RegionToFit)==0){
                                          text <- paste("ATTENTION: NO fit found for ",SpectName, " of the Source XPS-Sample. Operation stopped!", sep="")
                                          gmessage(msg=text , title = "BAD SOURCE CORELINE SELECTION",  icon = "warning")
                                          return()
                                      }
                                      if (length(DestFName[[destIndx]]@Components)>0) {  #a fit is present for the selected coreline
                                          text <- paste("ATTENTION: constraints of ",SpectName, " fit will be kept! Continue?", sep="")
                                          answ <- gconfirm(msg=text , title = "WARNING!",  icon = "warning")
                                          if (answ) {
                                              DestFName[[destIndx]] <<- XPSremove(DestFName[[destIndx]],"all")
                                          } else {
                                              return()
                                          }
                                      }

                                      DestFName[[destIndx]] <<- MakeBaseLine(SourceFile, SpectIndx, DestFName, destIndx)
                                      DestFName[[destIndx]]@Components <<- SourceFile[[SpectIndx]]@Components
                                      RescaleH <- max(DestFName[[destIndx]]@RegionToFit$y)/max(SourceFile[[SpectIndx]]@RegionToFit$y) #scale factor between source and destination corelines
                                      LL=length(DestFName[[destIndx]]@Components)
                                      tmp <- NULL
                                      if (LL > 0) {
                                          for(ii in 1:LL) {
                                              varmu <- getParam(DestFName[[destIndx]]@Components[[ii]],variable="mu")
                                              DestFName[[destIndx]]@Components[[ii]] <<- setParam(DestFName[[destIndx]]@Components[[ii]], parameter=NULL, variable="mu", value=varmu)
                                              varh <- getParam(DestFName[[destIndx]]@Components[[ii]],variable="h")
                                              DestFName[[destIndx]]@Components[[ii]] <<- Ycomponent(DestFName[[destIndx]]@Components[[ii]], x=DestFName[[destIndx]]@RegionToFit$x, y=DestFName[[destIndx]]@Baseline$y/RescaleH) #Rescale Baseline Y values
                                              DestFName[[destIndx]]@Components[[ii]]@ycoor <<- RescaleH*DestFName[[destIndx]]@Components[[ii]]@ycoor  #Rescale Values of Y Components respect DestFName Y data
                                              varh$start <- max(DestFName[[destIndx]]@Components[[ii]]@ycoor-DestFName[[destIndx]]@Baseline$y)  #Substract Baseline
                                              varh$min <- 0
                                              varh$max <- 5*varh$start
                                              DestFName[[destIndx]]@Components[[ii]] <<- setParam(DestFName[[destIndx]]@Components[[ii]], parameter=NULL, variable="h", value=varh)
                                          }
                                          tmp <- sapply(DestFName[[destIndx]]@Components, function(z) matrix(data=z@ycoor))
                                          DestFName[[destIndx]]@Fit$y <<- (colSums(t(tmp)) - length(DestFName[[destIndx]]@Components)*(DestFName[[destIndx]]@Baseline$y))
                                      }

                                      plot(DestFName[[destIndx]])
                                      activeSpectIndx <<- destIndx
                                      activeSpectName <<- SpectName
                                      msg <- paste(" Fit data added to ", SpectName, "PLASE SAVE DATA!", sep="")
                                      cat("\n ==> ", msg)
#                                      svalue(infoWin) <<- msg
                                 }
                                 enabled(SaveSpect) <- TRUE    #Enable SaveSpectrum
                                 enabled(SaveAndExit) <- TRUE
                             }
                         }, container=AddFrame5)


      layout1T1[7,1] <- AddFrame6 <- gframe(text=" ADD FIT ONLY ", horizontal=FALSE, spacing=3, container=layout1T1)
      Addbutton2 <- gbutton("Add Fit",  handler=function(h, ...){
                               SourceFile <- svalue(SourceFile1)
                               SourceCoreline <- svalue(SourceCoreline1)
                               DestFile <- svalue(DestFileName)
                               svalue(warning) <- ""
                               if (length(SourceCoreline)==0) {   #No coreline has been selected in source File
                                  gmessage(msg="Please select the source coreline please!" , title = "SOURCE CORELINE SELECTION",  icon = "warning")
                               } else if (length(DestFile)==0) {  #No Destination File has been selected
                                  gmessage(msg="Please select the destination file please!" , title = "DESTINATION FILE SELECTION",  icon = "warning")
                               } else {
                                  SourceFile <- get(SourceFile, envir=.GlobalEnv)
                                  DestSpectList <- XPSSpectList(DestFile)
                                  DestFName <<- get(DestFile, envir=.GlobalEnv)
                                  SourceCoreline <- unlist(strsplit(SourceCoreline, "\\."))   #split string at char "."
                                  SpectName <- SourceCoreline[2]
                                  SpectIndx <- as.integer(SourceCoreline[1])
                                  CoreLineList <- names(DestFName)
                                  destIndx <- grep(SpectName, CoreLineList)  #The selected CoreLine name could be in any posiiton in the Destination XPSSample => source Samp Index could be different from Dest Samp index
                                  if (length(destIndx) == 0) {
                                      text <- paste(SpectName, "not present in the Destination File: Fit copy stopped")
                                      gmessage(msg=text , title = "FIT COPY TO DESTINATION STOPPED",  icon = "warning")
                                  } else {

                                     if (length(destIndx) > 1){                #The same coreline can be present more than one time
                                        CtrlRepCL(destIndx, SpectName, DestSpectList)
                                        destIndx <- RepCLidx
                                     }

                                     if (length(DestFName[[destIndx]]@Components) > 0) {  #Core line already fitted
                                          text <- paste("ATTENTION: fit present on ",SpectName, " New fit constaints will be kept! Continue?")
                                          answ <- gconfirm(msg=text , title = "WARNING!",  icon = "warning")
                                          if (answ) {
                                              DestFName[[destIndx]] <<- XPSremove(DestFName[[destIndx]],"all")
                                          } else {
                                              return()
                                          }
                                      }
                                      if (length(DestFName[[destIndx]]@Baseline) == 0) {  #baseline not present
                                         text <- paste("Baseline not present in the Destination File: Fit copy stopped")
                                         gmessage(msg=text , title = "FIT COPY TO DESTINATION STOPPED",  icon = "warning")
                                         return()
                                      }

                                      DestFName[[destIndx]] <<- XPSsetRSF(DestFName[[destIndx]])
                                      DestFName[[destIndx]]@Components <<- SourceFile[[SpectIndx]]@Components
                                      RescaleH <- max(DestFName[[destIndx]]@RegionToFit$y)/max(SourceFile[[SpectIndx]]@RegionToFit$y) #fattore di scala tra source coreline e destination coreline
                                      LL=length(DestFName[[destIndx]]@Components)

                                      for(ii in 1:LL) {
                                          varmu <- getParam(DestFName[[destIndx]]@Components[[ii]],variable="mu")
                                          DestFName[[destIndx]]@Components[[ii]] <<- setParam(DestFName[[destIndx]]@Components[[ii]], parameter=NULL, variable="mu", value=varmu)
                                          varh <- getParam(DestFName[[destIndx]]@Components[[ii]],variable="h")
                                          DestFName[[destIndx]]@Components[[ii]] <<- Ycomponent(DestFName[[destIndx]]@Components[[ii]], x=DestFName[[destIndx]]@RegionToFit$x, y=DestFName[[destIndx]]@Baseline$y/RescaleH) #calcola la Y ed aggiunge la baseline
                                          DestFName[[destIndx]]@Components[[ii]]@ycoor <<- RescaleH*DestFName[[destIndx]]@Components[[ii]]@ycoor
                                          varh$start <- max(DestFName[[destIndx]]@Components[[ii]]@ycoor-DestFName[[destIndx]]@Baseline$y)
                                          varh$min <- 0
                                          varh$max <- 5*varh$start
                                          DestFName[[destIndx]]@Components[[ii]] <<- setParam(DestFName[[destIndx]]@Components[[ii]], parameter=NULL, variable="h", value=varh)
                                      }
                                      tmp <- sapply(DestFName[[destIndx]]@Components, function(z) matrix(data=z@ycoor))
                                      DestFName[[destIndx]]@Fit$y <<- (colSums(t(tmp)) - length(DestFName[[destIndx]]@Components)*(DestFName[[destIndx]]@Baseline$y))
                                      plot(DestFName[[destIndx]])
                                      activeSpectIndx <<- destIndx
                                      activeSpectName <<- SpectName
                                      msg <- paste(" Fit data added to ", SpectName, "PLASE SAVE DATA!", sep="")
                                      cat("\n ==> ", msg)
#                                      svalue(infoWin) <- msg
                                      plot(DestFName)
                                      enabled(SaveSpect) <- TRUE    #Enable SaveSpectrum
                                      enabled(SaveAndExit) <- TRUE
                                 }
                             }
                         }, container=AddFrame6)


#      layout1T1[8,1] <- warn2 <- glabel(text="                                                " , container=layout1T1)

      layout1T1[9,1] <- AddFrame7 <- gframe(text="OVERWRITE A CORELINE ", horizontal=FALSE, spacing=3, container=layout1T1)
      overbutton <- gbutton("Overwrite",  handler=function(h, ...){
                               SourceFile <- svalue(SourceFile1)
                               SourceCoreline <- svalue(SourceCoreline1)
                               DestFile <- svalue(DestFileName)
                               if (length(SourceCoreline)==0) {   #No coreline selected
                                   gmessage("Please select the source coreline please!" , title = "SOURCE CORELINE SELECTION",  icon = "warning")
                               } else if (length(DestFile)==0) {   #No dest File selected
                                   gmessage("Please select the destination file please!" , title = "DESTINATION FILE SELECTION",  icon = "warning")
                               } else {
                                  SourceFile <- get(SourceFile, envir=.GlobalEnv)
                                  DestFName <<- get(DestFile, envir=.GlobalEnv)
                                  DestSpectList <- XPSSpectList(DestFile)
                                  SourceCoreline <- unlist(strsplit(SourceCoreline, "\\."))
                                  SpectName <- SourceCoreline[2]
                                  SpectIndx <- as.integer(SourceCoreline[1])
                                  CoreLineList <- names(DestFName)
                                  destIndx <- grep(SpectName, CoreLineList)  #The selected CoreLine name could be in any posiiton in the Destination XPSSample => source Samp Index could be different from Dest Samp index
                                  if (length(destIndx) > 1){                 #The same coreline can be present more than one time
                                     CtrlRepCL(destIndx, SpectName, DestSpectList)
                                     destIndx <- RepCLidx
                                  }


                                  DestFName[[destIndx]] <<- SourceFile[[SpectIndx]]
                                  plot(DestFName[[destIndx]])
                                  activeSpectIndx <<- destIndx
                                  activeSpectName <<- SpectName
                                  msg <- paste("Coreline ", SpectName, "overwritten. PLASE SAVE DATA!")
                                  cat("\n ==> ", msg)
#                                  svalue(infoWin) <- msg
                                  plot(DestFName)
                                  enabled(SaveSpect) <- TRUE    #Enable SaveSpectrum
                                  enabled(SaveAndExit) <- TRUE
                              }
                         }, container=AddFrame7)
#      layout1T1[10,1] <- warn3 <- glabel(text="                                                " , container=layout1T1)

      layout1T1[11,1] <- AddFrame8 <- gframe(text=" REMOVE/DUPLICATE CORELINE ", horizontal=TRUE, spacing=3, container=layout1T1)
      delbutton <- gbutton("         Remove          ",  handler=function(h, ...){
                               SourceFile <- svalue(SourceFile1)
                               SourceCoreline <- svalue(SourceCoreline1)
                               DestFile <- svalue(DestFileName)
                               SpectList <<- XPSSpectList(SourceFile)
                               if (length(SourceCoreline)==0) {   #No coreline selected
                                   gmessage("Please select the source coreline please!" , title = "SOURCE CORELINE SELECTION",  icon = "warning")
                               }
                               if (length(SourceFile)==0) {   #No Source File selected
                                   gmessage("Please select the Source File please!" , title = "SOURCE FILE SELECTION",  icon = "warning")
                               } else {
                                  enabled(SaveSpect) <- TRUE # enabling the Save options if the destination file name is OK
                                  enabled(SaveAndExit) <- TRUE
                                  svalue(DestFileName) <<- SourceFile
                                  SourceFile <- get(SourceFile, envir=.GlobalEnv)
                                  SourceCoreline <- unlist(strsplit(SourceCoreline, "\\."))
                                  SpectIndx <- as.integer(SourceCoreline[1])
                                  SpectName <- SourceCoreline[2]
                                  plot(SourceFile[[SpectIndx]])
                                  text <- paste("ATTENTION: are you sure to delete the",svalue(SourceCoreline1), "core line ?", sep="")
                                  answ <- gconfirm(msg=text, title="WARNING", icon="warning")
                                  if (answ) {
                                     SourceFile[[SpectIndx]] <- NULL #this eliminates the coreline
                                     DestFName <<- SourceFile      #move updated XPSSample in the destinatin file for saving
                                     activeSpectIndx <<- 1
                                     activeSpectName <<- names(SourceFile)[1]
                                     msg <- paste("Coreline ", SpectName, "deleted. PLEASE SAVE DATA!")
                                     cat("\n ==> ", msg)
#                                     svalue(infoWin) <- msg
                                     plot(SourceFile)
                                  } else {
                                     return()  #do nothing
                                  }

                                  SpectList <<- names(SourceFile) #update the SpectrumList
                                  idx <- seq(1:length(SpectList)) #cannot use XPSSpectList because loads the SourceFile from the .GlobalEnv
                                  SpectList <<- paste(idx, SpectList, sep=".")
                                  delete(AddFrame2,SourceCoreline1)
                                  SourceCoreline1 <<- gcombobox(SpectList, selected=-1, handler=function(h, ...){ #update the selection of corelines
                                                           SourceCoreline <- svalue(SourceCoreline1)
                                                           SourceCoreline <- unlist(strsplit(SourceCoreline, "\\."))
                                                           SpectIndx <- as.integer(SourceCoreline[1])
                                                           plot(SourceFile[[SpectIndx]])
                                                           enabled(DestFileName) <- TRUE #enable selection of the destination file
                                                          }, editable=FALSE, container=AddFrame2)
                                  add(AddFrame2,SourceCoreline1)
                                  enabled(SaveSpect) <- TRUE # Enable Save Spectrum
                                  enabled(SaveAndExit) <- TRUE
                              }
                         }, container=AddFrame8)

      duplibutton <- gbutton("         Duplicate        ",  handler=function(h, ...){
                               SourceFile <- svalue(SourceFile1)
                               SourceCoreline <- svalue(SourceCoreline1)
                               DestFile <- svalue(DestFileName)
                               if (length(SourceCoreline)==0) {   #No Coreline selected
                                   gmessage("Please select the source coreline please!" , title = "SOURCE CORELINE SELECTION",  icon = "warning")
                               }
                               if (length(DestFile)==0) {   #Destination file not selected
                                   DestFile <- SourceFile
                               }
                               if (SourceFile != DestFile) {#Source different from Destination file
                                   DestFile <- SourceFile
                               }
                               if (length(SourceFile)==0) { #Source File not selected
                                   gmessage("Please select the Source File please!" , title = "SOURCE FILE SELECTION",  icon = "warning")
                               } else {
                                  enabled(SaveSpect) <- TRUE #Save Destination file enabled if SourceFile OK
                                  enabled(SaveAndExit) <- TRUE
                                  svalue(DestFileName) <<- SourceFile
                                  DestFName <<- get(SourceFile, envir=.GlobalEnv)
                                  SourceCoreline <- unlist(strsplit(SourceCoreline, "\\."))   #Skip the CoreLine index
                                  SpectIndx <- as.integer(SourceCoreline[1])
                                  SpectName <- SourceCoreline[2]
                                  CoreLineList <- names(DestFName)
                                  destIndx <- length(DestFName)+1
                                  DestFName[[destIndx]] <<- DestFName[[SpectIndx]]
                                  DestFName[[destIndx]]@Symbol <<- SpectName
                                  DestFName@names <<- c(CoreLineList,SpectName)
                                  activeSpectIndx <<- destIndx
                                  activeSpectName <<- SpectName
                                  msg <- paste(" Coreline: ",  SpectName, "added to XPS Sample", DestFName@Filename, ". PLEASE SAVE DATA!", sep="")
                                  cat("\n ==> ", msg)
#                                  svalue(infoWin) <- msg
                                  plot(DestFName)
                                  enabled(SaveSpect) <- TRUE #Enable SaveSpectrum
                                  enabled(SaveAndExit) <- TRUE
                              }
                         }, container=AddFrame8)

#---

      layout2T1 <- glayout(homogeneous=FALSE, spacing=3, container=T1group)

      layout2T1[1,1] <- AddFrame9 <- gframe(text=" PICK UP DATA ", horizontal=FALSE, spacing=3, container=layout2T1)
      posgroup <- ggroup(horizontal=FALSE, container=AddFrame9)
      CullFrom <- gedit(initial.msg ="From?", handler=function(h, ...){
                               CullData[1] <<- as.numeric(svalue(CullFrom))
                         }, container=posgroup)

      CullTo <- gedit(initial.msg ="To?", handler=function(h, ...){
                               CullData[2] <<- as.numeric(svalue(CullTo))
                         }, container=posgroup)

      gbutton("MOUSE", handler=function(h, ...){
                               SourceFile <- svalue(SourceFile1)
                               SourceCoreline <- svalue(SourceCoreline1)
                               DestFile <- svalue(DestFileName)
                               DestFile <- SourceFile
                               svalue(CullFrom) <<- ""
                               svalue(CullTo) <<- ""
                               SourceCoreline <- unlist(strsplit(SourceCoreline, "\\."))   #tolgo il "NUMERO." all'inizio del nome coreline
                               SpectIndx <- as.integer(SourceCoreline[1])
                               SpectName <- SourceCoreline[2]
                               svalue(CullFrom) <- NULL
                               svalue(CullTo) <- NULL
                               CullData <- NULL
                               if (length(SourceFile)==0) {   #non e' stata selezionata alcuna coreline
                                   gmessage("Please select the Source File please!" , title = "SOURCE FILE SELECTION",  icon = "warning")
                                   return()
                               }
                               if (length(SourceCoreline)==0) {   #non e' stata selezionata alcuna coreline
                                   gmessage("Please select the source coreline please!" , title = "SOURCE CORELINE SELECTION",  icon = "warning")
                                   return()
                               } else {
                                  svalue(DestFileName) <<- SourceFile
                                  DestFName <<- get(SourceFile, envir=.GlobalEnv)
                                  SourceCoreline <- unlist(strsplit(SourceCoreline, "\\."))   #tolgo il "NUMERO." all'inizio del nome coreline
                                  SpectIndx <- as.integer(SourceCoreline[1])
                                  SpectName <- SourceCoreline[2]
                                  if (length(DestFName[[SpectIndx]]@RegionToFit)>0) {   #non e' stata selezionata alcuna coreline
                                     gmessage("Baseline or Fit present. Reset analysis before culling data!" , title = " RESET ANALYSIS ",  icon = "warning")
                                     return()
                                  }
                                  plot(DestFName[[SpectIndx]])
                                  pos <- locator(n=2, type="p", col="red", lwd=2, pch=3)
                                  CullData[1] <<- pos$x[1]
                                  CullData[2] <<- pos$x[2]
                               }
                         }, container=posgroup)

      cullgroup <- ggroup(horizontal=TRUE, container=AddFrame9)
      cullButton <- gbutton("Select Data", handler=function(h, ...){
                               SourceFile <- svalue(SourceFile1)
                               SourceCoreline <- svalue(SourceCoreline1)
                               svalue(DestFileName) <<- SourceFile
                               DestFile <- SourceFile
                               if (length(SourceFile)==0) {   #non e' stata selezionata alcuna coreline
                                  gmessage("Please select the Source File please!" , title = "SOURCE FILE SELECTION",  icon = "warning")
                                  return()
                               }
                               if (length(SourceCoreline)==0) {   #non e' stata selezionata alcuna coreline
                                  gmessage("Please select the source coreline please!" , title = "SOURCE CORELINE SELECTION",  icon = "warning")
                               } else {
                                  SourceCoreline <- unlist(strsplit(SourceCoreline, "\\."))   #tolgo il "NUMERO." all'inizio del nome coreline
                                  SpectIndx <- as.integer(SourceCoreline[1])
                                  SpectName <- SourceCoreline[2]
                                  DestFName <<- get(SourceFile, envir=.GlobalEnv)
                                  if (length(DestFName[[SpectIndx]]@RegionToFit)>0) {   #non e' stata selezionata alcuna coreline
                                      gmessage("Baseline or Fit present. Reset analysis before culling data!" , title = " RESET ANALYSIS ",  icon = "warning")
                                      return()
                                  }
                                  if (is.null(CullData[1]) || is.null(CullData[2])) {
                                      gmessage("From/To points NOT defined!" , title = " DEFINE FROM/TO ",  icon = "warning")
                                      return()
                                  }
                                  from <- findXIndex(DestFName[[SpectIndx]]@.Data[[1]], CullData[1])
                                  to <- findXIndex(DestFName[[SpectIndx]]@.Data[[1]], CullData[2])
                                  if (from > to) {
                                     tmp <- from
                                     from <- to
                                     to <- tmp
                                  }
                                  DestFName[[SpectIndx]]@.Data[[1]] <<- DestFName[[SpectIndx]]@.Data[[1]][from:to]
                                  DestFName[[SpectIndx]]@.Data[[2]] <<- DestFName[[SpectIndx]]@.Data[[2]][from:to]
                                  DestFName[[SpectIndx]]@.Data[[3]] <<- DestFName[[SpectIndx]]@.Data[[3]][from:to]
                                  activeSpectIndx <<- SpectIndx
                                  activeSpectName <<- SpectName
                                  plot(DestFName[[SpectIndx]])
                                  msg <- "Please SAVE Selected Data"
                                  cat("\n ==> ", msg)
#                                  svalue(infoWin) <- msg
                                  enabled(SaveNewSpect) <- TRUE
                                  enabled(SaveAndExit) <- TRUE
                               }
                         }, container=AddFrame9)

      resetButton <- gbutton("RESET", handler=function(h, ...){
                               SourceFile <- svalue(SourceFile1)
                               SourceCoreline <- svalue(SourceCoreline1)
                               svalue(DestFileName) <<- SourceFile
                               svalue(CullFrom) <<- ""
                               svalue(CullTo) <<- ""
                               DestFName <<- get(SourceFile, envir=.GlobalEnv)
                               SourceCoreline <- unlist(strsplit(SourceCoreline, "\\."))
                               SpectIndx <- as.integer(SourceCoreline[1])
                               plot(DestFName[[SpectIndx]])
                         }, container=AddFrame9)


###---TAB2

      T2group <- ggroup(label="CORELINE MATH", horizontal=FALSE, spacing=3, container=NoteBK)
      layoutT2 <- glayout(homogeneous=FALSE, spacing=3, container=T2group)


      layoutT2[1,1] <- MathFrame1 <- gframe("SELECT XPS-SAMPLE 1", spacing=3,container=layoutT2)
      SourceFile11 <- gcombobox(FNameList, selected=-1, editable=FALSE, handler=function(h,...){
#                                  svalue(infoWin) <- ""
                                  CoreLine11 <<- updateObj(SourceFile11,MathGroup2,CoreLine11)
                                  addHandlerChanged(CoreLine11, handler=function(h,...) {
                                              SourceFile <- svalue(SourceFile11)
                                              SourceFile <- get(SourceFile, envir=.GlobalEnv)
                                              Coreline <- svalue(CoreLine11)
                                              Coreline <- unlist(strsplit(Coreline, "\\."))   #Split string at point
                                              SpectIndx <- as.integer(Coreline[1])
                                              Range <- range(SourceFile[[SpectIndx]]@.Data[1])
                                              Range <- round(Range, 1)
                                              svalue(CLRange11) <<- paste("Range", "  ", Coreline[2], ": ", Range[1], ",  ", Range[2], sep="")
                                              enabled(DestFileName) <- TRUE # Enable the choice of the destination file
                                  })
                             }, container = MathFrame1)

      layoutT2[1,2] <- MathFrame11 <- gframe("SELECT XPS-SAMPLE 2", spacing=3,container=layoutT2)
      SourceFile22 <- gcombobox(FNameList, selected=-1, editable=FALSE, handler=function(h,...){
#                                  svalue(infoWin) <- ""
                                  if (length(svalue(CoreLine11))==0) {
                                     gmessage(msg="Please select the CORELINE1 first", title = "WARNING: CORELINE1 NOT DEFINED ",  icon = "warning")
                                     svalue(SourceFile22) <- NULL
                                  } else {
                                     CoreLine22 <<- updateObj(SourceFile22,MathGroup22,CoreLine22)
                                     addHandlerChanged(CoreLine22, handler=function(h,...) {
                                              SourceFile <- svalue(SourceFile22)
                                              SourceFile <- get(SourceFile, envir=.GlobalEnv)
                                              Coreline <- svalue(CoreLine22)
                                              Coreline <- unlist(strsplit(Coreline, "\\."))
                                              SpectIndx <- as.integer(Coreline[1])
                                              Range <- range(SourceFile[[SpectIndx]]@.Data[1])
                                              Range <- round(Range, 1)
                                              svalue(CLRange22) <<- paste("Range", "  ", Coreline[2], ": ", Range[1], ",  ", Range[2], sep="")

                                     })
                                  }
                             }, container = MathFrame11)

      layoutT2[2,1] <- MathFrame2 <- gframe(text="SELECT CORELINE 1", spacing=3, container=layoutT2)
      MathGroup2 <- ggroup(horizontal=FALSE, spacing=3, container=MathFrame2)
      CLRange11 <- glabel(text="Range CoreLine1: ", container=MathGroup2)
      CoreLine11 <- gcombobox(SpectList, selected=-1, editable=FALSE, container=MathGroup2)
      enabled(CoreLine11) <- FALSE

      layoutT2[2,2] <- MathFrame22 <- gframe(text="SELECT CORELINE 2", spacing=3, container=layoutT2)
      MathGroup22 <- ggroup(horizontal=FALSE, spacing=3, container=MathFrame22)
      CLRange22 <- glabel(text="Range CoreLine2: ", container=MathGroup22)
      CoreLine22 <- gcombobox(SpectList, selected=-1, editable=FALSE, container=MathGroup22)
      enabled(CoreLine22) <- FALSE

      layoutT2[3,1] <- MathGroup3 <- ggroup(horizontal=FALSE, container=layoutT2)
      MathFrame3 <- gframe("ADD A CONSTANT VALUE TO CORELINE1", horizontal=TRUE,spacing=3,container=MathGroup3)
      AddValue <- gedit(selected=-1, editable=TRUE, handler=function(h, ...){
#                                   svalue(infoWin) <- ""
                                   value <- as.numeric(svalue(AddValue))
                                   SourceFile <- svalue(DestFileName) <- svalue(SourceFile11)
                                   SourceCoreline <- svalue(CoreLine11)
                                   SourceFile <- get(SourceFile, envir=.GlobalEnv)
                                   SourceCoreline <- unlist(strsplit(SourceCoreline, "\\."))   #tolgo il "NUMERO." all'inizio del nome coreline
                                   SpectIndx <- as.integer(SourceCoreline[1])
                                   SpectName <- SourceCoreline[2]
                                   SourceFile[[SpectIndx]]@.Data[[2]] <- SourceFile[[SpectIndx]]@.Data[[2]]+value
                                   if (length(SourceFile[[SpectIndx]]@RegionToFit)>0)
                                         SourceFile[[SpectIndx]]@RegionToFit[[2]] <- SourceFile[[SpectIndx]]@RegionToFit[[2]]+value
                                   if (length(SourceFile[[SpectIndx]]@Baseline)>0)
                                         SourceFile[[SpectIndx]]@Baseline[[2]] <- SourceFile[[SpectIndx]]@Baseline[[2]]+value
                                   LL <- length(SourceFile[[SpectIndx]]@Components)
                                   if (LL>0) {
                                      for (jj in c(1:LL)) {
				                              SourceFile[[SpectIndx]]@Components[[jj]]@ycoor <- SourceFile[[SpectIndx]]@Components[[jj]]@ycoor+value
			                             }
                                   }
                                   plot(SourceFile[[SpectIndx]])
                                   DestFName <<- SourceFile
                                   activeSpectIndx <<- SpectIndx
                                   activeSpectName <<- SpectName
                                   prefix <<- "M."
                                   enabled(SaveSpect) <- TRUE
                                   enabled(SaveNewSpect) <- TRUE
                                   enabled(SaveAndExit) <- TRUE
                                   msg <- paste("Constant ",value, " added to CoreLine ", SpectName)
                                   cat("\n ==> ", msg)
#                                   svalue(infoWin) <- msg
                             }, container=MathFrame3)

      MathFrame33 <- gframe("MULTIPLY CORELINE1 BY A CONSTANT", horizontal=TRUE, spacing=3, container=MathGroup3)
      MultValue <- gedit(selected=-1, editable=TRUE, handler=function(h, ...){
#                                   svalue(infoWin) <- ""
                                   value <- as.numeric(svalue(MultValue))
                                   SourceFile <- svalue(DestFileName) <- svalue(SourceFile11)
                                   SourceCoreline <- svalue(CoreLine11)
                                   SourceFile <- get(SourceFile, envir=.GlobalEnv)
                                   SourceCoreline <- unlist(strsplit(SourceCoreline, "\\."))   #tolgo il "NUMERO." all'inizio del nome coreline
                                   SpectIndx <- as.integer(SourceCoreline[1])
                                   SpectName <- SourceCoreline[2]
                                   SourceFile[[SpectIndx]]@.Data[[2]] <- SourceFile[[SpectIndx]]@.Data[[2]]*value
                                   if (length(SourceFile[[SpectIndx]]@RegionToFit)>0)
                                         SourceFile[[SpectIndx]]@RegionToFit[[2]] <- SourceFile[[SpectIndx]]@RegionToFit[[2]]*value
                                   if (length(SourceFile[[SpectIndx]]@Baseline)>0)
                                         SourceFile[[SpectIndx]]@Baseline[[2]] <- SourceFile[[SpectIndx]]@Baseline[[2]]*value
                                   LL <- length(SourceFile[[SpectIndx]]@Components)
                                   if (LL > 0) {
                                      for (jj in c(1:LL)) {
				                              SourceFile[[SpectIndx]]@Components[[jj]]@ycoor <- SourceFile[[SpectIndx]]@Components[[jj]]@ycoor*value
			                             }
			                             SourceFile[[SpectIndx]]@Fit$y <- SourceFile[[SpectIndx]]@Fit$y*value
                                   }
                                   plot(SourceFile[[SpectIndx]])
                                   DestFName <<- SourceFile
                                   activeSpectIndx <<- SpectIndx
                                   activeSpectName <<- SpectName
                                   enabled(SaveSpect) <- TRUE
                                   enabled(SaveNewSpect) <- TRUE
                                   enabled(SaveAndExit) <- TRUE
                                   msg <- paste("XX CoreLine ", SpectName, " multiplied by ", value)
                                   cat("\n ==> ", msg)
#                                   svalue(infoWin) <- msg
                             }, container=MathFrame33)



      layoutT2[3,2] <- MathFrame5 <- gframe(" COMBINE CORELINE1 and CORELINE2 \n Coreline2 added to Coreline1", horizontal=FALSE, spacing=3, container=layoutT2)
      CombineButton <- gbutton("COMBINE", handler=function(h, ...){
#                                   svalue(infoWin) <- ""
                                   SourceFile1 <- svalue(SourceFile11)
                                   CoreLine1 <- svalue(CoreLine11)
                                   SourceFile1 <- get(SourceFile1, envir=.GlobalEnv)
                                   CoreLine1 <- unlist(strsplit(CoreLine1, "\\."))   #tolgo il "NUMERO." all'inizio del nome coreline
                                   SpectIndx1 <- as.integer(CoreLine1[1])
                                   SpectName <- CoreLine1[2]

                                   SourceFile2 <- svalue(SourceFile22)
                                   CoreLine2 <- svalue(CoreLine22)
                                   SourceFile2 <- get(SourceFile2, envir=.GlobalEnv)
                                   CoreLine2 <- unlist(strsplit(CoreLine2, "\\."))   #tolgo il "NUMERO." all'inizio del nome coreline
                                   SpectIndx2 <- as.integer(CoreLine2[1])
                                   if (length(SourceFile1[[SpectIndx1]]@RegionToFit)>0 || length(SourceFile2[[SpectIndx2]]@RegionToFit)){
                                      txt <- paste("Combine CoreLines  allowed ONLY on RAW data! RESET ANALYSIS")
                                      gmessage(msg=txt, title = "WARNING: FIT PRESENT ",  icon = "warning")
                                   } else { 
                                      if (SourceFile1[[SpectIndx1]]@Flags[1] == TRUE){
                                          SourceFile1[[SpectIndx1]]@.Data[[1]] <- rev(SourceFile1[[SpectIndx1]]@.Data[[1]])#scale in ascending order (KE)
                                          SourceFile1[[SpectIndx1]]@.Data[[2]] <- rev(SourceFile1[[SpectIndx1]]@.Data[[2]])#order also Y data
                                      }
                                      if (SourceFile2[[SpectIndx2]]@Flags[1] == TRUE){
                                          SourceFile2[[SpectIndx2]]@.Data[[1]] <- rev(SourceFile2[[SpectIndx2]]@.Data[[1]])
                                          SourceFile2[[SpectIndx2]]@.Data[[2]] <- rev(SourceFile2[[SpectIndx2]]@.Data[[2]])
                                      }
                                      Rng1 <- range(SourceFile1[[SpectIndx1]]@.Data[[1]])
                                      Rng2 <- range(SourceFile2[[SpectIndx2]]@.Data[[1]])
                                      MaxX1 <- max(Rng1)
                                      MaxX2 <- max(Rng2)
                                      MinX1 <- min(Rng1)
                                      MinX2 <- min(Rng2)
                                      LL1 <- length(SourceFile1[[SpectIndx1]]@.Data[[1]])
                                      LL2 <- length(SourceFile2[[SpectIndx2]]@.Data[[1]])
                                      if ( MaxX1 > MaxX2) {
                                         if (MaxX2 < MinX1) {
                                             gmessage("XPS-SAMPLE1 and XPS-SAMPLE2 are separated! \n
                                                       XPS-SAMPLE1 must overlap XPS-SAMPLE2 at least in ONE point", icon="warning")
                                             return()
                                         } else {
                                             Idx2 <- findXIndex(SourceFile2[[SpectIndx2]]@.Data[[1]], MinX1)
                                             DataToCombineX <- SourceFile2[[SpectIndx2]]@.Data[[1]][1:Idx2]
                                             DataToCombineY <- SourceFile2[[SpectIndx2]]@.Data[[2]][1:Idx2]
                                             yy1 <- SourceFile1[[SpectIndx1]]@.Data[[2]][1]    #average intensity of second block of data in the combination point
                                             yy2 <- SourceFile2[[SpectIndx2]]@.Data[[2]][Idx2] #average intensity of first block of data in the combination point
                                             Dy <- yy2-yy1 #difference of the intensities between block1, block2 od data in the combination point
                                             SourceFile1[[SpectIndx1]]@.Data[[1]] <- c(DataToCombineX, SourceFile1[[SpectIndx1]]@.Data[[1]])
                                             SourceFile1[[SpectIndx1]]@.Data[[2]] <- c(DataToCombineY-Dy, SourceFile1[[SpectIndx1]]@.Data[[2]])
                                             if (SourceFile1[[SpectIndx1]]@Flags[1] == TRUE){ #Restore original order
                                                SourceFile1[[SpectIndx1]]@.Data[[1]] <- rev(SourceFile1[[SpectIndx1]]@.Data[[1]])#scale in ascending order (KE)
                                                SourceFile1[[SpectIndx1]]@.Data[[2]] <- rev(SourceFile1[[SpectIndx1]]@.Data[[2]])#order also Y data
                                             }

                                             plot(SourceFile1[[SpectIndx1]])
                                             svalue(DestFileName) <- svalue(SourceFile11)
                                             DestFName <<- SourceFile1
                                             activeSpectIndx <<- SpectIndx1
                                             activeSpectName <<- SpectName
                                         }
                                      } else if ( MaxX2 > MaxX1) {
                                          if (MaxX1 < MinX2) {
                                             gmessage("XPS-SAMPLE1 and XPS-SAMPLE2 are separated! \n
                                                       XPS-SAMPLE1 must overlap XPS-SAMPLE2 at least in ONE point", icon="warning")
                                             return()
                                         } else {
                                             MaxX1 <- SourceFile1[[SpectIndx1]]@.Data[[1]][LL1]
                                             Idx2 <- findXIndex(SourceFile2[[SpectIndx2]]@.Data[[1]], MaxX1)+1 #+1 needed to not start with the data next to MaxX1
                                             DataToCombineX <- SourceFile2[[SpectIndx2]]@.Data[[1]][Idx2:LL2]
                                             DataToCombineY <- SourceFile2[[SpectIndx2]]@.Data[[2]][Idx2:LL2]
                                             yy1 <- SourceFile1[[SpectIndx1]]@.Data[[2]][LL1] #average intensity of first block of data in the combination point
                                             yy2 <- SourceFile2[[SpectIndx2]]@.Data[[2]][Idx2] #average intensity of second block of data in the combination point
                                             Dy <- yy2-yy1 #difference of the intensities between block1, block2 od data in the combination point
                                             SourceFile1[[SpectIndx1]]@.Data[[1]] <- c(SourceFile1[[SpectIndx1]]@.Data[[1]], DataToCombineX)
                                             SourceFile1[[SpectIndx1]]@.Data[[2]] <- c(SourceFile1[[SpectIndx1]]@.Data[[2]], DataToCombineY-Dy)
                                             if (SourceFile1[[SpectIndx1]]@Flags[1] == TRUE){ #Restore original order
                                                SourceFile1[[SpectIndx1]]@.Data[[1]] <- rev(SourceFile1[[SpectIndx1]]@.Data[[1]])#scale in ascending order (KE)
                                                SourceFile1[[SpectIndx1]]@.Data[[2]] <- rev(SourceFile1[[SpectIndx1]]@.Data[[2]])#order also Y data
                                             }

                                             plot(SourceFile1[[SpectIndx1]])
                                             svalue(DestFileName) <- svalue(SourceFile11)
                                             DestFName <<- SourceFile1
                                             activeSpectIndx <<- SpectIndx1
                                             activeSpectName <<- SpectName
                                         }
                                      }
                                   }

                                   enabled(SaveSpect) <- TRUE
                                   enabled(SaveNewSpect) <- TRUE
                                   enabled(SaveAndExit) <- TRUE
                                   prefix <<- "M."
                                   msg <- paste(svalue(CoreLine11), "-", svalue(SourceFile11)," combined with ", svalue(CoreLine22), "-", svalue(SourceFile22), " Combined")
                                   cat("\n ==> ", msg)
#                                   svalue(infoWin) <- msg
                             }, container=MathFrame5)


      layoutT2[4,1] <- MathFrame4 <- gframe("DIFFERENTIATE", horizontal=FALSE, spacing=3, container=layoutT2)
      DiffButt <- gbutton("DIFFERENTIATE CORELINE1", handler=function(h, ...){
#                                   svalue(infoWin) <- ""
                                   SourceFile <- svalue(DestFileName) <- svalue(SourceFile11)
                                   CoreLine1 <- svalue(CoreLine11)
                                   SourceFile <- get(SourceFile, envir=.GlobalEnv)
                                   CoreLine1 <- unlist(strsplit(CoreLine1, "\\."))   #tolgo il "NUMERO." all'inizio del nome coreline
                                   SpectIndx <- as.integer(CoreLine1[1])
                                   SpectName <- CoreLine1[2]

                                   LL <- length(SourceFile[[SpectIndx]]@.Data[[1]])
                                   SourceFile[[SpectIndx]]@RegionToFit <- list(x=NULL, y=NULL)
                                   SourceFile[[SpectIndx]]@Baseline <- list(x=NULL, y=NULL, type=NULL)
                                   SourceFile[[SpectIndx]]@Baseline$x <- SourceFile[[SpectIndx]]@RegionToFit$x <- SourceFile[[SpectIndx]]@.Data[[1]]
                                   SourceFile[[SpectIndx]]@Baseline$y <- rep(0, LL) #creo una baseline fittizia nulla
                                   SourceFile[[SpectIndx]]@Baseline$type <- "linear"
                                   k <- max(SourceFile[[SpectIndx]]@.Data[[2]])-min(SourceFile[[SpectIndx]]@.Data[[2]])
                                   for (ii in 2:LL){
                                       SourceFile[[SpectIndx]]@RegionToFit$y[ii] <- SourceFile[[SpectIndx]]@.Data[[2]][ii]-SourceFile[[SpectIndx]]@.Data[[2]][ii-1]
                                   }
                                   SourceFile[[SpectIndx]]@RegionToFit$y[1] <- SourceFile[[SpectIndx]]@RegionToFit$y[2]
                                   Dmin <- min(SourceFile[[SpectIndx]]@RegionToFit$y)
                                   Dmax <- max(SourceFile[[SpectIndx]]@RegionToFit$y)
                                   SourceFile[[SpectIndx]]@RegionToFit$y <- k*(SourceFile[[SpectIndx]]@RegionToFit$y-Dmin)/(Dmax-Dmin)
                                   matplot(x=matrix(c(SourceFile[[SpectIndx]]@.Data[[1]], SourceFile[[SpectIndx]]@RegionToFit$x), nrow=LL, ncol=2),
                                           y=matrix(c(SourceFile[[SpectIndx]]@.Data[[2]], SourceFile[[SpectIndx]]@RegionToFit$y), nrow=LL, ncol=2),
                                           type="l", lty=c(1,1), col=c("black","blue"),
                                           xlab=SourceFile[[SpectIndx]]@units[1], ylab=SourceFile[[SpectIndx]]@units[2])
                                   DestFName <<- SourceFile
                                   activeSpectIndx <<- SpectIndx
                                   activeSpectName <<- SpectName
                                   enabled(SaveNewSpect) <- TRUE
                                   enabled(SaveAndExit) <- TRUE
                                   prefix <<- "D1."
                                   msg <- paste(CoreLine1[2], " differentiated", sep="")
                                   cat("\n ==> ", msg)
#                                   svalue(infoWin) <- msg
                             }, container=MathFrame4)



      layoutT2[4,2] <- MathFrame44 <- gframe("ADD CORELINE2 TO CORELINE1", horizontal=FALSE, spacing=3, container=layoutT2)
      AddButt <- gbutton("ADD SPECTRA", handler=function(h, ...){
#                                   svalue(infoWin) <- ""
                                   SourceFile1 <- svalue(DestFileName) <- svalue(SourceFile11)
                                   SourceFile1 <- get(SourceFile1, envir=.GlobalEnv)
                                   CoreLine1 <- svalue(CoreLine11)
                                   CoreLine1 <- unlist(strsplit(CoreLine1, "\\."))   #split string at point
                                   SpectIndx1 <- as.integer(CoreLine1[1])             #keep first part
                                   SpectName <- CoreLine1[2]
                                   
                                   SourceFile2 <- svalue(SourceFile22)
                                   SourceFile2 <- get(SourceFile2, envir=.GlobalEnv)
                                   CoreLine2 <- svalue(CoreLine22)
                                   CoreLine2 <- unlist(strsplit(CoreLine2, "\\."))
                                   SpectIndx2 <- as.integer(CoreLine2[1])

                                   if (hasBaseline(SourceFile1[[SpectIndx1]]) || hasBaseline(SourceFile2[[SpectIndx2]]) ){
                                      answ <- gmessage(msg="ADDITION CAN BE DONE ONLY ON RAW DATA: PLEASE REMOVE ANALYSIS" , title = "WARNING!",  icon = "warning")
                                   } else {
                                      txt="          ==> ADDITION WILL BE PERFOMED IN THE COMMON ENERGY RANGE\nSAVE AND SAVE&EXIT WILL OVERWRITE THE RESULT TO CORE LINE 1"

                                      answ <- gconfirm(msg=txt , title = "WARNING!",  icon = "warning")
                                      if (answ) {
                                         Range1 <- range(SourceFile1[[SpectIndx1]]@.Data[1])
                                         Range2 <- range(SourceFile2[[SpectIndx2]]@.Data[1])
                                         lim1 <- max(Range1[1], Range2[1])  #the greater of the lower limits
                                         lim2 <- min(Range1[2], Range2[2])  #the smaller of the higher limits
                                         idx1 <- findXIndex(SourceFile1[[SpectIndx1]]@.Data[[1]], lim1)
                                         idx2 <- findXIndex(SourceFile1[[SpectIndx1]]@.Data[[1]], lim2)
                                         if (SourceFile1[[SpectIndx1]]@Flags[1]) {  #Binding energy scale
                                            CoreLineX1 <- SourceFile1[[SpectIndx1]]@.Data[[1]][idx2:idx1]
                                            CoreLineY1 <- SourceFile1[[SpectIndx1]]@.Data[[2]][idx2:idx1]
                                            CoreLineSF1 <- SourceFile1[[SpectIndx1]]@.Data[[3]][idx2:idx1]    #Analyzer Transf Funct.
                                         } else {                               #Kinetic energy scale
                                            CoreLineX1 <- SourceFile1[[SpectIndx1]]@.Data[[1]][idx1:idx2]
                                            CoreLineY1 <- SourceFile1[[SpectIndx1]]@.Data[[2]][idx1:idx2]
                                            CoreLineSF1 <- SourceFile1[[SpectIndx1]]@.Data[[3]][idx1:idx2]    #Analyzer Transf Funct.
                                         }
                                         LL1 <- length(CoreLineY1)
                                         idx1 <- findXIndex(SourceFile2[[SpectIndx2]]@.Data[[1]], lim1)
                                         idx2 <- findXIndex(SourceFile2[[SpectIndx2]]@.Data[[1]], lim2)
                                         if (SourceFile2[[SpectIndx2]]@Flags[1]) {  #Binding energy scale
                                            CoreLineY2 <- SourceFile2[[SpectIndx2]]@.Data[[2]][idx2:idx1]
                                         } else {                               #Kinetic energy scale
                                            CoreLineY2 <- SourceFile2[[SpectIndx2]]@.Data[[2]][idx1:idx2]
                                         }
                                         LL2 <- length(CoreLineY2)
                                         if(LL1 > LL2){  #verify CoreLine1, CoreLine2 have same length
                                           kk <- LL1-LL2   #if CoreLine1 is longer
                                           CoreLineY2[LL2:(LL2+kk)] <- CoreLineY2[LL2]  #add lacking values to CoreLine2
                                         } else if (LL1 < LL2){  #if CoreLine2 longer
                                           kk <- LL2-LL1
                                           CoreLineY2 <- CoreLineY2[1:(LL2-kk)]  #drop exceeding values
                                         }

                                         Nswps1 <- strsplit(SourceFile1[[SpectIndx1]]@Info[2], "Sweeps:")     #split str in 2 parts: before and after  "sweeps:"
                                         Nsw1 <- as.numeric(substr(Nswps1[[1]][2], 1, 2))                 #keep the first 2 characters of the second part = Nsweeps and transform in integer
                                         Nswps2 <- strsplit(SourceFile1[[SpectIndx2]]@Info[2], "Sweeps:")
                                         Nsw2 <- as.numeric(substr(Nswps2[[1]][2], 1, 2))
                                         SourceFile1[[SpectIndx1]]@.Data[[1]] <- CoreLineX1
                                         SourceFile1[[SpectIndx1]]@.Data[[2]] <- (CoreLineY1*Nsw1+CoreLineY2*Nsw2)/(Nsw1+Nsw2) #renormalize data
                                         SourceFile1[[SpectIndx1]]@.Data[[3]] <- CoreLineSF1

                                         Nsw <- as.character(Nsw1+Nsw2)
                                         Nswps1 <- paste(Nswps1[[1]][1], "Sweeps: ",Nsw,substr(Nswps1[[1]][2], 3, nchar(Nswps1[[1]][2])), sep="")     #compose the new second part of the string with correct number of sweeps
                                         SourceFile1[[SpectIndx1]]@Info[2] <- Nswps1
                                         plot(SourceFile1[[SpectIndx1]])


                                         DestFName <<- SourceFile1
                                         activeSpectIndx <<- SpectIndx1
                                         activeSpectName <<- SpectName
                                         enabled(SaveSpect) <- TRUE
                                         enabled(SaveNewSpect) <- TRUE
                                         enabled(SaveAndExit) <- TRUE
                                         prefix <<- "M."
                                         msg <- paste(CoreLine2[2],"-", svalue(SourceFile22), " added to ", CoreLine1[2], "-", svalue(SourceFile11))
                                         cat("\n ==> ", msg)
#                                         svalue(infoWin) <- msg
                                      }
                                   }
                             }, container=MathFrame44)

      layoutT2[5,1] <- MathFrame5 <- gframe("NORMALIZE", horizontal=FALSE, spacing=3, container=layoutT2)
      DiffButt <- gbutton("NORMALIZE CORELINE1", handler=function(h, ...){
#                                   svalue(infoWin) <- ""
                                   SourceFile <- svalue(DestFileName) <- svalue(SourceFile11)
                                   CoreLine1 <- svalue(CoreLine11)
                                   SourceFile <- get(SourceFile, envir=.GlobalEnv)
                                   CoreLine1 <- unlist(strsplit(CoreLine1, "\\."))   #tolgo il "NUMERO." all'inizio del nome coreline
                                   SpectIndx <- as.integer(CoreLine1[1])
                                   SpectName <- CoreLine1[2]
	                                  maxY <- max(SourceFile[[SpectIndx]]@.Data[[2]])
	                                  minY <- min(SourceFile[[SpectIndx]]@.Data[[2]])
	                                  SourceFile[[SpectIndx]]@.Data[[2]] <- (SourceFile[[SpectIndx]]@.Data[[2]]-minY)/(maxY-minY)
                                   if(length(SourceFile[[SpectIndx]]@RegionToFit) > 0){
                                      SourceFile[[SpectIndx]]@RegionToFit$y <- (SourceFile[[SpectIndx]]@RegionToFit$y-minY)/(maxY-minY)
                                   }
                                   if(length(SourceFile[[SpectIndx]]@Baseline) > 0){
                                      SourceFile[[SpectIndx]]@Baseline$y <- (SourceFile[[SpectIndx]]@Baseline$y-minY)/(maxY-minY)
                                   }
                                   if(LL <- length(SourceFile[[SpectIndx]]@Components) > 0){
                                      for(ii in 1:LL){
                                         SourceFile[[SpectIndx]]@Components[[ii]]$ycoor <- (SourceFile[[SpectIndx]]@Components[[ii]]$ycoor-minY)/(maxY-minY)
                                      }
                                   }
                                   plot(SourceFile[[SpectIndx]])
                                   DestFName <<- SourceFile
                                   activeSpectIndx <<- SpectIndx
                                   activeSpectName <<- SpectName
                                   enabled(SaveSpect) <- TRUE
                                   enabled(SaveNewSpect) <- TRUE
                                   enabled(SaveAndExit) <- TRUE
                                   prefix <<- ""
                                   msg <- paste(CoreLine1[2], " normalized", sep="")
                                   cat("\n ==> ", msg)
#                                   svalue(infoWin) <- msg
                             }, container=MathFrame5)


      layoutT2[6,1] <- MathFrame6 <- gframe("SUBTRACT THE BASELINE", horizontal=FALSE, spacing=3, container=layoutT2)
      glabel(">>> WARNING: Data OUTSIDE the \n    BaseLine Range will be lost!", container=MathFrame6)
      SubtrBlButt <- gbutton("SUBTRACT", handler=function(h, ...){
#                                  svalue(infoWin) <- ""
                                  if (gconfirm(msg="ATTENTION: original data will be IRREVERSIBLY modified!" , title = "BASELINE SUBTRACTION",  icon = "warning")) {
                                      SourceFile <- svalue(DestFileName) <- svalue(SourceFile11)
                                      SourceFile <- get(SourceFile, envir=.GlobalEnv)
                                      CoreLine1 <- svalue(CoreLine11)
                                      CoreLine1 <- unlist(strsplit(CoreLine1, "\\."))   #tolgo il "NUMERO." all'inizio del nome coreline
                                      SpectIndx <- as.integer(CoreLine1[1])
                                      SpectName <- CoreLine1[2]
                                      if (length(SourceFile[[SpectIndx]]@RegionToFit)>0 && length(SourceFile[[SpectIndx]]@Baseline)>0) {
                                            SourceFile[[SpectIndx]]@RegionToFit[[2]] <- SourceFile[[SpectIndx]]@RegionToFit[[2]]-SourceFile[[SpectIndx]]@Baseline[[2]]
                                      }
                                      SourceFile[[SpectIndx]]@.Data <- SourceFile[[SpectIndx]]@RegionToFit
                                      LL <- length(SourceFile[[SpectIndx]]@Components)
                                      if (LL>0) {
                                         for (ii in 1:LL) {
				                                 SourceFile[[SpectIndx]]@Components[[ii]]@ycoor <- SourceFile[[SpectIndx]]@Components[[ii]]@ycoor - SourceFile[[SpectIndx]]@Baseline[[2]]
			                                }
                                      }
                                      SourceFile[[SpectIndx]]@Baseline[[2]] <- SourceFile[[SpectIndx]]@Baseline[[2]]*0 #annullo la baseline

                                      plot(SourceFile[[SpectIndx]])
                                      DestFName <<- SourceFile
                                      activeSpectIndx <<- SpectIndx
                                      activeSpectName <<- SpectName
                                      enabled(SaveSpect) <- TRUE
                                      enabled(SaveNewSpect) <- TRUE
                                      enabled(SaveAndExit) <- TRUE
                                      msg <- paste("Baseline substracted")
                                      cat("\n ==> ", msg)
#                                      svalue(infoWin) <- msg
                                   }
                             }, container=MathFrame6)

      layoutT2[6,2] <- MathFrame66 <- gframe("SUBTRACT CORELINE2 FROM CORELINE1", horizontal=FALSE, spacing=3, container=layoutT2)
      SubtrSpectButt <- gbutton("SUBTRACT SPECTRA", handler=function(h, ...){
#                                   svalue(infoWin) <- ""
                                   SourceFile1 <- svalue(DestFileName) <- svalue(SourceFile11)
                                   SourceFile1 <- get(SourceFile1, envir=.GlobalEnv)
                                   CoreLine1 <- svalue(CoreLine11)
                                   CoreLine1 <- unlist(strsplit(CoreLine1, "\\."))   #tolgo il "NUMERO." all'inizio del nome coreline
                                   SpectIndx1 <- as.integer(CoreLine1[1])
                                   SpectName <- CoreLine1[2]

                                   SourceFile2 <- svalue(SourceFile22)
                                   SourceFile2 <- get(SourceFile2, envir=.GlobalEnv)
                                   CoreLine2 <- svalue(CoreLine22)
                                   CoreLine2 <- unlist(strsplit(CoreLine2, "\\."))   #tolgo il "NUMERO." all'inizio del nome coreline
                                   SpectIndx2 <- as.integer(CoreLine2[1])

                                   if (hasBaseline(SourceFile1[[SpectIndx1]]) || hasBaseline(SourceFile2[[SpectIndx2]]) ){
                                      answ <- gmessage(msg="SUBTRACTION CAN BE DONE ONLY ON RAW DATA: PLEASE REMOVE ANALYSIS" , title = "WARNING!",  icon = "warning")
                                   } else {
                                      txt="          ==> SUBTRACTION WILL BE PERFOMED IN THE COMMON ENERGY RANGE\nSAVE AND SAVE&EXIT WILL OVERWRITE SUBTRACTION TO THE ORIGINAL CORE LINE"
                                      answ <- gconfirm(msg=txt , title = "WARNING!",  icon = "warning")
                                      if (answ) {
                                          Range1 <- range(SourceFile1[[SpectIndx1]]@.Data[1])
                                          Range2 <- range(SourceFile2[[SpectIndx2]]@.Data[1])

                                          lim1 <- max(Range1[1], Range2[1])  #il piu' grande dei limiti inferiori
                                          lim2 <- min(Range1[2], Range2[2])  #il piu' piccolo dei limiti superiori
                                          idx1 <- findXIndex(SourceFile1[[SpectIndx1]]@.Data[[1]], lim1)
                                          idx2 <- findXIndex(SourceFile1[[SpectIndx1]]@.Data[[1]], lim2)
                                          if (SourceFile1[[SpectIndx1]]@Flags[1]) {  #Binding energy scale
                                             CoreLineX1 <- SourceFile1[[SpectIndx1]]@.Data[[1]][idx2:idx1]
                                             CoreLineY1 <- SourceFile1[[SpectIndx1]]@.Data[[2]][idx2:idx1]
                                             CoreLineSF1 <- SourceFile1[[SpectIndx1]]@.Data[[3]][idx2:idx1]    #Analyzer Transf Funct.
                                          } else {                               #Kinetic energy scale
                                             CoreLineX1 <- SourceFile1[[SpectIndx1]]@.Data[[1]][idx1:idx2]
                                             CoreLineY1 <- SourceFile1[[SpectIndx1]]@.Data[[2]][idx1:idx2]
                                             CoreLineSF1 <- SourceFile1[[SpectIndx1]]@.Data[[3]][idx1:idx2]    #Analyzer Transf Funct.
                                          }
                                          LL1 <- length(CoreLineY1)
                                          idx1 <- findXIndex(SourceFile2[[SpectIndx2]]@.Data[[1]], lim1)
                                          idx2 <- findXIndex(SourceFile2[[SpectIndx2]]@.Data[[1]], lim2)
                                          if (SourceFile2[[SpectIndx2]]@Flags[1]) {  #Binding energy scale
                                             CoreLineY2 <- SourceFile2[[SpectIndx2]]@.Data[[2]][idx2:idx1]
                                          } else {                               #Kinetic energy scale
                                             CoreLineY2 <- SourceFile2[[SpectIndx2]]@.Data[[2]][idx1:idx2]
                                          }
                                          LL2 <- length(CoreLineY2)
                                          if(LL1 > LL2){  #verify CoreLine1, CoreLine2 have same length
                                            kk <- LL1-LL2   #if CoreLine1 is longer
                                            CoreLineY2[LL2:(LL2+kk)] <- CoreLineY2[LL2]  #add lacking values to CoreLine2
                                          } else if (LL1 < LL2){  #if CoreLine2 longer
                                            kk <- LL2-LL1
                                            CoreLineY2 <- CoreLineY2[1:(LL2-kk)]  #drop exceeding values
                                          }

                                          SourceFile1[[SpectIndx1]]@.Data[[1]] <- CoreLineX1
                                          SourceFile1[[SpectIndx1]]@.Data[[3]] <- CoreLineSF1

                                          Nswps1 <- strsplit(SourceFile1[[SpectIndx1]]@Info[2], "Sweeps:")    #taglia la str in 2 parti: ql prima e ql dopo "Sweeps:"
                                          Nsw1 <- as.numeric(substr(Nswps1[[1]][2], 1, 2))                 #della seconda parte prendo i primi due caratteri che indicano Nswps e trasforno in intero
                                          Nswps2 <- strsplit(SourceFile1[[SpectIndx2]]@Info[2], "Sweeps:")
                                          Nsw2 <- as.numeric(substr(Nswps2[[1]][2], 1, 2))
                                          if (Nsw1 >= Nsw2) {
                                             SourceFile1[[SpectIndx1]]@.Data[[2]] <- CoreLineY1-CoreLineY2
                                             Nsw1 <- as.character(Nsw1)
                                             Nswps1 <- paste(Nswps1[[1]][1], "Sweeps: ",Nsw1,substr(Nswps1[[1]][2], 3, nchar(Nswps1[[1]][2])), sep="")     #ricostruisco la secopnda parte della stinga
                                             SourceFile1[[SpectIndx1]]@Info[2] <- Nswps1

                                             plot(SourceFile1[[SpectIndx1]])
                                             DestFName <<- SourceFile1
                                             activeSpectIndx <<- SpectIndx1
                                             activeSpectName <<- SpectName
                                             enabled(SaveSpect) <- TRUE
                                             enabled(SaveNewSpect) <- TRUE
                                             enabled(SaveAndExit) <- TRUE
                                             prefix <<- "M."
                                             msg <- paste(CoreLine2[2],"-", svalue(SourceFile22), " subtracted from ", CoreLine1[2], "-", svalue(SourceFile11))
                                             cat("\n ==> ", msg)
#                                             svalue(infoWin) <- msg
                                          } else {
                                             answ <- gmessage(msg="Nsweeps CoreLine2 > Nsweeps CoreLine1: SUBTRACTION NOT ALLOWED" , title = "WARNING!",  icon = "warning")
                                          }
                                      }
                                  }
                             }, container=MathFrame66)


#---InfoWin
#      InfoFrame <- gframe("PROCESSING INFO", horizontal=FALSE, spacing=3, container=Addgroup)
#      infoWin <- glabel(text="                                                " , container=InfoFrame)

#---COMMON BUTTONS
      ButtGroup1 <- ggroup(horizontal=TRUE, spacing=5, container=Addgroup)

      SaveSpect <- gbutton("       SAVE       ", handler=function(h, ...) { SaveSpectrum() }, container=ButtGroup1)
      enabled(SaveSpect) <- FALSE # Saving data blocked: ctrls on Dest file needed

      SaveNewSpect <- gbutton("  SAVE AS A NEW CORELINE  ", handler=function(h, ...){ 
                                  SaveNewSpectrum()
                         }, container=ButtGroup1)
      enabled(SaveNewSpect) <- FALSE # Saving data blocked: ctrls on Dest file needed

      SaveAndExit <- gbutton("     SAVE & EXIT     ", handler=function(h, ...){
                                  if (prefix=="D1." || length(CullData>0)){  #differentiation was performed
                                      CullData <<- NULL
                                      SaveNewSpectrum()   #for spectral differentiation or culled data saving is forced in a new coreline
                                  } else {
                                      SaveSpectrum()
                                  }
                                  dispose(AddWin)
                                  XPSSaveRetrieveBkp("save")
                         }, container=ButtGroup1)
      enabled(SaveAndExit) <- FALSE # Saving data blocked: ctrls on Dest file needed

      gbutton("       EXIT       ",  handler=function(h, ...) {
                                  dispose(AddWin)
                                  XPSSaveRetrieveBkp("save")
                         },container=ButtGroup1)
      
      visible(AddWin) <- TRUE
      svalue(NoteBK) <- 2
      svalue(NoteBK) <- 1     #set the first page

}

