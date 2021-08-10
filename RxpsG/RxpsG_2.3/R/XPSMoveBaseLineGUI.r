#Macro to change Baseline extremes

#'Modifies the BaseLine level and limits of a given Coreline
#'
#'Function to modify the BaseLine ends and level for a given CoreLine.
#'No parameters are passed to this function.
#'
#'
#'@examples
#'
#'\dontrun{
#'	XPSMoveBaseLine()
#'}
#'
#'@export
#'


XPSMoveBaseLine <- function(){

   ReDraw <- function(){
       Xrange <- range(XPSSample[[indx]]@.Data[[1]])
       if (XPSSample[[indx]]@Flags[1]) {   #reverse if BE scale
           Xrange <- rev(Xrange)
       }

#--- Here the coreline and Baseline+Fit has to be displayed separately
       SampData <- as.matrix(XPSSample[[indx]]@.Data) #create spectrum data matrix for plot
       plot(x=SampData[[1]], y=SampData[[2]], xlim=Xrange, ylim=Yrange, type="l", lty="solid", lwd=1, col="black")
       SampData <- as(XPSSample[[indx]], "matrix") #create Baseline+Fit data matrix for plot
       NC <- ncol(SampData)
       if (NC > 2) { #there is a Baseline
          BaseLine <- SampData[,3]
          matlines(x=SampData[,1], y=BaseLine, xlim=Xrange, type="l", lty="solid", lwd=1, col="sienna")
       }
       if (NC > 3){ #there is a fit
          FitComp <- SampData[,4:NC-1]  #Only components and fit
          SpectFit <- SampData[,NC]  #fit
          matlines(x=SampData[,1], y=FitComp, xlim=Xrange, type="l", lty="solid", lwd=1, col="blue")
          matlines(x=SampData[,1], y=SpectFit, xlim=Xrange, type="l", lty="solid", lwd=1, col="red")
       }
       if (SetZoom == TRUE){   #set zoom area corners
           points(Corners, type="p", col=3, pch=3, cex=1.2, lwd=1)
           rect(Corners$x[1], Corners$y[1], Corners$x[4], Corners$y[4])
       }
   }



   updateObj <- function(h,...){
      SelectedXPSSample <- svalue(SourceFile)
      XPSSample <<- get(SelectedXPSSample,envir=.GlobalEnv)
      SpectList <<- XPSSpectList(SelectedXPSSample)
      delete(MBLFrame2,SourceCoreline)
      SourceCoreline <<- gcombobox(SpectList, selected=-1, handler=function(h, ...){
                                  SourceFile <- svalue(SourceFile)
                                  SourceCoreline <- svalue(SourceCoreline)
                                  SourceCoreline <- unlist(strsplit(SourceCoreline, "\\."))   #extract the spectrum idx
                                  indx <<- as.numeric(SourceCoreline[1])
                                  BLtype <- XPSSample[[indx]]@Baseline$type
                                  XPSSampleBkp <<- XPSSample #save for reset plot
                                  Yrange <<- range(XPSSample[[indx]]@.Data[[2]])
                                  ReDraw()
                             }, editable=FALSE, container=MBLFrame2)
      add(MBLFrame2,SourceCoreline)
      plot(XPSSample)
      enabled(SourceCoreline) <- TRUE #enable the selection of the coreline
   }


   do.baseline <- function(){
      deg <- NULL
      splinePoints <- NULL
      BLtype <- XPSSample[[indx]]@Baseline$type
      tmp <- NULL
      if ( indx != 0 && hasBoundaries(XPSSample[[indx]]) ) {
          XPSSample[[indx]] <<- XPSsetRegionToFit(XPSSample[[indx]])
          XPSSample[[indx]] <<- XPSbaseline(XPSSample[[indx]], BLtype, deg, splinePoints )
          LL <- length(XPSSample[[indx]]@Components)
          if (LL > 0) {
             for(ii in 1:LL){
                XPSSample[[indx]]@Components[[ii]] <<- Ycomponent(XPSSample[[indx]]@Components[[ii]], x=XPSSample[[indx]]@RegionToFit$x, y=XPSSample[[indx]]@Baseline$y) #calcola la Y eed aggiunge la baseline
             }
# update fit$y with sum of components
             tmp <- sapply(XPSSample[[indx]]@Components, function(z) matrix(data=z@ycoor))
             XPSSample[[indx]]@Fit$y <<- ( colSums(t(tmp)) - length(XPSSample[[indx]]@Components)*(XPSSample[[indx]]@Baseline$y))
          }
          plot(XPSSample[[indx]])
      }
   }


   MakeBaseLine <- function(){
        BLinfo <- XPSSample[[indx]]@Baseline$type
        BasLinType <- BLinfo[1]
        splinePoints <- NULL
        deg <- NULL
        Wgt <- NULL

#        BasLinType <- tolower(BasLinType)
        if (BasLinType == "shirley"){BasLinType <- "Shirley"}       #different names for old/new RXPSG packages
        if (BasLinType == "2P.shirley"){BasLinType <- "2P.Shirley"} #transform to default BaseLineNames of new RXPSG.
        if (BasLinType == "3P.Shirley"){BasLinType <- "3P.Shirley"}
        if (BasLinType == "LP.shirley"){BasLinType <- "LP.Shirley"}
        if (BasLinType == "2P.tougaard"){BasLinType <- "2P.Tougaard"}
        if (BasLinType == "3P.tougaard"){BasLinType <- "3P.Tougaard"}
        if (BasLinType == "4P.tougaard"){BasLinType <- "4P.Tougaard"}

        if (BasLinType == "linear" || BasLinType == "Shirley" || BasLinType == "2P.Shirley" || BasLinType == "2P.Tougaard" || BasLinType == "3P.Tougaard") {
           XPSSample[[indx]] <<- XPSsetRegionToFit(XPSSample[[indx]])
           XPSSample[[indx]] <<- XPSbaseline(XPSSample[[indx]], BLinfo[1], deg, Wgt, splinePoints )
           XPSSample[[indx]] <<- XPSsetRSF(XPSSample[[indx]])
        } else if (BasLinType == "polynomial") {
           deg <- as.numeric(BLinfo[2])
           XPSSample[[indx]] <<- XPSbaseline(XPSSample[[indx]], BLinfo[1], deg, Wgt, splinePoints )
        } else if (BasLinType == "spline") {
            splinePoints <- list(x=NULL, y=NULL)
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
            if (XPSSample[[indx]]@Flags[1] == TRUE) { decr <- TRUE }
            idx <- order(splinePoints$x, decreasing=decr)
            splinePoints$x <- splinePoints$x[idx] #splinePoints$x in ascending order
            splinePoints$y <- splinePoints$y[idx] #following same order select the correspondent splinePoints$y
            LL <- length(splinePoints$x)
            XPSSample[[indx]]@Boundaries$x <<- c(splinePoints$x[1],splinePoints$x[LL]) #set the boundaries of the baseline
            XPSSample[[indx]]@Boundaries$y <<- c(splinePoints$y[1],splinePoints$y[LL])
            XPSSample[[indx]] <<- XPSsetRegionToFit(XPSSample[[indx]])
            XPSSample[[indx]] <<- XPSbaseline(XPSSample[[indx]], BasLinType, deg, Wgt, splinePoints )
            XPSSample[[indx]] <<- XPSsetRSF(XPSSample[[indx]])
        } else if (BasLinType == "3P.Shirley" || BasLinType == "LP.Shirley" || BasLinType == "4P.Tougaard") {
            Wgt <- as.numeric(BLinfo[2])
            XPSSample[[indx]] <<- XPSsetRegionToFit(XPSSample[[indx]])
            XPSSample[[indx]] <<- XPSbaseline(XPSSample[[indx]], BasLinType, deg, Wgt, splinePoints )
            XPSSample[[indx]] <<- XPSsetRSF(XPSSample[[indx]])
        }
        LL <- length(XPSSample[[indx]]@Components)
        if (LL > 0) {
           for(ii in 1:LL){
              XPSSample[[indx]]@Components[[ii]] <<- Ycomponent(XPSSample[[indx]]@Components[[ii]], x=XPSSample[[indx]]@RegionToFit$x, y=XPSSample[[indx]]@Baseline$y) #calcola la Y eed aggiunge la baseline
           }
# update fit$y with sum of components
           tmp <- sapply(XPSSample[[indx]]@Components, function(z) matrix(data=z@ycoor))
           XPSSample[[indx]]@Fit$y <<- ( colSums(t(tmp)) - length(XPSSample[[indx]]@Components)*(XPSSample[[indx]]@Baseline$y))
        }
   }




#--- Variables ---

   rm(list = ls(envir = MyEnv, all.names = TRUE), envir = MyEnv)
   if (is.na(activeFName)){
       gmessage("No data present: please load and XPS Sample", title="XPS SAMPLES MISSING", icon="error")
       return()
   }

   XPSSample <- NULL
   XPSSampleBkp <- NULL
   XPSSampleList <- XPSFNameList()
   SpectList <- ""
   index <- NULL
   Yrange <- NULL
   Corners <- list(x=NULL, y=NULL)
   ExitWhile <- NULL
   SetZoom <- FALSE
   plot.new() #reset the graphic window

#--- GUI ---

   MBLwin <- gwindow("MOVE BASELINE", visible=FALSE)
   MBLgroup <- ggroup(label="", horizontal=FALSE, container=MBLwin)

   glabel("XPS SAMPLE", container=MBLgroup)
   MBLFrame1 <- gframe("Select the XPSSample data file", horizontal=TRUE,container=MBLgroup)
   SourceFile <- gcombobox(XPSSampleList, selected=-1, editable=FALSE, expand=FALSE, handler=updateObj, container = MBLFrame1)

   MBLFrame2 <- gframe(text="CORELINE", spacing=5, container=MBLgroup)
   SourceCoreline <- gcombobox(SpectList, selected=-1, editable=FALSE, container=MBLFrame2)

   MBLFrame3 <- gframe(text="Move Baseline", horizontal=FALSE, spacing=5, container=MBLgroup)
   MBLbutton <- gbutton(text="SET BASELINE BOUNDARIES", horizontal=FALSE, spacing=5, handler=function(h,...){
                    BLinfo <- XPSSample[[indx]]@Baseline$type
                    if (BLinfo[1] == "spline") {
                       txt <- "Spline background: \n ==> LEFT click to set spline points; RIGHT to exit"
                       gmessage(msg=txt, title="HELP INFO", icon="info")
                    } else {
                       txt <- paste(BLinfo[1], " background found!\n  ==> Set the Baseline Limits")
                       gmessage(msg=txt, title="HELP INFO", icon="info")
                       pos <- locator(n=2, type="p", col="red", lwd=2)
                       decr <- FALSE #Kinetic energy set
                       if (XPSSample[[indx]]@Flags[1] == TRUE) { decr <- TRUE } #Binding Energy set
                       idx <- order(pos$x, decreasing=decr)
                       pos$x <- pos$x[idx]
                       pos$y <- pos$y[idx]
                       XPSSample[[indx]]@Boundaries <<- pos
                    }
                    MakeBaseLine()
                    #if change the Baseline extension the fit components and the fit have to be recomputed
                    LL=length(XPSSample[[indx]]@Components)
                    tmp <- NULL
                    if (LL > 0) {
                       for(ii in 1:LL) {
                          XPSSample[[indx]]@Components[[ii]] <- Ycomponent(XPSSample[[indx]]@Components[[ii]], x=XPSSample[[indx]]@RegionToFit$x, y=XPSSample[[indx]]@Baseline$y) #computes the Y and add the Baseline
                       }
                       tmp <- sapply(XPSSample[[indx]]@Components, function(z) matrix(data=z@ycoor))
                       XPSSample[[indx]]@Fit$y <- (colSums(t(tmp)) - length(XPSSample[[indx]]@Components)*(XPSSample[[indx]]@Baseline$y))
                    }
                    plot(XPSSample[[indx]])
                 }, container=MBLFrame3)
#   MBLlabel <- glabel(text="Boundaries selection with Left mouse button ", container=MBLFrame3)

   gbutton("    SET ZOOM     ", handler = function(h, ...){
                    row1 <- " => Left click to define the 2 ZOOM REGION CORNERS (opposite in diagonal)"
                    row2 <- "\n => Click near corners to adjust Zoom Region Dimensions"
                    row3 <- "\n => Right click to Make Zoom when Zoom Region OK"
                    msg<-paste(row1, row2, row3, sep="")
                    gmessage( msg, icon="warning")
                    SetZoom <<- TRUE
                    pos <- locator(n=2, type="p", pch=3, cex=1, col=3, lwd=2)

                    Corners$x <<- c(pos$x[1],pos$x[1],pos$x[2],pos$x[2])
                    Corners$y <<- c(pos$y[1],pos$y[2],pos$y[1],pos$y[2])
                    Xrange <- range(XPSSample[[indx]]@.Data[[1]])
                    Yrange <<- range(XPSSample[[indx]]@.Data[[2]])
                    ReDraw()  #plots the zoom region
                    ExitWhile <- 1
                    while(ExitWhile > 0 ){
                          pos <- locator(n=1)
                          if (is.null(pos) ){
                              ExitWhile <- -1
                              SetZoom <<- FALSE
                              if (XPSSample[[indx]]@Flags[1]) { #Binding energy set
                                  Corners$x <- sort(Corners$x, decreasing=TRUE)
                              } else {                      #Kinetic energy set
                                  Corners$x <- sort(Corners$x, decreasing=FALSE)
                              }
                              Yrange <<- sort(c(Corners$y[1], Corners$y[2]), decreasing=FALSE)
                              ReDraw()
                              break()  #stops the while loop and exit
                          }
                          if (pos$x < Xrange[1]) {pos$x <- Xrange[1]}
                          if (pos$x > Xrange[2]) {pos$x <- Xrange[2]}
                          if (pos$y < Yrange[1]) {pos$y <- Yrange[1]}
                          if (pos$y > Yrange[2]) {pos$y <- Yrange[2]}

                          Dist <- NULL
                          Dmin <- ((pos$x-Corners$x[1])^2 + (pos$y-Corners$y[1])^2)^0.5  #valore di inizializzazione
                          for (ii in 1:4) {
                               Dist <- ((pos$x-Corners$x[ii])^2 + (pos$y-Corners$y[ii])^2)^0.5  #dist P0 P1
                               if(Dist <= Dmin){
                                  Dmin <- Dist
                                  idx <- ii
                               }
                          }
                          if (idx == 1){
                              Corners$x[1] <<- Corners$x[2] <<- pos$x
                              Corners$y[1] <<- Corners$y[3] <<- pos$y
                          } else if (idx == 2){
                              Corners$x[1] <<- Corners$x[2] <<- pos$x
                              Corners$y[2] <<- Corners$y[4] <<- pos$y
                          } else if (idx == 3){
                              Corners$x[3] <<- Corners$x[4] <<- pos$x
                              Corners$y[1] <<- Corners$y[3] <<- pos$y
                          } else if (idx == 4){
                              Corners$x[3] <<- Corners$x[4] <<- pos$x
                              Corners$y[2] <<- Corners$y[4] <<- pos$y
                          }
#                          if (Object@Flags[1]) { #Binding energy set
#                              point.coords$x <<- sort(c(Corners$x[1],Corners$x[3]), decreasing=TRUE) #pos$x in decreasing order
#                              point.coords$y <<- sort(c(Corners$y[1],Corners$y[4]), decreasing=FALSE)
#                          } else {
#                              point.coords$x <<- sort(c(Corners$x[1],Corners$x[3]), decreasing=FALSE) #pos$x in increasing order
#                              point.coords$y <<- sort(c(Corners$y[1],Corners$y[4]), decreasing=FALSE)
#                          }
                          ReDraw()
                    } #while loop end
                 }, container = MBLFrame3)


   gbutton("      RESET      ", handler=function(h,...){
                    Yrange <<- range(XPSSample[[indx]]@.Data[[2]])
                    SetZoom <<- FALSE

                    LL <- length(XPSSample[[indx]]@Components)
                    XPSSample[[indx]] <<- XPSSampleBkp[[indx]]
                    ReDraw()
                 }, container = MBLgroup)

   gbutton("       SAVE      ", handler=function(h,...){
    	              assign(activeFName, XPSSample, envir=.GlobalEnv)
                    XPSSaveRetrieveBkp("save")
                 }, container = MBLgroup)


   gbutton("   Save & EXIT   ", handler=function(h,...){
                    dispose(MBLwin)
     	              assign(activeFName, XPSSample, envir=.GlobalEnv)
                    XPSSaveRetrieveBkp("save")
                 }, container = MBLgroup)

   visible(MBLwin) <- TRUE

}
