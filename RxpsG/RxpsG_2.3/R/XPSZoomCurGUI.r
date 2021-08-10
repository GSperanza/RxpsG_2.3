## zoom function
#'To zoom a plotted spectrum
#'
#'To zoom a plotted spectrum.
#'No parameters are passed to this function
#'
#'
#'@examples
#'
#'\dontrun{
#'	XPSZoomCur()
#'}
#'
#'@export 
#'



XPSZoomCur <- function(){

ReDraw <- function(){   #redraw all spectrum with no restrictions to RegionToFit
   #XYrange <- get("ZM.XYrange", envir=MyEnv)
   SampData <- as.matrix(FName[[SpectIndx]]@.Data) #Whole coreline in SampData
   plot(x=SampData[[1]], y=SampData[[2]], xlim=XYrange$x, ylim=XYrange$y, type="l", lty="solid", lwd=1, col="black")
   SampData <- as(FName[[SpectIndx]], "matrix") # Regiontofit, Baseline, ecc in a matrix
   NC <- ncol(SampData)
   if (NC>2) { #a Baseline is present
       BaseLine <- SampData[,3]
       matlines(x=SampData[,1], y=BaseLine, xlim=XYrange$x, ylim=XYrange$y, type="l", lty="solid", lwd=1, col="sienna")
   }
   if (NC>3){ #c'e' un fit
      FitComp <- SampData[,4:NC-1]  #first three column skipped
      SpectFit <- SampData[,NC]  #fit
      matlines(x=SampData[,1], y=FitComp, xlim=XYrange$x, ylim=XYrange$y, type="l", lty="solid", lwd=1, col="blue")
      matlines(x=SampData[,1], y=SpectFit, xlim=XYrange$x, ylim=XYrange$y, type="l", lty="solid", lwd=1, col="red")
   }
}



FindNearest <- function(LocPos, Corners){

   D <- NULL
   Dmin <- ((LocPos$x-Corners$x[1])^2 + (LocPos$y-Corners$y[1])^2)^0.5  #init value
   for (ii in 1:4) {
       D[ii] <- ((LocPos$x-Corners$x[ii])^2 + (LocPos$y-Corners$y[ii])^2)^0.5  #dist P0 P1
       if(D[ii] <= Dmin){
          Dmin <- D[ii]
          idx=ii
       }
   }
   if (idx==1){
       Corners$x[1] <- Corners$x[2] <- LocPos$x
       Corners$y[1] <- Corners$y[3] <- LocPos$y
   } else if (idx==2){
       Corners$x[1] <- Corners$x[2] <- LocPos$x
       Corners$y[2] <- Corners$y[4] <- LocPos$y
   } else if (idx==3){
       Corners$x[3] <- Corners$x[4] <- LocPos$x
       Corners$y[1] <- Corners$y[3] <- LocPos$y
   } else if (idx==4){
       Corners$x[3] <- Corners$x[4] <- LocPos$x
       Corners$y[2] <- Corners$y[4] <- LocPos$y
   }
   Corners$x[idx] <- LocPos$x
   Corners$y[idx] <- LocPos$y
   return(Corners)
}



#----- Variabili -----
   coord <- list()
   LocPos <- list()
   Corners <- list(x=NULL, y=NULL)
   if (is.na(activeFName)){
       gmessage("No data present: please load and XPS Sample", title="XPS SAMPLES MISSING", icon="error")
       return()
   }
   FNameList <- XPSFNameList() #list of XPSSamples
   FName <- get(activeFName,envir=.GlobalEnv)
   SpectList <- XPSSpectList(activeFName)
   SpectIndx <- get("activeSpectIndx",envir=.GlobalEnv)
   activeSpectName <- get("activeSpectName",envir=.GlobalEnv)
   XYrange <- list(x=range(FName[[SpectIndx]]@.Data[1]), y=range(FName[[SpectIndx]]@.Data[2]))
   if (FName[[SpectIndx]]@Flags[1]) {   #reverse if BE scale
      XYrange$x <- rev(XYrange$x)
   }
   ResetXYrange <- XYrange

#===== Graphic Library Cntrl ===
   MatPlotMode <- get("MatPlotMode", envir=.GlobalEnv)
   if (MatPlotMode==FALSE){
      gmessage(msg="Overlay or CustomPlot active: DoubleClick on your XPSsample", title = "BAD GRAPHIC MODE",  icon = "error")
      return()
   }


#--- Redraw is used because plot is limited to the RegionToFit

   ReDraw()

#--- Zoom Cursor window ---

   ZMwin <- gwindow("ZOOM/CURSOR OPTION", parent=c(100,0) ,visible=FALSE)
   size(ZMwin) <- c(250, 270)

   ZMgroup <- ggroup(label="", horizontal=FALSE, container=ZMwin)

#   INFOframe <- gframe(text="INFO", horizontal=FALSE, spacing=1, container=ZMgroup)

   ZMframe0 <- gframe(text="XPS Sample and Core line Selection", horizontal=FALSE, spacing=5, container=ZMgroup)
   Fidx <- grep(activeFName, FNameList)
   XPS.Sample <- gcombobox(FNameList, selected=Fidx, editable=FALSE, handler=function(h,...){
                                 activeFName <- svalue(XPS.Sample)
                                 FName <<- get(activeFName, envir=.GlobalEnv)
                                 SpectList <<- XPSSpectList(activeFName)
                                 delete(ZMframe0, Core.Lines)
                                 Core.Lines <<- gcombobox(SpectList, selected=-1, editable=FALSE, handler=function(h,...){
                                                          SpectName <- svalue(Core.Lines)
                                                          SpectName <- unlist(strsplit(SpectName, "\\."))   #tolgo il N. all'inizio del nome coreline
                                                          SpectIndx <<- as.integer(SpectName[1])
                                                          XYrange <<- list(x=range(FName[[SpectIndx]]@.Data[1]), y=range(FName[[SpectIndx]]@.Data[2]))
                                                          if (length(FName[[SpectIndx]]@RegionToFit)>0) {   #reverse if BE scale
                                                              XYrange <<- list(x=range(FName[[SpectIndx]]@RegionToFit[1]), y=range(FName[[SpectIndx]]@RegionToFit[2]))
                                                          }
                                                          if (FName[[SpectIndx]]@Flags[1]) {   #reverse if BE scale
                                                              XYrange$x <<- rev(XYrange$x)
                                                          }
                                                          ReDraw()
                                                }, container=ZMframe0)
                       }, container = ZMframe0)
#   svalue(XPS.Sample) <- activeFName

   Core.Lines <- gcombobox(SpectList, selected=-1, editable=FALSE, handler=function(h,...){
                     SpectName <- svalue(Core.Lines)
                     SpectName <- unlist(strsplit(SpectName, "\\."))   #tolgo il N. all'inizio del nome coreline
                     SpectIndx <<- as.integer(SpectName[1])
                     XYrange <<- list(x=range(FName[[SpectIndx]]@.Data[1]), y=range(FName[[SpectIndx]]@.Data[2]))
                     if (length(FName[[SpectIndx]]@RegionToFit)>0) {   #reverse if BE scale
                        XYrange <<- list(x=range(FName[[SpectIndx]]@RegionToFit[1]), y=range(FName[[SpectIndx]]@RegionToFit[2]))
                     }
                     if (FName[[SpectIndx]]@Flags[1]) {   #reverse if BE scale
                        XYrange$x <<- rev(XYrange$x)
                     }
                     ReDraw()
                 }, container=ZMframe0)
   svalue(Core.Lines) <- paste(SpectIndx, ".",activeSpectName, sep="")

   ZMframe1 <- gframe(text="Set the zoom region", horizontal=FALSE, spacing=5, container=ZMgroup)

   ZMobj1 <- gbutton(" Set Zoom Limits ", handler=function(h,...){
                     txt <- "LEFT button to set the zoom area; RIGHT to exit \n Click near markers to modify the zoom area"
                     gmessage(msg=txt , title = "WARNING",  icon = "warning")
                     pos <- locator(n=2, type="p", pch=3, col="red", lwd=1.5) #first the two corners are drawn
                     rect(pos$x[1], min(pos$y), pos$x[2], max(pos$y))  #marker-Corners are ordered with ymin on Left and ymax on Right
                     Corners$x <- c(pos$x[1],pos$x[1],pos$x[2],pos$x[2])
                     Corners$y <- c(pos$y[1],pos$y[2],pos$y[1],pos$y[2])
                     points(Corners, type="p", pch=3, col="red", lwd=1.5)

                     LocPos <- list(x=0, y=0)
                     while (length(LocPos) > 0) {  #if pos1 not NULL a mouse butto was pressed
                        LocPos <- locator(n=1, type="p", pch=3, col="red", lwd=2) #to modify the zoom limits
                        if (length(LocPos$x) > 0) { #if the right mouse button NOT pressed
                           Corners <- FindNearest(LocPos,Corners)
                           if (FName[[SpectIndx]]@Flags[1]) { #Binding energy set
                              pos$x <- sort(c(Corners$x[1],Corners$x[3]), decreasing=TRUE) #pos$x in decrescent ordered => Corners$x[1]==Corners$x[2]
                           } else {
                              pos$x <- sort(c(Corners$x[1],Corners$x[3]), decreasing=FALSE) #pos$x in ascending order
                           }
                           pos$y <- c(Corners$y[1],Corners$y[2])
                           ReDraw()  #refresh graph
                           rect(pos$x[1], pos$y[1], pos$x[2], pos$y[2])
                           points(Corners, type="p", pch=3, col="red", lwd=1.5)
                        }
                     }
                     pos$y <- sort(c(Corners$y[1],Corners$y[2]), decreasing=FALSE) #pos$y in ascending order => Corners$y[1]==Corners$y[3]
                     XYrange <<- pos
                     plot(FName[[SpectIndx]], xlim=pos$x, ylim=pos$y) #zoom

                 }, container = ZMframe1)

   ZMobj2 <- gbutton("      RESET      ", handler=function(h,...){
                    XYrange <<- ResetXYrange
                    ReDraw()
                 }, container = ZMframe1)

   ZMframe2 <- gframe(text="Manual zoom values", horizontal=FALSE, spacing=5, container=ZMgroup)
   glabel(text="Exact Range Values:", spacing=1, container=ZMframe2)
   ZMgroup1 <- ggroup(horizontal=TRUE, container=ZMframe2)
   x1 <- gedit("", width=15, initial.msg = "Xmin= ", container=ZMgroup1)
   x2 <- gedit("", width=15, initial.msg = "Xmax= ", container=ZMgroup1)
   ZMgroup2 <- ggroup(horizontal=TRUE, container=ZMframe2)
   y1 <- gedit("", width=15, initial.msg = "Ymin= ", container=ZMgroup2)
   y2 <- gedit("", width=15, initial.msg = "Ymax= ", container=ZMgroup2)
   tkconfigure(x1$widget, width=18)
   tkconfigure(x2$widget, width=18)
   tkconfigure(y1$widget, width=18)
   tkconfigure(y2$widget, width=18)
   gbutton("  OK  ", handler=function(h,...){
                     x1 <- as.numeric(svalue(x1))
                     x2 <- as.numeric(svalue(x2))
                     y1 <- as.numeric(svalue(y1))
                     y2 <- as.numeric(svalue(y2))
                     if (FName[[SpectIndx]]@Flags) { #Binding energy set
                         XYrange$x <<- sort(c(x1, x2), decreasing=TRUE)
                         XYrange$y <<- sort(c(y1, y2))
                     } else {
                         XYrange$x <<- sort(c(x1, x2))
                         XYrange$y <<- sort(c(y1, y2))
                     }
                     ReDraw()
                 }, container = ZMframe2)


   ZMlabel <- glabel(text="    ", container=ZMgroup)
   ZMframe3 <- gframe(text="CURSOR POSITION", horizontal=FALSE, spacing=7, container=ZMgroup)

   ZMobj3 <- gbutton("Cursor Position", spacing=5, handler=function(h,...){
                     gmessage(msg="LEFT click to move marker's position; RIGHT to exit and zoom" , title = "WARNING",  icon = "warning")
                     pos <- c(1,1) # only to enter in  the loop
                     while (length(pos) > 0) {  #pos != NULL => mouse right button not pressed
                         pos <- locator(n=1, type="p", pch=3, cex=1.5, lwd=1.8, col="red")
                         if (length(pos) > 0) { #right mouse button not pressed
                            ReDraw() #refresh graph  to cancel previous cursor markers
                            points(pos, type="p", pch=3, cex=1.5, lwd=1.8, col="red")
                            pos <- round(x=c(pos$x, pos$y), digits=2)
                            txt <- paste("X: ", as.character(pos[1]), ", Y: ", as.character(pos[2]), sep="")
                            svalue(ZMobj4) <- txt
                            tcl("update", "idletasks")  #force writing cursor position in the glabel
                         }
                     }
                 }, container = ZMframe3)

   ZMobj4 <- glabel("Cursor position: ", container = ZMframe3)


   ZMobj6 <- gbutton("      CLOSE      ", handler=function(h,...){
                    rm(list = ls(envir = MyEnv, all.names = TRUE), envir = MyEnv)
                    dispose(ZMwin)
                 }, container = ZMgroup)

   visible(ZMwin) <- TRUE

}
