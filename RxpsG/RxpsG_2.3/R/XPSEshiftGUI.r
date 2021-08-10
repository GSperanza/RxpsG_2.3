#'XPSEShift function
#'
#'Function to alligne the energy scale through a reference peak
#'As reference peaks C1s from hydrocarbons at BE=285eV or Au4f 7/2 at BE=84eV
#'are generally chosen.
#'No parameters are passed to this function
#'
#'@examples
#'
#'\dontrun{
#'	XPSEshift()
#'}
#' 
#'@export
#'

XPSEshift <- function(){


#--- Variables ---
   if (is.na(activeFName)){
       gmessage("No data present: please load and XPS Sample", title="XPS SAMPLES MISSING", icon="error")
       return()
   }
   FName <- get(activeFName,envir=.GlobalEnv)  #this is the XPS Sample
   ActiveFName <- get("activeFName",envir=.GlobalEnv)  #this is the name of the XPS Sample
   SpectIndx <- get("activeSpectIndx",envir=.GlobalEnv)
   SpectList <- XPSSpectList(ActiveFName)
   FNameList <- XPSFNameList()
   FNidx <- grep(ActiveFName, FNameList)
   FitCompList <- ""
   position <- NULL
   Eshift <- NULL


   Ewin <- gwindow("ENERGY SHIFT", visible=FALSE)
   Egroup1 <- ggroup(container = Ewin, horizontal=FALSE, label = "ENERGY SHIFT")

   Eframe0 <- gframe(" Select XPSSample ", horizontal=FALSE, spacing=5, container=Egroup1)
   Eobj0 <- gcombobox(FNameList, selected=FNidx, editable=FALSE, handler=function(h,...){
                                  ActiveFName <<- svalue(Eobj0)
                                  FName <<- get(ActiveFName,envir=.GlobalEnv)
                                  SpectList <<- XPSSpectList(ActiveFName)
                                  FNameList <<- XPSFNameList()
                                  FNidx <<- grep(ActiveFName, FNameList)
                                  svalue(Eobj2) <<- ""
                                  delete(Eframe1, Eobj1)
                                  Eobj1 <<- gcombobox(SpectList, selected=-1, editable=FALSE, handler=function(h, ...){
                                                    SpectName <- svalue(Eobj1)
                                                    SpectName <- unlist(strsplit(SpectName, "\\."))   #skip the number at beginning of string
                                                    SpectIndx <<- as.numeric(SpectName[1])
                                                    plot(FName[[SpectIndx]])
                                                    tmp <- names(FName[[SpectIndx]]@Components)
                                                    FitCompList <- c(FitCompList,tmp)
                                                    delete(Eframe2,Eobj2)
                                                    Eobj2 <<- gcombobox(FitCompList, selected=-1, editable=FALSE, index=TRUE, handler=function(h, ...){
                                                                       Component <- svalue(Eobj2)
                                                           	           FName <- get(ActiveFName, envir=.GlobalEnv)  #reload FName if a new EShift will be set
                                                                       OldValue <- FName[[SpectIndx]]@Components[[Component]]@param["mu","start"] # actual value of the fit component position
                                                                       OldValue <- round(as.numeric(OldValue), digits=2)
                                                                       svalue(Eobj5) <- OldValue   #display the actual value of the component position in GEdit
                                                                       }, container=Eframe2)
                                                    enabled(Eobj2) <- TRUE
                                                    }, container=Eframe1)
                                  }, container=Eframe0)


   Eframe1 <- gframe(" Select CoreLine ", horizontal=FALSE, spacing=5, container=Egroup1)
   Eobj1 <- gcombobox(SpectList, selected=-1, editable=FALSE, handler=function(h, ...){
                                  SpectName <- svalue(Eobj1)
                                  SpectName <- unlist(strsplit(SpectName, "\\."))
                                  SpectIndx <<- as.numeric(SpectName[1])
                                  plot(FName[[SpectIndx]])
                                  tmp <- names(FName[[SpectIndx]]@Components)
                                  FitCompList <- c(FitCompList,tmp)
                                  delete(Eframe2,Eobj2)
                                  Eobj2 <<- gcombobox(FitCompList, selected=-1, editable=FALSE, index=TRUE, handler=function(h, ...){
                                                     Component <- svalue(Eobj2)
                                       	           FName <- get(ActiveFName, envir=.GlobalEnv)  #reload FName if a new EShift will be set
                                                     OldValue <- FName[[SpectIndx]]@Components[[Component]]@param["mu","start"] # actual value of the fit component position
                                                     OldValue <- round(as.numeric(OldValue), digits=2)
                                                     svalue(Eobj5) <- OldValue   #display the actual value of the component position in GEdit
                                                  }, container=Eframe2)
                                  enabled(Eobj2) <- TRUE
                                  }, container=Eframe1)

   Eframe2 <- gframe(" Select Component ", horizontal=FALSE, spacing=5, container=Egroup1)
   Eobj2 <- gcombobox(FitCompList, selected=-1, editable=FALSE, index=TRUE, container=Eframe2)
   enabled(Eobj2) <- FALSE

   Eframe3 <- gframe("Apply Shift" , spacing=5, container=Egroup1)
   ShiftObj <- gradio(c("All Corelines", "Selected Coreline Only"), selected=1, horizontal=TRUE, container=Eframe3)

   Eframe4 <- gframe("Set Zoom Region" , horizontal=FALSE, spacing=5, container=Egroup1)
   glabel("Set the zooming region corners", container=Eframe4)
   gbutton("Zoom Corners",  handler=function(h, ...){
                                  SpectName <- svalue(Eobj1)
                                  if(length(SpectName)==0) {
                                     gmessage(msg="WARNING: no coreline selected, zoom stopped", title = "WARNING",icon = "warning" )
                                     return()
                                  }
                                  pos <- locator(n=2, type="p", pch=3, col="red", lwd=1.5) #first the two conrners are drawn
                                  if (FName[[SpectIndx]]@Flags[1]) { #Binding energy set
                                     pos$x <- c(max(pos$x), min(pos$x))
                                     pos$y <- c(min(pos$y), max(pos$y))
                                  } else {
                                     pos$x <- c(min(pos$x), max(pos$x))
                                     pos$y <- c(min(pos$y), max(pos$y))
                                  }
                                  plot(FName[[SpectIndx]], xlim=pos$x, ylim=pos$y) #zoom
                                  }, container = Eframe4)


   Eframe5 <- gframe("Peak Position", horizontal=FALSE, spacing=5, container=Egroup1)
   Eobj3 <- gbutton("Cursor",  handler=function(h, ...){
                                  SpectName <- svalue(Eobj1)
                                  if(length(SpectName)==0) {
                                     gmessage(msg="WARNING: no coreline selected", title = "WARNING",icon = "warning" )
                                     return()
                                  }
                                  pos <- locator(n=1, type="p", col="red", lwd=2)
                                  pos$x <- round(pos$x, digits=2)
                                  svalue(Eobj5) <- pos$x
                                  position <<- pos$x
                                  }, container=Eframe5)

   Eobj5 <- gedit(initial.msg ="E?", handler=function(h,...){
                        NewE <- svalue(Eobj5)
                        Component <- svalue(Eobj2)
                        if (is.null(position) && length(Component)==0){
                            gmessage(msg="Cursor-Peak-Position or Fit Component lacking! Please select", title="LACKING REFERENCE POSITION", icon="warning")
                            svalue(Eobj5) <- "E?"
                            return()
                        }
                        if (NewE == "" ){
                           gmessage("New energy position NULL. Please insert a correct value", title="ENERGY POSITION NULL", icon="warning")
                           svalue(Eobj5) <- "E?"
                           return()
                        }
                        if (NewE != "E?"){
                        NewE <- as.numeric(NewE)
                           Component <- svalue(Eobj2)
  	                        if (length(Component)==0){ #no fit present, no components, position read from cursor
    	                         Eshift <<- as.numeric(NewE-position)
  	                        } else {
                               Component <- as.numeric(substr(Component, 2,nchar(Component))) #exttract the component name
      	                      CompPos <- FName[[SpectIndx]]@Components[[Component]]@param["mu","start"]
    	                         Eshift <<- as.numeric(NewE-CompPos)
                           }
                           Escale <- FName[[SpectIndx]]@units[1]
                           if (svalue(ShiftObj)=="All Corelines") {
                              NCoreLines <- length(FName)
                              for (ii in 1:NCoreLines){
                                   if (Escale == FName[[ii]]@units[1]) { #Eshift calculated on a BE scale and CoreLine[[ii]] sa same Energy Units
    		                              FName[[ii]] <<- XPSapplyshift(FName[[ii]], Eshift)
                                   } else {
    		                              FName[[ii]] <<- XPSapplyshift(FName[[ii]], -Eshift)  # Eshift calculated on a BE scale while CoreLine[[ii]] is in Kinetic (or viceversa)
                                   }
                              }
                              plot(FName[[SpectIndx]])
                           } else if (svalue(ShiftObj)=="Selected Coreline Only"){ #Apply Eshift only on the selected coreline
                              if (Escale == FName[[SpectIndx]]@units[1]) {         #Eshift calculated on a BE scale and CoreLine[[ii]] sa same Energy Units
    		                         FName[[SpectIndx]] <<- XPSapplyshift(FName[[SpectIndx]], Eshift)
                               } else {
    		                         FName[[SpectIndx]] <<- XPSapplyshift(FName[[SpectIndx]],-Eshift)  # Eshift calculated on a BE scale while CoreLine[[ii]] is in Kinetic (or viceversa)
                               }
                               plot(FName[[SpectIndx]])
 		                     }
                        }
                        svalue(Eobj5) <- "E?"
                    }, container=Eframe5)

    gbutton("RESET", container = Egroup1, handler = function(h, ...){
		                FName <- XPSapplyshift(FName)
		                svalue(Eobj5) <- "E?"
                      if (length(FName[[SpectIndx]]@RegionToFit)>0) {
                         rangeX <- range(FName[[SpectIndx]]@RegionToFit$x)
                      } else {
                         rangeX <- range(FName[[SpectIndx]][1])
                      }
                      plot(FName[[SpectIndx]], xlim=rangeX)
    	                assign(ActiveFName, FName, envir=.GlobalEnv)
           })


   gbutton("SAVE", container=Egroup1, handler=function(...){
    	                assign(ActiveFName, FName, envir=.GlobalEnv)
    	                XPSSaveRetrieveBkp("save")
           })


   gbutton("SAVE & EXIT", container=Egroup1, handler=function(...){
    	                dispose(Ewin)
    	                assign(ActiveFName, FName, envir=.GlobalEnv)
    	                XPSSaveRetrieveBkp("save")
          })
   visible(Ewin) <- TRUE
}
