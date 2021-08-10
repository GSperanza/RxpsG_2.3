#-----------------------------------------
# XPS Sprucing with gWidgets2 and tcltk
#-----------------------------------------
#'Sprucing XPSSample data
#'
#'GUI to correct original XPS spectral data
#'
#'@return Returns the \code{Object} with corrected Survey spectrum.
#'
#'@examples
#'
#'\dontrun{
#'	XPSSprucingGUI()
#'}
#'
#'@export
#'

XPSSprucingGUI <- function() {

  import::here(tktoplevel, tkwm.geometry, tkgrid, tkconfigure,
                tkframe, tktitle, tkwinfo, tkbind, tkfocus, .from=tcltk )
  import::here(tkrplot, tkrreplot, .from=tkrplot )

  my.coords <- function(x, y){
      coords <- NULL
      #xx, yy are given in pixel. Now normalize
      xx <- as.numeric(x)/ww
      yy <- 1 - as.numeric(y)/hh
      #now bottom left corner is (0,0), top right corner is (1,1)

      XpltRng <- parplt[2]-parplt[1]    #this is the fraction of the window used for plot
      YpltRng <- parplt[4]-parplt[3]

      XusrRng <- parusr[2] - parusr[1]  #this is x box extension containing the spectrum
      YusrRng <- parusr[4] - parusr[3]

      #     Eliminate borders: shifts the window (0,0) in the plot (0,0) i.e. eliminates figure borders
      #           |        Normalize the fraction of the window used for plot i.e. plot bottom left corner is (0,0), top right corner is (1,1)
      #           |            |     Transforms the X, Y [0, 1] range in the user range
      #           |            |        |       Shifts the scale in the correct position adding the Xmin, Ymin
      #           |            |        |           |
      #      -----------   --------  -------    ---------
      xx <- (xx-parplt[1]) /XpltRng  *XusrRng + parusr[1]
      yy <- (yy-parplt[3]) /YpltRng  *YusrRng + parusr[3]

      Xlim1 <- min(c(parusr[1], parusr[2]))   #limits coordinates in the plot area
      Xlim2 <- max(c(parusr[1], parusr[2]))
      Ylim1 <- min(c(parusr[3], parusr[4]))
      Ylim2 <- max(c(parusr[3], parusr[4]))

      if (xx < Xlim1 ) {xx <- Xlim1}
      if (xx > Xlim2 ) {xx <- Xlim2}
      if (yy < Ylim1 ) {yy <- Ylim1}
      if (yy > Ylim2 ) {yy <- Ylim2}
      coords <<- c(xx, yy)
      return()
  }


  LBmousedown <- function(x, y) {
		    my.coords(x, y)
     	point.coords$x[point.index] <<- coords[1]   #abscissa
     	point.coords$y[point.index] <<- coords[2]   #ordinate
     	if (point.index==1) {
     	   point.index<<-2    #to modify the second edge of the selected area
         Corners$x<<-c(point.coords$x[1],point.coords$x[1],point.coords$x[2],point.coords$x[2])
         Corners$y<<-c(point.coords$y[1],point.coords$y[2],point.coords$y[1],point.coords$y[2])
  	   } else if (point.index==2) {
         Corners$x<<-c(point.coords$x[1],point.coords$x[1],point.coords$x[2],point.coords$x[2])
         Corners$y<<-c(point.coords$y[1],point.coords$y[2],point.coords$y[1],point.coords$y[2])
         point.index<<-3
  	   } else if (point.index==3) {
         D<-vector("numeric", 4)
         Dmin<-((point.coords$x[3]-Corners$x[1])^2 + (point.coords$y[3]-Corners$y[1])^2)^0.5  #valore di inizializzazione
         for (ii in 1:4) {
             D[ii]<-((point.coords$x[3]-Corners$x[ii])^2 + (point.coords$y[3]-Corners$y[ii])^2)^0.5  #dist P0 P1
             if(D[ii] <= Dmin){
                Dmin<-D[ii]
                idx=ii
             }
         }
         if (idx==1){
            Corners$x[1]<<-Corners$x[2]<<-point.coords$x[3]
            Corners$y[1]<<-Corners$y[3]<<-point.coords$y[3]
         } else if (idx==2){
             Corners$x[1]<<-Corners$x[2]<<-point.coords$x[3]
             Corners$y[2]<<-Corners$y[4]<<-point.coords$y[3]
         } else if (idx==3){
             Corners$x[3]<<-Corners$x[4]<<-point.coords$x[3]
             Corners$y[1]<<-Corners$y[3]<<-point.coords$y[3]
         } else if (idx==4){
            Corners$x[3]<<-Corners$x[4]<<-point.coords$x[3]
            Corners$y[2]<<-Corners$y[4]<<-point.coords$y[3]
         }
         if (Object[[coreline]]@Flags[1]) { #Binding energy set
            point.coords$x<<-sort(c(Corners$x[1],Corners$x[3]), decreasing=TRUE) #ordina pos$x in ordine decrescente Corners$x[1]==Corners$x[2]
            point.coords$y<<-sort(c(Corners$y[1],Corners$y[4]), decreasing=FALSE)
         } else {
            point.coords$x<<-sort(c(Corners$x[1],Corners$x[3]), decreasing=FALSE) #ordina pos$x in ordine crescente
            point.coords$y<<-sort(c(Corners$y[1],Corners$y[4]), decreasing=FALSE)
         }
      }
      replot()
  }

  undo.plot <- function(...){
      if (SelReg==1) {
         reset.boundaries()
     	   replot()
      } else if (SelReg>1) {
	      Object[[coreline]]@Boundaries$x <<- OldCoords$x
	      Ylimits <<- OldCoords$y
     	   replot()
  	   }
  }

  draw.plot <- function(...) {
      Xlimits <- Object[[coreline]]@Boundaries$x
	     if (point.index <= 2) {
	  	      plot(Object[[coreline]], xlim=Xlimits)
          points(point.coords, col="red", cex=1, lwd=1.5, pch=3)
 	    } else if (point.index>2){
	  	       plot(Object[[coreline]], xlim=Xlimits, ylim=Ylimits)
           points(Corners, type="p", col="red", cex=1, lwd=1.5, pch=3)
          rect(point.coords$x[1], point.coords$y[1], point.coords$x[2], point.coords$y[2])
	     }
      parusr <<- par("usr")
      parplt <<- par("plt")
  }


  replot <- function(...) { tkrreplot(Graph) }


  reset.boundaries <- function(h, ...) {
	     Object[[coreline]] <<- XPSremove(Object[[coreline]], "all")
      LL<-length(Object[[coreline]]@.Data[[1]])
      point.coords$x[1] <<- Object[[coreline]]@.Data[[1]][1] #ascissa primo estremo 1 del survey
      point.coords$y[1] <<- Object[[coreline]]@.Data[[2]][1] #ordinata primo estremo 1 del survey
      point.coords$x[2] <<- Object[[coreline]]@.Data[[1]][LL] #ascissa  secondo estremo del survey
      point.coords$y[2] <<- Object[[coreline]]@.Data[[2]][LL] #ordinata secondo estremo del survey
      slot(Object[[coreline]],"Boundaries") <<- point.coords
      Ylimits<<-c(min(Object[[coreline]]@.Data[[2]]), max(Object[[coreline]]@.Data[[2]]))
      OldCoords <<- point.coords #for undo
      Corners <- point.coords
      point.index <<- 1
      parusr <<- par("usr")
      parplt <<- par("plt")
      replot()
  }

  do.editRegion <- function(h, ...){
              idx1 <<- point.coords$x[1]
              idx2 <<- point.coords$x[2]
              newcoreline <- Object[[coreline]]
              idx1 <<- findXIndex(unlist(newcoreline@.Data[1]), point.coords$x[1]) #indice relativo al limite superiore (inferiore KE?) della RegionToFit
              idx2 <<- findXIndex(unlist(newcoreline@.Data[1]), point.coords$x[2]) #indice relativo al limite inferiore (superiore KE?) della RegionToFit
              tmp <- unlist(Object[[coreline]]@.Data[1]) #estraggo le ascisse regione selezionata
              newcoreline@.Data[[1]] <- tmp[idx1:idx2]     #aggiungo e ascisse regione selezionata
              tmp <- unlist(Object[[coreline]]@.Data[2]) #estraggo le ordinate regione selezionata
              newcoreline@.Data[[2]] <- tmp[idx1:idx2]     #aggiungo le ordinate regione selezionata
              DataTable <<- as.data.frame(cbind(newcoreline@.Data[[1]], newcoreline@.Data[[2]]))
              names(DataTable) <<- c("  X  ", "  Y  ")
              delete(EditGroup, DataToChange) #the pointer DataToChange is still alive. I add the pointer to a still unknown widget to EditGroup
              DataToChange <<- gdf(items=DataTable, container=EditGroup) #Now I define the widget associated to the pointer DataToChange
              size(DataToChange) <<- c(200,180)
              addHandlerChanged(DataToChange, handler=function(h,...){ #addHandlerChanged scarica il dataframe modificato in NewFirParam che e' salvato al di fuori di saveFitParam attraverso la <<-
                                     DataTable <<- h$obj[]
                               })
              enabled(ButtGroup) <- TRUE
  }


#====== VARIABLES DEFINITION=======

  Object <- get(activeFName,envir=.GlobalEnv)
  SpectList <- XPSSpectList(activeFName)
  coreline <- get("activeSpectIndx",envir=.GlobalEnv)
  if (length(indx <- grep("survey",SpectList))==0){
      indx <- grep("Survey",SpectList)
  }
  if (length(indx)==0){
      msg <- paste("Sorry Core Line ", coreLine, " NOT present in ", activeFName, " XPSSample", sep="")
      tkmessageBox(message = msg, icon = "warning", type = "ok")
  }
  XPSSettings <- get("XPSSettings",envir=.GlobalEnv)
  WinSize <- as.numeric(XPSSettings$General[4])
  DataTable <- NULL
  idx1 <- NULL
  idx2 <- NULL

  point.coords <- list(x=NA,y=NA)
  point.index <- 1
  coords <- NA # for printing mouse coordinates on the plot
  xx <- NULL
  yy <- NULL
#Coreline boundaries
  LL <- length(Object[[coreline]]@.Data[[1]])

  point.coords$x[1] <- Object[[coreline]]@.Data[[1]][1] #ascissa primo estremo 1 del survey
  point.coords$y[1] <- Object[[coreline]]@.Data[[2]][1] #ordinata primo estremo 1 del survey
  point.coords$x[2] <- Object[[coreline]]@.Data[[1]][LL] #ascissa  secondo estremo del survey
  point.coords$y[2] <- Object[[coreline]]@.Data[[2]][LL] #ordinata secondo estremo del survey
  Ylimits<-c(min(Object[[coreline]]@.Data[[2]]), max(Object[[coreline]]@.Data[[2]]))
  OldCoords <- point.coords #for undo
  Corners <- point.coords
  SelReg <- 0

  Object[[coreline]]@Boundaries$x <- c(point.coords$x)
  Object[[coreline]]@Boundaries$y <- c(point.coords$y)
  parusr <- NA
  parplt <- NA

  DataToChange <- NULL
#====== Widget Definition =======
  SPwin <- gwindow("XPS SPRUCING GUI", parent=c(100,-10), toolkit = "tcltk", visible = FALSE)
  addHandlerDestroy(SPwin, handler=function(h, ...){  #if MainWindow unproperly closed with X
                               tkdestroy(TKwinGraph)    #closes the TK window
                               X11()                    #restore the X11() graphic window
                               plot.new()
                               plot(Object) #replot the CoreLine
                           })
  size(SPwin) <- c(250, 600)
  MainGroup <- ggroup(horizontal = FALSE, spacing=3, container = SPwin)



#  SPlayout <- glayout(homogeneous=FALSE, spacing=3, container=MainGroup)

#====== OPTIONS SECTION =======
#  SPlayout[1, 1] <- OptGroup <- ggroup(horizontal = FALSE, spacing = 3, container = SPlayout)
#  OptGroup <- ggroup(horizontal = FALSE, spacing = 3, container = MainGroup)

  gframe10 <- gframe(text="SELECT CORELINE", expand=TRUE, spacing=3, container=MainGroup)
  CLobj10 <- gcombobox(SpectList, selected=coreline, editable=FALSE, handler=function(h,...){
              XPSCL <- svalue(CLobj10)
              XPSCL <- unlist(strsplit(XPSCL, "\\."))   #remove number at CoreLine beginning
              coreline <<- as.numeric(XPSCL[1])
              #Coreline boundaries
              LL <- length(Object[[coreline]]@.Data[[1]])
              point.coords$x[1] <<- Object[[coreline]]@.Data[[1]][1] #ascissa primo estremo 1 del survey
              point.coords$y[1] <<- Object[[coreline]]@.Data[[2]][1] #ordinata primo estremo 1 del survey
              point.coords$x[2] <<- Object[[coreline]]@.Data[[1]][LL] #ascissa  secondo estremo del survey
              point.coords$y[2] <<- Object[[coreline]]@.Data[[2]][LL] #ordinata secondo estremo del survey
              Ylimits <<- c(min(Object[[coreline]]@.Data[[2]]), max(Object[[coreline]]@.Data[[2]]))
              OldCoords <<- point.coords #for undo
              Corners <<- point.coords
              SelReg <<- 0
              Object[[coreline]]@Boundaries$x <<- c(point.coords$x)
              Object[[coreline]]@Boundaries$y <<- c(point.coords$y)
              parusr <<- NA
              parplt <<- NA
              replot()
         }, container=gframe10)

  gframe22 <- gframe(text = " SELECT DATA ", spacing=3, horizontal = FALSE, container = MainGroup)

  gbutton(" SELECT REGION ", container = gframe22, handler = function(h, ...){
              OldCoords <<- Object[[coreline]]@Boundaries
              SelReg <<- SelReg+1
              rngX <- range(point.coords$x)
              rngX <- (rngX[2]-rngX[1])/20
              rngY <- range(point.coords$y)
              rngY <- (rngY[2]-rngY[1])/20

              if (Object[[coreline]]@Flags[1]) { #Binding energy set
                 point.coords$x<-sort(point.coords$x, decreasing=TRUE) #ordina pos$x in ordine decrescente Corners$x[1]==Corners$x[2]
                 point.coords$x[1] <- point.coords$x[1]+rngX/20
                 point.coords$x[2] <- point.coords$x[2]-rngX/20
              } else {
                 point.coords$x<-sort(point.coords$x, decreasing=FALSE) #ordina pos$x in ordine decrescente Corners$x[1]==Corners$x[2]
                 point.coords$x[1] <- point.coords$x[1]-rngX/20
                 point.coords$x[2] <- point.coords$x[2]+rngX/20
              }
              point.coords$y <- sort(point.coords$y, decreasing=FALSE)
              Ylimits <<- c(point.coords$y[1]-rngY/10, point.coords$y[2]+rngY/10)
 	            slot(Object[[coreline]],"Boundaries") <<- point.coords
              replot()
         } )

  gbutton(" EDIT REGION ", handler = function(h, ...){
              do.editRegion()
  	      }, container = gframe22 )

  gbutton(" UNDO ", handler = function(h, ...) {
  	           undo.plot()
  	      }, container = gframe22 )

  gbutton(" RESET BOUNDARIES ", handler = function(h, ...) {
  	           reset.boundaries()
  	      }, container = gframe22 )

  gframe23 <- gframe(text = " Sprucing ", horizontal=TRUE, container = MainGroup)

  EditGroup <- ggroup(container=gframe23)
  DataToChange <- gtext(" Data to correct:  ", container=EditGroup)
  size(DataToChange) <- c(200,180)

  ButtGroup <- ggroup(horizontal=FALSE, container=gframe23)
  buttOK<-gbutton("OK", handler=function(h, ...){
              Object[[coreline]]@.Data[[1]][idx1:idx2] <<- DataTable[[1]]
              Object[[coreline]]@.Data[[2]][idx1:idx2] <<- DataTable[[2]]
              delete(EditGroup, DataToChange)
              DataToChange <<- gtext(" Data to correct:  ", container=EditGroup)
              size(DataToChange) <<- c(200,180)
              add(EditGroup, DataToChange)
              replot(Object[[coreline]])
            }, container=ButtGroup)
  buttCanc<-gbutton("Cancel", handler=function(...) {
              delete(EditGroup, DataToChange)
              DataToChange <<- gtext(" Data to correct:  ", container=EditGroup)
              size(DataToChange) <<- c(200,180)
              add(EditGroup, DataToChange)
              reset.boundaries()
            }, container=ButtGroup)
  enabled(ButtGroup) <- FALSE

#=== TKplot size ===
  Egroup3 <- ggroup(horizontal=TRUE, spacing=3, container=MainGroup)
  glabel("Graphical Window size: ", container=Egroup3)
  WSvalue <- glabel(text=as.character(WinSize), container=Egroup3)

  WSize <- gslider(from = 1, to = 2.5, by = 0.1, value = WinSize, horizontal=TRUE, handler=function(h,...){
                        WinSize <<- svalue(WSize)
                        if (Sys.info()["sysname"] == "Windows") { #sets the interactive TK graphical window dimensions
                            WinScale <<- 1.3*WinSize
                        } else if (exists("X11", env=.GlobalEnv)) {
                            WinScale <<- 0.75*WinSize
                        } else {
                            stop("tkrplot only supports Windows and X11")
                        }
                        Graph$hscale <<- WinScale
                        Graph$vscale <<- WinScale
                        tcl("update", "idletasks") #needed to get the real width and heigth of TKwinGraph
                        delete(Egroup3, WSvalue)
                        WSvalue <<- glabel(text=as.character(WinSize), container=Egroup3)
                        replot()
                        tcl("update", "idletasks") #needed to get the real width and heigth of TKwinGraph
                        hh <<- as.numeric(tkwinfo("reqheight", TKwinGraph)) #store the height of the TKRplot window
                        ww <<- as.numeric(tkwinfo("reqwidth", TKwinGraph))  #store the width of the TKRplot window
         }, container=MainGroup)


#=== CLOSE button ===
  gbutton("Save", handler = function(h, ...){
	             assign(activeFName, Object, envir = .GlobalEnv)
  	           reset.boundaries()
  	           XPSSaveRetrieveBkp("save")
         }, container = MainGroup)


  gbutton("Save & Close", handler = function(h, ...) {
	             assign(activeFName, Object, envir = .GlobalEnv)
              dispose(SPwin)
              tkdestroy(TKwinGraph)
              XPSSaveRetrieveBkp("save")
         }, container = MainGroup)

  gbutton("Exit", handler = function(h, ...) {
              dispose(SPwin)
              tkdestroy(TKwinGraph)
              XPSSaveRetrieveBkp("save")
         }, container = MainGroup)

  ## status bar
  visible(SPwin) <- TRUE


#---- TKRPLOT interactive graphics------
     graphics.off() #switch off the X11() device
     #open a TK window to replace the X11 graphical device

     if (Sys.info()["sysname"] == "Windows") { #sets the interactive TK graphical window dimensions
          WinScale <- 1.3*WinSize
     } else if (exists("X11", env=.GlobalEnv)) {
          WinScale <- 0.75*WinSize
     } else {
          stop("tkrplot only supports Windows and X11")
     }

     TKwinGraph <- tktoplevel()
     tkwm.geometry(TKwinGraph,"+565+20")   #+xxx+yyy indicate the pixels from screen TopLeftCorner(screen) to TopLeftCorner(TopLevel)
     tktitle(TKwinGraph) <- paste("TKplot: ", activeFName)
     Graph <- tkrplot(TKwinGraph, draw.plot, hscale=WinScale, vscale=WinScale)
     tkgrid(Graph, padx=0, pady=0) #load the Graph image into the frame
     tcl("update", "idletasks") #needed to get the real width and heigth of TKwinGraph
     hh <- as.numeric(tkwinfo("reqheight", Graph)) #store the height of the TKRplot window
     ww <- as.numeric(tkwinfo("reqwidth", Graph))  #store the width of the TKRplot window
     tkconfigure(Graph, bg="white") #set to white the tkrplot background = grey
     tkconfigure(Graph, cursor = "crosshair")
     tkbind(Graph, "<Button-1>", LBmousedown)   #left mouse button
     tkfocus(force=TKwinGraph)
     replot()
cat("\n Please select the portion of the spectrum to control")


#====== Markers at Coreline extremes =======

  refresh<<-FALSE #ora plotto anche il marker componente
  replot()   #riplotta spettro e marker componente selezionata

}
