#-----------------------------------------
# XPS processing with gWidgets2 and tcltk
#-----------------------------------------
#'Extract a portion of spectrum from a XPS survey
#'
#'Extract a portion (abscissa) from a XPS survey in a XPSSample. It use a GUI
#'to interact with the x-axis and to get the \code{Symbol} of the new
#'XPSCoreLine. XPSExtract() does not need any parameter.
#'
#'@return returns the \code{object} with a coreline added.
#'@examples
#'
#'\dontrun{
#' XPSextractGUI()
#'}
#'
#'@export
#'

XPSExtract <- function() {

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


  dragmouse <- function(x, y) {
      my.coords(x, y)
      point.coords$x[point.index] <<- coords[1]
      point.coords$y[point.index] <<- coords[2]
      replot()
  }


  LBmousedown <- function(x, y) {
         my.coords(x, y)
      point.coords$x[point.index] <<- coords[1]   #abscissa
      point.coords$y[point.index] <<- coords[2]   #ordinate
      if (point.index == 1) {
         point.index<<-2    #to modify the second edge of the selected area
         Corners$x<<-c(point.coords$x[1],point.coords$x[1],point.coords$x[2],point.coords$x[2])
         Corners$y<<-c(point.coords$y[1],point.coords$y[2],point.coords$y[1],point.coords$y[2])
      } else if (point.index == 2) {
         Corners$x<<-c(point.coords$x[1],point.coords$x[1],point.coords$x[2],point.coords$x[2])
         Corners$y<<-c(point.coords$y[1],point.coords$y[2],point.coords$y[1],point.coords$y[2])
         point.index<<-3
      } else if (point.index == 3) {
         D<-vector("numeric", 4)
         Dmin<-((point.coords$x[3]-Corners$x[1])^2 + (point.coords$y[3]-Corners$y[1])^2)^0.5  #initialization value
         for (ii in 1:4) {
             D[ii]<-((point.coords$x[3]-Corners$x[ii])^2 + (point.coords$y[3]-Corners$y[ii])^2)^0.5  #distance P0 - P1
             if(D[ii] <= Dmin){
                Dmin<-D[ii]
                idx=ii
             }
         }
         if (idx == 1){
            Corners$x[1]<<-Corners$x[2]<<-point.coords$x[3]
            Corners$y[1]<<-Corners$y[3]<<-point.coords$y[3]
         } else if (idx == 2){
             Corners$x[1]<<-Corners$x[2]<<-point.coords$x[3]
             Corners$y[2]<<-Corners$y[4]<<-point.coords$y[3]
         } else if (idx == 3){
             Corners$x[3]<<-Corners$x[4]<<-point.coords$x[3]
             Corners$y[1]<<-Corners$y[3]<<-point.coords$y[3]
         } else if (idx == 4){
            Corners$x[3]<<-Corners$x[4]<<-point.coords$x[3]
            Corners$y[2]<<-Corners$y[4]<<-point.coords$y[3]
         }
         if (Object@Flags[1]) { #Binding energy set
            point.coords$x<<-sort(c(Corners$x[1],Corners$x[3]), decreasing=TRUE)
            point.coords$y<<-sort(c(Corners$y[1],Corners$y[4]), decreasing=FALSE)
         } else {
            point.coords$x<<-sort(c(Corners$x[1],Corners$x[3]), decreasing=FALSE)
            point.coords$y<<-sort(c(Corners$y[1],Corners$y[4]), decreasing=FALSE)
         }
      }
      replot()
  }


  undo.plot <- function(...){
      if (SelReg == 1) {
         reset.boundaries()
         replot()
      } else if (SelReg > 1) {
         Object@Boundaries$x <<- OldCoords$x
         Ylimits <<- OldCoords$y
         replot()
      }
  }


  draw.plot <- function(...) {
      Xlimits <- Object@Boundaries$x
      if (point.index == 0) {
          plot(Object, xlim=Xlimits, ylim=Ylimits)
      } else if (point.index == 1 || point.index == 2) {
          plot(Object, xlim=Xlimits)
          points(point.coords, col="red", cex=1, lwd=1.5, pch=3)
      } else if (point.index > 2){
          plot(Object, xlim=Xlimits, ylim=Ylimits)
          points(Corners, type="p", col="red", cex=1, lwd=1.5, pch=3)
          rect(point.coords$x[1], point.coords$y[1], point.coords$x[2], point.coords$y[2])
      }
      svalue(statbar) <- sprintf(paste("x =",round(coords[1],1), " y =",round(coords[2]), sep=" "))
      parusr <<- par("usr")
      parplt <<- par("plt")
  }


  replot <- function(...) { tkrreplot(Graph) }


  reset.boundaries <- function(h, ...) {
     Object <<- XPSremove(Object, "all")
     LL<-length(Object@.Data[[1]])
     point.coords$x[1] <<- Object@.Data[[1]][1]  #abscissa of the first survey edge
     point.coords$y[1] <<- Object@.Data[[2]][1]  #ordinate of the first survey edge
     point.coords$x[2] <<- Object@.Data[[1]][LL] #abscissa of the second survey edge
     point.coords$y[2] <<- Object@.Data[[2]][LL] #ordinate of the second survey edge
     slot(Object,"Boundaries") <<- point.coords
     Ylimits<<-c(min(Object@.Data[[2]]), max(Object@.Data[[2]]))
     OldCoords <<- point.coords #for undo
     Corners <- point.coords
     point.index <<- 1
     parusr <<- par("usr")
     parplt <<- par("plt")
     replot()
  }


  do.extract <- function(h, ...){
         slot(Object,"Boundaries") <<- point.coords
         Xlimits <<- Object@Boundaries$x   #visualize selected region
         Ylimits <<- sort(Object@Boundaries$y, decreasing=FALSE)
         point.index <<- 0
         replot()

         winExt <- gwindow("Extract Spectral Feature", visible=FALSE, parent=window) #
         gBox <- gvbox(container=winExt)
         flyt <- gformlayout(container=gBox)
         elesymbol <- gedit("", label="Element Name:", container=flyt)
         gseparator(container=gBox)
         bg <- ggroup(container=gBox)
         gbutton("OK", container=bg, handler=function(...){  #input spectrum name
                 Symbol <- svalue(elesymbol)
                 Symbol <- gsub(" ", "", Symbol)  #eliminates white spaces from Symbol
                 pattern <- c("[[:alpha:]]{1,2}")     #matches only the first two char
                 mpat <- regexpr(pattern, Symbol)
                 ## symbol element
                 Element <- regmatches(Symbol, mpat)
                 if (ElementCheck(Element)==FALSE ) {    #see XPSelement.r
                     yesno<-gconfirm(msg=" ATTENTION: element Name NOT found in Element Table! \n Proceed anyway?", icon="warning")
                     if (yesno==FALSE){
                         dispose(winExt)
                         return()
                    }
                 }
                 dispose(winExt)

                 newcoreline <- Object
                 Xmax <- max(range(newcoreline@.Data[1]))
                 Xmin <- min(range(newcoreline@.Data[1]))
                 #is selected region out of limits?
                 if (point.coords$x[1] > Xmax) {point.coords$x[1] <<- Xmax}
                 if (point.coords$x[1] < Xmin) {point.coords$x[1] <<- Xmin}
                 if (point.coords$x[2] > Xmax) {point.coords$x[2] <<- Xmax}
                 if (point.coords$x[2] < Xmin) {point.coords$x[2] <<- Xmin}

                 idx1 <- findXIndex(unlist(newcoreline@.Data[1]), point.coords$x[1]) #index corresponding to the selected BE1 (or KE1 value) of RegionToFit
                 idx2 <- findXIndex(unlist(newcoreline@.Data[1]), point.coords$x[2]) #index corresponding to the selected BE2 (or KE2 value) of RegionToFit
                 tmp <- unlist(Object@.Data[1])  #extract correspondent X values for the selected region
                 newcoreline@.Data[[1]] <- tmp[idx1:idx2]    #save the X values in the new coreline
                 OldEnergyScale <<- newcoreline@.Data[[1]]
                 newcoreline@Boundaries$x <- c(tmp[idx1], tmp[idx2])
                 tmp <- unlist(Object@.Data[2])  #extract correspondent Y values for the selected region
                 newcoreline@.Data[[2]] <- tmp[idx1:idx2]    #save the Y values in the new coreline
                 newcoreline@Boundaries$y <- range(tmp)
                 tmp <- unlist(Object@.Data[3])  #extract correspondent transmission Factor values for the selected region
                 newcoreline@.Data[[3]] <- tmp[idx1:idx2]    #save the transmission Factor values in the new coreline
                 slot(newcoreline,"Symbol") <- Symbol
                 ## add extracted coreline to original XPSSample
                 idx <- length(XPSSample) + 1
                 XPSSample[[idx]] <<- newcoreline
                 names(XPSSample) <<- unname(sapply(XPSSample, slot, "Symbol"))
                 point.coords$x[1] <<- Object@.Data[[1]][1]
                 point.coords$y[1] <<- Object@.Data[[2]][1]
                 point.coords$x[2] <<- Object@.Data[[1]][LL]
                 point.coords$y[2] <<- Object@.Data[[2]][LL]
                 svalue(statbar) <- sprintf("New %s added.", svalue(elesymbol))
                 Object <<- XPSSample[[idx]]
         })
         gbutton("Cancel", container=bg, handler=function(...) dispose(winExt))
         visible(winExt) <- TRUE
  }


#----- Variables -----

  XPSSample <- get(activeFName, envir=.GlobalEnv)       #load XPSdata values from main memory
  coreline <- get("activeSpectIndx", envir=.GlobalEnv)
  Object <- XPSSample[[coreline]]
  point.coords <- list(x=NA,y=NA)
  point.index <- 1
  coords <- NA # for printing mouse coordinates on the plot
  xx <- NULL
  yy <- NULL
  NO.Fit <- FALSE

#Coreline boundaries
  LL<-length(Object@.Data[[1]])
  point.coords$x[1] <- Object@.Data[[1]][1]
  point.coords$y[1] <- Object@.Data[[2]][1]
  point.coords$x[2] <- Object@.Data[[1]][LL]
  point.coords$y[2] <- Object@.Data[[2]][LL]
  Xlimits <- c(min(Object@.Data[[1]]), max(Object@.Data[[1]]))
  Ylimits <- c(min(Object@.Data[[2]]), max(Object@.Data[[2]]))
  Corners <- point.coords
  Object@Boundaries$x <- c(point.coords$x)
  Object@Boundaries$y <- c(point.coords$y)

  OldCoords <- point.coords #for undo
  OldEnergyScale <- Object@.Data[[1]]
  OldFlag <- Object@Flags[1]
  OldUnits <- Object@units[1]
  SelReg <- 0
  FNameList <- XPSFNameList()
  SpectList <- XPSSpectList(activeFName)

  WinSize <- as.numeric(XPSSettings$General[4])
  hh <- NULL
  ww <- NULL
  WinScale  <- NULL
  parusr <- NULL
  parplt <- NULL

#====== Widget definition =======
  Ewindow <- gwindow("XPS extract GUI", parent=c(100, 0), toolkit = "tcltk", visible = FALSE)
  addHandlerDestroy(Ewindow, handler=function(h, ...){  #if MainWindow unproperly closed with X
                               tkdestroy(TKwinGraph)    #closes the TK window
                               Gdev <- XPSSettings$General[6] #retrieve the graphic device used
                               eval(parse(text=Gdev),envir=.GlobalEnv) #restore the graphic window
                               plot.new()
                               plot(XPSSample) #replot the CoreLine
                           })
  Egroup1 <- ggroup(container = Ewindow, horizontal = TRUE)

  ## XPSSample and Core lines
  Egroup2 <- ggroup(expand = FALSE, horizontal = FALSE, spacing = 5, container = Egroup1)

  gframe20 <- gframe(text = " XPS Sample and Core line Selection ", container = Egroup2)
  XPS.Sample <- gcombobox(FNameList, selected=-1, editable=FALSE, handler=function(h,...){
                                 activeFName<-svalue(XPS.Sample)
                                 XPSSample <<- get(activeFName, envir=.GlobalEnv)
                                 SpectList <<- XPSSpectList(activeFName)
                                 indx <- grep("Survey", SpectList, value=FALSE)
                                 Object <<- XPSSample[[indx[1]]]
                                 point.coords$x[1] <<- Object@.Data[[1]][1]
                                 point.coords$y[1] <<- Object@.Data[[2]][1]
                                 point.coords$x[2] <<- Object@.Data[[1]][LL]
                                 point.coords$y[2] <<- Object@.Data[[2]][LL]
                                 delete(gframe20, Core.Lines)
                                 Core.Lines <<- gcombobox(SpectList, selected=1, expand = FALSE, handler = function(h, ...){
                                                CLine <- svalue(Core.Lines)
                                                CLine <- unlist(strsplit(CLine, "\\."))   #"number." and "CL name" are separated
                                                if (CLine[2] != "survey" && CLine[2] != "Survey"){
#                                                   gmessage("Wrong coreline: please select a survey spectrum", title="WRONG CORELINE", icon="warning")
                                                   yesno<-gconfirm(msg=" ATTENTION: This is NOT a Survey spectrum! \n Proceed anyway?", icon="warning")
                                                   if (yesno==FALSE){
                                                       return()
                                                   }
                                                }
                                                coreline <<- as.integer(CLine[1])   # Coreline index
                                                Object <<- XPSSample[[coreline]]
                                                LL<-length(Object@.Data[[1]])
                                                point.coords$x[1] <<- Object@.Data[[1]][1]
                                                point.coords$y[1] <<- Object@.Data[[2]][1]
                                                point.coords$x[2] <<- Object@.Data[[1]][LL]
                                                point.coords$y[2] <<- Object@.Data[[2]][LL]
                                                Xlimits <<- c(min(Object@.Data[[1]]), max(Object@.Data[[1]]))
                                                Ylimits <<- c(min(Object@.Data[[2]]), max(Object@.Data[[2]]))
                                                Corners <<- point.coords
                                                Object@Boundaries$x <<- c(point.coords$x)
                                                Object@Boundaries$y <<- c(point.coords$y)
                                                parusr <<- NA
                                                parplt <<- NA
                                                OldCoords <<- point.coords #for undo
                                                OldEnergyScale <<- Object@.Data[[1]]
                                                OldFlag <<- Object@Flags[1]
                                                OldUnits <<- Object@units[1]
                                                SelReg <<- 0
                                                replot()
                                 }, container = gframe20)
                                 replot()
                       }, container = gframe20)
  svalue(XPS.Sample) <- activeFName

  Core.Lines <- gcombobox(SpectList, selected=-1, expand = FALSE, handler = function(h, ...){
              CLine <- svalue(Core.Lines)
              CLine <- unlist(strsplit(CLine, "\\."))   #"number." and "CL name" are separated
              if (CLine[2] != "survey" && CLine[2] != "Survey"){
#                 gmessage("Wrong coreline: please select a survey spectrum", title="WRONG CORELINE", icon="warning")
                  yesno<-gconfirm(msg=" ATTENTION: This is NOT a Survey spectrum! \n Proceed anyway?", icon="warning")
                  if (yesno==FALSE){
                      return()
                  }
              }
              coreline <<- as.integer(CLine[1])   # Coreline index
              Object <<- XPSSample[[coreline]]
              LL <- length(Object@.Data[[1]])
              point.coords$x[1] <<- Object@.Data[[1]][1]
              point.coords$y[1] <<- Object@.Data[[2]][1]
              point.coords$x[2] <<- Object@.Data[[1]][LL]
              point.coords$y[2] <<- Object@.Data[[2]][LL]
              Xlimits <<- c(min(Object@.Data[[1]]), max(Object@.Data[[1]]))
              Ylimits <<- c(min(Object@.Data[[2]]), max(Object@.Data[[2]]))
              Corners <<- point.coords
              Object@Boundaries$x <<- c(point.coords$x)
              Object@Boundaries$y <<- c(point.coords$y)
              parusr <<- NA
              parplt <<- NA
              OldCoords <<- point.coords #for undo
              OldEnergyScale <<- Object@.Data[[1]]
              OldFlag <<- Object@Flags[1]
              OldUnits <<- Object@units[1]
              SelReg <<- 0
              replot()
         }, container = gframe20)

  gframe22 <- gframe(text = " Options ", horizontal = FALSE, container = Egroup2)
  gbutton("SELECT REGION", handler = function(h, ...){
              OldCoords <<- Object@Boundaries
              SelReg <<- SelReg+1
              rngX <- range(point.coords$x)
              rngX <- (rngX[2]-rngX[1])/20
              rngY <- range(point.coords$y)
              rngY <- (rngY[2]-rngY[1])/20
              if (Object@Flags[1]) { #Binding energy set
                 point.coords$x <- sort(point.coords$x, decreasing=TRUE)  #pos$x in decreasing order
                 point.coords$x[1] <- point.coords$x[1]+rngX/20
                 point.coords$x[2] <- point.coords$x[2]-rngX/20
              } else {
                 point.coords$x <- sort(point.coords$x, decreasing=FALSE) #pos$x in increasing order
                 point.coords$x[1] <- point.coords$x[1]-rngX/20
                 point.coords$x[2] <- point.coords$x[2]+rngX/20
              }
              point.coords$y <- sort(point.coords$y, decreasing=FALSE)
              Ylimits <<- c(point.coords$y[1]-rngY/10, point.coords$y[2]+rngY/10)
              slot(Object,"Boundaries") <<- point.coords
              replot()
         }, container = gframe22 )

  gbutton("EXTRACT REGION", handler = function(h, ...){
              do.extract()
         }, container = gframe22 )

  gbutton("UNDO", handler = function(h, ...) {
              undo.plot()
         }, container = gframe22 )

  gbutton("RESET BOUNDARIES", handler = function(h, ...) {
              Object@.Data[[1]]<<-OldEnergyScale
              Object@Flags[1]<<-OldFlag
              reset.boundaries()
         }, container = gframe22 )

  gframe23 <- gframe(text = " Plot ", container = Egroup2)
  SwitchE <- gcheckbox(" SWITCH BINDING/KINETIC ENERGY SCALE", checked=FALSE, handler = function(h, ...) {
                     XEnergy<-get("XPSSettings", envir=.GlobalEnv)$General[5] #the fifth element of the first column of XPSSettings
                     XEnergy<-as.numeric(XEnergy)
                     if (svalue(SwitchE)==TRUE){
                        idx1<-findXIndex(unlist(Object@.Data[[1]]), point.coords$x[1]) #index relative to the upper limit (lower if KE scale) of the X scale
                        idx2<-findXIndex(unlist(Object@.Data[[1]]), point.coords$x[2]) #index relative to the lower limit (upper if KE scale) of the X scale
                        Object@.Data[[1]]<<-XEnergy-OldEnergyScale
                        Object@Boundaries$x <<- c(Object@.Data[[1]][idx1],Object@.Data[[1]][idx2])
                        point.coords$x<<-XEnergy-point.coords$x
                        Corners$x<<-XEnergy-Corners$x
                        if (Object@Flags[1]==TRUE) { #original scaale is BE
                           Object@Flags[1]<<-FALSE   #set KE scale
                           Object@units[1]<<-"Kinetic Energy [eV]"
                        } else if (Object@Flags[1]==FALSE) { #original scaale is KE
                           Object@Flags[1]<<-TRUE   #set BE scale
                           Object@units[1]<<-"Binding Energy [eV]"
                        }
                     }
                     if (svalue(SwitchE)==FALSE){
                        idx1<-findXIndex(unlist(Object@.Data[[1]]), point.coords$x[1]) #index relative to the upper limit (lower if KE scale) of the X scale
                        idx2<-findXIndex(unlist(Object@.Data[[1]]), point.coords$x[2]) #index relative to the lower limit (upper if KE scale) of the X scale
                        Object@.Data[[1]]<<-OldEnergyScale
                        Object@Boundaries$x <<- c(Object@.Data[[1]][idx1],Object@.Data[[1]][idx2])
                        point.coords$x<<-XEnergy-point.coords$x
                        Corners$x<<-XEnergy-Corners$x
                        Object@Flags[1]<<-OldFlag
                        Object@units[1]<<-OldUnits
                     }
                     replot()
         }, container=gframe23 )

  gframe24 <- gframe(text = "TK Graphic Window Dimensions", horizontal=FALSE, container = Egroup2)
  Egroup3 <- ggroup(horizontal=TRUE, spacing=3, container=gframe24)
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
         }, container=gframe24)

#---- Buttons ----
  gbutton("SAVE", expand=FALSE, handler = function(h, ...){
              Object@.Data[[1]] <<- OldEnergyScale
              Object@Flags[1] <<- OldFlag
              Object@units[1] <<- OldUnits
              assign(activeFName, XPSSample, envir = .GlobalEnv)
              reset.boundaries()
              replot()
              XPSSaveRetrieveBkp("save")
         }, container = Egroup2 )


  gbutton("SAVE & EXIT", expand=FALSE, handler = function(h, ...){
              Object@.Data[[1]] <<- OldEnergyScale
              Object@Flags[1] <<- OldFlag
              Object@units[1] <<- OldUnits
              assign(activeFName, XPSSample, envir = .GlobalEnv)
              reset.boundaries()
              dispose(Ewindow)  #this calls the handlerdestroy(Ewindow...)
              tkdestroy(TKwinGraph)
              XPSSaveRetrieveBkp("save")
         }, container = Egroup2 )

  gbutton("EXIT", expand=FALSE, handler = function(h, ...){
              dispose(Ewindow)
              tkdestroy(TKwinGraph)
              XPSSaveRetrieveBkp("save")
         }, container = Egroup2 )

#---- Status bar
  statbar <- gstatusbar("status", container = Ewindow)

  visible(Ewindow) <- TRUE

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
cat("\n Please select the portion of the Survey containing the CoreLine to extract")
}
