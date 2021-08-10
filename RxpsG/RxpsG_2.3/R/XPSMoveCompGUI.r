#   Manual Move Components with tcltk and rpanel
#   FUNCTION CALL:  XPSMoveCompGUI(c2.pxt,CoreLine=2,hscale=1.8)

#'Function to modify component position and intensity in a fit
#'
#'Provides a userfriendly interface change position and intensity of each
#'individual fitting component of a selected XPSCoreline. Changes are saved in the
#'.GlobalEnv main software memory
#'
#'
#'@param XPSSample XPSSample object
#'@param Indx  index identifying the CoreLine to work on
#'@param hscale  dimension of the plotting window
#'@examples
#'
#'\dontrun{
#'	XPSMoveComp(SampData, 2, "medium")
#'}
#'
#'@export
#'


XPSMoveComp <- function(){

   import::here(tktoplevel, tkwm.geometry, tkgrid, tkconfigure,
                tkframe, tktitle, tkwinfo, tkbind, tkfocus, .from=tcltk )
   import::here(tkrplot, tkrreplot, .from=tkrplot )

   LoadCoreLine<-function(){
      XPSSample <<- get(activeFName, envir=.GlobalEnv)     #load the XPSSample data from main memory
      Indx <<- get("activeSpectIndx", envir=.GlobalEnv)    #get active Spectrum index
      Object <<- XPSSample[[Indx]]
      Xlimits <<- range(Object@RegionToFit$x)
      Ylimits <<- range(Object@RegionToFit$y)
      ComponentList <<- names(slot(Object,"Components"))
      if (length(ComponentList)==0) {
          gmessage("ATTENTION NO FIT FOUND: change coreline please!" , title = "WARNING",  icon = "warning")
          return()
      }
      delete(MCFrame2, FitComp) # selecting a new core line needs a new gradio running on the coreline fit-components

      if (length(ComponentList)>1){    #gradio works with at least 2 items Less than 2 items gcheckbox will be used
          FitComp <<- gradio(ComponentList, selected=1, handler = function(h,...){   #gradio handler has to be redefined
                   refresh <<- TRUE
                   replot()   #replot spectrum without marker
                   FComp <- svalue(FitComp)
                   FComp <- as.numeric(unlist(strsplit(FComp, split="C")))  #index of the selected component
                   FComp <- FComp[2]
                   xx <- Object@Components[[FComp]]@param[2,1] #component position mu
                   Rng <- range(Object@RegionToFit[[1]])
                   if (xx < Rng[1]) {xx <- Rng[1]}
                   if (xx > Rng[2]) {xx <- Rng[2]}
                   FuncName <- Object@Components[[FComp]]@funcName
                   yy <- Object@Components[[FComp]]@param[1,1] #component height h
                   yy <- yy/GetHvalue(Object,FComp, FuncName, 1)
                   Estep <- abs(Object@RegionToFit[[1]][1]-Object@RegionToFit[[1]][2])
                   Xindx <- which(Object@RegionToFit[[1]] > xx-Estep/2 & Object@RegionToFit[[1]] < xx+Estep/2)
                   yy <- yy+Object@Baseline$y[Xindx]  #spectral intensity + baseline at point xx
                   coords[1] <<- xx
                   coords[2] <<- yy
                   refresh <<- FALSE #now plot also the component marker
                   replot()   #replot all
               },  container = MCFrame2)
      } else {
         FitComp <<- gcheckboxgroup(ComponentList, checked=TRUE, handler = function(h,...){ #gradio handler has to be redefined
                   refresh<<-TRUE
                   replot()   #replot spectrum without marker
                   FComp <- svalue(FitComp)
                   FComp <- as.numeric(unlist(strsplit(FComp, split="C")))
                   FComp <- FComp[2]
                   xx <- Object@Components[[FComp]]@param[2,1] #component position mu
                   Rng <- range(Object@RegionToFit[[1]])
                   if (xx < Rng[1]) {xx <- Rng[1]}
                   if (xx > Rng[2]) {xx <- Rng[2]}
                   FuncName <- Object@Components[[FComp]]@funcName
                   yy <- Object@Components[[FComp]]@param[1,1] #component height h
                   yy <- yy/GetHvalue(Object,FComp, FuncName, 1)
                   Estep <- abs(Object@RegionToFit[[1]][1]-Object@RegionToFit[[1]][2])
                   Xindx <- which(Object@RegionToFit[[1]] > xx-Estep/2 & Object@RegionToFit[[1]] < xx+Estep/2)
                   yy <- yy+Object@Baseline$y[Xindx]  #spectral intensity + baseline at point xx
                   coords[1] <<- xx
                   coords[2] <<- yy
                   refresh <<- FALSE #now plot also the component marker
                   replot()   #plot all
              },  container = MCFrame2)
      }
      xx <- Object@Components[[1]]@param[2,1] #component position mu
      yy <- Object@Components[[1]]@param[1,1] #component height h
      Estep <- abs(Object@RegionToFit[[1]][1]-Object@RegionToFit[[1]][2])
      Xindx <- which(Object@RegionToFit[[1]] > xx-Estep/2 & Object@RegionToFit[[1]] < xx+Estep/2)
      yy <- yy+Object@Baseline$y[Xindx]  #spectral intensity + baseline at point xx
      coords[1] <<- xx
      coords[2] <<- yy
      refresh <<- FALSE #now plot the component marker
      replot()   #replot spectrum and marker of selected fit component
  }

  my.coord <- function(xx, yy){     #giving cursor position in pixels here compute cursor coordinates
      coords <- NULL

      #xx, yy are given in pixel. Now normalize
      xx <- as.numeric(xx)/ww
      yy <- 1 - as.numeric(yy)/hh
      #now bottom left corner is (0,0), top right corner is (1,1)

      XpltRng <- parplt[2]-parplt[1]    #this is the fraction of the window used for plot
      YpltRng <- parplt[4]-parplt[3]

      XusrRng <- parusr[2] - parusr[1]  #this is x box extension containing the spectrum
      YUsrRng <- parusr[4] - parusr[3]

     #      Eliminate borders: shifts the window (0,0) in the plot (0,0)
     #          |           Normalize the fraction of the window used for plot i.e. plot bottom left corner is (0,0), top right corner is (1,1)
     #          |              |      transforms the X, Y [0 - 1] range in the user range
     #          |              |        |          Shifts the scale in the correct positgion adding the Xmin, Ymin
     #          |              |        |          |
     #       -----------   --------  -------   ---------
      xx <- (xx-parplt[1]) /XpltRng  *XusrRng + parusr[1]
      yy <- (yy-parplt[3]) /YpltRng  *YUsrRng + parusr[3]

      Xlim1 <- min(Object@RegionToFit[[1]])   #limits coordinates in the Spectrum Range
      Xlim2 <- max(Object@RegionToFit[[1]])
      Ylim1 <- min(Object@RegionToFit[[2]])
      Ylim2 <- max(Object@RegionToFit[[2]])
      if (xx < Xlim1 ) {xx <- Xlim1}
      if (xx > Xlim2 ) {xx <- Xlim2}
      if (yy < Ylim1 ) {yy <- Ylim1}
      if (yy > Ylim2 ) {yy <- Ylim2}

      Estep <- abs(Object@RegionToFit[[1]][1]-Object@RegionToFit[[1]][2])
      Xindx <- which(Object@RegionToFit[[1]]>xx-Estep/2 & Object@RegionToFit[[1]]<xx+Estep/2)
      yy_BasLin <- yy-Object@Baseline$y[Xindx]  #spectral intensity at xx without Baseline
      coords <<- c(xx, yy, yy_BasLin)
      return()
  }




  LBmousedown <- function(x, y) {   #Left mouse button down
     my.coord(x, y)
	    xx <- coords[1]
	    yy <- coords[2]
     if (Object@Flags[1]) { #Binding energy set
        if (xx > parusr[1] || xx < parusr[2] || yy < parusr[3] || yy > parusr[4]){ return() } #if you click ouside XY box of the plot do noting
     } else {
        Xlimits<<-sort(c(point.coords$x[1], point.coords$x[2]), decreasing=FALSE) #pos$x in increasing order
        if (xx < parusr[1] || xx > parusr[2] || yy < parusr[3] || yy > parusr[4]){ return() }
     }
     if (SetZoom==FALSE) { #left button works only when SET ZOOM REGION inactive
        do.move()
        XPSquantify(XPSSample)
        refresh <<- FALSE
     }
     replot()  
  }


  RBmousedown <- function(x, y) {   #Right mouse button down
     my.coord(x, y)
	    yy <- coords[2]
	    xx <- coords[1]
     if (Object@Flags[1]) { #Binding energy set
        if (xx > parusr[1] || xx < parusr[2] || yy < parusr[3] || yy > parusr[4]){ return() } #if you click ouside XY box of the plot do noting
     } else {
        Xlimits <<- sort(c(point.coords$x[1], point.coords$x[2]), decreasing=FALSE) #pos$x in increasing order
        if (xx < parusr[1] || xx > parusr[2] || yy < parusr[3] || yy > parusr[4]){ return() }
     }
     if (SetZoom==TRUE) { #left button works only when SET ZOOM REGION button pressed
     	  point.coords$x[point.index] <<- coords[1]   #abscissa
     	  point.coords$y[point.index] <<- coords[2]   #ordinate
     	  if (point.index==1) {
     	     point.index <<- 2    #to modify the second edge of the selected area
           Corners$x <<- c(point.coords$x[1],point.coords$x[1],point.coords$x[2],point.coords$x[2])
           Corners$y <<- c(point.coords$y[1],point.coords$y[2],point.coords$y[1],point.coords$y[2])
  	     }  else if (point.index==2) {
           Corners$x <<- c(point.coords$x[1],point.coords$x[1],point.coords$x[2],point.coords$x[2])
           Corners$y <<- c(point.coords$y[1],point.coords$y[2],point.coords$y[1],point.coords$y[2])
           point.index <<- 3
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
              Corners$x[1] <<- Corners$x[2] <<- point.coords$x[3]
              Corners$y[1] <<- Corners$y[3] <<- point.coords$y[3]
           } else if (idx==2){
              Corners$x[1] <<- Corners$x[2] <<- point.coords$x[3]
              Corners$y[2] <<- Corners$y[4] <<- point.coords$y[3]
           } else if (idx==3){
              Corners$x[3] <<- Corners$x[4] <<- point.coords$x[3]
              Corners$y[1] <<- Corners$y[3] <<- point.coords$y[3]
           } else if (idx==4){
              Corners$x[3] <<- Corners$x[4] <<- point.coords$x[3]
              Corners$y[2] <<- Corners$y[4] <<- point.coords$y[3]
           }
           if (Object@Flags[1]) { #Binding energy set
              point.coords$x <<- sort(c(Corners$x[1],Corners$x[3]), decreasing=TRUE) #pos$x in decreasing order
              point.coords$y <<- sort(c(Corners$y[1],Corners$y[4]), decreasing=FALSE)
           } else {
              point.coords$x <<- sort(c(Corners$x[1],Corners$x[3]), decreasing=FALSE) #pos$x in increasing order
              point.coords$y <<- sort(c(Corners$y[1],Corners$y[4]), decreasing=FALSE)
           }
        }
        replot()
      }
  }


  do.move <- function(...) {
     FComp <- svalue(FitComp)
     if (length(FComp)==0) {
         gmessage(msg="Select Component Please", title="WARNING", icon = "warning")
     } else {
        FComp<-as.numeric(unlist(strsplit(FComp, split="C")))   #index of the selected component
        FComp<-FComp[2]
        xx <- coords[1]
        yy <- coords[2]  #Component max value with baseline
        zz <- coords[3]  #Component max value without baseline
        FitFunct<-Object@Components[[FComp]]@funcName
        newh<-GetHvalue(Object, FComp, FitFunct, zz)  #Get value computes the Component value given the fit parameters and the Ymax value

        #range limits for mu
        varmu <- getParam(Object@Components[[FComp]],variable="mu")
        minmu <- varmu$start-varmu$min
        maxmu <- varmu$max-varmu$start
        newmu <- c(xx, xx-minmu, xx+maxmu)
        #range limits for h
        varh <- getParam(Object@Components[[FComp]],variable="h")
        minh <- varh$start-varh$min
        maxh <- varh$max-varh$start

        if (maxh > 0) {
            newh <- c(newh, 0, newh*5)    # No constraints on h
        }
        if (maxh==0){
            newh <- c(newh, newh, newh)   # h is fixed
        }
        if (maxh<0){
            newh <- c(newh, 0, newh*5)    # maxh cannot be <0: => force newH to correct values
        }
        if (varh$start <0) {
            newh <- c(0, 0, 1e5)   #set a positive value for an hypotheic fit
        }
        Object@Components[[FComp]] <<- setParam(Object@Components[[FComp]], parameter=NULL, variable="mu", value=newmu)
        Object@Components[[FComp]] <<- setParam(Object@Components[[FComp]], parameter=NULL, variable="h", value=newh)
        Object@Components[[FComp]] <<- Ycomponent(Object@Components[[FComp]], x=Object@RegionToFit$x, y=Object@Baseline$y) #eomputes the Y value and add baseline
# Fit computed addind fit components with the modified ones
        tmp <- sapply(Object@Components, function(z) matrix(data=z@ycoor))
        Object@Fit$y <<- ( colSums(t(tmp)) - length(Object@Components)*(Object@Baseline$y))

        Object <<- sortComponents(Object)
#if component order changed then re-number them
        LL<-length(Object@Components) #N. fit components
        for (ii in 1:LL){
           if (xx == Object@Components[[ii]]@param["mu",1]) { #compare marker position with component positions
             indx <- ii
             break()
           }
        }
        svalue(FitComp) <- paste("C", indx, sep="")  #update component gradio
        XPSSample[[Indx]] <<- Object
     }
  }


  draw.plot <- function(...) {
     XPSSample <- get(activeFName, envir=.GlobalEnv) #Load XPSSample because changes in XPSSample[[activeSpectIndx]]
                                                     #may derive from changes of fit-parameters in XPSConstraintsGui()
                                                     #which are saved in the global variable activeFName[[activeSpectIndx]]
     Object <<- XPSSample[[activeSpectIndx]]
     if (point.index==1 && refresh==FALSE) {  #point.index==1 when moving mcomponent
         plot(Object, xlim=Xlimits, ylim=Ylimits)
         points(x=coords[1], y=coords[2], col=2, cex=1.2, lwd=2, pch=1)  # if refresh==FALSE plot spectrum with component marker
     } else if (SetZoom == TRUE){   #set zoom area corners
	        if (point.index <= 2) {  #define zoom area corners
 	           plot(Object)
             points(point.coords, type="p", col=3, cex=1.2, lwd=2.5, pch=3)
  	      } else if (point.index > 2){  #plot zoom area corners
 	           plot(Object, xlim=Xlimits, ylim=Ylimits)
             points(Corners, type="p", col=3, cex=1.2, lwd=2.5, pch=3)
             rect(point.coords$x[1], point.coords$y[1], point.coords$x[2], point.coords$y[2])
         }
     } else {
         plot(Object, xlim=Xlimits, ylim=Ylimits)
     }
     parusr <<- par('usr')
     parplt <<- par('plt')
     svalue(StatBar) <- paste("x =",round(coords[1],2), "y =",round(coords[2],2), sep="   ")
  }

  reset.plot <- function(h, ...) {
       point.coords$x <<- range(Object@RegionToFit$x) #set original X range
       point.coords$y <<- range(Object@RegionToFit$y) #set original Y range
       Object@Boundaries <<- point.coords
       Xlimits <<- point.coords$x
       Ylimits <<- sort(point.coords$y, decreasing=FALSE)
       Corners <<- point.coords
       parusr <<- par("usr")
       parplt <<- par("plt")
       tkconfigure(Graph, cursor = "crosshair")
       replot()
  }


  replot <- function(...) {
#XPSMoveComp() may be used in combination with XPSConstraintsGUI()
     XPSSample[[activeSpectIndx]] <<- Object #save changes done on OBJECT in .GlobalEnv
     assign(activeFName, XPSSample, envir = .GlobalEnv)  #in Draw.Plot we recall changes from .GlobalEnv
     tkrreplot(Graph)                        # changes may derive also from XPSFitConstraints()
  }


  refresh.plot <- function(...) {
     plot(Object)
     tkrreplot(Graph)
  }


  reset.vars <- function(){
     XPSSample <<- get(activeFName, envir=.GlobalEnv)       #load XPSdata values from main memory
     SpectName <<- get("activeSpectName", envir=.GlobalEnv)
     Indx <- grep(SpectName, XPSSpectList(activeFName))     #set the CoreLine Indx
     Object <<- XPSSample[[Indx]]
     ComponentList <<- names(slot(Object,"Components"))
     if (length(ComponentList)==0) {
        gmessage("ATTENTION NO FIT FOUND: change coreline please!" , title = "WARNING",  icon = "warning")
        Xlimits <<- range(Object@.Data[1])
        Ylimits <<- range(Object@.Data[2])
        NoFit <<- TRUE
        return()
     }

     FComp <<- svalue(FitComp)
     FComp <<- as.numeric(unlist(strsplit(FComp, split="C")))   #index selected component
     FComp <<- FComp[2]

     coords <<- c(xx=NA, yy=NA, yy_BasLin=NA)
     CompCoords <<- c(xx=NA, yy=NA, yy_BasLin=NA)
     FuncName <- Object@Components[[FComp]]@funcName
     vscale <<- GetHvalue(Object,FComp, FuncName, 1)
     xx <<- Object@Components[[FComp]]@param[2,1] #component position mu
     Rng <- range(Object@RegionToFit[[1]])
     if (xx < Rng[1]) {xx <- Rng[1]}
     if (xx > Rng[2]) {xx <- Rng[2]}
     Estep <- abs(Object@RegionToFit[[1]][1]-Object@RegionToFit[[1]][2])
     Xindx <- which(Object@RegionToFit[[1]] > xx-Estep/2 & Object@RegionToFit[[1]] < xx+Estep/2)
     yy <<- Object@Components[[FComp]]@param[1,1] #component height h
     yy <<- yy/hscale+Object@Baseline$y[Xindx]  #spectral intensity + baseline at point xx
     coords[1] <<- xx
     coords[2] <<- yy
     point.coords$x <<- range(Object@RegionToFit$x) #ordinata secondo estremo del survey
     point.coords$y <<- range(Object@RegionToFit$y) #ordinata secondo estremo del survey
     Xlimits <<- range(Object@RegionToFit$x)
     Ylimits <<- range(Object@RegionToFit$y)
     Object@Boundaries$x <<- Xlimits
     Object@Boundaries$y <<- Ylimits
     Corners <<- point.coords
     point.index <<- 1
     refresh <<- TRUE
     SetZoom <<- FALSE
     NoFit <<- FALSE
  }



# --- Varables definition ---
     XPSSample <- get(activeFName, envir=.GlobalEnv)       #load XPSdata values from main memory
     SpectName <- get("activeSpectName", envir=.GlobalEnv)
     Indx <- get("activeSpectIndx", envir=.GlobalEnv)
     Object <- XPSSample[[Indx]]
     ComponentList <-names(slot(Object,"Components"))
     FNameList <- XPSFNameList()
     SpectList <- XPSSpectList(activeFName)
     FComp <- NULL

     WinSize <- as.numeric(XPSSettings$General[4])
     hscale <- hscale <- as.numeric(WinSize)
     hh <- NULL
     ww <- NULL
     WinScale  <- NULL
     parusr <- par("usr")
     parplt <- par("plt")

     coords <- c(xx=NA, yy=NA, yy_BasLin=NA)
     CompCoords <- c(xx=NA, yy=NA, yy_BasLin=NA)
     point.coords <- list(x=range(Object@RegionToFit$x), #set original X range
                          y=range(Object@RegionToFit$y)) #set original Y range
     xx <- NULL
     yy <- NULL
     Corners <- point.coords
     point.index <- 1

     refresh <- TRUE
     SetZoom <- FALSE
     NoFit <- FALSE
#Coreline boundaries
     if (length(ComponentList)==0) {
        gmessage("ATTENTION NO FIT FOUND: change coreline please!" , title = "WARNING",  icon = "warning")
        Xlimits <- range(Object@.Data[1])
        Ylimits <- range(Object@.Data[2])
        NoFit <- TRUE
     } else {
        LL<-length(Object@.Data[[1]])
        point.coords$x <- range(Object@RegionToFit$x) #set the X window extension == x range
        point.coords$y <- range(Object@RegionToFit$y) #set the Y window extension == y range
        Xlimits <- range(Object@RegionToFit$x)
        Ylimits <- range(Object@RegionToFit$y)
        Object@Boundaries$x <- Xlimits
        Object@Boundaries$y <- Ylimits
     }


#--- GWidget definition ---

     MCWindow <- gwindow("XPS MOVE COMPONENT", parent=c(100,0), toolkit = "tcltk", visible = FALSE)
     addHandlerDestroy(MCWindow, handler=function(h, ...){  #if MainWindow unproperly closed with X
                               tkdestroy(TKwinGraph)           #closes the TK window
                               X11()                           #restore the X11() graphic window
                               plot.new()
                               plot(XPSSample[[Indx]])         #replot the CoreLine
                               rm(list=c("TKwinGraph", "Graph"), envir=.GlobalEnv) #clean .GlobalEnv
                           })

#--- Selection Group ---
     SelectGroup <- ggroup(horizontal=FALSE, spacing = 5, container=MCWindow)

     MCFrame1 <- gframe(text = " XPS Sample ", container = SelectGroup)

     XPS.Sample <- gcombobox(FNameList, selected=-1, editable=FALSE, handler=function(h,...){
                               ActFName <- svalue(XPS.Sample)
                               assign("activeFName", ActFName, envir=.GlobalEnv)
                               reset.vars()
                               if (NoFit == TRUE) { return() }
                               LoadCoreLine()
                               refresh <<- FALSE #now plot also the component marker
                               Object <<- XPSSample[[Indx]]
                               replot()
                           }, container = MCFrame1)
     svalue(XPS.Sample) <- activeFName

     CCLbutton <- gbutton(" Change Core Line ", handler=function(h,...){
                               CLmainwin <- gwindow("CORELINE SELECTION", parent=MCWindow, visible=FALSE)
                               size(CLmainwin) <- c(200, 100)
                               CLgroup <- ggroup(label="", container=CLmainwin)

                               CLframe <- gframe(text="SELECT THE CORELINE", spacing=5, container=CLgroup)
                               SpectList <- XPSSpectList(activeFName)
                               CLobj <- gradio(SpectList, selected=Indx, horizontal=FALSE, handler=function(h,...){
                                    XPS.CL <- svalue(CLobj)
                                    XPS.CL <- unlist(strsplit(XPS.CL, "\\."))   #drop "NUMber." in component name
                                    Indx <- as.integer(XPS.CL[1])
                                    SpectName <- XPS.CL[2]
                                    assign("activeSpectName", SpectName,envir=.GlobalEnv) #set activeSpectName == last selected spectrum
                                    assign("activeSpectIndx", Indx,envir=.GlobalEnv) #set the activeIndex == last selected spectrum
                                    LoadCoreLine()
                                    dispose(CLmainwin)
                               }, container=CLframe)
                               visible(CLmainwin)<-TRUE
                           }, container = SelectGroup)

     MCFrame2 <- gframe(text = " COMPONENTS ", spacing=5, container = SelectGroup)
     FitComplyt <- glayout(spacing=5, container=MCFrame2)
     if (length(ComponentList) > 1){    #gradio works with at least 2 items
         FitComp <- gradio(ComponentList, selected=1, handler = function(h,...){
                               refresh <<- TRUE    #cancel previous selections
                               replot()   #plot spectrum without marker
                               FComp <- svalue(FitComp)
                               FComp <- as.numeric(unlist(strsplit(FComp, split="C")))   #index selected component
                               FComp <- FComp[2]
                               xx <- Object@Components[[FComp]]@param[2,1] #component position mu
                               Rng <- range(Object@RegionToFit[[1]])
                               if (xx < Rng[1]) {xx <- Rng[1]}
                               if (xx > Rng[2]) {xx <- Rng[2]}
                               yy <- Object@Components[[FComp]]@param[1,1] #component height h
                               FuncName <- Object@Components[[FComp]]@funcName
                               yy <- yy/GetHvalue(Object,FComp, FuncName, 1)
                               Estep <- abs(Object@RegionToFit[[1]][1]-Object@RegionToFit[[1]][2])
                               Xindx <- which(Object@RegionToFit[[1]]>xx-Estep/2 & Object@RegionToFit[[1]] < xx+Estep/2)
                               yy <- yy+Object@Baseline$y[Xindx]  #spectral intensity + baseline at point xx
                               coords[1] <<- xx
                               coords[2] <<- yy
                               refresh <<- FALSE #now plot also the component marker
                               replot()   #replot spectrum and marker
                           },  container = MCFrame2)
               LL <- length(ComponentList)
               NCol <- ceiling(LL/5)   #gradio will be split in solumns of 5 elements
               for(ii in 1:LL){
                   tkpack.forget(FitComp$widgets[[ii]])  # unparent widgets (uses library call)
               }
               for(kk in 1:NCol){
                   NN <- (kk-1)*5
                   for (ii in 1:5) {
                        if((ii+NN) > LL) {break}
                        FitComplyt[ii,kk] <- FitComp$widgets[[(ii+NN)]]
                   }
               }
     } else {
         FitComp <- gcheckboxgroup(ComponentList, checked=TRUE, handler = function(h,...){
                               refresh<<-TRUE    #cancel previous component markers
                               replot()   #plot spectrum only
                               FComp <- svalue(FitComp)
                               FComp <- as.numeric(unlist(strsplit(FComp, split="C")))   #index selected compoent
                               FComp <- FComp[2]
                               xx <- Object@Components[[FComp]]@param[2,1] #component position mu
                               Rng <- range(Object@RegionToFit[[1]])
                               if (xx < Rng[1]) {xx <- Rng[1]}
                               if (xx > Rng[2]) {xx <- Rng[2]}
                               Estep <- abs(Object@RegionToFit[[1]][1]-Object@RegionToFit[[1]][2])
                               Xindx <- which(Object@RegionToFit[[1]] > xx-Estep/2 & Object@RegionToFit[[1]] < xx+Estep/2)
                               yy <- yy+Object@Baseline$y[Xindx]  #spectral intensity + Baseline at point xx
                               coords[1] <<- xx
                               coords[2] <<- yy
                               refresh <<- FALSE #now plot spectrum + component marker
                               replot()
                           },  container = MCFrame2)
     }

     MCFrame3 <- gframe(text = " OPTIONS ", spacing=5, container = SelectGroup)
     Buttlyt <- glayout(spacing=5, container=MCFrame3)

     Buttlyt[1,1] <- LBFitbutton <- gbutton("          FIT Lev.Marq.        ", handler=function(h,...){
                               FComp <- svalue(FitComp)
                               Object <<- XPSFitLM(Object, plt=FALSE)
                               xx <- Object@Components[[FComp]]@param[2,1] #component position mu
                               yy <- Object@Components[[FComp]]@param[1,1] #component height h
                               Estep <- abs(Object@RegionToFit[[1]][1]-Object@RegionToFit[[1]][2])
                               Xindx <- which(Object@RegionToFit[[1]] > xx-Estep/2 & Object@RegionToFit[[1]] < xx+Estep/2)
                               yy <- yy+Object@Baseline$y[Xindx]  #spectral intensity + Baseline at point xx
                               coords[1] <<- xx #coords of marker of the first fit component
                               coords[2] <<- yy
                               Object <<- sortComponents(Object)
                               refresh <<- FALSE  #now plot also the component marker
                               replot()
                           }, container = Buttlyt)

     Buttlyt[1,2] <- MFFitbutton <- gbutton("           FIT Modfit          ", handler=function(h,...){
                               FComp <- svalue(FitComp)
                               Object <<- XPSModFit(Object, plt=FALSE)
                               xx <- Object@Components[[FComp]]@param[2,1] #component position mu
                               yy <- Object@Components[[FComp]]@param[1,1] #component height h
                               Estep <- abs(Object@RegionToFit[[1]][1]-Object@RegionToFit[[1]][2])
                               Xindx <- which(Object@RegionToFit[[1]] > xx-Estep/2 & Object@RegionToFit[[1]] < xx+Estep/2)
                               yy <- yy+Object@Baseline$y[Xindx]  #spectral intensity + Baseline at point xx
                               coords[1] <<- xx
                               coords[2] <<- yy
                               Object <<- sortComponents(Object)
                               refresh <<- FALSE #now plot also the component marker
                               replot()
                           }, container = Buttlyt)

     Buttlyt[2,1] <- ZRbutton<-gbutton("         SET ZOOM REGION       ", handler = function(h, ...){
                               CompCoords <<- coords   #save the of position component_marker
                               point.coords <<- NULL   #point.coords contain the X, Y data ranges
                               enabled(CCLbutton) <- FALSE
                               enabled(LBFitbutton) <- FALSE
                               enabled(MFFitbutton) <- FALSE
                               enabled(RSTbutton) <- FALSE
                               row1<-" => Right Clicks to define the 2 ZOOM REGION CORNERS (opposite in diagonal)"
                               row2<-"\n => Right click near corners to adjust Zoom Region Dimensions"
                               row3<-"\n => When Zoom Region OK, press MAKE ZOOM"
                               msg<-paste(row1, row2, row3, sep="")
                               gmessage( msg, icon="warning")
                               SetZoom <<- TRUE
                               tkconfigure(Graph, cursor = "tcross")
                           }, container = Buttlyt)

     Buttlyt[2,2] <- MZbutton<-gbutton("           MAKE ZOOM           ", handler = function(h, ...){
                               if (Object@Flags[1]) { #Binding energy set
                                  point.coords$x <- sort(point.coords$x, decreasing=TRUE) #pos$x in decreasing order
                                  point.coords$x[1] <- point.coords$x[1]
                                  point.coords$x[2] <- point.coords$x[2]
                               } else {
                                  point.coords$x<-sort(point.coords$x, decreasing=FALSE) #pos$x in increasing order
                                  point.coords$x[1] <- point.coords$x[1]
                                  point.coords$x[2] <- point.coords$x[2]
                               }
                               Xlimits <<- point.coords$x
                               Ylimits <<- sort(point.coords$y, decreasing=FALSE)
 	                             slot(Object,"Boundaries") <<- point.coords
 	                             point.index <<- 1
 	                             coords <<- CompCoords #restore of position component_marker
 	                             refresh <<- FALSE
                               SetZoom <<- FALSE
                               tkconfigure(Graph, cursor = "crosshair")
                               replot()
                               enabled(CCLbutton) <- TRUE
                               enabled(LBFitbutton) <- TRUE
                               enabled(MFFitbutton) <- TRUE
                               enabled(RSTbutton) <- TRUE
                          }, container = Buttlyt)

     Buttlyt[3,1] <- RSTbutton<-gbutton("           RESET PLOT          ", handler = function(h, ...) {
                               SetZoom <<- FALSE
                               refresh <<- FALSE
  	                            point.index <<- 1
  	                            reset.plot()
                          }, container = Buttlyt)


     Buttlyt[3,2] <- gbutton("          RE-LOAD DATA         ", handler=function(h,...){
                               LoadCoreLine()
                          }, container = Buttlyt)

     Buttlyt[4,1] <- gbutton("              SAVE             ", handler=function(h,...){
#    With button SAVE the Component parameters are updated and are now available for FiTConstraints
                               Indx <- get("activeSpectIndx", envir=.GlobalEnv)
                               XPSSample[[Indx]] <<- Object
                               assign(activeFName, XPSSample, envir = .GlobalEnv)
                               replot()
                               XPSSaveRetrieveBkp("save")
                          }, container = Buttlyt)

     Buttlyt[4,2] <- gbutton("               EXIT            ", handler=function(h,...){
                               dispose(MCWindow)
                               tkdestroy(TKwinGraph)
                               XPSSaveRetrieveBkp("save")
                               plot.new()
                               plot(XPSSample[[Indx]])
                          }, container = Buttlyt)

     LabBox <- ggroup(spacing=1, horizontal=TRUE, container=SelectGroup)
     glabel("Graphical Window size: ", container=LabBox)
     WSvalue <- glabel(text=as.character(hscale), container=LabBox)
     WSize <- gslider(from = 1, to = 2.5, by = 0.1, value = hscale, horizontal=TRUE, handler=function(h,...){
                               WinSize <- svalue(WSize)
                               if (Sys.info()["sysname"] == "Windows") { #sets the interactive TK graphical window dimensions
                                   WinScale <- 1.3*WinSize
                               } else if (exists("X11", env=.GlobalEnv)) {
                                   WinScale <- 0.75*WinSize
                               } else {
                                   stop("tkrplot only supports Windows and X11")
                               }
                               Graph$hscale <<- WinScale
                               Graph$vscale <<- WinScale
                               delete(LabBox, WSvalue)
                               WSvalue <<- glabel(text=as.character(WinSize), container=LabBox)
                               refresh <<- FALSE #now plot also the component marker
                               replot(Object)
                               tcl("update", "idletasks") #needed to get the real width and heigth of TKwinGraph
                               hh <<- as.numeric(tkwinfo("reqheight", TKwinGraph)) #store the height of the TKRplot window
                               ww <<- as.numeric(tkwinfo("reqwidth", TKwinGraph))  #store the width of the TKRplot window
                          }, container=SelectGroup)

     StatBar <- gstatusbar("status", container = MCWindow)



     enabled(CCLbutton) <- TRUE
     enabled(LBFitbutton) <- TRUE
     enabled(MFFitbutton) <- TRUE
     enabled(RSTbutton) <- TRUE
     visible(MCWindow) <- TRUE

#--- TKRPLOT interactive graphics------
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
     assign("TKwinGraph", TKwinGraph, envir=.GlobalEnv)
     tkwm.geometry(TKwinGraph,"+565+20")   #+xxx+yyy indicate the pixels from screen TopLeftCorner(screen) to TopLeftCorner(TopLevel)
     tktitle(TKwinGraph) <- paste("TKplot: ", activeFName)
     TKfrmGraph <- tkframe(TKwinGraph, borderwidth = 0, relief = "groove")  #raised, sunken, flat, ridge, solid, and groove
     tkgrid(TKfrmGraph, padx=0, pady=0)
     Graph <- tkrplot(TKfrmGraph, draw.plot, hscale=WinScale, vscale=WinScale)
     tkgrid(Graph, padx=0, pady=0) #load the Graph image into the frame
     assign("Graph", Graph, envir=.GlobalEnv)
     tcl("update", "idletasks") #needed to get the real width and heigth of TKwinGraph
     hh <- as.numeric(tkwinfo("reqheight", Graph)) #store the height of the TKRplot window
     ww <- as.numeric(tkwinfo("reqwidth", Graph))  #store the width of the TKRplot window
     tkconfigure(Graph, bg="white") #set to white the tkrplot background = grey
     tkconfigure(Graph, cursor = "crosshair")
     tkbind(Graph, "<Button-1>", LBmousedown)   #left mouse button
     tkbind(Graph, "<Button-3>", RBmousedown)   #right mouse button
     tkfocus(Graph)
     assign("TKwinGraph", TKwinGraph, envir=.GlobalEnv) #to operate together with XPS FitConstraints
     assign("Graph", Graph, envir=.GlobalEnv)


#--- Marker-----
     if (NoFit==FALSE){
        coords[1] <- Object@Components[[1]]@param[2,1] #component position mu
        coords[2] <- Object@Components[[1]]@param[1,1] #component1 height h
        FuncName <- Object@Components[[1]]@funcName
        coords[2] <- coords[2]/GetHvalue(Object,1, FuncName, 1)
        Estep <- abs(Object@RegionToFit[[1]][1]-Object@RegionToFit[[1]][2])
        Xindx <- which(Object@RegionToFit[[1]] > coords[1]-Estep/2 & Object@RegionToFit[[1]] < coords[1]+Estep/2) #indice del vettore X corrispondente alla posizione della componente
        coords[2] <- coords[2]+Object@Baseline$y[Xindx]
        refresh <- FALSE
        replot()
        refresh <- TRUE
     }
}





