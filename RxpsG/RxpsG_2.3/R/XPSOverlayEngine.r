## XPSoverPlot   engine to perform overlay of XPSspectra

#'XPSOverlay Engine is the software powering the XPSOverlayGUI.
#'
#'@param PlotParameters the plot parameters asociated to the XPSOverlayGUI options;
#'@param Plot_Args list of plot options;
#'@param AutoKey_Args list of options for annotation;
#'@param SelectedNames list containing the XPSSample names and the Corelines to be plotted;
#'@param Xlim Xrange of the data to be plotted;
#'@param Ylim Yrange of the data to be plotted;
#'@return Return c(Xlim, Ylim) eventually modified
#'
#'@export
#'

XPSovEngine <-  function(PlotParameters, Plot_Args, AutoKey_Args, SelectedNames, Xlim, Ylim) {

#---  SetPltArgs sets the Plot_Arg list following selections in OverlayGUI
   SetPltArgs <- function(LType,SType , palette, FitStyle) {

        if (PlotParameters$OverlayMode=="Multi-Panel") {
            palette <- rep("black", 20)
        }
        Ylength <- lapply(Y, sapply, length)
        idx <- 1
        cx<<-list()
        levx<<-list()
        for (ii in seq_along(Ylength) ) {               #corro sulle CoreLines dei vari XPSSamples
             tmp1<-NULL
             tmp2<-NULL
             Cex<-Plot_Args$cex
             for ( jj in seq_along(Ylength[[ii]]) ) {    #jj corre sulle componenti Corelines
                   if (attr(Ylength[[ii]][jj], "names") == "MAIN"     #holds when just the spectrum is plotted
                       || attr(Ylength[[ii]][jj], "names") =="RTF"){  #holds when spectrum + Baseline or fit are plotted
                       Plot_Args$col[idx]<<-palette[ii]
                       Plot_Args$lty[idx]<<-LType[ii]
                       Plot_Args$pch[idx]<<-SType[ii]
                       Plot_Args$cex[idx]<<-Plot_Args$cex
                   }
                   if (attr(Ylength[[ii]][jj], "names") == "BASE"){
                       Plot_Args$col[idx]<<-FitStyle$Col[1]
                       Plot_Args$lty[idx]<<-"dashed"
                       Plot_Args$pch[idx]<<-3 #"Cross"
                       Plot_Args$cex[idx]<<-0.3*Plot_Args$cex
                   }
                   if (attr(Ylength[[ii]][jj], "names") == "COMPONENTS" ){
                       Plot_Args$col[idx]<<-FitStyle$Col[2]
                       Plot_Args$lty[idx]<<-FitStyle$Lty[1]
                       Plot_Args$pch[idx]<<-8 #"Star"  2 #"VoidTriangleUp"
                       Plot_Args$cex[idx]<<-0.4*Plot_Args$cex
                   }
                   if (attr(Ylength[[ii]][jj], "names") == "FIT" ){
                       Plot_Args$col[idx]<<-FitStyle$Col[3]
                       Plot_Args$lty[idx]<<-"solid"
                       Plot_Args$pch[idx]<<-8 #"Star" 2 #"VoidTriangleUp"
                       Plot_Args$cex[idx]<<-Plot_Args$cex
                   }
                   tmp1<-c(tmp1, rep(ii, times=as.integer(Ylength[[ii]][jj])))   #indice associato ai vari XPSSample
                   tmp2<-c(tmp2, rep(idx, times=as.integer(Ylength[[ii]][jj])))  #indice associato alle componenti della CoreLine (baseline o componenti del fit)
                   idx<-idx+1
             }
             levx[[ii]]<<-tmp1 #required to distinguish multiple panels
             cx[[ii]]<<-tmp2   #required to distinguish different curves
             Xlimits[[ii]]<<-rev(range(X[[ii]]))  #costruisco una lista di xlim invertiti
  	     }
    }


#---  rescale a vector so that it lies between the specified minimum and maximum
    rescale <- function(x, newrange=c(0,1)) {

	       if (!is.numeric(x) || !is.numeric(newrange)){
	          stop("Must supply numerics for the x and the new scale")
        }
        if (length(newrange) != 2) {
           stop("newrange must be a numeric vector with 2 elements")
        }
        newmin <- min(newrange)
        newmax <- max(newrange)
        oldrange <- range(x)
        if (oldrange[1] == oldrange[2]) {
           if (newmin==0) {
              return(x-oldrange[1])
           } else {
             warning("The supplied vector is a constant. Cannot rescale")
             return(x)
           }
	       } else {
	          ratio <- (newmax - newmin) / (oldrange[2] - oldrange[1])
	          return( newmin + (x - oldrange[1]) * ratio )
        }
    }

#----- OverlayEngine FUNCTIONS -----

#--- Define a NEW XPSSample containing all the spectra to be plotted with related fits (if present)
    if (length(SelectedNames$XPSSample) == 0) { return() } #no selected XPSSamples to plot
    overlayXPSSample <- new("XPSSample")
    LL<-length(SelectedNames$XPSSample)

    SpectLengths<-NULL
    idx<-1
    select<-list()
    for(ii in 1:LL){
        if (length(SelectedNames$XPSSample[ii]) > 0 && SelectedNames$XPSSample[ii] != "-----") { #Set XPSSample==FName only if SelectedNames$XPSSample != -----
            FName<-SelectedNames$XPSSample[ii]  #this is a string
            FName<-get(FName, envir=.GlobalEnv) #this is an XPSSample datafile
        }
        SpectName<-SelectedNames$CoreLines[ii]  #Load all the selected  indicated in SelectedNames$XPSSample
        SpectName<-unlist(strsplit(SpectName, "\\."))   #skip the initial number at beginning of the CoreLine name
        SpectIdx<-as.numeric(SpectName[1])
        FName[[SpectIdx]]@.Data[[2]]<-FName[[SpectIdx]]@.Data[[2]]*SelectedNames$Ampli[ii]    #if an amplification factor was seleted, AmpliFact != 1

        SpectLengths[ii]<-length(FName[[SpectIdx]]@.Data[[2]]) #it is needed for 3Dplot
        overlayXPSSample[[idx]] <- FName[[SpectIdx]]
        SpectName <- FName[[SpectIdx]]@Symbol
        names(overlayXPSSample)[idx] <- SpectName

#--- selection of corelines format: simple spectrum, spectrum+Baseline, spectrum+CompleteFit
        if (PlotParameters$OverlayType == "Spectrum") select[[idx]] <- "MAIN"
        if (PlotParameters$OverlayType == "Spectrum+Baseline") {  # CTRL if baseline and fit components are present
           if (length(FName[[SpectIdx]]@RegionToFit) > 0) {       # if not only baseline or main spectrum are ploted
              select[[idx]] <- c("RTF", "BASE")
           } else {
              select[[idx]] <- "MAIN"
           }
        }

        if (PlotParameters$OverlayType == "Spectrum+Fit") {
           if (length(FName[[SpectIdx]]@Components) > 0) {
              select[[idx]] <- c("RTF", "BASE", "COMPONENTS", "FIT")
           } else if (length(FName[[SpectIdx]]@RegionToFit) > 0) {
              select[[idx]] <- c("RTF", "BASE")
           } else {
              select[[idx]] <- "MAIN"
           }
        }
        idx <- idx+1
    }
    NXPSSamp <- idx-1

# set Titles and axis labels
    FName<-get(activeFName, envir=.GlobalEnv)               #load the active XPSSample
    SpectIndx<-get("activeSpectIndx", envir=.GlobalEnv)     #load the active spectrum index
    SpectName<-get("activeSpectName", envir=.GlobalEnv)     #load the active spectrum name
    if (length(Plot_Args$xlab$label) == 0) Plot_Args$xlab$label<-FName[[SpectName]]@units[1] #set the axis labels if not defined
    if (length(Plot_Args$ylab$label) == 0) Plot_Args$ylab$label<-FName[[SpectName]]@units[2]

#--- Now transform XPSSample into a list
#--- The asList function allows including/skipping fit components, baseline etc. following the select options
#--- NOTE USE of sapply instead of lapply!!!
    XPSSampLen <- length(overlayXPSSample)
    XPSSampNames <- names(overlayXPSSample)
    Nspettri<-length(XPSSampNames)
    X<-NULL
    Y<-NULL
    for (ii in 1:NXPSSamp){
        tmp <- as.matrix(asList(overlayXPSSample[[ii]], select = select[[ii]]))
        X <- c(X, tmp["x", ])   #X coords of the selected spectra
        Y <- c(Y, tmp["y", ])   #Y coords of the selected spectra
    }
    Xlim0 <- range(sapply(X, sapply, range))
    Ylim0 <- range(sapply(Y, sapply, range))

#--- Now all the data manipulation options
#--- X offset
    if ( ! is.null(PlotParameters$XOffset) ) {
	      if ( length(PlotParameters$XOffset)!=XPSSampLen ) {
           offset_sequence <- seq(from = 0, by = PlotParameters$XOffset, length.out = XPSSampLen)  #costruisco un vettore con gli Xoffset necessari a Xshiftare i vari spettri
       } else {
           offset_sequence <- PlotParameters$XOffset
       }
	      for (idx in 1:XPSSampLen) {
	          X[[idx]] <- lapply(X[[idx]], "+", offset_sequence[idx])
	      }
    }
#--- set xlim, and reverseX if BE
    if (is.null(Plot_Args$xlim)) {  #non ho fissato xlim per fare lo zoom
	      Plot_Args$xlim <- range(sapply(X, sapply, range))
	      wdth<-Plot_Args$xlim[2]-Plot_Args$xlim[1]
	      Plot_Args$xlim[1] <- Plot_Args$xlim[1]-wdth/15
	      Plot_Args$xlim[2] <- Plot_Args$xlim[2]+wdth/15
       Xlim <- Plot_Args$xlim    #Xlim iniziali senza alcuna operazione: vanno mantenute
    }
    if (PlotParameters$Reverse) Plot_Args$xlim <- sort(Plot_Args$xlim, decreasing=TRUE)

#--- Switch BE/KE scale
    if (PlotParameters$SwitchE) {
       XEnergy<-get("XPSSettings", envir=.GlobalEnv)$General[5] #the fifth element of the first column of XPSSettings
       XEnergy<-as.numeric(XEnergy)
       Plot_Args$xlim<-XEnergy-Plot_Args$xlim  #Transform X BE limits in X KE limits
       for (idx in 1:XPSSampLen) {
	          X[[idx]] <- lapply(X[[idx]], function(z, XEnergy){ XEnergy-z }, XEnergy) #Binding to Kinetic energy abscissa
       }
       if (FName[[SpectName]]@Flags[1]==TRUE){ #The original spectra has BE scale
          Plot_Args$xlab$label<-"Kinetic Energy [eV]"
       } else if (FName[[SpectName]]@Flags[1]==FALSE){ #The original spectra has KE scale
          Plot_Args$xlab$label<-"Binding Energy [eV]"
       }
    }

#--- Here Y alignment
    if (PlotParameters$Aligne) {
       LL<-length(Y)
       if ( all(sapply(Y, function(x) !is.na(charmatch("BASE", names(x))))) ) {
			       minybase <- sapply(Y, function(x) min(x$BASE))
			       for (idx in c(1:LL)) {
			          	Y[[idx]] <- lapply(Y[[idx]], "-", minybase[idx])
		       	}
       } else {
          for (idx in c(1:LL)) {
              Y[[idx]] <- lapply(Y[[idx]], function(j) {
	                               return( rescale(j, newrange = c(0, diff(range(j))) ) )
	                         })
          }
       }
    }

#--- Y normalization == scale c(0,1)

    if (PlotParameters$Normalize) {
			    maxy <- sapply(Y, function(x) max(sapply(x, max))) #here Y is the list of XPSSamples with baseline fitComp...
			    for (idx in c(1:XPSSampLen)) {
				       Y[[idx]] <- lapply(Y[[idx]], "/", maxy[idx])
			    }
    }

#--- Y offset
    if ( ! is.null(PlotParameters$YOffset) ) {
	      if ( length(PlotParameters$YOffset)!=XPSSampLen ) {
           offset_sequence <- seq(from = 0, by = PlotParameters$YOffset, length.out = XPSSampLen)  #costruisco un vettore con gli Xoffset necessari a Xshiftare i vari spettri
       } else {
           offset_sequence <- PlotParameters$YOffset
       }
		     for (idx in c(1:XPSSampLen)) {
		         Y[[idx]] <- lapply(Y[[idx]], "+", offset_sequence[idx])
		     }
    }

#--- After processing set Ylim
    if (is.null(Plot_Args$ylim)) {  #non ho fissato ylim per fare lo zoom
       	Plot_Args$ylim <- range(sapply(Y, sapply, range))
	       wdth<-Plot_Args$ylim[2]-Plot_Args$ylim[1]
	       Plot_Args$ylim[1] <- Plot_Args$ylim[1]-wdth/15
	       Plot_Args$ylim[2] <- Plot_Args$ylim[2]+wdth/15
    }
    Ylim <- Plot_Args$ylim

#------- APPLY GRAPHIC OPTION TO PLOTTING XYplot() ARGS -----------------
    Ylength <- lapply(Y, sapply, length)
    cx <- list()
	   levx <- list()
	   FitStyle <- list(Lty=NULL, Col=NULL)
	   panel <- sapply(Ylength, sum)
    PanelTitles <- NULL
    Xlimits<-list() # costruisco una lista di limiti nel caso Xaxis invertito revers=TRUE

    if ( Plot_Args$type=="l") { #lines are selected for plot
         AutoKey_Args$lines<-TRUE
         AutoKey_Args$points<-FALSE
         if (length(PlotParameters$Colors)==1) {   # B/W LINES
             LType <- Plot_Args$lty                # "solid", "dashed", "dotted" ....
             SType <- rep(NA, 20)
             palette <-  rep("black", 20)          # "Black","black","black",....
             FitStyle$Lty <- PlotParameters$CompLty
             FitStyle$Col <- c("black", "grey45", "black")
             SetPltArgs(LType, SType, palette, FitStyle)
         } else if (length(PlotParameters$Colors) > 1) {   # RainBow LINES
             LType <- rep(Plot_Args$lty[1], 20)    # "solid", "solid", "solid", ....
             SType <- rep(NA, 20)
             palette <- PlotParameters$Colors      #"black", "red", "green"....
             FitStyle$Lty <- PlotParameters$CompLty
             FitStyle$Col <- PlotParameters$FitCol         #Default colors("cadetblue", "grey45", "orangered3")
             SetPltArgs(LType, SType, palette, FitStyle)
         }
    } else if (Plot_Args$type=="p") { #symbols are selected for plot
         AutoKey_Args$lines<-FALSE
         AutoKey_Args$points<-TRUE
         if (length(PlotParameters$Colors)==1) {   # B/W  SYMBOLS
             LType <- rep(NA, 20)
             SType <- Plot_Args$pch                # VoidCircle", "VoidSquare", "VoidTriangleUp" ....
             palette <-  rep("black", 20)          # "black","black","black",....
             FitStyle$Lty <- PlotParameters$CompLty
             FitStyle$Col <- c("black", "grey45", "black")
             SetPltArgs(LType, SType, palette, FitStyle)
         } else if (length(PlotParameters$Colors) > 1) {   # RainBow SYMBOLS
             LType <- rep(NA, 20)
             SType <- rep(Plot_Args$pch[1], 20)    # "VoidCircle", "VoidCircle", "VoidCircle", ....
             palette <- PlotParameters$Colors      # "black", "red", "green"....
             FitStyle$Lty <- PlotParameters$CompLty
             FitStyle$Col <- PlotParameters$FitCol #c("cadetblue", "grey45", "orangered3")
             SetPltArgs(LType, SType, palette, FitStyle)
         }
    } else if (Plot_Args$type=="b") { #Lines + symbols are selected for plot
         AutoKey_Args$lines<-TRUE
         AutoKey_Args$points<-TRUE
         if (length(PlotParameters$Colors)==1) {   # B/W LINES & SYMBOLS
             LType <- Plot_Args$lty                # "solid", "dashed", "dotted" ....
             SType <- Plot_Args$pch                # "VoidCircle", "VoidSquare", "VoidTriangleUp" ....
             palette <-  rep("black", 20)          # "black","black","black",....
             FitStyle$Lty <- PlotParameters$CompLty
             FitStyle$Col <- c("black", "grey45", "black")
             SetPltArgs(LType, SType, palette, FitStyle)
         } else if (length(PlotParameters$Colors) > 1) {   # RainBow LINES & SYMBOLS
             LType <- rep(Plot_Args$lty[1], 20)    #"solid", "solid", "solid", ....
             SType <- rep(Plot_Args$pch[1], 20)    # "VoidCircle", "VoidCircle", "VoidCircle", ....
             palette <- PlotParameters$Colors      #"black", "red", "green"....
             FitStyle$Lty <- PlotParameters$CompLty
             FitStyle$Col <- PlotParameters$FitCol #c("cadetblue", "grey45", "orangered3")
             SetPltArgs(LType, SType, palette, FitStyle)
         }
    }

##--- SINGLE PANEL---
    if (PlotParameters$OverlayMode=="Single-Panel") {
       if (length(Plot_Args$main$label) == 0) Plot_Args$main$label<-SpectName
   	   df <- data.frame(x = unname(unlist(X)), y = unname(unlist(Y)) )
	      Plot_Args$x	<- formula("y ~ x")
       Plot_Args$data	<- df
       Plot_Args$groups	<- unlist(cx)
       graph <- do.call(xyplot, args = Plot_Args)
       plot(graph)
    }

##--- MULTI PANEL---
    if (PlotParameters$OverlayMode=="Multi-Panel") {
       #define row and columns of the panel matrix
       Nspettri<-length(XPSSampNames)
       Ncol<-1
       Nrow<-1
       rr<-FALSE
       while(Nspettri>Ncol*Nrow) {
          if (rr){
             Nrow<-Nrow+1
             rr<-FALSE
          } else {
             Ncol<-Ncol+1
             rr<-TRUE
          }
       }

       Plot_Args$xlim<-NULL  #X range is defined inside xyplot
       Plot_Args$ylim<-NULL  #Y range is defined inside xyplot
       if (PlotParameters$Reverse) { #If reverse is TRUE limits must be given through a list containing:
          Plot_Args$xlim<-Xlimits    #Xlimits[[1]]=X1max, X1min
       }                             #Xlimits[[2]]=X2max, X2min
                                    #   ...
       cx<-unlist(cx)
       levx<-unlist(levx)
       df <- data.frame(x = unname(unlist(X)), y = unname(unlist(Y)))

#in DF raggruppo curve secondo la loro categoria (spect, base, comp, fit)
       PanelTitles<-Plot_Args$PanelTitles  #recover Panel Titles from Plot_Args$PanelTitles. Plot_Args$PanelTitles ia a personal argument ignored by xyplot
       if (length(Plot_Args$main$label) > 0) { PanelTitles <- Plot_Args$main$label } #instead of the default MainLabel uses the title set by the user in OverlayGUI
#in formula y~x is plotted against levx: produces panels showing single XPSSamples
	      Plot_Args$x	<- formula("y ~ x| factor(levx, labels=PanelTitles)")
	      Plot_Args$data	<- df
       Plot_Args$par.settings$strip<-TRUE

	      Plot_Args$groups	<- cx
	      Plot_Args$layout	<- c(Nrow, Ncol)
       Plot_Args$main	<- NULL

	      graph <- do.call(xyplot, args = Plot_Args)
       plot(graph)
    }

##--- Pseudo  TreD ---
    if (PlotParameters$XOffset*PlotParameters$YOffset != 0  #both XOffset che YOffset different from zero
        && PlotParameters$OverlayMode=="Single-Panel" ){    #This means no 3D, no Multipanel was set
#if Yoffset is positive
#plot() draws the spectrum LL as the last => LL spectrum is in front. To make the first spectum to be
#in front plot() is applied to a reversed list of spectra: the LL spectrum is the first_plotted the 1 sectrum is the last_plotted
#the effect is the 1th spectrum in from, the LL spectrum on the back
#if Yoffset negative original order is OK
       tmp<-list()
       LL<-length(X)
       for(ii in 1:LL){
          tmp[[ii]]<-X[[(LL-ii+1)]]   #reverse X column order
       }
       X<-tmp
       if (PlotParameters$YOffset > 0){
           tmpLT<-array()
           tmpST<-array()
           tmpPal<-array()
           tmp<-list()
           for(ii in 1:LL){
              tmp[[ii]]<-Y[[(LL-ii+1)]]   #reverse Y column order
              tmpLT[ii]<-LType[(LL-ii+1)]
              tmpST[ii]<-SType[(LL-ii+1)]
              tmpPal[ii]<-palette[(LL-ii+1)]
           }
           Y<-tmp
#group properties cx and levx will be reversed in SetPltArgs
           SetPltArgs(tmpLT, tmpST, tmpPal, FitStyle)

           df <- data.frame(x = unname(unlist(X)), y = unname(unlist(Y)) )
           Plot_Args$data	<- df
	          Plot_Args$x	<- formula("y ~ x")
           Plot_Args$groups	<- unlist(cx)
           graph <- do.call(xyplot, args = Plot_Args)
           plot(graph)
       }
       if (PlotParameters$Pseudo3D==TRUE){
           Xrng<-range(X)
           Yrng<-range(Y)
           if (PlotParameters$XOffset >0) {
              SgmntX1<-c(Xrng[1], Xrng[1]+(Nspettri-1)*PlotParameters$XOffset) #bottom diag. segment1 on the right Xcoord  row1= x1,x2
              SgmntX2<-c(Xrng[1]+(Nspettri-1)*PlotParameters$XOffset, Xrng[2]) #bottom orizz. segment2 in front Xcoord  row2= x3,x4
              SgmntX3<-c(Xrng[1]+(Nspettri-1)*PlotParameters$XOffset, Xrng[1]+(Nspettri-1)*PlotParameters$XOffset) #segment3 vertical on the right Xcoord  row3= x5,x6
              SgmntY1<-c(Yrng[1], Yrng[1]+(Nspettri-1)*PlotParameters$YOffset) #segment1 Ycoord  row1= y1,y2
              SgmntY2<-c(Yrng[1]+(Nspettri-1)*PlotParameters$YOffset, Yrng[1]+(Nspettri-1)*PlotParameters$YOffset) #segment2 Ycoord  row2= y3,y4
              SgmntY3<-c(Yrng[1]+(Nspettri-1)*PlotParameters$YOffset, Yrng[2]) #segment3 Ycoord  row3= y5,y6
           }
           if (PlotParameters$XOffset <0) {
              SgmntX1<-c(Xrng[2], Xrng[2]+(Nspettri-1)*PlotParameters$XOffset) #segment1 diag. on the left Xcoord  row1= x1,x2
              SgmntX2<-c(Xrng[2]+(Nspettri-1)*PlotParameters$XOffset, Xrng[1]) #segment2 orizz. Xcoord  row2= x3,x4
              SgmntX3<-c(Xrng[2]+(Nspettri-1)*PlotParameters$XOffset, Xrng[2]+(Nspettri-1)*PlotParameters$XOffset) #segment3 vertical on the left Xcoord  row3= x5,x6
              SgmntY1<-c(Yrng[1], Yrng[1]+(Nspettri-1)*PlotParameters$YOffset) #segment1 Ycoord  row1= y1,y2
              SgmntY2<-c(Yrng[1]+(Nspettri-1)*PlotParameters$YOffset, Yrng[1]+(Nspettri-1)*PlotParameters$YOffset) #segment2 Ycoord  row2= y3,y4
              SgmntY3<-c(Yrng[1]+(Nspettri-1)*PlotParameters$YOffset, Yrng[2]) #segment3 Ycoord  row3= y5,y6
           }
           Plot_Args$x	<- formula("y ~ x")
           Plot_Args$groups	<- cx <-c(1,1)
           Plot_Args$col[1]<-"gray25"
           Plot_Args$lty[1]<-"dashed"

           Plot_Args$data <- data.frame(x = SgmntX1, y = SgmntY1, groups = factor(cx))
           segm1 <- do.call(xyplot, args = Plot_Args)

           Plot_Args$data <- data.frame(x = SgmntX2, y = SgmntY2, groups = factor(cx))
           segm2 <- do.call(xyplot, args = Plot_Args)
           segm1 <- segm1+as.layer(segm2)      #overlay segment 2 to previous plot stored in graph

           Plot_Args$data <- data.frame(x = SgmntX3, y = SgmntY3, groups = factor(cx))
           segm3 <- do.call(xyplot, args = Plot_Args)
           segm1 <- segm1+as.layer(segm3)      #overlay segment 3 to previous plot stored in graph

           segm1 <- segm1+as.layer(graph)
           plot(segm1)
       }
    }

##--- Real TreD ---
    if (PlotParameters$OverlayMode=="TreD") {
       Z<-NULL
       Cloud_Args<-list()
       for (ii in 1:XPSSampLen){
            Z<-c(Z, rep(ii, SpectLengths[ii]))
       }

   	   df <- data.frame(x =unname(unlist(X)), y=as.vector(Z), z=unname(unlist(Y)) )

#---Set Cloud Style parameters---
       LType <- Plot_Args$lty   # "solid", "dashed", "dotted" ....
       SType <- Plot_Args$pch   # "VoidCircle", "VoidSquare", "VoidTriangleUp" ....
       LW <- Plot_Args$lwd
       CX <- Plot_Args$cex
       if(Plot_Args$type !="l" && Plot_Args$type != "p" && Plot_Args$type != "b") { return() } #if style not defined return!

       if (length(PlotParameters$Colors)==1) { # B/W Lines
          if (Plot_Args$type=="l"){ #lines
             Cloud_Args<-list(lty=LType,cex=CX,lwd=LW,type="l")    #Plot_Args$lty = "solid", "dashed", "dotted" ....
          } else if(Plot_Args$type=="p"){ #symbols
             Cloud_Args<-list(pch=SType,cex=CX,lwd=LW,type="p")
          } else if(Plot_Args$type=="b"){ #lines & symbols
             Cloud_Args<-list(lty="solid", pch=SType,cex=CX,lwd=LW,type="b")
          }
          Cloud_Args$col <- Plot_Args$col <- rep("black", LL)
       } else if (length(PlotParameters$Colors) > 1) {   # Rainbow Lines
          if (Plot_Args$type=="l"){ # lines
             Cloud_Args<-list(lty="solid",cex=CX,lwd=LW,type="l")
          } else if(Plot_Args$type=="p"){  #symbols
             Cloud_Args<-list(pch=1,cex=CX,lwd=LW,type="p")  #pch=1  voidcircle
          } else if(Plot_Args$type=="b"){  # lines & symbols
             Cloud_Args<-list(lty="solid", pch=1,cex=CX,lwd=LW,type="b")
          }
          Cloud_Args$col <- Plot_Args$col <- PlotParameters$Colors
       }

#---axis options---
       if (length(Plot_Args$main$label) == 0) { Plot_Args$main$label <- SpectName }
       if (length(Plot_Args$xlab$label) == 0) { Plot_Args$xlab$label <- FName[[SpectName]]@units[1] }
       if (length(Plot_Args$ylab$label) == 0) { Plot_Args$ylab$label <- "Sample" }
       if (length(Plot_Args$zlab$label) == 0) { Plot_Args$zlab$label <- FName[[SpectName]]@units[2] }
       if (PlotParameters$Normalize == TRUE ) { Plot_Args$zlab$label <- "Intensity [a.u.]" }
       if (PlotParameters$Reverse) { #If reverse==TRUE compute limits and reverse
          Xmax <- max(sapply(Xlimits, max))
          Xmin <- min(sapply(Xlimits, min))
          Cloud_Args$xlim<-c(Xmax,Xmin)  #Xlimits[[1]]=X1max, X1min
       }

       Cloud_Args$ylim<-as.character(c(1:XPSSampLen))  #in Y ho messo i channels Z
       LogOnOff<-Plot_Args$scales$x$log #if x ax log TRUE all axes TRUE
       Cloud_Args$scales<-list(cex=Plot_Args$scales$cex, tck=c(1,0), alternating=c(1), relation="free",
                               x=list(log=LogOnOff), y=list(log=LogOnOff),z=list(log=LogOnOff),
                               arrows=FALSE)

#---3D rendering---
       Cloud_Args$aspect<-as.numeric(PlotParameters$TreDAspect)
       Cloud_Args$screen<-list(x=-60,
                               y= PlotParameters$AzymuthRot,
                               z= PlotParameters$ZenithRot)
       Cloud_Args$main<-list(label=Plot_Args$main$label,cex=Plot_Args$main$cex)
       Cloud_Args$xlab<-list(label=Plot_Args$xlab$label, rot=PlotParameters$AzymuthRot-10, cex=Plot_Args$xlab$cex)
       Cloud_Args$ylab<-list(label=Plot_Args$ylab$label, rot=PlotParameters$AzymuthRot-80, cex=Plot_Args$xlab$cex)
       Cloud_Args$zlab<-list(label=Plot_Args$zlab$label, rot=90, cex=Plot_Args$xlab$cex)
       LL<-length(Y)

#---legend options---
       if (length(Plot_Args$auto.key) > 1) { #auto.key TRUE
          Plot_Args$auto.key$space <- NULL   #Inside top right position
          Plot_Args$auto.key$corner<-c(1,1)
          Plot_Args$auto.key$x<- 0.95
          Plot_Args$auto.key$y<- 0.95

          Cloud_Args$auto.key<-Plot_Args$auto.key
          Cloud_Args$par.settings<-Plot_Args$par.settings
       }
#---plot commands---
	      Cloud_Args$x <- formula("z ~ x*y")
       Cloud_Args$data <- df
       Cloud_Args$groups	<- unlist(cx)

       graph <- do.call(cloud, args=Cloud_Args)
       plot(graph)
    }

    return(c(Xlim, Ylim))
}

