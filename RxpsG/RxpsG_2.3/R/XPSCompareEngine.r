## XPScomparePlot: engine to compare corelines in multi-panel mode

#'XPScompare Engine is the software powering the XPSCompareGUI.
#'This function compares corelines in multi-panel mode following
#'options selected in XPScompareGUI
#'
#'@param PlotParameters the plot parameters asociated to the XPSOverlayGUI options;
#'@param Plot_Args list of plot options;
#'@param SelectedNames list containing the XPSSample names and the Corelines to be plotted;
#'@param Xlim Xrange of the data to be plotted;
#'@param Ylim Yrange of the data to be plotted;
#'@return Returns c(Xlim, Ylim) eventually modified
#'
#'@export
#'

XPScompEngine <-  function(PlotParameters, Plot_Args, SelectedNames, Xlim, Ylim) {

#---  SetPltArgs sets the Plot_Arg list following selections in OverlayGUI
   SetPltArgs <- function(LType,SType , palette) {
            Ylength <- lapply(Y, sapply, length)
            N.XS <- length(SelectedNames$XPSSample)
            N.CL <- length(SelectedNames$CoreLines)
            idx <- 1
            cx <<- list()
            levx <<- list()
#now extract data: for ii, for jj in agreement with the compareXPSSample structure
#see the reading data section
#if C1s, O1s of XPSSamples X1, X2, X3 are compared then compareXPSSample is organized as follow:
#
#   C1s(X1) C1s(X2) C1s(X3) O1s(X1) O1s(X2) O1s(X3)
#
#Then the first for ii runs on the corelines, the second for jj runs on the XPSSamples

            for (ii in 1:N.CL) {        #ii runs on the CoreLines of the XPSSamples
                tmp1 <- NULL
                tmp2 <- NULL
                Xrng <- list()
                Yrng <- list()
                Cex <- Plot_Args$cex
 	              for (jj in 1:N.XS) {    #jj runs on the fit components Corelines
                   if (attr(Ylength[[idx]][1], "names") == "MAIN"){   #holds when just the spectrum is plotted
                       Plot_Args$col[idx] <<- palette[jj]
                       Plot_Args$lty[idx] <<- LType[jj]
                       Plot_Args$pch[idx] <<- SType[jj]
                       Plot_Args$cex[idx] <<- Plot_Args$cex
                       Plot_Args$auto.key$col[idx] <<- palette[jj]
                   }
                   tmp1 <- c( tmp1, rep(ii, times=as.integer(Ylength[[idx]][1])) ) #tmp1 contains indexes associated to the Coreline (group of spectra)
                   tmp2 <- c( tmp2, rep(idx, times=as.integer(Ylength[[idx]][1])) ) #tmp0 contains indexes associated to the XPSSamples... (distinguish different spectra inside the group)
                   Xrng[[jj]] <- range(X[[idx]]) #range computed on group C1s ( O1s ) of different XPSSamples
                   Yrng[[jj]] <- range(Y[[idx]]) #range computed on group C1s ( O1s ) of different XPSSamples
                   idx <- idx+1
                }
                levx[[ii]] <<- tmp1 #required to distinguish multiple panels
                cx[[ii]] <<- tmp2   #required to distinguish different curves
                Xlimits[[ii]] <<- range(Xrng)  #build a list of xlimits  for groups of C1s, O1s, ... CoreLines
                Ylimits[[ii]] <<- range(Yrng)  #build a list of xlimits  for groups of C1s, O1s, ... CoreLines
                if (length(PlotParameters$CustomXY) > 0 && ii==PlotParameters$CustomXY[1]){ #PlotParameters$CustomXY[1] indicated the CLine with custom XY scale
                    Xlimits[[ii]] <<- c(PlotParameters$CustomXY[2], PlotParameters$CustomXY[3]) #customX range
                    Ylimits[[ii]] <<- c(PlotParameters$CustomXY[4], PlotParameters$CustomXY[5]) #custom Yrange
                }
  	         }
   }

#---  rescale a vector so that it lies between the specified minimum and maximum
     rescale <- function(x, newrange=c(0,1)) {
	    if (!is.numeric(x) || !is.numeric(newrange)){
	        stop("Please supply numerics for the x and the new scale")
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


#----- Rading data  -----
#--- compareXPSSample will contain the selected XPSSamples and the selected CoreLines to compare

    compareXPSSample <- new("XPSSample")
    N.XS <- length(SelectedNames$XPSSample)
    N.CL <- length(SelectedNames$CoreLines)
    SpectLengths <- NULL
    idx <- 1
    select <- list()
    if (PlotParameters$OverlayType == "Compare.CoreLines") { #Compare CoreLines in multi.spectrum mode
       for(jj in 1:N.CL){
          SpectName <- SelectedNames$CoreLines[jj] #load all the selected corelines always when SelectedNames$XPSSample == '-----'
          for(ii in 1:N.XS){
             FName <- SelectedNames$XPSSample[ii]
             FName <- get(FName, envir=.GlobalEnv)
             FName[[SpectName]]@.Data[[2]] <- FName[[SpectName]]@.Data[[2]]*SelectedNames$Ampli[ii] #if an amplification factor != 1 was selected
             compareXPSSample[[idx]] <- new("XPSCoreLine")
             compareXPSSample[[idx]]@.Data[1] <- FName[[SpectName]]@.Data[1]  #SpectName not SpectIdx: is possible that Corelines are acquired in different order in different XPSSamples
             compareXPSSample[[idx]]@.Data[2] <- FName[[SpectName]]@.Data[2]  #SpectName not SpectIdx: is possible that Corelines are acquired in different order in different XPSSamples
             compareXPSSample[[idx]]@Flags <- FName[[SpectName]]@Flags
             names(compareXPSSample)[idx] <- SpectName
             select[[idx]] <- "MAIN"
             idx <- idx+1
          }
       }
    } else {
       return()
    }
    NSpect <- idx-1

#----- OverlayEngine FUNCTIONS -----


# set Titles and axis labels
    SpectIndx <- get("activeSpectIndx", envir=.GlobalEnv)     #carico nome spettro attivo
    SpectName <- get("activeSpectName", envir=.GlobalEnv)     #carico nome spettro attivo
    if (length(Plot_Args$xlab$label) == 0) Plot_Args$xlab$label <- compareXPSSample[[1]]@units[1]
    if (length(Plot_Args$ylab$label) == 0) Plot_Args$ylab$label <- compareXPSSample[[1]]@units[2]

#--- Now transform XPSSample into a list
#--- The asList function allows including/skipping fit components, baseline etc. following the select options
#--- NOTE USE of sapply instead of lapply!!!

    XPSSampLen <- length(compareXPSSample)
    XPSSampNames <- names(compareXPSSample)
    Nspettri <- length(XPSSampNames)
    X <- NULL
    Y <- NULL
    for (ii in 1:NSpect){
        tmp <- as.matrix(asList(compareXPSSample[[ii]], select = select[[ii]]))
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
	      wdth <- Plot_Args$xlim[2]-Plot_Args$xlim[1]
	      Plot_Args$xlim[1] <- Plot_Args$xlim[1]-wdth/15
	      Plot_Args$xlim[2] <- Plot_Args$xlim[2]+wdth/15
       Xlim <- Plot_Args$xlim    #Xlim iniziali senza alcuna operazione: vanno mantenute
    }
    if (PlotParameters$Reverse) Plot_Args$xlim <- sort(Plot_Args$xlim, decreasing=TRUE)

#--- Switch BE/KE scale
    if (PlotParameters$SwitchE) {
       XEnergy <- get("XPSSettings", envir=.GlobalEnv)$General[5] #the fifth element of the first column of XPSSettings
       XEnergy <- as.numeric(XEnergy)
       Plot_Args$xlim <- XEnergy-Plot_Args$xlim  #Transform X BE limits in X KE limits
       for (idx in 1:XPSSampLen) {
	        X[[idx]] <- lapply(X[[idx]], function(z, XEnergy){ XEnergy-z }, XEnergy) #Binding to Kinetic energy abscissa
       }
       if (FName[[SpectName]]@Flags[1]==TRUE){ #The original spectra has BE scale
          Plot_Args$xlab$label <- "Kinetic Energy [eV]"
       } else if (FName[[SpectName]]@Flags[1]==FALSE){ #The original spectra has KE scale
          Plot_Args$xlab$label <- "Binding Energy [eV]"
       }
    }

#--- Here Y alignment
    if (PlotParameters$Aligne) {
       LL <- length(Y)
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
    if (is.null(Plot_Args$ylim)) {  #No y limit fixed for zoom
       	Plot_Args$ylim <- range(sapply(Y, sapply, range))
	       wdth <- Plot_Args$ylim[2]-Plot_Args$ylim[1]
	       Plot_Args$ylim[1] <- Plot_Args$ylim[1]-wdth/15
	       Plot_Args$ylim[2] <- Plot_Args$ylim[2]+wdth/15
    }
    Ylim <- Plot_Args$ylim

#------- APPLY GRAPHIC OPTION TO PLOTTING XYplot() ARGS -----------------
    Ylength <- lapply(Y, sapply, length)
    cx <- list()
	   levx <- list()
	   panel <- sapply(Ylength, sum)
    PanelTitles <- NULL
    Xlimits <- list() # costruisco una lista di limiti nel caso Xaxis invertito revers=TRUE
    Ylimits <- list() # costruisco una lista di limiti nel caso Xaxis invertito revers=TRUE

    if ( Plot_Args$type=="l") { #lines are selected for plot
          Plot_Args$auto.key$lines <- TRUE
          Plot_Args$auto.key$points <- FALSE
          if (length(PlotParameters$Colors)==1) {   # B/W LINES
              LType <- Plot_Args$lty                # "solid", "dashed", "dotted" ....
              SType <- rep(NA, 20)
              palette <-  rep("black", 20)          # "Black","black","black",....
              SetPltArgs(LType, SType, palette)
          } else if (length(PlotParameters$Colors) > 1) {   # RainBow LINES
              LType <- rep(Plot_Args$lty[1], 20)    # "solid", "solid", "solid", ....
              SType <- rep(NA, 20)
              palette <- PlotParameters$Colors      #"black", "red", "green"....
              SetPltArgs(LType, SType, palette)
          }
    } else if (Plot_Args$type=="p") { #symbols are selected for plot
          Plot_Args$auto.key$lines <- FALSE
          Plot_Args$auto.key$points <- TRUE
          if (length(PlotParameters$Colors)==1) {   # B/W  SYMBOLS
              LType <- rep(NA, 20)
              SType <- Plot_Args$pch                # VoidCircle", "VoidSquare", "VoidTriangleUp" ....
              palette <-  rep("black", 20)          # "black","black","black",....
              SetPltArgs(LType, SType, palette)
          } else if (length(PlotParameters$Colors) > 1) {   # RainBow SYMBOLS
              LType <- rep(NA, 20)
              SType <- rep(Plot_Args$pch[1], 20)    # "VoidCircle", "VoidCircle", "VoidCircle", ....
              palette <- PlotParameters$Colors      # "black", "red", "green"....
              SetPltArgs(LType, SType, palette)
          }
    } else if (Plot_Args$type=="b") { #Lines + symbols are selected for plot
          Plot_Args$auto.key$lines <- TRUE
          Plot_Args$auto.key$points <- TRUE
          if (length(PlotParameters$Colors)==1) {   # B/W LINES & SYMBOLS
              LType <- Plot_Args$lty                # "solid", "dashed", "dotted" ....
              SType <- Plot_Args$pch                # "VoidCircle", "VoidSquare", "VoidTriangleUp" ....
              palette <-  rep("black", 20)          # "black","black","black",....
              SetPltArgs(LType, SType, palette)
          } else if (length(PlotParameters$Colors) > 1) {   # RainBow LINES & SYMBOLS
              LType <- rep(Plot_Args$lty[1], 20)    #"solid", "solid", "solid", ....
              SType <- rep(Plot_Args$pch[1], 20)    # "VoidCircle", "VoidCircle", "VoidCircle", ....
              palette <- PlotParameters$Colors      #"black", "red", "green"....
              SetPltArgs(LType, SType, palette)
          }
    }

##--- MULTI PANEL---

    if (PlotParameters$OverlayMode=="Multi-Panel") {
       #define row and columns of the panel matrix
       Nspect <- length(SelectedNames$CoreLines)
       Ncol <- 1
       Nrow <- 1
       rr <- FALSE
       while(Nspect > Ncol*Nrow) {
          if (rr){
             Nrow <- Nrow+1
             rr <- FALSE
          } else {
             Ncol <- Ncol+1
             rr <- TRUE
          }
       }

       if (PlotParameters$Reverse) {
          Xlimits <- lapply(Xlimits, sort, decreasing=TRUE) #reverse energy scale
       }

       Plot_Args$xlim <- Xlimits
       Plot_Args$ylim <- Ylimits
       cx <- unlist(cx)
       levx <- unlist(levx)
       df <- data.frame(x = unname(unlist(X)), y = unname(unlist(Y)))

#in df spectra are organized following their category (spect, base, comp, fit)
       PanelTitles <- Plot_Args$PanelTitles  #recover Panel Titles from Plot_Args$PanelTitles. Plot_Args$PanelTitles is a personal argument ignored by xyplot
#      if (length(Plot_Args$main$label) > 0) { PanelTitles <- Plot_Args$main$label } #instead of the default MainLabel uses the title set by the user in OverlayGUI

#in formula y~x is plotted against levx: produces panels showing single XPSSamples
	      Plot_Args$x	<- formula("y ~ x| factor(levx, labels=PanelTitles)")
	      Plot_Args$data	<- df
       Plot_Args$par.settings$strip <- TRUE

	      Plot_Args$groups <- cx
	      Plot_Args$layout <- c(Nrow, Ncol)
       Plot_Args$main	<- NULL
       if (Plot_Args$auto.key[[1]][1] == FALSE) {
          Plot_Args$auto.key <- list()
          Plot_Args$auto.key <- FALSE
       }

	      graph <- do.call(xyplot, args = Plot_Args)
       plot(graph)
    }
    return(c(Xlimits, Ylimits))
}

