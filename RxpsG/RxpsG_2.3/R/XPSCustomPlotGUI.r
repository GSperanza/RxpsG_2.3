#'Function to generate personalized plot
#'
#'This GUI allows setting the various parameters genereally used for
#'replotting data. Through a user friendly interface it is possible to
#'set colors, lines or symbols, their weight, and annotate the plot.
#'This Gui is based on the Lattice package.
#'No parameters are passed to this function.
#'
#'@seealso \link{Lattice}
#'@examples
#'
#'\dontrun{
#'	XPSCustomPlot()
#'}
#'
#'@export
#'

XPSCustomPlot <- function(){


   CtrlPlot <- function(){

            logTicks <- function (lim, loc = c(1, 5)) {
               ii <- floor(log10(range(lim))) + c(-1, 2)
               main <- 10^(ii[1]:ii[2])
               r <- as.numeric(outer(loc, main, "*"))
               r[lim[1] <= r & r <= lim[2]]
            }

            #Load the selected oiptions
            if ( svalue(T2obj4) == "ON") { # lines ON
               AutoKey_Args$lines <<- TRUE
               enabled(T3obj2) <- TRUE
               enabled(T3obj3) <- TRUE
               enabled(T4obj2) <- TRUE
               enabled(T4obj3) <- TRUE
               enabled(T5obj2) <- TRUE
               enabled(T5obj3) <- TRUE
            }
            if ( svalue(T2obj4) == "OFF") { # lines ON
               AutoKey_Args$lines <<- FALSE
               enabled(T3obj2) <- FALSE
               enabled(T3obj3) <- FALSE
               enabled(T4obj2) <- FALSE
               enabled(T4obj3) <- FALSE
               enabled(T5obj2) <- FALSE
               enabled(T5obj3) <- FALSE
            }
            if ( svalue(T2obj7) == "ON") {  # symbols ON
               AutoKey_Args$points <<- TRUE
               enabled(T3obj4) <- TRUE
               enabled(T3obj5) <- TRUE
               enabled(T4obj4) <- TRUE
               enabled(T4obj5) <- TRUE
               enabled(T5obj4) <- TRUE
               enabled(T5obj5) <- TRUE
            }
            if ( svalue(T2obj7) == "OFF") {  # symbols OFF
               AutoKey_Args$points <<- FALSE
               enabled(T3obj4) <- FALSE
               enabled(T3obj5) <- FALSE
               enabled(T4obj4) <- FALSE
               enabled(T4obj5) <- FALSE
               enabled(T5obj4) <- FALSE
               enabled(T5obj5) <- FALSE
            }

            if ( svalue(T2obj4) == "ON" && svalue(T2obj7) == "ON") {  # symbols OFF
               Plot_Args$type <<- "b"  # both: line and symbols
            }  #conditions on lines and symbols see above (T2obj4==ON   T2obj7==ON)

            idx <- svalue(T1obj7,index=TRUE)

            LL=length(SpectName)
            if (LL > 0){
	              graph <- do.call(xyplot, args = Plot_Args)
	              plot(graph)
            }

            if ( length(svalue(T2obj8)) > 0 ) {  #Request to plot Error bars
               PlotErrorBar()
            }
   }

   PlotErrorBar <- function(){

            ErrOptions <- as.numeric(svalue(T2obj8, index=TRUE))
            NOpt <- length(ErrOptions)
            xx <- FName[[SpectIndx]]@.Data[[1]]
            yy <- FName[[SpectIndx]]@.Data[[2]]
            Err <- FName[[SpectIndx]]@.Data[[4]]

            NData <- length(xx)
            EndAmpli <- as.numeric(svalue(T2obj9)) #Amplitude of the limiting ending bars
            EndLngth <- EndAmpli*(range(xx)[2]-range(xx)[1])/50
            colr <- svalue(T2obj1)

            MaxErr <- max(Err)
            MinErr <- min(Err)
            Plot_Args$ylim <<- c(Ylim[1]-MinErr, Ylim[2]+MaxErr)  #include error bars in the Y plot scale
            graph <- do.call(xyplot, args = Plot_Args)
            plot(graph)

            yy0 <- NULL
            yy1 <- NULL
            trellis.focus("panel", 1, 1, clip.off=TRUE, highlight=FALSE)  #to enable modify the lattice graph
            for(ii in 1:NOpt){
                if (ErrOptions[ii] == 1){ #plot the upper part of the stanbdard error
                    for(jj in 1:NData){
                        yy0[jj] <- yy[jj]
                        yy1[jj] <- yy[jj]+Err[jj]
                        if (svalue(T1obj7,index=TRUE) == 3) { # Y-log scale selected
                            yy0[jj] <- log10(yy[jj])
                            yy1[jj] <- log10(yy1[jj])
                        }
                        lsegments(x0=xx[jj], y0=yy0[jj], x1=xx[jj], y1=yy1[jj], col=colr) #plots the upper part of the segment
                    }
                }
                if (ErrOptions[ii] == 2){ #plot the lower part of the stanbdard error
                    for(jj in 1:NData){
                        yy0[jj] <- yy[jj]
                        yy1[jj] <- yy[jj]-Err[jj]
                        if (svalue(T1obj7,index=TRUE) == 3) { # Y-log scale selected
                            yy0[jj] <- log10(yy[jj])
                            yy1[jj] <- log10(yy1[jj])
                        }
                        lsegments(x0=xx[jj], y0=yy0[jj], x1=xx[jj], y1=yy1[jj], col=colr) #plots the upper part of the segment
                    }
                }
                if (ErrOptions[ii] == 3){ #plot the upper limiting bar
                    for(jj in 1:NData){
                        yy1[jj] <- yy[jj]+Err[jj]
                        if (svalue(T1obj7,index=TRUE) == 3) { # Y-log scale selected
                            yy1[jj] <- log10(yy1[jj])
                        }
                        lsegments(x0=xx[jj]-EndLngth, y0=yy1[jj], x1=xx[jj]+EndLngth, y1=yy1[jj], col=colr) #plots the upper part of the segment
                    }
                }
                if (ErrOptions[ii] == 4){ #plot the lower limiting bar
                    for(jj in 1:NData){
                        yy1[jj] <- yy[jj]-Err[jj]
                        if (svalue(T1obj7,index=TRUE) == 3) { # Y-log scale selected
                            yy1[jj] <- log10(yy1[jj])
                        }
                        lsegments(x0=xx[jj]-EndLngth, y0=yy1[jj], x1=xx[jj]+EndLngth, y1=yy1[jj], col=colr) #plots the upper part of the segment
                    }
                }
            }
            trellis.unfocus()
   }

   setRange <- function(){
               x1 <- as.numeric(svalue(XX1))
               y1 <- as.numeric(svalue(YY1))
               x2 <- as.numeric(svalue(XX2))
               y2 <- as.numeric(svalue(YY2))
               if (!is.na(x1) && !is.na(x2)){
                   if (svalue(RevAxis)) { #Binding energy set
                       Plot_Args$xlim  <<- Xlim <<- sort(c(x1, x2), decreasing=TRUE)
                   } else {
                       Plot_Args$xlim  <<- Xlim <<- sort(c(x1, x2))
                   }
               } else {
                   Plot_Args$xlim <<- range(FName[[SpectIndx]]@.Data[[1]])
               }
               if (!is.na(y1) && !is.na(y2)){
                   Plot_Args$ylim  <<- Ylim <<- sort(c(y1, y2))
               } else {
                   Plot_Args$ylim <<- range(FName[[SpectIndx]]@.Data[[2]])
               }
               CtrlPlot()
   }

#--- Routine for drawing Custom Axis
   CustomAx <- function(CustomDta){

               AxWin <- gwindow(title="CUSTOM AXIS", visible=FALSE)
               AxGroup1 <- ggroup(horizontal=FALSE, container=AxWin)
               txt1="1) Set your Min and Max scale value and the number of ticks Sub Ticks"
               glabel(txt1, container=AxGroup1)
               AxFrame <- gframe("SET SCALE RANGE", horizontal=FALSE, container=AxGroup1)
               AxLayout <- glayout(homogeneous=FALSE, spacing=3, container=AxFrame)
               X1 <- as.character(round(CustomDta[[1]], 2))
               X2 <- as.character(round(CustomDta[[2]], 2))
               msg <- paste("Xmin (min value=", X1, "):", sep="")
               AxLayout[1,1] <- EditXmin <- gedit(initial.msg =msg, width=40, container=AxLayout)
               msg <- paste("Xmax (max value=", X2, "):", sep="")
               AxLayout[1,2] <- EditXmax <- gedit(initial.msg =msg, width=40, container=AxLayout)
               AxLayout[2,1] <- EditNTicks <- gedit(initial.msg ="N. Ticks", container=AxLayout)
               AxLayout[2,2] <- EditStart <- gedit(initial.msg ="Starting at", container=AxLayout)
               AxLayout[3,1] <- EditStep <- gedit(initial.msg ="Step", container=AxLayout)

               gbutton("     SAVE & EXIT      ", handler=function(h,...){
                        AXmin <- as.numeric(svalue(EditXmin))
                        AXmax <- as.numeric(svalue(EditXmax))
                        NTicks <- as.numeric(svalue(EditNTicks))
                        Start <- as.numeric(svalue(EditStart))
                        Step <- as.numeric(svalue(EditStep))
                        if (is.null(AXmin)){
                            gmessage("Please Xmin value required!", icon="warning")
                        }
                        if (is.null(NTicks)){
                            gmessage("Please Xmax value required!", icon="warning")
                        }
                        if (is.null(NTicks)){
                            gmessage("Please N. Major Ticks  required!", icon="warning")
                        }
                        if (is.null(NTicks)){
                            gmessage("Please Start value required!", icon="warning")
                        }
                        if (is.null(NTicks)){
                            gmessage("Please Step value required!", icon="warning")
                        } else {
#                           dx <- (AXmax-AXmin)/NTicks
                           End <- Start+Step*NTicks
                           AXstp <- seq(from=Start, to=End, by=Step)
                           Ticklabels <- as.character(round(AXstp,digits=1))
                           if (CustomDta[[3]] == "X") {
                              Plot_Args$scales$x <<- list(at=AXstp, labels=Ticklabels)
                           } else if (CustomDta[[3]] == "Y") {
                              Plot_Args$scales$y <<- list(at=AXstp, labels=Ticklabels)
                           }
                           dispose(AxWin)
                           CtrlPlot()
                           Plot_Args$scales$relation <<- "same"
                        }
               }, container = AxFrame)
               visible(AxWin) <- TRUE
   }

   SetXYplotData <- function() {

               NComp <- length(FName[[SpectIndx]]@Components)
               idx <- svalue(T1obj8, index=TRUE)
               select <- ""
               code <- vector()
               if (idx==1) {
                  select <- "MAIN"   #plot raw data
               } else {
                  select <- "RTF"    #plot RegionToFit
               }
               code <- 1     #code is the label identifying the group of data possessing
#                          #same properties (linetype, lwd, color etc...)
#                          #and at the same time in which order the style options have to be applied
               if (svalue(BaseLineCK)==TRUE) {
                   select <- c(select, "BASE")
                   code <- c(code, 2)
               }
               if (svalue(ComponentCK)==TRUE) {
                   select <- c(select, "COMPONENTS")
                   code <- c(code, (3:(NComp+2)))
               }
               if (svalue(FitLineCK)==TRUE) {
                   select <- c(select, "FIT")
                   code <- c(code, (NComp+3))
               }
               tmp <- asList(FName[[SpectIndx]],select=select) #from coreline FName[[SpectIndx]] extract the selecteed regions
               X <- tmp$x # x list
               Y <- tmp$y # y list
               revAx <- svalue(RevAxis)
               if (length(Xlim)==0 || length(Ylim)==0){
                  Xlim  <<-sort(range(X, na.rm=TRUE), decreasing=revAx)
                  Ylim  <<-sort(range(Y, na.rm=TRUE))
                  Plot_Args$xlim <<- Xlim
                  Plot_Args$ylim <<- Ylim
               }
               Ylength <- lapply(Y, length)
               Ylength <- as.array(as.integer(Ylength))
               labX <- list()
	              levelX <- list()
               NN <- length(Ylength)
               for (ii in 1:NN){
                   labX[[ii]] <- rep(code[ii], times=Ylength[ii])
               }

               df <- data.frame(x = unname(unlist(X)), y = unname(unlist(Y)) )
               Plot_Args$data <<- df
               Plot_Args$groups <<- unlist(labX)

               if ( svalue(T2obj4) == "ON") {
                  Plot_Args$type <<-"l"
                  Plot_Args$col <<- svalue(T2obj1)
                  Plot_Args$lty <<- LType[svalue(T2obj2, index=TRUE)]
                  Plot_Args$lwd <<- svalue(T2obj3)
                  if (svalue(BaseLineCK)==TRUE) {
                     Plot_Args$col <<- c(Plot_Args$col, svalue(T3obj1))
                     Plot_Args$lty <<- c(Plot_Args$lty, LType[svalue(T3obj2, index=TRUE)])
                     Plot_Args$lwd <<- c(Plot_Args$lwd, svalue(T3obj3))
                  }
                  if (svalue(ComponentCK)==TRUE) {
                     Plot_Args$col <<- c(Plot_Args$col, rep(svalue(T4obj1), NComp))
                     Plot_Args$lty <<- c(Plot_Args$lty, rep(LType[svalue(T4obj2, index=TRUE)], NComp))
                     Plot_Args$lwd <<- c(Plot_Args$lwd, rep(svalue(T4obj3), NComp))
                  }
                  if (svalue(FitLineCK)==TRUE) {
                     Plot_Args$col <<- c(Plot_Args$col, svalue(T5obj1))
                     idx <- svalue(T5obj2, index=TRUE)
                     Plot_Args$lty <<- c(Plot_Args$lty, LType[idx])
                     Plot_Args$lwd <<- c(Plot_Args$lwd, svalue(T5obj3))
                  }
               }
               if ( svalue(T2obj7) == "ON") {
                  Plot_Args$type  <<-"p"
                  Plot_Args$col <<- svalue(T2obj1)
                  Plot_Args$pch <<- STypeIndx[svalue(T2obj5, index=TRUE)]
                  Plot_Args$cex <<- svalue(T2obj6)
                  if (svalue(BaseLineCK)==TRUE) {
                     Plot_Args$col <<- c(Plot_Args$col, svalue(T3obj1))
                     Plot_Args$pch <<- c(Plot_Args$pch, STypeIndx[svalue(T3obj4, index=TRUE)])
                     Plot_Args$cex <<- c(Plot_Args$cex, svalue(T3obj5))
                  }
                  if (svalue(ComponentCK)==TRUE) {
                     Plot_Args$col <<- c(Plot_Args$col, rep(svalue(T4obj1), NComp))
                     Plot_Args$pch <<- c(Plot_Args$pch, rep(STypeIndx[svalue(T4obj4, index=TRUE)], NComp))
                     Plot_Args$cex <<- c(Plot_Args$cex, rep(svalue(T4obj5), NComp))
                  }
                  if (svalue(FitLineCK)==TRUE) {
                     Plot_Args$col <<- c(Plot_Args$col, svalue(T5obj1))
                     Plot_Args$pch <<- c(Plot_Args$pch, STypeIndx[svalue(T5obj4, index=TRUE)])
                     Plot_Args$cex <<- c(Plot_Args$cex, svalue(T5obj5))
                  }
               }
               if ( svalue(T2obj4) == "ON" && svalue(T2obj7) == "ON") {
                  Plot_Args$type <<- "b"
               }
               if ( svalue(T2obj4) == "OFF" && svalue(T2obj7) == "OFF") {
                     Plot_Args$type <<- "x"   #linetype not defined figure cancelled!
               }
               CtrlPlot()
   }
   

   ResetPlot <- function(){
               svalue(RevAxis) <<- FALSE
               svalue(T1obj8) <<- 1 #original XY range
               svalue(BaseLineCK) <<- FALSE
               svalue(ComponentCK) <<- FALSE
               svalue(FitLineCK) <<- FALSE
               svalue(legendCK) <<- "FALSE"
               NComp <<- length(FName[[SpectIndx]]@Components)
               SampData <<- as(FName[[SpectIndx]],"matrix") #put spectrum, baseline, etc... in a matrix
               NColS <<- ncol(SampData)
               Xlim <<- sort(range(SampData[,1]))
               Ylim <<- sort(range(SampData[,2]))

               Plot_Args$xlim <<- Xlim
               Plot_Args$ylim <<- Ylim
               Plot_Args$data <<- data.frame(x=SampData[,1], y=SampData[,2])
               Plot_Args$pch <<- 1
               Plot_Args$cex <<- 1
               Plot_Args$lty <<- "solid"
               Plot_Args$lwd <<- 1
               Plot_Args$type <<- "l"
               Plot_Args$background <<- "transparent"
               Plot_Args$main <<- list(label <<- SpectName,cex=1.4)
               Plot_Args$xlab <<- list(label <<- FName[[SpectIndx]]@units[1], rot=0, cex=1.2)
               Plot_Args$ylab <<- list(label <<- FName[[SpectIndx]]@units[2], rot=90, cex=1.2)
               Plot_Args$scales <<- list(cex=1, tck=c(1,0), alternating=c(1), x=list(log=FALSE), y=list(log=FALSE))
               Plot_Args$xscale.components <<- xscale.components.subticks
               Plot_Args$yscale.components <<- yscale.components.subticks
               Plot_Args$col <<- "black"
               Plot_Args$groups <<- NULL
               Plot_Args$auto.key <<- FALSE
               Xlabel <<- FName[[SpectIndx]]@units[1]
               Ylabel <<- FName[[SpectIndx]]@units[2]

               AutoKey_Args <<- list(space="top",
                                    text=SpectName,
                                    cex = 1,
                                    type= "l",
                                    lines=TRUE,
                                    points=FALSE,
                                    border=FALSE,
                                    list(corner=NULL,x=NULL,y=NULL)
                                   )
               svalue(T1obj9) <<- "1.4"      #Titile size
               svalue(T1obj12) <<- "1"       #axis scale size
               svalue(T1obj13) <<- "1.2"     #axis label size
  }


#===== VARIABLES =====

   FName <- get(activeFName, envir=.GlobalEnv)
   ActiveFName <- get("activeFName", envir=.GlobalEnv)
   FNameList <- XPSFNameList()                 #list of all the XPSSamples
   FNameIdx <- grep(ActiveFName,FNameList)
   SpectIndx <- get("activeSpectIndx", envir=.GlobalEnv)
   SpectName <- get("activeSpectName", envir=.GlobalEnv)
   SpectList <- XPSSpectList(ActiveFName)      #list of all the corelines of the activeXPSSample
   OrigData <- as(FName[[SpectIndx]],"matrix")
   SampData <- as(FName[[SpectIndx]],"matrix")

   NComp <- length(FName[[SpectIndx]]@Components)
   NColS <- ncol(SampData)
   Xlim <- sort(range(SampData[,1]))
   Ylim <- sort(range(SampData[,2]))
   Xlabel <- FName[[SpectIndx]]@units[1]
   Ylabel <- FName[[SpectIndx]]@units[2]

   FitComp1 <- ""  #FitComp1 == vector containing the names of the fit components of the active spectrum
   for (ii in 1:NComp){
      FitComp1[ii] <- paste("C",ii, sep="")
   }
   Colors <- c("black", "red", "limegreen", "blue", "magenta", "orange", "cadetblue", "sienna", "darkgrey", "darkgreen", "gold", "darkviolet", "yellow", "cyan", "lightblue", "turquoise", "pink", "wheat", "thistle", "grey40")
   LType <- c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash", "F8", "431313", "22848222")
   LineTypes <- c("Solid", "Dashed", "Dotted", "Dotdash", "Longdash", "Twodash", "ExtraLongDash", "DashTwoDots", "DashTwoDotsDash")
   SType <- c("VoidSquare", "VoidCircle", "VoidTriangle", "Cross", "Diamond", "SolidSquare", "SolidCircle", "SolidTriangle", "SolidDiamond")
   STypeIndx <- c(0, 1, 2, 3, 5, 15, 16, 17, 18)
   LWidth <- c(0,1,1.5,2,2.5,3,3.5,4)
   SymSize <- c(0,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2)
   LCol <- 24    #initialize the color line to black
   LW <- 1       #initialize the linewidth to 1
   FontSize <- c(0.6, 0.8,1,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6,2.8,3)
   LegPos <- c("OutsideCenterTop", "OutsideTopRight", "OutsideTopLeft", "OutsideCenterRight", "OutsideCenterLeft",
             "OutsideCenterBottom", "InsideTopRight", "InsideTopLeft", "InsideBottomRight", "InsideBottomLeft")
   Orient <- c("Vertical", "Horizontal")
   LineWdh <- c(1,1.5,2,2.5,3,3.5,4,4.5,5)
   TxtCol <- c("Color", "Black")
   TxtSize <- c(0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6,2.8,3)
   Dist <- c(0.01,0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.2)

   Plot_Args <- list(x=formula("y ~ x"),data=NULL,xlim=NULL,ylim=NULL,
                     pch=1,cex=1,lty="solid",lwd=1,type="l",               #default settings for point and lines
                     background="transparent",
                     xlab=list(label=NULL, rot=0, cex=1.2),
                     ylab=list(label=NULL, rot=90, cex=1.2),
                     scales=list(cex=1, tck=c(1,0), alternating=c(1), x=list(log=FALSE), y=list(log=FALSE)),
                     xscale.components = xscale.components.subticks,
                     yscale.components = yscale.components.subticks,
                     main=list(label=NULL,cex=1.4),
                     col="black", groups=NULL,
                     par.settings = list(superpose.line=list(col=Colors)), #needed to set colors
                     auto.key = FALSE
                   )

   AutoKey_Args <- list( space="top",
                         text=SpectName,
                         cex = 1,
                         type= "l",
                         lines=TRUE,
                         points=FALSE,
                         border=FALSE,
                         list(corner=NULL,x=NULL,y=NULL)
                       )

#===== Reset graphic window =====

   plot.new()
   assign("MatPlotMode", FALSE, envir=.GlobalEnv)  #basic matplot function used to plot data


#===== NoteBook =====

   win <- gwindow("CUSTOM PLOT", visible=FALSE)
   maingp <- ggroup(horizontal=FALSE,container=win)
   NoteBk <- gnotebook(expand=TRUE, container = maingp)

# --- Tab1: Axes Options ---
       T1group1 <- ggroup(label="AXES", container=NoteBk)

       layoutAxis <- glayout(homogeneous=FALSE, spacing=3, container=T1group1)

       layoutAxis[1,1] <- T1frame1 <- gframe("XPSdata SELECTION", spacing=5, container=layoutAxis)
       T1obj1 <- gcombobox(FNameList, selected=FNameIdx, editable=FALSE, handler=function(h,...){
                             SelectedFName <- svalue(T1obj1)
                             FName <<- get(SelectedFName,envir=.GlobalEnv)  #load the XPSSample
                             SpectList <<- XPSSpectList(SelectedFName)
                             SpectIndx <<- 1
                             delete(T1frame2, T1obj2)
                             T1obj2  <<-gcombobox(SpectList, selected=-1, editable=FALSE, handler=function(h,...){
                                                 SpName <- svalue(T1obj2)
                                                 SpName <- unlist(strsplit(SpName, "\\."))
                                                 SpectIndx <<- as.numeric(SpName[1])
                                                 SpectName <<- SpName[2]
                                                 assign("activeSpectName",SpectName,envir=.GlobalEnv)
                                                 assign("activeSpectIndx",SpectIndx,envir=.GlobalEnv)
                                                 ResetPlot()

#                                                 AutoKey_Args$text <<- SpectName
#                                                 Plot_Args$main$label <<- SpectName
#                                                 Plot_Args$xlab$label <<- FName[[SpectIndx]]@units[1]
#                                                 Plot_Args$ylab$label <<- FName[[SpectIndx]]@units[2]
#                                                 AutoKey_Args$text <<- SpectName
#                                                 SampData <<- as(FName[[SpectIndx]],"matrix") #save X, spectrum, baseline, Components and Fit in a matrix
#                                                 Xlim  <<-sort(range(SampData[,1]))
#                                                 Ylim  <<-sort(range(SampData[,2]))
                                                 NColS <<- ncol(SampData)
	                                                wdth <- Xlim[2]-Xlim[1]
	                                                Xlim[1] <<- Xlim[1]-wdth/15
	                                                Xlim[2] <<- Xlim[2]+wdth/15
	                                                wdth <- Ylim[2]-Ylim[1]
	                                                Ylim[1] <<- Ylim[1]-wdth/15
	                                                Ylim[2] <<- Ylim[2]+wdth/15

                                                 if ( FName[[SpectIndx]]@Flags[1]) {   #reverse if BE is set
                                                    Xlim <<- rev(Xlim)
                                                 }
                                                 Plot_Args$xlim <<- Xlim
                                                 Plot_Args$ylim <<- Ylim
                                                 Plot_Args$data <<- data.frame(x=SampData[,1], y=SampData[,2])
                                                 if (NColS > 2) {
                                                    BaseLine <- SampData[,3]
                                                 }
                                                 if (NColS > 3){
                                                    FitComp <- SampData[,4:NColS-1] #skip the first three column
                                                    SpectFit <- SampData[,NColS]
                                                 }
                                                 CtrlPlot()
                             }, container=T1frame2)
                             add(T1frame2, T1obj2)
                             CtrlPlot() }, container=T1frame1)


       layoutAxis[1,2] <- T1frame2 <- gframe("CORE LINE SELECTION", spacing=5, container=layoutAxis)
       T1obj2 <- gcombobox(SpectList, selected=-1, editable=FALSE, handler=function(h,...){
                             SpName <- svalue(T1obj2)
                             SpName <- unlist(strsplit(SpName, "\\."))
                             SpectIndx <<- as.numeric(SpName[1])
                             SpectName <<- SpName[2]
                             assign("activeSpectName",SpectName,envir=.GlobalEnv)
                             assign("activeSpectIndx",SpectIndx,envir=.GlobalEnv)
                             ResetPlot()
#                             AutoKey_Args$text <<- SpectName
#                             Plot_Args$main$label <<- SpectName
#                             Plot_Args$xlab$label <<- FName[[SpectIndx]]@units[1]
#                             Plot_Args$ylab$label <<- FName[[SpectIndx]]@units[2]
#                             SampData <<- as(FName[[SpectIndx]],"matrix")
#                             Xlim <<- sort(range(SampData[,1]))
#                             Ylim <<- sort(range(SampData[,2]))
                             NColS <<- ncol(SampData)
	                            wdth <- Xlim[2]-Xlim[1]
	                            Xlim[1] <<- Xlim[1]-wdth/15
	                            Xlim[2] <<- Xlim[2]+wdth/15
	                            wdth <- Ylim[2]-Ylim[1]
	                            Ylim[1] <<- Ylim[1]-wdth/15
	                            Ylim[2] <<- Ylim[2]+wdth/15

                             Plot_Args$xlim <<- Xlim
                             Plot_Args$ylim <<- Ylim
                             Plot_Args$data <<- data.frame(x=SampData[,1], y=SampData[,2])
                             if (NColS > 2) {
                                BaseLine <- SampData[,3]
                             }
                             if (NColS > 3){
                                FitComp <- SampData[,4:NColS-1]
                                SpectFit <- SampData[,NColS]  #fit
                             }
#                             CtrlPlot() }, container=T1frame2)
                             CtrlPlot() 
                             }, container=T1frame2)


      layoutAxis[1,3] <- T1frame3 <- gframe("REVERSE X-axis", spacing=5, container=layoutAxis)
      RevAxis <- gcheckbox("Reverse X axis",checked=FALSE, handler=function(h,...){
                             if ( svalue(RevAxis)) {   #reverse scale if checkbox TRUE
                                Xlim <<- rev(Xlim)
                                Plot_Args$xlim <<- Xlim
                             } else {
                                Xlim <<- sort(Xlim, decreasing=FALSE)
                                Plot_Args$xlim <<- Xlim
                             }
                             CtrlPlot()
                 }, container=T1frame3)


      layoutAxis[2,1] <- T1frame6 <- gframe("TICKS", spacing=5, container=layoutAxis)
      T1obj6 <- gcombobox(c("Left & Bottom", "Top & Right", "Both", "Custom X", "Custom Y"), selected=1, editable=FALSE, handler= function(h,...){
                             idx <- svalue(T1obj6,index=TRUE)
                             if (idx==1) {
                                Plot_Args$scales$tck <<- c(1,0)
                                Plot_Args$scales$alternating <<- c(1)
                             } else if (idx==2) {
                                Plot_Args$scales$tck <<- c(0,1)
                                Plot_Args$scales$alternating <<- c(2)
                             } else if (idx==3) {
                                Plot_Args$scales$tck <<- c(1,1)
                                Plot_Args$scales$alternating <<- c(3)
                             } else if (idx==4 || idx==5) {
                                Plot_Args$scales$relation <<- "free"
                                if (svalue(T1obj6)=="Custom X") {
                                   CustomDta <- list(Xlim[1], Xlim[2], "X")
                                   CustomAx(CustomDta)
                                }
                                if (svalue(T1obj6)=="Custom Y") {
                                   txt1="1) Ymin, Ymax and the number of ticks on the Y axis: es. Ymin=0, Ymax=35, Nticks=7"
                                   txt2="2) Set Tick-Labels (as many labels as the ticks): es. Tick Labels= 0,5, ,15,20, ,30"
                                   CustomDta <- list(Ylim[1], Ylim[2], "Y")
                                   CustomAx(CustomDta)
                                }
                             }
                             CtrlPlot() }, container=T1frame6)

      layoutAxis[2,2] <- T1frame7 <- gframe("SCALE", spacing=5, container=layoutAxis)
      T1obj7 <- gcombobox(c("Standard", "X-Log10", "Y-Log10", "XY-Log10", "X-^10", "Y-^10", "XY-^10"), selected=1, editable=FALSE, handler= function(h,...){
                             idx <- svalue(T1obj7,index=TRUE)
                             Xlim <- sort(range(SampData[,1]))
                             Ylim <- sort(range(SampData[,2]))
                             if (idx == 1) {
                                Plot_Args$scales$x <<- list(log=FALSE)
                                Plot_Args$xscale.components <<- xscale.components.subticks
                                Plot_Args$scales$y <<- list(log=FALSE)
                                Plot_Args$yscale.components <<- yscale.components.subticks
                             } else if (idx == 2) {  #X-log10 scale
                                Xlim <- sort(range(SampData[,1]))
                                if (Xlim[1] < 0) {
                                    gmessage("Cannot plot negatige X-values !", title="WRONG X VALUES", icon="warning")
                                    return()
                                }
                                Plot_Args$scales$y <<- list(log=FALSE)
                                Plot_Args$yscale.components <<- yscale.components.subticks
                                Plot_Args$scales$x <<- list(log = 10)
                                Plot_Args$xscale.components <<- xscale.components.log10ticks
                             } else if (idx == 3) {  #Y-log10 scale
                                Ylim <- sort(range(SampData[,2]))
                                if (Ylim[1] < 0) {
                                    gmessage("Cannot plot negatige Y-values !", title="WRONG Y VALUES", icon="warning")
                                    return()
                                }
                                Plot_Args$scales$x <<- list(log=FALSE)
                                Plot_Args$xscale.components <<- xscale.components.subticks
                                Plot_Args$scales$y <<- list(log = 10)
                                Plot_Args$yscale.components <<- yscale.components.log10ticks
                             } else if (idx == 4) {  #XY-log10 scale
                                Xlim <- sort(range(SampData[,1]))
                                Ylim <- sort(range(SampData[,2]))
                                if (Xlim[1] < 0 || Ylim[1] < 0) {
                                    gmessage("Cannot plot negatige X, Y values !", title="WRONG X, Y VALUES", icon="warning")
                                    return()
                                }
                                Plot_Args$scales$x <<- list(log = 10)
                                Plot_Args$xscale.components <<- xscale.components.log10ticks
                                Plot_Args$scales$y <<- list(log = 10)
                                Plot_Args$yscale.components <<- yscale.components.log10ticks
                            } else if (idx == 5) {
                                Plot_Args$scales$x$log <<- TRUE
                                Plot_Args$scales$y$log <<- FALSE
                             } else if (idx == 6) {
                                Plot_Args$scales$x$log <<- FALSE
                                Plot_Args$scales$y$log <<- TRUE
                             } else if (idx == 7) {
                                Plot_Args$scales$x$log <<- TRUE
                                Plot_Args$scales$y$log <<- TRUE
                             }
                             CtrlPlot() }, container=T1frame7)

      layoutAxis[2,3] <- T1frame8 <- gframe("XY range", spacing=5, container=layoutAxis)
      T1obj8 <- gradio(items=c("Original_XYrange", "Fitted XY range"), selected=1, horizontal=FALSE,  handler=function(h,...){
                             if (length(svalue(T1obj2)) == 0) {
                                gmessage(msg="Please Select the Core Line!" , title = "No spectral Data selected",  icon = "warning")
                                return()
                             } else {
                                svalue(T1obj8) <- "Fitted XY range"
                                SetXYplotData()
                             }
                             CtrlPlot() }, container=T1frame8)


      layoutAxis[3,1] <- T1frame9 <- gframe("TITLE SIZE", spacing=5, container=layoutAxis)
      T1obj9 <- gcombobox(FontSize, selected=5, editable=FALSE, handler= function(h,...){
                             Plot_Args$main$cex <<- svalue(T1obj9)
                             CtrlPlot() }, container=T1frame9)

      layoutAxis[3,2] <- T1frame10 <- gframe("CHANGE TITLE", spacing=5, container=layoutAxis)
      T1obj10 <- gedit("", container=T1frame10, handler=function(h,...){
                             Plot_Args$main$label <<- svalue(T1obj10)
                             CtrlPlot() })

      layoutAxis[3,3] <- T1frame11 <- gframe("Exact Range Values", horizontal=FALSE, spacing=5, container=layoutAxis)
      T1group11a <- ggroup(horizontal=TRUE, container=T1frame11)
      XX1 <- gedit("", initial.msg = "Xmin= ", handler=function(h, ...){setRange()}, container=T1group11a)
      XX2 <- gedit("", initial.msg = "Xmax= ", handler=function(h, ...){setRange()}, container=T1group11a)
      T1group11b <- ggroup(horizontal=TRUE, container=T1frame11)
      YY1 <- gedit("", initial.msg = "Ymin= ", handler=function(h, ...){setRange()}, container=T1group11b)
      YY2 <- gedit("", initial.msg = "Ymax= ", handler=function(h, ...){setRange()}, container=T1group11b)
      tkconfigure(XX1$widget, width=10)
      tkconfigure(XX2$widget, width=10)
      tkconfigure(YY1$widget, width=10)
      tkconfigure(YY2$widget, width=10)


      layoutAxis[4,1] <- T1frame12 <- gframe("AXIS SCALE SIZE", spacing=5, container=layoutAxis)
      T1obj12 <- gcombobox(FontSize, selected=3, editable=FALSE, handler= function(h,...){
                             Plot_Args$scales$cex <<- svalue(T1obj12)
                             CtrlPlot() }, container=T1frame12)

      layoutAxis[4,2] <- T1frame13 <- gframe("AXIS LABEL SIZE", spacing=5, container=layoutAxis)
      T1obj13 <- gcombobox(FontSize, selected=3, editable=FALSE, handler= function(h,...){
                             Plot_Args$xlab$cex <<- svalue(T1obj13)
                             Plot_Args$ylab$cex <<- svalue(T1obj13)
                             CtrlPlot() }, container=T1frame13)


      layoutAxis[5,1] <- T1frame14 <- gframe("CHANGE X-LABEL", spacing=5, container=layoutAxis)
      T1obj14 <- gedit("", container=T1frame14, handler=function(h,...){
                             Plot_Args$xlab$label <<- svalue(T1obj14)
                             Xlabel <<- svalue(T1obj14)
                             CtrlPlot() })

      layoutAxis[5,2] <- T1frame15 <- gframe("CHANGE Y-LABEL", spacing=5, container=layoutAxis)
      T1obj15 <- gedit("", container=T1frame15,handler=function(h,...){
                             Plot_Args$ylab$label <<- svalue(T1obj15)
                             Ylabel <<- svalue(T1obj15)
                             CtrlPlot() })
                             

# --- Tab2: Spectrum Options ---

      T2group1 <- ggroup(label="SPECTRUM OPTIONS", horizontal=FALSE, container=NoteBk)

      layoutT2 <- glayout(homogeneous=FALSE, spacing=3, container=T2group1)

      layoutT2[1,1] <- T2frame1 <- gframe("COLOR", spacing=5, container=layoutT2)
      T2obj1 <- gcombobox(Colors, selected=1, editable=FALSE, handler=function(h,...){
                            Plot_Args$col <<- svalue(T2obj1)
                            SetXYplotData()
                          }, container=T2frame1)

      layoutT2[2,1] <- T2frame2 <- gframe("LINE TYPE", spacing=5, container=layoutT2)
      T2obj2 <- gcombobox(LineTypes, selected=1, editable=FALSE, handler=function(h,...){
                              Plot_Args$type <<- "l"
                              idx <- as.numeric(svalue(T2obj2, index=TRUE))
                              Plot_Args$lty <<- LType[idx]
                              Plot_Args$par.settings = list(superpose.line=list(col=Colors, lwd=1)) #needed to set legend colors
                              AutoKey_Args$lines <<- TRUE
                              AutoKey_Args$points <<- FALSE
                              if (svalue(T2obj4)=="ON") SetXYplotData()
                            }, container=T2frame2)

      layoutT2[2,2] <- T2frame3 <- gframe("LINE WIDTH", spacing=5, container=layoutT2)
      T2obj3 <- gcombobox(LWidth, selected=2, editable=FALSE, handler= function(h,...){
                              Plot_Args$lwd <<- as.numeric(svalue(T2obj3))
                              if (svalue(T2obj4)=="ON") SetXYplotData()
                            }, container=T2frame3)

      layoutT2[2,3] <- T2frame4 <- gframe("SET LINES", spacing=5, container=layoutT2)
      T2obj4 <- gradio(c("ON", "OFF"), selected=1, horizontal = TRUE, handler=function(h,...){
                              SetXYplotData()
                            }, container=T2frame4)


      layoutT2[3,1] <- T2frame5 <- gframe("SYMBOL", spacing=5, container=layoutT2)
      T2obj5 <- gcombobox(SType, selected=1, editable=FALSE, handler=function(h,...){
                              Plot_Args$type <<- "p"
                              Plot_Args$pch <<- STypeIndx[svalue(T2obj5, index=TRUE)]
                              Plot_Args$par.settings = list(superpose.symbol=list(col=Colors))
                              AutoKey_Args$lines <<- FALSE
                              AutoKey_Args$points <<- TRUE
                              if (svalue(T2obj7)=="ON") SetXYplotData()
                            }, container=T2frame5)

      layoutT2[3,2] <- T2frame6 <- gframe("SYMSIZE", spacing=5, container=layoutT2)
      T2obj6 <- gcombobox(SymSize, selected=5, editable=FALSE, handler= function(h,...){
                              Plot_Args$cex <<- svalue(T2obj6)
                              if (svalue(T2obj7)=="ON") SetXYplotData()
                            }, container=T2frame6)


      layoutT2[3,3] <- T2frame7 <- gframe("SET SYMBOLS", horizontal=TRUE, spacing=5, container=layoutT2)
      T2obj7 <- gradio(c("ON", "OFF"), selected=2, horizontal=TRUE, handler=function(h,...){
                              SetXYplotData()
                            }, container=T2frame7)

      T2frame8 <- gframe("ERROR BARS", horizontal=FALSE, spacing=5, container=T2group1)
      txt1 <- "Only generic X, Y data with Standard Errors loaded with XPSImport.asciiGUI()\n"
      txt2 <- "can be plotted using these options"
      glabel(paste(txt1, txt2), sep="", container=T2frame8)
      T2group2 <- ggroup(horizontal=TRUE, container=T2frame8)
      T2obj8 <- gcheckboxgroup(c("Upper Bar", "Lower Bar", "Upper End", "Lower End"), selected=-1, horizontal=TRUE, handler=function(h,...){
                              CtrlPlot()
                            }, container=T2group2)

      glabel("End Amplitude", container=T2group2)
      T2obj9 <- gedit(text = "1", initial.msg = "End Amplitude", handler=function(h,...){
                              CtrlPlot()
                            }, container=T2group2)



# --- Tab3: BaseLine Options ---

      T3group1 <- ggroup(label="BASELINE OPTIONS", container=NoteBk)

      layoutT3 <- glayout(homogeneous=FALSE, spacing=3, container=T3group1)

      layoutT3[1,1] <- T3CKframe <- gframe(text="Set BaseLine", spacing=5, container=layoutT3)
      BaseLineCK <- gcheckbox("Base Line ON/OFF", checked=FALSE,handler=function(h,...){
                            if (svalue(BaseLineCK)==TRUE) {
                               LL <- length(FName[[SpectIndx]]@Baseline)
                               if (LL == 0) {
                                   gmessage(msg="SORRY, NO BASELINE FOUND!" , title = "BASELINE PLOTTING ABORTED",  icon = "warning")
                                   svalue(BaseLineCK) <- "FALSE"
                               }
                            }
                            SetXYplotData()
                            CtrlPlot()
                          }, container=T3CKframe)

      layoutT3[2,1] <- T3frame1 <- gframe("COLOR", spacing=5, container=layoutT3)
      T3obj1 <- gcombobox(Colors, selected=1, editable=FALSE, handler=function(h,...){
                            SetXYplotData()
                            CtrlPlot() },
                            container=T3frame1)

      layoutT3[3,1] <- T3frame2 <- gframe("LINE TYPE", spacing=5, container=layoutT3)
      T3obj2 <- gcombobox(LineTypes, selected=1, editable=FALSE, handler=function(h,...){
                            SetXYplotData()
                            CtrlPlot()
                          }, container=T3frame2)

      layoutT3[3,2] <- T3frame3 <- gframe("LINE WIDTH", spacing=5, container=layoutT3)
      T3obj3 <- gcombobox(LWidth, selected=2, editable=FALSE, handler=function(h,...){
                            SetXYplotData()
                            CtrlPlot()
                          }, container=T3frame3)

      layoutT3[4,1] <- T3frame4 <- gframe("SYMBOL", spacing=5, container=layoutT3)
      T3obj4 <- gcombobox(SType, selected=1, editable=FALSE, handler=function(h,...){
                            Plot_Args$type <<- "p"
                            SetXYplotData()
                            CtrlPlot()
                          }, container=T3frame4)

      layoutT3[4,2] <- T3frame5 <- gframe("SYMSIZE", spacing=5, container=layoutT3)
      T3obj5 <- gcombobox(SymSize, selected=1, editable=FALSE, handler=function(h,...){
                            SSize <- c(svalue(T3obj5), svalue(T2obj5)) #vector of CoreLine colors + ComponentiFit
                            Plot_Args$cex <<- SSize
                            SetXYplotData()
                            CtrlPlot() },
                            container=T3frame5)


# --- Tab4: Fit Components Options ---

      T4group1 <- ggroup(label="FIT COMPONENT OPTIONS", container=NoteBk)

      layoutT4 <- glayout(homogeneous=FALSE, spacing=3, container=T4group1)

      layoutT4[1,1] <- T4CKframe <- gframe(text="Set Components", spacing=5, container=layoutT4)
      ComponentCK <- gcheckbox("Components ON/OFF", checked=FALSE,handler=function(h,...){
                              if (svalue(ComponentCK)==TRUE) {
                                 if (NComp == 0) {
                                    gmessage(msg="SORRY, NO FIT FOUND!" , title = "FIT PLOTTING ABORTED",  icon = "warning")
                                    svalue(ComponentCK) <- "FALSE"
                                 }
                              }
                              SetXYplotData()
                              CtrlPlot()
                            }, container=T4CKframe)

      layoutT4[2,1] <- T4frame1 <- gframe("COLOR", spacing=5, container=layoutT4)
      T4obj1 <- gcombobox(Colors, selected=1, editable=FALSE, handler=function(h,...){
                              SetXYplotData()
                              CtrlPlot()
                            }, container=T4frame1)


      layoutT4[3,1] <- T4frame2 <- gframe("LINE TYPE", spacing=5, container=layoutT4)
      T4obj2 <- gcombobox(LineTypes, selected=1, editable=FALSE, handler=function(h,...){
                              SetXYplotData()
                              CtrlPlot()
                            }, container=T4frame2)

      layoutT4[3,2] <- T4frame3 <- gframe("LIINE WIDTH", spacing=5, container=layoutT4)
      T4obj3 <- gcombobox(LWidth, selected=2, editable=FALSE, handler=function(h,...){
                              SetXYplotData()
                              CtrlPlot()
                            }, container=T4frame3)

      layoutT4[4,1] <- T4frame4 <- gframe("SYMBOL", spacing=5, container=layoutT4)
      T4obj4 <- gcombobox(SType, selected=1, editable=FALSE, handler=function(h,...){
                              Plot_Args$type <<- "p"
                              SetXYplotData()
                              CtrlPlot()
                            }, container=T4frame4)

      layoutT4[4,2] <- T4frame5 <- gframe("SYMSIZE", spacing=5, container=layoutT4)
      T4obj5 <- gcombobox(SymSize, selected=1, editable=FALSE, handler=function(h,...){
                              SSize <- c(svalue(T3obj5), svalue(T2obj5)) #vector of CoreLine colors + ComponentiFit
                              Plot_Args$cex <<- SSize
                              SetXYplotData()
                              CtrlPlot()
                            },  container=T4frame5)


# --- Tab5: Fit Options ---

      T5group1 <- ggroup(label="FIT OPTIONS", container=NoteBk)

      layoutT5 <- glayout(homogeneous=FALSE, spacing=3, container=T5group1)

      layoutT5[1,1] <- T5CKframe <- gframe(text="Set Fit", spacing=5, container=layoutT5)
      FitLineCK <- gcheckbox("Fit Line ON/OFF", checked=FALSE,handler=function(h,...){
                            if (svalue(FitLineCK)==TRUE) {
                               if (NComp == 0) {
                                  gmessage(msg="SORRY, NO FIT FOUND!" , title = "FIT PLOTTING ABORTED",  icon = "warning")
                                  svalue(FitLineCK) <- "FALSE"
                               }
                            }
                            SetXYplotData()
                            CtrlPlot()
                          }, container=T5CKframe)

      layoutT5[2,1] <- T5frame1 <- gframe("COLOR", spacing=5, container=layoutT5)
      T5obj1 <- gcombobox(Colors, selected=1, editable=FALSE, handler=function(h,...){
                            SetXYplotData()
                            CtrlPlot()
                          }, container=T5frame1)

      layoutT5[3,1] <- T5frame2 <- gframe("LINE TYPE", spacing=5, container=layoutT5)
      T5obj2 <- gcombobox(LineTypes, selected=1, editable=FALSE, handler=function(h,...){
                            SetXYplotData()
                            CtrlPlot()
                          }, container=T5frame2)

      layoutT5[3,2] <- T5frame3 <- gframe("LINE WIDTH", spacing=5, container=layoutT5)
      T5obj3 <- gcombobox(LWidth, selected=2, editable=FALSE, handler=function(h,...){
                            SetXYplotData()
                            CtrlPlot() }, container=T5frame3)

      layoutT5[4,1] <- T5frame4 <- gframe("SYMBOL", spacing=5, container=layoutT5)
      T5obj4 <- gcombobox(SType, selected=1, editable=FALSE, handler=function(h,...){
                            Plot_Args$type <<- "p"
                            Symbol <- c(STypeIndx[svalue(T5obj4, index=TRUE)],
                                      STypeIndx[svalue(T4obj4, index=TRUE)],
                                      STypeIndx[svalue(T3obj4, index=TRUE)],
                                      STypeIndx[svalue(T2obj4, index=TRUE)])
                            Plot_Args$pch <<- Symbol
                            SetXYplotData()
                            CtrlPlot()
                          }, container=T5frame4)

      layoutT5[4,2] <- T5frame5 <- gframe("SYMSIZE", spacing=5, container=layoutT5)
      T5obj5 <- gcombobox(SymSize, selected=1, editable=FALSE, handler=function(h,...){
                            SSize <- c(svalue(T5obj5), svalue(T4obj5), svalue(T3obj5), svalue(T2obj5)) #vector of Coreline and fitComp color
                            Plot_Args$cex <<- SSize
                            SetXYplotData()
                            CtrlPlot()
                          }, container=T5frame5)


# --- Tab6: Legend Options ---

   T6group1 <- ggroup(label="LEGEND", container=NoteBk)

   layoutLeg <- glayout(homogeneous=FALSE, spacing=3, container=T6group1)

   layoutLeg[1,1] <- Lframe1 <- gframe(text="Set Legend", spacing=5, container=layoutLeg)
   legendCK <- gcheckbox("Legend ON/OFF", checked=FALSE,handler=function(h,...){
                          if (svalue(legendCK) == TRUE) {
		           	                   Plot_Args$auto.key <<-AutoKey_Args
                             if ( svalue(T2obj4) == "ON") {
                                Plot_Args$par.settings$superpose.line$col <<- svalue(T2obj1)
                                Plot_Args$par.settings$superpose.line$lty <<- LType[svalue(T2obj2, index=TRUE)]
                             }
                             if ( svalue(T2obj7) == "ON") {
                                Plot_Args$par.settings$superpose.symbol$col <<- svalue(T2obj1)
                                Plot_Args$par.settings$superpose.symbol$pch <<- STypeIndx[svalue(T2obj5, index=TRUE)]
                             }
                          } else {
		                           Plot_Args$auto.key <<- FALSE
	           	             }
                          CtrlPlot()
                       }, container=Lframe1)


   layoutLeg[2,1] <- Lframe2 <- gframe(text="Legend Position", spacing=5, container=layoutLeg)
   LegPosCK <- gcombobox(LegPos, selected = -1, toolkit = guiToolkit(), handler=function(h,...){
			                       switch(svalue(LegPosCK),
                              "OutsideCenterTop"    = { Plot_Args$auto.key$space <<-"top" },
				                          "OutsideTopRight"     = { Plot_Args$auto.key$space <<-NULL
                                                        Plot_Args$auto.key$corner <<- c(1,1)
                                                        Plot_Args$auto.key$x <<-0.95
                                                        Plot_Args$auto.key$y <<-1.05 },
				                          "OutsideTopLeft"      = { Plot_Args$auto.key$space <<-NULL
                                                        Plot_Args$auto.key$corner <<- c(0,1)
                                                        Plot_Args$auto.key$x <<-0.05
                                                        Plot_Args$auto.key$y <<-1.05 },
				                          "OutsideCenterRight"  = { Plot_Args$auto.key$space <<-"right" },
				                          "OutsideCenterLeft"   = { Plot_Args$auto.key$space <<-"left" },
			                           "OutsideCenterBottom" = { Plot_Args$auto.key$space <<-"bottom" },
				                          "InsideTopRight"      = { Plot_Args$auto.key$space <<-NULL
                                                        Plot_Args$auto.key$corner <<- c(1,1)
                                                        Plot_Args$auto.key$x <<-0.95
                                                        Plot_Args$auto.key$y <<-0.95 },
				                          "InsideTopLeft"       = { Plot_Args$auto.key$space <<-NULL
                                                        Plot_Args$auto.key$corner <<- c(0,1)
                                                        Plot_Args$auto.key$x <<-0.05
                                                        Plot_Args$auto.key$y <<-0.95 },
                              "InsideBottomRight"   = { Plot_Args$auto.key$space <<-NULL
                                                        Plot_Args$auto.key$corner <<- c(1,0)
                                                        Plot_Args$auto.key$x <<-0.95
                                                        Plot_Args$auto.key$y <<-0.05 },
				                          "InsideBottomLeft"    = {	Plot_Args$auto.key$space <<-NULL
                                                        Plot_Args$auto.key$corner <<- c(0,0)
                                                        Plot_Args$auto.key$x <<-0.05
                                                        Plot_Args$auto.key$y <<-0.05 },
                          )
                          CtrlPlot()
                       }, container=Lframe2)

   layoutLeg[2,2] <- Lframe3 <- gframe(text="Legend Border", spacing=5, container=layoutLeg)
   BorderCK <- gcombobox(c("No", "Yes"),selected=-1, toolkit = guiToolkit(), handler=function(h,...){
                          BorderCK <- as.numeric(svalue(BorderCK, index=TRUE))
                          if (BorderCK == 1){
                             Plot_Args$auto.key$border <<- FALSE
                          } else if (BorderCK == 2){
                             Plot_Args$auto.key$border <<- TRUE
                          }
                          CtrlPlot()
                       }, container=Lframe3)

   layoutLeg[3,1] <- Lframe4 <- gframe(text="Line/Symbol weight", spacing=5, container=layoutLeg)
   LineWdhCK <- gcombobox(LineWdh,selected=1, toolkit = guiToolkit(), handler=function(h,...){
                          weight <- as.numeric(svalue(LineWdhCK))
                          if (length(svalue(T2obj4)) == 0) {   #ho scelto di plottare per linee

                             Plot_Args$par.settings$superpose.line$lty <<- LType[svalue(T2obj2, index=TRUE)]
                             Plot_Args$par.settings$superpose.line$lwd <<- weight
                          } else {
                             Plot_Args$par.settings$superpose.symbol$pch <<- STypeIndx[svalue(T2obj4, index=TRUE)]
                             Plot_Args$par.settings$superpose.symbol$cex <<- weight
                          }
                          CtrlPlot()
                       }, container=Lframe4)

   layoutLeg[3,2] <- Lframe5 <- gframe(text="Distance from Margin", spacing=5, container=layoutLeg)
   DistCK <- gcombobox(Dist,selected=5, toolkit = guiToolkit(), handler=function(h,...){
                          LegDist <- as.numeric(svalue(DistCK))
			                       switch(svalue(LegPosCK),
                              "OutsideTop"         = { Plot_Args$auto.key$space <<-"top"
                                                       Plot_Args$auto.key$y <<-1+LegDist },
				                          "OutsideTopRight"    = { Plot_Args$auto.key$space <<-NULL
                                                       Plot_Args$auto.key$corner <<- c(1,1)
                                                       Plot_Args$auto.key$x <<-0.95
                                                       Plot_Args$auto.key$y <<-1+LegDist },
				                          "OutsideTopLeft"     = { Plot_Args$auto.key$space <<-NULL
                                                       Plot_Args$auto.key$corner <<- c(0,1)
                                                       Plot_Args$auto.key$x <<-0.05
                                                       Plot_Args$auto.key$y <<-1+LegDist },
				                          "OutsideCenterRight" = { Plot_Args$auto.key$space <<-"right"
                                                       Plot_Args$par.settings$layout.widths$right.padding <<- 8-LegDist*40
                                                       Plot_Args$par.settings$layout.widths$key.right <<- LegDist*10 },
				                          "OutsideCenterLeft"  = { Plot_Args$auto.key$space <<-"left"
                                                       Plot_Args$par.settings$layout.widths$left.padding <<- 8-LegDist*40
                                                       Plot_Args$par.settings$layout.widths$key.left <<- LegDist*10 },
			                           "OutsideBottom"      = { Plot_Args$auto.key$space <<-"bottom"
                                                       Plot_Args$auto.key$y <<-1-LegDist },
				                          "InsideTopRight"     = { Plot_Args$auto.key$space <<-NULL
                                                       Plot_Args$auto.key$corner <<- c(1,1)
                                                       Plot_Args$auto.key$x <<-1-LegDist
                                                       Plot_Args$auto.key$y <<-1-LegDist },
				                          "InsideTopLeft"      = { Plot_Args$auto.key$space <<-NULL
                                                       Plot_Args$auto.key$corner <<- c(0,1)
                                                       Plot_Args$auto.key$x <<-LegDist
                                                       Plot_Args$auto.key$y <<-1-LegDist },
                              "InsideBottomRight"  = { Plot_Args$auto.key$space <<-NULL
                                                       Plot_Args$auto.key$corner <<- c(1,0)
                                                       Plot_Args$auto.key$x <<-1-LegDist
                                                       Plot_Args$auto.key$y <<-LegDist },
				                          "InsideBottomLeft"   = {	Plot_Args$auto.key$space <<-NULL
                                                       Plot_Args$auto.key$corner <<- c(0,0)
                                                       Plot_Args$auto.key$x <<-LegDist
                                                       Plot_Args$auto.key$y <<-LegDist },
                          )
                          CtrlPlot()
                       }, container=Lframe5)

   layoutLeg[4,1] <- Lframe6 <- gframe(text="Text Size", spacing=5, container=layoutLeg)
   TSizeCK <- gcombobox(TxtSize,selected=1, toolkit = guiToolkit(), handler=function(h,...){
		           	            Plot_Args$auto.key$cex <<-as.numeric(svalue(TSizeCK))
                          CtrlPlot()
                       }, container=Lframe6)

   layoutLeg[4,2] <- Lframe7 <- gframe(text="Legend text Color", spacing=5, container=layoutLeg)
   TxtColCK <- gcombobox(c("B/W", "Color"), selected=1, toolkit = guiToolkit(), handler=function(h,...){
                          if  (svalue(TxtColCK)=="B/W"){
                              Plot_Args$auto.key$col <<- "black"
                          } else {
                              Plot_Args$auto.key$col <<- svalue(T2obj1)
                          }
                          CtrlPlot()
                       }, container=Lframe7)

   layoutLeg[5,1] <- Lframe8 <- gframe(text="Change Legend", spacing=5, container=layoutLeg)
   NewLegend <- gedit(text="", initial.msg="New Label", container=Lframe8)
   addHandlerChanged(NewLegend,handler=function(h,...){
		           	            Plot_Args$auto.key$text <<- svalue(NewLegend)
                          CtrlPlot()
                       }, container=Lframe8)

   layoutLeg[5,2] <- Annotate <- gbutton(text=" Annotate ", handler=function(h,...){
                          XPSLattAnnotate(Xlim, Ylim)
                       }, container=layoutLeg)

#--- Common buttons

   groupBtn <- ggroup(label="BUTTONS", horizontal=TRUE, container=maingp)

   gbutton("               REFRESH                 ", expand=TRUE, handler=function(h,...){
                               SetXYplotData()
                            }, container = groupBtn)

   gbutton("               RESET PLOT              ", expand=TRUE, handler=function(h,...){
                               ResetPlot()
                               CtrlPlot()
                            }, container = groupBtn)

   gbutton("                 EXIT                  ", expand=TRUE, handler=function(h,...){
				                           dispose(win)
                            }, container = groupBtn)

   enabled(T3obj2) <- FALSE
   enabled(T3obj3) <- FALSE
   enabled(T4obj2) <- FALSE
   enabled(T4obj3) <- FALSE
   enabled(T5obj2) <- FALSE
   enabled(T5obj3) <- FALSE
   enabled(T3obj4) <- FALSE
   enabled(T3obj5) <- FALSE
   enabled(T4obj4) <- FALSE
   enabled(T4obj5) <- FALSE
   enabled(T5obj4) <- FALSE
   enabled(T5obj5) <- FALSE

   visible(win) <- TRUE
   for(ii in  6:1){
      svalue(NoteBk) <- ii #refresh notebook pages
   }
}



