#function to perform plots in overlay mode

#function to perform plots in overlay mode

#'Performs overlay of XPS-Spectra
#'
#'Provides a userfriendly interface to select XPS-Corelines to overlay
#'and a selection of plotting options for a personalized data representation
#'This function is based on the (/code{Lattice}) Package.  No parameters passed
#'No parameters passed to this function
#'
#'
#'@examples
#'
#'\dontrun{
#'	XPSOverlay()
#'}
#'
#'@export
#'


XPSCompare <- function(){

#---- calls the function to plot spectra folloowing option selection -----
   CtrlPlot <- function(){
            N.XS <- length(FNamesCoreLines$XPSSample)
            N.CL <- length(FNamesCoreLines$CoreLines)
            if (N.XS == 0) {
                gmessage(msg="Please select XPS Samples!", title="NO XPS SAMPLES SELECTED", icon="warning")
                return()
            }

            if (length(FNamesCoreLines$Ampli) == 0){
               N <- N.XS * N.CL
               FNamesCoreLines$Ampli <<- rep(1,N)
            }
            SelectedNames <- FNamesCoreLines
            SelectedNames$CoreLines <- svalue(T1CLineListCK)
            LL <- length(SelectedNames$CoreLines)
            if (LL == 0) {
                gmessage(msg="Please select The Core Lines to plot!", title="NO CORE-LINES SELECTED", icon="warning")
                return()
            }
            if (length(PanelTitles) == 0) {
               for (ii in 1:LL){
                    PanelTitles <<- c(PanelTitles, SelectedNames$CoreLines[ii]) #List of Titles for the Multipanel
               }
            }
            Plot_Args$PanelTitles <<- PanelTitles
            Limits <- XPScompEngine(PlotParameters, Plot_Args, SelectedNames, Xlim, Ylim)
   }


#----- Update table containing the selected XPSSample names and corelines
   RefreshTab <- function(FNamesCoreLines){
       delete(T1group1, NameTable)
       #If Compare.Coreline is selected likely the number of elements of FNamesCoreLines$XPSSample != FNamesCoreLines$CoreLines
       LL.X <- length(FNamesCoreLines$XPSSample)
       LL.CL <- length(FNamesCoreLines$CoreLines)
       if (LL.X < LL.CL){
           tmp <- rep("  ", (LL.CL-LL.X))
           FNamesCoreLines$XPSSample <- c(FNamesCoreLines$XPSSample, tmp)
       }
       if (LL.CL < LL.X){
           tmp <- rep("  ", (LL.X-LL.CL))
           FNamesCoreLines$CoreLines <- c(FNamesCoreLines$CoreLines, tmp)
       }
       if (LL.X == 0) {
          dummy <- list(XPSSample=c("   ", "  "),CoreLines=c("   ", "  "))   #dummy list to begin: NB each column has 2 initial element otherwise error...
          dummy$XPSSample <- encodeString(dummy$XPSSample, width=40, justify="right")
          dummy$CoreLines <- encodeString(dummy$CoreLines, width=40, justify="right")
          NameTable <<- gtable(dummy, expand=TRUE, fill=TRUE, container=T1group1) #table with the selected FNames
       } else {
          TabNames <- list(XPSSample=c(FNamesCoreLines$XPSSample, " "), CoreLines=c(FNamesCoreLines$CoreLines, " "))
          LL <- length(TabNames$XPSSample)
          for(ii in 1:LL){
             TabNames$XPSSample[ii] <- encodeString(TabNames$XPSSample[ii], width=40, justify="right")
             TabNames$CoreLines[ii] <- encodeString(TabNames$CoreLines[ii], width=40, justify="right")
          }
          NameTable <<- gtable(TabNames, expand=TRUE, container=T1group1) #Table with the list of selected names
       }
       names(NameTable) <<- c("XPSSample", "CoreLines")
   }


#--- Routine for drawing Custom Axis
   CustomAx <- function(CustomDta){
               AxWin <- gwindow(title="CUSTOM AXIS", visible=FALSE)
               AxGroup1 <- ggroup(horizontal=FALSE, container=AxWin)
               txt1="1) Set Axis min, Axis max values: es. min=0, max=10, Nticks=5"
               txt2="2) Set desired Number of Ticks: es. 5"
               glabel(txt1, container=AxGroup1)
               glabel(txt2, container=AxGroup1)
               AxFrame <- gframe("Set Axis Elements", horizontal=FALSE, container=AxGroup1)
               AxLayout <- glayout(homogeneous=FALSE, spacing=3, container=AxFrame)
               axMin <- as.character(round(CustomDta[[1]], 2))
               axMax <- as.character(round(CustomDta[[2]], 2))
               msg1 <- paste("Xmin (min value=", axMin, "):", sep="")
               msg2 <- paste("Xmax (max value=", axMax, "):", sep="")
               if (CustomDta[[3]] == "Y") {
                  msg1 <- paste("Ymin (min value=", axMin, "):", sep="")
                  msg2 <- paste("Ymax (max value=", axMax, "):", sep="")
               }
               AxLayout[1,1] <- EditXmin <- gedit(initial.msg =msg1, width=40, container=AxLayout)
               AxLayout[1,2] <- EditXmax <- gedit(initial.msg =msg2, width=40, container=AxLayout)
               AxLayout[3,1] <- EditNTicks <- gedit(initial.msg ="N.Ticks", container=AxLayout)

               gbutton("     SAVE & EXIT      ", handler=function(h,...){
                        axMin <- as.numeric(svalue(EditXmin))     #X or Y scale min value
                        axMax <- as.numeric(svalue(EditXmax))     #X or Y scale max value
                        axRange <- sort(c(axMin, axMax))          #X or R scale range
                        NTicks <- as.numeric(svalue(EditNTicks))
                        if (is.na(axMin*axMax)) {
                           gmessage("ATTENTION: plase set all the min, max values!", title = "CHANGE X Y RANGE", icon = "error")
                        }
                        if (is.null(NTicks)){
                            gmessage("Please N. Major Ticks  required!", icon="warning")
                        } else {
                            dx <- (axMax-axMin)/NTicks
                            axStp <- seq(from=axMin, to=axMax, by=dx)
                            Ticklabels <- as.character(round(axStp,digits=1))
                            if (CustomDta[[3]] == "X") {
                               if (FName[[SpectIndx]]@Flags) {  #Binding energy set reverse X axis
                                  axRange <- sort(c(axMin, axMax), decreasing=TRUE)
                               } else {
                                  axRange <- sort(c(axMin, axMax))
                               }
                               Plot_Args$scales$x <<- list(at=axStp, labels=Ticklabels)
                               Plot_Args$xlim <<- axRange
                               Xlim <<- axRange
                            } else if (CustomDta[[3]] == "Y") {
                               Plot_Args$scales$y <<- list(at=axStp, labels=Ticklabels)
                               Plot_Args$ylim <<- axRange
                               Ylim <<- axRange
                            }
                            dispose(AxWin)
                            CtrlPlot()
                            Plot_Args$scales$relation <<- "free"
                        }
               }, container = AxFrame)
               visible(AxWin) <- TRUE
   }


   SetLinesPoints <- function(){
         if ( svalue(T3_SetLines) == "OFF" && svalue(T3_SetSymbols) == "OFF") {
            Plot_Args$type <<- " "  #both: line and symbols
            AutoKey_Args$lines <<- FALSE
            AutoKey_Args$points <<- FALSE
            AutoKey_Args$col <<- "white"
            PlotParameters$Colors <<- "white"
            Plot_Args$par.settings$superpose.symbol$col <<- "white"
            Plot_Args$par.settings$superpose.symbol$pch <<- STypeIndx[1]
         }

         if ( svalue(T3_SetLines) == "ON" && svalue(T3_SetSymbols) == "OFF") {
            Plot_Args$type <<- "l"
            AutoKey_Args$lines <<- TRUE
            AutoKey_Args$points <<- FALSE
            AutoKey_Args$col <<- Colors
            PlotParameters$Colors <<- Colors
            Plot_Args$lty <<- LType
            Plot_Args$par.settings$superpose.line$col <<- Colors #Rainbow plot
            Plot_Args$par.settings$superpose.line$lty <<- "solid"
            if (svalue(T3_BW_Col)=="B/W") {
               AutoKey_Args$col <<- "black"
               PlotParameters$Colors <<- "black"
               Plot_Args$par.settings$superpose.line$col <<- "black" #B/W plot
               Plot_Args$par.settings$superpose.line$lty <<- LType
            }
         }

         if ( svalue(T3_SetLines) == "OFF" && svalue(T3_SetSymbols) == "ON") {
            Plot_Args$type <<- "p"  #both: line and symbols
            AutoKey_Args$lines <<- FALSE
            AutoKey_Args$points <<- TRUE
            AutoKey_Args$col <<- Colors
            PlotParameters$Colors <<- Colors
            Plot_Args$pch <<- STypeIndx
            Plot_Args$par.settings$superpose.symbol$col <<- Colors
            Plot_Args$par.settings$superpose.symbol$pch <<- STypeIndx[1]
            if (svalue(T3_BW_Col)=="B/W") {
               AutoKey_Args$col <<- "black"
               PlotParameters$Colors <<- "black"
               Plot_Args$par.settings$superpose.symbol$pch <<- STypeIndx
               Plot_Args$par.settings$superpose.symbol$col <<- "black"
            }
         }

         if ( svalue(T3_SetLines) == "ON" && svalue(T3_SetSymbols) == "ON") {
            Plot_Args$type <<- "b"  #both: line and symbols
            AutoKey_Args$lines <<- TRUE
            AutoKey_Args$points <<- TRUE
            Plot_Args$lty <<- LType
            Plot_Args$pch <<- STypeIndx
            if (svalue(T3_BW_Col)=="B/W") {
               AutoKey_Args$col <<- "black"
               PlotParameters$Colors <<- "black"
               Plot_Args$par.settings$superpose.line$lty <<- LType
               Plot_Args$par.settings$superpose.line$col <<- "black" #B/W plot
               Plot_Args$par.settings$superpose.symbol$pch <<- STypeIndx
               Plot_Args$par.settings$superpose.symbol$col <<- "black"
            } else {
               AutoKey_Args$col <<- Colors
               PlotParameters$Colors <<- Colors
               Plot_Args$par.settings$superpose.line$lty <<- "solid"
               Plot_Args$par.settings$superpose.line$col <<- Colors #Rainbow plot
               Plot_Args$par.settings$superpose.symbol$pch <<- STypeIndx[1]
               Plot_Args$par.settings$superpose.symbol$col <<- Colors
            }
         }
         CtrlPlot()
   }

   CheckCL <- function(){
       CommonCL <- sapply(CLlist[[1]], function(x) unlist(strsplit(x, "\\."))[2] )
       N.CL <- length(CommonCL)                        #N. Corelines in the reference XPS-Sample (first selected sample)
       N.XS <- length(FNamesCoreLines$XPSSample)       #Number of N. XPSSpectra selected
       XpSamp <- get(FNamesCoreLines$XPSSample[1], envir=.GlobalEnv)
       RngX <- RngY <- list(CLName=NULL, min=NULL, max=NULL)
       for(jj in 1:N.CL){ #Range of CoreLines of XPSSample1
          RngX[[1]][jj] <- CommonCL[jj]
          RngX[[2]][jj] <- range(XpSamp[[CommonCL[jj]]]@.Data[1])[1] #First col RngX contains range min value
          RngX[[3]][jj] <- range(XpSamp[[CommonCL[jj]]]@.Data[1])[2] #Second col RngX contains range max value
          RngY[[1]][jj] <- CommonCL[jj]
          RngY[[2]][jj] <- range(XpSamp[[CommonCL[jj]]]@.Data[2])[1]
          RngY[[3]][jj] <- range(XpSamp[[CommonCL[jj]]]@.Data[2])[2]
       }
#Here we load selected XPSSample to see which CoreLines contains.
#The X and Y range OK for all the common corelines of all XPSSamples
       if (N.XS > 1) {
           for (ii in 2:N.XS){
                XpSamp <- get(FNamesCoreLines$XPSSample[ii], envir=.GlobalEnv)
                CoreLines <- names(XpSamp) #set the CL of the first XPSSample as a reference
                for (jj in 1:N.CL){
                    xx <- grep(CommonCL[jj], CoreLines) #is CL string present in CLlist?
                    if (length(xx)==0) {          #pattern CommonCL[jj] not present in CoreLines vector
                        CommonCL <- CommonCL[-jj] #drop the elements non common element
                    }
                    #work out the X-range Y-range common to the selected Core Lines
                    RngXmin <- range(XpSamp[[CommonCL[jj]]]@.Data[1])[1]
                    RngXmax <- range(XpSamp[[CommonCL[jj]]]@.Data[1])[2]
                    if (RngXmin < RngX[[2]][jj]) {RngX[[2]][jj] <- RngXmin}
                    if (RngXmax > RngX[[3]][jj]) {RngX[[3]][jj] <- RngXmax}
                    RngYmin <- range(XpSamp[[CommonCL[jj]]]@.Data[2])[1]
                    RngYmax <- range(XpSamp[[CommonCL[jj]]]@.Data[2])[2]
                    if (RngYmin < RngY[[2]][jj]) {RngY[[2]][jj] <- RngYmin}
                    if (RngYmax > RngY[[3]][jj]) {RngY[[3]][jj] <- RngYmax}
                }
           }
       }
       FNamesCoreLines$CoreLines <<- CommonCL
       Xlim <<- RngX
       Ylim <<- RngY
   }

   CtrlRepCL <- function(ii){    #CTRL for repeated CL: search for core-lines with same name
       LL <- length(CLlist[[ii]])
       jj <- 1
       while(jj < LL){
           SpectName <- CLlist[[ii]][jj]
           SpectName <- unlist(strsplit(SpectName, "\\."))
           SpectName <- SpectName[2]
           Indx <- grep(SpectName, CLlist[[ii]])  #The selected CoreLine name could be in any posiiton in the Destination XPSSample => source Samp Index could be different from Dest Samp index
           if (length(Indx) > 1){                 #The same coreline can be present more than one time
              winCL <- gwindow("SELECT CORELINE", visible=FALSE)
              size(winCL) <- c(100, 200)
              groupCL <- ggroup(horizontal=FALSE, container=winCL)
              N.CL <- length(Indx)
              msg <- paste(" Found ", N.CL," ",SpectName, "spectra.\n Please select the coreline to compare")
              txt <- glabel(text=msg, container=groupCL)
              font(txt) <- list(family="sans",size=12)
              gseparator(horizontal=TRUE, container=groupCL)
              selectCL <- gradio(CLlist[[ii]][Indx], selected=1, horizontal=TRUE, container=groupCL)
              gbutton("    OK     ", handler=function(h, ...){
                                       zz <- as.numeric(svalue(selectCL, index=TRUE))
                                       Indx <<- Indx[-zz] #in Indx remain the component to eliminiate
                                       CLlist[[ii]] <<- CLlist[[ii]][-Indx] #eliminate the repeated spectra
                                       LL <<- length(CLlist[[ii]]) #update CLlist length
                                       dispose(winCL)
              }, container=groupCL)
              visible(winCL) <- TRUE
              winCL$set_modal(TRUE)  #nothing can be done while running this macro
           }                         #modal mode takes the control of the 'return' value which CANNOT be used to return variable values
           jj <- jj+1
       }
   }

#----- reset parameters to the initial values -----
   ResetPlot <- function(){
            svalue(T1FNameListCK) <<- NULL
            CLlist <<- list()
            delete(T1frameCLine, T1CLineListCK)
            T1CLineListCK <<- gcheckboxgroup(CLlist,checked=FALSE, container=T1frameCLine) #at beginning CLlist==list()
            NamesList <<- list(XPSSample=NULL, CoreLines=NULL)
            FNamesCoreLines <<- list(XPSSample=c(" ", " "),CoreLines=c(" "," "), Ampli=NULL)  #dummy lists to begin: NB each lcolumn contains 2 element otherwise error
            RefreshTab(FNamesCoreLines)   #update the table with the name of the selected FNames
            FNamesCoreLines <<- list(XPSSample=NULL,CoreLines=NULL,Ampli=NULL )   #dummy lists to begin: NB each lcolumn contains 2 element otherwise error
            delete(layoutT2, layoutT2[1,2])
            layoutT2[1,2] <<- objFunctAmpli <<- gcombobox(c("        "), selected=-1, editable=FALSE, handler=function(h,...){enabled(objFunctFact) <- TRUE}, container=layoutT2)
            SaveSelection <<- TRUE

            svalue(T1FNameListCK) <<- FALSE
            svalue(objFunctNorm) <<- FALSE
            svalue(objFunctAlign) <<- FALSE
            svalue(objFunctRev) <<- TRUE
            svalue(objFunctSwitch) <<- FALSE
            svalue(objFunctAmpli) <<- -1
            svalue(objFunctFact) <<- ""
            enabled(objFunctFact) <<- FALSE
            svalue(XOffsetobj) <<- 0
            svalue(YOffsetobj) <<- 0
            svalue(xx1) <<- ""
            svalue(xx2) <<- ""
            svalue(yy1) <<- ""
            svalue(yy2) <<- ""
            svalue(T3_BW_Col) <<- "B/W"
            svalue(T3_Grid) <<- "Grid OFF"
            svalue(T3_SetLines) <<- 1
            svalue(T3_SetSymbols) <<- 2
            svalue(T3_LineType) <<- "patterns"
            svalue(T3_LinWidth) <<- 1
            svalue(T3_SymType) <<- "single-symbol"
            svalue(T3_SymSize) <<- 0.8
            svalue(T3_PanStripCol) <<- ""
            svalue(T4_LBTR) <<- "LeftBottom"
            svalue(T4_XScale) <<- "Regular"
            svalue(T4_YScale) <<- "Regular"
            svalue(T4_TitSize) <<- 1.4
            svalue(T4_AxNumSize) <<- 1
            svalue(T4_AxLabSize) <<- 1
            svalue(T4_XAxNameChange) <<- ""
            svalue(T4_YAxNameChange) <<- ""
            svalue(T4_XStep) <<- ""
            svalue(T4_YStep) <<- ""
            svalue(legendCK) <<- TRUE
            svalue(LegColCK) <<- 1
            svalue(TSizeCK) <<- 1
            svalue(LineWdhCK) <<- 1
            svalue(TxtColCK) <<- "B/W"

            XPSSettings <<- get("XPSSettings", envir=.GlobalEnv)
            Colors <<- XPSSettings$Colors
            LType <<- XPSSettings$LType
            SType <<- XPSSettings$Symbols
            STypeIndx <<- XPSSettings$SymIndx
            FitColors <<- c(XPSSettings$BaseColor[1], XPSSettings$ComponentsColor[1], XPSSettings$FitColor[1])
            CLPalette <<- data.frame(Colors=Colors, stringsAsFactors=FALSE)
            FitPalette <<- data.frame(FitColors=FitColors, stringsAsFactors=FALSE)

            PlotParameters <<- DefaultPlotParameters

            Plot_Args <<- list( x=formula("y ~ x"), data=NULL, PanelTitles=list(), groups=NULL,layout=NULL,
                                  xlim=NULL,ylim=NULL,
                                  pch=STypeIndx,cex=1,lty=LType,lwd=1,type="l",
                                  background="transparent",  col="black",
                                  main=list(label=NULL,cex=1.5),
                                  xlab=list(label=NULL, rot=0, cex=1.2),
                                  ylab=list(label=NULL, rot=90, cex=1.2),
                                  zlab=NULL,
                                  scales=list(cex=1, tck=c(1,0), alternating=c(1), relation="free",
                                              x=list(log=FALSE), y=list(log=FALSE), axs="i"),
                                  xscale.components = xscale.components.subticks,
                                  yscale.components = yscale.components.subticks,
                                  las=0,
                                  par.settings = list(superpose.symbol=list(pch=STypeIndx, fill="black"), #set the symbol fill color
                                        superpose.line=list(lty=LType, col="black"), #needed to set the legend colors
                                        par.strip.text=list(cex=1),
                                        strip.background=list(col="grey90") ),
                                  auto.key = TRUE,
                                  grid = FALSE
                             )

            AutoKey_Args <<- list(space="top",
                                  text=" ",
                                  cex = 1,
                                  type= "l",
                                  lines=TRUE,
                                  points=FALSE,
                                  col="black",
                                  columns=1,   #leggends organized in a column
                                  list(corner=NULL,x=NULL,y=NULL)
                             )
            Xlim <<- NULL #reset Xlim
            Ylim <<- NULL #reset Ylim
            plot.new()
   }


#----- Variables -----
   if (is.na(activeFName)){
       gmessage("No data present: please load and XPS Sample", title="XPS SAMPLES MISSING", icon="error")
       return()
   }
   FName <- get(activeFName, envir=.GlobalEnv)   #load the active FName
   ActiveFName <- get("activeFName", envir=.GlobalEnv)  #load the name of the active FNamw (string)
   SpectIndx <- get("activeSpectIndx", envir=.GlobalEnv)#index of the active coreline
   SpectList <- XPSSpectList(ActiveFName)   #list of the active XPSSample core lines
   NComp=length(FName[[SpectIndx]]@Components)
   NCorelines <- NULL
   FitComp1 <- ""
   for (ii in 1:NComp){
      FitComp1[ii] <- paste("C",ii, sep="")
   }
   FNameListTot <- as.array(XPSFNameList())     #List of all XPSSamples loaded in the workspace
   LL=length(FNameListTot)
   jj <- 1
   FNamesCoreLines <- list(XPSSample=NULL, CoreLines=NULL, Ampli=NULL)
   NamesList <- list(XPSSample=NULL, CoreLines=NULL)
   CLlist <- list()
   SpectName <- ""

#--- list of graphical variables
   PatternList <- NULL
   FontSize <- c(0.6,0.8,1,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6,2.8,3)
   AxLabOrient <- c("Horizontal", "Rot-20", "Rot-45", "Rot-70", "Vertical")
   XPSSettings <- get("XPSSettings", envir=.GlobalEnv)
   Colors <- XPSSettings$Colors
   LType <- XPSSettings$LType
   SType <- XPSSettings$Symbols
   STypeIndx <- XPSSettings$SymIndx
   CLPalette <- as.matrix(Colors)
   CLPalette <- data.frame(Colors=CLPalette, stringsAsFactors=FALSE)
   FitColors <- as.matrix(c(XPSSettings$BaseColor[1], XPSSettings$ComponentsColor[1], XPSSettings$FitColor[1]))
   VarNames <- c("BasLnCol", "CompCol", "FitCol")
   FitPalette <- data.frame(Object=VarNames, Colors=FitColors, stringsAsFactors=FALSE)
#-------------------------------------------------------------------------------------------------
#   LType <- c("solid", "dashed", "dotted", "dotdash", "longdash",     #definisco 20 tipi divesi di line pattern
#            "twodash", "F8", "431313", "22848222", "12126262",
#            "12121262", "12626262", "52721272", "B454B222", "F313F313",
#            "71717313", "93213321", "66116611", "23111111", "222222A2" )
#
#   SType <- c("VoidCircle", "VoidSquare", "VoidTriangleUp", "VoidTriangleDwn",  "Diamond",
#            "X", "Star", "CrossSquare", "CrossCircle", "CrossDiamond",
#            "SolidSquare", "SolidCircle", "SolidTriangleUp", "SolidTriangleDwn", "SolidDiamond",
#            "DavidStar", "SquareCross", "SquareTriang", "CircleCross", "Cross")
#   STypeIndx <- c(1,  0,  2,  6,  5,
#                4,  8,  7,  10, 9,
#                15, 16, 17, 25, 18,
#                11, 12, 14, 13, 3)
#
#   Colors <- c("black", "red", "limegreen", "blue", "magenta", "orange", "cadetblue", "sienna",
#             "darkgrey", "forestgreen", "gold", "darkviolet", "greenyellow", "cyan", "lightblue",
#             "turquoise", "deeppink3", "wheat", "thistle", "grey40")
#-------------------------------------------------------------------------------------------------
   LWidth <- c(1,1.25,1.5,1.75,2,2.25,2.5,3, 3.5,4)
   SymSize <- c(0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2) #lattice prende indici simboli piuttosto che nomesimbolo
   PanelTitles <- NULL
   LegPos <- c("OutsideTop","OutsideRight","OutsideLeft", "OutsideBottom",
             "InsideTopRight","InsideTopLeft","InsideBottomRight","InsideBottomLeft")
   LegOrient <- c("Vertical", "Horizontal")
   LegLineWdh <- c(1,1.5,2,2.5,3,3.5,4,4.5,5)
   LegTxtCol <- c("RainBow", "Black")
   LegTxtSize <- c(0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6,2.8,3)
   LegDist <- c(0,0.01,0.02,0.04,0.08,0.1,0.12,0.14,0.16,0.18,0.2)
   ColorList <- NULL
   exit <- NULL
   Xlim <- NULL
   Ylim <- NULL

#--- general options
   PlotParameters <- list()
   PlotParameters$Aligne <- FALSE
   PlotParameters$RTFLtd <- FALSE #restrict plot to RTF region
   PlotParameters$Normalize <- FALSE
   PlotParameters$Reverse <- TRUE #reversed X axes for Bind. Energy
   PlotParameters$SwitchE <- FALSE
   PlotParameters$XOffset <- 0
   PlotParameters$YOffset <- 0
   PlotParameters$CustomXY <- NULL
   PlotParameters$OverlayType <- "Compare.CoreLines" #Compare.Corelines  and  Multi-Panel are fixed options
   PlotParameters$OverlayMode <- "Multi-Panel"
   PlotParameters$Colors <- "B/W"
#--- legend options
   PlotParameters$Labels <- NULL
   PlotParameters$Legenda <- FALSE
   PlotParameters$LegPos <- "left"  #Out side left legend position
   PlotParameters$LegLineWdh <- 1
   PlotParameters$LegTxtCol <- "RainBow"
   PlotParameters$LegTxtSize <- 1
   PlotParameters$LegDist <- 0

   DefaultPlotParameters <- PlotParameters

#--- comands for Lattice options
   Plot_Args <- list( x=formula("y ~ x"), data=NULL, PanelTitles=list(), groups=NULL,layout=NULL,
                    xlim=NULL, ylim=NULL,
                    pch=STypeIndx,cex=1,lty=LType,lwd=1,type="l",
                    background="transparent", col="black",
                    main=list(label=NULL,cex=1.5),
                    xlab=list(label=NULL, rot=0, cex=1.2),
                    ylab=list(label=NULL, rot=90, cex=1.2),
                    zlab=NULL,
                    scales=list(cex=1, tck=c(1,0), alternating=c(1), tick.number=NULL, relation="free",
                                x=list(log=FALSE), y=list(log=FALSE), axs="i"),
                    xscale.components = xscale.components.subticks,
                    yscale.components = yscale.components.subticks,
                    las=0,
                    par.settings = list(superpose.symbol=list(pch=STypeIndx,fill="black"), #set symbol filling color
                                        superpose.line=list(lty=LType, col="black"), #needed to set legend colors
                                        par.strip.text=list(cex=1),
                                        strip.background=list(col="grey90") ),
                    auto.key = TRUE,
                    grid = FALSE
                  )


   AutoKey_Args <- list( space="top",
                         text=" ",
                         cex = 1,
                         type= "l",
                         lines=TRUE,
                         points=FALSE,
                         col="black",
                         columns=1,  #legends organized in a column
                         list(corner=NULL,x=NULL,y=NULL)
                       )

   SaveSelection <- TRUE #at beginning force the control of the selection to TRUE to avoid error messages

#--- Reset graphical window
   plot.new()
   assign("MatPlotMode", FALSE, envir=.GlobalEnv)  #basic matplot function used to plot data


#===== NoteBook =====

   win <- gwindow(" COMPARE SPECTRA ", parent=c(0,50), visible=FALSE)
   size(win) <- c(400,400)
   maingroup <- ggroup(horizontal=FALSE, container=win)
   nb <- gnotebook(expand=TRUE, container = maingroup)

# --- TAB1 ---
#XPS Sample/Coreline selection

     T1group1 <- ggroup(label="XPS SAMPLE SELECTION", spacing=5, horizontal=FALSE, container=nb)
     layoutT1 <- glayout(homogeneous=FALSE, spacing=5, container=T1group1)



     layoutT1[1,1] <- T1frameButtT1 <- gframe(text="PLOT", spacing=5, container=layoutT1)
     T1groupButtons <- ggroup(horizontal=FALSE, container = T1frameButtT1)
     gbutton("PLOT", handler=function(h,...){
                            CtrlPlot() #plot selected XPS-SAmples
                   }, container=T1groupButtons)

     gbutton("RESET PLOT", handler=function(h,...){
                           ResetPlot()
     }, container=T1groupButtons)

     gbutton("UPDATE XPS-SAMPLE LIST", handler=function(h,...){
                           svalue(T1FNameListCK) <<- NULL
                           FName <<- get(activeFName, envir=.GlobalEnv)
                           ActiveFName <<- get("activeFName", envir=.GlobalEnv)
                           SpectIndx <<- get("activeSpectIndx", envir=.GlobalEnv)
                           SpectList <<- XPSSpectList(ActiveFName)   #sCoreLine list of the XPSSample
                           NComp <<- length(FName[[SpectIndx]]@Components)
                           NCorelines <<- NULL
                           FitComp1 <<- ""  #build vector containing names of the fit components on the Active Spectrum
                           for (ii in 1:NComp){
                               FitComp1[ii] <- paste("C",ii, sep="")
                           }
                           FNameListTot <- as.array(XPSFNameList())     #list of all XPSSample in Envir=.GlobalEnv
                           LL=length(FNameListTot)
                           jj <- 1
                           FNamesCoreLines <<- list(XPSSample=NULL, CoreLines=NULL, Ampli=NULL)
                           NamesList <<- list(XPSSample=NULL, CoreLines=NULL)
                           SaveSelection <<- TRUE

                           delete(T1frameFName, T1FNameListCK)       #update panel
                           T1FNameListCK <<- gcheckboxgroup(FNameListTot,checked=FALSE, handler=function(h,...){
                                                  FNamesCoreLines$XPSSample <<- svalue(T1FNameListCK)
                                                  FNamesCoreLines$CoreLines <<- list()
                                                  LL <- length(FNamesCoreLines$XPSSample)
                                                  if (LL > 3) LL <- 3
                                                  Plot_Args$auto.key$columns <<- LL
                                                  AutoKey_Args$columns <<- LL
                                                  Plot_Args$auto.key$text <<- unlist(FNamesCoreLines$XPSSample)
                                                  for (ii in 1:LL){
                                                      CLlist[[ii]] <<- XPSSpectList(FNamesCoreLines$XPSSample[ii])
                                                  }
                                                  CheckCL()   #checks for corelines common to selected XPSSamples
                                                  delete(T1frameCLine,T1CLineListCK)
                                                  T1CLineListCK <<- gcheckboxgroup(FNamesCoreLines$CoreLines, container=T1frameCLine) #refresh gcheckboxgroup for coreline selection
                                                  RefreshTab(FNamesCoreLines)   #update the table with the name of the selected FNames
                                                  delete(T2CLgroup, CLPanel)
                                                  CLPanel <<- gcombobox(FNamesCoreLines$CoreLines, selected = -1, editable = FALSE, container=T2CLgroup) #refresh combobox for custom XY scale
                           }, container=T1frameFName)
                           ResetPlot()
                           plot.new()
     }, container=T1groupButtons)

     gbutton("  EXIT  ", handler=function(h,...){
                     dispose(win)
                   }, container=T1groupButtons)

     layoutT1[1,2] <- T1frameFName <- gframe(text="SELECT XPS-SAMPLE", spacing=5, container=layoutT1)
     T1FNameListCK <- gcheckboxgroup(FNameListTot,checked=FALSE, handler=function(h,...){
                            FNamesCoreLines$XPSSample <<- svalue(T1FNameListCK)
                            if (length(FNamesCoreLines$XPSSample) == 0 ){
                               FNamesCoreLines$XPSSample <<- list()  #checkbox deselection
                               FNamesCoreLines$CoreLines <<- list()
                               Plot_Args$PanelTitles <<- list()
                               CLlist <<- list()
                            } else {
                               FNamesCoreLines$CoreLines <<- list()
                               LL <- length(FNamesCoreLines$XPSSample)
                               if (LL > 3) {
                                   Plot_Args$auto.key$columns <<- 3
                                   AutoKey_Args$columns <<- 3
                               }
                               Plot_Args$auto.key$text <<- unlist(FNamesCoreLines$XPSSample)
#                              Define only the CLlist for the new selected XPSSample
                               CLlist[[LL]] <<- XPSSpectList(FNamesCoreLines$XPSSample[LL])

                               CtrlRepCL(LL)  #controls if same spectra are repeated in CLlist[[LL]]
                               CheckCL()   #check for corelines common to selected XPSSamples
                            }
                            delete(T1frameCLine,T1CLineListCK)
                            T1CLineListCK <<- gcheckboxgroup(FNamesCoreLines$CoreLines, container=T1frameCLine) #refresh gcheckboxgroup for coreline selection
                            RefreshTab(FNamesCoreLines)   #update the table with the name of the selected FNames
                            delete(T2CLgroup, CLPanel)
                            CLPanel <<- gcombobox(FNamesCoreLines$CoreLines, selected = -1, editable = FALSE, container=T2CLgroup) #refresh combobox for custom XY scale
                   }, container=T1frameFName)

     layoutT1[1,3] <- T1frameCLine <- gframe(text="SELECT SPECTRA", spacing=5, container=layoutT1)
     T1CLineListCK <- gcheckboxgroup(CLlist,checked=FALSE, container=T1frameCLine) #at beginning CLlist==list()

     glabel("SELECTED XPS-SAMPLES AND COMMON CORE LINES", container=T1group1)
     dummy <- list(XPSSample=c("   ", "  "),CoreLines=c("   ", "  "))   #dummy list to begin: NB each column has 2 initial element otherwise error...
     dummy$XPSSample <- encodeString(dummy$XPSSample, width=40, justify="right")
     dummy$CoreLines <- encodeString(dummy$CoreLines, width=40, justify="right")
     NameTable <<- gtable(dummy, expand=TRUE, fill=TRUE, container=T1group1) #table with the selected FNames
     names(NameTable) <<- c("XPSSample", "CoreLines")


# --- TAB2 ---

###Funct1: NORMALIZE

   T2group1 <- ggroup(label="PLOT FUNCTIONS",horizontal=FALSE, container=nb)

   T2frame1 <- gframe(" FUNCTIONS ", horizontal=FALSE, spacing=5, container=T2group1)
   T2group2 <- ggroup(horizontal=TRUE, container=T2frame1)

   objFunctNorm <- gcheckbox("Normalize",checked=FALSE, handler=function(h,...){
                    PlotParameters$Normalize <<- svalue(objFunctNorm)
                    FName <- get(FNamesCoreLines$XPSSample[1], envir=.GlobalEnv) #retrieve a generic XPSSample from the selected ones
                    SpectName <- unlist(strsplit(FNamesCoreLines$CoreLines[1], "\\."))  #retrieve a generic coreline from the list of selected ones
                    indx <- as.numeric(SpectName[1])
                    Plot_Args$ylab$label <<- FName[[indx]]@units[2]   #retrieve the Y axis label
                    if ( svalue(objFunctNorm)) {   #Normalize option TRUE
                       Plot_Args$ylab$label <<- "Intensity [a.u.]"
                    }
                    CtrlPlot() ####PLOT FOLLOWING SELECTIONS
                 }, container=T2group2)

###Funct2: Y-Aligne

   objFunctAlign <- gcheckbox("Aligne bkg to 0",checked=FALSE, handler=function(h,...){
                    PlotParameters$Aligne <<- svalue(objFunctAlign)
                    CtrlPlot() ####PLOT FOLLOWING SELECTIONS
                 }, container=T2group2)

###Funct3: Reverse X axis

   objFunctRev <- gcheckbox("Reverse X axis",checked=TRUE, handler=function(h,...){
                    PlotParameters$Reverse <<- svalue(objFunctRev)
                    CtrlPlot() ####PLOT FOLLOWING SELECTIONS
                 }, container=T2group2)

###Funct4: Switch Binding to Kinetic Energy scale

   objFunctSwitch <- gcheckbox("Switch BE/KE energy scale",checked=FALSE, handler=function(h,...){
                    PlotParameters$SwitchE <<- svalue(objFunctSwitch)
                    CtrlPlot() ####PLOT FOLLOWING SELECTIONS
                 }, container=T2group2)

###Funct5: Amplify

   layoutT2 <- glayout(homogeneous=FALSE, spacing=5, container=T2frame1)
   layoutT2[1,1] <- glabel("XPSSamp.", container=layoutT2)
   layoutT2[1,2] <- objFunctAmpli <- gcombobox(c("   "), selected=-1, editable=FALSE, handler=function(h,...){
                    enabled(objFunctFact) <- TRUE
                 }, container=layoutT2)

   layoutT2[1,3] <- glabel("ScaleFact.", container=layoutT2)
   layoutT2[1,4] <- objFunctFact <- gedit("", handler=function(h,...){
                    indx <- as.numeric(svalue(objFunctAmpli, index=TRUE))
                    FNamesCoreLines$Ampli[indx] <<- as.numeric(svalue(objFunctFact))
                    CtrlPlot() ####PLOT FOLLOWING SELECTIONS
                 }, container=layoutT2)


###Funct6: X, Y offset
   layoutT2[2,1] <- glabel("X-Offset", container=layoutT2)
   layoutT2[2,2] <- XOffsetobj <- gedit("", initial.msg = "X_Off= ", container=layoutT2)
   addHandlerChanged(XOffsetobj, handler=function(h,...){
                        PlotParameters$XOffset <<- as.numeric(svalue(XOffsetobj))
                        CtrlPlot()
                    })
   layoutT2[2,3] <- glabel("Y-Offset", container=layoutT2)
   layoutT2[2,4] <- YOffsetobj <- gedit("", initial.msg = "Y_Off= ",container=layoutT2)
   addHandlerChanged(YOffsetobj, handler=function(h,...){
                        PlotParameters$YOffset <<- as.numeric(svalue(YOffsetobj))
                        CtrlPlot()
                    })

###Funct8: Custom XY scale
   T2frame2 <- gframe(text="EXACT X, Y RANGE", horizontal=FALSE, spacing=5, container=T2group1)
   glabel(text="Select the Spectrum", spacing=3, container=T2frame2)
   T2CLgroup <- ggroup(horizontal=TRUE, container=T2frame2)
   CLPanel <- gcombobox(" ", selected = -1, editable = FALSE, container=T2CLgroup)

   T2XYgroup <- ggroup(horizontal=TRUE, container=T2frame2)
   xx1 <- gedit("", initial.msg = "Xmin= ", container=T2XYgroup)
   xx2 <- gedit("", initial.msg = "Xmax= ", container=T2XYgroup)
   yy1 <- gedit("", initial.msg = "Ymin= ", container=T2XYgroup)
   yy2 <- gedit("", initial.msg = "Ymax= ", container=T2XYgroup)
   tkconfigure(xx1$widget, width=15)
   tkconfigure(xx2$widget, width=15)
   tkconfigure(yy1$widget, width=15)
   tkconfigure(yy2$widget, width=15)


   T2_ButtGroup <- ggroup(horizontal=TRUE, container=T2frame2)  #needed only to have a small OK button
   gbutton("  OK  ", width=25, handler=function(h,...){
                   panel <- svalue(CLPanel, index=TRUE)
                   xx1 <- as.numeric(svalue(xx1))
                   xx2 <- as.numeric(svalue(xx2))
                   yy1 <- as.numeric(svalue(yy1))
                   yy2 <- as.numeric(svalue(yy2))
                   if (is.na(xx1*xx2*yy1*yy2)) {
                       gmessage("ATTENTION: plase set all the xmin, xmax, ymin, ymax values!", title = "CHANGE X Y RANGE", icon = "error")
                   }
                   PlotParameters$CustomXY <<- c(panel, xx1, xx2, yy1, yy2)
                   CtrlPlot() 
                 }, container = T2_ButtGroup)

   gbutton("  RESET PLOT  ", handler=function(h,...){
                   ResetPlot()
                   CtrlPlot()
                 }, container = T2_ButtGroup)

   gbutton(" EXIT ", handler=function(h,...){
				       dispose(win)
                 }, container = T2_ButtGroup)


# --- TAB3 ---

# Rendering options
   T3group1 <- ggroup(label="RENDERING", horizontal=FALSE, container=nb)
   T3group2 <- ggroup(horizontal=TRUE, container=T3group1)
   T3group3 <- ggroup(horizontal=FALSE, container=T3group2)

   T3F_CL_Colors <- gframe("SET CORELINE PALETTE", spacing=5, container=T3group3)
   T3_CL_Colors <- gdf(CLPalette, container=T3F_CL_Colors)
   size(T3_CL_Colors) <- c(110, 250)
   addHandlerChanged(T3_CL_Colors, handler=function(h,...){   #edit Palette preferences
                             CLPalette$Colors <<- Colors <<- h$obj[]
                             PlotParameters$Colors <<- Colors
                             Plot_Args$par.settings$superpose.symbol$fill <<- Colors
                             Plot_Args$par.settings$superpose.line$col <<- Colors
                             AutoKey_Args$col <<- Colors
                             Plot_Args$par.settings$superpose.symbol$col <<- Colors
                             CtrlPlot()
                       } )

   layoutRend <- glayout(homogeneous=FALSE, spacing=3, container=T3group2)

   layoutRend[1,1] <- T3F_BW_Col <- gframe("COLOR", spacing=5, container=layoutRend)
   T3_BW_Col <- gcombobox(c("B/W", "RainBow"), selected=1, editable=FALSE, handler=function(h,...){
                             if(svalue(T3_BW_Col)=="B/W") {
                                svalue(T3_LineType) <<- "patterns"
                                svalue(TxtColCK) <<- "B/W"
                                PlotParameters$Colors <<- "black"
                                Plot_Args$lty <<- LType
                                Plot_Args$pch <<- STypeIndx
                                if (length(svalue(T3_LineType))==0) svalue(T3_SymType) <- "multi-symbols"
                                if (length(svalue(T3_SymType))==0) svalue(T3_LineType) <- "patterns"
                                Plot_Args$par.settings$superpose.symbol$col <<- "black"
                                Plot_Args$par.settings$superpose.symbol$pch <<- STypeIndx
                                Plot_Args$par.settings$superpose.line$col <<- "black"
                                Plot_Args$par.settings$superpose.line$lty <<- LType
                                Plot_Args$par.settings$strip.background$col <<- "grey90"
                                AutoKey_Args$col <<- "black"
                             } else {
                                svalue(T3_LineType) <<- "solid"
                                svalue(TxtColCK) <<- "RainBow"
                                Plot_Args$lty <<- "solid"
                                Plot_Args$pch <<- STypeIndx[1]
                                PlotParameters$Colors <<- Colors
                                Plot_Args$par.settings$superpose.symbol$fill <<- Colors
                                Plot_Args$par.settings$superpose.line$col <<- Colors
                                Plot_Args$par.settings$superpose.line$lty <<- "solid"
                                Plot_Args$par.settings$strip.background$col <<- "lightskyblue1"
                                AutoKey_Args$col <<- Colors
                             }
                             CtrlPlot() }, container=T3F_BW_Col)


   layoutRend[1,2] <- T3F_Grid <- gframe("GRID", spacing=5, container=layoutRend)
   T3_Grid <- gcombobox(c("Grid ON", "Grid OFF"), selected=-1, editable=FALSE, handler=function(h,...){
                             if(svalue(T3_Grid)=="Grid ON") {
                                Plot_Args$grid <<- TRUE
                             } else {
                                Plot_Args$grid <<- FALSE
                             }
                             CtrlPlot() }, container=T3F_Grid)

   layoutRend[2,1] <- T3F_SetLines <- gframe("SET LINES", spacing=5, container=layoutRend)
   T3_SetLines <- gradio(c("ON", "OFF"), selected=1, horizontal = TRUE, handler=function(h,...){
                               SetLinesPoints()
                           }, container=T3F_SetLines)

   layoutRend[2,2] <- T3F_SetSymbols <- gframe("SET SYMBOLS", horizontal=TRUE, spacing=5, container=layoutRend)
   T3_SetSymbols <- gradio(c("ON", "OFF"), selected=2, horizontal=TRUE, handler=function(h,...){
                               SetLinesPoints()
                            }, container=T3F_SetSymbols)

   layoutRend[3,1] <- T3F_SetLines <- gframe("LINE TYPE", spacing=5, container=layoutRend)
   T3_LineType <- gcombobox(c("solid", "patterns"), selected=2, editable=FALSE, handler=function(h,...){
                             Plot_Args$type <<- "l"
                             palette <- svalue(T3_BW_Col)
                             if (svalue(T3_LineType)=="solid") {
                                svalue(T3_BW_Col) <<- "RainBow"
                                Plot_Args$lty <<- "solid"
                                Plot_Args$pch <<- STypeIndx[1]
                                PlotParameters$Colors <<- Colors
                                Plot_Args$par.settings$superpose.symbol$fill <<- Colors
                                Plot_Args$par.settings$superpose.line$col <<- Colors
                                Plot_Args$par.settings$superpose.line$lty <<- "solid"
                                Plot_Args$par.settings$strip.background$col <<- "lightskyblue"
                                AutoKey_Args$col <<- Colors
                             }
                             if (svalue(T3_LineType)=="patterns") {
                                svalue(T3_BW_Col) <<- "B/W"
                                PlotParameters$Colors <<- "black"
                                Plot_Args$lty <<- LType
                                Plot_Args$pch <<- STypeIndx
                                Plot_Args$par.settings$superpose.symbol$col <<- "black"
                                Plot_Args$par.settings$superpose.symbol$pch <<- STypeIndx
                                Plot_Args$par.settings$superpose.line$col <<- "black"
                                Plot_Args$par.settings$superpose.line$lty <<- LType
                                Plot_Args$par.settings$strip.background$col <<- "gray90"
                                AutoKey_Args$col <<- "black"
                             }
                             CtrlPlot()
                           }, container=T3F_SetLines)

   layoutRend[3,2] <- T3F_LinWidth <- gframe("LINE WIDTH", spacing=5, container=layoutRend)
   T3_LinWidth <- gcombobox(LWidth, selected=1, editable=FALSE, handler= function(h,...){
                              Plot_Args$lwd <<- as.numeric(svalue(T3_LinWidth))
                              CtrlPlot()
                           }, container=T3F_LinWidth)


   layoutRend[4,1] <- T3F_SetSymbols <- gframe("SYMBOLS", spacing=5, container=layoutRend)
   T3_SymType <- gcombobox(c("single-symbol", "multi-symbols"), selected=2, editable=FALSE, handler=function(h,...){
                              if (svalue(T3_SymType)=="single-symbol") {
                                 svalue(T3_BW_Col) <<- "RainBow"
                                 Plot_Args$lty <<- "solid"
                                 Plot_Args$pch <<- STypeIndx[1]
                                 PlotParameters$Colors <<- Colors
                                 Plot_Args$par.settings$superpose.symbol$fill <<- Colors
                                 Plot_Args$par.settings$superpose.line$col <<- Colors
                                 Plot_Args$par.settings$superpose.line$lty <<- "solid"
                                 Plot_Args$par.settings$strip.background$col <<- "lightskyblue"
                                 AutoKey_Args$col <<- Colors
                              }
                              if (svalue(T3_SymType)=="multi-symbols") {
                                 svalue(T3_BW_Col) <<- "B/W"
                                 PlotParameters$Colors <<- "black"
                                 Plot_Args$lty <<- LType
                                 Plot_Args$pch <<- STypeIndx
                                 Plot_Args$par.settings$superpose.symbol$col <<- "black"
                                 Plot_Args$par.settings$superpose.symbol$pch <<- STypeIndx
                                 Plot_Args$par.settings$superpose.line$col <<- "black"
                                 Plot_Args$par.settings$superpose.line$lty <<- LType
                                 Plot_Args$par.settings$strip.background$col <<- "gray90"
                                 AutoKey_Args$col <<- "black"
                              }
                              CtrlPlot()
                            }, container=T3F_SetSymbols)

   layoutRend[4,2] <- T3F_SymSize <- gframe("SYMSIZE", spacing=5, container=layoutRend)
   T3_SymSize <- gcombobox(SymSize, selected=4, editable=FALSE, handler= function(h,...){
                              Plot_Args$cex <<- as.numeric(svalue(T3_SymSize))
                              CtrlPlot()
                            }, container=T3F_SymSize)

   layoutRend[5,1] <- T3F_PanStripCol <- gframe("PANEL STRIP COLOR", spacing=5, container=layoutRend)
   T3_PanStripCol <- gcombobox(c("white","grey", "darkgrey","lightblue","blue","darkblue","deepskyblue","lightbeige","beige","darkbeige","lightpink","pink","darkpink","lightgreen","green","darkgreen"), selected=-1, editable=FALSE, handler= function(h,...){
                             StripCol <- svalue(T3_PanStripCol)
                             if(StripCol=="grey")               { StripCol <- "grey90"
                             } else if (StripCol=="darkgrey")   { StripCol <- "gray60"

                             } else if (StripCol=="lightblue")  { StripCol <- "lightskyblue1"
                             } else if(StripCol=="blue")        { StripCol <- "lightskyblue3"
                             } else if(StripCol=="darkblue")    { StripCol <- "steelblue3"

                             } else if (StripCol=="lightbeige") { StripCol <- "beige"
                             } else if(StripCol=="beige")       { StripCol <- "bisque2"
                             } else if(StripCol=="darkbeige")   { StripCol <- "navajowhite4"

                             } else if (StripCol=="pink")       { StripCol <- "lightpink2"
                             } else if(StripCol=="darkpink")    { StripCol <- "lightpink4"

                             } else if (StripCol=="lightgreen") { StripCol <- "darkseagreen1"
                             } else if(StripCol=="green")       { StripCol <- "darkseagreen2"
                             } else if(StripCol=="darkgreen")    { StripCol <- "mediumseagreen"
                             }
                             Plot_Args$par.settings$strip.background$col <<- StripCol
                             CtrlPlot() }, container=T3F_PanStripCol)

   gbutton(" RESET PLOT ", handler=function(h,...){
                             ResetPlot()
                             CtrlPlot()
                            }, container=T3group1)

   gbutton(" EXIT ", handler=function(h,...){
				                  dispose(win)
                       }, container = T3group1)


# --- TAB4 ---

# Axis Rendering options

   T4group1 <- ggroup(label="AXES", horizontal=FALSE, container=nb)
   layoutAxis <- glayout(homogeneous=FALSE, spacing=3, container=T4group1)

   layoutAxis[1,1] <- T4F_LBTR <- gframe("TICKS", spacing=5, container=layoutAxis)
   T4_LBTR <- gcombobox(c("LeftBottom", "TopRight", "Both"), selected=1, editable=FALSE, handler= function(h,...){
                             if (svalue(T4_LBTR,index=TRUE)==1) {
                                Plot_Args$scales$tck <<- c(1,0)
                                Plot_Args$scales$alternating <<- c(1)
                             } else if (svalue(T4_LBTR,index=TRUE)==2) {
                                Plot_Args$scales$tck <<- c(0,1)
                                Plot_Args$scales$alternating <<- c(2)
                             } else if (svalue(T4_LBTR,index=TRUE)==3) {
                                Plot_Args$scales$tck <<- c(1,1)
                                Plot_Args$scales$alternating <<- c(3)
                             }
                             CtrlPlot()
                             }, container=T4F_LBTR)


   layoutAxis[1,2] <- T4F_XScale <- gframe("X SCALE", spacing=5, container=layoutAxis)
   T4_XScale <- gcombobox(c("Regular", "Power", "Log.10", "Log.e"), selected=1, editable=FALSE, handler= function(h,...){
                             if (svalue(T4_XScale,index=TRUE)==1) {
                                Plot_Args$scales$x$log <<- FALSE
                                Plot_Args$xscale.components <<- xscale.components.subticks
                             } else if (svalue(T4_XScale,index=TRUE)==2) {
                                Plot_Args$scales$X$log <<- 10    # 10^ power scale
                                Plot_Args$xscale.components <<- xscale.components.logpower
                             } else if (svalue(T4_XScale,index=TRUE)==3) {
                                Plot_Args$scales$X$log <<- 10    # log10 scale
                                Plot_Args$xscale.components <<- xscale.components.log10ticks
                             } else if (svalue(T4_XScale,index=TRUE)==4) {
                                Plot_Args$scales$X$log <<- "e"   # log e scale
                                Plot_Args$xscale.components <<- xscale.components.subticks
                             }
                             CtrlPlot() }, container=T4F_XScale)

   layoutAxis[1,3] <- T4F_YScale <- gframe("Y SCALE", spacing=5, container=layoutAxis)
   T4_YScale <- gcombobox(c("Regular", "Power", "Log.10", "Log.e"), selected=1, editable=FALSE, handler= function(h,...){
                             if (svalue(T4_YScale,index=TRUE)==1) {
                                Plot_Args$scales$y$log <<- FALSE
                                Plot_Args$yscale.components <<- yscale.components.subticks
                             } else if (svalue(T4_YScale,index=TRUE)==2) {
                                Plot_Args$scales$y$log <<- 10
                                Plot_Args$yscale.components <<- yscale.components.logpower
                             } else if (svalue(T4_YScale,index=TRUE)==3) {
                                Plot_Args$scales$y$log <<- 10
                                Plot_Args$yscale.components <<- yscale.components.log10ticks
                             } else if (svalue(T4_YScale,index=TRUE)==4) {
                                Plot_Args$scales$y$log <<- "e"
                                Plot_Args$yscale.components <<- yscale.components.subticks
                             }
                             CtrlPlot() }, container=T4F_YScale)

   layoutAxis[2,1] <- T4F_TitSize <- gframe("TITLE SIZE", spacing=5, container=layoutAxis)
   T4_TitSize <- gcombobox(FontSize, selected=5, editable=FALSE, handler= function(h,...){
                             if (PlotParameters$OverlayMode=="Single-Panel" || PlotParameters$OverlayMode=="TreD") {
                                 Plot_Args$main$cex <<- svalue(T4_TitSize)
                             } else if (PlotParameters$OverlayMode=="Multi-Panel") {
                                 Plot_Args$par.strip.text$cex <<- as.numeric(svalue(T4_TitSize))
                             }
                             CtrlPlot() }, container=T4F_TitSize)

   layoutAxis[2,2] <- T4F_PanelTitles <- gframe("CHANGE MULTI-PANEL TITLES", spacing=5, container=layoutAxis)
   T4_PanelTitles <- gbutton(text="Change Titles", spacing=5, handler=function(h,...){
                                TitleWin <- gwindow(title="MultiPanel Labels", visible=FALSE) #open a new window to contain a gdf() to change the titles of the panels
                                TitleGroup <- ggroup(horizontal=FALSE, container=TitleWin)
                                glabel("                     EDIT TITLES                           ", container=TitleGroup) #long lable to obtain a reasonable window dimension

                                LL=length(PanelTitles)
                                PTitles <- data.frame(TITLES=PanelTitles, stringsAsFactors=FALSE)
                                TitleDFrame <- gdf(items=PTitles, container=TitleGroup)  #here no handler it does not work in linux
                                size(TitleDFrame) <- c(100,200)   #size needed to obtain a non-null size for the gdf()
                                addHandlerChanged(TitleDFrame, handler=function(h,...){  #addHandlerChanged to add a handler to gdf() working also in linux
                                      PanelTitles <<- h$obj[]
                                })
                                gbutton("     SAVE TITLES AND EXIT      ", handler=function(h,...){
                                      Plot_Args$PanelTitles <<- PanelTitles
                                      dispose(TitleWin)
                                      CtrlPlot()
                                }, container = TitleGroup)
                                visible(TitleWin) <- TRUE
                       }, container=T4F_PanelTitles)



   layoutAxis[3,1] <- T4F_AxNumSize <- gframe("AXIS NUMBER SIZE", spacing=5, container=layoutAxis)
   T4_AxNumSize <- gcombobox(FontSize, selected=3, editable=FALSE, handler= function(h,...){
                             Plot_Args$scales$cex <<- svalue(T4_AxNumSize)
                             CtrlPlot() }, container=T4F_AxNumSize)

   layoutAxis[3,2] <- T4F_AxLabSize <- gframe("AXIS LABEL SIZE", spacing=5, container=layoutAxis)
   T4_AxLabSize <- gcombobox(FontSize, selected=3, editable=FALSE, handler= function(h,...){
                             Plot_Args$xlab$cex <<- svalue(T4_AxLabSize)
                             Plot_Args$ylab$cex <<- svalue(T4_AxLabSize)
                             CtrlPlot() }, container=T4F_AxLabSize)
                             
   layoutAxis[3,3] <- T4F_AxLabOrient <- gframe("AXIS NHUMBERg ORIENTATION", spacing=5, container=layoutAxis)
   T4_AxLabOrient <- gcombobox(AxLabOrient, selected=1, editable=FALSE, handler= function(h,...){
                             LabOrient <- svalue(T4_AxLabOrient)
                             if (LabOrient == "Horizontal"){LabOrient <- 0}
                             if (LabOrient == "Rot-20"){LabOrient <- 20}
                             if (LabOrient == "Rot-45"){LabOrient <- 45}
                             if (LabOrient == "Rot-70"){LabOrient <- 70}
                             if (LabOrient == "Vertical"){LabOrient <- 90}
                             Plot_Args$scales$rot <<- LabOrient
                             CtrlPlot() }, container=T4F_AxLabOrient)


   layoutAxis[4,1] <- T4F_XAxNameChange <- gframe("CHANGE X-LABEL", spacing=5, container=layoutAxis)
   T4_XAxNameChange <- gedit("", handler=function(h,...){
                             if(svalue(T4_XAxNameChange)==""){return()}
                             Plot_Args$xlab$label <<- svalue(T4_XAxNameChange)
                             CtrlPlot() } , container=T4F_XAxNameChange)

   layoutAxis[4,2] <- T4F_YAxNameChange <- gframe("CHANGE Y-LABEL", spacing=5, container=layoutAxis)
   T4_YAxNameChange <- gedit("",handler=function(h,...){
                             if(svalue(T4_YAxNameChange)==""){return()}
                             Plot_Args$ylab$label <<- svalue(T4_YAxNameChange) # in 2D Y is the vertical axis
                             CtrlPlot() }, container=T4F_YAxNameChange)

   layoutAxis[6,1] <- T4F_XStep <- gframe("X STEP", spacing=5, container=layoutAxis)
   T4_XStep <- gcheckbox("Custom X ticks", checked=FALSE, handler=function(h,...){
                             if(svalue(T4_XStep)==FALSE) {return()}
                             Core.Line <- svalue(T1CLineListCK)
                             NCL <- length(Core.Line)
                             RngXmin <- round(unlist(Xlim[[2]]), digits=0) #I need identify XRange components using the CoreLine names
                             RngXmax <- round(unlist(Xlim[[3]]), digits=0)
                             names(RngXmin) <- Xlim[[1]]
                             names(RngXmax) <- Xlim[[1]]
                             winTick <- gwindow(" X AXIS TICK INCREMENT ", parent=c(50,0), visible=FALSE)
                             DFgroup <- ggroup(horizontal=FALSE, container=winTick)
                             glabel("Please, give the increment between ticks", container=DFgroup)
                             Tick.Increment <- rep("?",NCL) #this is needed to construct correctly the data.frame
                             Tick.Increment <- data.frame(Core.Line,Tick.Increment, stringsAsFactors=FALSE) #in the dataframe add a column with variable names
                             DFrame <- gdf(items=Tick.Increment, container=DFgroup)
                             addHandlerChanged(DFrame, handler=function(h,...){ #addHandlerChanged dowload the dataFrame with modified parameters in NewFirParam (global variable)
                                                    Tick.Increment <<- h$obj[]
                                               })
                             gbutton(text=" EXIT ", handler= function(h, ...){
                                                    dispose(winTick)
                                                    Tick.Increment <- as.numeric(unlist(Tick.Increment[2])) #first element of Tick.Increment is the CL-names
                                                    AT <- list()
                                                    for(ii in 1:NCL){
                                                       AT[[ii]] <- seq(from=RngXmin[[Core.Line[ii]]], to=RngXmax[[Core.Line[ii]]], by=Tick.Increment[ii])
                                                    }
                                                    Plot_Args$scales$x$at <<- AT
                                                    CtrlPlot()
                                               }, container=DFgroup)
                                               visible(winTick) <- TRUE
                             }, container=T4F_XStep)

   layoutAxis[6,2] <- T4F_YStep <- gframe("Y STEP", spacing=5, container=layoutAxis)
   T4_YStep <- gcheckbox("Custom Y ticks ", checked=FALSE, handler=function(h,...){
                             if(svalue(T4_YStep)==FALSE) {return()}
                             Core.Line <- svalue(T1CLineListCK)
                             NCL <- length(Core.Line)
                             RngYmin <- round(unlist(Ylim[[2]]), digits=0) #I need identify XRange components using the CoreLine names
                             RngYmax <- round(unlist(Ylim[[3]]), digits=0)
                             names(RngYmin) <- Xlim[[1]]
                             names(RngYmax) <- Xlim[[1]]
                             winTick <- gwindow(" Y AXIS TICK INCREMENT ", parent=c(50,0), visible=FALSE)
                             DFgroup <- ggroup(horizontal=FALSE, container=winTick)
                             glabel("Please, give the increment between ticks", container=DFgroup)
                             Tick.Increment <- rep("?",NCL) #this is needed to construct correctly the data.frame
                             Tick.Increment <- data.frame(Core.Line,Tick.Increment, stringsAsFactors=FALSE) #in the dataframe add a column with variable names
                             DFrame <- gdf(items=Tick.Increment, container=DFgroup)
                             addHandlerChanged(DFrame, handler=function(h,...){ #addHandlerChanged dowload the dataFrame with modified parameters in NewFirParam (global variable)
                                                    Tick.Increment <<- h$obj[]
                                               })
                             gbutton(text=" EXIT ", handler= function(h, ...){
                                                    dispose(winTick)
                                                    Tick.Increment <- as.numeric(unlist(Tick.Increment[2])) #first element of Tick.Increment is the CL-names
                                                    AT <- list()
                                                    for(ii in 1:NCL){  #reshape the Y limits using rounded values
                                                       Ndigit <- nchar(RngYmin[[Core.Line[ii]]])
                                                       Y1 <- as.integer(RngYmin[[Core.Line[ii]]]/(10^(Ndigit-1)))* 10^(Ndigit-1)
                                                       Ndigit <- nchar(RngYmax[[Core.Line[ii]]])
                                                       Y2 <- as.integer(RngYmax[[Core.Line[ii]]]/(10^(Ndigit-2)))* 10^(Ndigit-2)+10^(Ndigit-2)
                                                       AT[[ii]] <- seq(from=Y1, to=Y2, by=Tick.Increment[ii])
                                                    }
                                                    Plot_Args$scales$y$at <<- AT
                                                    CtrlPlot()
                                               }, container=DFgroup)
                                               visible(winTick) <- TRUE
                             }, container=T4F_YStep)

   gbutton(" RESET PLOT ", handler=function(h,...){
                             ResetPlot()
                             CtrlPlot()
                             }, container=T4group1)

   gbutton(" EXIT ", handler=function(h,...){
				                         dispose(win)
                             }, container = T4group1)



# --- TAB5 ---

### LEGEND SETTINGS

   T5group1 <- ggroup(label="LEGEND", horizontal=FALSE, container=nb)

   layoutLeg <- glayout(homogeneous=FALSE, spacing=3, container=T5group1)

   layoutLeg[1,1] <- T5F_legendCK <- gframe(text="Enable Legend", spacing=5, container=layoutLeg)
   legendCK <- gcheckbox("Enable Legend ON/OFF", checked=TRUE,handler=function(h,...){
                          AutoKey_Args$text <<- unlist(FNamesCoreLines$XPSSample)  #load the Legends in the slot of the AutoKey_Args = List of parameters defining legend properties
	           	           Plot_Args$auto.key <<- AutoKey_Args #Save the AutoKey_Args list of par in Plot_Args$auto.key
                          if (svalue(legendCK)==TRUE) {
                             if (svalue(T3_SetLines)=="ON") {    #selezionate LINEE
                                Plot_Args$par.settings$superpose.line$col <<- "black" #B/W plot
                                Plot_Args$par.settings$superpose.line$lty <<- LType
                                Plot_Args$scales$relation <<- "free"
 		           	              if (svalue(T3_BW_Col)=="RainBow") {                    #COLOR plot
                                   Plot_Args$par.settings$superpose.line$col <<- Colors
                                   Plot_Args$par.settings$superpose.line$lty <<- "solid"
                                }
                             }
                             if (svalue(T3_SetSymbols)=="ON") {   #selezionate SIMBOLI
                                Plot_Args$par.settings$superpose.symbol$col <<- "black"  #B/W plot
                                Plot_Args$par.settings$superpose.symbol$pch <<- STypeIndx
                                Plot_Args$par.settings$superpose.symbol$pch <<- 1
                                Plot_Args$scales$relation <<- "free"
 		           	              if (svalue(T3_BW_Col)=="RainBow") {                    #COLOR plot
                                   Plot_Args$par.settings$superpose.symbol$col <<- Colors
                                }
                             }
                          } else {
		           	           Plot_Args$auto.key <<- FALSE
	           	           }
                          CtrlPlot()
                       }, container=T5F_legendCK)

   layoutLeg[1,2] <- T5F_LegFNameCK <- gframe(text="Add XPSSamp Name", spacing=5, container=layoutLeg)
   LegFNameCK <- gcheckbox("XPSSamp.Name ON/OFF", checked=TRUE,handler=function(h,...){
                          if (is.logical(Plot_Args$auto.key)){
                             gmessage("PLEASE ENABLE LEGENDS", icon="warning")
                             svalue(LegFNameCK) <- FALSE
                          } else {
                             if (svalue(LegFNameCK)==TRUE) {
                                Legends <- FNamesCoreLines$CoreLines
                                for (ii in seq_along(FNamesCoreLines$XPSSample)){
                                    tmp <- unlist(strsplit(Legends[ii], "\\."))  #skip the number at beginning coreline name
                                    Legends[ii] <- paste(tmp[2], "_", FNamesCoreLines$XPSSample[ii], sep="")
                                }
                                Plot_Args$auto.key$text <<- as.vector(Legends)
                             } else {
                                Legends <- FNamesCoreLines$CoreLines
                                for(ii in seq_along(Legends)){
                                   tmp <- unlist(strsplit(Legends[ii], "\\."))   #skip the number at beginning coreline name
                                   Legends[ii] <- tmp[2]
                                }
                                Plot_Args$auto.key$text <<- as.vector(Legends)
                             }
                          }
                          CtrlPlot()
                       }, container=T5F_LegFNameCK)

   layoutLeg[2,1] <- T5F_TSizeCK <- gframe(text="Text Size", spacing=5, container=layoutLeg)
   TSizeCK <- gcombobox(LegTxtSize,selected=1, toolkit = guiToolkit(), handler=function(h,...){
		           	        Plot_Args$auto.key$cex <<- as.numeric(svalue(TSizeCK))
                          CtrlPlot()
                       }, container=T5F_TSizeCK)



   layoutLeg[2,2] <- T5F_LegColCK <- gframe(text="Group Legend and organize in columns", spacing=5, container=layoutLeg)
   LegColCK <- gedit(initial.msg ="Col. numb.", selected=1, editable=FALSE, handler=function(h,...){
                          columns <- svalue(LegColCK)
                          Plot_Args$auto.key$columns <<- as.numeric(svalue(LegColCK))
                          CtrlPlot()
                       }, container=T5F_LegColCK)

   layoutLeg[3,1] <- T5F_LineWdhCK <- gframe(text="Line/Symbol weight", spacing=5, container=layoutLeg)
   LineWdhCK <- gcombobox(LWidth,selected=1, toolkit = guiToolkit(), handler=function(h,...){
                          weight <- as.numeric(svalue(LineWdhCK))
                          if (svalue(T3_SetLines)=="ON") {   #Lines selected
                             Plot_Args$par.settings$superpose.line$lwd <<- weight
                          }
                          if (svalue(T3_SetSymbols)=="ON") {   #Symbol selected
                             Plot_Args$par.settings$superpose.symbol$cex <<- weight
                          }
                          CtrlPlot()
                       }, container=T5F_LineWdhCK)

   layoutLeg[3,2] <- T5F_TxtColCK <- gframe(text="Legend text Color", spacing=5, container=layoutLeg)
   TxtColCK <- gcombobox(c("B/W", "RainBow"),selected=1, toolkit = guiToolkit(), handler=function(h,...){
                          if  (svalue(TxtColCK)=="B/W"){
                              Plot_Args$auto.key$col <<- "black"
                          } else {
                              Plot_Args$auto.key$col <<- Colors
                          }
                          CtrlPlot()
                       }, container=T5F_TxtColCK)



   layoutLeg[4,1] <- T5F_ChangLeg <- gbutton(text="Change Legend", spacing=5, handler=function(h,...){
                                LegWin <- gwindow(title="XPS Sample Legends", visible=FALSE) #open a new window to contain the list of new labels
                                LegGroup <- ggroup(horizontal=FALSE, container=LegWin)
                                glabel("                           EDIT LEGENDS                           ", container=LegGroup) #This label long to get a reasonable window dimension
                                LL=length(FNamesCoreLines$XPSSample)
                                Legends <- data.frame(LEGENDS=rep("-", LL), stringsAsFactors=FALSE)
                                LegDFrame <- gdf(items=Legends, container=LegGroup) #here no handler: it does not work in linux
                                size(LegDFrame) <- c(150,150)                           #size needed to generate a window with a non null-size window
                                addHandlerChanged(LegDFrame, handler=function(h,...){ #addHandlerChanged to add the handler to gdf() working also in linux.
                                      Legends <<- h$obj[]
                                })
                                gbutton("     SAVE LEGENDS & EXIT      ", handler=function(h,...){
                                      Plot_Args$auto.key$text <<- as.vector(Legends)
                                      dispose(LegWin)
                                      unblockHandler(Annotate)
                                      CtrlPlot()
                                }, container = LegGroup)
                                visible(LegWin) <- TRUE
                                CtrlPlot()
                       }, container=layoutLeg)

   layoutLeg[4,2] <- Annotate <- gbutton(text=" Annotate ", handler=function(h,...){
                                xx <- Plot_Args$xlim   #in the case of zoom Xlim, Ylim are not null
                                yy <<- Plot_Args$ylim
                                if (is.null(xx)){    #no zoom is present
                                   xx <- Xlim  #get the X range from OverlayEngine
                                   yy <- Ylim  #get the Y range from OverlayEngine
                                }
                                XPSLattAnnotate(xx, yy)
                       }, container=layoutLeg)


   gbutton(" RESET PLOT ", handler=function(h,...){
                                ResetPlot()
                                CtrlPlot()
                       }, container=T5group1)


   gbutton(" EXIT ", handler=function(h,...){
				                    dispose(win)
                       }, container = T5group1)


#----- END NOTEBOOK -----

   enabled(objFunctFact) <- FALSE
   svalue(nb) <- 5
   svalue(nb) <- 4
   svalue(nb) <- 3
   svalue(nb) <- 2
   svalue(nb) <- 1
   tcl("update", "idletasks")
   visible(win) <- TRUE
}
