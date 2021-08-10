#function to perform plots in overlay mode

#function to perform plots in overlay mode

#'Performs plots overlapping two different object with different Y scales
#'
#'Provides a userfriendly interface to select Objects to overlay
#'and a selection of major plotting options for a personalized data representation
#'Data are represented using different Left and Right Y scales
#'No parameters passed to this function
#'
#'
#'@examples
#'
#'\dontrun{
#'	XPSTwoScalePlot()
#'}
#'
#'@export 
#'


XPSTwoScalePlot <- function(){

   CtrlPlot <- function(){
#---- calls the macro which makes the plot following options
#            FName1 <- get(SelectedNames$XPSSample[1], envir=.GlobalEnv)

#            FName2 <- get(SelectedNames$XPSSample[2], envir=.GlobalEnv)


            Xlim <- range(c(Xlim1,Xlim2))
            if (svalue(T1_ReverseX)) {   #reverse X scale==TRUE
               Xlim <- sort(Xlim, decreasing=TRUE)
               Plot_Args$xlim <<- Xlim
            } else {
               Xlim <- sort(Xlim, decreasing=FALSE)
               Plot_Args$xlim <<- Xlim
            }

#--- plot XPSSamp1 ----
            Plot_Args$ylim <<- sort(Ylim1, decreasing=FALSE)
            Plot_Args$col <<- svalue(T2_Col1)
            Plot_Args$lty <<- svalue(T2_LineType1)
            Plot_Args$pch <<- STypeIndx[svalue(T2_SymType1, index=TRUE)]
            Plot_Args$data	<<- df1

            par(mar=c(5,5,5,5)) #set plot margins first
            with(df1, plot(x, y, xlim=Plot_Args$xlim, ylim=Plot_Args$ylim, axes=FALSE, ann=FALSE,
                           pch=Plot_Args$pch, cex=Plot_Args$cex, col=Plot_Args$col, #using dataframe df1 plot following parameters
                           lty=Plot_Args$lty, lwd=Plot_Args$lwd, type=Plot_Args$type,
#                           main=Plot_Args$main$label, cex.main=Plot_Args$main$cex,
#                           lines=Plot_Args$lines, point=Plot_Args$point, background=Plot_Args$background, col=Plot_Args$col,
#                           xlab=Plot_Args$xlab$label, cex.lab=Plot_Args$xlab$cex,
#                           ylab=Plot_Args$ylab$label, cex.lab=Plot_Args$ylab$cex,
#                           cex.axis=Plot_Args$scales$cex,
#                           xscale.components=Plot_Args$xscale.components, xscale.components=Plot_Args$xscale.components,
                           par.settings=Plot_Args$par.settings, grid=FALSE)
                )
            box() #frame around the plot
            mtext(side=3, line=1.5, text=Plot_Args$main$label, cex=Plot_Args$main$cex, col="black") # Now write the title

#--- draw X axis
            axis(side=1, labels=FALSE, col="black") #now plot X axis with ticks
            tks <- axTicks(1) # Use axTicks() to get default tick posiitons and numbers
            mtext(side=1, line=1, text=tks, at=tks, cex=Plot_Args$scales$cex, col="black") #now plot X numbers at distance line=1 (margins distance is in lines)
            mtext(side=1, line=2.5, text=Plot_Args$xlab$label, cex=Plot_Args$xlab$cex, col="black") # Now draw x label at distance = line=3

#--- draw Y1 axis
            axis(side=2, labels=FALSE, col=Plot_Args$col) #now plot left axis with color of dataset1
            tks <- axTicks(2) #Use axTicks() to get default tick posiitons and numbers
            mtext(side=2, line=1, text=tks, at=tks, cex=Plot_Args$scales$cex,     #now plot Y1 numbers at distance line=1 (margins distance is in lines)
                  col=Plot_Args$col, col.ticks=Plot_Args$col, col.axis=Plot_Args$col)
            mtext(side=2, line=2.7, text=Plot_Args$ylab$label, cex=Plot_Args$ylab$cex, col=Plot_Args$col) # Now draw Y1 label at distance = line=3

#--- plot XPSSamp2 ----
            Plot_Args$ylim <<- sort(Ylim2, decreasing=FALSE)
            Plot_Args$col <<- svalue(T2_Col2)
            Plot_Args$lty <<- svalue(T2_LineType2)
            Plot_Args$pch <<- STypeIndx[svalue(T2_SymType2, index=TRUE)]
            Plot_Args$data <<- df2

            par(new=T) #superpose new plot
            with(df2, plot(x, y, xlim=Plot_Args$xlim, ylim=Plot_Args$ylim, axes=FALSE, ann=FALSE,
                           pch=Plot_Args$pch, cex=Plot_Args$cex, col=Plot_Args$col, #using dataframe df1 plot following parameters
                           lty=Plot_Args$lty, lwd=Plot_Args$lwd, type=Plot_Args$type,
#                           main=Plot_Args$main$label, cex.main=Plot_Args$main$cex,
#                           lines=Plot_Args$lines, point=Plot_Args$point, background=Plot_Args$background,
#                           ylab=Plot_Args$ylab.right$label, cex.lab=Plot_Args$ylab.right$cex,
#                           main=Plot_Args$main, xlab=Plot_Args$xlab, ylab=Plot_Args$ylab,
#                           scales=Plot_Args$scales,
#                           xscale.components=Plot_Args$xscale.components, xscale.components=Plot_Args$xscale.components,
                           par.settings=Plot_Args$par.settings, grid=FALSE)
                )
            axis(side=4, labels=FALSE, col=Plot_Args$col) #now plot right axis with color of dataset2
            tks <- axTicks(2) # Use axTicks() to get default tick posiitons and numbers
            mtext(side=4, line=1, text=tks, at=tks, cex=Plot_Args$scales$cex,  #now plot Y2 numbers at distance line=1 (margins distance is in lines)
                  col=Plot_Args$col, col.ticks=Plot_Args$col, col.axis=Plot_Args$col)
            mtext(side=4, line=2.7, text=Plot_Args$ylab.right$label, cex=Plot_Args$ylab.right$cex, col=Plot_Args$col) # Now draw Y2 label at distance = line=3

            if(Plot_Args$grid) { grid() } #add grid to plot

            LTyp <- Sym <- NULL
            LW <- 0
            LegPos <- svalue(T3_LegPos1)
            if(LegPos != "OFF"){ #enable legend upon selction
               if (Plot_Args$type=="l" || Plot_Args$type=="b") { 
                  LTyp <- c(svalue(T2_LineType1), svalue(T2_LineType2))
                  LW <- Plot_Args$lwd
               }
               if (Plot_Args$type=="p" || Plot_Args$type=="b") { 
                  Sym <- c(STypeIndx[svalue(T2_SymType1, index=TRUE)], STypeIndx[svalue(T2_SymType2, index=TRUE)]) 
               }
               legend (LegPos, ncol=1, bty="n",    #position center, no box
                      legend=c(svalue(T3_Legend2), svalue(T3_Legend3)),
                      lty=LTyp, lwd=LW, pt.lwd=Plot_Args$lwd, pch=Sym, pt.cex=Plot_Args$cex,
                      col=c(svalue(T2_Col1), svalue(T2_Col2)))
            }
   }

#----- reset parametri ai valori iniziali -----

   ResetPlot <- function(){
            svalue(T1FNameListCK1) <- FALSE
            svalue(T1_CoreLineCK1) <- FALSE
            svalue(T1FNameListCK2) <- FALSE
            svalue(T1_CoreLineCK2) <- FALSE
            svalue(T1_ReverseX) <- FALSE
            svalue(T2_Col1, index=TRUE) <- 1
            svalue(T2_Col2, index=TRUE) <- 1
            svalue(T2_Grid, index=TRUE) <-  1
            svalue(T2_SetLines, index=TRUE) <- 1
            svalue(T2_SetSymbols, index=TRUE) <- 2
            svalue(T2_LineType1, index=TRUE) <- 1
            svalue(T2_LineType2, index=TRUE) <- 2
            svalue(T2_LinWidth, index=TRUE) <- 1
            svalue(T2_SymType1, index=TRUE) <- 1
            svalue(T2_SymType2, index=TRUE) <- 2
            svalue(T2_SymSize, index=TRUE) <- 4
            svalue(T3_TitSize, index=TRUE) <- 5
            svalue(T3_MainTitChange) <- ""
            svalue(T3_AxNumSize, index=TRUE) <- 3
            svalue(T3_AxLabSize, index=TRUE) <- 4
            svalue(T3_XAxNameChange) <- ""
            svalue(T3_LYAxNameChange) <- ""
            svalue(T3_RYAxNameChange) <- ""
            svalue(T3_Legend1) <- ""
            svalue(T3_Legend2) <- ""
            Xlim1<-NULL
            Xlim2<-NULL
            Ylim1<-NULL
            Ylim2<-NULL

            Plot_Args<<-list( x=formula("y ~ x"), data=NULL, groups=NULL,layout=NULL,
                    xlim=NULL, ylim=NULL,
                    pch=STypeIndx,cex=1,lty=LType,lwd=1,type="l",
                    lines=TRUE, point=FALSE,
                    background="transparent", col="black",
                    main=list(label=NULL,cex=1.4),
                    xlab=list(label="X", rot=0, cex=1),
                    ylab=list(label="Y1", rot=90, cex=1),
                    ylab.right=list(label="Y2", rot=90, cex=1),
                    par.settings = list(superpose.symbol=list(pch=STypeIndx,fill="black"), #setta colore riempimento simboli
                                        superpose.line=list(lty=LType, col="black")), #necessario per settare colori legende
                    grid = FALSE
                  )
            LegTxt<<-NULL
            CtrlPlot()
   }

#----- Variabili -----
   plot.new()                               #reset fineztra grafica
   if (is.na(activeFName)){
       gmessage("No data present: please load and XPS Sample", title="XPS SAMPLES MISSING", icon="error")
       return()
   }
   FNameListTot <- as.array(XPSFNameList())     #lista di tutti gli XPSSample caricati nel WorkSpace
   df1 <- data.frame(x=NULL, y=NULL)
   df2 <- data.frame(x=NULL, y=NULL)
   NCorelines <- NULL
   Xlim1 <- NULL
   Xlim2 <- NULL
   Ylim1 <- NULL
   Ylim2 <- NULL

   SelectedNames <- list(XPSSample=NULL, CoreLines=NULL)
   NamesList=list(XPSSample=NULL, CoreLines=NULL)
   FontSize <- c(0.6,0.8,1,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6,2.8,3)
   LWidth <- c(1,1.25,1.5,1.75,2,2.25,2.5,3, 3.5,4)
   SymSize <- c(0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2)
   LegTxt <- NULL
   LegPos <- c("OFF", "topleft", "top", "topright", "left",  "center", "right", "bottomleft", "bottom", "bottomright")

   XPSSettings <- get("XPSSettings", envir=.GlobalEnv)
   Colors <- XPSSettings$Colors
   LType <- XPSSettings$LType
   SType <- XPSSettings$Symbols
   STypeIndx <- XPSSettings$SymIndx

   Plot_Args<-list( x=formula("y ~ x"), data=NULL, groups=NULL,layout=NULL,
                    xlim=NULL, ylim=NULL,
                    pch=STypeIndx,cex=1,lty=LType,lwd=1,type="l",
                    lines=TRUE, point=FALSE,
                    background="transparent", col="black",
                    main=list(label=NULL,cex=1.2),
                    xlab=list(label="X", rot=0, cex=1),
                    ylab=list(label="Y1", rot=90, cex=1),
                    ylab.right=list(label="Y2", rot=90, cex=1),
                    par.settings = list(superpose.symbol=list(pch=STypeIndx,fill="black"), #set symbol fill color
                                        superpose.line=list(lty=LType, col="black")), #needed to set colors and legends
                    grid = FALSE
                  )

#----- GUI -----

   win <- gwindow(" TWO Y-SCALE PLOT ", visible=FALSE)
   maingroup <- ggroup(horizontal=FALSE, container=win)

   nb <- gnotebook(expand=TRUE, container = maingroup)

# --- TAB1 ---
#XPS Sample selecion & representation options

     T1group1 <- ggroup(label="XPS SAMPLE SELECTION", spacing=5, horizontal=FALSE, container=nb)

     layoutT1 <- glayout(homogeneous=FALSE, spacing=5, container=T1group1)

     layoutT1[1,1] <- T1F_FName1 <- gframe(text="SELECT XPS-SAMPLE 1", spacing=5, container=layoutT1)
     T1FNameListCK1 <- gcombobox(FNameListTot,selected=-1, handler=function(h,...){
                                    SelectedNames$XPSSample[1] <<- svalue(T1FNameListCK1)
                                    SpectList1 <- XPSSpectList(SelectedNames$XPSSample[1])   #sequenze delle CoreLines XPSSample attivo
                                    delete(T1F_CoreLines1,T1_CoreLineCK1) #cancello la lista delle coreline quando ho fatto la selezione
                                    T1_CoreLineCK1 <<- gcombobox(SpectList1,selected=-1, editable=FALSE, handler=function(h,...){
                                                       SelectedNames$CoreLines[1] <<- svalue(T1_CoreLineCK1)
                                                       FName1 <- get(SelectedNames$XPSSample[1], envir=.GlobalEnv)
                                                       idx <- unlist(strsplit(SelectedNames$CoreLines[1], "\\."))
                                                       idx <- as.numeric(idx[1])
                                                       Xlim1 <<- range(FName1[[idx]]@.Data[[1]])
                                                       Ylim1 <<- range(FName1[[idx]]@.Data[[2]])
                                                       df1 <<- data.frame(x=unlist(FName1[[idx]]@.Data[1]), y=unlist(FName1[[idx]]@.Data[2]))
                                    }, container=T1F_CoreLines1)
                                    enabled(T1_CoreLineCK1) <- TRUE
                               }, container=T1F_FName1)

     layoutT1[1,2] <- T1F_CoreLines1 <- gframe(text="SELECT CORE LINE",  spacing=5, container=layoutT1)
     T1_CoreLineCK1 <- gcombobox("       ",selected=-1,  handler=NULL, container=T1F_CoreLines1)


     layoutT1[2,1] <- T1F_FName2 <- gframe(text="SELECT XPS-SAMPLE 2", spacing=5, container=layoutT1)
     T1FNameListCK2 <- gcombobox(FNameListTot,selected=-1, handler=function(h,...){
                                    SelectedNames$XPSSample[2] <<- svalue(T1FNameListCK2)
                                    SpectList2 <- XPSSpectList(SelectedNames$XPSSample[2])   #sequenze delle CoreLines XPSSample attivo
                                    delete(T1F_CoreLines2 ,T1_CoreLineCK2) #cancello la lista delle coreline quando ho fatto la selezione
                                    T1_CoreLineCK2 <<- gcombobox(SpectList2,selected=-1, handler=function(h,...){
                                                       SelectedNames$CoreLines[2] <<- svalue(T1_CoreLineCK2)
                                                       FName2 <- get(SelectedNames$XPSSample[2], envir=.GlobalEnv)
                                                       idx <- unlist(strsplit(SelectedNames$CoreLines[2], "\\."))
                                                       idx <- as.numeric(idx[1])
                                                       Xlim2 <<- range(FName2[[idx]]@.Data[[1]])
                                                       Ylim2 <<- range(FName2[[idx]]@.Data[[2]])
                                                       df2 <<- data.frame(x=unlist(FName2[[idx]]@.Data[1]), y=unlist(FName2[[idx]]@.Data[2]))
                                   }, container=T1F_CoreLines2)
                                    enabled(T1_CoreLineCK2) <- TRUE
                               }, container=T1F_FName2)

     layoutT1[2,2] <- T1F_CoreLines2 <- gframe(text="SELECT CORE LINE",  spacing=5, container=layoutT1)
     T1_CoreLineCK2 <- gcombobox("       ",selected=-1, handler=NULL, container=T1F_CoreLines2)

     layoutT1[3,1] <- T1F_ReverseX <- gframe("REVERSE X-axis", spacing=5, container=layoutT1)
     T1_ReverseX <- gcheckbox("Reverse X axis",checked=FALSE, handler=function(h,...){
                                    CtrlPlot()
                               }, container=T1F_ReverseX)


     gbutton("PLOT", handler=function(h,...){
                                     if (SelectedNames$XPSSample[1]==SelectedNames$XPSSample[2]){
                                        answ <- gconfirm("XPSSample1 = XPSSample2! Proceed anyway?" , title = "WARNING: SAME XPSSAMPLE DATA SELECTED",  icon = "warning")
                                        if (! answ) { 
                                           return() 
                                        } else {
                                           CtrlPlot()
                                        }
                                     } else {
                                        CtrlPlot()
                                     }
                               }, container=maingroup)

     gbutton("EXIT", handler=function(h,...){
				                         dispose(win)
                               }, container=maingroup)

# --- TAB2 ---
# Rendering options

     T2group1 <- ggroup(label="RENDERING", horizontal=FALSE, container=nb)
     layoutRend <- glayout(homogeneous=FALSE, spacing=3, container=T2group1)

     layoutRend[1,1]  <- T2F_Col1 <- gframe("COLOR XPSSample1", spacing=5, container=layoutRend)
     T2_Col1 <- gcombobox(Colors, selected=1, editable=FALSE, handler=function(h,...){
                             Plot_Args$col <<- c(svalue(T2_Col1),svalue(T2_Col2))
                             Plot_Args$par.settings$superpose.line$col <<- c(svalue(T2_Col1),svalue(T2_Col2))
                             Plot_Args$par.settings$superpose.symbol$col <<- c(svalue(T2_Col1),svalue(T2_Col2))
                            CtrlPlot() }, container=T2F_Col1)

     layoutRend[1,2]  <- T2F_Col2 <- gframe("COLOR XPSSample2", spacing=5, container=layoutRend)
     T2_Col2 <- gcombobox(Colors, selected=1, editable=FALSE, handler=function(h,...){
                             Plot_Args$col <<- c(svalue(T2_Col1),svalue(T2_Col2))
                             Plot_Args$par.settings$superpose.line$col <<- c(svalue(T2_Col1),svalue(T2_Col2))
                             Plot_Args$par.settings$superpose.symbol$col <<- c(svalue(T2_Col1),svalue(T2_Col2))
                             CtrlPlot() }, container=T2F_Col2)

     layoutRend[1,3]  <- T2F_Grid <- gframe("GRID", spacing=5, container=layoutRend)
     T2_Grid <- gcombobox(c("Grid ON", "Grid OFF"), selected=-1, editable=FALSE, handler=function(h,...){
                             if(svalue(T2_Grid)=="Grid ON") {
                                Plot_Args$grid <<- TRUE
                             } else {
                                Plot_Args$grid <<- FALSE
                             }
                             CtrlPlot() }, container=T2F_Grid)

     layoutRend[2,1] <- T2F_SetLines <- gframe("SET LINES", spacing=5, container=layoutRend)
     T2_SetLines <-  gradio(c("ON", "OFF"), selected=1, horizontal = TRUE, handler=function(h,...){
                             if (svalue(T2_SetLines)=="ON" && svalue(T2_SetSymbols)=="OFF"){
                                 Plot_Args$type <<- "l"
                             } else if (svalue(T2_SetLines)=="OFF" && svalue(T2_SetSymbols)=="OFF"){
                                 Plot_Args$type <<- NA
                             } else if (svalue(T2_SetLines)=="OFF" && svalue(T2_SetSymbols)=="ON"){
                                 Plot_Args$type <<- "p"
                             } else if (svalue(T2_SetLines)=="ON" && svalue(T2_SetSymbols)=="ON"){
                                 Plot_Args$type <<- "b"
                             }
                             CtrlPlot() }, container=T2F_SetLines)

     layoutRend[2,2] <- T2F_SetSymbols <- gframe("SET SYMBOLS", horizontal=TRUE, spacing=5, container=layoutRend)
     T2_SetSymbols <-  gradio(c("ON", "OFF"), selected=2, horizontal=TRUE, handler=function(h,...){
                             if (svalue(T2_SetSymbols)=="ON" && svalue(T2_SetLines)=="OFF"){
                                 Plot_Args$type <<- "p"
                             } else if (svalue(T2_SetSymbols)=="OFF" && svalue(T2_SetLines)=="OFF"){
                                 Plot_Args$type <<- NA
                             } else if (svalue(T2_SetSymbols)=="OFF" && svalue(T2_SetLines)=="ON"){
                                 Plot_Args$type <<- "l"
                             } else if (svalue(T2_SetSymbols)=="ON" && svalue(T2_SetLines)=="ON"){
                                 Plot_Args$type <<- "b"
                             }
                             CtrlPlot() }, container=T2F_SetSymbols)

     layoutRend[3,1] <- T2F_SetLines1 <- gframe("LINE TYPE XPSSamp1", spacing=5, container=layoutRend)
     T2_LineType1 <- gcombobox(LType, selected=1, editable=FALSE, handler=function(h,...){
                             Plot_Args$par.settings$superpose.line$lty <<- c(svalue(T2_LineType1),svalue(T2_LineType2))
                             CtrlPlot()
                           }, container=T2F_SetLines1)

     layoutRend[3,2] <- T2F_SetLines2 <- gframe("LINE TYPE XPSSamp2", spacing=5, container=layoutRend)
     T2_LineType2 <- gcombobox(LType, selected=2, editable=FALSE, handler=function(h,...){
                             Plot_Args$par.settings$superpose.line$lty <<- c(svalue(T2_LineType1),svalue(T2_LineType2))
                             CtrlPlot()
                           }, container=T2F_SetLines2)

     layoutRend[3,3] <- T2F_LinWidth <- gframe("LINE WIDTH", spacing=5, container=layoutRend)
     T2_LinWidth <- gcombobox(LWidth, selected=1, editable=FALSE, handler= function(h,...){
                             Plot_Args$lwd <<- as.numeric(svalue(T2_LinWidth))
                             CtrlPlot()
                           }, container=T2F_LinWidth)


     layoutRend[4,1] <- T2F_SymType1 <- gframe("SYMBOLS XPSSamp1", spacing=5, container=layoutRend)
     T2_SymType1 <- gcombobox(SType, selected=1, editable=FALSE, handler=function(h,...){
                             Plot_Args$par.settings$superpose.symbol$pch <<- c(svalue(T2_SymType1),svalue(T2_SymType2))
                             CtrlPlot()
                             CtrlPlot()   #this is needed to refresh legend symbols 
                           }, container=T2F_SymType1)

     layoutRend[4,2] <- T2F_SymType2 <- gframe("SYMBOLS XPSSamp2", spacing=5, container=layoutRend)
     T2_SymType2 <- gcombobox(SType, selected=2, editable=FALSE, handler=function(h,...){
                             Plot_Args$par.settings$superpose.symbol$pch <<- c(svalue(T2_SymType1),svalue(T2_SymType2))
                             CtrlPlot()
                             CtrlPlot()   #this is needed to refresh legend symbols
                           }, container=T2F_SymType2)

     layoutRend[4,3] <- T2F_SymSize <- gframe("SYMSIZE", spacing=5, container=layoutRend)
     T2_SymSize <- gcombobox(SymSize, selected=4, editable=FALSE, handler= function(h,...){
                             Plot_Args$cex <<- as.numeric(svalue(T2_SymSize))
                             CtrlPlot()
                           }, container=T2F_SymSize)

     gbutton("RESET", handler=function(h,...){
                             ResetPlot()
                             CtrlPlot()
                            }, container=T2group1)


# --- TAB3 ---
# Axis Rendering options

   T3group1 <- ggroup(label="AXES", horizontal=FALSE, container=nb)
   layoutAxis <- glayout(homogeneous=FALSE, spacing=3, container=T3group1)


   layoutAxis[1,1] <- T3F_TitSize <- gframe("TITLE SIZE", spacing=5, container=layoutAxis)
   T3_TitSize <- gcombobox(FontSize, selected=5, editable=FALSE, handler= function(h,...){
                             Plot_Args$main$cex <<- as.numeric(svalue(T3_TitSize))
                             CtrlPlot() }, container=T3F_TitSize)

   layoutAxis[1,2] <- T3F_MainTitChange <- gframe("CHANGE MAIN TITLE", spacing=5, container=layoutAxis)
   T3_MainTitChange <- gedit("", handler=function(h,...){
                             Plot_Args$main$label <<- svalue(T3_MainTitChange)
                             CtrlPlot() }, container=T3F_MainTitChange)

   layoutAxis[2,1] <- T3F_AxNumSize <- gframe("AXIS NUMBER SIZE", spacing=3, container=layoutAxis)
   T3_AxNumSize <- gcombobox(FontSize, selected=3, editable=FALSE, handler= function(h,...){
                             Plot_Args$scales$cex <<- as.numeric(svalue(T3_AxNumSize))  #controls the size X axis only
                             CtrlPlot() }, container=T3F_AxNumSize)

   layoutAxis[2,2] <- T3F_AxLabSize <- gframe("AXIS LABEL SIZE", spacing=4, container=layoutAxis)
   T3_AxLabSize <- gcombobox(FontSize, selected=4, editable=FALSE, handler= function(h,...){
                             Plot_Args$xlab$cex <<- as.numeric(svalue(T3_AxLabSize))
                             Plot_Args$ylab$cex <<- as.numeric(svalue(T3_AxLabSize))
                             Plot_Args$ylab.right$cex <<- as.numeric(svalue(T3_AxLabSize))
                             CtrlPlot() }, container=T3F_AxLabSize)

   layoutAxis[3,1] <- T3F_XAxNameChange <- gframe("CHANGE X-LABEL", spacing=5, container=layoutAxis)
   T3_XAxNameChange <- gedit("", handler=function(h,...){
                             Plot_Args$xlab$label <<- svalue(T3_XAxNameChange)
                             CtrlPlot() } , container=T3F_XAxNameChange)

   layoutAxis[3,2] <- T3F_LYAxNameChange <- gframe("CHANGE LEFT Y-LABEL", spacing=5, container=layoutAxis)
   T3_LYAxNameChange <- gedit("",handler=function(h,...){
                             Plot_Args$ylab$label <<- svalue(T3_LYAxNameChange)
                             CtrlPlot() }, container=T3F_LYAxNameChange)

   layoutAxis[3,3] <- T3F_RYAxNameChange <- gframe("CHANGE RIGHT Y-LABEL", spacing=5, container=layoutAxis)
   T3_RYAxNameChange <- gedit("",handler=function(h,...){
                             Plot_Args$ylab.right$label <<- svalue(T3_RYAxNameChange)
                             CtrlPlot() }, container=T3F_RYAxNameChange)
                             
   layoutAxis[4,1] <- T3F_Legend1 <- gframe("PLOT LEGEND", spacing=5, container=layoutAxis)
   T3_LegPos1 <- gcombobox(LegPos, slected=1, handler=function(h,...){
                                    CtrlPlot()
                               }, container=T3F_Legend1)

   layoutAxis[4,2] <- T3F_Legend2 <- gframe("XPSSample1 LEGEND", spacing=5, container=layoutAxis)
   T3_Legend2 <- gedit("",handler=function(h,...){
                             LegTxt <<- c(svalue(T3_Legend2), svalue(T3_Legend3))
                             CtrlPlot() }, container=T3F_Legend2)

   layoutAxis[4,3] <- T3F_Legend3 <- gframe("XPSSample2 LEGEND", spacing=5, container=layoutAxis)
   T3_Legend3 <- gedit("",handler=function(h,...){
                             LegTxt <<- c(svalue(T3_Legend2), svalue(T3_Legend3))
                             CtrlPlot() }, container=T3F_Legend3)

#--- to modify the default X, Y1, Y2 scales
   layoutAxis[5,1] <- T3F_X <- gframe(" X range ", horizontal=FALSE, spacing=5, container=layoutAxis)
   T3_Xmin <- gedit(initial.msg="X min value:",handler=function(h,...){
                             Xmin <- as.numeric(svalue(T3_Xmin))
                             Xlim1 <<- c(Xmin, max(Xlim1))
                             CtrlPlot() }, container=T3F_X)
   T3_Xmax <- gedit(initial.msg="X max value:",handler=function(h,...){
                             Xmax <- as.numeric(svalue(T3_Xmax))
                             Xlim1 <<- c(min(Xlim1), Xmax)
                             CtrlPlot() }, container= T3F_X)

   layoutAxis[5,2] <- T3F_Y1 <- gframe(" Y1 range ", horizontal=FALSE, spacing=5, container=layoutAxis)
   T3_Ymin1 <- gedit(initial.msg="Y1 min value:",handler=function(h,...){
                             Ymin1 <- as.numeric(svalue(T3_Ymin1))
                             Ylim1 <<- c(Ymin1, max(Ylim1))
                             CtrlPlot() }, container= T3F_Y1)
   T3_Ymax1 <- gedit(initial.msg="Y1 max value:",handler=function(h,...){
                             Ymax1 <- as.numeric(svalue(T3_Ymax1))
                             Ylim1 <<- c(min(Ylim1), Ymax1)
                             CtrlPlot() }, container= T3F_Y1)

   layoutAxis[5,3] <- T3F_Y2 <- gframe(" Y2 range ", horizontal=FALSE, spacing=5, container=layoutAxis)
   T3_Ymin2 <- gedit(initial.msg="Y2 min value:",handler=function(h,...){
                             Ymin2 <- as.numeric(svalue(T3_Ymin2))
                             Ylim2 <<- c(Ymin2, max(Ylim1))
                             CtrlPlot() }, container= T3F_Y2)
   T3_Ymax2 <- gedit(initial.msg="Y2 max value:",handler=function(h,...){
                             Ymax2 <- as.numeric(svalue(T3_Ymax2))
                             Ylim2 <<- c(min(Ylim2), Ymax2)
                             CtrlPlot() }, container= T3F_Y2)

   gbutton("RESET", handler=function(h,...){
                             ResetPlot()
                             CtrlPlot()
                            }, container=T3group1)


   enabled(T1_CoreLineCK1) <- FALSE
   enabled(T1_CoreLineCK2) <- FALSE
   svalue(nb) <- 3
   svalue(nb) <- 2
   svalue(nb) <- 1
   visible(win) <- TRUE
   
}
