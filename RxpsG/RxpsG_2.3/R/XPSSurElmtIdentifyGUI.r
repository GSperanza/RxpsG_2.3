#'Identify Elements on the Survey
#'
#'GUI to Identify Elements on the Survey.
#'
#'No parameters are passed to this function
#'
#'@examples
#'
#'\dontrun{
#'	XPSSurveyElementIdentify()
#'}
#' 
#'@export
#'


XPSSurveyElementIdentify <- function(){

replot <- function(show.points = FALSE, ...) {
	  plot(newcoreline, xlim = point.coords$x, ylim = point.coords$y)
	     if ( show.points ) {
		       points(point.coords, col=2, cex=2, lwd=1.5, pch=10)
      }
   }

reset.boundaries <- function(h, ...) {
	  newcoreline <<- XPSremove(newcoreline, "all") # reset all
   point.coords <<- IniBoundaries
	  slot(newcoreline,"Boundaries") <<- point.coords
	  newcoreline <<- XPSsetRegionToFit(newcoreline)
}





#---- Var-Initialization

   Object <- get(activeFName, envir = .GlobalEnv)
   CurrentCoreLine <- activeSpectIndx
   nameCtrl <- Object[[CurrentCoreLine]]@Symbol
cat("\n xxxxx", nameCtrl)
   if (nameCtrl != "Survey" && nameCtrl != "survey") {
      gmessage("Please choose a survey spectrum" , title = "BAD CORELINE SELECTION",  icon = "warning") #controlla che lo spettro sia un survey
      return()
   }
   if ( is.numeric(CurrentCoreLine) ) {
   	newcoreline <- Object[[CurrentCoreLine]] # duplicate CoreLine
   } else if ( is.character(CurrentCoreLine) ) {
   	indx <- which(CurrentCoreLine == names(Object), arr.ind=TRUE)
   	newcoreline <- Object[[indx]] # duplicate CoreLine
   }
   if (newcoreline@Flags[3]) {
       ftype<-"scienta" #scienta filetype
   } else {
       ftype<-"kratos"  #kratos filetype
   }
   point.coords <<- list(x=range(newcoreline[1]), y=range(newcoreline[2]))
   IniBoundaries <<- point.coords       #save the original boundaries
   coord <- NULL
   Peaks <- NULL
   IdPeaks <- NULL
   rangeY <- range(Object[[CurrentCoreLine]]@.Data[2])
   ElmtList<-ReadElmtList("CoreLines") #reads the CoreLine/Auger List Table and save it in MyEnv  (macro in XPSSurveyUtility)
   ElmtList <- format(ElmtList, justify="centre", width=10)
   ElmtList <- as.data.frame(ElmtList,  stringsAsFactors = FALSE)
   RecPlot <- recordPlot()   #save graph for Refresh


#----- GUI -----

	  mainWin <- gwindow("ELEMENT IDENTIFICATION", visible = FALSE)
	  mainGroup <- ggroup(horizontal=FALSE, container = mainWin)

#--- Notebook

	  nb <- gnotebook(container = mainGroup, expand = FALSE)

#--- TAB1
	gtab1 <- ggroup(label = "Survey", horizontal=FALSE, spacing = 5, container = nb)
   gframe11 <- gframe(text = " Peak Detection ", horizontal = TRUE, spacing = 5, container = gtab1)
   glyt1<-glayout(homogeneous=FALSE, spacing=1, container=gframe11)

   glyt1[1,1]<-glabel(text = " Noise Level:", container=glyt1)
   glyt1[1,2]<-glabel(text = " Energy Window:", container=glyt1)
   glyt1[2,1]<-noiseLevel<-gedit(text = "7", width = 5, coerce.with=as.numeric, container=glyt1)
   glyt1[2,2]<-Erange<-gedit(text = "5", width = 5, coerce.with=as.numeric, container=glyt1)

   glyt1[3,1]<-detectionButton <- gbutton("  Detection  ", handler = function(h, ...) {
                      snmin<-svalue(noiseLevel)
                      Ewin<-svalue(Erange)
	                   Peaks <<- peakDetection(newcoreline, snmin, Ewin)
	                   plotPeaks(newcoreline, Peaks)
	                   svalue(sb) <- sprintf("Element detection done.")
	                   enabled(identificationButton) <- TRUE
	                   enabled(grfame13) <- TRUE
                      RecPlot <<- recordPlot()   #save graph for Refresh
               }, container = glyt1)

   gframe12 <- gframe(text = " Peak Identification ", horizontal = TRUE, spacing = 5, container = gtab1)
   glyt2<-glayout(homogeneous=FALSE, spacing=1, container=gframe12)

   glyt2[1,1]<-glabel(text = "  Precision (eV)  ", container = glyt2)
   glyt2[2,1]<-DEnergy<-gedit(text = "1", width = 5, container=glyt2)
	  glyt2[2,2]<-identificationButton <- gbutton("   Identification   ", handler = function(h, ...) {
                    DEnergy <- as.numeric(svalue(DEnergy))
                    IdPeaks <<- peakIdentification(Peaks, DEnergy, newcoreline, svalue(SpectType1), ElmtList)
                    ShowTablePeaks(IdPeaks, svalue(SpectType1))   #show table on the consolle
	              }, container = glyt2)
   enabled(identificationButton) <- FALSE

	  grfame13 <- gframe(text = " Plot ", container = gtab1, horizontal = TRUE, spacing = 5)
	  plotType <- gradio(items=c("normal", "corrected"), selected=1, handler = function(h, ...) {
		                  if (! is.null(Peaks)) {
		   	                 plotPeaks(newcoreline, Peaks, type = svalue(plotType))
		                  }
               }, container = grfame13)
	  enabled(grfame13) <- FALSE

  	gframe14 <- gframe(text = " Spectral Type selection ", container = gtab1, horizontal = FALSE)
   SpectType1 <- gradio(c("CoreLines", "AugerTransitions"), selected=1, horizontal = TRUE, spacing=5,  handler = function(h,...){
                     SType <- svalue(SpectType1)
                     ElmtList<<-ReadElmtList(SType) #reads the CoreLine/Auger List Table and save it in MyEnv  (macro in XPSSurveyUtility)
                     ElmtList<<- format(ElmtList, justify="centre", width=10)
                     ElmtList<<- as.data.frame(ElmtList,  stringsAsFactors = FALSE)
                     svalue(SpectType2)<-SType
                     svalue(SpectType3)<-SType
               }, container = gframe14)


#--- TAB2

	  gtab2 <- ggroup( label = " Single peak ", horizontal=FALSE, expand = FALSE, spacing = 5, container = nb)
	  gframe21 <- gframe(text = " Peak BE value  ", container = gtab2, horizontal = FALSE, expand = FALSE, spacing = 5)
	  glyt3 <- gformlayout(container = gframe21)
   BEbutton <- gbutton(label=" Peak Position: ", text="Press to locate", handler=function(h,...){
                    coord<<-locator(n=1, type="p", col="red", lwd=2)
                    xx<-as.character(round(as.numeric(coord$x),digits=2))
                    svalue(be_value)<-xx
               }, container = glyt3)

   be_value<- gedit(label="Energy value: ", initial.msg="BE?", handler = function(h,...) {
                    BE<-svalue(coord$x)
                    if (!is.null(BE)) { svalue(be_value)<-BE } #svuoto il valore del locator nel widget cursor position in modo che su qs valore siano attive tutte le funzioni
                    replayPlot(RecPlot)
#	 	              replot()
               }, container = glyt3)

	  be_win <- gedit("2", label="Precision (eV):", coerce.with = as.numeric, handler = function(h,...) {
                    replayPlot(RecPlot)
               }, container= glyt3)

	  gframe22 <- gframe(text = " Spectrum Type Selection ", container = gtab2, horizontal = FALSE)
   SpectType2 <- gradio(c("CoreLines", "AugerTransitions"), selected=1, horizontal = TRUE, spacing=5,  handler = function(h,...){
                     SType <- svalue(SpectType2)
                     ElmtList<<-ReadElmtList(SType) #reads the CoreLine/Auger List Table and save it in MyEnv  (macro in XPSSurveyUtility)
                     ElmtList<<- format(ElmtList, justify="centre", width=10)
                     ElmtList<<- as.data.frame(ElmtList,  stringsAsFactors = FALSE)
                     svalue(SpectType1)<-SType
                     svalue(SpectType3)<-SType
               }, container = gframe22)

	  gframe23 <- gframe(text = " Core Line Search Mode ", container = gtab2, horizontal = FALSE, expand = FALSE, spacing = 5)
   glyt4<-glayout(homogeneous=FALSE, spacing=1, container=gframe23)

	  glyt4[1,1]<-gbutton("  Find Elements Nearer to the Selected Energy  ", handler = function(h, ...) {
                     replot()
		                   if ( !is.na(svalue(be_value)) && !is.na(svalue(be_win)) ) {
			                     IdPeaks<<-NearerElement(svalue(be_value), svalue(be_win), newcoreline, ElmtList)
			                     if (length(IdPeaks)>0){
   		                      ShowTablePeaks(IdPeaks, svalue(SpectType2))
                        } else {
                           cat("\n No element found with the selected precision!")
                        }
		                   }
	              }, container = glyt4)
	  glyt4[2,1]<-gbutton(" Core Lines with Max RSF in the Selected Energy Range ", handler = function(h, ...) {
                     replot()
		                   if ( !is.na(svalue(be_value)) && !is.na(svalue(be_win)) ) {
			                     IdPeaks<<-CoreLinesMaxRSF(svalue(be_value), svalue(be_win), newcoreline, ElmtList)
			                     if (length(IdPeaks)>0){
   		                      ShowTablePeaks(IdPeaks, svalue(SpectType2))
                        } else {
                           cat("\n No element found with the selected precision!")
                        }
		                   }
	              }, container = glyt4)
	  glyt4[3,1]<-gbutton(" Find Any Element any Orbital in the Selected Energy Range ", handler = function(h, ...) {
                     replot()
	                    if ( !is.na(svalue(be_value)) && !is.na(svalue(be_win)) ) {
			                     IdPeaks<<-AllElements(svalue(be_value), svalue(be_win), newcoreline, ElmtList)
			                     if (length(IdPeaks)>0){
   		                      ShowTablePeaks(IdPeaks, svalue(SpectType2))
                        } else {
                           cat("\n No element found with the selected precision!")
                        }
		                   }
	              }, container = glyt4)

#--- TAB3

	gtab3 <- ggroup(label="Peak Table", horizontal=FALSE, spacing=5, container = nb)
   gframe31 <- gframe(text = "Peak Table", horizontal = FALSE, spacing=5, container = gtab3)
   gtable3 <- gtable(ElmtList, chosen.col=1, container = gframe31)
   size(gtable3) <- c(420,230)
   addHandlerDoubleclick(gtable3, handler = function(h, ...){
                     elmt<-svalue(gtable3)
                     HldPlt<-svalue(HoldPlot)
                     idx <- grep(elmt, ElmtList[,1])
                     elmtBE <- NULL
                     if (HldPlt==FALSE){
                        replot()    #refresh plot
                     }
                     for (ii in seq_along(idx)){
                         xx <- ElmtList[idx[ii],3]
#                         elmtBE <- c(elmtBE, ElmtList[idx[ii],3])
                         lines(x=c(xx, xx), y=rangeY, col="red")   #plot corelines of the selected elements
                     }

               })

   gframe32 <- gframe(text = "Spectral Type Selection", horizontal = TRUE, spacing=5, container = gtab3)
   SpectType3 <- gradio(c("CoreLines", "AugerTransitions"), selected=1, horizontal = TRUE, spacing=5,  handler = function(h,...){
                     SType <- svalue(SpectType3)
                     HldPlt<-svalue(HoldPlot)
                     ElmtList<<-ReadElmtList(SType) #reads the CoreLine/Auger List Table and save it in MyEnv  (macro in XPSSurveyUtility)
                     ElmtList<<-format(ElmtList, justify="centre", width=10)
                     ElmtList<<-as.data.frame(ElmtList,  stringsAsFactors = FALSE)
                     svalue(SpectType1)<-SType
                     svalue(SpectType2)<-SType
                     delete(gframe31, gtable3)
                     gtable3 <<- gtable(ElmtList, chosen.col=1, container = gframe31)
                     size(gtable3) <- c(420,230)
                     if(SType=="CoreLines") {
                        color="red"
                     } else if(SType=="AugerTransitions"){
                        color="blue"
                     }
                     addHandlerDoubleclick(gtable3, handler = function(h, ...){
                                       elmt<-svalue(gtable3)
                                       HldPlt<-svalue(HoldPlot)
                                       idx <- grep(elmt, ElmtList[,1])
                                       elmtBE <- NULL
                                       if (HldPlt==FALSE){
                                          replot()    #refresh plot
                                       }
                                       for (ii in seq_along(idx)){
                                           xx <- ElmtList[idx[ii],3]
                                           lines(x=c(xx, xx), y=rangeY, col=color)   #plot corelines of the selected elements
                                       }
                     })
                     add(gframe31, gtable3)
                     svalue(nb) <- 1     #force refresh of the notebook page
                     svalue(nb) <- 3
               }, container = gframe32)

    Search <-gedit(initial.msg="Element?", spacing=10, handler=function(h, ...){
                     elmt<-svalue(Search)
                     SType <- svalue(SpectType3)
                     HldPlt<-svalue(HoldPlot)
                     elmt<-paste(elmt, " ", sep="")
                     idx <- grep(elmt, ElmtList[,1])
                     if (HldPlt==FALSE){
                        replot()    #refresh plot
                     }
                     if(SType=="CoreLines") {
                        color="red"
                     } else if(SType=="AugerTransitions"){
                        color="blue"
                     }
                     for (ii in seq_along(idx)){
                         xx <- ElmtList[idx[ii],3]
                         lines(x=c(xx, xx), y=rangeY, col=color)   #plot corelines of the selected elements
                     }
               }, container=gframe32)
    tkconfigure(Search$widget, width=8)


    HoldPlot <-gcheckbox("Hold plot", checked = FALSE, spacing=10, container=gframe32)

#--- COMUNI

	CommFrame1 <- gframe(text = " Cursor and Zoom ", container = mainGroup, horizontal = FALSE)
    glabel("SX click to get position/define zoom area; DX to exit zoom", container = CommFrame1)

    glyt5<-glayout(homogeneous=FALSE, spacing=1, container=CommFrame1)
    glyt5[1,1]<-SetPos<-gbutton("Cursor Position", container = glyt5)
    addHandlerClicked(SetPos, handler=function(h,...){
                     pos<-locator(n=1, type="p", pch=3, cex=1.5, lwd=1.5, col="red") #to modify zoom limits
                     if (length(pos) > 0) { #non ho premuto tasto DX
                         replot()
                         points(pos, type="p", pch=3, cex=1.5, lwd=1.8, col="red")
                         pos<-round(x=c(pos$x, pos$y), digits=2)
                         txt<-paste("X: ", as.character(pos[1]), ", Y: ", as.character(pos[2]), sep="")
                         svalue(ShowPos)<-txt
                         tcl("update", "idletasks")  #force text to be shown in the glabel ShowPos
                     }
                  })

    glyt5[1,2]<-ShowPos<-glabel("Cursor position: ", container = glyt5)

    glyt6<-glayout(homogeneous=FALSE, spacing=1, container=CommFrame1)
    glyt6[1,1] <- gbutton("Set Zoom Boundaries", handler = function(h,...) {
	                 replot()
                  point.coords<<-locator(n=2, type="p", pch=3, cex=1.5, lwd=1.5, col="red")
                  DE<-abs(newcoreline[[1]][1]-newcoreline[[1]][2]) #Get the energy scale step: newcoreline[[1]] == newcoreline$x
	                 idx1 <- which(abs(newcoreline[[1]]-point.coords$x[1]) < DE) #get the element index corresponding to coords[1]
	                 idx1<-idx1[1] #in idx there could be more than one value: select the first
	                 idx2 <- which(abs(newcoreline[[1]]-point.coords$x[2]) < DE) #get the element index corresponding to coords[2]
	                 idx2<-idx2[1]
	                 point.coords$y<<-c(min(newcoreline[[2]][idx1:idx2]), max(newcoreline[[2]][idx1:idx2]))
                  slot(newcoreline,"Boundaries") <<- point.coords
	                 newcoreline <<- XPSsetRegionToFit(newcoreline)
		                replot()
             }, container = glyt6)

 	glyt6[1,2] <-gbutton("Reset Boundaries", handler = function(h, ...) {
  		              reset.boundaries()
  		              replot()
  	          }, container = glyt6)

   glyt6[1,3] <- gbutton("Refresh plot", handler = function(h,...) {
                  IdPeaks<<-NULL
                  enabled(identificationButton) <- TRUE
                  replayPlot(RecPlot)
             }, container = glyt6)

   glyt6[1,4] <- gbutton("Reset Analysis", handler = function(h,...) {
                  Peaks<<-NULL
                  IdPeaks<<-NULL
                  enabled(identificationButton) <- FALSE
                  replot()
                  RecPlot <<- recordPlot()   #save graph for Refresh
             }, container = glyt6)


 ## CLOSE button
  gseparator(container = mainGroup) # separator
  gbutton("Close", container = mainGroup, expand=FALSE, handler = function(h, ...) dispose(mainWin) )
  
  ## status
  sb <- gstatusbar("status", container = mainWin)
  svalue(sb) <- paste("Current Core Line :", as.character(CurrentCoreLine))
  plot(newcoreline)

  visible(mainWin) <- TRUE
  svalue(nb)<-3       #refresh notebook pages
  svalue(nb)<-2
  svalue(nb)<-1

}
