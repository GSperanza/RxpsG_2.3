# function to select the object to remove from the CoreLine spectral analysis
# regionToFit
# baseline
# components
# the whole analysis

# this function uses the function XPSremove() to remove the selected object

#'Resets a XPS-Coreline processing
#'
#'Function to delete individual elements of a Coreline analysis:
#'resets only the Best Fit;
#'resets an single or all the Fitting Components;
#'resets the BaseLine (i.e. eliminates alsio the Fit);
#'resets ALL alle the Coreline analysis.
#'no parameters are passed to this function
#'
#'
#'@examples
#'
#'\dontrun{
#'	XPSResetAnalysis()
#'}
#'
#'@export
#'


XPSResetAnalysis <- function(){

   if (is.na(activeFName)){
       gmessage("No data present: please load and XPS Sample", title="XPS SAMPLES MISSING", icon="error")
       return()
   }
   FNameList <- XPSFNameList()
   ActiveFName <- get("activeFName", envir=.GlobalEnv)  #load the XPSSample
   FName <- get(activeFName, envir=.GlobalEnv)
   SpectIndx <- get("activeSpectIndx", envir=.GlobalEnv)
   ActiveSpectName <- get("activeSpectName", envir=.GlobalEnv)
   SpectList <- XPSSpectList(ActiveFName)
   NComp <- length(FName[[SpectIndx]]@Components)
   FitComp <- ""
   if (NComp > 0){
      FitComp <- names(FName[[SpectIndx]]@Components)  #Define a vector containing the Component names of the Active Coreline Fit
   }


#===== NoteBook =====

   RstWin <- gwindow("RESET ANALYSIS", visible=FALSE)
   size(RstWin) <- c(500, 230)
   RstGroup <- ggroup(container = RstWin, horizontal=FALSE, label = "ENERGY SHIFT")

   nb <- gnotebook(expand=TRUE, container = RstGroup)

# --- Tab0 ---
   T0group1 <- ggroup(label=" XPS-SAMPLE & CORELINE SELECTION ", spacing=5, horizontal=FALSE, container=nb)
   T0group2 <- ggroup(label="  ",horizontal=TRUE, container=T0group1)
   T0frame2 <- gframe(" Select XPS-Sample", spacing=5, container=T0group2)
   T0XPSSample <- gcombobox(FNameList, selected=-1, editable=FALSE, handler=function(h,...){
                                 ActiveFName <<- svalue(T0XPSSample)
                                 FName <<- get(ActiveFName, envir=.GlobalEnv)
                                 SpectList <- XPSSpectList(activeFName)
                                 delete(T0frame3, T0CoreLine)
                                 T0CoreLine <<- gcombobox(SpectList, selected=-1, editable=FALSE, handler=function(h,...){
                                                          XPSComponent <- svalue(T0CoreLine)
                                                          XPSComponent <- unlist(strsplit(XPSComponent, "\\."))   #skip the  "NUMBER." at the beginning of the Coreline Name
                                                          SpectIndx <<- as.integer(XPSComponent[1])
                                                          SpectName <- XPSComponent[2]
                                                          if (length(FName[[SpectIndx]]@RegionToFit)==0){
                                                              gmessage(msg="ATTENTION: no defined fit region on this spectrum", title="WARNING!", icon="warning")
                                                              return()
                                                          }

                                                          assign("activeSpectName", SpectName,.GlobalEnv) #set active filename == last loaded XPSSample
                                                          assign("activeSpectIndx", SpectIndx,.GlobalEnv)
                                                          delete(T1frame1, T1obj1)
                                                          FitComp <- names(FName[[SpectIndx]]@Components)   #update component list of the combobox Notebook first page
                                                          if (length(FitComp) > 0){
                                                              FitComp <<- names(FName[[SpectIndx]]@Components)  #define the vector with the component names of the Active Coreline Fit
                                                              T1obj1 <<- gcombobox(FitComp, selected=-1, editable=FALSE, handler=function(h,...){
                                                                                   if (NComp == 0) {
                                                                                      gmessage(masg = "ATTENTION: NO FIT FOUND. CANNOT REMOVE COMPONENTS", icon = "warning", title="WARNING!")
                                                                                   }
                                                                         }, container=T1frame1)
                                                          } else {
                                                              FitComp <<- ""
                                                          }
                                                          add(T1frame1, T1obj1)
                                                          plot(FName[[SpectIndx]])
                                 }, container = T0frame3)
                                 plot(FName)
                       }, container = T0frame2)
   svalue(T0XPSSample) <- ActiveFName

   T0frame3 <- gframe(" Select Coreline ", spacing=5, container=T0group2)
   T0CoreLine <- gcombobox(SpectList, selected=-1, editable=FALSE, handler=function(h,...){
                                 XPSComponent <- svalue(T0CoreLine)
                                 XPSComponent <- unlist(strsplit(XPSComponent, "\\."))   #skip the  "NUMBER." at the beginning of the Coreline Name
                                 SpectIndx <<- as.integer(XPSComponent[1])
                                 SpectName <- XPSComponent[2]
                                 if (length(FName[[SpectIndx]]@RegionToFit)==0){
                                    gmessage(msg="ATTENTION: no defined fit region on this spectrum", title="WARNING!", icon="warning")
                                    return()
                                 }

                                 assign("activeSpectName", SpectName,.GlobalEnv) #set active filename == last loaded XPSSample
                                 assign("activeSpectIndx", SpectIndx,.GlobalEnv)
                                 delete(T1frame1, T1obj1)
                                 FitComp <- names(FName[[SpectIndx]]@Components)   #update component list of the combobox Notebook first page
                                 if (length(FitComp) > 0){
                                    T1obj1 <<- gcombobox(FitComp, selected=-1, editable=FALSE, handler=function(h,...){
                                                        if (NComp == 0) {
                                                          gmessage(masg = "ATTENTION: NO FIT FOUND. CANNOT REMOVE COMPONENTS", icon = "warning", title="WARNING!")
                                                        }
                                               }, container=T1frame1)
                                 }
                                 add(T1frame1, T1obj1)
                                 plot(FName[[SpectIndx]])
                                 NComp <<- length(FName[[SpectIndx]]@Components)
                                 FitComp <<- ""
                                 if (NComp > 0){
                                    FitComp <<- names(FName[[SpectIndx]]@Components)  #define the vector with the component names of the Active Coreline Fit
                                 }

                       }, container = T0frame3)
   svalue(T0CoreLine) <- paste(SpectIndx,".",ActiveSpectName, sep="")


# --- Tab1 ---
   T1group1 <- ggroup(label="  REMOVE FIT COMPONENT  ", selected=-1, horizontal=FALSE,container=nb)
   layout1 <- glayout(homogeneous=FALSE, spacing=3, container=T1group1)
   layout1[1,1] <- T1frame1 <- gframe(text="  SELECT FIT COMPONENT  ", spacing=5, container=layout1)

   T1obj1 <- gcombobox(FitComp, selected=-1, editable=FALSE, handler=function(h,...){
                            if (NComp == 0) {
                               gmessage(masg = "ATTENTION: NO FIT FOUND. CANNOT REMOVE COMPONENTS", icon = "warning", title="WARNING!")
                            }
                       },container=T1frame1)

   layout1[2,1] <- gbutton("     OK      ", expand=FALSE, handler=function(h,...){
                               answ <- gconfirm(msg = "ATTENTION: all the fit constraints will be lost! Proceed anyway?", icon = "warning", title="WARNING!")
                               if (answ) {
                                  comp <- svalue(T1obj1)
                                  object=as.integer(substr(comp,2,3)) #drops the first character (C) of the component name and take the following one or two
                                  FName[[SpectIndx]] <<- XPSremove(FName[[SpectIndx]],"components",object)
                                  NComp <<- length(FName[[SpectIndx]]@Components)
                                  if (NComp > 0){   #at least one fit component is still present
                                     tmp <- sapply(FName[[SpectIndx]]@Components, function(z) matrix(data=z@ycoor)) #compute the Fit without the selected component
	                                  FName[[SpectIndx]]@Fit$y <<- ( colSums(t(tmp)) - length(FName[[SpectIndx]]@Components)*(FName[[SpectIndx]]@Baseline$y)) #remove baseline
	                                  FName[[SpectIndx]]@Fit$fit <<- list()      #remove fit info
                                     FitComp <- names(FName[[SpectIndx]]@Components)   #update the component list
                                  } else {
                                     FName[[SpectIndx]]@Components <<- list()
                                     FName[[SpectIndx]]@Fit <<- list()
                                  }
                                  delete(T1frame1,T1obj1)  #update the coreline combobox
                                  T1obj1 <<- gcombobox(FitComp, selected=-1, editable=FALSE, handler=function(h,...){
                                                 if (NComp == 0) {
                                                    gmessage(msg = "ATTENTION: NO FIT FOUND. CANNOT REMOVE COMPONENTS", icon = "warning", title="WARNING!")
                                                 }
                                              },container=T1frame1)
                                  plot(FName[[SpectIndx]])
                               } else {
                                  return()
                               }
                  }, container = layout1)


# --- Tab2 ---
   T2group1 <- ggroup(label=" REMOVE FIT  ", horizontal=FALSE, container=nb)
   layout2 <- glayout(homogeneous=FALSE, spacing=3, container=T2group1)
   layout2[1,1] <- gbutton("     OK      ", handler=function(h,...){
                     if (NComp == 0) {
                        tkmessageBox(message = "WARNING: NO FIT FOUND. CANNOT REMOVE ANY FIT", icon = "warning", type = "ok")
                     } else {
                        FName[[SpectIndx]] <<- XPSremove(FName[[SpectIndx]],"fit")
                        FName[[SpectIndx]] <<- XPSremove(FName[[SpectIndx]],"components")
                        plot(FName[[SpectIndx]])
                     }
                  }, container = layout2)

# --- Tab3 ---
   T3group1 <- ggroup(label=" RESET ANALYSIS ", horizontal=FALSE, container=nb)
   layout3 <- glayout(homogeneous=FALSE, spacing=3, container=T3group1)
   layout3[1,1] <- gbutton("     OK      ", handler=function(h,...){
                     if (NComp == 0 && length(FName[[SpectIndx]]@RegionToFit)==0 && length(FName[[SpectIndx]]@Baseline)==0) {
                        gmessage(title="WARNING!", msg = "ATTENTION: NO BASELINE FOUND. Cannot remove anything...", icon = "warning")
                     } else {
                        FName[[SpectIndx]] <<- XPSremove(FName[[SpectIndx]],"all")
                        plot(FName[[SpectIndx]])
                     }
                  }, container = layout3)

# --- Common Buttons ---
   layout4 <- glayout(homogeneous=FALSE, spacing=3, container=RstGroup)

   layout4[1,1] <- gbutton("  SAVE ", handler=function(h,...){
                               assign(ActiveFName, FName, envir=.GlobalEnv)
                               XPSSaveRetrieveBkp("save")
                           }, container=layout4)

   layout4[1,2] <- gbutton("  REFRESH ", handler=function(h,...){
                               FName <<- get(activeFName, envir=.GlobalEnv)  #reload the CoreLine
                               FitComp <- names(FName[[SpectIndx]]@Components)   #Update component list in the combobox I Notebook page
                               delete(T1frame1,T1obj1)
                               T1obj1 <<- gcombobox(FitComp, selected=-1, editable=FALSE, handler=function(h,...){
                                                 if (NComp == 0) {
                                                    gmessage(masg = "ATTENTION: NO FIT FOUND. CANNOT REMOVE COMPONENTS", icon = "warning", title="WARNING!")
                                                 }
                                              },container=T1frame1)
                               add(T1frame1,T1obj1)
                               plot(FName[[SpectIndx]])
                           }, container=layout4)

   layout4[1,3] <- gbutton("  SAVE and EXIT ", handler=function(h,...){
                               assign(ActiveFName, FName, envir=.GlobalEnv)
                               dispose(RstWin)
                               XPSSaveRetrieveBkp("save")
                           }, container=layout4)



   visible(RstWin) <- TRUE
   for(ii in 4:1){
      svalue(nb) <- ii
   }
}
