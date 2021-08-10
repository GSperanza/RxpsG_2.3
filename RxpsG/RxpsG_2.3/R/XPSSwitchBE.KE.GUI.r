#'function to set the energy scale i <<-n a selected XPSSample
#'
#'The energy scale in a XPSSample may be set as Binding Energy or Kinetic Energy
#'This function converts Binding energy into Kinetic and viceversa. The selected energy scale
#'for a given core line will be used  to plot the correspondent spectral data.
#'No parameters are passed to this function
#'
#'
#'@examples
#'
#'\dontrun{
#'	XPSSetWinSize()
#'}
#'
#'@export
#'


XPSSwitch.BE.KE <- function() {

#--- Variabili
   if (is.na(activeFName)){
       gmessage("No data present: please load and XPS Sample", title="XPS SAMPLES MISSING", icon="error")
       return()
   }
   EnergyScale <- c("Binding", "Kinetic")
   FNameList <- XPSFNameList()  #list of the XPSSample Names loaded in the .GlobalEnv
   SpectList <- "" #List of Corelines in the XPSSample
   FName <- NULL #FName represents the selected XPSSample
   XPSSample <- NULL
   SpectIndx <- NULL
   SpectName <- NULL
   Escale <<- NULL


#--- GUI
   Ewin <- gwindow("ENERGY SCALE GUI", visible=FALSE)
   size(Ewin) <- c(300, 320)
   Egroup <- ggroup(horizontal=FALSE, container=Ewin)

   Eframe1 <- gframe(" Select XPS-Sample ", horizontal=FALSE, spacing=5, container=Egroup)
   SelectXPSSamp <- gcombobox(FNameList, selected=-1, editable=FALSE, handler=function(h, ...){
                        FName <<- svalue(SelectXPSSamp)
                        XPSSample <<- get(FName,envir=.GlobalEnv)  #load the XPSSample datafile
                        SpectList <<- c(XPSSpectList(FName), "All")
                        SpectName <<- NULL
                        SpectIndx <<- NULL
                        plot(XPSSample)

                        LL <- length(XPSSample)
                        delete(Eframe2, SelectCL)  #now the CoreLines of FName are known. Redefine the selectCL object
                        SelectCL <<- gcombobox(SpectList, selected=-1, editable=FALSE, handler=function(h, ...){ #ridefinisco il combobox con i nome delle corelines
                               CoreLine <- svalue(SelectCL)
                               if (CoreLine=="All") {
                                  SpectName <<- "All"
                                  svalue(SetE) <<- XPSSample[[1]]@units[1]  #Energy unit for the generic coreline
                                  SpectIndx <<- NULL
                                  return()
                               }
                               CoreLine <- unlist(strsplit(CoreLine, "\\."))   #Skip the X. number prior to the CoreLine Name
                               SpectIndx <<- as.integer(CoreLine[1])
                               SpectName <<- CoreLine[2]
                               Escale <- XPSSample[[SpectIndx]]@units[1]  #Energy unit for the selected coreline
                               Escale <- substr(Escale,1, 7)     #extract Binding (Kinetic) from "Binding Energy eV" (from "Binding Energy eV")
                               svalue(SetE) <- Escale
                               plot(XPSSample[[SpectIndx]])
                        }, container = Eframe2)
                    }, container=Eframe1)

   Eframe2 <- gframe("Select the Core Line", horizontal=FALSE, spacing=5, container=Egroup)
   SelectCL <- gcombobox(SpectList, selected=-1, editable=FALSE,container=Eframe2)

   Eframe3 <- gframe(text="Set the Energy Scale", spacing=5, container=Egroup)
   SetE <- gradio(items=EnergyScale, horizontal=TRUE, selected=1, handler=function(h,...){
                        XEnergy <- get("XPSSettings", envir=.GlobalEnv)$General[5] #the fifth element of the first column of XPSSettings: X radiation energy
                        XEnergy <- as.numeric(XEnergy)
                        if (is.null(FName)){
                           gmessage("Please select the XPS Sample", title="SELECTION OF XPS SAMPLE LACKING", icon="warning")
                        }
                        if (is.null(SpectName)){
                           gmessage("Please select a Core Line", title="SELECTION OF CORE LINE LACKING", icon="warning")
                        }
                        Escale <<- svalue(SetE)  #read Energy scale
                        if (SpectName=="All"){
                            LL <- length(XPSSample)
                            for(ii in 1:LL){     #run on all the XPSSample spectra
                               if (Escale=="Binding"){
                                  if (XPSSample[[ii]]@Flags[1]==FALSE){ #Kinetic energy set in the original XPSSample
                                      XPSSample[[ii]]@Flags[1] <<- TRUE
                                      XPSSample[[ii]]@units[1] <<- "Binding Energy [eV]"
                                  } else {
                                      return()
                                  }
                               }
                               if (Escale=="Kinetic"){
                                  if (XPSSample[[ii]]@Flags[1]==TRUE){ #Binding energy set in the original XPSSample
                                      XPSSample[[ii]]@Flags[1] <<- FALSE
                                      XPSSample[[ii]]@units[1] <<- "Kinetic Energy [eV]"
                                  } else {
                                      return()
                                  }
                               }
                               XPSSample[[ii]]@.Data[[1]] <<- XEnergy-XPSSample[[ii]]@.Data[[1]] #transform kinetic in binding
                               XPSSample[[ii]]@.Data[[2]] <<- rev(XPSSample[[ii]]@.Data[[2]])
                               XPSSample[[ii]]@Boundaries$x <<- XEnergy-XPSSample[[ii]]@Boundaries$x
                               if (length(XPSSample[[ii]]@RegionToFit)>0){
                                   XPSSample[[ii]]@RegionToFit$x <<- XEnergy-XPSSample[[ii]]@RegionToFit$x
                                   XPSSample[[ii]]@RegionToFit$y <<- rev(XPSSample[[ii]]@RegionToFit$y)
                               }
                               if (length(XPSSample[[ii]]@Baseline)>0){
                                   XPSSample[[ii]]@Baseline$x <<- XEnergy-XPSSample[[ii]]@Baseline$x
                                   XPSSample[[ii]]@Baseline$y <<- rev(XPSSample[[ii]]@Baseline$y)
                               }
                               if (length(XPSSample[[ii]]@Components)>0){
                                   LL <- length(XPSSample[[ii]]@Components)
                                   for(jj in 1:LL){    #transforms BE in KE for all the coreline of XPSSample parameter "mu"
                                       varmu <- getParam(XPSSample[[ii]]@Components[[jj]],variable="mu")
                                       varmu <- XEnergy-varmu
                                       XPSSample[[ii]]@Components[[jj]] <<- setParam(XPSSample[[ii]]@Components[[jj]], parameter=NULL, variable="mu", value=varmu)
                                       XPSSample[[ii]]@Components[[jj]]@ycoor <<- rev(XPSSample[[ii]]@Components[[jj]]@ycoor)
                                   }
                                   XPSSample[[ii]]@Fit$y <<- rev(XPSSample[[SpectIndx]]@Fit$y)
                               }
                            }
                            plot(XPSSample)
                        } else if (is.integer(SpectIndx)) {  #if a single Core Line was selected
                            if (Escale=="Binding"){
                               if (XPSSample[[SpectIndx]]@Flags[1]==FALSE){ #Kinetic energy set in the original XPSSample
                                   XPSSample[[SpectIndx]]@Flags[1] <<- TRUE
                                   XPSSample[[SpectIndx]]@units[1] <<- "Binding Energy [eV]"
                               } else {  #original Escale == Binding and option Binding selected
                                   return()
                               }
                            }
                            if (Escale=="Kinetic"){
                               if (XPSSample[[SpectIndx]]@Flags[1]==TRUE){ #Binding energy set in the original XPSSample
                                   XPSSample[[SpectIndx]]@Flags[1] <<- FALSE
                                   XPSSample[[SpectIndx]]@units[1] <<- "Kinetic Energy [eV]"
                               } else {  #original Escale == Kinetic and option Kinetic selected
                                   return()
                               }
                            }
                            XPSSample[[SpectIndx]]@.Data[[1]] <<- XEnergy-XPSSample[[SpectIndx]]@.Data[[1]] #transform kinetic in binding  and viceversa
                            XPSSample[[SpectIndx]]@.Data[[2]] <<- rev(XPSSample[[SpectIndx]]@.Data[[2]])
                            XPSSample[[SpectIndx]]@Boundaries$x <<- XEnergy-XPSSample[[SpectIndx]]@Boundaries$x
                            if (length(XPSSample[[SpectIndx]]@RegionToFit)>0){
                                XPSSample[[SpectIndx]]@RegionToFit$x <<- XEnergy-XPSSample[[SpectIndx]]@RegionToFit$x
                                XPSSample[[SpectIndx]]@RegionToFit$y <<- rev(XPSSample[[SpectIndx]]@RegionToFit$y)
                            }
                            if (length(XPSSample[[SpectIndx]]@Baseline)>0){
                                XPSSample[[SpectIndx]]@Baseline$x <<- XEnergy-XPSSample[[SpectIndx]]@Baseline$x
                                XPSSample[[SpectIndx]]@Baseline$y <<- rev(XPSSample[[SpectIndx]]@Baseline$y)
                            }
                            if (length(XPSSample[[SpectIndx]]@Components)>0){
                                LL <- length(XPSSample[[SpectIndx]]@Components)
                                for(jj in 1:LL){    #transforms BE in KE for all the coreline of XPSSample parameter "mu"
                                    varmu <- getParam(XPSSample[[SpectIndx]]@Components[[jj]],variable="mu")
                                    varmu <- XEnergy-varmu
                                    XPSSample[[SpectIndx]]@Components[[jj]] <<- setParam(XPSSample[[SpectIndx]]@Components[[jj]], parameter=NULL, variable="mu", value=varmu)
                                    XPSSample[[SpectIndx]]@Components[[jj]]@ycoor <<- rev(XPSSample[[SpectIndx]]@Components[[jj]]@ycoor)
                                }
                                XPSSample[[SpectIndx]]@Fit$y <<- rev(XPSSample[[SpectIndx]]@Fit$y)
                            }
                            plot(XPSSample[[SpectIndx]])
                        }
                    }, container=Eframe3)

   glabel(" \n \n", container=Egroup)   #some space before buttons

   gbutton("SAVE", handler=function(h,...){
                       assign(FName, XPSSample, envir = .GlobalEnv)  #save changes in the originasl datafile
                    }, container = Egroup)

   gbutton("SAVE and EXIT", handler=function(h,...){
                       assign(FName, XPSSample, envir = .GlobalEnv)  #save changes in the originasl datafile
                       dispose(Ewin)
                    }, container = Egroup)


   visible(Ewin) <- TRUE
}
