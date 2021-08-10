#'GUI to change the name of the XPSSample and/or the name associated to corelines
#'
#'
#'@examples
#'
#'\dontrun{
#'	XPSSpectNameChange()
#'}
#'
#'@export
#'


XPSSpectNameChange <- function(){


#--- Global variables definition ---
   if (is.na(activeFName)){
       gmessage("No data present: please load and XPS Sample", title="XPS SAMPLES MISSING", icon="error")
       return()
   }
   FNameList <- XPSFNameList()  #list of the XPSSample loaded in the Global Env
   SpectList <- "" #list of corelines present in the selected XPSSample
   FName <- NULL #FName represents the selected XPSSample(class XPSSample)
   activeFName <- NULL #activeFName is the name associated to the selected XPSSample (class character)

#--- Widget  
      
      Labwin <<- gwindow("CHANGE LABELS", visible=FALSE)

      Labgroup1 <- ggroup(container = Labwin, horizontal=FALSE)
      Labframe1 <- gframe(" Select XPS-Sample ", horizontal=FALSE, spacing=5, container=Labgroup1)
      Labobj1 <- gcombobox(FNameList, selected=-1, editable=FALSE, handler=function(h, ...){
                          activeFName<<- svalue(Labobj1)  #save the XPSSample name
                          FName <<- get(activeFName, envir=.GlobalEnv)  #load in FName the XPSSample data
                          OldSpectName <- names(FName)
                          SpectList <- data.frame("Names"=OldSpectName, stringsAsFactors=FALSE)

                          plot(FName)

                          svalue(Labobj5) <- FName@Filename
                          LL <- length(FName)
                          delete(Labframe2, LabGDF)
                          LabGDF <<- gdf(items=SpectList, container=Labframe2) #ridefine the combobox with the coreline names

                          size(LabGDF) <- c(80,100)
                          addHandlerChanged(LabGDF, handler=function(h, ...){
                                              idx <- svalue(LabGDF, drop=TRUE) <- ""
                                              NewSpectName <- as.character(h$obj[])
                                              names(FName) <<- NewSpectName
                                              for (ii in 1:LL){
                                                  if (NewSpectName[ii] != OldSpectName[ii]){
                                                     FName[[ii]]@Symbol <<- NewSpectName[ii]
                                                     FName[[ii]]@RSF <<- 0 #otherwise XPSClass does not set RSF (see next row)
                                                     FName[[ii]] <<- XPSsetRSF(FName[[ii]])
                                                  }
                                              }
                                              svalue(Labobj3) <- " Spectrum Label and RSF changed!"
                                          })
                     }, container=Labframe1)

      Labframe2 <- gframe(" Change Spectrum Label Change ", horizontal=FALSE, spacing=5, container=Labgroup1)
      Labobj3 <- glabel(" ", container=Labframe2)
      LabGDF <- gdf(items=SpectList, container=Labframe2)
      size(LabGDF) <- c(80,100)


      Labframe4 <- gframe(" Set the New XPS-Sample Name ", horizontal=FALSE, spacing=5, container=Labgroup1)
      Labobj5 <- gedit("", handler=function(h,...){
                          XPSSampName <- svalue(Labobj1)
                          if (length(XPSSampName)==0){
                             gmessage(msg="Please Select XPS-Sample name", title="WARNING", icon="warning")
                          }
                          # see rows  45 - 67: Labobj5 redefinition
                     }, container=Labframe4)
      Labobj6 <- glabel("", container=Labframe4)

      gbutton(" SAVE ", container=Labgroup1, handler=function(...){
       	                 assign(activeFName, FName, envir=.GlobalEnv)
                          svalue(Labobj3) <- "  "
                          plot(FName)
                          XPSSaveRetrieveBkp("save")
                     })

      gbutton(" SAVE and EXIT ", container=Labgroup1, handler=function(...){
       	                 assign(activeFName, FName, envir=.GlobalEnv)
                          plot(FName)
      	                 dispose(Labwin)
      	                 XPSSaveRetrieveBkp("save")
                     })

      visible(Labwin) <- TRUE
      Labwin$set_modal(TRUE) 

}
