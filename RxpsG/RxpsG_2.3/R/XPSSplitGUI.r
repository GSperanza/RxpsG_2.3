# GUI to split an XPSSample in groups its corelines
# The GUI shows a list of the XPSSamples loaded
# For the selected XPSSample the GUI shows the list of Corelines acquired
# A group of corelines can be chosen through a checkbox
# The chosen corelines may be saved in a new XPSSample.

#'To select split multiple acquisitions performed with the ScientaGammadata instrument
#'
#'Acquisitions on multiple samples may be included in a single .PXT file
#'This function allows splitting spectra corresponding to each sample in 
#'individual files.
#'No parameters are passed to this function
#'
#'
#'@examples
#'
#'\dontrun{
#'	XPSSplit()
#'}
#'
#'@export 
#'




XPSSplit <- function() {

      CKSlctdCL <- function(){ #Check which are the selected corelines in a column of checkboxes
             CheckedCL <- list()
             SelectedCL <<- NULL
             for(ii in 1:N.CKboxGrps){
                 CheckedCL[[ii]] <- svalue(CoreLineCK[[ii]])
             }
             for(ii in 1:N.CKboxGrps){
                SelectedCL <<- c(SelectedCL, unlist(CheckedCL[[ii]]))
             }
      }

      MakeSplitCL <- function(){  #construct the list of checkboxes describing all the XPSSample corelines
             CoreLines <- XPSSpectList(activeFName)
             LL <- length(CoreLines)
             jj <- 1
             for(ii in seq(from=1, to=LL, by=15) ){
                 NCL <- ii+15-1  #defines a set of 15 corelines (from ii to ii+15-1) for the jjth. checkbox group
                 if(LL < NCL ){  #if LL < ii+15 limits the CorelineList to LL
                    NCL <- LL
                 }
                 CoreLineList[[jj]] <<- CoreLines[ii:NCL]   #store CoreLines names in a list of vectors composed by 20-corelines each
                 CoreLineCK[[jj]] <<- gcheckboxgroup(CoreLineList[[jj]],checked=FALSE, handler=function(h, ...){
                                                        CKSlctdCL()
                                                     }, container=T1groupCoreLines)
                 jj <- jj+1
             }
             N.CKboxGrps <<- jj-1
      }



#--- Variables
     CoreLineList <- list()   #define a list for the XPSSample corelines
     CoreLineCK <- list()     #define a list of the Gwidget
     N.CKboxGrps <- NULL      #N. of checkboxgroups used to descibe the XPSSample corelines
     SelectedCL <- NULL
     CLidx <- NULL

     if (is.na(activeFName)){
         gmessage("No data present: please load and XPS Sample", title="XPS SAMPLES MISSING", icon="error")
         return()
     }

     FName <- get(activeFName, envir=.GlobalEnv)   #load the active XPSSample in memory
     activeFName <- get("activeFName", envir=.GlobalEnv)  #carico il nome XPSSample (stringa)
     FNameList <- XPSFNameList()     #list of all loaded XPSSamples in .GlobalEnv
     FNameIdx <- grep(activeFName, FNameList)

#--- GUI
     Swin <- gwindow(" SPLIT PXT FILEDATA ", visible=FALSE)
     maingroup <- ggroup(horizontal=TRUE, container=Swin)

     T1group1 <- ggroup(label="XPS SAMPLE SELECTION", spacing=5, horizontal=FALSE, container=maingroup)
     layoutT1 <- glayout(homogeneous=FALSE, spacing=5, container=T1group1)
     layoutT1[1,1] <- T1frameFName <- gframe(text="Select the XPS-SAMPLE", spacing=5, container=layoutT1)
     T1FNameListCK <- gradio(FNameList,selected=FNameIdx, handler = function(h, ...){
                                activeFName <<- svalue(T1FNameListCK)
                                assign("activeFName", activeFName, envir=.GlobalEnv)
                                FName <<- get(activeFName, envir=.GlobalEnv)
                                activeFName <- get("activeFName", envir=.GlobalEnv)  #load the active XPSSample in memory
                                FNameList <<- XPSFNameList()     #list of all loaded XPSSamples in .GlobalEnv
                                FNameIdx <<- grep(activeFName, FNameList)
                                delete(T1frameCoreLines,T1groupCoreLines) #reset the frame containinf the checkbox-corelines
                                T1groupCoreLines <<- ggroup(horizontal=TRUE,container = T1frameCoreLines)
                                MakeSplitCL()
                                plot(FName)
                             }, container=T1frameFName)

     layoutT1[1,2] <- T1frameCoreLines <- gframe(text="Select the CORE LINES to export",  spacing=5, container=layoutT1)
     T1groupCoreLines <- ggroup(horizontal=TRUE,container = T1frameCoreLines)
     MakeSplitCL()

     layoutT1[3,1] <- Export <- gbutton("SELECT AND SAVE SPECTRA", handler=function(h,...){
#--- building the NewXPSSample datafile
                            spectNames <- NULL
                            CLname <- NULL
                            CLidx <<- NULL
	                           NewFName <- new("XPSSample")
                            LL <- length(SelectedCL)
                            mm <- length(FName)
                            for (ii in 1:LL){
	                               NewFName[[ii]] <- new("XPSCoreLine")
                                CLname <- unlist(strsplit(SelectedCL[ii], "\\."))   #skip the number at coreline name beginning
                                CLidx[ii] <<- as.integer(CLname[1])
                                spectNames[ii] <- CLname[2]
                                NewFName[[ii]] <- FName[[CLidx[ii]]]
                            }
                            NewFName@Project <- FName@Project
                            NewFName@Sample <- FName@Sample
                            NewFName@Comments <- FName@Comments
                            NewFName@User <- FName@User
                            NewFName@names <- spectNames
                            plot(NewFName)

                            PathFile <- gfile(type=c("save"), initial.filename="*.RData", initial.dir=getwd()) #select folder and name to save data
                            if( length(PathFile) == 0 ){ return() }
                            PathName <- dirname(PathFile)
                            FileName <- basename(PathFile) #extract the filename from complete
                            FileName <- unlist(strsplit(FileName, "\\."))
                            if (is.na(FileName[2]) || FileName[2] != "RData"){
                               gmessage("Extension of the destination file forced to .RData!" , title = "DESTINATION FILE EXTENSION",  icon = "warning")
                            }
                            PathFile <- paste(PathName, "/", FileName[1], ".", "RData", sep="")
                            NewFName@Sample <- PathFile
                            FileName <- paste(FileName[1], ".", "RData", sep="")
                            NewFName@Filename <- FileName
                            assign(FileName,NewFName, envir=.GlobalEnv) #salvo il nuovo XPSSample nel GlobalEnvir
                            command=paste("save('", FileName,"', file='",PathFile, "', compress=TRUE)", sep="")
                            eval(parse(text=command),envir=.GlobalEnv)
                            cat("\n Data saved in: ", PathFile)
                            XPSSaveRetrieveBkp("save")

     }, container=layoutT1)

     layoutT1[3,2] <- Clear <- gbutton(" CLEAR SELECTIONS ", handler=function(h,...){  #Clear checkboxes and mark selected spectra with #
                           delete(T1frameCoreLines, T1groupCoreLines)
                           T1groupCoreLines <<- ggroup(horizontal=TRUE,container = T1frameCoreLines) #clear all checkboxes
                           for(jj in 1:N.CKboxGrps){
                               minIdx <- (jj-1)*15+1 #minIdx, maxIdx indicate the range of the CLindex (number preceding the spectname)
                               maxIdx <- jj*15
                               Idx <- CLidx-minIdx+1 #shift CLidx to match the CoreLineList[[jj]] components
                               LL <- length(Idx)
                               if (min(Idx) > 0 && min(Idx) < 15 ){ #all o part of the selected spectra belong to CoreLineList[[jj]]
                                  for(ii in 1:LL){
                                      if (Idx[ii] > 15) { break }  #we are across two CoreLine lists: consider only the part belonging to CoreLineList[[jj]]
                                      CoreLineList[[jj]][Idx[ii]] <<- paste("=>", CoreLineList[[jj]][Idx[ii]], sep="") #change the selected coreline name with the "#"
                                  }
                               }
                               if (min(Idx) < 0){ #second part of the selected spectra belong to CoreLineList[[jj+1]]
                                  for(ii in 1:LL){
                                      if (Idx[ii] > 15) { break }  #we are across two CoreLine lists: consider only the part belonging to CoreLineList[[jj]]
                                      if (Idx[ii] > 0){
                                          CoreLineList[[jj]][Idx[ii]] <<- paste("=>", CoreLineList[[jj]][Idx[ii]], sep="") #change the selected coreline name with the "#"
                                      }
                                  }
                               }
                               #Now regenerate the checkbox columns
                               CoreLineCK[[jj]] <<- gcheckboxgroup(CoreLineList[[jj]],checked=FALSE, handler=function(h, ...){
                                                                     CKSlctdCL()
                                                                  }, container=T1groupCoreLines)
                           }
     }, container=layoutT1)

     layoutT1[4,1] <- Exit <- gbutton(" EXIT ", handler=function(h,...){
    	                     dispose(Swin)
    	                     return(1)
     }, container=layoutT1)

     visible(Swin) <- TRUE
     Swin$set_modal(TRUE)

}


