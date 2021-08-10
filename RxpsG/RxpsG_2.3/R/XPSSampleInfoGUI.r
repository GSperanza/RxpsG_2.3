# Per recuperare e modificare le INFO dell'XPSSample
# revised October 2014

#'XPSSampleInfo to show/modify XPSSample INFOs is saved during acquisition
#'
#'@examples
#'
#'\dontrun{
#'	XPSSampleInfo()
#'}
#'
#'@export
#'


XPSSampleInfo <- function() {
      if (is.na(activeFName)){
          gmessage("No data present: please load and XPS Sample", title="XPS SAMPLES MISSING", icon="error")
          return()
      }

#--- Variable definition ---
      FName <- get(activeFName, envir=.GlobalEnv)
      CLineList <- XPSSpectList(activeFName)
      if (is.na(CLineList)){
         gmessage("ATTENTION NO CORELINES FOUND: please control your XPSSample datafile!" , title = "WARNING",  icon = "warning")
      }
      Data <- list()
      #read the XPS Sample information
      Data[1] <- FName@Project
      Data[2] <- FName@Sample
      Data[3] <- paste(FName@Comments, collapse=" ")
      Data[4] <- FName@User
      Data[5] <- paste(FName@names, collapse=" ")
      for(ii in 1:5){
         if ( is.na(Data[[ii]]) || is.null(Data[[ii]])) { Data[[ii]] <- "---" }
      }
      Data[[6]] <- "                                                  "  #add a row of 50 spaces to expand the GDF()
      VarNames <- c("Project", "Sample", "Comments", "User", "names", "  ")
      Data <- data.frame(INFO=cbind(VarNames,Data), stringsAsFactors=FALSE) #gdf() add a column to display the row names
      newData <- Data

#--- GUI ---
      DFwin <- gwindow(title="XPS-SAMPLE INFO", visible=FALSE) #define the main window to display the gdf()
      size(DFwin) <- c(600,500)
      DFgroup <- ggroup(horizontal=FALSE, container=DFwin)
      txt <- paste("   ", activeFName, " :   EDIT XPS-SAMPLE INFO", sep="")
      glabel(txt, container=DFgroup) #label to extend the window dimensions

      DFrame <- gdf(items=Data, container=DFgroup)
      addHandlerChanged(DFrame, handler=function(h,...){ #addHandlerChanged dowloads the dataFrame with modified parameters
                      newData <<- h$obj[]
                  })

      size(DFrame) <- c(600,150)

      if (length(CLineList)>1){      #gradio works with at least 2 items
          SelectCL <- gradio(CLineList, selected=-1, horizontal=TRUE, handler= function(h, ...){
                      idx <- svalue(SelectCL, index=TRUE)
                      txt <- NULL
                      txt[1] <- paste("Core Line : ",slot(FName[[idx]],"Symbol"),"\n", sep="")
                      txt[2] <- paste("N. data   : ",length(FName[[idx]]@.Data[[1]]))
                      txt[3] <- paste("baseline  : ",ifelse(hasBaseline(FName[[idx]]),"YES", "NO"),"\n", sep="")
                      txt[4] <- paste("fit       : ",ifelse(hasFit(FName[[idx]]),"YES", "NO"),"\n", sep="")
                      txt[5] <- paste("n. comp.  : ",ifelse(hasComponents(FName[[idx]]),length(FName[[idx]]@Components), "NONE"),"\n", sep="")
                      txt[6] <- (" Info:\n")
                      info<<- FName[[idx]]@Info
                      svalue(CLinfo) <- c(txt,info)
          }, container=DFgroup)
      } else {                       #if there is just 1 coreline in the XPSSample then use gcheckboxgroup()
          SelectCL <- gcheckboxgroup(CLineList, checked=FALSE, handler= function(h, ...){
                      txt <- NULL
                      txt[1] <- paste("Core Line : ",slot(FName[[1]],"Symbol"),"\n", sep="")
                      txt[2] <- paste("N. data   : ",length(FName[[1]]@.Data[[1]]))
                      txt[3] <- paste("baseline  : ",ifelse(hasBaseline(FName[[1]]),"YES", "NO"),"\n", sep="")
                      txt[4] <- paste("fit       : ",ifelse(hasFit(FName[[1]]),"YES", "NO"),"\n", sep="")
                      txt[5] <- paste("n. comp.  : ",ifelse(hasComponents(FName[[1]]),length(FName[[1]]@Components), "NONE"),"\n", sep="")
                      txt[6] <- (" Info:\n")
                      info<<- FName[[1]]@Info
                      svalue(CLinfo) <- c(txt,info)
          }, container=DFgroup)
      }

      CLinfo <- gtext(container=DFgroup)
      size(CLinfo) <- c(600,280)

      DFlayout <- glayout(homogeneous=FALSE, spacing=5, container=DFgroup)
      DFlayout[2,1] <- gbutton(" SAVE & EXIT ", handler=function(h,...){
                      FName@Project <<- newData[[2]][1]
                      FName@Sample <<- newData[[2]][2]
                      FName@Comments <<- newData[[2]][3]
                      FName@User <<- newData[[2]][4]
                      FName@names <<- unlist(strsplit(newData[[2]][5], " "))      #to correctly save the CoreLine Names
                      assign(activeFName, FName, envir=.GlobalEnv)
                      dispose(DFwin)
                      XPSSaveRetrieveBkp("save")
                  }, container=DFlayout)

      visible(DFwin) <- TRUE
}
