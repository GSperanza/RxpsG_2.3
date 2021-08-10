#Function to show the list of XPSSamples loaded in RXPSGUI and the list of their core lines

#'To select a XPS-Sample and a Core Line
#'
#'The list of XPS-Samples loaded is presented for selection.
#'After selection of the XPS-Sample the list of correspondent Corelines is available for selection.
#'No parameters are passed to this function
#'
#'
#'@examples
#'
#'\dontrun{
#'	XPSSetFNameCLine()
#'}
#'
#'@export
#'


XPSSetFNameCLine <- function() {

   updateObj <- function(h,...){
      SelectedFName<-svalue(T1obj1)
      FName=get(SelectedFName,envir=.GlobalEnv)  #carico in SampID il relativo XPSSAmple
      SpectList<<-XPSSpectList(SelectedFName)
      delete(T2frame1,T2obj1)
      T2obj1 <<- gcombobox(SpectList, selected=-1, editable=FALSE, handler=SetSpectrum, container=T2frame1)
      add(T2frame1,T2obj1)
      assign("activeFName", SelectedFName,envir=.GlobalEnv) #setto lo spettro attivo eguale ell'ultimo file caricato
      plot(FName)
   }



   SetSpectrum <- function(h,...){
      SpectName<-svalue(T2obj1)
      LL=length(SpectList)
      for (ii in 1:LL){   #dal nome del FileName riconosco lo SpectID ad esso associato
        if (SpectName==SpectList[ii]){
          indx<-ii
        }
      }

      assign("activeSpectName", SpectName,envir=.GlobalEnv) #setto lo spettro attivo eguale ell'ultimo file caricato
      assign("activeSpectIndx", indx, envir=.GlobalEnv) #setto lo spettro attivo eguale ell'ultimo file caricato

      #SampID e' un character devo caricare l'XPSSample dal GlobalEnv
      SampID=get(activeFName,envir=.GlobalEnv)  #carico in SampID il relativo XPSSAmple
      plot(SampID[[indx]])
   }




#carico la lista dei file ID e loro FileNames
   FNameList<-XPSFNameList()
   LL=length(FNameList)
   SampID<-""
   SpectList<-""

#===== NoteBook =====

   win <- gwindow("XPSsample & CORELINE SELECTION", visible=FALSE)
   nb <- gnotebook(expand=TRUE, container = win)

# --- Tab1 ---
   T1group1 <- ggroup(label=" SELECT XPSsample ", horizontal=FALSE, container=nb)
   T1layout<-glayout(container=T1group1)
   T1layout[1,1]<-T1frame1 <-gframe(" SELECT XPSsample ", spacing=5, container=T1layout)
   T1obj1 <- gcombobox(FNameList, selected=-1, editable=FALSE, handler=updateObj, container = T1frame1) # function(h,...){ }

   T1layout[2,1]<-T1obj2<-gbutton("     EXIT      ", handler=function(h,...){
                    dispose(win)
                    XPSSaveRetrieveBkp("save")
                 }, container = T1layout)

# --- Tab2 ---
   T2group1 <- ggroup(label=" SELECT CORELINE ",horizontal=FALSE, container=nb)
   T2layout<-glayout(container=T2group1)
   T2layout[1,1]<-T2frame1 <-gframe(text=" SELECT CORELINE ", spacing=5, container=T2layout)
   T2obj1 <- gcombobox(SpectList, selected=1, editable=FALSE, handler=SetSpectrum, container=T2frame1)

   T2layout[2,1]<-T2obj2<-gbutton("     EXIT      ", handler=function(h,...){
                    dispose(win)
                    XPSSaveRetrieveBkp("save")
                 }, container = T2layout)

   svalue(nb) <- 2     #refresh notebook pages
   svalue(nb) <- 1
   visible(win) <- TRUE

}