#Function to select Corelines for Fit Constraints

#'To select a Core Line
#'
#'The list of corelines of a give XPS-Sample are presented for selection
#'no parameters are passed to this function
#'
#'
#'@examples
#'
#'\dontrun{
#'	XPSSetCoreLine()
#'}
#'
#'@export
#'


XPSSetCoreLine <- function() {

#carico la lista dei file ID e loro FileNames
   ActiveFName<-get("activeFName", envir=.GlobalEnv)
   SpectList<-XPSSpectList(ActiveFName)


#===== CL SELECTION  FOR FIT CONSTRAINTS =====

   mainFCwin<-get("C.mainFCwin", mainFCwin, envir=MyEnv)
   CLgroup <- ggroup(label="", horizontal=FALSE, container=mainFCwin)

   CLframe <-gframe(text=" CORELINE SELECTION ", spacing=5, container=CLgroup)
   CLobj <- gcombobox(SpectList, selected=-1, editable=FALSE, handler=function(h,...){
                           XPSComponent<-svalue(CLobj)
                           XPSComponent<-unlist(strsplit(XPSComponent, "\\."))   #tolgo il "NUMERO." all'inizio del nome coreline
                           indx<-as.integer(XPSComponent[1])
                           SpectName<-XPSComponent[2]
                           assign("activeSpectName", SpectName,.GlobalEnv) #setto lo spettro attivo eguale ell'ultimo file caricato
                           assign("activeSpectIndx", indx,.GlobalEnv) #setto lo spettro attivo eguale ell'ultimo file caricato
                           #FName e' un character devo caricare l'XPSSample dal GlobalEnv
                           FName=get(activeFName,envir=.GlobalEnv)  #carico in FName il relativo XPSSAmple dataframe
                           plot(FName[[indx]])

                 }, container = CLframe)

   CLobj2<-gbutton("     SELECT    ", handler=function(h,...){
                    delete(mainFCwin,CLgroup)
                    XPSConstraints()
                 }, container = CLframe)
                 

#   visible(win) <- TRUE     NO GIA' ACCESA IN XPSSetMainFCWinGUI.r
}
