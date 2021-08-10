#Gui to set fit param MaXiTERATION, TOLERANCE, MINFACTOR, when NLSfit error occurs

#'XPSFitParam Gui to set fit param MaXiTERATION, TOLERANCE, MINFACTOR, when NLSfit error occurs
#'
#'@param ParamNames c("maxiteration", "tolerance", "minfactor") the set of parameter to work with
#'@param ParamValues the correspondent values
#'
#'@examples
#'\dontrun{
#'XPSFitParam(c("maxiteration", "tolerance", "minfactor"), c(1000, 1e-5, 1e-3))
#'}
#' 
#'@export
#'




XPSFitParam <- function(ParamNames, ParamValues){

            mainNLSwin<-gwindow("mainNLSwin", visible=FALSE)
            NLSgroup <- ggroup(label="", horizontal=FALSE, container=mainNLSwin)
            NLSframe1 <-gframe(text=" Select NLS parameters ", spacing=5, container=NLSgroup)
            NLSparam <- gcombobox(ParamNames, selected=-1, editable=FALSE, handler=function(h,...){
                           indx<-svalue(NLSparam, index=TRUE)
                           svalue(NLSobj1)<-ParamValues[indx]
                 }, container = NLSframe1)

            NLSframe2 <-gframe(text=" SET NLS parameter ", spacing=5, container=NLSgroup)
            NLSobj1 <- gedit(, selected=-1, editable=FALSE, handler=function(h,...){
                           indx<-svalue(NLSparam, index=TRUE)
                           ParamValues[indx] <<- as.numeric(svalue(NLSobj1))
                 }, container=NLSframe2)

            NLSframe3 <-gframe(text=" SET Parameters ", horizontal=FALSE, spacing=5, container=NLSgroup)
            NLSobj2<-gbutton("  SAVE Parameter   ", handler=function(h,...){
                    assign ("N.ParamValues", ParamValues, envir=MyEnv)
                 }, container = NLSframe3)

            NLSobj3<-gbutton("  Close and REFIT    ", handler=function(h,...){
                    fit<-NULL
                    tryAgain<-TRUE
                    assign("N.fit", fit, envir=MyEnv)   # set fit==NULL to repeat the fit if loop active i.e. tryAgain==FALSE
                    assign("N.tryAgain", tryAgain , envir=MyEnv)   # set tryAgain to loop
                    dispose(mainNLSwin)
                    return()
                 }, container = NLSframe3)

            NLSobj3<-gbutton("  EXIT FIT ROUTINE    ", handler=function(h,...){
                    tryAgain<-FALSE
                    assign("N.tryAgain", tryAgain , envir=MyEnv)   # set tryAgain==FALSE to exit the loop
                    dispose(mainNLSwin)
                    return()
                 }, container = NLSframe3)
            visible(mainNLSwin) <- TRUE
}
