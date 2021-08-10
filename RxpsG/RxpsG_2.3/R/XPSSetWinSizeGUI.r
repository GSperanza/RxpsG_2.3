#function to set the dimensions of the analysis window (XPSGUI.r)

#'To select a dimensions of the graphic window depending on the dimensions of the screen used
#'
#'To select a dimensions of the graphic window depending on the dimensions of the screen used
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


XPSSetWinSize <- function() {


   winsize<-NULL
   win <- gwindow("ANALYSIS WINDOW SIZE", visible=FALSE)
   size(win) <- c(300, 180)
   group1 <- ggroup(label="SELECT WINDOW DIMENSIONS", horizontal=FALSE, spacing = 5, container=win)
   frame1 <-gframe(text="  WINDOW DIMENSIONS  ", spacing=10, container=group1)
   obj1 <- gslider(from=1, to=3, by=0.2, value=1.6, expand=TRUE, handler=function(h,...){
                    winsize <<- as.numeric(svalue(obj1))
                    txt <- paste("Window Size  = ", winsize, sep="")
                    svalue(WSize) <- txt
                    }, container=frame1)

   WSize <- glabel("Window Size  = 1.6", spacing=10, container = group1)

   glabel("           ", spacing=10, container = group1) #just to add space
   gbutton("SET SIZE and EXIT", handler=function(h,...){
                    winsize
                    XPSSettings<-get("XPSSettings", envir=.GlobalEnv)   #set the new WinSize value in the XPSSettigs
                    XPSSettings$General[4] <- winsize
                    assign("XPSSettings", XPSSettings, envir=.GlobalEnv)
                    assign("winsize", winsize, envir=.GlobalEnv)
                    dispose(win)
                 }, container = group1)


   visible(win) <- TRUE
}
