#Setting the Graphic device for different Operating Systems

#'To select the kind of graphic device compatible with the operating system in use
#'
#'To select the kind of graphic device compatible with the operating system in use
#'Also a list of graphic formats is provided to save the content of the current
#'No parameters are passed to this function
#'
#'
#'@examples
#'
#'\dontrun{
#'	XPSSetGraphDev()
#'}
#'
#'@export
#'


XPSSetGraphDev <- function() {

   ChDir <- function(){
          workingDir <- getwd()
          PathName <- tk_choose.dir(default=workingDir)
          PathName <- paste(dirname(PathName), "/", basename(PathName), sep="") #changes from \\ to /
          return(PathName)
   }

   CutPathName <- function(PathName){  #taglia il pathname ad una lunghezza determinata
             splitPathName <- strsplit(PathName, "/")
             LL <- length(splitPathName[[1]])    
             headPathName <- paste(splitPathName[[1]][1],"/", splitPathName[[1]][2], "/ ... ", sep="") 
             ShortPathName <- paste(headPathName, splitPathName[[1]][LL], sep="")
             ii=0
             tailPathName <- ""
             while (nchar(ShortPathName) < 25) {
                   tailPathName <- paste("/",splitPathName[[1]][LL-ii],tailPathName, sep="")
                   ShortPathName <- paste(headPathName, tailPathName, sep="") 
                   ii <- ii+1
             }       
             return(ShortPathName)
   }



#--- variables
   OSList <- c("Windows", "MacOS-X", "Linux")
   FormatList <- c("png", "jpeg", "bmp","tiff","eps", "pdf")
   pathName <- getwd()
   XPSSettings <- get("XPSSettings", envir=.GlobalEnv)
   if (length(activeFName) == 0) activeFName <- " "   #no XPSSample are loaded

# --- Widget ---

   GDwin <- gwindow("GRAPHIC DEVICE ", parent=c(0, 200), visible=FALSE)
   GDgroup <- ggroup(label="GRAPHIC DEVICE OPTIONS", horizontal=FALSE, container=GDwin)

   frame1 <- gframe(" SELECT your Operating System ", spacing=5, container=GDgroup)

   obj1 <- gradio(OSList, selected=-1, horizontal=FALSE, handler=function(h,...){
                      OS <- svalue(obj1)
                      switch (OS,
                          "Linux" =   {Gdev <- "X11(xpos=600, ypos=5, title=' ')" },
                          "Windows" = {Gdev <- "X11(xpos=600, ypos=5, title=' ')"},
                          "MacOS-X" = {Gdev <- "quartz(title=' ')" },  #quartz() does allow setting the opening position
                          "Mac OS"  = {Gdev <- "quartz(title=' ')" },
                          "macOS"   = {Gdev <- "quartz(title=' ')" },
                          "Darwin"  = {Gdev <- "quartz(title=' ')" })
                      XPSSettings$General[6]<<- Gdev
                      graphics.off()
                      eval(parse(text=Gdev),envir=.GlobalEnv)
           }, container = frame1) 

   frame2 <- gframe("FILE FORMAT TO SAVE", horizontal=FALSE, spacing=5, container=GDgroup)

   obj2 <- gradio(FormatList, selected=-1, horizontal=FALSE, container=frame2)
   glabel(text="File Name (No extension):", container=frame2)
   obj3 <- gedit(text="", container=frame2)

   txt <- paste("Current dir: ", pathName, sep="")
   obj4 <- glabel(text=txt, container=frame2)

   gbutton("CHANGE DIRECTORY", handler=function(h,...){
                      pathName<<- ChDir()
                      pathName <- CutPathName(pathName)
                      txt <- paste("Current dir: ", pathName, sep="")
                      svalue(obj4)<<- txt
                   },container=frame2)

   gbutton("  EXPORT TO FILE   ", handler=function(h,...){
                      Format <- svalue(obj2)
                      PlotFileName <- svalue(obj3)
                      PlotFileName <- paste(pathName,"/",PlotFileName,".",Format,sep="")
                      if (Format == "png") dev.print(file=PlotFileName, device=png, bg="white", width=1024)
                      if (Format == "jpeg") dev.print(file=PlotFileName, device=jpeg,  bg="white", width=1024)
                      if (Format == "bmp") dev.print(file=PlotFileName, device=bmp,  bg="white",width=1024)
                      if (Format == "tiff") dev.print(file=PlotFileName, device=tiff, bg="white", width=1024)
                      if (Format == "eps") dev.print(file=PlotFileName, device=postscript, horizontal=FALSE, pointsize=1)
                      if (Format == "pdf") dev.print(file=PlotFileName,device=pdf)
                      cat("\n Graphic Device Exported to .", Format, " File")
                      Gdev <- dev.cur()
          }, container=frame2)

   gbutton(" SAVE&EXIT ", handler=function(h,...){
   #--- get System info and apply correspondent XPS Settings ---
                      OS <- Sys.info() #get system information
#Reading XPS settings which can be customized by users
                      Ini.pthName <- "/RxpsG/data/XPSSettings.ini"
                      LibPth <- .libPaths()  #path of all the R libraries where RxpsG could be located
                      LL <- length(LibPth)
                      Ini.pthName <- paste(LibPth, Ini.pthName, sep="") #paste acts on all the components of the LibPth vector
                      if (OS["sysname"] != "Linux") { #Windows and Mac OS systems
                         Ini.pthName <- gsub("/","\\", Ini.pthName, fixed=TRUE)   #path/filename for linux,  path\\filename for windows
                      }
                      for(ii in 1:LL){
                         fe <- file.exists(Ini.pthName[ii])
                         if (fe == TRUE) {
                            ColNames <- names(XPSSettings)
                            write.table(XPSSettings, file = Ini.pthName[ii], sep=" ", eol="\n", row.names=FALSE, col.names=ColNames)
                            break
                         }
                      }
                      if (fe == FALSE) {
                         gmessage(msg="ATTENTION: XPSSettings.ini file is lacking. Check RxpsG package", title = "WARNING",icon = "warning" )
                         return()
                      }
                      assign("XPSSettings", XPSSettings, envir=.GlobalEnv)
                      dispose(GDwin)
                 }, container = frame2)

   visible(GDwin) <- TRUE
}
