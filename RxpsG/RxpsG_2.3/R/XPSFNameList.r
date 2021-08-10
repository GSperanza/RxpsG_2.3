#crea la lista degli XPS-Sample Names caricati nel .GlobalEnv

#'Provide a list of the XPS-Sample names loaded
#'
#'Provide a list of names (strings) relative to the XPS-Samples
#'loaded in the XPS-Analysis software. No parameters are passed to this function.
#'
#'@param warn logical enables warning messages
#'
#'@examples
#'
#'\dontrun{
#'	XPSFNameList()
#'}
#'
#'@export
#'


XPSFNameList <- function(warn=TRUE){
   FNameList=NULL
   ClassFilter <- function(x) inherits(get(x), "XPSSample" ) #sets the "XPSSample" as the class to filter to select XPSSample Data files
   FNameList <- Filter(ClassFilter, ls(.GlobalEnv)) #list of XPSSample files loaded in .GlobalEnv
   if (length(FNameList)==0 && warn==TRUE){
     ErMsg <- gmessage(msg = "Cannot find XPS-Samples! Please load XPS data files", title="XPS SAMPLE LACKING", icon = "warning", type = "ok")
     focus(ErMsg) <- TRUE
     FNameList<-NULL
   }
   return(FNameList)
}
