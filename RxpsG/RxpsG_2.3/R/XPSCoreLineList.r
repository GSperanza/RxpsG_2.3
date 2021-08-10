#coreline list

#'Make the list of CoreLines of the loaded XPS-Sample
#'
#'The list of corelines which are contained in the loaded XPS-Sample
#'are displayed by this function. No parameters are passed to this function.
#'
#'@examples
#'
#'\dontrun{
#'	XPSCoreLineList()
#'}
#'
#'@export
#'

XPSCoreLineList <- function() {
  #load the list of file IDs
     ClassFilter <- function(x) inherits(get(x), "XPSSample" ) #set class==XPSSample for filtering
     IDList <- Filter(ClassFilter, ls(.GlobalEnv)) #filters the list of variables in the .GlobalEnv selecting those of class XPSSample
     if (length(IDList)==0){
       tkmessageBox(message = "NO SPECTRA LOADED", icon = "warning", type = "ok")
     } 
     LID=length(IDList)
     Spectra=""
     List=""
  #For the first XPSSample get the list of corelines: eval used to operate on XPSSamples in .GlobalEnv
     command<-paste("names(",IDList[1],")", sep="")
     Spectra<-eval(parse(text=command),envir =.GlobalEnv)  #Save the list of corelines in the return vector
     List<-Spectra
     LL=length(List)
  #Control if there is more than one XPSSample compare coreline_names
     if (LID>1){
        for(ii in 2:LID) {
          command<-paste("names(",IDList[ii],")", sep="")
          Spectra<-eval(parse(text=command),envir = .GlobalEnv)
          LS=length(Spectra)
          for (jj in 1:LS){
  #add core line names to the list of spectra if they are not present
            if(length(grep(Spectra[jj], List))==0) {
               List[LL+1]<-Spectra[jj]
               LL<-LL+1
            }
          }
        }
     }
     return(List)
}