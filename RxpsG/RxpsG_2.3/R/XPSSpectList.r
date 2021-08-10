# crea la lista degli spettri associati ad un dato Spectral identifier

#'Show the list of spctra related to an XPSSample identified by its Spectral ID
#'
#'@param activeFName the name (string) of the selected XPSSample
#'@param noIdx logical if TRUE a numerical index is added to any CoreLine name
#'
#'@examples
#'
#'\dontrun{
#'	XPSSampleInfo(activeFName, noIdx=FALSE)
#'}
#'
#'@export
#'


XPSSpectList <- function(activeFName, noIdx=FALSE) {
     if (is.na(activeFName)){
        gmessage("No data present: please load and XPS Sample", title="XPS SAMPLES MISSING", icon="error")
        return()
     }
     FName<-get(activeFName,envir=.GlobalEnv)
     CoreLineList=""
     CoreLineList<-FName@names
     LL<-length(CoreLineList)
     if (noIdx==FALSE){
        for (ii in 1:LL){
               CoreLineList[ii]<-paste(ii,".",CoreLineList[ii], sep="")   #Aggiungo N. all'inizio nome coreline per distinguere spettri con eguale nome
        }
     }
     return(CoreLineList)
}


