### $Id: readSPC.R,v 0.0.1
###
###             Read XPS spectra data (pxt or vamas)
###

#'Read '.pxt' and '.vms' format Data
#'
#'Reads a file and returns an XPSSpectra class item. In the case the file
#'extension is \code{.pxt} and the option \code{Genplot} is \code{TRUE}, then a
#'Genplot \code{RPL} file with the same name will be searched in a folder
#'\code{RPL} and it will be sourced if founded.
#'
#'@note \code{read.scienta} \code{read.vamas} and \code{readGenplot} 
#'must not used by the end-user
#'
#'@param file file name
#'@param Genplot flag to read genplot .RPL files for 'Scienta' format files
#'@param \dots further parameters to \code{read.scienta} and \code{read.vamas}
#'@return returns an object of class XPSSample
#'
#'@examples
#'
#'\dontrun{
#'SampData1 <- XPSread("SampData1.pxt")
#'}
#'
#'\dontrun{
#'SampData2 <- XPSread("SampData2.vms")
#'}
#'
#'@export
#'

XPSread <- function( file=NULL, Genplot=TRUE, ... )

{
	if ( is.null(file) ) file <- file.choose()

# check extension
	  FName <- basename(file)
	  dirName <- dirname(file)
   setwd(dirName)
   WorkingDir<-getwd()

#check filename extension to call relative reading function
	  fext <- unlist(strsplit(FName, "\\."))[2]
	  filewext <- unlist(strsplit(FName, "\\."))[1]

	  f.idx <- grep(fext, c("pxt", "PXT", "vms", "VMS"))
   if (length(f.idx)==0) { #any extension different from pxt or vms => old scienta file
	      object <- XPSRead.Oldscienta(file, ...)
   } else if ( f.idx == 1 || f.idx == 2 ) {  # pxt or PXT  extension
	     	object <- read.scienta(file, ...)
	     	if (Genplot) {       #Genplot==TRUE => .pxt+RPL option chosen. RPL folder contains old Genplot analysis information
	        		rplFile <- paste(filewext,".RPL", sep="")
        			checkfile <- list.files(dirName, pattern=rplFile, ignore.case=TRUE, recursive=TRUE, full.names=TRUE)
        			if (length(checkfile)) {
		           		object <- readGenplot(object, file=checkfile, ...)
			       	}
		     }
	  } else if ( f.idx == 3 || f.idx == 4) {  # vms or VMS extension
		   object <- read.vamas(file, ...)
	  }

	  slot(object, "Filename") <- FName
	  return(object)
}
