## Initialize Elements

ElementCheck <- function(element) {  # chiamato da XPSsetRSF
	PeriodicTableElements <- c("H", "He", "Li", "Be", "B", "C", "N", "O", "F",
         "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc",
         "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As",
         "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh",
         "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La",
         "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm",
         "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl",
         "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th")

	return(ifelse(element %in% PeriodicTableElements, TRUE, FALSE))
}


## Initialize Elements for Kratos or Scienta dataFiles

IniElement <- function(element) {

# -- search for the element table

   CL.pthName <- "/RxpsG/data/CoreLinesTable.lib"
   LibPth <- .libPaths()  #path of all the R libraries where RxpsG could be located
   LL <- length(LibPth)
   CL.pthName <- paste(LibPth, CL.pthName, sep="") #paste acts on all the components of the LibPth vector
   OS<-Sys.info() #get system information
   if (OS["sysname"] != "Linux") {  #Windows and Mac Os Systems
      CL.pthName <- gsub("/","\\", CL.pthName, fixed=TRUE)   #path/filename for linux,  path\\filename for windows
   }
   for(ii in 1:LL){
      fe <- file.exists(CL.pthName[ii])
      if (fe == TRUE) {
         fp <- CL.pthName[ii]  #load the path and filename in the file_pointer fp
         break                 #file found break loop
      }
   }
   if (fe == FALSE) {
       gmessage(msg="ATTENTION: CoreLinesTable file not Found. Check RxpsG package", title = "WARNING",icon = "warning" )
       return()
   }
# -- now read element table
	ElmtList <- scan(fp, skip=0, sep="", what = list("character", "character", "numeric", "numeric", "numeric", "numeric"), quiet=TRUE)
	ElmtList[[3]] <- as.numeric(ElmtList[[3]])
	ElmtList[[4]] <- as.numeric(ElmtList[[4]])
	ElmtList[[5]] <- as.numeric(ElmtList[[5]])
	ElmtList[[6]] <- as.numeric(ElmtList[[6]])

	ElmtList <- as.data.frame(ElmtList, stringsAsFactors = FALSE) # set it to FALSE
   names(ElmtList) <- c("Element", "Orbital", "BE", "KE", "RSF_K", "RSF_S")
   
	## search for element in the Element Table
	idx <- which(ElmtList[,"Element"] == element)
   return(ElmtList[idx,]) #Observe the comma to return all elements of ElmtList at rows==idx
}

#'Elements functions utility
#'
#'Unified functions to access and get values from the references files supplied
#'for Scienta and Kratos instruments. These function are internal functions and
#'not intended for the user.
#'
#'Unified functions to acces and get values from the references files supplied
#'for Scienta and Kratos instruments.
#'
#'@param element element symbol as character
#'@param orbital atom orbital
#'@param analyzer name of reference instrument dependent table of values 
#'@param what name of value to get. Default \code{"RSF"}
#'@return Return values named by \code{what}.
#'
#'@export
#'

getElementValue <- function(element, orbital, analyzer=c("scienta", "kratos"), what="RSF") {  # chiamato da XPSsetRSF

   tmp <- IniElement(element) #extract the desired element from the CoreLine table
	if ( length(tmp) == 0){
      txt <- paste("Warning: element ",element," or relative orbital not found!")
	   gmessage(msg=txt, title="WARNING", icon="warning")
		return(NULL)
   }
	idx <- grep(orbital, tmp[,"Orbital"]) #select the desired orbital
	if ( length(idx) == 0 ) {
	   txt <- paste("Orbital ", orbital, " not found in element ", element, sep="")
		gmessage(msg=txt, title="WARNING", icon="warning")
		return(NULL)
	} else {
      tmp <- tmp[idx, ]             #select the desired orbital among the other of the chosen element
	   if (analyzer == "scienta") {  #in the CoreLine Table select the line corresponding to Scienta
	       idx<-which(is.na(tmp[,"RSF_S"])==TRUE)
          if (length(idx) > 0 ){tmp <- tmp[-idx ,]}  #drop rows with NotAssigned Scienta RSF (in this line Kratos values are set)
	       tmp <- tmp[,-5]           #drop the column with Kratos RSF_K
	   } else if (analyzer == "kratos") {   #in the CoreLine Table select the line corresponding to Kratos
	       idx<-which(is.na(tmp[,"RSF_K"])==TRUE)
          if (length(idx) > 0 ){tmp <- tmp[-idx ,]}  #drop rows with NotAssigned Kratos RSF  in this line Scienta values are set)
		    tmp <- tmp[,-6]           #drop the column with Scienta RSF_S
		}
   }
   names(tmp) <- c("Element", "Orbital", "BE", "KE", "RSF")
	what <- match.arg(what, colnames(tmp))
	return(tmp[,what])
}

#'Elements functions utility
#'
#'Unified functions to acces and get values from the references files supplied
#'for Scienta and Kratos instruments. These function are internal functions and
#'not intended for the user.
#'
#'Unified functions to acces and get values from the references files supplied
#'for Scienta and Kratos instruments.
#'
#'@param element element symbol as character
#'@param analyzer name of reference list
#'@return Return values named by \code{what}.
#'
#'@export
#'

showTableElement <- function(element, analyzer=c("scienta", "kratos")) {   # chiamato da XPSsetRSF
   tmp <- IniElement(element)           #extract the desired element from the CoreLine table
	analyzer <- match.arg(analyzer)
   switch(analyzer,                     #in the CoreLine analyzer select the lline corresponding to Scienta or Kratos
	       "scienta" = tmp <- tmp[,-5],  #drop the column with Kratos RSF
		    "kratos"  = tmp <- tmp[,-6]   #drop the column with Scienta RSF
	)
   names(tmp) <- c("Element", "Orbital", "BE", "KE", "RSF")
	cat("\n Element:  ", tmp[1,"Element"])
	orbit <- sapply(tmp[,"Orbital"], function(x) sprintf("%-7s",x), USE.NAMES = FALSE)
	cat("\n Orbitals: ", orbit)
	BE <- sapply(tmp[,"BE"], function(x) sprintf("%-7s",x))
	cat("\n BE:       ", BE)
	KE <- sapply(tmp[,"KE"], function(x) sprintf("%-7s",x))
	cat("\n KE:       ", KE)
	rsf <- sapply(tmp[,"RSF"], function(x) sprintf("%-7s",x))
	cat("\n RSF:      ", rsf)
   table<-list()
   table$Element<-tmp[,"Element"]
   table$Orbital<-orbit
   table$BE<-BE
   table$KE<-KE
   table$RSF<-rsf
   return(table)
}
