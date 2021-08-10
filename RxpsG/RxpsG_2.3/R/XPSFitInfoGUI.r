#function to provide the list of fitting components and relative abundance for the selected coreline

#'Provides information about a Core Line fit
#'
#'Makes the computation of the integral intensity of each of the fitting components
#'of a given coreline. The sum of the component integral intensities =100% (it is the
#'best fit integral intensity). No parameters are passed to this function.
#'
#'@examples
#'
#'\dontrun{
#'	XPSFitInfo()
#'}
#'
#'@export
#'


XPSFitInfo <- function(){

   info <- function(CoreLine){

      CoreLineComp <- names(CoreLine@Components)
      maxNchar <- 0   
      sumCoreLine <- 0
      N_comp=length(CoreLine@Components) #number of fit components
      sumComp <- array(0,dim=N_comp)       #define a dummy vector of zeros
      RSF  <-  FName[[Indx]]@RSF 
      E_stp  <-  abs(FName[[Indx]]@.Data[[1]][2]-FName[[Indx]]@.Data[[1]][1]) #energy step

      if (N_comp==0) {
#--- solo Baseline No Fit ---
         sumCoreLine <- sum(FName[[Indx]]@RegionToFit$y-FName[[Indx]]@Baseline$y)*E_stp #Area of the Coreline NOT normalized for the sensitivity factor
         txt <- as.character(sumCoreLine)
         maxNchar <- nchar(txt)
         if (length(FName[[Indx]]@Baseline)>0 ){ #Baseline is defined
            txt <- paste("\n   File Name:", FName@Filename)  #Filename
            CellBord <- ""
            cell <- printCell("label",txt,CellBord,nchar(txt),"left")  #call celprint in LABEL mode
            area <- sprintf("%1.2f",sumCoreLine)  #print area Coreline in formatted way (2 decimals)
            txt <- paste(FName[[Indx]]@Symbol,"peak area: ", area, sep=" ")
            totLgth <- nchar(txt)
            cell <- printCell("separator", "-",CellBord,totLgth-5,"left")  #call celprint in SEPARATOR mode
            cat("\n",cell)
            CellBord <- "|"
            cell <- printCell("label",txt,CellBord,totLgth,"center")
            cat("\n",cell)
            return()
         }

      } else {
#--- Quantificazione sul FIT ---
         sumCoreLine <- sum(FName[[Indx]]@Fit$y)*E_stp  #Fit Contribution not corrected fro the sensitivity factor
         txt <- as.character(sumCoreLine)
         maxNchar <- nchar(txt)
         for(jj in 1:N_comp){    #jj runs on the fit components
            sumComp[jj] <- sum(FName[[Indx]]@Components[[jj]]@ycoor-FName[[Indx]]@Baseline$y)*E_stp #Contributo Componente Fit as acquired non corretto per RSF
         }

         if (11 > maxNchar) { maxNchar <- 11 }  #maxNchar = number among area values is computed
                                            #if maxChar > 11 then 11 is selected.
         lgth <- c(12, maxNchar, 8, 8, 7, 9)  #number of char reserved for "Components", "Area", ", FWHM", "BE(eV)", "RSF", "TOT%"
         totLgth <- sum(lgth)
         cat("\n")
         CellBord <- "|"
         txt <- paste("   File Name:", FName@Filename)  #Filename
         cell <- printCell("label",txt,CellBord, totLgth,"left")
         cat("\n",cell)

         txt <- "-"
         cell <- printCell("separator",txt,CellBord,totLgth,"left")
         cat("\n",cell)

         if (FName[[Indx]]@Flags[1]==TRUE){
            txt <- c("Components", "Area(cps)", "FWHM", "RSF", "BE(eV)", "TOT.(%)")    #Table names
         } else {
            txt <- c("Components", "Area(cps)", "FWHM", "RSF", "KE(eV)", "TOT.(%)")    #Voci Tabella
         }
         lgth <- c(12, maxNchar, 8, 7, 8, 9)
         cell <- printCell("tableRow",txt,CellBord,lgth, "center")  #call celprint in tableRow mode
         cat("\n",cell)

         txt <- "-"  #separatore
         cell <- printCell("separator", txt,CellBord,totLgth,left)
         cat("\n",cell)

         ClName <- FName[[Indx]]@Symbol
         Area <- sprintf("%1.2f", sumCoreLine)
         RSF <- sprintf("%1.3f",FName[[Indx]]@RSF)
         Conc <- sprintf("%1.2f",100)
         txt <- c(ClName, Area, " ", RSF, " ", Conc )  #total concentration associated to the CoreLine. FWHM e BE are not printed
         lgth <- c(12, maxNchar, 8, 7, 8, 9)
         cell <- printCell("tableRow", txt,CellBord,lgth, "center")
         cat("\n",cell)

         for(jj in 1:N_comp){ #jj runs on the fit components
            Area <- sprintf("%1.2f",sumComp[jj]) #area componente jj linea di core ii
            FWHM <- sprintf("%1.2f",FName[[Indx]]@Components[[jj]]@param["sigma", "start"]) #FWHM componente ii
            RSF <- sprintf("%1.3f",FName[[Indx]]@Components[[jj]]@rsf) #RSF componente ii
            BE <- sprintf("%1.2f",FName[[Indx]]@Components[[jj]]@param["mu","start"]) #BE componente ii
            Conc <- sprintf("%1.2f",100*sumComp[jj]/sumCoreLine)  #Concentrazione componente ii
            txt <- c(CoreLineComp[jj], Area, FWHM, RSF, BE, Conc) #compongo la stringa da stampare
            lgth <- c(12, maxNchar, 8, 7, 8, 9) #lunghezza varie celle riga tabella
            cell <- printCell("tableRow",txt,CellBord,lgth, "center") #stampo nella modalita' TABLEROW
            cat("\n",cell)
         }

         txt <- " "  #separatore
         cell <- printCell("separator",txt,CellBord,totLgth,left)
         cat("\n",cell)
      }
      return()
   }



#----- variabili -----

   FName <- NULL
   ActiveFName <- ""
   FName <- get(activeFName, envir = .GlobalEnv)   #carico l'XPSSample dataFrame attivo
   ActiveFName <- get("activeFName", envir = .GlobalEnv)  #carico il nome XPSSample (stringa)
   Indx <- get("activeSpectIndx", envir = .GlobalEnv)  #carico il nome XPSSample (stringa)
   activeSpectName <- get("activeSpectName", envir = .GlobalEnv)  #carico il nome XPSSample (stringa)
   if (length(FName[[Indx]]@RegionToFit) == 0){  #no information se il Baseline non presente
      message <- paste("No Baseline/Fit found for", ActiveFName, sep=" ")
      return()
   }
   message <- paste(activeSpectName, "Coreline Fit Info: ", sep=" ")
   cat("\n", message)
   info(FName[[Indx]])
}
