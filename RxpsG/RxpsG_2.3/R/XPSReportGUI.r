#function to provide the list of fitting components and relative abundance for the selected coreline

#'Provides information about a Core Line fit
#'
#'Makes the computation of the integral intensity of each of the fitting components
#'of a given coreline. The sum of the component integral intensities =100% (it is the
#'best fit integral intensity). No parameters are passed to this function. 
#'
#'
#'@examples
#'
#'\dontrun{
#'	XPSReport()
#'}
#'
#'@export
#'


XPSReport <- function(){


   doReport <- function(CoreLine){

      CoreLineComp<-names(CoreLine@Components)
      maxNchar<-0
      sumCoreLine<-0
      N_comp=length(CoreLine@Components) #e' il numero di componenti del fit
      sumComp<-array(0,dim=N_comp)  #definisco un vettore di zeri
      RSF <- FName[[Indx]]@RSF
      E_stp <- abs(FName[[Indx]]@.Data[[1]][2]-FName[[Indx]]@.Data[[1]][1]) #energy step

      if (length(CoreLine@Baseline)==0) {
#--- No Baseline No Fit ---
         txt<-paste("\n*** ", CoreLine@Symbol, ": no fit present", sep=" ")
         CellLength<-nchar(txt)
         TabTxt<-c(TabTxt,printCell("label",txt,cellB,CellLength,"left"))  #chiamo cellprint in modalita' LABEL
         return(TabTxt)

      } else if (length(CoreLine@Baseline) > 0 && N_comp==0) {
#--- solo Baseline No Fit ---
         txt<-paste("\n*** ", CoreLine@Symbol, "Coreline only Baseline present: ", sep=" ")
         CellLength<-nchar(txt)
         TabTxt<-c(TabTxt,printCell("label",txt,cellB,CellLength,"left"))  #chiamo cellprint in modalita' LABEL

         if (RSF==0){
            sumCoreLine<-sum(FName[[Indx]]@RegionToFit$y-FName[[Indx]]@Baseline$y*E_stp) #Contributo della linea di core integrale RSF==1
         } else {
            sumCoreLine<-sum(FName[[Indx]]@RegionToFit$y-FName[[Indx]]@Baseline$y)*E_stp/RSF #Contributo della linea di core integrale
         }
         area<-sprintf("%1.2f",sumCoreLine)  #converto
         txt<-paste("  ", FName[[Indx]]@Symbol,"peak area: ", area, sep=" ")
         CellLength<-nchar(txt)
         TabTxt<-c(TabTxt,printCell("label",txt,cellB,CellLength,"left"))  #chiamo cellprint in modalita' LABEL
         return(TabTxt)

      } else if (N_comp > 0) {
#--- Baseline + Fit ---
#--- Quantificazione sul FIT ---
         sumCoreLine<-0
         for(jj in 1:N_comp){    #jj corre sulle componenti del fit della coreline
            RSF<-FName[[Indx]]@Components[[jj]]@rsf
            if (RSF==0) { #non e' stato definito un RSF(es. Auger, VB spectra...) : non posso dividere
               sumComp[jj]<-sum(FName[[Indx]]@Components[[jj]]@ycoor-FName[[Indx]]@Baseline$y)*E_stp
            } else {
               sumComp[jj]<-sum(FName[[Indx]]@Components[[jj]]@ycoor-FName[[Indx]]@Baseline$y)*E_stp/RSF  #contributo della singola componente
            }
            sumCoreLine <- sumCoreLine + sumComp[jj]
         }
#                  sumCoreLine<-sum(FName[[Indx]]@Fit$y)/RSF #Contributo del Fit

# larghezza colonne("Components", "FitFunct.", "Area(cps)", "Intensity", "FWHM", "BE(eV)", "TOT.(%)")
         CellLength<-c(12, 15, 15, 13, 8, 8, 9)
         cellB<-" "
         txt<-paste("\n*** ", CoreLine@Symbol, "Coreline Fit Info: ", sep=" ")
         lgth<-nchar(txt)
         TabTxt<-c(TabTxt,printCell("label",txt,cellB,lgth,"left"))  #chiamo cellprint in modalita' LABEL
#Nomi Colonne
#         TabTxt<-c(TabTxt,"\n")
         txt<-c("Components", "Fit Funct.", "Area (cps)", "Intensity", "FWHM", "Position", "TOT.(%)")
         TabTxt<-c(TabTxt,printCell("tableRow", txt, cellB, CellLength, "center"))
#PrimaRiga

         ClName<-FName[[Indx]]@Symbol
         Area<-sprintf("%1.2f", sumCoreLine)  #tronco a 2 decimali e trasformo in stringa
         Conc<-sprintf("%1.2f",100)
         txt<-c(ClName, " ", Area, " ", " ", " ", Conc )  #concentrazioni totali della coreline: FitFunct, FWHM e BE non vengono stampate
         TabTxt<-c(TabTxt,printCell("tableRow", txt, cellB, CellLength, "center"))
#Componenti
         for(jj in 1:N_comp){ #ii corre sulle CoreLines, jj corre sulle componenti del fit della coreline
            Function<-sprintf("%s",FName[[Indx]]@Components[[jj]]@funcName) # Fit Funct name
            Area<-sprintf("%1.2f",sumComp[jj]) #area componente jj linea di core ii
            Intensity<-sprintf("%1.2f",FName[[Indx]]@Components[[jj]]@param[1,1]) # Intensita' componente
            FWHM<-sprintf("%1.2f",FName[[Indx]]@Components[[jj]]@param[3,1]) #FWHM componente ii
            BE<-sprintf("%1.2f",FName[[Indx]]@Components[[jj]]@param[2,1]) #BE componente ii
            Conc<-sprintf("%1.2f",100*sumComp[jj]/sumCoreLine)  #Concentrazione componente ii
            txt<-c(CoreLineComp[jj], Function, Area, Intensity, FWHM, BE, Conc) #compongo la stringa da stampare
            TabTxt<-c(TabTxt,printCell("tableRow", txt, cellB, CellLength, "center")) #stampo nella modalita' TABLEROW
         }
#Bordo inferiore
#         TabTxt<-c(TabTxt,"\n")
         return(TabTxt)
      }
   }



#----- variabili -----

   FName<-NULL
   ActiveFName<-""
   FName<-get(activeFName, envir = .GlobalEnv)   #carico l'XPSSample dataFrame attivo
   ActiveFName<-get("activeFName", envir = .GlobalEnv)  #carico il nome XPSSample (stringa)
   Indx<-get("activeSpectIndx", envir = .GlobalEnv)  #carico il nome XPSSample (stringa)
   activeSpectName<-get("activeSpectName", envir = .GlobalEnv)  #carico il nome XPSSample (stringa)
   TabTxt <- ""
   get("XPSSettings", envir=.GlobalEnv)
   Font <- XPSSettings$General[1]
   FStyle <- XPSSettings$General[2]
   FSize <- XPSSettings$General[3]

   NCorelines <- length(FName)
   txt<-paste("\n ==> File Name:", FName@Filename)  #Filename
   CellLength<-nchar(txt)
   TabTxt<-c(TabTxt,printCell("label",txt,cellB,CellLength,"left"))  #chiamo cellprint in modalita' LABEL

   for (Indx in 1:NCorelines){
      TabTxt<-doReport(FName[[Indx]])
      TabTxt<-c(TabTxt,sprintf("%s", "\n"))
   }

#---Widget
      FontAttr <- list(weight="normal", style="normal", family="normal", size="small")
      txtWin <- gwindow(title="XPS Sample Report", visible=FALSE)
      size(txtWin) <- c(750, 500)
      tabWin <- gtext(text=TabTxt, font.attr=FontAttr, container = txtWin)
      font(tabWin) <- list(weight="light", family=Font, style=FStyle, size=FSize)
      visible(txtWin)<-TRUE
}
