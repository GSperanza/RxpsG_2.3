## ===============================
## Identification of survey peaks
## ===============================

#--- function to read the Coreline table and identify elements
#--- A CoreLineTable is read containing a list of elements OK for Scienta and Kratos spectra
ReadElmtList<-function(SpectrumType){

   CL.pthName <- "/RxpsG/data/CoreLinesTable.lib"
   Aug.pthName <- "/RxpsG/data/AugerTransitionsTable.lib"
   LibPth <- .libPaths()  #path of all the R libraries where RxpsG could be located
   LL <- length(LibPth)
   CL.pthName <- paste(LibPth, CL.pthName, sep="") #paste acts on all the components of the LibPth vector
   Aug.pthName <- paste(LibPth, Aug.pthName, sep="") #paste acts on all the components of the LibPth vector
   OS<-Sys.info() #get system information
   if (OS["sysname"] != "Linux") {  # Windows and Mac OS systems
      CL.pthName <- gsub("/","\\", CL.pthName, fixed=TRUE)   #path/filename for linux,  path\\filename for windows
      Aug.pthName <- gsub("/","\\", Aug.pthName, fixed=TRUE) #also gsub acts on all the components of  .pthName
   }
   for(ii in 1:LL){
      if (SpectrumType == "CoreLines") {
         fe <- file.exists(CL.pthName[ii])
         if (fe == TRUE) {
            fp <- CL.pthName[ii]  #load the path and filename in the file_pointer fp
            break
         }
      } else if (SpectrumType == "AugerTransitions") {
         fe <- file.exists(Aug.pthName[ii])
         if (fe == TRUE) {
            fp <- Aug.pthName[ii] #load the path and filename in the file_pointer fp
            break
         }
      }
   }
   if (fe == FALSE) {
       gmessage(msg="ATTENTION: CoreLinesTable OR AugerTransitionsTable file not Found. Check RxpsG package", title = "WARNING",icon = "warning" )
       return()
   }

# -- now read the Coreline or Auger table
  	pp <- scan(fp, skip=0, sep="", what = list("character", "character", "numeric", "numeric", "numeric", "numeric"), quiet=TRUE)
   cat("\n ==>", SpectrumType, " Table loaded \n")
	  pp[[3]] <- as.numeric(pp[[3]])
	  pp[[4]] <- as.numeric(pp[[4]])
	  pp[[5]] <- as.numeric(pp[[5]])
	  pp[[6]] <- as.numeric(pp[[6]])

	  ElmtList <- as.data.frame(pp, stringsAsFactors = FALSE) # set it to FALSE
	  if (SpectrumType == "CoreLines") {
	    names(ElmtList) <- c("Element", "Orbital", "BE", "KE", "RSF_K", "RSF_S")
	  } else if (SpectrumType == "AugerTransitions") {
	    names(ElmtList) <- c("Element", "Transition", "BE", "KE", "RSF_K", "RSF_S")
   }
   return(ElmtList)
}


peakDetection <- function(object, snmin , Ewin) {
   import::here(baseline.peakDetection, .from=baseline)
 	 #- find peaks: DO NOT use baseline(object, method = "peakDetection")
 	 #- because it does not return the peaks index
 	 peaks <- baseline.peakDetection(matrix(data=object[[2]], nrow=1),
 									snminimum = snmin,
 									left.right = Ewin,
 									lwin.rwin = Ewin)
 	 #- subset peaks
 	 peaks <- peaks[c("baseline","corrected", "peaks")]
 	 idx <- unlist(peaks$peaks) ## peaks index
 	 #- coordinate dei picchi
 	 positionsCorr <- list(x=object$x[idx], y=peaks$corrected[1,idx])
 	 positionsOrig <- list(x=object$x[idx], y=object[[2]][idx])

 	 #- sort positions per intensity
   #- idx <- order(positionsCorr$y, decreasing=TRUE)

   #- peaks ordered from High BE to low BE
   idx <- c(1:length(positionsCorr$x))
  	peaks$table <- list(BE=positionsCorr$x[idx],
 		    					corr=positionsCorr$y[idx],
 						    	orig=positionsOrig$y[idx]
 						   )

	  names(peaks$table$BE) <- rep(NA, length(peaks$table$BE))
	  return(invisible(peaks))
}


#---Identify peaks in the survey spectrum
peakIdentification <- function(peaks, DeltaEnergy, object, SpectrumType, ElmtList) {
   RecPlot<-recordPlot()   #save graph for UNDO option

#-- init

   ElmtList1 <- ReadElmtList("CoreLines") #reads the CoreLine/Auger List Table and save it in MyEnv  (macro in XPSSurveyUtility)
   ElmtList1 <- format(ElmtList1, justify="centre", width=10)
   ElmtList1 <- as.data.frame(ElmtList1,  stringsAsFactors = FALSE)
   ElmtList2 <- ReadElmtList("AugerTransitions") #reads the CoreLine/Auger List Table and save it in MyEnv  (macro in XPSSurveyUtility)
   ElmtList2 <- format(ElmtList2, justify="centre", width=10)
   ElmtList2 <- as.data.frame(ElmtList2,  stringsAsFactors = FALSE)
   ElemLbl1 <- gsub(" ", "", ElmtList1[,1], fixed=TRUE) #suppress white spaces
   ElemLbl2 <- gsub(" ", "", ElmtList2[,1], fixed=TRUE) #suppress white spaces

	  peaks$table$BE <- unname(peaks$table$BE) # length of the list of identified peaks
	  ElmntTbl1 <- data.frame()
	  ElmntTbl2 <- data.frame()
	  FoundElmnt <- NULL
	  NElmnt <- 0
	  NOrbit <- NULL
	  BE<-FALSE  #set KE scale
 	 if (object@Flags[1]) {  # if object@Flags[1] == TRUE spectra are on a BE scale
	      BE <- TRUE
       Eidx=3  #BE column in the ElemenList tables
   } else {
       Eidx=4  #KE column in the ElemenList tables
   }
   peakIdx <- 1
#-- loop on peaks
   for (peakIdx in seq_along(peaks$table$BE)){
#-- if the peak is not assigned find for the element having higher RSF to that BE value
       NewElmnt<-NULL
       energy <- peaks$table$BE[peakIdx]
       NewElmnt <- FindElmtList(energy, DeltaEnergy, BE, ElmtList1) #list of elements matching the energy value
       if (length(NewElmnt)==0) {
           txt <- paste("No elements found at energy: ", energy)
           gmessage(msg=txt, title="NO ELEMENTS FOUND", icon="warning")
       } else {
           for (ii in seq_along(NewElmnt$Element)){ # Control if non void elements are present in the FoundElements list
               Elmnt <- NewElmnt$Element[ii]
               Elmnt <- gsub(" ", "", Elmnt, fixed=TRUE) #suppress white spaces
               idx <- grep(Elmnt, FoundElmnt) #check if the found element was already identified or discarded
               if (length(idx)==0) {  #the element found is new and not already present in list of FoundElements
                  if (length(Elmnt)==0){ #finds peak falling at energy position in the survey
                     idx<-findXIndx(object[[1]], energy)  #find index corresponding to x=energy
                     intensity<-max(object[[2]][(idx-5):(idx+5)])
                     newlab<-"??"  #display ?? for non identified spectral feature
                     DispLab(list(x=energy,y=intensity), label=newlab)
                     EE<-as.character(peaks$table$BE[peakIdx])
                     cat("\n Warning: No Element found at energy ", EE, "Skip value")
                     peakIdx <- peakIdx+1  #skip the peak
                  }
                  cat("\n Found Element: ", Elmnt)

                  replayPlot(RecPlot)  #restore original plot
                  idx1 <- which(Elmnt == ElemLbl1)    #greep does not work: if Elmnt=="O" extracts  rows corresponding to "O" and "Os"
                  Yrng<-range(object[[2]])
                  for (jj in seq_along(ElmtList1[idx1,Eidx])){
                     xx <- ElmtList1[idx1[jj],Eidx]
                     lines(x=c(xx, xx), y=Yrng, col="red") #plot CORELINES of the selected elements
                  }

                  idx2 <- which(Elmnt == ElemLbl2)    #greep does not work: if Elmnt=="O" extracts  rows corresponding to "O" and "Os"
                  for (jj in seq_along(ElmtList2[idx2,Eidx])){
                     xx <- ElmtList2[idx2[jj],Eidx]
                     lines(x=c(xx, xx), y=Yrng, col="blue") #plot AUGER TRANSITIONS of the selected elements
                  }

                  txt <- paste("Found correlation with ", Elmnt, " Core Lines. \nCorrect identification?", sep="")
                  answ=gconfirm(txt, title="FOUND ELEMENT", icon="warning")
                  if (answ) {
                     ElmntTbl1 <- rbind(ElmntTbl1, ElmtList1[idx1, ])
                     NElmnt <- NElmnt + 1
                     NOrbit[NElmnt] <-length(ElmtList1[idx1,2])  #number of orbitals of the identified element
                     ElmntTbl2 <- rbind(ElmntTbl2, ElmtList2[idx2, ])
                  }
                  FoundElmnt <- paste(FoundElmnt, Elmnt, sep=" ")
               }
           }
       }
    }  #end external for on detected peak

    replayPlot(RecPlot)  #restore original plot
    #now plot all the identified corelines and auger transition to show how many peaks are assigned
    for (jj in seq_along(ElmntTbl1$Element)){
        xx <- ElmntTbl1[jj,Eidx]
        lines(x=c(xx, xx), y=Yrng, col="red") #plot CORELINES of the selected elements
    }
    for (jj in seq_along(ElmntTbl2$Element)){
        xx <- ElmntTbl2[jj,Eidx]
        lines(x=c(xx, xx), y=Yrng, col="blue") #plot CORELINES of the selected elements
    }
    gmessage(msg=" Please look at identified peaks if they describes all the spectral features.\n Otherwise Refresh plot, change identification precision and repeat the analysis.", title="IDENTIFICATION RESULTS", icon="warning")

    if (length(ElmntTbl1$Element) == 0 ){
       return()
    } else {
       replayPlot(RecPlot)  #restore original plot
       Yrng<-range(object[[2]])
       DeltaY<-(Yrng[2]-Yrng[1])/25  #vertical shift to not superpose labels
       posLbl <- NULL
       arrow <- TRUE
       kk <- 1
       for(ii in 1:NElmnt){
          for(jj in 1:NOrbit[ii]){   #this for runs on the element orbitals which may be at same energy
             energy <- as.numeric(ElmntTbl1$BE[kk])  #kk runs on the ElmntTable rows
             idx <- findXIndex(object[[1]], energy)
             intensity <- object[[2]][idx]
             posLbl<-intensity+DeltaY*jj #jj increments for differen orbital same element
             if (SpectrumType =="CoreLines"){
                newlab <- paste(ElmntTbl1$Element[kk], ElmntTbl1$Orbital[kk], sep="")
                newlab <- gsub(" ", "", newlab, fixed=TRUE)  #drop white spaces
                DispLab(list(x=energy,y=posLbl), label=newlab, arrow, direction="right", color="blue")
             } else if (SpectrumType == "AugerTransitions") {
                newlab <- paste(ElmntTbl1[ii,"Element"], ElmntTbl1[ii,"Transition"], sep="")
                newlab <- gsub(" ", "", newlab, fixed=TRUE)  #drop white spaces
                DispLab(list(x=energy,y=posLbl[ii]), label=newlab, arrow, direction="right", color="blue")
             }
             kk <- kk+1 #kk runs on the ElmntTbl1 rows: it is not simply ii+jj because every new ii jj restarts from 1
          }
       }
    }
    return(ElmntTbl1)
}



#---Identify peaks in the survey spectrum
peakIdentificationOld <- function(peaks, DeltaEnergy, object, SpectrumType, ElmtList) {
   ## init
	PeaksLen <- length(peaks$table$BE) # length identified peaks
	## loop on peaks starting from that having highest intensity
	EleTable <- data.frame()
	BE<-FALSE  #set KE scale
	if (object@Flags[1]) { BE=TRUE } # if BE TRUE if spectra are on a BE scale

	for (peakIndex in 1:PeaksLen) {
       posLbl <- NULL
       elmLbl <- NULL
#	    DEn = 0.25
	    ## if the peak is not assigned find for the element having higher RSF to that BE value
       NewElement<-NULL
       energy <- peaks$table$BE[peakIndex]
       NewElements <- FindElmtList(energy, DeltaEnergy, BE, ElmtList)
       if (is.null(NewElements)){
          idx<-max(which(object[[1]] > energy))  #find all elements > energy
          intensity<-max(object[[2]][(idx-5):(idx+5)])
          newlab<-"??"
          DispLab(list(x=energy,y=intensity), label=newlab)
		    EE<-as.character(peaks$table$BE[peakIndex])
          cat("\n Warning: No Element found at energy ", EE, " Please check Energy scale Alignment!")
       }
       EleTable<-rbind(EleTable, NewElements) #for a given energy a list of element falling in the range: energy-DEn, energy+DEn is produced
       elmLbl<-rbind(elmLbl, NewElements)     #here only the elements found nearby the energy value
       if (is.null(elmLbl)) {
          # do nothing!
       } else {
          idx<-max(which(object[[1]]> energy))  # which gives the index of all the elements of object$x > energy
          intensity<-max(object[[2]][(idx-5):(idx+5)])
          ##print LBL elementi trovati
          Nrw<-nrow(elmLbl)
          Yrng<-range(object[[2]])
          DeltaY<-(Yrng[2]-Yrng[1])/20          #vertical shift to not superpose labels
          if ((intensity+Nrw*DeltaY)>Yrng[2]){  #if TRUE => out of Y boundaries
             intensity <- Yrng[2] - Nrw*DeltaY  #shift to lower y
          }
          for (ii in 1:Nrw){
             posLbl[ii]<-intensity+DeltaY*ii
          }
          arrow=TRUE
          for(ii in 1:Nrw){
             if (length(elmLbl)>0 && SpectrumType =="CoreLines"){
                newlab <- paste(elmLbl[ii,"Element"], elmLbl[ii,"Orbital"], sep="")
                newlab <- gsub(" ", "", newlab, fixed=TRUE)  #drop white spaces
                DispLab(list(x=energy,y=posLbl[ii]), label=newlab, arrow, direction="right", color="blue")
             } else if (length(elmLbl)>0 && SpectrumType == "AugerTransitions") {
                newlab <- paste(elmLbl[ii,"Element"], elmLbl[ii,"Transition"], sep="")
                newlab <- gsub(" ", "", newlab, fixed=TRUE)  #drop white spaces
                DispLab(list(x=energy,y=posLbl[ii]), label=newlab, arrow, direction="right", color="blue")
             }
             arrow=FALSE
          }
       }
   }
	  names(EleTable)<-c("Element", "Orbital", "BE", "KE", "RSF_K", "RSF_S")
   return(EleTable)
}


#---Function to extract list of elements in a given energy interval
FindElmtList<-function(energy, DeltaEnergy, BE, ElmtList){

   EleTable<-NULL
   if (BE) {
      Escale="BE"
   } else {
      Escale="KE"
   }
   ElmtList[,Escale] <- as.numeric(ElmtList[,Escale])
   EleTable <- subset(ElmtList, ElmtList[,Escale] >= energy-DeltaEnergy) #extract all element with energy >= energy+DeltaEnergy
	  EleTable <- subset(EleTable, EleTable[,Escale] <= energy+DeltaEnergy) #from that list select those with energy <= energy+DeltaEnergy
   if (nrow(EleTable) == 0){ EleTable<-NULL }
   return(EleTable)
}



#--- find first element in list whose position is close to 'energy' less than DeltaEnergy
NearerElement <- function(energy, DeltaEnergy, object, ElmtList) {


	energy<-as.numeric(energy)
	BE<-FALSE  #imposto KE
	if (object@Flags[1]) { BE=TRUE } # energy scale is BINDING ENERGY
   EleTable<-NULL
   DEn = 0.3
   while(is.null(EleTable)){
       EleTable<-FindElmtList(energy, DEn, BE, ElmtList)
       DEn <- DEn + 0.2
       if (DEn >  DeltaEnergy) {
          cat("\n No elements found with the selected precision!")
          break
       }
   }

   idx<-max(which(object[[1]]> energy)) # which gives the index of all the elements of object$x > energy
   intensity<-max(object[[2]][(idx-5):(idx+5)])

   if (length(EleTable)>0 && names(ElmtList)[2] =="Orbital"){  #list of corelines was selected
       newlab <- paste(EleTable[1,"Element"], EleTable[1,"Orbital"], sep="")
       newlab <- gsub(" ", "", newlab, fixed=TRUE)  #drop white spaces
       DispLab(list(x=energy,y=intensity), label=newlab, arrow=TRUE, direction="right")
   } else if (length(EleTable)>0 && names(ElmtList)[2] == "Transitions") { #list of Auger Transitions was selected
       newlab <- paste(EleTable[1,"Element"], EleTable[1,"Transition"], sep=" ")
       newlab <- gsub(" ", "", newlab, fixed=TRUE)  #drop white spaces
       DispLab(list(x=energy,y=intensity), label=newlab, arrow=TRUE, direction="right")
   } else {
      newlab<-"??"
      DispLab(list(x=energy,y=intensity), label=newlab)
   }
   if (!is.null(EleTable)) {
      space<-data.frame(Element=" ", Orbital=" ", BE=" ", KE=" ", RSF_K=" ", RSF_S=" ")
      EleTable<-rbind(EleTable, space)
   }
   return(EleTable)
}


#--- Check All elements which RSFmax distance is less than DeltaEnergy
CoreLinesMaxRSF <- function(energy, DeltaEnergy=3, object, ElmtList) { #, index=1
   posLbl <- NULL
	  energy<-as.numeric(energy)
	  if (object@Flags[1]) { BE=TRUE } # energy scale is BINDING ENERGY

   EleTable<-FindElmtList(energy, DeltaEnergy, BE, ElmtList)
   if (object@Flags[3] == TRUE){
      idx<-which(EleTable[ ,6] == max(EleTable[ ,6]), arr.ind = TRUE) #search for the max in the RSF column Scienta
   } else {
      idx<-which(EleTable[ ,5] == max(EleTable[ ,5]), arr.ind = TRUE) #search for the max in the RSF column Kratos
   }
   EleTable<-EleTable[idx,] #select the EleTable-row with max RSF
   Nrw<-nrow(EleTable)
   if (is.null(Nrw)) {
#--- object[[1]] $x values,   object[[2]] $y values
      idx<-max(which(object[[1]] > energy)) # which gives the index of all the elements of object$x > energy
      intensity<-max(object[[2]][(idx-5):(idx+5)])
      newlab<-"??"
      DispLab(list(x=energy,y=intensity), label=newlab)
   } else if (Nrw > 0) {
      ii=1
#--- If same element repeats itself in the EleTable only that wit max RSF is kept
      while(ii < Nrw){ #cannot use a FOR-loop since it is not sensitive to Nrw changes
           element<-EleTable[ii,1] #
           RepIndx<-grep(element, EleTable[ ,1]) #how many times an element is present in the EleTable
           MaxIndx<-which(EleTable[RepIndx,3]==max(EleTable[RepIndx,3]), arr.ind=TRUE)  #select the max RSF among EleTable-rows
           LL<-length(RepIndx)   #N. rows of tmp
           for (jj in LL:1){ #ATTENZIONE: regressive FOR because rows are cancelled from the EleTable
               if (RepIndx[jj] != RepIndx[MaxIndx]){  #if the EleTable-row hasn't max RSF
                   EleTable<-EleTable[-RepIndx[jj], ] #drop the EleTable-row
                   Nrw<-Nrw-1
               }
           }
           ii<-ii+1
      }

      Nrw<-nrow(EleTable)
      Yrng<-range(object[[2]])
      DeltaY<-(Yrng[2]-Yrng[1])/20  #vertical shift to not superpose labels
      idx<-max(which(object[[1]] > energy)) # which gives the index of all the elements of object$x > energy
      intensity<-object[[2]][idx]
      if ((intensity+Nrw*DeltaY)>Yrng[2]){  #if TRUE => out of Ygraph boundaries
          intensity <- Yrng[2] - Nrw*DeltaY
      }
      for (ii in 1:Nrw){
         posLbl[ii]<-intensity+DeltaY*ii
      }
      arrow=TRUE

      for (ii in 1:Nrw){
          if (length(EleTable)>0 && names(ElmtList)[2] =="Orbital"){  #list of corelines was selected
             newlab <- paste(EleTable[ii,"Element"], EleTable[ii,"Orbital"], sep="")
             newlab <- gsub(" ", "", newlab, fixed=TRUE)  #drop white spaces
             DispLab(list(x=energy,y=posLbl[ii]), label=newlab, arrow, direction="right", color="blue")
          } else if (length(EleTable)>0 && names(ElmtList)[2] == "Transitions") { #list of Auger Transitions was selected
             newlab <- paste(EleTable[ii,"Element"], EleTable[ii,"Transition"], sep=" ")
             newlab <- gsub(" ", "", newlab, fixed=TRUE)  #drop white spaces
             DispLab(list(x=energy,y=posLbl[ii]), label=newlab, arrow, direction="right", color="blue")
          }
          arrow=FALSE
      }
   } 
	  return(EleTable)
}


#--- Find ANY element with energy
AllElements <- function(energy, DeltaEnergy=2, object, ElmtList) {

   posLbl <- NULL
   posY <- NULL
	  energy<-as.numeric(energy)
	  if (object@Flags[1]) { BE=TRUE } # energy scale is BINDING ENERGY
   EleTable<-FindElmtList(energy, DeltaEnergy, BE, ElmtList)

   Nrw<-nrow(EleTable)
   if (is.null(Nrw)) {
       idx<-max(which(object[[1]] > energy)) # which gives the index of all the elements of object$x > energy
       intensity<-max(object[[2]][(idx-5):(idx+5)])
       newlab<-"??"
       DispLab(list(x=energy,y=intensity), label=newlab, arrow=TRUE, direction="right")
   } else if (Nrw > 0) {
       idx<-max(which(object[[1]] > energy))
       intensity<-max(object[[2]][(idx-5):(idx+5)])
       Yrng<-range(object[[2]])
       DeltaY<-(Yrng[2]-Yrng[1])/20  #vertical shift to not superpose labels
       if ((intensity+Nrw*DeltaY)>Yrng[2]){  #if TRUE => out of Ygraph boundaries
           intensity <- Yrng[2] - Nrw*DeltaY
       }
       for (ii in 1:Nrw){
          posLbl[ii]<-intensity+DeltaY*ii
       }
       arrow=TRUE
       for (ii in 1:Nrw){
          idx<-max(which(object[[1]] > energy))
          intensity<-object[[2]][idx]
          if (length(EleTable)>0 && names(ElmtList)[2] =="Orbital"){  #list of corelines was selected
             newlab <- paste(EleTable[ii,"Element"], EleTable[ii,"Orbital"], sep="")
             newlab <- gsub(" ", "", newlab, fixed=TRUE)  #drop white spaces
             DispLab(list(x=energy,y=posLbl[ii]), label=newlab, arrow, direction="right", color="blue")
          } else if (length(EleTable)>0 && names(ElmtList)[2] == "Transitions") { #list of Auger Transitions was selected
            newlab <- paste(EleTable[ii,"Element"], EleTable[ii,"Transition"], sep=" ")
            newlab <- gsub(" ", "", newlab, fixed=TRUE)  #drop white spaces
            DispLab(list(x=energy,y=posLbl[ii]), label=newlab, arrow, direction="right", color="blue")
          }
          arrow<-FALSE
       }
   }
	return(EleTable)
}

#--- prints the ElementTable in a formatted fashion
ShowTablePeaks <- function(EleTable, SpectrumType) {
   if(SpectrumType == "CoreLines") { Wdth=10 }
   if(SpectrumType == "AugerTransitions") { Wdth=12 }
   cat("\n\n\n                   *****  IDENTIFIED ELEMENTS ***** \n")

   LL<-nrow(EleTable)
   sss<-format(names(EleTable), justify="right", width=Wdth)
   cat("\n", sss)
   for(ii in 1:LL){
       sss<-as.character(gsub(" ", "", EleTable[ii, ], fixed=TRUE))
       sss<-format(sss, justify="right", width=Wdth)
       cat("\n",sss)
   }
}

plotPeaks <- function(object, peaks=NULL, type = c("normal","corrected"), ...) {
	  if (is.null(peaks) ) stop("\npeaks is NULL")
	  type <- match.arg(type)
	  ## plot curve
	  if (type == "normal") {
		    plot(object)
		    lines(object$x, peaks$baseline[1,], col=2) ## baseline
		    points(peaks$table$BE, peaks$table$orig, col=3, cex=1.0) ## peaks
	  } else {
		    plot(object$x, peaks$corrected[1,], type="l", xlim=rev(c(range(object$x))), xlab=object@units[1], ylab=object@units[2], main = object@Symbol  )
		    points(peaks$table$BE, peaks$table$corr, col=3, cex=1.0) ## peaks
	  }
	  ## Peaks Label
	  if ( ! all(is.na(names(peaks$table$BE)) ) ) {
		    for (PeakBEIndx in seq_along(peaks$table$BE)) {
		 	     if (type == "normal") {
		 		       DispLab(list(x=peaks$table$BE[PeakBEIndx],
		 					            y=peaks$table$orig[PeakBEIndx]),
		 					            label=names(peaks$table$BE[PeakBEIndx]),
		 					            arrow=TRUE,
                    direction="right",
                    ...)
		 	     } else {
		 	        DispLab(list(x=peaks$table$BE[PeakBEIndx],
		 				             y=peaks$table$corr[PeakBEIndx]),
		 				             label=names(peaks$table$BE[PeakBEIndx]),
		 				             arrow=TRUE,
                    direction="right",
                    ...)
		 				  }
		    }
	  }
}


DispLab <- function(position.list,
						              label,
						              color="blue",
						              arrow = TRUE,
						              xspanf = 1, yspanf = 3, lenarrow = 0.08,
						              direction = c("right", "left", "center"), 
                    ...) {

	  # position.list = list(x,y)
	  direction <- match.arg(direction)

	  X0 <- position.list$x
	  # orizontal segment xspanf % of xrange
	  xspan <- abs(par("usr")[2]-par("usr")[1]) * xspanf/100
	  # vertical segment 3% of yrange
	  yspan <- abs(par("usr")[4]-par("usr")[3]) * yspanf/100

	  Yspan <- rep(yspan, length(X0))
	  # reverse yspan if y+2*yspan > ymax
	  if (any((position.list$y + 2*yspan) > par("usr")[4])) {
	    	Yspan[which((position.list$y + 2*yspan) > par("usr")[4])] <- -yspan
		 }

	  Y0 <- position.list$y + 2*Yspan
	  Y1 <- position.list$y + Yspan
 	 Y2 <- position.list$y + 3*Yspan
	  # arrow from(x0,y0) to (x1,y1)
	  if (arrow) {
		 arrows(X0, Y0, X0, Y1, length = lenarrow, col = color)
		 switch(direction,
		        "right" = {
			            X1 <- X0 + xspan
               #segments(X0,Y0, X1, Y0, col=color)
			            text(X1, Y0, labels = label, col=color, pos=4, offset=0.5, ...)## etichetta
			       },
		        "left" = {
			            X1 <- X0 - xspan
               #segments(X0,Y0, X1, Y0, col=color)
			            text(X1, Y0, labels = label, col=color, pos=2, offset=0.5, ...)## etichetta
			       },
		        "center" = {
			            text(X0, Y0, labels = label, col=color, pos=3, offset=0.5, ...)## etichetta
		       	}
		       )
	   }
	   else text(X0, Y0, labels = label, pos=4, col=color, ...)
}

