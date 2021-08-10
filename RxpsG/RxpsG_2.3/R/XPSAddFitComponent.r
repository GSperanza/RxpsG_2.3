##  =================================================
#  			XPSaddFitComponent
#  Impostazione delle variabili: numero componenti
#  posizione, Sigma e altezza delle componenti  
##  =================================================

#'Add a fit component to the XPSCoreLine
#'
#'Add a fit component to the XPSCoreLine object. This function is called after
#'the definition of the baseline (\code{XPSbaseline}). Different types of
#'function could be added as component, see \code{\link{XPSfitAlgorithms}} for
#'a list of the implemented functions.
#'
#'@param Object XPSCoreLine object
#'@param type name of the fitting function used. Default "Gauss"
#'@param range_h the value used to determine \code{min} and \code{max} range
#'for \code{Intensity}
#'@param range_mu the value used to determine \code{min} and \code{max} range
#'for \code{xCenter}
#'@param peakPosition a list with x,y value corresponding the center and height
#'of the new component. Default is \code{NULL} then the position is asked
#'trough the plot with the cursor
#'@return Usually the values \code{xCenter} and \code{Amplitude} are get trough
#'a click with the mouse on the plot of the XPSCoreLine. All the parameters
#'will be added to the slot \code{Components} of the XPSCoreLine.
#'@seealso \link{XPSfitAlgorithms}, \link{XPSbaseline}, \link{findXIndex}
#'
#'@examples
#'
#'\dontrun{
#'	XPSdata[["C1s"]] <- XPSaddFitComponent(XPSdata[["C1s"]], type="Voigt", range)
#'}
#'
#'@export
#'

XPSaddFitComponent <- function(Object, type="Gauss", range_h=5, range_mu=1.5, peakPosition=NULL)
{
	Object <- XPSaddComponent(Object, type, range_h, range_mu, peakPosition)
	plot(Object)
	return(Object)
}

#-----------
XPSaddComponent <- function(Object, type, range_h=5, range_mu=1.5, peakPosition=NULL)
{
	# check if Baseline is already present
	if ( ! hasBaseline(Object) ) {
		stop("Baseline not defined for Core Line ", slot(Object,"Symbol") )
   }
#modified Giorgio 17-1-2017

	fit_functions <- sapply(fitAlgorithms, slot, "funcName") #fit_functions=list of all the fitting function-names
	funct <- match.arg(type, unname(fit_functions) )         #selection of the fit function

#modified Giorgio 5-2-2019
	if (type=="Initialize"){                                 #creation of an empty component
	   num <- length(names(Object@Components)) + 1
           slot(Object,"Components")[[num]] <- fitAlgorithms[[funct]]
	   ## Now set h, mu for first call

	   ## set Amplitude:start,min,max
	   slot(Object,"Components")[[num]] <- setParam(Object@Components[[num]],variable="h", value = c(0, 0, 0)) #l'intensita' della componente puo' muoversi tra 0 e 5 volte l'intensita' settata col mouse in XPSprocess add_component

	   ## set xCenter: start,min,max
	   slot(Object,"Components")[[num]] <- setParam(Object@Components[[num]],variable="mu",value = c(0, 0, 0))

	   ## set sigma: start,min,max
	   slot(Object,"Components")[[num]] <- setParam(Object@Components[[num]],variable="sigma",value = c(0, 0, 0))

	   # label (used in the plot) # #1, #2, ...
	   slot(Object@Components[[num]],"label") <- paste("#", as.character(num),sep="")

	   # name of this component
	   names(Object@Components)[num] <- paste("C",as.character(num),sep="") # C1, C2, ...

   	   ## add the y values for this component
	   Object@Components[[num]] <- Ycomponent(Object@Components[[num]], x=Object@RegionToFit$x, y=Object@Baseline$y)

	   ## now set the rsf for this component
	   ## if Object@RSF is not set then call XPSsetRSF
	   if ( slot(Object,"RSF") != 0 ) {slot(Object@Components[[num]],"rsf") <- slot(Object,"RSF")}

	   return(Object)
   }

#modified Giorgio 5-2-2019
	if (type=="Generic"){  #creation of an empty component
	   num <- length(names(Object@Components)) + 1
   	   slot(Object,"Components")[[num]] <- fitAlgorithms[[funct]]
	   return(Object)
   }


	if (type=="Linear"){   #creation of an empty component of type Linear
	   num <- length(names(Object@Components)) + 1
   	   slot(Object,"Components")[[num]] <- fitAlgorithms[[funct]]
	   ## Now set h, mu for first call

	   ## set Amplitude:start,min,max
	   m <-(peakPosition$y[2]-peakPosition$y[1])/(peakPosition$x[1]-peakPosition$x[2]) #line slope
	   slot(Object,"Components")[[num]] <- setParam(Object@Components[[num]],variable="m", value = c(m, -Inf, Inf)) #l'intensita' della componente puo' muoversi tra 0 e 5 volte l'intensita' settata col mouse in XPSprocess add_component

	   ## set xCenter: start,min,max
	   slot(Object,"Components")[[num]] <- setParam(Object@Components[[num]],variable="c",value = c(0, -Inf, Inf))

	   # label (used in the plot) # #1, #2, ...
	   slot(Object@Components[[num]],"label") <-  paste("#", as.character(num),sep="")  #labels of the components in the plot

	   # name of this component
	   names(Object@Components)[num] <- paste("C",as.character(num),sep="") # C1, C2, ...

   	   ## add the y values for this component
	   Object@Components[[num]] <- Ycomponent(Object@Components[[num]], x=Object@RegionToFit$x, y=Object@Baseline$y)

	   ## now set the rsf for this component
	   ## if Object@RSF is not set then call XPSsetRSF
	   if ( slot(Object,"RSF") != 0 ) {slot(Object@Components[[num]],"rsf") <- slot(Object,"RSF")}


	   return(Object)
   }


#modified Giorgio 10-11-2018
   if (type=="HillSigmoid"){   #creation of an empty component of type VBtop: needed to store VBtop Position

      num <- length(names(Object@Components)) + 1

      ## Set the x,y coord
      X <- peakPosition$x[2]   #this corresponds to the mu parameter = position of the Sigmoid flex point

      # substract the baseline value correspondent to the X position
      idx <- findXIndex(Object@Baseline$x, X)
      Y <- peakPosition$y[2] - Object@Baseline$y[idx]
      RTF_x<-Object@RegionToFit$x
      #HillSigmoid needs positive abscissas defined from 0 to LL*dx  dx=energy step
      dx<-abs(RTF_x[2]-RTF_x[1])
      LL<-length(RTF_x)
      XXX<-NULL
      for (ii in 1:LL){        #new X coords: Step and number of absissas are equal to the original ones
          XXX[ii]<-dx*(ii-1)
      }
      newX<-XXX[idx]

      slot(Object,"Components")[[num]] <- fitAlgorithms[[funct]]

      ## set Amplitude:start,min,max
      slot(Object,"Components")[[num]] <- setParam(Object@Components[[num]],variable="h", value = c(Y, 0, Y*1.5)) #l'intensita' della componente puo' muoversi tra 0 e 5 volte l'intensita' settata col mouse in XPSprocess add_component

      ## set xCenter: start,min,max
      slot(Object,"Components")[[num]] <- setParam(Object@Components[[num]],variable="mu",value = c(newX, 0, newX+5))

      ## set Powr: start,min,max
      slot(Object,"Components")[[num]] <- setParam(Object@Components[[num]],variable="pow",value = c(2, 1, 10))

      ## set A: start,min,max
      slot(Object,"Components")[[num]] <- setParam(Object@Components[[num]],variable="A",value = c(peakPosition$y[1], 0.5*peakPosition$y[1], 1.5*peakPosition$y[1]))

      ## set B: start,min,max
      slot(Object,"Components")[[num]] <- setParam(Object@Components[[num]],variable="B",value = c(peakPosition$y[3], -1, 1.5*peakPosition$y[3]))

	   ## Now just set the component identifier
	   # label (used in the plot)
	   slot(Object@Components[[num]],"label") <-"HS"   #Label indicating the Derivate in the plot

	   # name of this component
	   names(Object@Components)[num] <- "C1"  # The name of the HillSigmoid component

	   Object@Components[[num]] <- Ycomponent(Object@Components[[num]], x=XXX, y=Object@Baseline$y) #Ycomponent saves the fitting component but does not modify the RegionToFit$x
      Object@Fit <- list(x=XXX, y=NULL, idx=idx)  #save new absissas and index of flexpoint
	   return(Object)
   }


#modified Giorgio 17-1-2017
	if (type=="VBtop"){   #creation of an empty component of type VBtop: needed to store VBtop Position
	   num <- length(names(Object@Components)) + 1
   	slot(Object,"Components")[[num]] <- fitAlgorithms[[funct]]

   	## add the y values for this component: in this case a series of NA (see XPSFitAlgorithms.r)
	   Object@Components[[num]] <- Ycomponent(Object@Components[[num]], x=Object@RegionToFit$x, y=Object@Baseline$y)

	   ## Now just set the component identifier
	   # label (used in the plot) # #1, #2, ...
	   slot(Object@Components[[num]],"label") <-"VBtop"   #Label indicating the VBtop in the plot

	   # name of this component
	   names(Object@Components)[num] <- paste("VBt",as.character(num),sep="") # V name of the VBtop component

	   return(Object)
   }

	if (type=="Derivative"){   #creation of an empty component of type VBtop: needed to store VBtop Position

	   num <- length(names(Object@Components)) + 1
   	slot(Object,"Components")[[num]] <- fitAlgorithms[[funct]]

   	## add the y values for this component: in this case a series of NA (see XPSFitAlgorithms.r)
	   Object@Components[[num]] <- Ycomponent(Object@Components[[num]], x=Object@RegionToFit$x, y=Object@Baseline$y)

	   ## Now just set the component identifier
	   # label (used in the plot)
	   slot(Object@Components[[num]],"label") <-"D1"   #Label indicating the Derivate in the plot

	   # name of this component
	   names(Object@Components)[num] <- "D1"  # The name of the Derivate component

	   return(Object)
   }



	if (is.null(peakPosition)) {
		## get the position from the plot
		plot(Object)
		peakPosition <- locator(1)
   }
	## Set the x,y coord
	X <- peakPosition$x   #this corresponds to the mu parameter = position of the fitting comp.
	# substract the baseline value correspondent to the X position
   idx<-findXIndex(Object@Baseline$x, X)
	Y <- peakPosition$y - Object@Baseline$y[idx]

	# increase the number of components
	num <- length(names(Object@Components)) + 1

	## add the new component: recall the fitAlgorithm with all the preset parameters (start, min, max) values
	#  and fill the slots Component@param,  Component@rsf, Component@ycoor, Component@link
	slot(Object,"Components")[[num]] <- fitAlgorithms[[funct]]
	## Now set h, mu for first call

	## set Amplitude:start,min,max
	slot(Object,"Components")[[num]] <- setParam(Object@Components[[num]],variable="h", value = c(Y, 0, Y*range_h)) #l'intensita' della componente puo' muoversi tra 0 e 5 volte l'intensita' settata col mouse in XPSprocess add_component

	## set xCenter: start,min,max
	slot(Object,"Components")[[num]] <- setParam(Object@Components[[num]],variable="mu",value = c(X, X-range_mu, X+range_mu))

	# label (used in the plot) # #1, #2, ...
	slot(Object@Components[[num]],"label") <- paste("#", as.character(num),sep="")

#--- Calcolo fattore conversione per ottenere initensita' componente coincidente con posizione mouse
#    ATTENZIONE: questa operazione puo' esser eseguita solo quando tutti gli slot della componente aggiunta sono stati definiti

   Yfactor<-GetHvalue(Object, num, type, 1)  #calcola il valore della corrispondente funzione di fit considerati i suoi parametri
#--- Re-set Amplitude   dopo aggiustamento intensita'

	slot(Object,"Components")[[num]] <- setParam(Object@Components[[num]],variable="h", value = c(Yfactor*Y, 0, Yfactor*Y*range_h)) #l'intensita' della componente puo' muoversi tra 0 e 5 volte l'intensita' settata col mouse in XPSprocess add_component

	# name of this component

	names(Object@Components)[num] <- paste("C",as.character(num),sep="") # C1, C2, ...

	## add the y values for this component
	Object@Components[[num]] <- Ycomponent(Object@Components[[num]], x=Object@RegionToFit$x, y=Object@Baseline$y)
	
	## now set the rsf for this component
	## if Object@RSF is not set then call XPSsetRSF
	if ( slot(Object,"RSF") != 0 )
		slot(Object@Components[[num]],"rsf") <- slot(Object,"RSF")
  

	## sort the Components slot
	Object <- sortComponents(Object)
	return(Object)
}

