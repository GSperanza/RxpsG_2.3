###
### Classes and methods for describing fit components
### 

#'Class "fitComponents"
#'
#'A class that describes a fit component.  The class has the definition of all
#'the necessary informations needed to use a component algorithm with the
#'optimisation framework and the graphical user interface (but see Notes
#'below).
#'
#'
#'@docType class
#'@note The goal is that the optimisation framework and the GUI code should get
#'all information about available baseline algorithms through a list of
#'\code{fitComponents} objects.  This will make it relatively simple to add new
#'baseline algorithms.
#'@section Objects from the Class: Objects can be created by calls of the form
#'\code{new("fitComponents", ...)}.
#'@author Roberto Canteri
#'@keywords classes
#'@examples
#'
#'\dontrun{
#'showClass("fitComponents")
#'}
#'@keywords classes
#'@docType class
#'
#'@export
#'

setClass("fitComponents",
         representation(funcName = "character",
                        description = "character",
                        label = "character",
                        param = "data.frame",
                        rsf = "numeric",
                        ycoor = "numeric",
                        link = "list"
                        ),
         prototype(param = data.frame( start = NA, min = NA, max = NA)[0,] ),
         validity = function(object) {
             if (!identical(names(object@param),
                        c("start","min","max")))
                 return("The param slot does not have the correct column names")
                 return(TRUE)
             }
         )

## Accessor functions:
#setGeneric("label", function(object) standardGeneric("label"))
#setMethod("label", "fitComponents", function(object) object@label)
#setGeneric("funcName", function(object) standardGeneric("funcName"))
#setMethod("funcName", "fitComponents", function(object) object@funcName)

setGeneric("setParam", function(object, parameter=NULL, variable=NULL, value= NULL) standardGeneric("setParam"))
setMethod("setParam", "fitComponents", function(object, parameter=NULL, variable=NULL, value=NULL) 
  {
     if ( is.null(value) ) stop(" value is required.")
     rownamesParam <- rownames(slot(object,"param"))
     if ( ! is.null(parameter) ) {
        parameter <- match.arg(parameter,c("start","min", "max"))
        if ( ! is.null(variable) ) {
           variable <- match.arg(variable, rownamesParam)
           slot(object,"param")[variable,parameter] <- value
        } else {
           if (length(value) != length(rownamesParam) )
              stop("wrong value length")
           else
              slot(object,"param")[,parameter] <- value
        }
     } else {
        if ( is.null(variable) ) {
             warning(" if 'parameter' = NULL then 'variable' must be != NULL.")
        } else {
             variable <- match.arg(variable, rownamesParam)
             if (length(value) != 3 ) { stop("wrong value length") }
             else { slot(object,"param")[variable,] <- value }
        }
    }
    invisible(object)
    }
)

setGeneric("getParam", function(object, parameter=NULL, variable=NULL) standardGeneric("getParam"))
setMethod("getParam", "fitComponents", function(object, parameter=NULL, variable=NULL) 
  {
     if ( ! is.null(parameter) ) {
        parameter <- match.arg(parameter,c("start","min", "max"))
        if ( is.null(variable) ) { value <- slot(object,"param")[,parameter] }
        else {
           variable <- match.arg(variable, rownames(object@param))
           value <- slot(object,"param")[variable,parameter]
        }
     } else {
        if ( is.null(variable) ) {
           warning(" if 'parameter' = NULL then 'variable' must be != NULL.")
        } else {
           variable <- match.arg(variable, rownames(object@param))
           value <- slot(object,"param")[variable,]
        }
     }
     return(value)
  }
)

##==============================================================================
# returns a list with the y values for the component n
##==============================================================================
setGeneric("Ycomponent", function(object, x, y)   standardGeneric("Ycomponent"))

setMethod("Ycomponent", signature(object="fitComponents"),
	function(object, x, y)
	{

		## let's take the starting values
		starter <- getParam(object, parameter = "start")
		# combine x with start values
		parm <- list(x=x) # object@RegionToFit$x
		parm <- c(parm,starter)
		# formula e parameters to do the call
		fmla <- slot(object,"funcName")
		names(parm) <- formalArgs(fmla)
		## calculate y values
		ycomponent <- do.call(fmla, parm)
		## exit with y + Baseline$y
		slot(object, "ycoor") <- ( ycomponent + y ) # object@Baseline$y
		return(object)
	}
)


## A list with objects for each implemented algorithm:
fitAlgorithms <- list(

#Modified Giorgio 17-1-2017 to insert empty component structure (needed for VBTop procedure

Initialize = new("fitComponents",
    funcName = "Initialize",
    description = "Symmetric Function(h, mu, sigma)",
    label = "",
    param = data.frame(
        row.names = c("h", "mu", "sigma"),     # min and max values are set in XPSaddFitComponent()
        start = c(NA, NA, NA),
        min = c(NA, NA, NA),
        max = c(NA, NA, NA) ),
	 ycoor = 0,
	 rsf = 0,
	 link = list()
    ),

Generic = new("fitComponents",
    funcName = "Generic",
    description = "Generic Function(x)",
    label = "",
    param = data.frame(
        start = c(NA, NA, NA),
        min = c(NA, NA, NA),
        max = c(NA, NA, NA) ),
	 ycoor = 0,
	 rsf = 0,
	 link = list()
    ),

Gauss = new("fitComponents",
    funcName = "Gauss",
    description = "Symmetric Gaussian shape",
    label = "Gauss",
    param = data.frame(
        row.names = c("h", "mu", "sigma"),
        start = c(1, NA, 1),
        min = c(0, 0, 0.1),
        max = c(Inf, Inf, 10) ),
    ycoor = 0,
	 rsf = 0,
	 link = list()
    ),

Lorentz = new("fitComponents",
    funcName = "Lorentz",
    description = "Symmetric Lorentz shape",
    label = "Lorentz",
   	param = data.frame(
        row.names = c("h", "mu", "sigma"),
        start = c(1, NA, 1),
        min = c(0, 0, 0.1),
        max = c(Inf, Inf, 10) ),
    ycoor = 0,
	 rsf = 0,
	 link = list()
    ),

Voigt = new("fitComponents",
    funcName = "Voigt",
    description = "Symmetric Voigt shape",
    label = "Voigt",
    param = data.frame(
        row.names = c("h", "mu", "sigma", "lg"),
        start = c(1, NA, 0.3, 0.2),
        min = c(0, 0, 0.05, 0.01),
        max = c(Inf, Inf, 10, 1) 	),
    ycoor = 0,
	 rsf = 0,
	 link = list()
    ),

Sech2 = new("fitComponents",
    funcName = "Sech2",
    description = "Symmetric Sech2 shape",
    label = "Sech2",
    param = data.frame(
        row.names = c("h", "mu", "sigma"),
        start = c(1, NA, 0.5),
        min = c(0, 0, 0.1),
        max = c(Inf, Inf, 10) ),
    ycoor = 0,
	 rsf = 0,
	 link = list()
    ),

GaussLorentzProd = new("fitComponents",
    funcName = "GaussLorentzProd",
    description = "Symmetric Gaussian Lorentz cross product form",
    label = "GaussLorentzProd",
    param = data.frame(
        row.names = c("h", "mu", "sigma", "lg"),
        start = c(1, NA, 1.2, 0.75),
        min = c(0, 0, 0.1, 0.01),
        max = c(Inf, Inf, 10, 1) ),
    ycoor = 0,
	 rsf = 0,
	 link = list()
    ),

GaussLorentzSum = new("fitComponents",
    funcName = "GaussLorentzSum",
    description = "Symmetric Gaussian Lorentz Sum Form",
    label = "GaussLorentzSum",
    param = data.frame(
        row.names = c("h", "mu", "sigma", "lg"),
        start = c(1, NA, 1.2, 0.95),
        min = c(0, 0, 0.1, 0),
        max = c(Inf, Inf, 10, 1) ),
    ycoor = 0,
	 rsf = 0,
	 link = list()
    ),

AsymmGauss = new("fitComponents",
    funcName = "AsymmGauss",
    description = "Simple Asymmetric Gaussian shape",
    label = "AsymmGauss",
    param = data.frame(
        row.names = c("h", "mu", "sigma", "asym"),           #min e max values are set in XPSaddFitComponent()
        start = c(1, NA, 1, 0.3),
        min = c(0, 0, 0.05, 0.01),
        max = c(Inf, Inf, 10, 1) 	),
	 ycoor = 0,
	 rsf = 0,
	 link = list()
    ),

AsymmVoigt = new("fitComponents",
    funcName = "AsymmVoigt",
    description = "Asymmetric Voigt shape",
    label = "AsymmVoigt",
    param = data.frame(
        row.names = c("h", "mu", "sigma", "lg", "asym"),     #min e max values are set in XPSaddFitComponent()
        start = c(1, NA, 0.4, 0.1, 0.1),
        min = c(0, 0, 0.1, 0.01, 0.01),
        max = c(Inf, Inf, 10, 1, 1) ),
	 ycoor = 0,
	 rsf = 0,
	 link = list()
    ),


AsymmGaussLorentz = new("fitComponents",
    funcName = "AsymmGaussLorentz",
    description = "Asymmetric Gaussian Lorentz used in Genplot",
    label = "AsymmGaussLorentz",
    param = data.frame(
        row.names = c("h", "mu", "sigma", "lg", "asym"),
        start = c(1, NA, 1, 0.3, 0.45),
        min = c(0, 0, 0.05, 0.4, 1e-5),
        max = c(Inf, Inf, 5, 1, 0.4)  ),
    ycoor = 0,
	 rsf = 0,
	 link = list()
    ),

AsymmGaussVoigt = new("fitComponents",
    funcName = "AsymmGaussVoigt",
    description = "Asymmetric Line-Shape from CasaXPS",
    label = "AsymmGaussVoigt",
    param = data.frame(
        row.names = c("h", "mu", "sigma", "lg", "asym", "gv"),
        start = c(1, NA, 0.6, 0.1, 0.3, 0.4),
        min = c(0, 0, 0.01, 0.01, 0.01, 0.01),
        max = c(Inf, Inf, 5, 1, 1, 1)  ),
    ycoor = 0,
	 rsf = 0,
	 link = list()
    ),

AsymmGaussLorentzProd = new("fitComponents",
    funcName = "AsymmGaussLorentzProd",
    description = "Asym Gaussian Lorentz cross product from UNIFIT Publication",
    label = "AsymmGaussLorentzProd",
    param = data.frame(
        row.names = c("h", "mu", "sigma", "asym", "lg"),
        start = c(1, NA, 1, 0.3, 0.8),
        min = c(0, 0, 0.05, 0.1, 1e-5),
        max = c(Inf, Inf, 5, 1, 0.4)  ),
    ycoor = 0,
	 rsf = 0,
	 link = list()
    ),


DoniachSunjic = new("fitComponents",
    funcName = "DoniachSunjic",
    description = "Asymm. DoniachSunjic shape",
    label = "DoniachSunjic",
    param = data.frame(
        row.names = c("h", "mu", "sigmaDS", "asym"),
        start = c(1, NA, 1, 0.05),
        min = c(0, 0, 0.01, 0.01),
        max = c(Inf, Inf, 10, 1)),
    ycoor = 0,
	 rsf = 0,
	 link = list()
    ),

DoniachSunjicTail = new("fitComponents",
    funcName = "DoniachSunjicTail",
    description = "Asymm. DoniachSunjic shape Tail corrected",
    label = "DoniachSunjicTail",
    param = data.frame(
        row.names = c("h", "mu", "sigmaDS", "asym", "tail"),
        start = c(1, NA, 1, 0.05, 0.6),
        min = c(0, 0, 0.01, 0.01, 0.01),
        max = c(Inf, Inf, 10, 1, 5)  ),
    ycoor = 0,
	 rsf = 0,
	 link = list()
    ),

DoniachSunjicGauss = new("fitComponents",
    funcName = "DoniachSunjicGauss",
    description = "Asymm. DoniachSunjic shape + GaussBroadening",
    label = "DoniachSunjicGauss",
    param = data.frame(
        row.names = c("h", "mu", "sigmaDS", "sigmaG", "asym"),
        start = c(1, NA, 1, 0.8, 0.15),
        min = c(0, 0, 0.01, 0.01, 0.01),
        max = c(Inf, Inf, 5, 5, 1) ),
    ycoor = 0,
	 rsf = 0,
	 link = list()
    ),

DoniachSunjicGaussTail = new("fitComponents",
    funcName = "DoniachSunjicGaussTail",
    description = "Asymm. DoniachSunjic shape + GaussBroadening and Tail correction",
    label = "DoniachSunjicGaussTail",
    param = data.frame(
        row.names = c("h", "mu", "sigmaDS", "sigmaG", "asym", "tail"),
        start = c(1, NA, 1, 0.8, 0.05, 0.8),
        min = c(0, 0, 0.01, 0.01, 0.01, 0.01),
        max = c(Inf, Inf, 5, 5, 1, 5) ),
    ycoor = 0,
	 rsf = 0,
	 link = list()
    ),

SimplfiedDoniachSunjic = new("fitComponents",
    funcName = "SimplfiedDoniachSunjic",
    description = "Asymmetric DoniachSunjic shape",
    label = "SimplfiedDoniachSunjic",
    param = data.frame(
        row.names = c("h", "mu", "sigma", "asym"),
        start = c(1, NA, 0.6, 0.01),
        min = c(0, 0, 0.01, 0.01),
        max = c(Inf, Inf, 5, 1)	),
    ycoor = 0,
	 rsf = 0,
	 link = list()
    ),

Linear = new("fitComponents",
    funcName = "Linear",
    description = "Linear Fit",
    label = "Linear",
    param = data.frame(
        row.names = c("m", "c", "mu"),          #min e max values are set in XPSaddFitComponent()
        start = c(NA, 0, NA),
        min = c(-Inf, -Inf, NA),
        max = c(Inf, Inf, NA) ),
	 ycoor = 0,
	 rsf = 0,
	 link = list()
    ),

ExpDecay = new("fitComponents",
    funcName = "ExpDecay",
    description = "Exponential Decay",
    label = "ExpDecay",
    param = data.frame(
        row.names = c("h", "mu", "k", "c"),     #min e max values are set in XPSaddFitComponent()
        start = c(1, NA, 1, 0),
        min = c(0, 0, 0, 0),
        max = c(Inf, Inf, 5, Inf) ),
	 ycoor = 0,
	 rsf = 0,
	 link = list()
    ),

PowerDecay = new("fitComponents",
    funcName = "PowerDecay",
    description = "Power Law Decay",
    label = "PowerDecay",
    param = data.frame(
        row.names = c("h", "mu", "pow", "c"),   #min e max values are set in XPSaddFitComponent()
        start = c(1, NA, 2, 0),
        min = c(0, 0, 0.1, 0),
        max = c(Inf, Inf, 5, Inf) ),
	 ycoor = 0,
	 rsf = 0,
	 link = list()
    ),

Sigmoid = new("fitComponents",
    funcName = "Sigmoid",
    description = "Sigmoid function",
    label = "Sigmoid",
    param = data.frame(
        row.names = c("h", "mu", "k", "c"),     #min e max values are set in XPSaddFitComponent()
        start = c(1, NA, 1, 0),
        min = c(0, 0, 0, 0),
        max = c(Inf, Inf, 5, Inf) ),
	 ycoor = 0,
	 rsf = 0,
	 link = list()
    ),

HillSigmoid = new("fitComponents",
    funcName = "HillSigmoid",
    description = "Hill Sigmoid function",
    label = "HillSigmoid",
    param = data.frame(
        row.names = c("h", "mu", "pow", "A", "B"),     #mu, A1, A2 values ar set by do.fit() in XPSVBTopGUI()
        start = c(NA, NA, 1, NA, NA),
        min = c(0, 0, 0, 0, 0),
        max = c(Inf, Inf, 5, Inf, Inf) ),
	 ycoor = 0,
	 rsf = 0,
	 link = list()
    ),


VBtop = new("fitComponents",     #Together with FitAlgorithms this is needeed to create a slot
    funcName = "VBtop",          #in the ...@Components[[ii]]@param location of the XPSSample
    description = "VBtop fit",   #where to save the VBtop position.
    label = "VBtop",             #This slot is sensitive to the EnergyShift function
    param = data.frame(
        row.names = c("mu"),
        start = c(NA),
        min = c(NA),
        max = c(NA) ),
	 ycoor = 0,
	 rsf = 0,
	 link = list()
    ),

Derivative = new("fitComponents",     #Together with FitAlgorithms this is needeed to create a slot
    funcName = "Derivative",          #in the ...@Components[[ii]]@param location of the XPSSample
    description = "First Derivative", #where to save the core line first derivative.
    label = "D1",                     #This slot is sensitive to the EnergyShift function
    param = data.frame(
        row.names = c("mu"),
        start = c(NA),
        min = c(NA),
        max = c(NA) ),
	 ycoor = 0,
	 rsf = 0,
	 link = list()
    ),

FitProfile = new("fitComponents",
    funcName = "FitProfile",
    description = "Depth Profile",
    label = "DepthP.",
    param = data.frame(
        row.names = c("R.I", "d", "La", "Lb"),     #min e max values are set in XPSaddFitComponent()
        start = c(NA, 10, NA, NA),
        min = c(NA, 0, NA, NA),
        max = c(NA, 50, NA, NA) ),
	 ycoor = 0,
	 rsf = 0,
	 link = list()
    )
)
