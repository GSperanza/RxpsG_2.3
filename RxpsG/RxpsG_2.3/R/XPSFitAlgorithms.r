## =======================================================================
## Fit function
## =======================================================================


XPSfitAlgorithms <- function() {
     bho <- sapply(fitAlgorithms, function(x) {
  	       cat(sprintf("%15s  %-30s\n", slot(x,"funcName"), slot(x,"description") ))
  	  })
}
## =======================================================================

#>>>> IMPORTANT: IF NEW FUNCTIONS ARE ADDED ALSO  XPSFitCompClass() and getHvalue() HAS TO BE MODIFIED CORRESPONDINGLY !!!!!!




## =======================================================================
## Symmetric functions
## =======================================================================

Initialize <- function(x, h, mu, sigma) { 
#'Initialize: initalizes parameter slot for a generic function
#'
#'@param x numeric vector
#'@param h function amplitude
#'@param mu position of function center
#'@param sigma function full width at half maximum
#'
#'@export

    return(x*NA)     #initialize the component with NA values
}


## =======================================================================
#'Generic function to create the new Fit Component Slots
#'
#'@param x numeric vector
#'@export
Generic <- function(x) { 

    return(x*NA)     #initialize the component with NA values
}


## =======================================================================
#'Gauss function
#'
#'@param x numeric vector
#'@param h function amplitude
#'@param mu position of function center
#'@param sigma function full width at half maximum
#'
#'@export
Gauss <- function(x, h, mu, sigma) { 

    return( h*exp(-(2.77258872*((x-mu)/sigma)^2)) ) 
}


## =======================================================================
#'Lorentz function
#'
#'@param x numeric vector
#'@param h function amplitude
#'@param mu position of function center
#'@param sigma function full width at half maximum
#'
#'@export
Lorentz <- function(x, h, mu, sigma) { 

     return( h/(1+(4*(((x-mu)/sigma)^2))) ) 
}


## =======================================================================
## Voigt require NORMT3 package
## http://en.wikipedia.org/wiki/Voigt_profile

#'Voigt function
#'
#'@param x numeric vector
#'@param h function amplitude
#'@param mu position of function center
#'@param sigma function full width at half maximum
#'@param lg mix Gauss-Lorentz
#'
#'@export
Voigt <- function(x, h, mu, sigma, lg) {

    import::here(wofz, .from=NORMT3)
    xx <- (x - mu)/sigma/1.41421356     # sqrt(2)=1.41421356
    yy <- lg/sigma/1.41421356
    return( h * Re(wofz( complex(real = xx, imaginary = yy) ))/sigma/2.50662827 )  #sqrt(2*PI)=2.50662827
}


## =======================================================================
## Definition: pulses with a temporal intensity profile which has the shape of a sech2 function
## http://www.rp-photonics.com/sech2_shaped_pulses.html
#'Sech2 function
#'
#'@param x numeric vector
#'@param h function amplitude
#'@param mu position of function center
#'@param sigma function full width at half maximum

#'@export
Sech2 <- function(x, h, mu, sigma) { 

    return( h/( cosh((x-mu)/sigma)*cosh((x-mu)/sigma) ) ) 
}


## =======================================================================
## Gaussian Lorentz cross product
# http://www.casaxps.com/help_manual/line_shapes.htm

#'GaussLorentzProd function
#'
#'@param x numeric vector
#'@param h function amplitude
#'@param mu position of function center
#'@param sigma function full width at half maximum
#'@param lg mix Gauss-Lorentz
#'
#'@export
GaussLorentzProd <- function(x, h, mu, sigma, lg) {

 	  return( h / (1+(4*lg*((x-mu)/sigma)^2)) * exp(-2.77258872*(1-lg)*((x-mu)/sigma)^2) )
}


## =======================================================================
## Gaussian Lorentz Sum form
# http://www.casaxps.com/help_manual/line_shapes.htm

#'GaussLorentzSum function
#'
#'@param x numeric vector
#'@param h function amplitude
#'@param mu position of function center
#'@param sigma function full width at half maximum
#'@param lg mix Gauss-Lorentz
#'
#'@export
GaussLorentzSum <- function(x, h, mu, sigma, lg) {

		  return(h*(lg*exp(-2.77258872*((x-mu)/sigma)^2) + (1-lg)*(1/(1+4*(((x-mu)/sigma)^2)))))
}



## =======================================================================
## Asym functions
## =======================================================================



## =======================================================================
# Tail modified Gauss used in Genplot

#'AsymmGauss function
#'Tail modified Gauss by an exponential function
#'
#'@param x numeric vector
#'@param h function amplitude
#'@param mu position of function center
#'@param sigma function full width at half maximum
#'@param asym function asymmetry value
#'
#'@export
AsymmGauss <- function(x, h, mu, sigma, asym) {

	   Tail <- function(x, mu, asym) {
	   	   Y <- vector("numeric", length=length(x))
	   	   Ex <- (x-mu)
	   	   index <- which(Ex >= 0)
	   	   Y[index] <- exp(-abs(Ex[index])*(1-asym)/asym)
	   	   return(Y)
	   }
	   return( Gauss(x, h, mu, sigma) + (h - Gauss(x, h, mu, sigma))*Tail(x, mu, asym) )
}

## =======================================================================
#'AsymmVoigt
#'
#'@param x numeric vector
#'@param h function amplitude
#'@param mu position of function center
#'@param sigma function full width at half maximum
#'@param lg mix Gauss-Lorentz
#'@param asym function asymmetry value
#'
#'@export
AsymmVoigt <- function(x, h, mu, sigma, lg, asym) {

	   Tail <- function(x, mu, sigma, asym) {
            Y<- asym / (sigma^2 + (x - mu)^2)^(1-asym/4)
            return(Y)
    }

	   Y <- vector("numeric", length=length(x))
	   Ex <- (x-mu)
	   indx1 <- which(Ex < 0)                          #DX part: normal Voigt
	   Y[indx1] <- Voigt(x[indx1], h, mu, sigma, lg)
    MM<-max(Y[indx1])

	   indx2 <- which(Ex >= 0)                        #SX part: Voigt + voigt*tail
#   Y[indx2] <- Voigt(x[indx2], h, mu, sigma, lg) +  (h- Voigt(x[indx2], h, mu, sigma, lg))*Tail(x[indx2], mu, sigma, asym)
    Y[indx2] <- Voigt(x[indx2], h, mu, sigma, lg)+h*Tail(x[indx2], mu, sigma, asym)
#   Y[indx2] <- h* Tail(x[indx2], mu, sigma, asym)
	   Y[indx2]<-MM * Y[indx2]/max(Y[indx2])
	   return (Y)
}


## =======================================================================
# Tail modified Gauss-Lorentz used in Genplot

#'AsymmGaussLorentz function
#'Tail modified Gauss-Lorentz Sum modifed by using an exponential function
#'
#'@param x numeric vector
#'@param h function amplitude
#'@param mu position of function center
#'@param sigma function full width at half maximum
#'@param lg mix Gauss-Lorentz
#'@param asym function asymmetry value
#'
#'@export
AsymmGaussLorentz <- function(x, h, mu, sigma, lg, asym) {

	   Tail <- function(x, mu, asym) {
		          X <- vector("numeric", length=length(x))
		          Ex <- (x-mu)
		          index <- which(Ex >= 0)
		          X[index] <- exp(-abs(Ex[index])*(1-asym)/asym)
		          return(X)
		  }
	   return( GaussLorentzSum(x, h, mu, sigma, lg) + (h - GaussLorentzSum(x, h, mu, sigma, lg))*Tail(x, mu, asym) )
}



## =======================================================================
#'AsymmGaussVoigt Asymmetric Gaussian-Voigt function see CasaXPS
#'
#'@param x numeric vector
#'@param h function amplitude
#'@param mu position of function center
#'@param sigma function full width at half maximum
#'@param lg mix Gauss-Lorentz
#'@param asym function asymmetry value
#'@param gv mix Gauss-Voigt
#'
#'@export
AsymmGaussVoigt <- function(x, h, mu, sigma, lg, asym, gv) {
	   Tail <- function(x, mu, asym) {
		          Y <- vector("numeric", length=length(x))
	          	Ex <- (x-mu)
          		index <- which(Ex >= 0)
          		Y[index] <- exp(-abs(Ex[index])*(1-asym)/asym)
          		return(Y)
	   }
	   return((gv*Gauss(x, h, mu, sigma)+(h - Gauss(x, h, mu, sigma))*Tail(x, mu, asym)) + (1-gv)*Voigt(x, h, mu, sigma, lg) )
}


## =======================================================================
## Asym Gaussian Lorentz cross product from UNIFIT Publication

#'AsymmGaussLorentzProd function
#'Asym Gaussian Lorentz cross product from UNIFIT Publication
#'
#'@param x numeric vector
#'@param h function amplitude
#'@param mu position of function center
#'@param sigma function full width at half maximum
#'@param asym function asymmetry value
#'@param lg mix Gauss-Lorentz
#'
#'@export
AsymmGaussLorentzProd <- function(x, h, mu, sigma, asym, lg) {

	   z <- (x-mu)/(sigma + asym * (x-mu))
   	return(h / (1+(4*lg*(z)^2)) * exp(-2.77258872*(1-lg)*(z)^2) )
}

## =======================================================================
## DoniachSunjic(x, h, mu, sigmaDS, asym)
#Lineshape Doniach -Sunjic corretta pura
#(vedi anche Leiro et al. J. El. Spectr. Rel. Phen. 128, 205, (2003).

#'DoniachSunjic function
#'Correct version of the Doniach-Sunjic function (see Wertheim PRB 25(3), 1987, (1982))
#'
#'
#'@param x numeric vector
#'@param h function amplitude
#'@param mu position of function center
#'@param sigmaDS function full width at half maximum
#'@param asym function asymmetry value
#'
#'@export
DoniachSunjic <- function(x, h, mu, sigmaDS, asym) {

 	DS <- h/4*( (gamma(1-asym)/((mu-x)^2+(sigmaDS/2)^2)^((1-asym)/2) ) * cos((pi*asym)/2+(1-asym)*atan((mu-x)*2/sigmaDS)) )  #DoniachSunjic
        return (DS)
}


## =======================================================================
## DoniachSunjic(x, h, mu, sigmaDS, asym, tail)
#corrected version of  Doniach -Sunjic Lineshape adding a Tail
#on the low BE (high KE) side following the model of Wertheim PRB 25(3), 1987, (1982)
#(see also Leiro et al. J. El. Spectr. Rel. Phen. 128, 205, (2003).
#mu-x instead of x-mu to work on energies < mu
#tail damps the lineshape at low BE

#'DoniachSunjic function
#'Version of the Doniach-Sunjic function (see Wertheim PRB 25(3), 1987, (1982))
#'Corrected for an exponential decay (tail) on the low BE side
#'
#'@param x numeric vector
#'@param h function amplitude
#'@param mu position of function center
#'@param sigmaDS function full width at half maximum
#'@param asym function asymmetry value
#'@param tail amplitude of the spectral tail on the low BE (high KE) side
#'
#'@export
#'
DoniachSunjicTail <- function(x, h, mu, sigmaDS, asym, tail) {

 	  Tail <- function(x, mu, tail) {
		          X <- rep(1, length(x))  #unity vettore
		          Ex <- (mu-x)    #tail at the right of mu (component position)
	          	index <- which(Ex >= 0)
	          	X[index] <- exp(-abs(Ex[index])*(1-tail)/tail)
		          return(X)
	   }

 	  DS <- h/4*( (gamma(1-asym)/((mu-x)^2+(sigmaDS/2)^2)^((1-asym)/2) ) * cos((pi*asym)/2+(1-asym)*atan((mu-x)*2/sigmaDS)) )  #DoniachSunjic
    return (DS * Tail(x, mu, tail) )
}


## =======================================================================
## DSunjicGauss(x, h, mu, sigmaDS, sigmaG, asym)
#a gaussian broadeniing is added by multiplying the DS function with a Gauss

#'DoniachSunjicGauss function
#'Doniach Sunjic function multiplied for a Gaussian broadening
#'
#'@param x numeric vector
#'@param h function amplitude
#'@param mu position of function center
#'@param sigmaDS function full width at half maximum
#'@param sigmaG full width at half maximum of superimpossed Gaussian broadening
#'@param asym function asymmetry value
#'
#'@export
DoniachSunjicGauss <- function(x, h, mu, sigmaDS, sigmaG, asym) {

    DS = ( (gamma(1-asym)/((mu-x)^2+(sigmaDS/2)^2)^((1-asym)/2) ) * cos((pi*asym)/2+(1-asym)*atan((mu-x)*2/sigmaDS)) )  #DoniachSunjic
    maxDS<-max(DS)
    minDS<-min(DS)
    DS<-(DS-minDS)/(maxDS-minDS)                   #normalize
    Gs = exp(-(2.77258872*(0.5*(x-mu)/sigmaG)^2))  #Gauss function

    LL <- length(DS)
    stp <- (DS[LL]-DS[1])/LL
    bkg <- stp*seq(1,LL,1)                         #DS background generation
    DS <- DS-bkg-DS[1]                             #BKG subtraction: recovers convolution distortions inroduced by the non-zer values of the asymmetric tail of the DS

    ConvDSGs<-convolve(DS, rev(Gs), type = "open")
    ConvDSGs<-h/4*maxDS*ConvDSGs/max(ConvDSGs)     #normalize and moltiply by h/4 to get the proper height

    LL=length(ConvDSGs)               #the convolution doubles the original number of data
    X <- ConvDSGs[seq(1,LL,2)]        #decimation 1 each two values
    return(X)
}

## =======================================================================
## DSunjicGauss(x, h, mu, sigmaDS, sigmaG, asym, tail)
#corrected version of  Doniach -Sunjic Lineshape adding a Tail
#on the low BE (high KE) side following the model of Wertheim PRB 25(3), 1987, (1982)
#(see also Leiro et al. J. El. Spectr. Rel. Phen. 128, 205, (2003).
#mu-x instead of x-mu to work on energies < mu
#tail damps the lineshape at low BE
#a gaussian broadeniing is added by multiplying the DS function with a Gauss

#'DoniachSunjicGauss function
#'Doniach Sunjic function multiplied for a Gaussian broadening
#'corrected for an exponential decay on the low BE side
#'
#'@param x numeric vector
#'@param h function amplitude
#'@param mu position of function center
#'@param sigmaDS function full width at half maximum
#'@param sigmaG full width at half maximum of superimpossed Gaussian broadening
#'@param asym function asymmetry value
#'@param tail amplitude of the spectral tail on the low BE (high KE) side
#'
#'@export
DoniachSunjicGaussTail <- function(x, h, mu, sigmaDS, sigmaG, asym, tail) {

 	  Tail <- function(x, mu, tail) {
		          X <- rep(1, length(x))  #vettore contenete tutti 1
		          Ex <- (mu-x)    #coda a DX di mu (posizione componente)
	          	index <- which(Ex >= 0)
	          	X[index] <- exp(-abs(Ex[index])*(1-tail)/tail)
	          	return(X)
	   }

    DS = ( (gamma(1-asym)/((mu-x)^2+(sigmaDS/2)^2)^((1-asym)/2) ) * cos((pi*asym)/2+(1-asym)*atan((mu-x)*2/sigmaDS)) )  #DoniachSunjic
    maxDS<-max(DS)
    minDS<-min(DS)
    DS<-(DS-minDS)/(maxDS-minDS)                    #normalize
    Gs = exp(-(2.77258872*(0.5*(x-mu)/sigmaG)^2))   #Gauss function

    LL <- length(DS)
    stp <- (DS[LL]-DS[1])/LL
    bkg <- stp*seq(1,LL,1)                         #DS background generation
    DS <- DS-bkg-DS[1]                             #BKG subtraction: recovers convolution distortions inroduced by the non-zer values of the asymmetric tail of the DS

    ConvDSGs<-convolve(DS, rev(Gs), type = "open")
    ConvDSGs<-h/4*maxDS*ConvDSGs/max(ConvDSGs)     #normalize and moltiply by h/4 to get the proper height

    LL=length(ConvDSGs)               #convoluzione provides the double of elements
    X <- ConvDSGs[seq(1,LL,2)]        #since DS and Gs have the same number of elements: decimation 1 each two elements
    return (X * Tail(x, mu, tail) )   #generation of DS gaussian broadened + exponential decay on the low BE side
}

## =======================================================================
## SimplfiedDoniachSunjic(x, h, mu, sigma, asym)
# http://www.casaxps.com/help_manual/line_shapes.htm 
# mu-x invece di x-mu per avere asym ad energie superiori a mu

#'SimplfiedDoniachSunjic function
#'see http://www.casaxps.com/help_manual/line_shapes.htm
#'
#'@param x numeric vector
#'@param h function amplitude
#'@param mu position of function center
#'@param sigma function full width at half maximum
#'@param asym function asymmetry value
#'
#'@export
SimplfiedDoniachSunjic <- function(x, h, mu, sigma, asym) {

	   return( h/2 * cos(pi*asym/2 + (1-asym)*atan((mu-x)/sigma)) / (sigma^2 + (mu-x)^2)^((1-asym)/2) )
}

## =======================================================================
#'linear function
#'
#'@param x numeric vector
#'@param m slope
#'@param c constant value
#'@param mu function position
#'
#'@export
Linear <- function(x, m, c, mu) {

    return( m*x+c )
}

## =======================================================================
## Functions useful for VB top estimation
## =======================================================================

## =======================================================================
#'Eponential decay function
#'
#'@param x numeric vector
#'@param h function amplitude
#'@param mu position of function center
#'@param k decay constant
#'@param c constant value
#'
#'@export
ExpDecay <- function(x, h, mu, k, c) {

    return( h*exp(-k*(mu-x))+c)
}

## =======================================================================
#'Power decay function
#'
#'@param x numeric vector
#'@param h function amplitude
#'@param mu position of function center
#'@param pow value of the power exponent
#'@param c constant value
#'
#'@export
PowerDecay <- function(x, h, mu, pow, c) {

    return( h/((mu-x+1)^pow+1)+c)
}

## =======================================================================
#'Sigmoid function
#'
#'@param x numeric vector
#'@param h function amplitude
#'@param mu position of function center
#'@param k sigmoid decay rate
#'@param c constant value
#'
#'@export
Sigmoid <- function(x, h, mu, k, c) {

    return( h/(1+exp(-k*(x-mu))+c))
}

## =======================================================================
#'HillSigmoid function
#'
#'@param x numeric vector
#'@param h function amplitude
#'@param mu position of function center
#'@param pow sigmoid decay rate
#'@param A Sigmoid upper limit
#'@param B Sigmoid lower limit
#'
#'@export
HillSigmoid <- function(x, h, mu, pow, A, B) {

    return( A+(B-A)*x^pow/(mu^pow+x^pow) )
}

## =======================================================================
#VBtop function needed to store VBtop position
#evaluated either using Linear Fit or Decay Fit
#
VBtop <- function(x, mu) {

   return(x*NA)
}

## =======================================================================
#Derivate function needed to store the Core Line derivate
#
Derivative <- function(x, mu) {

   return(x*NA)
}



## =======================================================================
## ETG
# Empirically Transformed Gaussian function (presa da xcms)
#etg <- function(x, H, t1, tt, k1, kt, lambda1, lambdat, alpha, beta) {
#    2*H*exp(0.5)/((1+lambda1*exp(k1*(t1-x)))^alpha + (1+lambdat*exp(kt*(x-tt)))^beta - 1) }
## =======================================================================
# Empirically Transformed Gaussian function
# Journal of Chromatography A, 952 (2002) 63-70
# Comparison of the capability of peak functions in describing real chromatographic peaks
# Jianwei Li
# http://scholar.google.com/scholar?as_q=Comparison+of+the+capability+of+peak+functions+in+describing+real+chromatographic+peaks.
# t denotes time; substitute with x
# H is related to peak amplitude;
# lambdaL and lambdaT are pre-exponential parameters,
# kL and kT are the parameters related to the speeds of the rise and fall of the leading and trailing edges, respectively;
# xL and xT are the inflection times of the leading and trailing edges, respectively - The selection of xL and xT does not have to be accurate, and can be replaced by the times at half height.

# alpha and beta are the parameters to further modify the shapes of the leading and trailing edges, respectively.

#etg <- function(x, h, xL, xT, kL, kT, lambdaL, lambdaT, alpha, beta) {
#    h/((1+lambdaL*exp(kL*(xL-x)))^alpha + (1+lambdaT*exp(kT*(x-xT)))^beta - 1)
#}
## =======================================================================





