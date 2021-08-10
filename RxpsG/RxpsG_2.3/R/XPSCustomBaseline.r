# definition of CUSTOM baseline functions



#====================================
#    Spline background
#====================================
#'Spline Baseline
#'
#'Definition of Spline Algorythms
#'
#'Computes the Spline curve trough a set of splinePoints
#'
#'@param spectra coreline on which apply spline function
#'@param splinePoints set of points where the spline has to pass through
#'@param limits the abscissa range where the spline is defined
#'
#'@export
#'

baseSpline <- function(spectra, splinePoints, limits) {
	       object <- list(y=spectra[1,]) #spectra contiene solo l'ordinata
          npt <- length(object$y)
          BGND<-spline(x=splinePoints$x, y=splinePoints$y, method="natural", n=npt)
          if (limits$x[1] > limits$x[2]) {
             BGND$y<-rev(BGND$y)   #se ho scala in BE inverto la spline
          }
          baseline <- matrix(data=BGND$y, nrow=1)
#here a list following the baseline package format is returned
          return(list(baseline = baseline, corrected = spectra - baseline, spectra=spectra))
      }


#============================================
# 2 Parameter CrossSection Shirley background
#============================================
# Shirley2P, 2 parameter Shirley background derived from J. Vegh, Surface Science 563, (2004), 183
# This background is based on Aproximated Shirley Cross section function defined by
# SC = Integral[ L(E)* K(E'-E) ]
# The energy loss (background) is then defined by
#
# Shirley(E) = Integral[ L(E)*K(E'-E) * J(E)dE' ]    J(E) = acquired spectrum
#
# which may be expressed by a two paramenter function:
#
# Shirley(E) = Integral[ BsT/(Cs + T^2) * J(E)dT ]  where T=E'-E
#
#'Shirley2P  is a two parameter Bs, Cs expressely evaluated for the Shirley background
#'Shirley Cross Section described by: B*T/(C + T^2)     T=energy loss
#'see J. Vegh, Surface Science 563 (2004) 183,    J. Vegh  J.Electr.Spect.Rel.Phen. 151 (2006) 159
#'
#'@param object XPSSample object
#'@param limits limits of the XPSSample object
#'@return returns the Two_parameter Shirley background.
#'
#'@export
#'


Shirley2P<-function(object, limits){
   LL <- length(object@RegionToFit$x)
#KE scale: reverse vectors to make the routine compatible with BE scale
   if (object@Flags[1] == FALSE){  #KE scale
       object@RegionToFit$x<-rev(object@RegionToFit$x)
   }
   Lim1<-limits$y[1]
   Lim2<-limits$y[2]
   spect<-object@RegionToFit$y - Lim2

   Cs<-750
   shirley2P<-vector(mode="numeric", length=LL)
   lK<-vector(mode="numeric", length=LL)
# Evaluation of the Shirley Cross Section lK() the Bs parameter
   tmp<-0
   dE <-(object@RegionToFit$x[2] - object@RegionToFit$x[1]) # Estep
   for (ii in 1:LL){
       deltaE<-(ii-1)*dE
       num<-dE                    #dE
       denom<-(Cs+deltaE^2)       #(Cs+[E'-E]^2)
       lK[ii]<-num/denom          #dE/{Cs+[E'-E]^2)   the Shirley Cross Section/B
       tmp<-tmp+lK[ii]*spect[ii]       #j(E')dE/{Cs+[E'-E]^2}
   }
   Bs<-(Lim1-Lim2)/tmp

#--- Computation shirley2P background: OK for both KE and BE scale

   for (ii in 1:LL){
       shirley2P[ii]<-Bs*sum(lK[1:(LL-ii+1)]*spect[ii:LL])
   }
   shirley2P<-shirley2P+Lim2    #add the Low-Energy-Side value of the spectrum
   return(shirley2P)
}

#============================================
# 3 Parameter CrossSection Shirley background
#============================================
# 3 parameter Shirley background derived from J. Vegh  J.Electr.Spect.Rel.Phen. 151 (2006) 159
# This background is based on Aproximated Shirley Cross section function defined by
# SC = Integral( L(E)* K(E'-E) )
# The energy loss (background) is then defined by
#
# Shirley(E) = Integral( L(E)*K(E'-E) * J(E)dE'     J(E) = acquired spectrum
#
# which may be expressed by a two paramenter function:
#
# Shirley3p = Integral[ BsT/(Cs + T^2)*(1-exp(-Ds*T)) * J(E)dT ]   where T=E'-E
#
#'Shirley3p is a three parameters Bs=30eV^2, Cs=750eV^2, Ds=0.75 are expressely evaluated for the Shirley background
#'Shirley Cross Section described by: B*T/(C + T^2)*(1-exp(-Ds*T))     T=energy loss
#'see J. Vegh  J.Electr.Spect.Rel.Phen. 151 (2006) 159
#'
#'@param object XPSSample object
#'@param Wgt XPSSample object
#'@param limits limits of the XPSSample object
#'@return returns the 3-parameter Shirley background.
#'
#'@export
#'

Shirley3P<-function(object, Wgt, limits){
   LL <- length(object@RegionToFit$x)
#KE scale: reverse vectors to make the routine compatible with BE scale
   if (object@Flags[1] == FALSE){
       object@RegionToFit$x<-rev(object@RegionToFit$x)
   }
   Lim1<-limits$y[1]
   Lim2<-limits$y[2]
   spect<-object@RegionToFit$y-Lim2

   Cs<-2500
   shirley3P<-vector(mode="numeric", length=LL)
   lK<-vector(mode="numeric", length=LL)
# Evaluation of the Shirley Cross Section lK() the Bs parameter
   tmp<-0
   dE <-(object@RegionToFit$x[2] - object@RegionToFit$x[1]) # Estep
   for (ii in 1:LL){
       deltaE<-(ii-1)*dE
       num<-1-exp(-Wgt*abs(deltaE))     #(1-exp(-Ds*T)*dE
       denom<-(Cs+deltaE^2)           #(Cs+[E'-E]^2)
       lK[ii]<-dE*num/denom           #(1-exp(-Ds*T)*dE/(Cs+[E'-E]^2)   the Shirley Cross Section/B
       tmp<-tmp+lK[ii]*spect[ii]      #j(E')dE/{Cs+[E'-E]^2}
   }
   Bs<-(Lim1-Lim2)/tmp

#--- Computation shirley3P background: OK for both KE and BE scale

   for (ii in 1:LL){
       shirley3P[ii]<-Bs*sum(lK[1:(LL-ii+1)]*spect[ii:LL])
   }
   shirley3P<-shirley3P+Lim2    #add the Low-Energy-Side value of the spectrum
   return(shirley3P)
}


#====================================
#    Polynomial+Shirley background
#====================================
#'LPShirley is a LinearPolynomial * Shirley Baseline see Practical Surface Analysis Briggs and Seah Wiley Ed.
#'
#'Applies in all the cases where the Shirley background crosses the spectrum
#'causing the Shirley algorith to diverge. 
#'In LPShirley firstly a linear background subtraction is performed to recognize presence of multiple peaks
#'Modified Bishop polynom weakens the Shirley cross section
#'Bishop polynom:  PP[ii]<-1-m*ii*abs(dX)   dX = energy step, m=coeff manually selected
#'in LPShirley the classical Shirley expression is multiplied by PP:
#'S(Ei) =  Shirley(E) = Integral( L(E)*K(E'-E)* PP(E)* J(E)dE'     J(E) = acquired spectrum
#'where the integral is extended between Ei and infinite and j=i-1
#'
#'@param object  coreline region where to apply LPshirley function
#'@param m  coefficient to
#'@param limits limits of the XPSSample object
#'@return returns the LPShirley background
#'
#'@export
#'

LPShirley <- function(object, m, limits) {
   LL<-length(object@RegionToFit$y)
   objX <- object@RegionToFit$x
   objY <- object@RegionToFit$y
   if (object@Flags[1] == FALSE){  #KE scale
       objX<-rev(objX)
   }
   Lim1<-limits$y[1]
   Lim2<-limits$y[2]
#   objY<-objY-Lim2

   MinLy<-min(limits$y)
   MaxLy<-max(limits$y)
   dX<-abs(objX[1]-objX[2])  #energy step

	BGND <- rep.int(MinLy, LL)
	SumBGND <- sum(BGND)
	SumRTF <- sum(objY)
	RangeY <- diff(range(limits$y))
   maxit <- 50
   err <- 1e-6

   dY<-(limits$y[2]-limits$y[1])/LL
   LinBkg<-vector(mode="numeric", length=LL)
   for (ii in 1:LL){ #generating a linear Bkg under the RegionToFit
       LinBkg[ii]<-limits$y[1]+(ii-1)*dY
   }
   object@Baseline$x<-objX
   object@Baseline$y<-LinBkg

   plot(object)

#LinBKg subtraction
   objY_L<-objY-LinBkg
   maxY_L<-max(objY_L)
#search for multiple peaks: no derivative used because noise causes uncontrollable derivative oscillation
#are there multiple peaks (spin orbit splitting, oxidized components...)
   idx<-which(objY_L > maxY_L/2)  #regions where peaks are > maxY_L/2
   Lidx<-length(idx)
#idx contains the indexes of the spectral points where the intensity I > maxY/2
#in presence of spin orbit splitting there are two regions satisfying this condition
#in presence of oxidized components the spin orbit leads to 4 peaks: 
#two for the pure element two for the oxidized element.
#In general more than one peak could be present => more than one region where I > maxY/2
#How to identify the number of these regions? 
#idx2[ii] = idx[ii]-idx[ii-1]
#if idx[1]==0 this difference > 1 in correspondence of the beginning of each region:
#
   idx[1]<-0                        #idx==   0 225 226  ...  344 345  655 656 657  ...   943 944   1164 1165 1166 ...
   idx2<-idx[2:Lidx]-idx[1:Lidx-1]  #idx2==  225 1  1  1   ...  1  1  310  1  1  1   ...    1  1   220  1   1   1  ....
                                    #        ----------REG1---------  ++++++++++++REG2++++++++++   ----------REG3---------  ...
                                    #     255-0=255                   655-345==310                 1164-944=220
   Regs<-which(idx2>1)
   Nr<-length(Regs)

#now working on the original non subtracted data
#position  of the last peak if BE scale, of the first peak if KE scale
   PeakMax<-max(objY[idx[Regs[Nr]+1]:LL])    #max of the last region
   posPeak<-which(objY==PeakMax)
#Is there any flat region at the edge of the coreline?
   tailTresh<-(max(objY)-MinLy)/50+MinLy  #define a treshold on the basis of peak and edge intensities
	for (Tpos in LL:1){   #Tpos = beginning of the flat region
       if(objY[Tpos] > tailTresh) { break }
   }
   PP<-vector(mode="numeric", length=LL)
   for (ii in 1:posPeak) {
       PP[ii]<-1-m*ii*abs(dX)  #Bishop weakening polynom see Practical Surface Analysis Briggs and Seah Wiley Ed.
   }
#Bishop polynom modified to better describe the flat region at high BE (lowKE)
#PP decreases linearly until the max of peak at higher BE (lower KE) then linearly
#increases till PP==1 the value assumed in the flat region at the end (beginning) of the
#core line. PP==1 means background == non weakened classical Shirley
   dPP<-(1-PP[posPeak])/abs(Tpos-posPeak)
   for (ii in (posPeak+1):Tpos) {
       PP[ii]<-PP[posPeak]+dPP*(ii-posPeak)
   }
   PP[Tpos:LL]<-1
   for (nloop in seq_len(maxit)) {
		for (ii in 1:LL) {
			BGND[ii] = ((RangeY/(SumRTF-sum(BGND)))*(sum(objY[ii:LL]*PP[ii:LL])-sum(BGND[ii:LL])))+MinLy
		}
		if ( abs( (sum(BGND)-SumBGND)/SumBGND ) < err ) { break }
		SumBGND <- abs(sum(BGND))
	}

## reverse y-baseline in case of reverse bgnd
#	if ( limits$y[which.max(limits$x)] < limits$y[which.min(limits$x)]  ) {
#	baseline <- matrix(data=rev(BGND), nrow=1) }
#	else { baseline <- matrix(data=BGND, nrow=1) }
	baseline <- matrix(data=BGND, nrow=1)
	MinBGND <- min(c(BGND[1], BGND[LL]))
	BGND<-(limits$y[1]-limits$y[2])/(BGND[1]-BGND[LL]) * (BGND-MinBGND)+MinLy   #Adjust BGND on the limits$y values
   return(BGND)
}


#====================================
#    Classic Shirley background
#====================================

IterShirley <- function(object) {   #Iterative Shirley following the Sherwood method
	obj <- object@RegionToFit$y
	limits<-list(x=object@Boundaries$x, y=object@Boundaries$y)
   LL<-length(object@RegionToFit$y)
   if (object@Flags[1] == FALSE){   #KE scale
       object@RegionToFit$x<-rev(object@RegionToFit$x)
   }
   Lim1<-limits$y[1]
   Lim2<-limits$y[2]
   spect<-object@RegionToFit$y-Lim2
   LL<-length(obj)
   Ymin <- min(limits$y)
	BGND <- rep.int(Ymin, LL)
	SumBGND <- sum(BGND)
	SumRTF <- sum(obj)
	RangeY <- diff(range(limits$y))
   maxit=50
   err=1e-6

   for (nloop in seq_len(maxit)) {
		for ( idx in seq_len(LL)) {
			BGND[idx] = ((RangeY/(SumRTF-sum(BGND)))*(sum(obj[idx:LL])-sum(BGND[idx:LL])))+Ymin
		}
		if ( abs( (sum(BGND)-SumBGND)/SumBGND ) < err ) { break }
		SumBGND <- abs(sum(BGND))
	}
   return(BGND)
}


#====================================
#    2 Parameter Tougaard background
#====================================

# 2 parameter Tougaard background for metallic like samples
# This background is based on the Universal Cross section function defined by
# UC = Integral( L(E)* K(E'-E) )
# The energy loss (background) is then defined by
#
# Tougaard(E) = Integral( L(E)*K(E'-E) * J(E)dE'     J(E) = acquired spectrum
#
# which may be expressed by a two paramenter function:
#
# Tougaard2p = Integral( BT/(C + T^2)^2 * J(E)dT   where T=E'-E
#
#'Two parameter B, C Tougaard background
#'Universal Cross Section described by: B*T/(C + T^2)^2     T=energy loss
#'see S. Tougaard, Surf. Sci. (1989), 216, 343;  S. Tougaard, Sol. Stat. Comm.(1987), 61(9), 547
#'
#'@param object XPSSample object
#'@param limits limits of the XPSSample object
#'@return returns the 2parameter Tougaard background.
#'
#'@export
#'
Tougaard2P<-function(object, limits){
   LL <- length(object@RegionToFit$x)
   if (object@Flags[1] == FALSE){  #KE scale
       object@RegionToFit$x<-rev(object@RegionToFit$x)
   }
   avg1<-limits$y[1]
   avg2<-sum(object@RegionToFit$y[(LL-4):LL])/5
   object@RegionToFit$y<-object@RegionToFit$y-avg2

   C<-1643
   tougaard<-vector(mode="numeric", length=LL)
   lK<-vector(mode="numeric", length=LL)
# Evaluation of the Universal Cross Section lK() and the B parameter
   tmp<-0
   dE <-(object@RegionToFit$x[2] - object@RegionToFit$x[1]) # Estep
   for (ii in 1:LL){
       deltaE<-(ii-1)*dE
       num<-deltaE*dE             #(E'-E)*dE
       denom<-(C+deltaE^2)^2      #(C+[E'-E]^2}^2
       lK[ii]<-num/denom          #(E'-E)/{C+[E'-E]^2}^2  the uniiversal Cross Section/B
       tmp<-tmp+lK[ii]*object@RegionToFit$y[ii]       #j(E')(E'-E)/{C+[E'-E]^2}^2
   }
   B<-(avg1-avg2)/tmp

#--- Computation Tougaard background: OK for both KE and BE scale

   for (ii in 1:LL){
       tougaard[ii]<-B*sum(lK[1:(LL-ii+1)]*object@RegionToFit$y[ii:LL])
   }

   tougaard<-tougaard+avg2        #add the Low-Energy-Side value of the spectrum
   return(tougaard)
}


#====================================
#    3 Parameter Tougaard background
#====================================

# 3 parameter Tougaard background for non metallic  matter such as polymers
# This background is based on the Universal Cross section function defined by
# UC = Integral( L(E)* K(E'-E) )
# The energy loss (background) is then defined by
#
# Tougaard(E) = Integral( L(E)*K(E'-E) * J(E)dE'     J(E) = acquired spectrum
#
# which may be expressed by a three paramenter B, C, D function:
#
# Tougaard3p = Integral( BT/[(C + T^2)^2 + D*T^2] * J(E)dT   where T=E'-E
#
#'Tougaard3P is the three parameter B, C and D Tougaard background
#'Universal Cross Section described by: B*T/[(C + T^2)^2 + D*T^2]    T=energy loss
#'see S. Hajati, Surf. Sci. (2006), 600, 3015
#'
#'@param object XPSSample object
#'@param limits limits of the XPSSample object
#'@return returns the 3parameter Tougaard background.
#'
#'@export

Tougaard3P<-function(object, limits){
   LL <- length(object@RegionToFit$x)
   if (object@Flags[1] == FALSE){  #KE scale
       object@RegionToFit$x<-rev(object@RegionToFit$x)
   }
   avg1<-limits$y[1]
   avg2<-sum(object@RegionToFit$y[(LL-4):LL])/5
   object@RegionToFit$y<-object@RegionToFit$y-avg2

   C<-551
   D<-436
   tougaard<-vector(mode="numeric", length=LL)
   lK<-vector(mode="numeric", length=LL)
# Evaluation of the Universal Cross Section lK() and the B parameter
   tmp<-0
   dE <-(object@RegionToFit$x[2] - object@RegionToFit$x[1]) # Estep
   for (ii in 1:LL){
       deltaE<-(ii-1)*dE
       num<-deltaE*dE                      #(E'-E)*dE
       denom<-(C+deltaE^2)^2 + D*deltaE^2  #(C+[E'-E]^2}^2 + D(E'-E)^2
       lK[ii]<-num/denom                   #(E'-E)/{C+[E'-E]^2}^2  the modified uniiversal Cross Section L(E)*K(E,T)  T=E'-E
       tmp<-tmp+lK[ii]*object@RegionToFit$y[ii]       #j(E')(E'-E)/{C+[E'-E]^2}^2
   }
   B<-(avg1-avg2)/tmp

#--- Computation Tougaard background: OK for both KE and BE scale

   for (ii in 1:LL){
       tougaard[ii]<-B*sum(lK[1:(LL-ii+1)]*object@RegionToFit$y[ii:LL])
   }

   tougaard<-tougaard+avg2    #add the Low-Energy-Side value of the spectrum
   return(tougaard)
}



#====================================
#    4 Parameter Tougaard background
#====================================

# 4 parameter Tougaard background for non metallic homogeneous matter such as polymers
# This background is based on the Universal Cross section function defined by
# UC = Integral( L(E)* K(E'-E) )
# The energy loss (background) is then defined by
#
# Tougaard(E) = Integral( L(E)*K(E'-E) * J(E)dE'     J(E) = acquired spectrum
#
# which may be expressed by a three paramenter B, C, C', D function:
#
# Tougaard4p = Integral( BT/[(C + C'*T^2)^2 + D*T^2] * J(E)dT   where T=E'-E
#
#'Tougaard4P is the four parameters B, C, C' and D Tougaard background
#'Universal Cross Section described by: B*T/[(C + C'*T^2)^2 + D*T^2]    T=energy loss
#'see R. Hesse et al., Surf. Interface Anal. (2011), 43, 1514
#'
#'@param object XPSSample object
#'@param C1 numeric is a parameter value
#'@param limits limits of the XPSSample object
#'@return Return the 4parameter Tougaard background.
#'
#'@export

Tougaard4P<-function(object, C1, limits){
   LL <- length(object@RegionToFit$x)
   if (object@Flags[1] == FALSE){  #KE scale
       object@RegionToFit$x<-rev(object@RegionToFit$x)
   }
   avg1<-limits$y[1]
   avg2<-sum(object@RegionToFit$y[(LL-4):LL])/5
   object@RegionToFit$y<-object@RegionToFit$y-avg2

   C<-551
   D<-436
   tougaard<-vector(mode="numeric", length=LL)
   lK<-vector(mode="numeric", length=LL)
# Evaluation of the Universal Cross Section lK() and the B parameter
   tmp<-0
   dE <-(object@RegionToFit$x[2] - object@RegionToFit$x[1]) # Estep
   for (ii in 1:LL){
       deltaE<-(ii-1)*dE
       num<-deltaE*dE                      #(E'-E)*dE
       denom<-(C+C1*deltaE^2)^2 + D*deltaE^2  #(C+[E'-E]^2}^2 + D(E'-E)^2
       lK[ii]<-num/denom                   #(E'-E)/{C+[E'-E]^2}^2  the modified uniiversal Cross Section L(E)*K(E,T)  T=E'-E
       tmp<-tmp+lK[ii]*object@RegionToFit$y[ii]       #j(E')(E'-E)/{C+[E'-E]^2}^2
   }
   B<-(avg1-avg2)/tmp

#--- Computation Tougaard background: OK for both KE and BE scale

   for (ii in 1:LL){
       tougaard[ii]<-B*sum(lK[1:(LL-ii+1)]*object@RegionToFit$y[ii:LL])
   }

   tougaard<-tougaard+avg2    #add the Low-Energy-Side value of the spectrum
   return(tougaard)
}

