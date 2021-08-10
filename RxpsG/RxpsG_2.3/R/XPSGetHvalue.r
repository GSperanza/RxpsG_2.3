# This function computes the which f(y) value corresponds to a y obtained by mouse coord in
# a graphic window. This routine considers the function type following the algorithms described
# in XPSFitAlgorithm.r
#

#'In interactive graphic panels the ordinate of the cursor positin is used to
#'compute the intensity of the selected fit function to add or to move.
#'see analysis or Move Component options.
#'
#'@param Object the fitted CoreLine;
#'@param Comp the fitting component which is added or moved;
#'@param FitFunct thee fitting function corresponding to the component;
#'@param A amplitude of the fitting component obtained from the mouse locator
#'@return Returns the component intensity.
#'
#'@export
#'

#for complex functions computes the amplitue to match the component height modified
#with mouse in XPSMoveComponent() or XPSProcessing() in interactive graphical windows
#Essentially here we estimate which value would be returned by a given function
#We estimate the function value if A were equal to 1 as see if the function gives expected values
#In some cases the amplitude has to be normalized (Voigt, AsymmGaussVoigt, DoniachSunjic...)

GetHvalue<-function(Object,Comp, FitFunct, A){
         Param<-Object@Components[[Comp]]@param
         h<-NULL
       	switch(FitFunct,
	       "Gauss" ={
#Gauss = A*exp(-(2.77258872*((x-mu)/sigma)^2))   #if A==1 the
                    h<-A      #returned value is simply the amplitude of Gauss
	       },
	       "Lorentz" ={
#Lorentz =  A/(1+(4*(((x-mu)/sigma)^2)))
                    h<-A      #returned value is simply the amplitude of Lorentz
	       },
	       "Voigt" ={
#xx <- (x - mu)/sigma/1.41421356
#yy <- lg/sigma/1.41421356
#Voigt =  A * Re(wofz( complex(real = xx, imaginary = yy) ))/sigma/2.50662827
#If A==1 which value is returned by Voigt function? Computes the reciprocal of this value
                    import::here(wofz, .from=NORMT3)
                    sigma<-Param[[3,1]] #sigma is stored in Param[[3,1]]
                    lg<-Param[[4,1]]
                    yy <- lg/sigma/1.41421356
                    h<-A*sigma*2.50662827/Re(wofz( complex(real = 0, imaginary = yy) ))
	       },
	       "Sech2" ={
#Sech2 = A/( cosh((x-mu)/sigma)*cosh((x-mu)/sigma) )  Se A==1 quale valore mi ritorna la Sech2
                    h<-A       #returned value is simply the amplitude of Sech2
	       },
	       "GaussLorentzProd" ={
#GLProd = A / (1+(4*mix*((x-mu)/sigma)^2)) * exp(-2.77258872*(1-mix)*((x-mu)/sigma)^2) )
                    h<-A
	       },
	       "GaussLorentzSum" ={
#GLSum = A*(mix*exp(-2.77258872*((x-mu)/sigma)^2) + (1-mix)*(1/(1+4*(((x-mu)/sigma)^2))))
                    h<-A
	       },
	       "AsymmGauss" ={
#SimpAsG =  Gauss(x, A, mu, sigma) + (A - Gauss(x, A, mu, sigma))*Tail(x, mu, asym)
#per x=mu Tail=1, AsymmGauss = Gauss
                    h<-A
	       },
	       "AsymmVoigt" ={
#xx <- (x - mu)/sigma/1.41421356
#yy <- lg/sigma/1.41421356
#Voigt =  A * Re(wofz( complex(real = xx, imaginary = yy) ))/sigma/2.50662827   Se A==1 quale valore mi ritorna la Voigt
                    import::here(wofz, .from=NORMT3)
                    sigma<-Param[[3,1]] #sigma is stored in Param[[3,1]]
                    lg<-Param[[4,1]]
                    yy <- lg/sigma/1.41421356
                    h<-A*sigma*2.50662827/Re(wofz( complex(real = 0, imaginary = yy) ))
	       },

	       "AsymmGaussLorentz" ={
#AsymGLSum =  GaussLorentzSum(x, A, mu, sigma, mix) + (A - GaussLorentzSum(x, A, mu, sigma, mix))*Tail(x, mu, asym) )
#per x=mu Tail=1, AsymGLSum= GLSum
#                    mix<-Param[[5,1]]
#                    h<-(A-1+mix)/mix
                    h<-A
	       },
        "AsymmGaussLorentzProd" ={
# z <- (x-mu)/(sigma + asym * (x-mu))
# AsGLProd = A / (1+(4*mix*(z)^2)) * exp(-2.77258872*(1-mix)*(z)^2)
                    h<-A
        },
	       "AsymmGaussVoigt" ={
#AsymmGaussVoigt = A * (gv / (sigma^2 + (x - mu)^2)^((1-asym)/2)) + (1-gv)*Voigt(x, A, mu, sigma, lg)
#AsymmGaussVoigt = A * (  (gv / (sigma^2 + (x - mu)^2)^((1-asym)/2)) + (1-gv) * Re(wofz( complex(real = 0, imaginary = yy) ))/sigma/2.50662827 )

                    import::here(wofz, .from=NORMT3)
                    mu<-Param[[2,1]]
                    sigma<-Param[[3,1]]
                    lg<-Param[[4,1]]
                    asym<-Param[[5,1]]
                    gv<-Param[[6,1]]
                    yy <- lg/sigma/1.41421356

                    h<-A/( (gv / (sigma^2)^((1-asym)/2)) + (1-gv)*Re(wofz( complex(real = 0, imaginary = yy) ))/sigma/2.50662827)
	       },
        "DoniachSunjic" = {
#DS = A/4*(gamma(1-asym)/((mu-x)^2+(sigmaDS/2)^2)^((1-asym)/2)*cos((pi*asym)/2+(1-asym)*atan((mu-x)*2/sigmaDS)))  #DoniachSunjic
                    sigmaDS<-Param[[3,1]]
                    asym<-Param[[4,1]]
 	                  h <- A*4/(gamma(1-asym)/((sigmaDS/2)^2)^((1-asym)/2)*cos(pi*asym/2) )  #DoniachSunjic per x=mu
        },
        "DoniachSunjicTail" = {
#DS = A/4*(gamma(1-asym)/((mu-x)^2+(sigmaDS/2)^2)^((1-asym)/2)*cos((pi*asym)/2+(1-asym)*atan((mu-x)*2/sigmaDS)))  #DoniachSunjic
                    sigmaDS<-Param[[3,1]]
                    asym<-Param[[4,1]]
 	                  h <- A*4/(gamma(1-asym)/((sigmaDS/2)^2)^((1-asym)/2)*cos(pi*asym/2) )  #DoniachSunjic per x=mu
        },
	       "DoniachSunjicGauss" ={

#DS = A/2*(gamma(1-asym)/((mu-x)^2+(sigmaDS/2)^2)^((1-asym)/2)*cos((pi*asym)/2+(1-asym)*atan((mu-x)*2/sigmaDS)))  #DoniachSunjic
#DSmax<-max(DS)
#GS <- DSmax*exp(-(2.77258872*(0.5*(x-mu)/sigmaG)^2))      #Gauss function
#DoniachSunjicGauss = 0.5*(DS  + Gs)* Tail(x, mu, tail) )  #Sum of Gaussian broadening Gs and multiply for a damping factor for low-BE tail
#per x=mu Tail=1, Gauss=1, DoniachSunjicGauss = DS
                    sigmaDS<-Param[[3,1]]
                    asym<-Param[[5,1]]
 	                  h<-A*4/( (gamma(1-asym)/((sigmaDS/2)^2)^((1-asym)/2) ) * cos((pi*asym)/2) )   #DoniachSunjic per x=mu
	       },
	       "DoniachSunjicGaussTail" ={

#DS = A/2*(gamma(1-asym)/((mu-x)^2+(sigmaDS/2)^2)^((1-asym)/2)*cos((pi*asym)/2+(1-asym)*atan((mu-x)*2/sigmaDS)))  #DoniachSunjic
#DSmax<-max(DS)
#GS <- DSmax*exp(-(2.77258872*(0.5*(x-mu)/sigmaG)^2))      #Gauss function
#DoniachSunjicGauss = 0.5*(DS  + Gs)* Tail(x, mu, tail) )  #Sum of Gaussian broadening Gs and multiply for a damping factor for low-BE tail
#per x=mu Tail=1, Gauss=1, DoniachSunjicGauss = DS
                    sigmaDS<-Param[[3,1]]
                    asym<-Param[[5,1]]
 	                 h<-A*4/( (gamma(1-asym)/((sigmaDS/2)^2)^((1-asym)/2) ) * cos((pi*asym)/2) )   #DoniachSunjic for x=mu
	       },
	       "SimplfiedDoniachSunjic" ={
#SimplfiedDoniachSunjic =  A/2 * cos(pi*asym/2 + (1-asym)*atan((mu-x)/sigma)) / (sigma^2 + (mu-x)^2)^((1-asym)/2)
                    sigma<-Param[[3,1]]
                    asym<-Param[[4,1]]
                    h<-A*2*(sigma^2)^((1-asym)/2) /cos(pi*asym/2)
	       },

#---- Definitions for VBtop macro
	       "ExpDecay" ={
#ExpDecay = A*exp(-k*(mu-x))+c
                    h<-A
	       },
	       "PowerDecay" ={
#PowerDecay = A/((mu-x)^pow+0.1)+c
#if x=mu e A==1 mi viene ritornato
                    h<-A
	       },
	       "Sigmoid" ={
#ExpDecay = A/(1+exp-k*(mu-x))+c
                    h<-A*2
	       },


        )
return(h)
}
