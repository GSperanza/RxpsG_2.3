### $Id: readScienta.R,v 0.0.1
###
###             Read pxt files data from XPS Scienta
###

#'read.scienta
#'
#'Read Scienta format
#'
#'Read Scienta format
#'
#'@param file File
#'@param debug debug switch
#'
#'@export
#'

read.scienta <- function(file=NULL, debug = FALSE) {

#funzioni per leggere i vari tipi di dati: char, short, fload, double.
## i char sono separati dal carattere space == 20 oppure da 00 ????
## per leggere i character meglio usare readChar
        get_char <- function(con, POS=NULL, n=1) {
            if ( ! is.null(POS) ) seek(con, where=POS, "start")
            return(readBin(con, "char", n))
        	}

        get_short  <- function(con, POS=NULL) {
            if ( ! is.null(POS) ) seek(con, where=POS, "start")
            return(readBin(con, "int", 1, 2, signed=FALSE))
        	}

        get_float  <- function(con, POS=NULL) {
            if ( ! is.null(POS) ) seek(con, where=POS, "start") #
            return(readBin(con, "double", 1, 4))
        	}

        get_double  <- function(con, POS=NULL) {
            if ( ! is.null(POS) ) seek(con, where=POS, "start") #
            return(readBin(con, "double", 1, 8))
        	}

        check_pxt_version <- function(con) {
        	version <- try( get_short(con), silent = TRUE )
            return( ifelse(identical(class(version), "try-error") || length(version) == 0, 0, 1))
        	}
        
#--- getInfo for XPSCoreLine
# returns the value from info slot of XPSCoreLine
		getInfo <- function(object, label, type=c("numeric","character")) {
			type <- match.arg(type)
			idx <- charmatch(label, object)
			out <- unlist(strsplit(object[idx], "\\="))[2]
			if (type == "numeric") out <- as.numeric(out)
			return(out)
		}


#--- read.pxt: da file Ses.h
# define MAX_SES_NAME2 40
# Maximum length of wave name in version 1 and 2 files. Does not include the trailing null.
		
		#typedef struct SesHeader {
		#	short version;						// Version number for backwards compatibility.
		#	short hole1;						// 
		#	short offset_EndSES;				// .
		#	char  hole2[78];					// .
		#	short npts;						    // .
		#	char  hole3[14];				    // .
		#	char  curve_name[MAX_SES_NAME2];	// .
		#	short npts2;		                // .
		#	char  hole4[14];					// .
		#   double step_E;                      // +156 = Step Energy.
		#   double d1;
		#    double d2;
		#    double d3;
		#    double x_1;                         // +188 = first X value.
		#    char  hole5[196];
		#}sesHeader;
		
		#struct SesCurve {
		#    sesHeader sesH;                     // header 
		#    double *data;                       // point values
		#    char *labels;                       // labels terminated width "EndSES"
		#}sesCurve[1024];	
		#define SES_BYTES   (sizeof(sesHeader))

		read.pxt <- function(con) {
			# struct SesHeader
			unused <- readBin(fp, "raw", 2)
			offset_EndSES <- get_short(fp)#,4)
			labelpos <- seek(fp, where=NA) ## offset dove siamo x ENDSES
			unused <- readBin(fp, "raw", 78) #[78] char
			npts <- get_short(fp) # ,84) ## npts
			unused <- readBin(fp, "raw", 14) #[14] char
			curve_name <- readChar(fp, 40) # con questa leggo comunque 40 char e non ho bisogno di seek
			npts2 <- get_short(fp) #,140) ## npts
			unused <- readBin(fp, "raw", 14) #[14] char
			step_E <- get_double(fp) #,156)
			d1 <- get_double(fp)  # ,164)
			d2 <- get_double(fp)  # ,172)
			d3 <- get_double(fp)  # ,180)
			x_1 <- get_double(fp) # ,188)
			unused <- readBin(fp, "raw", 196) #[196] char

			if (debug) {
				#cat("\nversion = ",version)
				cat("\nNPT = ",npts)
				#cat("\nNPT2 = ",npts2)
				cat("\noffset_EndSES = ",offset_EndSES)
				cat("\nEndSES = ",(labelpos+offset_EndSES+1))
				cat("\ncurve_name = ",curve_name)
				cat("\nSteps = ", step_E)
				cat("\nd1 = ",d1)
				cat("\nd2 = ",d2)
				cat("\nd3 = ",d3)
				cat("\nx_1 = ",x_1)
				cat("\n")
			}
			## posizionamento inizio dati
			if (debug) cat(" Posizione partenza dati",seek(fp, where=NA),"\n")
			#seek(fp, where=196, "current") 
			vect <- readBin(fp, "double", npts, 8) 

			startSESlabel <- seek(fp, where=NA)
			if (debug) cat(" Posizione partenza sezione SES ",startSESlabel,"\n")
			endSESlabel <- (labelpos+offset_EndSES+1)
			SESlabellength <- endSESlabel-startSESlabel+1

			if (debug) cat(" lunghezza sezione Label",SESlabellength,"\n")
			## read labels terminated width "EndSES"	
			lab <- readChar(fp, SESlabellength)
			m <- gregexpr("\r", lab)

			pxtlabel <- regmatches(lab, m, invert=TRUE)
			if (debug) cat("")
			#print(pxtlabel)

			return(list(data=vect, info=unlist(pxtlabel)))
		}
		
## ================================================
# End of utility functions
## ================================================



#--- Start of main code

	if (is.null(file)) file <- file.choose()
	fp <- file(file, open="rb")
	on.exit(close(fp))
	fName <- basename(file)   #this is the filename without the pathname
	idx <- 0
	ans <- new("XPSSample")
	# Version number for backwards compatibility
	while (	check_pxt_version(fp) ) {
			idx <- idx + 1
         SpectInfo <- NULL  #scienta spetrum info
         T<-NULL
			tmp <- read.pxt(fp)
			
			symbol <- unlist(strsplit(getInfo(tmp$info, "Region Name", "char"), "\\#"))[1]

			# Kinetic Energy or Binding Energy
			X <- seq(from=getInfo(tmp$info, "Low Energy"),
						to=getInfo(tmp$info, "High Energy"),
						by=getInfo(tmp$info, "Energy Step"))

			# Acquisizione come Binding?
			binding <- ifelse( getInfo(tmp$info, "Energy Scale", "char") == "Binding", TRUE, FALSE)
			if (binding) { 
				X <- getInfo(tmp$info, "Excitation Energy") - X  #convert X from kinetic to binding energy
				units <- c("Binding Energy [eV]","Intensity [cps]")
				}
			else {
				units <- c("Kinetic Energy [eV]","Intensity [cps]")
			}
			## (y scale) Intensity [counts] vs. Counts per Second (cps) conversion	
			## ALWAYS TRUE		
			nsweep <- getInfo(tmp$info, "Number of Sweeps")
			acq_time <- getInfo(tmp$info, "Step Time")
			## acq_time could be ms
			if ( acq_time > 90 ) acq_time <- acq_time/1000  #time in seconds
			steps <- getInfo(tmp$info, "Energy Step")
			passEnergy <- getInfo(tmp$info, "Pass Energy")
			
			## Conversion factor
			ConvFactor <- abs( nsweep * acq_time * passEnergy / steps / 10.33 )
			## Convert the intensity in CPS
			Y <- tmp$data / ConvFactor
			
			## 18/9/2012 Add isScienta Flags[3] == TRUE for pxt files
         LL1<-length(X)
         LL2<-length(Y)
         if (LL1>LL2) {
            X<-X[-(LL2+1:LL1)] #ascisse X piu' lunghe delle ordinate Y: elimino dati X in piu'
         }
         if (LL2>LL1) {
            Y<-Y[-(LL1+1:LL2)] #ascisse X piu' corte delle ordinate Y:  elimino  dati Y in piu' '
         }
         LL1<-length(X)
         LL2<-length(Y)
         T<-rep(1, LL2) #data for the analyzer transfer function == unity for the Scienta Analyzers
			TotTime<-LL2*acq_time*nsweep
         SpectInfo[1]<-paste("   XPS     Spectrum    Lens Mode:default   Resolution:Pass energy ", as.character(passEnergy),"   Iris(Aper):manual", sep="")
         SpectInfo[2]<-paste("   Acqn. Time(s): ", as.character(TotTime),"  Sweeps: ", as.character(nsweep),"   Anode:Mono Al Ka=1486.6eV   Step(meV): ",as.character(steps), sep="")
         SpectInfo[3]<-paste("   Dwell Time(ms):", as.character(acq_time*1000),"  Charge Neutraliser :manual   Acquired On :unspecified")
			SEScurve <- new("XPSCoreLine",
							.Data = list(x = X, y = Y, t=T, dummy=NULL),   #X, Y, T=Analizer, transfer function, and dummy column for future use
							units = units,
							Flags = c(binding, TRUE, TRUE),
							Info = SpectInfo,
							Symbol = symbol
						)
			ans[[idx]] <- SEScurve
	}
	slot(ans, "Sample") <- getInfo(ans[[1]]@Info, "Sample", "char")
	slot(ans, "Comments") <- getInfo(ans[[1]]@Info, "Comments", "char")
	slot(ans, "User") <- getInfo(ans[[1]]@Info, "User", "char")
	slot(ans, "Filename") <- fName
	slot(ans, "Project") <- getInfo(ans[[1]]@Info, "Location", "char") 
	names(ans) <- sapply(ans, function(x) x@Symbol)
	return(ans)
}
