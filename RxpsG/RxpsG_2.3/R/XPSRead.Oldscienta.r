#'Read XPS data from old-DOS-Scienta format
#'
#'Reads XPS data 'old-Scienta' format Data.
#'
#'
#'@param filename file name or list of file names
#'@param project  comment
#'@param experiment  comment
#'@param out output class type default list()
#'@return list object class
#'@export
#'@examples
#'
#'\dontrun{
#'pk1 <- read.xps("/Volumes/GoogleDrive/Il mio Drive/R/OldScientaFiles/Specific/DLC.48")
#'}
#'
XPSRead.Oldscienta <- function(filename  = NULL,
					                project    = "",
				                   experiment = "",
					                out        = c("list")) {


#----------------------------------------------------------------------
# Read information in Specific folder
#----------------------------------------------------------------------
  readXPS.Specific <- function(filename, debug){
     tmp <- list()
     tmp$Filename <-  basename(filename)

     fp <- file(filename, open="rb")
	  #on.exit(close(fp))

	  seek(fp, where = 459, "start")
	  tmp$nr_curve <- readBin(fp, "integer", n=1, size=2, signed = FALSE)
     NSpectra<<-tmp$nr_curve
# check_spc_header file format

	  seek(fp, where = 0, "start") # skip first byte
	  buffer <- readBin(fp, "integer", n=1, size=2, signed = FALSE)
	  if (buffer != 256) warning("Wrong SPECIFIC file format.\n")
# read labels
	  seek(fp, where = 0x0D, "current") ## skip 13 bytes
	  label <- list()
     for (idx in c(1:4)) {
		   npt 	<- readBin(fp, "integer", n=1, size=1)
		   if (npt != 0) {
		      txt <- readChar(fp, nchars=npt)
		      skip <- (40 - npt)
		   }  else {
		      txt <- readChar(fp, nchars=40)
		      skip <- 0
		   }

		   if ( nchar(tmp) > 0 ) {
		      label[[idx]] <- txt
		   }
		   seek(fp, where = skip, "current")
     }
     tmp$comments <- label
     seek(fp, where = 189, "start")
	  year  <- readBin(fp, "integer", n=1, size=2)
	  month <- readBin(fp, "integer", n=1, size=2)
	  day   <- readBin(fp, "integer", n=1, size=2)
	  hour  <- readBin(fp, "integer", n=1, size=2)
	  min   <- readBin(fp, "integer", n=1, size=2)
	  tmp$Date <- paste(year, "/", month, "/", day, "-", hour, ":", min, sep="")

	  curve <- list()
	  for (idx in c(1:tmp$nr_curve)) {
		  npt 	<- readBin(fp, "integer", n=1, size=1)
		  curve[[idx]] <- readChar(fp, nchars=npt)
		  SpectNames[idx]<<-curve[[idx]]
		  if (SpectNames[idx]=="SURVEY" || SpectNames[idx]=="SUR"){SpectNames[idx]<<-"survey"}
		  seek(fp, where = (12 - npt), "current")

	  }
	  tmp$curves <- curve

# now close file
     close(fp)
# save information into the XPSSample
     OldScientaObj@Project<<-paste("Acquisition date: ", tmp$Date, sep=" ")
     OldScientaObj@Sample<<-filename
     OldScientaObj@Comments<<-c(tmp$comments[[1]], tmp$comments[[2]])
     OldScientaObj@User<<-"Unknown"
     OldScientaObj@Filename<<-tmp$Filename

#cat(str(tmp)) ## DEBUG
#cat("\n ******\n", str(OldScientaObj)) ## DEBUG
#scan(n=1, quiet=TRUE)
    }

  
#----------------------------------------------------------------------
# Read information in Regions folder
#----------------------------------------------------------------------
  readXPS.Regions <- function(filename, nr_curve){
     fp <- file(filename, open="rb")
     #on.exit(close(fp))

     tmp <- list()
     for (idx in c(1:nr_curve)) {
        SpectInfo<-NULL
        buffer <- readBin(fp, "integer", n=1, size=2, signed = FALSE)
        if (buffer != 256) warning("Wrong REGIONS file format.\n")

        seek(fp, where = 177, "current")

        tmp$mode 	<- readBin(fp, "integer", n=1, size=1)
        tmp$pass_energy <- readBin(fp, "integer", n=1, size=1)
        if (tmp$pass_energy==5) {tmp$pass_energy<-1000}
        if (tmp$pass_energy==4) {tmp$pass_energy<-500}
        if (tmp$pass_energy==3) {tmp$pass_energy<-300}
        if (tmp$pass_energy==2) {tmp$pass_energy<-150}
        if (tmp$pass_energy==1) {tmp$pass_energy<-75}

        seek(fp, where = 54, "current") ## skip 54 bytes
        tmp$mode_type <- readBin(fp, "integer", n=1, size=1)
        tmp$xmin      <- readBin(fp, "double", n=1, size=8)
        tmp$xmax      <- readBin(fp, "double", n=1, size=8)
        tmp$xstep     <- readBin(fp, "double", n=1, size=8)
        tmp$t_acq     <- readBin(fp, "double", n=1, size=8)
        tmp$sweeps    <- readBin(fp, "integer", n=1, size=2)
        tmp$half_range <- readBin(fp, "double", n=1, size=8)
        tmp$t_acq_2    <- readBin(fp, "double", n=1, size=8)
        seek(fp, where = 8, "current")

# generate the XPSCoreline structure
 	    		OldScientaObj[[idx]]<<- new("XPSCoreLine")

# save information into the XPSSample Coreline
        if(tmp$xmin > tmp$xmax) { #Acquisition made using BE scale 
           OldScientaObj[[idx]]@units<<-c("Binding Energy [eV]","Intensity [cps]")
           OldScientaObj[[idx]]@Flags<<-c(TRUE, FALSE, TRUE, FALSE)    #BE scale, CPS, ScientaData, TransmissionCorrection
           sgn <- -1
        } else if(tmp$xmin < tmp$xmax) { #Acquisition made using KE scale
           OldScientaObj[[idx]]@units<<-c("Kinetic Energy [eV]","Intensity [cps]")
           OldScientaObj[[idx]]@Flags<<-c(FALSE, FALSE, TRUE, FALSE)
           sgn <- 1
        }
        OldScientaObj[[idx]]@.Data[[1]]<<-seq(from=tmp$xmin, to=tmp$xmax, by=(sgn*tmp$xstep)) #X axis has to be genrated from range and  Xstep
        TotTime<-tmp$t_acq*(abs(tmp$xmax-tmp$xmin)/tmp$xstep)*tmp$sweeps
        SpectInfo[1]<-paste("   XPS     Spectrum    Lens Mode:",tmp$mode,"   Resolution:Pass energy ", as.character(tmp$pass_energy),"   Iris(Aper):manual", sep="")
        SpectInfo[2]<-paste("   Acqn. Time(s): ", as.character(TotTime),"  Sweeps: ", as.character(tmp$sweeps),"   Anode:Mono Al Ka=1486.6eV   Step(meV): ",as.character(tmp$xstep), sep="")
        SpectInfo[3]<-paste("   Dwell Time(ms):", as.character(tmp$t_acq*1000),"  Charge Neutraliser :manual   Acquired On :unspecified")

        OldScientaObj[[idx]]@Boundaries[[1]]<<-c(tmp$xmin,tmp$xmax)
        OldScientaObj[[idx]]@Info<<-SpectInfo
        OldScientaObj[[idx]]@Symbol<<-SpectNames[idx]
#save acquisition info to change COUNTS in CPS
        Tint[idx]<<-tmp$t_acq
        Nsweeps[idx]<<-tmp$sweeps
        PE[idx]<<-tmp$pass_energy
        Estep[idx]<<-tmp$xstep
      }
      OldScientaObj@names<<-SpectNames #All regions defined, now is possible to save SpectNames
      close(fp) # now close file

#cat(str(tmp)) ## DEBUG
#cat("\n ******\n", str(OldScientaObj)) ## DEBUG
#scan(n=1, quiet=TRUE)
   }

#----------------------------------------------------------------------
# Read information in Results folder
#----------------------------------------------------------------------
  readXPS.Results <- function(filename, nr_curve){
     fp <- file(filename, open="rb")
     #on.exit(close(fp))

     ymin <- NULL
     ymax <- NULL
     tmp <- list()
     for (idx in c(1:nr_curve)) {

        buffer <- readBin(fp, "integer", n=1, size=2, signed = FALSE)
        if (buffer != 256) warning("Wrong REGIONS file format.\n")

        seek(fp, where = 180, "current")
        tmp$n_182     <- readBin(fp, "double", n=1, size=8)
        tmp$npt       <- readBin(fp, "integer", n=1, size=2)
        tmp$half_range <- readBin(fp, "double", n=1, size=8)
        tmp$xrange    <- readBin(fp, "double", n=1, size=8)
# cat(seek(fp, where = NA),"\n")
        tmp$data      <- readBin(fp, "integer", n=tmp$npt, size=4)
        seek(fp, where = (128 - ((4*tmp$npt + 208)%%128)), "current")

# save information into the XPSSample Coreline
        if ( Tint[idx] > 99 ) Tint[idx] <- Tint[idx]/1000  #acq time could be expressed in msec.
        tmp$data<-tmp$data/abs( Nsweeps[idx] * Tint[idx] * PE[idx] / Estep[idx] / 10.33) #counts to CPS conversion and normalization for different PE
        ymin<-min(tmp$data)
        ymax<-max(tmp$data)
        OldScientaObj[[idx]]@Boundaries[[2]]<<-c(ymin,ymax)
        OldScientaObj[[idx]]@.Data[[2]]<<-tmp$data
        LL <- length(OldScientaObj[[idx]]@.Data[[2]])
        OldScientaObj[[idx]]@.Data[[3]]<<-rep(1, LL) #data for the analyzer transfer function == unity for the Scienta Analyzers 
     }
#cat(str(tmp)) ## DEBUG
#cat("\n ******\n", str(OldScientaObj)) ## DEBUG
#scan(n=1, quiet=TRUE)
     close(fp) # now close file
  }


#----------------------------------------------------------------------
# Main Variables
#----------------------------------------------------------------------
  NSpectra<-NULL
  SpectNames<-NULL
  Tint<-NULL
  Nsweeps<-NULL
  PE<-NULL
  Estep<-NULL
  OldScientaObj <- new("XPSSample") #, Project="", Sample="", Comments="", User="", Filename="")
#---
  if ( is.null(filename) ) stop("No Filename null?\n")
  FName <-  basename(filename)
  rootdir <- dirname(dirname(filename))
cat("DIR:", rootdir,"\n")
cat("FILE:", FName,"\n")

# read Specific
cat("\n Read Specific")
  readXPS.Specific(paste(rootdir, "SPECIFIC", FName, sep="/"))

# read Regions
cat("\n Read Regions")
	 readXPS.Regions(paste(rootdir, "REGIONS", FName, sep="/"), NSpectra)

# read Results
cat("\n Read Results")
  res <- readXPS.Results(paste(rootdir, "RESULTS", FName, sep="/"), NSpectra)

  cat("\n ==> ",OldScientaObj@Sample, OldScientaObj@Project)
  cat("\n Sample: ",OldScientaObj@Comments)
  txt<-c("Element", "Xmin", "Xmax", "Ymin", "Ymax", "N. sweeps", "t int.")
  txt<-format(txt, digits=2, justify="right", width=10)
  cat("\n", txt)
  for (idx in 1:NSpectra) {
       xmin<-min(OldScientaObj[[idx]]@Boundaries[[1]])
       xmax<-max(OldScientaObj[[idx]]@Boundaries[[1]])
       OldScientaObj[[idx]]@Boundaries[[2]]<-round(OldScientaObj[[idx]]@Boundaries[[2]], 2)
       txt<-c(OldScientaObj[[idx]]@Symbol, xmin, xmax, OldScientaObj[[idx]]@Boundaries[[2]], Nsweeps[idx], Tint[idx])
       txt<-format(txt, justify="right", width=10)
       cat("\n", txt)
  }
  return(OldScientaObj)
}
