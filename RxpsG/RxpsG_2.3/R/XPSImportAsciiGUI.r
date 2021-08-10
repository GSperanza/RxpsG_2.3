#'XPSImport.Ascii function
#'
#'Function allowing import of textual (ascii) file data
#'Options are available to account for header and data separator
#'Options are also provided to store data in a XPS-Sample data-structure.
#'No parameters are passed to this function
#'
#'@return Returns the processed \code{object}.
#'
#'\dontrun{
#'XPSImport.Ascii()
#'}
#'
#'@export
#'

XPSImport.Ascii <- function() {
    options(guiToolkit = "tcltk")

#--- warnings if required selections are lacking
	   check_selection <- function(){
         if (nchar(svalue(ColX)) * nchar(svalue(ColY))==0) {  #X and Y columns must be indicated
              gmessage(msg="Please give ColX and ColY to be imported", title="WARNING: column lacking", icon="warning")
              return(FALSE)
         }
         if ( ! (svalue(reverseX) || svalue(NOreverseX))) {   #one of reverseX or NOreverseX should be SELECTED
              gmessage(msg="Reverse X axis? Plsease check!", title="WARNING X axis!", icon="warning")
              return(FALSE)
         }
         if ( nchar(svalue(CLname))==0) {   #Column name must be indicated
              gmessage(msg="Core Line Name please!", title="WARNING Core Line Name", icon="warning")
              return(FALSE)
         }
         return(TRUE)
      }

#--- read data from file
      update_output <- function(...) {
          opt <- svalue(OptLayout)
          Nrws <- as.numeric(svalue(NRowHeader))
          scf <- svalue(scanFile)
          if (scf==FALSE){
             out <- read.table(file=FNameIN,sep=seps[opt$Separator],dec=decs[opt$Decimal],
                          skip=Nrws, colClasses="numeric" )
          } else {
             fp <- file(FNameIN, open="r")
             Ncol <- as.numeric(svalue(dataCol))
             tmp <- NULL
             out <- NULL   #ora leggo i dati
             tmp <- scan(fp, what="character",nlines=1, skip=(Nrws-1), quiet=TRUE) #skip header
             while (length(tmp)>0) {
                 tmp <- scan(fp, what="character", n=Ncol, quiet=TRUE)
                 tmp <- sub(", ", "  ", tmp)   #changes separation "," with " ": for data  1, 2,143, 5,723  generates  1  2,143  5,723
                 tmp <- sub(",", ".", tmp)     #changes decimal "," with ".": for data  1  2,143  5,723  generates  1  2.143  5.723
                 if (is.na(as.numeric(tmp))) break #stop reading if there are characters which cannot translated in numbers
                 out <-  rbind(out,as.numeric(tmp))
             }

          }
          output[] <- out
          invisible(out)
      }

#--- Add a new XY data in a New CoreLine in an existing XPSSample
      addCoreLine <- function(){

              Xidx <- as.numeric(svalue(ColX))
              Yidx <- as.numeric(svalue(ColY))
              DF <- update_output()
              LL <- length(DF[[Xidx]])
              if (svalue(reverseX) && DF[[Xidx]][1] < DF[[Xidx]][LL]) { #reverse X axis selected but X is ascending ordered
                 answ <- gconfirm(msg="X is in ascending order. Do you want to reverse X axis? ", title="CONFIRM REVERSE AXIS", AICON="WARNING")
                 if (answ == TRUE ){
                    DF[[Xidx]] <- rev(DF[[Xidx]]) #reverse X in descending order
                    DF[[Yidx]] <- rev(DF[[Yidx]]) #reverse X in descending order
                 } else {
                    svalue(reverseX) <- FALSE
                    svalue(NOreverseX) <- TRUE
                 }
              }
#	             FName	<- get(activeFName, envir = .GlobalEnv)
	             LL <- length(XPSSample)+1
 	            mynewcoreline <- new("XPSCoreLine",
				                .Data = list(x = DF[[Xidx]], y = DF[[Yidx]], t=NULL, err=NULL),   #err is dedicated to standard errors on Y data
				                units = c(svalue(unitX), svalue(unitY)),
				                Flags = c(svalue(reverseX), TRUE, FALSE, FALSE),
				                Symbol= svalue(CLname)
                   )
	             CLnames <- names(XPSSample)
	             XPSSample[[LL]] <<- mynewcoreline
	             names(XPSSample) <<- c(CLnames, as.character(svalue(CLname)))
              plot(XPSSample)
       }

#--- Add a new XY data in a New CoreLine in an existing XPSSample
      addErrors <- function(){
              Xidx <- as.numeric(svalue(ColX))
              Yidx <- as.numeric(svalue(ColY))
              DF <- update_output()
              LL <- length(DF[[Xidx]])
              if (svalue(reverseX) && DF[[Xidx]][1] < DF[[Xidx]][LL]) { #reverse X axis selected but X is ascending ordered
                 answ <- gconfirm(msg="X is in ascending order. Do you want to reverse X axis? ", title="CONFIRM REVERSE AXIS", AICON="WARNING")
                 if (answ == TRUE ){
                    DF[[Xidx]] <- rev(DF[[Xidx]]) #reverse X in descending order
                    DF[[Yidx]] <- rev(DF[[Yidx]]) #reverse X in descending order
                 } else {
                    svalue(reverseX) <- FALSE
                    svalue(NOreverseX) <- TRUE
                 }
              }
	             LL <- length(XPSSample)
 	            XPSSample[[LL]]@.Data[[4]] <<- DF[[Yidx]]
              plot(XPSSample, col="blue")
       }


#--- Variabili ---

       FNameIN <- NULL
       FName <- NULL
       XPSSample <- NULL
       activeFName <- NULL

#--- Widget definition ---
       ImportWin <- gwindow("Import Ascii Data", parent(100, 0), visible=FALSE)
       size(ImportWin) <- c(300, 400)
       MainGroup <- ggroup(horizontal=FALSE, spacing = 1, container = ImportWin)
       ImportGroup <- ggroup(horizontal=TRUE, container = MainGroup)
       OptGroup <- ggroup(horizontal=FALSE, spacing = 1, container = ImportGroup)
       WoutGroup <- ggroup(horizontal=FALSE, spacing = 1, container = ImportGroup)


#--- Import options ---
       read.opt <- gframe(text=" Import Options ", horizontal = FALSE, spacing=5,  container = OptGroup)
       LoadButt <- gbutton("Open Data File", handler=function(h,...){
                           FNameIN <<- gfile(text = "Select a file ...", type = "open",
				                              ,filter = list("Ascii files" = list(patterns = c(".asc",".txt", ".prn", ".dat"))),
					                      ,multi = FALSE)
                           activeFName <<- basename(FNameIN)
                           pathFile <- dirname(FNameIN)
	                   XPSSample <<- new("XPSSample",
			                     Project = " ",
					     Comments = " ",
					     User=Sys.getenv('USER'),
					     Filename=activeFName )
                           setwd(pathFile)
                           svalue(raw_input) <-  paste(readLines(FNameIN), collapse="\n")
                           enabled(import_btn) <- TRUE
                           enabled(addErr_btn) <- FALSE
                           enabled(save_btn) <- FALSE
                           enabled(AddToXPSSamp) <- FALSE
                           enabled(exit_btn) <- FALSE
                        }, container=read.opt)

#--- read options
       HeaderLayout <- gformlayout(container = read.opt)

       heading <- gradio(c("Yes", "No"), horizontal=TRUE, label="Header", selected = 2,  handler=function(h, ...){
                             answ <- svalue(heading)
                             if (answ == "Yes") {
                                svalue(NRowHeader) <- "1"
                               enabled(NRowHeader)  <-  TRUE
                              } else if (answ == "No") {
                                svalue(NRowHeader) <- "0"
                                enabled(NRowHeader) <- FALSE
                             }
                        }, container = HeaderLayout)
       NRowHeader <- gedit(text="0", label="Header Rows", container=HeaderLayout)
       enabled(NRowHeader) <- FALSE
       gseparator(container = read.opt) # separator

       OptLayout <- gformlayout(container = read.opt)

       seps <- c("Tab"="\t", "Whitespace"=" ", "Comma" = ",", Semicolon=";", "Unspecified"="")
       gcombobox(names(seps), selected=1, label="Separator", container = OptLayout)

       decs <- c("Period"=".", "Comma"=",")
       gcombobox(names(decs), label="Decimal", container = OptLayout)

       quotes 	<- c("No quote" = "", "Double quote (\")" = '"', "Single quote (')" = "'")
       gcombobox(names(quotes), label="Quote", container = OptLayout)

       gseparator(container = read.opt) # separator
       FmtGroup <- gpanedgroup(expand=TRUE, horizontal=TRUE,container = read.opt)
       scanFile <- gcheckbox(text="Scan DataFile", checked = FALSE, handler=function(h, ...){
                            scf <- svalue(scanFile)
                            if (scf==TRUE) {
                               enabled(dataCol) <- TRUE
                            } else {
                               enabled(dataCol) <- FALSE
                            }
                        }, container = FmtGroup)
       dataCol <- gedit(initial.msg="Data Ncol", container=FmtGroup)
       enabled(dataCol) <- FALSE

#--- Try to read data with selected options and unpdate OUTPUT data gtable to see if reading options are correct
       gbutton("Try to read data", spacing=5, handler=function(h, ...) update_output(), container = read.opt)

#---  Asci data Information to store Ascii_Data in a XPSSample DataFrame
       Elementframe <- gframe(text=" XPS Core Line ", spacing=5, horizontal = FALSE, container = OptGroup)
       CLparam	<- gformlayout(container = Elementframe)
       CLname <- gedit("", initial.msg=" Core Line Name", container = CLparam, label="Core Line Name")
       unitX <- gedit("Binding Energy [eV]", container = CLparam, label="X Units")
       unitY <- gedit("Intensity [counts]", container = CLparam, label="Y Units")
       CLyesno <- ggroup( horizontal=TRUE, container = Elementframe)
       glabel(text="Reverse X axis?", container = CLyesno)
       reverseX <- gcheckbox("Yes", selected=FALSE, handler=function(h, ...){
                    svalue(NOreverseX) <- FALSE
		                  enabled(import_btn) <- TRUE
                 }, container = CLyesno)
       NOreverseX <- gcheckbox("No", selected=FALSE, handler=function(h, ...){
                    svalue(reverseX) <- FALSE
		                  enabled(import_btn) <- TRUE
                 }, container = CLyesno)


#--- define INPUT Window
       InputGroup <- ggroup(horizontal=FALSE, container = WoutGroup)
       WinIn <- gvbox(expand=TRUE, fill=TRUE, spacing=5, container = InputGroup)
       labtxt	<- glabel(gettext("Input data:"), container = WinIn, anchor=c(-1,0))
       font(labtxt) <- list(weight="bold")
       raw_input <- gtext("", wrap=FALSE,  #
                   font.attr=list(family="monospace"),
                   container=WinIn, expand=TRUE, fill=TRUE)
       size(raw_input) <- c(200,150)

#--- define Load Window
       WinOUT <- ggroup(horizontal=FALSE, spacing=5, container =InputGroup )
       labtxt <- glabel("Loaded data:", container = WinOUT)
       font(labtxt) <- list(weight="bold")
       output <- gtable("", container = WinOUT)
       size(output) <- c(200,150)

#--- Which column to read?
       Colframe <- gformlayout(container = WoutGroup)
       ColX <- gedit(text="", container = Colframe, label="ColX to read (1 ?)")
       ColY <- gedit(text="", container = Colframe, label="ColY to read (2 or 3 or...)")


##--- BUTTONS
       ButtLT <- glayout(spacing = 5, container = MainGroup)
       import_btn <- ButtLT[1,1] <- gbutton(" IMPORT ", handler=function(h,...) {
                         if (! check_selection()){return()}  #controls all the needed information are given
                         addCoreLine()    #add a new coreline
                         LL <- length(XPSSample)
                         cat("\n ----- Data Info -----")
                         cat("\n ===> Data File: ", activeFName, ", Core Line: ", XPSSample[[LL]]@Symbol)
                         cat("\n ===> Xmin= ", min(XPSSample[[LL]]@.Data[[1]]), "Xmax: ", max(XPSSample[[1]]@.Data[[1]]))
                         cat("\n ===> Ymin= ", min(XPSSample[[LL]]@.Data[[2]]), "Ymax: ", max(XPSSample[[1]]@.Data[[2]]))
                         cat("\n")
                         enabled(addErr_btn) <- TRUE
                         enabled(save_btn) <- TRUE
                         enabled(AddToXPSSamp) <- TRUE
                         enabled(exit_btn) <- TRUE
       }, container = ButtLT)

       addErr_btn <- ButtLT[1,2] <- gbutton(" ADD ERRORS ", handler=function(h, ...){
                         if (! check_selection()){return()}  #controls all the needed information are given
                         addErrors()    #add a new coreline
                         LL <- length(XPSSample)
                         cat("\n ----- Data Info -----")
                         cat("\n ===> Data File: ", activeFName, ", Core Line: ", XPSSample[[LL]]@Symbol)
                         cat("\n ===> Standard Deviation Error Added to Last Saved Data", min(XPSSample[[LL]]@.Data[[1]]), "Xmax: ", max(XPSSample[[1]]@.Data[[1]]))
                         cat("\n")
                     }, container = ButtLT)    #save without exiting ImportAscii(): possibility to import other data from ascii file

       save_btn <- ButtLT[1,3] <- gbutton(" SAVE ", handler=function(h, ...){
                         LL <- length(XPSSample) #number of Corelines of the source XPSSample
                         assign(activeFName, XPSSample, envir=.GlobalEnv)  #save the XPSSample in the .GlobalEnv
                         assign("activeFName", activeFName, envir=.GlobalEnv)  #Set the activeSpectName to the last name of imported data
                         assign("activeSpectName", XPSSample[[LL]]@Symbol, envir=.GlobalEnv)
                         assign("activeSpectIndx", LL, envir=.GlobalEnv)   #set the activeSpectIndx to the last imported data
                         XPSSaveRetrieveBkp("save")   
                         enabled(addErr_btn) <- FALSE
                         enabled(save_btn) <- FALSE
                         enabled(AddToXPSSamp) <- FALSE
                     }, container = ButtLT)    #save without exiting ImportAscii(): possibility to import other data from ascii file

       AddToXPSSamp <- ButtLT[1,4] <- gbutton(text = " Save in an existing XPS Sample ", label = " ", checked = FALSE, handler = function(h, ...) {
                                   XPSSamplesList <- XPSFNameList()
		                   if (length(XPSSamplesList) > 0 ) {
			               SelectWin <- gwindow(" SELECT SAMPLE ", visible=FALSE)
			               size(SelectWin) <- c(150,250)
			               gwinsave	<- ggroup(container = SelectWin)
                                       samplesfr <- gframe(" XPS Samples ", spacing=5, horizontal=FALSE, container= gwinsave)
                                       SampleName <- gtable(items=XPSSamplesList, container=samplesfr)
                                       size(SampleName) <- c(120,200)
                                       gbutton("Select", handler = function(...) {
                                            if (length(svalue(SampleName)) > 0 ) {
		 	                        activeFName	<<- svalue(SampleName)
                                                dispose(SelectWin)
                                                FName <<- get(activeFName, envir=.GlobalEnv)
                                                LL <- length(FName)      #Number of corelines in the destination XPSSample
                                                CLnames <- names(FName)
                                                LLL <- length(XPSSample) #Number of CoreLines in the source XPSSample
                                                FName[[LL+1]] <<- XPSSample[[LLL]] #save last imported Corelines in the destinaton XPSSample
                                                FName@names <<- c(CLnames, svalue(CLname)) #add names of new CoreLines
                                                XPSSample<<-FName      #set the source XPSSample == destination file with all spectra
                                                assign(activeFName, FName, envir=.GlobalEnv)  #Save the destination XPSSample in GlobalEnv
                                                assign("activeFName", activeFName, envir=.GlobalEnv)
                                                assign("activeSpectName", FName[[LL+1]]@Symbol, envir=.GlobalEnv)
                                                assign("activeSpectIndx", 1, envir=.GlobalEnv)
                                                plot(FName)
                                                cat("\n Data saved in ", activeFName)
                                                XPSSaveRetrieveBkp("save")
                                                enabled(addErr_btn) <- FALSE
                                                enabled(save_btn) <- FALSE
                                                enabled(AddToXPSSamp) <- FALSE
                                            } 
                                       }, container = samplesfr)
                                       visible(SelectWin) <- TRUE
                                   }   
                     }, container = ButtLT)

       exit_btn	<- ButtLT[1,5] <- gbutton(" SAVE and EXIT ", handler=function(h, ...){
                         dispose(ImportWin)
                         LL <- length(XPSSample)
                         assign(activeFName, XPSSample, envir=.GlobalEnv)  #save XPSSample in GlobalEnv
                         assign("activeFName", activeFName, envir=.GlobalEnv)  #Set the activeSpectName to the last name of imported data
                         assign("activeSpectName", XPSSample[[LL]]@Symbol, envir=.GlobalEnv)
                         assign("activeSpectIndx", LL, envir=.GlobalEnv)   #set the activeSpectIndx to the last imported data
                         XPSSaveRetrieveBkp("save")
                         return(XPSSample)
                     }, container = ButtLT)




       enabled(import_btn) <- FALSE
       enabled(addErr_btn) <- FALSE
       enabled(save_btn) <- FALSE
       enabled(AddToXPSSamp) <- FALSE
       enabled(exit_btn) <- FALSE

       visible(ImportWin) <- TRUE
       ImportWin$set_modal(TRUE)

       return(XPSSample)
}     




