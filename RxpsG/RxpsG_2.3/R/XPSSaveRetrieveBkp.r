XPSSaveRetrieveBkp<-function(opt){

#--- get System info and apply correspondent XPS Settings ---
   OS<-Sys.info() #get system information
#Reading XPS settings which can be customized by users
   Bkp.pthName <- "/RxpsG/data/BkpData.Rdata"
   LibPth <- .libPaths()  #path of all the R libraries where RxpsG could be located
   LL <- length(LibPth)
   Bkp.pthName <- paste(LibPth, Bkp.pthName, sep="") #paste acts on all the components of the LibPth vector
   if (OS["sysname"] != "Linux") {  #Windows and Mac OS systems
      Bkp.pthName <- gsub("/","\\", Bkp.pthName, fixed=TRUE)   #path/filename for linux,  path\\filename for windows
   }   #gsub acts on a vector of path names
   for(ii in 1:LL){
      fe <- file.exists(Bkp.pthName[ii])  #search for the BKpFile in all the R-Library locations
      if (fe == TRUE) {     # in Windows two library location could be present: (1) under c:/program Files, (2) under c:/user/documents/
         Bkp.pthName <- Bkp.pthName[ii] #File pointer to the BKP-file location
         break
      }
   }
   if (fe == FALSE) {
      gmessage("File BkpData.Rdata not found in the folder /Data of the RxpsG library!",title="BACKUP FILE MISSING", icon="error")
      return()
   }

   switch(opt,
       "save"={
          FNameList<-XPSFNameList() #read the list of XPSSample loaded in the .GlobalEnv
          save(list=FNameList, file=Bkp.pthName, compress=TRUE)   #save all the XPSSamples

#          txt<-paste("\n Bkp-Data saved in: ", SaveFile, sep="")
#          cat("\n", txt)
       },
       "retrieve"={
          FNameList<-load(Bkp.pthName,envir=.GlobalEnv)   #load the data in the .GlobalEnv not in the local memory
                                                          #assign not necessary
       } )

}
