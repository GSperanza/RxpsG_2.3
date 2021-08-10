 q<- function(save="no", ...){    #to exit without requiring to save workspace
     quit(save=save,...)
 }


options(stringsAsFactors=FALSE)  # the default for DataSets in RxpsG


setHook(packageEvent("RxpsG", "attach"), function(...) {   # setHook forces the execution of xps() only after loading the package RxpsG
                       cat("\014")  #clears R consolle
                       xps()            
                       cat("\n Welcome to Rxps!")
                       setwd("~") #sets WD in ..users/NomeUtente/Documents (windows)     .../home/  (linux)
                      })

 DefPackages<-getOption("defaultPackages")   #loads default packages
 opt<-options(defaultPackages = c(DefPackages, "RxpsG")) #attach RxpsG


.Last <- function()  cat("\n   Goodbye!\n\n")
