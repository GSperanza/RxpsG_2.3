#Tougaard bacgrkund test macro
#non optimized version
#See XPSCustomBaseline.r

TougaardBkg<-function(object, limits){

   LL <- length(object@RegionToFit[[1]])
   tougaard<-vector(mode = "numeric", length = LL)   #to store the backgound
   lK<-vector(mode = "numeric", length = LL)         #vector conttaining the Universal Cross Section


   C<-1643
   cat("\n New C: ")
   C<-scan(n=1, quiet=TRUE)

   tmp<-0
   dE <-(object@RegionToFit[[1]][2] - object@RegionToFit[[1]][1]) # Estep
   
    
   for (ii in 1:LL){
       deltaE1<-ii*dE^2
       deltaE2<-(C+deltaE1^2)^2
       lK[ii]<-deltaE1/deltaE2
       tmp<-tmp+deltaE1*object@RegionToFit[[2]][ii]/deltaE2   #j(E')(E'-E)/{C+[E'-E]^2}^2
   }

   if (object@Flags[1] == TRUE){  #BE scale
      avg1<-sum(object@RegionToFit[[2]][1:5])/5
      avg2<-sum(object@RegionToFit[[2]][(LL-4):LL])/5
      level<-limits$y[1]
   } else {
      avg1<-sum(object@RegionToFit[[2]][(LL-4):LL])/5
      avg2<-sum(object@RegionToFit[[2]][1:5])/5
      level<-limits$y[2]
   }
   object@RegionToFit[[2]]<-object@RegionToFit[[2]]-avg2    #set the less intense part of the spectrum to 0
   B<-level/tmp
   cat("\n >>> Parameter B= ", B)

#--- Computation Tougaard background: OK for both KE and BE scale

   tougaard<-vector(mode = "numeric", length = LL)
cat("\n Background Data points:", LL)
#scan(n=1, quiet=TRUE)

#   k<-1
#   for (ii in 1:LL){
#       tmp<-0
#       for (jj in ii:LL){
#           deltaE<-object@.Data[[1]][jj]-object@.Data[[1]][ii]  #attenzione: Ekin(j)-Ekin(i)=BE(i)-BE(j) la formula di tougaard funziona sulle Ekin!!!!!!
#    	     denom<-(C+deltaE^2)^2                                #{C+[E'-E]^2}^2
#    	     tmp<-tmp+deltaE*object@.Data[[2]][jj]*dE/denom     # B*j(E')(E'-E)/{C+[E'-E]^2}^2
#       }
#
#
#       tougaard[k]<-B*tmp
#cat("\n Processing data ", k, tougaard[k])
#       k<-k+1
#   }
   

   for (ii in 1:LL){
       tougaard[ii]<-B*sum(lK[1:(LL-ii)]*object@RegionToFit[[2]][ii:LL])
   }


   tougaard<-tougaard+avg2
   return(tougaard)
}
