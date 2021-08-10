# General purpose functions
#
# Note on FLAGS in XPSSamples:
# Flag[1] == Binding Energy --> If TRUE then x scale is Binding Energy
# Flag[2] == cps 	--> If TRUE then y scale is cps. NB!! With vamas and scienta files is always TRUE!!
# Flag[3] == Scienta --> If TRUE then sample is Scienta sample
# Flag[4] == Vamas correction for transmission factor
##==============================================================================

##==============================================================================
# findXIndex: finds the index corresponding to a given value in an array of abscisse
##==============================================================================
#'findXIndex function utility
#'
#'finds the index of X vector corresponding to a given abscissa value
#'
#'@param X numeric vector
#'@param Value numeric value
#'@return Returns X index
#'
#'@examples
#'
#'\dontrun{
#'	findXIndex(X, Value)
#'}
#' 
#'@export
#'

findXIndex <- function(X, Value) {

   if (is.na(Value)){
       gmessage("WARNING: value to search for index == NA", title="BAD SEARCH VALUE", icon="error")
       return()
   }
   RngX <- range(X)
   LL <- length(X)
   DE <- diff(X[1:2])
   if (abs(Value-X[1]) < abs(DE/10^3)) return(1)      #to avoid probems at X edges due to decimal errors in double precision numbers
   if (abs(Value-X[LL]) < abs(DE/10^3)) return(LL)
   if (Value > max(RngX) || Value < min(RngX)){
      gmessage(msg="X-Value out of Core Line range!", title = "WARNING: WRONG LIMITS",  icon = "warning")
      return()
   }
   if ( DE > 0 ) {  ## if data are increasing sorted
      return(max(which(X<=Value)))
   } else {
      return(max(which(X>=Value)))
   }
}

##==============================================================================
# findYIndex: finds the index corresponding to a given value in an array of ordinate
##==============================================================================
#'findYIndex function utility
#'
#'finds the index of Y vector corresponding to a given Value
#'
#'@param Y numeric vector
#'@param Value numeric value
#'@param tolerance numeric value, precision required to find out Y
#'@return Returns Y index
#'
#'@examples
#'
#'\dontrun{
#'	findYIndex(Y, Value, tolerance=0.01)
#'}
#' 
#'@export
#'

findYIndex <- function(Y, Value, tolerance=0.01) {

   if (is.na(Value)){
       gmessage("WARNING: value to search for index == NA", title="BAD SEARCH VALUE", icon="error")
       return()
   }
   Ymax <- max(Y)
   Ymin <- min(Y)
   LL<-length(Y)
   DE<- (Ymax-Ymin)*tolerance #tolerance defines the range around Value where to search for Y index
   Yindx <- NULL
   if (Value > Ymax || Value < Ymin){
      gmessage(msg="Y-Value out of Core Line range!", title = "WARNING: WRONG LIMITS",  icon = "warning")
      return()
   }
   jj <- 1
   for (ii in 1:LL) {
       if (Y[ii] < Value+DE && Y[ii] > Value-DE){
          Yindx[jj] <- ii
          jj <- jj+1
       }
   }
   return(Yindx)
}




##==============================================================================
# findY: finds the Y value corresponding to a given X value
##==============================================================================
#'findY function utility
#'
#'finds the Y value corresponding to a given X value
#'
#'@param XY numeric XY matrix
#'@param Xvalue numeric value
#'@return Returns Y value
#'
#'@examples
#'
#'\dontrun{
#'	findYIndex(XY, Xvalue)
#'}
#' 
#'@export
#'
findY <- function(XY, Xvalue) {    #XY is a list of X and Y numbers for example FName@RegionToFit
   idx<-findXIndex(XY$x, Xvalue)   #XY$x is the array ob abscissas, Value the Xvalue to look for the correspondent Yvalue
   Yvalue<-XY$y[idx]
   return(Yvalue)
}


##==============================================================================
# findMaxPos: finds the abscissa corresponding to the Y max
##==============================================================================
#'findMaxPos function utility
#'
#'finds the abscissa corresponding to the y-vector maximum
#'
#'@param XY numeric XY matrix
#'@return Returns c(x, max) values
#'
#'@examples
#'
#'\dontrun{
#'	findYIndex(XY)
#'}
#' 
#'@export
#'

findMaxPos <- function(XY) {    #XY is a list of X and Y numbers for example FName@RegionToFit
   MM<-max(XY$y)
   idx<-which(XY$y==MM)
   BE<-XY$x[idx]                #BE corresponding to the spectrum max
   return(c(BE,MM))
}


##==============================================================================
# labformula used in plot but useful also for label any type of plot
##==============================================================================
labformula <- function(position.list, label, color="blue", ...) {
  ## position.list = list(x,y)
  X0 <- position.list$x
  RngX<-par("usr")[1:2]
  X1 <- X0 + diff(par("usr")[1:2])*0.020 ## orizzontal dim 2.0% of xrange
  yspan <- diff(par("usr")[3:4])*0.030   ## vertical dim 3% of yrange
  Yspan <- rep(yspan, length(X0))
  # reverse yspan if y+2*yspan > ymax 
  if (any((position.list$y + 2*yspan) > par("usr")[4])) { 
    Yspan[which((position.list$y + 2*yspan) > par("usr")[4])] <- -yspan
  }
  Y0 <- position.list$y + 2*Yspan
  Y1 <- position.list$y + Yspan
  arrows(X0, Y0, X0, Y1, length=0.08, col=color) # arraw from(x0,y0) to (x1,y1)
  segments(X0,Y0, X1, Y0, col=color) ##  horizontal segment
  text(X1, Y0, labels = label, col=color, ...) ## label
}


##==============================================================================
# MakeTable to write tables using ASCII text parameters: Font=GCourierNew
##==============================================================================
#'Font=GCourierNew is the CourierNew modified with symbols for drawing tab borders.
#'Provides commands to print borders and cell text alignment (left, center, right)
#'Borders are made using ascii characters
#'
#'@param Type string = BorderUp, BorderIn, BorderBottom, TableRow: prints border or table rows
#'@param txt array containing the txt of the table row
#'@param CellLength array containing the number of character of each table column
#'@param align string = left, center, right: text alignment of the table cells
#'@return Returns the processed \code{object}.
#'
#'@examples
#'
#'\dontrun{
#'   txt<-c("Components", "Area(cps)", "FWHM", "RSF", "BE(eV)", "TOT.(%)")    #Voci Tabella
#'   CellLgth<-c(12, 15, 8, 7, 8, 9)
#'   cell<-MakeTable("MultiCell", txt,CellLength, "center")
#'}
#'
#'@export
#'
    MakeTable <- function(Type,txt,CellLength, align){
         cell<-""
         if(Type=="BordTopM"){      # Bordo Top MultiCell tipo:   |```|`````|```````|```|
            HorBar<-rawToChar(as.raw(151))  #HorizontalBar
            Tdown<-rawToChar(as.raw(236))   # Carattere a forma di T
            cell<-rawToChar(as.raw(204))    #Curve left down
            CellLength[1]<-CellLength[1]-1  #ho aggiunto il carattere iniziale del bordo
            LL<-length(CellLength)          #txt puo' essere  una singola stringa o un vettore
            if (LL-1>0){
               for(ii in 1:(LL-1)){
                  Nbar<-CellLength[ii]-1    # -1 perche' aggiunge il Cross
                  Border<-paste(rep(HorBar,Nbar), collapse="")
                  cell<-sprintf("%s%s%s", cell, Border, Tdown)
               }
            }

            Nbar<-CellLength[LL]-1          # -1 perche' aggiunge il termine bordo
            Border<-paste(rep(HorBar,Nbar), collapse="")
            cell<-sprintf("%s%s%s", cell, Border, rawToChar(as.raw(208)))  #termine bordo curva right down
         }
         if(Type=="BordTopS"){      # Bordo Top SingleCell tipo:  |``````````````````````|
            HorBar<-rawToChar(as.raw(151))  #HorizontalBar
            Border<-paste(rep(HorBar,CellLength-2), collapse="")
            cell<-sprintf("%s%s%s", rawToChar(as.raw(204)), Border, rawToChar(as.raw(208)))  #termine bordo curva right down
         }
         if(Type=="BordInM"){        #Bordo interno MultiCell tipo |---+-----+-------+---|
            HorBar<-rawToChar(as.raw(151))  #HorizontalBar
            Cross<-rawToChar(as.raw(252))   #carattere a froma di croce
            cell<-rawToChar(as.raw(220))    #T orizzontale
            CellLength[1]<-CellLength[1]-1  #ho aggiunto il carattere iniziale del bordo
            LL<-length(CellLength)          #txt puo' essere  una singola stringa o un vettore
            if (LL-1>0){
               for(ii in 1:(LL-1)){
                   Nbar<-CellLength[ii]-1   # -1 perche' aggiunge il Cross
                   Border<-paste(rep(HorBar,Nbar), collapse="")
                   cell<-sprintf("%s%s%s", cell, Border, Cross)
               }
            }
            Nbar<-CellLength[LL]-1   # -1 perche' aggiunge il termine bordo
            Border<-paste(rep(HorBar,Nbar), collapse="")
            cell<-sprintf("%s%s%s", cell, Border, rawToChar(as.raw(228))) #termine bordo T orizzontale
         }
         if(Type=="BordInS"){        #Bordo interno SingleCell tipo |----------------------|
            HorBar<-rawToChar(as.raw(151)) #HorizontalBar
            Cross<-rawToChar(as.raw(252))  #carattere a froma di croce
            Border<-paste(rep(HorBar,CellLength-2), collapse="")
            cell<-sprintf("%s%s%s", rawToChar(as.raw(220)), Border, rawToChar(as.raw(228))) #termine bordo T orizzontale SX, DX
         }
         if(Type=="BordInMU"){        #Bordo interno MultiCellUP tipo |---'-----'-------'---| separatore solo sopra
            HorBar<-rawToChar(as.raw(151)) #HorizontalBar
            Tup<-rawToChar(as.raw(244))    #carattere a forma di T rovesciata in alto
            cell<-rawToChar(as.raw(220))   #T orizzontale
            CellLength[1]<-CellLength[1]-1 #ho aggiunto il carattere iniziale del bordo
            LL<-length(CellLength)         #txt puo' essere  una singola stringa o un vettore
            for(ii in 1:(LL-1)){
                Nbar<-CellLength[ii]-1     # -1 perche' aggiunge il Cross
                Border<-paste(rep(HorBar,Nbar), collapse="")
                cell<-sprintf("%s%s%s", cell, Border, Tup)
            }
            Nbar<-CellLength[LL]-1         # -1 perche' aggiunge il termine bordo
            Border<-paste(rep(HorBar,Nbar), collapse="")
            cell<-sprintf("%s%s%s", cell, Border, rawToChar(as.raw(228))) #termine bordo T orizzontale
         }
          if(Type=="BordInMD"){        #Bordo interno MultiCellDown separatore: |---,-----,-------,---|  separatore solo sotto
            HorBar<-rawToChar(as.raw(151)) #HorizontalBar
            Tdown<-rawToChar(as.raw(236))  # Carattere a forma di T verso basso
            cell<-rawToChar(as.raw(220))   #T orizzontale  SX
            CellLength[1]<-CellLength[1]-1 #ho aggiunto il carattere iniziale del bordo
            LL<-length(CellLength)         #txt puo' essere  una singola stringa o un vettore
            for(ii in 1:(LL-1)){
                Nbar<-CellLength[ii]-1     # -1 perche' aggiunge il Cross
                Border<-paste(rep(HorBar,Nbar), collapse="")
                cell<-sprintf("%s%s%s", cell, Border, Tdown)
            }
            Nbar<-CellLength[LL]-1         # -1 perche' aggiunge il termine bordo
            Border<-paste(rep(HorBar,Nbar), collapse="")
            cell<-sprintf("%s%s%s", cell, Border, rawToChar(as.raw(228))) #termine bordo T orizzontale DX
         }
        if(Type=="BordBotM"){          #Bordo Bottom MultiCellDown separatore: |___,_____,_______,___|
            HorBar<-rawToChar(as.raw(151)) #HorizontalBar
            Tup<-rawToChar(as.raw(244))    #carattere a forma di T rovesciata in alto
            cell<-rawToChar(as.raw(212))   #Curve left up
            CellLength[1]<-CellLength[1]-1 #ho aggiunto il carattere iniziale del bordo
            LL<-length(CellLength)         #txt puo' essere  una singola stringa o un vettore
            for(ii in 1:(LL-1)){
               Nbar<-CellLength[ii]-1      # -1 perche' aggiunge il Cross
               Border<-paste(rep(HorBar,Nbar), collapse="")
               cell<-sprintf("%s%s%s", cell, Border, Tup)
            }
            Nbar<-CellLength[LL]-1         # -1 perche' aggiunge il termine bordo
            Border<-paste(rep(HorBar,Nbar), collapse="")
            cell<-sprintf("%s%s%s", cell, Border, rawToChar(as.raw(216))) #termine bordo curva right up
         }
        if(Type=="BordBotS"){          #Bordo Bottom MultiCellDown separatore: |_____________________|
            HorBar<-rawToChar(as.raw(151)) #HorizontalBar
            Border<-paste(rep(HorBar,CellLength-2), collapse="")
            cell<-sprintf("%s%s%s", rawToChar(as.raw(212)), Border, rawToChar(as.raw(216)))  #termine bordo curva right down
         }
         if(Type=="SingleCell"){           #cella singola con testo :              | titolo:             |
            VertBar<-rawToChar(as.raw(194))
            Nspaces<-CellLength-nchar(txt)-2
            spaces<-paste(rep(" ",Nspaces), collapse="")
            cell<-sprintf("%s%s%s%s",VertBar,txt,spaces,VertBar)
         }
         if(Type=="MultiCell"){            #celle multiple separate + testo        |xxx| yy  | zzzzz | t |
           LL=length(txt)                  #txt e' un vettore
           VertBar<-rawToChar(as.raw(194))
           cell<-VertBar
           CellLength[1]<-CellLength[1]-1  #ho aggiunto il carattere iniziale del bordo
           for (ii in 1:LL){
               if (align=="right"){
                  Nspaces<-CellLength[ii]-nchar(txt[ii])-1   # -1 perche' aggiunge il VertBar
                  spaces<-paste(rep(" ",Nspaces), collapse="")
                  cell<-sprintf("%s%s%s%s", cell, spaces, txt[ii], VertBar)
               }
               if (align=="left"){
                  Nspaces<-CellLength[ii]-nchar(txt[ii])-1
                  spaces<-paste(rep(" ",Nspaces), collapse="")
                  cell<-sprintf("%s%s%s%s", cell, txt[ii], spaces, VertBar)
               }
               if (align=="center"){
                  Nspaces<-CellLength[ii]-nchar(txt[ii])-1
                  Nspaces1<-as.integer(Nspaces/2)
                  Nspaces2<-Nspaces-Nspaces1 #recupero spazio divisione non intera
                  spaces1<-paste(rep(" ",Nspaces1), collapse="")
                  spaces2<-paste(rep(" ",Nspaces2), collapse="")
                  cell<-sprintf("%s%s%s%s%s", cell, spaces1, txt[ii], spaces2, VertBar)
               }
            }
         }
         return(cell)
    }

##==============================================================================
# printCell to print data in table fashion using text just text chacters
##==============================================================================
#'printCell prints data in a table constructed by MakeTable()
#'Borders are made using ascii characters  - and |
#'
#'@param Type label for labels such as the title, 
#'         separator to print a ---------- like separator,
#'         tableRow to print a sequence of data in a table row
#'@param txt the text to print, single string in label mode, an array of strungs in tableRow
#'@param CellB cell border can be "|", " " or ""
#'@param CellLength = in tableRow mode is an array containing the width (nchar) of each table column,
#'                    in label mode it is the witdth of the tableRow (sum of column widths)
#'@param align text alignment used in tableRow mode: left, center, right
#'
#'@examples
#'
#'\dontrun{
#'   CellLgth<-c(12, 15, 8, 7, 8, 9)
#'   cell<-printCell(Type="separator", txt="-", cellB="", CellLength=sum(CellLgth), align="")
#'   txt<-"XPS Au24.vms"
#'   cell<-printCell("label", txt, "|", sum(CellLgth), "")
#'   txt<-c("Components", "Area(cps)", "FWHM", "RSF", "BE(eV)", "TOT.(%)")    #Voci Tabella
#'   cell<-printCell("tableRow", txt, "|", CellLgth, "center")
#'   cell<-printCell("separator", "-", "|", sum(CellLgth), "")
#'}
#'
#'@export
#'

    printCell <- function(Type,txt,CellB, CellLength, align){
         cell<-""
         LL=length(txt)   #txt may be a single string or a vector

         if(Type=="label"){                       #if txt is a vector select the max of nchar
            Nspaces<-max(CellLength-nchar(txt)+6) #6 are the character used for the table vertical separators
            spaces<-paste(rep(" ",Nspaces), collapse="")
            cell<-paste(txt, spaces, sep="")
         }
         if(Type=="separator"){
            Nspaces<-0
            txt<-paste(rep(txt,CellLength+6), collapse="")  #6 sono i caratteri aggiuntivi per fare i separatori verticali nella tabella
            cell<-paste(txt, sep="")
         }
         if(Type=="tableRow"){
            for (ii in 1:LL){
               if (align=="right"){
                  Nspaces<-CellLength[ii]-nchar(txt[ii])
                spaces<-paste(rep(" ",Nspaces), collapse="")
                  cell<-paste(cell,CellB, spaces, txt[ii], sep="")
               }
               if (align=="left"){
                  Nspaces<-CellLength[ii]-nchar(txt[ii])
                  spaces<-paste(rep(" ",Nspaces), collapse="")
                  cell<-paste(cell,CellB, txt[ii], spaces, sep="")
               }
               if (align=="center"){
                  Nspaces<-CellLength[ii]-nchar(txt[ii])
                  Nspaces1<-as.integer(Nspaces/2)
                  Nspaces2<-Nspaces-Nspaces1 #recupero spazio divisione non intera
                  spaces1<-paste(rep(" ",Nspaces1), collapse="")
                  spaces2<-paste(rep(" ",Nspaces2), collapse="")
                  cell<-paste(cell,CellB,spaces1, txt[ii], spaces2, sep="")
               }
            }
            cell<-paste(cell, CellB, sep="")
         }
         return(cell)
    }


##==============================================================================
# MinDist Min distance between point P0 and points P1, P2
##==============================================================================
#MinDist function utility
#'
#'finds the minimum distance between point P0 and points P1, P2
#'
#'@param posP0 numeric position of probe point
#'@param posP1 numeric position of reference point 1
#'@param posP2 numeric position of reference point 2
#'@return Returns index 1 or 2 indicating min distance from P1 or P2
#'
#'@examples
#'
#'\dontrun{
#'	MinDist(posP0, posP1, posP2)
#'}
#' 
#'@export
#'
MinDist<-function(posP0, posP1, posP2){   #minima distanza tra un punto e altri due punti nel piano XY
         D1<-((posP0$x-posP1$x)^2 + (posP0$y-posP1$y)^2)^0.5  #dist P0 P1
         D2<-((posP0$x-posP2$x)^2 + (posP0$y-posP2$y)^2)^0.5  #dist P0 P2
         if (D1 < D2) return(1)
         if (D2 < D1) return(2)
}

##==============================================================================
# FitLin: linear fitting functin
##==============================================================================
#'FitLin function utility
#'
#'linear fit function
#'
#'@param X abscissa
#'@param Y ordinate
#'@return Returns c(m, c) slope m and intercept c
#'
#'@examples
#'
#'\dontrun{
#'	MinDist(posP0, posP1, posP2)
#'}
#' 
#'@export
#'
FitLin<-function(X, Y){   #make the linear fit on a set of data

     SumXi<-0         # y=mx + c
     SumXiXi<-0       # m = [ N*Sum(XiYi) - Sum(Xi)Sum(Yi)]/[N*Sum(XiXi)-(Sum(Xi))^2]
     SumXiYi<-0       # c = [Sum(Yi)*Sum(XiXi)-Sum(XiYi)*Sum(Xi)]/[N*Sum(XiXi)-(Sum(Xi))^2]
     SumYi<-0

     N<-length(X)
     for(ii in 1:N){
         SumXi<-SumXi+X[ii]
         SumXiXi<-SumXiXi+X[ii]*X[ii]
         SumXiYi<-SumXiYi+X[ii]*Y[ii]
         SumYi<-SumYi+Y[ii]
     }
     m <- (N*SumXiYi - SumXi*SumYi)/(N*SumXiXi - SumXi^2)
     c <- (SumYi*SumXiXi - SumXiYi*SumXi)/(N*SumXiXi - SumXi^2)
     return(c(m, c))
}














