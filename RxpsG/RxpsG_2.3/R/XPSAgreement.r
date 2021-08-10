#'
#'Agreement to use the RxpsG package
#'
#'@examples
#'
#'\dontrun{
#'XPSAgreement()
#'}
#'
#'@export
#'


XPSAgreement<-function(){
       MainWin<-gwindow("SOFTWARE AGREEMENT", visible=FALSE)
#       size(MainWin) <- c(20, 50)
       MainGroup <- ggroup(horizontal=FALSE, container = MainWin)

       msg<-"           Welcome to RXPSG!

       This program is free software. You can use it under the
       terms of the  GNU Affero General Public License.
       http://www.gnu.org/licenses/
       and licenses linked to the use of the R, Rstudio platforms.

       - Authors decline any responsibility deriving from
         the use of the software or from software bugs.

       - Users are kindly requested to cite the software
         authors in their publications:
       
         Giorgio Speranza, Roberto Canteri
         Fondazione Bruno Kessler, Sommarive str. 18
         38123 Trento Italy.  "



       txt <- gtext(text=msg, wrap=FALSE, container=MainGroup)
       font(txt)<-list(weight="normal",style="normal", family="sans", size="small")
       size(txt) <- c(400, 300)
       ButtGroup <- ggroup(horizontal=TRUE, container = MainGroup)
       Butt1<-gbutton("Agree", handler=function(h, ...){
                          xps()
                          dispose(MainWin)
                      }, container=ButtGroup)
       Butt2<-gbutton("Decline", handler=function(h, ...){
                          dispose(MainWin)
                      }, container=ButtGroup)
       visible(MainWin)<-TRUE

#--- FONT ATTRIBUTES
#weight
#in c("light", "normal", "bold", "heavy")
#style
#in c("normal", "oblique", "italic")
#family
#in c("sans", "helvetica", "times", "monospace")
#size
#in c("xx-small", "x-small", "small", "medium", "large", "x-large", "xx-large")

}
