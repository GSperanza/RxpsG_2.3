#Depth profile sinthesizer

XPSDPsynth <- function(){

    Find.Elmt <- function(jj){
        ElemTab$Element <- paste(ElemTab$Element, ElemTab$Orbital, sep="")
        idx <- which(ElemTab$Element == ElmtName[jj], arr.ind = TRUE)
        idx <- idx[1]
        RSF[jj] <<- ElemTab$RSF_K[idx]
        KE[jj] <<- ElemTab$KE[idx]
        AttLth[jj] <<- (1000/Density)*( 49/KE[jj]^2 + 0.11*KE[jj]^0.5) #attenuation length at energy correspondent to jj coreline position.
    }

    Calc.TF <- function(){
    #Function calculating the Transfer function TF
    #TF is a function of the element Att.Length, of the depth (i.e. which layer) at which the element is and of the tilt angle
    #for this reason values of TF are contained in a 3 dimensional matrix.
    #The function computes also the sensitivity factor Rjt[j,tilt] for eac element j and tilt angle
    #Finally the function computes the normalization factor to compute the concentrations i.e. SUM_j( Intensity[j,tt]/Rjt[j,tt] )

        LBeer <- function(lambda, angle){  # Lambert Beer decay function for the generic layer ii: at depth ii*thickness the decay is TF^(ii-1) see eqs 3,4,5 Livesey
            return(exp(-LThickness/(lambda*sin(angle))))   #Attention here the tilt is the angle between the sample surface and the vertical analyzer axis => sin(Tilt)
        }                                                  #in Livesey tilt is the angle of photoelectron emission measured respect to the analyzer axis => cos(Tilt)

        TF <<- array(dim=c(NElmnts, NLayers, NAngles))     #defines a 3D matrix with the indicated dimensions containing the value of the transfer function
                                                           #for all the tilt angle all the samples and all the depth=ii*LThickness of the Nlayers
        Rjt <<- matrix(nrow=NElmnts, ncol=NAngles)         #initialize the matrix containing the transfer function associated to the element jj through the layers for given tilt
        for(tt in 1:NAngles){ #this loop is on the tilt since we want to model the intensity of element jj as a function of the tilt tt
            for(jj in 1:NElmnts){     #this loop is on the different elements of the layers
                for(ii in 1:NLayers){
                    TF[jj,ii,tt] <<- LBeer(AttLth[jj], Tilt[tt])^(ii-1)  #contribution to the intensity of element ii in layer jj at tilt tt
                }
                Rjt[jj,tt] <<- sum(TF[jj, ,tt]) #1/kk * sensitivity factor for element jj (see eqq.6 Livesey) for given tilt
#cat("\n\n xxx", jj, tt, Rjt[jj,tt])
            }
            cat("\n *** Tilt",Tilt[tt]*180/3.1415, "Transf.Funct. ",Rjt[ ,tt])
        }
        cat("\n *** Transfer Function,and Normalyzing Factor R ready")
#scan(n=1, quiet=TRUE)
    }

    SynthDP <- function(){
        dpth <- NULL

        cat("\n NElmnts", NElmnts)
        cat("\n Elmnt.Name", ElmtName)
        cat("\n Att.Length", AttLth)
        cat("\n NAngles", NAngles)
        cat("\n Dteta", Dteta)
        cat("\n NLayers", NLayers)
        cat("\n Lthickness", LThickness)
        cat("\n From", Depth[ ,1])
        cat("\n To", Depth[ ,2])
        cat("\n Concentration", ElmtConc)

#Set Depth, Tilt, Nji, TF, Rjt arrays
        for(jj in 1:NElmnts){
                Depth[jj, ] <<- sort(Depth[jj, ], decreasing=FALSE) #Depth goes from surface to NLayers
        }
        for(tt in 1:NAngles){
            Tilt[tt] <<- (90-(tt-1)*Dteta)*pi/180 #series of tilt angles in radiants
        }
        for(ii in 1:NLayers){ #Depth increases from surface ii=0 to MaxDepth ii=NLayers
            dpth[ii] <- (ii-1)*LThickness
            Ntot <- 0
            for(jj in 1:NElmnts){
                if(Depth[jj,1] <= dpth[ii] && Depth[jj,2] >= dpth[ii]){  #element jj is present at depth rangung from Depth[jj,1] to Depth[jj,2]
                   Nji[jj,ii] <<- ElmtConc[jj]
                }
                Ntot <- Ntot+Nji[jj,ii]   #Total element concentration of layer ii
            }
            Nji[ ,ii] <<- Nji[ ,ii]/Ntot    #Concentration normalization condition for layer ii
        }
        Calc.TF()
#Estimate element concentrations and intensities
        for(tt in 1:NAngles){
            Itot <- 0
            for(jj in 1:NElmnts){
                Ijt[jj,tt] <<- sum(Nji[jj, ]*TF[jj, ,tt])/Rjt[jj,tt]
                Itot <- Itot + Ijt[jj,tt]
            }
            Xjt[ ,tt] <<- Ijt[ ,tt]/Itot #simulated element concentration for tilt angle tt
            for(jj in 1:NElmnts){
                Variance[jj,tt] <<- 10-Xjt[jj,tt]*10  #simulated variance affecting the element concentration Xjt
            }
        }
#Plot Data
        xx <- yy <- NULL
        for(jj in 1:NElmnts){
            xx <- cbind(xx, dpth, deparse.level=0)
            yy <- cbind(yy, Nji[jj, ])
        }
        txt <- paste(ElmtName[jj], " conc. (%)", sep="")
        clr <- seq(1:NElmnts)
        matplot(x=xx, y=yy, ylim=c(0,1.2), main="Element conc. profile",
                type="b", lty=c(1, 1), pch=sym, cex=0.9, col=clr, xlab="Depth nm.", ylab=txt)
        legend("top", legend=ElmtName, pch=sym, lty=rep(1,NElmnts), ncol=3, col=clr)
        scan(n=1, quiet=TRUE)
#----------
        xx <- yy <- NULL

        for(jj in 1:NElmnts){
            xx <- cbind(xx, Tilt*180/3.1415, deparse.level=0)
            yy <- cbind(yy, Ijt[jj, ])
        }
        txt <- paste(ElmtName[jj], "Intensity [A.U.]", sep="")
        matplot(x=xx, y=yy, ylim=c(0, 1.2), main="Spectral Intensities v.s Tilt",
                type="b", lty=c(1, 1), pch=sym, cex=0.9, col=clr, xlab="Tilt deg.", ylab=txt)
        legend("top", legend=ElmtName, pch=sym, lty=rep(1,NElmnts), ncol=3, col=clr)

    }



    NElmnts <- 1     # N. elements used to build the final spectral intensity
    ElmtName <- NULL #element names
    ElName <- list() #list of pointers to widgets for Elmt Name input
    SensFact <- NULL #Eelement sensitivity factors
    RSF <- list()    #list of pointers to widget for RSF input
    KE <- NULL       #element BE
    AttLth <- NULL   #element attenuation length
    From <- list()   #list of pointers to widget for element depth input
    To <- list()     #list of pointers to widget for element depth input
    Depth <- NULL    #element depth distribution
    Conc <- list()   #list of pointers to widget for element real concentration
    ElmtConc <- NULL #element concentrations
    NAngles <- NULL  #N. tilt
    Dteta <- NULL    #tilt increment
    NLayers <- NULL  #number of layers to discretize the profile
    LThickness <- NULL #Layer thickness
    Density <- NULL  #Material density
    Ijt <- NULL      #spectral intensity of element j
    Xjt <- NULL      #simulated apparent concentration as measured on the surface by XPS
    Nji <- NULL      #depth distribution of element concentration
    Variance <- NULL #Variance associated to the Xjt concentrations
    TF <- NULL       #material transfer function
    Rjt <- NULL      #R=sum(TF[jj, ,tt])
    sym <- c(1, 2, 5, 0, 6, 16, 17, 18, 15, 25)
    
    cat("
    DPSynth() requires a set of information to generate spectral intensities and concentrations
    as a function of the tilt angle (tilt=90 <=> surface normal to the analyzer axis).
    In particular it is supposed that the material is formed by a series of layers of different
    composition. The program requires:
    - the name of elements composing the material;
    - for each element the depth distribution: i.e. ElementName=Al2p  from=3, to=7 means that Al
       is distributed in a layer starting at a depth 3nm form the surface with thickness 4nm
    - for each element its concentration in the layer.
    DPSynth() requires the CoreLinesTable.txt to get element RSF and Kinetic energy
    needed to compute the element attenuation length.

    DPSynth() saves a FileName.RData in the  Z:/X/LAVORI/R/Analysis/test/DepthProf/   folder:
    data are ready to be loaded by XPSMaxEnt()")
    scan(n=1, quiet=TRUE)

    if (exists("ElemTab")==FALSE){ #if MaxEntData file not loaded
       cat("\n LOAD CORELINE TABLE")
       PathFile <- gfile(text = "CoreLinesTable.lib", type = c("open"),
               				      filter = list("RData files" = list(patterns= c("*.lib"))),
          	              multi = FALSE)
       OS<-Sys.info()
       if (OS["sysname"] != "Linux") {  # Windows and Mac OS systems
           PathFile <- gsub("/","\\", PathFile, fixed=TRUE)   #path/filename for linux,  path\\filename for windows
       }
       pp <- scan(PathFile, skip=0, sep="", what = list("character", "character", "numeric", "numeric", "numeric", "numeric"), quiet=TRUE)
       cat("\n ==> Coreline Table loaded \n")
       pp[[3]] <- as.numeric(pp[[3]])
       pp[[4]] <- as.numeric(pp[[4]])
       pp[[5]] <- as.numeric(pp[[5]])
       pp[[6]] <- as.numeric(pp[[6]])
       ElemTab <- as.data.frame(pp, stringsAsFactors = FALSE) # set it to FALSE
  	    names(ElemTab) <- c("Element", "Orbital", "BE", "KE", "RSF_K", "RSF_S")
    } else {
       ElemTab <- get("ElemTab", envir=.GlobalEnv)
    }


    SynthWin <- gwindow("ARXPS: DEPTH PROFILE ANALYSIS", visible=FALSE)
    size(SynthWin) <- c(400, 300)
    Synthgroup <- ggroup(horizontal=FALSE, container=SynthWin)
    NoteBK<-gnotebook(expand=TRUE, container = Synthgroup)
    DpInit <- ggroup(label="DepthPro INIT.", horizontal=FALSE, spacing=5, container=NoteBK)

    LayoutDP <- glayout(homogeneous=FALSE, spacing=5, container=DpInit)
    LayoutDP[1,1] <- glabel("N. elements", container=LayoutDP)
    Nelem <- LayoutDP[1,2] <- gedit("", initial.msg="N. elements", container = LayoutDP)
    LayoutDP[1,3] <- glabel("Mater. Density", container=LayoutDP)
    Density <- LayoutDP[1,4] <- gedit("", initial.msg="Mater. Density", container = LayoutDP)
    LayoutDP[2,1] <- glabel("N. angles", container=LayoutDP)
    Nang <- LayoutDP[2,2] <- gedit("", initial.msg="N. angles", container = LayoutDP)
    LayoutDP[2,3] <- glabel("Tilt increment", container=LayoutDP)
    TilInc <- LayoutDP[2,4] <- gedit("", initial.msg="Tilt increment", container = LayoutDP)
    LayoutDP[3,1] <- glabel("N. layers", container=LayoutDP)
    Nlay<- LayoutDP[3,2] <- gedit("", initial.msg="N. layers", container = LayoutDP)
    LayoutDP[3,3] <- glabel("Layer Thickness", container=LayoutDP)
    Lthick <- LayoutDP[3,4] <- gedit("", initial.msg="Layer Thickness", container = LayoutDP)
    ButtInit1 <- LayoutDP[4,1] <- gbutton(text = "   OK   ", handler = function(h, ...){
                                              NElmnts <<- 4
                                              NAngles <<- 10+1   #from 90 to 20 included
                                              Dteta <<- 7
                                              NLayers <<- 40        #40    20    10
                                              LThickness <<- 0.25    #0.25  0.5   1
                                              Density <<- 1080

                                              #NElmnts <<- as.double(svalue(Nelem))
                                              #NAngles <<- as.double(svalue(Nang))
                                              #Dteta <<- as.double(svalue(TilInc))
                                              #NLayers <<- as.double(svalue(Nlay))
                                              #LThickness <<- as.double(svalue(Lthick))
                                              Tilt <<- array(data=NA, dim=NAngles) #array of tilt angles
                                              Xjt <<- matrix(data=0, nrow=NElmnts, ncol=NAngles) #Spectral intensity of the elements
                                              Variance <<- matrix(data=0, nrow=NElmnts, ncol=NAngles) #Spectral intensity of the elements
                                              Nji <<- matrix(data=0, nrow=NElmnts, ncol=NLayers) #trend of the synthesized elements conc along the NLayers
                                              Ijt <<- matrix(data=0, nrow=NElmnts, ncol=NAngles) #trend of the synthesized elements intensity with tilt
                                              Depth <<- matrix(data=NA, nrow=NElmnts, ncol=2) #contains the depth distribution of elements
                                              for(ii in 1:NElmnts){
                                                  ElName[[ii]] <<- LayoutEl[ii,1] <- gedit("", initial.msg="C1s,Si2p...", container = LayoutEl)
                                                  tkconfigure(ElName[[ii]]$widget, width=12)
                                                  From[[ii]] <<- LayoutEl[ii,2] <- gedit("", initial.msg="From depth", container = LayoutEl)
                                                  tkconfigure(From[[ii]]$widget, width=12)
                                                  To[[ii]] <<- LayoutEl[ii,3] <- gedit("", initial.msg="To depth", container = LayoutEl)
                                                  tkconfigure(To[[ii]]$widget, width=12)
                                                  Conc[[ii]] <<- LayoutEl[ii,4] <- gedit("", initial.msg="Conc. %", container = LayoutEl)
                                                  tkconfigure(Conc[[ii]]$widget, width=12)
                                              }
                                              ButtEl1 <- LayoutDP[ii+1,1] <- gbutton(text = "SYNTHESIZE", handler = function(h, ...){
                                                                                      ElmtName <<- c("C1s", "N1s", "O1s", "Si2p")
                                                                                      Depth[ ,1] <<- c(0, 2, 5, 8)
                                                                                      Depth[ ,2] <<- c(2, 5, 8, 10)
                                                                                      ElmtConc <<- c(20, 30, 30, 20)  #in percentuale
                                                                                      for(jj in 1:NElmnts){
                                                                                      #    ElmtName[ii] <<- svalue(ElName[[ii]])
                                                                                      #    Depth[ii,1] <<- svalue(From[[ii]])
                                                                                      #    Depth[ii,2] <<- svalue(To[[ii]])
                                                                                      #    ElmtConc[ii] <<- svalue(Conc[[ii]])
                                                                                          Find.Elmt(jj)
                                                                                      }
                                                                                      SynthDP()
                                                                                   }, container = LayoutEl)
                                              ButtEl2 <- LayoutDP[ii+1,2] <- gbutton(text = "SAVE DATA ", handler = function(h, ...){
                                                                                      PathFileName <- gfile(text = " ", type = "save", initial.dir = getwd(),
                                                                                        multi = FALSE)
                                                                                      MaxEntData <- list(Variance=Variance, Concent=Xjt,
                                                                                        Intensity=Ijt, RSF=RSF, Att.L=AttLth,
                                                                                        Tilt=Tilt, N.Layers=NLayers, L.Thickness=LThickness, 
                                                                                        SelectedFName="XPS_DPsynth", SelectedCL=ElmtName)
cat("\n *************** MaxEnt List **************\n")
print(MaxEntData)
                                                                                      assign("MaxEntData", MaxEntData, envir=.GlobalEnv)  #save the xxx.Rdata XPSSample in the .GlobalEnv
                                                                                      save(list="MaxEntData", file=PathFileName, compress=TRUE)

                                                                                   }, container = LayoutEl)

                                              ButtEl3 <- LayoutDP[ii+1,3] <- gbutton(text = "   EXIT   ", handler = function(h, ...){
                                                                                      dispose(SynthWin)
                                                                                   }, container = LayoutEl)
                                              svalue(NoteBK) <<- 2
                                           }, container = LayoutDP)

    ButtInit2 <- LayoutDP[4,2] <- gbutton(text = "  EXIT  ", handler = function(h, ...){
                                              dispose(SynthWin)
                                           }, container = LayoutDP)

    ElemInit <- ggroup(label="Element INIT.", horizontal=FALSE, spacing=5, container=NoteBK)
    glabel("Element name must like C1s, O1s, Si2p... with orbital indicated", container=ElemInit)
    LayoutEl <- glayout(homogeneous=FALSE, spacing=5, container=ElemInit)



    svalue(NoteBK) <- 1
    visible(SynthWin) <- TRUE
}