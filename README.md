# RxpsG
Processing tool for X-ray Photoelectron Spectroscopy Data

RxpsG installation instructions
    1. Install R and Rstudio on your computer;
    
    2. Run Rstudio and open the Tools menu -> Global Options;
    
    3. For Windows Systems set the R version path (es. C:\Program Files\R\R-3.6.0);
    
    4. Install the devtools package from CRAN;
    
    5. Load the devtools package using the command: library(devtools);
    
    6. Use install_github("GSperanza/RxpsG");
    
    
Manual Installation:

    1. Click on the RxpsG_xx.xx.tar.gz package and download. Exit the unzipping procedure if it starts automatically.
    
    2. Control in the Dowloads folder the RxpsG_xx.xx.tar.gz package is present (it could be the .gz extension is lacking do not worry).
    
    3. Run RStudio
    
    4. Install all the RxpsG dependencies using: install.packages(pkgs=c("baseline", "digest", "FME", "gWidgets2 (>= 1.0-7)", "gWidgets2tcltk", "latticeExtra",     "memoise", "minpack.lm", "NORMT3", "signal", "sm", "SparseM", "tkrplot", "wavelets"))
    
    5. Under RStudio open the Tools menu -> Install Packages -> Package Archive File ( *.tar.gz )
    
    6. Install RxpsG browsing the downloaded Rxpsg.xx.xx.tar.gz file
