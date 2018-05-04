# Installing Libraries
   Libraries in R are packages that have functions written to serve specific purposes. We must have these libraries installed in the system as well as loaded in the R session in order to be able to use them. I have installed R-3.4.1 for win10. And then I study how to install R packages in R. For more details on the R program and its installation, visit [R](https://www.r-project.org/).
### Three ways to install packages:
   1、By selecting the menu, of course, it needs a internet connection:
  * From the `Packages` menu in the toolbar, select `Install package(s)....`
  * We can choose the packages we want in the pop-up dialog box, and then OK.
   
   2、By using the command: 
  
  `> install.packages("packagename","dir")`
  
   If we want to know more about how to use the command, we can type `?install.packages` in R console.  
   But in the jupyter notebook, I cannot directly run `install.packages("biomaRt")`, we can find a list CRAN in [R-CRAN](https://cran.r-project.org/mirrors.html) and then I try the following command to choose a CRAN mirror from Tsinghua University before installing packages:
  
  > options(repos=structure(c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")))  
  > install.packages("dplyr")

   During studying the bioinformatics, I always choose the [Bioconductor](http://www.bioconductor.org/) to install packages:
  
  > source("http://www.bioconductor.org/biocLite.R") <br>
  > biocLite("biomaRt") <br>
  > library(biomaRt) <br>
  
   library(packagename) is used to load the installed libraries/packages in R. The BioCmirror can be changed by:
  
  > options(BioC_mirror="http://mirrors.ustc.edu.cn/bioc/")
  
  3、By local files to install:
  * Download the corresponding package compressed file: In windows, unix, macOS operating system installation file extension is not the same: 
  1) linux environment compile and run: tar.gz file
  2) windows environment compile and run: .zip file
  3) MacOSg environment compile and run: .tgz file
  
  * And then we can type the following commmand to install:
  
  > install.packages("path/to/mypackage.tar.gz", repos = NULL, type="source")
  
  * Or from the R toolbar, choose `Packages` menu, select`Install package(s) from local files...`
  
  ## Practice in Jupyter notebok(R)
  The related commmand is in [practice-installing-libraries](https://github.com/Chengshu21/bio-start-with-R/blob/master/R/practice-installing%20libraries.r).
