
# reading CEL files
To have a better understanding of how to analyze microarrays data, we will start with the reading and loading of microarray data, and then follow by its preprocessing, analysis, mining, and finally, its visualization using R. Lastly, we learn about the solution related to the biological use of the data.<br>
The CEL file contains all the intensity-related information of the pixels(including the intensity itself, the standard deviation, the number of pixels, and other meta information) on the array. Every experiment usually has more than one sample and replicates, so there will be one CEL file present for each sample or replicated file. The CEL files must be read to get the data in a workable format in the R workspace. To do the first step of microarray data analysis, as an example, we will use the NCBI GEO database for a simple dataset on breast cancer. Furthermore, we will need some R packages(`affy`).<br>
为了更好地理解如何分析微阵列数据，我们将从微阵列数据的读取和加载开始，然后进行预处理，分析，挖掘，最后使用R进行可视化。最后，我们了解与数据的生物学使用有关的解决方案。 <br>
“CEL”文件包含像素的所有强度相关信息（包括强度本身，标准偏差，像素的数量和其他元信息）。每个实验通常都有多个样本并进行复制，所以每个样本或复制文件都会有一个“CEL”文件。必须读取“CEL”文件才能在R工作区中以可行的格式获取数据。为了完成微阵列数据分析的第一步，我们将使用NCBI GEO数据库作为乳腺癌的简单数据集作为例子。此外，我们需要一些R包（`affy`）。

At first, download `GSE24460` from [NCBI GEO](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE24460). This will give a file named `GSE24460_RAW.tar` in desired directory. But I run in windows operating sysytem, then use the `type = "win.binary"`. 
Unzip this file to get CEL files. The files have also been provided with the code files the sub directory `<GSE24460_RAW>`.<br>
首先在NCBI GEO数据库里下载‘GSE24460’文件，下载完成，在目标路径下获得的是文件压缩包，需要解压缩，获得‘.CEL’格式的文件。


```R
source("http://www.bioconductor.org/biocLite.R")
```

    Bioconductor version 3.6 (BiocInstaller 1.28.0), ?biocLite for help
    


```R
biocLite("affy")
```

    BioC_mirror: https://bioconductor.org
    Using Bioconductor 3.6 (BiocInstaller 1.28.0), R 3.4.3 (2017-11-30).
    Installing package(s) 'affy'
    

    package 'affy' successfully unpacked and MD5 sums checked
    
    The downloaded binary packages are in
    	C:\Users\Master1\AppData\Local\Temp\RtmpcHoKSH\downloaded_packages
    

    Old packages: 'GenomicRanges'
    


```R
biocLite("affydata")
```

    BioC_mirror: https://bioconductor.org
    Using Bioconductor 3.6 (BiocInstaller 1.28.0), R 3.4.3 (2017-11-30).
    Installing package(s) 'affydata'
    installing the source package 'affydata'
    
    Old packages: 'GenomicRanges'
    


```R
biocLite("hgu133a2cdf")
```

    BioC_mirror: https://bioconductor.org
    Using Bioconductor 3.6 (BiocInstaller 1.28.0), R 3.4.3 (2017-11-30).
    Installing package(s) 'hgu133a2cdf'
    installing the source package 'hgu133a2cdf'
    
    Old packages: 'GenomicRanges'
    


```R
library(affy)
```

    Loading required package: BiocGenerics
    Loading required package: parallel
    
    Attaching package: 'BiocGenerics'
    
    The following objects are masked from 'package:parallel':
    
        clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
        clusterExport, clusterMap, parApply, parCapply, parLapply,
        parLapplyLB, parRapply, parSapply, parSapplyLB
    
    The following objects are masked from 'package:stats':
    
        IQR, mad, sd, var, xtabs
    
    The following objects are masked from 'package:base':
    
        anyDuplicated, append, as.data.frame, cbind, colMeans, colnames,
        colSums, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
        grepl, intersect, is.unsorted, lapply, lengths, Map, mapply, match,
        mget, order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
        rbind, Reduce, rowMeans, rownames, rowSums, sapply, setdiff, sort,
        table, tapply, union, unique, unsplit, which, which.max, which.min
    
    Loading required package: Biobase
    Welcome to Bioconductor
    
        Vignettes contain introductory material; view with
        'browseVignettes()'. To cite Bioconductor, see
        'citation("Biobase")', and for packages 'citation("pkgname")'.
    
    


```R
library(affydata)
```

         Package    LibPath                                        Item      
    [1,] "affydata" "C:/Users/Master1/Documents/R/win-library/3.4" "Dilution"
         Title                        
    [1,] "AffyBatch instance Dilution"
    


```R
library(hgu133a2cdf)
```


```R
mydata <- ReadAffy(celfile.path= "D:/Try-practice/Chapter 5/GSE24460_RAW/")
```


```R
mydata <- ReadAffy(filenames="D:/Try-practice/Chapter 5/GSE24460_RAW/GSM602658_MCF71.CEL")
```


```R
mydata
```


    AffyBatch object
    size of arrays=732x732 features (17 kb)
    cdf=HG-U133A_2 (22277 affyids)
    number of samples=1
    number of genes=22277
    annotation=hgu133a2
    notes=


The `ReadAffy` function is a wrapper for another function called `read.affybatch` that allows the reading of the CEL file along with other relevant data into an `AffyBatch` object. The `ReadAffy` function gets the list of all of the CEL files (even the compressed TAR files) in the current working directory and passes it to the `read.affybatch` function to create the `AffyBatch` object from these files. It is possible to provide arguments for other information such as `phenoData` and specific CEL filenames and paths to another directory with the same function. When no arguments are passed into the `ReadAffy` function (as `ReadAffy()`), the function reads all the CEL files in the working directory.<br>
ReadAffy函数是另一个称为read.affybatch的函数的包装器，它允许将CEL文件和其他相关数据读入到AffyBatch对象中。ReadAffy函数获取当前工作目录中所有CEL文件（甚至是压缩的TAR文件）的列表，并将其传递给read.affybatch函数从而以这些文件创建AffyBatch对象。可以为其他信息提供参数，例如phenoData和特定的CEL文件名以及具有相同功能的另一个目录的路径。当没有参数传入ReadAffy函数时(如ReadAffy())，该函数将读取工作目录中的所有CEL文件。

There is an example:


```R
 if(require(affydata)){
    celpath <- system.file("celfiles", package="affydata")
     fns <- list.celfiles(path=celpath,full.names=TRUE)
     cat("Reading files:\n",paste(fns,collapse="\n"),"\n")
     abatch <- ReadAffy(filenames=fns[1])                # read a binary celfile
     abatch <- ReadAffy(filenames=fns[2])                # read a text celfile
     abatch <- ReadAffy(celfile.path=celpath)            # read all files in that dir 
}
```

    Reading files:
     C:/Users/Master1/Documents/R/win-library/3.4/affydata/celfiles/binary.cel
    C:/Users/Master1/Documents/R/win-library/3.4/affydata/celfiles/text.cel 
    
