
# Checking the quality of data
The quality of data can be affected at each step of the microarray experiment pipeline. Quality-related problems could be produced from hybridization due to uneven fluorescence on the chip that causes variable intensity distributions, and a nonspecific binding or other biological/technical reasons creating background noise in the data. Another possible situation can be an inappropriate experimental design that may affect the dataset as a whole. Using such data will
result in wrongful or inconclusive inference during data analysis. Therefore, as analysts we must ensure the data quality before we begin the data analysis. This is achieved by looking for outlying arrays, distributions within arrays, batch effects, and so on. There are various analyses and diagnostic plots that can be used to compute these measures that explain the quality of array data under analysis.<br>

数据的质量可能会影响微阵列实验过程的每一步。质量相关问题可能源于芯片上荧光不均导致的杂交，从而导致可变强度分布。非特异性结合或其他生物/技术原因可能会在数据中产生背景噪音。另一种可能的情况可能是不恰当的实验设计，可能会影响整个数据集。使用这些数据会导致数据分析过程中的错误或不确定的推断。因此，作为分析师，我们必须在开始数据分析之前确保数据质量。可以通过查找外围数组，数组中的分布，批处理效果等来实现。有多种分析和诊断图可用于计算这些解释分析中阵列数据质量的度量。配方将解释数据质量检查的各种诊断步骤。


1、Install and load the `arrayQualityMetrics` library from the Bioconductor repository, and load the `mydata`.


```R
source("http://www.bioconductor.org/biocLite.R")
```

    Bioconductor version 3.6 (BiocInstaller 1.28.0), ?biocLite for help
    


```R
biocLite("arrayQualityMetrics")
```

    BioC_mirror: https://bioconductor.org
    Using Bioconductor 3.6 (BiocInstaller 1.28.0), R 3.4.3 (2017-11-30).
    Installing package(s) 'arrayQualityMetrics'
    

    package 'arrayQualityMetrics' successfully unpacked and MD5 sums checked
    
    The downloaded binary packages are in
    	C:\Users\Master1\AppData\Local\Temp\RtmpURrkBW\downloaded_packages
    

    Old packages: 'GenomicRanges'
    


```R
library(arrayQualityMetrics)
```


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
mydata <- ReadAffy(celfile.path= "D:/Try-practice/Chapter 5/GSE24460_RAW/")
mydata
```


    AffyBatch object
    size of arrays=732x732 features (18 kb)
    cdf=HG-U133A_2 (22277 affyids)
    number of samples=4
    number of genes=22277
    annotation=hgu133a2
    notes=


2、Use the `arrayQualityMetrics` function to create plots to assess the data quality. `ArrayQualityMetrics` function takes the data and performs different types of checks including measuring between arrays distances, Principal Component Analysis(PCA), density plots, MA plots, and RNA degradation plots. The outlier detection among the arrays is reflected by measurements such as the distance between the arrays, boxplots, and MA plots. A detailed description is available in the HTML file produced.<br>
ArrayQualityMetrics函数获取数据并执行不同类型的检查，包括在阵列距离的测量，主成分分析（PCA），密度图，MA图和RNA降解图。阵列之间的距离，箱形图和MA图之类的测量等能反映阵列之间的异常值。在生成的HTML文件中有详细的描述。


```R
arrayQualityMetrics(mydata, outdir="D:/Try-practice/Chapter 5/quality_assesment")
```

    The directory 'D:/Try-practice/Chapter 5/quality_assesment' has been created.
    Warning message in svgStyleAttributes(style, svgdev):
    "Removing non-SVG style attribute name(s): subscripts, group.number, group.value"Warning message in svgStyleAttributes(style, svgdev):
    "Removing non-SVG style attribute name(s): subscripts, group.number, group.value"Warning message in svgStyleAttributes(style, svgdev):
    "Removing non-SVG style attribute name(s): subscripts, group.number, group.value"Warning message in svgStyleAttributes(style, svgdev):
    "Removing non-SVG style attribute name(s): subscripts, group.number, group.value"Warning message in KernSmooth::bkde2D(x, gridsize = nbin, bandwidth = bandwidth):
    "Binning grid too coarse for current (small) bandwidth: consider increasing 'gridsize'"Warning message in KernSmooth::bkde2D(x, gridsize = nbin, bandwidth = bandwidth):
    "Binning grid too coarse for current (small) bandwidth: consider increasing 'gridsize'"Warning message in KernSmooth::bkde2D(x, gridsize = nbin, bandwidth = bandwidth):
    "Binning grid too coarse for current (small) bandwidth: consider increasing 'gridsize'"Warning message in KernSmooth::bkde2D(x, gridsize = nbin, bandwidth = bandwidth):
    "Binning grid too coarse for current (small) bandwidth: consider increasing 'gridsize'"Warning message in KernSmooth::bkde2D(x, gridsize = nbin, bandwidth = bandwidth):
    "Binning grid too coarse for current (small) bandwidth: consider increasing 'gridsize'"Warning message in KernSmooth::bkde2D(x, gridsize = nbin, bandwidth = bandwidth):
    "Binning grid too coarse for current (small) bandwidth: consider increasing 'gridsize'"Warning message in KernSmooth::bkde2D(x, gridsize = nbin, bandwidth = bandwidth):
    "Binning grid too coarse for current (small) bandwidth: consider increasing 'gridsize'"Warning message in KernSmooth::bkde2D(x, gridsize = nbin, bandwidth = bandwidth):
    "Binning grid too coarse for current (small) bandwidth: consider increasing 'gridsize'"


```R
biocLite("hexbin")
```

    BioC_mirror: https://bioconductor.org
    Using Bioconductor 3.6 (BiocInstaller 1.28.0), R 3.4.3 (2017-11-30).
    Installing package(s) 'hexbin'
    

    package 'hexbin' successfully unpacked and MD5 sums checked
    
    The downloaded binary packages are in
    	C:\Users\Master1\AppData\Local\Temp\RtmpURrkBW\downloaded_packages
    

    Old packages: 'GenomicRanges'
    


```R
library(hexbin)
```

3、Open the HTML file in a browser


```R
browseURL("D:/Try-practice/Chapter 5/quality_assesment/index.html")
```

4、create these plots and assessments individually, like an MA plot by `MAplot()`. 


```R
MAplot(mydata, pairs=TRUE, plot.method="smoothScatter")
```


[MAplot](https://github.com/Chengshu21/Chapter-5-Analyzing-Microarray-Data-with-R/blob/master/md/pic/MAplot(Checking%20the%20quality%20of%20data).png)


5、plot the log densities by `plotDensity.AffyBatch`


```R
plotDensity.AffyBatch(mydata)
```


[indensity plot](https://github.com/Chengshu21/Chapter-5-Analyzing-Microarray-Data-with-R/blob/master/md/pic/indensity%20plot(Checking%20the%20quality%20of%20data).png)


6、create the boxplots by using the `boxplot` function in the `AffyBatch` object


```R
boxplot(mydata)
```


[boxplot](https://github.com/Chengshu21/Chapter-5-Analyzing-Microarray-Data-with-R/blob/master/md/pic/boxplot（checkingthequalityof%20data）.png)


7、get the RNA degradation plot by using the `AffyRNAdeg` function and then `plotAffyRNAdeg`


```R
rnaDeg <- AffyRNAdeg(mydata)
plotAffyRNAdeg(rnaDeg)
```


[RNA degradation](https://github.com/Chengshu21/Chapter-5-Analyzing-Microarray-Data-with-R/blob/master/md/pic/RNA%20degradation%20plot(Checking%20the%20quality%20of%20data).png)


8、Check the details of the `rnaDeg` object


```R
summaryAffyRNAdeg(rnaDeg)
```


<table>
<thead><tr><th></th><th scope=col>GSM602658_MCF71.CEL</th><th scope=col>GSM602659_MCF72.CEL</th><th scope=col>GSM602660_MCF7226ng.CEL</th><th scope=col>GSM602661_MCF7262ng.CEL</th></tr></thead>
<tbody>
	<tr><th scope=row>slope</th><td>2.19e+00</td><td>2.22e+00</td><td>2.73e+00</td><td>1.74e+00</td></tr>
	<tr><th scope=row>pvalue</th><td>2.87e-13</td><td>5.01e-13</td><td>3.40e-11</td><td>1.51e-06</td></tr>
</tbody>
</table>



More introduction of `ArrayQualityMetrics` function by typing `?ArrayQualityMetrics`.<br>
Turning to [ArrayQualityMetrics](), understanding of MA plot, intensity plot, boxplots and RNA degradation plot.
