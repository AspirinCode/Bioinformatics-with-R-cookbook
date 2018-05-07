
# Overcoming batch effects in expression data
As processing samples in different batches, the systematic errors will be caused, which is named batch effects. Systematic errors can be reduced, to some extent, by careful, experimental design but cannot be eliminated completely unless the study is performed under a single batch. It is difficult to combine data from different batches due to batch effects, which makes the statistical analysis of the data unsuccessful. Therefore, if we do preprocessing before the batches are combined, the analysis will be easy. <br>
处理不同批次的样品会导致系统误差，这种称为批量效应。 在一定程度上，通过仔细严谨的实验设计可以减少系统误差，但是除非单一批次进行研究，否则不能完全消除。由于批处理效应而难以合并不同批次的数据，从而使得数据的统计分析不顺利。因此，如果我们在批次合并之前进行预处理，分析将变得简单许多。<br>

In this section, we need a dataset that shows the batch effect without any preprocessing. Take the `bladderbatch` data as an example, which consists of five batches. The data is a part of the `bladderbatch` package. We will use `sva` library as well.<br>
在本节中，我们需要一个显示批处理效果的数据集，而无需任何预处理。 以`bladderbatch`数据为例，它由五批次组成。这些数据是`bladderbatch`软件包的一部分。 我们也会使用`sva`库。

1、Install and load all the required libraries.


```R
source("http://bioconductor.org/biocLite.R")
```

    Bioconductor version 3.6 (BiocInstaller 1.28.0), ?biocLite for help
    


```R
biocLite(c("sva", "bladderbatch"))
```

    BioC_mirror: https://bioconductor.org
    Using Bioconductor 3.6 (BiocInstaller 1.28.0), R 3.4.3 (2017-11-30).
    Installing package(s) 'sva', 'bladderbatch'
    

    package 'sva' successfully unpacked and MD5 sums checked
    
    The downloaded binary packages are in
    	C:\Users\Master1\AppData\Local\Temp\Rtmpy4D4wz\downloaded_packages
    

    installing the source package 'bladderbatch'
    
    Old packages: 'GenomicRanges'
    


```R
library(sva)
```

    Loading required package: mgcv
    Loading required package: nlme
    This is mgcv 1.8-23. For overview type 'help("mgcv-package")'.
    Loading required package: genefilter
    Loading required package: BiocParallel
    


```R
library(bladderbatch)
```

    Loading required package: Biobase
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
    
    Welcome to Bioconductor
    
        Vignettes contain introductory material; view with
        'browseVignettes()'. To cite Bioconductor, see
        'citation("Biobase")', and for packages 'citation("pkgname")'.
    
    

2、Load `bladderdata`. And extract the expression matrix and pheno data from it.


```R
data(bladderdata)
```


```R
pheno <- pData(bladderEset)
edata <- exprs(bladderEset)
pheno
```


<table>
<thead><tr><th></th><th scope=col>sample</th><th scope=col>outcome</th><th scope=col>batch</th><th scope=col>cancer</th></tr></thead>
<tbody>
	<tr><th scope=row>GSM71019.CEL</th><td> 1      </td><td>Normal  </td><td>3       </td><td>Normal  </td></tr>
	<tr><th scope=row>GSM71020.CEL</th><td> 2      </td><td>Normal  </td><td>2       </td><td>Normal  </td></tr>
	<tr><th scope=row>GSM71021.CEL</th><td> 3      </td><td>Normal  </td><td>2       </td><td>Normal  </td></tr>
	<tr><th scope=row>GSM71022.CEL</th><td> 4      </td><td>Normal  </td><td>3       </td><td>Normal  </td></tr>
	<tr><th scope=row>GSM71023.CEL</th><td> 5      </td><td>Normal  </td><td>3       </td><td>Normal  </td></tr>
	<tr><th scope=row>GSM71024.CEL</th><td> 6      </td><td>Normal  </td><td>3       </td><td>Normal  </td></tr>
	<tr><th scope=row>GSM71025.CEL</th><td> 7      </td><td>Normal  </td><td>2       </td><td>Normal  </td></tr>
	<tr><th scope=row>GSM71026.CEL</th><td> 8      </td><td>Normal  </td><td>2       </td><td>Normal  </td></tr>
	<tr><th scope=row>GSM71028.CEL</th><td> 9      </td><td>sTCC+CIS</td><td>5       </td><td>Cancer  </td></tr>
	<tr><th scope=row>GSM71029.CEL</th><td>10      </td><td>sTCC-CIS</td><td>2       </td><td>Cancer  </td></tr>
	<tr><th scope=row>GSM71030.CEL</th><td>11      </td><td>sTCC-CIS</td><td>5       </td><td>Cancer  </td></tr>
	<tr><th scope=row>GSM71031.CEL</th><td>12      </td><td>sTCC-CIS</td><td>2       </td><td>Cancer  </td></tr>
	<tr><th scope=row>GSM71032.CEL</th><td>13      </td><td>sTCC+CIS</td><td>5       </td><td>Cancer  </td></tr>
	<tr><th scope=row>GSM71033.CEL</th><td>14      </td><td>sTCC-CIS</td><td>2       </td><td>Cancer  </td></tr>
	<tr><th scope=row>GSM71034.CEL</th><td>15      </td><td>sTCC+CIS</td><td>5       </td><td>Cancer  </td></tr>
	<tr><th scope=row>GSM71035.CEL</th><td>16      </td><td>sTCC+CIS</td><td>5       </td><td>Cancer  </td></tr>
	<tr><th scope=row>GSM71036.CEL</th><td>17      </td><td>sTCC-CIS</td><td>2       </td><td>Cancer  </td></tr>
	<tr><th scope=row>GSM71037.CEL</th><td>18      </td><td>mTCC    </td><td>1       </td><td>Cancer  </td></tr>
	<tr><th scope=row>GSM71038.CEL</th><td>19      </td><td>sTCC+CIS</td><td>5       </td><td>Cancer  </td></tr>
	<tr><th scope=row>GSM71039.CEL</th><td>20      </td><td>mTCC    </td><td>1       </td><td>Cancer  </td></tr>
	<tr><th scope=row>GSM71040.CEL</th><td>21      </td><td>mTCC    </td><td>2       </td><td>Cancer  </td></tr>
	<tr><th scope=row>GSM71041.CEL</th><td>22      </td><td>mTCC    </td><td>1       </td><td>Cancer  </td></tr>
	<tr><th scope=row>GSM71042.CEL</th><td>23      </td><td>sTCC-CIS</td><td>2       </td><td>Cancer  </td></tr>
	<tr><th scope=row>GSM71043.CEL</th><td>24      </td><td>sTCC+CIS</td><td>5       </td><td>Cancer  </td></tr>
	<tr><th scope=row>GSM71044.CEL</th><td>25      </td><td>sTCC-CIS</td><td>2       </td><td>Cancer  </td></tr>
	<tr><th scope=row>GSM71045.CEL</th><td>26      </td><td>sTCC-CIS</td><td>2       </td><td>Cancer  </td></tr>
	<tr><th scope=row>GSM71046.CEL</th><td>27      </td><td>sTCC+CIS</td><td>5       </td><td>Cancer  </td></tr>
	<tr><th scope=row>GSM71047.CEL</th><td>28      </td><td>mTCC    </td><td>1       </td><td>Cancer  </td></tr>
	<tr><th scope=row>GSM71048.CEL</th><td>29      </td><td>mTCC    </td><td>1       </td><td>Cancer  </td></tr>
	<tr><th scope=row>GSM71049.CEL</th><td>30      </td><td>sTCC-CIS</td><td>2       </td><td>Cancer  </td></tr>
	<tr><th scope=row>GSM71050.CEL</th><td>31      </td><td>mTCC    </td><td>1       </td><td>Cancer  </td></tr>
	<tr><th scope=row>GSM71051.CEL</th><td>32      </td><td>mTCC    </td><td>1       </td><td>Cancer  </td></tr>
	<tr><th scope=row>GSM71052.CEL</th><td>33      </td><td>mTCC    </td><td>1       </td><td>Cancer  </td></tr>
	<tr><th scope=row>GSM71053.CEL</th><td>34      </td><td>sTCC+CIS</td><td>5       </td><td>Cancer  </td></tr>
	<tr><th scope=row>GSM71054.CEL</th><td>35      </td><td>mTCC    </td><td>1       </td><td>Cancer  </td></tr>
	<tr><th scope=row>GSM71055.CEL</th><td>36      </td><td>sTCC-CIS</td><td>2       </td><td>Cancer  </td></tr>
	<tr><th scope=row>GSM71056.CEL</th><td>37      </td><td>sTCC-CIS</td><td>2       </td><td>Cancer  </td></tr>
	<tr><th scope=row>GSM71058.CEL</th><td>38      </td><td>sTCC-CIS</td><td>2       </td><td>Cancer  </td></tr>
	<tr><th scope=row>GSM71059.CEL</th><td>39      </td><td>sTCC-CIS</td><td>2       </td><td>Cancer  </td></tr>
	<tr><th scope=row>GSM71060.CEL</th><td>40      </td><td>mTCC    </td><td>1       </td><td>Cancer  </td></tr>
	<tr><th scope=row>GSM71061.CEL</th><td>41      </td><td>sTCC+CIS</td><td>5       </td><td>Cancer  </td></tr>
	<tr><th scope=row>GSM71062.CEL</th><td>42      </td><td>sTCC+CIS</td><td>5       </td><td>Cancer  </td></tr>
	<tr><th scope=row>GSM71063.CEL</th><td>43      </td><td>sTCC+CIS</td><td>5       </td><td>Cancer  </td></tr>
	<tr><th scope=row>GSM71064.CEL</th><td>44      </td><td>sTCC-CIS</td><td>2       </td><td>Cancer  </td></tr>
	<tr><th scope=row>GSM71065.CEL</th><td>45      </td><td>sTCC-CIS</td><td>5       </td><td>Cancer  </td></tr>
	<tr><th scope=row>GSM71066.CEL</th><td>46      </td><td>mTCC    </td><td>1       </td><td>Cancer  </td></tr>
	<tr><th scope=row>GSM71067.CEL</th><td>47      </td><td>sTCC-CIS</td><td>5       </td><td>Cancer  </td></tr>
	<tr><th scope=row>GSM71068.CEL</th><td>48      </td><td>sTCC+CIS</td><td>5       </td><td>Cancer  </td></tr>
	<tr><th scope=row>GSM71069.CEL</th><td>49      </td><td>Biopsy  </td><td>4       </td><td>Biopsy  </td></tr>
	<tr><th scope=row>GSM71070.CEL</th><td>50      </td><td>Biopsy  </td><td>4       </td><td>Biopsy  </td></tr>
	<tr><th scope=row>GSM71071.CEL</th><td>51      </td><td>Biopsy  </td><td>5       </td><td>Biopsy  </td></tr>
	<tr><th scope=row>GSM71072.CEL</th><td>52      </td><td>Biopsy  </td><td>5       </td><td>Biopsy  </td></tr>
	<tr><th scope=row>GSM71073.CEL</th><td>53      </td><td>Biopsy  </td><td>5       </td><td>Biopsy  </td></tr>
	<tr><th scope=row>GSM71074.CEL</th><td>54      </td><td>Biopsy  </td><td>5       </td><td>Biopsy  </td></tr>
	<tr><th scope=row>GSM71075.CEL</th><td>55      </td><td>Biopsy  </td><td>4       </td><td>Biopsy  </td></tr>
	<tr><th scope=row>GSM71076.CEL</th><td>56      </td><td>Biopsy  </td><td>4       </td><td>Biopsy  </td></tr>
	<tr><th scope=row>GSM71077.CEL</th><td>57      </td><td>Biopsy  </td><td>4       </td><td>Biopsy  </td></tr>
</tbody>
</table>



3、We can see that the first eight samples are normal cells but split into two batches(batch number 2 and 3). Use these eight samples to demonstrate how to remove the batch effects.


```R
mydata <- bladderEset[,sampleNames(bladderEset)[1:8]]
mydata
```


    ExpressionSet (storageMode: lockedEnvironment)
    assayData: 22283 features, 8 samples 
      element names: exprs, se.exprs 
    protocolData: none
    phenoData
      sampleNames: GSM71019.CEL GSM71020.CEL ... GSM71026.CEL (8 total)
      varLabels: sample outcome batch cancer
      varMetadata: labelDescription
    featureData: none
    experimentData: use 'experimentData(object)'
    Annotation: hgu133a 


4、Use `arrayQualityMetrics` to perform a quality check on the data, and then we will have a look at the batch effects.


```R
library(arrayQualityMetrics)
```


```R
arrayQualityMetrics(mydata, outdir="D:/Try-practice/Chapter 5/qc_be")
```

    The directory 'D:/Try-practice/Chapter 5/qc_be' has been created.
    Warning message in svgStyleAttributes(style, svgdev):
    "Removing non-SVG style attribute name(s): subscripts, group.number, group.value"Warning message in svgStyleAttributes(style, svgdev):
    "Removing non-SVG style attribute name(s): subscripts, group.number, group.value"Warning message in svgStyleAttributes(style, svgdev):
    "Removing non-SVG style attribute name(s): subscripts, group.number, group.value"Warning message in svgStyleAttributes(style, svgdev):
    "Removing non-SVG style attribute name(s): subscripts, group.number, group.value"Warning message in svgStyleAttributes(style, svgdev):
    "Removing non-SVG style attribute name(s): subscripts, group.number, group.value"Warning message in svgStyleAttributes(style, svgdev):
    "Removing non-SVG style attribute name(s): subscripts, group.number, group.value"Warning message in svgStyleAttributes(style, svgdev):
    "Removing non-SVG style attribute name(s): subscripts, group.number, group.value"Warning message in svgStyleAttributes(style, svgdev):
    "Removing non-SVG style attribute name(s): subscripts, group.number, group.value"(loaded the KernSmooth namespace)
    

5、Take a look at the heatmap and other plots produced for the samples to check for the batch effect. In the filefolder, there are files called 'hm' and 'ma'.


```R
library(hexbin)
```


```R
browseURL("D:/Try-practice/Chapter 5/qc_be/index.html")
```

5、Create the model matrix for the dataset (note that only the first and third columns have been used from the model matrix as the data has only one condition).<br>为数据集创建模型矩阵(请注意，从模型矩阵中只使用第一列和第三列，因为数据只有一个条件).


```R
mod1 <- model.matrix(~as.factor(cancer), data=pData(mydata))[,c(1,3)]
mod1
```


<table>
<thead><tr><th></th><th scope=col>(Intercept)</th><th scope=col>as.factor(cancer)Normal</th></tr></thead>
<tbody>
	<tr><th scope=row>GSM71019.CEL</th><td>1</td><td>1</td></tr>
	<tr><th scope=row>GSM71020.CEL</th><td>1</td><td>1</td></tr>
	<tr><th scope=row>GSM71021.CEL</th><td>1</td><td>1</td></tr>
	<tr><th scope=row>GSM71022.CEL</th><td>1</td><td>1</td></tr>
	<tr><th scope=row>GSM71023.CEL</th><td>1</td><td>1</td></tr>
	<tr><th scope=row>GSM71024.CEL</th><td>1</td><td>1</td></tr>
	<tr><th scope=row>GSM71025.CEL</th><td>1</td><td>1</td></tr>
	<tr><th scope=row>GSM71026.CEL</th><td>1</td><td>1</td></tr>
</tbody>
</table>



6、define the batches


```R
batch <- pData(mydata)$batch
```

7、Extract the expression matrix from the expression set object myData, where the batch effect has to be removed.<br>(从表达式集对象mydata中提取表达式矩阵，其中必须删除批处理效果)


```R
edata <- exprs(mydata)
```

8、When we have all the objects ready, run the `ComBat` function.


```R
combat_edata <- ComBat(dat=edata, batch=batch, mod=mod1, par.prior=TRUE) # numCovs=NULL?
```

    Found2batches
    Adjusting for0covariate(s) or covariate level(s)
    

    Standardizing Data across genes
    

    Fitting L/S model and finding priors
    Finding parametric adjustments
    Adjusting the Data
    
    

9、Create an ExpressionSet object with everything as the original input data, except the expression matrix—which is replaced by the matrix received as a result of the `ComBat` function in the last step.


```R
mydata2 <- mydata
exprs(mydata2) <- combat_edata
```

10、Rerun the `arrayQualityMetrics` function to check for the elimination of batch effects with this new object as the input.


```R
arrayQualityMetrics(mydata2, outdir="D:/Try-practice/Chapter 5/qc_nbe")
```

    The report will be written into directory 'D:/Try-practice/Chapter 5/qc_nbe'. 
    Warning message in svgStyleAttributes(style, svgdev):
    "Removing non-SVG style attribute name(s): subscripts, group.number, group.value"Warning message in svgStyleAttributes(style, svgdev):
    "Removing non-SVG style attribute name(s): subscripts, group.number, group.value"Warning message in svgStyleAttributes(style, svgdev):
    "Removing non-SVG style attribute name(s): subscripts, group.number, group.value"Warning message in svgStyleAttributes(style, svgdev):
    "Removing non-SVG style attribute name(s): subscripts, group.number, group.value"Warning message in svgStyleAttributes(style, svgdev):
    "Removing non-SVG style attribute name(s): subscripts, group.number, group.value"Warning message in svgStyleAttributes(style, svgdev):
    "Removing non-SVG style attribute name(s): subscripts, group.number, group.value"Warning message in svgStyleAttributes(style, svgdev):
    "Removing non-SVG style attribute name(s): subscripts, group.number, group.value"Warning message in svgStyleAttributes(style, svgdev):
    "Removing non-SVG style attribute name(s): subscripts, group.number, group.value"


```R
browseURL("D:/Try-practice/Chapter 5/qc_nbe/index.html")
```

For other optional input arguments, typing ?ComBat.
