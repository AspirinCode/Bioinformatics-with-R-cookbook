
# Finding the differentially expressed genes(查找差异表达的基因)
At the genome level, we know cells have similar genes. But what makes them different? As we all known, only a fraction of a genome is expressed in each cell, and this phenomenon of selective expression of genes based on cell types is the baseline behind the concept of differential gene expression. Thus, it is important to find which genes show differential expression in a particular cell. This will be performed by comparing the cell under study with a reference, usually called control. We will present how to find the DE genes for a cell based on the expression levels of the control and treatment cells.<br>
在基因组水平上，我们知道细胞具有相似的基因。但是因为什么才使得他们有所不同？众所周知，每个细胞只有一小部分基因组表达，这种基于细胞类型选择性表达基因的现象是差异基因表达概念的基础。因此，找出哪些基因在特定细胞中表现出差异表达很重要。这将通过比较实验组细胞和对照组细胞来进行参考。我们将根据对照组和实验组细胞的表达水平来演示如何找到细胞的DE基因。<br>

Required:<br>
1、The normalized expression data for treatment and control samples. Take the `quantile` normalized data as an example.实验组和对照组样本的标准化表达数据.以`quantile`标准化数据为例.<br>
2、The experiment and phenotype details, which are part of the `affyBatch` or `ExpressionSet` object.需要`affyBatch`或`ExpressionSet`对象的实验和表型细节部分.<br>
4、The R library, `limma`, that houses one of the most popular methods in R for differential gene expression analysis.R 库,limma是R中用于差异基因表达分析的最流行的方法之一.<br>
5、For demonstration, use normal colon cancer preprocessed `affy` data from the `antiProfilesData` package.为了演示,使用来自`antiProfilesData`软件包的正常结肠癌的预处理`affy`数据.<br>

The `limma` stands for linear models for microarray dat, which is used for analyzing gene expression microarray data, especially the use of linear models for analyzing gene expression data. It implements several methods of linear modeling for microarray data that can be used to identify DE genes. At first, it fits a linear model for each gene in the data given a set of arrays. Thereafter, it uses an empirical Bayes method to assess differential expression. This computes the statistical test and corresponding score in the form of p-values, log fold change, and so on.<br>
“limma”代表微阵列数据的线性模型，用于分析基因表达微阵列数据，尤其是使用线性模型分析基因表达数据。它实现了可用于识别DE基因的微阵列数据的多种线性建模方法。首先，它对给定一组数组的数据中的每个基因进行线性模型的拟合。此后，它使用经验贝叶斯方法来评估差异表达。这以p值，log倍数变化等形式计算统计测试和相应分数。

1、Install and load the `limma` library together with the `affy` package and `antiProfilesData`.安装并载入limma, affy和antiProfilesData.


```R
source("http://bioconductor.org/biocLite.R")
```

    Bioconductor version 3.6 (BiocInstaller 1.28.0), ?biocLite for help
    


```R
biocLite(c("limma", "antiProfilesData"))
```

    BioC_mirror: http://bioconductor.org
    Using Bioconductor 3.6 (BiocInstaller 1.28.0), R 3.4.3 (2017-11-30).
    Installing package(s) 'limma', 'antiProfilesData'
    

    package 'limma' successfully unpacked and MD5 sums checked
    
    The downloaded binary packages are in
    	C:\Users\Master1\AppData\Local\Temp\Rtmp4CNv4Q\downloaded_packages
    

    installing the source package 'antiProfilesData'
    
    Old packages: 'GenomicRanges', 'yaml'
    


```R
library(affy) # Package for affy data handling
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
library(antiProfilesData) # Package containing input data
```


```R
library(affyPLM) # Normalization package for eSet
```

    Loading required package: gcrma
    Loading required package: preprocessCore
    


```R
library(limma) # limma analysis package
```

    
    Attaching package: 'limma'
    
    The following object is masked from 'package:BiocGenerics':
    
        plotMA
    
    

2、Get the data (the colon cancer data) from the antiProfilesData package and subset the first 16 samples that represent the normal and tumor samples (eight each).<br>
从antiProfilesData软件包中获取数据（结肠癌数据），并将代表正常和肿瘤样本的前16个样本（各8个）的子集合。


```R
data(apColonData)
mydata <- apColonData[, sampleNames(apColonData)[1:16]]
mydata
```


    ExpressionSet (storageMode: lockedEnvironment)
    assayData: 5339 features, 16 samples 
      element names: exprs 
    protocolData: none
    phenoData
      sampleNames: GSM95473 GSM95474 ... GSM95488 (16 total)
      varLabels: filename DB_ID ... Status (7 total)
      varMetadata: labelDescription
    featureData: none
    experimentData: use 'experimentData(object)'
    Annotation: hgu133plus2 



```R
mydata_quantile <- normalize.ExpressionSet.quantiles(mydata)
mydata_quantile
```


    ExpressionSet (storageMode: lockedEnvironment)
    assayData: 5339 features, 16 samples 
      element names: exprs 
    protocolData: none
    phenoData
      sampleNames: GSM95473 GSM95474 ... GSM95488 (16 total)
      varLabels: filename DB_ID ... Status (7 total)
      varMetadata: labelDescription
    featureData: none
    experimentData: use 'experimentData(object)'
    Annotation: hgu133plus2 



```R
?apColonData
```

3、Prepare a design matrix based on the experiment details and `phenoData`.根据实验细节和`phenoData`准备一个设计矩阵.<br>
The design matrix describes the experiment condition in each of its column. In our case, we have only two conditions, and hence, the single column design matrix with appropriate experiment indication works as well.<br>
设计矩阵描述了每个列中的实验条件。在我们的例子中，我们只有两个条件，因此，具有适当实验指示的单列设计矩阵也适用.


```R
design <- model.matrix(~0 + pData(mydata)$Status)   # design <- model.matrix(~-1 + factor(pData(myData_quantile)$type))
design
```


<table>
<thead><tr><th></th><th scope=col>pData(mydata)$Status</th></tr></thead>
<tbody>
	<tr><th scope=row>1</th><td>0</td></tr>
	<tr><th scope=row>2</th><td>0</td></tr>
	<tr><th scope=row>3</th><td>0</td></tr>
	<tr><th scope=row>4</th><td>0</td></tr>
	<tr><th scope=row>5</th><td>0</td></tr>
	<tr><th scope=row>6</th><td>0</td></tr>
	<tr><th scope=row>7</th><td>0</td></tr>
	<tr><th scope=row>8</th><td>0</td></tr>
	<tr><th scope=row>9</th><td>1</td></tr>
	<tr><th scope=row>10</th><td>1</td></tr>
	<tr><th scope=row>11</th><td>1</td></tr>
	<tr><th scope=row>12</th><td>1</td></tr>
	<tr><th scope=row>13</th><td>1</td></tr>
	<tr><th scope=row>14</th><td>1</td></tr>
	<tr><th scope=row>15</th><td>1</td></tr>
	<tr><th scope=row>16</th><td>1</td></tr>
</tbody>
</table>



4、Fit a linear model using the expression data and design matrix.使用表达式数据和设计矩阵拟合线性模型.


```R
fit <- lmFit(mydata_quantile,design)
fit
```


    An object of class "MArrayLM"
    $coefficients
                 pData(mydata)$Status
    1555078_at             -0.3891325
    238493_at               0.3318299
    1562133_x_at           -0.3822317
    1559616_x_at           -0.9623383
    235687_at              -0.6461947
    5334 more rows ...
    
    $rank
    [1] 1
    
    $assign
    [1] 1
    
    $qr
    $qr
      pData(mydata)$Status
    1            -2.828427
    2             0.000000
    3             0.000000
    4             0.000000
    5             0.000000
    11 more rows ...
    
    $qraux
    [1] 1
    
    $pivot
    [1] 1
    
    $tol
    [1] 1e-07
    
    $rank
    [1] 1
    
    
    $df.residual
    [1] 15 15 15 15 15
    5334 more elements ...
    
    $sigma
      1555078_at    238493_at 1562133_x_at 1559616_x_at    235687_at 
       0.3893793    1.3899476    0.4058717    0.6304099    0.6146187 
    5334 more elements ...
    
    $cov.coefficients
                         pData(mydata)$Status
    pData(mydata)$Status                0.125
    
    $stdev.unscaled
                 pData(mydata)$Status
    1555078_at              0.3535534
    238493_at               0.3535534
    1562133_x_at            0.3535534
    1559616_x_at            0.3535534
    235687_at               0.3535534
    5334 more rows ...
    
    $pivot
    [1] 1
    
    $Amean
      1555078_at    238493_at 1562133_x_at 1559616_x_at    235687_at 
      -0.2382232    0.6460257   -0.4434705   -0.8408450   -0.6992017 
    5334 more elements ...
    
    $method
    [1] "ls"
    
    $design
      pData(mydata)$Status
    1                    0
    2                    0
    3                    0
    4                    0
    5                    0
    11 more rows ...
    


5、After having a linear model fit, compute the moderated statistics for it by `eBayes` function.<br>
在进行线性模型拟合之后，通过eBayes函数计算它的调节统计量.


```R
fitE <- eBayes(fit)
fitE
```


    An object of class "MArrayLM"
    $coefficients
                 pData(mydata)$Status
    1555078_at             -0.3891325
    238493_at               0.3318299
    1562133_x_at           -0.3822317
    1559616_x_at           -0.9623383
    235687_at              -0.6461947
    5334 more rows ...
    
    $rank
    [1] 1
    
    $assign
    [1] 1
    
    $qr
    $qr
      pData(mydata)$Status
    1            -2.828427
    2             0.000000
    3             0.000000
    4             0.000000
    5             0.000000
    11 more rows ...
    
    $qraux
    [1] 1
    
    $pivot
    [1] 1
    
    $tol
    [1] 1e-07
    
    $rank
    [1] 1
    
    
    $df.residual
    [1] 15 15 15 15 15
    5334 more elements ...
    
    $sigma
      1555078_at    238493_at 1562133_x_at 1559616_x_at    235687_at 
       0.3893793    1.3899476    0.4058717    0.6304099    0.6146187 
    5334 more elements ...
    
    $cov.coefficients
                         pData(mydata)$Status
    pData(mydata)$Status                0.125
    
    $stdev.unscaled
                 pData(mydata)$Status
    1555078_at              0.3535534
    238493_at               0.3535534
    1562133_x_at            0.3535534
    1559616_x_at            0.3535534
    235687_at               0.3535534
    5334 more rows ...
    
    $pivot
    [1] 1
    
    $Amean
      1555078_at    238493_at 1562133_x_at 1559616_x_at    235687_at 
      -0.2382232    0.6460257   -0.4434705   -0.8408450   -0.6992017 
    5334 more elements ...
    
    $method
    [1] "ls"
    
    $design
      pData(mydata)$Status
    1                    0
    2                    0
    3                    0
    4                    0
    5                    0
    11 more rows ...
    
    $df.prior
    [1] 1.316078
    
    $s2.prior
    [1] 0.2580634
    
    $var.prior
    [1] 9.156446
    
    $proportion
    [1] 0.01
    
    $s2.post
      1555078_at    238493_at 1562133_x_at 1559616_x_at    235687_at 
       0.1602024    1.7969358    0.1722601    0.3861762    0.3681016 
    5334 more elements ...
    
    $t
                 pData(mydata)$Status
    1555078_at              -2.749843
    238493_at                0.700155
    1562133_x_at            -2.604830
    1559616_x_at            -4.380060
    235687_at               -3.012482
    5334 more rows ...
    
    $df.total
    [1] 16.31608 16.31608 16.31608 16.31608 16.31608
    5334 more elements ...
    
    $p.value
                 pData(mydata)$Status
    1555078_at           0.0140512222
    238493_at            0.4936939029
    1562133_x_at         0.0189345321
    1559616_x_at         0.0004465372
    235687_at            0.0081226177
    5334 more rows ...
    
    $lods
                 pData(mydata)$Status
    1555078_at             -3.5057863
    238493_at              -6.4960522
    1562133_x_at           -3.7865081
    1559616_x_at           -0.1540363
    235687_at              -2.9844557
    5334 more rows ...
    
    $F
    [1]  7.561638  0.490217  6.785137 19.184923  9.075050
    5334 more elements ...
    
    $F.p.value
    [1] 0.0140512222 0.4936939029 0.0189345321 0.0004465372 0.0081226177
    5334 more elements ...
    


6、The output can be ranked and top-ranking genes extracted from it.<br>
可以对输出内容进行排序并从中选择排在前面的基因，


```R
tested <- topTable(fitE, adjust="fdr", sort.by="B", number=Inf)
tested
```


<table>
<thead><tr><th></th><th scope=col>logFC</th><th scope=col>AveExpr</th><th scope=col>t</th><th scope=col>P.Value</th><th scope=col>adj.P.Val</th><th scope=col>B</th></tr></thead>
<tbody>
	<tr><th scope=row>210372_s_at</th><td> 5.200189   </td><td> 2.9369060  </td><td> 18.586820  </td><td>2.082126e-12</td><td>1.111647e-08</td><td>17.909406   </td></tr>
	<tr><th scope=row>205470_s_at</th><td> 6.255458   </td><td> 3.4624353  </td><td> 14.523453  </td><td>9.365043e-11</td><td>2.499998e-07</td><td>14.665745   </td></tr>
	<tr><th scope=row>204855_at</th><td>13.342429   </td><td> 8.3325984  </td><td> 12.511677  </td><td>8.820153e-10</td><td>1.260731e-06</td><td>12.634755   </td></tr>
	<tr><th scope=row>219795_at</th><td> 6.773997   </td><td> 3.2696607  </td><td> 12.454113  </td><td>9.445448e-10</td><td>1.260731e-06</td><td>12.571572   </td></tr>
	<tr><th scope=row>220133_at</th><td> 3.223162   </td><td> 1.5617431  </td><td> 12.127986  </td><td>1.399296e-09</td><td>1.494168e-06</td><td>12.207825   </td></tr>
	<tr><th scope=row>204259_at</th><td> 7.111360   </td><td> 3.6682683  </td><td> 10.809794  </td><td>7.511799e-09</td><td>6.684249e-06</td><td>10.631112   </td></tr>
	<tr><th scope=row>203585_at</th><td> 3.814912   </td><td> 2.3970818  </td><td> 10.575189  </td><td>1.030032e-08</td><td>7.014783e-06</td><td>10.331384   </td></tr>
	<tr><th scope=row>203485_at</th><td>-1.477512   </td><td>-0.7942994  </td><td>-10.526023  </td><td>1.101229e-08</td><td>7.014783e-06</td><td>10.267797   </td></tr>
	<tr><th scope=row>219752_at</th><td> 2.948165   </td><td> 1.8787889  </td><td> 10.473843  </td><td>1.182488e-08</td><td>7.014783e-06</td><td>10.200014   </td></tr>
	<tr><th scope=row>232151_at</th><td> 5.490196   </td><td> 3.4301028  </td><td> 10.321526  </td><td>1.457882e-08</td><td>7.783631e-06</td><td>10.000392   </td></tr>
	<tr><th scope=row>209309_at</th><td> 6.996555   </td><td> 4.1278207  </td><td>  9.848596  </td><td>2.835113e-08</td><td>1.376061e-05</td><td> 9.363467   </td></tr>
	<tr><th scope=row>230398_at</th><td> 5.665631   </td><td> 3.3003175  </td><td>  9.766528  </td><td>3.189540e-08</td><td>1.419079e-05</td><td> 9.250240   </td></tr>
	<tr><th scope=row>218412_s_at</th><td> 7.497466   </td><td> 5.2023793  </td><td>  9.423892  </td><td>5.256847e-08</td><td>2.158947e-05</td><td> 8.768637   </td></tr>
	<tr><th scope=row>203962_s_at</th><td>10.751804   </td><td> 7.5600927  </td><td>  9.370139  </td><td>5.692160e-08</td><td>2.166592e-05</td><td> 8.691761   </td></tr>
	<tr><th scope=row>212531_at</th><td>17.064917   </td><td>11.9881353  </td><td>  9.317129  </td><td>6.158685e-08</td><td>2.166592e-05</td><td> 8.615594   </td></tr>
	<tr><th scope=row>241782_at</th><td> 2.195624   </td><td> 1.1963186  </td><td>  9.253373  </td><td>6.773548e-08</td><td>2.166592e-05</td><td> 8.523515   </td></tr>
	<tr><th scope=row>235245_at</th><td> 2.334599   </td><td> 1.4756777  </td><td>  9.227164  </td><td>7.044712e-08</td><td>2.166592e-05</td><td> 8.485513   </td></tr>
	<tr><th scope=row>208209_s_at</th><td> 3.822005   </td><td> 2.4020938  </td><td>  9.203031  </td><td>7.304486e-08</td><td>2.166592e-05</td><td> 8.450445   </td></tr>
	<tr><th scope=row>1566766_a_at</th><td> 4.146974   </td><td> 2.6918098  </td><td>  9.136375  </td><td>8.075648e-08</td><td>2.269257e-05</td><td> 8.353196   </td></tr>
	<tr><th scope=row>244180_at</th><td>-1.344145   </td><td>-0.8307531  </td><td> -9.061854  </td><td>9.040027e-08</td><td>2.413235e-05</td><td> 8.243797   </td></tr>
	<tr><th scope=row>205566_at</th><td> 3.226118   </td><td> 1.9002566  </td><td>  9.012534  </td><td>9.744220e-08</td><td>2.477352e-05</td><td> 8.171000   </td></tr>
	<tr><th scope=row>206941_x_at</th><td>-1.005474   </td><td>-0.6635295  </td><td> -8.962673  </td><td>1.051493e-07</td><td>2.551782e-05</td><td> 8.097083   </td></tr>
	<tr><th scope=row>205349_at</th><td> 3.212065   </td><td> 1.6030105  </td><td>  8.791665  </td><td>1.368203e-07</td><td>3.176015e-05</td><td> 7.841098   </td></tr>
	<tr><th scope=row>202917_s_at</th><td> 9.110590   </td><td> 4.6962185  </td><td>  8.695275  </td><td>1.589515e-07</td><td>3.340168e-05</td><td> 7.695107   </td></tr>
	<tr><th scope=row>202357_s_at</th><td> 8.271597   </td><td> 5.6286043  </td><td>  8.689031  </td><td>1.605088e-07</td><td>3.340168e-05</td><td> 7.685609   </td></tr>
	<tr><th scope=row>205977_s_at</th><td> 4.119680   </td><td> 2.8835908  </td><td>  8.680509  </td><td>1.626604e-07</td><td>3.340168e-05</td><td> 7.672635   </td></tr>
	<tr><th scope=row>209409_at</th><td> 5.113576   </td><td> 3.5591380  </td><td>  8.444412  </td><td>2.360515e-07</td><td>4.667700e-05</td><td> 7.309324   </td></tr>
	<tr><th scope=row>218704_at</th><td> 7.976579   </td><td> 5.7665953  </td><td>  8.404005  </td><td>2.517568e-07</td><td>4.800462e-05</td><td> 7.246389   </td></tr>
	<tr><th scope=row>221815_at</th><td> 4.468376   </td><td> 3.0759494  </td><td>  8.351318  </td><td>2.738977e-07</td><td>5.042551e-05</td><td> 7.163994   </td></tr>
	<tr><th scope=row>205476_at</th><td>14.466160   </td><td>10.4372247  </td><td>  8.248990  </td><td>3.229343e-07</td><td>5.747153e-05</td><td> 7.002877   </td></tr>
	<tr><th scope=row>...</th><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td></tr>
	<tr><th scope=row>1562211_a_at</th><td>-3.627688e-03</td><td> 0.033859078 </td><td>-0.0281154771</td><td>0.9779111    </td><td>0.9832518    </td><td>-6.748436    </td></tr>
	<tr><th scope=row>220373_at</th><td> 5.740751e-03</td><td>-0.258401883 </td><td> 0.0269088363</td><td>0.9788588    </td><td>0.9840194    </td><td>-6.748470    </td></tr>
	<tr><th scope=row>231626_at</th><td> 2.757279e-02</td><td> 1.857394629 </td><td> 0.0253280978</td><td>0.9801005    </td><td>0.9850070    </td><td>-6.748514    </td></tr>
	<tr><th scope=row>220343_at</th><td>-2.956402e-03</td><td> 0.004189963 </td><td>-0.0251884424</td><td>0.9802102    </td><td>0.9850070    </td><td>-6.748517    </td></tr>
	<tr><th scope=row>1559429_a_at</th><td> 3.890462e-03</td><td>-0.040891499 </td><td> 0.0217618683</td><td>0.9829018    </td><td>0.9875259    </td><td>-6.748601    </td></tr>
	<tr><th scope=row>221690_s_at</th><td> 1.831913e-02</td><td> 0.659423984 </td><td> 0.0214927376</td><td>0.9831133    </td><td>0.9875525    </td><td>-6.748608    </td></tr>
	<tr><th scope=row>1553059_at</th><td>-3.434531e-03</td><td>-0.163317401 </td><td>-0.0210531485</td><td>0.9834586    </td><td>0.9877136    </td><td>-6.748617    </td></tr>
	<tr><th scope=row>206834_at</th><td>-4.011657e-03</td><td>-0.167439369 </td><td>-0.0198192016</td><td>0.9844279    </td><td>0.9885012    </td><td>-6.748644    </td></tr>
	<tr><th scope=row>234361_at</th><td> 3.016104e-03</td><td> 0.044745530 </td><td> 0.0192995924</td><td>0.9848362    </td><td>0.9887251    </td><td>-6.748654    </td></tr>
	<tr><th scope=row>224390_s_at</th><td>-1.756399e-03</td><td>-0.001132483 </td><td>-0.0188067738</td><td>0.9852233    </td><td>0.9889247    </td><td>-6.748664    </td></tr>
	<tr><th scope=row>240967_at</th><td> 2.842829e-03</td><td>-0.164153264 </td><td> 0.0185750184</td><td>0.9854054    </td><td>0.9889247    </td><td>-6.748669    </td></tr>
	<tr><th scope=row>243231_at</th><td>-4.517521e-03</td><td> 0.252822841 </td><td>-0.0168659091</td><td>0.9867481    </td><td>0.9900861    </td><td>-6.748700    </td></tr>
	<tr><th scope=row>1552526_at</th><td>-1.509388e-03</td><td> 0.042989550 </td><td>-0.0147194442</td><td>0.9884345    </td><td>0.9915918    </td><td>-6.748736    </td></tr>
	<tr><th scope=row>213780_at</th><td>-1.830945e-03</td><td>-0.012559300 </td><td>-0.0120276478</td><td>0.9905494    </td><td>0.9935268    </td><td>-6.748774    </td></tr>
	<tr><th scope=row>229152_at</th><td> 7.274305e-03</td><td> 0.785249728 </td><td> 0.0094557174</td><td>0.9925702    </td><td>0.9953667    </td><td>-6.748803    </td></tr>
	<tr><th scope=row>242883_at</th><td> 1.284877e-03</td><td>-0.202721125 </td><td> 0.0086463667</td><td>0.9932061    </td><td>0.9958174    </td><td>-6.748810    </td></tr>
	<tr><th scope=row>221658_s_at</th><td>-1.067008e-03</td><td>-0.026685009 </td><td>-0.0081065141</td><td>0.9936303    </td><td>0.9960556    </td><td>-6.748815    </td></tr>
	<tr><th scope=row>231386_at</th><td> 8.382728e-04</td><td> 0.027066565 </td><td> 0.0077754172</td><td>0.9938905    </td><td>0.9961049    </td><td>-6.748818    </td></tr>
	<tr><th scope=row>1555643_s_at</th><td>-1.162331e-03</td><td> 0.030274093 </td><td>-0.0075690004</td><td>0.9940526    </td><td>0.9961049    </td><td>-6.748819    </td></tr>
	<tr><th scope=row>207972_at</th><td>-9.074408e-04</td><td>-0.120894936 </td><td>-0.0060321130</td><td>0.9952602    </td><td>0.9971208    </td><td>-6.748830    </td></tr>
	<tr><th scope=row>210402_at</th><td>-8.584665e-04</td><td> 0.164386348 </td><td>-0.0058033702</td><td>0.9954400    </td><td>0.9971208    </td><td>-6.748832    </td></tr>
	<tr><th scope=row>1552948_at</th><td>-9.712917e-04</td><td> 0.037384033 </td><td>-0.0055466975</td><td>0.9956417    </td><td>0.9971358    </td><td>-6.748833    </td></tr>
	<tr><th scope=row>205685_at</th><td>-8.604690e-04</td><td> 0.230873696 </td><td>-0.0047693364</td><td>0.9962525    </td><td>0.9975604    </td><td>-6.748837    </td></tr>
	<tr><th scope=row>231794_at</th><td>-6.487468e-04</td><td> 0.098780981 </td><td>-0.0034424825</td><td>0.9972950    </td><td>0.9984065    </td><td>-6.748843    </td></tr>
	<tr><th scope=row>242809_at</th><td>-5.042769e-04</td><td>-0.022999856 </td><td>-0.0029907888</td><td>0.9976500    </td><td>0.9984065    </td><td>-6.748845    </td></tr>
	<tr><th scope=row>203603_s_at</th><td> 2.719045e-03</td><td> 1.850278399 </td><td> 0.0027853887</td><td>0.9978114    </td><td>0.9984065    </td><td>-6.748845    </td></tr>
	<tr><th scope=row>208402_at</th><td>-3.644787e-04</td><td>-0.131042892 </td><td>-0.0026030083</td><td>0.9979547    </td><td>0.9984065    </td><td>-6.748846    </td></tr>
	<tr><th scope=row>1553534_at</th><td>-2.792301e-04</td><td>-0.044046503 </td><td>-0.0025039235</td><td>0.9980325    </td><td>0.9984065    </td><td>-6.748846    </td></tr>
	<tr><th scope=row>237320_at</th><td> 4.826854e-05</td><td>-0.213631604 </td><td> 0.0002714129</td><td>0.9997867    </td><td>0.9999727    </td><td>-6.748849    </td></tr>
	<tr><th scope=row>207067_s_at</th><td>-7.696790e-06</td><td> 0.334753187 </td><td>-0.0000347936</td><td>0.9999727    </td><td>0.9999727    </td><td>-6.748849    </td></tr>
</tbody>
</table>



7、To add the conditions of the p-values or other conditions to get the DE genes.添加p值或其他条件的条件以获得DE基因.


```R
DE <- tested[tested$adj.P.Val<0.01,]
```


```R
dim(DE)
```


<ol class=list-inline>
	<li>1211</li>
	<li>6</li>
</ol>




```R
DE <- tested[tested$adj.P.Val< 0.01 & abs(tested$logFC) >2, ]
```


```R
dim(DE)
```


<ol class=list-inline>
	<li>558</li>
	<li>6</li>
</ol>



There are several other packages for differential gene expression analysis. One of them is EMA, which can be used for many purposes. We will use the package for clustering and heatmap generation. To learn more about `EMA` package for differential gene analysis, check the help file.<br>
还有几个其他的差异基因表达分析软件包.其中之一是EMA,可用于多种用途.我们将使用该包进行聚类和热图生成.要了解有关使用EMA包进行差异基因分析的更多信息，请查看相应的帮助文件。
> library(EMA)<br}>
> help(package='EMA')
