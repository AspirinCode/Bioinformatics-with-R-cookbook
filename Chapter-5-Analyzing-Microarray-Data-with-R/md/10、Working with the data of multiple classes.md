
# Working with the data of multiple classes
We know how to make a analysis with two experimental groups, namely, treatment and control. But if we faces three or more, how can we compare them systematically against each other?<br>
我们能够将对照组和实验组数据进行分析，但是如果我们遇到三个甚至更多的组别的时候该怎么将他们进行对比呢？<br>

Take another dataset from the `leukemiasEset` package as an example. The data is from 60 bone marrow samples of patients with one of the four main types of leukemia (ALL, AML, CLL, and CML) and non-leukemia controls. However, for demonstration purposes, we just take three samples from each of these categories. <br>

We need to principally focus on creating the design matrix for the comparison, which is the key difference. And then present it on imaginary data with three conditions and pairwise comparisons.主要需要专注于创建用于比较的设计矩阵，这是关键区别。 用三个条件和成对比较在虚数数据上显示，

1、Install and load the required libraries `leukemiasEset`.安装并载入库。


```R
source("http://www.bioconductor.org/biocLite.R")
```

    Bioconductor version 3.6 (BiocInstaller 1.28.0), ?biocLite for help
    


```R
biocLite("leukemiasEset")
```

    BioC_mirror: https://bioconductor.org
    Using Bioconductor 3.6 (BiocInstaller 1.28.0), R 3.4.3 (2017-11-30).
    Installing package(s) 'leukemiasEset'
    installing the source package 'leukemiasEset'
    
    Old packages: 'GenomicRanges', 'yaml'
    


```R
library(leukemiasEset)
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
    
    

2、Load the dataset wanted to start with from the `leukemiasEset` library.


```R
data(leukemiasEset)
```

3、Define the three conditions with two replicates for each.


```R
pheno <- pData(leukemiasEset)
```

4、select three samples from each set (use the corresponding indexes here).


```R
mydata <- leukemiasEset[, sampleNames(leukemiasEset)[c(1:3, 13:15, 25:27, 49:51)]]
```

5、Create the design matrix based on the condition variables.


```R
design <- model.matrix(~0 + factor(pData(mydata)$LeukemiaType))
```

6、Rename the columns of the design matrix.


```R
colnames(design) <- unique(as.character(pData(mydata)$LeukemiaType))
```

7、See the created design matrix.


```R
design # First three columns being different types of leukemiaand the fourth one being the control non leukemia samples
```


<table>
<thead><tr><th></th><th scope=col>ALL</th><th scope=col>AML</th><th scope=col>CLL</th><th scope=col>NoL</th></tr></thead>
<tbody>
	<tr><th scope=row>1</th><td>1</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>2</th><td>1</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>3</th><td>1</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>4</th><td>0</td><td>1</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>5</th><td>0</td><td>1</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>6</th><td>0</td><td>1</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>7</th><td>0</td><td>0</td><td>1</td><td>0</td></tr>
	<tr><th scope=row>8</th><td>0</td><td>0</td><td>1</td><td>0</td></tr>
	<tr><th scope=row>9</th><td>0</td><td>0</td><td>1</td><td>0</td></tr>
	<tr><th scope=row>10</th><td>0</td><td>0</td><td>0</td><td>1</td></tr>
	<tr><th scope=row>11</th><td>0</td><td>0</td><td>0</td><td>1</td></tr>
	<tr><th scope=row>12</th><td>0</td><td>0</td><td>0</td><td>1</td></tr>
</tbody>
</table>



8、Proceed by fitting a linear model to the data.


```R
library(limma)
```

    
    Attaching package: 'limma'
    
    The following object is masked from 'package:BiocGenerics':
    
        plotMA
    
    


```R
fit <- lmFit(mydata, design)
fit
```


    An object of class "MArrayLM"
    $coefficients
                         ALL      AML      CLL      NoL
    ENSG00000000003 3.478096 3.435959 3.608932 3.860124
    ENSG00000000005 3.540522 3.332321 3.379161 3.278787
    ENSG00000000419 9.083140 9.124083 9.639916 9.551689
    ENSG00000000457 4.865089 6.217431 5.838287 5.503014
    ENSG00000000460 4.294320 3.979242 3.612462 5.218094
    20167 more rows ...
    
    $rank
    [1] 4
    
    $assign
    [1] 1 1 1 1
    
    $qr
    $qr
             ALL        AML       CLL       NoL
    1 -1.7320508  0.0000000  0.000000  0.000000
    2  0.5773503 -1.7320508  0.000000  0.000000
    3  0.5773503  0.0000000 -1.732051  0.000000
    4  0.0000000  0.5773503  0.000000 -1.732051
    5  0.0000000  0.5773503  0.000000  0.000000
    7 more rows ...
    
    $qraux
    [1] 1.57735 1.00000 1.00000 1.00000
    
    $pivot
    [1] 1 2 3 4
    
    $tol
    [1] 1e-07
    
    $rank
    [1] 4
    
    
    $df.residual
    [1] 8 8 8 8 8
    20167 more elements ...
    
    $sigma
    ENSG00000000003 ENSG00000000005 ENSG00000000419 ENSG00000000457 ENSG00000000460 
          0.2263650       0.1942866       0.6732108       0.4570501       0.6726086 
    20167 more elements ...
    
    $cov.coefficients
              ALL       AML       CLL       NoL
    ALL 0.3333333 0.0000000 0.0000000 0.0000000
    AML 0.0000000 0.3333333 0.0000000 0.0000000
    CLL 0.0000000 0.0000000 0.3333333 0.0000000
    NoL 0.0000000 0.0000000 0.0000000 0.3333333
    
    $stdev.unscaled
                          ALL       AML       CLL       NoL
    ENSG00000000003 0.5773503 0.5773503 0.5773503 0.5773503
    ENSG00000000005 0.5773503 0.5773503 0.5773503 0.5773503
    ENSG00000000419 0.5773503 0.5773503 0.5773503 0.5773503
    ENSG00000000457 0.5773503 0.5773503 0.5773503 0.5773503
    ENSG00000000460 0.5773503 0.5773503 0.5773503 0.5773503
    20167 more rows ...
    
    $pivot
    [1] 1 2 3 4
    
    $Amean
    ENSG00000000003 ENSG00000000005 ENSG00000000419 ENSG00000000457 ENSG00000000460 
           3.595778        3.382697        9.349707        5.605955        4.276029 
    20167 more elements ...
    
    $method
    [1] "ls"
    
    $design
      ALL AML CLL NoL
    1   1   0   0   0
    2   1   0   0   0
    3   1   0   0   0
    4   0   1   0   0
    5   0   1   0   0
    7 more rows ...
    


9、Create a contrast matrix for pairwise comparisons.


```R
contrast.matrix <- makeContrasts(NoL- ALL, NoL- AML, NoL- CLL, levels = design)
contrast.matrix
```


<table>
<thead><tr><th></th><th scope=col>NoL - ALL</th><th scope=col>NoL - AML</th><th scope=col>NoL - CLL</th></tr></thead>
<tbody>
	<tr><th scope=row>ALL</th><td>-1</td><td> 0</td><td> 0</td></tr>
	<tr><th scope=row>AML</th><td> 0</td><td>-1</td><td> 0</td></tr>
	<tr><th scope=row>CLL</th><td> 0</td><td> 0</td><td>-1</td></tr>
	<tr><th scope=row>NoL</th><td> 1</td><td> 1</td><td> 1</td></tr>
</tbody>
</table>



10、Fit the linear model using this contrast matrix.


```R
fit2 <- contrasts.fit(fit, contrast.matrix)
```

11、Perform the empirical Bayes analysis of the model.


```R
fit2 <- eBayes(fit2)
```

12、Extract the differentially expressed gene for each of the pairwise comparisons using the `coef` argument in the `topTable` function. For the first pairwise comparison, set `coef=1` to compare non-lukemia control with Acute Lymphoblastic Leukemia (ALL) leukemia.


```R
tested2 <- topTable(fit2,adjust="fdr",sort.by="B",number=Inf, coef=1)
tested2
```


<table>
<thead><tr><th></th><th scope=col>logFC</th><th scope=col>AveExpr</th><th scope=col>t</th><th scope=col>P.Value</th><th scope=col>adj.P.Val</th><th scope=col>B</th></tr></thead>
<tbody>
	<tr><th scope=row>ENSG00000152078</th><td> 4.510507   </td><td>4.856523    </td><td> 28.13988   </td><td>4.463747e-11</td><td>9.004270e-07</td><td>14.014724   </td></tr>
	<tr><th scope=row>ENSG00000117519</th><td>-4.185175   </td><td>4.791585    </td><td>-22.73888   </td><td>3.878292e-10</td><td>3.911645e-06</td><td>12.697381   </td></tr>
	<tr><th scope=row>ENSG00000145850</th><td> 4.142236   </td><td>4.507655    </td><td> 17.38636   </td><td>5.759942e-09</td><td>2.925048e-05</td><td>10.727823   </td></tr>
	<tr><th scope=row>ENSG00000170180</th><td> 5.681327   </td><td>5.734169    </td><td> 17.37423   </td><td>5.800214e-09</td><td>2.925048e-05</td><td>10.722311   </td></tr>
	<tr><th scope=row>ENSG00000087586</th><td> 3.952183   </td><td>5.720789    </td><td> 16.45393   </td><td>9.977396e-09</td><td>3.111188e-05</td><td>10.287051   </td></tr>
	<tr><th scope=row>ENSG00000047597</th><td> 5.362419   </td><td>5.108415    </td><td> 16.32474   </td><td>1.079114e-08</td><td>3.111188e-05</td><td>10.223153   </td></tr>
	<tr><th scope=row>ENSG00000175449</th><td> 3.954293   </td><td>4.667288    </td><td> 16.32395   </td><td>1.079631e-08</td><td>3.111188e-05</td><td>10.222762   </td></tr>
	<tr><th scope=row>ENSG00000104043</th><td> 4.594551   </td><td>5.877417    </td><td> 15.91546   </td><td>1.388748e-08</td><td>3.501728e-05</td><td>10.015920   </td></tr>
	<tr><th scope=row>ENSG00000024526</th><td> 3.048430   </td><td>4.657230    </td><td> 15.70978   </td><td>1.580050e-08</td><td>3.541419e-05</td><td> 9.908949   </td></tr>
	<tr><th scope=row>ENSG00000115641</th><td> 3.660332   </td><td>5.460531    </td><td> 15.49887   </td><td>1.806560e-08</td><td>3.644192e-05</td><td> 9.797237   </td></tr>
	<tr><th scope=row>ENSG00000188672</th><td> 4.614139   </td><td>5.380265    </td><td> 15.26785   </td><td>2.096196e-08</td><td>3.844043e-05</td><td> 9.672452   </td></tr>
	<tr><th scope=row>ENSG00000133063</th><td> 4.857495   </td><td>5.630711    </td><td> 14.76852   </td><td>2.911671e-08</td><td>4.894519e-05</td><td> 9.393827   </td></tr>
	<tr><th scope=row>ENSG00000164330</th><td>-4.548742   </td><td>5.139590    </td><td>-14.54086   </td><td>3.393791e-08</td><td>4.941928e-05</td><td> 9.262590   </td></tr>
	<tr><th scope=row>ENSG00000173372</th><td> 2.784815   </td><td>5.523598    </td><td> 14.52528   </td><td>3.429853e-08</td><td>4.941928e-05</td><td> 9.253507   </td></tr>
	<tr><th scope=row>ENSG00000234955</th><td>-3.181503   </td><td>4.290691    </td><td>-14.21654   </td><td>4.237846e-08</td><td>5.699055e-05</td><td> 9.070902   </td></tr>
	<tr><th scope=row>ENSG00000102935</th><td>-5.341545   </td><td>4.829078    </td><td>-13.98406   </td><td>4.983512e-08</td><td>6.282962e-05</td><td> 8.929961   </td></tr>
	<tr><th scope=row>ENSG00000159189</th><td> 3.512535   </td><td>4.917536    </td><td> 13.86537   </td><td>5.418509e-08</td><td>6.429540e-05</td><td> 8.856847   </td></tr>
	<tr><th scope=row>ENSG00000107562</th><td> 3.359836   </td><td>4.984015    </td><td> 13.50747   </td><td>7.001837e-08</td><td>7.434224e-05</td><td> 8.631454   </td></tr>
	<tr><th scope=row>ENSG00000086506</th><td> 4.011569   </td><td>5.859682    </td><td> 13.50738   </td><td>7.002293e-08</td><td>7.434224e-05</td><td> 8.631396   </td></tr>
	<tr><th scope=row>ENSG00000130635</th><td>-4.395195   </td><td>5.362425    </td><td>-13.27887   </td><td>8.274167e-08</td><td>8.345325e-05</td><td> 8.483522   </td></tr>
	<tr><th scope=row>ENSG00000066923</th><td>-4.520295   </td><td>6.494213    </td><td>-12.93093   </td><td>1.072199e-07</td><td>1.029924e-04</td><td> 8.252195   </td></tr>
	<tr><th scope=row>ENSG00000162692</th><td> 4.844911   </td><td>4.630924    </td><td> 12.65173   </td><td>1.326070e-07</td><td>1.215885e-04</td><td> 8.060998   </td></tr>
	<tr><th scope=row>ENSG00000198336</th><td> 3.883021   </td><td>6.768568    </td><td> 12.52486   </td><td>1.462526e-07</td><td>1.232467e-04</td><td> 7.972428   </td></tr>
	<tr><th scope=row>ENSG00000109501</th><td>-4.056787   </td><td>5.692457    </td><td>-12.51679   </td><td>1.471700e-07</td><td>1.232467e-04</td><td> 7.966765   </td></tr>
	<tr><th scope=row>ENSG00000168497</th><td> 3.108121   </td><td>4.741989    </td><td> 12.40452   </td><td>1.606225e-07</td><td>1.232467e-04</td><td> 7.887420   </td></tr>
	<tr><th scope=row>ENSG00000080819</th><td> 4.365147   </td><td>5.626793    </td><td> 12.34964   </td><td>1.676802e-07</td><td>1.232467e-04</td><td> 7.848332   </td></tr>
	<tr><th scope=row>ENSG00000130821</th><td> 3.111452   </td><td>6.040674    </td><td> 12.31654   </td><td>1.721015e-07</td><td>1.232467e-04</td><td> 7.824650   </td></tr>
	<tr><th scope=row>ENSG00000162599</th><td> 2.402039   </td><td>4.391735    </td><td> 12.29547   </td><td>1.749813e-07</td><td>1.232467e-04</td><td> 7.809539   </td></tr>
	<tr><th scope=row>ENSG00000173391</th><td> 4.511073   </td><td>4.289218    </td><td> 12.27961   </td><td>1.771840e-07</td><td>1.232467e-04</td><td> 7.798144   </td></tr>
	<tr><th scope=row>ENSG00000071967</th><td> 3.474916   </td><td>5.842135    </td><td> 11.93721   </td><td>2.329409e-07</td><td>1.524925e-04</td><td> 7.547824   </td></tr>
	<tr><th scope=row>...</th><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td></tr>
	<tr><th scope=row>ENSG00000184898</th><td>-6.679917e-04</td><td> 4.385725    </td><td>-2.006986e-03</td><td>0.9984370    </td><td>0.9997165    </td><td>-6.787415    </td></tr>
	<tr><th scope=row>ENSG00000162526</th><td> 5.333565e-04</td><td> 6.845486    </td><td> 2.006777e-03</td><td>0.9984372    </td><td>0.9997165    </td><td>-6.787415    </td></tr>
	<tr><th scope=row>ENSG00000117174</th><td> 1.246917e-03</td><td> 7.133992    </td><td> 1.978443e-03</td><td>0.9984592    </td><td>0.9997165    </td><td>-6.787415    </td></tr>
	<tr><th scope=row>ENSG00000205038</th><td> 4.642149e-04</td><td> 3.707259    </td><td> 1.977952e-03</td><td>0.9984596    </td><td>0.9997165    </td><td>-6.787415    </td></tr>
	<tr><th scope=row>ENSG00000245927</th><td> 3.255004e-04</td><td> 3.583726    </td><td> 1.954928e-03</td><td>0.9984776    </td><td>0.9997165    </td><td>-6.787415    </td></tr>
	<tr><th scope=row>ENSG00000113327</th><td>-2.485544e-04</td><td> 3.238793    </td><td>-1.678736e-03</td><td>0.9986926    </td><td>0.9997335    </td><td>-6.787416    </td></tr>
	<tr><th scope=row>ENSG00000163479</th><td> 5.506234e-04</td><td>10.494864    </td><td> 1.663374e-03</td><td>0.9987046    </td><td>0.9997335    </td><td>-6.787416    </td></tr>
	<tr><th scope=row>ENSG00000011260</th><td>-5.297753e-04</td><td> 9.563794    </td><td>-1.649236e-03</td><td>0.9987156    </td><td>0.9997335    </td><td>-6.787416    </td></tr>
	<tr><th scope=row>ENSG00000121743</th><td> 2.468503e-04</td><td> 3.792622    </td><td> 1.588047e-03</td><td>0.9987633    </td><td>0.9997335    </td><td>-6.787416    </td></tr>
	<tr><th scope=row>ENSG00000148335</th><td>-3.661422e-04</td><td> 6.499836    </td><td>-1.554274e-03</td><td>0.9987896    </td><td>0.9997335    </td><td>-6.787416    </td></tr>
	<tr><th scope=row>ENSG00000092345</th><td>-2.874116e-04</td><td> 2.950834    </td><td>-1.551385e-03</td><td>0.9987918    </td><td>0.9997335    </td><td>-6.787416    </td></tr>
	<tr><th scope=row>ENSG00000166348</th><td> 1.128627e-03</td><td> 5.727922    </td><td> 1.482353e-03</td><td>0.9988456    </td><td>0.9997377    </td><td>-6.787416    </td></tr>
	<tr><th scope=row>ENSG00000018236</th><td> 3.532907e-04</td><td> 3.782055    </td><td> 1.403421e-03</td><td>0.9989071    </td><td>0.9997426    </td><td>-6.787416    </td></tr>
	<tr><th scope=row>ENSG00000164604</th><td>-2.339540e-04</td><td> 3.328498    </td><td>-1.327779e-03</td><td>0.9989660    </td><td>0.9997426    </td><td>-6.787416    </td></tr>
	<tr><th scope=row>ENSG00000118596</th><td> 5.541809e-04</td><td> 4.865807    </td><td> 1.285061e-03</td><td>0.9989992    </td><td>0.9997426    </td><td>-6.787416    </td></tr>
	<tr><th scope=row>ENSG00000230743</th><td> 3.698132e-04</td><td> 9.729000    </td><td> 1.052519e-03</td><td>0.9991803    </td><td>0.9997845    </td><td>-6.787417    </td></tr>
	<tr><th scope=row>ENSG00000101350</th><td> 5.139717e-04</td><td> 6.506761    </td><td> 1.017790e-03</td><td>0.9992074    </td><td>0.9997845    </td><td>-6.787417    </td></tr>
	<tr><th scope=row>ENSG00000168591</th><td>-3.643930e-04</td><td> 7.228878    </td><td>-9.431005e-04</td><td>0.9992655    </td><td>0.9997845    </td><td>-6.787417    </td></tr>
	<tr><th scope=row>ENSG00000101441</th><td>-2.494937e-04</td><td> 4.053601    </td><td>-8.461279e-04</td><td>0.9993411    </td><td>0.9997845    </td><td>-6.787417    </td></tr>
	<tr><th scope=row>ENSG00000197472</th><td> 1.479224e-04</td><td> 2.947163    </td><td> 8.299368e-04</td><td>0.9993537    </td><td>0.9997845    </td><td>-6.787417    </td></tr>
	<tr><th scope=row>ENSG00000139656</th><td>-1.622774e-04</td><td> 3.339458    </td><td>-8.232701e-04</td><td>0.9993589    </td><td>0.9997845    </td><td>-6.787417    </td></tr>
	<tr><th scope=row>ENSG00000147439</th><td>-3.089195e-04</td><td> 8.237133    </td><td>-7.858004e-04</td><td>0.9993880    </td><td>0.9997845    </td><td>-6.787417    </td></tr>
	<tr><th scope=row>ENSG00000002745</th><td> 1.185971e-04</td><td> 3.225855    </td><td> 5.995940e-04</td><td>0.9995331    </td><td>0.9998287    </td><td>-6.787417    </td></tr>
	<tr><th scope=row>ENSG00000134138</th><td>-2.667634e-04</td><td> 3.451379    </td><td>-5.995502e-04</td><td>0.9995331    </td><td>0.9998287    </td><td>-6.787417    </td></tr>
	<tr><th scope=row>ENSG00000122884</th><td> 2.610773e-04</td><td> 7.224659    </td><td> 5.358619e-04</td><td>0.9995827    </td><td>0.9998287    </td><td>-6.787417    </td></tr>
	<tr><th scope=row>ENSG00000182187</th><td>-7.211313e-05</td><td> 5.232893    </td><td>-4.745563e-04</td><td>0.9996304    </td><td>0.9998287    </td><td>-6.787417    </td></tr>
	<tr><th scope=row>ENSG00000239961</th><td> 2.724331e-04</td><td> 6.020185    </td><td> 3.815179e-04</td><td>0.9997029    </td><td>0.9998516    </td><td>-6.787417    </td></tr>
	<tr><th scope=row>ENSG00000163449</th><td> 5.805249e-05</td><td> 3.751309    </td><td> 3.007044e-04</td><td>0.9997658    </td><td>0.9998650    </td><td>-6.787417    </td></tr>
	<tr><th scope=row>ENSG00000151718</th><td> 8.595550e-06</td><td> 3.492013    </td><td> 5.490059e-05</td><td>0.9999572    </td><td>0.9999758    </td><td>-6.787417    </td></tr>
	<tr><th scope=row>ENSG00000138768</th><td>-1.614100e-05</td><td> 7.786855    </td><td>-3.104010e-05</td><td>0.9999758    </td><td>0.9999758    </td><td>-6.787417    </td></tr>
</tbody>
</table>




```R
DE2 <- tested2[tested2$adj.P.Val < 0.01,]
dim(DE2)
```


<ol class=list-inline>
	<li>252</li>
	<li>6</li>
</ol>



The major difference compared to the two conditions(Finding the differentially expressed genes) is the use of the contrast matrix and the `coef` argument in the `topTable` function. The contrast matrix enables the pairwise comparison for the computation of p-values. The model fitted returns a set of p-values for each comparison. The desired set of DE genes for a specific comparison can then be extracted using the proper `coef` value for comparisons.<br>
与两种条件（寻找差异表达基因）相比，主要区别在于使用对照矩阵和`topTable`函数中的`coef`参数。对比矩阵可以用于计算p值的成对比较。拟合的模型会为每次比较返回一组p值。然后可以使用适当的“coef”值来提取用于特定比较的期望的一组DE基因进行比较。
