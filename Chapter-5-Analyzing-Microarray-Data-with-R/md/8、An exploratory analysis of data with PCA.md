
# An exploratory analysis of data with PCA
If we measure too many human genes(>2000) in many samples, we will get a multiple column-row matrix. That will makes it difficult to build a pattern. So it is better to smplify the multidimensional pattern to a lower dimensions in order to explain and graphically represent the data.<br>
Principal Components Analysis (PCA) is a method that achieves this by performing a covariance analysis between factors. This finds the orthogonal components that represent the data and each component (called principal components) that represents the dimension where the features are more extended. Thus, PCA projects data onto a lower dimensional space and can be used as an exploratory method to serve purposes such as finding patterns in data and noise reduction. In this section, we will deal with PCA for microarray data.<br>
如果我们在很多样本中测量了太多的人类基因（>2000）,我们将得到一个多列-多行矩阵。这会让人难以构建出一种模式。因此，将多维模式缩小到更低的维度，能够很好的解释并图形化这些数据。<br>
主成分分析（PCA）是一种通过执行因素之间的协方差分析来实现这一点的方法。这将会找到表现数据的正交分量和每个组件（称为主要组件），这些组件都有着被极大扩展的维度特征。因此，PCA将数据投影到较低维的空间上，又可以被用作探索性方法来提供诸如查找数据模式和降低噪声等的用途。 在本节中，我们将处理PCA的微阵列数据。


```R
library(affy)
mydata <- ReadAffy(celfile.path= "D:/Try-practice/Chapter 5/GSE24460_RAW/")
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
    
    

1、select your data for the analysis.


```R
mydata.pca <- exprs(mydata)
```

2、Transpose the matrix to get the genes (features) as columns。


```R
mypca <- prcomp(mydata.pca, scale=TRUE)
```

3、check the principal components, take a look at the summary of created objects.


```R
summary(mypca)
```


    Importance of components:
                              PC1     PC2     PC3     PC4
    Standard deviation     1.9305 0.46011 0.17821 0.17182
    Proportion of Variance 0.9317 0.05293 0.00794 0.00738
    Cumulative Proportion  0.9317 0.98468 0.99262 1.00000


5、Create a vector of colors for every sample in the data。


```R
colors <- c("green","cyan","violet","magenta")
```

6、plot the principal components, use the `pairs` function.(note that it might take some time depending on the size of the data)


```R
pairs(mypca$x, col=colors)
```


[output](https://github.com/Chengshu21/Chapter-5-Analyzing-Microarray-Data-with-R/blob/master/md/pic/output-An%20exploratory%20analysis%20of%20data%20with%20PCA.png)


The PCA computation via the `prcomp` function performs the principal component analysis on the data matrix. It returns the principal components, their standard deviations (the square roots of Eigen values), and the rotation (containing the Eigen vectors).
