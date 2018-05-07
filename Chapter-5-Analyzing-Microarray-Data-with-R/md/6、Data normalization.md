
# Data normalization
Microarrays can measure the expression levels of thousands of genes simultaneously. To make the data be comparable, normalization of data is essential. In this section, there are three methods(the `vsn`, `loess`, and `quantile` normalizations) of many normalization methods developed for data normalization.<br>
To start with this section, we need to define our data as an `ExpressionSet` or `AffyBatch` object. Here, take the breast cancer data object that we created in the first part as an example. It is an `AffyBatch` object.

1、Get the data as an `AffyBatch` object or a matrix.


```R
library(affy)
```


```R
library(hgu133a2cdf)
```


```R
mydata <- ReadAffy(celfile.path= "D:/Try-practice/Chapter 5/GSE24460_RAW")
mydata
```


    AffyBatch object
    size of arrays=732x732 features (18 kb)
    cdf=HG-U133A_2 (22277 affyids)
    number of samples=4
    number of genes=22277
    annotation=hgu133a2
    notes=


2、The `loess` normalization uses the `affy` library. To do the `loess` normalization,use the `normalize.AffyBatch.loess` function.


```R
mydata.loess <- normalize.AffyBatch.loess(mydata)
mydata.loess
boxplot(mydata.loess)
```

    Done with 1 vs 2 in iteration 1 
    Done with 1 vs 3 in iteration 1 
    Done with 1 vs 4 in iteration 1 
    Done with 2 vs 3 in iteration 1 
    Done with 2 vs 4 in iteration 1 
    Done with 3 vs 4 in iteration 1 
    1 0.5262025 
    


    AffyBatch object
    size of arrays=732x732 features (18 kb)
    cdf=HG-U133A_2 (22277 affyids)
    number of samples=4
    number of genes=22277
    annotation=hgu133a2
    notes=



[output-loess](https://github.com/Chengshu21/Chapter-5-Analyzing-Microarray-Data-with-R/blob/master/md/pic/output-loess.png)


3、The `quantile` normalization uses the `affy` library. To do the `quantile` normalization, use the `normalize.AffyBatch.quantile` function.


```R
mydata.quantile <- normalize.AffyBatch.quantiles(mydata)
mydata.quantile
boxplot(mydata.quantile)
```


    AffyBatch object
    size of arrays=732x732 features (18 kb)
    cdf=HG-U133A_2 (22277 affyids)
    number of samples=4
    number of genes=22277
    annotation=hgu133a2
    notes=



[output-quantile](https://github.com/Chengshu21/Chapter-5-Analyzing-Microarray-Data-with-R/blob/master/md/pic/output-quantile.png)


Normalization can also be performed using the `normalize.loess` and `normalize.quantile` functions with the data matrix as the input argument if the data is only a matrix of intensities.<br>
After that, do a second round of quality check for the normalized data and observe the effect of normalization in the same way as before using `arrayQualityMetrics`.
