
# Handling time series data

In fact, a cell sample is given a certain treatment, and its expression can change along the course of time. Consider that during stem cell or embryonic development, the expression of genes at different time points will vary. And handling such time course expression data, though not much different from the standard protocol described earlier, needs small modifications here.<br>
事实上，细胞样本被给予一定的处理，其表达会随着时间的推移而变化。考虑到在干细胞或胚胎发育期间，不同时间点的基因表达会有所不同。处理这样的时间过程表达数据虽然与前面描述的标准协议没有多大区别，但在这里需要很小的修改。

Required:<br>
1、We need time course data. Take the yeast (Saccharomyces cerevisae) dataset from the `Mfuzz` package as an example.们需要时间相关数据。以“Mfuzz”包中的酵母（酿酒酵母）数据集为例。<br> 
2、The installation of the `Mfuzz` package. 安装`Mfuzz`包。<br>

1、Install and load the `Mfuzz` library.安装并载入`Mfuzz`库。


```R
library(Mfuzz)
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
    
    Loading required package: e1071
    

2、As the dataset is `ExpressionSet`, which needs some other packages for direct handling (the `affy` package needs the `AffyBatch` object), use the `affyPLM` library.


```R
library(affy)
```


```R
library(antiProfilesData)
```


```R
library(affyPLM)
```

    Loading required package: gcrma
    Loading required package: preprocessCore
    

3、Load the data using the following `data` function.


```R
data(yeast)
```

4、Check the quality of the data using the density plots, boxplots, and so on.


```R
class(yeast)
```


'ExpressionSet'



```R
plotDensity(yeast)
```


[plotDensity](https://github.com/Chengshu21/Chapter-5-Analyzing-Microarray-Data-with-R/blob/master/md/pic/output-Handling%20time%20series%20data.png)



```R
boxplot(yeast)
```


[boxplot](https://github.com/Chengshu21/Chapter-5-Analyzing-Microarray-Data-with-R/blob/master/md/pic/output-2Handling%20time%20series%20data.png)


5、To normalize the data, use `normalize.ExpressionSet.quantiles` from the `affyPLM` package(here, use `quantile` normalization).对数据进行标准化，使用`affyPLM`中的`normalize.ExpressionSet.quantiles`函数。


```R
library(affy)
```


```R
yeast_norm <- normalize.ExpressionSet.quantiles(yeast)
```

6、Perform the quality assessment for the normalized data again.对标准化数据做质量评价。


```R
plotDensity(yeast_norm)
```


[plotDensity](https://github.com/Chengshu21/Chapter-5-Analyzing-Microarray-Data-with-R/blob/master/md/pic/output-3Handling%20time%20series%20data.png)



```R
boxplot(yeast_norm)
```


[boxplot](https://github.com/Chengshu21/Chapter-5-Analyzing-Microarray-Data-with-R/blob/master/md/pic/output-4Handling%20time%20series%20data.png)


7、To check the attributes of the data, get the component details of the `ExpressionSet` object.


```R
pData(yeast_norm) 
#有16个样本（从0点往后每隔10个单位直到160）
```


<table>
<thead><tr><th></th><th scope=col>index</th><th scope=col>label</th><th scope=col>time</th></tr></thead>
<tbody>
	<tr><th scope=row>cdc28_0</th><td> 0       </td><td>cdc28_0  </td><td>  0      </td></tr>
	<tr><th scope=row>cdc28_10</th><td> 1       </td><td>cdc28_10 </td><td> 10      </td></tr>
	<tr><th scope=row>cdc28_20</th><td> 2       </td><td>cdc28_20 </td><td> 20      </td></tr>
	<tr><th scope=row>cdc28_30</th><td> 3       </td><td>cdc28_30 </td><td> 30      </td></tr>
	<tr><th scope=row>cdc28_40</th><td> 4       </td><td>cdc28_40 </td><td> 40      </td></tr>
	<tr><th scope=row>cdc28_50</th><td> 5       </td><td>cdc28_50 </td><td> 50      </td></tr>
	<tr><th scope=row>cdc28_60</th><td> 6       </td><td>cdc28_60 </td><td> 60      </td></tr>
	<tr><th scope=row>cdc28_70</th><td> 7       </td><td>cdc28_70 </td><td> 70      </td></tr>
	<tr><th scope=row>cdc28_80</th><td> 8       </td><td>cdc28_80 </td><td> 80      </td></tr>
	<tr><th scope=row>cdc28_90</th><td> 9       </td><td>cdc28_90 </td><td> 90      </td></tr>
	<tr><th scope=row>cdc28_100</th><td>10       </td><td>cdc28_100</td><td>100      </td></tr>
	<tr><th scope=row>cdc28_110</th><td>11       </td><td>cdc28_110</td><td>110      </td></tr>
	<tr><th scope=row>cdc28_120</th><td>12       </td><td>cdc28_120</td><td>120      </td></tr>
	<tr><th scope=row>cdc28_130</th><td>13       </td><td>cdc28_130</td><td>130      </td></tr>
	<tr><th scope=row>cdc28_140</th><td>14       </td><td>cdc28_140</td><td>140      </td></tr>
	<tr><th scope=row>cdc28_150</th><td>15       </td><td>cdc28_150</td><td>150      </td></tr>
	<tr><th scope=row>cdc28_160</th><td>16       </td><td>cdc28_160</td><td>160      </td></tr>
</tbody>
</table>



8、The design matrix can accordingly be created for a time series where there will be controls and time points. For example, if you have two replicates, C, for each control and time points, T1 and T2,


```R
times <- pData(yeast_norm)$time
times <- as.factor(times)
times
```


<ol class=list-inline>
	<li>0</li>
	<li>10</li>
	<li>20</li>
	<li>30</li>
	<li>40</li>
	<li>50</li>
	<li>60</li>
	<li>70</li>
	<li>80</li>
	<li>90</li>
	<li>100</li>
	<li>110</li>
	<li>120</li>
	<li>130</li>
	<li>140</li>
	<li>150</li>
	<li>160</li>
</ol>




```R
design <- model.matrix(~0 +factor(pData(yeast_norm)$time))
```


```R
colnames(design)[1:17] <- c("C", paste("T", 1:16, sep=""))
```

9、Create a contrast matrix using the `0th` point as a reference and all the other time points as the treatment (usually done for samples from the same culture).


```R
cont <- makeContrasts(C-T1, C-T2, C-T3, C-T4, C-T5, C-T6, C-T7, C-T8, C-T9, C-T10,
                      C-T11, C-T12, C-T13, C-T14, C-T15, C-T16, levels=design)
```

10、Use this matrix to fit the linear model, followed by the `eBayes` function to compute statistics.


```R
library(limma)
```


```R
fit <- lmFit(yeast_norm, cont)
fitE <- eBayes(fit)
```

    Warning message:
    "Partial NA coefficients for 332 probe(s)"

11、Filter the top-ranking genes using the `topTable` function.


```R
x <- topTable(fitE, adjust="fdr", sort.by="F", number=100)
y <- x[x$adj.P.Val< 0.05,]
head(y)
```


<table>
<thead><tr><th></th><th scope=col>C...T1</th><th scope=col>C...T2</th><th scope=col>C...T3</th><th scope=col>C...T4</th><th scope=col>C...T5</th><th scope=col>C...T6</th><th scope=col>C...T7</th><th scope=col>C...T8</th><th scope=col>C...T9</th><th scope=col>C...T10</th><th scope=col>C...T11</th><th scope=col>C...T12</th><th scope=col>C...T13</th><th scope=col>C...T14</th><th scope=col>C...T15</th><th scope=col>C...T16</th><th scope=col>AveExpr</th><th scope=col>F</th><th scope=col>P.Value</th><th scope=col>adj.P.Val</th></tr></thead>
<tbody>
	<tr><th scope=row>YDR535C</th><td> 0.045373350 </td><td>-0.3547546   </td><td> 3.0812557   </td><td> 3.0812557   </td><td>-0.9260739   </td><td> 0.4494472   </td><td>-1.056314    </td><td>-0.6702775   </td><td>-1.5162388   </td><td>-0.41473668  </td><td>-0.9354636   </td><td>-0.529467891 </td><td>-1.0627284   </td><td> 0.5007095   </td><td> 0.96456085  </td><td>-0.8230849   </td><td>-0.0005090025</td><td>256.5624     </td><td>2.244303e-05 </td><td>0.02482905   </td></tr>
	<tr><th scope=row>YMR317W</th><td>-1.361618986 </td><td>-0.6120542   </td><td> 0.1716029   </td><td>-1.2922435   </td><td>-1.2019479   </td><td>-2.2924229   </td><td>-1.580176    </td><td>-0.8586329   </td><td>-0.4953870   </td><td> 0.66528863  </td><td> 0.2092902   </td><td> 1.339882439 </td><td> 1.4571198   </td><td> 1.5013949   </td><td> 2.47249069  </td><td> 2.1328350   </td><td>-0.0283052587</td><td>193.5176     </td><td>4.037131e-05 </td><td>0.02482905   </td></tr>
	<tr><th scope=row>YHR137W</th><td> 0.085269790 </td><td> 0.9183755   </td><td> 1.1131313   </td><td> 0.9793744   </td><td> 1.6141288   </td><td> 1.9857307   </td><td> 1.291442    </td><td>-0.2478004   </td><td>-0.3620259   </td><td>-0.72751145  </td><td>-0.8521049   </td><td>-0.705844221 </td><td>-1.0954868   </td><td>-1.5246807   </td><td>-1.27434537  </td><td>-1.1303175   </td><td> 0.0120720297</td><td>163.5566     </td><td>5.728278e-05 </td><td>0.02482905   </td></tr>
	<tr><th scope=row>YML047C</th><td>-0.459388638 </td><td>-0.7911001   </td><td> 0.2786879   </td><td> 1.8677038   </td><td> 0.2557558   </td><td>-0.6688607   </td><td> 1.407027    </td><td> 1.3683844   </td><td> 1.4283458   </td><td>-0.43175983  </td><td>-0.1842532   </td><td> 0.651554550 </td><td>-1.3090766   </td><td>-1.2034157   </td><td>-1.33061493  </td><td>-0.9217954   </td><td> 0.0045702974</td><td>153.9416     </td><td>6.497142e-05 </td><td>0.02482905   </td></tr>
	<tr><th scope=row>YMR032W</th><td> 1.801686846 </td><td> 0.9413805   </td><td> 1.3044926   </td><td> 1.0622749   </td><td>-0.4679252   </td><td>-1.9259800   </td><td>-2.229041    </td><td>-1.4127112   </td><td>-0.4110722   </td><td>-0.02723012  </td><td> 0.7586565   </td><td> 1.155280123 </td><td>-0.1675567   </td><td>-1.0783823   </td><td>-0.87879003  </td><td>-0.7665813   </td><td> 0.0350770571</td><td>146.6830     </td><td>7.183158e-05 </td><td>0.02482905   </td></tr>
	<tr><th scope=row>YGL184C</th><td>-0.001825009 </td><td>-0.2441705   </td><td>-1.6977556   </td><td>-1.9718563   </td><td>-1.4590626   </td><td> 0.4270147   </td><td> 1.178113    </td><td> 1.3716358   </td><td> 0.7443744   </td><td> 0.94560458  </td><td> 0.9217415   </td><td> 0.007869278 </td><td> 0.1914265   </td><td> 0.3503931   </td><td> 0.01746599  </td><td>-0.6565194   </td><td>-0.0103381969</td><td>130.2669     </td><td>9.191468e-05 </td><td>0.02482905   </td></tr>
</tbody>
</table>


