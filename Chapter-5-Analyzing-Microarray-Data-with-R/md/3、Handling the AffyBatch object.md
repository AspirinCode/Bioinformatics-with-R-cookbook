
# Handling the AffyBatch object
During analyzing the microarray data, the `AffyBatch` object will be always used. As we have seen, it can be created
by reading the CEL files for an experiment together with the other allocated information. In this part, the main is to
look at the various components of such an object. Use the `AffyBatch` object created in the `mydata` object from the part 'Reading CEL files' as an example.<br>
在分析微阵列数据的过程中，总是会使用AffyBatch对象。 正如我们所看到的，它可以通过与其他分配的信息一起读取“CEL”文件进行实验来创建。在这一部分中，主要是看这种对象的各个组成部分。 以'Reading CEL files'部分的`mydata`对象中创建的`AffyBatch`对象为例。


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
library(hgu133a2cdf)
```

1、Check an `AffyBatch` object by simply typing in the object name, and check the structure of the object by `str()`.


```R
mydata <- ReadAffy(filenames="D:/Try-practice/Chapter 5/GSE24460_RAW/GSM602658_MCF71.CEL")
mydata
```


    AffyBatch object
    size of arrays=732x732 features (17 kb)
    cdf=HG-U133A_2 (22277 affyids)
    number of samples=1
    number of genes=22277
    annotation=hgu133a2
    notes=



```R
str(mydata)
```

    Formal class 'AffyBatch' [package "affy"] with 10 slots
      ..@ cdfName          : chr "HG-U133A_2"
      ..@ nrow             : Named int 732
      .. ..- attr(*, "names")= chr "Rows"
      ..@ ncol             : Named int 732
      .. ..- attr(*, "names")= chr "Cols"
      ..@ assayData        :<environment: 0x000000001929e910> 
      ..@ phenoData        :Formal class 'AnnotatedDataFrame' [package "Biobase"] with 4 slots
      .. .. ..@ varMetadata      :'data.frame':	1 obs. of  1 variable:
      .. .. .. ..$ labelDescription: chr "arbitrary numbering"
      .. .. ..@ data             :'data.frame':	1 obs. of  1 variable:
      .. .. .. ..$ sample: int 1
      .. .. ..@ dimLabels        : chr [1:2] "sampleNames" "sampleColumns"
      .. .. ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slot
      .. .. .. .. ..@ .Data:List of 1
      .. .. .. .. .. ..$ : int [1:3] 1 1 0
      ..@ featureData      :Formal class 'AnnotatedDataFrame' [package "Biobase"] with 4 slots
      .. .. ..@ varMetadata      :'data.frame':	0 obs. of  1 variable:
      .. .. .. ..$ labelDescription: chr(0) 
      .. .. ..@ data             :'data.frame':	535824 obs. of  0 variables
      .. .. ..@ dimLabels        : chr [1:2] "featureNames" "featureColumns"
      .. .. ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slot
      .. .. .. .. ..@ .Data:List of 1
      .. .. .. .. .. ..$ : int [1:3] 1 1 0
      ..@ experimentData   :Formal class 'MIAME' [package "Biobase"] with 13 slots
      .. .. ..@ name             : chr ""
      .. .. ..@ lab              : chr ""
      .. .. ..@ contact          : chr ""
      .. .. ..@ title            : chr ""
      .. .. ..@ abstract         : chr ""
      .. .. ..@ url              : chr ""
      .. .. ..@ pubMedIds        : chr ""
      .. .. ..@ samples          : list()
      .. .. ..@ hybridizations   : list()
      .. .. ..@ normControls     : list()
      .. .. ..@ preprocessing    :List of 2
      .. .. .. ..$ filenames  : chr "D:/Try-practice/Chapter 5/GSE24460_RAW/GSM602658_MCF71.CEL"
      .. .. .. ..$ affyversion: chr NA
      .. .. ..@ other            :List of 1
      .. .. .. ..$ : chr ""
      .. .. ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slot
      .. .. .. .. ..@ .Data:List of 2
      .. .. .. .. .. ..$ : int [1:3] 1 0 0
      .. .. .. .. .. ..$ : int [1:3] 1 1 0
      ..@ annotation       : chr "hgu133a2"
      ..@ protocolData     :Formal class 'AnnotatedDataFrame' [package "Biobase"] with 4 slots
      .. .. ..@ varMetadata      :'data.frame':	1 obs. of  1 variable:
      .. .. .. ..$ labelDescription: chr NA
      .. .. ..@ data             :'data.frame':	1 obs. of  1 variable:
      .. .. .. ..$ ScanDate: chr "01/10/06 11:14:07"
      .. .. ..@ dimLabels        : chr [1:2] "sampleNames" "sampleColumns"
      .. .. ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slot
      .. .. .. .. ..@ .Data:List of 1
      .. .. .. .. .. ..$ : int [1:3] 1 1 0
      ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slot
      .. .. ..@ .Data:List of 4
      .. .. .. ..$ : int [1:3] 3 4 2
      .. .. .. ..$ : int [1:3] 2 38 0
      .. .. .. ..$ : int [1:3] 1 3 0
      .. .. .. ..$ : int [1:3] 1 2 0
    

3、Check the phenotype data of the object using the `pData` and `phenoData` functions.


```R
pData(mydata)
```


<table>
<thead><tr><th></th><th scope=col>sample</th></tr></thead>
<tbody>
	<tr><th scope=row>GSM602658_MCF71.CEL</th><td>1</td></tr>
</tbody>
</table>




```R
phenoData(mydata)
```


    An object of class 'AnnotatedDataFrame'
      sampleNames: GSM602658_MCF71.CEL
      varLabels: sample
      varMetadata: labelDescription


4、Use `exprs` function to get the expression data as a matrix(These values can be written to a separate file using the `write.csv` function.<br>
使用`exprs`函数将表达式数据作为矩阵(可以使用`write.csv`函数将这些值写入单独的文件中)


```R
exprs(mydata)
```


<table>
<thead><tr><th></th><th scope=col>GSM602658_MCF71.CEL</th></tr></thead>
<tbody>
	<tr><th scope=row>1</th><td>  88</td></tr>
	<tr><th scope=row>2</th><td>7068</td></tr>
	<tr><th scope=row>3</th><td>  99</td></tr>
	<tr><th scope=row>4</th><td>7452</td></tr>
	<tr><th scope=row>5</th><td>  87</td></tr>
	<tr><th scope=row>6</th><td>  80</td></tr>
	<tr><th scope=row>7</th><td>7794</td></tr>
	<tr><th scope=row>8</th><td> 101</td></tr>
	<tr><th scope=row>9</th><td>7493</td></tr>
	<tr><th scope=row>10</th><td>  95</td></tr>
	<tr><th scope=row>11</th><td>7230</td></tr>
	<tr><th scope=row>12</th><td>  81</td></tr>
	<tr><th scope=row>13</th><td>7530</td></tr>
	<tr><th scope=row>14</th><td>  88</td></tr>
	<tr><th scope=row>15</th><td>7577</td></tr>
	<tr><th scope=row>16</th><td>  85</td></tr>
	<tr><th scope=row>17</th><td>7699</td></tr>
	<tr><th scope=row>18</th><td>  88</td></tr>
	<tr><th scope=row>19</th><td>7759</td></tr>
	<tr><th scope=row>20</th><td>  81</td></tr>
	<tr><th scope=row>21</th><td>7911</td></tr>
	<tr><th scope=row>22</th><td>  75</td></tr>
	<tr><th scope=row>23</th><td>7770</td></tr>
	<tr><th scope=row>24</th><td>  76</td></tr>
	<tr><th scope=row>25</th><td>7615</td></tr>
	<tr><th scope=row>26</th><td>  90</td></tr>
	<tr><th scope=row>27</th><td>7669</td></tr>
	<tr><th scope=row>28</th><td>  85</td></tr>
	<tr><th scope=row>29</th><td>7724</td></tr>
	<tr><th scope=row>30</th><td> 101</td></tr>
	<tr><th scope=row>...</th><td>...</td></tr>
	<tr><th scope=row>535795</th><td>11251</td></tr>
	<tr><th scope=row>535796</th><td>  136</td></tr>
	<tr><th scope=row>535797</th><td>11221</td></tr>
	<tr><th scope=row>535798</th><td>  107</td></tr>
	<tr><th scope=row>535799</th><td> 9126</td></tr>
	<tr><th scope=row>535800</th><td>  108</td></tr>
	<tr><th scope=row>535801</th><td>11946</td></tr>
	<tr><th scope=row>535802</th><td>  110</td></tr>
	<tr><th scope=row>535803</th><td>11589</td></tr>
	<tr><th scope=row>535804</th><td>  116</td></tr>
	<tr><th scope=row>535805</th><td>11741</td></tr>
	<tr><th scope=row>535806</th><td>  121</td></tr>
	<tr><th scope=row>535807</th><td>11715</td></tr>
	<tr><th scope=row>535808</th><td>  112</td></tr>
	<tr><th scope=row>535809</th><td>11069</td></tr>
	<tr><th scope=row>535810</th><td>  117</td></tr>
	<tr><th scope=row>535811</th><td>11213</td></tr>
	<tr><th scope=row>535812</th><td>  108</td></tr>
	<tr><th scope=row>535813</th><td>10196</td></tr>
	<tr><th scope=row>535814</th><td>  117</td></tr>
	<tr><th scope=row>535815</th><td>10649</td></tr>
	<tr><th scope=row>535816</th><td>  121</td></tr>
	<tr><th scope=row>535817</th><td>10106</td></tr>
	<tr><th scope=row>535818</th><td>  120</td></tr>
	<tr><th scope=row>535819</th><td>10496</td></tr>
	<tr><th scope=row>535820</th><td>  105</td></tr>
	<tr><th scope=row>535821</th><td>11432</td></tr>
	<tr><th scope=row>535822</th><td>  134</td></tr>
	<tr><th scope=row>535823</th><td>10979</td></tr>
	<tr><th scope=row>535824</th><td>  107</td></tr>
</tbody>
</table>



5、Use the `annotation` function to get the annotation name of the object.<br>
（使用`annotation`函数来获取对象的注释名称）


```R
annotation(mydata)
```


'hgu133a2'


6、Use `probNames` and `sampleNames` functions to get the probe names or sample names in the data.<br>
（使用`probNames`和`sampleNames`函数来获取数据中的探测名称或样本名称）


```R
probeNames(mydata)
```

    IOPub data rate exceeded.
    The notebook server will temporarily stop sending output
    to the client in order to avoid crashing it.
    To change this limit, set the config variable
    `--NotebookApp.iopub_data_rate_limit`.
    


```R
write.csv(probeNames(mydata), file = "D:/Try-practice/Chapter 5/probenames.csv")  # 读出数据
```

    Your code contains a unicode char which cannot be displayed in your
    current locale and R will silently convert it to an escaped form when the
    R kernel executes this code. This can lead to subtle errors if you use
    such chars to do comparisons. For more information, please see
    https://github.com/IRkernel/repr/wiki/Problems-with-unicode-on-windows


```R
 sampleNames(mydata)
```


'GSM602658_MCF71.CEL'


The `AffyBatch` object has a complex structure that contains many components and subcomponents. At every step, we check individual components of the complete `AffyBatch` object. The components of the `AffyBatch` object, such as expression data, can be extracted and written into a separate file using the `write.csv` or similar functions.<br>
AffyBatch对象具有复杂的结构，其中包含许多组件和子组件。 在每一步中，我们检查`AffyBatch`对象的各个组件。可以使用`write.csv`或类似函数将`AffyBatch`对象的组件（如表达式数据）读出来并写入单独的文件中。
