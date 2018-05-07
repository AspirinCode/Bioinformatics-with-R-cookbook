
# 10、Analyzing methylation data

DNA methylation data can help to detect an epigenetic marker for which a detailed mechanism of mitotic inheritance has been described. It is usually done via methylation-specific DNA sequencing or via microarray data epigenetic regulation. Recent advances in NGS and microarray technology make it possible to map genome-wide DNA methylation at a high resolution and in a large number of samples. One of the common goals in analyzing methylation data is to identify differentially methylated regions (DMRs) when comparing control and treatment conditions. This recipe aims to address this issue.<br>
Here use the `methyAnalysis` package and the built-in dataset as an example.<br>
DNA甲基化数据可以帮助检测已经描述了有丝分裂遗传的详细机制的表观遗传标记。它通常通过甲基化特异性DNA测序或通过微阵列数据表观遗传调节完成。 NGS和微阵列技术的最新进展使得以高分辨率和大量样本绘制全基因组DNA甲基化成为可能。分析甲基化数据的共同目标之一是在比较对照和治疗条件时鉴别差异甲基化区域（DMR）。这个配方旨在解决这个问题。
这里以methyAnalysis包和内置数据集为例。

1、Install and load the library and dataset.


```R
source("http://www.bioconductor.org/biocLite.R")
```

    Bioconductor version 3.6 (BiocInstaller 1.28.0), ?biocLite for help
    


```R
biocLite(c("methyAnalysis", "TxDb.Hsapiens.UCSC.hg19.knownGene"))
```

    BioC_mirror: https://bioconductor.org
    Using Bioconductor 3.6 (BiocInstaller 1.28.0), R 3.4.4 (2018-03-15).
    Installing package(s) 'methyAnalysis', 'TxDb.Hsapiens.UCSC.hg19.knownGene'
    

    package 'methyAnalysis' successfully unpacked and MD5 sums checked
    
    The downloaded binary packages are in
    	C:\Users\Master1\AppData\Local\Temp\RtmpUziMsd\downloaded_packages
    

    installing the source package 'TxDb.Hsapiens.UCSC.hg19.knownGene'
    
    Old packages: 'foreign', 'futile.options', 'htmlwidgets', 'httpuv', 'lambda.r',
      'party', 'R.oo', 'readxl', 'XML'
    


```R
library(methyAnalysis)
data(exampleMethyGenoSet)
```

2、Look at different components to get some information, such as samples and location for the data.
查看数据中不同组成。


```R
slotNames(exampleMethyGenoSet)  # 对象数据的根目录类别
dim(exprs(exampleMethyGenoSet)) # 对象数据二维行列
str(exampleMethyGenoSet)        # 对象数据结构
# head(locData(exampleMethyGenoSet))
# pData(exampleMethyGenoSet)
class(exampleMethyGenoSet)      # 对象数据类别
colnames(exampleMethyGenoSet)   # 对象数据的行名
exampleMethyGenoSet             # 对象数据内容
```

    Your code contains a unicode char which cannot be displayed in your
    current locale and R will silently convert it to an escaped form when the
    R kernel executes this code. This can lead to subtle errors if you use
    such chars to do comparisons. For more information, please see
    https://github.com/IRkernel/repr/wiki/Problems-with-unicode-on-windows


<ol class=list-inline>
	<li>'history'</li>
	<li>'annotation'</li>
	<li>'rowRanges'</li>
	<li>'colData'</li>
	<li>'assays'</li>
	<li>'NAMES'</li>
	<li>'elementMetadata'</li>
	<li>'metadata'</li>
</ol>




<ol class=list-inline>
	<li>4243</li>
	<li>8</li>
</ol>



    Formal class 'MethyGenoSet' [package "methyAnalysis"] with 8 slots
      ..@ history        :'data.frame':	0 obs. of  4 variables:
      .. ..$ submitted  :Class 'AsIs'  logi(0) 
      .. ..$ finished   :Class 'AsIs'  logi(0) 
      .. ..$ command    :Class 'AsIs'  logi(0) 
      .. ..$ lumiVersion:Class 'AsIs'  logi(0) 
      ..@ annotation     : chr "IlluminaHumanMethylation450k.db"
      ..@ rowRanges      :Formal class 'GRanges' [package "GenomicRanges"] with 6 slots
      .. .. ..@ seqnames       :Formal class 'Rle' [package "IRanges"] with 4 slots
      .. .. .. .. ..@ values         : Factor w/ 1 level "chr21": 1
      .. .. .. .. ..@ lengths        : int 4243
      .. .. .. .. ..@ elementMetadata: NULL
      .. .. .. .. ..@ metadata       : list()
      .. .. ..@ ranges         :Formal class 'IRanges' [package "IRanges"] with 6 slots
      .. .. .. .. ..@ start          : int [1:4243] 10882029 10883548 10884748 10884967 10884969 10885409 10888524 10891858 10895603 10895780 ...
      .. .. .. .. ..@ width          : int [1:4243] 1 1 1 1 1 1 1 1 1 1 ...
      .. .. .. .. ..@ NAMES          : chr [1:4243] "cg17035109" "cg06187584" "cg12459059" "cg25450479" ...
      .. .. .. .. ..@ elementType    : chr "integer"
      .. .. .. .. ..@ elementMetadata: NULL
      .. .. .. .. ..@ metadata       : list()
      .. .. ..@ strand         :Formal class 'Rle' [package "IRanges"] with 4 slots
      .. .. .. .. ..@ values         : Factor w/ 3 levels "+","-","*": 3
      .. .. .. .. ..@ lengths        : int 4243
      .. .. .. .. ..@ elementMetadata: NULL
      .. .. .. .. ..@ metadata       : list()
      .. .. ..@ elementMetadata:Formal class 'DataFrame' [package "S4Vectors"] with 6 slots
      .. .. .. .. ..@ rownames       : NULL
      .. .. .. .. ..@ nrows          : int 4243
      .. .. .. .. ..@ listData       :List of 1
      .. .. .. .. .. ..$ ID: Factor w/ 4243 levels "cg00002080","cg00004533",..: 2788 1120 2103 3925 3598 708 3075 1762 375 2572 ...
      .. .. .. .. ..@ elementType    : chr "ANY"
      .. .. .. .. ..@ elementMetadata: NULL
      .. .. .. .. ..@ metadata       : list()
      .. .. ..@ seqinfo        :Formal class 'Seqinfo' [package "GenomicRanges"] with 4 slots
      .. .. .. .. ..@ seqnames   : chr "chr21"
      .. .. .. .. ..@ seqlengths : int NA
      .. .. .. .. ..@ is_circular: logi NA
      .. .. .. .. ..@ genome     : chr "hg19"
      .. .. ..@ metadata       : list()
      ..@ colData        :Formal class 'DataFrame' [package "S4Vectors"] with 6 slots
      .. .. ..@ rownames       : chr [1:8] "Sample1" "Sample2" "Sample3" "Sample4" ...
      .. .. ..@ nrows          : int 8
      .. .. ..@ listData       :List of 1
      .. .. .. ..$ SampleType: chr [1:8] "Type1" "Type1" "Type1" "Type1" ...
      .. .. ..@ elementType    : chr "ANY"
      .. .. ..@ elementMetadata: NULL
      .. .. ..@ metadata       : list()
      ..@ assays         :Reference class 'ShallowSimpleListAssays' [package "SummarizedExperiment"] with 1 field
      .. ..$ data: NULL
      .. ..and 14 methods.
      ..@ NAMES          : NULL
      ..@ elementMetadata:Formal class 'DataFrame' [package "S4Vectors"] with 6 slots
      .. .. ..@ rownames       : NULL
      .. .. ..@ nrows          : int 4243
      .. .. ..@ listData       : Named list()
      .. .. ..@ elementType    : chr "ANY"
      .. .. ..@ elementMetadata: NULL
      .. .. ..@ metadata       : list()
      ..@ metadata       : list()
    


'MethyGenoSet'



<ol class=list-inline>
	<li>'Sample1'</li>
	<li>'Sample2'</li>
	<li>'Sample3'</li>
	<li>'Sample4'</li>
	<li>'Sample5'</li>
	<li>'Sample6'</li>
	<li>'Sample7'</li>
	<li>'Sample8'</li>
</ol>




    class: MethyGenoSet 
    dim: 4243 8 
    metadata(0):
    assays(4): unmethylated methylated detection exprs
    rownames(4243): cg17035109 cg06187584 ... cg07468397 cg08821909
    rowData names(1): ID
    colnames(8): Sample1 Sample2 ... Sample7 Sample8
    colData names(1): SampleType


3、The data has two sample conditions, which are Type1 and Type2(we can see it is in the `colData` component). Smoothen the input data with a desired window size (here, it is 200).<br>
数据有两个样本条件，即Type1和Type2（我们可以看到它在`colData`组件中）。用所需的窗口大小平滑输入数据(这里是200)。


```R
methylSmooth <- smoothMethyData(exampleMethyGenoSet, winSize =200) #Might get some warning messages
attr(methylSmooth,'windowSize')
```

    Smoothing Chromosome chr21 ...
    


200


4、Extract the sample conditions from the `colData` component, and take the input data and sample type from the preceding step and detect the `DMRs` in the data. And then, take a look at the object created.<br>
从`colData`部分提取样本条件，并从上一步获取输入数据和样本类型，检测数据中的`DMRs`。然后，看看创建的对象。


```R
conditons <- colData(exampleMethyGenoSet)$SampleType                    #Extract the sample conditions from the colData
myDMR <- detectDMR.slideWin(exampleMethyGenoSet, sampleType = conditons)#Might get some warning messages
head(myDMR)
#The following screenshot shows the first six entries in the MyDMR object:

```

    Smoothing Chromosome chr21 ...
    


    GRanges object with 6 ranges and 11 metadata columns:
                 seqnames               ranges strand |    PROBEID difference
                    <Rle>            <IRanges>  <Rle> |   <factor>  <numeric>
      cg17035109    chr21 [10882029, 10882029]      * | cg17035109 -1.8411605
      cg06187584    chr21 [10883548, 10883548]      * | cg06187584 -0.4566059
      cg12459059    chr21 [10884748, 10884748]      * | cg12459059 -0.3591179
      cg25450479    chr21 [10884967, 10884967]      * | cg25450479 -0.3591179
      cg23347501    chr21 [10884969, 10884969]      * | cg23347501 -0.3591179
      cg03661019    chr21 [10885409, 10885409]      * | cg03661019 -0.3532662
                    p.value  p.adjust     tscore startWinIndex endWinIndex
                  <numeric> <numeric>  <numeric>     <numeric>   <numeric>
      cg17035109 0.06276449 0.1888488 -2.2804091             1           1
      cg06187584 0.41601486 0.6149596 -0.8734252             2           2
      cg12459059 0.36542152 0.5627904 -0.9789270             3           5
      cg25450479 0.36542152 0.5627904 -0.9789270             3           5
      cg23347501 0.36542152 0.5627904 -0.9789270             3           5
      cg03661019 0.38065600 0.5782099 -0.9460314             6           6
                 startLocation endLocation mean_Type1  mean_Type2
                     <integer>   <integer>  <numeric>   <numeric>
      cg17035109      10882029    10882029 -2.4183775 -0.57721699
      cg06187584      10883548    10883548 -2.2297567 -1.77315084
      cg12459059      10884748    10884969  0.2594151  0.61853304
      cg25450479      10884748    10884969  0.2594151  0.61853304
      cg23347501      10884748    10884969  0.2594151  0.61853304
      cg03661019      10885409    10885409 -0.4170363 -0.06377013
      -------
      seqinfo: 1 sequence from hg19 genome; no seqlengths


5、Finally, to identify the significant DMRs, run `identifySigDMR` in the output.最后，为了确定重要的DMR，在输出中运行identifySigDMR。


```R
mySigDMR <- identifySigDMR(myDMR)
mySigDMR
```


    $sigDMRInfo
    GRanges object with 14 ranges and 7 metadata columns:
           seqnames               ranges strand | NumOfProbe  min_p.value
              <Rle>            <IRanges>  <Rle> |  <integer>    <numeric>
       [1]    chr21 [19191045, 19191320]      * |          3 0.0001235924
       [2]    chr21 [34522588, 34522588]      * |          1 7.677523e-06
       [3]    chr21 [37851847, 37851847]      * |          1 9.416121e-05
       [4]    chr21 [38066047, 38066047]      * |          1 1.131734e-05
       [5]    chr21 [38075599, 38077042]      * |          4 2.609326e-07
       ...      ...                  ...    ... .        ...          ...
      [10]    chr21 [42217001, 42217001]      * |          1  7.22781e-05
      [11]    chr21 [43652704, 43652704]      * |          1 0.0002190184
      [12]    chr21 [45139229, 45139379]      * |          2 0.0002073719
      [13]    chr21 [47876058, 47876058]      * |          1 4.405417e-06
      [14]    chr21 [47878552, 47878975]      * |          5  5.00709e-07
           max_difference min_p.adjust    max_tscore        mean_Type1
                <numeric>    <numeric>     <numeric>         <numeric>
       [1]  -1.1907873545  0.007061388   -8.74702112 -4.02704733856667
       [2]  -1.2848047397  0.003093883 -14.182749104     -4.5115552816
       [3]  -1.4922516167  0.006768435  -9.179626306      0.6201243857
       [4]  -2.5480810202  0.003692065 -13.270634374      0.5772609227
       [5]  -3.5166289767 0.0005533076 -25.138160458    0.697308788975
       ...            ...          ...           ...               ...
      [10]  -2.0397570049  0.005676508  -9.618344098      0.4670764024
      [11]  -2.0082802789  0.009900444  -7.894394819      0.3099782311
      [12]  -1.0405871865  0.009900444   -7.97261521     -4.9126547368
      [13]  -2.7791642879  0.002335422  -15.59262922      0.1149510409
      [14]  -1.1235181909 0.0005676619 -22.527841219    -4.98131326204
                  mean_Type2
                   <numeric>
       [1] -2.92264438366667
       [2]      -3.226750542
       [3]       2.112376002
       [4]       3.125341943
       [5]     3.54771257025
       ...               ...
      [10]       2.506833407
      [11]        2.31825851
      [12]       -3.87206755
      [13]       2.894115329
      [14]     -3.9667573434
      -------
      seqinfo: 1 sequence from an unspecified genome; no seqlengths
    
    $sigDataInfo
    GRanges object with 21 ranges and 11 metadata columns:
                      seqnames               ranges strand |         PROBEID
                         <Rle>            <IRanges>  <Rle> |        <factor>
           cg12430776    chr21 [19191096, 19191096]      * |      cg12430776
      ch.21.33444458F    chr21 [34522588, 34522588]      * | ch.21.33444458F
           cg02417033    chr21 [37851847, 37851847]      * |      cg02417033
           cg10445315    chr21 [38066047, 38066047]      * |      cg10445315
           cg22711869    chr21 [38075599, 38075599]      * |      cg22711869
                  ...      ...                  ...    ... .             ...
           cg14522549    chr21 [45139379, 45139379]      * |      cg14522549
           cg09387528    chr21 [47876058, 47876058]      * |      cg09387528
           cg19247551    chr21 [47878552, 47878552]      * |      cg19247551
           cg15775835    chr21 [47878727, 47878727]      * |      cg15775835
           cg12533308    chr21 [47878746, 47878746]      * |      cg12533308
                      difference      p.value     p.adjust     tscore startWinIndex
                       <numeric>    <numeric>    <numeric>  <numeric>     <numeric>
           cg12430776  -1.071523 1.235924e-04  0.007061388  -8.747021           279
      ch.21.33444458F  -1.284805 7.677523e-06  0.003093883 -14.182749           998
           cg02417033  -1.492252 9.416121e-05  0.006768435  -9.179626          1486
           cg10445315  -2.548081 1.131734e-05  0.003692065 -13.270634          1514
           cg22711869  -1.552162 2.332110e-04  0.009900444  -7.805318          1555
                  ...        ...          ...          ...        ...           ...
           cg14522549  -1.040587 2.073719e-04 0.0099004442  -7.972615          2848
           cg09387528  -2.779164 4.405417e-06 0.0023354219 -15.592629          4167
           cg19247551  -1.123518 5.007090e-07 0.0005676619 -22.527841          4172
           cg15775835  -1.021457 8.754208e-06 0.0030938831 -13.867824          4172
           cg12533308  -1.021457 8.754208e-06 0.0030938831 -13.867824          4172
                      endWinIndex startLocation endLocation mean_Type1 mean_Type2
                        <numeric>     <integer>   <integer>  <numeric>  <numeric>
           cg12430776         281      19191045    19191320 -4.0922713  -3.020748
      ch.21.33444458F         998      34522588    34522588 -4.5115553  -3.226751
           cg02417033        1486      37851847    37851847  0.6201244   2.112376
           cg10445315        1514      38066047    38066047  0.5772609   3.125342
           cg22711869        1555      38075599    38075599  0.9672058   2.519367
                  ...         ...           ...         ...        ...        ...
           cg14522549        2849      45139229    45139379  -4.912655  -3.872068
           cg09387528        4167      47876058    47876058   0.114951   2.894115
           cg19247551        4174      47878552    47878746  -5.482018  -4.358500
           cg15775835        4176      47878552    47878975  -4.969093  -3.947636
           cg12533308        4176      47878552    47878975  -4.969093  -3.947636
      -------
      seqinfo: 1 sequence from hg19 genome; no seqlengths
    


6、Get the annotation information for the regions from UCSC, and finally, export the results of analysis.<br>
从UCSC获取区域的注释信息，最后导出分析结果.


```R
dmr_anno <- annotateDMRInfo(mySigDMR,'TxDb.Hsapiens.UCSC.hg19.knownGene')
export.DMRInfo(dmr_anno, savePrefix='testExample')
```

The example data we used in this part is an extension of the `eSet` class, bearing slots for the different types of information such as samples, chromosomes, and genomic ranges. The analysis begins with data smoothing. The function considers the correlation in the nearby CpG sites in terms of the defined window (here, it is 200), thus reducing the noise. With the less noisy data in place, `detectDMR.slideWin` checks whether the region (smoothed) is differentially methylated based on the `t test` or `wilcox test`, which is then merged to get the continuous regions by calling the `getContinuousRegion` function. The results annotated with the genome data are finally exported as `CSV` files in step 6. <br>

A differentially methylated region, which we aimed to find in this part, is a region where most of the CpG sites are methylated. What we did here was the detection of probe methylation status from the `exampleMethyGenoSet` object. Finally, we mapped these regions onto the chromosomal maps. In this part, we saw the SIM2, ERG, DIP2, and so on, regions being substantially methylated on chromosome 21, where the p-value was less than 0.05. The input data for the part was an `eSet` object. However, it can be in other formats such as a `CSV` or `txt` file. The approach to create `eSet` from these files has been already explained in Chapter 5, Analyzing Microarray Data with R.<br>

在这部分中使用的示例数据是eSet类的扩展，它为不同类型的信息（如样本，染色体和基因组范围）提供插槽。分析从数据处理（标准化）开始。该函数根据所定义的窗口（这里是200）考虑附近CpG站点的相关性，从而减少噪声。使用较少噪声的数据，detectDMR.slideWin会根据t检验或wilcox检验来检查区域（已平滑）是否有差异甲基化，然后通过调用getContinuousRegion函数合并以获得连续区域。基因组数据注释的结果在步骤6中最终导出为CSV文件。

我们在本部分中寻找的差异甲基化区域是大多数CpG位点被甲基化的区域。我们在这里做的是从示例MethyGenoSet对象中检测探针甲基化状态。最后，我们将这些区域映射到染色体图上。在这部分，我们看到了SIM2，ERG，DIP2等等，染色体21上的区域基本上是甲基化的，其中p值小于0.05。该零件的输入数据是一个eSet对象。但是，它可以是其他格式，例如CSV或txt文件。从这些文件中创建eSet的方法已经在第5章用R分析微阵列数据中进行了解释。


```R

```
