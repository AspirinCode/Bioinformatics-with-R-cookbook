
# 6、Analyzing RNAseq data with the edgeR package

To determine whether the count for a transcript is significantly different or differentially expressed under the treatment condition, we need to do a differential count analysis for the data.<br> 

为了确定转录本的计数是否在处理条件下显着不同或差异表达，我们需要对数据进行差异计数分析。<br>

To perform differential analysis for the RNAseq data, we will need two R packages, namely, `edgeR` and `goseq`, for the recipe, and use the dataset: the default data for these packages. The data comes from an experiment that examined the effect of androgen stimulation on a human prostate cancer cell line, LNCaP. The data has four control and three treated samples.<br>

在本节内容里，为了实现对RNA-seq数据的差异分析，我们需要两个R包，即`edgeR`和`goseq`，并使用数据集：这些包的默认数据。数据来自一项实验，该实验检测雄激素刺激对人前列腺癌细胞系LNCaP的影响。 数据有四个对照样品和三个处理过的样品。

The `edgeR` package uses the count data by modeling it via an overdispersed Poisson model, assuming it to be a negative binomial distribution. Then, it follows an Empirical Bayes procedure to moderate the degree of overdispersion across genes by conditional maximum likelihood, conditioned on the total count for that gene (step 5). To compute the differential expression of a tag/gene, a Fisher's exact test (step 6) is performed, yielding the corresponding statistical scores. Finally, the top ranking tags (according to the p-values) are returned in steps 7. The `topTags` function, by default, returns the top 10 results; we can return desired number of top ranked tags by setting the argument n to this number in the function.<br>

`edgeR`库使用计数数据，通过对超分散泊松模型进行建模，假设它是一个负二项式分布。然后，它遵循经验贝叶斯公式，通过条件最大可能性来调节基因的过度分散程度，以该基因的总计数为条件（步骤5）。为了计算标签/基因的差异表达，进行Fisher精确检验（步骤6），得到相应的统计分数。最后，排名最高的标签（根据p值）在步骤7中返回。默认情况下，`topTags`函数返回前10个结果; 我们可以通过在函数中将参数n设置为该数字来返回所需数量的排名最靠前的标签。

1、Install and load required libraries.安装并载入所需库。


```R
source("http://bioconductor.org/biocLite.R")
biocLite("edgeR")
biocLite("goseq")
```


```R
library(edgeR)
library(goseq)
```

    Loading required package: limma
    Loading required package: BiasedUrn
    Loading required package: geneLenDataBase
    
    

2、Read in the input data from the `goseq` library data directory，and take a look at the content of a part of the data.<br> 从`goseq`库中读取内置数据并查看部分数据。


```R
myData <- read.table(system.file("extdata", "Li_sum.txt", package='goseq'), sep = '\t', 
                     header = TRUE, stringsAsFactors = FALSE,row.names=1)
head(myData)
```


<table>
<thead><tr><th></th><th scope=col>lane1</th><th scope=col>lane2</th><th scope=col>lane3</th><th scope=col>lane4</th><th scope=col>lane5</th><th scope=col>lane6</th><th scope=col>lane8</th></tr></thead>
<tbody>
	<tr><th scope=row>ENSG00000215688</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000215689</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000220823</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000242499</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000224938</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000239242</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
</tbody>
</table>



4、The first four columns in the data are controls and the last three are the treatment samples. Assign these attributes to the data.<br>数据中前四行是对照组，后三行是实验组样本。给这些数据设置属性值。


```R
myTreat <- factor(rep(c("Control","Treatment"),times = c(4,3)))
```

5、Create a `DGElist` object using all the count data and treatment information.<br>
用所有计数数据和治疗信息创建一个DGElist对象。


```R
myDG <- DGEList(myData,lib.size = colSums(myData),group = myTreat)
# The `DGElist` object is a list with two components: counts and sample (treatment information).
myDG
```


<dl>
	<dt>$counts</dt>
		<dd><table>
<thead><tr><th></th><th scope=col>lane1</th><th scope=col>lane2</th><th scope=col>lane3</th><th scope=col>lane4</th><th scope=col>lane5</th><th scope=col>lane6</th><th scope=col>lane8</th></tr></thead>
<tbody>
	<tr><th scope=row>ENSG00000215688</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000215689</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000220823</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000242499</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000224938</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000239242</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000243140</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000240187</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000241444</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000242468</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000241826</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000239150</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000220023</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000238303</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000239060</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000238664</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000238600</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000234317</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000225412</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000233023</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000231129</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000233149</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000225459</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000226732</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000223844</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000229493</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000234100</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000234372</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000224230</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000225859</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>...</th><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td></tr>
	<tr><th scope=row>ENSG00000229284</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000206437</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000206308</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000232902</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000235223</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000236934</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000206439</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000235757</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000237559</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000224044</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000228104</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000183574</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000206510</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000243594</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000235814</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000206458</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000215425</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000206366</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000229353</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000242685</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000221974</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000235630</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000235952</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000229199</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000232345</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000231030</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000206296</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000212866</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000228904</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000226795</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
</tbody>
</table>
</dd>
	<dt>$samples</dt>
		<dd><table>
<thead><tr><th></th><th scope=col>group</th><th scope=col>lib.size</th><th scope=col>norm.factors</th></tr></thead>
<tbody>
	<tr><th scope=row>lane1</th><td>Control  </td><td>1178832  </td><td>1        </td></tr>
	<tr><th scope=row>lane2</th><td>Control  </td><td>1384945  </td><td>1        </td></tr>
	<tr><th scope=row>lane3</th><td>Control  </td><td>1716355  </td><td>1        </td></tr>
	<tr><th scope=row>lane4</th><td>Control  </td><td>1767927  </td><td>1        </td></tr>
	<tr><th scope=row>lane5</th><td>Treatment</td><td>2127868  </td><td>1        </td></tr>
	<tr><th scope=row>lane6</th><td>Treatment</td><td>2142158  </td><td>1        </td></tr>
	<tr><th scope=row>lane8</th><td>Treatment</td><td> 816171  </td><td>1        </td></tr>
</tbody>
</table>
</dd>
</dl>



6、Estimate the dispersion in the data, and then use Fisher's exact test.<br>
评估数据中的分散性（正态性），然后使用Fisher精确检验。


```R
myDisp <- estimateCommonDisp(myDG)
mytest <- exactTest(myDisp)
```

7、Extract the top DE tags ranked by the p-value (or the absolute log fold change) using `topTags` function, and to see the results by `as.data.frame` function.<br>
使用`topTags`函数提取按p值排列的顶级DE标记（或绝对对数倍数变化），并通过`as.data.frame`函数查看结果。


```R
myRes <- topTags(mytest, sort.by = "PValue")
summary(myRes)
as.data.frame(myRes)
```


                  Length Class      Mode     
    table         4      data.frame list     
    adjust.method 1      -none-     character
    comparison    2      -none-     character
    test          1      -none-     character



<table>
<thead><tr><th></th><th scope=col>logFC</th><th scope=col>logCPM</th><th scope=col>PValue</th><th scope=col>FDR</th></tr></thead>
<tbody>
	<tr><th scope=row>ENSG00000127954</th><td>11.557868   </td><td>6.680748    </td><td>2.574972e-80</td><td>1.274766e-75</td></tr>
	<tr><th scope=row>ENSG00000151503</th><td> 5.398963   </td><td>8.499530    </td><td>1.781732e-65</td><td>4.410321e-61</td></tr>
	<tr><th scope=row>ENSG00000096060</th><td> 4.897600   </td><td>9.446705    </td><td>7.983756e-60</td><td>1.317479e-55</td></tr>
	<tr><th scope=row>ENSG00000091879</th><td> 5.737627   </td><td>6.282646    </td><td>1.207655e-54</td><td>1.494654e-50</td></tr>
	<tr><th scope=row>ENSG00000132437</th><td>-5.880436   </td><td>7.951910    </td><td>2.950042e-52</td><td>2.920896e-48</td></tr>
	<tr><th scope=row>ENSG00000166451</th><td> 4.564246   </td><td>8.458467    </td><td>7.126763e-52</td><td>5.880292e-48</td></tr>
	<tr><th scope=row>ENSG00000131016</th><td> 5.254737   </td><td>6.607957    </td><td>1.066807e-51</td><td>7.544766e-48</td></tr>
	<tr><th scope=row>ENSG00000163492</th><td> 7.085400   </td><td>5.128514    </td><td>2.716461e-45</td><td>1.681014e-41</td></tr>
	<tr><th scope=row>ENSG00000113594</th><td> 4.051053   </td><td>8.603264    </td><td>9.272066e-44</td><td>5.100255e-40</td></tr>
	<tr><th scope=row>ENSG00000116285</th><td> 4.108522   </td><td>7.864773    </td><td>6.422468e-43</td><td>3.179507e-39</td></tr>
</tbody>
</table>




```R

```
