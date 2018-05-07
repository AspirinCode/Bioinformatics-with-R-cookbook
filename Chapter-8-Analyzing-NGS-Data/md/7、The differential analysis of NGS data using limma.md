
# 7、The differential analysis of NGS data using limma

In previous part, we have learned how to find differentially expressed genes in analyzing microarray data by `limma`. It can handle multiple experiments via Empirical Bayes statistical methods and uses normalized read counts for each gene. This section will explain the use of the `limma` package for differential gene analysis with NGS data.<br>

在前面的部分中，我们已经学习如何通过`limma`分析微阵列数据中的差异表达基因。它可以通过Empirical Bayes统计方法处理多个实验，并使用每个基因的片段标准化计数。本节将解释使用`limma`软件包与NGS数据进行差异基因分析的方法。

Besides using the `limma` library, we will use the `Pasilla` dataset here. The dataset can be obtained from Bioconductor and consists of sequence counts from a perturbation experiment in Drosophila. To know more about the data, refer to [the Conservation of an RNA regulatory map between Drosophila and mammals article by Brooks and others](http://genome.cshlp.org/content/early/2010/10/04/gr.108662.110).<br>

除了使用`limma`库，我们将在这里使用`Pasilla`数据集。数据集可以从Bioconductor获得，并且由来自果蝇扰动实验的序列计数组成。要了解更多关于这些数据的信息，请参阅[Brooks等人在果蝇和哺乳动物之间的RNA调控图的保存](http://genome.cshlp.org/content/early/2010/10/04/gr.108662.110)。

1、Install and load the packages and data required, and check the data.


```R
source("http://bioconductor.org/biocLite.R")
biocLite(c("DESeq","pasilla"))               # install
biocLite("survival")                         # install
```


```R
library(limma)                #load 
library(DESeq)                #load
library(pasilla)              #load
library(survival)             #load
data(pasillaGenes)            #check the data
```

2、Create an expression set using the counts from the `pasillaGenes` dataset and take a look at it, we will find seven samples in the dataset: the first three are the treatment samples and the last four are controls. Assign this 7 samples to our `eset` data, and check it.<br>
使用来自“pasillaGenes”数据集的计数创建一个表达集，并查看它，我们将在数据集中找到七个样本：前三个是处理样本，后四个是对照。把它分配给我们的`eset`数据，并检查它。


```R
eset <- counts(pasillaGenes)      #create
eset                              #take a look
colnames(eset) <- c(paste("T", 1:3, sep="_"), paste("C", 1:4, sep="_"))   #assign the data
```


<table>
<thead><tr><th></th><th scope=col>treated1fb</th><th scope=col>treated2fb</th><th scope=col>treated3fb</th><th scope=col>untreated1fb</th><th scope=col>untreated2fb</th><th scope=col>untreated3fb</th><th scope=col>untreated4fb</th></tr></thead>
<tbody>
	<tr><th scope=row>FBgn0000003</th><td>    0</td><td>    0</td><td>    1</td><td>    0</td><td>    0</td><td>    0</td><td>    0</td></tr>
	<tr><th scope=row>FBgn0000008</th><td>   78</td><td>   46</td><td>   43</td><td>   47</td><td>   89</td><td>   53</td><td>   27</td></tr>
	<tr><th scope=row>FBgn0000014</th><td>    2</td><td>    0</td><td>    0</td><td>    0</td><td>    0</td><td>    1</td><td>    0</td></tr>
	<tr><th scope=row>FBgn0000015</th><td>    1</td><td>    0</td><td>    1</td><td>    0</td><td>    1</td><td>    1</td><td>    2</td></tr>
	<tr><th scope=row>FBgn0000017</th><td> 3187</td><td> 1672</td><td> 1859</td><td> 2445</td><td> 4615</td><td> 2063</td><td> 1711</td></tr>
	<tr><th scope=row>FBgn0000018</th><td>  369</td><td>  150</td><td>  176</td><td>  288</td><td>  383</td><td>  135</td><td>  174</td></tr>
	<tr><th scope=row>FBgn0000022</th><td>    0</td><td>    0</td><td>    0</td><td>    0</td><td>    1</td><td>    0</td><td>    0</td></tr>
	<tr><th scope=row>FBgn0000024</th><td>    4</td><td>    5</td><td>    3</td><td>    4</td><td>    7</td><td>    1</td><td>    0</td></tr>
	<tr><th scope=row>FBgn0000028</th><td>    0</td><td>    1</td><td>    1</td><td>    0</td><td>    1</td><td>    0</td><td>    0</td></tr>
	<tr><th scope=row>FBgn0000032</th><td>  942</td><td>  465</td><td>  536</td><td>  767</td><td>  956</td><td>  464</td><td>  471</td></tr>
	<tr><th scope=row>FBgn0000036</th><td>    1</td><td>    0</td><td>    0</td><td>    1</td><td>    1</td><td>    0</td><td>    0</td></tr>
	<tr><th scope=row>FBgn0000037</th><td>   10</td><td>   10</td><td>   11</td><td>    7</td><td>   13</td><td>    3</td><td>    5</td></tr>
	<tr><th scope=row>FBgn0000038</th><td>    0</td><td>    0</td><td>    0</td><td>    0</td><td>    0</td><td>    0</td><td>    0</td></tr>
	<tr><th scope=row>FBgn0000039</th><td>    0</td><td>    0</td><td>    0</td><td>    0</td><td>    1</td><td>    0</td><td>    0</td></tr>
	<tr><th scope=row>FBgn0000042</th><td>64357</td><td>39356</td><td>43480</td><td>52791</td><td>62319</td><td>24399</td><td>27687</td></tr>
	<tr><th scope=row>FBgn0000043</th><td>28458</td><td>15671</td><td>17180</td><td>16485</td><td>21077</td><td> 8376</td><td> 9284</td></tr>
	<tr><th scope=row>FBgn0000044+FBgn0065059</th><td>   30</td><td>   17</td><td>   17</td><td>   11</td><td>   43</td><td>   10</td><td>   11</td></tr>
	<tr><th scope=row>FBgn0000045</th><td>    8</td><td>    4</td><td>    8</td><td>    5</td><td>   10</td><td>    1</td><td>    1</td></tr>
	<tr><th scope=row>FBgn0000046</th><td>   36</td><td>   17</td><td>   17</td><td>    5</td><td>   33</td><td>   14</td><td>   12</td></tr>
	<tr><th scope=row>FBgn0000047</th><td>    8</td><td>    7</td><td>    7</td><td>    6</td><td>    9</td><td>    3</td><td>    2</td></tr>
	<tr><th scope=row>FBgn0000052</th><td> 1647</td><td>  996</td><td> 1043</td><td> 1491</td><td> 1811</td><td>  899</td><td> 1039</td></tr>
	<tr><th scope=row>FBgn0000053</th><td> 1603</td><td> 1110</td><td> 1222</td><td> 1108</td><td> 1531</td><td>  768</td><td>  877</td></tr>
	<tr><th scope=row>FBgn0000054</th><td>  542</td><td>  314</td><td>  382</td><td>  479</td><td>  588</td><td>  302</td><td>  317</td></tr>
	<tr><th scope=row>FBgn0000055+FBgn0000056</th><td>    6</td><td>    2</td><td>    4</td><td>    4</td><td>    7</td><td>    0</td><td>    6</td></tr>
	<tr><th scope=row>FBgn0000057</th><td>  513</td><td>  321</td><td>  311</td><td>  389</td><td>  596</td><td>  266</td><td>  287</td></tr>
	<tr><th scope=row>FBgn0000061</th><td>    2</td><td>    0</td><td>    0</td><td>    0</td><td>    1</td><td>    1</td><td>    4</td></tr>
	<tr><th scope=row>FBgn0000063</th><td>  383</td><td>  208</td><td>  203</td><td>  268</td><td>  395</td><td>  153</td><td>  159</td></tr>
	<tr><th scope=row>FBgn0000064</th><td> 6701</td><td> 4199</td><td> 5148</td><td> 4139</td><td> 6265</td><td> 3360</td><td> 3342</td></tr>
	<tr><th scope=row>FBgn0000071</th><td>  533</td><td>  266</td><td>  311</td><td>   32</td><td>  105</td><td>   48</td><td>   47</td></tr>
	<tr><th scope=row>FBgn0000075</th><td>    2</td><td>    0</td><td>    1</td><td>    0</td><td>    2</td><td>    2</td><td>    1</td></tr>
	<tr><th scope=row>...</th><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td></tr>
	<tr><th scope=row>FBgn0261535</th><td> 470</td><td> 261</td><td> 299</td><td> 428</td><td> 562</td><td> 232</td><td> 263</td></tr>
	<tr><th scope=row>FBgn0261538+FBgn0050080</th><td>   0</td><td>   0</td><td>   0</td><td>   0</td><td>   1</td><td>   1</td><td>   1</td></tr>
	<tr><th scope=row>FBgn0261545</th><td>1278</td><td> 581</td><td> 766</td><td> 942</td><td>1939</td><td> 946</td><td> 904</td></tr>
	<tr><th scope=row>FBgn0261546</th><td>  12</td><td>   5</td><td>   7</td><td>  11</td><td>  43</td><td>  22</td><td>  20</td></tr>
	<tr><th scope=row>FBgn0261547</th><td> 437</td><td> 252</td><td> 242</td><td> 402</td><td> 760</td><td> 370</td><td> 291</td></tr>
	<tr><th scope=row>FBgn0261548</th><td>1563</td><td> 742</td><td> 859</td><td> 788</td><td>1655</td><td> 679</td><td> 657</td></tr>
	<tr><th scope=row>FBgn0261549</th><td> 528</td><td> 253</td><td> 378</td><td> 225</td><td> 561</td><td> 283</td><td> 254</td></tr>
	<tr><th scope=row>FBgn0261550</th><td>6179</td><td>4185</td><td>4853</td><td>6276</td><td>8639</td><td>5194</td><td>4922</td></tr>
	<tr><th scope=row>FBgn0261551</th><td>3985</td><td>2556</td><td>2770</td><td>2346</td><td>5256</td><td>2449</td><td>2158</td></tr>
	<tr><th scope=row>FBgn0261552</th><td>2206</td><td> 739</td><td> 804</td><td>3662</td><td>8839</td><td>3285</td><td>2751</td></tr>
	<tr><th scope=row>FBgn0261553</th><td> 934</td><td> 573</td><td> 571</td><td> 675</td><td>1433</td><td> 678</td><td> 650</td></tr>
	<tr><th scope=row>FBgn0261554</th><td> 917</td><td> 523</td><td> 574</td><td> 749</td><td>1003</td><td> 548</td><td> 473</td></tr>
	<tr><th scope=row>FBgn0261555</th><td>  63</td><td>  51</td><td>  60</td><td>  38</td><td>  63</td><td>  20</td><td>  19</td></tr>
	<tr><th scope=row>FBgn0261556</th><td>2218</td><td>1564</td><td>1796</td><td>1610</td><td>2332</td><td>1305</td><td>1121</td></tr>
	<tr><th scope=row>FBgn0261560</th><td> 997</td><td> 488</td><td> 588</td><td>1217</td><td>2501</td><td> 743</td><td> 915</td></tr>
	<tr><th scope=row>FBgn0261561</th><td>  95</td><td>  40</td><td>  55</td><td>  40</td><td>  49</td><td>  29</td><td>  36</td></tr>
	<tr><th scope=row>FBgn0261562</th><td> 298</td><td> 185</td><td> 241</td><td> 251</td><td> 400</td><td> 180</td><td> 200</td></tr>
	<tr><th scope=row>FBgn0261563</th><td>1918</td><td>1121</td><td>1117</td><td> 596</td><td>1353</td><td> 687</td><td> 650</td></tr>
	<tr><th scope=row>FBgn0261564</th><td>1033</td><td> 564</td><td> 768</td><td> 724</td><td>1233</td><td> 754</td><td> 687</td></tr>
	<tr><th scope=row>FBgn0261565</th><td> 514</td><td> 335</td><td> 332</td><td> 453</td><td> 590</td><td> 345</td><td> 328</td></tr>
	<tr><th scope=row>FBgn0261566</th><td>   0</td><td>   0</td><td>   0</td><td>   0</td><td>   0</td><td>   0</td><td>   0</td></tr>
	<tr><th scope=row>FBgn0261567</th><td>   0</td><td>   0</td><td>   0</td><td>   0</td><td>   0</td><td>   0</td><td>   0</td></tr>
	<tr><th scope=row>FBgn0261568</th><td>   0</td><td>   0</td><td>   1</td><td>   0</td><td>   0</td><td>   0</td><td>   0</td></tr>
	<tr><th scope=row>FBgn0261569</th><td>   1</td><td>   0</td><td>   0</td><td>   0</td><td>   0</td><td>   0</td><td>   0</td></tr>
	<tr><th scope=row>FBgn0261570</th><td>2531</td><td>1739</td><td>1880</td><td>1712</td><td>2566</td><td>1316</td><td>1203</td></tr>
	<tr><th scope=row>FBgn0261571</th><td>   1</td><td>   0</td><td>   0</td><td>   0</td><td>   0</td><td>   0</td><td>   0</td></tr>
	<tr><th scope=row>FBgn0261572</th><td>   1</td><td>   1</td><td>   3</td><td>   2</td><td>   6</td><td>   3</td><td>   6</td></tr>
	<tr><th scope=row>FBgn0261573</th><td>1835</td><td>1171</td><td>1336</td><td>1405</td><td>2003</td><td>1199</td><td>1081</td></tr>
	<tr><th scope=row>FBgn0261574</th><td>5412</td><td>2046</td><td>1850</td><td>3262</td><td>4962</td><td>1876</td><td>1603</td></tr>
	<tr><th scope=row>FBgn0261575</th><td>  12</td><td>   1</td><td>   3</td><td>   5</td><td>   8</td><td>   1</td><td>   0</td></tr>
</tbody>
</table>




```R
head(eset)   # check the assigned data
```


<table>
<thead><tr><th></th><th scope=col>T_1</th><th scope=col>T_2</th><th scope=col>T_3</th><th scope=col>C_1</th><th scope=col>C_2</th><th scope=col>C_3</th><th scope=col>C_4</th></tr></thead>
<tbody>
	<tr><th scope=row>FBgn0000003</th><td>   0</td><td>   0</td><td>   1</td><td>   0</td><td>   0</td><td>   0</td><td>   0</td></tr>
	<tr><th scope=row>FBgn0000008</th><td>  78</td><td>  46</td><td>  43</td><td>  47</td><td>  89</td><td>  53</td><td>  27</td></tr>
	<tr><th scope=row>FBgn0000014</th><td>   2</td><td>   0</td><td>   0</td><td>   0</td><td>   0</td><td>   1</td><td>   0</td></tr>
	<tr><th scope=row>FBgn0000015</th><td>   1</td><td>   0</td><td>   1</td><td>   0</td><td>   1</td><td>   1</td><td>   2</td></tr>
	<tr><th scope=row>FBgn0000017</th><td>3187</td><td>1672</td><td>1859</td><td>2445</td><td>4615</td><td>2063</td><td>1711</td></tr>
	<tr><th scope=row>FBgn0000018</th><td> 369</td><td> 150</td><td> 176</td><td> 288</td><td> 383</td><td> 135</td><td> 174</td></tr>
</tbody>
</table>



3、Create a design matrix，and perform a `voom` transformation using the experiment design matrix.<br>
创建一个设计矩阵，并使用实验设计矩阵执行`voom`变换.


```R
design <- cbind(Intercept=1,trt=c(1,1,1,0,0,0,0))    # add a new column to the `eset` data
eset_voom <- voom(eset, design, plot=FALSE)          # perform voom transformation
eset_voom
```


<dl>
	<dt>$E</dt>
		<dd><table>
<thead><tr><th></th><th scope=col>T_1</th><th scope=col>T_2</th><th scope=col>T_3</th><th scope=col>C_1</th><th scope=col>C_2</th><th scope=col>C_3</th><th scope=col>C_4</th></tr></thead>
<tbody>
	<tr><th scope=row>FBgn0000003</th><td>-4.30792027</td><td>-3.50962618</td><td>-2.1190901 </td><td>-3.87180837</td><td>-4.5531169 </td><td>-3.4137300 </td><td>-3.5722974 </td></tr>
	<tr><th scope=row>FBgn0000008</th><td> 2.98670048</td><td> 3.02953264</td><td> 2.7388909 </td><td> 2.69804724</td><td> 2.9306989 </td><td> 3.3277370 </td><td> 2.2090623 </td></tr>
	<tr><th scope=row>FBgn0000014</th><td>-1.98599217</td><td>-3.50962618</td><td>-3.7040526 </td><td>-3.87180837</td><td>-4.5531169 </td><td>-1.8287675 </td><td>-3.5722974 </td></tr>
	<tr><th scope=row>FBgn0000015</th><td>-2.72295777</td><td>-3.50962618</td><td>-2.1190901 </td><td>-3.87180837</td><td>-2.9681544 </td><td>-1.8287675 </td><td>-1.2503693 </td></tr>
	<tr><th scope=row>FBgn0000017</th><td> 8.33028936</td><td> 8.19816432</td><td> 8.1566465 </td><td> 8.38410538</td><td> 8.6191543 </td><td> 8.5971477 </td><td> 8.1687481 </td></tr>
	<tr><th scope=row>FBgn0000018</th><td> 5.22151028</td><td> 4.72399350</td><td> 4.7594718 </td><td> 5.30061914</td><td> 5.0299658 </td><td> 4.6684190 </td><td> 4.8747858 </td></tr>
	<tr><th scope=row>FBgn0000022</th><td>-4.30792027</td><td>-3.50962618</td><td>-3.7040526 </td><td>-3.87180837</td><td>-2.9681544 </td><td>-3.4137300 </td><td>-3.5722974 </td></tr>
	<tr><th scope=row>FBgn0000024</th><td>-1.13799527</td><td>-0.05019456</td><td>-0.8966977 </td><td>-0.70188337</td><td>-0.6462263 </td><td>-1.8287675 </td><td>-3.5722974 </td></tr>
	<tr><th scope=row>FBgn0000028</th><td>-4.30792027</td><td>-1.92466367</td><td>-2.1190901 </td><td>-3.87180837</td><td>-2.9681544 </td><td>-3.4137300 </td><td>-3.5722974 </td></tr>
	<tr><th scope=row>FBgn0000032</th><td> 6.57242854</td><td> 6.35301118</td><td> 6.3633818 </td><td> 6.71221457</td><td> 6.3485042 </td><td> 6.4458048 </td><td> 6.3088165 </td></tr>
	<tr><th scope=row>FBgn0000036</th><td>-2.72295777</td><td>-3.50962618</td><td>-3.7040526 </td><td>-2.28684587</td><td>-2.9681544 </td><td>-3.4137300 </td><td>-3.5722974 </td></tr>
	<tr><th scope=row>FBgn0000037</th><td> 0.08439715</td><td> 0.88269125</td><td> 0.8195094 </td><td> 0.03508223</td><td> 0.2017706 </td><td>-0.6063751 </td><td>-0.1128658 </td></tr>
	<tr><th scope=row>FBgn0000038</th><td>-4.30792027</td><td>-3.50962618</td><td>-3.7040526 </td><td>-3.87180837</td><td>-4.5531169 </td><td>-3.4137300 </td><td>-3.5722974 </td></tr>
	<tr><th scope=row>FBgn0000039</th><td>-4.30792027</td><td>-3.50962618</td><td>-3.7040526 </td><td>-3.87180837</td><td>-2.9681544 </td><td>-3.4137300 </td><td>-3.5722974 </td></tr>
	<tr><th scope=row>FBgn0000042</th><td>12.66590040</td><td>12.75468813</td><td>12.7040283 </td><td>12.81620967</td><td>12.3742791 </td><td>12.1608340 </td><td>12.1846497 </td></tr>
	<tr><th scope=row>FBgn0000043</th><td>11.48865173</td><td>11.42622948</td><td>11.3644318 </td><td>11.13710165</td><td>10.8102992 </td><td>10.6184018 </td><td>10.6083111 </td></tr>
	<tr><th scope=row>FBgn0000044+FBgn0065059</th><td> 1.62281707</td><td> 1.61965684</td><td> 1.4252304 </td><td> 0.65175359</td><td> 1.8898266 </td><td> 0.9785874 </td><td> 0.9512645 </td></tr>
	<tr><th scope=row>FBgn0000045</th><td>-0.22045743</td><td>-0.33970117</td><td> 0.3834103 </td><td>-0.41237675</td><td>-0.1607995 </td><td>-1.8287675 </td><td>-1.9873349 </td></tr>
	<tr><th scope=row>FBgn0000046</th><td> 1.88190429</td><td> 1.61965684</td><td> 1.4252304 </td><td>-0.41237675</td><td> 1.5129723 </td><td> 1.4442510 </td><td> 1.0715588 </td></tr>
	<tr><th scope=row>FBgn0000047</th><td>-0.22045743</td><td> 0.39726442</td><td> 0.2028380 </td><td>-0.17136865</td><td>-0.3051894 </td><td>-0.6063751 </td><td>-1.2503693 </td></tr>
	<tr><th scope=row>FBgn0000052</th><td> 7.37814248</td><td> 7.45109982</td><td> 7.3231623 </td><td> 7.67073989</td><td> 7.2698522 </td><td> 7.3992495 </td><td> 7.4493766 </td></tr>
	<tr><th scope=row>FBgn0000053</th><td> 7.33908837</td><td> 7.60736750</td><td> 7.5515662 </td><td> 7.24258469</td><td> 7.0276127 </td><td> 7.1721714 </td><td> 7.2049579 </td></tr>
	<tr><th scope=row>FBgn0000054</th><td> 5.77555906</td><td> 5.78729003</td><td> 5.8752634 </td><td> 6.03357864</td><td> 5.6477817 </td><td> 5.8270613 </td><td> 5.7383153 </td></tr>
	<tr><th scope=row>FBgn0000055+FBgn0000056</th><td>-0.60748055</td><td>-1.18769808</td><td>-0.5341276 </td><td>-0.70188337</td><td>-0.6462263 </td><td>-3.4137300 </td><td> 0.1281423 </td></tr>
	<tr><th scope=row>FBgn0000057</th><td> 5.69630020</td><td> 5.81904875</td><td> 5.5790358 </td><td> 5.73367115</td><td> 5.6672614 </td><td> 5.6442617 </td><td> 5.5951207 </td></tr>
	<tr><th scope=row>FBgn0000061</th><td>-1.98599217</td><td>-3.50962618</td><td>-3.7040526 </td><td>-3.87180837</td><td>-2.9681544 </td><td>-1.8287675 </td><td>-0.4023724 </td></tr>
	<tr><th scope=row>FBgn0000063</th><td> 5.27516250</td><td> 5.19427740</td><td> 4.9648324 </td><td> 5.19696991</td><td> 5.0744170 </td><td> 4.8483648 </td><td> 4.7451152 </td></tr>
	<tr><th scope=row>FBgn0000064</th><td> 9.40234807</td><td> 9.52637568</td><td> 9.6258839 </td><td> 9.14343243</td><td> 9.0600970 </td><td> 9.3007302 </td><td> 9.1344144 </td></tr>
	<tr><th scope=row>FBgn0000071</th><td> 5.75142419</td><td> 5.54836555</td><td> 5.5790358 </td><td> 2.15055944</td><td> 3.1679823 </td><td> 3.1861828 </td><td> 2.9975582 </td></tr>
	<tr><th scope=row>FBgn0000075</th><td>-1.98599217</td><td>-3.50962618</td><td>-2.1190901 </td><td>-3.87180837</td><td>-2.2311888 </td><td>-1.0918019 </td><td>-1.9873349 </td></tr>
	<tr><th scope=row>...</th><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td></tr>
	<tr><th scope=row>FBgn0261535</th><td> 5.5701306 </td><td> 5.52104096</td><td> 5.5223596 </td><td> 5.8713430 </td><td> 5.5825924 </td><td> 5.4473569 </td><td> 5.4693617 </td></tr>
	<tr><th scope=row>FBgn0261538+FBgn0050080</th><td>-4.3079203 </td><td>-3.50962618</td><td>-3.7040526 </td><td>-3.8718084 </td><td>-2.9681544 </td><td>-1.8287675 </td><td>-1.9873349 </td></tr>
	<tr><th scope=row>FBgn0261545</th><td> 7.0123162 </td><td> 6.67400921</td><td> 6.8780894 </td><td> 7.0085404 </td><td> 7.3683521 </td><td> 7.4727287 </td><td> 7.2486793 </td></tr>
	<tr><th scope=row>FBgn0261546</th><td> 0.3359359 </td><td>-0.05019456</td><td> 0.2028380 </td><td> 0.6517536 </td><td> 1.8898266 </td><td> 2.0781231 </td><td> 1.7852546 </td></tr>
	<tr><th scope=row>FBgn0261547</th><td> 5.4652189 </td><td> 5.47051340</td><td> 5.2177884 </td><td> 5.7810366 </td><td> 6.0176875 </td><td> 6.1195997 </td><td> 5.6150546 </td></tr>
	<tr><th scope=row>FBgn0261548</th><td> 7.3026432 </td><td> 7.02662104</td><td> 7.0433012 </td><td> 6.7511586 </td><td> 7.1399344 </td><td> 6.9945997 </td><td> 6.7885496 </td></tr>
	<tr><th scope=row>FBgn0261549</th><td> 5.7378394 </td><td> 5.47621576</td><td> 5.8600969 </td><td> 4.9451753 </td><td> 5.5800253 </td><td> 5.7334749 </td><td> 5.4192244 </td></tr>
	<tr><th scope=row>FBgn0261550</th><td> 9.2853541 </td><td> 9.52155809</td><td> 9.5407572 </td><td> 9.7439362 </td><td> 9.5236152 </td><td> 9.9290392 </td><td> 9.6928781 </td></tr>
	<tr><th scope=row>FBgn0261551</th><td> 8.6526247 </td><td> 8.81032814</td><td> 8.7318781 </td><td> 8.3244864 </td><td> 8.8067699 </td><td> 8.8445416 </td><td> 8.5035159 </td></tr>
	<tr><th scope=row>FBgn0261552</th><td> 7.7996238 </td><td> 7.02078016</td><td> 6.9478960 </td><td> 8.9668047 </td><td> 9.5566321 </td><td> 9.2681672 </td><td> 8.8537052 </td></tr>
	<tr><th scope=row>FBgn0261553</th><td> 6.5601306 </td><td> 6.65402350</td><td> 6.4545571 </td><td> 6.5280036 </td><td> 6.9322093 </td><td> 6.9924750 </td><td> 6.7731078 </td></tr>
	<tr><th scope=row>FBgn0261554</th><td> 6.5336441 </td><td> 6.52241955</td><td> 6.4621105 </td><td> 6.6779763 </td><td> 6.4177080 </td><td> 6.6856178 </td><td> 6.3149232 </td></tr>
	<tr><th scope=row>FBgn0261555</th><td> 2.6807644 </td><td> 3.17687435</td><td> 3.2148107 </td><td> 2.3949782 </td><td> 2.4355678 </td><td> 1.9438220 </td><td> 1.7131048 </td></tr>
	<tr><th scope=row>FBgn0261556</th><td> 7.8074486 </td><td> 8.10185977</td><td> 8.1069206 </td><td> 7.7814846 </td><td> 7.6345444 </td><td> 7.9366567 </td><td> 7.5589165 </td></tr>
	<tr><th scope=row>FBgn0261560</th><td> 6.6542528 </td><td> 6.42258858</td><td> 6.4968460 </td><td> 7.3778977 </td><td> 7.7354608 </td><td> 7.1244589 </td><td> 7.2661186 </td></tr>
	<tr><th scope=row>FBgn0261561</th><td> 3.2695086 </td><td> 2.83022383</td><td> 3.0903633 </td><td> 2.4680416 </td><td> 2.0762397 </td><td> 2.4689130 </td><td> 2.6175271 </td></tr>
	<tr><th scope=row>FBgn0261562</th><td> 4.9136669 </td><td> 5.02564920</td><td> 5.2118268 </td><td> 5.1026062 </td><td> 5.0925415 </td><td> 5.0821250 </td><td> 5.0751610 </td></tr>
	<tr><th scope=row>FBgn0261563</th><td> 7.5978428 </td><td> 7.62158773</td><td> 7.4220065 </td><td> 6.3485700 </td><td> 6.8493623 </td><td> 7.0114859 </td><td> 6.7731078 </td></tr>
	<tr><th scope=row>FBgn0261564</th><td> 6.7054024 </td><td> 6.63120360</td><td> 6.8818489 </td><td> 6.6290335 </td><td> 6.7154251 </td><td> 7.1456471 </td><td> 6.8529185 </td></tr>
	<tr><th scope=row>FBgn0261565</th><td> 5.6991070 </td><td> 5.88054278</td><td> 5.6731579 </td><td> 5.9531504 </td><td> 5.6526763 </td><td> 6.0188119 </td><td> 5.7874521 </td></tr>
	<tr><th scope=row>FBgn0261566</th><td>-4.3079203 </td><td>-3.50962618</td><td>-3.7040526 </td><td>-3.8718084 </td><td>-4.5531169 </td><td>-3.4137300 </td><td>-3.5722974 </td></tr>
	<tr><th scope=row>FBgn0261567</th><td>-4.3079203 </td><td>-3.50962618</td><td>-3.7040526 </td><td>-3.8718084 </td><td>-4.5531169 </td><td>-3.4137300 </td><td>-3.5722974 </td></tr>
	<tr><th scope=row>FBgn0261568</th><td>-4.3079203 </td><td>-3.50962618</td><td>-2.1190901 </td><td>-3.8718084 </td><td>-4.5531169 </td><td>-3.4137300 </td><td>-3.5722974 </td></tr>
	<tr><th scope=row>FBgn0261569</th><td>-2.7229578 </td><td>-3.50962618</td><td>-3.7040526 </td><td>-3.8718084 </td><td>-4.5531169 </td><td>-3.4137300 </td><td>-3.5722974 </td></tr>
	<tr><th scope=row>FBgn0261570</th><td> 7.9978565 </td><td> 8.25483079</td><td> 8.1728480 </td><td> 7.8700799 </td><td> 7.7724696 </td><td> 7.9487618 </td><td> 7.6607230 </td></tr>
	<tr><th scope=row>FBgn0261571</th><td>-2.7229578 </td><td>-3.50962618</td><td>-3.7040526 </td><td>-3.8718084 </td><td>-4.5531169 </td><td>-3.4137300 </td><td>-3.5722974 </td></tr>
	<tr><th scope=row>FBgn0261572</th><td>-2.7229578 </td><td>-1.92466367</td><td>-0.8966977 </td><td>-1.5498803 </td><td>-0.8526772 </td><td>-0.6063751 </td><td> 0.1281423 </td></tr>
	<tr><th scope=row>FBgn0261573</th><td> 7.5340371 </td><td> 7.68451506</td><td> 7.6801915 </td><td> 7.5850594 </td><td> 7.4151899 </td><td> 7.8144874 </td><td> 7.5065205 </td></tr>
	<tr><th scope=row>FBgn0261574</th><td> 9.0941591 </td><td> 8.48931678</td><td> 8.1496468 </td><td> 8.7999538 </td><td> 8.7237345 </td><td> 8.4600986 </td><td> 8.0747112 </td></tr>
	<tr><th scope=row>FBgn0261575</th><td> 0.3359359 </td><td>-1.92466367</td><td>-0.8966977 </td><td>-0.4123768 </td><td>-0.4656541 </td><td>-1.8287675 </td><td>-3.5722974 </td></tr>
</tbody>
</table>
</dd>
	<dt>$weights</dt>
		<dd><table>
<tbody>
	<tr><td> 1.366814</td><td> 1.235393</td><td> 1.260553</td><td> 1.235393</td><td> 1.311863</td><td> 1.235393</td><td> 1.235393</td></tr>
	<tr><td>10.828358</td><td> 7.349417</td><td> 8.086271</td><td> 8.247937</td><td>11.457612</td><td> 6.594662</td><td> 7.122025</td></tr>
	<tr><td> 1.417297</td><td> 1.268646</td><td> 1.300879</td><td> 1.264253</td><td> 1.387106</td><td> 1.235393</td><td> 1.235393</td></tr>
	<tr><td> 1.480700</td><td> 1.316390</td><td> 1.352282</td><td> 1.450354</td><td> 1.622582</td><td> 1.354920</td><td> 1.386275</td></tr>
	<tr><td>40.774463</td><td>41.077587</td><td>41.173969</td><td>40.990275</td><td>40.295368</td><td>41.160309</td><td>41.171015</td></tr>
	<tr><td>25.975332</td><td>18.848497</td><td>20.518924</td><td>22.606543</td><td>28.748978</td><td>18.606583</td><td>19.962593</td></tr>
	<tr><td> 1.272510</td><td> 1.235393</td><td> 1.235393</td><td> 1.264253</td><td> 1.387106</td><td> 1.235393</td><td> 1.235393</td></tr>
	<tr><td> 2.209454</td><td> 1.863157</td><td> 1.938323</td><td> 1.654535</td><td> 1.882236</td><td> 1.529991</td><td> 1.570957</td></tr>
	<tr><td> 1.480700</td><td> 1.316390</td><td> 1.352282</td><td> 1.264253</td><td> 1.387106</td><td> 1.235393</td><td> 1.235393</td></tr>
	<tr><td>37.508934</td><td>32.245278</td><td>33.698287</td><td>35.014641</td><td>38.881477</td><td>31.676404</td><td>32.908108</td></tr>
	<tr><td> 1.366814</td><td> 1.235393</td><td> 1.260553</td><td> 1.331859</td><td> 1.473341</td><td> 1.254535</td><td> 1.279756</td></tr>
	<tr><td> 3.607259</td><td> 2.602963</td><td> 2.811078</td><td> 2.287815</td><td> 2.961572</td><td> 2.058331</td><td> 2.130329</td></tr>
	<tr><td> 1.272510</td><td> 1.235393</td><td> 1.235393</td><td> 1.235393</td><td> 1.311863</td><td> 1.235393</td><td> 1.235393</td></tr>
	<tr><td> 1.272510</td><td> 1.235393</td><td> 1.235393</td><td> 1.264253</td><td> 1.387106</td><td> 1.235393</td><td> 1.235393</td></tr>
	<tr><td>33.968460</td><td>35.761589</td><td>35.340076</td><td>35.680305</td><td>34.153778</td><td>36.634792</td><td>36.311559</td></tr>
	<tr><td>36.759193</td><td>38.164831</td><td>37.857925</td><td>38.556839</td><td>37.487475</td><td>39.118364</td><td>38.936766</td></tr>
	<tr><td> 5.599681</td><td> 3.872632</td><td> 4.228295</td><td> 3.745543</td><td> 5.113741</td><td> 3.083514</td><td> 3.292745</td></tr>
	<tr><td> 2.756531</td><td> 2.129867</td><td> 2.226547</td><td> 1.848344</td><td> 2.132241</td><td> 1.694311</td><td> 1.744646</td></tr>
	<tr><td> 5.835252</td><td> 4.025241</td><td> 4.398600</td><td> 3.413867</td><td> 4.629524</td><td> 2.831086</td><td> 3.014824</td></tr>
	<tr><td> 2.963912</td><td> 2.221770</td><td> 2.346556</td><td> 2.056273</td><td> 2.466419</td><td> 1.869001</td><td> 1.930321</td></tr>
	<tr><td>41.030835</td><td>38.329708</td><td>39.258457</td><td>40.167808</td><td>41.176852</td><td>38.160421</td><td>38.939882</td></tr>
	<tr><td>41.148209</td><td>38.890772</td><td>39.734147</td><td>39.004824</td><td>41.054126</td><td>36.589718</td><td>37.489703</td></tr>
	<tr><td>33.607318</td><td>26.986503</td><td>28.706386</td><td>30.148560</td><td>35.271436</td><td>26.116175</td><td>27.559056</td></tr>
	<tr><td> 2.167949</td><td> 1.833270</td><td> 1.906169</td><td> 1.826095</td><td> 2.103947</td><td> 1.675737</td><td> 1.724935</td></tr>
	<tr><td>32.753414</td><td>25.958087</td><td>27.714290</td><td>28.840695</td><td>34.260057</td><td>24.730117</td><td>26.177514</td></tr>
	<tr><td> 1.417297</td><td> 1.268646</td><td> 1.300879</td><td> 1.499801</td><td> 1.684835</td><td> 1.397163</td><td> 1.430896</td></tr>
	<tr><td>28.146166</td><td>20.948468</td><td>22.685609</td><td>22.585981</td><td>28.729206</td><td>18.588203</td><td>19.943827</td></tr>
	<tr><td>39.511214</td><td>40.263108</td><td>40.075903</td><td>40.259566</td><td>39.619452</td><td>40.734356</td><td>40.565262</td></tr>
	<tr><td>32.207365</td><td>25.290121</td><td>27.059930</td><td> 8.593260</td><td>11.925130</td><td> 6.869756</td><td> 7.421487</td></tr>
	<tr><td> 1.540952</td><td> 1.362113</td><td> 1.401243</td><td> 1.493139</td><td> 1.676425</td><td> 1.391493</td><td> 1.424900</td></tr>
	<tr><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td></tr>
	<tr><td>31.496273</td><td>24.497257</td><td>26.261684</td><td>28.248245</td><td>33.783149</td><td>24.115339</td><td>25.557649</td></tr>
	<tr><td> 1.272510</td><td> 1.235393</td><td> 1.235393</td><td> 1.410135</td><td> 1.572049</td><td> 1.320690</td><td> 1.350015</td></tr>
	<tr><td>39.574016</td><td>35.268998</td><td>36.486347</td><td>39.509897</td><td>41.155246</td><td>37.231590</td><td>38.087453</td></tr>
	<tr><td> 3.007939</td><td> 2.241205</td><td> 2.374051</td><td> 4.663489</td><td> 6.435432</td><td> 3.787345</td><td> 4.064923</td></tr>
	<tr><td>30.226854</td><td>23.094076</td><td>24.861736</td><td>30.761178</td><td>35.739955</td><td>26.756252</td><td>28.171779</td></tr>
	<tr><td>40.515905</td><td>36.927545</td><td>37.991723</td><td>37.798058</td><td>40.630987</td><td>35.058068</td><td>36.069295</td></tr>
	<tr><td>32.701671</td><td>25.894460</td><td>27.657239</td><td>26.703800</td><td>32.497732</td><td>22.541631</td><td>23.979173</td></tr>
	<tr><td>39.577046</td><td>40.332221</td><td>40.141540</td><td>39.730245</td><td>39.041081</td><td>40.158160</td><td>40.007983</td></tr>
	<tr><td>40.251749</td><td>41.047022</td><td>40.875490</td><td>40.821382</td><td>40.122520</td><td>41.163199</td><td>41.085327</td></tr>
	<tr><td>40.822220</td><td>37.668984</td><td>38.666354</td><td>40.257956</td><td>39.617865</td><td>40.732565</td><td>40.563526</td></tr>
	<tr><td>38.176534</td><td>33.216829</td><td>34.591236</td><td>37.175290</td><td>40.291997</td><td>34.309938</td><td>35.363456</td></tr>
	<tr><td>37.918643</td><td>32.827518</td><td>34.238050</td><td>35.476826</td><td>39.207078</td><td>32.234383</td><td>33.438759</td></tr>
	<tr><td>11.390519</td><td> 7.740024</td><td> 8.513032</td><td> 5.957338</td><td> 8.295242</td><td> 4.800789</td><td> 5.169442</td></tr>
	<tr><td>40.990872</td><td>40.720638</td><td>41.048543</td><td>40.896429</td><td>41.016026</td><td>39.491273</td><td>40.109904</td></tr>
	<tr><td>38.013491</td><td>32.970416</td><td>34.370317</td><td>39.909896</td><td>41.173376</td><td>37.794720</td><td>38.612022</td></tr>
	<tr><td>11.603821</td><td> 7.890855</td><td> 8.678108</td><td> 6.836963</td><td> 9.530466</td><td> 5.488349</td><td> 5.918553</td></tr>
	<tr><td>27.318736</td><td>20.125555</td><td>21.832042</td><td>23.685815</td><td>29.781590</td><td>19.626678</td><td>21.002394</td></tr>
	<tr><td>41.160139</td><td>39.112131</td><td>39.923418</td><td>36.836784</td><td>40.105471</td><td>33.878553</td><td>34.964325</td></tr>
	<tr><td>39.068635</td><td>34.517567</td><td>35.782583</td><td>37.339800</td><td>40.382268</td><td>34.520268</td><td>35.557826</td></tr>
	<tr><td>33.162120</td><td>26.427698</td><td>28.166129</td><td>30.512950</td><td>35.545092</td><td>26.483047</td><td>27.911154</td></tr>
	<tr><td> 1.272510</td><td> 1.235393</td><td> 1.235393</td><td> 1.235393</td><td> 1.311863</td><td> 1.235393</td><td> 1.235393</td></tr>
	<tr><td> 1.272510</td><td> 1.235393</td><td> 1.235393</td><td> 1.235393</td><td> 1.311863</td><td> 1.235393</td><td> 1.235393</td></tr>
	<tr><td> 1.366814</td><td> 1.235393</td><td> 1.260553</td><td> 1.235393</td><td> 1.311863</td><td> 1.235393</td><td> 1.235393</td></tr>
	<tr><td> 1.366814</td><td> 1.235393</td><td> 1.260553</td><td> 1.235393</td><td> 1.311863</td><td> 1.235393</td><td> 1.235393</td></tr>
	<tr><td>40.861849</td><td>40.989587</td><td>41.158500</td><td>41.023475</td><td>40.942835</td><td>39.826436</td><td>40.371268</td></tr>
	<tr><td> 1.366814</td><td> 1.235393</td><td> 1.260553</td><td> 1.235393</td><td> 1.311863</td><td> 1.235393</td><td> 1.235393</td></tr>
	<tr><td> 1.739728</td><td> 1.513720</td><td> 1.563236</td><td> 1.997115</td><td> 2.348299</td><td> 1.818942</td><td> 1.877094</td></tr>
	<tr><td>41.175474</td><td>39.494839</td><td>40.222560</td><td>40.577953</td><td>41.122084</td><td>38.821995</td><td>39.534319</td></tr>
	<tr><td>40.406164</td><td>41.141055</td><td>41.015448</td><td>40.923153</td><td>40.224380</td><td>41.173242</td><td>41.141224</td></tr>
	<tr><td> 2.143318</td><td> 1.814607</td><td> 1.885979</td><td> 1.689918</td><td> 1.928089</td><td> 1.560064</td><td> 1.602673</td></tr>
</tbody>
</table>
</dd>
	<dt>$design</dt>
		<dd><table>
<thead><tr><th scope=col>Intercept</th><th scope=col>trt</th></tr></thead>
<tbody>
	<tr><td>1</td><td>1</td></tr>
	<tr><td>1</td><td>1</td></tr>
	<tr><td>1</td><td>1</td></tr>
	<tr><td>1</td><td>0</td></tr>
	<tr><td>1</td><td>0</td></tr>
	<tr><td>1</td><td>0</td></tr>
	<tr><td>1</td><td>0</td></tr>
</tbody>
</table>
</dd>
	<dt>$targets</dt>
		<dd><table>
<thead><tr><th></th><th scope=col>lib.size</th></tr></thead>
<tbody>
	<tr><th scope=row>T_1</th><td> 9903374</td></tr>
	<tr><th scope=row>T_2</th><td> 5694724</td></tr>
	<tr><th scope=row>T_3</th><td> 6516297</td></tr>
	<tr><th scope=row>C_1</th><td> 7319820</td></tr>
	<tr><th scope=row>C_2</th><td>11738017</td></tr>
	<tr><th scope=row>C_3</th><td> 5328501</td></tr>
	<tr><th scope=row>C_4</th><td> 5947557</td></tr>
</tbody>
</table>
</dd>
</dl>



4、Fit a linear model on the `eset` data and design matrix， and use `eBayes` to perform the computation of statistics.<br>
在`eset`数据和设计矩阵中拟合一个线性模型，并使用eBayes来执行统计计算。


```R
fit <- lmFit(eset_voom,design)
fitE <- eBayes(fit)
```

5、 Use `topTable` function to find the top genes and filter them based on the corresponding p-values.<br>
使用topTable函数查找top基因并根据相应的p值对它们进行过滤。


```R
topAll <- topTable(fitE, n=nrow(eset), coef = 2, adjust = "BH")
topAll
DEgenes <- rownames(topAll[which(topAll$adj.P.Val<0.05),])
DEgenes
```


<table>
<thead><tr><th></th><th scope=col>logFC</th><th scope=col>AveExpr</th><th scope=col>t</th><th scope=col>P.Value</th><th scope=col>adj.P.Val</th><th scope=col>B</th></tr></thead>
<tbody>
	<tr><th scope=row>FBgn0029167</th><td>-2.1473232  </td><td> 7.8828594  </td><td>-22.057481  </td><td>3.268794e-10</td><td>2.875966e-06</td><td>13.859933   </td></tr>
	<tr><th scope=row>FBgn0035085</th><td>-2.4294695  </td><td> 5.2437331  </td><td>-19.828037  </td><td>9.885392e-10</td><td>4.768054e-06</td><td>12.308149   </td></tr>
	<tr><th scope=row>FBgn0039155</th><td>-4.3966425  </td><td> 4.8505796  </td><td>-21.646659  </td><td>3.975074e-10</td><td>2.875966e-06</td><td>11.973515   </td></tr>
	<tr><th scope=row>FBgn0001226</th><td> 1.7553191  </td><td> 6.3299240  </td><td> 14.065912  </td><td>3.340458e-08</td><td>7.614614e-05</td><td> 9.525138   </td></tr>
	<tr><th scope=row>FBgn0011260</th><td> 2.4179643  </td><td> 3.8409988  </td><td> 14.108427  </td><td>3.239964e-08</td><td>7.614614e-05</td><td> 9.203474   </td></tr>
	<tr><th scope=row>FBgn0034736</th><td>-3.2489931  </td><td> 3.3469966  </td><td>-15.562978  </td><td>1.194885e-08</td><td>4.322496e-05</td><td> 9.079040   </td></tr>
	<tr><th scope=row>FBgn0029896</th><td>-2.5112555  </td><td> 4.6935590  </td><td>-13.623414  </td><td>4.613932e-08</td><td>7.614614e-05</td><td> 8.938973   </td></tr>
	<tr><th scope=row>FBgn0000071</th><td> 2.7514504  </td><td> 4.0544440  </td><td> 13.172847  </td><td>6.474037e-08</td><td>9.367931e-05</td><td> 8.650658   </td></tr>
	<tr><th scope=row>FBgn0051092</th><td> 2.4912495  </td><td> 3.1985276  </td><td> 13.588162  </td><td>4.736111e-08</td><td>7.614614e-05</td><td> 8.632674   </td></tr>
	<tr><th scope=row>FBgn0040091</th><td>-1.4807943  </td><td> 6.2259329  </td><td>-12.566810  </td><td>1.038228e-07</td><td>1.251930e-04</td><td> 8.409547   </td></tr>
	<tr><th scope=row>FBgn0026562</th><td>-2.3074332  </td><td>11.2395357  </td><td>-12.322033  </td><td>1.263671e-07</td><td>1.406563e-04</td><td> 8.206601   </td></tr>
	<tr><th scope=row>FBgn0023479</th><td>-1.4362778  </td><td> 7.8071862  </td><td>-11.755540  </td><td>2.018355e-07</td><td>1.958230e-04</td><td> 7.730486   </td></tr>
	<tr><th scope=row>FBgn0003501</th><td> 2.1879952  </td><td> 4.0953602  </td><td> 11.748750  </td><td>2.029955e-07</td><td>1.958230e-04</td><td> 7.654270   </td></tr>
	<tr><th scope=row>FBgn0033764</th><td> 3.5672217  </td><td> 1.6028610  </td><td> 13.591134  </td><td>4.725675e-08</td><td>7.614614e-05</td><td> 7.539485   </td></tr>
	<tr><th scope=row>FBgn0003137</th><td> 0.9626954  </td><td> 8.4186859  </td><td> 11.019348  </td><td>3.822832e-07</td><td>3.048920e-04</td><td> 7.073661   </td></tr>
	<tr><th scope=row>FBgn0003748</th><td> 1.0724705  </td><td> 6.3881966  </td><td> 10.967714  </td><td>4.003419e-07</td><td>3.048920e-04</td><td> 7.053516   </td></tr>
	<tr><th scope=row>FBgn0035189</th><td> 3.0412367  </td><td> 3.2592484  </td><td> 11.097435  </td><td>3.566323e-07</td><td>3.035571e-04</td><td> 6.943769   </td></tr>
	<tr><th scope=row>FBgn0001225</th><td> 1.3268358  </td><td> 4.7947688  </td><td> 10.815053  </td><td>4.593802e-07</td><td>3.165348e-04</td><td> 6.941237   </td></tr>
	<tr><th scope=row>FBgn0031191+FBgn0027279</th><td>-0.9535197  </td><td> 8.0737077  </td><td>-10.443465  </td><td>6.465812e-07</td><td>4.070329e-04</td><td> 6.536008   </td></tr>
	<tr><th scope=row>FBgn0001258</th><td> 0.9694595  </td><td> 6.8767036  </td><td> 10.412625  </td><td>6.654907e-07</td><td>4.070329e-04</td><td> 6.522279   </td></tr>
	<tr><th scope=row>FBgn0039479</th><td> 1.6856529  </td><td> 3.5402736  </td><td> 10.397306  </td><td>6.751064e-07</td><td>4.070329e-04</td><td> 6.478714   </td></tr>
	<tr><th scope=row>FBgn0034897</th><td>-1.7227324  </td><td> 5.9723986  </td><td>-10.221359  </td><td>7.970469e-07</td><td>4.435873e-04</td><td> 6.388645   </td></tr>
	<tr><th scope=row>FBgn0034434</th><td>-3.6109506  </td><td> 2.5597985  </td><td>-11.165049  </td><td>3.359297e-07</td><td>3.035571e-04</td><td> 6.212774   </td></tr>
	<tr><th scope=row>FBgn0001224</th><td> 1.5299878  </td><td> 5.2645851  </td><td> 10.007762  </td><td>9.781946e-07</td><td>5.055170e-04</td><td> 6.189658   </td></tr>
	<tr><th scope=row>FBgn0033913</th><td>-1.1113301  </td><td> 8.5071266  </td><td>-10.062709  </td><td>9.276766e-07</td><td>4.971659e-04</td><td> 6.162981   </td></tr>
	<tr><th scope=row>FBgn0037290</th><td> 3.6225526  </td><td> 1.7009716  </td><td> 10.832129  </td><td>4.523291e-07</td><td>3.165348e-04</td><td> 6.105169   </td></tr>
	<tr><th scope=row>FBgn0034751</th><td>-1.0174301  </td><td>10.1675244  </td><td> -9.943721  </td><td>1.040873e-06</td><td>5.193596e-04</td><td> 6.051593   </td></tr>
	<tr><th scope=row>FBgn0260011</th><td> 2.2894646  </td><td> 2.2430213  </td><td> 10.254590  </td><td>7.722992e-07</td><td>4.435873e-04</td><td> 6.011483   </td></tr>
	<tr><th scope=row>FBgn0037153</th><td> 1.9007262  </td><td> 2.8750861  </td><td>  9.852796  </td><td>1.137466e-06</td><td>5.486378e-04</td><td> 5.883048   </td></tr>
	<tr><th scope=row>FBgn0024288</th><td>-4.7192479  </td><td> 0.8696533  </td><td>-12.580464  </td><td>1.027012e-07</td><td>1.251930e-04</td><td> 5.791253   </td></tr>
	<tr><th scope=row>...</th><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td></tr>
	<tr><th scope=row>FBgn0035121</th><td> 0.0002311982</td><td>9.001532     </td><td> 0.002363025 </td><td>0.9981584    </td><td>0.9989661    </td><td>-7.244785    </td></tr>
	<tr><th scope=row>FBgn0033688</th><td> 0.0051489254</td><td>8.874872     </td><td> 0.052110494 </td><td>0.9594081    </td><td>0.9867547    </td><td>-7.244833    </td></tr>
	<tr><th scope=row>FBgn0011217</th><td>-0.0169947318</td><td>8.637299     </td><td>-0.086366332 </td><td>0.9327823    </td><td>0.9867547    </td><td>-7.244914    </td></tr>
	<tr><th scope=row>FBgn0011837+FBgn0053963</th><td> 0.0080708276</td><td>7.818126     </td><td> 0.057190834 </td><td>0.9554552    </td><td>0.9867547    </td><td>-7.245085    </td></tr>
	<tr><th scope=row>FBgn0010341</th><td>-0.0064128966</td><td>8.775797     </td><td>-0.059842214 </td><td>0.9533927    </td><td>0.9867547    </td><td>-7.245540    </td></tr>
	<tr><th scope=row>FBgn0013984</th><td>-0.0190160579</td><td>8.484710     </td><td>-0.093084377 </td><td>0.9275697    </td><td>0.9867547    </td><td>-7.245563    </td></tr>
	<tr><th scope=row>FBgn0028982</th><td>-0.0080417611</td><td>8.140941     </td><td>-0.096219005 </td><td>0.9251386    </td><td>0.9867547    </td><td>-7.245908    </td></tr>
	<tr><th scope=row>FBgn0022959</th><td>-0.0098744118</td><td>8.294319     </td><td>-0.094825264 </td><td>0.9262194    </td><td>0.9867547    </td><td>-7.246254    </td></tr>
	<tr><th scope=row>FBgn0005536</th><td>-0.0097699694</td><td>8.469107     </td><td>-0.086578333 </td><td>0.9326178    </td><td>0.9867547    </td><td>-7.246294    </td></tr>
	<tr><th scope=row>FBgn0051716</th><td>-0.0041347383</td><td>7.855832     </td><td>-0.043609608 </td><td>0.9660249    </td><td>0.9867547    </td><td>-7.246569    </td></tr>
	<tr><th scope=row>FBgn0002781</th><td>-0.0066747665</td><td>8.475092     </td><td>-0.078481548 </td><td>0.9389045    </td><td>0.9867547    </td><td>-7.246979    </td></tr>
	<tr><th scope=row>FBgn0037912+FBgn0037913</th><td>-0.0108428346</td><td>8.133232     </td><td>-0.083402905 </td><td>0.9350828    </td><td>0.9867547    </td><td>-7.247139    </td></tr>
	<tr><th scope=row>FBgn0020443</th><td> 0.0019068986</td><td>8.781464     </td><td> 0.020886651 </td><td>0.9837234    </td><td>0.9923307    </td><td>-7.247173    </td></tr>
	<tr><th scope=row>FBgn0020497</th><td>-0.0025093031</td><td>8.745295     </td><td>-0.019945568 </td><td>0.9844567    </td><td>0.9924816    </td><td>-7.247605    </td></tr>
	<tr><th scope=row>FBgn0030504</th><td> 0.0028675532</td><td>7.890754     </td><td> 0.030611244 </td><td>0.9761474    </td><td>0.9883048    </td><td>-7.247778    </td></tr>
	<tr><th scope=row>FBgn0030740</th><td>-0.0117377188</td><td>7.989674     </td><td>-0.058791020 </td><td>0.9542104    </td><td>0.9867547    </td><td>-7.247793    </td></tr>
	<tr><th scope=row>FBgn0013981</th><td>-0.0081550532</td><td>8.274796     </td><td>-0.075827751 </td><td>0.9409660    </td><td>0.9867547    </td><td>-7.248036    </td></tr>
	<tr><th scope=row>FBgn0003717</th><td>-0.0115861823</td><td>8.238971     </td><td>-0.074016632 </td><td>0.9423731    </td><td>0.9867547    </td><td>-7.248185    </td></tr>
	<tr><th scope=row>FBgn0026577</th><td>-0.0113688428</td><td>8.337305     </td><td>-0.070315210 </td><td>0.9452496    </td><td>0.9867547    </td><td>-7.248320    </td></tr>
	<tr><th scope=row>FBgn0036448</th><td> 0.0063642000</td><td>7.995062     </td><td> 0.050530450 </td><td>0.9606377    </td><td>0.9867547    </td><td>-7.248354    </td></tr>
	<tr><th scope=row>FBgn0002524</th><td>-0.0026342480</td><td>8.655167     </td><td>-0.024532292 </td><td>0.9808830    </td><td>0.9907425    </td><td>-7.248447    </td></tr>
	<tr><th scope=row>FBgn0010609</th><td> 0.0032335056</td><td>8.565167     </td><td> 0.032835762 </td><td>0.9744147    </td><td>0.9878005    </td><td>-7.249011    </td></tr>
	<tr><th scope=row>FBgn0031030</th><td> 0.0004140565</td><td>8.627531     </td><td> 0.003205618 </td><td>0.9975017    </td><td>0.9986547    </td><td>-7.249039    </td></tr>
	<tr><th scope=row>FBgn0003392</th><td>-0.0020112331</td><td>7.991076     </td><td>-0.019567173 </td><td>0.9847515    </td><td>0.9927096    </td><td>-7.249486    </td></tr>
	<tr><th scope=row>FBgn0005674</th><td>-0.0029027731</td><td>8.506550     </td><td>-0.027173549 </td><td>0.9788253    </td><td>0.9898380    </td><td>-7.249676    </td></tr>
	<tr><th scope=row>FBgn0035872</th><td> 0.0011257442</td><td>8.030984     </td><td> 0.011770947 </td><td>0.9908266    </td><td>0.9959720    </td><td>-7.250053    </td></tr>
	<tr><th scope=row>FBgn0004924</th><td>-0.0008096367</td><td>8.143069     </td><td>-0.008304370 </td><td>0.9935281    </td><td>0.9967430    </td><td>-7.250909    </td></tr>
	<tr><th scope=row>FBgn0040335</th><td> 0.0002434155</td><td>8.132772     </td><td> 0.001684001 </td><td>0.9986876    </td><td>0.9991019    </td><td>-7.250917    </td></tr>
	<tr><th scope=row>FBgn0027338</th><td> 0.0006004041</td><td>8.156461     </td><td> 0.005700789 </td><td>0.9955571    </td><td>0.9975564    </td><td>-7.250964    </td></tr>
	<tr><th scope=row>FBgn0031408</th><td>-0.0008151292</td><td>8.211378     </td><td>-0.009702040 </td><td>0.9924389    </td><td>0.9967337    </td><td>-7.251074    </td></tr>
</tbody>
</table>




<ol class=list-inline>
	<li>'FBgn0029167'</li>
	<li>'FBgn0035085'</li>
	<li>'FBgn0039155'</li>
	<li>'FBgn0001226'</li>
	<li>'FBgn0011260'</li>
	<li>'FBgn0034736'</li>
	<li>'FBgn0029896'</li>
	<li>'FBgn0000071'</li>
	<li>'FBgn0051092'</li>
	<li>'FBgn0040091'</li>
	<li>'FBgn0026562'</li>
	<li>'FBgn0023479'</li>
	<li>'FBgn0003501'</li>
	<li>'FBgn0033764'</li>
	<li>'FBgn0003137'</li>
	<li>'FBgn0003748'</li>
	<li>'FBgn0035189'</li>
	<li>'FBgn0001225'</li>
	<li>'FBgn0031191+FBgn0027279'</li>
	<li>'FBgn0001258'</li>
	<li>'FBgn0039479'</li>
	<li>'FBgn0034897'</li>
	<li>'FBgn0034434'</li>
	<li>'FBgn0001224'</li>
	<li>'FBgn0033913'</li>
	<li>'FBgn0037290'</li>
	<li>'FBgn0034751'</li>
	<li>'FBgn0260011'</li>
	<li>'FBgn0037153'</li>
	<li>'FBgn0024288'</li>
	<li>'FBgn0051642'</li>
	<li>'FBgn0034438'</li>
	<li>'FBgn0038149'</li>
	<li>'FBgn0039109'</li>
	<li>'FBgn0063649'</li>
	<li>'FBgn0038832'</li>
	<li>'FBgn0039419'</li>
	<li>'FBgn0050147'</li>
	<li>'FBgn0016715'</li>
	<li>'FBgn0261362'</li>
	<li>'FBgn0035968'</li>
	<li>'FBgn0261238'</li>
	<li>'FBgn0040271'</li>
	<li>'FBgn0032405'</li>
	<li>'FBgn0034405'</li>
	<li>'FBgn0003076'</li>
	<li>'FBgn0261552'</li>
	<li>'FBgn0004108'</li>
	<li>'FBgn0039113'</li>
	<li>'FBgn0051038'</li>
	<li>'FBgn0036299'</li>
	<li>'FBgn0025111+FBgn0003360'</li>
	<li>'FBgn0001124'</li>
	<li>'FBgn0038528'</li>
	<li>'FBgn0027515'</li>
	<li>'FBgn0031150'</li>
	<li>'FBgn0031327'</li>
	<li>'FBgn0020248'</li>
	<li>'FBgn0052407'</li>
	<li>'FBgn0015568'</li>
	<li>'FBgn0030805'</li>
	<li>'FBgn0024315'</li>
	<li>'FBgn0051195'</li>
	<li>'FBgn0038293'</li>
	<li>'FBgn0261284'</li>
	<li>'FBgn0033367'</li>
	<li>'FBgn0036968'</li>
	<li>'FBgn0040827'</li>
	<li>'FBgn0038341'</li>
	<li>'FBgn0038198'</li>
	<li>'FBgn0000116'</li>
	<li>'FBgn0038805'</li>
	<li>'FBgn0035147'</li>
	<li>'FBgn0036684'</li>
	<li>'FBgn0003502'</li>
	<li>'FBgn0034010'</li>
	<li>'FBgn0030598'</li>
	<li>'FBgn0035765'</li>
	<li>'FBgn0029801'</li>
	<li>'FBgn0033724'</li>
	<li>'FBgn0086910'</li>
	<li>'FBgn0039827'</li>
	<li>'FBgn0051363'</li>
	<li>'FBgn0036007'</li>
	<li>'FBgn0035403'</li>
	<li>'FBgn0050463'</li>
	<li>'FBgn0033095'</li>
	<li>'FBgn0024984'</li>
	<li>'FBgn0003317'</li>
	<li>'FBgn0031805'</li>
	<li>'FBgn0051523'</li>
	<li>'FBgn0035811'</li>
	<li>'FBgn0039593'</li>
	<li>'FBgn0030037'</li>
	<li>'FBgn0034885'</li>
	<li>'FBgn0029002'</li>
	<li>'FBgn0035344'</li>
	<li>'FBgn0027580'</li>
	<li>'FBgn0001137'</li>
	<li>'FBgn0031689'</li>
	<li>'FBgn0040398'</li>
	<li>'FBgn0001228+FBgn0001223'</li>
	<li>'FBgn0261015+FBgn0053191+FBgn0046874'</li>
	<li>'FBgn0028988+FBgn0033113'</li>
	<li>'FBgn0260933'</li>
	<li>'FBgn0038436'</li>
	<li>'FBgn0025678'</li>
	<li>'FBgn0004507'</li>
	<li>'FBgn0032004'</li>
	<li>'FBgn0086355'</li>
	<li>'FBgn0037635'</li>
	<li>'FBgn0003074'</li>
	<li>'FBgn0034391'</li>
	<li>'FBgn0035266'</li>
	<li>'FBgn0031307'</li>
	<li>'FBgn0037680'</li>
	<li>'FBgn0001091'</li>
	<li>'FBgn0039640'</li>
	<li>'FBgn0051431'</li>
	<li>'FBgn0053307'</li>
	<li>'FBgn0033775'</li>
	<li>'FBgn0023549'</li>
	<li>'FBgn0034389'</li>
	<li>'FBgn0037754'</li>
	<li>'FBgn0032105'</li>
	<li>'FBgn0036985'</li>
	<li>'FBgn0000299'</li>
	<li>'FBgn0040099'</li>
	<li>'FBgn0039257'</li>
	<li>'FBgn0032089'</li>
	<li>'FBgn0034067'</li>
	<li>'FBgn0010591'</li>
	<li>'FBgn0039098'</li>
	<li>'FBgn0015541'</li>
	<li>'FBgn0015777'</li>
	<li>'FBgn0030763'</li>
	<li>'FBgn0032088'</li>
	<li>'FBgn0031449'</li>
	<li>'FBgn0011674'</li>
	<li>'FBgn0085359'</li>
	<li>'FBgn0028327'</li>
	<li>'FBgn0034718'</li>
	<li>'FBgn0027091'</li>
	<li>'FBgn0031516'</li>
	<li>'FBgn0052625'</li>
	<li>'FBgn0063498'</li>
	<li>'FBgn0026415'</li>
	<li>'FBgn0000146'</li>
	<li>'FBgn0000567'</li>
	<li>'FBgn0031374'</li>
	<li>'FBgn0014869'</li>
	<li>'FBgn0033782'</li>
	<li>'FBgn0259715'</li>
	<li>'FBgn0040752'</li>
	<li>'FBgn0039464'</li>
	<li>'FBgn0029090'</li>
	<li>'FBgn0038545'</li>
	<li>'FBgn0028990'</li>
	<li>'FBgn0051116'</li>
	<li>'FBgn0032820'</li>
	<li>'FBgn0004369'</li>
	<li>'FBgn0039485'</li>
	<li>'FBgn0033065'</li>
	<li>'FBgn0052677'</li>
	<li>'FBgn0020240'</li>
	<li>'FBgn0037683'</li>
	<li>'FBgn0031268'</li>
	<li>'FBgn0024491'</li>
	<li>'FBgn0261560'</li>
	<li>'FBgn0085446'</li>
	<li>'FBgn0044047'</li>
	<li>'FBgn0033368'</li>
	<li>'FBgn0030145'</li>
	<li>'FBgn0087007'</li>
	<li>'FBgn0030041'</li>
	<li>'FBgn0034793'</li>
	<li>'FBgn0040388'</li>
	<li>'FBgn0030452'</li>
	<li>'FBgn0261278+FBgn0002652'</li>
	<li>'FBgn0010387'</li>
	<li>'FBgn0250906'</li>
	<li>'FBgn0027836'</li>
	<li>'FBgn0011710'</li>
	<li>'FBgn0027949'</li>
	<li>'FBgn0039525'</li>
	<li>'FBgn0038237'</li>
	<li>'FBgn0035761'</li>
	<li>'FBgn0031703'</li>
	<li>'FBgn0052021'</li>
	<li>'FBgn0000406'</li>
	<li>'FBgn0032420'</li>
	<li>'FBgn0033635'</li>
	<li>'FBgn0031313'</li>
	<li>'FBgn0032436'</li>
	<li>'FBgn0011722'</li>
	<li>'FBgn0010223'</li>
	<li>'FBgn0053318'</li>
	<li>'FBgn0250904'</li>
	<li>'FBgn0035232'</li>
	<li>'FBgn0038012'</li>
	<li>'FBgn0039637'</li>
	<li>'FBgn0261563'</li>
	<li>'FBgn0038877'</li>
	<li>'FBgn0086906'</li>
	<li>'FBgn0031322'</li>
	<li>'FBgn0053516'</li>
	<li>'FBgn0038720'</li>
	<li>'FBgn0050035'</li>
	<li>'FBgn0000579'</li>
	<li>'FBgn0032036'</li>
	<li>'FBgn0041629'</li>
	<li>'FBgn0035937'</li>
	<li>'FBgn0032638'</li>
	<li>'FBgn0261444'</li>
	<li>'FBgn0031688'</li>
	<li>'FBgn0063667'</li>
	<li>'FBgn0030964'</li>
	<li>'FBgn0034354'</li>
	<li>'FBgn0032775'</li>
	<li>'FBgn0085402'</li>
	<li>'FBgn0039735'</li>
	<li>'FBgn0038865'</li>
	<li>'FBgn0037739'</li>
	<li>'FBgn0037239'</li>
	<li>'FBgn0010470'</li>
	<li>'FBgn0031912'</li>
	<li>'FBgn0031461'</li>
	<li>'FBgn0050185'</li>
	<li>'FBgn0015371'</li>
	<li>'FBgn0045852'</li>
	<li>'FBgn0004893'</li>
	<li>'FBgn0259111'</li>
	<li>'FBgn0039538'</li>
	<li>'FBgn0033426'</li>
	<li>'FBgn0034638'</li>
	<li>'FBgn0037607'</li>
	<li>'FBgn0027596'</li>
	<li>'FBgn0036688'</li>
	<li>'FBgn0030318'</li>
	<li>'FBgn0033188'</li>
	<li>'FBgn0032451'</li>
	<li>'FBgn0053508'</li>
	<li>'FBgn0011823'</li>
	<li>'FBgn0038815'</li>
	<li>'FBgn0001281+FBgn0001280'</li>
	<li>'FBgn0051856'</li>
	<li>'FBgn0016075'</li>
	<li>'FBgn0027843'</li>
	<li>'FBgn0030529'</li>
	<li>'FBgn0033760'</li>
	<li>'FBgn0013811'</li>
	<li>'FBgn0040308'</li>
	<li>'FBgn0001122'</li>
	<li>'FBgn0032421'</li>
	<li>'FBgn0259224'</li>
	<li>'FBgn0031888'</li>
	<li>'FBgn0000405'</li>
	<li>'FBgn0035539'</li>
	<li>'FBgn0023170'</li>
	<li>'FBgn0250757'</li>
	<li>'FBgn0037896'</li>
	<li>'FBgn0038351'</li>
	<li>'FBgn0036641'</li>
	<li>'FBgn0050069'</li>
	<li>'FBgn0032598'</li>
	<li>'FBgn0052813'</li>
	<li>'FBgn0027657'</li>
	<li>'FBgn0259749'</li>
	<li>'FBgn0051555'</li>
	<li>'FBgn0035150'</li>
	<li>'FBgn0031816'</li>
	<li>'FBgn0036576'</li>
	<li>'FBgn0086365'</li>
	<li>'FBgn0033079'</li>
	<li>'FBgn0001197'</li>
	<li>'FBgn0027597'</li>
	<li>'FBgn0002527'</li>
	<li>'FBgn0031117'</li>
	<li>'FBgn0037872'</li>
	<li>'FBgn0029664'</li>
	<li>'FBgn0031538'</li>
	<li>'FBgn0046258'</li>
	<li>'FBgn0013763'</li>
	<li>'FBgn0038660'</li>
	<li>'FBgn0041604'</li>
	<li>'FBgn0260966'</li>
	<li>'FBgn0030343'</li>
	<li>'FBgn0032029'</li>
	<li>'FBgn0037901'</li>
	<li>'FBgn0031489'</li>
	<li>'FBgn0002868'</li>
	<li>'FBgn0035091'</li>
	<li>'FBgn0037143'</li>
	<li>'FBgn0021795'</li>
	<li>'FBgn0036663'</li>
	<li>'FBgn0037850'</li>
	<li>'FBgn0026189'</li>
	<li>'FBgn0013765'</li>
	<li>'FBgn0004396'</li>
	<li>'FBgn0033268'</li>
	<li>'FBgn0032635'</li>
	<li>'FBgn0036896'</li>
	<li>'FBgn0259236'</li>
	<li>'FBgn0034564'</li>
	<li>'FBgn0005593'</li>
	<li>'FBgn0031245'</li>
	<li>'FBgn0050000'</li>
	<li>'FBgn0259935'</li>
	<li>'FBgn0031769'</li>
	<li>'FBgn0052681'</li>
	<li>'FBgn0037223'</li>
	<li>'FBgn0037191'</li>
	<li>'FBgn0022338'</li>
	<li>'FBgn0035879'</li>
	<li>'FBgn0032078'</li>
	<li>'FBgn0003366'</li>
	<li>'FBgn0033777'</li>
	<li>'FBgn0015522'</li>
	<li>'FBgn0000928'</li>
	<li>'FBgn0051614'</li>
	<li>'FBgn0033205'</li>
	<li>'FBgn0002578'</li>
	<li>'FBgn0001114'</li>
	<li>'FBgn0028514'</li>
	<li>'FBgn0030349'</li>
	<li>'FBgn0031055'</li>
	<li>'FBgn0030966'</li>
	<li>'FBgn0051663'</li>
	<li>'FBgn0000043'</li>
	<li>'FBgn0010225'</li>
	<li>'FBgn0020376'</li>
	<li>'FBgn0034694'</li>
	<li>'FBgn0033926'</li>
	<li>'FBgn0085339+FBgn0000392'</li>
	<li>'FBgn0051337'</li>
	<li>'FBgn0011606'</li>
	<li>'FBgn0085403'</li>
	<li>'FBgn0025628'</li>
	<li>'FBgn0052700'</li>
	<li>'FBgn0259740'</li>
	<li>'FBgn0010786'</li>
	<li>'FBgn0014029'</li>
	<li>'FBgn0032785'</li>
	<li>'FBgn0036862'</li>
	<li>'FBgn0063492'</li>
	<li>'FBgn0035076'</li>
	<li>'FBgn0030322'</li>
	<li>'FBgn0023091'</li>
	<li>'FBgn0051361'</li>
	<li>'FBgn0085407'</li>
	<li>'FBgn0085434'</li>
	<li>'FBgn0011259'</li>
	<li>'FBgn0039937'</li>
	<li>'FBgn0033733'</li>
	<li>'FBgn0050324'</li>
	<li>'FBgn0040813'</li>
	<li>'FBgn0033478'</li>
	<li>'FBgn0028473'</li>
	<li>'FBgn0005660'</li>
	<li>'FBgn0039052'</li>
	<li>'FBgn0041342'</li>
	<li>'FBgn0015576'</li>
	<li>'FBgn0260768'</li>
	<li>'FBgn0034898'</li>
	<li>'FBgn0028939'</li>
	<li>'FBgn0030954'</li>
	<li>'FBgn0010768'</li>
	<li>'FBgn0026585'</li>
	<li>'FBgn0053556'</li>
	<li>'FBgn0031996'</li>
	<li>'FBgn0250837'</li>
	<li>'FBgn0033288'</li>
	<li>'FBgn0001257'</li>
	<li>'FBgn0010228'</li>
	<li>'FBgn0037672'</li>
	<li>'FBgn0038381'</li>
	<li>'FBgn0259714'</li>
	<li>'FBgn0038476'</li>
	<li>'FBgn0039788'</li>
	<li>'FBgn0037646'</li>
	<li>'FBgn0259217'</li>
	<li>'FBgn0034911'</li>
	<li>'FBgn0030237'</li>
	<li>'FBgn0026376'</li>
	<li>'FBgn0085354'</li>
	<li>'FBgn0250829'</li>
	<li>'FBgn0038049'</li>
	<li>'FBgn0052135'</li>
	<li>'FBgn0051776'</li>
	<li>'FBgn0001332'</li>
	<li>'FBgn0031547'</li>
	<li>'FBgn0053193'</li>
	<li>'FBgn0032770'</li>
	<li>'FBgn0020639'</li>
	<li>'FBgn0053681+FBgn0027348'</li>
	<li>'FBgn0032881'</li>
	<li>'FBgn0259734'</li>
	<li>'FBgn0039000'</li>
	<li>'FBgn0004512'</li>
	<li>'FBgn0031070'</li>
	<li>'FBgn0005640'</li>
	<li>'FBgn0086358'</li>
	<li>'FBgn0028433'</li>
	<li>'FBgn0037717'</li>
	<li>'FBgn0037468'</li>
	<li>'FBgn0036290'</li>
	<li>'FBgn0037844'</li>
	<li>'FBgn0003886'</li>
	<li>'FBgn0036373'</li>
	<li>'FBgn0033677'</li>
	<li>'FBgn0027945'</li>
	<li>'FBgn0052521'</li>
	<li>'FBgn0034249'</li>
	<li>'FBgn0035390'</li>
	<li>'FBgn0031148'</li>
	<li>'FBgn0036732'</li>
	<li>'FBgn0025681'</li>
	<li>'FBgn0024891'</li>
	<li>'FBgn0001206'</li>
	<li>'FBgn0029861'</li>
	<li>'FBgn0033502+FBgn0033504'</li>
	<li>'FBgn0020257'</li>
	<li>'FBgn0259246'</li>
	<li>'FBgn0003041'</li>
	<li>'FBgn0261109'</li>
	<li>'FBgn0038347'</li>
	<li>'FBgn0028542'</li>
	<li>'FBgn0034400'</li>
	<li>'FBgn0034013'</li>
	<li>'FBgn0031738'</li>
	<li>'FBgn0053995+FBgn0031674'</li>
	<li>'FBgn0026179'</li>
	<li>'FBgn0035211'</li>
	<li>'FBgn0033357'</li>
	<li>'FBgn0013548'</li>
	<li>'FBgn0002778'</li>
	<li>'FBgn0034390'</li>
	<li>'FBgn0028490'</li>
	<li>'FBgn0033539'</li>
	<li>'FBgn0086450'</li>
	<li>'FBgn0039644'</li>
	<li>'FBgn0035500'</li>
	<li>'FBgn0053138'</li>
	<li>'FBgn0053558'</li>
	<li>'FBgn0029095'</li>
</ol>



The `limma` analysis for a data has already been described in Chapter 5, Analyzing Microarray Data with R. In this part, we used the `count data` in place of `normalized expression values`. Another difference is the use of the `voom transformation` that estimates the mean-variance relationship for the log counts, which robustly generates a precision weight for each individual normalized observation. The plot shown in the following plots shows a histogram of the tags with different ranges of p-values and a volcano plot for the `limma` analysis.<br>
对数据的“limma”分析已经在第5章用R分析微阵列数据中进行了描述。在这一部分中，我们使用了“计数数据”来代替“归一化表达式值”。另一个区别是使用'voom变换'来估计对数计数的均值-方差关系，这为每个单独的归一化观测值生成了一个精确的权重。下面的图显示了具有不同p值范围的标签直方图和“limma”分析的火山图。


```R
hist(topAll$adj.P.Val, xlab="adjusted P values",main="(A) Histogram for P values")
clr <- rep("black",nrow(topAll)) # creates a vector for color
clr[which(topAll$adj.P.Val<0.05)] <- "red"# sets color for DE to red
plot(x = topAll$logFC, y = -log10(topAll$adj.P.Val), col = clr, xlab = "log fold change", 
     ylab = "-log P",main = "(B) Volcano Plot") # Do a volcano plot
abline(h=-log10(0.05), col="blue")# Draw a horizontal line marking a P value threshold for 0.05

# The following plots show the distribution of p-values in the analysis results and the number of tags with low p-values.
```


[output_1](https://github.com/Chengshu21/Chapter-8-Analyzing-NGS-Data/blob/master/MD/pic/output_The%20differential%20analysis%20of%20NGS%20data%20using%20limma.png)



[output_2](https://github.com/Chengshu21/Chapter-8-Analyzing-NGS-Data/blob/master/MD/pic/output_1The%20differential%20analysis%20of%20NGS%20data%20using%20limma.png)


The results show that we have 445 tags that show DE tags, that is, tags that were significantly more expressive between the two conditions. We can check how many of these show a higher fold change.<br>


```R
rownames(topAll[which(topAll$adj.P.Val<0.05&abs(topAll$logFC) > 2),])
# We can see that 59 tags show a change of more than two folds (positive or negative).
```


<ol class=list-inline>
	<li>'FBgn0029167'</li>
	<li>'FBgn0035085'</li>
	<li>'FBgn0039155'</li>
	<li>'FBgn0011260'</li>
	<li>'FBgn0034736'</li>
	<li>'FBgn0029896'</li>
	<li>'FBgn0000071'</li>
	<li>'FBgn0051092'</li>
	<li>'FBgn0026562'</li>
	<li>'FBgn0003501'</li>
	<li>'FBgn0033764'</li>
	<li>'FBgn0035189'</li>
	<li>'FBgn0034434'</li>
	<li>'FBgn0037290'</li>
	<li>'FBgn0260011'</li>
	<li>'FBgn0024288'</li>
	<li>'FBgn0051642'</li>
	<li>'FBgn0034438'</li>
	<li>'FBgn0038832'</li>
	<li>'FBgn0032405'</li>
	<li>'FBgn0020248'</li>
	<li>'FBgn0052407'</li>
	<li>'FBgn0261284'</li>
	<li>'FBgn0040827'</li>
	<li>'FBgn0038198'</li>
	<li>'FBgn0030598'</li>
	<li>'FBgn0039827'</li>
	<li>'FBgn0050463'</li>
	<li>'FBgn0039593'</li>
	<li>'FBgn0037754'</li>
	<li>'FBgn0030763'</li>
	<li>'FBgn0085359'</li>
	<li>'FBgn0033065'</li>
	<li>'FBgn0030041'</li>
	<li>'FBgn0010387'</li>
	<li>'FBgn0038237'</li>
	<li>'FBgn0032436'</li>
	<li>'FBgn0053318'</li>
	<li>'FBgn0038012'</li>
	<li>'FBgn0063667'</li>
	<li>'FBgn0030964'</li>
	<li>'FBgn0033760'</li>
	<li>'FBgn0051555'</li>
	<li>'FBgn0046258'</li>
	<li>'FBgn0037143'</li>
	<li>'FBgn0259236'</li>
	<li>'FBgn0037223'</li>
	<li>'FBgn0037191'</li>
	<li>'FBgn0002578'</li>
	<li>'FBgn0051663'</li>
	<li>'FBgn0052700'</li>
	<li>'FBgn0039937'</li>
	<li>'FBgn0033733'</li>
	<li>'FBgn0050324'</li>
	<li>'FBgn0034898'</li>
	<li>'FBgn0028939'</li>
	<li>'FBgn0051776'</li>
	<li>'FBgn0032770'</li>
	<li>'FBgn0020639'</li>
</ol>




```R

```
