
# Fold changes in microarray data
 
Fold change refers to the ratio of final value to initial value. In terms of gene expression, it can be defined as the ratio of the final quantification of mRNA to the initial content. The initial and final stages can be the time points or treatment and control conditions. It represents the change rather than an ambiguous absolute quantity. It has been suggested that while extracting DE genes from a dataset, fold changes can serve as more reproducible identifiers. This recipe will explain the use of fold changes for such purposes.<br>
倍数变化是指最终值与初始值的比率。就基因表达而言,它可以定义为mRNA的最终定量与初始含量的比率。初始阶段和最终阶段可以是时间点或治疗和控制条件。它代表了变化而不是一个模棱两可的绝对数量。从数据集中提取DE基因的同时,认为倍数变化可以作为更可重复的标识符。这个部分将解释用于这种目的的倍数变化的使用。<br>

Required:<br>
1、The leukemia dataset.<br>
2、Directly use the results from the Finding the differentially expressed genes.<br>
3、The results only from the ALL type of leukemia.<br>
Fold changes are a part of tables that are derived from `limma`. Some other interesting operations on (or about) fold change are explained.

1、Use the result from the `limma` analysis to get the fold changes. The table generated has a column for the fold change associated with the probes (an example table has been shown here). Refer to the Working with the data of multiple classes part in this chapter to create the `tested2` object. Take a look.



```R
library(leukemiasEset)
data(leukemiasEset)
pheno <- pData(leukemiasEset)
mydata <- leukemiasEset[, sampleNames(leukemiasEset)[c(1:3, 13:15, 25:27, 49:51)]]
design <- model.matrix(~0 + factor(pData(mydata)$LeukemiaType))
colnames(design) <- unique(as.character(pData(mydata)$LeukemiaType))
design
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




```R
library(limma)
fit <- lmFit(mydata, design)
contrast.matrix <- makeContrasts(NoL- ALL, NoL- AML, NoL- CLL, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
tested2 <- topTable(fit2,adjust="fdr",sort.by="B",number=Inf, coef=1)
DE2 <- tested2[tested2$adj.P.Val < 0.01,]
```


```R
head(DE2)
```


<table>
<thead><tr><th></th><th scope=col>logFC</th><th scope=col>AveExpr</th><th scope=col>t</th><th scope=col>P.Value</th><th scope=col>adj.P.Val</th><th scope=col>B</th></tr></thead>
<tbody>
	<tr><th scope=row>ENSG00000152078</th><td> 4.510507   </td><td>4.856523    </td><td> 28.13988   </td><td>4.463747e-11</td><td>9.004270e-07</td><td>14.01472    </td></tr>
	<tr><th scope=row>ENSG00000117519</th><td>-4.185175   </td><td>4.791585    </td><td>-22.73888   </td><td>3.878292e-10</td><td>3.911645e-06</td><td>12.69738    </td></tr>
	<tr><th scope=row>ENSG00000145850</th><td> 4.142236   </td><td>4.507655    </td><td> 17.38636   </td><td>5.759942e-09</td><td>2.925048e-05</td><td>10.72782    </td></tr>
	<tr><th scope=row>ENSG00000170180</th><td> 5.681327   </td><td>5.734169    </td><td> 17.37423   </td><td>5.800214e-09</td><td>2.925048e-05</td><td>10.72231    </td></tr>
	<tr><th scope=row>ENSG00000087586</th><td> 3.952183   </td><td>5.720789    </td><td> 16.45393   </td><td>9.977396e-09</td><td>3.111188e-05</td><td>10.28705    </td></tr>
	<tr><th scope=row>ENSG00000047597</th><td> 5.362419   </td><td>5.108415    </td><td> 16.32474   </td><td>1.079114e-08</td><td>3.111188e-05</td><td>10.22315    </td></tr>
</tbody>
</table>



2、Extract the relevant columns into a separate data frame for the top 10,000 probes.


```R
fit <- eBayes(fit)
myTable <- topTable(fit, number=10000)
logratio <- tested2$logFC
```

3、The fold change to log values and vice versa by the `gtools` library.


```R
library(gtools)
```


```R
LR <- foldchange2logratio(1, base=2)
```


```R
FC <- logratio2foldchange(logratio, base=2)
```


```R
??foldchange2logratio
```

    starting httpd help server ... done
    



4、Now, visualize the log fold change and p-value relations in a volcano plot.


```R
plot(tested2$logFC, -log10(tested2$P.Value),xlim=c(-10, 10), 
     ylim=c(0, 15), xlab="log2 fold change", ylab="-log10 p-value")
```


[plot](https://github.com/Chengshu21/Chapter-5-Analyzing-Microarray-Data-with-R/blob/master/md/pic/output-Fold%20changes%20in%20microarray%20data.png)


5、While selecting significant genes from the `limma` generated table, use the log fold change column as another criterion to select DE genes in combination with the p-values.


```R
myTable[tested2$P.Val < 0.05&tested2$logFC > 1.5,]
```


<table>
<thead><tr><th></th><th scope=col>ALL</th><th scope=col>AML</th><th scope=col>CLL</th><th scope=col>NoL</th><th scope=col>AveExpr</th><th scope=col>F</th><th scope=col>P.Value</th><th scope=col>adj.P.Val</th></tr></thead>
<tbody>
	<tr><th scope=row>ENSG00000177954</th><td>14.48894    </td><td>14.39784    </td><td>14.61716    </td><td>14.27750    </td><td>14.44536    </td><td>26127.989   </td><td>1.456821e-20</td><td>1.430395e-16</td></tr>
	<tr><th scope=row>ENSG00000142534</th><td>13.93954    </td><td>13.75007    </td><td>13.95506    </td><td>13.72552    </td><td>13.84255    </td><td>23133.285   </td><td>2.726058e-20</td><td>1.430395e-16</td></tr>
	<tr><th scope=row>ENSG00000198637</th><td>14.13433    </td><td>14.11080    </td><td>14.23767    </td><td>14.03697    </td><td>14.12994    </td><td>22955.645   </td><td>2.836398e-20</td><td>1.430395e-16</td></tr>
	<tr><th scope=row>ENSG00000133112</th><td>13.78061    </td><td>13.78450    </td><td>13.77267    </td><td>13.56294    </td><td>13.72518    </td><td>20350.440   </td><td>5.272850e-20</td><td>2.014582e-16</td></tr>
	<tr><th scope=row>ENSG00000198034</th><td>13.51314    </td><td>13.80660    </td><td>13.71761    </td><td>13.38727    </td><td>13.60616    </td><td>19851.026   </td><td>5.992212e-20</td><td>2.014582e-16</td></tr>
	<tr><th scope=row>ENSG00000140988</th><td>14.10058    </td><td>14.14974    </td><td>14.22823    </td><td>14.07113    </td><td>14.13742    </td><td>16599.698   </td><td>1.504609e-19</td><td>4.185933e-16</td></tr>
	<tr><th scope=row>ENSG00000233927</th><td>13.32112    </td><td>13.23596    </td><td>13.42129    </td><td>13.17491    </td><td>13.28832    </td><td>16016.015   </td><td>1.808992e-19</td><td>4.185933e-16</td></tr>
	<tr><th scope=row>ENSG00000166710</th><td>12.99352    </td><td>13.13625    </td><td>13.27690    </td><td>13.31242    </td><td>13.17977    </td><td>15917.092   </td><td>1.867608e-19</td><td>4.185933e-16</td></tr>
	<tr><th scope=row>ENSG00000109475</th><td>14.08755    </td><td>13.98899    </td><td>14.11584    </td><td>13.77873    </td><td>13.99278    </td><td>14880.393   </td><td>2.641368e-19</td><td>4.813734e-16</td></tr>
	<tr><th scope=row>ENSG00000105193</th><td>13.29223    </td><td>13.12494    </td><td>13.47154    </td><td>13.10956    </td><td>13.24957    </td><td>14837.774   </td><td>2.680649e-19</td><td>4.813734e-16</td></tr>
	<tr><th scope=row>ENSG00000198918</th><td>13.72933    </td><td>13.80642    </td><td>13.73162    </td><td>13.66407    </td><td>13.73286    </td><td>14639.246   </td><td>2.873091e-19</td><td>4.813734e-16</td></tr>
	<tr><th scope=row>ENSG00000156482</th><td>13.45121    </td><td>13.36451    </td><td>13.64896    </td><td>13.24524    </td><td>13.42748    </td><td>13028.933   </td><td>5.233919e-19</td><td>7.541329e-16</td></tr>
	<tr><th scope=row>ENSG00000108298</th><td>13.52101    </td><td>13.21202    </td><td>13.68027    </td><td>13.15634    </td><td>13.39241    </td><td>12122.278   </td><td>7.586549e-19</td><td>9.002110e-16</td></tr>
	<tr><th scope=row>ENSG00000170315</th><td>13.04716    </td><td>13.48633    </td><td>13.18056    </td><td>13.64973    </td><td>13.34094    </td><td>11905.139   </td><td>8.326147e-19</td><td>9.297433e-16</td></tr>
	<tr><th scope=row>ENSG00000135486</th><td>13.33097    </td><td>13.41526    </td><td>13.29131    </td><td>13.26601    </td><td>13.32589    </td><td>11788.936   </td><td>8.757249e-19</td><td>9.297433e-16</td></tr>
	<tr><th scope=row>ENSG00000110700</th><td>12.99999    </td><td>13.11217    </td><td>13.08990    </td><td>12.98353    </td><td>13.04640    </td><td>11063.208   </td><td>1.214449e-18</td><td>1.113539e-15</td></tr>
	<tr><th scope=row>ENSG00000149806</th><td>12.99400    </td><td>12.99275    </td><td>13.31644    </td><td>12.85388    </td><td>13.03927    </td><td>10875.890   </td><td>1.326010e-18</td><td>1.162968e-15</td></tr>
	<tr><th scope=row>ENSG00000169567</th><td>11.49506    </td><td>12.05842    </td><td>12.06379    </td><td>11.50760    </td><td>11.78122    </td><td> 9895.511   </td><td>2.156186e-18</td><td>1.739784e-15</td></tr>
	<tr><th scope=row>ENSG00000143543</th><td>10.81122    </td><td>11.22164    </td><td>11.51480    </td><td>11.23056    </td><td>11.19456    </td><td> 9391.415   </td><td>2.821917e-18</td><td>2.089925e-15</td></tr>
	<tr><th scope=row>ENSG00000137154</th><td>13.37647    </td><td>13.46902    </td><td>13.53522    </td><td>13.32310    </td><td>13.42595    </td><td> 9369.702   </td><td>2.855733e-18</td><td>2.089925e-15</td></tr>
	<tr><th scope=row>ENSG00000232112</th><td>12.23607    </td><td>12.49263    </td><td>12.71735    </td><td>12.66173    </td><td>12.52694    </td><td> 9341.144   </td><td>2.900947e-18</td><td>2.089925e-15</td></tr>
	<tr><th scope=row>ENSG00000187514</th><td>13.36346    </td><td>13.26131    </td><td>13.01263    </td><td>12.66995    </td><td>13.07684    </td><td> 9267.683   </td><td>3.021244e-18</td><td>2.101536e-15</td></tr>
	<tr><th scope=row>ENSG00000118816</th><td>11.68461    </td><td>11.75722    </td><td>11.81722    </td><td>11.50098    </td><td>11.69001    </td><td> 9069.223   </td><td>3.377271e-18</td><td>2.270877e-15</td></tr>
	<tr><th scope=row>ENSG00000174444</th><td>13.43060    </td><td>13.03562    </td><td>13.49196    </td><td>13.10332    </td><td>13.26538    </td><td> 8599.458   </td><td>4.440549e-18</td><td>2.889508e-15</td></tr>
	<tr><th scope=row>ENSG00000080824</th><td>12.16780    </td><td>12.09354    </td><td>12.35360    </td><td>12.02030    </td><td>12.15881    </td><td> 8450.570   </td><td>4.858143e-18</td><td>2.971319e-15</td></tr>
	<tr><th scope=row>ENSG00000122406</th><td>13.23121    </td><td>13.26863    </td><td>13.40527    </td><td>12.91735    </td><td>13.20562    </td><td> 8417.489   </td><td>4.957197e-18</td><td>2.971319e-15</td></tr>
	<tr><th scope=row>ENSG00000130255</th><td>12.88832    </td><td>12.83524    </td><td>13.18521    </td><td>12.72745    </td><td>12.90905    </td><td> 8386.328   </td><td>5.052718e-18</td><td>2.971319e-15</td></tr>
	<tr><th scope=row>ENSG00000108654</th><td>12.50890    </td><td>11.83196    </td><td>12.46445    </td><td>12.21406    </td><td>12.25484    </td><td> 8353.583   </td><td>5.155472e-18</td><td>2.971319e-15</td></tr>
	<tr><th scope=row>ENSG00000087460</th><td>12.24961    </td><td>11.91147    </td><td>12.51694    </td><td>12.07412    </td><td>12.18803    </td><td> 8231.482   </td><td>5.561299e-18</td><td>2.980396e-15</td></tr>
	<tr><th scope=row>ENSG00000221983</th><td>12.79697    </td><td>12.70643    </td><td>13.09305    </td><td>12.82913    </td><td>12.85639    </td><td> 8160.056   </td><td>5.816385e-18</td><td>2.980396e-15</td></tr>
	<tr><th scope=row>...</th><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td></tr>
	<tr><th scope=row>ENSG00000130821</th><td>5.214284    </td><td>4.974164    </td><td>5.648512    </td><td>8.325736    </td><td>6.040674    </td><td>1199.917    </td><td>1.111065e-13</td><td>7.206560e-13</td></tr>
	<tr><th scope=row>ENSG00000130518</th><td>6.986330    </td><td>5.407438    </td><td>6.937936    </td><td>4.750937    </td><td>6.020661    </td><td>1199.166    </td><td>1.114643e-13</td><td>7.215846e-13</td></tr>
	<tr><th scope=row>ENSG00000163521</th><td>4.912902    </td><td>5.218787    </td><td>5.275215    </td><td>5.390328    </td><td>5.199308    </td><td>1195.602    </td><td>1.131810e-13</td><td>7.284900e-13</td></tr>
	<tr><th scope=row>ENSG00000212992</th><td>4.823410    </td><td>4.225729    </td><td>4.196469    </td><td>4.097073    </td><td>4.335670    </td><td>1192.855    </td><td>1.145258e-13</td><td>7.307768e-13</td></tr>
	<tr><th scope=row>ENSG00000157890</th><td>3.678097    </td><td>3.586679    </td><td>3.865001    </td><td>3.647675    </td><td>3.694363    </td><td>1188.597    </td><td>1.166481e-13</td><td>7.371633e-13</td></tr>
	<tr><th scope=row>ENSG00000061936</th><td>8.061970    </td><td>7.394516    </td><td>8.229012    </td><td>7.432428    </td><td>7.779482    </td><td>1187.675    </td><td>1.171138e-13</td><td>7.378703e-13</td></tr>
	<tr><th scope=row>ENSG00000247357</th><td>5.139469    </td><td>5.839189    </td><td>4.948686    </td><td>5.061485    </td><td>5.247207    </td><td>1187.485    </td><td>1.172101e-13</td><td>7.378703e-13</td></tr>
	<tr><th scope=row>ENSG00000125831</th><td>3.624877    </td><td>3.533432    </td><td>3.682283    </td><td>3.397352    </td><td>3.559486    </td><td>1187.044    </td><td>1.174340e-13</td><td>7.378703e-13</td></tr>
	<tr><th scope=row>ENSG00000247412</th><td>3.577447    </td><td>3.399866    </td><td>3.701217    </td><td>3.341264    </td><td>3.504948    </td><td>1183.147    </td><td>1.194337e-13</td><td>7.438471e-13</td></tr>
	<tr><th scope=row>ENSG00000179562</th><td>6.386263    </td><td>5.811791    </td><td>6.411045    </td><td>5.690265    </td><td>6.074841    </td><td>1182.730    </td><td>1.196500e-13</td><td>7.442430e-13</td></tr>
	<tr><th scope=row>ENSG00000115902</th><td>6.367353    </td><td>7.003458    </td><td>6.293532    </td><td>6.773750    </td><td>6.609524    </td><td>1177.576    </td><td>1.223634e-13</td><td>7.515389e-13</td></tr>
	<tr><th scope=row>ENSG00000095139</th><td>7.880052    </td><td>9.166382    </td><td>9.011994    </td><td>8.991784    </td><td>8.762553    </td><td>1175.992    </td><td>1.232126e-13</td><td>7.533931e-13</td></tr>
	<tr><th scope=row>ENSG00000122718</th><td>3.829331    </td><td>3.879242    </td><td>4.253310    </td><td>4.102226    </td><td>4.016027    </td><td>1174.753    </td><td>1.238809e-13</td><td>7.558757e-13</td></tr>
	<tr><th scope=row>ENSG00000007080</th><td>6.900640    </td><td>6.561585    </td><td>6.623781    </td><td>6.463479    </td><td>6.637371    </td><td>1165.568    </td><td>1.289762e-13</td><td>7.724784e-13</td></tr>
	<tr><th scope=row>ENSG00000206073</th><td>3.108975    </td><td>3.119414    </td><td>3.128686    </td><td>3.056851    </td><td>3.103482    </td><td>1161.641    </td><td>1.312303e-13</td><td>7.790379e-13</td></tr>
	<tr><th scope=row>ENSG00000058085</th><td>3.748232    </td><td>3.578609    </td><td>3.720734    </td><td>3.558487    </td><td>3.651516    </td><td>1159.281    </td><td>1.326082e-13</td><td>7.835303e-13</td></tr>
	<tr><th scope=row>ENSG00000143368</th><td>8.068664    </td><td>9.029597    </td><td>8.010834    </td><td>8.573280    </td><td>8.420594    </td><td>1156.463    </td><td>1.342755e-13</td><td>7.894732e-13</td></tr>
	<tr><th scope=row>ENSG00000164299</th><td>3.082035    </td><td>3.031434    </td><td>3.007472    </td><td>2.970033    </td><td>3.022743    </td><td>1153.881    </td><td>1.358255e-13</td><td>7.933014e-13</td></tr>
	<tr><th scope=row>ENSG00000073670</th><td>4.968278    </td><td>4.803982    </td><td>5.093394    </td><td>4.572368    </td><td>4.859505    </td><td>1152.884    </td><td>1.364296e-13</td><td>7.933014e-13</td></tr>
	<tr><th scope=row>ENSG00000179335</th><td>8.615910    </td><td>7.672596    </td><td>8.495549    </td><td>8.095217    </td><td>8.219818    </td><td>1151.646    </td><td>1.371843e-13</td><td>7.958818e-13</td></tr>
	<tr><th scope=row>ENSG00000243489</th><td>6.172070    </td><td>6.262636    </td><td>6.711030    </td><td>6.752572    </td><td>6.474577    </td><td>1151.134    </td><td>1.374982e-13</td><td>7.974736e-13</td></tr>
	<tr><th scope=row>ENSG00000166888</th><td>7.985756    </td><td>8.265127    </td><td>8.811501    </td><td>7.776192    </td><td>8.209644    </td><td>1148.875    </td><td>1.388920e-13</td><td>8.017814e-13</td></tr>
	<tr><th scope=row>ENSG00000163719</th><td>8.576895    </td><td>8.257451    </td><td>8.679623    </td><td>8.515046    </td><td>8.507254    </td><td>1147.992    </td><td>1.394410e-13</td><td>8.034289e-13</td></tr>
	<tr><th scope=row>ENSG00000137869</th><td>3.779607    </td><td>3.532846    </td><td>3.691546    </td><td>3.652205    </td><td>3.664051    </td><td>1147.578    </td><td>1.396995e-13</td><td>8.046882e-13</td></tr>
	<tr><th scope=row>ENSG00000233398</th><td>3.874477    </td><td>3.706202    </td><td>3.656185    </td><td>3.668896    </td><td>3.726440    </td><td>1147.227    </td><td>1.399191e-13</td><td>8.047740e-13</td></tr>
	<tr><th scope=row>ENSG00000105701</th><td>6.903015    </td><td>6.954748    </td><td>7.146769    </td><td>6.827252    </td><td>6.957946    </td><td>1145.809    </td><td>1.408107e-13</td><td>8.064829e-13</td></tr>
	<tr><th scope=row>ENSG00000143207</th><td>7.004623    </td><td>7.892597    </td><td>7.313433    </td><td>7.807295    </td><td>7.504487    </td><td>1144.101    </td><td>1.418934e-13</td><td>8.113023e-13</td></tr>
	<tr><th scope=row>ENSG00000123737</th><td>8.388098    </td><td>7.678710    </td><td>8.530180    </td><td>8.109720    </td><td>8.176677    </td><td>1138.635    </td><td>1.454253e-13</td><td>8.221206e-13</td></tr>
	<tr><th scope=row>ENSG00000214160</th><td>6.994047    </td><td>7.913299    </td><td>7.110934    </td><td>7.412355    </td><td>7.357659    </td><td>1138.421    </td><td>1.455661e-13</td><td>8.223719e-13</td></tr>
	<tr><th scope=row>ENSG00000109846</th><td>4.249335    </td><td>4.032090    </td><td>4.224073    </td><td>4.154921    </td><td>4.165105    </td><td>1131.348    </td><td>1.502994e-13</td><td>8.356778e-13</td></tr>
</tbody>
</table>



The working of the preceding code is straightforward and self-explanatory. The log-fold changes are computed based on the final (treatment) and initial (control) values. The log used is to the base of 2. The volcano plot simply creates the plot of log fold change to log odds in the data. The plotting is a simple scatter plot with log fold changes along the x axis and –log(p-values) along the y axis. Transforming the p-values into a log scale gives better resolution for visualization.
