
# Getting a co-expression network from microarray data
   Recently, biologists are prone to generate some interactive networks from gene expression data to represent the relations between the genes based on the data. There are many possible ways to plot these relationships from the data. This section will explore the relations based on the correlation among the genes.<br>
   最近，生物学家倾向于从基因表达数据中产生一些交互式网络，以表示基于数据的基因之间的关系。 从数据中绘制这些关系有很多种可能的方法。 本节将根据基因之间的相关性探讨它们之间的关系。<br>
 
Requires:<br>
1、A dataset().<br>
2、`WGCNA` and `RBGL` packages.<br>

The `eset` has annotation of `hgu133a2`. 


```R
library(leukemiasEset)
data(leukemiasEset)
pheno <- pData(leukemiasEset)
mydata <- leukemiasEset[, sampleNames(leukemiasEset)[c(1:3, 13:15, 25:27, 49:51)]]
mydata
design <- model.matrix(~0 + factor(pData(mydata)$LeukemiaType))
colnames(design) <- unique(as.character(pData(mydata)$LeukemiaType))
library(limma)
fit <- lmFit(mydata, design)
contrast.matrix <- makeContrasts(NoL- ALL, NoL- AML, NoL- CLL, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
tested2 <- topTable(fit2,adjust="fdr",sort.by="B",number=Inf, coef=1)
DE2 <- tested2[tested2$adj.P.Val < 0.01,]
```


    ExpressionSet (storageMode: lockedEnvironment)
    assayData: 20172 features, 12 samples 
      element names: exprs, se.exprs 
    protocolData
      sampleNames: GSM330151.CEL GSM330153.CEL ... GSM331663.CEL (12 total)
      varLabels: ScanDate
      varMetadata: labelDescription
    phenoData
      sampleNames: GSM330151.CEL GSM330153.CEL ... GSM331663.CEL (12 total)
      varLabels: Project Tissue ... Subtype (5 total)
      varMetadata: labelDescription
    featureData: none
    experimentData: use 'experimentData(object)'
    Annotation: genemapperhgu133plus2 



```R
library(affy)
library(affydata)
library(hgu133a2cdf)
mydata <- ReadAffy(celfile.path= "D:/Try-practice/Chapter 5/GSE24460_RAW/")
# mydata <- ReadAffy(filenames="D:/Try-practice/Chapter 5/GSE24460_RAW/GSM602658_MCF71.CEL")
library(Biobase)
DIR <- system.file("extdata", package="Biobase")
exprsLoc <- file.path(DIR, "exprsData.txt") 
pDataLoc <- file.path(DIR, "pData.txt")
DIR
exprsLoc
pDataLoc
exprs <- as.matrix(read.csv(exprsLoc, header = TRUE, sep = "\t", row.names = 1, as.is = TRUE))
pData <- read.table(pDataLoc, row.names = 1, header = TRUE, sep = "\t")
pData <- new("AnnotatedDataFrame", data = pData)
exData <- new("MIAME", name="ABCabc", lab="XYZ Lab", contact="abc@xyz", title="", abstract="", url="www.xyz")
eset <- new("ExpressionSet", exprs = exprs, phenoData = pData, experimentData = exData, annotation = "hgu133a2")
eset
```


'C:/Users/Master1/Documents/R/win-library/3.4/Biobase/extdata'



'C:/Users/Master1/Documents/R/win-library/3.4/Biobase/extdata/exprsData.txt'



'C:/Users/Master1/Documents/R/win-library/3.4/Biobase/extdata/pData.txt'



    ExpressionSet (storageMode: lockedEnvironment)
    assayData: 500 features, 26 samples 
      element names: exprs 
    protocolData: none
    phenoData
      sampleNames: A B ... Z (26 total)
      varLabels: gender type score
      varMetadata: labelDescription
    featureData: none
    experimentData: use 'experimentData(object)'
    Annotation: hgu133a2 



```R
library(WGCNA)
library(RBGL)
```

1、Take the dataset, and use only significant genes at this stage to reduce the noise and simultaneously consume less time. Only 50 genes have been used.


```R
a <- rownames(DE2)[1:25]
a
eset[1:25, 1:3]
```


<ol class=list-inline>
	<li>'ENSG00000152078'</li>
	<li>'ENSG00000117519'</li>
	<li>'ENSG00000145850'</li>
	<li>'ENSG00000170180'</li>
	<li>'ENSG00000087586'</li>
	<li>'ENSG00000047597'</li>
	<li>'ENSG00000175449'</li>
	<li>'ENSG00000104043'</li>
	<li>'ENSG00000024526'</li>
	<li>'ENSG00000115641'</li>
	<li>'ENSG00000188672'</li>
	<li>'ENSG00000133063'</li>
	<li>'ENSG00000164330'</li>
	<li>'ENSG00000173372'</li>
	<li>'ENSG00000234955'</li>
	<li>'ENSG00000102935'</li>
	<li>'ENSG00000159189'</li>
	<li>'ENSG00000107562'</li>
	<li>'ENSG00000086506'</li>
	<li>'ENSG00000130635'</li>
	<li>'ENSG00000066923'</li>
	<li>'ENSG00000162692'</li>
	<li>'ENSG00000198336'</li>
	<li>'ENSG00000109501'</li>
	<li>'ENSG00000168497'</li>
</ol>




    ExpressionSet (storageMode: lockedEnvironment)
    assayData: 25 features, 3 samples 
      element names: exprs 
    protocolData: none
    phenoData
      sampleNames: A B C
      varLabels: gender type score
      varMetadata: labelDescription
    featureData: none
    experimentData: use 'experimentData(object)'
    Annotation: hgu133a2 



```R
myData_Sel <- exprs(eset[1:25, 1:3])
myData_Sel
```


<table>
<thead><tr><th></th><th scope=col>A</th><th scope=col>B</th><th scope=col>C</th></tr></thead>
<tbody>
	<tr><th scope=row>AFFX-MurIL2_at</th><td> 192.7420 </td><td>  85.75330</td><td> 176.7570 </td></tr>
	<tr><th scope=row>AFFX-MurIL10_at</th><td>  97.1370 </td><td> 126.19600</td><td>  77.9216 </td></tr>
	<tr><th scope=row>AFFX-MurIL4_at</th><td>  45.8192 </td><td>   8.83135</td><td>  33.0632 </td></tr>
	<tr><th scope=row>AFFX-MurFAS_at</th><td>  22.5445 </td><td>   3.60093</td><td>  14.6883 </td></tr>
	<tr><th scope=row>AFFX-BioB-5_at</th><td>  96.7875 </td><td>  30.43800</td><td>  46.1271 </td></tr>
	<tr><th scope=row>AFFX-BioB-M_at</th><td>  89.0730 </td><td>  25.84610</td><td>  57.2033 </td></tr>
	<tr><th scope=row>AFFX-BioB-3_at</th><td> 265.9640 </td><td> 181.08000</td><td> 164.9260 </td></tr>
	<tr><th scope=row>AFFX-BioC-5_at</th><td> 110.1360 </td><td>  57.28890</td><td>  67.3980 </td></tr>
	<tr><th scope=row>AFFX-BioC-3_at</th><td>  43.0794 </td><td>  16.80060</td><td>  37.6002 </td></tr>
	<tr><th scope=row>AFFX-BioDn-5_at</th><td>  10.9187 </td><td>  16.17890</td><td>  10.1495 </td></tr>
	<tr><th scope=row>AFFX-BioDn-3_at</th><td> 751.2270 </td><td> 515.00400</td><td> 622.9010 </td></tr>
	<tr><th scope=row>AFFX-CreX-5_at</th><td>  76.9437 </td><td>  40.90700</td><td>  62.0314 </td></tr>
	<tr><th scope=row>AFFX-CreX-3_at</th><td> 105.3780 </td><td>  97.49320</td><td>  74.0299 </td></tr>
	<tr><th scope=row>AFFX-BioB-5_st</th><td>  40.4826 </td><td>   7.45801</td><td>  19.4069 </td></tr>
	<tr><th scope=row>AFFX-BioB-M_st</th><td>  58.1706 </td><td>  15.79260</td><td>  25.1962 </td></tr>
	<tr><th scope=row>AFFX-BioB-3_st</th><td> 257.6190 </td><td> 113.69000</td><td> 187.7960 </td></tr>
	<tr><th scope=row>AFFX-BioC-5_st</th><td> 129.0560 </td><td>  74.60950</td><td>  82.8271 </td></tr>
	<tr><th scope=row>AFFX-BioC-3_st</th><td>  61.7251 </td><td>  50.23720</td><td>  61.6710 </td></tr>
	<tr><th scope=row>AFFX-BioDn-5_st</th><td> -40.9349 </td><td> -83.93020</td><td> -28.7050 </td></tr>
	<tr><th scope=row>AFFX-BioDn-3_st</th><td> 284.4070 </td><td> 208.09900</td><td> 239.0390 </td></tr>
	<tr><th scope=row>AFFX-CreX-5_st</th><td> 178.7450 </td><td> 101.30000</td><td> 118.6990 </td></tr>
	<tr><th scope=row>AFFX-CreX-3_st</th><td>  79.7368 </td><td>  55.56320</td><td>  68.5976 </td></tr>
	<tr><th scope=row>AFFX-hum_alu_at</th><td>9903.1900 </td><td>8501.62000</td><td>9453.0000 </td></tr>
	<tr><th scope=row>AFFX-DapX-5_at</th><td>  61.2671 </td><td>  37.47400</td><td>  44.7525 </td></tr>
	<tr><th scope=row>AFFX-DapX-M_at</th><td> 120.5440 </td><td>  75.98540</td><td> 126.3740 </td></tr>
</tbody>
</table>



2、The data has sample names in the columns and genes in the rows. So transpose the data.


```R
myData_Sel <- t(myData_Sel)
myData_Sel
```


<table>
<thead><tr><th></th><th scope=col>AFFX-MurIL2_at</th><th scope=col>AFFX-MurIL10_at</th><th scope=col>AFFX-MurIL4_at</th><th scope=col>AFFX-MurFAS_at</th><th scope=col>AFFX-BioB-5_at</th><th scope=col>AFFX-BioB-M_at</th><th scope=col>AFFX-BioB-3_at</th><th scope=col>AFFX-BioC-5_at</th><th scope=col>AFFX-BioC-3_at</th><th scope=col>AFFX-BioDn-5_at</th><th scope=col>...</th><th scope=col>AFFX-BioB-3_st</th><th scope=col>AFFX-BioC-5_st</th><th scope=col>AFFX-BioC-3_st</th><th scope=col>AFFX-BioDn-5_st</th><th scope=col>AFFX-BioDn-3_st</th><th scope=col>AFFX-CreX-5_st</th><th scope=col>AFFX-CreX-3_st</th><th scope=col>AFFX-hum_alu_at</th><th scope=col>AFFX-DapX-5_at</th><th scope=col>AFFX-DapX-M_at</th></tr></thead>
<tbody>
	<tr><th scope=row>A</th><td>192.7420</td><td> 97.1370</td><td>45.81920</td><td>22.54450</td><td>96.7875 </td><td>89.0730 </td><td>265.964 </td><td>110.1360</td><td>43.0794 </td><td>10.9187 </td><td>...     </td><td>257.619 </td><td>129.0560</td><td>61.7251 </td><td>-40.9349</td><td>284.407 </td><td>178.745 </td><td>79.7368 </td><td>9903.19 </td><td>61.2671 </td><td>120.5440</td></tr>
	<tr><th scope=row>B</th><td> 85.7533</td><td>126.1960</td><td> 8.83135</td><td> 3.60093</td><td>30.4380 </td><td>25.8461 </td><td>181.080 </td><td> 57.2889</td><td>16.8006 </td><td>16.1789 </td><td>...     </td><td>113.690 </td><td> 74.6095</td><td>50.2372 </td><td>-83.9302</td><td>208.099 </td><td>101.300 </td><td>55.5632 </td><td>8501.62 </td><td>37.4740 </td><td> 75.9854</td></tr>
	<tr><th scope=row>C</th><td>176.7570</td><td> 77.9216</td><td>33.06320</td><td>14.68830</td><td>46.1271 </td><td>57.2033 </td><td>164.926 </td><td> 67.3980</td><td>37.6002 </td><td>10.1495 </td><td>...     </td><td>187.796 </td><td> 82.8271</td><td>61.6710 </td><td>-28.7050</td><td>239.039 </td><td>118.699 </td><td>68.5976 </td><td>9453.00 </td><td>44.7525 </td><td>126.3740</td></tr>
</tbody>
</table>



3、Either run the correlation computation (it takes a lot of time) or use the adjacency function from `WGCNA` to compute the adjacency matrix.<br>计算矩阵.


```R
myMat <- adjacency(myData_Sel, type="signed")
myMat
# The results give a square matrix of size equal to the number of genes. 
# Each entry in the matrix is the connectivity between the genes.
```


<table>
<thead><tr><th></th><th scope=col>AFFX-MurIL2_at</th><th scope=col>AFFX-MurIL10_at</th><th scope=col>AFFX-MurIL4_at</th><th scope=col>AFFX-MurFAS_at</th><th scope=col>AFFX-BioB-5_at</th><th scope=col>AFFX-BioB-M_at</th><th scope=col>AFFX-BioB-3_at</th><th scope=col>AFFX-BioC-5_at</th><th scope=col>AFFX-BioC-3_at</th><th scope=col>AFFX-BioDn-5_at</th><th scope=col>...</th><th scope=col>AFFX-BioB-3_st</th><th scope=col>AFFX-BioC-5_st</th><th scope=col>AFFX-BioC-3_st</th><th scope=col>AFFX-BioDn-5_st</th><th scope=col>AFFX-BioDn-3_st</th><th scope=col>AFFX-CreX-5_st</th><th scope=col>AFFX-CreX-3_st</th><th scope=col>AFFX-hum_alu_at</th><th scope=col>AFFX-DapX-5_at</th><th scope=col>AFFX-DapX-M_at</th></tr></thead>
<tbody>
	<tr><th scope=row>AFFX-MurIL2_at</th><td>1.000000e+00</td><td>1.455403e-07</td><td>9.373990e-01</td><td>8.837797e-01</td><td>4.928595e-01</td><td>7.954628e-01</td><td>0.171623608 </td><td>0.444372324 </td><td>9.946171e-01</td><td>1.907480e-11</td><td>...         </td><td>8.156942e-01</td><td>0.4037651381</td><td>9.730722e-01</td><td>8.301657e-01</td><td>6.914793e-01</td><td>4.799091e-01</td><td>8.405878e-01</td><td>9.519479e-01</td><td>5.724229e-01</td><td>9.137639e-01</td></tr>
	<tr><th scope=row>AFFX-MurIL10_at</th><td>1.455403e-07</td><td>1.000000e+00</td><td>6.079850e-06</td><td>1.892503e-05</td><td>1.308513e-03</td><td>6.990998e-05</td><td>0.018894884 </td><td>0.001935609 </td><td>4.917313e-07</td><td>8.819052e-01</td><td>...         </td><td>5.371687e-05</td><td>0.0026759431</td><td>5.137398e-09</td><td>6.808438e-13</td><td>2.235073e-04</td><td>1.453855e-03</td><td>3.789065e-05</td><td>4.043552e-06</td><td>6.722836e-04</td><td>1.263022e-10</td></tr>
	<tr><th scope=row>AFFX-MurIL4_at</th><td>9.373990e-01</td><td>6.079850e-06</td><td>1.000000e+00</td><td>9.906591e-01</td><td>7.130040e-01</td><td>9.514999e-01</td><td>0.326727380 </td><td>0.663826782 </td><td>9.678930e-01</td><td>2.172729e-08</td><td>...         </td><td>9.622561e-01</td><td>0.6204134991</td><td>8.381323e-01</td><td>6.221678e-01</td><td>8.841767e-01</td><td>7.001490e-01</td><td>9.741780e-01</td><td>9.989593e-01</td><td>7.875276e-01</td><td>7.337988e-01</td></tr>
	<tr><th scope=row>AFFX-MurFAS_at</th><td>8.837797e-01</td><td>1.892503e-05</td><td>9.906591e-01</td><td>1.000000e+00</td><td>7.914758e-01</td><td>9.842539e-01</td><td>0.400524968 </td><td>0.745604017 </td><td>9.257679e-01</td><td>1.380911e-07</td><td>...         </td><td>9.902084e-01</td><td>0.7040544782</td><td>7.646807e-01</td><td>5.374527e-01</td><td>9.377028e-01</td><td>7.796149e-01</td><td>9.958063e-01</td><td>9.834561e-01</td><td>8.581795e-01</td><td>6.513910e-01</td></tr>
	<tr><th scope=row>AFFX-BioB-5_at</th><td>4.928595e-01</td><td>1.308513e-03</td><td>7.130040e-01</td><td>7.914758e-01</td><td>1.000000e+00</td><td>8.805309e-01</td><td>0.806449161 </td><td>0.996688741 </td><td>5.564141e-01</td><td>7.140448e-05</td><td>...         </td><td>8.631531e-01</td><td>0.9885134748</td><td>3.588041e-01</td><td>1.878861e-01</td><td>9.490007e-01</td><td>9.997686e-01</td><td>8.395418e-01</td><td>6.856132e-01</td><td>9.916117e-01</td><td>2.639864e-01</td></tr>
	<tr><th scope=row>AFFX-BioB-M_at</th><td>7.954628e-01</td><td>6.990998e-05</td><td>9.514999e-01</td><td>9.842539e-01</td><td>8.805309e-01</td><td>1.000000e+00</td><td>0.505217865 </td><td>0.842040128 </td><td>8.492020e-01</td><td>1.041157e-06</td><td>...         </td><td>9.992835e-01</td><td>0.8056096836</td><td>6.590401e-01</td><td>4.306478e-01</td><td>9.838835e-01</td><td>8.707699e-01</td><td>9.962696e-01</td><td>9.368935e-01</td><td>9.324057e-01</td><td>5.413067e-01</td></tr>
	<tr><th scope=row>AFFX-BioB-3_at</th><td>1.716236e-01</td><td>1.889488e-02</td><td>3.267274e-01</td><td>4.005250e-01</td><td>8.064492e-01</td><td>5.052179e-01</td><td>1.000000000 </td><td>0.848199089 </td><td>2.102798e-01</td><td>2.694197e-03</td><td>...         </td><td>4.823420e-01</td><td>0.8812318972</td><td>1.030196e-01</td><td>3.828800e-02</td><td>6.163032e-01</td><td>8.178186e-01</td><td>4.534161e-01</td><td>3.037795e-01</td><td>7.336065e-01</td><td>6.404777e-02</td></tr>
	<tr><th scope=row>AFFX-BioC-5_at</th><td>4.443723e-01</td><td>1.935609e-03</td><td>6.638268e-01</td><td>7.456040e-01</td><td>9.966887e-01</td><td>8.420401e-01</td><td>0.848199089 </td><td>1.000000000 </td><td>5.064631e-01</td><td>1.227810e-04</td><td>...         </td><td>8.228123e-01</td><td>0.9975150007</td><td>3.164124e-01</td><td>1.594107e-01</td><td>9.211746e-01</td><td>9.982059e-01</td><td>7.970447e-01</td><td>6.357969e-01</td><td>9.779256e-01</td><td>2.283942e-01</td></tr>
	<tr><th scope=row>AFFX-BioC-3_at</th><td>9.946171e-01</td><td>4.917313e-07</td><td>9.678930e-01</td><td>9.257679e-01</td><td>5.564141e-01</td><td>8.492020e-01</td><td>0.210279772 </td><td>0.506463052 </td><td>1.000000e+00</td><td>2.336069e-10</td><td>...         </td><td>8.672758e-01</td><td>0.4640452053</td><td>9.445672e-01</td><td>7.744525e-01</td><td>7.527638e-01</td><td>5.431469e-01</td><td>8.891232e-01</td><td>9.782344e-01</td><td>6.367208e-01</td><td>8.693896e-01</td></tr>
	<tr><th scope=row>AFFX-BioDn-5_at</th><td>1.907480e-11</td><td>8.819052e-01</td><td>2.172729e-08</td><td>1.380911e-07</td><td>7.140448e-05</td><td>1.041157e-06</td><td>0.002694197 </td><td>0.000122781 </td><td>2.336069e-10</td><td>1.000000e+00</td><td>...         </td><td>6.984189e-07</td><td>0.0001916592</td><td>2.521861e-15</td><td>1.304229e-16</td><td>5.811910e-06</td><td>8.264866e-05</td><td>4.090860e-07</td><td>1.089949e-08</td><td>2.810378e-05</td><td>1.324426e-27</td></tr>
	<tr><th scope=row>AFFX-BioDn-3_at</th><td>7.516318e-01</td><td>1.179275e-04</td><td>9.253768e-01</td><td>9.676273e-01</td><td>9.133420e-01</td><td>9.969386e-01</td><td>0.553212017 </td><td>0.879126270 </td><td>8.092112e-01</td><td>2.274845e-06</td><td>...         </td><td>9.932743e-01</td><td>0.8458843696</td><td>6.107672e-01</td><td>3.862265e-01</td><td>9.948189e-01</td><td>9.047693e-01</td><td>9.865143e-01</td><td>9.078824e-01</td><td>9.572388e-01</td><td>4.935131e-01</td></tr>
	<tr><th scope=row>AFFX-CreX-5_at</th><td>8.845720e-01</td><td>1.866003e-05</td><td>9.909015e-01</td><td>9.999984e-01</td><td>7.904930e-01</td><td>9.839371e-01</td><td>0.399513629 </td><td>0.744564251 </td><td>9.264220e-01</td><td>1.350341e-07</td><td>...         </td><td>9.899571e-01</td><td>0.7029786262</td><td>7.656979e-01</td><td>5.385575e-01</td><td>9.370955e-01</td><td>7.786154e-01</td><td>9.956406e-01</td><td>9.837776e-01</td><td>8.573202e-01</td><td>6.524941e-01</td></tr>
	<tr><th scope=row>AFFX-CreX-3_at</th><td>6.320865e-03</td><td>2.966397e-01</td><td>2.304297e-02</td><td>3.530305e-02</td><td>1.880222e-01</td><td>5.855607e-02</td><td>0.520685249 </td><td>0.219665783 </td><td>9.449321e-03</td><td>1.110234e-01</td><td>...         </td><td>5.281722e-02</td><td>0.2496153439</td><td>2.328743e-03</td><td>3.376555e-04</td><td>9.276726e-02</td><td>1.960743e-01</td><td>4.610953e-02</td><td>1.983716e-02</td><td>1.440959e-01</td><td>9.233380e-04</td></tr>
	<tr><th scope=row>AFFX-BioB-5_st</th><td>6.391427e-01</td><td>3.703120e-04</td><td>8.439924e-01</td><td>9.058586e-01</td><td>9.721459e-01</td><td>9.652202e-01</td><td>0.669042109 </td><td>0.950310301 </td><td>7.024019e-01</td><td>1.204001e-05</td><td>...         </td><td>9.548422e-01</td><td>0.9267714474</td><td>4.952990e-01</td><td>2.888861e-01</td><td>9.963321e-01</td><td>9.669563e-01</td><td>9.396947e-01</td><td>8.208199e-01</td><td>9.942034e-01</td><td>3.844246e-01</td></tr>
	<tr><th scope=row>AFFX-BioB-M_st</th><td>4.769028e-01</td><td>1.489708e-03</td><td>6.971358e-01</td><td>7.768210e-01</td><td>9.996481e-01</td><td>8.684503e-01</td><td>0.820436198 </td><td>0.998493467 </td><td>5.400593e-01</td><td>8.548812e-05</td><td>...         </td><td>8.504354e-01</td><td>0.9921602490</td><td>3.447051e-01</td><td>1.782646e-01</td><td>9.405370e-01</td><td>9.999874e-01</td><td>8.260759e-01</td><td>6.694920e-01</td><td>9.878533e-01</td><td>2.520524e-01</td></tr>
	<tr><th scope=row>AFFX-BioB-3_st</th><td>8.156942e-01</td><td>5.371687e-05</td><td>9.622561e-01</td><td>9.902084e-01</td><td>8.631531e-01</td><td>9.992835e-01</td><td>0.482342009 </td><td>0.822812327 </td><td>8.672758e-01</td><td>6.984189e-07</td><td>...         </td><td>1.000000e+00</td><td>0.7850410401</td><td>6.821172e-01</td><td>4.527724e-01</td><td>9.764847e-01</td><td>8.528730e-01</td><td>9.988201e-01</td><td>9.491369e-01</td><td>9.185866e-01</td><td>5.646676e-01</td></tr>
	<tr><th scope=row>AFFX-BioC-5_st</th><td>4.037651e-01</td><td>2.675943e-03</td><td>6.204135e-01</td><td>7.040545e-01</td><td>9.885135e-01</td><td>8.056097e-01</td><td>0.881231897 </td><td>0.997515001 </td><td>4.640452e-01</td><td>1.916592e-04</td><td>...         </td><td>7.850410e-01</td><td>1.0000000000</td><td>2.819303e-01</td><td>1.372580e-01</td><td>8.928695e-01</td><td>9.915233e-01</td><td>7.577506e-01</td><td>5.921533e-01</td><td>9.610399e-01</td><td>2.000940e-01</td></tr>
	<tr><th scope=row>AFFX-BioC-3_st</th><td>9.730722e-01</td><td>5.137398e-09</td><td>8.381323e-01</td><td>7.646807e-01</td><td>3.588041e-01</td><td>6.590401e-01</td><td>0.103019606 </td><td>0.316412362 </td><td>9.445672e-01</td><td>2.521861e-15</td><td>...         </td><td>6.821172e-01</td><td>0.2819302798</td><td>1.000000e+00</td><td>9.319777e-01</td><td>5.476961e-01</td><td>3.473502e-01</td><td>7.113295e-01</td><td>8.604234e-01</td><td>4.313557e-01</td><td>9.819911e-01</td></tr>
	<tr><th scope=row>AFFX-BioDn-5_st</th><td>8.301657e-01</td><td>6.808438e-13</td><td>6.221678e-01</td><td>5.374527e-01</td><td>1.878861e-01</td><td>4.306478e-01</td><td>0.038288002 </td><td>0.159410681 </td><td>7.744525e-01</td><td>1.304229e-16</td><td>...         </td><td>4.527724e-01</td><td>0.1372579694</td><td>9.319777e-01</td><td>1.000000e+00</td><td>3.316036e-01</td><td>1.800582e-01</td><td>4.816862e-01</td><td>6.503014e-01</td><td>2.397785e-01</td><td>9.831357e-01</td></tr>
	<tr><th scope=row>AFFX-BioDn-3_st</th><td>6.914793e-01</td><td>2.235073e-04</td><td>8.841767e-01</td><td>9.377028e-01</td><td>9.490007e-01</td><td>9.838835e-01</td><td>0.616303198 </td><td>0.921174626 </td><td>7.527638e-01</td><td>5.811910e-06</td><td>...         </td><td>9.764847e-01</td><td>0.8928694517</td><td>5.476961e-01</td><td>3.316036e-01</td><td>1.000000e+00</td><td>9.421829e-01</td><td>9.650327e-01</td><td>8.633617e-01</td><td>9.814352e-01</td><td>4.330638e-01</td></tr>
	<tr><th scope=row>AFFX-CreX-5_st</th><td>4.799091e-01</td><td>1.453855e-03</td><td>7.001490e-01</td><td>7.796149e-01</td><td>9.997686e-01</td><td>8.707699e-01</td><td>0.817818616 </td><td>0.998205855 </td><td>5.431469e-01</td><td>8.264866e-05</td><td>...         </td><td>8.528730e-01</td><td>0.9915232507</td><td>3.473502e-01</td><td>1.800582e-01</td><td>9.421829e-01</td><td>1.000000e+00</td><td>8.286516e-01</td><td>6.725497e-01</td><td>9.886163e-01</td><td>2.542841e-01</td></tr>
	<tr><th scope=row>AFFX-CreX-3_st</th><td>8.405878e-01</td><td>3.789065e-05</td><td>9.741780e-01</td><td>9.958063e-01</td><td>8.395418e-01</td><td>9.962696e-01</td><td>0.453416126 </td><td>0.797044736 </td><td>8.891232e-01</td><td>4.090860e-07</td><td>...         </td><td>9.988201e-01</td><td>0.7577506037</td><td>7.113295e-01</td><td>4.816862e-01</td><td>9.650327e-01</td><td>8.286516e-01</td><td>1.000000e+00</td><td>9.630430e-01</td><td>8.992321e-01</td><td>5.947582e-01</td></tr>
	<tr><th scope=row>AFFX-hum_alu_at</th><td>9.519479e-01</td><td>4.043552e-06</td><td>9.989593e-01</td><td>9.834561e-01</td><td>6.856132e-01</td><td>9.368935e-01</td><td>0.303779521 </td><td>0.635796884 </td><td>9.782344e-01</td><td>1.089949e-08</td><td>...         </td><td>9.491369e-01</td><td>0.5921532755</td><td>8.604234e-01</td><td>6.503014e-01</td><td>8.633617e-01</td><td>6.725497e-01</td><td>9.630430e-01</td><td>1.000000e+00</td><td>7.620073e-01</td><td>7.601477e-01</td></tr>
	<tr><th scope=row>AFFX-DapX-5_at</th><td>5.724229e-01</td><td>6.722836e-04</td><td>7.875276e-01</td><td>8.581795e-01</td><td>9.916117e-01</td><td>9.324057e-01</td><td>0.733606541 </td><td>0.977925574 </td><td>6.367208e-01</td><td>2.810378e-05</td><td>...         </td><td>9.185866e-01</td><td>0.9610399151</td><td>4.313557e-01</td><td>2.397785e-01</td><td>9.814352e-01</td><td>9.886163e-01</td><td>8.992321e-01</td><td>7.620073e-01</td><td>1.000000e+00</td><td>3.268958e-01</td></tr>
	<tr><th scope=row>AFFX-DapX-M_at</th><td>9.137639e-01</td><td>1.263022e-10</td><td>7.337988e-01</td><td>6.513910e-01</td><td>2.639864e-01</td><td>5.413067e-01</td><td>0.064047766 </td><td>0.228394189 </td><td>8.693896e-01</td><td>1.324426e-27</td><td>...         </td><td>5.646676e-01</td><td>0.2000939591</td><td>9.819911e-01</td><td>9.831357e-01</td><td>4.330638e-01</td><td>2.542841e-01</td><td>5.947582e-01</td><td>7.601477e-01</td><td>3.268958e-01</td><td>1.000000e+00</td></tr>
</tbody>
</table>



4、The values in the resulting adjacency matrix can be set to either 0 (edge absent) or 1 (edge present) via dichotomization in different ways to get the final adjacency matrix. Simply use a threshold of 0.90 (a high value to get the most correlated vertices), as one of the simplest methods, but not the optimal one.<br>
可以用不同的方式使用二分法将得到的邻接矩阵中的值设置为0（边缘不存在）或1（边缘存在）以获得最终的邻接矩阵。 用0.90的阈值（高值来获得最相关的顶点）只是最简单的方法之一，但不是最优的方法.


```R
adjMat <- myMat
adjMat[abs(adjMat)>0.90] <- 1
adjMat[abs(adjMat)<=0.90] <- 0
diag(adjMat) <- 0
```

5、The final adjacency matrix can be converted into a `graphNEL` object to be rendered as a graph with nodes and edges(note that plotting the graph might take some time depending on the network size).<br>
最终的邻接矩阵可以被转换成一个`graphNEL`对象来渲染为一个带有节点和边的图（注意绘制图可能需要一些时间，这取决于网络的大小）.


```R
myGraph <- as(adjMat, "graphNEL")
myGraph
plot(myGraph, nodeAttrs=makeNodeAttrs(myGraph,fontsize=28, fillcolor="grey"))
```


    A graphNEL graph with undirected edges
    Number of Nodes = 25 
    Number of Edges = 86 



[plot](https://github.com/Chengshu21/Chapter-5-Analyzing-Microarray-Data-with-R/blob/master/md/pic/output-Getting%20a%20co-expression%20network%20from%20microarray%20data.png)


The method in this part is based on the computation of the relationship between the genes in terms of correlation or similarity measures. The function computes the pairwise similarity or correlation among the genes based on the expression data and returns this as a matrix. The threshold set defines only highly correlated or similar genes connected by an edge in the network otherwise no connection between the genes. This adjacency matrix is then
used to get the graph object.<br>
这个部分使用的方法是基于根据相关性或相似性度量计算基因之间的相互关系。 基于表达数据该函数计算基因之间的成对相似性或相关性，并将其作为矩阵返回。 阈值集仅定义了高度相关或相似的基因，通过网络中的边缘连接，否则基因之间没有连接。 然后这个邻接矩阵用来获得图形对象。


```R

```
