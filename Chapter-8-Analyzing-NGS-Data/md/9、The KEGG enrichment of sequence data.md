
# 9、The KEGG enrichment of sequence data

We have done the GO enrichment of sequence data, and we also can do a similar analysis in terms of KEGG annotations.
Take the same data and packages used in the GO analysis.

1、Start with loading the required packages, `goseq`, `edgeR` and `org.Hs.eg.db`. And load the data used for analysis from the `goseq` package.


```R
library(goseq)
library(edgeR)
library(org.Hs.eg.db)   #load the library
```

2、The first four columns in the data are controls and the last three are the treatment samples. Assign these attributes to the data and perform the differential tag computation.<br>
数据的前四列是对照，最后三列是实验样本。将这些属性分配给数据并执行差分标签计算。


```R
myData <- read.table(system.file("extdata", "Li_sum.txt",package = 'goseq'), 
                     sep = '\t', header = TRUE, stringsAsFactors = FALSE,row.names = 1)   #assign the data
# myData <- read.table("path/to/code/files/Li_sum.txt", sep = '\t', header = TRUE, stringsAsFactors = FALSE,row.names=1)
```


```R
myTreat <- factor(rep(c("Control","Treatment"),times = c(4,3)))
myDG <- DGEList(myData,lib.size = colSums(myData),group = myTreat)    #create `DGElist` object
myDisp <- estimateCommonDisp(myDG)
mytest <- exactTest(myDisp)
```

3、Use the genes from this analysis for enrichment, so extract the genes with the desired p-value and log fold change condition with their corresponding gene names, creating a named vector.<br>
使用该分析中的基因进行富集，因此用相应的基因名称提取具有期望的p值和对数倍数变化条件的基因，从而创建命名的载体。


```R
myTags <- as.integer(p.adjust(mytest$table$PValue[mytest$table$logFC!=0], method = "BH") < 0.05)
names(myTags) <- row.names(mytest$table[mytest$table$logFC!=0,])
```

4、Compile the KEGG data for enrichment, starting with the conversion of ENSEMBL IDs to Entrez.<br>
为KEGG富集整理数据，首先将ensembl id 转化为Entrez id.


```R
en2eg <- as.list(org.Hs.egENSEMBL2EG)
```

5、Get all the KEGG IDs for the compiled Entrez IDs, and get the KEGG and Entrez IDs mapped together.<br>
将所有的KEGG ids和entrez id 进行映射。


```R
eg2kegg <- as.list(org.Hs.egPATH)
grepKEGG <- function(id,mapkeys){unique(unlist(mapkeys[id], use.names = FALSE))}
kegg <- lapply(en2eg,grepKEGG,eg2kegg)
```

6、Compute the probability weighting function, and use the `goseq` function with KEGG mappings for the enrichment of tags.<br>
计算概率加权函数，并使用带有KEGG映射的goseq函数来丰富标签。


```R
pwf <- nullp(myTags,"hg19","ensGene")
KEGG <- goseq(pwf,gene2cat = kegg)
head(KEGG)
```

    Loading hg19 length data...
    Warning message in pcls(G):
    "initial point very close to some inequality constraints"Using manually entered categories.
    For 18659 genes, we could not find any categories. These genes will be excluded.
    To force their use, please run with use_genes_without_cat=TRUE (see documentation).
    This was the default behavior for version 1.15.1 and earlier.
    Calculating the p-values...
    


<table>
<thead><tr><th></th><th scope=col>category</th><th scope=col>over_represented_pvalue</th><th scope=col>under_represented_pvalue</th><th scope=col>numDEInCat</th><th scope=col>numInCat</th></tr></thead>
<tbody>
	<tr><th scope=row>88</th><td>03010       </td><td>1.222285e-05</td><td>0.9999960   </td><td>29          </td><td>89          </td></tr>
	<tr><th scope=row>77</th><td>00900       </td><td>2.395639e-04</td><td>0.9999709   </td><td>10          </td><td>15          </td></tr>
	<tr><th scope=row>113</th><td>04115       </td><td>8.190518e-04</td><td>0.9996824   </td><td>26          </td><td>64          </td></tr>
	<tr><th scope=row>175</th><td>04964       </td><td>2.154090e-03</td><td>0.9995918   </td><td>10          </td><td>17          </td></tr>
	<tr><th scope=row>27</th><td>00330       </td><td>3.677049e-03</td><td>0.9986561   </td><td>18          </td><td>44          </td></tr>
	<tr><th scope=row>20</th><td>00250       </td><td>5.209339e-03</td><td>0.9984326   </td><td>13          </td><td>28          </td></tr>
</tbody>
</table>




[output](https://github.com/Chengshu21/Chapter-8-Analyzing-NGS-Data/blob/master/MD/pic/output_9、The%20KEGG%20enrichment%20of%20sequence%20data.png)


The above plot shows the KEGG enrichment of the data.

This part is very similar to the previous part. However, in this part, its aim is to find the KEGG annotation in place of GO categories. Until step 3, we fetch the differentially expressed tags from data as we did in previous part followed by the extraction of genes of interest(meeting the p-value and log fold change values). Step 4 gets the Entrez mapping of the ENSEMBL genes, which is followed by the corresponding KEGG mapping (of Entrez genes)later. Finally, we extract the KEGG annotation of the interesting genes from the entire list and compute the enrichment scores.<br>
这部分与前一部分非常相似。但是，在这一部分中，其目标是找到KEGG注释来代替GO类别。直到第3步，我们从数据中获取差异表达的标签，就像我们在前一部分中所做的那样，然后提取感兴趣的基因（满足p值和log倍数变化值）。步骤4获得ENSEMBL基因的Entrez图谱，随后是相应的（Entrez基因的）KEGG图谱。最后，我们从整个列表中提取有趣基因的KEGG注释并计算富集分数。


```R

```
