
# The functional enrichment of data

We already know the DE gene from the array data. We can understand their role in biological function by analyzing the GO categories in the collection, so as to learn more about this group of genes at the biological level. This section is about the enrichment of the gene set with the GO term.<br>
我们已经从阵列数据中知道了DE基因，可以通过分析集合中的GO类别来了解它们在生物学功能方面的作用，从而在生物学层面更多地了解这组基因。这个部分是关于GO术语的基因集的富集。<br>

In this part, it requires the gene set data that comes from our analysis and the annotation package indicated in the `ExpressionSet` object (`hgu95av2.db`).<br>
The following steps is to describe the enrichment of genes with GO terms via a hypergeometric test.

1、Install and load the annotation database and `GOstats` library.


```R
library(hgu95av2.db)
```

    Loading required package: AnnotationDbi
    Loading required package: stats4
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
    
    Loading required package: IRanges
    Loading required package: S4Vectors
    
    Attaching package: 'S4Vectors'
    
    The following object is masked from 'package:base':
    
        expand.grid
    
    Loading required package: org.Hs.eg.db
    
    
    


```R
library(GOstats)
```

    Loading required package: Category
    Loading required package: Matrix
    
    Attaching package: 'Matrix'
    
    The following object is masked from 'package:S4Vectors':
    
        expand
    
    Loading required package: graph
    
    
    Attaching package: 'GOstats'
    
    The following object is masked from 'package:AnnotationDbi':
    
        makeGOGraph
    
    


```R
library(biomaRt)
```

2、Prepare the input data from the results of the leukemia data analysis (Working with the data of multiple classes recipe). Create two sets, one that consists of all the genes in the data and the other that consists of DE genes.<br>
根据白血病数据分析的结果准备输入数据（Working with the data of multiple classes）。 创建两个数据集，一个由数据中的所有基因组成，另一个由差异表达基因组成。


```R
# Working with the data of multiple classes 
library(leukemiasEset)
data(leukemiasEset)
pheno <- pData(leukemiasEset)
mydata <- leukemiasEset[, sampleNames(leukemiasEset)[c(1:3, 13:15, 25:27, 49:51)]]
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

    
    Attaching package: 'limma'
    
    The following object is masked from 'package:BiocGenerics':
    
        plotMA
    
    


```R
all_genes <- rownames(tested2)
head(all_genes)
```


<ol class=list-inline>
	<li>'ENSG00000152078'</li>
	<li>'ENSG00000117519'</li>
	<li>'ENSG00000145850'</li>
	<li>'ENSG00000170180'</li>
	<li>'ENSG00000087586'</li>
	<li>'ENSG00000047597'</li>
</ol>




```R
sel_genes <- rownames(DE2)
head(sel_genes)
```


<ol class=list-inline>
	<li>'ENSG00000152078'</li>
	<li>'ENSG00000117519'</li>
	<li>'ENSG00000145850'</li>
	<li>'ENSG00000170180'</li>
	<li>'ENSG00000087586'</li>
	<li>'ENSG00000047597'</li>
</ol>



3、Map these sets to their Entrez IDs.进行配对


```R
# set the mart
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

# get entrez ids for all genes
all_genes <- c(getBM(filters= "ensembl_gene_id", attributes=c("entrezgene"),
                     values= all_genes, mart= mart)) 

# get entrez ids for DE
sel_genes <- c(getBM(filters= "ensembl_gene_id", attributes=c("entrezgene"),
                     values= sel_genes, mart= mart)) 
```

4、Define a cutoff for the test statistics.为检验统计定义一个临界值（截止点）.


```R
hgCutoff <- 0.05
```

5、Use a `GOHyperGParams` object as an input parameter for the enrichment computations. <br>使用`GOHyperGParams`对象作为富集计算的输入参数。


```R
params <- new("GOHyperGParams", geneIds = sel_genes, universeGeneIds = all_genes, 
              annotation = "hgu95av2.db", ontology = "BP", pvalueCutoff = hgCutoff, 
              conditional = FALSE, testDirection = "over")
```

    Warning message in makeValidParams(.Object):
    "converting geneIds from list to atomic vector via unlist"Warning message in makeValidParams(.Object):
    "converting univ from list to atomic vector via unlist"Warning message in makeValidParams(.Object):
    "removing duplicate IDs in universeGeneIds"

6、Once having `GOHyperGParams` object, perform a hypergeometric test to get the p-value for the GO annotations.<br>
一旦拥有`GOHyperGParams`对象，执行超几何检验以获取GO注释.


```R
hgOver <- hyperGTest(params)
hgOver
```


    Gene to GO BP  test for over-representation 
    3178 GO BP ids tested (445 have p < 0.05)
    Selected gene set size: 126 
        Gene universe size: 7898 
        Annotation package: hgu95av2 


7、Check the summary of the object that was created.


```R
summary(hgOver)
```


<table>
<thead><tr><th scope=col>GOBPID</th><th scope=col>Pvalue</th><th scope=col>OddsRatio</th><th scope=col>ExpCount</th><th scope=col>Count</th><th scope=col>Size</th><th scope=col>Term</th></tr></thead>
<tbody>
	<tr><td>GO:0046501                                                          </td><td>1.182128e-09                                                        </td><td>129.483333                                                          </td><td>0.14358065                                                          </td><td> 6                                                                  </td><td>  9                                                                 </td><td>protoporphyrinogen IX metabolic process                             </td></tr>
	<tr><td>GO:0006779                                                          </td><td>1.287903e-09                                                        </td><td> 57.088235                                                          </td><td>0.23930109                                                          </td><td> 7                                                                  </td><td> 15                                                                 </td><td>porphyrin-containing compound biosynthetic process                  </td></tr>
	<tr><td>GO:0006778                                                          </td><td>1.992940e-09                                                        </td><td> 32.864407                                                          </td><td>0.38288174                                                          </td><td> 8                                                                  </td><td> 24                                                                 </td><td>porphyrin-containing compound metabolic process                     </td></tr>
	<tr><td>GO:0033014                                                          </td><td>2.259418e-09                                                        </td><td> 50.738562                                                          </td><td>0.25525449                                                          </td><td> 7                                                                  </td><td> 16                                                                 </td><td>tetrapyrrole biosynthetic process                                   </td></tr>
	<tr><td>GO:0042168                                                          </td><td>9.563327e-09                                                        </td><td> 38.039216                                                          </td><td>0.30311471                                                          </td><td> 7                                                                  </td><td> 19                                                                 </td><td>heme metabolic process                                              </td></tr>
	<tr><td>GO:0006782                                                          </td><td>1.954270e-08                                                        </td><td>160.537190                                                          </td><td>0.11167384                                                          </td><td> 5                                                                  </td><td>  7                                                                 </td><td>protoporphyrinogen IX biosynthetic process                          </td></tr>
	<tr><td>GO:0006783                                                          </td><td>2.291555e-08                                                        </td><td> 55.464286                                                          </td><td>0.20739428                                                          </td><td> 6                                                                  </td><td> 13                                                                 </td><td>heme biosynthetic process                                           </td></tr>
	<tr><td>GO:0033013                                                          </td><td>6.984665e-08                                                        </td><td> 18.750605                                                          </td><td>0.57432261                                                          </td><td> 8                                                                  </td><td> 36                                                                 </td><td>tetrapyrrole metabolic process                                      </td></tr>
	<tr><td>GO:0051188                                                          </td><td>1.968627e-06                                                        </td><td>  7.986082                                                          </td><td>1.48366675                                                          </td><td>10                                                                  </td><td> 93                                                                 </td><td>cofactor biosynthetic process                                       </td></tr>
	<tr><td>GO:0042440                                                          </td><td>6.103108e-06                                                        </td><td> 11.972136                                                          </td><td>0.71790327                                                          </td><td> 7                                                                  </td><td> 45                                                                 </td><td>pigment metabolic process                                           </td></tr>
	<tr><td>GO:0019755                                                          </td><td>7.315551e-06                                                        </td><td> 50.931148                                                          </td><td>0.14358065                                                          </td><td> 4                                                                  </td><td>  9                                                                 </td><td>one-carbon compound transport                                       </td></tr>
	<tr><td>GO:0046148                                                          </td><td>2.658675e-05                                                        </td><td> 12.093750                                                          </td><td>0.60622943                                                          </td><td> 6                                                                  </td><td> 38                                                                 </td><td>pigment biosynthetic process                                        </td></tr>
	<tr><td>GO:1905446                                                          </td><td>3.873562e-05                                                        </td><td> 94.756098                                                          </td><td>0.07976703                                                          </td><td> 3                                                                  </td><td>  5                                                                 </td><td>regulation of mitochondrial ATP synthesis coupled electron transport</td></tr>
	<tr><td>GO:0055065                                                          </td><td>5.142005e-05                                                        </td><td>  3.260141                                                          </td><td>6.31754875                                                          </td><td>18                                                                  </td><td>396                                                                 </td><td>metal ion homeostasis                                               </td></tr>
	<tr><td>GO:0055080                                                          </td><td>6.662986e-05                                                        </td><td>  3.085018                                                          </td><td>7.05140542                                                          </td><td>19                                                                  </td><td>442                                                                 </td><td>cation homeostasis                                                  </td></tr>
	<tr><td>GO:0098771                                                          </td><td>8.729942e-05                                                        </td><td>  3.017047                                                          </td><td>7.19498607                                                          </td><td>19                                                                  </td><td>451                                                                 </td><td>inorganic ion homeostasis                                           </td></tr>
	<tr><td>GO:0050801                                                          </td><td>9.153868e-05                                                        </td><td>  2.918132                                                          </td><td>7.84907572                                                          </td><td>20                                                                  </td><td>492                                                                 </td><td>ion homeostasis                                                     </td></tr>
	<tr><td>GO:0051597                                                          </td><td>1.324303e-04                                                        </td><td> 47.365854                                                          </td><td>0.11167384                                                          </td><td> 3                                                                  </td><td>  7                                                                 </td><td>response to methylmercury                                           </td></tr>
	<tr><td>GO:0055072                                                          </td><td>1.818521e-04                                                        </td><td>  8.218085                                                          </td><td>0.84553051                                                          </td><td> 6                                                                  </td><td> 53                                                                 </td><td>iron ion homeostasis                                                </td></tr>
	<tr><td>GO:0048821                                                          </td><td>1.988150e-04                                                        </td><td> 16.955191                                                          </td><td>0.30311471                                                          </td><td> 4                                                                  </td><td> 19                                                                 </td><td>erythrocyte development                                             </td></tr>
	<tr><td>GO:0032536                                                          </td><td>2.094195e-04                                                        </td><td> 37.887805                                                          </td><td>0.12762725                                                          </td><td> 3                                                                  </td><td>  8                                                                 </td><td>regulation of cell projection size                                  </td></tr>
	<tr><td>GO:0010038                                                          </td><td>2.262548e-04                                                        </td><td>  3.534441                                                          </td><td>4.11597873                                                          </td><td>13                                                                  </td><td>258                                                                 </td><td>response to metal ion                                               </td></tr>
	<tr><td>GO:0055076                                                          </td><td>2.298500e-04                                                        </td><td>  6.380282                                                          </td><td>1.24436566                                                          </td><td> 7                                                                  </td><td> 78                                                                 </td><td>transition metal ion homeostasis                                    </td></tr>
	<tr><td>GO:0030218                                                          </td><td>2.488975e-04                                                        </td><td>  6.290850                                                          </td><td>1.26031907                                                          </td><td> 7                                                                  </td><td> 79                                                                 </td><td>erythrocyte differentiation                                         </td></tr>
	<tr><td>GO:0035378                                                          </td><td>2.525232e-04                                                        </td><td>       Inf                                                          </td><td>0.03190681                                                          </td><td> 2                                                                  </td><td>  2                                                                 </td><td>carbon dioxide transmembrane transport                              </td></tr>
	<tr><td>GO:1901880                                                          </td><td>2.763075e-04                                                        </td><td>  9.994835                                                          </td><td>0.59027602                                                          </td><td> 5                                                                  </td><td> 37                                                                 </td><td>negative regulation of protein depolymerization                     </td></tr>
	<tr><td>GO:0002262                                                          </td><td>2.921525e-04                                                        </td><td>  5.254580                                                          </td><td>1.70701443                                                          </td><td> 8                                                                  </td><td>107                                                                 </td><td>myeloid cell homeostasis                                            </td></tr>
	<tr><td>GO:0051261                                                          </td><td>3.004397e-04                                                        </td><td>  7.423077                                                          </td><td>0.92529754                                                          </td><td> 6                                                                  </td><td> 58                                                                 </td><td>protein depolymerization                                            </td></tr>
	<tr><td>GO:0017085                                                          </td><td>3.104714e-04                                                        </td><td> 31.569106                                                          </td><td>0.14358065                                                          </td><td> 3                                                                  </td><td>  9                                                                 </td><td>response to insecticide                                             </td></tr>
	<tr><td>GO:0051494                                                          </td><td>3.382146e-04                                                        </td><td>  5.956656                                                          </td><td>1.32413269                                                          </td><td> 7                                                                  </td><td> 83                                                                 </td><td>negative regulation of cytoskeleton organization                    </td></tr>
	<tr><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td></tr>
	<tr><td>GO:0070495                                                          </td><td>0.04710661                                                          </td><td>31.080000                                                           </td><td>0.04786022                                                          </td><td>1                                                                   </td><td>  3                                                                 </td><td>negative regulation of thrombin-activated receptor signaling pathway</td></tr>
	<tr><td>GO:0070837                                                          </td><td>0.04710661                                                          </td><td>31.080000                                                           </td><td>0.04786022                                                          </td><td>1                                                                   </td><td>  3                                                                 </td><td>dehydroascorbic acid transport                                      </td></tr>
	<tr><td>GO:0071284                                                          </td><td>0.04710661                                                          </td><td>31.080000                                                           </td><td>0.04786022                                                          </td><td>1                                                                   </td><td>  3                                                                 </td><td>cellular response to lead ion                                       </td></tr>
	<tr><td>GO:0071288                                                          </td><td>0.04710661                                                          </td><td>31.080000                                                           </td><td>0.04786022                                                          </td><td>1                                                                   </td><td>  3                                                                 </td><td>cellular response to mercury ion                                    </td></tr>
	<tr><td>GO:0071692                                                          </td><td>0.04710661                                                          </td><td>31.080000                                                           </td><td>0.04786022                                                          </td><td>1                                                                   </td><td>  3                                                                 </td><td>protein localization to extracellular region                        </td></tr>
	<tr><td>GO:0071694                                                          </td><td>0.04710661                                                          </td><td>31.080000                                                           </td><td>0.04786022                                                          </td><td>1                                                                   </td><td>  3                                                                 </td><td>maintenance of protein location in extracellular region             </td></tr>
	<tr><td>GO:0071918                                                          </td><td>0.04710661                                                          </td><td>31.080000                                                           </td><td>0.04786022                                                          </td><td>1                                                                   </td><td>  3                                                                 </td><td>urea transmembrane transport                                        </td></tr>
	<tr><td>GO:0072236                                                          </td><td>0.04710661                                                          </td><td>31.080000                                                           </td><td>0.04786022                                                          </td><td>1                                                                   </td><td>  3                                                                 </td><td>metanephric loop of Henle development                               </td></tr>
	<tr><td>GO:0090306                                                          </td><td>0.04710661                                                          </td><td>31.080000                                                           </td><td>0.04786022                                                          </td><td>1                                                                   </td><td>  3                                                                 </td><td>spindle assembly involved in meiosis                                </td></tr>
	<tr><td>GO:0090324                                                          </td><td>0.04710661                                                          </td><td>31.080000                                                           </td><td>0.04786022                                                          </td><td>1                                                                   </td><td>  3                                                                 </td><td>negative regulation of oxidative phosphorylation                    </td></tr>
	<tr><td>GO:0098935                                                          </td><td>0.04710661                                                          </td><td>31.080000                                                           </td><td>0.04786022                                                          </td><td>1                                                                   </td><td>  3                                                                 </td><td>dendritic transport                                                 </td></tr>
	<tr><td>GO:0098937                                                          </td><td>0.04710661                                                          </td><td>31.080000                                                           </td><td>0.04786022                                                          </td><td>1                                                                   </td><td>  3                                                                 </td><td>anterograde dendritic transport                                     </td></tr>
	<tr><td>GO:1902302                                                          </td><td>0.04710661                                                          </td><td>31.080000                                                           </td><td>0.04786022                                                          </td><td>1                                                                   </td><td>  3                                                                 </td><td>regulation of potassium ion export                                  </td></tr>
	<tr><td>GO:1902603                                                          </td><td>0.04710661                                                          </td><td>31.080000                                                           </td><td>0.04786022                                                          </td><td>1                                                                   </td><td>  3                                                                 </td><td>carnitine transmembrane transport                                   </td></tr>
	<tr><td>GO:1902896                                                          </td><td>0.04710661                                                          </td><td>31.080000                                                           </td><td>0.04786022                                                          </td><td>1                                                                   </td><td>  3                                                                 </td><td>terminal web assembly                                               </td></tr>
	<tr><td>GO:1902956                                                          </td><td>0.04710661                                                          </td><td>31.080000                                                           </td><td>0.04786022                                                          </td><td>1                                                                   </td><td>  3                                                                 </td><td>regulation of mitochondrial electron transport, NADH to ubiquinone  </td></tr>
	<tr><td>GO:1903232                                                          </td><td>0.04710661                                                          </td><td>31.080000                                                           </td><td>0.04786022                                                          </td><td>1                                                                   </td><td>  3                                                                 </td><td>melanosome assembly                                                 </td></tr>
	<tr><td>GO:1903237                                                          </td><td>0.04710661                                                          </td><td>31.080000                                                           </td><td>0.04786022                                                          </td><td>1                                                                   </td><td>  3                                                                 </td><td>negative regulation of leukocyte tethering or rolling               </td></tr>
	<tr><td>GO:1904016                                                          </td><td>0.04710661                                                          </td><td>31.080000                                                           </td><td>0.04786022                                                          </td><td>1                                                                   </td><td>  3                                                                 </td><td>response to Thyroglobulin triiodothyronine                          </td></tr>
	<tr><td>GO:1904715                                                          </td><td>0.04710661                                                          </td><td>31.080000                                                           </td><td>0.04786022                                                          </td><td>1                                                                   </td><td>  3                                                                 </td><td>negative regulation of chaperone-mediated autophagy                 </td></tr>
	<tr><td>GO:1905606                                                          </td><td>0.04710661                                                          </td><td>31.080000                                                           </td><td>0.04786022                                                          </td><td>1                                                                   </td><td>  3                                                                 </td><td>regulation of presynapse assembly                                   </td></tr>
	<tr><td>GO:1905881                                                          </td><td>0.04710661                                                          </td><td>31.080000                                                           </td><td>0.04786022                                                          </td><td>1                                                                   </td><td>  3                                                                 </td><td>positive regulation of oogenesis                                    </td></tr>
	<tr><td>GO:2000821                                                          </td><td>0.04710661                                                          </td><td>31.080000                                                           </td><td>0.04786022                                                          </td><td>1                                                                   </td><td>  3                                                                 </td><td>regulation of grooming behavior                                     </td></tr>
	<tr><td>GO:0010591                                                          </td><td>0.04737911                                                          </td><td> 6.251613                                                           </td><td>0.35097493                                                          </td><td>2                                                                   </td><td> 22                                                                 </td><td>regulation of lamellipodium assembly                                </td></tr>
	<tr><td>GO:0030261                                                          </td><td>0.04737911                                                          </td><td> 6.251613                                                           </td><td>0.35097493                                                          </td><td>2                                                                   </td><td> 22                                                                 </td><td>chromosome condensation                                             </td></tr>
	<tr><td>GO:0043267                                                          </td><td>0.04737911                                                          </td><td> 6.251613                                                           </td><td>0.35097493                                                          </td><td>2                                                                   </td><td> 22                                                                 </td><td>negative regulation of potassium ion transport                      </td></tr>
	<tr><td>GO:0045841                                                          </td><td>0.04737911                                                          </td><td> 6.251613                                                           </td><td>0.35097493                                                          </td><td>2                                                                   </td><td> 22                                                                 </td><td>negative regulation of mitotic metaphase/anaphase transition        </td></tr>
	<tr><td>GO:0030104                                                          </td><td>0.04743354                                                          </td><td> 3.924797                                                           </td><td>0.81362370                                                          </td><td>3                                                                   </td><td> 51                                                                 </td><td>water homeostasis                                                   </td></tr>
	<tr><td>GO:0008016                                                          </td><td>0.04931757                                                          </td><td> 2.378750                                                           </td><td>2.64826538                                                          </td><td>6                                                                   </td><td>166                                                                 </td><td>regulation of heart contraction                                     </td></tr>
	<tr><td>GO:0070252                                                          </td><td>0.04994437                                                          </td><td> 3.037330                                                           </td><td>1.38794632                                                          </td><td>4                                                                   </td><td> 87                                                                 </td><td>actin-mediated cell contraction                                     </td></tr>
</tbody>
</table>



8、Get the number of genes associated with the different categories.<br>获取与不同类别相关的基因数量.


```R
head(geneCounts(hgOver))
head(universeCounts(hgOver))
```


<dl class=dl-horizontal>
	<dt>GO:0046501</dt>
		<dd>6</dd>
	<dt>GO:0006779</dt>
		<dd>7</dd>
	<dt>GO:0006778</dt>
		<dd>8</dd>
	<dt>GO:0033014</dt>
		<dd>7</dd>
	<dt>GO:0042168</dt>
		<dd>7</dd>
	<dt>GO:0006782</dt>
		<dd>5</dd>
</dl>




<dl class=dl-horizontal>
	<dt>GO:0046501</dt>
		<dd>9</dd>
	<dt>GO:0006779</dt>
		<dd>15</dd>
	<dt>GO:0006778</dt>
		<dd>24</dd>
	<dt>GO:0033014</dt>
		<dd>16</dd>
	<dt>GO:0042168</dt>
		<dd>19</dd>
	<dt>GO:0006782</dt>
		<dd>7</dd>
</dl>



9、Plot the GO directed acyclic graph (DAG). 绘制GO有向无环图（DAG）.<br>
The plot is too much, so I just take some data of them to make a plot.


```R
library(Rgraphviz)    # `Rgraphviz` is used to show the DAG plot of GO enrichment
```


```R
subhgOver <- subGraph(snodes=as.character(summary(hgOver)[1:9,1]), graph = goDag(hgOver))
plot(subhgOver)
```


[plot](https://github.com/Chengshu21/Chapter-5-Analyzing-Microarray-Data-with-R/blob/master/md/pic/output-The%20functional%20enrichment%20of%20data.png)


10、Generate the report as an HTML file that can be read using any browser.生成可以使用任何浏览器读取的HTML文件


```R
htmlReport(hgOver, file="D:/Try-practice/Chapter 5/ALL_hgco.html")
```
