# Analyzing GWAS data
## 全基因组关联分析

Genome-wide association studies (`GWAS`) is a method of identifying susceptibility loci for complex diseases, which is based on the technique of scanning genomes of many subjects in order to identify the genetic variation possibly responsible for a disease through statistical tests. To identify the variation reliably, a large number of subjects is needed. Due to the high-density microarray technologies for Single Nucleotide Polymorphisms (`SNP`) detection and the [`HapMap`](http://hapmap.ncbi.nlm.nih.gov) project together with the Human Genome Project (`HGP`), this scanning technique has been made possible. And then through analysing GWAS data, we can make some biological inference.<br>
   
全基因组关联研究（`GWAS`）是一种鉴定复杂疾病易感基因座的方法，该方法基于扫描许多受试者基因组的技术，以便通过统计检验来鉴定可能导致疾病的遗传变异。 为了使识别变体的结果具有可靠性，需要大量的受试对象。 由于单核苷酸多态性（SNP）检测的高密度微阵列技术和[`HapMap`](http://hapmap.ncbi.nlm.nih.gov)项目与人类基因组计划（`HGP `），这种扫描技术已成为可能。分析GWAS数据可以进行生物学推断。<br>
   
Here are some resources:<br>
* A number of SNP studies are available on [SNPedia](http://snpedia.com/index.php/SNPedia), which shares information about the effects of variations in DNA based on peer-reviewed literature. 
* A popular database for SNPs is `dbSNP`, which is available on the [Entrez page](http://www.ncbi.nlm.nih.gov/SNP/). 
* [Online Mendelian Inheritence in Man (OMIM)](http://omim.org) is another  source of information on human genes and genetic phenotypes.<br>

GWAS uses the control and the disease groups(for example,the former without disease or without the desired trait, while the latter with desired traits). SNPs are detected from the subjects. They refer to a variation in the DNA sequence where a single nucleotide differs between two phenotypes or chromosomes. To illustrate, in a sequence `ATCGTACG`, the variant can be `ATCGCACG`, where `T` has been changed to `C`. To find the association of these SNPs with certain phenotypes, the presence of such SNPs is statistically tested with samples from these phenotypes. For example, the SNPs that are significantly more frequent in people who are suffering from the disease, compared to those who aren't, are said to be associated with the disease. We can illustrate this with an example of a disease-control GWA study. For control cases and disease cases, the allele count of SNPs is computed and then a statistical test, such as a Chi-squared test, is performed to identify variants associated with the disease under consideration. This involves the analysis of a huge amount of genotypic data. The GWAS data can be the size of hundreds of MBs up to GBs.<br>
    
GWAS使用对照组和实验组（例如，对照组没有疾病或没有期望的性状，而实验组具有期望的性状）。从受试者中检测到SNP。它们是指在两个表型或染色体之间单个核苷酸不同的DNA序列中的变异。为了说明，在序列ATCGTACG中，变体可以是ATCGCACG，其中T已经变为C.为了找到这些SNP与某些表型的关联，用来自这些表型的样品统计测试这些SNP的存在。例如，患有该疾病的人比那些不患病的人显着更频繁的SNP据说与该疾病有关。我们可以用疾病控制GWA研究的一个例子来说明这一点。对于对照病例和疾病病例，计算SNP的等位基因计数，然后进行统计检验，如卡方检验，以鉴定与所考虑的疾病相关的变体。这涉及对大量基因型数据的分析。GWAS数据可以是数百MB到GB的大小。<br>
    
GWAS data in any format is obtained from genome-wide studies. It includes information about the variants as well as their mappings onto chromosomes for every individual. This information is structured in different ways depending on the data format, such as tables or flat files. Some of these data formats will be discussed during the course of this chapter. To get a detailed idea of GWAS, refer to the Genome wide association studies article by Bush and Moore at http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1002822.<br>
   
任何格式的GWAS数据都是从全基因组研究中获得的。它包括变体的信息以及它们在每个个体染色体上的映射。由于数据格式不同（例如表格或平面文件），此信息的结构也有所不同。想详细了解GWAS，请参阅Bush和Moore撰写的Genome wide association studies文章，网址为http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1002822。<br>

GWAS是在某一特定人群中研究遗传突变和表型之间的相关性。GWAS的理论基础是连锁不平衡定律（linkage disequilibrium， LD），既假设观察到的SNP与真正的致病突变（causal variant）之间存在很强的LD。基于基因芯片设计的GWAS目前着重关注人群中的常见变异（common SNPs, 通常指最小等位基因频率MAF > 0.01），因此，通过GWAS发现的疾病易感位点主要集中在常见变异上，也即是通常所说的common disease，common variants理论。GWAS的概念是在群体遗传学的概念下发展起来的，其统计效应受到样本量的直接影响。

### 参考文章：<br>
1、Klein RJ, Zeiss C, Chew EY,ect. Complement factor H polymorphism in age-related macular degeneration. Science, 2005, 308(5720): 385−389.<br>
2、Samani NJ, Erdmann J, Hall AS, ect. Genomewide associationanalysis of coronary artery disease. N Engl J Med, 2007,357(5): 443−453.<br>
3、Herbert A, Gerry NP, McQueen MB, ect. A common geneticvariant is associated with adult and childhood obesity.Science, 2006, 312(5771): 279−283.<br>
4、Rosskopf D, Bornhorst A, Rimmbach C, ect. Comment on "A common genetic variant is associatedwith adult and childhood obesity". Science, 2007,315(5809): 187: author reply 187.<br>
5、Samani NJ, Erdmann J, Hall AS, Hengstenberg C,ect. Genomewide association analysis of coronary artery disease. N Engl J Med, 2007, 357(5): 443−453.<br>

## Steps:<br>
In  this chapter, try the following steps to perform GWAS data analysis:
### [1、The SNP association analysis](https://github.com/Chengshu21/Chapter-6-Analyzing-GWAS-Data/blob/master/md/1、The%20SNP%20association%20analysis.md)
### [2、Running association scans for SNPs](https://github.com/Chengshu21/Chapter-6-Analyzing-GWAS-Data/blob/master/md/2、Running%20association%20scans%20for%20SNPs.md)
### [3、The whole genome SNP association analysis](https://github.com/Chengshu21/Chapter-6-Analyzing-GWAS-Data/blob/master/md/3%E3%80%81The%20whole%20genome%20SNP%20association%20analysis.md)
