# Analyzing NGS Data
HGP(human genome project) aimed to determine the sequences that make up human DNA. It was completed in 2003, but it isn't an end. The project not only seeded the many sequencing projects, but also encouraged the development of technologies that could enable faster and economical sequencing of genomes. A single human genome is not enough to understand genetic information on health; rather, many genomes are required for such studies. Due to the demand for cheaper and faster sequencing methods, the development of NGS(Next Generation Sequencing) should be driven quickly. The NGS platforms perform massively parallel sequencing, where millions of fragments of DNA from a single sample are sequenced in parallel, facilitating high-throughput sequencing. This allows an entire genome to be sequenced faster and at a lower cost.<br>
HGP（人类基因组计划）旨在确定组成人类DNA的序列。 它在2003年完成，但它还没有结束。 该项目不仅延伸了许多测序项目，而且还鼓励开发能够对基因组进行更快更经济的测序的技术。 单一的人类基因组不足以了解健康相关的遗传信息; 相反，这样的研究需要许多基因组。 由于对便宜和更快的测序方法的需求，下一代测序（NGS）的发展迫在眉睫。 NGS平台执行大规模并行测序，其中来自单个样品的数百万个DNA片段被并行测序，促进高通量测序。 这使得整个基因组的测序速度更快，成本更低。<br>

NGS(Next Generation Sequencing) is a term used to describe a number of different modern sequencing technologies.<br>
下一代测序技术是用来描述不同测序技术的一个术语。<br>
It includes some popular technologies, explored as follows:<br>
* Illumina (Solexa) sequencing （Illumina (Solexa)测序）
* Roche 454 sequencing（罗氏454测序）
* Ion torrent (proton and PGM sequencing)
* SOLiD sequencing

Sanger sequencing is the first method to do sequencing. But it is slower, more expensive and so on. While these technologies above allow faster and cheaper sequencing of DNA and RNA. Chapter 8 will deal with the results obtained from these technologies and some related technologies. To know more on NGS technologies and their application, refer to [Metzker's Sequencing technologies-the next generation](http://www.nature.com/nrg/journal/v11/n1/full/nrg2626.html).<br>
Sanger测序相当于是第一代测序技术的方法。但他的速度较慢，成本较高。然而上述技术能对DNA和RNA进行更快速和更便宜的测序。第8章将讨论从这些技术和相关技术中获得的结果。 欲了解NGS技术及其应用的更多信息，请参阅“Metzker测序技术-下一代”。<br>

Before we do data analysis, we should have a understanding of the data. The popular formats of sequence data in NGS are the `BAM`, `FASTA`, `FASTQ` format, etc. But now, bioinformaticians are perferrable to use the `FASTQ`format. The `FASTQ` data format consists of four lines. The first line is for the sequence name, the second is for the sequence itself, the third is for optional information about the sequence, and the fourth is for a confidence or accuracy measurement of bases. <br>
在进行数据分析之前，我们应该了解所拥有的数据。 NGS中常用的序列数据格式是'BAM'，'FASTA'，'FASTQ'等格式。但是现在，生物信息学家们更喜欢用`FASTQ`格式。`FASTQ`数据格式由四行组成。第一行是序列名称，第二行是序列本身，第三行是关于序列的可选信息，第四行是碱基置信度或准确性测量。<br>

The data quality in NGS is measured in terms of a metrics called the `Phred score`. A `Phred score` is assigned to each base during the sequencing process; therefore, we have a corresponding character for each base in FASTQ data. Mathematically, the `Phred score` (Q):<br>
                              Q = -10(lgP)<br>
 `P` is the estimated error probability in base call (process of assigning bases to peaks). This establishes a logarithmic relationship between the quality and base calling error, which allows you to work with very small errors (close to zero) and deal in high accuracies numerically. Thus, 99.999 percent (1 in 100,000) accuracy in base calling yields a score equal to 50. 
The FASTQ format displays this quality measure in terms of ASCII characters (usually, the (33+Q)th ASCII character is used to represent a Q value). That's why, we see such characters throughout every fourth line of every FASTQ file. The FASTQ data allows the storage of the sequence and quality information for each read in a compact text-based format. To know more about the FASTQ format, refer to The Sanger FASTQ file format for sequences with quality scores, and [the Solexa/Illumina FASTQ variants article
by Cock and others](http://nar.oxfordjournals.org/content/38/6/1767.full).<br>
NGS中的数据质量以称为“Phred分数”的度量标准衡量。在测序过程中将“Phred评分”分配给每个碱基; 因此，我们为FASTQ数据中的每个基数都有相应的字符。“P”是基地调用中的估计错误概率（将基地分配给高峰的过程）。这建立了质量和基本调用错误之间的对数关系，它允许您使用非常小的错误（接近于零）并以数字方式处理高精度。因此，99.999％（1万分之一）的基本调用精确度得分等于50。FASTQ格式以ASCII字符显示此质量度量（通常，第（33 + Q）个ASCII字符用于表示Q值）。这就是为什么我们在每个FASTQ文件的每第四行都会看到这样的字符。 FASTQ数据允许以紧凑的基于文本的格式存储每次读取的序列和质量信息。要了解有关FASTQ格式的更多信息，请参阅Sanger FASTQ文件格式以获得具有质量得分的序列及相关文章。<br>


This Chapter introduce the following steps to do NGS analysis.<br>
* [1、Querying the SRA databases](https://github.com/Chengshu21/Chapter-8-Analyzing-NGS-Data/blob/master/MD/1、Querying%20the%20SRA%20database.md)
* [2、Downloading data from the SRA database](https://github.com/Chengshu21/Chapter-8-Analyzing-NGS-Data/blob/master/MD/2、Downloading%20data%20from%20the%20SRA%20database.md)
* [3、Reading FASTQ files in R](https://github.com/Chengshu21/Chapter-8-Analyzing-NGS-Data/blob/master/MD/3、Reading%20FASTQ%20files%20in%20R.md#3reading-fastq-files-in-r)
* [4、Reading alignment data](https://github.com/Chengshu21/Chapter-8-Analyzing-NGS-Data/blob/master/MD/4、Reading%20alignment%20data.md)
* [5、Preprocessing the raw NGS data](https://github.com/Chengshu21/Chapter-8-Analyzing-NGS-Data/blob/master/MD/5、Preprocessing%20the%20raw%20NGS%20data.md)
* [6、Analyzing RNAseq data with the edgeR package](https://github.com/Chengshu21/Chapter-8-Analyzing-NGS-Data/blob/master/MD/6、Analyzing%20RNAseq%20data%20with%20the%20edgeR%20package.md)
* [7、The differential analysis of NGS data using limma](https://github.com/Chengshu21/Chapter-8-Analyzing-NGS-Data/blob/master/MD/7、The%20differential%20analysis%20of%20NGS%20data%20using%20limma.md)
* [8、Enriching RNAseq data with GO terms](https://github.com/Chengshu21/Chapter-8-Analyzing-NGS-Data/blob/master/MD/8、Enriching%20RNAseq%20data%20with%20GO%20terms.md)
* [9、The KEGG enrichment of sequence data](https://github.com/Chengshu21/Chapter-8-Analyzing-NGS-Data/blob/master/MD/9、The%20KEGG%20enrichment%20of%20sequence%20data.md)
* [10、Analyzing methylation data](https://github.com/Chengshu21/Chapter-8-Analyzing-NGS-Data/blob/master/MD/10、Analyzing%20methylation%20data.md)
* [11、Analyzing ChipSeq data](https://github.com/Chengshu21/Chapter-8-Analyzing-NGS-Data/blob/master/MD/11、Analyzing%20ChipSeq%20data.md)
