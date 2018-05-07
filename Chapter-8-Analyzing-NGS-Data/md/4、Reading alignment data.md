
# 4、Reading alignment data

Generally, what we get from the NGS experiments are sequence reads, and we need align and map them to a reference genome. At first, to do NGS data analysis, we should align our reads to the reference genome. There are tools beyond R to do this. Most commonly used alignment tools include BWA and Bowtie.<br>

In R, we use 'Biostrings' and 'muscle' packages to do pairwise or multiple sequence alignment.<br> 

After mapping the reads in the FASTQ file to the reference genome, what we get is a sequence alignment map (SAM) or BAM (binary version of SAM) file. This part will discuss how to read/load these files in the R workspace.<br>

Here, we use the `Rsamtools` package to read the `BAM` files. The example data used is a sequence set sample (HG00107) from the 1000 Genome Project. The file `HG00107.chrom11.ILLUMINA.bwa.GBR.low_coverage.20130415.bam` is used.

一般来说，我们从NGS实验中得到的是序列片段，我们需要对齐并将它们映射到参考基因组。首先，为了进行NGS数据分析，我们应该将我们的片段序列与参考基因组进行比对。除了R之外，还有其他工具可以做到这一点，最常用的对齐工具包括BWA和Bowtie。

在R中，我们使用'Biostrings'和'muscle'包进行成对或多重序列比对。

将FASTQ文件中的片段映射到参考基因组后，我们得到的是序列比对图（SAM）或BAM（SAM的二进制版本）文件。本部分将讨论如何在R工作区中读取/加载这些文件

在这里，我们使用`Rsamtools`包来读取`BAM`文件。
使用的示例数据是来自1000 Genome Project的序列组样本（[HG00107](ftp://ftp.ncbi.nlm.nih.gov/1000genomes/ftp/data/HG00107/alignment/))。
使用文件`HG00107.chrom11.ILLUMINA.bwa.GBR.low_coverage.20130415.bam`。<br>
HG00107:GBR(population：British individuals from England and Scotland), SRS006848

1、Download an example file from within R.下载一个bam文件举例。


```R
download.file(url="ftp://ftp.ncbi.nlm.nih.gov/1000genomes/ftp/data/HG00107/alignment/HG00107.chrom11.ILLUMINA.bwa.GBR.low_coverage.20130415.bam",
              destfile = "D:/Try-practice/Chapter-8/4-read alignment file/HG00107.chrom11.ILLUMINA.bwa.GBR.low_coverage.20130415.bam")
```

2、To read the `BAM` file, use the `scanBam` function from the `Rsamtools` package.<br>
用`Rsamtools`包中的`scanBam`来读取`BAM`文件。<br>


```R
library(Rsamtools)
```


```R
bam <- scanBam("D:/Try-practice/Chapter-8/4-read alignment file/HG00107.chrom11.ILLUMINA.bwa.GBR.low_coverage.20130415.bam")
```

3、Take a look at the attributes for the first list element of the read data.<br>
查看数据列表第一列元素的属性特征。


```R
names(bam[[1]])
```


<ol class=list-inline>
	<li>'qname'</li>
	<li>'flag'</li>
	<li>'rname'</li>
	<li>'strand'</li>
	<li>'pos'</li>
	<li>'qwidth'</li>
	<li>'mapq'</li>
	<li>'cigar'</li>
	<li>'mrnm'</li>
	<li>'mpos'</li>
	<li>'isize'</li>
	<li>'seq'</li>
	<li>'qual'</li>
</ol>



4、Check the count of the records in the data.查看数据中的记录数。


```R
countBam("D:/Try-practice/Chapter-8/4-read alignment file/HG00107.chrom11.ILLUMINA.bwa.GBR.low_coverage.20130415.bam")
```


<table>
<thead><tr><th scope=col>space</th><th scope=col>start</th><th scope=col>end</th><th scope=col>width</th><th scope=col>file</th><th scope=col>records</th><th scope=col>nucleotides</th></tr></thead>
<tbody>
	<tr><td>NA                                                        </td><td>NA                                                        </td><td>NA                                                        </td><td>NA                                                        </td><td>HG00107.chrom11.ILLUMINA.bwa.GBR.low_coverage.20130415.bam</td><td>10237012                                                  </td><td>1033938212                                                </td></tr>
</tbody>
</table>



`scanBam` function together with another input function called `countBam` imports the binary alignment maps in an NGS analysis. The `what` argument defines the attributes that have to be imported from the data.

5、If want to read only selected attributes, set them as parameters.<br>
只看所需要的属性特征参数，将其设置为参数。


```R
what <- c("rname", "strand", "pos", "qwidth", "seq")
param <- ScanBamParam(what=what)
bam2 <- scanBam("D:/Try-practice/Chapter-8/4-read alignment file/HG00107.chrom11.ILLUMINA.bwa.GBR.low_coverage.20130415.bam",
                param=param)
names(bam2[[1]])
```


<ol class=list-inline>
	<li>'rname'</li>
	<li>'strand'</li>
	<li>'pos'</li>
	<li>'qwidth'</li>
	<li>'seq'</li>
</ol>



6、Read the data as a `DataFrame` object using the `do.call` function.<br>
使用do.call函数将数据作为DataFrame对象读取。


```R
bam_df <- do.call("DataFrame", bam[[1]])
head(bam_df)
```


    DataFrame with 6 rows and 13 columns
                    qname      flag    rname   strand       pos    qwidth      mapq
              <character> <integer> <factor> <factor> <integer> <integer> <integer>
    1  ERR229778.27253914        73       11        +     60009       101         0
    2  ERR229778.27253914       133       11       NA        NA        NA        NA
    3 ERR229778.105643880       163       11        +     60063       101         0
    4  ERR229778.63233022       163       11        +     60068       101         0
    5  ERR229778.76311833        99       11        +     60114       101         0
    6 ERR229778.116261665       163       11        +     60148       101         0
            cigar     mrnm      mpos     isize                     seq
      <character> <factor> <integer> <integer>          <DNAStringSet>
    1        101M       11        NA        NA CATTAGAAAA...CTTCCTCATA
    2          NA       11     60009        NA GAATTTGTTT...AACTTCTGGA
    3        101M       11     60354       394 CCAAAACAAT...ACAAAAATCA
    4       99M2S       11     60328       360 ACAATCTATC...AATCAAGGCG
    5       96M5S       11     60370       355 ACTGGGAGAC...ATTAATGCAA
    6        101M       11     60394       346 ATAGCCACAA...AATACAGATC
                         qual
               <PhredQuality>
    1 BEFEDDGHHH...GLEHHKFGEE
    2 CGHFFGFEGG...GJGKJGJDDE
    3 CFGIHIFGIG...JGIGFFEDEE
    4 @DCDCCA?CC...FCFHFHD###
    5 B8==9>GCG)...EGDH######
    6 BFFFDCGFFJ...IIFCFHFGDD


Each element in the file is an alignment. The list can be easily transformed into a `data.frame` object as shown in step 6. As a `data.frame` object, each entry in the list is a row and the attributes of the alignment are shown in terms of columns. Finally, we use the `table` function to determine how many of the sequences in the file meet the desired conditions.<br>
列表中的每个元素都是一个对齐。该列表可以很容易地转换为data.frame对象，如步骤6所示。作为data.frame对象，列表中的每个条目都是一行，并且列的属性以列的形式显示。最后，我们使用`table`函数来确定文件中有多少个序列符合所需的条件。

7、From this `DataFrame` object, extract the sequences that fulfill certain conditions.
从数据框对象中筛选满足特定条件的序列。


```R
table(bam_df$rname == '21' & bam_df$flag == 16)
```


    
       FALSE 
    10237012 



```R
subset(bam_df, bam_df$rname == '21' & bam_df$flag == 16)
```


    DataFrame with 0 rows and 13 columns



```R

```
