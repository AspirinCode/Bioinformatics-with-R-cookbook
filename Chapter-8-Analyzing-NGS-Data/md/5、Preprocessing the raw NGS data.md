
# 5、Preprocessing the raw NGS data

`FASTQ` data has the sequences (the bases) as the corresponding quality scores (Phred) in terms of ASCII characters, as explained in the introductory part of the chapter. Once read into the R workspace, the data is ready to be analyzed. However, it needs some preprocessing to meet the desired conditions on quality and data instance according to our interest. For example, we need higher Phred scores and a particular strand. This preprocessing involves quality assessment and filtering. This part will deal with these aspects, specifically filtering and quality checks.<br>

we will use the data downloaded from [chipseqBioc2012](http://biocluster.ucr.edu/~tgirke/HTML_Presentations/Manuals/Rngsapps/chipseqBioc2012/data.zip). We will also continue to use the `ShortRead` library.<br>

FASTQ数据具有序列（基数）作为相应的质量分数（Phred），以ASCII字符表示。一旦读入R工作区，就可以分析数据了。但是，有时候需要对他进行一些预处理来满足质量和数据实例的所需条件。例如，我们需要更高的Phred分数和一个特定的链。这个预处理涉及质量评估和过滤。 这部分将处理这些方面，特别是过滤和质量检查方面。

1、Download the required files, and unzip the downloaded files from within R.


```R
download.file(url = "http://biocluster.ucr.edu/~tgirke/HTML_Presentations/Manuals/Rngsapps/chipseqBioc2012/data.zip",
                destfile = "data.zip")
unzip("data.zip")           # Unzips and removes the original bunzip file
install.packages("R.utils") # install the R.utils from CRAN
library(R.utils)
```

2、To assess the quality, use the `FastQuality` function from the `ShortRead` library.


```R
library(ShortRead)
myFiles <- list.files("D:/Try-practice/Chapter-8/data/fastq", "fastq", full=TRUE)
myFQ <- lapply(myFiles, readFastq)
myQual <- FastqQuality(quality(quality(myFQ[[1]])))
```

3、Convert the quality measure in terms of a matrix.按照矩阵转换质量度量。


```R
readM <- as(myQual, "matrix")
```

4、Visualize the results as a boxplot.


```R
boxplot(as.data.frame(readM), outline = FALSE, main="Per Cycle Read Quality", xlab="Cycle", ylab="Phred Quality")
# the x axis shows the cycles and the y axis represents the Phred quality score available within the FASTQ file.
```


[output](https://github.com/Chengshu21/Chapter-8-Analyzing-NGS-Data/blob/master/MD/pic/output_Preprocessing%20the%20raw%20NGS%20data.png)


5、Another interesting preprocessing step involves filtering the sequences while reading alignment data. For this, you first need alignment data. And check the object.另一预处理步骤涉及在读取对齐数据的同时过滤序列。为此，首先需要对齐数据。然后检查对象。


```R
myData <- readAligned("D:/Try-practice/Chapter-8/myBowtie.txt", type = "Bowtie")
myData
```


    class: AlignedRead
    length: 1000 reads; width: 35 cycles
    chromosome: chr5 chr10 ... chr15 chr16 
    position: 151311502 35505989 ... 25552487 33204792 
    strand: - - ... + + 
    alignQuality: NumericQuality 
    alignData varLabels: similar mismatch 


6、There are different filters that can be used, such as chromosomes, sequence length, strands, and so on. To use these, first define the required filter (in this case, define a filter for the + strand). And use the filter.<br>
可以使用不同的筛选方式，例如染色体，序列长度，链等等。要使用这些，首先定义所需的过滤器（在这种情况下，为+链定义一个过滤器）并使用。


```R
strand <- strandFilter("+")
myRead_strand <- readAligned("D:/Try-practice/Chapter-8/myBowtie.txt", filter=strand, type = "Bowtie")
```

7、Combine more than one filters and then use with the `compose` function. And then take a look at the filtered data (`myRead_filt`).


```R
chromosome <- chromosomeFilter("3")
myFilt <- compose(strand, chromosome)
myRead_filt <- readAligned("D:/Try-practice/Chapter-8/myBowtie.txt", filter=myFilt, type = "Bowtie")
myRead_filt
```


    class: AlignedRead
    length: 55 reads; width: 35 cycles
    chromosome: chr3 chr13 ... chr3 chr13 
    position: 149312235 27275680 ... 136314155 52437523 
    strand: + + ... + + 
    alignQuality: NumericQuality 
    alignData varLabels: similar mismatch 


The quality check function described in this recipe uses the Phred scores assigned in terms of ASCII characters in the FASTQ files, and computes the actual quality score for each cycle. This is reflected in the plot we saw in the previous boxplot. The filtering function simply checks the data and extracts the alignment reads that meet the filtering conditions. The other filters that
can be used include idFilter, positionFilter, and so on. We can see the difference
between the myData object that has 1,000 reads with both the + and – strands, whereas
the filtered object, myRead_filt, has only 55 reads with the + strand.
There's more...
