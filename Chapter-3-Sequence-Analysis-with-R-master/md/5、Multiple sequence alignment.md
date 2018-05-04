
# Multiple sequence alignment
Last time, what I have learned is the alignment between two sequences. Now it is the alignment among many sequences. This helps us to know how similar or dissimilar the sequences are against each other. <br>
A new package, `muscle`, is needed.<br>
I will use the protein sequences of Cytochrome oxidase (COX-2) from five different species(download from [Uniprot](http://www.uniprot.org/).`MUSCLE` algorithms that are very well established and frequently used by bioinformaticians.


```R
source("http://www.bioconductor.org/biocLite.R")
```

    Bioconductor version 3.6 (BiocInstaller 1.28.0), ?biocLite for help
    


```R
biocLite("muscle")
```

    BioC_mirror: https://bioconductor.org
    Using Bioconductor 3.6 (BiocInstaller 1.28.0), R 3.4.2 (2017-09-28).
    Installing package(s) 'muscle'
    

    package 'muscle' successfully unpacked and MD5 sums checked
    
    The downloaded binary packages are in
    	C:\Users\Master1\AppData\Local\Temp\Rtmp6Fd49p\downloaded_packages
    

    Old packages: 'GenomicFeatures', 'GenomicRanges', 'matrixStats', 'tibble',
      'digest', 'pbdZMQ', 'stringi'
    


```R
library(muscle)
library(Biostrings)
```


```R
library(seqinr)
```

    
    Attaching package: 'seqinr'
    
    The following object is masked from 'package:Biostrings':
    
        translate
    
    


```R
myseq <- readAAStringSet("D:/Try-practice/Chapter 3/fastaMSA.fasta")
myseq
myMSA <- muscle(myseq, out=NULL, quiet = FALSE)
myMSA   # muscle(myseq, out=NULL, quiet = FALSE)
```


      A AAStringSet instance of length 10
         width seq                                              names               
     [1]   227 MAHAAQVGLQDATSPIMEELITF...FMPIVLELIPLKIFEMGPVFTL sp|P00403|COX2_HUMAN
     [2]   227 MAYPMQLGFQDATSPIMEELLHF...FMPIVLELVPLKYFEKWSASML sp|P68530|COX2_BOVIN
     [3]   227 MAYPFQLGLQDATSPIMEELTNF...FMPIVLEMVPLKYFENWSASMI sp|P00406|COX2_RAT
     [4]   231 MNNFFQGYNLLFQHSLFASYMDW...MPIALEVTLLDNFKSWCFGTME sp|P24894|COX2_CAEEL
     [5]   227 MAYPFQLGFQDATSPIMEELLHF...FMPIVLELVPLKHFEEWSASML sp|P48660|COX2_HORSE
     [6]   227 MAYPLQLGFQDATSPVMEELLHF...FMPIVLELVPLKYFESWSASLA sp|Q38PR9|COX2_MAMPR
     [7]   228 MAYPLQLGFQDATSPVMEELLHF...MPIVLELVPLKYFENWSASLAQ sp|Q9TA26|COX2_LOXAF
     [8]   227 MAHAAQVGLQDATSPIMEELITF...FMPIVLELIPLKIFEMGPVFAL sp|P26456|COX2_GORGO
     [9]   227 MAYPFQLGFQDATSPIMEELLHF...FMPIVLELVPLTYFEKWSASML sp|P48890|COX2_FELCA
    [10]   227 MAYPFQLGLQDATSPIMEELMNF...FMPIVLEMVPLKYFENWSASMI sp|P00405|COX2_MOUSE


    Option -out must have value
    Invalid option "out"
    
    Basic usage
    
        muscle -in <inputfile> -out <outputfile>
    
    Common options (for a complete list please see the User Guide):
    
        -in <inputfile>    Input file in FASTA format (default stdin)
        -out <outputfile>  Output alignment in FASTA format (default stdout)
        -diags             Find diagonals (faster for similar sequences)
        -maxiters <n>      Maximum number of iterations (integer, default 16)
        -maxhours <h>      Maximum time to iterate in hours (default no limit)
        -html              Write output in HTML format (default FASTA)
        -msf               Write output in GCG MSF format (default FASTA)
        -clw               Write output in CLUSTALW format (default FASTA)
        -clwstrict         As -clw, with 'CLUSTAL W (1.81)' header
        -log[a] <logfile>  Log to file (append if -loga, overwrite if -log)
        -quiet             Do not write progress messages to the screen
        -version           Display version information and exit
    
    Without refinement (very fast, avg accuracy similar to T-Coffee): -maxiters 2
    Fastest possible (amino acids): -maxiters 1 -diags -sv -distance1 kbit20_3
    Fastest possible (nucleotides): -maxiters 1 -diags
    file239c6b6b5361 10 seqs, max length 231, avg  length 227
    1807 MB(14%)00:00:00                Iter   1  100.00%  K-mer dist pass 1
    1807 MB(14%)00:00:00                Iter   1  100.00%  K-mer dist pass 2
    1807 MB(14%)00:00:00                Iter   1  100.00%  Align node
    1807 MB(14%)00:00:00                Iter   1  100.00%  Root alignment
    1807 MB(14%)00:00:00                Iter   2  100.00%  Refine tree
    1807 MB(14%)00:00:00                Iter   2  100.00%  Root alignment
    1807 MB(14%)00:00:00                Iter   2  100.00%  Root alignment
    1807 MB(14%)00:00:00                Iter   3  100.00%  Refine biparts
    1807 MB(14%)00:00:00                Iter   4  100.00%  Refine biparts
    1807 MB(14%)00:00:00                Iter   5  100.00%  Refine biparts
    1807 MB(14%)00:00:00                Iter   6  100.00%  Refine biparts
    1807 MB(14%)00:00:00                Iter   7  100.00%  Refine biparts
    1807 MB(14%)00:00:00                Iter   8  100.00%  Refine biparts
    1807 MB(14%)00:00:00                Iter   9  100.00%  Refine biparts



    AAMultipleAlignment with 10 rows and 232 columns
          aln                                                   names               
     [1] MAHAAQ---VGLQDATSPIMEELIT...HSFMPIVLELIPLKIFEMGPVFTL- sp|P00403|COX2_HUMAN
     [2] MAYPMQ---LGFQDATSPIMEELLH...HSFMPIVLELVPLKYFEKWSASML- sp|P68530|COX2_BOVIN
     [3] MAYPFQ---LGLQDATSPIMEELTN...HSFMPIVLEMVPLKYFENWSASMI- sp|P00406|COX2_RAT
     [4] MNNFFQGYNLLFQHSLFASYMDWFH...HSFMPIALEVTLLDNFKSWCFGTME sp|P24894|COX2_CAEEL
     [5] MAYPFQ---LGFQDATSPIMEELLH...HSFMPIVLELVPLKHFEEWSASML- sp|P48660|COX2_HORSE
     [6] MAYPLQ---LGFQDATSPVMEELLH...HSFMPIVLELVPLKYFESWSASLA- sp|Q38PR9|COX2_MAMPR
     [7] MAYPLQ---LGFQDATSPVMEELLH...HSFMPIVLELVPLKYFENWSASLAQ sp|Q9TA26|COX2_LOXAF
     [8] MAHAAQ---VGLQDATSPIMEELIT...HSFMPIVLELIPLKIFEMGPVFAL- sp|P26456|COX2_GORGO
     [9] MAYPFQ---LGFQDATSPIMEELLH...HSFMPIVLELVPLTYFEKWSASML- sp|P48890|COX2_FELCA
    [10] MAYPFQ---LGLQDATSPIMEELMN...HSFMPIVLEMVPLKYFENWSASMI- sp|P00405|COX2_MOUSE



```R
library(muscle)
mySeq<- muscle::read.fasta("D:/Try-practice/Chapter 3/fastaMSA.fasta")
# Please assign the path to the file fastaMSA.fasta
MyMSA <- muscle(mySeq)
print(MyMSA, from=1, to=51)
```


    Error: 'read.fasta' is not an exported object from 'namespace:muscle'
    Traceback:
    

    1. muscle::read.fasta

    2. getExportedValue(pkg, name)

    3. stop(gettextf("'%s' is not an exported object from 'namespace:%s'", 
     .     name, getNamespaceName(ns)), call. = FALSE, domain = NA)



```R

```
