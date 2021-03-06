
# Reading and writing the FASTA file

The `FASTA` format is a simple and widely used format for storing biological (DNA or protein) sequences.<br> 
It begins with a single-line description that starts with a `>` symbol. The description consists of virtually anything regarding the sequence but usually carries the sequence name, ID, name of the species, name of the author, and so on. The line that follows carries the sequences (nucleotide or protein).<br>
The following is an example of a FASTA file for a protein sequence taken from Protein Data Bank (PDB) or (ID- 2BQ0):<br>
`>2BQ0:A|PDBID|CHAIN|SEQUENCE
MTMDKSELVQKAKLAEQAERYDDMAAAMKAVTEQGHELSNEERNLLSVAYKNVVGARRSSWRVISSIEQK
TERNEKKQQMGKEYREKIEAELQDICNDVLELLDKYLIPNATQPESKVFYLKMKGDYFRYLSEVASGDNK
QTTVSNQQAYQEAFEISKKEMQPTHPIRLGLALNFSVFYYEILNSPEKACSLAKTAFDEAIAELDTLNEE
SYKDSTLIMQLLRDNLTWTSENQGDEGENLYFQ`<br>
`FASTA` is a standard format for storing sequence data in databases such as `GenBank`. It is better to learn how to import and export this file format into R. 

copy and paste a FASTA sequence into a text editor and save it as a `fasta` file (even `.txt` should serve the purpose). Besides this, the `seqinr` library will be required, which we discussed in the previous `Retrieving a sequence in R`.


```R
library('seqinr')
```


```R
mysequence <- read.fasta(file = "D:/Try-practice/Chapter 3/myfasta.fasta")
mysequence
# read a 'fasta' file
```


    $`sp|P68871|HBB_HUMAN`
      [1] "m" "v" "h" "l" "t" "p" "e" "e" "k" "s" "a" "v" "t" "a" "l" "w" "g" "k"
     [19] "v" "n" "v" "d" "e" "v" "g" "g" "e" "a" "l" "g" "r" "l" "l" "v" "v" "y"
     [37] "p" "w" "t" "q" "r" "f" "f" "e" "s" "f" "g" "d" "l" "s" "t" "p" "d" "a"
     [55] "v" "m" "g" "n" "p" "k" "v" "k" "a" "h" "g" "k" "k" "v" "l" "g" "a" "f"
     [73] "s" "d" "g" "l" "a" "h" "l" "d" "n" "l" "k" "g" "t" "f" "a" "t" "l" "s"
     [91] "e" "l" "h" "c" "d" "k" "l" "h" "v" "d" "p" "e" "n" "f" "r" "l" "l" "g"
    [109] "n" "v" "l" "v" "c" "v" "l" "a" "h" "h" "f" "g" "k" "e" "f" "t" "p" "p"
    [127] "v" "q" "a" "a" "y" "q" "k" "v" "v" "a" "g" "v" "a" "n" "a" "l" "a" "h"
    [145] "k" "y" "h"
    attr(,"name")
    [1] "sp|P68871|HBB_HUMAN"
    attr(,"Annot")
    [1] ">sp|P68871|HBB_HUMAN"
    attr(,"class")
    [1] "SeqFastadna"
    
    $`sp|P68873|HBB_PANTR`
      [1] "m" "v" "h" "l" "t" "p" "e" "e" "k" "s" "a" "v" "t" "a" "l" "w" "g" "k"
     [19] "v" "n" "v" "d" "e" "v" "g" "g" "e" "a" "l" "g" "r" "l" "l" "v" "v" "y"
     [37] "p" "w" "t" "q" "r" "f" "f" "e" "s" "f" "g" "d" "l" "s" "t" "p" "d" "a"
     [55] "v" "m" "g" "n" "p" "k" "v" "k" "a" "h" "g" "k" "k" "v" "l" "g" "a" "f"
     [73] "s" "d" "g" "l" "a" "h" "l" "d" "n" "l" "k" "g" "t" "f" "a" "t" "l" "s"
     [91] "e" "l" "h" "c" "d" "k" "l" "h" "v" "d" "p" "e" "n" "f" "r" "l" "l" "g"
    [109] "n" "v" "l" "v" "c" "v" "l" "a" "h" "h" "f" "g" "k" "e" "f" "t" "p" "p"
    [127] "v" "q" "a" "a" "y" "q" "k" "v" "v" "a" "g" "v" "a" "n" "a" "l" "a" "h"
    [145] "k" "y" "h"
    attr(,"name")
    [1] "sp|P68873|HBB_PANTR"
    attr(,"Annot")
    [1] ">sp|P68873|HBB_PANTR"
    attr(,"class")
    [1] "SeqFastadna"
    
    $`sp|P02112|HBB_CHICK`
      [1] "m" "v" "h" "w" "t" "a" "e" "e" "k" "q" "l" "i" "t" "g" "l" "w" "g" "k"
     [19] "v" "n" "v" "a" "e" "c" "g" "a" "e" "a" "l" "a" "r" "l" "l" "i" "v" "y"
     [37] "p" "w" "t" "q" "r" "f" "f" "a" "s" "f" "g" "n" "l" "s" "s" "p" "t" "a"
     [55] "i" "l" "g" "n" "p" "m" "v" "r" "a" "h" "g" "k" "k" "v" "l" "t" "s" "f"
     [73] "g" "d" "a" "v" "k" "n" "l" "d" "n" "i" "k" "n" "t" "f" "s" "q" "l" "s"
     [91] "e" "l" "h" "c" "d" "k" "l" "h" "v" "d" "p" "e" "n" "f" "r" "l" "l" "g"
    [109] "d" "i" "l" "i" "i" "v" "l" "a" "a" "h" "f" "s" "k" "d" "f" "t" "p" "e"
    [127] "c" "q" "a" "a" "w" "q" "k" "l" "v" "r" "v" "v" "a" "h" "a" "l" "a" "r"
    [145] "k" "y" "h"
    attr(,"name")
    [1] "sp|P02112|HBB_CHICK"
    attr(,"Annot")
    [1] ">sp|P02112|HBB_CHICK"
    attr(,"class")
    [1] "SeqFastadna"
    
    $`sp|Q90486|HBB1_DANRE`
      [1] "m" "v" "e" "w" "t" "d" "a" "e" "r" "t" "a" "i" "l" "g" "l" "w" "g" "k"
     [19] "l" "n" "i" "d" "e" "i" "g" "p" "q" "a" "l" "s" "r" "c" "l" "i" "v" "y"
     [37] "p" "w" "t" "q" "r" "y" "f" "a" "t" "f" "g" "n" "l" "s" "s" "p" "a" "a"
     [55] "i" "m" "g" "n" "p" "k" "v" "a" "a" "h" "g" "r" "t" "v" "m" "g" "g" "l"
     [73] "e" "r" "a" "i" "k" "n" "m" "d" "n" "v" "k" "n" "t" "y" "a" "a" "l" "s"
     [91] "v" "m" "h" "s" "e" "k" "l" "h" "v" "d" "p" "d" "n" "f" "r" "l" "l" "a"
    [109] "d" "c" "i" "t" "v" "c" "a" "a" "m" "k" "f" "g" "q" "a" "g" "f" "n" "a"
    [127] "d" "v" "q" "e" "a" "w" "q" "k" "f" "l" "a" "v" "v" "v" "s" "a" "l" "c"
    [145] "r" "q" "y" "h"
    attr(,"name")
    [1] "sp|Q90486|HBB1_DANRE"
    attr(,"Annot")
    [1] ">sp|Q90486|HBB1_DANRE"
    attr(,"class")
    [1] "SeqFastadna"
    
    $`sp|P02062|HBB_HORSE`
      [1] "v" "q" "l" "s" "g" "e" "e" "k" "a" "a" "v" "l" "a" "l" "w" "d" "k" "v"
     [19] "n" "e" "e" "e" "v" "g" "g" "e" "a" "l" "g" "r" "l" "l" "v" "v" "y" "p"
     [37] "w" "t" "q" "r" "f" "f" "d" "s" "f" "g" "d" "l" "s" "n" "p" "g" "a" "v"
     [55] "m" "g" "n" "p" "k" "v" "k" "a" "h" "g" "k" "k" "v" "l" "h" "s" "f" "g"
     [73] "e" "g" "v" "h" "h" "l" "d" "n" "l" "k" "g" "t" "f" "a" "a" "l" "s" "e"
     [91] "l" "h" "c" "d" "k" "l" "h" "v" "d" "p" "e" "n" "f" "r" "l" "l" "g" "n"
    [109] "v" "l" "v" "v" "v" "l" "a" "r" "h" "f" "g" "k" "d" "f" "t" "p" "e" "l"
    [127] "q" "a" "s" "y" "q" "k" "v" "v" "a" "g" "v" "a" "n" "a" "l" "a" "h" "k"
    [145] "y" "h"
    attr(,"name")
    [1] "sp|P02062|HBB_HORSE"
    attr(,"Annot")
    [1] ">sp|P02062|HBB_HORSE"
    attr(,"class")
    [1] "SeqFastadna"
    



```R
choosebank("genbank")
```


```R
q <- query(listname = "BRCA1", query="SP=homo sapiens AND K=BRCA1")
q
# the search keyword 'K' and the species we are interested in as 'SP'
```


    176 SQ for SP=homo sapiens AND K=BRCA1



```R
# To look for the sequence identifiers, extract the name attribute of the sequences
myseqs <- getSequence(q)
head(myseqs)
```


<ol>
	<li><ol class=list-inline>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
</ol>
</li>
	<li><ol class=list-inline>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
</ol>
</li>
	<li><ol class=list-inline>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
</ol>
</li>
	<li><ol class=list-inline>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
</ol>
</li>
	<li><ol class=list-inline>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
</ol>
</li>
	<li><ol class=list-inline>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
</ol>
</li>
</ol>




```R
mynames <- getName(q)
mynames
```


<ol class=list-inline>
	<li>'AB621825.BRCA1'</li>
	<li>'AF005068.BRCA1'</li>
	<li>'AF284812.BRCA1'</li>
	<li>'AF507076.BRCA1'</li>
	<li>'AF507077.BRCA1'</li>
	<li>'AF507078.BRCA1'</li>
	<li>'AY093484.BRCA1'</li>
	<li>'AY093485.BRCA1'</li>
	<li>'AY093486.BRCA1'</li>
	<li>'AY093487.BRCA1'</li>
	<li>'AY093489.BRCA1'</li>
	<li>'AY093490.BRCA1'</li>
	<li>'AY093491.BRCA1'</li>
	<li>'AY093492.BRCA1'</li>
	<li>'AY093493.BRCA1'</li>
	<li>'AY144588'</li>
	<li>'AY150865'</li>
	<li>'AY273801.BRCA1'</li>
	<li>'AY304547.BRCA1'</li>
	<li>'AY438030'</li>
	<li>'AY438031'</li>
	<li>'AY706911'</li>
	<li>'AY706912.BRCA1'</li>
	<li>'AY706913'</li>
	<li>'AY751490'</li>
	<li>'BC030969.BRCA1'</li>
	<li>'BC062429'</li>
	<li>'BC072418.BRCA1'</li>
	<li>'BC085615.BRCA1'</li>
	<li>'BC106745.BRCA1'</li>
	<li>'BC106746'</li>
	<li>'BC115037.BRCA1'</li>
	<li>'DQ075361'</li>
	<li>'DQ116737'</li>
	<li>'DQ145822'</li>
	<li>'DQ145823'</li>
	<li>'DQ145824'</li>
	<li>'DQ145825.BRCA1'</li>
	<li>'DQ145826'</li>
	<li>'DQ190450.BRCA1'</li>
	<li>'DQ190451.BRCA1'</li>
	<li>'DQ190452.BRCA1'</li>
	<li>'DQ190453.BRCA1'</li>
	<li>'DQ190454.BRCA1'</li>
	<li>'DQ190455.BRCA1'</li>
	<li>'DQ190456.BRCA1'</li>
	<li>'DQ190457.BRCA1'</li>
	<li>'DQ299305'</li>
	<li>'DQ299306'</li>
	<li>'DQ299307'</li>
	<li>'DQ299308'</li>
	<li>'DQ299309'</li>
	<li>'DQ299310'</li>
	<li>'DQ299311'</li>
	<li>'DQ299312'</li>
	<li>'DQ299313'</li>
	<li>'DQ299314'</li>
	<li>'DQ299315'</li>
	<li>'DQ299316'</li>
	<li>'DQ299317'</li>
	<li>'DQ299318'</li>
	<li>'DQ299319'</li>
	<li>'DQ299320'</li>
	<li>'DQ299321'</li>
	<li>'DQ299322'</li>
	<li>'DQ299323'</li>
	<li>'DQ299324'</li>
	<li>'DQ299325'</li>
	<li>'DQ299326'</li>
	<li>'DQ299327'</li>
	<li>'DQ299328'</li>
	<li>'DQ299330'</li>
	<li>'DQ299331'</li>
	<li>'DQ363751.BRCA1'</li>
	<li>'DQ478408.BRCA1'</li>
	<li>'FJ940752.BRCA1'</li>
	<li>'HE600032'</li>
	<li>'HE600033.BRCA1'</li>
	<li>'HE600034'</li>
	<li>'HE600035'</li>
	<li>'HE600036'</li>
	<li>'HE600037'</li>
	<li>'HE600038'</li>
	<li>'HSU14680.BRCA1'</li>
	<li>'HSU18009'</li>
	<li>'HSU18018'</li>
	<li>'HSU37574.BRCA1'</li>
	<li>'HSU61268.BRCA1'</li>
	<li>'HSU64805'</li>
	<li>'HSU68041.BRCA1'</li>
	<li>'JN384124.BRCA1'</li>
	<li>'JN686490'</li>
	<li>'JX480460.BRCA1'</li>
	<li>'JX480461'</li>
	<li>'JX480462.BRCA1'</li>
	<li>'JX480463.BRCA1'</li>
	<li>'JX480464.BRCA1'</li>
	<li>'JX480465.BRCA1'</li>
	<li>'JX480466.BRCA1'</li>
	<li>'JX480467'</li>
	<li>'KJ625149.BRCA1'</li>
	<li>'KJ625150.BRCA1'</li>
	<li>'KJ625151.BRCA1'</li>
	<li>'KJ625152.BRCA1'</li>
	<li>'KJ625153.BRCA1'</li>
	<li>'KJ625154.BRCA1'</li>
	<li>'KJ625155.BRCA1'</li>
	<li>'KJ625156.BRCA1'</li>
	<li>'KJ625157.BRCA1'</li>
	<li>'KJ625158.BRCA1'</li>
	<li>'KJ625159.BRCA1'</li>
	<li>'KJ625160.BRCA1'</li>
	<li>'KJ625161.BRCA1'</li>
	<li>'KJ625162'</li>
	<li>'KJ625163.BRCA1'</li>
	<li>'KJ625164.BRCA1'</li>
	<li>'KJ625165.BRCA1'</li>
	<li>'KJ625166.BRCA1'</li>
	<li>'KJ625167.BRCA1'</li>
	<li>'KJ625168.BRCA1'</li>
	<li>'KJ625169.BRCA1'</li>
	<li>'KJ625170.BRCA1'</li>
	<li>'KJ625171.BRCA1'</li>
	<li>'KJ625172.BRCA1'</li>
	<li>'KJ625173.BRCA1'</li>
	<li>'KJ625174.BRCA1'</li>
	<li>'KJ625175.BRCA1'</li>
	<li>'KJ625176.BRCA1'</li>
	<li>'KJ625176.PE2'</li>
	<li>'KJ625177.BRCA1'</li>
	<li>'KJ625178.BRCA1'</li>
	<li>'KJ625179.BRCA1'</li>
	<li>'KM434065'</li>
	<li>'KP255396.BRCA1'</li>
	<li>'KP255397.BRCA1'</li>
	<li>'KP255398.BRCA1'</li>
	<li>'KP255399.BRCA1'</li>
	<li>'KP255400.BRCA1'</li>
	<li>'KP255401.BRCA1'</li>
	<li>'KP255402.BRCA1'</li>
	<li>'KP255403.BRCA1'</li>
	<li>'KP272102.BRCA1'</li>
	<li>'KP272103.BRCA1'</li>
	<li>'KP272104.BRCA1'</li>
	<li>'KP272105.BRCA1'</li>
	<li>'KP272106.BRCA1'</li>
	<li>'KP404097'</li>
	<li>'KP455327.BRCA1'</li>
	<li>'KP701015'</li>
	<li>'KP701016.BRCA1'</li>
	<li>'KP729136'</li>
	<li>'KP729137'</li>
	<li>'KP744861'</li>
	<li>'KP753383.BRCA1'</li>
	<li>'KT120061.BRCA1'</li>
	<li>'KT152888'</li>
	<li>'KT152889.BRCA1'</li>
	<li>'KT152890.BRCA1'</li>
	<li>'KT844468.BRCA1'</li>
	<li>'KT844469.BRCA1'</li>
	<li>'KU359055'</li>
	<li>'KU359056'</li>
	<li>'KU359057'</li>
	<li>'KU359058'</li>
	<li>'KU359059'</li>
	<li>'KU359060'</li>
	<li>'KU359061'</li>
	<li>'KU359062'</li>
	<li>'KU359063'</li>
	<li>'KU359064'</li>
	<li>'KX580312.BRCA1'</li>
	<li>'KX944478.BRCA1'</li>
	<li>'L78833.BRCA1'</li>
	<li>'MG969402'</li>
	<li>'S78558.BRCA1'</li>
	<li>'Y08757.BRCA1'</li>
</ol>




```R
length(mynames)
# Check the number of hits by simply checking the length of the mynames object that contains the names of the sequences found in the hits as a vector
```


176



```R
write.fasta(myseqs, mynames, file.out = "D:/Try-practice/Chapter 3/myBRCA.fasta")
# writing the retrieved sequence data carrying information in different fields as multiple sequences in a 'FASTA' file named 'MyBRCA.fasta'
```


```R
mybrca <- read.fasta(file = "D:/Try-practice/Chapter 3/myBRCA.fasta")
mybrca
# read the 'myBRCA.fasta' 
```


    $AB621825.BRCA1
     [1] "a" "t" "g" "g" "a" "t" "t" "t" "a" "t" "c" "t" "g" "c" "t" "c" "t" "t" "c"
    [20] "g" "c" "g" "t" "t" "g" "a" "a" "g" "a" "a" "g" "t" "a" "c" "a" "a" "a" "a"
    [39] "t" "g" "t" "c" "a" "t" "t" "a" "a" "t" "g" "c" "t" "a" "t" "g" "c" "a" "g"
    [58] "a" "a" "a" "a" "t" "c" "t" "t" "a" "g" "a" "g" "t" "g"
    attr(,"name")
    [1] "AB621825.BRCA1"
    attr(,"Annot")
    [1] ">AB621825.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $AF005068.BRCA1
       [1] "a" "t" "g" "a" "g" "c" "c" "t" "a" "c" "a" "a" "g" "a" "a" "a" "g" "t"
      [19] "a" "c" "g" "a" "g" "a" "t" "t" "c" "a" "g" "t" "c" "a" "a" "c" "t" "t"
      [37] "g" "t" "t" "g" "a" "a" "g" "a" "g" "c" "t" "a" "t" "t" "g" "a" "a" "a"
      [55] "a" "c" "c" "a" "t" "t" "t" "g" "t" "g" "c" "t" "t" "t" "t" "c" "a" "g"
      [73] "c" "t" "t" "g" "a" "c" "a" "c" "a" "g" "g" "t" "t" "t" "g" "g" "a" "g"
      [91] "t" "a" "t" "g" "c" "a" "a" "a" "c" "a" "g" "c" "t" "a" "t" "a" "a" "t"
     [109] "t" "t" "t" "g" "c" "a" "a" "a" "a" "a" "a" "g" "g" "a" "a" "a" "a" "t"
     [127] "a" "a" "c" "t" "c" "t" "c" "c" "t" "g" "a" "a" "c" "a" "t" "c" "t" "a"
     [145] "a" "a" "a" "g" "a" "t" "g" "a" "a" "g" "t" "t" "t" "c" "t" "a" "t" "c"
     [163] "a" "t" "c" "c" "a" "a" "a" "g" "t" "a" "t" "g" "g" "g" "c" "t" "a" "c"
     [181] "a" "g" "a" "a" "a" "c" "c" "g" "t" "g" "c" "c" "a" "a" "a" "a" "g" "a"
     [199] "c" "t" "t" "c" "t" "a" "c" "a" "g" "a" "g" "t" "g" "a" "a" "c" "c" "c"
     [217] "g" "a" "a" "a" "a" "t" "c" "c" "t" "t" "c" "c" "t" "t" "g" "g" "a" "a"
     [235] "a" "c" "c" "a" "g" "t" "c" "t" "c" "a" "g" "t" "g" "t" "c" "c" "a" "a"
     [253] "c" "t" "c" "t" "c" "t" "a" "a" "c" "c" "t" "t" "g" "g" "a" "a" "c" "t"
     [271] "g" "t" "g" "a" "g" "a" "a" "c" "t" "c" "t" "g" "a" "g" "g" "a" "c" "a"
     [289] "a" "a" "g" "c" "a" "g" "c" "g" "g" "a" "t" "a" "c" "a" "a" "c" "c" "t"
     [307] "c" "a" "a" "a" "a" "g" "a" "c" "g" "t" "c" "t" "g" "t" "c" "t" "a" "c"
     [325] "a" "t" "t" "g" "a" "a" "t" "t" "g" "g" "g" "a" "t" "c" "t" "g" "a" "t"
     [343] "t" "c" "t" "t" "c" "t" "g" "a" "a" "g" "a" "t" "a" "c" "c" "g" "t" "t"
     [361] "a" "a" "t" "a" "a" "g" "g" "c" "a" "a" "c" "t" "t" "a" "t" "t" "g" "c"
     [379] "a" "g" "t" "g" "t" "g" "g" "g" "a" "g" "a" "t" "c" "a" "a" "g" "a" "a"
     [397] "t" "t" "g" "t" "t" "a" "c" "a" "a" "a" "t" "c" "a" "c" "c" "c" "c" "t"
     [415] "c" "a" "a" "g" "g" "a" "a" "c" "c" "a" "g" "g" "g" "a" "t" "g" "a" "a"
     [433] "a" "t" "c" "a" "g" "t" "t" "t" "g" "g" "a" "c" "t" "c" "t" "g" "c" "a"
     [451] "a" "a" "a" "a" "a" "g" "g" "c" "t" "g" "c" "t" "t" "g" "t" "g" "a" "a"
     [469] "t" "t" "t" "t" "c" "t" "g" "a" "g" "a" "c" "g" "g" "a" "t" "g" "t" "a"
     [487] "a" "c" "a" "a" "a" "t" "a" "c" "t" "g" "a" "a" "c" "a" "t" "c" "a" "t"
     [505] "c" "a" "a" "c" "c" "c" "a" "g" "t" "a" "a" "t" "a" "a" "t" "g" "a" "t"
     [523] "t" "t" "g" "a" "a" "c" "a" "c" "c" "a" "c" "t" "g" "a" "g" "a" "a" "g"
     [541] "c" "g" "t" "g" "c" "a" "g" "c" "t" "g" "a" "g" "a" "g" "g" "c" "a" "t"
     [559] "c" "c" "a" "g" "a" "a" "a" "a" "g" "t" "a" "t" "c" "a" "g" "g" "g" "t"
     [577] "a" "g" "t" "t" "c" "t" "g" "t" "t" "t" "c" "a" "a" "a" "c" "t" "t" "g"
     [595] "c" "a" "t" "g" "t" "g" "g" "a" "g" "c" "c" "a" "t" "g" "t" "g" "g" "c"
     [613] "a" "c" "a" "a" "a" "t" "a" "c" "t" "c" "a" "t" "g" "c" "c" "a" "g" "c"
     [631] "t" "c" "a" "t" "t" "a" "c" "a" "g" "c" "a" "t" "g" "a" "g" "a" "a" "c"
     [649] "a" "g" "c" "a" "g" "t" "t" "t" "a" "t" "t" "a" "c" "t" "c" "a" "c" "t"
     [667] "a" "a" "a" "g" "a" "c" "a" "g" "a" "a" "t" "g" "a" "a" "t" "g" "t" "a"
     [685] "g" "a" "a" "a" "a" "g" "g" "c" "t" "g" "a" "a" "t" "t" "c" "t" "g" "t"
     [703] "a" "a" "t" "a" "a" "a" "a" "g" "c" "a" "a" "a" "c" "a" "g" "c" "c" "t"
     [721] "g" "g" "c" "t" "t" "a" "g" "c" "a" "a" "g" "g" "a" "g" "c" "c" "a" "a"
     [739] "c" "a" "t" "a" "a" "c" "a" "g" "a" "t" "g" "g" "g" "c" "t" "g" "g" "a"
     [757] "a" "g" "t" "a" "a" "g" "g" "a" "a" "a" "c" "a" "t" "g" "t" "a" "a" "t"
     [775] "g" "a" "t" "a" "g" "g" "c" "g" "g" "a" "c" "t" "c" "c" "c" "a" "g" "c"
     [793] "a" "c" "a" "g" "a" "a" "a" "a" "a" "a" "a" "g" "g" "t" "a" "g" "a" "t"
     [811] "c" "t" "g" "a" "a" "t" "g" "c" "t" "g" "a" "t" "c" "c" "c" "c" "t" "g"
     [829] "t" "g" "t" "g" "a" "g" "a" "g" "a" "a" "a" "a" "g" "a" "a" "t" "g" "g"
     [847] "a" "a" "t" "a" "a" "g" "c" "a" "g" "a" "a" "a" "c" "t" "g" "c" "c" "a"
     [865] "t" "g" "c" "t" "c" "a" "g" "a" "g" "a" "a" "t" "c" "c" "t" "a" "g" "a"
     [883] "g" "a" "t" "a" "c" "t" "g" "a" "a" "g" "a" "t" "g" "t" "t" "c" "c" "t"
     [901] "t" "g" "g" "a" "t" "a" "a" "c" "a" "c" "t" "a" "a" "a" "t" "a" "g" "c"
     [919] "a" "g" "c" "a" "t" "t" "c" "a" "g" "a" "a" "a" "g" "t" "t" "a" "a" "t"
     [937] "g" "a" "g" "t" "g" "g" "t" "t" "t" "t" "c" "c" "a" "g" "a" "a" "g" "t"
     [955] "g" "a" "t" "g" "a" "a" "c" "t" "g" "t" "t" "a" "g" "g" "t" "t" "c" "t"
     [973] "g" "a" "t" "g" "a" "c" "t" "c" "a" "c" "a" "t" "g" "a" "t" "g" "g" "g"
     [991] "g" "a" "g" "t" "c" "t" "g" "a" "a" "t" "c" "a" "a" "a" "t" "g" "c" "c"
    [1009] "a" "a" "a" "g" "t" "a" "g" "c" "t" "g" "a" "t" "g" "t" "a" "t" "t" "g"
    [1027] "g" "a" "c" "g" "t" "t" "c" "t" "a" "a" "a" "t" "g" "a" "g" "g" "t" "a"
    [1045] "g" "a" "t" "g" "a" "a" "t" "a" "t" "t" "c" "t" "g" "g" "t" "t" "c" "t"
    [1063] "t" "c" "a" "g" "a" "g" "a" "a" "a" "a" "t" "a" "g" "a" "c" "t" "t" "a"
    [1081] "c" "t" "g" "g" "c" "c" "a" "g" "t" "g" "a" "t" "c" "c" "t" "c" "a" "t"
    [1099] "g" "a" "g" "g" "c" "t" "t" "t" "a" "a" "t" "a" "t" "g" "t" "a" "a" "a"
    [1117] "a" "g" "t" "g" "a" "a" "a" "g" "a" "g" "t" "t" "c" "a" "c" "t" "c" "c"
    [1135] "a" "a" "a" "t" "c" "a" "g" "t" "a" "g" "a" "g" "a" "g" "t" "a" "a" "t"
    [1153] "a" "t" "t" "g" "a" "a" "g" "a" "c" "a" "a" "a" "a" "t" "a" "t" "t" "t"
    [1171] "g" "g" "g" "a" "a" "a" "a" "c" "c" "t" "a" "t" "c" "g" "g" "a" "a" "g"
    [1189] "a" "a" "g" "g" "c" "a" "a" "g" "c" "c" "t" "c" "c" "c" "c" "a" "a" "c"
    [1207] "t" "t" "a" "a" "g" "c" "c" "a" "t" "g" "t" "a" "a" "c" "t" "g" "a" "a"
    [1225] "a" "a" "t" "c" "t" "a" "a" "t" "t" "a" "t" "a" "g" "g" "a" "g" "c" "a"
    [1243] "t" "t" "t" "g" "t" "t" "a" "c" "t" "g" "a" "g" "c" "c" "a" "c" "a" "g"
    [1261] "a" "t" "a" "a" "t" "a" "c" "a" "a" "g" "a" "g" "c" "g" "t" "c" "c" "c"
    [1279] "c" "t" "c" "a" "c" "a" "a" "a" "t" "a" "a" "a" "t" "t" "a" "a" "a" "g"
    [1297] "c" "g" "t" "a" "a" "a" "a" "g" "g" "a" "g" "a" "c" "c" "t" "a" "c" "a"
    [1315] "t" "c" "a" "g" "g" "c" "c" "t" "t" "c" "a" "t" "c" "c" "t" "g" "a" "g"
    [1333] "g" "a" "t" "t" "t" "t" "a" "t" "c" "a" "a" "g" "a" "a" "a" "g" "c" "a"
    [1351] "g" "a" "t" "t" "t" "g" "g" "c" "a" "g" "t" "t" "c" "a" "a" "a" "a" "g"
    [1369] "a" "c" "t" "c" "c" "t" "g" "a" "a" "a" "t" "g" "a" "t" "a" "a" "a" "t"
    [1387] "c" "a" "g" "g" "g" "a" "a" "c" "t" "a" "a" "c" "c" "a" "a" "a" "c" "g"
    [1405] "g" "a" "g" "c" "a" "g" "a" "a" "t" "g" "g" "t" "c" "a" "a" "g" "t" "g"
    [1423] "a" "t" "g" "a" "a" "t" "a" "t" "t" "a" "c" "t" "a" "a" "t" "a" "g" "t"
    [1441] "g" "g" "t" "c" "a" "t" "g" "a" "g" "a" "a" "t" "a" "a" "a" "a" "c" "a"
    [1459] "a" "a" "a" "g" "g" "t" "g" "a" "t" "t" "c" "t" "a" "t" "t" "c" "a" "g"
    [1477] "a" "a" "t" "g" "a" "g" "a" "a" "a" "a" "a" "t" "c" "c" "t" "a" "a" "c"
    [1495] "c" "c" "a" "a" "t" "a" "g" "a" "a" "t" "c" "a" "c" "t" "c" "g" "a" "a"
    [1513] "a" "a" "a" "g" "a" "a" "t" "c" "t" "g" "c" "t" "t" "t" "c" "a" "a" "a"
    [1531] "a" "c" "g" "a" "a" "a" "g" "c" "t" "g" "a" "a" "c" "c" "t" "a" "t" "a"
    [1549] "a" "g" "c" "a" "g" "c" "a" "g" "t" "a" "t" "a" "a" "g" "c" "a" "a" "t"
    [1567] "a" "t" "g" "g" "a" "a" "c" "t" "c" "g" "a" "a" "t" "t" "a" "a" "a" "t"
    [1585] "a" "t" "c" "c" "a" "c" "a" "a" "t" "t" "c" "a" "a" "a" "a" "g" "c" "a"
    [1603] "c" "c" "t" "a" "a" "a" "a" "a" "g" "a" "a" "t" "a" "g" "g" "c" "t" "g"
    [1621] "a" "g" "g" "a" "g" "g" "a" "a" "g" "t" "c" "t" "t" "c" "t" "a" "c" "c"
    [1639] "a" "g" "g" "c" "a" "t" "a" "t" "t" "c" "a" "t" "g" "c" "g" "c" "t" "t"
    [1657] "g" "a" "a" "c" "t" "a" "g" "t" "a" "g" "t" "c" "a" "g" "t" "a" "g" "a"
    [1675] "a" "a" "t" "c" "t" "a" "a" "g" "c" "c" "c" "a" "c" "c" "t" "a" "a" "t"
    [1693] "t" "g" "t" "a" "c" "t" "g" "a" "a" "t" "t" "g" "c" "a" "a" "a" "t" "t"
    [1711] "g" "a" "t" "a" "g" "t" "t" "g" "t" "t" "c" "t" "a" "g" "c" "a" "g" "t"
    [1729] "g" "a" "a" "g" "a" "g" "a" "t" "a" "a" "a" "g" "a" "a" "a" "a" "a" "a"
    [1747] "a" "a" "g" "t" "a" "c" "a" "a" "c" "c" "a" "a" "a" "t" "g" "c" "c" "a"
    [1765] "g" "t" "c" "a" "g" "g" "c" "a" "c" "a" "g" "c" "a" "g" "a" "a" "a" "c"
    [1783] "c" "t" "a" "c" "a" "a" "c" "t" "c" "a" "t" "g" "g" "a" "a" "g" "g" "t"
    [1801] "a" "a" "a" "g" "a" "a" "c" "c" "t" "g" "c" "a" "a" "c" "t" "g" "g" "a"
    [1819] "g" "c" "c" "a" "a" "g" "a" "a" "g" "a" "g" "t" "a" "a" "c" "a" "a" "g"
    [1837] "c" "c" "a" "a" "a" "t" "g" "a" "a" "c" "a" "g" "a" "c" "a" "a" "g" "t"
    [1855] "a" "a" "a" "a" "g" "a" "c" "a" "t" "g" "a" "c" "a" "g" "c" "g" "a" "t"
    [1873] "a" "c" "t" "t" "t" "c" "c" "c" "a" "g" "a" "g" "c" "t" "g" "a" "a" "g"
    [1891] "t" "t" "a" "a" "c" "a" "a" "a" "t" "g" "c" "a" "c" "c" "t" "g" "g" "t"
    [1909] "t" "c" "t" "t" "t" "t" "a" "c" "t" "a" "a" "g" "t" "g" "t" "t" "c" "a"
    [1927] "a" "a" "t" "a" "c" "c" "a" "g" "t" "g" "a" "a" "c" "t" "t" "a" "a" "a"
    [1945] "g" "a" "a" "t" "t" "t" "g" "t" "c" "a" "a" "t" "c" "c" "t" "a" "g" "c"
    [1963] "c" "t" "t" "c" "c" "a" "a" "g" "a" "g" "a" "a" "g" "a" "a" "a" "a" "a"
    [1981] "g" "a" "a" "g" "a" "g" "a" "a" "a" "c" "t" "a" "g" "a" "a" "a" "c" "a"
    [1999] "g" "t" "t" "a" "a" "a" "g" "t" "g" "t" "c" "t" "a" "a" "t" "a" "a" "t"
    [2017] "g" "c" "t" "g" "a" "a" "g" "a" "c" "c" "c" "c" "a" "a" "a" "g" "a" "t"
    [2035] "c" "t" "c" "a" "t" "g" "t" "t" "a" "a" "g" "t" "g" "g" "a" "g" "a" "a"
    [2053] "a" "g" "g" "g" "t" "t" "t" "t" "g" "c" "a" "a" "a" "c" "t" "g" "a" "a"
    [2071] "a" "g" "a" "t" "c" "t" "g" "t" "a" "g" "a" "g" "a" "g" "t" "a" "g" "c"
    [2089] "a" "g" "t" "a" "t" "t" "t" "c" "a" "t" "t" "g" "g" "t" "a" "c" "c" "t"
    [2107] "g" "g" "t" "a" "c" "t" "g" "a" "t" "t" "a" "t" "g" "g" "c" "a" "c" "t"
    [2125] "c" "a" "g" "g" "a" "a" "a" "g" "t" "a" "t" "c" "t" "c" "g" "t" "t" "a"
    [2143] "c" "t" "g" "g" "a" "a" "g" "t" "t" "a" "g" "c" "a" "c" "t" "c" "t" "a"
    [2161] "g" "g" "g" "a" "a" "g" "g" "c" "a" "a" "a" "a" "a" "c" "a" "g" "a" "a"
    [2179] "c" "c" "a" "a" "a" "t" "a" "a" "a" "t" "g" "t" "g" "t" "g" "a" "g" "t"
    [2197] "c" "a" "g" "t" "g" "t" "g" "c" "a" "g" "c" "a" "t" "t" "t" "g" "a" "a"
    [2215] "a" "a" "c" "c" "c" "c" "a" "a" "g" "g" "g" "a" "c" "t" "a" "a" "t" "t"
    [2233] "c" "a" "t" "g" "g" "t" "t" "g" "t" "t" "c" "c" "a" "a" "a" "g" "a" "t"
    [2251] "a" "a" "t" "a" "g" "a" "a" "a" "t" "g" "a" "c" "a" "c" "a" "g" "a" "a"
    [2269] "g" "g" "c" "t" "t" "t" "a" "a" "g" "t" "a" "t" "c" "c" "a" "t" "t" "g"
    [2287] "g" "g" "a" "c" "a" "t" "g" "a" "a" "g" "t" "t" "a" "a" "c" "c" "a" "c"
    [2305] "a" "g" "t" "c" "g" "g" "g" "a" "a" "a" "c" "a" "a" "g" "c" "a" "t" "a"
    [2323] "g" "a" "a" "a" "t" "g" "g" "a" "a" "g" "a" "a" "a" "g" "t" "g" "a" "a"
    [2341] "c" "t" "t" "g" "a" "t" "g" "c" "t" "c" "a" "g" "t" "a" "t" "t" "t" "g"
    [2359] "c" "a" "g" "a" "a" "t" "a" "c" "a" "t" "t" "c" "a" "a" "g" "g" "t" "t"
    [2377] "t" "c" "a" "a" "a" "g" "c" "g" "c" "c" "a" "g" "t" "c" "a" "t" "t" "t"
    [2395] "g" "c" "t" "c" "c" "g" "t" "t" "t" "t" "c" "a" "a" "a" "t" "c" "c" "a"
    [2413] "g" "g" "a" "a" "a" "t" "g" "c" "a" "g" "a" "a" "g" "a" "g" "g" "a" "a"
    [2431] "t" "g" "t" "g" "c" "a" "a" "c" "a" "t" "t" "c" "t" "c" "t" "g" "c" "c"
    [2449] "c" "a" "c" "t" "c" "t" "g" "g" "g" "t" "c" "c" "t" "t" "a" "a" "a" "g"
    [2467] "a" "a" "a" "c" "a" "a" "a" "g" "t" "c" "c" "a" "a" "a" "a" "g" "t" "c"
    [2485] "a" "c" "t" "t" "t" "t" "g" "a" "a" "t" "g" "t" "g" "a" "a" "c" "a" "a"
    [2503] "a" "a" "g" "g" "a" "a" "g" "a" "a" "a" "a" "t" "c" "a" "a" "g" "g" "a"
    [2521] "a" "a" "g" "a" "a" "t" "g" "a" "g" "t" "c" "t" "a" "a" "t" "a" "t" "c"
    [2539] "a" "a" "g" "c" "c" "t" "g" "t" "a" "c" "a" "g" "a" "c" "a" "g" "t" "t"
    [2557] "a" "a" "t" "a" "t" "c" "a" "c" "t" "g" "c" "a" "g" "g" "c" "t" "t" "t"
    [2575] "c" "c" "t" "g" "t" "g" "g" "t" "t" "g" "g" "t" "c" "a" "g" "a" "a" "a"
    [2593] "g" "a" "t" "a" "a" "g" "c" "c" "a" "g" "t" "t" "g" "a" "t" "a" "a" "t"
    [2611] "g" "c" "c" "a" "a" "a" "t" "g" "t" "a" "g" "t" "a" "t" "c" "a" "a" "a"
    [2629] "g" "g" "a" "g" "g" "c" "t" "c" "t" "a" "g" "g" "t" "t" "t" "t" "g" "t"
    [2647] "c" "t" "a" "t" "c" "a" "t" "c" "t" "c" "a" "g" "t" "t" "c" "a" "g" "a"
    [2665] "g" "g" "c" "a" "a" "c" "g" "a" "a" "a" "c" "t" "g" "g" "a" "c" "t" "c"
    [2683] "a" "t" "t" "a" "c" "t" "c" "c" "a" "a" "a" "t" "a" "a" "a" "c" "a" "t"
    [2701] "g" "g" "a" "c" "t" "t" "t" "t" "a" "c" "a" "a" "a" "a" "c" "c" "c" "a"
    [2719] "t" "a" "t" "c" "g" "t" "a" "t" "a" "c" "c" "a" "c" "c" "a" "c" "t" "t"
    [2737] "t" "t" "t" "c" "c" "c" "a" "t" "c" "a" "a" "g" "t" "c" "a" "t" "t" "t"
    [2755] "g" "t" "t" "a" "a" "a" "a" "c" "t" "a" "a" "a" "t" "g" "t" "a" "a" "g"
    [2773] "a" "a" "a" "a" "a" "t" "c" "t" "g" "c" "t" "a" "g" "a" "g" "g" "a" "a"
    [2791] "a" "a" "c" "t" "t" "t" "g" "a" "g" "g" "a" "a" "c" "a" "t" "t" "c" "a"
    [2809] "a" "t" "g" "t" "c" "a" "c" "c" "t" "g" "a" "a" "a" "g" "a" "g" "a" "a"
    [2827] "a" "t" "g" "g" "g" "a" "a" "a" "t" "g" "a" "g" "a" "a" "c" "a" "t" "t"
    [2845] "c" "c" "a" "a" "g" "t" "a" "c" "a" "g" "t" "g" "a" "g" "c" "a" "c" "a"
    [2863] "a" "t" "t" "a" "g" "c" "c" "g" "t" "a" "a" "t" "a" "a" "c" "a" "t" "t"
    [2881] "a" "g" "a" "g" "a" "a" "a" "a" "t" "g" "t" "t" "t" "t" "t" "a" "a" "a"
    [2899] "g" "a" "a" "g" "c" "c" "a" "g" "c" "t" "c" "a" "a" "g" "c" "a" "a" "t"
    [2917] "a" "t" "t" "a" "a" "t" "g" "a" "a" "g" "t" "a" "g" "g" "t" "t" "c" "c"
    [2935] "a" "g" "t" "a" "c" "t" "a" "a" "t" "g" "a" "a" "g" "t" "g" "g" "g" "c"
    [2953] "t" "c" "c" "a" "g" "t" "a" "t" "t" "a" "a" "t" "g" "a" "a" "a" "t" "a"
    [2971] "g" "g" "t" "t" "c" "c" "a" "g" "t" "g" "a" "t" "g" "a" "a" "a" "a" "c"
    [2989] "a" "t" "t" "c" "a" "a" "g" "c" "a" "g" "a" "a" "c" "t" "a" "g" "g" "t"
    [3007] "a" "g" "a" "a" "a" "c" "a" "g" "a" "a" "g" "g" "c" "c" "c" "a" "a" "a"
    [3025] "t" "t" "g" "a" "a" "t" "g" "c" "t" "a" "t" "g" "c" "t" "t" "a" "g" "a"
    [3043] "t" "t" "a" "g" "g" "g" "g" "t" "t" "t" "t" "g" "c" "a" "a" "c" "c" "t"
    [3061] "g" "a" "g" "g" "t" "c" "t" "a" "t" "a" "a" "a" "c" "a" "a" "a" "g" "t"
    [3079] "c" "t" "t" "c" "c" "t" "g" "g" "a" "a" "g" "t" "a" "a" "t" "t" "g" "t"
    [3097] "a" "a" "g" "c" "a" "t" "c" "c" "t" "g" "a" "a" "a" "t" "a" "a" "a" "a"
    [3115] "a" "a" "g" "c" "a" "a" "g" "a" "a" "t" "a" "t" "g" "a" "a" "g" "a" "a"
    [3133] "g" "t" "a" "g" "t" "t" "c" "a" "g" "a" "c" "t" "g" "t" "t" "a" "a" "t"
    [3151] "a" "c" "a" "g" "a" "t" "t" "t" "c" "t" "c" "t" "c" "c" "a" "t" "a" "t"
    [3169] "c" "t" "g" "a" "t" "t" "t" "c" "a" "g" "a" "t" "a" "a" "c" "t" "t" "a"
    [3187] "g" "a" "a" "c" "a" "g" "c" "c" "t" "a" "t" "g" "g" "g" "a" "a" "g" "t"
    [3205] "a" "g" "t" "c" "a" "t" "g" "c" "a" "t" "c" "t" "c" "a" "g" "g" "t" "t"
    [3223] "t" "g" "t" "t" "c" "t" "g" "a" "g" "a" "c" "a" "c" "c" "t" "g" "a" "t"
    [3241] "g" "a" "c" "c" "t" "g" "t" "t" "a" "g" "a" "t" "g" "a" "t" "g" "g" "t"
    [3259] "g" "a" "a" "a" "t" "a" "a" "a" "g" "g" "a" "a" "g" "a" "t" "a" "c" "t"
    [3277] "a" "g" "t" "t" "t" "t" "g" "c" "t" "g" "a" "a" "a" "a" "t" "g" "a" "c"
    [3295] "a" "t" "t" "a" "a" "g" "g" "a" "a" "a" "g" "t" "t" "c" "t" "g" "c" "t"
    [3313] "g" "t" "t" "t" "t" "t" "a" "g" "c" "a" "a" "a" "a" "g" "c" "g" "t" "c"
    [3331] "c" "a" "g" "a" "a" "a" "g" "g" "a" "g" "a" "g" "c" "t" "t" "a" "g" "c"
    [3349] "a" "g" "g" "a" "g" "t" "c" "c" "t" "a" "g" "c" "c" "c" "t" "t" "t" "c"
    [3367] "a" "c" "c" "c" "a" "t" "a" "c" "a" "c" "a" "t" "t" "t" "g" "g" "c" "t"
    [3385] "c" "a" "g" "g" "g" "t" "t" "a" "c" "c" "g" "a" "a" "g" "a" "g" "g" "g"
    [3403] "g" "c" "c" "a" "a" "g" "a" "a" "a" "t" "t" "a" "g" "a" "g" "t" "c" "c"
    [3421] "t" "c" "a" "g" "a" "a" "g" "a" "g" "a" "a" "c" "t" "t" "a" "t" "c" "t"
    [3439] "a" "g" "t" "g" "a" "g" "g" "a" "t" "g" "a" "a" "g" "a" "g" "c" "t" "t"
    [3457] "c" "c" "c" "t" "g" "c" "t" "t" "c" "c" "a" "a" "c" "a" "c" "t" "t" "g"
    [3475] "t" "t" "a" "t" "t" "t" "g" "g" "t" "a" "a" "a" "g" "t" "a" "a" "a" "c"
    [3493] "a" "a" "t" "a" "t" "a" "c" "c" "t" "t" "c" "t" "c" "a" "g" "t" "c" "t"
    [3511] "a" "c" "t" "a" "g" "g" "c" "a" "t" "a" "g" "c" "a" "c" "c" "g" "t" "t"
    [3529] "g" "c" "t" "a" "c" "c" "g" "a" "g" "t" "g" "t" "c" "t" "g" "t" "c" "t"
    [3547] "a" "a" "g" "a" "a" "c" "a" "c" "a" "g" "a" "g" "g" "a" "g" "a" "a" "t"
    [3565] "t" "t" "a" "t" "t" "a" "t" "c" "a" "t" "t" "g" "a" "a" "g" "a" "a" "t"
    [3583] "a" "g" "c" "t" "t" "a" "a" "a" "t" "g" "a" "c" "t" "g" "c" "a" "g" "t"
    [3601] "a" "a" "c" "c" "a" "g" "g" "t" "a" "a" "t" "a" "t" "t" "g" "g" "c" "a"
    [3619] "a" "a" "g" "g" "c" "a" "t" "c" "t" "c" "a" "g" "g" "a" "a" "c" "a" "t"
    [3637] "c" "a" "c" "c" "t" "t" "a" "g" "t" "g" "a" "g" "g" "a" "a" "a" "c" "a"
    [3655] "a" "a" "a" "t" "g" "t" "t" "c" "t" "g" "c" "t" "a" "g" "c" "t" "t" "g"
    [3673] "t" "t" "t" "t" "c" "t" "t" "c" "a" "c" "a" "g" "t" "g" "c" "a" "g" "t"
    [3691] "g" "a" "a" "t" "t" "g" "g" "a" "a" "g" "a" "c" "t" "t" "g" "a" "c" "t"
    [3709] "g" "c" "a" "a" "a" "t" "a" "c" "a" "a" "a" "c" "a" "c" "c" "c" "a" "g"
    [3727] "g" "a" "t" "c" "c" "t" "t" "t" "c" "t" "t" "g" "a" "t" "t" "g" "g" "t"
    [3745] "t" "c" "t" "t" "c" "c" "a" "a" "a" "c" "a" "a" "a" "t" "g" "a" "g" "g"
    [3763] "c" "a" "t" "c" "a" "g" "t" "c" "t" "g" "a" "a" "a" "g" "c" "c" "a" "g"
    [3781] "g" "g" "a" "g" "t" "t" "g" "g" "t" "c" "t" "g" "a" "g" "t" "g" "a" "c"
    [3799] "a" "a" "g" "g" "a" "a" "t" "t" "g" "g" "t" "t" "t" "c" "a" "g" "a" "t"
    [3817] "g" "a" "t" "g" "a" "a" "g" "a" "a" "a" "g" "a" "g" "g" "a" "a" "c" "g"
    [3835] "g" "g" "c" "t" "t" "g" "g" "a" "a" "g" "a" "a" "a" "a" "t" "a" "a" "t"
    [3853] "c" "a" "a" "g" "a" "a" "g" "a" "g" "c" "a" "a" "a" "g" "c" "a" "t" "g"
    [3871] "g" "a" "t" "t" "c" "a" "a" "a" "c" "t" "t" "a" "g" "g" "t" "g" "a" "a"
    [3889] "g" "c" "a" "g" "c" "a" "t" "c" "t" "g" "g" "g" "t" "g" "t" "g" "a" "g"
    [3907] "a" "g" "t" "g" "a" "a" "a" "c" "a" "a" "g" "c" "g" "t" "c" "t" "c" "t"
    [3925] "g" "a" "a" "g" "a" "c" "t" "g" "c" "t" "c" "a" "g" "g" "g" "c" "t" "a"
    [3943] "t" "c" "c" "t" "c" "t" "c" "a" "g" "a" "g" "t" "g" "a" "c" "a" "t" "t"
    [3961] "t" "t" "a" "a" "c" "c" "a" "c" "t" "c" "a" "g" "c" "a" "g" "a" "g" "g"
    [3979] "g" "a" "t" "a" "c" "c" "a" "t" "g" "c" "a" "a" "c" "a" "t" "a" "a" "c"
    [3997] "c" "t" "g" "a" "t" "a" "a" "a" "g" "c" "t" "c" "c" "a" "g" "c" "a" "g"
    [4015] "g" "a" "a" "a" "t" "g" "g" "c" "t" "g" "a" "a" "c" "t" "a" "g" "a" "a"
    [4033] "g" "c" "t" "g" "t" "g" "t" "t" "a" "g" "a" "a" "c" "a" "g" "c" "a" "t"
    [4051] "g" "g" "g" "a" "g" "c" "c" "a" "g" "c" "c" "t" "t" "c" "t" "a" "a" "c"
    [4069] "a" "g" "c" "t" "a" "c" "c" "c" "t" "t" "c" "c" "a" "t" "c" "a" "t" "a"
    [4087] "a" "g" "t" "g" "a" "c" "t" "c" "t" "t" "c" "t" "g" "c" "c" "c" "t" "t"
    [4105] "g" "a" "g" "g" "a" "c" "c" "t" "g" "c" "g" "a" "a" "a" "t" "c" "c" "a"
    [4123] "g" "a" "a" "c" "a" "a" "a" "g" "c" "a" "c" "a" "t" "c" "a" "g" "a" "a"
    [4141] "a" "a" "a" "g" "c" "a" "g" "t" "a" "t" "t" "a" "a" "c" "t" "t" "c" "a"
    [4159] "c" "a" "g" "a" "a" "a" "a" "g" "t" "a" "g" "t" "g" "a" "a" "t" "a" "c"
    [4177] "c" "c" "t" "a" "t" "a" "a" "g" "c" "c" "a" "g" "a" "a" "t" "c" "c" "a"
    [4195] "g" "a" "a" "g" "g" "c" "c" "t" "t" "t" "c" "t" "g" "c" "t" "g" "a" "c"
    [4213] "a" "a" "g" "t" "t" "t" "g" "a" "g" "g" "t" "g" "t" "c" "t" "g" "c" "a"
    [4231] "g" "a" "t" "a" "g" "t" "t" "c" "t" "a" "c" "c" "a" "g" "t" "a" "a" "a"
    [4249] "a" "a" "t" "a" "a" "a" "g" "a" "a" "c" "c" "a" "g" "g" "a" "g" "t" "g"
    [4267] "g" "a" "a" "a" "g" "g" "t" "c" "a" "t" "c" "c" "c" "c" "t" "t" "c" "t"
    [4285] "a" "a" "a" "t" "g" "c" "c" "c" "a" "t" "c" "a" "t" "t" "a" "g" "a" "t"
    [4303] "g" "a" "t" "a" "g" "g" "t" "g" "g" "t" "a" "c" "a" "t" "g" "c" "a" "c"
    [4321] "a" "g" "t" "t" "g" "c" "t" "c" "t" "g" "g" "g" "a" "g" "t" "c" "t" "t"
    [4339] "c" "a" "g" "a" "a" "t" "a" "g" "a" "a" "a" "c" "t" "a" "c" "c" "c" "a"
    [4357] "t" "c" "t" "c" "a" "a" "g" "a" "g" "g" "a" "g" "c" "t" "c" "a" "t" "t"
    [4375] "a" "a" "g" "g" "t" "t" "g" "t" "t" "g" "a" "t" "g" "t" "g" "g" "a" "g"
    [4393] "g" "a" "g" "c" "a" "a" "c" "a" "g" "c" "t" "g" "g" "a" "a" "g" "a" "g"
    [4411] "t" "c" "t" "g" "g" "g" "c" "c" "a" "c" "a" "c" "g" "a" "t" "t" "t" "g"
    [4429] "a" "c" "g" "g" "a" "a" "a" "c" "a" "t" "c" "t" "t" "a" "c" "t" "t" "g"
    [4447] "c" "c" "a" "a" "g" "g" "c" "a" "a" "g" "a" "t" "c" "t" "a" "g" "a" "g"
    [4465] "g" "g" "a" "a" "c" "c" "c" "c" "t" "t" "a" "c" "c" "t" "g" "g" "a" "a"
    [4483] "t" "c" "t" "g" "g" "a" "a" "t" "c" "a" "g" "c" "c" "t" "c" "t" "t" "c"
    [4501] "t" "c" "t" "g" "a" "t" "g" "a" "c" "c" "c" "t" "g" "a" "a" "t" "c" "t"
    [4519] "g" "a" "t" "c" "c" "t" "t" "c" "t" "g" "a" "a" "g" "a" "c" "a" "g" "a"
    [4537] "g" "c" "c" "c" "c" "a" "g" "a" "g" "t" "c" "a" "g" "c" "t" "c" "g" "t"
    [4555] "g" "t" "t" "g" "g" "c" "a" "a" "c" "a" "t" "a" "c" "c" "a" "t" "c" "t"
    [4573] "t" "c" "a" "a" "c" "c" "t" "c" "t" "g" "c" "a" "t" "t" "g" "a" "a" "a"
    [4591] "g" "t" "t" "c" "c" "c" "c" "a" "a" "t" "t" "g" "a" "a" "a" "g" "t" "t"
    [4609] "g" "c" "a" "g" "a" "a" "t" "c" "t" "g" "c" "c" "c" "a" "g" "a" "g" "t"
    [4627] "c" "c" "a" "g" "c" "t" "g" "c" "t" "g" "c" "t" "c" "a" "t" "a" "c" "t"
    [4645] "a" "c" "t" "g" "a" "t" "a" "c" "t" "g" "c" "t" "g" "g" "g" "t" "a" "t"
    [4663] "a" "a" "t" "g" "c" "a" "a" "t" "g" "g" "a" "a" "g" "a" "a" "a" "g" "t"
    [4681] "g" "t" "g" "a" "g" "c" "a" "g" "g" "g" "a" "g" "a" "a" "g" "c" "c" "a"
    [4699] "g" "a" "a" "t" "t" "g" "a" "c" "a" "g" "c" "t" "t" "c" "a" "a" "c" "a"
    [4717] "g" "a" "a" "a" "g" "g" "g" "t" "c" "a" "a" "c" "a" "a" "a" "a" "g" "a"
    [4735] "a" "t" "g" "t" "c" "c" "a" "t" "g" "g" "t" "g" "g" "t" "g" "t" "c" "t"
    [4753] "g" "g" "c" "c" "t" "g" "a" "c" "c" "c" "c" "a" "g" "a" "a" "g" "a" "a"
    [4771] "t" "t" "t" "a" "t" "g" "c" "t" "c" "g" "t" "g" "t" "a" "c" "a" "a" "g"
    [4789] "t" "t" "t" "g" "c" "c" "a" "g" "a" "a" "a" "a" "c" "a" "c" "c" "a" "c"
    [4807] "a" "t" "c" "a" "c" "t" "t" "t" "a" "a" "c" "t" "a" "a" "t" "c" "t" "a"
    [4825] "a" "t" "t" "a" "c" "t" "g" "a" "a" "g" "a" "g" "a" "c" "t" "a" "c" "t"
    [4843] "c" "a" "t" "g" "t" "t" "g" "t" "t" "a" "t" "g" "a" "a" "a" "a" "c" "a"
    [4861] "g" "a" "t" "g" "c" "t" "g" "a" "g" "t" "t" "t" "g" "t" "g" "t" "g" "t"
    [4879] "g" "a" "a" "c" "g" "g" "a" "c" "a" "c" "t" "g" "a" "a" "a" "t" "a" "t"
    [4897] "t" "t" "t" "c" "t" "a" "g" "g" "a" "a" "t" "t" "g" "c" "g" "g" "g" "a"
    [4915] "g" "g" "a" "a" "a" "a" "t" "g" "g" "g" "t" "a" "g" "t" "t" "a" "g" "c"
    [4933] "t" "a" "t" "t" "t" "c" "t" "g" "g" "g" "t" "g" "a" "c" "c" "c" "a" "g"
    [4951] "t" "c" "t" "a" "t" "t" "a" "a" "a" "g" "a" "a" "a" "g" "a" "a" "a" "a"
    [4969] "a" "t" "g" "c" "t" "g" "a" "a" "t" "g" "a" "g" "c" "a" "t" "g" "a" "t"
    [4987] "t" "t" "t" "g" "a" "a" "g" "t" "c" "a" "g" "a" "g" "g" "a" "g" "a" "t"
    [5005] "g" "t" "g" "g" "t" "c" "a" "a" "t" "g" "g" "a" "a" "g" "a" "a" "a" "c"
    [5023] "c" "a" "c" "c" "a" "a" "g" "g" "t" "c" "c" "a" "a" "a" "g" "c" "g" "a"
    [5041] "g" "c" "a" "a" "g" "a" "g" "a" "a" "t" "c" "c" "c" "a" "g" "g" "a" "c"
    [5059] "a" "g" "a" "a" "a" "g" "a" "t" "c" "t" "t" "c" "a" "g" "g" "g" "g" "g"
    [5077] "c" "t" "a" "g" "a" "a" "a" "t" "c" "t" "g" "t" "t" "g" "c" "t" "a" "t"
    [5095] "g" "g" "g" "c" "c" "c" "t" "t" "c" "a" "c" "c" "a" "a" "c" "a" "t" "g"
    [5113] "c" "c" "c" "a" "c" "a" "g" "a" "t" "c" "a" "a" "c" "t" "g" "g" "a" "a"
    [5131] "t" "g" "g" "a" "t" "g" "g" "t" "a" "c" "a" "g" "c" "t" "g" "t" "g" "t"
    [5149] "g" "g" "t" "g" "c" "t" "t" "c" "t" "g" "t" "g" "g" "t" "g" "a" "a" "g"
    [5167] "g" "a" "g" "c" "t" "t" "t" "c" "a" "t" "c" "a" "t" "t" "c" "a" "c" "c"
    [5185] "c" "t" "t" "g" "g" "c" "a" "c" "a" "g" "g" "t" "g" "t" "c" "c" "a" "c"
    [5203] "c" "c" "a" "a" "t" "t" "g" "t" "g" "g" "t" "t" "g" "t" "g" "c" "a" "g"
    [5221] "c" "c" "a" "g" "a" "t" "g" "c" "c" "t" "g" "g" "a" "c" "a" "g" "a" "g"
    [5239] "g" "a" "c" "a" "a" "t" "g" "g" "c" "t" "t" "c" "c" "a" "t" "g" "c" "a"
    [5257] "a" "t" "t" "g" "g" "g" "c" "a" "g" "a" "t" "g" "t" "g" "t" "g" "a" "g"
    [5275] "g" "c" "a" "c" "c" "t" "g" "t" "g" "g" "t" "g" "a" "c" "c" "c" "g" "a"
    [5293] "g" "a" "g" "t" "g" "g" "g" "t" "g" "t" "t" "g" "g" "a" "c" "a" "g" "t"
    [5311] "g" "t" "a" "g" "c" "a" "c" "t" "c" "t" "a" "c" "c" "a" "g" "t" "g" "c"
    [5329] "c" "a" "g" "g" "a" "g" "c" "t" "g" "g" "a" "c" "a" "c" "c" "t" "a" "c"
    [5347] "c" "t" "g" "a" "t" "a" "c" "c" "c" "c" "a" "g" "a" "t" "c" "c" "c" "c"
    [5365] "c" "a" "c" "a" "g" "c" "c" "a" "c" "t" "a" "c" "t" "g" "a"
    attr(,"name")
    [1] "AF005068.BRCA1"
    attr(,"Annot")
    [1] ">AF005068.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $AF284812.BRCA1
     [1] "c" "a" "t" "g" "a" "t" "t" "t" "t" "g" "a" "a" "g" "t" "c" "a" "g" "a" "g"
    [20] "g" "a" "g" "a" "t" "g" "t" "g" "g" "t" "c" "a" "a" "t" "g" "g" "a" "a" "g"
    [39] "a" "a" "a" "c" "c" "a" "c" "c" "a" "a" "g" "g" "t" "c" "c" "a" "a" "a" "g"
    [58] "c" "g" "a" "g" "c" "a" "a" "g" "a" "g" "a" "a" "t" "c" "c" "c" "a" "g" "g"
    [77] "a" "c" "a" "g" "a" "a" "a" "g"
    attr(,"name")
    [1] "AF284812.BRCA1"
    attr(,"Annot")
    [1] ">AF284812.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $AF507076.BRCA1
     [1] "c" "a" "t" "g" "a" "t" "t" "t" "t" "g" "a" "a" "g" "t" "c" "a" "g" "a" "g"
    [20] "g" "a" "g" "a" "t" "g" "t" "g" "g" "t" "c" "a" "a" "t" "g" "g" "a" "a" "g"
    [39] "a" "a" "a" "c" "c" "a" "c" "c" "a" "a" "g" "g" "t" "c" "c" "a" "a" "a" "g"
    [58] "c" "g" "a" "g" "c" "a" "a" "g" "a" "g" "a" "a" "t" "c" "c" "c" "a" "g" "g"
    [77] "a" "c" "a" "g" "a" "a" "a" "g"
    attr(,"name")
    [1] "AF507076.BRCA1"
    attr(,"Annot")
    [1] ">AF507076.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $AF507077.BRCA1
     [1] "c" "a" "t" "g" "a" "t" "t" "t" "t" "g" "a" "a" "g" "t" "c" "c" "g" "a" "g"
    [20] "g" "a" "g" "a" "t" "g" "t" "g" "g" "g" "c" "a" "a" "t" "g" "g" "a" "a" "g"
    [39] "a" "a" "a" "c" "c" "a" "c" "c" "a" "a" "g" "g" "t" "g" "c" "a" "a" "a" "g"
    [58] "c" "g" "a" "g" "c" "a" "a" "g" "a" "g" "a" "a" "t" "c" "c" "c" "a" "g" "g"
    [77] "a" "c" "a" "g" "a" "a" "a" "g"
    attr(,"name")
    [1] "AF507077.BRCA1"
    attr(,"Annot")
    [1] ">AF507077.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $AF507078.BRCA1
     [1] "c" "a" "t" "g" "a" "t" "t" "t" "t" "g" "a" "a" "g" "t" "c" "a" "g" "a" "g"
    [20] "g" "a" "g" "a" "t" "g" "t" "g" "g" "t" "c" "a" "a" "t" "g" "g" "a" "a" "g"
    [39] "a" "a" "a" "c" "c" "a" "c" "c" "a" "a" "g" "g" "t" "c" "c" "a" "a" "a" "g"
    [58] "c" "g" "a" "g" "c" "a" "a" "g" "a" "g" "a" "a" "t" "c" "c" "c" "a" "g" "g"
    [77] "a" "c" "a" "g" "a" "a" "a" "g"
    attr(,"name")
    [1] "AF507078.BRCA1"
    attr(,"Annot")
    [1] ">AF507078.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $AY093484.BRCA1
     [1] "a" "a" "a" "c" "c" "a" "c" "c" "a" "a" "g" "g" "t" "c" "c" "a" "a" "a" "g"
    [20] "c" "g" "a" "g" "c" "a" "a" "g" "a" "g" "a" "a" "t" "c" "c" "c" "a" "g" "g"
    [39] "a" "c" "a" "g" "a" "a" "a" "g"
    attr(,"name")
    [1] "AY093484.BRCA1"
    attr(,"Annot")
    [1] ">AY093484.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $AY093485.BRCA1
     [1] "c" "a" "t" "g" "a" "t" "t" "t" "t" "g" "a" "a" "g" "t" "c" "a" "g" "a" "g"
    [20] "g" "a" "g" "a" "t" "g" "t" "g" "g" "t" "c" "a" "a" "t" "g" "g" "a" "a" "g"
    [39] "a" "a" "a" "c" "c" "a" "c" "c" "a" "a" "g" "g" "t" "c" "c" "a" "a" "a" "g"
    [58] "c" "g" "a" "g" "c" "a" "a" "g" "a" "g" "a" "a" "t" "c" "c" "c" "a" "g" "g"
    [77] "a" "c" "a" "g" "a" "a" "a" "g"
    attr(,"name")
    [1] "AY093485.BRCA1"
    attr(,"Annot")
    [1] ">AY093485.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $AY093486.BRCA1
     [1] "c" "a" "t" "g" "a" "t" "t" "t" "t" "g" "a" "a" "g" "t" "c" "a" "g" "a" "g"
    [20] "g" "a" "g" "a" "t" "g" "t" "g" "g" "t" "c" "a" "a" "t" "g" "g" "a" "a" "g"
    [39] "a" "a" "a" "c" "c" "a" "c" "c" "a" "a" "g" "g" "t" "c" "c" "a" "a" "a" "g"
    [58] "c" "g" "a" "g" "c" "a" "a" "g" "a" "g" "a" "a" "t" "c" "c" "c" "a" "g" "g"
    [77] "a" "c" "a" "g" "a" "a" "a" "g"
    attr(,"name")
    [1] "AY093486.BRCA1"
    attr(,"Annot")
    [1] ">AY093486.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $AY093487.BRCA1
     [1] "c" "a" "t" "g" "a" "t" "t" "t" "t" "g" "a" "a" "g" "t" "c" "a" "g" "a" "g"
    [20] "g" "a" "g" "a" "t" "g" "t" "g" "g" "t" "c" "a" "a" "t" "g" "g" "a" "a" "g"
    [39] "a" "a" "a" "c" "c" "a" "c" "c" "a" "a" "g" "g" "t" "c" "c" "a" "a" "a" "g"
    [58] "c" "g" "a" "g" "c" "a" "a" "g" "a" "g" "a" "a" "t" "c" "c" "c" "a" "g" "g"
    [77] "a" "c" "a" "g" "a" "a" "a" "g"
    attr(,"name")
    [1] "AY093487.BRCA1"
    attr(,"Annot")
    [1] ">AY093487.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $AY093489.BRCA1
     [1] "c" "a" "t" "g" "a" "t" "t" "t" "t" "g" "a" "a" "g" "t" "c" "a" "g" "a" "g"
    [20] "g" "a" "g" "a" "t" "g" "t" "g" "g" "t" "c" "a" "a" "t" "g" "g" "a" "a" "g"
    [39] "a" "a" "a" "c" "c" "a" "c" "c" "a" "a" "g" "g" "t" "c" "c" "a" "a" "a" "g"
    [58] "c" "g" "a" "g" "c" "a" "a" "g" "a" "g" "a" "a" "t" "c" "c" "c" "a" "g" "g"
    [77] "a" "c" "a" "g" "a" "a" "a" "g"
    attr(,"name")
    [1] "AY093489.BRCA1"
    attr(,"Annot")
    [1] ">AY093489.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $AY093490.BRCA1
     [1] "c" "a" "c" "g" "a" "t" "t" "t" "t" "g" "a" "a" "g" "t" "c" "a" "g" "a" "g"
    [20] "g" "a" "g" "a" "t" "g" "t" "g" "g" "t" "t" "a" "a" "t" "g" "g" "a" "a" "g"
    [39] "a" "a" "a" "c" "c" "a" "c" "c" "a" "a" "g" "g" "t" "c" "c" "a" "a" "a" "g"
    [58] "c" "g" "a" "g" "c" "a" "a" "g" "a" "g" "a" "a" "t" "c" "c" "c" "a" "g" "g"
    [77] "a" "c" "a" "g" "a" "a" "a" "g"
    attr(,"name")
    [1] "AY093490.BRCA1"
    attr(,"Annot")
    [1] ">AY093490.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $AY093491.BRCA1
     [1] "c" "a" "t" "g" "a" "t" "t" "t" "t" "g" "a" "a" "g" "t" "c" "a" "g" "a" "g"
    [20] "g" "a" "t" "a" "t" "g" "t" "g" "g" "t" "c" "a" "a" "t" "g" "g" "a" "a" "g"
    [39] "a" "a" "a" "c" "c" "a" "c" "c" "a" "a" "g" "g" "t" "c" "c" "a" "a" "a" "g"
    [58] "c" "g" "a" "g" "c" "a" "a" "g" "a" "g" "a" "a" "t" "c" "c" "c" "a" "g" "g"
    [77] "a" "c" "a" "g" "a" "a" "a" "g"
    attr(,"name")
    [1] "AY093491.BRCA1"
    attr(,"Annot")
    [1] ">AY093491.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $AY093492.BRCA1
     [1] "c" "a" "t" "g" "a" "t" "t" "c" "t" "g" "a" "a" "g" "t" "c" "a" "g" "a" "g"
    [20] "g" "a" "g" "a" "t" "g" "t" "g" "g" "t" "c" "a" "a" "t" "g" "g" "a" "a" "g"
    [39] "a" "a" "a" "c" "c" "a" "c" "c" "a" "a" "g" "g" "t" "c" "c" "a" "a" "a" "g"
    [58] "c" "g" "a" "g" "c" "a" "a" "g" "a" "g" "a" "a" "t" "c" "c" "c" "a" "g" "g"
    [77] "a" "c" "a" "g" "a" "a" "a" "g"
    attr(,"name")
    [1] "AY093492.BRCA1"
    attr(,"Annot")
    [1] ">AY093492.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $AY093493.BRCA1
     [1] "c" "a" "t" "g" "a" "t" "t" "t" "t" "g" "a" "a" "g" "t" "c" "a" "g" "a" "g"
    [20] "g" "a" "g" "a" "t" "g" "t" "g" "g" "t" "c" "a" "a" "t" "g" "g" "a" "a" "g"
    [39] "a" "a" "a" "c" "c" "a" "c" "c" "a" "a" "g" "g" "t" "c" "c" "a" "a" "a" "g"
    [58] "c" "g" "a" "g" "c" "a" "a" "g" "a" "g" "a" "a" "t" "c" "c" "c" "a" "g" "g"
    [77] "a" "c" "a" "g" "a" "a" "a" "g"
    attr(,"name")
    [1] "AY093493.BRCA1"
    attr(,"Annot")
    [1] ">AY093493.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $AY144588
     [1] "g" "t" "g" "a" "a" "g" "c" "a" "g" "c" "a" "t" "c" "t" "g" "g" "g" "t" "g"
    [20] "t" "g" "a" "g" "a" "g" "t" "g" "a" "a" "a" "c" "a" "a" "g" "c" "g" "t" "c"
    [39] "t" "c" "t" "g" "a" "a" "g" "a" "c" "t" "g" "c" "t" "c" "a" "g" "g" "g" "c"
    [58] "t" "a" "t" "c" "a" "g" "a" "g" "t" "g" "a"
    attr(,"name")
    [1] "AY144588"
    attr(,"Annot")
    [1] ">AY144588"
    attr(,"class")
    [1] "SeqFastadna"
    
    $AY150865
     [1] "g" "t" "g" "a" "a" "g" "c" "a" "g" "c" "a" "t" "c" "t" "g" "g" "g" "t" "g"
    [20] "t" "g" "a" "g" "a" "g" "t" "g" "a" "a" "a" "c" "a" "a" "g" "c" "g" "t" "c"
    [39] "t" "c" "t" "g" "a" "a" "g" "a" "c" "t" "g" "c" "t" "c" "a" "g" "g" "g" "c"
    [58] "t" "a" "t" "c" "a" "g" "a" "g" "t" "g" "a"
    attr(,"name")
    [1] "AY150865"
    attr(,"Annot")
    [1] ">AY150865"
    attr(,"class")
    [1] "SeqFastadna"
    
    $AY273801.BRCA1
       [1] "a" "t" "g" "g" "a" "t" "t" "t" "a" "t" "c" "t" "g" "c" "t" "c" "t" "t"
      [19] "c" "g" "c" "g" "t" "t" "g" "a" "a" "g" "a" "a" "g" "t" "a" "c" "a" "a"
      [37] "a" "a" "t" "g" "t" "c" "a" "t" "t" "a" "a" "t" "g" "c" "t" "a" "t" "g"
      [55] "c" "a" "g" "a" "a" "a" "a" "t" "c" "t" "t" "a" "g" "a" "g" "t" "g" "t"
      [73] "c" "c" "c" "a" "t" "c" "t" "g" "t" "c" "t" "g" "g" "a" "g" "t" "t" "g"
      [91] "a" "t" "c" "a" "a" "g" "g" "a" "a" "c" "c" "t" "g" "t" "c" "t" "c" "c"
     [109] "a" "c" "a" "a" "a" "g" "t" "g" "t" "g" "a" "c" "c" "a" "c" "a" "t" "a"
     [127] "t" "t" "t" "t" "g" "c" "a" "a" "a" "t" "t" "t" "t" "g" "c" "a" "t" "g"
     [145] "c" "t" "g" "a" "a" "a" "c" "t" "t" "c" "t" "c" "a" "a" "c" "c" "a" "g"
     [163] "a" "a" "g" "a" "a" "a" "g" "g" "g" "c" "c" "t" "t" "c" "a" "c" "a" "g"
     [181] "t" "g" "t" "c" "c" "t" "t" "t" "a" "t" "g" "t" "a" "a" "g" "a" "a" "t"
     [199] "g" "a" "t" "a" "t" "a" "a" "c" "c" "a" "a" "a" "a" "g" "g" "a" "g" "c"
     [217] "c" "t" "a" "c" "a" "a" "g" "a" "a" "a" "g" "t" "a" "c" "g" "a" "g" "a"
     [235] "t" "t" "t" "a" "g" "t" "c" "a" "a" "c" "t" "t" "g" "t" "t" "g" "a" "a"
     [253] "g" "a" "g" "c" "t" "a" "t" "t" "g" "a" "a" "a" "a" "t" "c" "a" "t" "t"
     [271] "t" "g" "t" "g" "c" "t" "t" "t" "t" "c" "a" "g" "c" "t" "t" "g" "a" "c"
     [289] "a" "c" "a" "g" "g" "t" "t" "t" "g" "g" "a" "g" "t" "a" "t" "g" "c" "a"
     [307] "a" "a" "c" "a" "g" "c" "t" "a" "t" "a" "a" "t" "t" "t" "t" "g" "c" "a"
     [325] "a" "a" "a" "a" "a" "g" "g" "a" "a" "a" "a" "t" "a" "a" "c" "t" "c" "t"
     [343] "c" "c" "t" "g" "a" "a" "c" "a" "t" "c" "t" "a" "a" "a" "a" "g" "a" "t"
     [361] "g" "a" "a" "g" "t" "t" "t" "c" "t" "a" "t" "c" "a" "t" "c" "c" "a" "a"
     [379] "a" "g" "t" "a" "t" "g" "g" "g" "c" "t" "a" "c" "a" "g" "a" "a" "a" "c"
     [397] "c" "g" "t" "g" "c" "c" "a" "a" "a" "a" "g" "a" "c" "t" "t" "c" "t" "a"
     [415] "c" "a" "g" "a" "g" "t" "g" "a" "a" "c" "c" "c" "g" "a" "a" "a" "a" "t"
     [433] "c" "c" "t" "t" "c" "c" "t" "t" "g" "c" "a" "g" "g" "a" "a" "a" "c" "c"
     [451] "a" "g" "t" "c" "t" "c" "a" "g" "t" "g" "t" "c" "c" "a" "a" "c" "t" "c"
     [469] "t" "c" "t" "a" "a" "c" "c" "t" "t" "g" "g" "a" "a" "c" "t" "g" "t" "g"
     [487] "a" "g" "a" "a" "c" "t" "c" "t" "g" "a" "g" "g" "a" "c" "a" "a" "a" "g"
     [505] "c" "a" "g" "c" "g" "g" "a" "t" "a" "c" "a" "a" "c" "c" "t" "c" "a" "a"
     [523] "a" "a" "g" "a" "c" "g" "t" "c" "t" "g" "t" "c" "t" "a" "c" "a" "t" "t"
     [541] "g" "a" "a" "t" "t" "g" "g" "g" "a" "t" "c" "t" "g" "a" "t" "t" "c" "t"
     [559] "t" "c" "t" "g" "a" "a" "g" "a" "t" "a" "c" "c" "g" "t" "t" "a" "a" "t"
     [577] "a" "a" "g" "g" "c" "a" "a" "c" "t" "t" "a" "t" "t" "g" "c" "a" "g" "t"
     [595] "g" "t" "g" "g" "g" "a" "g" "a" "t" "c" "a" "a" "g" "a" "a" "t" "t" "g"
     [613] "t" "t" "a" "c" "a" "a" "a" "t" "c" "a" "c" "c" "c" "c" "t" "c" "a" "a"
     [631] "g" "g" "a" "a" "c" "c" "a" "g" "g" "g" "a" "t" "g" "a" "a" "a" "t" "c"
     [649] "a" "g" "t" "t" "t" "g" "g" "a" "t" "t" "c" "t" "g" "c" "a" "a" "a" "a"
     [667] "a" "a" "g" "g" "c" "t" "g" "c" "t" "t" "g" "t" "g" "a" "a" "t" "t" "t"
     [685] "t" "c" "t" "g" "a" "g" "a" "c" "g" "g" "a" "t" "g" "t" "a" "a" "c" "a"
     [703] "a" "a" "t" "a" "c" "t" "g" "a" "a" "c" "a" "t" "c" "a" "t" "c" "a" "a"
     [721] "c" "c" "c" "a" "g" "t" "a" "a" "t" "a" "a" "t" "g" "a" "t" "t" "t" "g"
     [739] "a" "a" "c" "a" "c" "c" "a" "c" "t" "g" "a" "g" "a" "a" "g" "c" "g" "t"
     [757] "g" "c" "a" "g" "c" "t" "g" "a" "g" "a" "g" "g" "c" "a" "t" "c" "c" "a"
     [775] "g" "a" "a" "a" "a" "g" "t" "a" "t" "c" "a" "g" "g" "g" "t" "a" "g" "t"
     [793] "t" "c" "t" "g" "t" "t" "t" "c" "a" "a" "a" "c" "t" "t" "g" "c" "a" "t"
     [811] "g" "t" "g" "g" "a" "g" "c" "c" "a" "t" "g" "t" "g" "g" "c" "a" "c" "a"
     [829] "a" "a" "t" "a" "c" "t" "c" "a" "t" "g" "c" "c" "a" "g" "c" "t" "c" "a"
     [847] "t" "t" "a" "c" "a" "g" "c" "a" "t" "g" "a" "g" "a" "a" "c" "a" "g" "c"
     [865] "a" "g" "t" "t" "t" "a" "t" "t" "a" "c" "t" "c" "a" "c" "t" "a" "a" "a"
     [883] "g" "a" "c" "a" "g" "a" "a" "t" "g" "a" "a" "t" "g" "t" "a" "g" "a" "a"
     [901] "a" "a" "g" "g" "c" "t" "g" "a" "a" "t" "t" "c" "t" "g" "t" "a" "a" "t"
     [919] "a" "a" "a" "a" "g" "c" "a" "a" "a" "c" "a" "g" "c" "c" "t" "g" "g" "c"
     [937] "t" "t" "a" "g" "c" "a" "a" "g" "g" "a" "g" "c" "c" "a" "a" "c" "a" "t"
     [955] "a" "a" "c" "a" "g" "a" "t" "g" "g" "g" "c" "t" "g" "g" "a" "a" "g" "t"
     [973] "a" "a" "g" "g" "a" "a" "a" "c" "a" "t" "g" "t" "a" "a" "t" "g" "a" "t"
     [991] "a" "g" "g" "c" "g" "g" "a" "c" "t" "c" "c" "c" "a" "g" "c" "a" "c" "a"
    [1009] "g" "a" "a" "a" "a" "a" "a" "a" "g" "g" "t" "a" "g" "a" "t" "c" "t" "g"
    [1027] "a" "a" "t" "g" "c" "t" "g" "a" "t" "c" "c" "c" "c" "t" "g" "t" "g" "t"
    [1045] "g" "a" "g" "a" "g" "a" "a" "a" "a" "g" "a" "a" "t" "g" "g" "a" "a" "t"
    [1063] "a" "a" "g" "c" "a" "g" "a" "a" "a" "c" "t" "g" "c" "c" "a" "t" "g" "c"
    [1081] "t" "c" "a" "g" "a" "g" "a" "a" "t" "c" "c" "t" "a" "g" "a" "g" "a" "t"
    [1099] "a" "c" "t" "g" "a" "a" "g" "a" "t" "g" "t" "t" "c" "c" "t" "t" "g" "g"
    [1117] "a" "t" "a" "a" "c" "a" "c" "t" "a" "a" "a" "t" "a" "g" "c" "a" "g" "c"
    [1135] "a" "t" "t" "c" "a" "g" "a" "a" "a" "g" "t" "t" "a" "a" "t" "g" "a" "g"
    [1153] "t" "g" "g" "t" "t" "t" "t" "c" "c" "a" "g" "a" "a" "g" "t" "g" "a" "t"
    [1171] "g" "a" "a" "c" "t" "g" "t" "t" "a" "g" "g" "t" "t" "c" "t" "g" "a" "t"
    [1189] "g" "a" "c" "t" "c" "a" "c" "a" "t" "g" "a" "t" "g" "g" "g" "g" "a" "g"
    [1207] "t" "c" "t" "g" "a" "a" "t" "c" "a" "a" "a" "t" "g" "c" "c" "a" "a" "a"
    [1225] "g" "t" "a" "g" "c" "t" "g" "a" "t" "g" "t" "a" "t" "t" "g" "g" "a" "c"
    [1243] "g" "t" "t" "c" "t" "a" "a" "a" "t" "g" "a" "g" "g" "t" "a" "g" "a" "t"
    [1261] "g" "a" "a" "t" "a" "t" "t" "c" "t" "g" "g" "t" "t" "c" "t" "t" "c" "a"
    [1279] "g" "a" "g" "a" "a" "a" "a" "t" "a" "g" "a" "c" "t" "t" "a" "c" "t" "g"
    [1297] "g" "c" "c" "a" "g" "t" "g" "a" "t" "c" "c" "t" "c" "a" "t" "g" "a" "g"
    [1315] "g" "c" "t" "t" "t" "a" "a" "t" "a" "t" "g" "t" "a" "a" "a" "a" "g" "t"
    [1333] "g" "a" "a" "a" "g" "a" "g" "t" "t" "c" "a" "c" "t" "c" "c" "a" "a" "a"
    [1351] "t" "c" "a" "g" "t" "a" "g" "a" "g" "a" "g" "t" "a" "a" "t" "a" "t" "t"
    [1369] "g" "a" "a" "g" "a" "c" "a" "a" "a" "a" "t" "a" "t" "t" "t" "g" "g" "g"
    [1387] "a" "a" "a" "a" "c" "c" "t" "a" "t" "c" "g" "g" "a" "a" "g" "a" "a" "g"
    [1405] "g" "c" "a" "a" "g" "c" "c" "t" "c" "c" "c" "c" "a" "a" "c" "t" "t" "a"
    [1423] "a" "g" "c" "c" "a" "t" "g" "t" "a" "a" "c" "t" "g" "a" "a" "a" "a" "t"
    [1441] "c" "t" "a" "a" "t" "t" "a" "t" "a" "g" "g" "a" "g" "c" "a" "t" "t" "t"
    [1459] "g" "t" "t" "a" "c" "t" "g" "a" "g" "c" "c" "a" "c" "a" "g" "a" "t" "a"
    [1477] "a" "t" "a" "c" "a" "a" "g" "a" "g" "c" "g" "t" "c" "c" "c" "c" "t" "c"
    [1495] "a" "c" "a" "a" "a" "t" "a" "a" "a" "t" "t" "a" "a" "a" "g" "c" "g" "t"
    [1513] "a" "a" "a" "a" "g" "g" "a" "g" "a" "c" "c" "t" "a" "c" "a" "t" "c" "a"
    [1531] "g" "g" "c" "c" "t" "t" "c" "a" "t" "c" "c" "t" "g" "a" "g" "g" "a" "t"
    [1549] "t" "t" "t" "a" "t" "c" "a" "a" "g" "a" "a" "a" "g" "c" "a" "g" "a" "t"
    [1567] "t" "t" "g" "g" "c" "a" "g" "t" "t" "c" "a" "a" "a" "a" "g" "a" "c" "t"
    [1585] "c" "c" "t" "g" "a" "a" "a" "t" "g" "a" "t" "a" "a" "a" "t" "c" "a" "g"
    [1603] "g" "g" "a" "a" "c" "t" "a" "a" "c" "c" "a" "a" "a" "c" "g" "g" "a" "g"
    [1621] "c" "a" "g" "a" "a" "t" "g" "g" "t" "c" "a" "a" "g" "t" "g" "a" "t" "g"
    [1639] "a" "a" "t" "a" "t" "t" "a" "c" "t" "a" "a" "t" "a" "g" "t" "g" "g" "t"
    [1657] "c" "a" "t" "g" "a" "g" "a" "a" "t" "a" "a" "a" "a" "c" "a" "a" "a" "a"
    [1675] "g" "g" "t" "g" "a" "t" "t" "c" "t" "a" "t" "t" "c" "a" "g" "a" "a" "t"
    [1693] "g" "a" "g" "a" "a" "a" "a" "a" "t" "c" "c" "t" "a" "a" "c" "c" "c" "a"
    [1711] "a" "t" "a" "g" "a" "a" "t" "c" "a" "c" "t" "c" "g" "a" "a" "a" "a" "a"
    [1729] "g" "a" "a" "t" "c" "t" "g" "c" "t" "t" "t" "c" "a" "a" "a" "a" "c" "g"
    [1747] "a" "a" "a" "g" "c" "t" "g" "a" "a" "c" "c" "t" "a" "t" "a" "a" "g" "c"
    [1765] "a" "g" "c" "a" "g" "t" "a" "t" "a" "a" "g" "c" "a" "a" "t" "a" "t" "g"
    [1783] "g" "a" "a" "c" "t" "c" "g" "a" "a" "t" "t" "a" "a" "a" "t" "a" "t" "c"
    [1801] "c" "a" "c" "a" "a" "t" "t" "c" "a" "a" "a" "a" "g" "c" "a" "c" "c" "t"
    [1819] "a" "a" "a" "a" "a" "g" "a" "a" "t" "a" "g" "g" "c" "t" "g" "a" "g" "g"
    [1837] "a" "g" "g" "a" "a" "g" "t" "c" "t" "t" "c" "t" "a" "c" "c" "a" "g" "g"
    [1855] "c" "a" "t" "a" "t" "t" "c" "a" "t" "g" "c" "g" "c" "t" "t" "g" "a" "a"
    [1873] "c" "t" "a" "g" "t" "a" "g" "t" "c" "a" "g" "t" "a" "g" "a" "a" "a" "t"
    [1891] "c" "t" "a" "a" "g" "c" "c" "c" "a" "c" "c" "t" "a" "a" "t" "t" "g" "t"
    [1909] "a" "c" "t" "g" "a" "a" "t" "t" "g" "c" "a" "a" "a" "t" "t" "g" "a" "t"
    [1927] "a" "g" "t" "t" "g" "t" "t" "c" "t" "a" "g" "c" "a" "g" "t" "g" "a" "a"
    [1945] "g" "a" "g" "a" "t" "a" "a" "a" "g" "a" "a" "a" "a" "a" "a" "a" "a" "g"
    [1963] "t" "a" "c" "a" "a" "c" "c" "a" "a" "a" "t" "g" "c" "c" "a" "g" "t" "c"
    [1981] "a" "g" "g" "c" "a" "c" "a" "g" "c" "a" "g" "a" "a" "a" "c" "c" "t" "a"
    [1999] "c" "a" "a" "c" "t" "c" "a" "t" "g" "g" "a" "a" "g" "g" "t" "a" "a" "a"
    [2017] "g" "a" "a" "c" "c" "t" "g" "c" "a" "a" "c" "t" "g" "g" "a" "g" "c" "c"
    [2035] "a" "a" "g" "a" "a" "g" "a" "g" "t" "a" "a" "c" "a" "a" "g" "c" "c" "a"
    [2053] "a" "a" "t" "g" "a" "a" "c" "a" "g" "a" "c" "a" "a" "g" "t" "a" "a" "a"
    [2071] "a" "g" "a" "c" "a" "t" "g" "a" "c" "a" "g" "c" "g" "a" "t" "a" "c" "t"
    [2089] "t" "t" "c" "c" "c" "a" "g" "a" "g" "c" "t" "g" "a" "a" "g" "t" "t" "a"
    [2107] "a" "c" "a" "a" "a" "t" "g" "c" "a" "c" "c" "t" "g" "g" "t" "t" "c" "t"
    [2125] "t" "t" "t" "a" "c" "t" "a" "a" "g" "t" "g" "t" "t" "c" "a" "a" "a" "t"
    [2143] "a" "c" "c" "a" "g" "t" "g" "a" "a" "c" "t" "t" "a" "a" "a" "g" "a" "a"
    [2161] "t" "t" "t" "g" "t" "c" "a" "a" "t" "c" "c" "t" "a" "g" "c" "c" "t" "t"
    [2179] "c" "c" "a" "a" "g" "a" "g" "a" "a" "g" "a" "a" "a" "a" "a" "g" "a" "a"
    [2197] "g" "a" "g" "a" "a" "a" "c" "t" "a" "g" "a" "a" "a" "c" "a" "g" "t" "t"
    [2215] "a" "a" "a" "g" "t" "g" "t" "c" "t" "a" "a" "t" "a" "a" "t" "g" "c" "t"
    [2233] "g" "a" "a" "g" "a" "c" "c" "c" "c" "a" "a" "a" "g" "a" "t" "c" "t" "c"
    [2251] "a" "t" "g" "t" "t" "a" "a" "g" "t" "g" "g" "a" "g" "a" "a" "a" "g" "g"
    [2269] "g" "t" "t" "t" "t" "g" "c" "a" "a" "a" "c" "t" "g" "a" "a" "a" "g" "a"
    [2287] "t" "c" "t" "g" "t" "a" "g" "a" "g" "a" "g" "t" "a" "g" "c" "a" "g" "t"
    [2305] "a" "t" "t" "t" "c" "a" "t" "t" "g" "g" "t" "a" "c" "c" "t" "g" "g" "t"
    [2323] "a" "c" "t" "g" "a" "t" "t" "a" "t" "g" "g" "c" "a" "c" "t" "c" "a" "g"
    [2341] "g" "a" "a" "a" "g" "t" "a" "t" "c" "t" "c" "g" "t" "t" "a" "c" "t" "g"
    [2359] "g" "a" "a" "g" "t" "t" "a" "g" "c" "a" "c" "t" "c" "t" "a" "g" "g" "g"
    [2377] "a" "a" "g" "g" "c" "a" "a" "a" "a" "a" "c" "a" "g" "a" "a" "c" "c" "a"
    [2395] "a" "a" "t" "a" "a" "a" "t" "g" "t" "g" "t" "g" "a" "g" "t" "c" "a" "g"
    [2413] "t" "g" "t" "g" "c" "a" "g" "c" "a" "t" "t" "t" "g" "a" "a" "a" "a" "c"
    [2431] "c" "c" "c" "a" "a" "g" "g" "g" "a" "c" "t" "a" "a" "t" "t" "c" "a" "t"
    [2449] "g" "g" "t" "t" "g" "t" "t" "c" "c" "a" "a" "a" "g" "a" "t" "a" "a" "t"
    [2467] "a" "g" "a" "a" "a" "t" "g" "a" "c" "a" "c" "a" "g" "a" "a" "g" "g" "c"
    [2485] "t" "t" "t" "a" "a" "g" "t" "a" "t" "c" "c" "a" "t" "t" "g" "g" "g" "a"
    [2503] "c" "a" "t" "g" "a" "a" "g" "t" "t" "a" "a" "c" "c" "a" "c" "a" "g" "t"
    [2521] "c" "g" "g" "g" "a" "a" "a" "c" "a" "a" "g" "c" "a" "t" "a" "g" "a" "a"
    [2539] "a" "t" "g" "g" "a" "a" "g" "a" "a" "a" "g" "t" "g" "a" "a" "c" "t" "t"
    [2557] "g" "a" "t" "g" "c" "t" "c" "a" "g" "t" "a" "t" "t" "t" "g" "c" "a" "g"
    [2575] "a" "a" "t" "a" "c" "a" "t" "t" "c" "a" "a" "g" "g" "t" "t" "t" "c" "a"
    [2593] "a" "a" "g" "c" "g" "c" "c" "a" "g" "t" "c" "a" "t" "t" "t" "g" "c" "t"
    [2611] "c" "c" "g" "t" "t" "t" "t" "c" "a" "a" "a" "t" "c" "c" "a" "g" "g" "a"
    [2629] "a" "a" "t" "g" "c" "a" "g" "a" "a" "g" "a" "g" "g" "a" "a" "t" "g" "t"
    [2647] "g" "c" "a" "a" "c" "a" "t" "t" "c" "t" "c" "t" "g" "c" "c" "c" "a" "c"
    [2665] "t" "c" "t" "g" "g" "g" "t" "c" "c" "t" "t" "a" "a" "a" "g" "a" "a" "a"
    [2683] "c" "a" "a" "a" "g" "t" "c" "c" "a" "a" "a" "a" "g" "t" "c" "a" "c" "t"
    [2701] "t" "t" "t" "g" "a" "a" "t" "g" "t" "g" "a" "a" "c" "a" "a" "a" "a" "g"
    [2719] "g" "a" "a" "g" "a" "a" "a" "a" "t" "c" "a" "a" "g" "g" "a" "a" "a" "g"
    [2737] "a" "a" "t" "g" "a" "g" "t" "c" "t" "a" "a" "t" "a" "t" "c" "a" "a" "g"
    [2755] "c" "c" "t" "g" "t" "a" "c" "a" "g" "a" "c" "a" "g" "t" "t" "a" "a" "t"
    [2773] "a" "t" "c" "a" "c" "t" "g" "c" "a" "g" "g" "c" "t" "t" "t" "c" "c" "t"
    [2791] "g" "t" "g" "g" "t" "t" "g" "g" "t" "c" "a" "g" "a" "a" "a" "g" "a" "t"
    [2809] "a" "a" "g" "c" "c" "a" "g" "t" "t" "g" "a" "t" "a" "a" "t" "g" "c" "c"
    [2827] "a" "a" "a" "t" "g" "t" "a" "g" "t" "a" "t" "c" "a" "a" "a" "g" "g" "a"
    [2845] "g" "g" "c" "t" "c" "t" "a" "g" "g" "t" "t" "t" "t" "g" "t" "c" "t" "a"
    [2863] "t" "c" "a" "t" "c" "t" "c" "a" "g" "t" "t" "c" "a" "g" "a" "g" "g" "c"
    [2881] "a" "a" "c" "g" "a" "a" "a" "c" "t" "g" "g" "a" "c" "t" "c" "a" "t" "t"
    [2899] "a" "c" "t" "c" "c" "a" "a" "a" "t" "a" "a" "a" "c" "a" "t" "g" "g" "a"
    [2917] "c" "t" "t" "t" "t" "a" "c" "a" "a" "a" "a" "c" "c" "c" "a" "t" "a" "t"
    [2935] "c" "g" "t" "a" "t" "a" "c" "c" "a" "c" "c" "a" "c" "t" "t" "t" "t" "t"
    [2953] "c" "c" "c" "a" "t" "c" "a" "a" "g" "t" "c" "a" "t" "t" "t" "g" "t" "t"
    [2971] "a" "a" "a" "a" "c" "t" "a" "a" "a" "t" "g" "t" "a" "a" "g" "a" "a" "a"
    [2989] "a" "a" "t" "c" "t" "g" "c" "t" "a" "g" "a" "g" "g" "a" "a" "a" "a" "c"
    [3007] "t" "t" "t" "g" "a" "g" "g" "a" "a" "c" "a" "t" "t" "c" "a" "a" "t" "g"
    [3025] "t" "c" "a" "c" "c" "t" "g" "a" "a" "a" "g" "a" "g" "a" "a" "a" "t" "g"
    [3043] "g" "g" "a" "a" "a" "t" "g" "a" "g" "a" "a" "c" "a" "t" "t" "c" "c" "a"
    [3061] "a" "g" "t" "a" "c" "a" "g" "t" "g" "a" "g" "c" "a" "c" "a" "a" "t" "t"
    [3079] "a" "g" "c" "c" "g" "t" "a" "a" "t" "a" "a" "c" "a" "t" "t" "a" "g" "a"
    [3097] "g" "a" "a" "a" "a" "t" "g" "t" "t" "t" "t" "t" "a" "a" "a" "g" "a" "a"
    [3115] "g" "c" "c" "a" "g" "c" "t" "c" "a" "a" "g" "c" "a" "a" "t" "a" "t" "t"
    [3133] "a" "a" "t" "g" "a" "a" "g" "t" "a" "g" "g" "t" "t" "c" "c" "a" "g" "t"
    [3151] "a" "c" "t" "a" "a" "t" "g" "a" "a" "g" "t" "g" "g" "g" "c" "t" "c" "c"
    [3169] "a" "g" "t" "a" "t" "t" "a" "a" "t" "g" "a" "a" "a" "t" "a" "g" "g" "t"
    [3187] "t" "c" "c" "a" "g" "t" "g" "a" "t" "g" "a" "a" "a" "a" "c" "a" "t" "t"
    [3205] "c" "a" "a" "g" "c" "a" "g" "a" "a" "c" "t" "a" "g" "g" "t" "a" "g" "a"
    [3223] "a" "a" "c" "a" "g" "a" "g" "g" "g" "c" "c" "a" "a" "a" "a" "t" "t" "g"
    [3241] "a" "a" "t" "g" "c" "t" "a" "t" "g" "c" "t" "t" "a" "g" "a" "t" "t" "a"
    [3259] "g" "g" "g" "g" "t" "t" "t" "t" "g" "c" "a" "a" "c" "c" "t" "g" "a" "g"
    [3277] "g" "t" "c" "t" "a" "t" "a" "a" "a" "c" "a" "a" "a" "g" "t" "c" "t" "t"
    [3295] "c" "c" "t" "g" "g" "a" "a" "g" "t" "a" "a" "t" "t" "g" "t" "a" "a" "g"
    [3313] "c" "a" "t" "c" "c" "t" "g" "a" "a" "a" "t" "a" "a" "a" "a" "a" "a" "g"
    [3331] "c" "a" "a" "g" "a" "a" "t" "a" "t" "g" "a" "a" "g" "a" "a" "g" "t" "a"
    [3349] "g" "t" "t" "c" "a" "g" "a" "c" "t" "g" "t" "t" "a" "a" "t" "a" "c" "a"
    [3367] "g" "a" "t" "t" "t" "c" "t" "c" "t" "c" "c" "a" "t" "a" "t" "c" "t" "g"
    [3385] "a" "t" "t" "t" "c" "a" "g" "a" "t" "a" "a" "c" "t" "t" "a" "g" "a" "a"
    [3403] "c" "a" "g" "c" "c" "t" "a" "t" "g" "g" "g" "a" "a" "g" "t" "a" "g" "t"
    [3421] "c" "a" "t" "g" "c" "a" "t" "c" "t" "c" "a" "g" "g" "t" "t" "t" "g" "t"
    [3439] "t" "c" "t" "g" "a" "g" "a" "c" "a" "c" "c" "t" "g" "a" "t" "g" "a" "c"
    [3457] "c" "t" "g" "t" "t" "a" "g" "a" "t" "g" "a" "t" "g" "g" "t" "g" "a" "a"
    [3475] "a" "t" "a" "a" "a" "g" "g" "a" "a" "g" "a" "t" "a" "c" "t" "a" "g" "t"
    [3493] "t" "t" "t" "g" "c" "t" "g" "a" "a" "a" "a" "t" "g" "a" "c" "a" "t" "t"
    [3511] "a" "a" "g" "g" "a" "a" "a" "g" "t" "t" "c" "t" "g" "c" "t" "g" "t" "t"
    [3529] "t" "t" "t" "a" "g" "c" "a" "a" "a" "a" "g" "c" "g" "t" "c" "c" "a" "g"
    [3547] "a" "a" "a" "g" "g" "a" "g" "a" "g" "c" "t" "t" "a" "g" "c" "a" "g" "g"
    [3565] "a" "g" "t" "c" "c" "t" "a" "g" "c" "c" "c" "t" "t" "t" "c" "a" "c" "c"
    [3583] "c" "a" "t" "a" "c" "a" "c" "a" "t" "t" "t" "g" "g" "c" "t" "c" "a" "g"
    [3601] "g" "g" "t" "t" "a" "c" "c" "g" "a" "a" "g" "a" "g" "g" "g" "g" "c" "c"
    [3619] "a" "a" "g" "a" "a" "a" "t" "t" "a" "g" "a" "g" "t" "c" "c" "t" "c" "a"
    [3637] "g" "a" "a" "g" "a" "g" "a" "a" "c" "t" "t" "a" "t" "c" "t" "a" "g" "t"
    [3655] "g" "a" "g" "g" "a" "t" "g" "a" "a" "g" "a" "g" "c" "t" "t" "c" "c" "c"
    [3673] "t" "g" "c" "t" "t" "c" "c" "a" "a" "c" "a" "c" "t" "t" "g" "t" "t" "a"
    [3691] "t" "t" "t" "g" "g" "t" "a" "a" "a" "g" "t" "a" "a" "a" "c" "a" "a" "t"
    [3709] "a" "t" "a" "c" "c" "t" "t" "c" "t" "c" "a" "g" "t" "c" "t" "a" "c" "t"
    [3727] "a" "g" "g" "c" "a" "t" "a" "g" "c" "a" "c" "c" "g" "t" "t" "g" "c" "t"
    [3745] "a" "c" "c" "g" "a" "g" "t" "g" "t" "c" "t" "g" "t" "c" "t" "a" "a" "g"
    [3763] "a" "a" "c" "a" "c" "a" "g" "a" "g" "g" "a" "g" "a" "a" "t" "t" "t" "a"
    [3781] "t" "t" "a" "t" "c" "a" "t" "t" "g" "a" "a" "g" "a" "a" "t" "a" "g" "c"
    [3799] "t" "t" "a" "a" "a" "t" "g" "a" "c" "t" "g" "c" "a" "g" "t" "a" "a" "c"
    [3817] "c" "a" "g" "g" "t" "a" "a" "t" "a" "t" "t" "g" "g" "c" "a" "a" "a" "g"
    [3835] "g" "c" "a" "t" "c" "t" "c" "a" "g" "g" "a" "a" "c" "a" "t" "c" "a" "c"
    [3853] "c" "t" "t" "a" "g" "t" "g" "a" "g" "g" "a" "a" "a" "c" "a" "a" "a" "a"
    [3871] "t" "g" "t" "t" "c" "t" "g" "c" "t" "a" "g" "c" "t" "t" "g" "t" "t" "t"
    [3889] "t" "c" "t" "t" "c" "a" "c" "a" "g" "t" "g" "c" "a" "g" "t" "g" "a" "a"
    [3907] "t" "t" "g" "g" "a" "a" "g" "a" "c" "t" "t" "g" "a" "c" "t" "g" "c" "a"
    [3925] "a" "a" "t" "a" "c" "a" "a" "a" "c" "a" "c" "c" "c" "a" "g" "g" "a" "t"
    [3943] "c" "c" "t" "t" "t" "c" "t" "t" "g" "a" "t" "t" "g" "g" "t" "t" "c" "t"
    [3961] "t" "c" "c" "a" "a" "a" "c" "a" "a" "a" "t" "g" "a" "g" "g" "c" "a" "t"
    [3979] "c" "a" "g" "t" "c" "t" "g" "a" "a" "a" "g" "c" "c" "a" "g" "g" "g" "a"
    [3997] "g" "t" "t" "g" "g" "t" "c" "t" "g" "a" "g" "t" "g" "a" "c" "a" "a" "g"
    [4015] "g" "a" "a" "t" "t" "g" "g" "t" "t" "t" "c" "a" "g" "a" "t" "g" "a" "t"
    [4033] "g" "a" "a" "g" "a" "a" "a" "g" "a" "g" "g" "a" "a" "c" "g" "g" "g" "c"
    [4051] "t" "t" "g" "g" "a" "a" "g" "a" "a" "a" "a" "t" "a" "a" "t" "c" "a" "a"
    [4069] "g" "a" "a" "g" "a" "g" "c" "a" "a" "a" "g" "c" "a" "t" "g" "g" "a" "t"
    [4087] "t" "c" "a" "a" "a" "c" "t" "t" "a" "g" "g" "t" "g" "a" "a" "g" "c" "a"
    [4105] "g" "c" "a" "t" "c" "t" "g" "g" "g" "t" "g" "t" "g" "a" "g" "a" "g" "t"
    [4123] "g" "a" "a" "a" "c" "a" "a" "g" "c" "g" "t" "c" "t" "c" "t" "g" "a" "a"
    [4141] "g" "a" "c" "t" "g" "c" "t" "c" "a" "g" "g" "g" "c" "t" "a" "t" "c" "c"
    [4159] "t" "c" "t" "c" "a" "g" "a" "g" "t" "g" "a" "c" "a" "t" "t" "t" "t" "a"
    [4177] "a" "c" "c" "a" "c" "t" "c" "a" "g" "c" "a" "g" "a" "g" "g" "g" "a" "t"
    [4195] "a" "c" "c" "a" "t" "g" "c" "a" "a" "c" "a" "t" "a" "a" "c" "c" "t" "g"
    [4213] "a" "t" "a" "a" "a" "g" "c" "t" "c" "c" "a" "g" "c" "a" "g" "g" "a" "a"
    [4231] "a" "t" "g" "g" "c" "t" "g" "a" "a" "c" "t" "a" "g" "a" "a" "g" "c" "t"
    [4249] "g" "t" "g" "t" "t" "a" "g" "a" "a" "c" "a" "g" "c" "a" "t" "g" "g" "g"
    [4267] "a" "g" "c" "c" "a" "g" "c" "c" "t" "t" "c" "t" "a" "a" "c" "a" "g" "c"
    [4285] "t" "a" "c" "c" "c" "t" "t" "c" "c" "a" "t" "c" "a" "t" "a" "a" "g" "t"
    [4303] "g" "a" "c" "t" "c" "t" "t" "c" "t" "g" "c" "c" "c" "t" "t" "g" "a" "g"
    [4321] "g" "a" "c" "c" "t" "g" "c" "g" "a" "a" "a" "t" "c" "c" "a" "g" "a" "a"
    [4339] "c" "a" "a" "a" "g" "c" "a" "c" "a" "t" "c" "a" "g" "a" "a" "a" "a" "a"
    [4357] "g" "c" "a" "g" "t" "a" "t" "t" "a" "a" "c" "t" "t" "c" "a" "c" "a" "g"
    [4375] "a" "a" "a" "a" "g" "t" "a" "g" "t" "g" "a" "a" "t" "a" "c" "c" "c" "t"
    [4393] "a" "t" "a" "a" "g" "c" "c" "a" "g" "a" "a" "t" "c" "c" "a" "g" "a" "a"
    [4411] "g" "g" "c" "c" "t" "t" "t" "c" "t" "g" "c" "t" "g" "a" "c" "a" "a" "g"
    [4429] "t" "t" "t" "g" "a" "g" "g" "t" "g" "t" "c" "t" "g" "c" "a" "g" "a" "t"
    [4447] "a" "g" "t" "t" "c" "t" "a" "c" "c" "a" "g" "t" "a" "a" "a" "a" "a" "t"
    [4465] "a" "a" "a" "g" "a" "a" "c" "c" "a" "g" "g" "a" "g" "t" "g" "g" "a" "a"
    [4483] "a" "g" "g" "t" "c" "a" "t" "c" "c" "c" "c" "t" "t" "c" "t" "a" "a" "a"
    [4501] "t" "g" "c" "c" "c" "a" "t" "c" "a" "t" "t" "a" "g" "a" "t" "g" "a" "t"
    [4519] "a" "g" "g" "t" "g" "g" "t" "a" "c" "a" "t" "g" "c" "a" "c" "a" "g" "t"
    [4537] "t" "g" "c" "t" "c" "t" "g" "g" "g" "a" "g" "t" "c" "t" "t" "c" "a" "g"
    [4555] "a" "a" "t" "a" "g" "a" "a" "a" "c" "t" "a" "c" "c" "c" "a" "t" "c" "t"
    [4573] "c" "a" "a" "g" "a" "g" "g" "a" "g" "c" "t" "c" "a" "t" "t" "a" "a" "g"
    [4591] "g" "t" "t" "g" "t" "t" "g" "a" "t" "g" "t" "g" "g" "a" "g" "g" "a" "g"
    [4609] "c" "a" "a" "c" "a" "g" "c" "t" "g" "g" "a" "a" "g" "a" "g" "t" "c" "t"
    [4627] "g" "g" "g" "c" "c" "a" "c" "a" "c" "g" "a" "t" "t" "t" "g" "a" "c" "g"
    [4645] "g" "a" "a" "a" "c" "a" "t" "c" "t" "t" "a" "c" "t" "t" "g" "c" "c" "a"
    [4663] "a" "g" "g" "c" "a" "a" "g" "a" "t" "c" "t" "a" "g" "a" "g" "g" "g" "a"
    [4681] "a" "c" "c" "c" "c" "t" "t" "a" "c" "c" "t" "g" "g" "a" "a" "t" "c" "t"
    [4699] "g" "g" "a" "a" "t" "c" "a" "g" "c" "c" "t" "c" "t" "t" "c" "t" "c" "t"
    [4717] "g" "a" "t" "g" "a" "c" "c" "c" "t" "g" "a" "a" "t" "c" "t" "g" "a" "t"
    [4735] "c" "c" "t" "t" "c" "t" "g" "a" "a" "g" "a" "c" "a" "g" "a" "g" "c" "c"
    [4753] "c" "c" "a" "g" "a" "g" "t" "c" "a" "g" "c" "t" "c" "g" "t" "g" "t" "t"
    [4771] "g" "g" "c" "a" "a" "c" "a" "t" "a" "c" "c" "a" "t" "c" "t" "t" "c" "a"
    [4789] "a" "c" "c" "t" "c" "t" "g" "c" "a" "t" "t" "g" "a" "a" "a" "g" "t" "t"
    [4807] "c" "c" "c" "c" "a" "a" "t" "t" "g" "a" "a" "a" "g" "t" "t" "g" "c" "a"
    [4825] "g" "a" "a" "t" "c" "t" "g" "c" "c" "c" "a" "g" "a" "g" "t" "c" "c" "a"
    [4843] "g" "c" "t" "g" "c" "t" "g" "c" "t" "c" "a" "t" "a" "c" "t" "a" "c" "t"
    [4861] "g" "a" "t" "a" "c" "t" "g" "c" "t" "g" "g" "g" "t" "a" "t" "a" "a" "t"
    [4879] "g" "c" "a" "a" "t" "g" "g" "a" "a" "g" "a" "a" "a" "g" "t" "g" "t" "g"
    [4897] "a" "g" "c" "a" "g" "g" "g" "a" "g" "a" "a" "g" "c" "c" "a" "g" "a" "a"
    [4915] "t" "t" "g" "a" "c" "a" "g" "c" "t" "t" "c" "a" "a" "c" "a" "g" "a" "a"
    [4933] "a" "g" "g" "g" "t" "c" "a" "a" "c" "a" "a" "a" "a" "g" "a" "a" "t" "g"
    [4951] "t" "c" "c" "a" "t" "g" "g" "t" "g" "g" "t" "g" "t" "c" "t" "g" "g" "c"
    [4969] "c" "t" "g" "a" "c" "c" "c" "c" "a" "g" "a" "a" "g" "a" "a" "t" "t" "t"
    [4987] "a" "t" "g" "c" "t" "c" "g" "t" "g" "t" "a" "c" "a" "a" "g" "t" "t" "t"
    [5005] "g" "c" "c" "a" "g" "a" "a" "a" "a" "c" "a" "c" "c" "a" "c" "a" "t" "c"
    [5023] "a" "c" "t" "t" "t" "a" "a" "c" "t" "a" "a" "t" "c" "t" "a" "a" "t" "t"
    [5041] "a" "c" "t" "g" "a" "a" "g" "a" "g" "a" "c" "t" "a" "c" "t" "c" "a" "t"
    [5059] "g" "t" "t" "g" "t" "t" "a" "t" "g" "a" "a" "a" "a" "c" "a" "g" "a" "t"
    [5077] "g" "c" "t" "g" "a" "g" "t" "t" "t" "g" "t" "g" "t" "g" "t" "g" "a" "a"
    [5095] "c" "g" "g" "a" "c" "a" "c" "t" "g" "a" "a" "a" "t" "a" "t" "t" "t" "t"
    [5113] "c" "t" "a" "g" "g" "a" "a" "t" "t" "g" "c" "g" "g" "g" "a" "g" "g" "a"
    [5131] "a" "a" "a" "t" "g" "g" "g" "t" "a" "g" "t" "t" "a" "g" "c" "t" "a" "t"
    [5149] "t" "t" "c" "t" "g" "g" "g" "t" "g" "a" "c" "c" "c" "a" "g" "t" "c" "t"
    [5167] "a" "t" "t" "a" "a" "a" "g" "a" "a" "a" "g" "a" "a" "a" "a" "a" "t" "g"
    [5185] "c" "t" "g" "a" "a" "t" "g" "a" "g" "c" "a" "t" "g" "a" "t" "t" "t" "t"
    [5203] "g" "a" "a" "g" "t" "c" "a" "g" "a" "g" "g" "a" "g" "a" "t" "g" "t" "g"
    [5221] "g" "t" "c" "a" "a" "t" "g" "g" "a" "a" "g" "a" "a" "a" "c" "c" "a" "c"
    [5239] "c" "a" "a" "g" "g" "t" "c" "c" "a" "a" "a" "g" "c" "g" "a" "g" "c" "a"
    [5257] "a" "g" "a" "g" "a" "a" "t" "c" "c" "c" "a" "g" "g" "a" "c" "a" "g" "a"
    [5275] "a" "a" "g" "a" "t" "c" "t" "t" "c" "a" "g" "g" "g" "g" "g" "c" "t" "a"
    [5293] "g" "a" "a" "a" "t" "c" "t" "g" "t" "t" "g" "c" "t" "a" "t" "g" "g" "g"
    [5311] "c" "c" "c" "t" "t" "c" "a" "c" "c" "a" "a" "c" "a" "t" "g" "c" "c" "c"
    [5329] "a" "c" "a" "g" "a" "t" "c" "a" "a" "c" "t" "g" "g" "a" "a" "t" "g" "g"
    [5347] "a" "t" "g" "g" "t" "a" "c" "a" "g" "c" "t" "g" "t" "g" "t" "g" "g" "t"
    [5365] "g" "c" "t" "t" "c" "t" "g" "t" "g" "g" "t" "g" "a" "a" "g" "g" "a" "g"
    [5383] "c" "t" "t" "t" "c" "a" "t" "c" "a" "t" "t" "c" "a" "c" "c" "c" "t" "t"
    [5401] "g" "g" "c" "a" "c" "a" "g" "g" "t" "g" "t" "c" "c" "a" "c" "c" "c" "a"
    [5419] "a" "t" "t" "g" "t" "g" "g" "t" "t" "g" "t" "g" "c" "a" "g" "c" "c" "a"
    [5437] "g" "a" "t" "g" "c" "c" "t" "g" "g" "a" "c" "a" "g" "a" "g" "g" "a" "c"
    [5455] "a" "a" "t" "g" "g" "c" "t" "t" "c" "c" "a" "t" "g" "c" "a" "a" "t" "t"
    [5473] "g" "g" "g" "c" "a" "g" "a" "t" "g" "t" "g" "t" "g" "a" "g" "g" "c" "a"
    [5491] "c" "c" "t" "g" "t" "g" "g" "t" "g" "a" "c" "c" "c" "g" "a" "g" "a" "g"
    [5509] "t" "g" "g" "g" "t" "g" "t" "t" "g" "g" "a" "c" "a" "g" "t" "g" "t" "a"
    [5527] "g" "c" "a" "c" "t" "c" "t" "a" "c" "c" "a" "g" "t" "g" "c" "c" "a" "g"
    [5545] "g" "a" "g" "c" "t" "g" "g" "a" "c" "a" "c" "c" "t" "a" "c" "c" "t" "g"
    [5563] "a" "t" "a" "c" "c" "c" "c" "a" "g" "a" "t" "c" "c" "c" "c" "c" "a" "c"
    [5581] "a" "g" "c" "c" "a" "c" "t" "a" "c" "t" "g" "a"
    attr(,"name")
    [1] "AY273801.BRCA1"
    attr(,"Annot")
    [1] ">AY273801.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $AY304547.BRCA1
       [1] "c" "t" "g" "c" "t" "t" "g" "t" "g" "a" "a" "t" "t" "t" "t" "c" "t" "g"
      [19] "a" "g" "a" "c" "g" "g" "a" "t" "g" "t" "a" "a" "c" "a" "a" "a" "t" "a"
      [37] "c" "t" "g" "a" "a" "c" "a" "t" "c" "a" "t" "c" "a" "a" "c" "c" "c" "a"
      [55] "g" "t" "a" "a" "t" "a" "a" "t" "g" "a" "t" "t" "t" "g" "a" "a" "c" "a"
      [73] "c" "c" "a" "c" "t" "g" "a" "g" "a" "a" "g" "c" "g" "t" "g" "c" "a" "g"
      [91] "c" "t" "g" "a" "g" "a" "g" "g" "c" "a" "t" "c" "c" "a" "g" "a" "a" "a"
     [109] "a" "g" "t" "a" "t" "c" "a" "g" "g" "g" "t" "a" "g" "t" "t" "c" "t" "g"
     [127] "t" "t" "t" "c" "a" "a" "a" "c" "t" "t" "g" "c" "a" "t" "g" "t" "g" "g"
     [145] "a" "g" "c" "c" "a" "t" "g" "t" "g" "g" "c" "a" "c" "a" "a" "a" "t" "a"
     [163] "c" "t" "c" "a" "t" "g" "c" "c" "a" "g" "c" "t" "c" "a" "t" "t" "a" "c"
     [181] "a" "g" "c" "a" "t" "g" "a" "g" "a" "a" "c" "a" "g" "c" "a" "g" "t" "t"
     [199] "t" "a" "t" "t" "a" "c" "t" "c" "a" "c" "t" "a" "a" "a" "g" "a" "c" "a"
     [217] "g" "a" "a" "t" "g" "a" "a" "t" "g" "t" "a" "g" "a" "a" "a" "a" "g" "g"
     [235] "c" "t" "g" "a" "a" "t" "t" "c" "t" "g" "t" "a" "a" "t" "a" "a" "a" "a"
     [253] "g" "c" "a" "a" "a" "c" "a" "g" "c" "c" "t" "g" "g" "c" "t" "t" "a" "g"
     [271] "c" "a" "a" "g" "g" "a" "g" "c" "c" "a" "a" "c" "a" "t" "a" "a" "c" "a"
     [289] "g" "a" "t" "g" "g" "g" "c" "t" "g" "g" "a" "a" "g" "t" "a" "a" "g" "g"
     [307] "a" "a" "a" "c" "a" "t" "g" "t" "a" "a" "t" "g" "a" "t" "a" "g" "g" "c"
     [325] "g" "g" "a" "c" "t" "c" "c" "c" "a" "g" "c" "a" "c" "a" "g" "a" "a" "a"
     [343] "a" "a" "a" "a" "g" "g" "t" "a" "g" "a" "t" "c" "t" "g" "a" "a" "t" "g"
     [361] "c" "t" "g" "a" "t" "c" "c" "c" "c" "t" "g" "t" "g" "t" "g" "a" "g" "a"
     [379] "g" "a" "a" "a" "a" "g" "a" "a" "t" "g" "g" "a" "a" "t" "a" "a" "g" "c"
     [397] "a" "g" "a" "a" "a" "c" "t" "g" "c" "c" "a" "t" "g" "c" "t" "c" "a" "g"
     [415] "a" "g" "a" "a" "t" "c" "c" "t" "a" "g" "a" "g" "a" "t" "a" "c" "t" "g"
     [433] "a" "a" "g" "a" "t" "g" "t" "t" "c" "c" "t" "t" "g" "g" "a" "t" "a" "a"
     [451] "c" "a" "c" "t" "a" "a" "a" "t" "a" "g" "c" "a" "g" "c" "a" "t" "t" "c"
     [469] "a" "g" "a" "a" "a" "g" "t" "t" "a" "a" "t" "g" "a" "g" "t" "g" "g" "t"
     [487] "t" "t" "t" "c" "c" "a" "g" "a" "a" "g" "t" "g" "a" "t" "g" "a" "a" "c"
     [505] "t" "g" "t" "t" "a" "g" "g" "t" "t" "c" "t" "g" "a" "t" "g" "a" "c" "t"
     [523] "c" "a" "c" "a" "t" "g" "a" "t" "g" "g" "g" "g" "a" "g" "t" "c" "t" "g"
     [541] "a" "a" "t" "c" "a" "a" "a" "t" "g" "c" "c" "a" "a" "a" "g" "t" "a" "g"
     [559] "c" "t" "g" "a" "t" "g" "t" "a" "t" "t" "g" "g" "a" "c" "g" "t" "t" "c"
     [577] "t" "a" "a" "a" "t" "g" "a" "g" "g" "t" "a" "g" "a" "t" "g" "a" "a" "t"
     [595] "a" "t" "t" "c" "t" "g" "g" "t" "t" "c" "t" "t" "c" "a" "g" "a" "g" "a"
     [613] "a" "a" "a" "t" "a" "g" "a" "c" "t" "t" "a" "c" "t" "g" "g" "c" "c" "a"
     [631] "g" "t" "g" "a" "t" "c" "c" "t" "c" "a" "t" "g" "a" "g" "g" "c" "t" "t"
     [649] "t" "a" "a" "t" "a" "t" "g" "t" "a" "a" "a" "a" "g" "t" "g" "a" "a" "a"
     [667] "g" "a" "g" "t" "t" "c" "a" "c" "t" "c" "c" "a" "a" "a" "t" "c" "a" "g"
     [685] "t" "a" "g" "a" "g" "a" "g" "t" "a" "a" "t" "a" "t" "t" "g" "a" "a" "g"
     [703] "a" "c" "a" "a" "a" "a" "t" "a" "t" "t" "t" "g" "g" "g" "a" "a" "a" "a"
     [721] "c" "c" "t" "a" "t" "c" "g" "g" "a" "a" "g" "a" "a" "g" "g" "c" "a" "a"
     [739] "g" "c" "c" "t" "c" "c" "c" "c" "a" "a" "c" "t" "t" "a" "a" "g" "c" "c"
     [757] "a" "t" "g" "t" "a" "a" "c" "t" "g" "a" "a" "a" "a" "t" "c" "t" "a" "a"
     [775] "t" "t" "a" "t" "a" "g" "g" "a" "g" "c" "a" "t" "t" "t" "g" "t" "t" "a"
     [793] "c" "t" "g" "a" "g" "c" "c" "a" "c" "a" "g" "a" "t" "a" "a" "t" "a" "c"
     [811] "a" "a" "g" "a" "g" "c" "g" "t" "c" "c" "c" "c" "t" "c" "a" "c" "a" "a"
     [829] "a" "t" "a" "a" "a" "t" "t" "a" "a" "a" "g" "c" "g" "t" "a" "a" "a" "a"
     [847] "g" "g" "a" "g" "a" "c" "c" "t" "a" "c" "a" "t" "c" "a" "g" "g" "c" "c"
     [865] "t" "t" "c" "a" "t" "c" "c" "t" "g" "a" "g" "g" "a" "t" "t" "t" "t" "a"
     [883] "t" "c" "a" "a" "g" "a" "a" "a" "g" "c" "a" "g" "a" "t" "t" "t" "g" "g"
     [901] "c" "a" "g" "t" "t" "c" "a" "a" "a" "a" "g" "a" "c" "t" "c" "c" "t" "g"
     [919] "a" "a" "a" "t" "g" "a" "t" "a" "a" "a" "t" "c" "a" "g" "g" "g" "a" "a"
     [937] "c" "t" "a" "a" "c" "c" "a" "a" "a" "c" "g" "g" "a" "g" "c" "a" "g" "a"
     [955] "a" "t" "g" "g" "t" "c" "a" "a" "g" "t" "g" "a" "t" "g" "a" "a" "t" "a"
     [973] "t" "t" "a" "c" "t" "a" "a" "t" "a" "g" "t" "g" "g" "t" "c" "a" "t" "g"
     [991] "a" "g" "a" "a" "t" "a" "a" "a" "a" "c" "a" "a" "a" "a" "g" "g" "t" "g"
    [1009] "a" "t" "t" "c" "t" "a" "t" "t" "c" "a" "g" "a" "a" "t" "g" "a" "g" "a"
    [1027] "a" "a" "a" "a" "t" "c" "c" "t" "a" "a" "c" "c" "c" "a" "a" "t" "a" "g"
    [1045] "a" "a" "t" "c" "a" "c" "t" "c" "g" "a" "a" "a" "a" "a" "g" "a" "a" "t"
    [1063] "c" "t" "g" "c" "t" "t" "t" "c" "a" "a" "a" "a" "c" "g" "a" "a" "a" "g"
    [1081] "c" "t" "g" "a" "a" "c" "c" "t" "a" "t" "a" "a" "g" "c" "a" "g" "c" "a"
    [1099] "g" "t" "a" "t" "a" "a" "g" "c" "a" "a" "t" "a" "t" "g" "g" "a" "a" "c"
    [1117] "t" "c" "g" "a" "a" "t" "t" "a" "a" "a" "t" "a" "t" "c" "c" "a" "c" "a"
    [1135] "a" "t" "t" "c" "a" "a" "a" "a" "g" "c" "a" "c" "c" "t" "a" "a" "a" "a"
    [1153] "a" "g" "a" "a" "t" "a" "g" "g" "c" "t" "g" "a" "g" "g" "a" "g" "g" "a"
    [1171] "a" "g" "t" "c" "t" "t" "c" "t" "a" "c" "c" "a" "g" "g" "c" "a" "t" "a"
    [1189] "t" "t" "c" "a" "t" "g" "c" "g" "c" "t" "t" "g" "a" "a" "c" "t" "a" "g"
    [1207] "t" "a" "g" "t" "c" "a" "g" "t" "a" "g" "a" "a" "a" "t" "c" "t" "a" "a"
    [1225] "g" "c" "c" "c" "a" "c" "c" "t" "a" "a" "t" "t" "g" "t" "a" "c" "t" "g"
    [1243] "a" "a" "t" "t" "g" "c" "a" "a" "a" "t" "t" "g" "a" "t" "a" "g" "t" "t"
    [1261] "g" "t" "t" "c" "t" "a" "g" "c" "a" "g" "t" "g" "a" "a" "g" "a" "g" "a"
    [1279] "t" "a" "a" "a" "g" "a" "a" "a" "a" "a" "a" "a" "a" "g" "t" "a" "c" "a"
    [1297] "a" "c" "c" "a" "a" "a" "t" "g" "c" "c" "a" "g" "t" "c" "a" "g" "g" "c"
    [1315] "a" "c" "a" "g" "c" "a" "g" "a" "a" "a" "c" "c" "t" "a" "c" "a" "a" "c"
    [1333] "t" "c" "a" "t" "g" "g" "a" "a" "g" "g" "t" "a" "a" "a" "g" "a" "a" "c"
    [1351] "c" "t" "g" "c" "a" "a" "c" "t" "g" "g" "a" "g" "c" "c" "a" "a" "g" "a"
    [1369] "a" "g" "a" "g" "t" "a" "a" "c" "a" "a" "g" "c" "c" "a" "a" "a" "t" "g"
    [1387] "a" "a" "c" "a" "g" "a" "c" "a" "a" "g" "t" "a" "a" "a" "a" "g" "a" "c"
    [1405] "a" "t" "g" "a" "c" "a" "g" "t" "g" "a" "t" "a" "c" "t" "t" "t" "c" "c"
    [1423] "c" "a" "g" "a" "g" "c" "t" "g" "a" "a" "g" "t" "t" "a" "a" "c" "a" "a"
    [1441] "a" "t" "g" "c" "a" "c" "c" "t" "g" "g" "t" "t" "c" "t" "t" "t" "t" "a"
    [1459] "c" "t" "a" "a" "g" "t" "g" "t" "t" "c" "a" "a" "a" "t" "a" "c" "c" "a"
    [1477] "g" "t" "g" "a" "a" "c" "t" "t" "a" "a" "a" "g" "a" "a" "t" "t" "t" "g"
    [1495] "t" "c" "a" "a" "t" "c" "c" "t" "a" "g" "c" "c" "t" "t" "c" "c" "a" "a"
    [1513] "g" "a" "g" "a" "a" "g" "a" "a" "a" "a" "a" "g" "a" "a" "g" "a" "g" "a"
    [1531] "a" "a" "c" "t" "a" "g" "a" "a" "a" "c" "a" "g" "t" "t" "a" "a" "a" "g"
    [1549] "t" "g" "t" "c" "t" "a" "a" "t" "a" "a" "t" "g" "c" "t" "g" "a" "a" "g"
    [1567] "a" "c" "c" "c" "c" "a" "a" "a" "g" "a" "t" "c" "t" "c" "a" "t" "g" "t"
    [1585] "t" "a" "a" "g" "t" "g" "g" "a" "g" "a" "a" "a" "g" "g" "g" "t" "t" "t"
    [1603] "t" "g" "c" "a" "a" "a" "c" "t" "g" "a" "a" "a" "g" "a" "t" "c" "t" "g"
    [1621] "t" "a" "g" "a" "g" "a" "g" "t" "a" "g" "c" "a" "g" "t" "a" "t" "t" "t"
    [1639] "c" "a" "c" "t" "g" "g" "t" "a" "c" "c" "t" "g" "g" "t" "a" "c" "t" "g"
    [1657] "a" "t" "t" "a" "t" "g" "g" "c" "a" "c" "t" "c" "a" "g" "g" "a" "a" "a"
    [1675] "g" "t" "a" "t" "c" "t" "c" "g" "t" "t" "a" "c" "t" "g" "g" "a" "a" "g"
    [1693] "t" "t" "a" "g" "c" "a" "c" "t" "c" "t" "a" "g" "g" "g" "a" "a" "g" "g"
    [1711] "c" "a" "a" "a" "a" "a" "c" "a" "g" "a" "a" "c" "c" "a" "a" "a" "t" "a"
    [1729] "a" "a" "t" "g" "t" "g" "t" "g" "a" "g" "t" "c" "a" "g" "t" "g" "t" "g"
    [1747] "c" "a" "g" "c" "a" "t" "t" "t" "g" "a" "a" "a" "a" "c" "c" "c" "c" "a"
    [1765] "a" "g" "g" "g" "a" "c" "t" "a" "a" "t" "t" "c" "a" "t" "g" "g" "t" "t"
    [1783] "g" "t" "t" "c" "c" "a" "a" "a" "g" "a" "t" "a" "a" "t" "a" "g" "a" "a"
    [1801] "a" "t" "g" "a" "c" "a" "c" "a" "g" "a" "a" "g" "g" "c" "t" "t" "t" "a"
    [1819] "a" "g" "t" "a" "t" "c" "c" "a" "t" "t" "g" "g" "g" "a" "c" "a" "t" "g"
    [1837] "a" "a" "g" "t" "t" "a" "a" "c" "c" "a" "c" "a" "g" "t" "c" "g" "g" "g"
    [1855] "a" "a" "a" "c" "a" "a" "g" "c" "a" "t" "a" "g" "a" "a" "a" "t" "g" "g"
    [1873] "a" "a" "g" "a" "a" "a" "g" "t" "g" "a" "a" "c" "t" "t" "g" "a" "t" "g"
    [1891] "c" "t" "c" "a" "g" "t" "a" "t" "t" "t" "g" "c" "a" "g" "a" "a" "t" "a"
    [1909] "c" "a" "t" "t" "c" "a" "a" "g" "g" "t" "t" "t" "c" "a" "a" "a" "g" "c"
    [1927] "g" "c" "c" "a" "g" "t" "c" "a" "t" "t" "t" "g" "c" "t" "c" "t" "g" "t"
    [1945] "t" "t" "t" "c" "a" "a" "a" "t" "c" "c" "a" "g" "g" "a" "a" "a" "t" "g"
    [1963] "c" "a" "g" "a" "a" "g" "a" "g" "g" "a" "a" "t" "g" "t" "g" "c" "a" "a"
    [1981] "c" "a" "t" "t" "c" "t" "c" "t" "g" "c" "c" "c" "a" "c" "t" "c" "t" "g"
    [1999] "g" "g" "t" "c" "c" "t" "t" "a" "a" "a" "g" "a" "a" "a" "c" "a" "a" "a"
    [2017] "g" "t" "c" "c" "a" "a" "a" "a" "g" "t" "c" "a" "c" "t" "t" "t" "t" "g"
    [2035] "a" "a" "t" "g" "t" "g" "a" "a" "c" "a" "a" "a" "a" "g" "g" "a" "a" "g"
    [2053] "a" "a" "a" "a" "t" "c" "a" "a" "g" "g" "a" "a" "a" "g" "a" "a" "t" "g"
    [2071] "a" "g" "t" "c" "t" "a" "a" "t" "a" "t" "c" "a" "a" "g" "c" "c" "t" "g"
    [2089] "t" "a" "c" "a" "g" "a" "c" "a" "g" "t" "t" "a" "a" "t" "a" "t" "c" "a"
    [2107] "c" "t" "g" "c" "a" "g" "g" "c" "t" "t" "t" "c" "c" "t" "g" "t" "g" "g"
    [2125] "t" "t" "g" "g" "t" "c" "a" "g" "a" "a" "a" "g" "a" "t" "a" "a" "g" "c"
    [2143] "c" "a" "g" "t" "t" "g" "a" "t" "a" "a" "t" "g" "c" "c" "a" "a" "a" "t"
    [2161] "g" "t" "a" "g" "t" "a" "t" "c" "a" "a" "a" "g" "g" "a" "g" "g" "c" "t"
    [2179] "c" "t" "a" "g" "g" "t" "t" "t" "t" "g" "t" "c" "t" "a" "t" "c" "a" "t"
    [2197] "c" "t" "c" "a" "g" "t" "t" "c" "a" "g" "a" "g" "g" "c" "a" "a" "c" "g"
    [2215] "a" "a" "a" "c" "t" "g" "g" "a" "c" "t" "c" "a" "t" "t" "a" "c" "t" "c"
    [2233] "c" "a" "a" "a" "t" "a" "a" "a" "c" "a" "t" "g" "g" "a" "c" "t" "t" "t"
    [2251] "t" "a" "c" "a" "a" "a" "a" "c" "c" "c" "a" "t" "a" "t" "c" "g" "t" "a"
    [2269] "t" "a" "c" "c" "a" "c" "c" "a" "c" "t" "t" "t" "t" "t" "c" "c" "c" "a"
    [2287] "t" "c" "a" "a" "g" "t" "c" "a" "t" "t" "t" "g" "t" "t" "a" "a" "a" "a"
    [2305] "c" "t" "a" "a" "a" "t" "g" "t" "a" "a" "g" "a" "a" "a" "a" "a" "t" "c"
    [2323] "t" "g" "c" "t" "a" "g" "a" "g" "g" "a" "a" "a" "a" "c" "t" "t" "t" "g"
    [2341] "a" "g" "g" "a" "a" "c" "a" "t" "t" "c" "a" "a" "t" "g" "t" "c" "a" "c"
    [2359] "c" "t" "g" "a" "a" "a" "g" "a" "g" "a" "a" "a" "t" "g" "g" "g" "a" "a"
    [2377] "a" "t" "g" "a" "g" "a" "a" "c" "a" "t" "t" "c" "c" "a" "a" "g" "t" "a"
    [2395] "c" "a" "g" "t" "g" "a" "g" "c" "a" "c" "a" "a" "t" "t" "a" "g" "c" "c"
    [2413] "g" "t" "a" "a" "t" "a" "a" "c" "a" "t" "t" "a" "g" "a" "g" "a" "a" "a"
    [2431] "a" "t" "g" "t" "t" "t" "t" "t" "a" "a" "a" "g" "g" "a" "g" "c" "c" "a"
    [2449] "g" "c" "t" "c" "a" "a" "g" "c" "a" "a" "t" "a" "t" "t" "a" "a" "t" "g"
    [2467] "a" "a" "g" "t" "a" "g" "g" "t" "t" "c" "c" "a" "g" "t" "a" "c" "t" "a"
    [2485] "a" "t" "g" "a" "a" "g" "t" "g" "g" "g" "c" "t" "c" "c" "a" "g" "t" "a"
    [2503] "t" "t" "a" "a" "t" "g" "a" "a" "a" "t" "a" "g" "g" "t" "t" "c" "c" "a"
    [2521] "g" "t" "g" "a" "t" "g" "a" "a" "a" "a" "c" "a" "t" "t" "c" "a" "a" "g"
    [2539] "c" "a" "g" "a" "a" "c" "t" "a" "g" "g" "t" "a" "g" "a" "a" "a" "c" "a"
    [2557] "g" "a" "g" "g" "g" "c" "c" "a" "a" "a" "a" "t" "t" "g" "a" "a" "t" "g"
    [2575] "c" "t" "a" "t" "g" "c" "t" "t" "a" "g" "a" "t" "t" "a" "g" "g" "g" "g"
    [2593] "t" "t" "t" "t" "g" "c" "a" "a" "c" "c" "t" "g" "a" "g" "g" "t" "c" "t"
    [2611] "a" "t" "a" "a" "a" "c" "a" "a" "a" "g" "t" "c" "t" "t" "c" "c" "t" "g"
    [2629] "g" "a" "a" "g" "t" "a" "a" "t" "t" "g" "t" "a" "a" "g" "c" "a" "t" "c"
    [2647] "c" "t" "g" "a" "a" "a" "t" "a" "a" "a" "a" "a" "a" "g" "c" "a" "a" "g"
    [2665] "a" "a" "t" "a" "t" "g" "a" "a" "g" "a" "a" "g" "t" "a" "g" "t" "t" "c"
    [2683] "a" "g" "a" "c" "t" "g" "t" "t" "a" "a" "t" "a" "c" "a" "g" "a" "t" "t"
    [2701] "t" "c" "t" "c" "t" "c" "c" "a" "t" "a" "t" "c" "t" "g" "a" "t" "t" "t"
    [2719] "c" "a" "g" "a" "t" "a" "a" "c" "t" "t" "a" "g" "a" "a" "c" "a" "g" "c"
    [2737] "c" "t" "a" "t" "g" "g" "g" "a" "a" "g" "t" "a" "g" "t" "c" "a" "t" "g"
    [2755] "c" "a" "t" "c" "t" "c" "a" "g" "g" "t" "t" "t" "g" "t" "t" "c" "t" "g"
    [2773] "a" "g" "a" "c" "a" "c" "c" "t" "g" "a" "t" "g" "a" "c" "c" "t" "g" "t"
    [2791] "t" "a" "g" "a" "t" "g" "a" "t" "g" "g" "t" "g" "a" "a" "a" "t" "a" "a"
    [2809] "a" "g" "g" "a" "a" "g" "a" "t" "a" "c" "t" "a" "g" "t" "t" "t" "t" "g"
    [2827] "c" "t" "g" "a" "a" "a" "a" "t" "g" "a" "c" "a" "t" "t" "a" "a" "g" "g"
    [2845] "a" "a" "a" "g" "t" "t" "c" "t" "g" "c" "t" "g" "t" "t" "t" "t" "t" "a"
    [2863] "g" "c" "a" "a" "a" "a" "g" "c" "g" "t" "c" "c" "a" "g" "a" "g" "a" "g"
    [2881] "g" "a" "g" "a" "g" "c" "t" "t" "a" "g" "c" "a" "g" "g" "a" "g" "t" "c"
    [2899] "c" "t" "a" "g" "c" "c" "c" "t" "t" "t" "c" "a" "c" "c" "c" "a" "t" "a"
    [2917] "c" "a" "c" "a" "t" "t" "t" "g" "g" "c" "t" "c" "a" "g" "g" "g" "t" "t"
    [2935] "a" "c" "c" "g" "a" "a" "g" "a" "g" "g" "g" "g" "c" "c" "a" "a" "g" "a"
    [2953] "a" "a" "t" "t" "a" "g" "a" "g" "t" "c" "c" "t" "c" "a" "g" "a" "a" "g"
    [2971] "a" "g" "a" "a" "c" "t" "t" "a" "t" "c" "t" "a" "g" "t" "g" "a" "g" "g"
    [2989] "a" "t" "g" "a" "a" "g" "a" "g" "c" "t" "t" "c" "c" "c" "t" "g" "c" "t"
    [3007] "t" "c" "c" "a" "a" "c" "a" "c" "t" "t" "g" "t" "t" "a" "t" "t" "t" "g"
    [3025] "g" "t" "a" "a" "a" "g" "t" "a" "a" "a" "c" "a" "a" "t" "a" "t" "a" "c"
    [3043] "c" "t" "t" "c" "t" "c" "a" "g" "t" "c" "t" "a" "c" "t" "a" "g" "g" "c"
    [3061] "a" "t" "a" "g" "c" "a" "c" "c" "g" "t" "t" "g" "c" "t" "a" "c" "c" "g"
    [3079] "a" "g" "t" "g" "t" "c" "t" "g" "t" "c" "t" "a" "a" "g" "a" "a" "c" "a"
    [3097] "c" "a" "g" "a" "g" "g" "a" "g" "a" "a" "t" "t" "t" "a" "t" "t" "a" "t"
    [3115] "c" "a" "t" "t" "g" "a" "a" "g" "a" "a" "t" "a" "g" "c" "t" "t" "a" "a"
    [3133] "a" "t" "g" "a" "c" "t" "g" "c" "a" "g" "t" "a" "a" "c" "c" "a" "g" "g"
    [3151] "t" "a" "a" "t" "a" "t" "t" "g" "g" "c" "a" "a" "a" "g" "g" "c" "a" "t"
    [3169] "c" "t" "c" "a" "g" "g" "a" "a" "c" "a" "t" "c" "a" "c" "c" "t" "t" "a"
    [3187] "g" "t" "g" "a" "g" "g" "a" "a" "a" "c" "a" "a" "a" "a" "t" "g" "t" "t"
    [3205] "c" "t" "g" "c" "t" "a" "g" "c" "t" "t" "g" "t" "t" "t" "t" "c" "t" "t"
    [3223] "c" "a" "c" "a" "g" "t" "g" "c" "a" "g" "t" "g" "a" "a" "t" "t" "g" "g"
    [3241] "a" "a" "g" "a" "c" "t" "t" "g" "a" "c" "t" "g" "c" "a" "a" "a" "t" "a"
    [3259] "c" "a" "a" "a" "c" "a" "c" "c" "c" "a" "g" "g" "a" "t" "c" "c" "t" "t"
    [3277] "t" "c" "t" "t" "g" "a" "t" "t" "g" "g" "t" "t" "c" "t" "t" "c" "c" "a"
    [3295] "a" "a" "c" "a" "a" "a" "t" "g" "a" "g" "g" "c" "a" "t" "c" "a" "g" "t"
    [3313] "c" "t" "g" "a" "a" "a" "g" "c" "c" "a" "g" "g" "g" "a" "g" "t" "t" "g"
    [3331] "g" "t" "c" "t" "g" "a" "g" "t" "g" "a" "c" "a" "a" "g" "g" "a" "a" "t"
    [3349] "t" "g" "g" "t" "t" "t" "c" "a" "g" "a" "t" "g" "a" "t" "g" "a" "a" "g"
    [3367] "a" "a" "a" "g" "a" "g" "g" "a" "a" "c" "g" "g" "g" "c" "t" "t" "g" "g"
    [3385] "a" "a" "g" "a" "a" "a" "a" "t" "a" "a" "t" "c" "a" "a" "g" "a" "a" "g"
    [3403] "a" "g" "c" "a" "a" "a" "g" "c" "a" "t" "g" "g" "a" "t" "t" "c" "a" "a"
    [3421] "a" "c" "t" "t" "a" "g"
    attr(,"name")
    [1] "AY304547.BRCA1"
    attr(,"Annot")
    [1] ">AY304547.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $AY438030
      [1] "a" "t" "g" "g" "c" "g" "g" "t" "g" "c" "a" "g" "g" "t" "g" "g" "t" "g"
     [19] "c" "a" "g" "g" "c" "g" "g" "t" "g" "c" "a" "g" "g" "c" "g" "g" "t" "t"
     [37] "c" "a" "t" "c" "t" "c" "g" "a" "g" "t" "c" "t" "g" "a" "c" "g" "c" "t"
     [55] "t" "t" "c" "c" "t" "c" "g" "t" "t" "t" "g" "t" "c" "t" "c" "a" "a" "c"
     [73] "c" "a" "c" "g" "c" "t" "c" "t" "g" "a" "g" "c" "a" "c" "a" "g" "a" "g"
     [91] "a" "a" "g" "g" "a" "g" "g" "a" "a" "g" "t" "a" "a" "t" "g" "g" "g" "g"
    [109] "c" "t" "g" "t" "g" "c" "a" "t" "a" "g" "g" "g" "g" "a" "g" "t" "t" "g"
    [127] "a" "a" "c" "g" "a" "t" "g" "a" "t" "a" "c" "a" "a" "g" "g" "a" "g" "t"
    [145] "g" "a" "c" "t" "c" "c" "a" "a" "a" "t" "t" "t" "g" "c" "a" "t" "a" "t"
    [163] "a" "c" "t" "g" "g" "a" "a" "c" "t" "g" "a" "a" "a" "t" "g" "c" "g" "c"
    [181] "a" "c" "a" "g" "t" "t" "g" "c" "t" "g" "a" "a" "a" "a" "g" "g" "t" "t"
    [199] "g" "a" "t" "g" "c" "c" "g" "t" "c" "a" "g" "a" "a" "t" "t" "g" "t" "t"
    [217] "c" "a" "c" "a" "t" "t" "c" "a" "t" "t" "c" "t" "g" "t" "c" "a" "t" "c"
    [235] "a" "t" "c" "t" "t" "a" "c" "g" "a" "c" "g" "t" "t" "c" "t" "g" "a" "t"
    [253] "a" "a" "g" "a" "g" "g" "a" "a" "g" "g" "a" "c" "c" "g" "a" "g" "t" "a"
    [271] "g" "a" "a" "a" "t" "t" "t" "c" "t" "c" "c" "a" "g" "a" "g" "c" "a" "g"
    [289] "c" "t" "g" "t" "c" "t" "g" "c" "a" "g" "c" "t" "t" "c" "a" "a" "c" "a"
    [307] "g" "a" "g" "g" "c" "a" "g" "a" "g" "a" "g" "g" "t" "t" "g" "g" "c" "t"
    [325] "g" "a" "a" "c" "t" "g" "a" "c" "a" "g" "g" "c" "c" "g" "c" "c" "c" "c"
    [343] "a" "t" "g" "a" "g" "a" "g" "t" "t" "g" "t" "g" "g" "g" "c" "t" "g" "g"
    [361] "t" "a" "t" "c" "a" "t" "t" "c" "c" "c" "a" "t" "c" "c" "t" "c" "a" "t"
    [379] "a" "t" "a" "a" "c" "t" "g" "t" "t" "t" "g" "g" "c" "c" "t" "t" "c" "a"
    [397] "c" "a" "t" "g" "t" "t" "g" "a" "t" "g" "t" "t" "c" "g" "c" "a" "c" "a"
    [415] "c" "a" "a" "g" "c" "c" "a" "t" "g" "t" "a" "c" "c" "a" "g" "a" "t" "g"
    [433] "a" "t" "g" "g" "a" "t" "c" "a" "a" "g" "g" "c" "t" "t" "t" "g" "t" "a"
    [451] "g" "g" "a" "c" "t" "t" "a" "t" "t" "t" "t" "t" "t" "c" "c" "t" "g" "t"
    [469] "t" "t" "c" "a" "t" "a" "g" "a" "a" "g" "a" "t" "a" "a" "g" "a" "a" "c"
    [487] "a" "c" "a" "a" "a" "g" "a" "c" "t" "g" "g" "c" "c" "g" "g" "g" "t" "a"
    [505] "c" "t" "c" "t" "a" "c" "a" "c" "t" "t" "g" "c" "t" "t" "c" "c" "a" "a"
    [523] "t" "c" "c" "a" "t" "a" "c" "a" "g" "g" "c" "c" "c" "a" "a" "a" "a" "g"
    [541] "a" "g" "t" "t" "c" "a" "g" "a" "g" "t" "c" "c" "c" "t" "t" "c" "a" "t"
    [559] "g" "g" "t" "c" "c" "a" "c" "g" "a" "g" "a" "c" "t" "t" "c" "t" "g" "g"
    [577] "a" "g" "c" "t" "c" "c" "a" "g" "c" "c" "a" "g" "c" "a" "c" "a" "t" "c"
    [595] "t" "c" "c" "a" "t" "t" "g" "a" "g" "g" "g" "c" "c" "a" "g" "a" "a" "g"
    [613] "g" "a" "a" "g" "a" "g" "g" "a" "a" "a" "g" "g" "t" "a" "t" "g" "a" "g"
    [631] "a" "g" "a" "a" "t" "c" "g" "a" "a" "a" "t" "c" "c" "c" "a" "a" "t" "c"
    [649] "c" "a" "t" "a" "t" "t" "g" "t" "a" "c" "c" "t" "c" "a" "t" "g" "t" "c"
    [667] "a" "c" "t" "a" "t" "c" "g" "g" "g" "a" "a" "a" "g" "t" "g" "t" "g" "c"
    [685] "c" "t" "t" "g" "a" "a" "t" "c" "a" "g" "c" "a" "g" "t" "a" "g" "a" "g"
    [703] "c" "t" "g" "c" "c" "c" "a" "a" "g" "a" "t" "c" "c" "t" "g" "t" "g" "c"
    [721] "c" "a" "g" "g" "a" "g" "g" "a" "g" "c" "a" "g" "g" "a" "t" "g" "c" "g"
    [739] "t" "a" "t" "a" "g" "g" "a" "g" "g" "a" "t" "c" "c" "a" "c" "a" "g" "c"
    [757] "c" "t" "t" "a" "c" "a" "c" "a" "t" "c" "t" "g" "g" "a" "c" "t" "c" "a"
    [775] "g" "t" "a" "a" "c" "c" "a" "a" "g" "a" "t" "c" "c" "a" "t" "a" "a" "t"
    [793] "g" "g" "c" "t" "c" "a" "g" "t" "g" "t" "t" "t" "a" "c" "c" "a" "a" "g"
    [811] "a" "a" "t" "c" "t" "g" "t" "g" "c" "a" "g" "t" "c" "a" "g" "a" "t" "g"
    [829] "t" "c" "g" "g" "c" "a" "g" "t" "c" "a" "g" "c" "g" "g" "g" "c" "c" "t"
    [847] "c" "t" "c" "c" "t" "a" "c" "a" "g" "t" "g" "g" "t" "t" "g" "g" "a" "g"
    [865] "g" "a" "c" "a" "g" "a" "c" "t" "g" "g" "a" "g" "c" "a" "a" "a" "a" "c"
    [883] "c" "a" "a" "c" "a" "g" "c" "a" "t" "t" "t" "g" "c" "a" "g" "g" "a" "a"
    [901] "t" "t" "a" "c" "a" "a" "c" "a" "a" "g" "a" "a" "a" "a" "g" "g" "a" "a"
    [919] "g" "a" "g" "c" "t" "t" "a" "t" "g" "c" "a" "a" "g" "a" "a" "c" "t" "t"
    [937] "t" "c" "t" "t" "c" "t" "c" "t" "a" "g" "a" "a" "t" "a" "a"
    attr(,"name")
    [1] "AY438030"
    attr(,"Annot")
    [1] ">AY438030"
    attr(,"class")
    [1] "SeqFastadna"
    
    $AY438031
       [1] "a" "t" "g" "t" "c" "c" "c" "c" "a" "g" "a" "a" "g" "t" "g" "g" "c" "c"
      [19] "t" "t" "g" "a" "a" "c" "c" "g" "a" "a" "t" "a" "t" "c" "t" "c" "c" "a"
      [37] "a" "t" "g" "c" "t" "c" "t" "c" "c" "c" "c" "t" "t" "t" "c" "a" "t" "a"
      [55] "t" "c" "t" "a" "g" "c" "g" "t" "g" "g" "t" "c" "c" "g" "g" "a" "a" "t"
      [73] "g" "g" "a" "a" "a" "a" "g" "t" "g" "g" "g" "a" "c" "t" "g" "g" "a" "t"
      [91] "g" "c" "t" "a" "c" "a" "a" "a" "c" "t" "g" "t" "t" "t" "g" "a" "g" "g"
     [109] "a" "t" "a" "a" "c" "t" "g" "a" "c" "t" "t" "a" "a" "a" "a" "t" "c" "t"
     [127] "g" "g" "c" "t" "g" "c" "a" "c" "a" "t" "c" "a" "t" "t" "g" "a" "c" "t"
     [145] "c" "c" "t" "g" "g" "g" "c" "c" "c" "a" "a" "c" "t" "g" "t" "g" "a" "c"
     [163] "c" "g" "a" "t" "t" "t" "a" "a" "a" "c" "t" "g" "c" "a" "c" "a" "t" "a"
     [181] "c" "c" "a" "t" "a" "t" "g" "c" "t" "g" "g" "a" "g" "a" "g" "a" "c" "a"
     [199] "t" "t" "a" "a" "a" "g" "t" "g" "g" "g" "a" "t" "a" "t" "c" "a" "t" "t"
     [217] "t" "t" "c" "a" "a" "t" "g" "c" "c" "c" "a" "a" "t" "a" "c" "c" "c" "a"
     [235] "g" "a" "a" "c" "t" "g" "c" "c" "t" "c" "c" "c" "g" "a" "t" "t" "t" "t"
     [253] "a" "t" "c" "t" "t" "t" "g" "g" "a" "g" "a" "a" "g" "a" "t" "g" "c" "t"
     [271] "g" "a" "a" "t" "t" "c" "c" "t" "g" "c" "c" "a" "g" "a" "c" "c" "c" "c"
     [289] "t" "c" "a" "g" "c" "t" "t" "t" "g" "c" "a" "g" "a" "a" "t" "c" "t" "t"
     [307] "g" "c" "c" "t" "c" "c" "t" "g" "g" "a" "a" "t" "c" "c" "t" "t" "c" "a"
     [325] "a" "a" "t" "c" "c" "t" "g" "a" "a" "t" "g" "t" "c" "t" "c" "t" "t" "a"
     [343] "c" "t" "t" "g" "t" "g" "g" "t" "g" "a" "a" "g" "g" "a" "a" "c" "t" "t"
     [361] "g" "t" "g" "c" "a" "a" "c" "a" "a" "t" "a" "t" "c" "a" "c" "c" "a" "a"
     [379] "t" "t" "c" "c" "a" "a" "t" "g" "t" "a" "g" "c" "c" "g" "c" "c" "t" "c"
     [397] "c" "g" "g" "g" "a" "g" "a" "g" "c" "t" "c" "c" "c" "g" "c" "c" "t" "c"
     [415] "a" "t" "g" "t" "t" "t" "g" "a" "a" "t" "a" "c" "c" "a" "g" "a" "c" "a"
     [433] "t" "t" "a" "c" "t" "g" "g" "a" "g" "g" "a" "g" "c" "c" "a" "c" "a" "g"
     [451] "t" "a" "t" "g" "g" "a" "g" "a" "g" "a" "a" "c" "a" "t" "g" "g" "a" "a"
     [469] "a" "t" "t" "t" "a" "t" "g" "c" "t" "g" "g" "g" "a" "a" "a" "a" "a" "a"
     [487] "a" "a" "c" "a" "a" "c" "t" "g" "g" "a" "c" "t" "g" "g" "t" "g" "a" "a"
     [505] "t" "t" "t" "t" "c" "a" "g" "c" "t" "c" "g" "t" "t" "t" "c" "c" "t" "t"
     [523] "t" "t" "g" "a" "a" "g" "c" "t" "g" "c" "c" "c" "g" "t" "a" "g" "a" "t"
     [541] "t" "t" "c" "a" "g" "c" "a" "a" "t" "a" "t" "c" "c" "c" "c" "a" "c" "a"
     [559] "t" "a" "c" "c" "t" "t" "c" "t" "c" "a" "a" "g" "g" "a" "t" "g" "t" "a"
     [577] "a" "a" "t" "g" "a" "a" "g" "a" "c" "c" "c" "t" "g" "g" "a" "g" "a" "a"
     [595] "g" "a" "t" "g" "t" "g" "g" "c" "c" "c" "t" "c" "c" "t" "c" "t" "c" "t"
     [613] "g" "t" "t" "a" "g" "t" "t" "t" "t" "g" "a" "g" "g" "a" "c" "a" "c" "t"
     [631] "g" "a" "a" "g" "c" "c" "a" "c" "c" "c" "a" "g" "g" "t" "g" "t" "a" "c"
     [649] "c" "c" "c" "a" "a" "g" "c" "t" "g" "t" "a" "c" "t" "t" "g" "t" "c" "a"
     [667] "c" "c" "t" "c" "g" "a" "a" "t" "t" "g" "a" "g" "c" "a" "t" "g" "c" "a"
     [685] "c" "t" "t" "g" "g" "a" "g" "g" "c" "t" "c" "c" "t" "c" "a" "g" "c" "t"
     [703] "c" "t" "t" "c" "a" "t" "a" "t" "c" "c" "c" "a" "g" "c" "t" "t" "t" "t"
     [721] "c" "c" "a" "g" "g" "a" "g" "g" "a" "g" "g" "a" "t" "g" "t" "c" "t" "c"
     [739] "a" "t" "t" "g" "a" "t" "t" "a" "c" "g" "t" "t" "c" "c" "t" "c" "a" "a"
     [757] "g" "t" "a" "t" "g" "c" "c" "a" "c" "c" "t" "g" "c" "t" "c" "a" "c" "c"
     [775] "a" "a" "c" "a" "a" "g" "g" "t" "g" "c" "a" "g" "t" "a" "c" "g" "t" "g"
     [793] "a" "t" "t" "c" "a" "a" "g" "g" "g" "t" "a" "t" "c" "a" "c" "a" "a" "a"
     [811] "a" "g" "a" "a" "g" "a" "g" "a" "g" "t" "a" "t" "a" "t" "t" "g" "c" "t"
     [829] "g" "c" "t" "t" "t" "t" "c" "t" "c" "a" "g" "t" "c" "a" "c" "t" "t" "t"
     [847] "g" "g" "c" "a" "c" "a" "g" "g" "t" "g" "t" "c" "g" "t" "g" "g" "a" "a"
     [865] "t" "a" "t" "g" "a" "t" "g" "c" "a" "g" "a" "a" "g" "g" "c" "t" "t" "t"
     [883] "a" "c" "a" "a" "a" "a" "c" "t" "c" "a" "c" "t" "c" "t" "g" "c" "t" "g"
     [901] "c" "t" "g" "a" "t" "g" "t" "g" "g" "a" "a" "a" "g" "a" "t" "t" "t" "t"
     [919] "t" "g" "t" "t" "t" "t" "c" "t" "t" "g" "t" "a" "c" "a" "c" "a" "t" "t"
     [937] "g" "a" "c" "c" "t" "g" "c" "c" "t" "c" "t" "g" "t" "t" "t" "t" "t" "c"
     [955] "c" "c" "t" "c" "g" "a" "g" "a" "c" "c" "a" "g" "c" "c" "a" "a" "c" "t"
     [973] "c" "t" "c" "a" "c" "a" "t" "t" "t" "c" "a" "g" "t" "c" "c" "g" "t" "t"
     [991] "t" "a" "t" "c" "a" "c" "t" "t" "t" "a" "c" "c" "a" "a" "c" "a" "g" "t"
    [1009] "g" "g" "a" "c" "a" "g" "c" "t" "t" "t" "a" "c" "t" "c" "c" "c" "a" "g"
    [1027] "g" "c" "c" "c" "a" "a" "a" "a" "a" "a" "a" "t" "t" "a" "t" "c" "c" "g"
    [1045] "t" "a" "c" "a" "g" "c" "c" "c" "c" "a" "g" "a" "t" "g" "g" "g" "a" "t"
    [1063] "g" "g" "a" "a" "a" "t" "g" "a" "a" "a" "t" "g" "g" "c" "c" "a" "a" "a"
    [1081] "a" "g" "a" "g" "c" "a" "a" "a" "g" "g" "c" "t" "t" "a" "t" "t" "t" "c"
    [1099] "a" "a" "a" "a" "c" "c" "t" "t" "t" "g" "t" "c" "c" "c" "t" "c" "a" "g"
    [1117] "t" "t" "c" "c" "a" "g" "g" "a" "g" "g" "c" "a" "g" "c" "a" "t" "t" "t"
    [1135] "g" "c" "c" "a" "a" "t" "g" "g" "a" "a" "a" "g" "c" "t" "c" "t" "a" "g"
    attr(,"name")
    [1] "AY438031"
    attr(,"Annot")
    [1] ">AY438031"
    attr(,"class")
    [1] "SeqFastadna"
    
    $AY706911
      [1] "g" "a" "g" "c" "c" "a" "c" "a" "g" "a" "t" "a" "a" "t" "a" "c" "a" "a"
     [19] "g" "a" "g" "c" "g" "t" "c" "c" "c" "c" "t" "c" "a" "c" "a" "a" "a" "t"
     [37] "a" "a" "a" "t" "t" "a" "a" "a" "g" "c" "a" "t" "a" "a" "a" "a" "g" "g"
     [55] "a" "g" "a" "c" "c" "t" "a" "c" "a" "t" "c" "a" "g" "g" "c" "c" "t" "t"
     [73] "c" "a" "t" "c" "c" "t" "g" "a" "g" "g" "a" "t" "t" "t" "t" "a" "t" "c"
     [91] "a" "a" "g" "a" "a" "a" "g" "c" "a" "g" "a" "t" "t" "t" "g" "g" "c" "a"
    [109] "g" "t" "t" "c" "a" "a" "a" "a" "g" "a" "c" "t" "c" "c" "t" "g" "a" "a"
    [127] "a" "t" "g" "a" "t" "a" "a" "a" "t" "c" "a" "g" "g" "g" "a" "a" "c" "t"
    [145] "a" "a" "c" "c" "a" "a" "a" "c" "g" "g" "a" "g" "c" "a" "g" "a" "a" "t"
    [163] "g" "g" "t" "c" "a" "a" "g" "t" "g" "a" "t" "g" "a" "a" "t" "a" "t" "t"
    [181] "a" "c" "t" "a" "a" "t" "a" "g" "t" "g" "g" "t" "c" "a" "t" "g" "a" "g"
    [199] "a" "a" "t" "a" "a" "a" "a" "c" "a" "a" "a" "a" "g" "g" "t" "g" "a" "t"
    [217] "t" "c" "t" "a" "t" "t" "c" "a" "g" "a" "a" "t" "g" "a" "g" "a" "a" "a"
    [235] "a" "a" "t" "c" "c" "t" "a" "a" "c" "c" "c" "a" "a" "t" "a" "g" "a" "a"
    [253] "t" "c" "a" "c" "t" "c" "g" "a" "a" "a" "a" "a" "g" "a" "a" "t" "c" "t"
    [271] "g" "c"
    attr(,"name")
    [1] "AY706911"
    attr(,"Annot")
    [1] ">AY706911"
    attr(,"class")
    [1] "SeqFastadna"
    
    $AY706912.BRCA1
     [1] "c" "a" "g" "t" "a" "t" "t" "a" "a" "c" "t" "t" "c" "a" "c" "a" "g" "a" "a"
    [20] "a" "a" "g" "t" "a" "g" "t" "g" "a" "a" "t" "a" "c" "c" "c" "t" "a" "t" "a"
    [39] "a" "g" "c" "t" "a" "g"
    attr(,"name")
    [1] "AY706912.BRCA1"
    attr(,"Annot")
    [1] ">AY706912.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $AY706913
     [1] "a" "g" "a" "t" "g" "c" "t" "c" "g" "t" "g" "t" "a" "c" "a" "a" "g" "t" "t"
    [20] "t" "g" "c" "c" "a" "g" "a" "a" "a" "a" "c" "a" "c" "c" "a" "c" "a" "t" "c"
    [39] "a" "c" "t" "t" "t" "t" "a" "a" "c" "t" "a" "a"
    attr(,"name")
    [1] "AY706913"
    attr(,"Annot")
    [1] ">AY706913"
    attr(,"class")
    [1] "SeqFastadna"
    
    $AY751490
       [1] "a" "t" "g" "g" "a" "t" "t" "t" "a" "t" "c" "t" "g" "c" "t" "c" "t" "t"
      [19] "c" "g" "c" "g" "t" "t" "g" "a" "a" "g" "a" "a" "g" "t" "a" "c" "a" "a"
      [37] "a" "a" "t" "g" "t" "c" "a" "t" "t" "a" "a" "t" "g" "c" "t" "a" "t" "g"
      [55] "c" "a" "g" "a" "a" "a" "a" "t" "c" "t" "t" "a" "g" "a" "g" "t" "g" "t"
      [73] "c" "c" "c" "a" "t" "c" "t" "g" "t" "c" "t" "g" "g" "a" "g" "t" "t" "g"
      [91] "a" "t" "c" "a" "a" "g" "g" "a" "a" "c" "c" "t" "g" "t" "c" "t" "c" "c"
     [109] "a" "c" "a" "a" "a" "g" "t" "g" "t" "g" "a" "c" "c" "a" "c" "a" "t" "a"
     [127] "t" "t" "t" "t" "g" "c" "a" "a" "a" "t" "t" "t" "t" "g" "c" "a" "t" "g"
     [145] "c" "t" "g" "a" "a" "a" "c" "t" "t" "c" "t" "c" "a" "a" "c" "c" "a" "g"
     [163] "a" "a" "g" "a" "a" "a" "g" "g" "g" "c" "c" "t" "t" "c" "a" "c" "a" "g"
     [181] "t" "g" "t" "c" "c" "t" "t" "t" "a" "t" "g" "t" "a" "a" "g" "a" "a" "t"
     [199] "g" "a" "t" "a" "t" "a" "a" "c" "c" "a" "a" "a" "a" "g" "g" "a" "g" "c"
     [217] "c" "t" "a" "c" "a" "a" "g" "a" "a" "a" "g" "t" "a" "c" "g" "a" "g" "a"
     [235] "t" "t" "t" "a" "g" "t" "c" "a" "a" "c" "t" "t" "g" "t" "t" "g" "a" "a"
     [253] "g" "a" "g" "c" "t" "a" "t" "t" "g" "a" "a" "a" "a" "t" "c" "a" "t" "t"
     [271] "t" "g" "t" "g" "c" "t" "t" "t" "t" "c" "a" "g" "c" "t" "t" "g" "a" "c"
     [289] "a" "c" "a" "g" "g" "t" "t" "t" "g" "g" "a" "g" "t" "a" "t" "g" "c" "a"
     [307] "a" "a" "c" "a" "g" "c" "t" "a" "t" "a" "a" "t" "t" "t" "t" "g" "c" "a"
     [325] "a" "a" "a" "a" "a" "g" "g" "a" "a" "a" "a" "t" "a" "a" "c" "t" "c" "t"
     [343] "c" "c" "t" "g" "a" "a" "c" "a" "t" "c" "t" "a" "a" "a" "a" "g" "a" "t"
     [361] "g" "a" "a" "g" "t" "t" "t" "c" "t" "a" "t" "c" "a" "t" "c" "c" "a" "a"
     [379] "a" "g" "t" "a" "t" "g" "g" "g" "c" "t" "a" "c" "a" "g" "a" "a" "a" "c"
     [397] "c" "g" "t" "g" "c" "c" "a" "a" "a" "a" "g" "a" "c" "t" "t" "c" "t" "a"
     [415] "c" "a" "g" "a" "g" "t" "g" "a" "a" "c" "c" "c" "g" "a" "a" "a" "a" "t"
     [433] "c" "c" "t" "t" "c" "c" "t" "t" "g" "c" "a" "g" "g" "a" "a" "a" "c" "c"
     [451] "a" "g" "t" "c" "t" "c" "a" "g" "t" "g" "t" "c" "c" "a" "a" "c" "t" "c"
     [469] "t" "c" "t" "a" "a" "c" "c" "t" "t" "g" "g" "a" "a" "c" "t" "g" "t" "g"
     [487] "a" "g" "a" "a" "c" "t" "c" "t" "g" "a" "g" "g" "a" "c" "a" "a" "a" "g"
     [505] "c" "a" "g" "c" "g" "g" "a" "t" "a" "c" "a" "a" "c" "c" "t" "c" "a" "a"
     [523] "a" "a" "g" "a" "c" "g" "t" "c" "t" "g" "t" "c" "t" "a" "c" "a" "t" "t"
     [541] "g" "a" "a" "t" "t" "g" "g" "g" "a" "t" "c" "t" "g" "a" "t" "t" "c" "t"
     [559] "t" "c" "t" "g" "a" "a" "g" "a" "t" "a" "c" "c" "g" "t" "t" "a" "a" "t"
     [577] "a" "a" "g" "g" "c" "a" "a" "c" "t" "t" "a" "t" "t" "g" "c" "a" "g" "t"
     [595] "g" "t" "g" "g" "g" "a" "g" "a" "t" "c" "a" "a" "g" "a" "a" "t" "t" "g"
     [613] "t" "t" "a" "c" "a" "a" "a" "t" "c" "a" "c" "c" "c" "c" "t" "c" "a" "a"
     [631] "g" "g" "a" "a" "c" "c" "a" "g" "g" "g" "a" "t" "g" "a" "a" "a" "t" "c"
     [649] "a" "g" "t" "t" "t" "g" "g" "a" "t" "t" "c" "t" "g" "c" "a" "a" "a" "a"
     [667] "a" "a" "g" "g" "c" "t" "g" "c" "t" "t" "g" "t" "g" "a" "a" "t" "t" "t"
     [685] "t" "c" "t" "g" "a" "g" "a" "c" "g" "g" "a" "t" "g" "t" "a" "a" "c" "a"
     [703] "a" "a" "t" "a" "c" "t" "g" "a" "a" "c" "a" "t" "c" "a" "t" "c" "a" "a"
     [721] "c" "c" "c" "a" "g" "t" "a" "a" "t" "a" "a" "t" "g" "a" "t" "t" "t" "g"
     [739] "a" "a" "c" "a" "c" "c" "a" "c" "t" "g" "a" "g" "a" "a" "g" "c" "g" "t"
     [757] "g" "c" "a" "g" "c" "t" "g" "a" "g" "a" "g" "g" "c" "a" "t" "c" "c" "a"
     [775] "g" "a" "a" "a" "a" "g" "t" "a" "t" "c" "a" "g" "g" "g" "t" "a" "g" "t"
     [793] "t" "c" "t" "g" "t" "t" "t" "c" "a" "a" "a" "c" "t" "t" "g" "c" "a" "t"
     [811] "g" "t" "g" "g" "a" "g" "c" "c" "a" "t" "g" "t" "g" "g" "c" "a" "c" "a"
     [829] "a" "a" "t" "a" "c" "t" "c" "a" "t" "g" "c" "c" "a" "g" "c" "t" "c" "a"
     [847] "t" "t" "a" "c" "a" "g" "c" "a" "t" "g" "a" "g" "a" "a" "c" "a" "g" "c"
     [865] "a" "g" "t" "t" "t" "a" "t" "t" "a" "c" "t" "c" "a" "c" "t" "a" "a" "a"
     [883] "g" "a" "c" "a" "g" "a" "a" "t" "g" "a" "a" "t" "g" "t" "a" "g" "a" "a"
     [901] "a" "a" "g" "g" "c" "t" "g" "a" "a" "t" "t" "c" "t" "g" "t" "a" "a" "t"
     [919] "a" "a" "a" "a" "g" "c" "a" "a" "a" "c" "a" "g" "c" "c" "t" "g" "g" "c"
     [937] "t" "t" "a" "g" "c" "a" "a" "g" "g" "a" "g" "c" "c" "a" "a" "c" "a" "t"
     [955] "a" "a" "c" "a" "g" "a" "t" "g" "g" "g" "c" "t" "g" "g" "a" "a" "g" "t"
     [973] "a" "a" "g" "g" "a" "a" "a" "c" "a" "t" "g" "t" "a" "a" "t" "g" "a" "t"
     [991] "a" "g" "g" "c" "g" "g" "a" "c" "t" "c" "c" "c" "a" "g" "c" "a" "c" "a"
    [1009] "g" "a" "a" "a" "a" "a" "a" "a" "g" "g" "t" "a" "g" "a" "t" "c" "t" "g"
    [1027] "a" "a" "t" "g" "c" "t" "g" "a" "t" "c" "c" "c" "c" "t" "g" "t" "g" "t"
    [1045] "g" "a" "g" "a" "g" "a" "a" "a" "a" "g" "a" "a" "t" "g" "g" "a" "a" "t"
    [1063] "a" "a" "g" "c" "a" "g" "a" "a" "a" "c" "t" "g" "c" "c" "a" "t" "g" "c"
    [1081] "t" "c" "a" "g" "a" "g" "a" "a" "t" "c" "c" "t" "a" "g" "a" "g" "a" "t"
    [1099] "a" "c" "t" "g" "a" "a" "g" "a" "t" "g" "t" "t" "c" "c" "t" "t" "g" "g"
    [1117] "a" "t" "a" "a" "c" "a" "c" "t" "a" "a" "a" "t" "a" "g" "c" "a" "g" "c"
    [1135] "a" "t" "t" "c" "a" "g" "a" "a" "a" "g" "t" "t" "a" "a" "t" "g" "a" "g"
    [1153] "t" "g" "g" "t" "t" "t" "t" "c" "c" "a" "g" "a" "a" "g" "t" "g" "a" "t"
    [1171] "g" "a" "a" "c" "t" "g" "t" "t" "a" "g" "g" "t" "t" "c" "t" "g" "a" "t"
    [1189] "g" "a" "c" "t" "c" "a" "c" "a" "t" "g" "a" "t" "g" "g" "g" "g" "a" "g"
    [1207] "t" "c" "t" "g" "a" "a" "t" "c" "a" "a" "a" "t" "g" "c" "c" "a" "a" "a"
    [1225] "g" "t" "a" "g" "c" "t" "g" "a" "t" "g" "t" "a" "t" "t" "g" "g" "a" "c"
    [1243] "g" "t" "t" "c" "t" "a" "a" "a" "t" "g" "a" "g" "g" "t" "a" "g" "a" "t"
    [1261] "g" "a" "a" "t" "a" "t" "t" "c" "t" "g" "g" "t" "t" "c" "t" "t" "c" "a"
    [1279] "g" "a" "g" "a" "a" "a" "a" "t" "a" "g" "a" "c" "t" "t" "a" "c" "t" "g"
    [1297] "g" "c" "c" "a" "g" "t" "g" "a" "t" "c" "c" "t" "c" "a" "t" "g" "a" "g"
    [1315] "g" "c" "t" "t" "t" "a" "a" "t" "a" "t" "g" "t" "a" "a" "a" "a" "g" "t"
    [1333] "g" "a" "a" "a" "g" "a" "g" "t" "t" "c" "a" "c" "t" "c" "c" "a" "a" "a"
    [1351] "t" "c" "a" "g" "t" "a" "g" "a" "g" "a" "g" "t" "a" "a" "t" "a" "t" "t"
    [1369] "g" "a" "a" "g" "a" "c" "a" "a" "a" "a" "t" "a" "t" "t" "t" "g" "g" "g"
    [1387] "a" "a" "a" "a" "c" "c" "t" "a" "t" "c" "g" "g" "a" "a" "g" "a" "a" "g"
    [1405] "g" "c" "a" "a" "g" "c" "c" "t" "c" "c" "c" "c" "a" "a" "c" "t" "t" "a"
    [1423] "a" "g" "c" "c" "a" "t" "g" "t" "a" "a" "c" "t" "g" "a" "a" "a" "a" "t"
    [1441] "c" "t" "a" "a" "t" "t" "a" "t" "a" "g" "g" "a" "g" "c" "a" "t" "t" "t"
    [1459] "g" "t" "t" "a" "c" "t" "g" "a" "g" "c" "c" "a" "c" "a" "g" "a" "t" "a"
    [1477] "a" "t" "a" "c" "a" "a" "g" "a" "g" "c" "g" "t" "c" "c" "c" "c" "t" "c"
    [1495] "a" "c" "a" "a" "a" "t" "a" "a" "a" "t" "t" "a" "a" "a" "g" "c" "g" "t"
    [1513] "a" "a" "a" "a" "g" "g" "a" "g" "a" "c" "c" "t" "a" "c" "a" "t" "c" "a"
    [1531] "g" "g" "c" "c" "t" "t" "c" "a" "t" "c" "c" "t" "g" "a" "g" "g" "a" "t"
    [1549] "t" "t" "t" "a" "t" "c" "a" "a" "g" "a" "a" "a" "g" "c" "a" "g" "a" "t"
    [1567] "t" "t" "g" "g" "c" "a" "g" "t" "t" "c" "a" "a" "a" "a" "g" "a" "c" "t"
    [1585] "c" "c" "t" "g" "a" "a" "a" "t" "g" "a" "t" "a" "a" "a" "t" "c" "a" "g"
    [1603] "g" "g" "a" "a" "c" "t" "a" "a" "c" "c" "a" "a" "a" "c" "g" "g" "a" "g"
    [1621] "c" "a" "g" "a" "a" "t" "g" "g" "t" "c" "a" "a" "g" "t" "g" "a" "t" "g"
    [1639] "a" "a" "t" "a" "t" "t" "a" "c" "t" "a" "a" "t" "a" "g" "t" "g" "g" "t"
    [1657] "c" "a" "t" "g" "a" "g" "a" "a" "t" "a" "a" "a" "a" "c" "a" "a" "a" "a"
    [1675] "g" "g" "t" "g" "a" "t" "t" "c" "t" "a" "t" "t" "c" "a" "g" "a" "a" "t"
    [1693] "g" "a" "g" "a" "a" "a" "a" "a" "t" "c" "c" "t" "a" "a" "c" "c" "c" "a"
    [1711] "a" "t" "a" "g" "a" "a" "t" "c" "a" "c" "t" "c" "g" "a" "a" "a" "a" "a"
    [1729] "g" "a" "a" "t" "c" "t" "g" "c" "t" "t" "t" "c" "a" "a" "a" "a" "c" "g"
    [1747] "a" "a" "a" "g" "c" "t" "g" "a" "a" "c" "c" "t" "a" "t" "a" "a" "g" "c"
    [1765] "a" "g" "c" "a" "g" "t" "a" "t" "a" "a" "g" "c" "a" "a" "t" "a" "t" "g"
    [1783] "g" "a" "a" "c" "t" "c" "g" "a" "a" "t" "t" "a" "a" "a" "t" "a" "t" "c"
    [1801] "c" "a" "c" "a" "a" "t" "t" "c" "a" "a" "a" "a" "g" "c" "a" "c" "c" "t"
    [1819] "a" "a" "a" "a" "a" "g" "a" "a" "t" "a" "g" "g" "c" "t" "g" "a" "g" "g"
    [1837] "a" "g" "g" "a" "a" "g" "t" "c" "t" "t" "c" "t" "a" "c" "c" "a" "g" "g"
    [1855] "c" "a" "t" "a" "t" "t" "c" "a" "t" "g" "c" "g" "c" "t" "t" "g" "a" "a"
    [1873] "c" "t" "a" "g" "t" "a" "g" "t" "c" "a" "g" "t" "a" "g" "a" "a" "a" "t"
    [1891] "c" "t" "a" "a" "g" "c" "c" "c" "a" "c" "c" "t" "a" "a" "t" "t" "g" "t"
    [1909] "a" "c" "t" "g" "a" "a" "t" "t" "g" "c" "a" "a" "a" "t" "t" "g" "a" "t"
    [1927] "a" "g" "t" "t" "g" "t" "t" "c" "t" "a" "g" "c" "a" "g" "t" "g" "a" "a"
    [1945] "g" "a" "g" "a" "t" "a" "a" "a" "g" "a" "a" "a" "a" "a" "a" "a" "a" "g"
    [1963] "t" "a" "c" "a" "a" "c" "c" "a" "a" "a" "t" "g" "c" "c" "a" "g" "t" "c"
    [1981] "a" "g" "g" "c" "a" "c" "a" "g" "c" "a" "g" "a" "a" "a" "c" "c" "t" "a"
    [1999] "c" "a" "a" "c" "t" "c" "a" "t" "g" "g" "a" "a" "g" "g" "t" "a" "a" "a"
    [2017] "g" "a" "a" "c" "c" "t" "g" "c" "a" "a" "c" "t" "g" "g" "a" "g" "c" "c"
    [2035] "a" "a" "g" "a" "a" "g" "a" "g" "t" "a" "a" "c" "a" "a" "g" "c" "c" "a"
    [2053] "a" "a" "t" "g" "a" "a" "c" "a" "g" "a" "c" "a" "a" "g" "t" "a" "a" "a"
    [2071] "a" "g" "a" "c" "a" "t" "g" "a" "c" "a" "g" "t" "g" "a" "t" "a" "c" "t"
    [2089] "t" "t" "c" "c" "c" "a" "g" "a" "g" "c" "t" "g" "a" "a" "g" "t" "t" "a"
    [2107] "a" "c" "a" "a" "a" "t" "g" "c" "a" "c" "c" "t" "g" "g" "t" "t" "c" "t"
    [2125] "t" "t" "t" "a" "c" "t" "a" "a" "g" "t" "g" "t" "t" "c" "a" "a" "a" "t"
    [2143] "a" "c" "c" "a" "g" "t" "g" "a" "a" "c" "t" "t" "a" "a" "a" "g" "a" "a"
    [2161] "t" "t" "t" "g" "t" "c" "a" "a" "t" "c" "c" "t" "a" "g" "c" "c" "t" "t"
    [2179] "c" "c" "a" "a" "g" "a" "g" "a" "a" "g" "a" "a" "a" "a" "a" "g" "a" "a"
    [2197] "g" "a" "g" "a" "a" "a" "c" "t" "a" "g" "a" "a" "a" "c" "a" "g" "t" "t"
    [2215] "a" "a" "a" "g" "t" "g" "t" "c" "t" "a" "a" "t" "a" "a" "t" "g" "c" "t"
    [2233] "g" "a" "a" "g" "a" "c" "c" "c" "c" "a" "a" "a" "g" "a" "t" "c" "t" "c"
    [2251] "a" "t" "g" "t" "t" "a" "a" "g" "t" "g" "g" "a" "g" "a" "a" "a" "g" "g"
    [2269] "g" "t" "t" "t" "t" "g" "c" "a" "a" "a" "c" "t" "g" "a" "a" "a" "g" "a"
    [2287] "t" "c" "t" "g" "t" "a" "g" "a" "g" "a" "g" "t" "a" "g" "c" "a" "g" "t"
    [2305] "a" "t" "t" "t" "c" "a" "c" "t" "g" "g" "t" "a" "c" "c" "t" "g" "g" "t"
    [2323] "a" "c" "t" "g" "a" "t" "t" "a" "t" "g" "g" "c" "a" "c" "t" "c" "a" "g"
    [2341] "g" "a" "a" "a" "g" "t" "a" "t" "c" "t" "c" "g" "t" "t" "a" "c" "t" "g"
    [2359] "g" "a" "a" "g" "t" "t" "a" "g" "c" "a" "c" "t" "c" "t" "a" "g" "g" "g"
    [2377] "a" "a" "g" "g" "c" "a" "a" "a" "a" "a" "c" "a" "g" "a" "a" "c" "c" "a"
    [2395] "a" "a" "t" "a" "a" "a" "t" "g" "t" "g" "t" "g" "a" "g" "t" "c" "a" "g"
    [2413] "t" "g" "t" "g" "c" "a" "g" "c" "a" "t" "t" "t" "g" "a" "a" "a" "a" "c"
    [2431] "c" "c" "c" "a" "a" "g" "g" "g" "a" "c" "t" "a" "a" "t" "t" "c" "a" "t"
    [2449] "g" "g" "t" "t" "g" "t" "t" "c" "c" "a" "a" "a" "g" "a" "t" "a" "a" "t"
    [2467] "a" "g" "a" "a" "a" "t" "g" "a" "c" "a" "c" "a" "g" "a" "a" "g" "g" "c"
    [2485] "t" "t" "t" "a" "a" "g" "t" "a" "t" "c" "c" "a" "t" "t" "g" "g" "g" "a"
    [2503] "c" "a" "t" "g" "a" "a" "g" "t" "t" "a" "a" "c" "c" "a" "c" "a" "g" "t"
    [2521] "c" "g" "g" "g" "a" "a" "a" "c" "a" "a" "g" "c" "a" "t" "a" "g" "a" "a"
    [2539] "a" "t" "g" "g" "a" "a" "g" "a" "a" "a" "g" "t" "g" "a" "a" "c" "t" "t"
    [2557] "g" "a" "t" "g" "c" "t" "c" "a" "g" "t" "a" "t" "t" "t" "g" "c" "a" "g"
    [2575] "a" "a" "t" "a" "c" "a" "t" "t" "c" "a" "a" "g" "g" "t" "t" "t" "c" "a"
    [2593] "a" "a" "g" "c" "g" "c" "c" "a" "g" "t" "c" "a" "t" "t" "t" "g" "c" "t"
    [2611] "c" "t" "g" "t" "t" "t" "t" "c" "a" "a" "a" "t" "c" "c" "a" "g" "g" "a"
    [2629] "a" "a" "t" "g" "c" "a" "g" "a" "a" "g" "a" "g" "g" "a" "a" "t" "g" "t"
    [2647] "g" "c" "a" "a" "c" "a" "t" "t" "c" "t" "c" "t" "g" "c" "c" "c" "a" "c"
    [2665] "t" "c" "t" "g" "g" "g" "t" "c" "c" "t" "t" "a" "a" "a" "g" "a" "a" "a"
    [2683] "c" "a" "a" "a" "g" "t" "c" "c" "a" "a" "a" "a" "g" "t" "c" "a" "c" "t"
    [2701] "t" "t" "t" "g" "a" "a" "t" "g" "t" "g" "a" "a" "c" "a" "a" "a" "a" "g"
    [2719] "g" "a" "a" "g" "a" "a" "a" "a" "t" "c" "a" "a" "g" "g" "a" "a" "a" "g"
    [2737] "a" "a" "t" "g" "a" "g" "t" "c" "t" "a" "a" "t" "a" "t" "c" "a" "a" "g"
    [2755] "c" "c" "t" "g" "t" "a" "c" "a" "g" "a" "c" "a" "g" "t" "t" "a" "a" "t"
    [2773] "a" "t" "c" "a" "c" "t" "g" "c" "a" "g" "g" "c" "t" "t" "t" "c" "c" "t"
    [2791] "g" "t" "g" "g" "t" "t" "g" "g" "t" "c" "a" "g" "a" "a" "a" "g" "a" "t"
    [2809] "a" "a" "g" "c" "c" "a" "g" "t" "t" "g" "a" "t" "a" "a" "t" "g" "c" "c"
    [2827] "a" "a" "a" "t" "g" "t" "a" "g" "t" "a" "t" "c" "a" "a" "a" "g" "g" "a"
    [2845] "g" "g" "c" "t" "c" "t" "a" "g" "g" "t" "t" "t" "t" "g" "t" "c" "t" "a"
    [2863] "t" "c" "a" "t" "c" "t" "c" "a" "g" "t" "t" "c" "a" "g" "a" "g" "g" "c"
    [2881] "a" "a" "c" "g" "a" "a" "a" "c" "t" "g" "g" "a" "c" "t" "c" "a" "t" "t"
    [2899] "a" "c" "t" "c" "c" "a" "a" "a" "t" "a" "a" "a" "c" "a" "t" "g" "g" "a"
    [2917] "c" "t" "t" "t" "t" "a" "c" "a" "a" "a" "a" "c" "c" "c" "a" "t" "a" "t"
    [2935] "c" "g" "t" "a" "t" "a" "c" "c" "a" "c" "c" "a" "c" "t" "t" "t" "t" "t"
    [2953] "c" "c" "c" "a" "t" "c" "a" "a" "g" "t" "c" "a" "t" "t" "t" "g" "t" "t"
    [2971] "a" "a" "a" "a" "c" "t" "a" "a" "a" "t" "g" "t" "a" "a" "g" "a" "a" "a"
    [2989] "a" "a" "t" "c" "t" "g" "c" "t" "a" "g" "a" "g" "g" "a" "a" "a" "a" "c"
    [3007] "t" "t" "t" "g" "a" "g" "g" "a" "a" "c" "a" "t" "t" "c" "a" "a" "t" "g"
    [3025] "t" "c" "a" "c" "c" "t" "g" "a" "a" "a" "g" "a" "g" "a" "a" "a" "t" "g"
    [3043] "g" "g" "a" "a" "a" "t" "g" "a" "g" "a" "a" "c" "a" "t" "t" "c" "c" "a"
    [3061] "a" "g" "t" "a" "c" "a" "g" "t" "g" "a" "g" "c" "a" "c" "a" "a" "t" "t"
    [3079] "a" "g" "c" "c" "g" "t" "a" "a" "t" "a" "a" "c" "a" "t" "t" "a" "g" "a"
    [3097] "g" "a" "a" "a" "a" "t" "g" "t" "t" "t" "t" "t" "a" "a" "a" "g" "g" "a"
    [3115] "g" "c" "c" "a" "g" "c" "t" "c" "a" "a" "g" "c" "a" "a" "t" "a" "t" "t"
    [3133] "a" "a" "t" "g" "a" "a" "g" "t" "a" "g" "g" "t" "t" "c" "c" "a" "g" "t"
    [3151] "a" "c" "t" "a" "a" "t" "g" "a" "a" "g" "t" "g" "g" "g" "c" "t" "c" "c"
    [3169] "a" "g" "t" "a" "t" "t" "a" "a" "t" "g" "a" "a" "a" "t" "a" "g" "g" "t"
    [3187] "t" "c" "c" "a" "g" "t" "g" "a" "t" "g" "a" "a" "a" "a" "c" "a" "t" "t"
    [3205] "c" "a" "a" "g" "c" "a" "g" "a" "a" "c" "t" "a" "g" "g" "t" "a" "g" "a"
    [3223] "a" "a" "c" "a" "g" "a" "g" "g" "g" "c" "c" "a" "a" "a" "a" "t" "t" "g"
    [3241] "a" "a" "t" "g" "c" "t" "a" "t" "g" "c" "t" "t" "a" "g" "a" "t" "t" "a"
    [3259] "g" "g" "g" "g" "t" "t" "t" "t" "g" "c" "a" "a" "c" "c" "t" "g" "a" "g"
    [3277] "g" "t" "c" "t" "a" "t" "a" "a" "a" "c" "a" "a" "a" "g" "t" "c" "t" "t"
    [3295] "c" "c" "t" "g" "g" "a" "a" "g" "t" "a" "a" "t" "t" "g" "t" "a" "a" "g"
    [3313] "c" "a" "t" "c" "c" "t" "g" "a" "a" "a" "t" "a" "a" "a" "a" "a" "a" "g"
    [3331] "c" "a" "a" "g" "a" "a" "t" "a" "t" "g" "a" "a" "g" "a" "a" "g" "t" "a"
    [3349] "g" "t" "t" "c" "a" "g" "a" "c" "t" "g" "t" "t" "a" "a" "t" "a" "c" "a"
    [3367] "g" "a" "t" "t" "t" "c" "t" "c" "t" "c" "c" "a" "t" "a" "t" "c" "t" "g"
    [3385] "a" "t" "t" "t" "c" "a" "g" "a" "t" "a" "a" "c" "t" "t" "a" "g" "a" "a"
    [3403] "c" "a" "g" "c" "c" "t" "a" "t" "g" "g" "g" "a" "a" "g" "t" "a" "g" "t"
    [3421] "c" "a" "t" "g" "c" "a" "t" "c" "t" "c" "a" "g" "g" "t" "t" "t" "g" "t"
    [3439] "t" "c" "t" "g" "a" "g" "a" "c" "a" "c" "c" "t" "g" "a" "t" "g" "a" "c"
    [3457] "c" "t" "g" "t" "t" "a" "g" "a" "t" "g" "a" "t" "g" "g" "t" "g" "a" "a"
    [3475] "a" "t" "a" "a" "a" "g" "g" "a" "a" "g" "a" "t" "a" "c" "t" "a" "g" "t"
    [3493] "t" "t" "t" "g" "c" "t" "g" "a" "a" "a" "a" "t" "g" "a" "c" "a" "t" "t"
    [3511] "a" "a" "g" "g" "a" "a" "a" "g" "t" "t" "c" "t" "g" "c" "t" "g" "t" "t"
    [3529] "t" "t" "t" "a" "g" "c" "a" "a" "a" "a" "g" "c" "g" "t" "c" "c" "a" "g"
    [3547] "a" "g" "a" "g" "g" "a" "g" "a" "g" "c" "t" "t" "a" "g" "c" "a" "g" "g"
    [3565] "a" "g" "t" "c" "c" "t" "a" "g" "c" "c" "c" "t" "t" "t" "c" "a" "c" "c"
    [3583] "c" "a" "t" "a" "c" "a" "c" "a" "t" "t" "t" "g" "g" "c" "t" "c" "a" "g"
    [3601] "g" "g" "t" "t" "a" "c" "c" "g" "a" "a" "g" "a" "g" "g" "g" "g" "c" "c"
    [3619] "a" "a" "g" "a" "a" "a" "t" "t" "a" "g" "a" "g" "t" "c" "c" "t" "c" "a"
    [3637] "g" "a" "a" "g" "a" "g" "a" "a" "c" "t" "t" "a" "t" "c" "t" "a" "g" "t"
    [3655] "g" "a" "g" "g" "a" "t" "g" "a" "a" "g" "a" "g" "c" "t" "t" "c" "c" "c"
    [3673] "t" "g" "c" "t" "t" "c" "c" "a" "a" "c" "a" "c" "t" "t" "g" "t" "t" "a"
    [3691] "t" "t" "t" "g" "g" "t" "a" "a" "a" "g" "t" "a" "a" "a" "c" "a" "a" "t"
    [3709] "a" "t" "a" "c" "c" "t" "t" "c" "t" "c" "a" "g" "t" "c" "t" "a" "c" "t"
    [3727] "a" "g" "g" "c" "a" "t" "a" "g" "c" "a" "c" "c" "g" "t" "t" "g" "c" "t"
    [3745] "a" "c" "c" "g" "a" "g" "t" "g" "t" "c" "t" "g" "t" "c" "t" "a" "a" "g"
    [3763] "a" "a" "c" "a" "c" "a" "g" "a" "g" "g" "a" "g" "a" "a" "t" "t" "t" "a"
    [3781] "t" "t" "a" "t" "c" "a" "t" "t" "g" "a" "a" "g" "a" "a" "t" "a" "g" "c"
    [3799] "t" "t" "a" "a" "a" "t" "g" "a" "c" "t" "g" "c" "a" "g" "t" "a" "a" "c"
    [3817] "c" "a" "g" "g" "t" "a" "a" "t" "a" "t" "t" "g" "g" "c" "a" "a" "a" "g"
    [3835] "g" "c" "a" "t" "c" "t" "c" "a" "g" "g" "a" "a" "c" "a" "t" "c" "a" "c"
    [3853] "c" "t" "t" "a" "g" "t" "g" "a" "g" "g" "a" "a" "a" "c" "a" "a" "a" "a"
    [3871] "t" "g" "t" "t" "c" "t" "g" "c" "t" "a" "g" "c" "t" "t" "g" "t" "t" "t"
    [3889] "t" "c" "t" "t" "c" "a" "c" "a" "g" "t" "g" "c" "a" "g" "t" "g" "a" "a"
    [3907] "t" "t" "g" "g" "a" "a" "g" "a" "c" "t" "t" "g" "a" "c" "t" "g" "c" "a"
    [3925] "a" "a" "t" "a" "c" "a" "a" "a" "c" "a" "c" "c" "c" "a" "g" "g" "a" "t"
    [3943] "c" "c" "t" "t" "t" "c" "t" "t" "g" "a" "t" "t" "g" "g" "t" "t" "c" "t"
    [3961] "t" "c" "c" "a" "a" "a" "c" "a" "a" "a" "t" "g" "a" "g" "g" "c" "a" "t"
    [3979] "c" "a" "g" "t" "c" "t" "g" "a" "a" "a" "g" "c" "c" "a" "g" "g" "g" "a"
    [3997] "g" "t" "t" "g" "g" "t" "c" "t" "g" "a" "g" "t" "g" "a" "c" "a" "a" "g"
    [4015] "g" "a" "a" "t" "t" "g" "g" "t" "t" "t" "c" "a" "g" "a" "t" "g" "a" "t"
    [4033] "g" "a" "a" "g" "a" "a" "a" "g" "a" "g" "g" "a" "a" "c" "g" "g" "g" "c"
    [4051] "t" "t" "g" "g" "a" "a" "g" "a" "a" "a" "a" "t" "a" "a" "t" "c" "a" "a"
    [4069] "g" "a" "a" "g" "a" "g" "c" "a" "a" "a" "g" "c" "a" "t" "g" "g" "a" "t"
    [4087] "t" "c" "a" "a" "a" "c" "t" "t" "a" "g" "g" "t" "g" "a" "a" "g" "c" "a"
    [4105] "g" "c" "a" "t" "c" "t" "g" "g" "g" "t" "g" "t" "g" "a" "g" "a" "g" "t"
    [4123] "g" "a" "a" "a" "c" "a" "a" "g" "c" "g" "t" "c" "t" "c" "t" "g" "a" "a"
    [4141] "g" "a" "c" "t" "g" "c" "t" "c" "a" "g" "g" "g" "c" "t" "a" "t" "c" "c"
    [4159] "t" "c" "t" "c" "a" "g" "a" "g" "t" "g" "a" "c" "a" "t" "t" "t" "t" "a"
    [4177] "a" "c" "c" "a" "c" "t" "c" "a" "g" "c" "a" "g" "a" "g" "g" "g" "a" "t"
    [4195] "a" "c" "c" "a" "t" "g" "c" "a" "a" "c" "a" "t" "a" "a" "c" "c" "t" "g"
    [4213] "a" "t" "a" "a" "a" "g" "c" "t" "c" "c" "a" "g" "c" "a" "g" "g" "a" "a"
    [4231] "a" "t" "g" "g" "c" "t" "g" "a" "a" "c" "t" "a" "g" "a" "a" "g" "c" "t"
    [4249] "g" "t" "g" "t" "t" "a" "g" "a" "a" "c" "a" "g" "c" "a" "t" "g" "g" "g"
    [4267] "a" "g" "c" "c" "a" "g" "c" "c" "t" "t" "c" "t" "a" "a" "c" "a" "g" "c"
    [4285] "t" "a" "c" "c" "c" "t" "t" "c" "c" "a" "t" "c" "a" "t" "a" "a" "g" "t"
    [4303] "g" "a" "c" "t" "c" "c" "t" "c" "t" "g" "c" "c" "c" "t" "t" "g" "a" "g"
    [4321] "g" "a" "c" "c" "t" "g" "c" "g" "a" "a" "a" "t" "c" "c" "a" "g" "a" "a"
    [4339] "c" "a" "a" "a" "g" "c" "a" "c" "a" "t" "c" "a" "g" "a" "a" "a" "a" "a"
    [4357] "g" "c" "a" "g" "t" "a" "t" "t" "a" "a" "c" "t" "t" "c" "a" "c" "a" "g"
    [4375] "a" "a" "a" "a" "g" "t" "a" "g" "t" "g" "a" "a" "t" "a" "c" "c" "c" "t"
    [4393] "a" "t" "a" "a" "g" "c" "c" "a" "g" "a" "a" "t" "c" "c" "a" "g" "a" "a"
    [4411] "g" "g" "c" "c" "t" "t" "t" "c" "t" "g" "c" "t" "g" "a" "c" "a" "a" "g"
    [4429] "t" "t" "t" "g" "a" "g" "g" "t" "g" "t" "c" "t" "g" "c" "a" "g" "a" "t"
    [4447] "a" "g" "t" "t" "c" "t" "a" "c" "c" "a" "g" "t" "a" "a" "a" "a" "a" "t"
    [4465] "a" "a" "a" "g" "a" "a" "c" "c" "a" "g" "g" "a" "g" "t" "g" "g" "a" "a"
    [4483] "a" "g" "g" "t" "c" "a" "t" "c" "c" "c" "c" "t" "t" "c" "t" "a" "a" "a"
    [4501] "t" "g" "c" "c" "c" "a" "t" "c" "a" "t" "t" "a" "g" "a" "t" "g" "a" "t"
    [4519] "a" "g" "g" "t" "g" "g" "t" "a" "c" "a" "t" "g" "c" "a" "c" "a" "g" "t"
    [4537] "t" "g" "c" "t" "c" "t" "g" "g" "g" "a" "g" "t" "c" "t" "t" "c" "a" "g"
    [4555] "a" "a" "t" "a" "g" "a" "a" "a" "c" "t" "a" "c" "c" "c" "a" "t" "c" "t"
    [4573] "c" "a" "a" "g" "a" "g" "g" "a" "g" "c" "t" "c" "a" "t" "t" "a" "a" "g"
    [4591] "g" "t" "t" "g" "t" "t" "g" "a" "t" "g" "t" "g" "g" "a" "g" "g" "a" "g"
    [4609] "c" "a" "a" "c" "a" "g" "c" "t" "g" "g" "a" "a" "g" "a" "g" "t" "c" "t"
    [4627] "g" "g" "g" "c" "c" "a" "c" "a" "c" "g" "a" "t" "t" "t" "g" "a" "c" "g"
    [4645] "g" "a" "a" "a" "c" "a" "t" "c" "t" "t" "a" "c" "t" "t" "g" "c" "c" "a"
    [4663] "a" "g" "g" "c" "a" "a" "g" "a" "t" "c" "t" "a" "g" "a" "g" "g" "g" "a"
    [4681] "a" "c" "c" "c" "c" "t" "t" "a" "c" "c" "t" "g" "g" "a" "a" "t" "c" "t"
    [4699] "g" "g" "a" "a" "t" "c" "a" "g" "c" "c" "t" "c" "t" "t" "c" "t" "c" "t"
    [4717] "g" "a" "t" "g" "a" "c" "c" "c" "t" "g" "a" "a" "t" "c" "t" "g" "a" "t"
    [4735] "c" "c" "t" "t" "c" "t" "g" "a" "a" "g" "a" "c" "a" "g" "a" "g" "c" "c"
    [4753] "c" "c" "a" "g" "a" "g" "t" "c" "a" "g" "c" "t" "c" "g" "t" "g" "t" "t"
    [4771] "g" "g" "c" "a" "a" "c" "a" "t" "a" "c" "c" "a" "t" "c" "t" "t" "c" "a"
    [4789] "a" "c" "c" "t" "c" "t" "g" "c" "a" "t" "t" "g" "a" "a" "a" "g" "t" "t"
    [4807] "c" "c" "c" "c" "a" "a" "t" "t" "g" "a" "a" "a" "g" "t" "t" "g" "c" "a"
    [4825] "g" "a" "a" "t" "c" "t" "g" "c" "c" "c" "a" "g" "g" "g" "t" "c" "c" "a"
    [4843] "g" "c" "t" "g" "c" "t" "g" "c" "t" "c" "a" "t" "a" "c" "t" "a" "c" "t"
    [4861] "g" "a" "t" "a" "c" "t" "g" "c" "t" "g" "g" "g" "t" "a" "t" "a" "a" "t"
    [4879] "g" "c" "a" "a" "t" "g" "g" "a" "a" "g" "a" "a" "a" "g" "t" "g" "t" "g"
    [4897] "a" "g" "c" "a" "g" "g" "g" "a" "g" "a" "a" "g" "c" "c" "a" "g" "a" "a"
    [4915] "t" "t" "g" "a" "c" "a" "g" "c" "t" "t" "c" "a" "a" "c" "a" "g" "a" "a"
    [4933] "a" "g" "g" "g" "t" "c" "a" "a" "c" "a" "a" "a" "a" "g" "a" "a" "t" "g"
    [4951] "t" "c" "c" "a" "t" "g" "g" "t" "g" "g" "t" "g" "t" "c" "t" "g" "g" "c"
    [4969] "c" "t" "g" "a" "c" "c" "c" "c" "a" "g" "a" "a" "g" "a" "a" "t" "t" "t"
    [4987] "a" "t" "g" "c" "t" "c" "g" "t" "g" "t" "a" "c" "a" "a" "g" "t" "t" "t"
    [5005] "g" "c" "c" "a" "g" "a" "a" "a" "a" "c" "a" "c" "c" "a" "c" "a" "t" "c"
    [5023] "a" "c" "t" "t" "t" "a" "a" "c" "t" "a" "a" "t" "c" "t" "a" "a" "t" "t"
    [5041] "a" "c" "t" "g" "a" "a" "g" "a" "g" "a" "c" "t" "a" "c" "t" "c" "a" "t"
    [5059] "g" "t" "t" "g" "t" "t" "a" "t" "g" "a" "a" "a" "a" "c" "a" "g" "a" "t"
    [5077] "g" "c" "t" "g" "a" "g" "t" "t" "t" "g" "t" "g" "t" "g" "t" "g" "a" "a"
    [5095] "c" "g" "g" "a" "c" "a" "c" "t" "g" "a" "a" "a" "t" "a" "t" "t" "t" "t"
    [5113] "c" "t" "a" "g" "g" "a" "a" "t" "t" "g" "c" "g" "g" "g" "a" "g" "g" "a"
    [5131] "a" "a" "a" "t" "g" "g" "g" "t" "a" "g" "t" "t" "a" "g" "c" "t" "a" "t"
    [5149] "t" "t" "c" "t" "g" "g" "g" "t" "g" "a" "c" "c" "c" "a" "g" "t" "c" "t"
    [5167] "a" "t" "t" "a" "a" "a" "g" "a" "a" "a" "g" "a" "a" "a" "a" "a" "t" "g"
    [5185] "c" "t" "g" "a" "a" "t" "g" "a" "g" "c" "a" "t" "g" "a" "t" "t" "t" "t"
    [5203] "g" "a" "a" "g" "t" "c" "a" "g" "a" "g" "g" "a" "g" "a" "t" "g" "t" "g"
    [5221] "g" "t" "c" "a" "a" "t" "g" "g" "a" "a" "g" "a" "a" "a" "c" "c" "a" "c"
    [5239] "c" "a" "a" "g" "g" "t" "c" "c" "a" "a" "a" "g" "c" "g" "a" "g" "c" "a"
    [5257] "a" "g" "a" "g" "a" "a" "t" "c" "c" "c" "a" "g" "g" "a" "c" "a" "g" "a"
    [5275] "a" "a" "g" "a" "t" "c" "t" "t" "c" "a" "g" "g" "g" "g" "g" "c" "t" "a"
    [5293] "g" "a" "a" "a" "t" "c" "t" "g" "t" "t" "g" "c" "t" "a" "t" "g" "g" "g"
    [5311] "c" "c" "c" "t" "t" "c" "a" "c" "c" "a" "a" "c" "a" "t" "g" "c" "c" "c"
    [5329] "a" "c" "a" "g" "a" "t" "c" "a" "a" "c" "t" "g" "g" "a" "a" "t" "g" "g"
    [5347] "a" "t" "g" "g" "t" "a" "c" "a" "g" "c" "t" "g" "t" "g" "t" "g" "g" "t"
    [5365] "g" "c" "t" "t" "c" "t" "g" "t" "g" "g" "t" "g" "a" "a" "g" "g" "a" "g"
    [5383] "c" "t" "t" "t" "c" "a" "t" "c" "a" "t" "t" "c" "a" "c" "c" "c" "t" "t"
    [5401] "g" "g" "c" "a" "c" "a" "g" "g" "t" "g" "t" "c" "c" "a" "c" "c" "c" "a"
    [5419] "a" "t" "t" "g" "t" "g" "g" "t" "t" "g" "t" "g" "c" "a" "g" "c" "c" "a"
    [5437] "g" "a" "t" "g" "c" "c" "t" "g" "g" "a" "c" "a" "g" "a" "g" "g" "a" "c"
    [5455] "a" "a" "t" "g" "g" "c" "t" "t" "c" "c" "a" "t" "g" "c" "a" "a" "t" "t"
    [5473] "g" "g" "g" "c" "a" "g" "a" "t" "g" "t" "g" "t" "g" "a" "g" "g" "c" "a"
    [5491] "c" "c" "t" "g" "t" "g" "g" "t" "g" "a" "c" "c" "c" "g" "a" "g" "a" "g"
    [5509] "t" "g" "g" "g" "t" "g" "t" "t" "g" "g" "a" "c" "a" "g" "t" "g"
    attr(,"name")
    [1] "AY751490"
    attr(,"Annot")
    [1] ">AY751490"
    attr(,"class")
    [1] "SeqFastadna"
    
    $BC030969.BRCA1
       [1] "a" "t" "g" "g" "a" "t" "t" "t" "a" "t" "c" "t" "g" "c" "t" "c" "t" "t"
      [19] "c" "g" "c" "g" "t" "t" "g" "a" "a" "g" "a" "a" "g" "t" "a" "c" "a" "a"
      [37] "a" "a" "t" "g" "t" "c" "a" "t" "t" "a" "a" "t" "g" "c" "t" "a" "t" "g"
      [55] "c" "a" "g" "a" "a" "a" "a" "t" "c" "t" "t" "a" "g" "a" "g" "t" "g" "t"
      [73] "c" "c" "c" "a" "t" "c" "t" "g" "t" "c" "t" "g" "g" "a" "g" "t" "t" "g"
      [91] "a" "t" "c" "a" "a" "g" "g" "a" "a" "c" "c" "t" "g" "t" "c" "t" "c" "c"
     [109] "a" "c" "a" "a" "a" "g" "t" "g" "t" "g" "a" "c" "c" "a" "c" "a" "t" "a"
     [127] "t" "t" "t" "t" "g" "c" "a" "a" "a" "t" "t" "t" "t" "g" "c" "a" "t" "g"
     [145] "c" "t" "g" "a" "a" "a" "c" "t" "t" "c" "t" "c" "a" "a" "c" "c" "a" "g"
     [163] "a" "a" "g" "a" "a" "a" "g" "g" "g" "c" "c" "t" "t" "c" "a" "c" "a" "g"
     [181] "t" "g" "t" "c" "c" "t" "t" "t" "a" "t" "g" "t" "a" "a" "g" "a" "a" "t"
     [199] "g" "a" "t" "a" "t" "a" "a" "c" "c" "a" "a" "a" "a" "g" "g" "a" "g" "c"
     [217] "c" "t" "a" "c" "a" "a" "g" "a" "a" "a" "g" "t" "a" "c" "g" "a" "g" "a"
     [235] "t" "t" "t" "a" "g" "t" "c" "a" "a" "c" "t" "t" "g" "t" "t" "g" "a" "a"
     [253] "g" "a" "g" "c" "t" "a" "t" "t" "g" "a" "a" "a" "a" "t" "c" "a" "t" "t"
     [271] "t" "g" "t" "g" "c" "t" "t" "t" "t" "c" "a" "g" "c" "t" "t" "g" "a" "c"
     [289] "a" "c" "a" "g" "g" "t" "t" "t" "g" "g" "a" "g" "t" "a" "t" "g" "c" "a"
     [307] "a" "a" "c" "a" "g" "c" "t" "a" "t" "a" "a" "t" "t" "t" "t" "g" "c" "a"
     [325] "a" "a" "a" "a" "a" "g" "g" "a" "a" "a" "a" "t" "a" "a" "c" "t" "c" "t"
     [343] "c" "c" "t" "g" "a" "a" "c" "a" "t" "c" "t" "a" "a" "a" "a" "g" "a" "t"
     [361] "g" "a" "a" "g" "t" "t" "t" "c" "t" "a" "t" "c" "a" "t" "c" "c" "a" "a"
     [379] "a" "g" "t" "a" "t" "g" "g" "g" "c" "t" "a" "c" "a" "g" "a" "a" "a" "c"
     [397] "c" "g" "t" "g" "c" "c" "a" "a" "a" "a" "g" "a" "c" "t" "t" "c" "t" "a"
     [415] "c" "a" "g" "a" "g" "t" "g" "a" "a" "c" "c" "c" "g" "a" "a" "a" "a" "t"
     [433] "c" "c" "t" "t" "c" "c" "t" "t" "g" "c" "a" "g" "g" "a" "a" "a" "c" "c"
     [451] "a" "g" "t" "c" "t" "c" "a" "g" "t" "g" "t" "c" "c" "a" "a" "c" "t" "c"
     [469] "t" "c" "t" "a" "a" "c" "c" "t" "t" "g" "g" "a" "a" "c" "t" "g" "t" "g"
     [487] "a" "g" "a" "a" "c" "t" "c" "t" "g" "a" "g" "g" "a" "c" "a" "a" "a" "g"
     [505] "c" "a" "g" "c" "g" "g" "a" "t" "a" "c" "a" "a" "c" "c" "t" "c" "a" "a"
     [523] "a" "a" "g" "a" "c" "g" "t" "c" "t" "g" "t" "c" "t" "a" "c" "a" "t" "t"
     [541] "g" "a" "a" "t" "t" "g" "g" "g" "a" "t" "c" "t" "g" "a" "t" "t" "c" "t"
     [559] "t" "c" "t" "g" "a" "a" "g" "a" "t" "a" "c" "c" "g" "t" "t" "a" "a" "t"
     [577] "a" "a" "g" "g" "c" "a" "a" "c" "t" "t" "a" "t" "t" "g" "c" "a" "g" "t"
     [595] "g" "t" "g" "g" "g" "a" "g" "a" "t" "c" "a" "a" "g" "a" "a" "t" "t" "g"
     [613] "t" "t" "a" "c" "a" "a" "a" "t" "c" "a" "c" "c" "c" "c" "t" "c" "a" "a"
     [631] "g" "g" "a" "a" "c" "c" "a" "g" "g" "g" "a" "t" "g" "a" "a" "a" "t" "c"
     [649] "a" "g" "t" "t" "t" "g" "g" "a" "t" "t" "c" "t" "g" "c" "a" "a" "a" "a"
     [667] "a" "a" "g" "g" "c" "t" "g" "c" "t" "t" "g" "t" "g" "a" "a" "t" "t" "t"
     [685] "t" "c" "t" "g" "a" "g" "a" "c" "g" "g" "a" "t" "g" "t" "a" "a" "c" "a"
     [703] "a" "a" "t" "a" "c" "t" "g" "a" "a" "c" "a" "t" "c" "a" "t" "c" "a" "a"
     [721] "c" "c" "c" "a" "g" "t" "a" "a" "t" "a" "a" "t" "g" "a" "t" "t" "t" "g"
     [739] "a" "a" "c" "a" "c" "c" "a" "c" "t" "g" "a" "g" "a" "a" "g" "c" "g" "t"
     [757] "g" "c" "a" "g" "c" "t" "g" "a" "g" "a" "g" "g" "c" "a" "t" "c" "c" "a"
     [775] "g" "a" "a" "a" "a" "g" "t" "a" "t" "c" "a" "g" "g" "g" "t" "a" "g" "t"
     [793] "t" "c" "t" "g" "t" "t" "t" "c" "a" "a" "a" "c" "t" "t" "g" "c" "a" "t"
     [811] "g" "t" "g" "g" "a" "g" "c" "c" "a" "t" "g" "t" "g" "g" "c" "a" "c" "a"
     [829] "a" "a" "t" "a" "c" "t" "c" "a" "t" "g" "c" "c" "a" "g" "c" "t" "c" "a"
     [847] "t" "t" "a" "c" "a" "g" "c" "a" "t" "g" "a" "g" "a" "a" "c" "a" "g" "c"
     [865] "a" "g" "t" "t" "t" "a" "t" "t" "a" "c" "t" "c" "a" "c" "t" "a" "a" "a"
     [883] "g" "a" "c" "a" "g" "a" "a" "t" "g" "a" "a" "t" "g" "t" "a" "g" "a" "a"
     [901] "a" "a" "g" "g" "c" "t" "g" "a" "a" "t" "t" "c" "t" "g" "t" "a" "a" "t"
     [919] "a" "a" "a" "a" "g" "c" "a" "a" "a" "c" "a" "g" "c" "c" "t" "g" "g" "c"
     [937] "t" "t" "a" "g" "c" "a" "a" "g" "g" "a" "g" "c" "c" "a" "a" "c" "a" "t"
     [955] "a" "a" "c" "a" "g" "a" "t" "g" "g" "g" "c" "t" "g" "g" "a" "a" "g" "t"
     [973] "a" "a" "g" "g" "a" "a" "a" "c" "a" "t" "g" "t" "a" "a" "t" "g" "a" "t"
     [991] "a" "g" "g" "c" "g" "g" "a" "c" "t" "c" "c" "c" "a" "g" "c" "a" "c" "a"
    [1009] "g" "a" "a" "a" "a" "a" "a" "a" "g" "g" "t" "a" "g" "a" "t" "c" "t" "g"
    [1027] "a" "a" "t" "g" "c" "t" "g" "a" "t" "c" "c" "c" "c" "t" "g" "t" "g" "t"
    [1045] "g" "a" "g" "a" "g" "a" "a" "a" "a" "g" "a" "a" "t" "g" "g" "a" "a" "t"
    [1063] "a" "a" "g" "c" "a" "g" "a" "a" "a" "c" "t" "g" "c" "c" "a" "t" "g" "c"
    [1081] "t" "c" "a" "g" "a" "g" "a" "a" "t" "c" "c" "t" "a" "g" "a" "g" "a" "t"
    [1099] "a" "c" "t" "g" "a" "a" "g" "a" "t" "g" "t" "t" "c" "c" "t" "t" "g" "g"
    [1117] "a" "t" "a" "a" "c" "a" "c" "t" "a" "a" "a" "t" "a" "g" "c" "a" "g" "c"
    [1135] "a" "t" "t" "c" "a" "g" "a" "a" "a" "g" "t" "t" "a" "a" "t" "g" "a" "g"
    [1153] "t" "g" "g" "t" "t" "t" "t" "c" "c" "a" "g" "a" "a" "g" "t" "g" "a" "t"
    [1171] "g" "a" "a" "c" "t" "g" "t" "t" "a" "g" "g" "t" "t" "c" "t" "g" "a" "t"
    [1189] "g" "a" "c" "t" "c" "a" "c" "a" "t" "g" "a" "t" "g" "g" "g" "g" "a" "g"
    [1207] "t" "c" "t" "g" "a" "a" "t" "c" "a" "a" "a" "t" "g" "c" "c" "a" "a" "a"
    [1225] "g" "t" "a" "g" "c" "t" "g" "a" "t" "g" "t" "a" "t" "t" "g" "g" "a" "c"
    [1243] "g" "t" "t" "c" "t" "a" "a" "a" "t" "g" "a" "g" "g" "t" "a" "g" "a" "t"
    [1261] "g" "a" "a" "t" "a" "t" "t" "c" "t" "g" "g" "t" "t" "c" "t" "t" "c" "a"
    [1279] "g" "a" "g" "a" "a" "a" "a" "t" "a" "g" "a" "c" "t" "t" "a" "c" "t" "g"
    [1297] "g" "c" "c" "a" "g" "t" "g" "a" "t" "c" "c" "t" "c" "a" "t" "g" "a" "g"
    [1315] "g" "c" "t" "t" "t" "a" "a" "t" "a" "t" "g" "t" "a" "a" "a" "a" "g" "t"
    [1333] "g" "a" "a" "a" "g" "a" "g" "t" "t" "c" "a" "c" "t" "c" "c" "a" "a" "a"
    [1351] "t" "c" "a" "g" "t" "a" "g" "a" "g" "a" "g" "t" "a" "a" "t" "a" "t" "t"
    [1369] "g" "a" "a" "g" "a" "c" "a" "a" "a" "a" "t" "a" "t" "t" "t" "g" "g" "g"
    [1387] "a" "a" "a" "a" "c" "c" "t" "a" "t" "c" "g" "g" "a" "a" "g" "a" "a" "g"
    [1405] "g" "c" "a" "a" "g" "c" "c" "t" "c" "c" "c" "c" "a" "a" "c" "t" "t" "a"
    [1423] "a" "g" "c" "c" "a" "t" "g" "t" "a" "a" "c" "t" "g" "a" "a" "a" "a" "t"
    [1441] "c" "t" "a" "a" "t" "t" "a" "t" "a" "g" "g" "a" "g" "c" "a" "t" "t" "t"
    [1459] "g" "t" "t" "a" "c" "t" "g" "a" "g" "c" "c" "a" "c" "a" "g" "a" "t" "a"
    [1477] "a" "t" "a" "c" "a" "a" "g" "a" "g" "c" "g" "t" "c" "c" "c" "c" "t" "c"
    [1495] "a" "c" "a" "a" "a" "t" "a" "a" "a" "t" "t" "a" "a" "a" "g" "c" "g" "t"
    [1513] "a" "a" "a" "a" "g" "g" "a" "g" "a" "c" "c" "t" "a" "c" "a" "t" "c" "a"
    [1531] "g" "g" "c" "c" "t" "t" "c" "a" "t" "c" "c" "t" "g" "a" "g" "g" "a" "t"
    [1549] "t" "t" "t" "a" "t" "c" "a" "a" "g" "a" "a" "a" "g" "c" "a" "g" "a" "t"
    [1567] "t" "t" "g" "g" "c" "a" "g" "t" "t" "c" "a" "a" "a" "a" "g" "a" "c" "t"
    [1585] "c" "c" "t" "g" "a" "a" "a" "t" "g" "a" "t" "a" "a" "a" "t" "c" "a" "g"
    [1603] "g" "g" "a" "a" "c" "t" "a" "a" "c" "c" "a" "a" "a" "c" "g" "g" "a" "g"
    [1621] "c" "a" "g" "a" "a" "t" "g" "g" "t" "c" "a" "a" "g" "t" "g" "a" "t" "g"
    [1639] "a" "a" "t" "a" "t" "t" "a" "c" "t" "a" "a" "t" "a" "g" "t" "g" "g" "t"
    [1657] "c" "a" "t" "g" "a" "g" "a" "a" "t" "a" "a" "a" "a" "c" "a" "a" "a" "a"
    [1675] "g" "g" "t" "g" "a" "t" "t" "c" "t" "a" "t" "t" "c" "a" "g" "a" "a" "t"
    [1693] "g" "a" "g" "a" "a" "a" "a" "a" "t" "c" "c" "t" "a" "a" "c" "c" "c" "a"
    [1711] "a" "t" "a" "g" "a" "a" "t" "c" "a" "c" "t" "c" "g" "a" "a" "a" "a" "a"
    [1729] "g" "a" "a" "t" "c" "t" "g" "c" "t" "t" "t" "c" "a" "a" "a" "a" "c" "g"
    [1747] "a" "a" "a" "g" "c" "t" "g" "a" "a" "c" "c" "t" "a" "t" "a" "a" "g" "c"
    [1765] "a" "g" "c" "a" "g" "t" "a" "t" "a" "a" "g" "c" "a" "a" "t" "a" "t" "g"
    [1783] "g" "a" "a" "c" "t" "c" "g" "a" "a" "t" "t" "a" "a" "a" "t" "a" "t" "c"
    [1801] "c" "a" "c" "a" "a" "t" "t" "c" "a" "a" "a" "a" "g" "c" "a" "c" "c" "t"
    [1819] "a" "a" "a" "a" "a" "g" "a" "a" "t" "a" "g" "g" "c" "t" "g" "a" "g" "g"
    [1837] "a" "g" "g" "a" "a" "g" "t" "c" "t" "t" "c" "t" "a" "c" "c" "a" "g" "g"
    [1855] "c" "a" "t" "a" "t" "t" "c" "a" "t" "g" "c" "g" "c" "t" "t" "g" "a" "a"
    [1873] "c" "t" "a" "g" "t" "a" "g" "t" "c" "a" "g" "t" "a" "g" "a" "a" "a" "t"
    [1891] "c" "t" "a" "a" "g" "c" "c" "c" "a" "c" "c" "t" "a" "a" "t" "t" "g" "t"
    [1909] "a" "c" "t" "g" "a" "a" "t" "t" "g" "c" "a" "a" "a" "t" "t" "g" "a" "t"
    [1927] "a" "g" "t" "t" "g" "t" "t" "c" "t" "a" "g" "c" "a" "g" "t" "a" "a" "a"
    [1945] "a" "a" "a" "a" "a" "a" "a" "a" "a" "a" "a" "a" "a" "a" "a" "a" "a" "a"
    [1963] "a" "a" "a" "a" "a" "a" "a" "a" "a"
    attr(,"name")
    [1] "BC030969.BRCA1"
    attr(,"Annot")
    [1] ">BC030969.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $BC062429
       [1] "g" "g" "g" "g" "t" "t" "t" "g" "g" "a" "t" "t" "c" "t" "g" "c" "a" "a"
      [19] "a" "a" "a" "a" "g" "g" "c" "t" "g" "c" "t" "t" "g" "t" "g" "a" "a" "t"
      [37] "t" "t" "t" "c" "t" "g" "a" "g" "a" "c" "g" "g" "a" "t" "g" "t" "a" "a"
      [55] "c" "a" "a" "a" "t" "a" "c" "t" "g" "a" "a" "c" "a" "t" "c" "a" "t" "c"
      [73] "a" "a" "c" "c" "c" "a" "g" "t" "a" "a" "t" "a" "a" "t" "g" "a" "t" "t"
      [91] "t" "g" "a" "a" "c" "a" "c" "c" "a" "c" "t" "g" "a" "g" "a" "a" "g" "c"
     [109] "g" "t" "g" "c" "a" "g" "c" "t" "g" "a" "g" "a" "g" "g" "c" "a" "t" "c"
     [127] "c" "a" "g" "a" "a" "a" "a" "g" "t" "a" "t" "c" "a" "g" "g" "g" "t" "a"
     [145] "g" "t" "t" "c" "t" "g" "t" "t" "t" "c" "a" "a" "a" "c" "t" "t" "g" "c"
     [163] "a" "t" "g" "t" "g" "g" "a" "g" "c" "c" "a" "t" "g" "t" "g" "g" "c" "a"
     [181] "c" "a" "a" "a" "t" "a" "c" "t" "c" "a" "t" "g" "c" "c" "a" "g" "c" "t"
     [199] "c" "a" "t" "t" "a" "c" "a" "g" "c" "a" "t" "g" "a" "g" "a" "a" "c" "a"
     [217] "g" "c" "a" "g" "t" "t" "t" "a" "t" "t" "a" "c" "t" "c" "a" "c" "t" "a"
     [235] "a" "a" "g" "a" "c" "a" "g" "a" "a" "t" "g" "a" "a" "t" "g" "t" "a" "g"
     [253] "a" "a" "a" "a" "g" "g" "c" "t" "g" "a" "a" "t" "t" "c" "t" "g" "t" "a"
     [271] "a" "t" "a" "a" "a" "a" "g" "c" "a" "a" "a" "c" "a" "g" "c" "c" "t" "g"
     [289] "g" "c" "t" "t" "a" "g" "c" "a" "a" "g" "g" "a" "g" "c" "c" "a" "a" "c"
     [307] "a" "t" "a" "a" "c" "a" "g" "a" "t" "g" "g" "g" "c" "t" "g" "g" "a" "a"
     [325] "g" "t" "a" "a" "g" "g" "a" "a" "a" "c" "a" "t" "g" "t" "a" "a" "t" "g"
     [343] "a" "t" "a" "g" "g" "c" "g" "g" "a" "c" "t" "c" "c" "c" "a" "g" "c" "a"
     [361] "c" "a" "g" "a" "a" "a" "a" "a" "a" "a" "g" "g" "t" "a" "g" "a" "t" "c"
     [379] "t" "g" "a" "a" "t" "g" "c" "t" "g" "a" "t" "c" "c" "c" "c" "t" "g" "t"
     [397] "g" "t" "g" "a" "g" "a" "g" "a" "a" "a" "a" "g" "a" "a" "t" "g" "g" "a"
     [415] "a" "t" "a" "a" "g" "c" "g" "g" "a" "a" "a" "c" "t" "g" "c" "c" "a" "t"
     [433] "g" "c" "t" "c" "a" "g" "a" "g" "a" "a" "t" "c" "c" "t" "a" "g" "a" "g"
     [451] "a" "t" "a" "c" "t" "g" "a" "a" "g" "a" "t" "g" "t" "t" "c" "c" "t" "t"
     [469] "g" "g" "a" "t" "a" "a" "c" "a" "c" "t" "a" "a" "a" "t" "a" "g" "c" "a"
     [487] "g" "c" "a" "t" "t" "c" "a" "g" "a" "a" "a" "g" "t" "t" "a" "a" "t" "g"
     [505] "a" "g" "t" "g" "g" "t" "t" "t" "t" "c" "c" "a" "g" "a" "a" "g" "t" "g"
     [523] "a" "t" "g" "a" "a" "c" "t" "g" "t" "t" "a" "g" "g" "t" "t" "c" "t" "g"
     [541] "a" "t" "g" "a" "c" "t" "c" "a" "c" "a" "t" "g" "a" "t" "g" "g" "g" "g"
     [559] "a" "g" "t" "c" "t" "g" "a" "a" "t" "c" "a" "a" "a" "t" "g" "c" "c" "a"
     [577] "a" "a" "g" "t" "a" "g" "c" "t" "g" "a" "t" "g" "t" "a" "t" "t" "g" "g"
     [595] "a" "c" "g" "t" "t" "c" "t" "a" "a" "a" "t" "g" "a" "g" "g" "t" "a" "g"
     [613] "a" "t" "g" "a" "a" "t" "a" "t" "t" "c" "t" "g" "g" "t" "t" "c" "t" "t"
     [631] "c" "a" "g" "a" "g" "a" "a" "a" "a" "t" "a" "g" "a" "c" "t" "t" "a" "c"
     [649] "t" "g" "g" "c" "c" "a" "g" "t" "g" "a" "t" "c" "c" "t" "c" "a" "t" "g"
     [667] "a" "g" "g" "c" "t" "t" "t" "a" "a" "t" "a" "t" "g" "t" "a" "a" "a" "a"
     [685] "g" "t" "g" "a" "a" "a" "g" "a" "g" "t" "t" "c" "a" "c" "t" "c" "c" "a"
     [703] "a" "a" "t" "c" "a" "g" "t" "a" "g" "a" "g" "a" "g" "t" "a" "a" "t" "a"
     [721] "t" "t" "g" "a" "a" "g" "a" "c" "a" "a" "a" "a" "t" "a" "t" "t" "t" "g"
     [739] "g" "g" "a" "a" "a" "a" "c" "c" "t" "a" "t" "c" "g" "g" "a" "a" "g" "a"
     [757] "a" "g" "g" "c" "a" "a" "g" "c" "c" "t" "c" "c" "c" "c" "a" "a" "c" "t"
     [775] "t" "a" "a" "g" "c" "c" "a" "t" "g" "t" "a" "a" "c" "t" "g" "a" "a" "a"
     [793] "a" "t" "c" "t" "a" "a" "t" "t" "a" "t" "a" "g" "g" "a" "g" "c" "a" "t"
     [811] "t" "t" "g" "t" "t" "a" "c" "t" "g" "a" "g" "c" "c" "a" "c" "a" "g" "a"
     [829] "t" "a" "a" "t" "a" "c" "a" "a" "g" "a" "g" "c" "g" "t" "c" "c" "c" "c"
     [847] "t" "c" "a" "c" "a" "a" "a" "t" "a" "a" "a" "t" "t" "a" "a" "a" "g" "c"
     [865] "g" "t" "a" "a" "a" "a" "g" "g" "a" "g" "a" "c" "c" "c" "a" "c" "a" "t"
     [883] "c" "a" "g" "g" "c" "c" "t" "t" "c" "a" "t" "c" "c" "t" "g" "a" "g" "g"
     [901] "a" "t" "t" "t" "t" "a" "t" "c" "a" "a" "g" "a" "a" "a" "g" "c" "a" "g"
     [919] "a" "t" "t" "t" "g" "g" "c" "a" "g" "t" "t" "c" "a" "a" "a" "a" "g" "a"
     [937] "c" "t" "c" "c" "t" "g" "a" "a" "a" "t" "g" "a" "t" "a" "a" "a" "t" "c"
     [955] "a" "g" "g" "g" "a" "a" "c" "t" "a" "a" "c" "c" "a" "a" "a" "c" "g" "g"
     [973] "a" "g" "c" "a" "g" "a" "a" "t" "g" "g" "t" "c" "a" "a" "g" "t" "g" "a"
     [991] "t" "g" "a" "a" "t" "a" "t" "t" "a" "c" "t" "a" "a" "t" "a" "g" "t" "g"
    [1009] "g" "t" "c" "a" "t" "g" "a" "g" "a" "a" "t" "a" "a" "a" "a" "c" "a" "a"
    [1027] "a" "a" "g" "g" "t" "g" "a" "t" "t" "c" "t" "a" "t" "t" "c" "a" "g" "a"
    [1045] "a" "t" "g" "a" "g" "a" "a" "a" "a" "a" "t" "c" "c" "t" "a" "a" "c" "c"
    [1063] "c" "a" "a" "t" "a" "g" "a" "a" "t" "c" "a" "c" "t" "c" "g" "a" "a" "a"
    [1081] "a" "a" "g" "a" "a" "t" "c" "t" "g" "c" "t" "t" "t" "c" "a" "a" "a" "a"
    [1099] "c" "g" "a" "a" "a" "g" "c" "t" "g" "a" "a" "c" "c" "t" "a" "t" "a" "a"
    [1117] "g" "c" "a" "g" "c" "a" "g" "t" "a" "t" "a" "a" "g" "c" "a" "a" "t" "a"
    [1135] "t" "g" "g" "a" "a" "c" "t" "c" "g" "a" "a" "t" "t" "a" "a" "a" "t" "a"
    [1153] "t" "c" "c" "a" "c" "a" "a" "t" "t" "c" "a" "a" "a" "a" "g" "c" "a" "c"
    [1171] "c" "t" "a" "a" "a" "a" "a" "g" "a" "a" "t" "a" "g" "g" "c" "t" "g" "a"
    [1189] "g" "g" "a" "g" "g" "a" "a" "g" "t" "c" "t" "t" "c" "t" "a" "c" "c" "a"
    [1207] "g" "g" "c" "a" "t" "a" "t" "t" "c" "a" "t" "g" "c" "g" "c" "t" "t" "g"
    [1225] "a" "a" "c" "t" "a" "g" "t" "a" "g" "t" "c" "a" "g" "t" "a" "g" "a" "a"
    [1243] "a" "t" "c" "t" "a" "a" "g" "c" "c" "c" "a" "c" "c" "t" "a" "a" "t" "t"
    [1261] "g" "t" "a" "c" "t" "g" "a" "a" "t" "t" "g" "c" "a" "a" "a" "t" "t" "g"
    [1279] "a" "t" "a" "g" "t" "t" "g" "t" "t" "c" "t" "a" "g" "c" "a" "a" "a" "a"
    [1297] "a" "a" "a" "a" "a" "a" "a" "a" "a" "a" "a" "a" "a" "a" "a" "a" "a" "a"
    [1315] "a" "a" "a" "a" "a" "a" "a" "a"
    attr(,"name")
    [1] "BC062429"
    attr(,"Annot")
    [1] ">BC062429"
    attr(,"class")
    [1] "SeqFastadna"
    
    $BC072418.BRCA1
       [1] "a" "t" "g" "g" "a" "t" "t" "t" "a" "t" "c" "t" "g" "c" "t" "c" "t" "t"
      [19] "c" "g" "c" "g" "t" "t" "g" "a" "a" "g" "a" "a" "g" "t" "a" "c" "a" "a"
      [37] "a" "a" "t" "g" "t" "c" "a" "t" "t" "a" "a" "t" "g" "c" "t" "a" "t" "g"
      [55] "c" "a" "g" "a" "a" "a" "a" "t" "c" "t" "t" "a" "g" "a" "g" "t" "g" "t"
      [73] "c" "c" "c" "a" "t" "c" "t" "g" "t" "c" "t" "g" "g" "a" "g" "t" "t" "g"
      [91] "a" "t" "c" "a" "a" "g" "g" "a" "a" "c" "c" "t" "g" "t" "c" "t" "c" "c"
     [109] "a" "c" "a" "a" "a" "g" "t" "g" "t" "g" "a" "c" "c" "a" "c" "a" "t" "a"
     [127] "t" "t" "t" "t" "g" "c" "a" "a" "a" "t" "t" "t" "t" "g" "c" "a" "t" "g"
     [145] "c" "t" "g" "a" "a" "a" "c" "t" "t" "c" "t" "c" "a" "a" "c" "c" "a" "g"
     [163] "a" "a" "g" "a" "a" "a" "g" "g" "g" "c" "c" "t" "t" "c" "a" "c" "a" "g"
     [181] "t" "g" "t" "c" "c" "t" "t" "t" "a" "t" "g" "t" "a" "a" "g" "a" "a" "t"
     [199] "g" "a" "t" "a" "t" "a" "a" "c" "c" "a" "a" "a" "a" "g" "g" "a" "g" "c"
     [217] "c" "t" "a" "c" "a" "a" "g" "a" "a" "a" "g" "t" "a" "c" "g" "a" "g" "a"
     [235] "t" "t" "t" "a" "g" "t" "c" "a" "a" "c" "t" "t" "g" "t" "t" "g" "a" "a"
     [253] "g" "a" "g" "c" "t" "a" "t" "t" "g" "a" "a" "a" "a" "t" "c" "a" "t" "t"
     [271] "t" "g" "t" "g" "c" "t" "t" "t" "t" "c" "a" "g" "c" "t" "t" "g" "a" "c"
     [289] "a" "c" "a" "g" "g" "t" "t" "t" "g" "g" "a" "g" "t" "a" "t" "g" "c" "a"
     [307] "a" "a" "c" "a" "g" "c" "t" "a" "t" "a" "a" "t" "t" "t" "t" "g" "c" "a"
     [325] "a" "a" "a" "a" "a" "g" "g" "a" "a" "a" "a" "t" "a" "a" "c" "t" "c" "t"
     [343] "c" "c" "t" "g" "a" "a" "c" "a" "t" "c" "t" "a" "a" "a" "a" "g" "a" "t"
     [361] "g" "a" "a" "g" "t" "t" "t" "c" "t" "a" "t" "c" "a" "t" "c" "c" "a" "a"
     [379] "a" "g" "t" "a" "t" "g" "g" "g" "c" "t" "a" "c" "a" "g" "a" "a" "a" "c"
     [397] "c" "g" "t" "g" "c" "c" "a" "a" "a" "a" "g" "a" "c" "t" "t" "c" "t" "a"
     [415] "c" "a" "g" "a" "g" "t" "g" "a" "a" "c" "c" "c" "g" "a" "a" "a" "a" "t"
     [433] "c" "c" "t" "t" "c" "c" "t" "t" "g" "c" "a" "g" "g" "a" "a" "a" "c" "c"
     [451] "a" "g" "t" "c" "t" "c" "a" "g" "t" "g" "t" "c" "c" "a" "a" "c" "t" "c"
     [469] "t" "c" "t" "a" "a" "c" "c" "t" "t" "g" "g" "a" "a" "c" "t" "g" "t" "g"
     [487] "a" "g" "a" "a" "c" "t" "c" "t" "g" "a" "g" "g" "a" "c" "a" "a" "a" "g"
     [505] "c" "a" "g" "c" "g" "g" "a" "t" "a" "c" "a" "a" "c" "c" "t" "c" "a" "a"
     [523] "a" "a" "g" "a" "c" "g" "t" "c" "t" "g" "t" "c" "t" "a" "c" "a" "t" "t"
     [541] "g" "a" "a" "t" "t" "g" "g" "g" "a" "t" "c" "t" "g" "a" "t" "t" "c" "t"
     [559] "t" "c" "t" "g" "a" "a" "g" "a" "t" "a" "c" "c" "g" "t" "t" "a" "a" "t"
     [577] "a" "a" "g" "g" "c" "a" "a" "c" "t" "t" "a" "t" "t" "g" "c" "a" "g" "t"
     [595] "g" "t" "g" "g" "g" "a" "g" "a" "t" "c" "a" "a" "g" "a" "a" "t" "t" "g"
     [613] "t" "t" "a" "c" "a" "a" "a" "t" "c" "a" "c" "c" "c" "c" "t" "c" "a" "a"
     [631] "g" "g" "a" "a" "c" "c" "a" "g" "g" "g" "a" "t" "g" "a" "a" "a" "t" "c"
     [649] "a" "g" "t" "t" "t" "g" "g" "a" "t" "t" "c" "t" "g" "c" "a" "a" "a" "a"
     [667] "a" "a" "g" "g" "c" "t" "g" "c" "t" "t" "g" "t" "g" "a" "a" "t" "t" "t"
     [685] "t" "c" "t" "g" "a" "g" "a" "c" "g" "g" "a" "t" "g" "t" "a" "a" "c" "a"
     [703] "a" "a" "t" "a" "c" "t" "g" "a" "a" "c" "a" "t" "c" "a" "t" "c" "a" "a"
     [721] "c" "c" "c" "a" "g" "t" "a" "a" "t" "a" "a" "t" "g" "a" "t" "t" "t" "g"
     [739] "a" "a" "c" "a" "c" "c" "a" "c" "t" "g" "a" "g" "a" "a" "g" "c" "g" "t"
     [757] "g" "c" "a" "g" "c" "t" "g" "a" "g" "a" "g" "g" "c" "a" "t" "c" "c" "a"
     [775] "g" "a" "a" "a" "a" "g" "t" "a" "t" "c" "a" "g" "g" "g" "t" "g" "a" "a"
     [793] "g" "c" "a" "g" "c" "a" "t" "c" "t" "g" "g" "g" "t" "g" "t" "g" "a" "g"
     [811] "a" "g" "t" "g" "a" "a" "a" "c" "a" "a" "g" "c" "g" "t" "c" "t" "c" "t"
     [829] "g" "a" "a" "g" "a" "c" "t" "g" "c" "t" "c" "a" "g" "g" "g" "c" "t" "a"
     [847] "t" "c" "c" "t" "c" "t" "c" "a" "g" "a" "g" "t" "g" "a" "c" "a" "t" "t"
     [865] "t" "t" "a" "a" "c" "c" "a" "c" "t" "c" "a" "g" "c" "a" "g" "a" "g" "g"
     [883] "g" "a" "t" "a" "c" "c" "a" "t" "g" "c" "a" "a" "c" "a" "t" "a" "a" "c"
     [901] "c" "t" "g" "a" "t" "a" "a" "a" "g" "c" "t" "c" "c" "a" "g" "c" "a" "g"
     [919] "g" "a" "a" "a" "t" "g" "g" "c" "t" "g" "a" "a" "c" "t" "a" "g" "a" "a"
     [937] "g" "c" "t" "g" "t" "g" "t" "t" "a" "g" "a" "a" "c" "a" "g" "c" "a" "t"
     [955] "g" "g" "g" "a" "g" "c" "c" "a" "g" "c" "c" "t" "t" "c" "t" "a" "a" "c"
     [973] "a" "g" "c" "t" "a" "c" "c" "c" "t" "t" "c" "c" "a" "t" "c" "a" "t" "a"
     [991] "a" "g" "t" "g" "a" "c" "t" "c" "t" "t" "c" "t" "g" "c" "c" "c" "t" "t"
    [1009] "g" "a" "g" "g" "a" "c" "c" "t" "g" "c" "g" "a" "a" "a" "t" "c" "c" "a"
    [1027] "g" "a" "a" "c" "a" "a" "a" "g" "c" "a" "c" "a" "t" "c" "a" "g" "a" "a"
    [1045] "a" "a" "a" "g" "t" "a" "t" "t" "a" "a" "c" "t" "t" "c" "a" "c" "a" "g"
    [1063] "a" "a" "a" "a" "g" "t" "a" "g" "t" "g" "a" "a" "t" "a" "c" "c" "c" "t"
    [1081] "a" "t" "a" "a" "g" "c" "c" "a" "g" "a" "a" "t" "c" "c" "a" "g" "a" "a"
    [1099] "g" "g" "c" "c" "t" "t" "t" "c" "t" "g" "c" "t" "g" "a" "c" "a" "a" "g"
    [1117] "t" "t" "t" "g" "a" "g" "g" "t" "g" "t" "c" "t" "g" "c" "a" "g" "a" "t"
    [1135] "a" "g" "t" "t" "c" "t" "a" "c" "c" "a" "g" "t" "a" "a" "a" "a" "a" "t"
    [1153] "a" "a" "a" "g" "a" "a" "c" "c" "a" "g" "g" "a" "g" "t" "g" "g" "a" "a"
    [1171] "a" "g" "g" "t" "c" "a" "t" "c" "c" "c" "c" "t" "t" "c" "t" "a" "a" "a"
    [1189] "t" "g" "c" "c" "c" "a" "t" "c" "a" "t" "t" "a" "g" "a" "t" "g" "a" "t"
    [1207] "a" "g" "g" "t" "g" "g" "t" "a" "c" "a" "t" "g" "c" "a" "c" "a" "g" "t"
    [1225] "t" "g" "c" "t" "c" "t" "g" "g" "g" "a" "g" "t" "c" "t" "t" "c" "a" "g"
    [1243] "a" "a" "t" "a" "g" "a" "a" "a" "c" "t" "a" "c" "c" "c" "a" "t" "c" "t"
    [1261] "c" "a" "a" "g" "a" "g" "g" "a" "g" "c" "t" "c" "a" "t" "t" "a" "a" "g"
    [1279] "g" "t" "t" "g" "t" "t" "g" "a" "t" "g" "t" "g" "g" "a" "g" "g" "a" "g"
    [1297] "c" "a" "a" "c" "a" "g" "c" "t" "g" "g" "a" "a" "g" "a" "g" "t" "c" "t"
    [1315] "g" "g" "g" "c" "c" "a" "c" "a" "c" "g" "a" "t" "t" "t" "g" "a" "c" "g"
    [1333] "g" "a" "a" "a" "c" "a" "t" "c" "t" "t" "a" "c" "t" "t" "g" "c" "c" "a"
    [1351] "a" "g" "g" "c" "a" "a" "g" "a" "t" "c" "t" "a" "g" "a" "g" "g" "g" "a"
    [1369] "a" "c" "c" "c" "c" "t" "t" "a" "c" "c" "t" "g" "g" "a" "a" "t" "c" "t"
    [1387] "g" "g" "a" "a" "t" "c" "a" "g" "c" "c" "t" "c" "t" "t" "c" "t" "c" "t"
    [1405] "g" "a" "t" "g" "a" "c" "c" "c" "t" "g" "a" "a" "t" "c" "t" "g" "a" "t"
    [1423] "c" "c" "t" "t" "c" "t" "g" "a" "a" "g" "a" "c" "a" "g" "a" "g" "c" "c"
    [1441] "c" "c" "a" "g" "a" "g" "t" "c" "a" "g" "c" "t" "c" "g" "t" "g" "t" "t"
    [1459] "g" "g" "c" "a" "a" "c" "a" "t" "a" "c" "c" "a" "t" "c" "t" "t" "c" "a"
    [1477] "a" "c" "c" "t" "c" "t" "g" "c" "a" "t" "t" "g" "a" "a" "a" "g" "t" "t"
    [1495] "c" "c" "c" "c" "a" "a" "t" "t" "g" "a" "a" "a" "g" "t" "t" "g" "c" "a"
    [1513] "g" "a" "a" "t" "c" "t" "g" "c" "c" "c" "a" "g" "a" "g" "t" "c" "c" "a"
    [1531] "g" "c" "t" "g" "c" "t" "g" "c" "t" "c" "a" "t" "a" "c" "t" "a" "c" "t"
    [1549] "g" "a" "t" "a" "c" "t" "g" "c" "t" "g" "g" "g" "t" "a" "t" "a" "a" "t"
    [1567] "g" "c" "a" "a" "t" "g" "g" "a" "a" "g" "a" "a" "a" "g" "t" "g" "t" "g"
    [1585] "a" "g" "c" "a" "g" "g" "g" "a" "g" "a" "a" "g" "c" "c" "a" "g" "a" "a"
    [1603] "t" "t" "g" "a" "c" "a" "g" "c" "t" "t" "c" "a" "a" "c" "a" "g" "a" "a"
    [1621] "a" "g" "g" "g" "t" "c" "a" "a" "c" "a" "a" "a" "a" "g" "a" "a" "t" "g"
    [1639] "t" "c" "c" "a" "t" "g" "g" "t" "g" "g" "t" "g" "t" "c" "t" "g" "g" "c"
    [1657] "c" "t" "g" "a" "c" "c" "c" "c" "a" "g" "a" "a" "g" "a" "a" "t" "t" "t"
    [1675] "a" "t" "g" "c" "t" "c" "g" "t" "g" "t" "a" "c" "a" "a" "g" "t" "t" "t"
    [1693] "g" "c" "c" "a" "g" "a" "a" "a" "a" "c" "a" "c" "c" "a" "c" "a" "t" "c"
    [1711] "a" "c" "t" "t" "t" "a" "a" "c" "t" "a" "a" "t" "c" "t" "a" "a" "t" "t"
    [1729] "a" "c" "t" "g" "a" "a" "g" "a" "g" "a" "c" "t" "a" "c" "t" "c" "a" "t"
    [1747] "g" "t" "t" "g" "t" "t" "a" "t" "g" "a" "a" "a" "a" "c" "a" "g" "a" "t"
    [1765] "g" "c" "t" "g" "a" "g" "t" "t" "t" "g" "t" "g" "t" "g" "t" "g" "a" "a"
    [1783] "c" "g" "g" "a" "c" "a" "c" "t" "g" "a" "a" "a" "t" "a" "t" "t" "t" "t"
    [1801] "c" "t" "a" "g" "g" "a" "a" "t" "t" "g" "c" "g" "g" "g" "a" "g" "g" "a"
    [1819] "a" "a" "a" "t" "g" "g" "g" "t" "a" "g" "t" "t" "a" "g" "c" "t" "a" "t"
    [1837] "t" "t" "c" "t" "g" "g" "g" "t" "g" "a" "c" "c" "c" "a" "g" "t" "c" "t"
    [1855] "a" "t" "t" "a" "a" "a" "g" "a" "a" "a" "g" "a" "a" "a" "a" "a" "t" "g"
    [1873] "c" "t" "g" "a" "a" "t" "g" "a" "g" "c" "a" "t" "g" "a" "t" "t" "t" "t"
    [1891] "g" "a" "a" "g" "t" "c" "a" "g" "a" "g" "g" "a" "g" "a" "t" "g" "t" "g"
    [1909] "g" "t" "c" "a" "a" "t" "g" "g" "a" "a" "g" "a" "a" "a" "c" "c" "a" "c"
    [1927] "c" "a" "a" "g" "g" "t" "c" "c" "a" "a" "a" "g" "c" "g" "a" "g" "c" "a"
    [1945] "a" "g" "a" "g" "a" "a" "t" "c" "c" "c" "a" "g" "g" "a" "c" "a" "g" "a"
    [1963] "a" "a" "g" "a" "t" "c" "t" "t" "c" "a" "g" "g" "g" "g" "g" "c" "t" "a"
    [1981] "g" "a" "a" "a" "t" "c" "t" "g" "t" "t" "g" "c" "t" "a" "t" "g" "g" "g"
    [1999] "c" "c" "c" "t" "t" "c" "a" "c" "c" "a" "a" "c" "a" "t" "g" "c" "c" "c"
    [2017] "a" "c" "a" "g" "g" "g" "t" "g" "t" "c" "c" "a" "c" "c" "c" "a" "a" "t"
    [2035] "t" "g" "t" "g" "g" "t" "t" "g" "t" "g" "c" "a" "g" "c" "c" "a" "g" "a"
    [2053] "t" "g" "c" "c" "t" "g" "g" "a" "c" "a" "g" "a" "g" "g" "a" "c" "a" "a"
    [2071] "t" "g" "g" "c" "t" "t" "c" "c" "a" "t" "g" "c" "a" "a" "t" "t" "g" "g"
    [2089] "g" "c" "a" "g" "a" "t" "g" "t" "g" "t" "g" "a"
    attr(,"name")
    [1] "BC072418.BRCA1"
    attr(,"Annot")
    [1] ">BC072418.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $BC085615.BRCA1
       [1] "a" "t" "g" "g" "a" "t" "t" "t" "a" "t" "c" "t" "g" "c" "t" "c" "t" "t"
      [19] "c" "g" "c" "g" "t" "t" "g" "a" "a" "g" "a" "a" "g" "t" "a" "c" "a" "a"
      [37] "a" "a" "t" "g" "t" "c" "a" "t" "t" "a" "a" "t" "g" "c" "t" "a" "t" "g"
      [55] "c" "a" "g" "a" "a" "a" "a" "t" "c" "t" "t" "a" "g" "a" "g" "t" "g" "t"
      [73] "c" "c" "c" "a" "t" "c" "t" "g" "t" "c" "t" "g" "g" "a" "g" "t" "t" "g"
      [91] "a" "t" "c" "a" "a" "g" "g" "a" "a" "c" "c" "t" "g" "t" "c" "t" "c" "c"
     [109] "a" "c" "a" "a" "a" "g" "t" "g" "t" "g" "a" "c" "c" "a" "c" "a" "t" "a"
     [127] "t" "t" "t" "t" "g" "c" "a" "a" "a" "t" "t" "t" "t" "g" "c" "a" "t" "g"
     [145] "c" "t" "g" "a" "a" "a" "c" "t" "t" "c" "t" "c" "a" "a" "c" "c" "a" "g"
     [163] "a" "a" "g" "a" "a" "a" "g" "g" "g" "c" "c" "t" "t" "c" "a" "c" "a" "g"
     [181] "t" "g" "t" "c" "c" "t" "t" "t" "a" "t" "g" "t" "a" "a" "g" "a" "a" "t"
     [199] "g" "a" "t" "a" "t" "a" "a" "c" "c" "a" "a" "a" "a" "g" "g" "a" "g" "c"
     [217] "c" "t" "a" "c" "a" "a" "g" "a" "a" "a" "g" "t" "a" "c" "g" "a" "g" "a"
     [235] "t" "t" "t" "a" "g" "t" "c" "a" "a" "c" "t" "t" "g" "t" "t" "g" "a" "a"
     [253] "g" "a" "g" "c" "t" "a" "t" "t" "g" "a" "a" "a" "a" "t" "c" "a" "t" "t"
     [271] "t" "g" "t" "g" "c" "t" "t" "t" "t" "c" "a" "g" "c" "t" "t" "g" "a" "c"
     [289] "a" "c" "a" "g" "g" "t" "t" "t" "g" "g" "a" "g" "t" "a" "t" "g" "c" "a"
     [307] "a" "a" "c" "a" "g" "c" "t" "a" "t" "a" "a" "t" "t" "t" "t" "g" "c" "a"
     [325] "a" "a" "a" "a" "a" "g" "g" "a" "a" "a" "a" "t" "a" "a" "c" "t" "c" "t"
     [343] "c" "c" "t" "g" "a" "a" "c" "a" "t" "c" "t" "a" "a" "a" "a" "g" "a" "t"
     [361] "g" "a" "a" "g" "t" "t" "t" "c" "t" "a" "t" "c" "a" "t" "c" "c" "a" "a"
     [379] "a" "g" "t" "a" "t" "g" "g" "g" "c" "t" "a" "c" "a" "g" "a" "a" "a" "c"
     [397] "c" "g" "t" "g" "c" "c" "a" "a" "a" "a" "g" "a" "c" "t" "t" "c" "t" "a"
     [415] "c" "a" "g" "a" "g" "t" "g" "a" "a" "c" "c" "c" "g" "a" "a" "a" "a" "t"
     [433] "c" "c" "t" "t" "c" "c" "t" "t" "g" "g" "a" "a" "a" "c" "c" "a" "g" "t"
     [451] "c" "t" "c" "a" "g" "t" "g" "t" "c" "c" "a" "a" "c" "t" "c" "t" "c" "t"
     [469] "a" "a" "c" "c" "t" "t" "g" "g" "a" "a" "c" "t" "g" "t" "g" "a" "g" "a"
     [487] "a" "c" "t" "c" "t" "g" "a" "g" "g" "a" "c" "a" "a" "a" "g" "c" "a" "g"
     [505] "c" "g" "g" "a" "t" "a" "c" "a" "a" "c" "c" "t" "c" "a" "a" "a" "a" "g"
     [523] "a" "c" "g" "t" "c" "t" "g" "t" "c" "t" "a" "c" "a" "t" "t" "g" "a" "a"
     [541] "t" "t" "g" "g" "g" "a" "t" "c" "t" "g" "a" "t" "t" "c" "t" "t" "c" "t"
     [559] "g" "a" "a" "g" "a" "t" "a" "c" "c" "g" "t" "t" "a" "a" "t" "a" "a" "g"
     [577] "g" "c" "a" "a" "c" "t" "t" "a" "t" "t" "g" "c" "a" "g" "t" "g" "t" "g"
     [595] "g" "g" "a" "g" "a" "t" "c" "a" "a" "g" "a" "a" "t" "t" "g" "t" "t" "a"
     [613] "c" "a" "a" "a" "t" "c" "a" "c" "c" "c" "c" "t" "c" "a" "a" "g" "g" "a"
     [631] "a" "c" "c" "a" "g" "g" "g" "a" "t" "g" "a" "a" "a" "t" "c" "a" "g" "t"
     [649] "t" "t" "g" "g" "a" "t" "t" "c" "t" "g" "c" "a" "a" "a" "a" "a" "a" "g"
     [667] "g" "c" "t" "g" "c" "t" "t" "g" "t" "g" "a" "a" "t" "t" "t" "t" "c" "t"
     [685] "g" "a" "g" "a" "c" "g" "g" "a" "t" "g" "t" "a" "a" "c" "a" "a" "a" "t"
     [703] "a" "c" "t" "g" "a" "a" "c" "a" "t" "c" "a" "t" "c" "a" "a" "c" "c" "c"
     [721] "a" "g" "t" "a" "a" "t" "a" "a" "t" "g" "a" "t" "t" "t" "g" "a" "a" "c"
     [739] "a" "c" "c" "a" "c" "t" "g" "a" "g" "a" "a" "g" "c" "g" "t" "g" "c" "a"
     [757] "g" "c" "t" "g" "a" "g" "a" "g" "g" "c" "a" "t" "c" "c" "a" "g" "a" "a"
     [775] "a" "a" "g" "t" "a" "t" "c" "a" "g" "g" "g" "t" "g" "a" "a" "g" "c" "a"
     [793] "g" "c" "a" "t" "c" "t" "g" "g" "g" "t" "g" "t" "g" "a" "g" "a" "g" "t"
     [811] "g" "a" "a" "a" "c" "a" "a" "g" "c" "g" "t" "c" "t" "c" "t" "g" "a" "a"
     [829] "g" "a" "c" "t" "g" "c" "t" "c" "a" "g" "g" "g" "c" "t" "a" "t" "c" "c"
     [847] "t" "c" "t" "c" "a" "g" "a" "g" "t" "g" "a" "c" "a" "t" "t" "t" "t" "a"
     [865] "a" "c" "c" "a" "c" "t" "c" "a" "g" "c" "a" "g" "a" "g" "g" "g" "a" "t"
     [883] "a" "c" "c" "a" "t" "g" "c" "a" "a" "c" "a" "t" "a" "a" "c" "c" "t" "g"
     [901] "a" "t" "a" "a" "a" "g" "c" "t" "c" "c" "a" "g" "c" "a" "g" "g" "a" "a"
     [919] "a" "t" "g" "g" "c" "t" "g" "a" "a" "c" "t" "a" "g" "a" "a" "g" "c" "t"
     [937] "g" "t" "g" "t" "t" "a" "g" "a" "a" "c" "a" "g" "c" "a" "t" "g" "g" "g"
     [955] "a" "g" "c" "c" "a" "g" "c" "c" "t" "t" "c" "t" "a" "a" "c" "a" "g" "c"
     [973] "t" "a" "c" "c" "c" "t" "t" "c" "c" "a" "t" "c" "a" "t" "a" "a" "g" "t"
     [991] "g" "a" "c" "t" "c" "t" "t" "c" "t" "g" "c" "c" "c" "t" "t" "g" "a" "g"
    [1009] "g" "a" "c" "c" "t" "g" "c" "g" "a" "a" "a" "t" "c" "c" "a" "g" "a" "a"
    [1027] "c" "a" "a" "a" "g" "c" "a" "c" "a" "t" "c" "a" "g" "a" "a" "a" "a" "a"
    [1045] "g" "c" "a" "g" "t" "a" "t" "t" "a" "a" "c" "t" "t" "c" "a" "c" "a" "g"
    [1063] "a" "a" "a" "a" "g" "t" "a" "g" "t" "g" "a" "a" "t" "a" "c" "c" "c" "t"
    [1081] "a" "t" "a" "a" "g" "c" "c" "a" "g" "a" "a" "t" "c" "c" "a" "g" "a" "a"
    [1099] "g" "g" "c" "c" "t" "t" "t" "c" "t" "g" "c" "t" "g" "a" "c" "a" "a" "g"
    [1117] "t" "t" "t" "g" "a" "g" "g" "t" "g" "t" "c" "t" "g" "c" "a" "g" "a" "t"
    [1135] "a" "g" "t" "t" "c" "t" "a" "c" "c" "a" "g" "t" "a" "a" "a" "a" "a" "t"
    [1153] "a" "a" "a" "g" "a" "a" "c" "c" "a" "g" "g" "a" "g" "t" "g" "g" "a" "a"
    [1171] "a" "g" "g" "t" "c" "a" "t" "c" "c" "c" "c" "t" "t" "c" "t" "a" "a" "a"
    [1189] "t" "g" "c" "c" "c" "a" "t" "c" "a" "t" "t" "a" "g" "a" "t" "g" "a" "t"
    [1207] "a" "g" "g" "t" "g" "g" "t" "a" "c" "a" "t" "g" "c" "a" "c" "a" "g" "t"
    [1225] "t" "g" "c" "t" "c" "t" "g" "g" "g" "a" "g" "t" "c" "t" "t" "c" "a" "g"
    [1243] "a" "a" "t" "a" "g" "a" "a" "a" "c" "t" "a" "c" "c" "c" "a" "t" "c" "t"
    [1261] "c" "a" "a" "g" "a" "g" "g" "a" "g" "c" "t" "c" "a" "t" "t" "a" "a" "g"
    [1279] "g" "t" "t" "g" "t" "t" "g" "a" "t" "g" "t" "g" "g" "a" "g" "g" "a" "g"
    [1297] "c" "a" "a" "c" "a" "g" "c" "t" "g" "g" "a" "a" "g" "a" "g" "t" "c" "t"
    [1315] "g" "g" "g" "c" "c" "a" "c" "a" "c" "g" "a" "t" "t" "t" "g" "a" "c" "g"
    [1333] "g" "a" "a" "a" "c" "a" "t" "c" "t" "t" "a" "c" "t" "t" "g" "c" "c" "a"
    [1351] "a" "g" "g" "c" "a" "a" "g" "a" "t" "c" "t" "a" "g" "a" "g" "g" "g" "a"
    [1369] "a" "c" "c" "c" "c" "t" "t" "a" "c" "c" "t" "g" "g" "a" "a" "t" "c" "t"
    [1387] "g" "g" "a" "a" "t" "c" "a" "g" "c" "c" "t" "c" "t" "t" "c" "t" "c" "t"
    [1405] "g" "a" "t" "g" "a" "c" "c" "c" "t" "g" "a" "a" "t" "c" "t" "g" "a" "t"
    [1423] "c" "c" "t" "t" "c" "t" "g" "a" "a" "g" "a" "c" "a" "g" "a" "g" "c" "c"
    [1441] "c" "c" "a" "g" "a" "g" "t" "c" "a" "g" "c" "t" "c" "g" "t" "g" "t" "t"
    [1459] "g" "g" "c" "a" "a" "c" "a" "t" "a" "c" "c" "a" "t" "c" "t" "t" "c" "a"
    [1477] "a" "c" "c" "t" "c" "t" "g" "c" "a" "t" "t" "g" "a" "a" "a" "g" "t" "t"
    [1495] "c" "c" "c" "c" "a" "a" "t" "t" "g" "a" "a" "a" "g" "t" "t" "g" "c" "a"
    [1513] "g" "a" "a" "t" "c" "t" "g" "c" "c" "c" "a" "g" "a" "g" "t" "c" "c" "a"
    [1531] "g" "c" "t" "g" "c" "t" "g" "c" "t" "c" "a" "t" "a" "c" "t" "a" "c" "t"
    [1549] "g" "a" "t" "a" "c" "t" "g" "c" "t" "g" "g" "g" "t" "a" "t" "a" "a" "t"
    [1567] "g" "c" "a" "a" "t" "g" "g" "a" "a" "g" "a" "a" "a" "g" "t" "g" "t" "g"
    [1585] "a" "g" "c" "a" "g" "g" "g" "a" "g" "a" "a" "g" "c" "c" "a" "g" "a" "a"
    [1603] "t" "t" "g" "a" "c" "a" "g" "c" "t" "t" "c" "a" "a" "c" "a" "g" "a" "a"
    [1621] "a" "g" "g" "g" "t" "c" "a" "a" "c" "a" "a" "a" "a" "g" "a" "a" "t" "g"
    [1639] "t" "c" "c" "a" "t" "g" "g" "t" "g" "g" "t" "g" "t" "c" "t" "g" "g" "c"
    [1657] "c" "t" "g" "a" "c" "c" "c" "c" "a" "g" "a" "a" "g" "a" "a" "t" "t" "t"
    [1675] "a" "t" "g" "c" "t" "c" "g" "t" "g" "t" "a" "c" "a" "a" "g" "t" "t" "t"
    [1693] "g" "c" "c" "a" "g" "a" "a" "a" "a" "c" "a" "c" "c" "a" "c" "a" "t" "c"
    [1711] "a" "c" "t" "t" "t" "a" "a" "c" "t" "a" "a" "t" "c" "t" "a" "a" "t" "t"
    [1729] "a" "c" "t" "g" "a" "a" "g" "a" "g" "a" "c" "t" "a" "c" "t" "c" "a" "t"
    [1747] "g" "t" "t" "g" "t" "t" "a" "t" "g" "a" "a" "a" "a" "c" "a" "g" "a" "t"
    [1765] "g" "c" "t" "g" "a" "g" "t" "t" "t" "g" "t" "g" "t" "g" "t" "g" "a" "a"
    [1783] "c" "g" "g" "a" "c" "a" "c" "t" "g" "a" "a" "a" "t" "a" "t" "t" "t" "t"
    [1801] "c" "t" "a" "g" "g" "a" "a" "t" "t" "g" "c" "g" "g" "g" "a" "g" "g" "a"
    [1819] "a" "a" "a" "t" "g" "g" "g" "t" "a" "g" "t" "t" "a" "g" "c" "t" "a" "t"
    [1837] "t" "t" "c" "t" "g" "g" "g" "t" "g" "a" "c" "c" "c" "a" "g" "t" "c" "t"
    [1855] "a" "t" "t" "a" "a" "a" "a" "a" "a" "a" "a" "a" "a" "a" "a" "a" "a" "a"
    attr(,"name")
    [1] "BC085615.BRCA1"
    attr(,"Annot")
    [1] ">BC085615.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $BC106745.BRCA1
       [1] "a" "t" "g" "g" "a" "t" "t" "t" "a" "t" "c" "t" "g" "c" "t" "c" "t" "t"
      [19] "c" "g" "c" "g" "t" "t" "g" "a" "a" "g" "a" "a" "g" "t" "a" "c" "a" "a"
      [37] "a" "a" "t" "g" "t" "c" "a" "t" "t" "a" "a" "t" "g" "c" "t" "a" "t" "g"
      [55] "c" "a" "g" "a" "a" "a" "a" "t" "c" "t" "t" "a" "g" "a" "g" "t" "g" "t"
      [73] "c" "c" "c" "a" "t" "c" "t" "g" "t" "c" "t" "g" "g" "a" "g" "t" "t" "g"
      [91] "a" "t" "c" "a" "a" "g" "g" "a" "a" "c" "c" "t" "g" "t" "c" "t" "c" "c"
     [109] "a" "c" "a" "a" "a" "g" "t" "g" "t" "g" "a" "c" "c" "a" "c" "a" "t" "a"
     [127] "t" "t" "t" "t" "g" "c" "a" "a" "a" "t" "t" "t" "t" "g" "c" "a" "t" "g"
     [145] "c" "t" "g" "a" "a" "a" "c" "t" "t" "c" "t" "c" "a" "a" "c" "c" "a" "g"
     [163] "a" "a" "g" "a" "a" "a" "g" "g" "g" "c" "c" "t" "t" "c" "a" "c" "a" "g"
     [181] "t" "g" "t" "c" "c" "t" "t" "t" "a" "t" "g" "t" "a" "a" "g" "a" "a" "t"
     [199] "g" "a" "t" "a" "t" "a" "a" "c" "c" "a" "a" "a" "a" "g" "g" "a" "g" "c"
     [217] "c" "t" "a" "c" "a" "a" "g" "a" "a" "a" "g" "t" "a" "c" "g" "a" "g" "a"
     [235] "t" "t" "t" "a" "g" "t" "c" "a" "a" "c" "t" "t" "g" "t" "t" "g" "a" "a"
     [253] "g" "a" "g" "c" "t" "a" "t" "t" "g" "a" "a" "a" "a" "t" "c" "a" "t" "t"
     [271] "t" "g" "t" "g" "c" "t" "t" "t" "t" "c" "a" "g" "c" "t" "t" "g" "a" "c"
     [289] "a" "c" "a" "g" "g" "t" "t" "t" "g" "g" "a" "g" "t" "a" "t" "g" "c" "a"
     [307] "a" "a" "c" "a" "g" "c" "t" "a" "t" "a" "a" "t" "t" "t" "t" "g" "c" "a"
     [325] "a" "a" "a" "a" "a" "g" "g" "a" "a" "a" "a" "t" "a" "a" "c" "t" "c" "t"
     [343] "c" "c" "t" "g" "a" "a" "c" "a" "t" "c" "t" "a" "a" "a" "a" "g" "a" "t"
     [361] "g" "a" "a" "g" "t" "t" "t" "c" "t" "a" "t" "c" "a" "t" "c" "c" "a" "a"
     [379] "a" "g" "t" "a" "t" "g" "g" "g" "c" "t" "a" "c" "a" "g" "a" "a" "a" "c"
     [397] "c" "g" "t" "g" "c" "c" "a" "a" "a" "a" "g" "a" "c" "t" "t" "c" "t" "a"
     [415] "c" "a" "g" "a" "g" "t" "g" "a" "a" "c" "c" "c" "g" "a" "a" "a" "a" "t"
     [433] "c" "c" "t" "t" "c" "c" "t" "t" "g" "c" "a" "g" "g" "a" "a" "a" "c" "c"
     [451] "a" "g" "t" "c" "t" "c" "a" "g" "t" "g" "t" "c" "c" "a" "a" "c" "t" "c"
     [469] "t" "c" "t" "a" "a" "c" "c" "t" "t" "g" "g" "a" "a" "c" "t" "g" "t" "g"
     [487] "a" "g" "a" "a" "c" "t" "c" "t" "g" "a" "g" "g" "a" "c" "a" "a" "a" "g"
     [505] "c" "a" "g" "c" "g" "g" "a" "t" "a" "c" "a" "a" "c" "c" "t" "c" "a" "a"
     [523] "a" "a" "g" "a" "c" "g" "t" "c" "t" "g" "t" "c" "t" "a" "c" "a" "t" "t"
     [541] "g" "a" "a" "t" "t" "g" "g" "g" "a" "t" "c" "t" "g" "a" "t" "t" "c" "t"
     [559] "t" "c" "t" "g" "a" "a" "g" "a" "t" "a" "c" "c" "g" "t" "t" "a" "a" "t"
     [577] "a" "a" "g" "g" "c" "a" "a" "c" "t" "t" "a" "t" "t" "g" "c" "a" "g" "t"
     [595] "g" "t" "g" "g" "g" "a" "g" "a" "t" "c" "a" "a" "g" "a" "a" "t" "t" "g"
     [613] "t" "t" "a" "c" "a" "a" "a" "t" "c" "a" "c" "c" "c" "c" "t" "c" "a" "a"
     [631] "g" "g" "a" "a" "c" "c" "a" "g" "g" "g" "a" "t" "g" "a" "a" "a" "t" "c"
     [649] "a" "g" "t" "t" "t" "g" "g" "a" "t" "t" "c" "t" "g" "c" "a" "a" "a" "a"
     [667] "a" "a" "g" "g" "c" "t" "g" "c" "t" "t" "g" "t" "g" "a" "a" "t" "t" "t"
     [685] "t" "c" "t" "g" "a" "g" "a" "c" "g" "g" "a" "t" "g" "t" "a" "a" "c" "a"
     [703] "a" "a" "t" "a" "c" "t" "g" "a" "a" "c" "a" "t" "c" "a" "t" "c" "a" "a"
     [721] "c" "c" "c" "a" "g" "t" "a" "a" "t" "a" "a" "t" "g" "a" "t" "t" "t" "g"
     [739] "a" "a" "c" "a" "c" "c" "a" "c" "t" "g" "a" "g" "a" "a" "g" "c" "g" "t"
     [757] "g" "c" "a" "g" "c" "t" "g" "a" "g" "a" "g" "g" "c" "a" "t" "c" "c" "a"
     [775] "g" "a" "a" "a" "a" "g" "t" "a" "t" "c" "a" "g" "g" "g" "t" "a" "g" "t"
     [793] "t" "c" "t" "g" "t" "t" "t" "c" "a" "a" "a" "c" "t" "t" "g" "c" "a" "t"
     [811] "g" "t" "g" "g" "a" "g" "c" "c" "a" "t" "g" "t" "g" "g" "c" "a" "c" "a"
     [829] "a" "a" "t" "a" "c" "t" "c" "a" "t" "g" "c" "c" "a" "g" "c" "t" "c" "a"
     [847] "t" "t" "a" "c" "a" "g" "c" "a" "t" "g" "a" "g" "a" "a" "c" "a" "g" "c"
     [865] "a" "g" "t" "t" "t" "a" "t" "t" "a" "c" "t" "c" "a" "c" "t" "a" "a" "a"
     [883] "g" "a" "c" "a" "g" "a" "a" "t" "g" "a" "a" "t" "g" "t" "a" "g" "a" "a"
     [901] "a" "a" "g" "g" "c" "t" "g" "a" "a" "t" "t" "c" "t" "g" "t" "a" "a" "t"
     [919] "a" "a" "a" "a" "g" "c" "a" "a" "a" "c" "a" "g" "c" "c" "t" "g" "g" "c"
     [937] "t" "t" "a" "g" "c" "a" "a" "g" "g" "a" "g" "c" "c" "a" "a" "c" "a" "t"
     [955] "a" "a" "c" "a" "g" "a" "t" "g" "g" "g" "c" "t" "g" "g" "a" "a" "g" "t"
     [973] "a" "a" "g" "g" "a" "a" "a" "c" "a" "t" "g" "t" "a" "a" "t" "g" "a" "t"
     [991] "a" "g" "g" "c" "g" "g" "a" "c" "t" "c" "c" "c" "a" "g" "c" "a" "c" "a"
    [1009] "g" "a" "a" "a" "a" "a" "a" "a" "g" "g" "t" "a" "g" "a" "t" "c" "t" "g"
    [1027] "a" "a" "t" "g" "c" "t" "g" "a" "t" "c" "c" "c" "c" "t" "g" "t" "g" "t"
    [1045] "g" "a" "g" "a" "g" "a" "a" "a" "a" "g" "a" "a" "t" "g" "g" "a" "a" "t"
    [1063] "a" "a" "g" "c" "a" "g" "a" "a" "a" "c" "t" "g" "c" "c" "a" "t" "g" "c"
    [1081] "t" "c" "a" "g" "a" "g" "a" "a" "t" "c" "c" "t" "a" "g" "a" "g" "a" "t"
    [1099] "a" "c" "t" "g" "a" "a" "g" "a" "t" "g" "t" "t" "c" "c" "t" "t" "g" "g"
    [1117] "a" "t" "a" "a" "c" "a" "c" "t" "a" "a" "a" "t" "a" "g" "c" "a" "g" "c"
    [1135] "a" "t" "t" "c" "a" "g" "a" "a" "a" "g" "t" "t" "a" "a" "t" "g" "a" "g"
    [1153] "t" "g" "g" "t" "t" "t" "t" "c" "c" "a" "g" "a" "a" "g" "t" "g" "a" "t"
    [1171] "g" "a" "a" "c" "t" "g" "t" "t" "a" "g" "g" "t" "t" "c" "t" "g" "a" "t"
    [1189] "g" "a" "c" "t" "c" "a" "c" "a" "t" "g" "a" "t" "g" "g" "g" "g" "a" "g"
    [1207] "t" "c" "t" "g" "a" "a" "t" "c" "a" "a" "a" "t" "g" "c" "c" "a" "a" "a"
    [1225] "g" "t" "a" "g" "c" "t" "g" "a" "t" "g" "t" "a" "t" "t" "g" "g" "a" "c"
    [1243] "g" "t" "t" "c" "t" "a" "a" "a" "t" "g" "a" "g" "g" "t" "a" "g" "a" "t"
    [1261] "g" "a" "a" "t" "a" "t" "t" "c" "t" "g" "g" "t" "t" "c" "t" "t" "c" "a"
    [1279] "g" "a" "g" "a" "a" "a" "a" "t" "a" "g" "a" "c" "t" "t" "a" "c" "t" "g"
    [1297] "g" "c" "c" "a" "g" "t" "g" "a" "t" "c" "c" "t" "c" "a" "t" "g" "a" "g"
    [1315] "g" "c" "t" "t" "t" "a" "a" "t" "a" "t" "g" "t" "a" "a" "a" "a" "g" "t"
    [1333] "g" "a" "a" "a" "g" "a" "g" "t" "t" "c" "a" "c" "t" "c" "c" "a" "a" "a"
    [1351] "t" "c" "a" "g" "t" "a" "g" "a" "g" "a" "g" "t" "a" "a" "t" "a" "t" "t"
    [1369] "g" "a" "a" "g" "a" "c" "a" "a" "a" "a" "t" "a" "t" "t" "t" "g" "g" "g"
    [1387] "a" "a" "a" "a" "c" "c" "t" "a" "t" "c" "g" "g" "a" "a" "g" "a" "a" "g"
    [1405] "g" "c" "a" "a" "g" "c" "c" "t" "c" "c" "c" "c" "a" "a" "c"
    attr(,"name")
    [1] "BC106745.BRCA1"
    attr(,"Annot")
    [1] ">BC106745.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $BC106746
      [1] "g" "g" "a" "g" "g" "c" "c" "t" "t" "c" "a" "c" "c" "c" "t" "c" "t" "g"
     [19] "c" "t" "c" "t" "g" "g" "g" "t" "a" "a" "a" "g" "c" "t" "g" "c" "t" "t"
     [37] "g" "t" "g" "a" "a" "t" "t" "t" "t" "c" "t" "g" "a" "g" "a" "c" "g" "g"
     [55] "a" "t" "g" "t" "a" "a" "c" "a" "a" "a" "t" "a" "c" "t" "g" "a" "a" "c"
     [73] "a" "t" "c" "a" "t" "c" "a" "a" "c" "c" "c" "a" "g" "t" "a" "a" "t" "a"
     [91] "a" "t" "g" "a" "t" "t" "t" "g" "a" "a" "c" "a" "c" "c" "a" "c" "t" "g"
    [109] "a" "g" "a" "a" "g" "c" "g" "t" "g" "c" "a" "g" "c" "t" "g" "a" "g" "a"
    [127] "g" "g" "c" "a" "t" "c" "c" "a" "g" "a" "a" "a" "a" "g" "t" "a" "t" "c"
    [145] "a" "g" "g" "g" "t" "a" "g" "t" "t" "c" "t" "g" "t" "t" "t" "c" "a" "a"
    [163] "a" "c" "t" "t" "g" "c" "a" "t" "g" "t" "g" "g" "a" "g" "c" "c" "a" "t"
    [181] "g" "t" "g" "g" "c" "a" "c" "a" "a" "a" "t" "a" "c" "t" "c" "a" "t" "g"
    [199] "c" "c" "a" "g" "c" "t" "c" "a" "t" "t" "a" "c" "a" "g" "c" "a" "t" "g"
    [217] "a" "g" "a" "a" "c" "a" "g" "c" "a" "g" "t" "t" "t" "a" "t" "t" "a" "c"
    [235] "t" "c" "a" "c" "t" "a" "a" "a" "g" "a" "c" "a" "g" "a" "a" "t" "g" "a"
    [253] "a" "t" "g" "t" "a" "g" "a" "a" "a" "a" "g" "g" "c" "t" "g" "a" "a" "t"
    [271] "t" "c" "t" "g" "t" "a" "a" "t" "a" "a" "a" "a" "g" "c" "a" "a" "a" "c"
    [289] "a" "g" "c" "c" "t" "g" "g" "c" "t" "t" "a" "g" "c" "a" "a" "g" "g" "a"
    [307] "g" "c" "c" "a" "a" "c" "a" "t" "a" "a" "c" "a" "g" "a" "t" "g" "g" "g"
    [325] "c" "t" "g" "g" "a" "a" "g" "t" "a" "a" "g" "g" "a" "a" "a" "c" "a" "t"
    [343] "g" "t" "a" "a" "t" "g" "a" "t" "a" "g" "g" "c" "g" "g" "a" "c" "t" "c"
    [361] "c" "c" "a" "g" "c" "a" "c" "a" "g" "a" "a" "a" "a" "a" "a" "a" "g" "g"
    [379] "t" "a" "g" "a" "t" "c" "t" "g" "a" "a" "t" "g" "c" "t" "g" "a" "t" "c"
    [397] "c" "c" "c" "t" "g" "t" "g" "t" "g" "a" "g" "a" "g" "a" "a" "a" "a" "g"
    [415] "a" "a" "t" "g" "g" "a" "a" "t" "a" "a" "g" "c" "a" "g" "a" "a" "a" "c"
    [433] "t" "g" "c" "c" "a" "t" "g" "c" "t" "c" "a" "g" "a" "g" "a" "a" "t" "c"
    [451] "c" "t" "a" "g" "a" "g" "a" "t" "a" "c" "t" "g" "a" "a" "g" "a" "t" "g"
    [469] "t" "t" "c" "c" "t" "t" "g" "g" "a" "t" "a" "a" "c" "a" "c" "t" "a" "a"
    [487] "a" "t" "a" "g" "c" "a" "g" "c" "a" "t" "t" "c" "a" "g" "a" "a" "a" "g"
    [505] "t" "t" "a" "a" "t" "g" "a" "g" "t" "g" "g" "t" "t" "t" "t" "c" "c" "a"
    [523] "g" "a" "a" "g" "t" "g" "a" "t" "g" "a" "a" "c" "t" "g" "t" "t" "a" "g"
    [541] "g" "t" "t" "c" "t" "g" "a" "t" "g" "a" "c" "t" "c" "a" "c" "a" "t" "g"
    [559] "a" "t" "g" "g" "g" "g" "a" "g" "t" "c" "t" "g" "a" "a" "t" "c" "a" "a"
    [577] "a" "t" "g" "c" "c" "a" "a" "a" "g" "t" "a" "g" "c" "t" "g" "a" "t" "g"
    [595] "t" "a" "t" "t" "g" "g" "a" "c" "g" "t" "t" "c" "t" "a" "a" "a" "t" "g"
    [613] "a" "g" "g" "t" "a" "g" "a" "t" "g" "a" "a" "t" "a" "t" "t" "c" "t" "g"
    [631] "g" "t" "t" "c" "t" "t" "c" "a" "g" "a" "g" "a" "a" "a" "a" "t" "a" "g"
    [649] "a" "c" "t" "t" "a" "c" "t" "g" "g" "c" "c" "a" "g" "t" "g" "a" "t" "c"
    [667] "c" "t" "c" "a" "t" "g" "a" "g" "g" "c" "t" "t" "t" "a" "a" "t" "a" "t"
    [685] "g" "t" "a" "a" "a" "a" "g" "t" "g" "a" "a" "a" "g" "a" "g" "t" "t" "c"
    [703] "a" "c" "t" "c" "c" "a" "a" "a" "t" "c" "a" "g" "t" "a" "g" "a" "g" "a"
    [721] "g" "t" "a" "a" "t" "a" "t" "t" "g" "a" "a" "g" "a" "c" "a" "a" "a" "a"
    [739] "t" "a" "t" "t" "t" "g" "g" "g" "a" "a" "a" "a" "c" "c" "t" "a" "t" "c"
    [757] "g" "g" "a" "a" "g" "a" "a" "g" "g" "c" "a" "a" "g" "c" "c" "t" "c" "c"
    [775] "c" "c" "a" "a" "c"
    attr(,"name")
    [1] "BC106746"
    attr(,"Annot")
    [1] ">BC106746"
    attr(,"class")
    [1] "SeqFastadna"
    
    $BC115037.BRCA1
       [1] "a" "t" "g" "a" "t" "a" "a" "a" "t" "c" "a" "g" "g" "g" "a" "a" "c" "t"
      [19] "a" "a" "c" "c" "a" "a" "a" "c" "g" "g" "a" "g" "c" "a" "g" "a" "a" "t"
      [37] "g" "g" "t" "c" "a" "a" "g" "t" "g" "a" "t" "g" "a" "a" "t" "a" "t" "t"
      [55] "a" "c" "t" "a" "a" "t" "a" "g" "t" "g" "g" "t" "c" "a" "t" "g" "a" "g"
      [73] "a" "a" "t" "a" "a" "a" "a" "c" "a" "a" "a" "a" "g" "g" "t" "g" "a" "t"
      [91] "t" "c" "t" "a" "t" "t" "c" "a" "g" "a" "a" "t" "g" "a" "g" "a" "a" "a"
     [109] "a" "a" "t" "c" "c" "t" "a" "a" "c" "c" "c" "a" "a" "t" "a" "g" "a" "a"
     [127] "t" "c" "a" "c" "t" "c" "g" "a" "a" "a" "a" "a" "g" "a" "a" "t" "c" "t"
     [145] "g" "c" "t" "t" "t" "c" "a" "a" "a" "a" "c" "g" "a" "a" "a" "g" "c" "t"
     [163] "g" "a" "a" "c" "c" "t" "a" "t" "a" "a" "g" "c" "a" "g" "c" "a" "g" "t"
     [181] "a" "t" "a" "a" "g" "c" "a" "a" "t" "a" "t" "g" "g" "a" "a" "c" "t" "c"
     [199] "g" "a" "a" "t" "t" "a" "a" "a" "t" "a" "t" "c" "c" "a" "c" "a" "a" "t"
     [217] "t" "c" "a" "a" "a" "a" "g" "c" "a" "c" "c" "t" "a" "a" "a" "a" "a" "g"
     [235] "a" "a" "t" "a" "g" "g" "c" "t" "g" "a" "g" "g" "a" "g" "g" "a" "a" "g"
     [253] "t" "c" "t" "t" "c" "t" "a" "c" "c" "a" "g" "g" "c" "a" "t" "a" "t" "t"
     [271] "c" "a" "t" "g" "c" "g" "c" "t" "t" "g" "a" "a" "c" "t" "a" "g" "t" "a"
     [289] "g" "t" "c" "a" "g" "t" "a" "g" "a" "a" "a" "t" "c" "t" "a" "a" "g" "c"
     [307] "c" "c" "a" "c" "c" "t" "a" "a" "t" "t" "g" "t" "a" "c" "t" "g" "a" "a"
     [325] "t" "t" "g" "c" "a" "a" "a" "t" "t" "g" "a" "t" "a" "g" "t" "t" "g" "t"
     [343] "t" "c" "t" "a" "g" "c" "a" "g" "t" "g" "a" "a" "g" "a" "g" "a" "t" "a"
     [361] "a" "a" "g" "a" "a" "a" "a" "a" "a" "a" "a" "g" "t" "a" "c" "a" "a" "c"
     [379] "c" "a" "a" "a" "t" "g" "c" "c" "a" "g" "t" "c" "a" "g" "g" "c" "a" "c"
     [397] "a" "g" "c" "a" "g" "a" "a" "a" "c" "c" "t" "a" "c" "a" "a" "c" "t" "c"
     [415] "a" "t" "g" "g" "a" "a" "g" "g" "t" "a" "a" "a" "g" "a" "a" "c" "c" "t"
     [433] "g" "c" "a" "a" "c" "t" "g" "g" "a" "g" "c" "c" "a" "a" "g" "a" "a" "g"
     [451] "a" "g" "t" "a" "a" "c" "a" "a" "g" "c" "c" "a" "a" "a" "t" "g" "a" "a"
     [469] "c" "a" "g" "a" "c" "a" "a" "g" "t" "a" "a" "a" "a" "g" "a" "c" "a" "t"
     [487] "g" "a" "c" "a" "g" "t" "g" "a" "t" "a" "c" "t" "t" "t" "c" "c" "c" "a"
     [505] "g" "a" "g" "c" "t" "g" "a" "a" "g" "t" "t" "a" "a" "c" "a" "a" "a" "t"
     [523] "g" "c" "a" "c" "c" "t" "g" "g" "t" "t" "c" "t" "t" "t" "t" "a" "c" "t"
     [541] "a" "a" "g" "t" "g" "t" "c" "c" "a" "a" "a" "t" "a" "c" "c" "a" "g" "t"
     [559] "g" "a" "a" "c" "t" "t" "a" "a" "a" "g" "a" "a" "t" "t" "t" "g" "t" "c"
     [577] "a" "a" "t" "c" "c" "t" "a" "g" "c" "c" "t" "t" "c" "c" "a" "a" "g" "a"
     [595] "g" "a" "a" "g" "a" "a" "a" "a" "a" "g" "a" "a" "g" "a" "g" "a" "a" "a"
     [613] "c" "t" "a" "g" "a" "a" "a" "c" "a" "g" "t" "t" "a" "a" "a" "g" "t" "g"
     [631] "t" "c" "t" "a" "a" "t" "a" "a" "t" "g" "c" "t" "g" "a" "a" "g" "a" "c"
     [649] "c" "c" "c" "a" "a" "a" "g" "a" "t" "c" "t" "c" "a" "t" "g" "t" "t" "a"
     [667] "a" "g" "t" "g" "g" "a" "g" "a" "a" "a" "g" "g" "g" "t" "t" "t" "t" "g"
     [685] "c" "a" "a" "a" "c" "t" "g" "a" "a" "a" "g" "a" "t" "c" "t" "g" "t" "a"
     [703] "g" "a" "g" "a" "g" "t" "a" "g" "c" "a" "g" "t" "a" "t" "t" "t" "c" "a"
     [721] "c" "t" "g" "g" "t" "a" "c" "c" "t" "g" "g" "t" "a" "c" "t" "g" "a" "t"
     [739] "t" "a" "t" "g" "g" "c" "a" "c" "t" "c" "a" "g" "g" "a" "a" "a" "g" "t"
     [757] "a" "t" "c" "t" "c" "g" "t" "t" "a" "c" "t" "g" "g" "a" "a" "g" "t" "t"
     [775] "a" "g" "c" "a" "c" "t" "c" "t" "a" "g" "g" "g" "a" "a" "g" "g" "c" "a"
     [793] "a" "a" "a" "a" "c" "a" "g" "a" "a" "c" "c" "a" "a" "a" "t" "a" "a" "a"
     [811] "t" "g" "t" "g" "t" "g" "a" "g" "t" "c" "a" "g" "t" "g" "t" "g" "c" "a"
     [829] "g" "c" "a" "t" "t" "t" "g" "a" "a" "a" "a" "c" "c" "c" "c" "a" "a" "g"
     [847] "g" "g" "a" "c" "t" "a" "a" "t" "t" "c" "a" "t" "g" "g" "t" "t" "g" "t"
     [865] "t" "c" "c" "a" "a" "a" "g" "a" "t" "a" "a" "t" "a" "g" "a" "a" "a" "t"
     [883] "g" "a" "c" "a" "c" "a" "g" "a" "a" "g" "g" "c" "t" "t" "t" "a" "a" "g"
     [901] "t" "a" "t" "c" "c" "a" "t" "t" "g" "g" "g" "a" "c" "a" "t" "g" "a" "a"
     [919] "g" "t" "t" "a" "a" "c" "c" "a" "c" "a" "g" "t" "c" "g" "g" "g" "a" "a"
     [937] "a" "c" "a" "a" "g" "c" "a" "t" "a" "g" "a" "a" "a" "t" "g" "g" "a" "a"
     [955] "g" "a" "a" "a" "g" "t" "g" "a" "a" "c" "t" "c" "g" "a" "t" "g" "c" "t"
     [973] "c" "a" "g" "t" "a" "t" "t" "t" "g" "c" "a" "g" "a" "a" "t" "a" "c" "a"
     [991] "t" "t" "c" "a" "a" "g" "g" "t" "t" "t" "c" "a" "a" "a" "g" "c" "g" "c"
    [1009] "c" "a" "g" "t" "c" "a" "t" "t" "t" "g" "c" "t" "c" "t" "g" "t" "t" "t"
    [1027] "t" "c" "a" "a" "a" "t" "c" "c" "a" "g" "g" "a" "a" "a" "t" "g" "c" "a"
    [1045] "g" "a" "a" "g" "a" "g" "g" "a" "a" "t" "g" "t" "g" "c" "a" "a" "c" "a"
    [1063] "t" "t" "c" "t" "c" "t" "g" "c" "c" "c" "a" "c" "t" "c" "t" "g" "g" "g"
    [1081] "t" "c" "c" "t" "t" "a" "a" "a" "g" "a" "a" "a" "c" "a" "a" "a" "g" "t"
    [1099] "c" "c" "a" "a" "a" "a" "g" "t" "c" "a" "c" "t" "t" "t" "t" "g" "a" "a"
    [1117] "t" "g" "t" "g" "a" "a" "c" "a" "a" "a" "a" "g" "g" "a" "a" "g" "a" "a"
    [1135] "a" "a" "t" "c" "a" "a" "g" "g" "a" "a" "a" "g" "a" "a" "t" "g" "a" "g"
    [1153] "t" "c" "t" "a" "a" "t" "a" "t" "c" "a" "a" "g" "c" "c" "t" "g" "t" "a"
    [1171] "c" "a" "g" "a" "c" "a" "g" "t" "t" "a" "a" "t" "a" "t" "c" "a" "c" "t"
    [1189] "g" "c" "a" "g" "g" "c" "t" "t" "t" "c" "c" "t" "g" "t" "g" "g" "t" "t"
    [1207] "g" "g" "t" "c" "a" "g" "a" "a" "a" "g" "a" "t" "a" "a" "g" "c" "c" "a"
    [1225] "g" "t" "t" "g" "a" "t" "a" "a" "t" "g" "c" "c" "a" "a" "a" "t" "g" "t"
    [1243] "a" "g" "t" "a" "t" "c" "a" "a" "a" "g" "g" "a" "g" "g" "c" "t" "c" "t"
    [1261] "a" "g" "g" "t" "t" "t" "t" "g" "t" "c" "t" "a" "t" "c" "a" "t" "c" "t"
    [1279] "c" "a" "g" "t" "t" "c" "a" "g" "a" "g" "g" "c" "a" "a" "c" "g" "a" "a"
    [1297] "a" "c" "t" "g" "g" "a" "c" "t" "c" "a" "t" "t" "a" "c" "t" "c" "c" "a"
    [1315] "a" "a" "t" "a" "a" "a" "c" "a" "t" "g" "g" "a" "c" "t" "t" "t" "t" "a"
    [1333] "c" "a" "a" "a" "a" "c" "c" "c" "a" "t" "a" "t" "c" "g" "t" "a" "t" "a"
    [1351] "c" "c" "a" "c" "c" "a" "c" "t" "t" "t" "t" "t" "c" "c" "c" "a" "t" "c"
    [1369] "a" "a" "g" "t" "c" "a" "t" "t" "t" "g" "t" "t" "a" "a" "a" "a" "c" "t"
    [1387] "a" "a" "a" "t" "g" "t" "a" "a" "g" "a" "a" "a" "a" "a" "t" "c" "t" "g"
    [1405] "c" "t" "a" "g" "a" "g" "g" "a" "a" "a" "a" "c" "t" "t" "t" "g" "a" "g"
    [1423] "g" "a" "a" "c" "a" "t" "t" "c" "a" "a" "t" "g" "t" "c" "a" "c" "c" "t"
    [1441] "g" "a" "a" "a" "g" "a" "g" "a" "a" "a" "t" "g" "g" "g" "a" "a" "a" "t"
    [1459] "g" "a" "g" "a" "a" "c" "a" "t" "t" "c" "c" "a" "a" "g" "t" "a" "c" "a"
    [1477] "g" "t" "g" "a" "g" "c" "a" "c" "a" "a" "t" "t" "a" "g" "c" "c" "g" "t"
    [1495] "a" "a" "t" "a" "a" "c" "a" "t" "t" "a" "g" "a" "g" "a" "a" "a" "a" "t"
    [1513] "g" "t" "t" "t" "t" "t" "a" "a" "a" "g" "g" "a" "g" "c" "c" "a" "g" "c"
    [1531] "t" "c" "a" "a" "g" "c" "a" "a" "t" "a" "t" "t" "a" "a" "t" "g" "a" "a"
    [1549] "g" "t" "a" "g" "g" "t" "t" "c" "c" "a" "g" "t" "a" "c" "t" "a" "a" "t"
    [1567] "g" "a" "a" "g" "t" "g" "g" "g" "c" "t" "c" "c" "a" "g" "t" "a" "t" "t"
    [1585] "a" "a" "t" "g" "a" "a" "a" "t" "a" "g" "g" "t" "t" "c" "c" "a" "g" "t"
    [1603] "g" "a" "t" "g" "a" "a" "a" "a" "c" "a" "t" "t" "c" "a" "a" "g" "c" "a"
    [1621] "g" "a" "a" "c" "t" "a" "g" "g" "t" "a" "g" "a" "a" "a" "c" "a" "g" "a"
    [1639] "g" "g" "g" "c" "c" "a" "a" "a" "a" "t" "t" "g" "a" "a" "t" "g" "c" "t"
    [1657] "a" "t" "g" "c" "t" "t" "a" "g" "a" "t" "t" "a" "g" "g" "g" "g" "t" "t"
    [1675] "t" "t" "g" "c" "a" "a" "c" "c" "t" "g" "a" "g" "g" "t" "c" "t" "a" "t"
    [1693] "a" "a" "a" "c" "a" "a" "a" "g" "t" "c" "t" "t" "c" "c" "t" "g" "g" "a"
    [1711] "a" "g" "t" "a" "a" "t" "t" "g" "t" "a" "a" "g" "c" "a" "t" "c" "c" "t"
    [1729] "g" "a" "a" "a" "t" "a" "a" "a" "a" "a" "a" "g" "c" "a" "a" "g" "a" "a"
    [1747] "t" "a" "t" "g" "a" "a" "g" "a" "a" "g" "t" "a" "g" "t" "t" "c" "a" "g"
    [1765] "a" "c" "t" "g" "t" "t" "a" "a" "t" "a" "c" "a" "g" "a" "t" "t" "t" "c"
    [1783] "t" "c" "t" "c" "c" "a" "t" "a" "t" "c" "t" "g" "a" "t" "t" "t" "c" "a"
    [1801] "g" "a" "t" "a" "a" "c" "t" "t" "a" "g" "a" "a" "c" "a" "g" "c" "c" "t"
    [1819] "a" "t" "g" "g" "g" "a" "a" "g" "t" "a" "g" "t" "c" "a" "t" "g" "c" "a"
    [1837] "t" "c" "t" "c" "a" "g" "g" "t" "t" "t" "g" "t" "t" "c" "t" "g" "a" "g"
    [1855] "a" "c" "a" "c" "c" "t" "g" "a" "t" "g" "a" "c" "c" "t" "g" "t" "t" "a"
    [1873] "g" "a" "t" "g" "a" "t" "g" "g" "t" "g" "a" "a" "a" "t" "a" "a" "a" "g"
    [1891] "g" "a" "a" "g" "a" "t" "a" "c" "t" "a" "g" "t" "t" "t" "t" "g" "c" "t"
    [1909] "g" "a" "a" "a" "a" "t" "g" "a" "c" "a" "t" "t" "a" "a" "g" "g" "a" "a"
    [1927] "a" "g" "t" "t" "c" "t" "g" "c" "t" "g" "t" "t" "t" "t" "t" "a" "g" "c"
    [1945] "a" "a" "a" "a" "g" "c" "g" "t" "c" "c" "a" "g" "a" "g" "a" "g" "g" "a"
    [1963] "g" "a" "g" "c" "t" "t" "a" "g" "c" "a" "g" "g" "a" "g" "t" "c" "c" "t"
    [1981] "a" "g" "c" "c" "c" "t" "t" "t" "c" "a" "c" "c" "c" "a" "t" "a" "c" "a"
    [1999] "c" "a" "t" "t" "t" "g" "g" "c" "t" "c" "a" "g" "g" "g" "t" "t" "a" "c"
    [2017] "c" "g" "a" "a" "g" "a" "g" "g" "g" "g" "c" "c" "a" "a" "g" "a" "a" "a"
    [2035] "t" "t" "a" "g" "a" "g" "t" "c" "c" "t" "c" "a" "g" "a" "a" "g" "a" "g"
    [2053] "a" "a" "c" "t" "t" "a" "t" "c" "t" "a" "g" "t" "g" "a" "g" "g" "a" "t"
    [2071] "g" "a" "a" "g" "a" "g" "c" "t" "t" "c" "c" "c" "t" "g" "c" "t" "t" "c"
    [2089] "c" "a" "a" "c" "a" "c" "t" "t" "g" "t" "t" "a" "t" "t" "t" "g" "g" "t"
    [2107] "a" "a" "a" "g" "t" "a" "a" "a" "c" "a" "a" "t" "a" "t" "a" "c" "c" "t"
    [2125] "t" "c" "t" "c" "a" "g" "t" "c" "t" "a" "c" "t" "a" "g" "g" "c" "a" "t"
    [2143] "a" "g" "c" "a" "c" "c" "g" "t" "t" "g" "c" "t" "a" "c" "c" "g" "a" "g"
    [2161] "t" "g" "t" "c" "t" "g" "t" "c" "t" "a" "a" "g" "a" "a" "c" "a" "c" "a"
    [2179] "g" "a" "g" "g" "a" "g" "a" "a" "t" "t" "t" "a" "t" "t" "a" "t" "c" "a"
    [2197] "t" "t" "g" "a" "a" "g" "a" "a" "t" "a" "g" "c" "t" "t" "a" "a" "a" "t"
    [2215] "g" "a" "c" "t" "g" "c" "a" "g" "t" "a" "a" "c" "c" "a" "g" "g" "t" "a"
    [2233] "a" "t" "a" "t" "t" "g" "g" "c" "a" "a" "a" "g" "g" "c" "a" "t" "c" "t"
    [2251] "c" "a" "g" "g" "a" "a" "c" "a" "t" "c" "a" "c" "c" "t" "t" "a" "g" "t"
    [2269] "g" "a" "g" "g" "a" "a" "a" "c" "a" "a" "a" "a" "t" "g" "t" "t" "c" "t"
    [2287] "g" "c" "t" "a" "g" "c" "t" "t" "g" "t" "t" "t" "t" "c" "t" "t" "c" "a"
    [2305] "c" "a" "g" "t" "g" "c" "a" "g" "t" "g" "a" "a" "t" "t" "g" "g" "a" "a"
    [2323] "g" "a" "c" "t" "t" "g" "a" "c" "t" "g" "c" "a" "a" "a" "t" "a" "c" "a"
    [2341] "a" "a" "c" "a" "c" "c" "c" "a" "g" "g" "a" "t" "c" "c" "t" "t" "t" "c"
    [2359] "t" "t" "g" "a" "t" "t" "g" "g" "t" "t" "c" "t" "t" "c" "c" "a" "a" "a"
    [2377] "c" "a" "a" "a" "t" "g" "a" "g" "g" "c" "a" "t" "c" "a" "g" "t" "c" "t"
    [2395] "g" "a" "a" "a" "g" "c" "c" "a" "g" "g" "g" "a" "g" "t" "t" "g" "g" "t"
    [2413] "c" "t" "g" "a" "g" "t" "g" "a" "c" "a" "a" "g" "g" "a" "a" "t" "t" "g"
    [2431] "g" "t" "t" "t" "c" "a" "g" "a" "t" "g" "a" "t" "g" "a" "a" "g" "a" "a"
    [2449] "a" "g" "a" "g" "g" "a" "a" "c" "g" "g" "g" "c" "t" "t" "g" "g" "a" "a"
    [2467] "g" "a" "a" "a" "a" "t" "a" "a" "t" "c" "a" "a" "g" "a" "a" "g" "a" "g"
    [2485] "c" "a" "a" "a" "g" "c" "a" "t" "g" "g" "a" "t" "t" "c" "a" "a" "a" "c"
    [2503] "t" "t" "a" "g" "g" "t" "g" "a" "a" "g" "c" "a" "g" "c" "a" "t" "c" "t"
    [2521] "g" "g" "g" "t" "g" "t" "g" "a" "g" "a" "g" "t" "g" "a" "a" "a" "c" "a"
    [2539] "a" "g" "c" "g" "t" "c" "t" "c" "t" "g" "a" "a" "g" "a" "c" "t" "g" "c"
    [2557] "t" "c" "a" "g" "g" "g" "c" "t" "a" "t" "c" "c" "t" "c" "t" "c" "a" "g"
    [2575] "a" "g" "t" "g" "a" "c" "a" "t" "t" "t" "t" "a" "a" "c" "c" "a" "c" "t"
    [2593] "c" "a" "g" "c" "a" "g" "a" "g" "g" "g" "a" "t" "a" "c" "c" "a" "t" "g"
    [2611] "c" "a" "a" "c" "a" "t" "a" "a" "c" "c" "t" "g" "a" "t" "a" "a" "a" "g"
    [2629] "c" "t" "c" "c" "a" "g" "c" "a" "g" "g" "a" "a" "a" "t" "g" "g" "c" "t"
    [2647] "g" "a" "a" "c" "t" "a" "g" "a" "a" "g" "c" "t" "g" "t" "g" "t" "t" "a"
    [2665] "g" "a" "a" "c" "a" "g" "c" "a" "t" "g" "g" "g" "a" "g" "c" "c" "a" "g"
    [2683] "c" "c" "t" "t" "c" "t" "a" "a" "c" "a" "g" "c" "t" "a" "c" "c" "c" "t"
    [2701] "t" "c" "c" "a" "t" "c" "a" "t" "a" "a" "g" "t" "g" "a" "c" "t" "c" "c"
    [2719] "t" "c" "t" "g" "c" "c" "c" "t" "t" "g" "a" "g" "g" "a" "c" "c" "t" "g"
    [2737] "c" "g" "a" "a" "a" "t" "c" "c" "a" "g" "a" "a" "c" "a" "a" "a" "g" "c"
    [2755] "a" "c" "a" "t" "c" "a" "g" "a" "a" "a" "a" "a" "g" "a" "t" "t" "c" "g"
    [2773] "c" "a" "t" "a" "t" "a" "c" "a" "t" "g" "g" "c" "c" "a" "a" "a" "g" "g"
    [2791] "g" "a" "c" "a" "a" "c" "t" "c" "c" "a" "t" "g" "t" "t" "t" "t" "c" "t"
    [2809] "a" "a" "a" "a" "g" "g" "c" "c" "t" "a" "g" "a" "g" "a" "a" "c" "a" "t"
    [2827] "a" "t" "a" "t" "c" "a" "g" "t" "a" "t" "t" "a" "a" "c" "t" "t" "c" "a"
    [2845] "c" "a" "g" "a" "a" "a" "a" "g" "t" "a" "g" "t" "g" "a" "a" "t" "a" "c"
    [2863] "c" "c" "t" "a" "t" "a" "a" "g" "c" "c" "a" "g" "a" "a" "t" "c" "c" "a"
    [2881] "g" "a" "a" "g" "g" "c" "c" "t" "t" "t" "c" "t" "g" "c" "t" "g" "a" "c"
    [2899] "a" "a" "g" "t" "t" "t" "g" "a" "g" "g" "t" "g" "t" "c" "t" "g" "c" "a"
    [2917] "g" "a" "t" "a" "g" "t" "t" "c" "t" "a" "c" "c" "a" "g" "t" "a" "a" "a"
    [2935] "a" "a" "t" "a" "a" "a" "g" "a" "a" "c" "c" "a" "g" "g" "a" "g" "t" "g"
    [2953] "g" "a" "a" "a" "g" "g" "t" "c" "a" "t" "c" "c" "c" "c" "t" "t" "c" "t"
    [2971] "a" "a" "a" "t" "g" "c" "c" "c" "a" "t" "c" "a" "t" "t" "a" "g" "a" "t"
    [2989] "g" "a" "t" "a" "g" "g" "t" "g" "g" "t" "a" "c" "a" "t" "g" "c" "a" "c"
    [3007] "a" "g" "t" "t" "g" "c" "t" "c" "t" "g" "g" "g" "a" "g" "t" "c" "t" "t"
    [3025] "c" "a" "g" "a" "a" "t" "a" "g" "a" "a" "a" "c" "t" "a" "c" "c" "c" "a"
    [3043] "t" "c" "t" "c" "a" "a" "g" "a" "g" "g" "g" "g" "c" "t" "c" "a" "t" "t"
    [3061] "a" "a" "g" "g" "t" "t" "g" "t" "t" "g" "a" "t" "g" "t" "g" "g" "a" "g"
    [3079] "g" "a" "g" "c" "a" "a" "c" "a" "g" "c" "t" "g" "g" "a" "a" "g" "a" "g"
    [3097] "t" "c" "t" "g" "g" "g" "c" "c" "a" "c" "a" "c" "g" "a" "t" "t" "t" "g"
    [3115] "a" "c" "g" "g" "a" "a" "a" "c" "a" "t" "c" "t" "t" "a" "c" "t" "t" "g"
    [3133] "c" "c" "a" "a" "g" "g" "c" "a" "a" "g" "a" "t" "c" "t" "a" "g" "a" "g"
    [3151] "g" "g" "a" "a" "c" "c" "c" "c" "t" "t" "a" "c" "c" "t" "g" "g" "a" "a"
    [3169] "t" "c" "t" "g" "g" "a" "a" "t" "c" "a" "g" "c" "c" "t" "c" "t" "t" "c"
    [3187] "t" "c" "t" "g" "a" "t" "g" "a" "c" "c" "c" "t" "g" "a" "a" "t" "c" "t"
    [3205] "g" "a" "t" "c" "c" "t" "t" "c" "t" "g" "a" "a" "g" "a" "c" "a" "g" "a"
    [3223] "g" "c" "c" "c" "c" "a" "g" "a" "g" "t" "c" "a" "g" "c" "t" "c" "g" "t"
    [3241] "g" "t" "t" "g" "g" "c" "a" "a" "c" "a" "t" "a" "c" "c" "a" "t" "c" "t"
    [3259] "t" "c" "a" "a" "c" "c" "t" "c" "t" "g" "c" "a" "t" "t" "g" "a" "a" "a"
    [3277] "g" "t" "t" "c" "c" "c" "c" "a" "a" "t" "t" "g" "a" "a" "a" "g" "t" "t"
    [3295] "g" "c" "a" "g" "a" "a" "t" "c" "t" "g" "c" "c" "c" "a" "g" "g" "g" "t"
    [3313] "c" "c" "a" "g" "c" "t" "g" "c" "t" "g" "c" "t" "c" "a" "t" "a" "c" "t"
    [3331] "a" "c" "t" "g" "a" "t" "a" "c" "t" "g" "c" "t" "g" "g" "g" "t" "a" "t"
    [3349] "a" "a" "t" "g" "c" "a" "a" "t" "g" "g" "a" "a" "g" "a" "a" "a" "g" "t"
    [3367] "g" "t" "g" "a" "g" "c" "a" "g" "g" "g" "a" "g" "a" "a" "g" "c" "c" "a"
    [3385] "g" "a" "a" "t" "t" "g" "a" "c" "a" "g" "c" "t" "t" "c" "a" "a" "c" "a"
    [3403] "g" "a" "a" "a" "g" "g" "g" "t" "c" "a" "a" "c" "a" "a" "a" "a" "g" "a"
    [3421] "a" "t" "g" "t" "c" "c" "a" "t" "a" "g" "t" "g" "g" "t" "g" "t" "c" "t"
    [3439] "g" "g" "c" "c" "t" "g" "a" "c" "c" "c" "c" "a" "g" "a" "a" "g" "a" "a"
    [3457] "t" "t" "t" "a" "t" "g" "c" "t" "c" "g" "t" "g" "t" "a" "c" "a" "a" "g"
    [3475] "t" "t" "t" "g" "c" "c" "a" "g" "a" "a" "a" "a" "c" "a" "c" "c" "a" "c"
    [3493] "a" "t" "c" "a" "c" "t" "t" "t" "a" "a" "c" "t" "a" "a" "t" "c" "t" "a"
    [3511] "a" "t" "t" "a" "c" "t" "g" "a" "a" "g" "a" "g" "a" "c" "t" "a" "c" "t"
    [3529] "c" "a" "t" "g" "t" "t" "g" "t" "t" "a" "t" "g" "a" "a" "a" "a" "c" "a"
    [3547] "g" "a" "t" "g" "c" "t" "g" "a" "g" "t" "t" "t" "g" "t" "g" "t" "g" "t"
    [3565] "g" "a" "a" "c" "g" "g" "a" "c" "a" "c" "t" "g" "a" "a" "a" "t" "a" "t"
    [3583] "t" "t" "t" "c" "t" "a" "g" "g" "a" "a" "t" "t" "g" "c" "g" "g" "g" "a"
    [3601] "g" "g" "a" "a" "a" "a" "t" "g" "g" "g" "t" "a" "g" "t" "t" "a" "g" "c"
    [3619] "t" "a" "t" "t" "t" "c" "t" "g" "g" "g" "t" "g" "a" "c" "c" "c" "a" "g"
    [3637] "t" "c" "t" "a" "t" "t" "a" "a" "a" "g" "a" "a" "a" "g" "a" "a" "a" "a"
    [3655] "a" "t" "g" "c" "t" "g" "a" "a" "t" "g" "a" "g" "c" "a" "t" "g" "a" "t"
    [3673] "t" "t" "t" "g" "a" "a" "g" "t" "c" "a" "g" "a" "g" "g" "a" "g" "a" "t"
    [3691] "g" "t" "g" "g" "t" "c" "a" "a" "t" "g" "g" "a" "a" "g" "a" "a" "a" "c"
    [3709] "c" "a" "c" "c" "a" "a" "g" "g" "t" "c" "c" "a" "a" "a" "g" "c" "g" "a"
    [3727] "g" "c" "a" "a" "g" "a" "g" "a" "a" "t" "c" "c" "c" "a" "g" "g" "a" "c"
    [3745] "a" "g" "a" "a" "a" "g" "a" "t" "c" "t" "t" "c" "a" "g" "g" "g" "g" "g"
    [3763] "c" "t" "a" "g" "a" "a" "a" "t" "c" "t" "g" "t" "t" "g" "c" "t" "a" "t"
    [3781] "g" "g" "g" "c" "c" "c" "t" "t" "c" "a" "c" "c" "a" "a" "c" "a" "t" "g"
    [3799] "c" "c" "c" "a" "c" "a" "g" "a" "t" "c" "a" "a" "c" "t" "g" "g" "a" "a"
    [3817] "t" "g" "g" "a" "t" "g" "g" "t" "a" "c" "a" "g" "c" "t" "g" "t" "g" "t"
    [3835] "g" "g" "t" "g" "c" "t" "t" "c" "t" "g" "t" "g" "g" "t" "g" "a" "a" "g"
    [3853] "g" "a" "g" "c" "t" "t" "t" "c" "a" "t" "c" "a" "t" "t" "c" "a" "c" "c"
    [3871] "c" "t" "t" "g" "g" "c" "a" "c" "a" "g" "g" "t" "g" "t" "c" "c" "a" "c"
    [3889] "c" "c" "a" "a" "t" "t" "g" "t" "g" "g" "t" "t" "g" "t" "g" "c" "a" "g"
    [3907] "c" "c" "a" "g" "a" "t" "g" "c" "c" "t" "g" "g" "a" "c" "a" "g" "a" "g"
    [3925] "g" "a" "c" "a" "a" "t" "g" "g" "c" "t" "t" "c" "c" "a" "t" "g" "c" "a"
    [3943] "a" "t" "t" "g" "g" "g" "c" "a" "g" "a" "t" "g" "t" "g" "t" "g" "a" "g"
    [3961] "g" "c" "a" "c" "c" "t" "g" "t" "g" "g" "t" "g" "a" "c" "c" "c" "g" "a"
    [3979] "g" "a" "g" "t" "g" "g" "g" "t" "g" "t" "t" "g" "g" "a" "c" "a" "g" "t"
    [3997] "g" "t" "a" "g" "c" "a" "c" "t" "c" "t" "a" "c" "c" "a" "g" "t" "g" "c"
    [4015] "c" "a" "g" "g" "a" "g" "c" "t" "g" "g" "a" "c" "a" "c" "c" "t" "a" "c"
    [4033] "c" "t" "g" "a" "t" "a" "c" "c" "c" "c" "a" "g" "a" "t" "c" "c" "c" "c"
    [4051] "c" "a" "c" "a" "g" "c" "c" "a" "c" "t" "a" "c" "t" "g" "a"
    attr(,"name")
    [1] "BC115037.BRCA1"
    attr(,"Annot")
    [1] ">BC115037.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $DQ075361
     [1] "a" "g" "g" "g" "a" "a" "c" "c" "c" "c" "t" "t" "a" "c" "c" "t" "g" "g" "a"
    [20] "a" "t" "c" "t" "g" "g" "a" "a" "t" "c" "t" "g" "g" "a" "a" "t" "c" "a" "g"
    [39] "c" "c" "t" "c" "t" "t" "c" "t" "c" "t" "g" "a"
    attr(,"name")
    [1] "DQ075361"
    attr(,"Annot")
    [1] ">DQ075361"
    attr(,"class")
    [1] "SeqFastadna"
    
    $DQ116737
      [1] "a" "a" "g" "c" "g" "t" "g" "c" "a" "g" "c" "t" "g" "a" "g" "a" "g" "g"
     [19] "c" "a" "t" "c" "c" "a" "g" "a" "a" "a" "a" "g" "t" "a" "t" "c" "a" "g"
     [37] "g" "g" "t" "a" "g" "t" "t" "c" "t" "g" "t" "t" "t" "c" "a" "a" "a" "c"
     [55] "t" "t" "g" "c" "a" "t" "g" "t" "g" "g" "a" "g" "c" "c" "a" "t" "g" "t"
     [73] "g" "g" "c" "a" "c" "a" "a" "a" "t" "a" "c" "t" "c" "a" "t" "g" "c" "c"
     [91] "a" "g" "c" "t" "c" "a" "t" "t" "a" "c" "a" "g" "c" "a" "t" "g" "a" "g"
    [109] "a" "a" "c" "a" "g" "c" "a" "g" "t" "t" "t" "a" "t" "t" "a" "c" "t" "c"
    [127] "a" "c" "t" "a" "a" "a" "g" "a" "c" "a" "g" "a" "a" "t" "g" "a" "a" "t"
    [145] "g" "t" "a" "g" "a" "a" "a" "a" "g" "g" "c" "t" "g" "a" "a" "t" "t" "c"
    [163] "t" "g" "t" "a" "a" "t" "a" "a" "a" "a" "g" "c" "a" "a" "a" "c" "a" "g"
    [181] "c" "c" "t" "g" "g" "c" "t" "t" "a" "g" "c" "a" "a" "g" "g" "a" "g" "c"
    [199] "c" "a" "a" "c" "a" "t" "a" "a" "c" "a" "g" "a" "t" "g" "g" "g" "c" "t"
    [217] "g" "g" "a" "a" "g" "t" "a" "a" "g" "g" "a" "a" "a" "c" "a" "t" "g" "t"
    [235] "a" "a" "t" "g" "a" "t" "a" "g" "g" "c" "g" "g" "a" "c" "t" "c" "c" "c"
    [253] "a" "g" "c" "a" "c" "a" "g" "a" "a" "a" "a" "a" "a" "a" "g" "g" "t" "a"
    [271] "g" "a" "t" "c" "t" "g" "a" "a" "t" "g" "c" "t" "g" "a" "t" "c" "c" "c"
    [289] "c" "t" "g" "t" "g" "t" "g" "a" "g" "a" "g" "a" "a" "a" "a" "g" "a" "a"
    [307] "t" "g" "g" "a" "a" "t" "a" "a" "a" "c" "a" "g" "a" "a" "a" "c" "t" "g"
    [325] "c" "c" "a" "t" "g" "c" "t" "c" "a" "g" "a" "g" "a" "a" "t" "c" "c" "t"
    [343] "a" "g" "a" "g" "a" "t" "a" "c" "t" "g" "a" "a" "g" "a" "t" "g" "t" "t"
    [361] "c" "c" "t" "t" "g" "g" "a" "t" "a" "a" "c" "a" "c" "t" "a" "a" "a" "t"
    [379] "a" "g" "c" "a" "g" "c" "a" "t" "t" "c" "a" "g" "a" "a" "a" "g" "t" "t"
    [397] "a" "a" "t" "g" "a" "g" "t" "g" "g" "t" "t" "t" "t" "c" "c" "a" "g" "a"
    [415] "a" "g" "t" "g" "a" "t" "g" "a" "a" "c" "t" "g" "t" "t" "a" "g" "g" "t"
    [433] "t" "c" "t" "g" "a" "t" "g" "a" "c" "t" "c" "a" "c" "a" "t" "g" "a" "t"
    [451] "g" "g" "g" "g" "a" "g" "t" "c" "t" "g" "a" "a" "t" "c" "a" "a" "a" "t"
    [469] "g" "c" "c" "a" "a" "a" "g" "t" "a" "g" "c" "t" "g" "a" "t" "g" "t" "a"
    [487] "t" "t" "g" "g" "a" "c" "g" "t" "t" "c" "t" "a"
    attr(,"name")
    [1] "DQ116737"
    attr(,"Annot")
    [1] ">DQ116737"
    attr(,"class")
    [1] "SeqFastadna"
    
    $DQ145822
      [1] "a" "a" "a" "t" "g" "a" "c" "a" "t" "t" "a" "a" "g" "g" "a" "a" "a" "g"
     [19] "t" "t" "c" "t" "g" "c" "t" "g" "t" "t" "t" "t" "t" "a" "g" "c" "a" "a"
     [37] "a" "a" "g" "c" "g" "t" "c" "c" "a" "g" "a" "a" "a" "g" "g" "a" "g" "a"
     [55] "g" "c" "t" "t" "a" "g" "c" "a" "g" "g" "a" "g" "t" "c" "c" "t" "a" "g"
     [73] "c" "c" "c" "t" "t" "t" "c" "a" "c" "c" "c" "a" "t" "a" "c" "a" "c" "a"
     [91] "t" "t" "t" "g" "g" "c" "t" "c" "a" "g" "g" "g" "t" "t" "a" "c" "c" "g"
    [109] "a" "a" "g" "a" "g" "g" "g" "g" "c" "c" "a" "a"
    attr(,"name")
    [1] "DQ145822"
    attr(,"Annot")
    [1] ">DQ145822"
    attr(,"class")
    [1] "SeqFastadna"
    
    $DQ145823
      [1] "a" "g" "c" "a" "c" "c" "g" "t" "t" "g" "c" "t" "a" "c" "c" "g" "a" "g"
     [19] "t" "g" "t" "c" "t" "g" "t" "c" "t" "a" "a" "g" "a" "a" "c" "a" "c" "a"
     [37] "g" "a" "g" "g" "a" "g" "a" "a" "t" "t" "t" "a" "t" "t" "a" "t" "c" "a"
     [55] "t" "t" "g" "a" "a" "g" "a" "a" "t" "a" "g" "c" "t" "t" "a" "a" "a" "t"
     [73] "g" "a" "c" "t" "g" "c" "g" "g" "t" "a" "a" "c" "c" "a" "g" "g" "t" "a"
     [91] "a" "t" "a" "t" "t" "g" "g" "c" "a" "a" "a"
    attr(,"name")
    [1] "DQ145823"
    attr(,"Annot")
    [1] ">DQ145823"
    attr(,"class")
    [1] "SeqFastadna"
    
    $DQ145824
      [1] "c" "a" "a" "t" "a" "t" "a" "c" "c" "t" "t" "c" "t" "c" "a" "g" "t" "c"
     [19] "t" "a" "c" "t" "a" "g" "g" "c" "a" "t" "a" "g" "c" "a" "c" "c" "g" "t"
     [37] "t" "g" "c" "t" "a" "c" "c" "g" "a" "g" "t" "g" "t" "c" "t" "g" "t" "c"
     [55] "t" "a" "a" "g" "a" "a" "c" "a" "c" "a" "g" "a" "g" "g" "a" "g" "a" "a"
     [73] "t" "t" "t" "a" "t" "t" "a" "t" "c" "a" "t" "t" "g" "a" "a" "g" "a" "a"
     [91] "t" "a" "g" "c" "t" "t" "a" "a" "a" "t" "g" "a" "g" "t" "g" "c" "g" "g"
    [109] "t" "a" "a" "c" "c" "a"
    attr(,"name")
    [1] "DQ145824"
    attr(,"Annot")
    [1] ">DQ145824"
    attr(,"class")
    [1] "SeqFastadna"
    
    $DQ145825.BRCA1
      [1] "t" "g" "c" "t" "g" "a" "a" "a" "a" "t" "g" "a" "c" "a" "t" "t" "a" "a"
     [19] "g" "g" "a" "a" "a" "g" "t" "t" "c" "t" "g" "c" "t" "g" "t" "t" "t" "t"
     [37] "t" "a" "g" "c" "a" "a" "a" "a" "g" "c" "g" "t" "c" "c" "a" "g" "a" "a"
     [55] "a" "g" "g" "a" "g" "a" "g" "c" "t" "t" "a" "g" "c" "a" "g" "g" "a" "g"
     [73] "t" "c" "c" "t" "a" "g" "c" "c" "c" "t" "t" "t" "c" "a" "c" "c" "c" "a"
     [91] "t" "a" "c" "a" "c" "a" "t" "t" "t" "g" "g" "c" "t" "c" "a" "g"
    attr(,"name")
    [1] "DQ145825.BRCA1"
    attr(,"Annot")
    [1] ">DQ145825.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $DQ145826
      [1] "t" "t" "a" "a" "a" "g" "t" "g" "t" "c" "t" "a" "a" "t" "a" "a" "t" "g"
     [19] "c" "t" "g" "a" "a" "g" "a" "c" "c" "c" "c" "a" "a" "a" "g" "a" "t" "c"
     [37] "t" "c" "a" "t" "g" "t" "t" "a" "a" "g" "t" "g" "g" "a" "g" "a" "a" "a"
     [55] "g" "g" "g" "t" "t" "t" "t" "g" "c" "a" "a" "a" "c" "t" "g" "a" "a" "a"
     [73] "g" "a" "t" "c" "t" "g" "t" "a" "g" "a" "g" "a" "g" "t" "a" "g" "c" "a"
     [91] "g" "t" "a" "t" "t" "t" "c" "a" "t" "t" "g" "g" "t" "a" "c" "c" "t" "g"
    [109] "g" "t" "a" "c" "t" "g" "a" "t" "t" "a" "t" "g"
    attr(,"name")
    [1] "DQ145826"
    attr(,"Annot")
    [1] ">DQ145826"
    attr(,"class")
    [1] "SeqFastadna"
    
    $DQ190450.BRCA1
       [1] "a" "t" "g" "g" "a" "t" "t" "t" "a" "t" "c" "t" "g" "c" "t" "c" "t" "t"
      [19] "c" "g" "c" "g" "t" "t" "g" "a" "a" "g" "a" "a" "g" "t" "a" "c" "a" "a"
      [37] "a" "a" "t" "g" "t" "c" "a" "t" "t" "a" "a" "t" "g" "c" "t" "a" "t" "g"
      [55] "c" "a" "g" "a" "a" "a" "a" "t" "c" "t" "t" "a" "g" "a" "g" "t" "g" "t"
      [73] "c" "c" "c" "a" "t" "c" "t" "g" "t" "c" "t" "g" "g" "a" "g" "t" "t" "g"
      [91] "a" "t" "c" "a" "a" "g" "g" "a" "a" "c" "c" "t" "g" "t" "c" "t" "c" "c"
     [109] "a" "c" "a" "a" "a" "g" "t" "g" "t" "g" "a" "c" "c" "a" "c" "a" "t" "a"
     [127] "t" "t" "t" "t" "g" "c" "a" "a" "a" "t" "t" "t" "t" "g" "c" "a" "t" "g"
     [145] "c" "t" "g" "a" "a" "a" "c" "t" "t" "c" "t" "c" "a" "a" "c" "c" "a" "g"
     [163] "a" "a" "g" "a" "a" "a" "g" "g" "g" "c" "c" "t" "t" "c" "a" "c" "a" "g"
     [181] "t" "g" "t" "c" "c" "t" "t" "t" "a" "t" "g" "t" "a" "a" "g" "a" "a" "t"
     [199] "g" "a" "t" "a" "t" "a" "a" "c" "c" "a" "a" "a" "a" "g" "g" "a" "g" "c"
     [217] "c" "t" "a" "c" "a" "a" "g" "a" "a" "a" "g" "t" "a" "c" "g" "a" "g" "a"
     [235] "t" "t" "t" "a" "g" "t" "c" "a" "a" "c" "t" "t" "g" "t" "t" "g" "a" "a"
     [253] "g" "a" "g" "c" "t" "a" "t" "t" "g" "a" "a" "a" "a" "t" "c" "a" "t" "t"
     [271] "t" "g" "t" "g" "c" "t" "t" "t" "t" "c" "a" "g" "c" "t" "t" "g" "a" "c"
     [289] "a" "c" "a" "g" "g" "t" "t" "t" "g" "g" "a" "g" "t" "a" "t" "g" "c" "a"
     [307] "a" "a" "c" "a" "g" "c" "t" "a" "t" "a" "a" "t" "t" "t" "t" "g" "c" "a"
     [325] "a" "a" "a" "a" "a" "g" "g" "a" "a" "a" "a" "t" "a" "a" "c" "t" "c" "t"
     [343] "c" "c" "t" "g" "a" "a" "c" "a" "t" "c" "t" "a" "a" "a" "a" "g" "a" "t"
     [361] "g" "a" "a" "g" "t" "t" "t" "c" "t" "a" "t" "c" "a" "t" "c" "c" "a" "a"
     [379] "a" "g" "t" "a" "t" "g" "g" "g" "c" "t" "a" "c" "a" "g" "a" "a" "a" "c"
     [397] "c" "g" "t" "g" "c" "c" "a" "a" "a" "a" "g" "a" "c" "t" "t" "c" "t" "a"
     [415] "c" "a" "g" "a" "g" "t" "g" "a" "a" "c" "c" "c" "g" "a" "a" "a" "a" "t"
     [433] "c" "c" "t" "t" "c" "c" "t" "t" "g" "c" "a" "g" "g" "a" "a" "a" "c" "c"
     [451] "a" "g" "t" "c" "t" "c" "a" "g" "t" "g" "t" "c" "c" "a" "a" "c" "t" "c"
     [469] "t" "c" "t" "a" "a" "c" "c" "t" "t" "g" "g" "a" "a" "c" "t" "g" "t" "g"
     [487] "a" "g" "a" "a" "c" "t" "c" "t" "g" "a" "g" "g" "a" "c" "a" "a" "a" "g"
     [505] "c" "a" "g" "c" "g" "g" "a" "t" "a" "c" "a" "a" "c" "c" "t" "c" "a" "a"
     [523] "a" "a" "g" "a" "c" "g" "t" "c" "t" "g" "t" "c" "t" "a" "c" "a" "t" "t"
     [541] "g" "a" "a" "t" "t" "g" "g" "g" "a" "t" "c" "t" "g" "a" "t" "t" "c" "t"
     [559] "t" "c" "t" "g" "a" "a" "g" "a" "t" "a" "c" "c" "g" "t" "t" "a" "a" "t"
     [577] "a" "a" "g" "g" "c" "a" "a" "c" "t" "t" "a" "t" "t" "g" "c" "a" "g" "t"
     [595] "g" "t" "g" "g" "g" "a" "g" "a" "t" "c" "a" "a" "g" "a" "a" "t" "t" "g"
     [613] "t" "t" "a" "c" "a" "a" "a" "t" "c" "a" "c" "c" "c" "c" "t" "c" "a" "a"
     [631] "g" "g" "a" "a" "c" "c" "a" "g" "g" "g" "a" "t" "g" "a" "a" "a" "t" "c"
     [649] "a" "g" "t" "t" "t" "g" "g" "a" "t" "t" "c" "t" "g" "c" "a" "a" "a" "a"
     [667] "a" "a" "g" "g" "c" "t" "g" "c" "t" "t" "g" "t" "g" "a" "a" "t" "t" "t"
     [685] "t" "c" "t" "g" "a" "g" "a" "c" "g" "g" "a" "t" "g" "t" "a" "a" "c" "a"
     [703] "a" "a" "t" "a" "c" "t" "g" "a" "a" "c" "a" "t" "c" "a" "t" "c" "a" "a"
     [721] "c" "c" "c" "a" "g" "t" "a" "a" "t" "a" "a" "t" "g" "a" "t" "t" "t" "g"
     [739] "a" "a" "c" "a" "c" "c" "a" "c" "t" "g" "a" "g" "a" "a" "g" "c" "g" "t"
     [757] "g" "c" "a" "g" "c" "t" "g" "a" "g" "a" "g" "g" "c" "a" "t" "c" "c" "a"
     [775] "g" "a" "a" "a" "a" "g" "t" "a" "t" "c" "a" "g" "g" "g" "t" "a" "g" "t"
     [793] "t" "c" "t" "g" "t" "t" "t" "c" "a" "a" "a" "c" "t" "t" "g" "c" "a" "t"
     [811] "g" "t" "g" "g" "a" "g" "c" "c" "a" "t" "g" "t" "g" "g" "c" "a" "c" "a"
     [829] "a" "a" "t" "a" "c" "t" "c" "a" "t" "g" "c" "c" "a" "g" "c" "t" "c" "a"
     [847] "t" "t" "a" "c" "a" "g" "c" "a" "t" "g" "a" "g" "a" "a" "c" "a" "g" "c"
     [865] "a" "g" "t" "t" "t" "a" "t" "t" "a" "c" "t" "c" "a" "c" "t" "a" "a" "a"
     [883] "g" "a" "c" "a" "g" "a" "a" "t" "g" "a" "a" "t" "g" "t" "a" "g" "a" "a"
     [901] "a" "a" "g" "g" "c" "t" "g" "a" "a" "t" "t" "c" "t" "g" "t" "a" "a" "t"
     [919] "a" "a" "a" "a" "g" "c" "a" "a" "a" "c" "a" "g" "c" "c" "t" "g" "g" "c"
     [937] "t" "t" "a" "g" "c" "a" "a" "g" "g" "a" "g" "c" "c" "a" "a" "c" "a" "t"
     [955] "a" "a" "c" "a" "g" "a" "t" "g" "g" "g" "c" "t" "g" "g" "a" "a" "g" "t"
     [973] "a" "a" "g" "g" "a" "a" "a" "c" "a" "t" "g" "t" "a" "a" "t" "g" "a" "t"
     [991] "a" "g" "g" "c" "g" "g" "a" "c" "t" "c" "c" "c" "a" "g" "c" "a" "c" "a"
    [1009] "g" "a" "a" "a" "a" "a" "a" "a" "g" "g" "t" "a" "g" "a" "t" "c" "t" "g"
    [1027] "a" "a" "t" "g" "c" "t" "g" "a" "t" "c" "c" "c" "c" "t" "g" "t" "g" "t"
    [1045] "g" "a" "g" "a" "g" "a" "a" "a" "a" "g" "a" "a" "t" "g" "g" "a" "a" "t"
    [1063] "a" "a" "g" "c" "a" "g" "a" "a" "a" "c" "t" "g" "c" "c" "a" "t" "g" "c"
    [1081] "t" "c" "a" "g" "a" "g" "a" "a" "t" "c" "c" "t" "a" "g" "a" "g" "a" "t"
    [1099] "a" "c" "t" "g" "a" "a" "g" "a" "t" "g" "t" "t" "c" "c" "t" "t" "g" "g"
    [1117] "a" "t" "a" "a" "c" "a" "c" "t" "a" "a" "a" "t" "a" "g" "c" "a" "g" "c"
    [1135] "a" "t" "t" "c" "a" "g" "a" "a" "a" "g" "t" "t" "a" "a" "t" "g" "a" "g"
    [1153] "t" "g" "g" "t" "t" "t" "t" "c" "c" "a" "g" "a" "a" "g" "t" "g" "a" "t"
    [1171] "g" "a" "a" "c" "t" "g" "t" "t" "a" "g" "g" "t" "t" "c" "t" "g" "a" "t"
    [1189] "g" "a" "c" "t" "c" "a" "c" "a" "t" "g" "a" "t" "g" "g" "g" "g" "a" "g"
    [1207] "t" "c" "t" "g" "a" "a" "t" "c" "a" "a" "a" "t" "g" "c" "c" "a" "a" "a"
    [1225] "g" "t" "a" "g" "c" "t" "g" "a" "t" "g" "t" "a" "t" "t" "g" "g" "a" "c"
    [1243] "g" "t" "t" "c" "t" "a" "a" "a" "t" "g" "a" "g" "g" "t" "a" "g" "a" "t"
    [1261] "g" "a" "a" "t" "a" "t" "t" "c" "t" "g" "g" "t" "t" "c" "t" "t" "c" "a"
    [1279] "g" "a" "g" "a" "a" "a" "a" "t" "a" "g" "a" "c" "t" "t" "a" "c" "t" "g"
    [1297] "g" "c" "c" "a" "g" "t" "g" "a" "t" "c" "c" "t" "c" "a" "t" "g" "a" "g"
    [1315] "g" "c" "t" "t" "t" "a" "a" "t" "a" "t" "g" "t" "a" "a" "a" "a" "g" "t"
    [1333] "g" "a" "a" "a" "g" "a" "g" "t" "t" "c" "a" "c" "t" "c" "c" "a" "a" "a"
    [1351] "t" "c" "a" "g" "t" "a" "g" "a" "g" "a" "g" "t" "a" "a" "t" "a" "t" "t"
    [1369] "g" "a" "a" "g" "a" "c" "a" "a" "a" "a" "t" "a" "t" "t" "t" "g" "g" "g"
    [1387] "a" "a" "a" "a" "c" "c" "t" "a" "t" "c" "g" "g" "a" "a" "g" "a" "a" "g"
    [1405] "g" "c" "a" "a" "g" "c" "c" "t" "c" "c" "c" "c" "a" "a" "c" "t" "t" "a"
    [1423] "a" "g" "c" "c" "a" "t" "g" "t" "a" "a" "c" "t" "g" "a" "a" "a" "a" "t"
    [1441] "c" "t" "a" "a" "t" "t" "a" "t" "a" "g" "g" "a" "g" "c" "a" "t" "t" "t"
    [1459] "g" "t" "t" "a" "c" "t" "g" "a" "g" "c" "c" "a" "c" "a" "g" "a" "t" "a"
    [1477] "a" "t" "a" "c" "a" "a" "g" "a" "g" "c" "g" "t" "c" "c" "c" "c" "t" "c"
    [1495] "a" "c" "a" "a" "a" "t" "a" "a" "a" "t" "t" "a" "a" "a" "g" "c" "g" "t"
    [1513] "a" "a" "a" "a" "g" "g" "a" "g" "a" "c" "c" "t" "a" "c" "a" "t" "c" "a"
    [1531] "g" "g" "c" "c" "t" "t" "c" "a" "t" "c" "c" "t" "g" "a" "g" "g" "a" "t"
    [1549] "t" "t" "t" "a" "t" "c" "a" "a" "g" "a" "a" "a" "g" "c" "a" "g" "a" "t"
    [1567] "t" "t" "g" "g" "c" "a" "g" "t" "t" "c" "a" "a" "a" "a" "g" "a" "c" "t"
    [1585] "c" "c" "t" "g" "a" "a" "a" "t" "g" "a" "t" "a" "a" "a" "t" "c" "a" "g"
    [1603] "g" "g" "a" "a" "c" "t" "a" "a" "c" "c" "a" "a" "a" "c" "g" "g" "a" "g"
    [1621] "c" "a" "g" "a" "a" "t" "g" "g" "t" "c" "a" "a" "g" "t" "g" "a" "t" "g"
    [1639] "a" "a" "t" "a" "t" "t" "a" "c" "t" "a" "a" "t" "a" "g" "t" "g" "g" "t"
    [1657] "c" "a" "t" "g" "a" "g" "a" "a" "t" "a" "a" "a" "a" "c" "a" "a" "a" "a"
    [1675] "g" "g" "t" "g" "a" "t" "t" "c" "t" "a" "t" "t" "c" "a" "g" "a" "a" "t"
    [1693] "g" "a" "g" "a" "a" "a" "a" "a" "t" "c" "c" "t" "a" "a" "c" "c" "c" "a"
    [1711] "a" "t" "a" "g" "a" "a" "t" "c" "a" "c" "t" "c" "g" "a" "a" "a" "a" "a"
    [1729] "g" "a" "a" "t" "c" "t" "g" "c" "t" "t" "t" "c" "a" "a" "a" "a" "c" "g"
    [1747] "a" "a" "a" "g" "c" "t" "g" "a" "a" "c" "c" "t" "a" "t" "a" "a" "g" "c"
    [1765] "a" "g" "c" "a" "g" "t" "a" "t" "a" "a" "g" "c" "a" "a" "t" "a" "t" "g"
    [1783] "g" "a" "a" "c" "t" "c" "g" "a" "a" "t" "t" "a" "a" "a" "t" "a" "t" "c"
    [1801] "c" "a" "c" "a" "a" "t" "t" "c" "a" "a" "a" "a" "g" "c" "a" "c" "c" "t"
    [1819] "a" "a" "a" "a" "a" "g" "a" "a" "t" "a" "g" "g" "c" "t" "g" "a" "g" "g"
    [1837] "a" "g" "g" "a" "a" "g" "t" "c" "t" "t" "c" "t" "a" "c" "c" "a" "g" "g"
    [1855] "c" "a" "t" "a" "t" "t" "c" "a" "t" "g" "c" "g" "c" "t" "t" "g" "a" "a"
    [1873] "c" "t" "a" "g" "t" "a" "g" "t" "c" "a" "g" "t" "a" "g" "a" "a" "a" "t"
    [1891] "c" "t" "a" "a" "g" "c" "c" "c" "a" "c" "c" "t" "a" "a" "t" "t" "g" "t"
    [1909] "a" "c" "t" "g" "a" "a" "t" "t" "g" "c" "a" "a" "a" "t" "t" "g" "a" "t"
    [1927] "a" "g" "t" "t" "g" "t" "t" "c" "t" "a" "g" "c" "a" "g" "t" "g" "a" "a"
    [1945] "g" "a" "g" "a" "t" "a" "a" "a" "g" "a" "a" "a" "a" "a" "a" "a" "a" "g"
    [1963] "t" "a" "c" "a" "a" "c" "c" "a" "a" "a" "t" "g" "c" "c" "a" "g" "t" "c"
    [1981] "a" "g" "g" "c" "a" "c" "a" "g" "c" "a" "g" "a" "a" "a" "c" "c" "t" "a"
    [1999] "c" "a" "a" "c" "t" "c" "a" "t" "g" "g" "a" "a" "g" "g" "t" "a" "a" "a"
    [2017] "g" "a" "a" "c" "c" "t" "g" "c" "a" "a" "c" "t" "g" "g" "a" "g" "c" "c"
    [2035] "a" "a" "g" "a" "a" "g" "a" "g" "t" "a" "a" "c" "a" "a" "g" "c" "c" "a"
    [2053] "a" "a" "t" "g" "a" "a" "c" "a" "g" "a" "c" "a" "a" "g" "t" "a" "a" "a"
    [2071] "a" "g" "a" "c" "a" "t" "g" "a" "c" "a" "g" "c" "g" "a" "t" "a" "c" "t"
    [2089] "t" "t" "c" "c" "c" "a" "g" "a" "g" "c" "t" "g" "a" "a" "g" "t" "t" "a"
    [2107] "a" "c" "a" "a" "a" "t" "g" "c" "a" "c" "c" "t" "g" "g" "t" "t" "c" "t"
    [2125] "t" "t" "t" "a" "c" "t" "a" "a" "g" "t" "g" "t" "t" "c" "a" "a" "a" "t"
    [2143] "a" "c" "c" "a" "g" "t" "g" "a" "a" "c" "t" "t" "a" "a" "a" "g" "a" "a"
    [2161] "t" "t" "t" "g" "t" "c" "a" "a" "t" "c" "c" "t" "a" "g" "c" "c" "t" "t"
    [2179] "c" "c" "a" "a" "g" "a" "g" "a" "a" "g" "a" "a" "a" "a" "a" "g" "a" "a"
    [2197] "g" "a" "g" "a" "a" "a" "c" "t" "a" "g" "a" "a" "a" "c" "a" "g" "t" "t"
    [2215] "a" "a" "a" "g" "t" "g" "t" "c" "t" "a" "a" "t" "a" "a" "t" "g" "c" "t"
    [2233] "g" "a" "a" "g" "a" "c" "c" "c" "c" "a" "a" "a" "g" "a" "t" "c" "t" "c"
    [2251] "a" "t" "g" "t" "t" "a" "a" "g" "t" "g" "g" "a" "g" "a" "a" "a" "g" "g"
    [2269] "g" "t" "t" "t" "t" "g" "c" "a" "a" "a" "c" "t" "g" "a" "a" "a" "g" "a"
    [2287] "t" "c" "t" "g" "t" "a" "g" "a" "g" "a" "g" "t" "a" "g" "c" "a" "g" "t"
    [2305] "a" "t" "t" "t" "c" "a" "t" "t" "g" "g" "t" "a" "c" "c" "t" "g" "g" "t"
    [2323] "a" "c" "t" "g" "a" "t" "t" "a" "t" "g" "g" "c" "a" "c" "t" "c" "a" "g"
    [2341] "g" "a" "a" "a" "g" "t" "a" "t" "c" "t" "c" "g" "t" "t" "a" "c" "t" "g"
    [2359] "g" "a" "a" "g" "t" "t" "a" "g" "c" "a" "c" "t" "c" "t" "a" "g" "g" "g"
    [2377] "a" "a" "g" "g" "c" "a" "a" "a" "a" "a" "c" "a" "g" "a" "a" "c" "c" "a"
    [2395] "a" "a" "t" "a" "a" "a" "t" "g" "t" "g" "t" "g" "a" "g" "t" "c" "a" "g"
    [2413] "t" "g" "t" "g" "c" "a" "g" "c" "a" "t" "t" "t" "g" "a" "a" "a" "a" "c"
    [2431] "c" "c" "c" "a" "a" "g" "g" "g" "a" "c" "t" "a" "a" "t" "t" "c" "a" "t"
    [2449] "g" "g" "t" "t" "g" "t" "t" "c" "c" "a" "a" "a" "g" "a" "t" "a" "a" "t"
    [2467] "a" "g" "a" "a" "a" "t" "g" "a" "c" "a" "c" "a" "g" "a" "a" "g" "g" "c"
    [2485] "t" "t" "t" "a" "a" "g" "t" "a" "t" "c" "c" "a" "t" "t" "g" "g" "g" "a"
    [2503] "c" "a" "t" "g" "a" "a" "g" "t" "t" "a" "a" "c" "c" "a" "c" "a" "g" "t"
    [2521] "c" "g" "g" "g" "a" "a" "a" "c" "a" "a" "g" "c" "a" "t" "a" "g" "a" "a"
    [2539] "a" "t" "g" "g" "a" "a" "g" "a" "a" "a" "g" "t" "g" "a" "a" "c" "t" "t"
    [2557] "g" "a" "t" "g" "c" "t" "c" "a" "g" "t" "a" "t" "t" "t" "g" "c" "a" "g"
    [2575] "a" "a" "t" "a" "c" "a" "t" "t" "c" "a" "a" "g" "g" "t" "t" "t" "c" "a"
    [2593] "a" "a" "g" "c" "g" "c" "c" "a" "g" "t" "c" "a" "t" "t" "t" "g" "c" "t"
    [2611] "c" "c" "g" "t" "t" "t" "t" "c" "a" "a" "a" "t" "c" "c" "a" "g" "g" "a"
    [2629] "a" "a" "t" "g" "c" "a" "g" "a" "a" "g" "a" "g" "g" "a" "a" "t" "g" "t"
    [2647] "g" "c" "a" "a" "c" "a" "t" "t" "c" "t" "c" "t" "g" "c" "c" "c" "a" "c"
    [2665] "t" "c" "t" "g" "g" "g" "t" "c" "c" "t" "t" "a" "a" "a" "g" "a" "a" "a"
    [2683] "c" "a" "a" "a" "g" "t" "c" "c" "a" "a" "a" "a" "g" "t" "c" "a" "c" "t"
    [2701] "t" "t" "t" "g" "a" "a" "t" "g" "t" "g" "a" "a" "c" "a" "a" "a" "a" "g"
    [2719] "g" "a" "a" "g" "a" "a" "a" "a" "t" "c" "a" "a" "g" "g" "a" "a" "a" "g"
    [2737] "a" "a" "t" "g" "a" "g" "t" "c" "t" "a" "a" "t" "a" "t" "c" "a" "a" "g"
    [2755] "c" "c" "t" "g" "t" "a" "c" "a" "g" "a" "c" "a" "g" "t" "t" "a" "a" "t"
    [2773] "a" "t" "c" "a" "c" "t" "g" "c" "a" "g" "g" "c" "t" "t" "t" "c" "c" "t"
    [2791] "g" "t" "g" "g" "t" "t" "g" "g" "t" "c" "a" "g" "a" "a" "a" "g" "a" "t"
    [2809] "a" "a" "g" "c" "c" "a" "g" "t" "t" "g" "a" "t" "a" "a" "t" "g" "c" "c"
    [2827] "a" "a" "a" "t" "g" "t" "a" "g" "t" "a" "t" "c" "a" "a" "a" "g" "g" "a"
    [2845] "g" "g" "c" "t" "c" "t" "a" "g" "g" "t" "t" "t" "t" "g" "t" "c" "t" "a"
    [2863] "t" "c" "a" "t" "c" "t" "c" "a" "g" "t" "t" "c" "a" "g" "a" "g" "g" "c"
    [2881] "a" "a" "c" "g" "a" "a" "a" "c" "t" "g" "g" "a" "c" "t" "c" "a" "t" "t"
    [2899] "a" "c" "t" "c" "c" "a" "a" "a" "t" "a" "a" "a" "c" "a" "t" "g" "g" "a"
    [2917] "c" "t" "t" "t" "t" "a" "c" "a" "a" "a" "a" "c" "c" "c" "a" "t" "a" "t"
    [2935] "c" "g" "t" "a" "t" "a" "c" "c" "a" "c" "c" "a" "c" "t" "t" "t" "t" "t"
    [2953] "c" "c" "c" "a" "t" "c" "a" "a" "g" "t" "c" "a" "t" "t" "t" "g" "t" "t"
    [2971] "a" "a" "a" "a" "c" "t" "a" "a" "a" "t" "g" "t" "a" "a" "g" "a" "a" "a"
    [2989] "a" "a" "t" "c" "t" "g" "c" "t" "a" "g" "a" "g" "g" "a" "a" "a" "a" "c"
    [3007] "t" "t" "t" "g" "a" "g" "g" "a" "a" "c" "a" "t" "t" "c" "a" "a" "t" "g"
    [3025] "t" "c" "a" "c" "c" "t" "g" "a" "a" "a" "g" "a" "g" "a" "a" "a" "t" "g"
    [3043] "g" "g" "a" "a" "a" "t" "g" "a" "g" "a" "a" "c" "a" "t" "t" "c" "c" "a"
    [3061] "a" "g" "t" "a" "c" "a" "g" "t" "g" "a" "g" "c" "a" "c" "a" "a" "t" "t"
    [3079] "a" "g" "c" "c" "g" "t" "a" "a" "t" "a" "a" "c" "a" "t" "t" "a" "g" "a"
    [3097] "g" "a" "a" "a" "a" "t" "g" "t" "t" "t" "t" "t" "a" "a" "a" "g" "a" "a"
    [3115] "g" "c" "c" "a" "g" "c" "t" "c" "a" "a" "g" "c" "a" "a" "t" "a" "t" "t"
    [3133] "a" "a" "t" "g" "a" "a" "g" "t" "a" "g" "g" "t" "t" "c" "c" "a" "g" "t"
    [3151] "a" "c" "t" "a" "a" "t" "g" "a" "a" "g" "t" "g" "g" "g" "c" "t" "c" "c"
    [3169] "a" "g" "t" "a" "t" "t" "a" "a" "t" "g" "a" "a" "a" "t" "a" "g" "g" "t"
    [3187] "t" "c" "c" "a" "g" "t" "g" "a" "t" "g" "a" "a" "a" "a" "c" "a" "t" "t"
    [3205] "c" "a" "a" "g" "c" "a" "g" "a" "a" "c" "t" "a" "g" "g" "t" "a" "g" "a"
    [3223] "a" "a" "c" "a" "g" "a" "g" "g" "g" "c" "c" "a" "a" "a" "a" "t" "t" "g"
    [3241] "a" "a" "t" "g" "c" "t" "a" "t" "g" "c" "t" "t" "a" "g" "a" "t" "t" "a"
    [3259] "g" "g" "g" "g" "t" "t" "t" "t" "g" "c" "a" "a" "c" "c" "t" "g" "a" "g"
    [3277] "g" "t" "c" "t" "a" "t" "a" "a" "a" "c" "a" "a" "a" "g" "t" "c" "t" "t"
    [3295] "c" "c" "t" "g" "g" "a" "a" "g" "t" "a" "a" "t" "t" "g" "t" "a" "a" "g"
    [3313] "c" "a" "t" "c" "c" "t" "g" "a" "a" "a" "t" "a" "a" "a" "a" "a" "a" "g"
    [3331] "c" "a" "a" "g" "a" "a" "t" "a" "t" "g" "a" "a" "g" "a" "a" "g" "t" "a"
    [3349] "g" "t" "t" "c" "a" "g" "a" "c" "t" "g" "t" "t" "a" "a" "t" "a" "c" "a"
    [3367] "g" "a" "t" "t" "t" "c" "t" "c" "t" "c" "c" "a" "t" "a" "t" "c" "t" "g"
    [3385] "a" "t" "t" "t" "c" "a" "g" "a" "t" "a" "a" "c" "t" "t" "a" "g" "a" "a"
    [3403] "c" "a" "g" "c" "c" "t" "a" "t" "g" "g" "g" "a" "a" "g" "t" "a" "g" "t"
    [3421] "c" "a" "t" "g" "c" "a" "t" "c" "t" "c" "a" "g" "g" "t" "t" "t" "g" "t"
    [3439] "t" "c" "t" "g" "a" "g" "a" "c" "a" "c" "c" "t" "g" "a" "t" "g" "a" "c"
    [3457] "c" "t" "g" "t" "t" "a" "g" "a" "t" "g" "a" "t" "g" "g" "t" "g" "a" "a"
    [3475] "a" "t" "a" "a" "a" "g" "g" "a" "a" "g" "a" "t" "a" "c" "t" "a" "g" "t"
    [3493] "t" "t" "t" "g" "c" "t" "g" "a" "a" "a" "a" "t" "g" "a" "c" "a" "t" "t"
    [3511] "a" "a" "g" "g" "a" "a" "a" "g" "t" "t" "c" "t" "g" "c" "t" "g" "t" "t"
    [3529] "t" "t" "t" "a" "g" "c" "a" "a" "a" "a" "g" "c" "g" "t" "c" "c" "a" "g"
    [3547] "a" "a" "a" "g" "g" "a" "g" "a" "g" "c" "t" "t" "a" "g" "c" "a" "g" "g"
    [3565] "a" "g" "t" "c" "c" "t" "a" "g" "c" "c" "c" "t" "t" "t" "c" "a" "c" "c"
    [3583] "c" "a" "t" "a" "c" "a" "c" "a" "t" "t" "t" "g" "g" "c" "t" "c" "a" "g"
    [3601] "g" "g" "t" "t" "a" "c" "c" "g" "a" "a" "g" "a" "g" "g" "g" "g" "c" "c"
    [3619] "a" "a" "g" "a" "a" "a" "t" "t" "a" "g" "a" "g" "t" "c" "c" "t" "c" "a"
    [3637] "g" "a" "a" "g" "a" "g" "a" "a" "c" "t" "t" "a" "t" "c" "t" "a" "g" "t"
    [3655] "g" "a" "g" "g" "a" "t" "g" "a" "a" "g" "a" "g" "c" "t" "t" "c" "c" "c"
    [3673] "t" "g" "c" "t" "t" "c" "c" "a" "a" "c" "a" "c" "t" "t" "g" "t" "t" "a"
    [3691] "t" "t" "t" "g" "g" "t" "a" "a" "a" "g" "t" "a" "a" "a" "c" "a" "a" "t"
    [3709] "a" "t" "a" "c" "c" "t" "t" "c" "t" "c" "a" "g" "t" "c" "t" "a" "c" "t"
    [3727] "a" "g" "g" "c" "a" "t" "a" "g" "c" "a" "c" "c" "g" "t" "t" "g" "c" "t"
    [3745] "a" "c" "c" "g" "a" "g" "t" "g" "t" "c" "t" "g" "t" "c" "t" "a" "a" "g"
    [3763] "a" "a" "c" "a" "c" "a" "g" "a" "g" "g" "a" "g" "a" "a" "t" "t" "t" "a"
    [3781] "t" "t" "a" "t" "c" "a" "t" "t" "g" "a" "a" "g" "a" "a" "t" "a" "g" "c"
    [3799] "t" "t" "a" "a" "a" "t" "g" "a" "c" "t" "g" "c" "a" "g" "t" "a" "a" "c"
    [3817] "c" "a" "g" "g" "t" "a" "a" "t" "a" "t" "t" "g" "g" "c" "a" "a" "a" "g"
    [3835] "g" "c" "a" "t" "c" "t" "c" "a" "g" "g" "a" "a" "c" "a" "t" "c" "a" "c"
    [3853] "c" "t" "t" "a" "g" "t" "g" "a" "g" "g" "a" "a" "a" "c" "a" "a" "a" "a"
    [3871] "t" "g" "t" "t" "c" "t" "g" "c" "t" "a" "g" "c" "t" "t" "g" "t" "t" "t"
    [3889] "t" "c" "t" "t" "c" "a" "c" "a" "g" "t" "g" "c" "a" "g" "t" "g" "a" "a"
    [3907] "t" "t" "g" "g" "a" "a" "g" "a" "c" "t" "t" "g" "a" "c" "t" "g" "c" "a"
    [3925] "a" "a" "t" "a" "c" "a" "a" "a" "c" "a" "c" "c" "c" "a" "g" "g" "a" "t"
    [3943] "c" "c" "t" "t" "t" "c" "t" "t" "g" "a" "t" "t" "g" "g" "t" "t" "c" "t"
    [3961] "t" "c" "c" "a" "a" "a" "c" "a" "a" "a" "t" "g" "a" "g" "g" "c" "a" "t"
    [3979] "c" "a" "g" "t" "c" "t" "g" "a" "a" "a" "g" "c" "c" "a" "g" "g" "g" "a"
    [3997] "g" "t" "t" "g" "g" "t" "c" "t" "g" "a" "g" "t" "g" "a" "c" "a" "a" "g"
    [4015] "g" "a" "a" "t" "t" "g" "g" "t" "t" "t" "c" "a" "g" "a" "t" "g" "a" "t"
    [4033] "g" "a" "a" "g" "a" "a" "a" "g" "a" "g" "g" "a" "a" "c" "g" "g" "g" "c"
    [4051] "t" "t" "g" "g" "a" "a" "g" "a" "a" "a" "a" "t" "a" "a" "t" "c" "a" "a"
    [4069] "g" "a" "a" "g" "a" "g" "c" "a" "a" "a" "g" "c" "a" "t" "g" "g" "a" "t"
    [4087] "t" "c" "a" "a" "a" "c" "t" "t" "a" "g" "g" "t" "g" "a" "a" "g" "c" "a"
    [4105] "g" "c" "a" "t" "c" "t" "g" "g" "g" "t" "g" "t" "g" "a" "g" "a" "g" "t"
    [4123] "g" "a" "a" "a" "c" "a" "a" "g" "c" "g" "t" "c" "t" "c" "t" "g" "a" "a"
    [4141] "g" "a" "c" "t" "g" "c" "t" "c" "a" "g" "g" "g" "c" "t" "a" "t" "c" "c"
    [4159] "t" "c" "t" "c" "a" "g" "a" "g" "t" "g" "a" "c" "a" "t" "t" "t" "t" "a"
    [4177] "a" "c" "c" "a" "c" "t" "c" "a" "g" "c" "a" "g" "a" "g" "g" "g" "a" "t"
    [4195] "a" "c" "c" "a" "t" "g" "c" "a" "a" "c" "a" "t" "a" "a" "c" "c" "t" "g"
    [4213] "a" "t" "a" "a" "a" "g" "c" "t" "c" "c" "a" "g" "c" "a" "g" "g" "a" "a"
    [4231] "a" "t" "g" "g" "c" "t" "g" "a" "a" "c" "t" "a" "g" "a" "a" "g" "c" "t"
    [4249] "g" "t" "g" "t" "t" "a" "g" "a" "a" "c" "a" "g" "c" "a" "t" "g" "g" "g"
    [4267] "a" "g" "c" "c" "a" "g" "c" "c" "t" "t" "c" "t" "a" "a" "c" "a" "g" "c"
    [4285] "t" "a" "c" "c" "c" "t" "t" "c" "c" "a" "t" "c" "a" "t" "a" "a" "g" "t"
    [4303] "g" "a" "c" "t" "c" "t" "t" "c" "t" "g" "c" "c" "c" "t" "t" "g" "a" "g"
    [4321] "g" "a" "c" "c" "t" "g" "c" "g" "a" "a" "a" "t" "c" "c" "a" "g" "a" "a"
    [4339] "c" "a" "a" "a" "g" "c" "a" "c" "a" "t" "c" "a" "g" "a" "a" "a" "a" "a"
    [4357] "g" "c" "a" "g" "t" "a" "t" "t" "a" "a" "c" "t" "t" "c" "a" "c" "a" "g"
    [4375] "a" "a" "a" "a" "g" "t" "a" "g" "t" "g" "a" "a" "t" "a" "c" "c" "c" "t"
    [4393] "a" "t" "a" "a" "g" "c" "c" "a" "g" "a" "a" "t" "c" "c" "a" "g" "a" "a"
    [4411] "g" "g" "c" "c" "t" "t" "t" "c" "t" "g" "c" "t" "g" "a" "c" "a" "a" "g"
    [4429] "t" "t" "t" "g" "a" "g" "g" "t" "g" "t" "c" "t" "g" "c" "a" "g" "a" "t"
    [4447] "a" "g" "t" "t" "c" "t" "a" "c" "c" "a" "g" "t" "a" "a" "a" "a" "a" "t"
    [4465] "a" "a" "a" "g" "a" "a" "c" "c" "a" "g" "g" "a" "g" "t" "g" "g" "a" "a"
    [4483] "a" "g" "g" "t" "c" "a" "t" "c" "c" "c" "c" "t" "t" "c" "t" "a" "a" "a"
    [4501] "t" "g" "c" "c" "c" "a" "t" "c" "a" "t" "t" "a" "g" "a" "t" "g" "a" "t"
    [4519] "a" "g" "g" "t" "g" "g" "t" "a" "c" "a" "t" "g" "c" "a" "c" "a" "g" "t"
    [4537] "t" "g" "c" "t" "c" "t" "g" "g" "g" "a" "g" "t" "c" "t" "t" "c" "a" "g"
    [4555] "a" "a" "t" "a" "g" "a" "a" "a" "c" "t" "a" "c" "c" "c" "a" "t" "c" "t"
    [4573] "c" "a" "a" "g" "a" "g" "g" "a" "g" "c" "t" "c" "a" "t" "t" "a" "a" "g"
    [4591] "g" "t" "t" "g" "t" "t" "g" "a" "t" "g" "t" "g" "g" "a" "g" "g" "a" "g"
    [4609] "c" "a" "a" "c" "a" "g" "c" "t" "g" "g" "a" "a" "g" "a" "g" "t" "c" "t"
    [4627] "g" "g" "g" "c" "c" "a" "c" "a" "c" "g" "a" "t" "t" "t" "g" "a" "c" "g"
    [4645] "g" "a" "a" "a" "c" "a" "t" "c" "t" "t" "a" "c" "t" "t" "g" "c" "c" "a"
    [4663] "a" "g" "g" "c" "a" "a" "g" "a" "t" "c" "t" "a" "g" "a" "g" "g" "g" "a"
    [4681] "a" "c" "c" "c" "c" "t" "t" "a" "c" "c" "t" "g" "g" "a" "a" "t" "c" "t"
    [4699] "g" "g" "a" "a" "t" "c" "a" "g" "c" "c" "t" "c" "t" "t" "c" "t" "c" "t"
    [4717] "g" "a" "t" "g" "a" "c" "c" "c" "t" "g" "a" "a" "t" "c" "t" "g" "a" "t"
    [4735] "c" "c" "t" "t" "c" "t" "g" "a" "a" "g" "a" "c" "a" "g" "a" "g" "c" "c"
    [4753] "c" "c" "a" "g" "a" "g" "t" "c" "a" "g" "c" "t" "c" "g" "t" "g" "t" "t"
    [4771] "g" "g" "c" "a" "a" "c" "a" "t" "a" "c" "c" "a" "t" "c" "t" "t" "c" "a"
    [4789] "a" "c" "c" "t" "c" "t" "g" "c" "a" "t" "t" "g" "a" "a" "a" "g" "t" "t"
    [4807] "c" "c" "c" "c" "a" "a" "t" "t" "g" "a" "a" "a" "g" "t" "t" "g" "c" "a"
    [4825] "g" "a" "a" "t" "c" "t" "g" "c" "c" "c" "a" "g" "a" "g" "t" "c" "c" "a"
    [4843] "g" "c" "t" "g" "c" "t" "g" "c" "t" "c" "a" "t" "a" "c" "t" "a" "c" "t"
    [4861] "g" "a" "t" "a" "c" "t" "g" "c" "t" "g" "g" "g" "t" "a" "t" "a" "a" "t"
    [4879] "g" "c" "a" "a" "t" "g" "g" "a" "a" "g" "a" "a" "a" "g" "t" "g" "t" "g"
    [4897] "a" "g" "c" "a" "g" "g" "g" "a" "g" "a" "a" "g" "c" "c" "a" "g" "a" "a"
    [4915] "t" "t" "g" "a" "c" "a" "g" "c" "t" "t" "c" "a" "a" "c" "a" "g" "a" "a"
    [4933] "a" "g" "g" "g" "t" "c" "a" "a" "c" "a" "a" "a" "a" "g" "a" "a" "t" "g"
    [4951] "t" "c" "c" "a" "t" "g" "g" "t" "g" "g" "t" "g" "t" "c" "t" "g" "g" "c"
    [4969] "c" "t" "g" "a" "c" "c" "c" "c" "a" "g" "a" "a" "g" "a" "a" "t" "t" "t"
    [4987] "a" "t" "g" "c" "t" "c" "g" "t" "g" "t" "a" "c" "a" "a" "g" "t" "t" "t"
    [5005] "g" "c" "c" "a" "g" "a" "a" "a" "a" "c" "a" "c" "c" "a" "c" "a" "t" "c"
    [5023] "a" "c" "t" "t" "t" "a" "a" "c" "t" "a" "a" "t" "c" "t" "a" "a" "t" "t"
    [5041] "a" "c" "t" "g" "a" "a" "g" "a" "g" "a" "c" "t" "a" "c" "t" "c" "a" "t"
    [5059] "g" "t" "t" "g" "t" "t" "a" "t" "g" "a" "a" "a" "a" "c" "a" "g" "a" "t"
    [5077] "g" "c" "t" "g" "a" "g" "t" "t" "t" "g" "t" "g" "t" "g" "t" "g" "a" "a"
    [5095] "c" "g" "g" "a" "c" "a" "c" "t" "g" "a" "a" "a" "t" "a" "t" "t" "t" "t"
    [5113] "c" "t" "a" "g" "g" "a" "a" "t" "t" "g" "c" "g" "g" "g" "a" "g" "g" "a"
    [5131] "a" "a" "a" "t" "g" "g" "g" "t" "a" "g" "t" "t" "a" "g" "c" "t" "a" "t"
    [5149] "t" "t" "c" "t" "g" "g" "g" "t" "g" "a" "c" "c" "c" "a" "g" "t" "c" "t"
    [5167] "a" "t" "t" "a" "a" "a" "g" "a" "a" "a" "g" "a" "a" "a" "a" "a" "t" "g"
    [5185] "c" "t" "g" "a" "a" "t" "g" "a" "g" "c" "a" "t" "g" "a" "t" "t" "t" "t"
    [5203] "g" "a" "a" "g" "t" "c" "a" "g" "a" "g" "g" "a" "g" "a" "t" "g" "t" "g"
    [5221] "g" "t" "c" "a" "a" "t" "g" "g" "a" "a" "g" "a" "a" "a" "c" "c" "a" "c"
    [5239] "c" "a" "a" "g" "g" "t" "c" "c" "a" "a" "a" "g" "c" "g" "a" "g" "c" "a"
    [5257] "a" "g" "a" "g" "a" "a" "t" "c" "c" "c" "a" "g" "g" "a" "c" "a" "g" "a"
    [5275] "a" "a" "g" "a" "t" "c" "t" "t" "c" "a" "g" "g" "g" "g" "g" "c" "t" "a"
    [5293] "g" "a" "a" "a" "t" "c" "t" "g" "t" "t" "g" "c" "t" "a" "t" "g" "g" "g"
    [5311] "c" "c" "c" "t" "t" "c" "a" "c" "c" "a" "a" "c" "a" "t" "g" "c" "c" "c"
    [5329] "a" "c" "a" "g" "a" "t" "c" "a" "a" "c" "t" "g" "g" "a" "a" "t" "g" "g"
    [5347] "a" "t" "g" "g" "t" "a" "c" "a" "g" "c" "t" "g" "t" "g" "t" "g" "g" "t"
    [5365] "g" "c" "t" "t" "c" "t" "g" "t" "g" "g" "t" "g" "a" "a" "g" "g" "a" "g"
    [5383] "c" "t" "t" "t" "c" "a" "t" "c" "a" "t" "t" "c" "a" "c" "c" "c" "t" "t"
    [5401] "g" "g" "c" "a" "c" "a" "g" "g" "t" "g" "t" "c" "c" "a" "c" "c" "c" "a"
    [5419] "a" "t" "t" "g" "t" "g" "g" "t" "t" "g" "t" "g" "c" "a" "g" "c" "c" "a"
    [5437] "g" "a" "t" "g" "c" "c" "t" "g" "g" "a" "c" "a" "g" "a" "g" "g" "a" "c"
    [5455] "a" "a" "t" "g" "g" "c" "t" "t" "c" "c" "a" "t" "g" "c" "a" "a" "t" "t"
    [5473] "g" "g" "g" "c" "a" "g" "a" "t" "g" "t" "g" "t" "g" "a" "g" "g" "c" "a"
    [5491] "c" "c" "t" "g" "t" "g" "g" "t" "g" "a" "c" "c" "c" "g" "a" "g" "a" "g"
    [5509] "t" "g" "g" "g" "t" "g" "t" "t" "g" "g" "a" "c" "a" "g" "t" "g" "t" "a"
    [5527] "g" "c" "a" "c" "t" "c" "t" "a" "c" "c" "a" "g" "t" "g" "c" "c" "a" "g"
    [5545] "g" "a" "g" "c" "t" "g" "g" "a" "c" "a" "c" "c" "t" "a" "c" "c" "t" "g"
    [5563] "a" "t" "a" "c" "c" "c" "c" "a" "g" "a" "t" "c" "c" "c" "c" "c" "a" "c"
    [5581] "a" "g" "c" "c" "a" "c" "t" "a" "c" "t" "g" "a"
    attr(,"name")
    [1] "DQ190450.BRCA1"
    attr(,"Annot")
    [1] ">DQ190450.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $DQ190451.BRCA1
       [1] "a" "t" "g" "g" "a" "t" "t" "t" "a" "t" "c" "t" "g" "c" "t" "c" "t" "t"
      [19] "c" "g" "c" "g" "t" "t" "g" "a" "a" "g" "a" "a" "g" "t" "a" "c" "a" "a"
      [37] "a" "a" "t" "g" "t" "c" "a" "t" "t" "a" "a" "t" "g" "c" "t" "a" "t" "g"
      [55] "c" "a" "g" "a" "a" "a" "a" "t" "c" "t" "t" "a" "g" "a" "g" "t" "g" "t"
      [73] "c" "c" "c" "a" "t" "c" "t" "g" "t" "c" "t" "g" "g" "a" "g" "t" "t" "g"
      [91] "a" "t" "c" "a" "a" "g" "g" "a" "a" "c" "c" "t" "g" "t" "c" "t" "c" "c"
     [109] "a" "c" "a" "a" "a" "g" "t" "g" "t" "g" "a" "c" "c" "a" "c" "a" "t" "a"
     [127] "t" "t" "t" "t" "g" "c" "a" "a" "a" "t" "t" "t" "t" "g" "c" "a" "t" "g"
     [145] "c" "t" "g" "a" "a" "a" "c" "t" "t" "c" "t" "c" "a" "a" "c" "c" "a" "g"
     [163] "a" "a" "g" "a" "a" "a" "g" "g" "g" "c" "c" "t" "t" "c" "a" "c" "a" "g"
     [181] "t" "g" "t" "c" "c" "t" "t" "t" "a" "t" "g" "t" "a" "a" "g" "a" "a" "t"
     [199] "g" "a" "t" "a" "t" "a" "a" "c" "c" "a" "a" "a" "a" "g" "g" "a" "g" "c"
     [217] "c" "t" "a" "c" "a" "a" "g" "a" "a" "a" "g" "t" "a" "c" "g" "a" "g" "a"
     [235] "t" "t" "t" "a" "g" "t" "c" "a" "a" "c" "t" "t" "g" "t" "t" "g" "a" "a"
     [253] "g" "a" "g" "c" "t" "a" "t" "t" "g" "a" "a" "a" "a" "t" "c" "a" "t" "t"
     [271] "t" "g" "t" "g" "c" "t" "t" "t" "t" "c" "a" "g" "c" "t" "t" "g" "a" "c"
     [289] "a" "c" "a" "g" "g" "t" "t" "t" "g" "g" "a" "g" "t" "a" "t" "g" "c" "a"
     [307] "a" "a" "c" "a" "g" "c" "t" "a" "t" "a" "a" "t" "t" "t" "t" "g" "c" "a"
     [325] "a" "a" "a" "a" "a" "g" "g" "a" "a" "a" "a" "t" "a" "a" "c" "t" "c" "t"
     [343] "c" "c" "t" "g" "a" "a" "c" "a" "t" "c" "t" "a" "a" "a" "a" "g" "a" "t"
     [361] "g" "a" "a" "g" "t" "t" "t" "c" "t" "a" "t" "c" "a" "t" "c" "c" "a" "a"
     [379] "a" "g" "t" "a" "t" "g" "g" "g" "c" "t" "a" "c" "a" "g" "a" "a" "a" "c"
     [397] "c" "g" "t" "g" "c" "c" "a" "a" "a" "a" "g" "a" "c" "t" "t" "c" "t" "a"
     [415] "c" "a" "g" "a" "g" "t" "g" "a" "a" "c" "c" "c" "g" "a" "a" "a" "a" "t"
     [433] "c" "c" "t" "t" "c" "c" "t" "t" "g" "c" "a" "g" "g" "a" "a" "a" "c" "c"
     [451] "a" "g" "t" "c" "t" "c" "a" "g" "t" "g" "t" "c" "c" "a" "a" "c" "t" "c"
     [469] "t" "c" "t" "a" "a" "c" "c" "t" "t" "g" "g" "a" "a" "c" "t" "g" "t" "g"
     [487] "a" "g" "a" "a" "c" "t" "c" "t" "g" "a" "g" "g" "a" "c" "a" "a" "a" "g"
     [505] "c" "a" "g" "c" "g" "g" "a" "t" "a" "c" "a" "a" "c" "c" "t" "c" "a" "a"
     [523] "a" "a" "g" "a" "c" "g" "t" "c" "t" "g" "t" "c" "t" "a" "c" "a" "t" "t"
     [541] "g" "a" "a" "t" "t" "g" "g" "g" "a" "t" "c" "t" "g" "a" "t" "t" "c" "t"
     [559] "t" "c" "t" "g" "a" "a" "g" "a" "t" "a" "c" "c" "g" "t" "t" "a" "a" "t"
     [577] "a" "a" "g" "g" "c" "a" "a" "c" "t" "t" "a" "t" "t" "g" "c" "a" "g" "t"
     [595] "g" "t" "g" "g" "g" "a" "g" "a" "t" "c" "a" "a" "g" "a" "a" "t" "t" "g"
     [613] "t" "t" "a" "c" "a" "a" "a" "t" "c" "a" "c" "c" "c" "c" "t" "c" "a" "a"
     [631] "g" "g" "a" "a" "c" "c" "a" "g" "g" "g" "a" "t" "g" "a" "a" "a" "t" "c"
     [649] "a" "g" "t" "t" "t" "g" "g" "a" "t" "t" "c" "t" "g" "c" "a" "a" "a" "a"
     [667] "a" "a" "g" "g" "c" "t" "g" "c" "t" "t" "g" "t" "g" "a" "a" "t" "t" "t"
     [685] "t" "c" "t" "g" "a" "g" "a" "c" "g" "g" "a" "t" "g" "t" "a" "a" "c" "a"
     [703] "a" "a" "t" "a" "c" "t" "g" "a" "a" "c" "a" "t" "c" "a" "t" "c" "a" "a"
     [721] "c" "c" "c" "a" "g" "t" "a" "a" "t" "a" "a" "t" "g" "a" "t" "t" "t" "g"
     [739] "a" "a" "c" "a" "c" "c" "a" "c" "t" "g" "a" "g" "a" "a" "g" "c" "g" "t"
     [757] "g" "c" "a" "g" "c" "t" "g" "a" "g" "a" "g" "g" "c" "a" "t" "c" "c" "a"
     [775] "g" "a" "a" "a" "a" "g" "t" "a" "t" "c" "a" "g" "g" "g" "t" "a" "g" "t"
     [793] "t" "c" "t" "g" "t" "t" "t" "c" "a" "a" "a" "c" "t" "t" "g" "c" "a" "t"
     [811] "g" "t" "g" "g" "a" "g" "c" "c" "a" "t" "g" "t" "g" "g" "c" "a" "c" "a"
     [829] "a" "a" "t" "a" "c" "t" "c" "a" "t" "g" "c" "c" "a" "g" "c" "t" "c" "a"
     [847] "t" "t" "a" "c" "a" "g" "c" "a" "t" "g" "a" "g" "a" "a" "c" "a" "g" "c"
     [865] "a" "g" "t" "t" "t" "a" "t" "t" "a" "c" "t" "c" "a" "c" "t" "a" "a" "a"
     [883] "g" "a" "c" "a" "g" "a" "a" "t" "g" "a" "a" "t" "g" "t" "a" "g" "a" "a"
     [901] "a" "a" "g" "g" "c" "t" "g" "a" "a" "t" "t" "c" "t" "g" "t" "a" "a" "t"
     [919] "a" "a" "a" "a" "g" "c" "a" "a" "a" "c" "a" "g" "c" "c" "t" "g" "g" "c"
     [937] "t" "t" "a" "g" "c" "a" "a" "g" "g" "a" "g" "c" "c" "a" "a" "c" "a" "t"
     [955] "a" "a" "c" "a" "g" "a" "t" "g" "g" "g" "c" "t" "g" "g" "a" "a" "g" "t"
     [973] "a" "a" "g" "g" "a" "a" "a" "c" "a" "t" "g" "t" "a" "a" "t" "g" "a" "t"
     [991] "a" "g" "g" "c" "g" "g" "a" "c" "t" "c" "c" "c" "a" "g" "c" "a" "c" "a"
    [1009] "g" "a" "a" "a" "a" "a" "a" "a" "g" "g" "t" "a" "g" "a" "t" "c" "t" "g"
    [1027] "a" "a" "t" "g" "c" "t" "g" "a" "t" "c" "c" "c" "c" "t" "g" "t" "g" "t"
    [1045] "g" "a" "g" "a" "g" "a" "a" "a" "a" "g" "a" "a" "t" "g" "g" "a" "a" "t"
    [1063] "a" "a" "g" "c" "a" "g" "a" "a" "a" "c" "t" "g" "c" "c" "a" "t" "g" "c"
    [1081] "t" "c" "a" "g" "a" "g" "a" "a" "t" "c" "c" "t" "a" "g" "a" "g" "a" "t"
    [1099] "a" "c" "t" "g" "a" "a" "g" "a" "t" "g" "t" "t" "c" "c" "t" "t" "g" "g"
    [1117] "a" "t" "a" "a" "c" "a" "c" "t" "a" "a" "a" "t" "a" "g" "c" "a" "g" "c"
    [1135] "a" "t" "t" "c" "a" "g" "a" "a" "a" "g" "t" "t" "a" "a" "t" "g" "a" "g"
    [1153] "t" "g" "g" "t" "t" "t" "t" "c" "c" "a" "g" "a" "a" "g" "t" "g" "a" "t"
    [1171] "g" "a" "a" "c" "t" "g" "t" "t" "a" "g" "g" "t" "t" "c" "t" "g" "a" "t"
    [1189] "g" "a" "c" "t" "c" "a" "c" "a" "t" "g" "a" "t" "g" "g" "g" "g" "a" "g"
    [1207] "t" "c" "t" "g" "a" "a" "t" "c" "a" "a" "a" "t" "g" "c" "c" "a" "a" "a"
    [1225] "g" "t" "a" "g" "c" "t" "g" "a" "t" "g" "t" "a" "t" "t" "g" "g" "a" "c"
    [1243] "g" "t" "t" "c" "t" "a" "a" "a" "t" "g" "a" "g" "g" "t" "a" "g" "a" "t"
    [1261] "g" "a" "a" "t" "a" "t" "t" "c" "t" "g" "g" "t" "t" "c" "t" "t" "c" "a"
    [1279] "g" "a" "g" "a" "a" "a" "a" "t" "a" "g" "a" "c" "t" "t" "a" "c" "t" "g"
    [1297] "g" "c" "c" "a" "g" "t" "g" "a" "t" "c" "c" "t" "c" "a" "t" "g" "a" "g"
    [1315] "g" "c" "t" "t" "t" "a" "a" "t" "a" "t" "g" "t" "a" "a" "a" "a" "g" "t"
    [1333] "g" "a" "a" "a" "g" "a" "g" "t" "t" "c" "a" "c" "t" "c" "c" "a" "a" "a"
    [1351] "t" "c" "a" "g" "t" "a" "g" "a" "g" "a" "g" "t" "a" "a" "t" "a" "t" "t"
    [1369] "g" "a" "a" "g" "a" "c" "a" "a" "a" "a" "t" "a" "t" "t" "t" "g" "g" "g"
    [1387] "a" "a" "a" "a" "c" "c" "t" "a" "t" "c" "g" "g" "a" "a" "g" "a" "a" "g"
    [1405] "g" "c" "a" "a" "g" "c" "c" "t" "c" "c" "c" "c" "a" "a" "c" "t" "t" "a"
    [1423] "a" "g" "c" "c" "a" "t" "g" "t" "a" "a" "c" "t" "g" "a" "a" "a" "a" "t"
    [1441] "c" "t" "a" "a" "t" "t" "a" "t" "a" "g" "g" "a" "g" "c" "a" "t" "t" "t"
    [1459] "g" "t" "t" "a" "c" "t" "g" "a" "g" "c" "c" "a" "c" "a" "g" "a" "t" "a"
    [1477] "a" "t" "a" "c" "a" "a" "g" "a" "g" "c" "g" "t" "c" "c" "c" "c" "t" "c"
    [1495] "a" "c" "a" "a" "a" "t" "a" "a" "a" "t" "t" "a" "a" "a" "g" "c" "g" "t"
    [1513] "a" "a" "a" "a" "g" "g" "a" "g" "a" "c" "c" "t" "a" "c" "a" "t" "c" "a"
    [1531] "g" "g" "c" "c" "t" "t" "c" "a" "t" "c" "c" "t" "g" "a" "g" "g" "a" "t"
    [1549] "t" "t" "t" "a" "t" "c" "a" "a" "g" "a" "a" "a" "g" "c" "a" "g" "a" "t"
    [1567] "t" "t" "g" "g" "c" "a" "g" "t" "t" "c" "a" "a" "a" "a" "g" "a" "c" "t"
    [1585] "c" "c" "t" "g" "a" "a" "a" "t" "g" "a" "t" "a" "a" "a" "t" "c" "a" "g"
    [1603] "g" "g" "a" "a" "c" "t" "a" "a" "c" "c" "a" "a" "a" "c" "g" "g" "a" "g"
    [1621] "c" "a" "g" "a" "a" "t" "g" "g" "t" "c" "a" "a" "g" "t" "g" "a" "t" "g"
    [1639] "a" "a" "t" "a" "t" "t" "a" "c" "t" "a" "a" "t" "a" "g" "t" "g" "g" "t"
    [1657] "c" "a" "t" "g" "a" "g" "a" "a" "t" "a" "a" "a" "a" "c" "a" "a" "a" "a"
    [1675] "g" "g" "t" "g" "a" "t" "t" "c" "t" "a" "t" "t" "c" "a" "g" "a" "a" "t"
    [1693] "g" "a" "g" "a" "a" "a" "a" "a" "t" "c" "c" "t" "a" "a" "c" "c" "c" "a"
    [1711] "a" "t" "a" "g" "a" "a" "t" "c" "a" "c" "t" "c" "g" "a" "a" "a" "a" "a"
    [1729] "g" "a" "a" "t" "c" "t" "g" "c" "t" "t" "t" "c" "a" "a" "a" "a" "c" "g"
    [1747] "a" "a" "a" "g" "c" "t" "g" "a" "a" "c" "c" "t" "a" "t" "a" "a" "g" "c"
    [1765] "a" "g" "c" "a" "g" "t" "a" "t" "a" "a" "g" "c" "a" "a" "t" "a" "t" "g"
    [1783] "g" "a" "a" "c" "t" "c" "g" "a" "a" "t" "t" "a" "a" "a" "t" "a" "t" "c"
    [1801] "c" "a" "c" "a" "a" "t" "t" "c" "a" "a" "a" "a" "g" "c" "a" "c" "c" "t"
    [1819] "a" "a" "a" "a" "a" "g" "a" "a" "t" "a" "g" "g" "c" "t" "g" "a" "g" "g"
    [1837] "a" "g" "g" "a" "a" "g" "t" "c" "t" "t" "c" "t" "a" "c" "c" "a" "g" "g"
    [1855] "c" "a" "t" "a" "t" "t" "c" "a" "t" "g" "c" "g" "c" "t" "t" "g" "a" "a"
    [1873] "c" "t" "a" "g" "t" "a" "g" "t" "c" "a" "g" "t" "a" "g" "a" "a" "a" "t"
    [1891] "c" "t" "a" "a" "g" "c" "c" "c" "a" "c" "c" "t" "a" "a" "t" "t" "g" "t"
    [1909] "a" "c" "t" "g" "a" "a" "t" "t" "g" "c" "a" "a" "a" "t" "t" "g" "a" "t"
    [1927] "a" "g" "t" "t" "g" "t" "t" "c" "t" "a" "g" "c" "a" "g" "t" "g" "a" "a"
    [1945] "g" "a" "g" "a" "t" "a" "a" "a" "g" "a" "a" "a" "a" "a" "a" "a" "a" "g"
    [1963] "t" "a" "c" "a" "a" "c" "c" "a" "a" "a" "t" "g" "c" "c" "a" "g" "t" "c"
    [1981] "a" "g" "g" "c" "a" "c" "a" "g" "c" "a" "g" "a" "a" "a" "c" "c" "t" "a"
    [1999] "c" "a" "a" "c" "t" "c" "a" "t" "g" "g" "a" "a" "g" "g" "t" "a" "a" "a"
    [2017] "g" "a" "a" "c" "c" "t" "g" "c" "a" "a" "c" "t" "g" "g" "a" "g" "c" "c"
    [2035] "a" "a" "g" "a" "a" "g" "a" "g" "t" "a" "a" "c" "a" "a" "g" "c" "c" "a"
    [2053] "a" "a" "t" "g" "a" "a" "c" "a" "g" "a" "c" "a" "a" "g" "t" "a" "a" "a"
    [2071] "a" "g" "a" "c" "a" "t" "g" "a" "c" "a" "g" "c" "g" "a" "t" "a" "c" "t"
    [2089] "t" "t" "c" "c" "c" "a" "g" "a" "g" "c" "t" "g" "a" "a" "g" "t" "t" "a"
    [2107] "a" "c" "a" "a" "a" "t" "g" "c" "a" "c" "c" "t" "g" "g" "t" "t" "c" "t"
    [2125] "t" "t" "t" "a" "c" "t" "a" "a" "g" "t" "g" "t" "t" "c" "a" "a" "a" "t"
    [2143] "a" "c" "c" "a" "g" "t" "g" "a" "a" "c" "t" "t" "a" "a" "a" "g" "a" "a"
    [2161] "t" "t" "t" "g" "t" "c" "a" "a" "t" "c" "c" "t" "a" "g" "c" "c" "t" "t"
    [2179] "c" "c" "a" "a" "g" "a" "g" "a" "a" "g" "a" "a" "a" "a" "a" "g" "a" "a"
    [2197] "g" "a" "g" "a" "a" "a" "c" "t" "a" "g" "a" "a" "a" "c" "a" "g" "t" "t"
    [2215] "a" "a" "a" "g" "t" "g" "t" "c" "t" "a" "a" "t" "a" "a" "t" "g" "c" "t"
    [2233] "g" "a" "a" "g" "a" "c" "c" "c" "c" "a" "a" "a" "g" "a" "t" "c" "t" "c"
    [2251] "a" "t" "g" "t" "t" "a" "a" "g" "t" "g" "g" "a" "g" "a" "a" "a" "g" "g"
    [2269] "g" "t" "t" "t" "t" "g" "c" "a" "a" "a" "c" "t" "g" "a" "a" "a" "g" "a"
    [2287] "t" "c" "t" "g" "t" "a" "g" "a" "g" "a" "g" "t" "a" "g" "c" "a" "g" "t"
    [2305] "a" "t" "t" "t" "c" "a" "t" "t" "g" "g" "t" "a" "c" "c" "t" "g" "g" "t"
    [2323] "a" "c" "t" "g" "a" "t" "t" "a" "t" "g" "g" "c" "a" "c" "t" "c" "a" "g"
    [2341] "g" "a" "a" "a" "g" "t" "a" "t" "c" "t" "c" "g" "t" "t" "a" "c" "t" "g"
    [2359] "g" "a" "a" "g" "t" "t" "a" "g" "c" "a" "c" "t" "c" "t" "a" "g" "g" "g"
    [2377] "a" "a" "g" "g" "c" "a" "a" "a" "a" "a" "c" "a" "g" "a" "a" "c" "c" "a"
    [2395] "a" "a" "t" "a" "a" "a" "t" "g" "t" "g" "t" "g" "a" "g" "t" "c" "a" "g"
    [2413] "t" "g" "t" "g" "c" "a" "g" "c" "a" "t" "t" "t" "g" "a" "a" "a" "a" "c"
    [2431] "c" "c" "c" "a" "a" "g" "g" "g" "a" "c" "t" "a" "a" "t" "t" "c" "a" "t"
    [2449] "g" "g" "t" "t" "g" "t" "t" "c" "c" "a" "a" "a" "g" "a" "t" "a" "a" "t"
    [2467] "a" "g" "a" "a" "a" "t" "g" "a" "c" "a" "c" "a" "g" "a" "a" "g" "g" "c"
    [2485] "t" "t" "t" "a" "a" "g" "t" "a" "t" "c" "c" "a" "t" "t" "g" "g" "g" "a"
    [2503] "c" "a" "t" "g" "a" "a" "g" "t" "t" "a" "a" "c" "c" "a" "c" "a" "g" "t"
    [2521] "c" "g" "g" "g" "a" "a" "a" "c" "a" "a" "g" "c" "a" "t" "a" "g" "a" "a"
    [2539] "a" "t" "g" "g" "a" "a" "g" "a" "a" "a" "g" "t" "g" "a" "a" "c" "t" "t"
    [2557] "g" "a" "t" "g" "c" "t" "c" "a" "g" "t" "a" "t" "t" "t" "g" "c" "a" "g"
    [2575] "a" "a" "t" "a" "c" "a" "t" "t" "c" "a" "a" "g" "g" "t" "t" "t" "c" "a"
    [2593] "a" "a" "g" "c" "g" "c" "c" "a" "g" "t" "c" "a" "t" "t" "t" "g" "c" "t"
    [2611] "c" "c" "g" "t" "t" "t" "t" "c" "a" "a" "a" "t" "c" "c" "a" "g" "g" "a"
    [2629] "a" "a" "t" "g" "c" "a" "g" "a" "a" "g" "a" "g" "g" "a" "a" "t" "g" "t"
    [2647] "g" "c" "a" "a" "c" "a" "t" "t" "c" "t" "c" "t" "g" "c" "c" "c" "a" "c"
    [2665] "t" "c" "t" "g" "g" "g" "t" "c" "c" "t" "t" "a" "a" "a" "g" "a" "a" "a"
    [2683] "c" "a" "a" "a" "g" "t" "c" "c" "a" "a" "a" "a" "g" "t" "c" "a" "c" "t"
    [2701] "t" "t" "t" "g" "a" "a" "t" "g" "t" "g" "a" "a" "c" "a" "a" "a" "a" "g"
    [2719] "g" "a" "a" "g" "a" "a" "a" "a" "t" "c" "a" "a" "g" "g" "a" "a" "a" "g"
    [2737] "a" "a" "t" "g" "a" "g" "t" "c" "t" "a" "a" "t" "a" "t" "c" "a" "a" "g"
    [2755] "c" "c" "t" "g" "t" "a" "c" "a" "g" "a" "c" "a" "g" "t" "t" "a" "a" "t"
    [2773] "a" "t" "c" "a" "c" "t" "g" "c" "a" "g" "g" "c" "t" "t" "t" "c" "c" "t"
    [2791] "g" "t" "g" "g" "t" "t" "g" "g" "t" "c" "a" "g" "a" "a" "a" "g" "a" "t"
    [2809] "a" "a" "g" "c" "c" "a" "g" "t" "t" "g" "a" "t" "a" "a" "t" "g" "c" "c"
    [2827] "a" "a" "a" "t" "g" "t" "a" "g" "t" "a" "t" "c" "a" "a" "a" "g" "g" "a"
    [2845] "g" "g" "c" "t" "c" "t" "a" "g" "g" "t" "t" "t" "t" "g" "t" "c" "t" "a"
    [2863] "t" "c" "a" "t" "c" "t" "c" "a" "g" "t" "t" "c" "a" "g" "a" "g" "g" "c"
    [2881] "a" "a" "c" "g" "a" "a" "a" "c" "t" "g" "g" "a" "c" "t" "c" "a" "t" "t"
    [2899] "a" "c" "t" "c" "c" "a" "a" "a" "t" "a" "a" "a" "c" "a" "t" "g" "g" "a"
    [2917] "c" "t" "t" "t" "t" "a" "c" "a" "a" "a" "a" "c" "c" "c" "a" "t" "a" "t"
    [2935] "c" "g" "t" "a" "t" "a" "c" "c" "a" "c" "c" "a" "c" "t" "t" "t" "t" "t"
    [2953] "c" "c" "c" "a" "t" "c" "a" "a" "g" "t" "c" "a" "t" "t" "t" "g" "t" "t"
    [2971] "a" "a" "a" "a" "c" "t" "a" "a" "a" "t" "g" "t" "a" "a" "g" "a" "a" "a"
    [2989] "a" "a" "t" "c" "t" "g" "c" "t" "a" "g" "a" "g" "g" "a" "a" "a" "a" "c"
    [3007] "t" "t" "t" "g" "a" "g" "g" "a" "a" "c" "a" "t" "t" "c" "a" "a" "t" "g"
    [3025] "t" "c" "a" "c" "c" "t" "g" "a" "a" "a" "g" "a" "g" "a" "a" "a" "t" "g"
    [3043] "g" "g" "a" "a" "a" "t" "g" "a" "g" "a" "a" "c" "a" "t" "t" "c" "c" "a"
    [3061] "a" "g" "t" "a" "c" "a" "g" "t" "g" "a" "g" "c" "a" "c" "a" "a" "t" "t"
    [3079] "a" "g" "c" "c" "g" "t" "a" "a" "t" "a" "a" "c" "a" "t" "t" "a" "g" "a"
    [3097] "g" "a" "a" "a" "a" "t" "g" "t" "t" "t" "t" "t" "a" "a" "a" "g" "a" "a"
    [3115] "g" "c" "c" "a" "g" "c" "t" "c" "a" "a" "g" "c" "a" "a" "t" "a" "t" "t"
    [3133] "a" "a" "t" "g" "a" "a" "g" "t" "a" "g" "g" "t" "t" "c" "c" "a" "g" "t"
    [3151] "a" "c" "t" "a" "a" "t" "g" "a" "a" "g" "t" "g" "g" "g" "c" "t" "c" "c"
    [3169] "a" "g" "t" "a" "t" "t" "a" "a" "t" "g" "a" "a" "a" "t" "a" "g" "g" "t"
    [3187] "t" "c" "c" "a" "g" "t" "g" "a" "t" "g" "a" "a" "a" "a" "c" "a" "t" "t"
    [3205] "c" "a" "a" "g" "c" "a" "g" "a" "a" "c" "t" "a" "g" "g" "t" "a" "g" "a"
    [3223] "a" "a" "c" "a" "g" "a" "g" "g" "g" "c" "c" "a" "a" "a" "a" "t" "t" "g"
    [3241] "a" "a" "t" "g" "c" "t" "a" "t" "g" "c" "t" "t" "a" "g" "a" "t" "t" "a"
    [3259] "g" "g" "g" "g" "t" "t" "t" "t" "g" "c" "a" "a" "c" "c" "t" "g" "a" "g"
    [3277] "g" "t" "c" "t" "a" "t" "a" "a" "a" "c" "a" "a" "a" "g" "t" "c" "t" "t"
    [3295] "c" "c" "t" "g" "g" "a" "a" "g" "t" "a" "a" "t" "t" "g" "t" "a" "a" "g"
    [3313] "c" "a" "t" "c" "c" "t" "g" "a" "a" "a" "t" "a" "a" "a" "a" "a" "a" "g"
    [3331] "c" "a" "a" "g" "a" "a" "t" "a" "t" "g" "a" "a" "g" "a" "a" "g" "t" "a"
    [3349] "g" "t" "t" "c" "a" "g" "a" "c" "t" "g" "t" "t" "a" "a" "t" "a" "c" "a"
    [3367] "g" "a" "t" "t" "t" "c" "t" "c" "t" "c" "c" "a" "t" "a" "t" "c" "t" "g"
    [3385] "a" "t" "t" "t" "c" "a" "g" "a" "t" "a" "a" "c" "t" "t" "a" "g" "a" "a"
    [3403] "c" "a" "g" "c" "c" "t" "a" "t" "g" "g" "g" "a" "a" "g" "t" "a" "g" "t"
    [3421] "c" "a" "t" "g" "c" "a" "t" "c" "t" "c" "a" "g" "g" "t" "t" "t" "g" "t"
    [3439] "t" "c" "t" "g" "a" "g" "a" "c" "a" "c" "c" "t" "g" "a" "t" "g" "a" "c"
    [3457] "c" "t" "g" "t" "t" "a" "g" "a" "t" "g" "a" "t" "g" "g" "t" "g" "a" "a"
    [3475] "a" "t" "a" "a" "a" "g" "g" "a" "a" "g" "a" "t" "a" "c" "t" "a" "g" "t"
    [3493] "t" "t" "t" "g" "c" "t" "g" "a" "a" "a" "a" "t" "g" "a" "c" "a" "t" "t"
    [3511] "a" "a" "g" "g" "a" "a" "a" "g" "t" "t" "c" "t" "g" "c" "t" "g" "t" "t"
    [3529] "t" "t" "t" "a" "g" "c" "a" "a" "a" "a" "g" "c" "g" "t" "c" "c" "a" "g"
    [3547] "a" "a" "a" "g" "g" "a" "g" "a" "g" "c" "t" "t" "a" "g" "c" "a" "g" "g"
    [3565] "a" "g" "t" "c" "c" "t" "a" "g" "c" "c" "c" "t" "t" "t" "c" "a" "c" "c"
    [3583] "c" "a" "t" "a" "c" "a" "c" "a" "t" "t" "t" "g" "g" "c" "t" "c" "a" "g"
    [3601] "g" "g" "t" "t" "a" "c" "c" "g" "a" "a" "g" "a" "g" "g" "g" "g" "c" "c"
    [3619] "a" "a" "g" "a" "a" "a" "t" "t" "a" "g" "a" "g" "t" "c" "c" "t" "c" "a"
    [3637] "g" "a" "a" "g" "a" "g" "a" "a" "c" "t" "t" "a" "t" "c" "t" "a" "g" "t"
    [3655] "g" "a" "g" "g" "a" "t" "g" "a" "a" "g" "a" "g" "c" "t" "t" "c" "c" "c"
    [3673] "t" "g" "c" "t" "t" "c" "c" "a" "a" "c" "a" "c" "t" "t" "g" "t" "t" "a"
    [3691] "t" "t" "t" "g" "g" "t" "a" "a" "a" "g" "t" "a" "a" "a" "c" "a" "a" "t"
    [3709] "a" "t" "a" "c" "c" "t" "t" "c" "t" "c" "a" "g" "t" "c" "t" "a" "c" "t"
    [3727] "a" "g" "g" "c" "a" "t" "a" "g" "c" "a" "c" "c" "g" "t" "t" "g" "c" "t"
    [3745] "a" "c" "c" "g" "a" "g" "t" "g" "t" "c" "t" "g" "t" "c" "t" "a" "a" "g"
    [3763] "a" "a" "c" "a" "c" "a" "g" "a" "g" "g" "a" "g" "a" "a" "t" "t" "t" "a"
    [3781] "t" "t" "a" "t" "c" "a" "t" "t" "g" "a" "a" "g" "a" "a" "t" "a" "g" "c"
    [3799] "t" "t" "a" "a" "a" "t" "g" "a" "c" "t" "g" "c" "a" "g" "t" "a" "a" "c"
    [3817] "c" "a" "g" "g" "t" "a" "a" "t" "a" "t" "t" "g" "g" "c" "a" "a" "a" "g"
    [3835] "g" "c" "a" "t" "c" "t" "c" "a" "g" "g" "a" "a" "c" "a" "t" "c" "a" "c"
    [3853] "c" "t" "t" "a" "g" "t" "g" "a" "g" "g" "a" "a" "a" "c" "a" "a" "a" "a"
    [3871] "t" "g" "t" "t" "c" "t" "g" "c" "t" "a" "g" "c" "t" "t" "g" "t" "t" "t"
    [3889] "t" "c" "t" "t" "c" "a" "c" "a" "g" "t" "g" "c" "a" "g" "t" "g" "a" "a"
    [3907] "t" "t" "g" "g" "a" "a" "g" "a" "c" "t" "t" "g" "a" "c" "t" "g" "c" "a"
    [3925] "a" "a" "t" "a" "c" "a" "a" "a" "c" "a" "c" "c" "c" "a" "g" "g" "a" "t"
    [3943] "c" "c" "t" "t" "t" "c" "t" "t" "g" "a" "t" "t" "g" "g" "t" "t" "c" "t"
    [3961] "t" "c" "c" "a" "a" "a" "c" "a" "a" "a" "t" "g" "a" "g" "g" "c" "a" "t"
    [3979] "c" "a" "g" "t" "c" "t" "g" "a" "a" "a" "g" "c" "c" "a" "g" "g" "g" "a"
    [3997] "g" "t" "t" "g" "g" "t" "c" "t" "g" "a" "g" "t" "g" "a" "c" "a" "a" "g"
    [4015] "g" "a" "a" "t" "t" "g" "g" "t" "t" "t" "c" "a" "g" "a" "t" "g" "a" "t"
    [4033] "g" "a" "a" "g" "a" "a" "a" "g" "a" "g" "g" "a" "a" "c" "g" "g" "g" "c"
    [4051] "t" "t" "g" "g" "a" "a" "g" "a" "a" "a" "a" "t" "a" "a" "t" "c" "a" "a"
    [4069] "g" "a" "a" "g" "a" "g" "c" "a" "a" "a" "g" "c" "a" "t" "g" "g" "a" "t"
    [4087] "t" "c" "a" "a" "a" "c" "t" "t" "a" "g" "g" "t" "g" "a" "a" "g" "c" "a"
    [4105] "g" "c" "a" "t" "c" "t" "g" "g" "g" "t" "g" "t" "g" "a" "g" "a" "g" "t"
    [4123] "g" "a" "a" "a" "c" "a" "a" "g" "c" "g" "t" "c" "t" "c" "t" "g" "a" "a"
    [4141] "g" "a" "c" "t" "g" "c" "t" "c" "a" "g" "g" "g" "c" "t" "a" "t" "c" "c"
    [4159] "t" "c" "t" "c" "a" "g" "a" "g" "t" "g" "a" "c" "a" "t" "t" "t" "t" "a"
    [4177] "a" "c" "c" "a" "c" "t" "c" "a" "g" "c" "a" "g" "a" "g" "g" "g" "a" "t"
    [4195] "a" "c" "c" "a" "t" "g" "c" "a" "a" "c" "a" "t" "a" "a" "c" "c" "t" "g"
    [4213] "a" "t" "a" "a" "a" "g" "c" "t" "c" "c" "a" "g" "c" "a" "g" "g" "a" "a"
    [4231] "a" "t" "g" "g" "c" "t" "g" "a" "a" "c" "t" "a" "g" "a" "a" "g" "c" "t"
    [4249] "g" "t" "g" "t" "t" "a" "g" "a" "a" "c" "a" "g" "c" "a" "t" "g" "g" "g"
    [4267] "a" "g" "c" "c" "a" "g" "c" "c" "t" "t" "c" "t" "a" "a" "c" "a" "g" "c"
    [4285] "t" "a" "c" "c" "c" "t" "t" "c" "c" "a" "t" "c" "a" "t" "a" "a" "g" "t"
    [4303] "g" "a" "c" "t" "c" "t" "t" "c" "t" "g" "c" "c" "c" "t" "t" "g" "a" "g"
    [4321] "g" "a" "c" "c" "t" "g" "c" "g" "a" "a" "a" "t" "c" "c" "a" "g" "a" "a"
    [4339] "c" "a" "a" "a" "g" "c" "a" "c" "a" "t" "c" "a" "g" "a" "a" "a" "a" "a"
    [4357] "g" "c" "a" "g" "t" "a" "t" "t" "a" "a" "c" "t" "t" "c" "a" "c" "a" "g"
    [4375] "a" "a" "a" "a" "g" "t" "a" "g" "t" "g" "a" "a" "t" "a" "c" "c" "c" "t"
    [4393] "a" "t" "a" "a" "g" "c" "c" "a" "g" "a" "a" "t" "c" "c" "a" "g" "a" "a"
    [4411] "g" "g" "c" "c" "t" "t" "t" "c" "t" "g" "c" "t" "g" "a" "c" "a" "a" "g"
    [4429] "t" "t" "t" "g" "a" "g" "g" "t" "g" "t" "c" "t" "g" "c" "a" "g" "a" "t"
    [4447] "a" "g" "t" "t" "c" "t" "a" "c" "c" "a" "g" "t" "a" "a" "a" "a" "a" "t"
    [4465] "a" "a" "a" "g" "a" "a" "c" "c" "a" "g" "g" "a" "g" "t" "g" "g" "a" "a"
    [4483] "a" "g" "g" "t" "c" "a" "t" "c" "c" "c" "c" "t" "t" "c" "t" "a" "a" "a"
    [4501] "t" "g" "c" "c" "c" "a" "t" "c" "a" "t" "t" "a" "g" "a" "t" "g" "a" "t"
    [4519] "a" "g" "g" "t" "g" "g" "t" "a" "c" "a" "t" "g" "c" "a" "c" "a" "g" "t"
    [4537] "t" "g" "c" "t" "c" "t" "g" "g" "g" "a" "g" "t" "c" "t" "t" "c" "a" "g"
    [4555] "a" "a" "t" "a" "g" "a" "a" "a" "c" "t" "a" "c" "c" "c" "a" "t" "c" "t"
    [4573] "c" "a" "a" "g" "a" "g" "g" "a" "g" "c" "t" "c" "a" "t" "t" "a" "a" "g"
    [4591] "g" "t" "t" "g" "t" "t" "g" "a" "t" "g" "t" "g" "g" "a" "g" "g" "a" "g"
    [4609] "c" "a" "a" "c" "a" "g" "c" "t" "g" "g" "a" "a" "g" "a" "g" "t" "c" "t"
    [4627] "g" "g" "g" "c" "c" "a" "c" "a" "c" "g" "a" "t" "t" "t" "g" "a" "c" "g"
    [4645] "g" "a" "a" "a" "c" "a" "t" "c" "t" "t" "a" "c" "t" "t" "g" "c" "c" "a"
    [4663] "a" "g" "g" "c" "a" "a" "g" "a" "t" "c" "t" "a" "g" "a" "g" "g" "g" "a"
    [4681] "a" "c" "c" "c" "c" "t" "t" "a" "c" "c" "t" "g" "g" "a" "a" "t" "c" "t"
    [4699] "g" "g" "a" "a" "t" "c" "a" "g" "c" "c" "t" "c" "t" "t" "c" "t" "c" "t"
    [4717] "g" "a" "t" "g" "a" "c" "c" "c" "t" "g" "a" "a" "t" "c" "t" "g" "a" "t"
    [4735] "c" "c" "t" "t" "c" "t" "g" "a" "a" "g" "a" "c" "a" "g" "a" "g" "c" "c"
    [4753] "c" "c" "a" "g" "a" "g" "t" "c" "a" "g" "c" "t" "c" "g" "t" "g" "t" "t"
    [4771] "g" "g" "c" "a" "a" "c" "a" "t" "a" "c" "c" "a" "t" "c" "t" "t" "c" "a"
    [4789] "a" "c" "c" "t" "c" "t" "g" "c" "a" "t" "t" "g" "a" "a" "a" "g" "t" "t"
    [4807] "c" "c" "c" "c" "a" "a" "t" "t" "g" "a" "a" "a" "g" "t" "t" "g" "c" "a"
    [4825] "g" "a" "a" "t" "c" "t" "g" "c" "c" "c" "a" "g" "a" "g" "t" "c" "c" "a"
    [4843] "g" "c" "t" "g" "c" "t" "g" "c" "t" "c" "a" "t" "a" "c" "t" "a" "c" "t"
    [4861] "g" "a" "t" "a" "c" "t" "g" "c" "t" "g" "g" "g" "t" "a" "t" "a" "a" "t"
    [4879] "g" "c" "a" "a" "t" "g" "g" "a" "a" "g" "a" "a" "a" "g" "t" "g" "t" "g"
    [4897] "a" "g" "c" "a" "g" "g" "g" "a" "g" "a" "a" "g" "c" "c" "a" "g" "a" "a"
    [4915] "t" "t" "g" "a" "c" "a" "g" "c" "t" "t" "c" "a" "a" "c" "a" "g" "a" "a"
    [4933] "a" "g" "g" "g" "t" "c" "a" "a" "c" "a" "a" "a" "a" "g" "a" "a" "t" "g"
    [4951] "t" "c" "c" "a" "t" "g" "g" "t" "g" "g" "t" "g" "t" "c" "t" "g" "g" "c"
    [4969] "c" "t" "g" "a" "c" "c" "c" "c" "a" "g" "a" "a" "g" "a" "a" "t" "t" "t"
    [4987] "a" "t" "g" "c" "t" "c" "g" "t" "g" "t" "a" "c" "a" "a" "g" "t" "t" "t"
    [5005] "g" "c" "c" "a" "g" "a" "a" "a" "a" "c" "a" "c" "c" "a" "c" "a" "t" "c"
    [5023] "a" "c" "t" "t" "t" "a" "a" "c" "t" "a" "a" "t" "c" "t" "a" "a" "t" "t"
    [5041] "a" "c" "t" "g" "a" "a" "g" "a" "g" "a" "c" "t" "a" "c" "t" "c" "a" "t"
    [5059] "g" "t" "t" "g" "t" "t" "a" "t" "g" "a" "a" "a" "a" "c" "a" "g" "a" "t"
    [5077] "g" "c" "t" "g" "a" "g" "t" "t" "t" "g" "t" "g" "t" "g" "t" "g" "a" "a"
    [5095] "c" "g" "g" "a" "c" "a" "c" "t" "g" "a" "a" "a" "t" "a" "t" "t" "t" "t"
    [5113] "c" "t" "a" "g" "g" "a" "a" "t" "t" "g" "c" "g" "g" "g" "a" "g" "g" "a"
    [5131] "a" "a" "a" "t" "g" "g" "g" "t" "a" "g" "t" "t" "a" "g" "c" "t" "a" "t"
    [5149] "t" "t" "c" "t" "g" "g" "g" "t" "g" "a" "c" "c" "c" "a" "g" "t" "c" "t"
    [5167] "a" "t" "t" "a" "a" "a" "g" "a" "a" "a" "g" "a" "a" "a" "a" "a" "t" "g"
    [5185] "c" "t" "g" "a" "a" "t" "g" "a" "g" "c" "a" "t" "g" "a" "t" "t" "t" "t"
    [5203] "g" "a" "a" "g" "t" "c" "a" "g" "a" "g" "g" "a" "g" "a" "t" "g" "t" "g"
    [5221] "g" "t" "c" "a" "a" "t" "g" "g" "a" "a" "g" "a" "a" "a" "c" "c" "a" "c"
    [5239] "c" "a" "a" "g" "g" "t" "c" "c" "a" "a" "a" "g" "c" "g" "a" "g" "c" "a"
    [5257] "a" "g" "a" "g" "a" "a" "t" "c" "c" "c" "a" "g" "g" "a" "c" "a" "g" "a"
    [5275] "a" "a" "g" "a" "t" "c" "t" "t" "c" "a" "g" "g" "g" "g" "g" "c" "t" "a"
    [5293] "g" "a" "a" "a" "t" "c" "t" "g" "t" "t" "g" "c" "t" "a" "t" "g" "g" "g"
    [5311] "c" "c" "c" "t" "t" "c" "a" "c" "c" "a" "a" "c" "a" "t" "g" "c" "c" "c"
    [5329] "a" "c" "a" "g" "a" "t" "c" "a" "a" "c" "t" "g" "g" "a" "a" "t" "g" "g"
    [5347] "a" "t" "g" "g" "t" "a" "c" "a" "g" "c" "t" "g" "t" "g" "t" "g" "g" "t"
    [5365] "g" "c" "t" "t" "c" "t" "g" "t" "g" "g" "t" "g" "a" "a" "g" "g" "a" "g"
    [5383] "c" "t" "t" "t" "c" "a" "t" "c" "a" "t" "t" "c" "a" "c" "c" "c" "t" "t"
    [5401] "g" "g" "c" "a" "c" "a" "g" "g" "t" "g" "t" "c" "c" "a" "c" "c" "c" "a"
    [5419] "a" "t" "t" "g" "t" "g" "g" "t" "t" "g" "t" "g" "c" "a" "g" "c" "c" "a"
    [5437] "g" "a" "t" "g" "c" "c" "t" "g" "g" "a" "c" "a" "g" "a" "g" "g" "a" "c"
    [5455] "a" "a" "t" "g" "g" "c" "t" "t" "c" "c" "a" "t" "g" "c" "a" "a" "t" "t"
    [5473] "g" "g" "g" "c" "a" "g" "a" "t" "g" "t" "g" "t" "g" "a" "g" "g" "c" "a"
    [5491] "c" "c" "t" "g" "t" "g" "g" "t" "g" "a" "c" "c" "c" "g" "a" "g" "a" "g"
    [5509] "t" "g" "g" "g" "t" "g" "t" "t" "g" "g" "a" "c" "a" "g" "t" "g" "t" "a"
    [5527] "g" "c" "a" "c" "t" "c" "t" "a" "c" "c" "a" "g" "t" "g" "c" "c" "a" "g"
    [5545] "g" "a" "g" "c" "t" "g" "g" "a" "c" "a" "c" "c" "t" "a" "c" "c" "t" "g"
    [5563] "a" "t" "a" "c" "c" "c" "c" "a" "g" "a" "t" "c" "c" "c" "c" "c" "a" "c"
    [5581] "a" "g" "c" "c" "a" "c" "t" "a" "c" "t" "g" "a"
    attr(,"name")
    [1] "DQ190451.BRCA1"
    attr(,"Annot")
    [1] ">DQ190451.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $DQ190452.BRCA1
       [1] "a" "t" "g" "g" "a" "t" "t" "t" "a" "t" "c" "t" "g" "c" "t" "c" "t" "t"
      [19] "c" "g" "c" "g" "t" "t" "g" "a" "a" "g" "a" "a" "g" "t" "a" "c" "a" "a"
      [37] "a" "a" "t" "g" "t" "c" "a" "t" "t" "a" "a" "t" "g" "c" "t" "a" "t" "g"
      [55] "c" "a" "g" "a" "a" "a" "a" "t" "c" "t" "t" "a" "g" "a" "g" "t" "g" "t"
      [73] "c" "c" "c" "a" "t" "c" "t" "g" "t" "c" "t" "g" "g" "a" "g" "t" "t" "g"
      [91] "a" "t" "c" "a" "a" "g" "g" "a" "a" "c" "c" "t" "g" "t" "c" "t" "c" "c"
     [109] "a" "c" "a" "a" "a" "g" "t" "g" "t" "g" "a" "c" "c" "a" "c" "a" "t" "a"
     [127] "t" "t" "t" "t" "g" "c" "a" "a" "a" "t" "t" "t" "t" "g" "c" "a" "t" "g"
     [145] "c" "t" "g" "a" "a" "a" "c" "t" "t" "c" "t" "c" "a" "a" "c" "c" "a" "g"
     [163] "a" "a" "g" "a" "a" "a" "g" "g" "g" "c" "c" "t" "t" "c" "a" "c" "a" "g"
     [181] "t" "g" "t" "c" "c" "t" "t" "t" "a" "t" "g" "t" "a" "a" "g" "a" "a" "t"
     [199] "g" "a" "t" "a" "t" "a" "a" "c" "c" "a" "a" "a" "a" "g" "g" "a" "g" "c"
     [217] "c" "t" "a" "c" "a" "a" "g" "a" "a" "a" "g" "t" "a" "c" "g" "a" "g" "a"
     [235] "t" "t" "t" "a" "g" "t" "c" "a" "a" "c" "t" "t" "g" "t" "t" "g" "a" "a"
     [253] "g" "a" "g" "c" "t" "a" "t" "t" "g" "a" "a" "a" "a" "t" "c" "a" "t" "t"
     [271] "t" "g" "t" "g" "c" "t" "t" "t" "t" "c" "a" "g" "c" "t" "t" "g" "a" "c"
     [289] "a" "c" "a" "g" "g" "t" "t" "t" "g" "g" "a" "g" "t" "a" "t" "g" "c" "a"
     [307] "a" "a" "c" "a" "g" "c" "t" "a" "t" "a" "a" "t" "t" "t" "t" "g" "c" "a"
     [325] "a" "a" "a" "a" "a" "g" "g" "a" "a" "a" "a" "t" "a" "a" "c" "t" "c" "t"
     [343] "c" "c" "t" "g" "a" "a" "c" "a" "t" "c" "t" "a" "a" "a" "a" "g" "a" "t"
     [361] "g" "a" "a" "g" "t" "t" "t" "c" "t" "a" "t" "c" "a" "t" "c" "c" "a" "a"
     [379] "a" "g" "t" "a" "t" "g" "g" "g" "c" "t" "a" "c" "a" "g" "a" "a" "a" "c"
     [397] "c" "g" "t" "g" "c" "c" "a" "a" "a" "a" "g" "a" "c" "t" "t" "c" "t" "a"
     [415] "c" "a" "g" "a" "g" "t" "g" "a" "a" "c" "c" "c" "g" "a" "a" "a" "a" "t"
     [433] "c" "c" "t" "t" "c" "c" "t" "t" "g" "c" "a" "g" "g" "a" "a" "a" "c" "c"
     [451] "a" "g" "t" "c" "t" "c" "a" "g" "t" "g" "t" "c" "c" "a" "a" "c" "t" "c"
     [469] "t" "c" "t" "a" "a" "c" "c" "t" "t" "g" "g" "a" "a" "c" "t" "g" "t" "g"
     [487] "a" "g" "a" "a" "c" "t" "c" "t" "g" "a" "g" "g" "a" "c" "a" "a" "a" "g"
     [505] "c" "a" "g" "c" "g" "g" "a" "t" "a" "c" "a" "a" "c" "c" "t" "c" "a" "a"
     [523] "a" "a" "g" "a" "c" "g" "t" "c" "t" "g" "t" "c" "t" "a" "c" "a" "t" "t"
     [541] "g" "a" "a" "t" "t" "g" "g" "g" "a" "t" "c" "t" "g" "a" "t" "t" "c" "t"
     [559] "t" "c" "t" "g" "a" "a" "g" "a" "t" "a" "c" "c" "g" "t" "t" "a" "a" "t"
     [577] "a" "a" "g" "g" "c" "a" "a" "c" "t" "t" "a" "t" "t" "g" "c" "a" "g" "t"
     [595] "g" "t" "g" "g" "g" "a" "g" "a" "t" "c" "a" "a" "g" "a" "a" "t" "t" "g"
     [613] "t" "t" "a" "c" "a" "a" "a" "t" "c" "a" "c" "c" "c" "c" "t" "c" "a" "a"
     [631] "g" "g" "a" "a" "c" "c" "a" "g" "g" "g" "a" "t" "g" "a" "a" "a" "t" "c"
     [649] "a" "g" "t" "t" "t" "g" "g" "a" "t" "t" "c" "t" "g" "c" "a" "a" "a" "a"
     [667] "a" "a" "g" "g" "c" "t" "g" "c" "t" "t" "g" "t" "g" "a" "a" "t" "t" "t"
     [685] "t" "c" "t" "g" "a" "g" "a" "c" "g" "g" "a" "t" "g" "t" "a" "a" "c" "a"
     [703] "a" "a" "t" "a" "c" "t" "g" "a" "a" "c" "a" "t" "c" "a" "t" "c" "a" "a"
     [721] "c" "c" "c" "a" "g" "t" "a" "a" "t" "a" "a" "t" "g" "a" "t" "t" "t" "g"
     [739] "a" "a" "c" "a" "c" "c" "a" "c" "t" "g" "a" "g" "a" "a" "g" "c" "g" "t"
     [757] "g" "c" "a" "g" "c" "t" "g" "a" "g" "a" "g" "g" "c" "a" "t" "c" "c" "a"
     [775] "g" "a" "a" "a" "a" "g" "t" "a" "t" "c" "a" "g" "g" "g" "t" "a" "g" "t"
     [793] "t" "c" "t" "g" "t" "t" "t" "c" "a" "a" "a" "c" "t" "t" "g" "c" "a" "t"
     [811] "g" "t" "g" "g" "a" "g" "c" "c" "a" "t" "g" "t" "g" "g" "c" "a" "c" "a"
     [829] "a" "a" "t" "a" "c" "t" "c" "a" "t" "g" "c" "c" "a" "g" "c" "t" "c" "a"
     [847] "t" "t" "a" "c" "a" "g" "c" "a" "t" "g" "a" "g" "a" "a" "c" "a" "g" "c"
     [865] "a" "g" "t" "t" "t" "a" "t" "t" "a" "c" "t" "c" "a" "c" "t" "a" "a" "a"
     [883] "g" "a" "c" "a" "g" "a" "a" "t" "g" "a" "a" "t" "g" "t" "a" "g" "a" "a"
     [901] "a" "a" "g" "g" "c" "t" "g" "a" "a" "t" "t" "c" "t" "g" "t" "a" "a" "t"
     [919] "a" "a" "a" "a" "g" "c" "a" "a" "a" "c" "a" "g" "c" "c" "t" "g" "g" "c"
     [937] "t" "t" "a" "g" "c" "a" "a" "g" "g" "a" "g" "c" "c" "a" "a" "c" "a" "t"
     [955] "a" "a" "c" "a" "g" "a" "t" "g" "g" "g" "c" "t" "g" "g" "a" "a" "g" "t"
     [973] "a" "a" "g" "g" "a" "a" "a" "c" "a" "t" "g" "t" "a" "a" "t" "g" "a" "t"
     [991] "a" "g" "g" "c" "g" "g" "a" "c" "t" "c" "c" "c" "a" "g" "c" "a" "c" "a"
    [1009] "g" "a" "a" "a" "a" "a" "a" "a" "g" "g" "t" "a" "g" "a" "t" "c" "t" "g"
    [1027] "a" "a" "t" "g" "c" "t" "g" "a" "t" "c" "c" "c" "c" "t" "g" "t" "g" "t"
    [1045] "g" "a" "g" "a" "g" "a" "a" "a" "a" "g" "a" "a" "t" "g" "g" "a" "a" "t"
    [1063] "a" "a" "g" "c" "a" "g" "a" "a" "a" "c" "t" "g" "c" "c" "a" "t" "g" "c"
    [1081] "t" "c" "a" "g" "a" "g" "a" "a" "t" "c" "c" "t" "a" "g" "a" "g" "a" "t"
    [1099] "a" "c" "t" "g" "a" "a" "g" "a" "t" "g" "t" "t" "c" "c" "t" "t" "g" "g"
    [1117] "a" "t" "a" "a" "c" "a" "c" "t" "a" "a" "a" "t" "a" "g" "c" "a" "g" "c"
    [1135] "a" "t" "t" "c" "a" "g" "a" "a" "a" "g" "t" "t" "a" "a" "t" "g" "a" "g"
    [1153] "t" "g" "g" "t" "t" "t" "t" "c" "c" "a" "g" "a" "a" "g" "t" "g" "a" "t"
    [1171] "g" "a" "a" "c" "t" "g" "t" "t" "a" "g" "g" "t" "t" "c" "t" "g" "a" "t"
    [1189] "g" "a" "c" "t" "c" "a" "c" "a" "t" "g" "a" "t" "g" "g" "g" "g" "a" "g"
    [1207] "t" "c" "t" "g" "a" "a" "t" "c" "a" "a" "a" "t" "g" "c" "c" "a" "a" "a"
    [1225] "g" "t" "a" "g" "c" "t" "g" "a" "t" "g" "t" "a" "t" "t" "g" "g" "a" "c"
    [1243] "g" "t" "t" "c" "t" "a" "a" "a" "t" "g" "a" "g" "g" "t" "a" "g" "a" "t"
    [1261] "g" "a" "a" "t" "a" "t" "t" "c" "t" "g" "g" "t" "t" "c" "t" "t" "c" "a"
    [1279] "g" "a" "g" "a" "a" "a" "a" "t" "a" "g" "a" "c" "t" "t" "a" "c" "t" "g"
    [1297] "g" "c" "c" "a" "g" "t" "g" "a" "t" "c" "c" "t" "c" "a" "t" "g" "a" "g"
    [1315] "g" "c" "t" "t" "t" "a" "a" "t" "a" "t" "g" "t" "a" "a" "a" "a" "g" "t"
    [1333] "g" "a" "a" "a" "g" "a" "g" "t" "t" "c" "a" "c" "t" "c" "c" "a" "a" "a"
    [1351] "t" "c" "a" "g" "t" "a" "g" "a" "g" "a" "g" "t" "a" "a" "t" "a" "t" "t"
    [1369] "g" "a" "a" "g" "a" "c" "a" "a" "a" "a" "t" "a" "t" "t" "t" "g" "g" "g"
    [1387] "a" "a" "a" "a" "c" "c" "t" "a" "t" "c" "g" "g" "a" "a" "g" "a" "a" "g"
    [1405] "g" "c" "a" "a" "g" "c" "c" "t" "c" "c" "c" "c" "a" "a" "c" "t" "t" "a"
    [1423] "a" "g" "c" "c" "a" "t" "g" "t" "a" "a" "c" "t" "g" "a" "a" "a" "a" "t"
    [1441] "c" "t" "a" "a" "t" "t" "a" "t" "a" "g" "g" "a" "g" "c" "a" "t" "t" "t"
    [1459] "g" "t" "t" "a" "c" "t" "g" "a" "g" "c" "c" "a" "c" "a" "g" "a" "t" "a"
    [1477] "a" "t" "a" "c" "a" "a" "g" "a" "g" "c" "g" "t" "c" "c" "c" "c" "t" "c"
    [1495] "a" "c" "a" "a" "a" "t" "a" "a" "a" "t" "t" "a" "a" "a" "g" "c" "g" "t"
    [1513] "a" "a" "a" "a" "g" "g" "a" "g" "a" "c" "c" "t" "a" "c" "a" "t" "c" "a"
    [1531] "g" "g" "c" "c" "t" "t" "c" "a" "t" "c" "c" "t" "g" "a" "g" "g" "a" "t"
    [1549] "t" "t" "t" "a" "t" "c" "a" "a" "g" "a" "a" "a" "g" "c" "a" "g" "a" "t"
    [1567] "t" "t" "g" "g" "c" "a" "g" "t" "t" "c" "a" "a" "a" "a" "g" "a" "c" "t"
    [1585] "c" "c" "t" "g" "a" "a" "a" "t" "g" "a" "t" "a" "a" "a" "t" "c" "a" "g"
    [1603] "g" "g" "a" "a" "c" "t" "a" "a" "c" "c" "a" "a" "a" "c" "g" "g" "a" "g"
    [1621] "c" "a" "g" "a" "a" "t" "g" "g" "t" "c" "a" "a" "g" "t" "g" "a" "t" "g"
    [1639] "a" "a" "t" "a" "t" "t" "a" "c" "t" "a" "a" "t" "a" "g" "t" "g" "g" "t"
    [1657] "c" "a" "t" "g" "a" "g" "a" "a" "t" "a" "a" "a" "a" "c" "a" "a" "a" "a"
    [1675] "g" "g" "t" "g" "a" "t" "t" "c" "t" "a" "t" "t" "c" "a" "g" "a" "a" "t"
    [1693] "g" "a" "g" "a" "a" "a" "a" "a" "t" "c" "c" "t" "a" "a" "c" "c" "c" "a"
    [1711] "a" "t" "a" "g" "a" "a" "t" "c" "a" "c" "t" "c" "g" "a" "a" "a" "a" "a"
    [1729] "g" "a" "a" "t" "c" "t" "g" "c" "t" "t" "t" "c" "a" "a" "a" "a" "c" "g"
    [1747] "a" "a" "a" "g" "c" "t" "g" "a" "a" "c" "c" "t" "a" "t" "a" "a" "g" "c"
    [1765] "a" "g" "c" "a" "g" "t" "a" "t" "a" "a" "g" "c" "a" "a" "t" "a" "t" "g"
    [1783] "g" "a" "a" "c" "t" "c" "g" "a" "a" "t" "t" "a" "a" "a" "t" "a" "t" "c"
    [1801] "c" "a" "c" "a" "a" "t" "t" "c" "a" "a" "a" "a" "g" "c" "a" "c" "c" "t"
    [1819] "a" "a" "a" "a" "a" "g" "a" "a" "t" "a" "g" "g" "c" "t" "g" "a" "g" "g"
    [1837] "a" "g" "g" "a" "a" "g" "t" "c" "t" "t" "c" "t" "a" "c" "c" "a" "g" "g"
    [1855] "c" "a" "t" "a" "t" "t" "c" "a" "t" "g" "c" "g" "c" "t" "t" "g" "a" "a"
    [1873] "c" "t" "a" "g" "t" "a" "g" "t" "c" "a" "g" "t" "a" "g" "a" "a" "a" "t"
    [1891] "c" "t" "a" "a" "g" "c" "c" "c" "a" "c" "c" "t" "a" "a" "t" "t" "g" "t"
    [1909] "a" "c" "t" "g" "a" "a" "t" "t" "g" "c" "a" "a" "a" "t" "t" "g" "a" "t"
    [1927] "a" "g" "t" "t" "g" "t" "t" "c" "t" "a" "g" "c" "a" "g" "t" "g" "a" "a"
    [1945] "g" "a" "g" "a" "t" "a" "a" "a" "g" "a" "a" "a" "a" "a" "a" "a" "a" "g"
    [1963] "t" "a" "c" "a" "a" "c" "c" "a" "a" "a" "t" "g" "c" "c" "a" "g" "t" "c"
    [1981] "a" "g" "g" "c" "a" "c" "a" "g" "c" "a" "g" "a" "a" "a" "c" "c" "t" "a"
    [1999] "c" "a" "a" "c" "t" "c" "a" "t" "g" "g" "a" "a" "g" "g" "t" "a" "a" "a"
    [2017] "g" "a" "a" "c" "c" "t" "g" "c" "a" "a" "c" "t" "g" "g" "a" "g" "c" "c"
    [2035] "a" "a" "g" "a" "a" "g" "a" "g" "t" "a" "a" "c" "a" "a" "g" "c" "c" "a"
    [2053] "a" "a" "t" "g" "a" "a" "c" "a" "g" "a" "c" "a" "a" "g" "t" "a" "a" "a"
    [2071] "a" "g" "a" "c" "a" "t" "g" "a" "c" "a" "g" "c" "g" "a" "t" "a" "c" "t"
    [2089] "t" "t" "c" "c" "c" "a" "g" "a" "g" "c" "t" "g" "a" "a" "g" "t" "t" "a"
    [2107] "a" "c" "a" "a" "a" "t" "g" "c" "a" "c" "c" "t" "g" "g" "t" "t" "c" "t"
    [2125] "t" "t" "t" "a" "c" "t" "a" "a" "g" "t" "g" "t" "t" "c" "a" "a" "a" "t"
    [2143] "a" "c" "c" "a" "g" "t" "g" "a" "a" "c" "t" "t" "a" "a" "a" "g" "a" "a"
    [2161] "t" "t" "t" "g" "t" "c" "a" "a" "t" "c" "c" "t" "a" "g" "c" "c" "t" "t"
    [2179] "c" "c" "a" "a" "g" "a" "g" "a" "a" "g" "a" "a" "a" "a" "a" "g" "a" "a"
    [2197] "g" "a" "g" "a" "a" "a" "c" "t" "a" "g" "a" "a" "a" "c" "a" "g" "t" "t"
    [2215] "a" "a" "a" "g" "t" "g" "t" "c" "t" "a" "a" "t" "a" "a" "t" "g" "c" "t"
    [2233] "g" "a" "a" "g" "a" "c" "c" "c" "c" "a" "a" "a" "g" "a" "t" "c" "t" "c"
    [2251] "a" "t" "g" "t" "t" "a" "a" "g" "t" "g" "g" "a" "g" "a" "a" "a" "g" "g"
    [2269] "g" "t" "t" "t" "t" "g" "c" "a" "a" "a" "c" "t" "g" "a" "a" "a" "g" "a"
    [2287] "t" "c" "t" "g" "t" "a" "g" "a" "g" "a" "g" "t" "a" "g" "c" "a" "g" "t"
    [2305] "a" "t" "t" "t" "c" "a" "t" "t" "g" "g" "t" "a" "c" "c" "t" "g" "g" "t"
    [2323] "a" "c" "t" "g" "a" "t" "t" "a" "t" "g" "g" "c" "a" "c" "t" "c" "a" "g"
    [2341] "g" "a" "a" "a" "g" "t" "a" "t" "c" "t" "c" "g" "t" "t" "a" "c" "t" "g"
    [2359] "g" "a" "a" "g" "t" "t" "a" "g" "c" "a" "c" "t" "c" "t" "a" "g" "g" "g"
    [2377] "a" "a" "g" "g" "c" "a" "a" "a" "a" "a" "c" "a" "g" "a" "a" "c" "c" "a"
    [2395] "a" "a" "t" "a" "a" "a" "t" "g" "t" "g" "t" "g" "a" "g" "t" "c" "a" "g"
    [2413] "t" "g" "t" "g" "c" "a" "g" "c" "a" "t" "t" "t" "g" "a" "a" "a" "a" "c"
    [2431] "c" "c" "c" "a" "a" "g" "g" "g" "a" "c" "t" "a" "a" "t" "t" "c" "a" "t"
    [2449] "g" "g" "t" "t" "g" "t" "t" "c" "c" "a" "a" "a" "g" "a" "t" "a" "a" "t"
    [2467] "a" "g" "a" "a" "a" "t" "g" "a" "c" "a" "c" "a" "g" "a" "a" "g" "g" "c"
    [2485] "t" "t" "t" "a" "a" "g" "t" "a" "t" "c" "c" "a" "t" "t" "g" "g" "g" "a"
    [2503] "c" "a" "t" "g" "a" "a" "g" "t" "t" "a" "a" "c" "c" "a" "c" "a" "g" "t"
    [2521] "c" "g" "g" "g" "a" "a" "a" "c" "a" "a" "g" "c" "a" "t" "a" "g" "a" "a"
    [2539] "a" "t" "g" "g" "a" "a" "g" "a" "a" "a" "g" "t" "g" "a" "a" "c" "t" "t"
    [2557] "g" "a" "t" "g" "c" "t" "c" "a" "g" "t" "a" "t" "t" "t" "g" "c" "a" "g"
    [2575] "a" "a" "t" "a" "c" "a" "t" "t" "c" "a" "a" "g" "g" "t" "t" "t" "c" "a"
    [2593] "a" "a" "g" "c" "g" "c" "c" "a" "g" "t" "c" "a" "t" "t" "t" "g" "c" "t"
    [2611] "c" "c" "g" "t" "t" "t" "t" "c" "a" "a" "a" "t" "c" "c" "a" "g" "g" "a"
    [2629] "a" "a" "t" "g" "c" "a" "g" "a" "a" "g" "a" "g" "g" "a" "a" "t" "g" "t"
    [2647] "g" "c" "a" "a" "c" "a" "t" "t" "c" "t" "c" "t" "g" "c" "c" "c" "a" "c"
    [2665] "t" "c" "t" "g" "g" "g" "t" "c" "c" "t" "t" "a" "a" "a" "g" "a" "a" "a"
    [2683] "c" "a" "a" "a" "g" "t" "c" "c" "a" "a" "a" "a" "g" "t" "c" "a" "c" "t"
    [2701] "t" "t" "t" "g" "a" "a" "t" "g" "t" "g" "a" "a" "c" "a" "a" "a" "a" "g"
    [2719] "g" "a" "a" "g" "a" "a" "a" "a" "t" "c" "a" "a" "g" "g" "a" "a" "a" "g"
    [2737] "a" "a" "t" "g" "a" "g" "t" "c" "t" "a" "a" "t" "a" "t" "c" "a" "a" "g"
    [2755] "c" "c" "t" "g" "t" "a" "c" "a" "g" "a" "c" "a" "g" "t" "t" "a" "a" "t"
    [2773] "a" "t" "c" "a" "c" "t" "g" "c" "a" "g" "g" "c" "t" "t" "t" "c" "c" "t"
    [2791] "g" "t" "g" "g" "t" "t" "g" "g" "t" "c" "a" "g" "a" "a" "a" "g" "a" "t"
    [2809] "a" "a" "g" "c" "c" "a" "g" "t" "t" "g" "a" "t" "a" "a" "t" "g" "c" "c"
    [2827] "a" "a" "a" "t" "g" "t" "a" "g" "t" "a" "t" "c" "a" "a" "a" "g" "g" "a"
    [2845] "g" "g" "c" "t" "c" "t" "a" "g" "g" "t" "t" "t" "t" "g" "t" "c" "t" "a"
    [2863] "t" "c" "a" "t" "c" "t" "c" "a" "g" "t" "t" "c" "a" "g" "a" "g" "g" "c"
    [2881] "a" "a" "c" "g" "a" "a" "a" "c" "t" "g" "g" "a" "c" "t" "c" "a" "t" "t"
    [2899] "a" "c" "t" "c" "c" "a" "a" "a" "t" "a" "a" "a" "c" "a" "t" "g" "g" "a"
    [2917] "c" "t" "t" "t" "t" "a" "c" "a" "a" "a" "a" "c" "c" "c" "a" "t" "a" "t"
    [2935] "c" "g" "t" "a" "t" "a" "c" "c" "a" "c" "c" "a" "c" "t" "t" "t" "t" "t"
    [2953] "c" "c" "c" "a" "t" "c" "a" "a" "g" "t" "c" "a" "t" "t" "t" "g" "t" "t"
    [2971] "a" "a" "a" "a" "c" "t" "a" "a" "a" "t" "g" "t" "a" "a" "g" "a" "a" "a"
    [2989] "a" "a" "t" "c" "t" "g" "c" "t" "a" "g" "a" "g" "g" "a" "a" "a" "a" "c"
    [3007] "t" "t" "t" "g" "a" "g" "g" "a" "a" "c" "a" "t" "t" "c" "a" "a" "t" "g"
    [3025] "t" "c" "a" "c" "c" "t" "g" "a" "a" "a" "g" "a" "g" "a" "a" "a" "t" "g"
    [3043] "g" "g" "a" "a" "a" "t" "g" "a" "g" "a" "a" "c" "a" "t" "t" "c" "c" "a"
    [3061] "a" "g" "t" "a" "c" "a" "g" "t" "g" "a" "g" "c" "a" "c" "a" "a" "t" "t"
    [3079] "a" "g" "c" "c" "g" "t" "a" "a" "t" "a" "a" "c" "a" "t" "t" "a" "g" "a"
    [3097] "g" "a" "a" "a" "a" "t" "g" "t" "t" "t" "t" "t" "a" "a" "a" "g" "a" "a"
    [3115] "g" "c" "c" "a" "g" "c" "t" "c" "a" "a" "g" "c" "a" "a" "t" "a" "t" "t"
    [3133] "a" "a" "t" "g" "a" "a" "g" "t" "a" "g" "g" "t" "t" "c" "c" "a" "g" "t"
    [3151] "a" "c" "t" "a" "a" "t" "g" "a" "a" "g" "t" "g" "g" "g" "c" "t" "c" "c"
    [3169] "a" "g" "t" "a" "t" "t" "a" "a" "t" "g" "a" "a" "a" "t" "a" "g" "g" "t"
    [3187] "t" "c" "c" "a" "g" "t" "g" "a" "t" "g" "a" "a" "a" "a" "c" "a" "t" "t"
    [3205] "c" "a" "a" "g" "c" "a" "g" "a" "a" "c" "t" "a" "g" "g" "t" "a" "g" "a"
    [3223] "a" "a" "c" "a" "g" "a" "g" "g" "g" "c" "c" "a" "a" "a" "a" "t" "t" "g"
    [3241] "a" "a" "t" "g" "c" "t" "a" "t" "g" "c" "t" "t" "a" "g" "a" "t" "t" "a"
    [3259] "g" "g" "g" "g" "t" "t" "t" "t" "g" "c" "a" "a" "c" "c" "t" "g" "a" "g"
    [3277] "g" "t" "c" "t" "a" "t" "a" "a" "a" "c" "a" "a" "a" "g" "t" "c" "t" "t"
    [3295] "c" "c" "t" "g" "g" "a" "a" "g" "t" "a" "a" "t" "t" "g" "t" "a" "a" "g"
    [3313] "c" "a" "t" "c" "c" "t" "g" "a" "a" "a" "t" "a" "a" "a" "a" "a" "a" "g"
    [3331] "c" "a" "a" "g" "a" "a" "t" "a" "t" "g" "a" "a" "g" "a" "a" "g" "t" "a"
    [3349] "g" "t" "t" "c" "a" "g" "a" "c" "t" "g" "t" "t" "a" "a" "t" "a" "c" "a"
    [3367] "g" "a" "t" "t" "t" "c" "t" "c" "t" "c" "c" "a" "t" "a" "t" "c" "t" "g"
    [3385] "a" "t" "t" "t" "c" "a" "g" "a" "t" "a" "a" "c" "t" "t" "a" "g" "a" "a"
    [3403] "c" "a" "g" "c" "c" "t" "a" "t" "g" "g" "g" "a" "a" "g" "t" "a" "g" "t"
    [3421] "c" "a" "t" "g" "c" "a" "t" "c" "t" "c" "a" "g" "g" "t" "t" "t" "g" "t"
    [3439] "t" "c" "t" "g" "a" "g" "a" "c" "a" "c" "c" "t" "g" "a" "t" "g" "a" "c"
    [3457] "c" "t" "g" "t" "t" "a" "g" "a" "t" "g" "a" "t" "g" "g" "t" "g" "a" "a"
    [3475] "a" "t" "a" "a" "a" "g" "g" "a" "a" "g" "a" "t" "a" "c" "t" "a" "g" "t"
    [3493] "t" "t" "t" "g" "c" "t" "g" "a" "a" "a" "a" "t" "g" "a" "c" "a" "t" "t"
    [3511] "a" "a" "g" "g" "a" "a" "a" "g" "t" "t" "c" "t" "g" "c" "t" "g" "t" "t"
    [3529] "t" "t" "t" "a" "g" "c" "a" "a" "a" "a" "g" "c" "g" "t" "c" "c" "a" "g"
    [3547] "a" "a" "a" "g" "g" "a" "g" "a" "g" "c" "t" "t" "a" "g" "c" "a" "g" "g"
    [3565] "a" "g" "t" "c" "c" "t" "a" "g" "c" "c" "c" "t" "t" "t" "c" "a" "c" "c"
    [3583] "c" "a" "t" "a" "c" "a" "c" "a" "t" "t" "t" "g" "g" "c" "t" "c" "a" "g"
    [3601] "g" "g" "t" "t" "a" "c" "c" "g" "a" "a" "g" "a" "g" "g" "g" "g" "c" "c"
    [3619] "a" "a" "g" "a" "a" "a" "t" "t" "a" "g" "a" "g" "t" "c" "c" "t" "c" "a"
    [3637] "g" "a" "a" "g" "a" "g" "a" "a" "c" "t" "t" "a" "t" "c" "t" "a" "g" "t"
    [3655] "g" "a" "g" "g" "a" "t" "g" "a" "a" "g" "a" "g" "c" "t" "t" "c" "c" "c"
    [3673] "t" "g" "c" "t" "t" "c" "c" "a" "a" "c" "a" "c" "t" "t" "g" "t" "t" "a"
    [3691] "t" "t" "t" "g" "g" "t" "a" "a" "a" "g" "t" "a" "a" "a" "c" "a" "a" "t"
    [3709] "a" "t" "a" "c" "c" "t" "t" "c" "t" "c" "a" "g" "t" "c" "t" "a" "c" "t"
    [3727] "a" "g" "g" "c" "a" "t" "a" "g" "c" "a" "c" "c" "g" "t" "t" "g" "c" "t"
    [3745] "a" "c" "c" "g" "a" "g" "t" "g" "t" "c" "t" "g" "t" "c" "t" "a" "a" "g"
    [3763] "a" "a" "c" "a" "c" "a" "g" "a" "g" "g" "a" "g" "a" "a" "t" "t" "t" "a"
    [3781] "t" "t" "a" "t" "c" "a" "t" "t" "g" "a" "a" "g" "a" "a" "t" "a" "g" "c"
    [3799] "t" "t" "a" "a" "a" "t" "g" "a" "c" "t" "g" "c" "a" "g" "t" "a" "a" "c"
    [3817] "c" "a" "g" "g" "t" "a" "a" "t" "a" "t" "t" "g" "g" "c" "a" "a" "a" "g"
    [3835] "g" "c" "a" "t" "c" "t" "c" "a" "g" "g" "a" "a" "c" "a" "t" "c" "a" "c"
    [3853] "c" "t" "t" "a" "g" "t" "g" "a" "g" "g" "a" "a" "a" "c" "a" "a" "a" "a"
    [3871] "t" "g" "t" "t" "c" "t" "g" "c" "t" "a" "g" "c" "t" "t" "g" "t" "t" "t"
    [3889] "t" "c" "t" "t" "c" "a" "c" "a" "g" "t" "g" "c" "a" "g" "t" "g" "a" "a"
    [3907] "t" "t" "g" "g" "a" "a" "g" "a" "c" "t" "t" "g" "a" "c" "t" "g" "c" "a"
    [3925] "a" "a" "t" "a" "c" "a" "a" "a" "c" "a" "c" "c" "c" "a" "g" "g" "a" "t"
    [3943] "c" "c" "t" "t" "t" "c" "t" "t" "g" "a" "t" "t" "g" "g" "t" "t" "c" "t"
    [3961] "t" "c" "c" "a" "a" "a" "c" "a" "a" "a" "t" "g" "a" "g" "g" "c" "a" "t"
    [3979] "c" "a" "g" "t" "c" "t" "g" "a" "a" "a" "g" "c" "c" "a" "g" "g" "g" "a"
    [3997] "g" "t" "t" "g" "g" "t" "c" "t" "g" "a" "g" "t" "g" "a" "c" "a" "a" "g"
    [4015] "g" "a" "a" "t" "t" "g" "g" "t" "t" "t" "c" "a" "g" "a" "t" "g" "a" "t"
    [4033] "g" "a" "a" "g" "a" "a" "a" "g" "a" "g" "g" "a" "a" "c" "g" "g" "g" "c"
    [4051] "t" "t" "g" "g" "a" "a" "g" "a" "a" "a" "a" "t" "a" "a" "t" "c" "a" "a"
    [4069] "g" "a" "a" "g" "a" "g" "c" "a" "a" "a" "g" "c" "a" "t" "g" "g" "a" "t"
    [4087] "t" "c" "a" "a" "a" "c" "t" "t" "a" "g" "g" "t" "g" "a" "a" "g" "c" "a"
    [4105] "g" "c" "a" "t" "c" "t" "g" "g" "g" "t" "g" "t" "g" "a" "g" "a" "g" "t"
    [4123] "g" "a" "a" "a" "c" "a" "a" "g" "c" "g" "t" "c" "t" "c" "t" "g" "a" "a"
    [4141] "g" "a" "c" "t" "g" "c" "t" "c" "a" "g" "g" "g" "c" "t" "a" "t" "c" "c"
    [4159] "t" "c" "t" "c" "a" "g" "a" "g" "t" "g" "a" "c" "a" "t" "t" "t" "t" "a"
    [4177] "a" "c" "c" "a" "c" "t" "c" "a" "g" "c" "a" "g" "a" "g" "g" "g" "a" "t"
    [4195] "a" "c" "c" "a" "t" "g" "c" "a" "a" "c" "a" "t" "a" "a" "c" "c" "t" "g"
    [4213] "a" "t" "a" "a" "a" "g" "c" "t" "c" "c" "a" "g" "c" "a" "g" "g" "a" "a"
    [4231] "a" "t" "g" "g" "c" "t" "g" "a" "a" "c" "t" "a" "g" "a" "a" "g" "c" "t"
    [4249] "g" "t" "g" "t" "t" "a" "g" "a" "a" "c" "a" "g" "c" "a" "t" "g" "g" "g"
    [4267] "a" "g" "c" "c" "a" "g" "c" "c" "t" "t" "c" "t" "a" "a" "c" "a" "g" "c"
    [4285] "t" "a" "c" "c" "c" "t" "t" "c" "c" "a" "t" "c" "a" "t" "a" "a" "g" "t"
    [4303] "g" "a" "c" "t" "c" "t" "t" "c" "t" "g" "c" "c" "c" "t" "t" "g" "a" "g"
    [4321] "g" "a" "c" "c" "t" "g" "c" "g" "a" "a" "a" "t" "c" "c" "a" "g" "a" "a"
    [4339] "c" "a" "a" "a" "g" "c" "a" "c" "a" "t" "c" "a" "g" "a" "a" "a" "a" "a"
    [4357] "g" "c" "a" "g" "t" "a" "t" "t" "a" "a" "c" "t" "t" "c" "a" "c" "a" "g"
    [4375] "a" "a" "a" "a" "g" "t" "a" "g" "t" "g" "a" "a" "t" "a" "c" "c" "c" "t"
    [4393] "a" "t" "a" "a" "g" "c" "c" "a" "g" "a" "a" "t" "c" "c" "a" "g" "a" "a"
    [4411] "g" "g" "c" "c" "t" "t" "t" "c" "t" "g" "c" "t" "g" "a" "c" "a" "a" "g"
    [4429] "t" "t" "t" "g" "a" "g" "g" "t" "g" "t" "c" "t" "g" "c" "a" "g" "a" "t"
    [4447] "a" "g" "t" "t" "c" "t" "a" "c" "c" "a" "g" "t" "a" "a" "a" "a" "a" "t"
    [4465] "a" "a" "a" "g" "a" "a" "c" "c" "a" "g" "g" "a" "g" "t" "g" "g" "a" "a"
    [4483] "a" "g" "g" "t" "c" "a" "t" "c" "c" "c" "c" "t" "t" "c" "t" "a" "a" "a"
    [4501] "t" "g" "c" "c" "c" "a" "t" "c" "a" "t" "t" "a" "g" "a" "t" "g" "a" "t"
    [4519] "a" "g" "g" "t" "g" "g" "t" "a" "c" "a" "t" "g" "c" "a" "c" "a" "g" "t"
    [4537] "t" "g" "c" "t" "c" "t" "g" "g" "g" "a" "g" "t" "c" "t" "t" "c" "a" "g"
    [4555] "a" "a" "t" "a" "g" "a" "a" "a" "c" "t" "a" "c" "c" "c" "a" "t" "c" "t"
    [4573] "c" "a" "a" "g" "a" "g" "g" "a" "g" "c" "t" "c" "a" "t" "t" "a" "a" "g"
    [4591] "g" "t" "t" "g" "t" "t" "g" "a" "t" "g" "t" "g" "g" "a" "g" "g" "a" "g"
    [4609] "c" "a" "a" "c" "a" "g" "c" "t" "g" "g" "a" "a" "g" "a" "g" "t" "c" "t"
    [4627] "g" "g" "g" "c" "c" "a" "c" "a" "c" "g" "a" "t" "t" "t" "g" "a" "c" "g"
    [4645] "g" "a" "a" "a" "c" "a" "t" "c" "t" "t" "a" "c" "t" "t" "g" "c" "c" "a"
    [4663] "a" "g" "g" "c" "a" "a" "g" "a" "t" "c" "t" "a" "g" "a" "g" "g" "g" "a"
    [4681] "a" "c" "c" "c" "c" "t" "t" "a" "c" "c" "t" "g" "g" "a" "a" "t" "c" "t"
    [4699] "g" "g" "a" "a" "t" "c" "a" "g" "c" "c" "t" "c" "t" "t" "c" "t" "c" "t"
    [4717] "g" "a" "t" "g" "a" "c" "c" "c" "t" "g" "a" "a" "t" "c" "t" "g" "a" "t"
    [4735] "c" "c" "t" "t" "c" "t" "g" "a" "a" "g" "a" "c" "a" "g" "a" "g" "c" "c"
    [4753] "c" "c" "a" "g" "a" "g" "t" "c" "a" "g" "c" "t" "c" "g" "t" "g" "t" "t"
    [4771] "g" "g" "c" "a" "a" "c" "a" "t" "a" "c" "c" "a" "t" "c" "t" "t" "c" "a"
    [4789] "a" "c" "c" "t" "c" "t" "g" "c" "a" "t" "t" "g" "a" "a" "a" "g" "t" "t"
    [4807] "c" "c" "c" "c" "a" "a" "t" "t" "g" "a" "a" "a" "g" "t" "t" "g" "c" "a"
    [4825] "g" "a" "a" "t" "c" "t" "g" "c" "c" "c" "a" "g" "a" "g" "t" "c" "c" "a"
    [4843] "g" "c" "t" "g" "c" "t" "g" "c" "t" "c" "a" "t" "a" "c" "t" "a" "c" "t"
    [4861] "g" "a" "t" "a" "c" "t" "g" "c" "t" "g" "g" "g" "t" "a" "t" "a" "a" "t"
    [4879] "g" "c" "a" "a" "t" "g" "g" "a" "a" "g" "a" "a" "a" "g" "t" "g" "t" "g"
    [4897] "a" "g" "c" "a" "g" "g" "g" "a" "g" "a" "a" "g" "c" "c" "a" "g" "a" "a"
    [4915] "t" "t" "g" "a" "c" "a" "g" "c" "t" "t" "c" "a" "a" "c" "a" "g" "a" "a"
    [4933] "a" "g" "g" "g" "t" "c" "a" "a" "c" "a" "a" "a" "a" "g" "a" "a" "t" "g"
    [4951] "t" "c" "c" "a" "t" "g" "g" "t" "g" "g" "t" "g" "t" "c" "t" "g" "g" "c"
    [4969] "c" "t" "g" "a" "c" "c" "c" "c" "a" "g" "a" "a" "g" "a" "a" "t" "t" "t"
    [4987] "a" "t" "g" "c" "t" "c" "g" "t" "g" "t" "a" "c" "a" "a" "g" "t" "t" "t"
    [5005] "g" "c" "c" "a" "g" "a" "a" "a" "a" "c" "a" "c" "c" "a" "c" "a" "t" "c"
    [5023] "a" "c" "t" "t" "t" "a" "a" "c" "t" "a" "a" "t" "c" "t" "a" "a" "t" "t"
    [5041] "a" "c" "t" "g" "a" "a" "g" "a" "g" "a" "c" "t" "a" "c" "t" "c" "a" "t"
    [5059] "g" "t" "t" "g" "t" "t" "a" "t" "g" "a" "a" "a" "a" "c" "a" "g" "a" "t"
    [5077] "g" "c" "t" "g" "a" "g" "t" "t" "t" "g" "t" "g" "t" "g" "t" "g" "a" "a"
    [5095] "c" "g" "g" "a" "c" "a" "c" "t" "g" "a" "a" "a" "t" "a" "t" "t" "t" "t"
    [5113] "c" "t" "a" "g" "g" "a" "a" "t" "t" "g" "c" "g" "g" "g" "a" "g" "g" "a"
    [5131] "a" "a" "a" "t" "g" "g" "g" "t" "a" "g" "t" "t" "a" "g" "c" "t" "a" "t"
    [5149] "t" "t" "c" "t" "g" "g" "g" "t" "g" "a" "c" "c" "c" "a" "g" "t" "c" "t"
    [5167] "a" "t" "t" "a" "a" "a" "g" "a" "a" "a" "g" "a" "a" "a" "a" "a" "t" "g"
    [5185] "c" "t" "g" "a" "a" "t" "g" "a" "g" "c" "a" "t" "g" "a" "t" "t" "t" "t"
    [5203] "g" "a" "a" "g" "t" "c" "a" "g" "a" "g" "g" "a" "g" "a" "t" "g" "t" "g"
    [5221] "g" "t" "c" "a" "a" "t" "g" "g" "a" "a" "g" "a" "a" "a" "c" "c" "a" "c"
    [5239] "c" "a" "a" "g" "g" "t" "c" "c" "a" "a" "a" "g" "c" "g" "a" "g" "c" "a"
    [5257] "a" "g" "a" "g" "a" "a" "t" "c" "c" "c" "a" "g" "g" "a" "c" "a" "g" "a"
    [5275] "a" "a" "g" "a" "t" "c" "t" "t" "c" "a" "g" "g" "g" "g" "g" "c" "t" "a"
    [5293] "g" "a" "a" "a" "t" "c" "t" "g" "t" "t" "g" "c" "t" "a" "t" "g" "g" "g"
    [5311] "c" "c" "c" "t" "t" "c" "a" "c" "c" "a" "a" "c" "a" "t" "g" "c" "c" "c"
    [5329] "a" "c" "a" "g" "a" "t" "c" "a" "a" "c" "t" "g" "g" "a" "a" "t" "g" "g"
    [5347] "a" "t" "g" "g" "t" "a" "c" "a" "g" "c" "t" "g" "t" "g" "t" "g" "g" "t"
    [5365] "g" "c" "t" "t" "c" "t" "g" "t" "g" "g" "t" "g" "a" "a" "g" "g" "a" "g"
    [5383] "c" "t" "t" "t" "c" "a" "t" "c" "a" "t" "t" "c" "a" "c" "c" "c" "t" "t"
    [5401] "g" "g" "c" "a" "c" "a" "g" "g" "t" "g" "t" "c" "c" "a" "c" "c" "c" "a"
    [5419] "a" "t" "t" "g" "t" "g" "g" "t" "t" "g" "t" "g" "c" "a" "g" "c" "c" "a"
    [5437] "g" "a" "t" "g" "c" "c" "t" "g" "g" "a" "c" "a" "g" "a" "g" "g" "a" "c"
    [5455] "a" "a" "t" "g" "g" "c" "t" "t" "c" "c" "a" "t" "g" "c" "a" "a" "t" "t"
    [5473] "g" "g" "g" "c" "a" "g" "a" "t" "g" "t" "g" "t" "g" "a" "g" "g" "c" "a"
    [5491] "c" "c" "t" "g" "t" "g" "g" "t" "g" "a" "c" "c" "c" "g" "a" "g" "a" "g"
    [5509] "t" "g" "g" "g" "t" "g" "t" "t" "g" "g" "a" "c" "a" "g" "t" "g" "t" "a"
    [5527] "g" "c" "a" "c" "t" "c" "t" "a" "c" "c" "a" "g" "t" "g" "c" "c" "a" "g"
    [5545] "g" "a" "g" "c" "t" "g" "g" "a" "c" "a" "c" "c" "t" "a" "c" "c" "t" "g"
    [5563] "a" "t" "a" "c" "c" "c" "c" "a" "g" "a" "t" "c" "c" "c" "c" "c" "a" "c"
    [5581] "a" "g" "c" "c" "a" "c" "t" "a" "c" "t" "g" "a"
    attr(,"name")
    [1] "DQ190452.BRCA1"
    attr(,"Annot")
    [1] ">DQ190452.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $DQ190453.BRCA1
       [1] "a" "t" "g" "g" "a" "t" "t" "t" "a" "t" "c" "t" "g" "c" "t" "c" "t" "t"
      [19] "c" "g" "c" "g" "t" "t" "g" "a" "a" "g" "a" "a" "g" "t" "a" "c" "a" "a"
      [37] "a" "a" "t" "g" "t" "c" "a" "t" "t" "a" "a" "t" "g" "c" "t" "a" "t" "g"
      [55] "c" "a" "g" "a" "a" "a" "a" "t" "c" "t" "t" "a" "g" "a" "g" "t" "g" "t"
      [73] "c" "c" "c" "a" "t" "c" "t" "g" "t" "c" "t" "g" "g" "a" "g" "t" "t" "g"
      [91] "a" "t" "c" "a" "a" "g" "g" "a" "a" "c" "c" "t" "g" "t" "c" "t" "c" "c"
     [109] "a" "c" "a" "a" "a" "g" "t" "g" "t" "g" "a" "c" "c" "a" "c" "a" "t" "a"
     [127] "t" "t" "t" "t" "g" "c" "a" "a" "a" "t" "t" "t" "t" "g" "c" "a" "t" "g"
     [145] "c" "t" "g" "a" "a" "a" "c" "t" "t" "c" "t" "c" "a" "a" "c" "c" "a" "g"
     [163] "a" "a" "g" "a" "a" "a" "g" "g" "g" "c" "c" "t" "t" "c" "a" "c" "a" "g"
     [181] "t" "g" "t" "c" "c" "t" "t" "t" "a" "t" "g" "t" "a" "a" "g" "a" "a" "t"
     [199] "g" "a" "t" "a" "t" "a" "a" "c" "c" "a" "a" "a" "a" "g" "g" "a" "g" "c"
     [217] "c" "t" "a" "c" "a" "a" "g" "a" "a" "a" "g" "t" "a" "c" "g" "a" "g" "a"
     [235] "t" "t" "t" "a" "g" "t" "c" "a" "a" "c" "t" "t" "g" "t" "t" "g" "a" "a"
     [253] "g" "a" "g" "c" "t" "a" "t" "t" "g" "a" "a" "a" "a" "t" "c" "a" "t" "t"
     [271] "t" "g" "t" "g" "c" "t" "t" "t" "t" "c" "a" "g" "c" "t" "t" "g" "a" "c"
     [289] "a" "c" "a" "g" "g" "t" "t" "t" "g" "g" "a" "g" "t" "a" "t" "g" "c" "a"
     [307] "a" "a" "c" "a" "g" "c" "t" "a" "t" "a" "a" "t" "t" "t" "t" "g" "c" "a"
     [325] "a" "a" "a" "a" "a" "g" "g" "a" "a" "a" "a" "t" "a" "a" "c" "t" "c" "t"
     [343] "c" "c" "t" "g" "a" "a" "c" "a" "t" "c" "t" "a" "a" "a" "a" "g" "a" "t"
     [361] "g" "a" "a" "g" "t" "t" "t" "c" "t" "a" "t" "c" "a" "t" "c" "c" "a" "a"
     [379] "a" "g" "t" "a" "t" "g" "g" "g" "c" "t" "a" "c" "a" "g" "a" "a" "a" "c"
     [397] "c" "g" "t" "g" "c" "c" "a" "a" "a" "a" "g" "a" "c" "t" "t" "c" "t" "a"
     [415] "c" "a" "g" "a" "g" "t" "g" "a" "a" "c" "c" "c" "g" "a" "a" "a" "a" "t"
     [433] "c" "c" "t" "t" "c" "c" "t" "t" "g" "c" "a" "g" "g" "a" "a" "a" "c" "c"
     [451] "a" "g" "t" "c" "t" "c" "a" "g" "t" "g" "t" "c" "c" "a" "a" "c" "t" "c"
     [469] "t" "c" "t" "a" "a" "c" "c" "t" "t" "g" "g" "a" "a" "c" "t" "g" "t" "g"
     [487] "a" "g" "a" "a" "c" "t" "c" "t" "g" "a" "g" "g" "a" "c" "a" "a" "a" "g"
     [505] "c" "a" "g" "c" "g" "g" "a" "t" "a" "c" "a" "a" "c" "c" "t" "c" "a" "a"
     [523] "a" "a" "g" "a" "c" "g" "t" "c" "t" "g" "t" "c" "t" "a" "c" "a" "t" "t"
     [541] "g" "a" "a" "t" "t" "g" "g" "g" "a" "t" "c" "t" "g" "a" "t" "t" "c" "t"
     [559] "t" "c" "t" "g" "a" "a" "g" "a" "t" "a" "c" "c" "g" "t" "t" "a" "a" "t"
     [577] "a" "a" "g" "g" "c" "a" "a" "c" "t" "t" "a" "t" "t" "g" "c" "a" "g" "t"
     [595] "g" "t" "g" "g" "g" "a" "g" "a" "t" "c" "a" "a" "g" "a" "a" "t" "t" "g"
     [613] "t" "t" "a" "c" "a" "a" "a" "t" "c" "a" "c" "c" "c" "c" "t" "c" "a" "a"
     [631] "g" "g" "a" "a" "c" "c" "a" "g" "g" "g" "a" "t" "g" "a" "a" "a" "t" "c"
     [649] "a" "g" "t" "t" "t" "g" "g" "a" "t" "t" "c" "t" "g" "c" "a" "a" "a" "a"
     [667] "a" "a" "g" "g" "c" "t" "g" "c" "t" "t" "g" "t" "g" "a" "a" "t" "t" "t"
     [685] "t" "c" "t" "g" "a" "g" "a" "c" "g" "g" "a" "t" "g" "t" "a" "a" "c" "a"
     [703] "a" "a" "t" "a" "c" "t" "g" "a" "a" "c" "a" "t" "c" "a" "t" "c" "a" "a"
     [721] "c" "c" "c" "a" "g" "t" "a" "a" "t" "a" "a" "t" "g" "a" "t" "t" "t" "g"
     [739] "a" "a" "c" "a" "c" "c" "a" "c" "t" "g" "a" "g" "a" "a" "g" "c" "g" "t"
     [757] "g" "c" "a" "g" "c" "t" "g" "a" "g" "a" "g" "g" "c" "a" "t" "c" "c" "a"
     [775] "g" "a" "a" "a" "a" "g" "t" "a" "t" "c" "a" "g" "g" "g" "t" "a" "g" "t"
     [793] "t" "c" "t" "g" "t" "t" "t" "c" "a" "a" "a" "c" "t" "t" "g" "c" "a" "t"
     [811] "g" "t" "g" "g" "a" "g" "c" "c" "a" "t" "g" "t" "g" "g" "c" "a" "c" "a"
     [829] "a" "a" "t" "a" "c" "t" "c" "a" "t" "g" "c" "c" "a" "g" "c" "t" "c" "a"
     [847] "t" "t" "a" "c" "a" "g" "c" "a" "t" "g" "a" "g" "a" "a" "c" "a" "g" "c"
     [865] "a" "g" "t" "t" "t" "a" "t" "t" "a" "c" "t" "c" "a" "c" "t" "a" "a" "a"
     [883] "g" "a" "c" "a" "g" "a" "a" "t" "g" "a" "a" "t" "g" "t" "a" "g" "a" "a"
     [901] "a" "a" "g" "g" "c" "t" "g" "a" "a" "t" "t" "c" "t" "g" "t" "a" "a" "t"
     [919] "a" "a" "a" "a" "g" "c" "a" "a" "a" "c" "a" "g" "c" "c" "t" "g" "g" "c"
     [937] "t" "t" "a" "g" "c" "a" "a" "g" "g" "a" "g" "c" "c" "a" "a" "c" "a" "t"
     [955] "a" "a" "c" "a" "g" "a" "t" "g" "g" "g" "c" "t" "g" "g" "a" "a" "g" "t"
     [973] "a" "a" "g" "g" "a" "a" "a" "c" "a" "t" "g" "t" "a" "a" "t" "g" "a" "t"
     [991] "a" "g" "g" "c" "g" "g" "a" "c" "t" "c" "c" "c" "a" "g" "c" "a" "c" "a"
    [1009] "g" "a" "a" "a" "a" "a" "a" "a" "g" "g" "t" "a" "g" "a" "t" "c" "t" "g"
    [1027] "a" "a" "t" "g" "c" "t" "g" "a" "t" "c" "c" "c" "c" "t" "g" "t" "g" "t"
    [1045] "g" "a" "g" "a" "g" "a" "a" "a" "a" "g" "a" "a" "t" "g" "g" "a" "a" "t"
    [1063] "a" "a" "g" "c" "g" "g" "a" "a" "a" "c" "t" "g" "c" "c" "a" "t" "g" "c"
    [1081] "t" "c" "a" "g" "a" "g" "a" "a" "t" "c" "c" "t" "a" "g" "a" "g" "a" "t"
    [1099] "a" "c" "t" "g" "a" "a" "g" "a" "t" "g" "t" "t" "c" "c" "t" "t" "g" "g"
    [1117] "a" "t" "a" "a" "c" "a" "c" "t" "a" "a" "a" "t" "a" "g" "c" "a" "g" "c"
    [1135] "a" "t" "t" "c" "a" "g" "a" "a" "a" "g" "t" "t" "a" "a" "t" "g" "a" "g"
    [1153] "t" "g" "g" "t" "t" "t" "t" "c" "c" "a" "g" "a" "a" "g" "t" "g" "a" "t"
    [1171] "g" "a" "a" "c" "t" "g" "t" "t" "a" "g" "g" "t" "t" "c" "t" "g" "a" "t"
    [1189] "g" "a" "c" "t" "c" "a" "c" "a" "t" "g" "a" "t" "g" "g" "g" "g" "a" "g"
    [1207] "t" "c" "t" "g" "a" "a" "t" "c" "a" "a" "a" "t" "g" "c" "c" "a" "a" "a"
    [1225] "g" "t" "a" "g" "c" "t" "g" "a" "t" "g" "t" "a" "t" "t" "g" "g" "a" "c"
    [1243] "g" "t" "t" "c" "t" "a" "a" "a" "t" "g" "a" "g" "g" "t" "a" "g" "a" "t"
    [1261] "g" "a" "a" "t" "a" "t" "t" "c" "t" "g" "g" "t" "t" "c" "t" "t" "c" "a"
    [1279] "g" "a" "g" "a" "a" "a" "a" "t" "a" "g" "a" "c" "t" "t" "a" "c" "t" "g"
    [1297] "g" "c" "c" "a" "g" "t" "g" "a" "t" "c" "c" "t" "c" "a" "t" "g" "a" "g"
    [1315] "g" "c" "t" "t" "t" "a" "a" "t" "a" "t" "g" "t" "a" "a" "a" "a" "g" "t"
    [1333] "g" "a" "a" "a" "g" "a" "g" "t" "t" "c" "a" "c" "t" "c" "c" "a" "a" "a"
    [1351] "t" "c" "a" "g" "t" "a" "g" "a" "g" "a" "g" "t" "a" "a" "t" "a" "t" "t"
    [1369] "g" "a" "a" "g" "a" "c" "a" "a" "a" "a" "t" "a" "t" "t" "t" "g" "g" "g"
    [1387] "a" "a" "a" "a" "c" "c" "t" "a" "t" "c" "g" "g" "a" "a" "g" "a" "a" "g"
    [1405] "g" "c" "a" "a" "g" "c" "c" "t" "c" "c" "c" "c" "a" "a" "c" "t" "t" "a"
    [1423] "a" "g" "c" "c" "a" "t" "g" "t" "a" "a" "c" "t" "g" "a" "a" "a" "a" "t"
    [1441] "c" "t" "a" "a" "t" "t" "a" "t" "a" "g" "g" "a" "g" "c" "a" "t" "t" "t"
    [1459] "g" "t" "t" "a" "c" "t" "g" "a" "g" "c" "c" "a" "c" "a" "g" "a" "t" "a"
    [1477] "a" "t" "a" "c" "a" "a" "g" "a" "g" "c" "g" "t" "c" "c" "c" "c" "t" "c"
    [1495] "a" "c" "a" "a" "a" "t" "a" "a" "a" "t" "t" "a" "a" "a" "g" "c" "g" "t"
    [1513] "a" "a" "a" "a" "g" "g" "a" "g" "a" "c" "c" "t" "a" "c" "a" "t" "c" "a"
    [1531] "g" "g" "c" "c" "t" "t" "c" "a" "t" "c" "c" "t" "g" "a" "g" "g" "a" "t"
    [1549] "t" "t" "t" "a" "t" "c" "a" "a" "g" "a" "a" "a" "g" "c" "a" "g" "a" "t"
    [1567] "t" "t" "g" "g" "c" "a" "g" "t" "t" "c" "a" "a" "a" "a" "g" "a" "c" "t"
    [1585] "c" "c" "t" "g" "a" "a" "a" "t" "g" "a" "t" "a" "a" "a" "t" "c" "a" "g"
    [1603] "g" "g" "a" "a" "c" "t" "a" "a" "c" "c" "a" "a" "a" "c" "g" "g" "a" "g"
    [1621] "c" "a" "g" "a" "a" "t" "g" "g" "t" "c" "a" "a" "g" "t" "g" "a" "t" "g"
    [1639] "a" "a" "t" "a" "t" "t" "a" "c" "t" "a" "a" "t" "a" "g" "t" "g" "g" "t"
    [1657] "c" "a" "t" "g" "a" "g" "a" "a" "t" "a" "a" "a" "a" "c" "a" "a" "a" "a"
    [1675] "g" "g" "t" "g" "a" "t" "t" "c" "t" "a" "t" "t" "c" "a" "g" "a" "a" "t"
    [1693] "g" "a" "g" "a" "a" "a" "a" "a" "t" "c" "c" "t" "a" "a" "c" "c" "c" "a"
    [1711] "a" "t" "a" "g" "a" "a" "t" "c" "a" "c" "t" "c" "g" "a" "a" "a" "a" "a"
    [1729] "g" "a" "a" "t" "c" "t" "g" "c" "t" "t" "t" "c" "a" "a" "a" "a" "c" "g"
    [1747] "a" "a" "a" "g" "c" "t" "g" "a" "a" "c" "c" "t" "a" "t" "a" "a" "g" "c"
    [1765] "a" "g" "c" "a" "g" "t" "a" "t" "a" "a" "g" "c" "a" "a" "t" "a" "t" "g"
    [1783] "g" "a" "a" "c" "t" "c" "g" "a" "a" "t" "t" "a" "a" "a" "t" "a" "t" "c"
    [1801] "c" "a" "c" "a" "a" "t" "t" "c" "a" "a" "a" "a" "g" "c" "a" "c" "c" "t"
    [1819] "a" "a" "a" "a" "a" "g" "a" "a" "t" "a" "g" "g" "c" "t" "g" "a" "g" "g"
    [1837] "a" "g" "g" "a" "a" "g" "t" "c" "t" "t" "c" "t" "a" "c" "c" "a" "g" "g"
    [1855] "c" "a" "t" "a" "t" "t" "c" "a" "t" "g" "c" "g" "c" "t" "t" "g" "a" "a"
    [1873] "c" "t" "a" "g" "t" "a" "g" "t" "c" "a" "g" "t" "a" "g" "a" "a" "a" "t"
    [1891] "c" "t" "a" "a" "g" "c" "c" "c" "a" "c" "c" "t" "a" "a" "t" "t" "g" "t"
    [1909] "a" "c" "t" "g" "a" "a" "t" "t" "g" "c" "a" "a" "a" "t" "t" "g" "a" "t"
    [1927] "a" "g" "t" "t" "g" "t" "t" "c" "t" "a" "g" "c" "a" "g" "t" "g" "a" "a"
    [1945] "g" "a" "g" "a" "t" "a" "a" "a" "g" "a" "a" "a" "a" "a" "a" "a" "a" "g"
    [1963] "t" "a" "c" "a" "a" "c" "c" "a" "a" "a" "t" "g" "c" "c" "a" "g" "t" "c"
    [1981] "a" "g" "g" "c" "a" "c" "a" "g" "c" "a" "g" "a" "a" "a" "c" "c" "t" "a"
    [1999] "c" "a" "a" "c" "t" "c" "a" "t" "g" "g" "a" "a" "g" "g" "t" "a" "a" "a"
    [2017] "g" "a" "a" "c" "c" "t" "g" "c" "a" "a" "c" "t" "g" "g" "a" "g" "c" "c"
    [2035] "a" "a" "g" "a" "a" "g" "a" "g" "t" "a" "a" "c" "a" "a" "g" "c" "c" "a"
    [2053] "a" "a" "t" "g" "a" "a" "c" "a" "g" "a" "c" "a" "a" "g" "t" "a" "a" "a"
    [2071] "a" "g" "a" "c" "a" "t" "g" "a" "c" "a" "g" "c" "g" "a" "t" "a" "c" "t"
    [2089] "t" "t" "c" "c" "c" "a" "g" "a" "g" "c" "t" "g" "a" "a" "g" "t" "t" "a"
    [2107] "a" "c" "a" "a" "a" "t" "g" "c" "a" "c" "c" "t" "g" "g" "t" "t" "c" "t"
    [2125] "t" "t" "t" "a" "c" "t" "a" "a" "g" "t" "g" "t" "t" "c" "a" "a" "a" "t"
    [2143] "a" "c" "c" "a" "g" "t" "g" "a" "a" "c" "t" "t" "a" "a" "a" "g" "a" "a"
    [2161] "t" "t" "t" "g" "t" "c" "a" "a" "t" "c" "c" "t" "a" "g" "c" "c" "t" "t"
    [2179] "c" "c" "a" "a" "g" "a" "g" "a" "a" "g" "a" "a" "a" "a" "a" "g" "a" "a"
    [2197] "g" "a" "g" "a" "a" "a" "c" "t" "a" "g" "a" "a" "a" "c" "a" "g" "t" "t"
    [2215] "a" "a" "a" "g" "t" "g" "t" "c" "t" "a" "a" "t" "a" "a" "t" "g" "c" "t"
    [2233] "g" "a" "a" "g" "a" "c" "c" "c" "c" "a" "a" "a" "g" "a" "t" "c" "t" "c"
    [2251] "a" "t" "g" "t" "t" "a" "a" "g" "t" "g" "g" "a" "g" "a" "a" "a" "g" "g"
    [2269] "g" "t" "t" "t" "t" "g" "c" "a" "a" "a" "c" "t" "g" "a" "a" "a" "g" "a"
    [2287] "t" "c" "t" "g" "t" "a" "g" "a" "g" "a" "g" "t" "a" "g" "c" "a" "g" "t"
    [2305] "a" "t" "t" "t" "c" "a" "t" "t" "g" "g" "t" "a" "c" "c" "t" "g" "g" "t"
    [2323] "a" "c" "t" "g" "a" "t" "t" "a" "t" "g" "g" "c" "a" "c" "t" "c" "a" "g"
    [2341] "g" "a" "a" "a" "g" "t" "a" "t" "c" "t" "c" "g" "t" "t" "a" "c" "t" "g"
    [2359] "g" "a" "a" "g" "t" "t" "a" "g" "c" "a" "c" "t" "c" "t" "a" "g" "g" "g"
    [2377] "a" "a" "g" "g" "c" "a" "a" "a" "a" "a" "c" "a" "g" "a" "a" "c" "c" "a"
    [2395] "a" "a" "t" "a" "a" "a" "t" "g" "t" "g" "t" "g" "a" "g" "t" "c" "a" "g"
    [2413] "t" "g" "t" "g" "c" "a" "g" "c" "a" "t" "t" "t" "g" "a" "a" "a" "a" "c"
    [2431] "c" "c" "c" "a" "a" "g" "g" "g" "a" "c" "t" "a" "a" "t" "t" "c" "a" "t"
    [2449] "g" "g" "t" "t" "g" "t" "t" "c" "c" "a" "a" "a" "g" "a" "t" "a" "a" "t"
    [2467] "a" "g" "a" "a" "a" "t" "g" "a" "c" "a" "c" "a" "g" "a" "a" "g" "g" "c"
    [2485] "t" "t" "t" "a" "a" "g" "t" "a" "t" "c" "c" "a" "t" "t" "g" "g" "g" "a"
    [2503] "c" "a" "t" "g" "a" "a" "g" "t" "t" "a" "a" "c" "c" "a" "c" "a" "g" "t"
    [2521] "c" "g" "g" "g" "a" "a" "a" "c" "a" "a" "g" "c" "a" "t" "a" "g" "a" "a"
    [2539] "a" "t" "g" "g" "a" "a" "g" "a" "a" "a" "g" "t" "g" "a" "a" "c" "t" "t"
    [2557] "g" "a" "t" "g" "c" "t" "c" "a" "g" "t" "a" "t" "t" "t" "g" "c" "a" "g"
    [2575] "a" "a" "t" "a" "c" "a" "t" "t" "c" "a" "a" "g" "g" "t" "t" "t" "c" "a"
    [2593] "a" "a" "g" "c" "g" "c" "c" "a" "g" "t" "c" "a" "t" "t" "t" "g" "c" "t"
    [2611] "c" "c" "g" "t" "t" "t" "t" "c" "a" "a" "a" "t" "c" "c" "a" "g" "g" "a"
    [2629] "a" "a" "t" "g" "c" "a" "g" "a" "a" "g" "a" "g" "g" "a" "a" "t" "g" "t"
    [2647] "g" "c" "a" "a" "c" "a" "t" "t" "c" "t" "c" "t" "g" "c" "c" "c" "a" "c"
    [2665] "t" "c" "t" "g" "g" "g" "t" "c" "c" "t" "t" "a" "a" "a" "g" "a" "a" "a"
    [2683] "c" "a" "a" "a" "g" "t" "c" "c" "a" "a" "a" "a" "g" "t" "c" "a" "c" "t"
    [2701] "t" "t" "t" "g" "a" "a" "t" "g" "t" "g" "a" "a" "c" "a" "a" "a" "a" "g"
    [2719] "g" "a" "a" "g" "a" "a" "a" "a" "t" "c" "a" "a" "g" "g" "a" "a" "a" "g"
    [2737] "a" "a" "t" "g" "a" "g" "t" "c" "t" "a" "a" "t" "a" "t" "c" "a" "a" "g"
    [2755] "c" "c" "t" "g" "t" "a" "c" "a" "g" "a" "c" "a" "g" "t" "t" "a" "a" "t"
    [2773] "a" "t" "c" "a" "c" "t" "g" "c" "a" "g" "g" "c" "t" "t" "t" "c" "c" "t"
    [2791] "g" "t" "g" "g" "t" "t" "g" "g" "t" "c" "a" "g" "a" "a" "a" "g" "a" "t"
    [2809] "a" "a" "g" "c" "c" "a" "g" "t" "t" "g" "a" "t" "a" "a" "t" "g" "c" "c"
    [2827] "a" "a" "a" "t" "g" "t" "a" "g" "t" "a" "t" "c" "a" "a" "a" "g" "g" "a"
    [2845] "g" "g" "c" "t" "c" "t" "a" "g" "g" "t" "t" "t" "t" "g" "t" "c" "t" "a"
    [2863] "t" "c" "a" "t" "c" "t" "c" "a" "g" "t" "t" "c" "a" "g" "a" "g" "g" "c"
    [2881] "a" "a" "c" "g" "a" "a" "a" "c" "t" "g" "g" "a" "c" "t" "c" "a" "t" "t"
    [2899] "a" "c" "t" "c" "c" "a" "a" "a" "t" "a" "a" "a" "c" "a" "t" "g" "g" "a"
    [2917] "c" "t" "t" "t" "t" "a" "c" "a" "a" "a" "a" "c" "c" "c" "a" "t" "a" "t"
    [2935] "c" "g" "t" "a" "t" "a" "c" "c" "a" "c" "c" "a" "c" "t" "t" "t" "t" "t"
    [2953] "c" "c" "c" "a" "t" "c" "a" "a" "g" "t" "c" "a" "t" "t" "t" "g" "t" "t"
    [2971] "a" "a" "a" "a" "c" "t" "a" "a" "a" "t" "g" "t" "a" "a" "g" "a" "a" "a"
    [2989] "a" "a" "t" "c" "t" "g" "c" "t" "a" "g" "a" "g" "g" "a" "a" "a" "a" "c"
    [3007] "t" "t" "t" "g" "a" "g" "g" "a" "a" "c" "a" "t" "t" "c" "a" "a" "t" "g"
    [3025] "t" "c" "a" "c" "c" "t" "g" "a" "a" "a" "g" "a" "g" "a" "a" "a" "t" "g"
    [3043] "g" "g" "a" "a" "a" "t" "g" "a" "g" "a" "a" "c" "a" "t" "t" "c" "c" "a"
    [3061] "a" "g" "t" "a" "c" "a" "g" "t" "g" "a" "g" "c" "a" "c" "a" "a" "t" "t"
    [3079] "a" "g" "c" "c" "g" "t" "a" "a" "t" "a" "a" "c" "a" "t" "t" "a" "g" "a"
    [3097] "g" "a" "a" "a" "a" "t" "g" "t" "t" "t" "t" "t" "a" "a" "a" "g" "a" "a"
    [3115] "g" "c" "c" "a" "g" "c" "t" "c" "a" "a" "g" "c" "a" "a" "t" "a" "t" "t"
    [3133] "a" "a" "t" "g" "a" "a" "g" "t" "a" "g" "g" "t" "t" "c" "c" "a" "g" "t"
    [3151] "a" "c" "t" "a" "a" "t" "g" "a" "a" "g" "t" "g" "g" "g" "c" "t" "c" "c"
    [3169] "a" "g" "t" "a" "t" "t" "a" "a" "t" "g" "a" "a" "a" "t" "a" "g" "g" "t"
    [3187] "t" "c" "c" "a" "g" "t" "g" "a" "t" "g" "a" "a" "a" "a" "c" "a" "t" "t"
    [3205] "c" "a" "a" "g" "c" "a" "g" "a" "a" "c" "t" "a" "g" "g" "t" "a" "g" "a"
    [3223] "a" "a" "c" "a" "g" "a" "g" "g" "g" "c" "c" "a" "a" "a" "a" "t" "t" "g"
    [3241] "a" "a" "t" "g" "c" "t" "a" "t" "g" "c" "t" "t" "a" "g" "a" "t" "t" "a"
    [3259] "g" "g" "g" "g" "t" "t" "t" "t" "g" "c" "a" "a" "c" "c" "t" "g" "a" "g"
    [3277] "g" "t" "c" "t" "a" "t" "a" "a" "a" "c" "a" "a" "a" "g" "t" "c" "t" "t"
    [3295] "c" "c" "t" "g" "g" "a" "a" "g" "t" "a" "a" "t" "t" "g" "t" "a" "a" "g"
    [3313] "c" "a" "t" "c" "c" "t" "g" "a" "a" "a" "t" "a" "a" "a" "a" "a" "a" "g"
    [3331] "c" "a" "a" "g" "a" "a" "t" "a" "t" "g" "a" "a" "g" "a" "a" "g" "t" "a"
    [3349] "g" "t" "t" "c" "a" "g" "a" "c" "t" "g" "t" "t" "a" "a" "t" "a" "c" "a"
    [3367] "g" "a" "t" "t" "t" "c" "t" "c" "t" "c" "c" "a" "t" "a" "t" "c" "t" "g"
    [3385] "a" "t" "t" "t" "c" "a" "g" "a" "t" "a" "a" "c" "t" "t" "a" "g" "a" "a"
    [3403] "c" "a" "g" "c" "c" "t" "a" "t" "g" "g" "g" "a" "a" "g" "t" "a" "g" "t"
    [3421] "c" "a" "t" "g" "c" "a" "t" "c" "t" "c" "a" "g" "g" "t" "t" "t" "g" "t"
    [3439] "t" "c" "t" "g" "a" "g" "a" "c" "a" "c" "c" "t" "g" "a" "t" "g" "a" "c"
    [3457] "c" "t" "g" "t" "t" "a" "g" "a" "t" "g" "a" "t" "g" "g" "t" "g" "a" "a"
    [3475] "a" "t" "a" "a" "a" "g" "g" "a" "a" "g" "a" "t" "a" "c" "t" "a" "g" "t"
    [3493] "t" "t" "t" "g" "c" "t" "g" "a" "a" "a" "a" "t" "g" "a" "c" "a" "t" "t"
    [3511] "a" "a" "g" "g" "a" "a" "a" "g" "t" "t" "c" "t" "g" "c" "t" "g" "t" "t"
    [3529] "t" "t" "t" "a" "g" "c" "a" "a" "a" "a" "g" "c" "g" "t" "c" "c" "a" "g"
    [3547] "a" "a" "a" "g" "g" "a" "g" "a" "g" "c" "t" "t" "a" "g" "c" "a" "g" "g"
    [3565] "a" "g" "t" "c" "c" "t" "a" "g" "c" "c" "c" "t" "t" "t" "c" "a" "c" "c"
    [3583] "c" "a" "t" "a" "c" "a" "c" "a" "t" "t" "t" "g" "g" "c" "t" "c" "a" "g"
    [3601] "g" "g" "t" "t" "a" "c" "c" "g" "a" "a" "g" "a" "g" "g" "g" "g" "c" "c"
    [3619] "a" "a" "g" "a" "a" "a" "t" "t" "a" "g" "a" "g" "t" "c" "c" "t" "c" "a"
    [3637] "g" "a" "a" "g" "a" "g" "a" "a" "c" "t" "t" "a" "t" "c" "t" "a" "g" "t"
    [3655] "g" "a" "g" "g" "a" "t" "g" "a" "a" "g" "a" "g" "c" "t" "t" "c" "c" "c"
    [3673] "t" "g" "c" "t" "t" "c" "c" "a" "a" "c" "a" "c" "t" "t" "g" "t" "t" "a"
    [3691] "t" "t" "t" "g" "g" "t" "a" "a" "a" "g" "t" "a" "a" "a" "c" "a" "a" "t"
    [3709] "a" "t" "a" "c" "c" "t" "t" "c" "t" "c" "a" "g" "t" "c" "t" "a" "c" "t"
    [3727] "a" "g" "g" "c" "a" "t" "a" "g" "c" "a" "c" "c" "g" "t" "t" "g" "c" "t"
    [3745] "a" "c" "c" "g" "a" "g" "t" "g" "t" "c" "t" "g" "t" "c" "t" "a" "a" "g"
    [3763] "a" "a" "c" "a" "c" "a" "g" "a" "g" "g" "a" "g" "a" "a" "t" "t" "t" "a"
    [3781] "t" "t" "a" "t" "c" "a" "t" "t" "g" "a" "a" "g" "a" "a" "t" "a" "g" "c"
    [3799] "t" "t" "a" "a" "a" "t" "g" "a" "c" "t" "g" "c" "a" "g" "t" "a" "a" "c"
    [3817] "c" "a" "g" "g" "t" "a" "a" "t" "a" "t" "t" "g" "g" "c" "a" "a" "a" "g"
    [3835] "g" "c" "a" "t" "c" "t" "c" "a" "g" "g" "a" "a" "c" "a" "t" "c" "a" "c"
    [3853] "c" "t" "t" "a" "g" "t" "g" "a" "g" "g" "a" "a" "a" "c" "a" "a" "a" "a"
    [3871] "t" "g" "t" "t" "c" "t" "g" "c" "t" "a" "g" "c" "t" "t" "g" "t" "t" "t"
    [3889] "t" "c" "t" "t" "c" "a" "c" "a" "g" "t" "g" "c" "a" "g" "t" "g" "a" "a"
    [3907] "t" "t" "g" "g" "a" "a" "g" "a" "c" "t" "t" "g" "a" "c" "t" "g" "c" "a"
    [3925] "a" "a" "t" "a" "c" "a" "a" "a" "c" "a" "c" "c" "c" "a" "g" "g" "a" "t"
    [3943] "c" "c" "t" "t" "t" "c" "t" "t" "g" "a" "t" "t" "g" "g" "t" "t" "c" "t"
    [3961] "t" "c" "c" "a" "a" "a" "c" "a" "a" "a" "t" "g" "a" "g" "g" "c" "a" "t"
    [3979] "c" "a" "g" "t" "c" "t" "g" "a" "a" "a" "g" "c" "c" "a" "g" "g" "g" "a"
    [3997] "g" "t" "t" "g" "g" "t" "c" "t" "g" "a" "g" "t" "g" "a" "c" "a" "a" "g"
    [4015] "g" "a" "a" "t" "t" "g" "g" "t" "t" "t" "c" "a" "g" "a" "t" "g" "a" "t"
    [4033] "g" "a" "a" "g" "a" "a" "a" "g" "a" "g" "g" "a" "a" "c" "g" "g" "g" "c"
    [4051] "t" "t" "g" "g" "a" "a" "g" "a" "a" "a" "a" "t" "a" "a" "t" "c" "a" "a"
    [4069] "g" "a" "a" "g" "a" "g" "c" "a" "a" "a" "g" "c" "a" "t" "g" "g" "a" "t"
    [4087] "t" "c" "a" "a" "a" "c" "t" "t" "a" "g" "g" "t" "g" "a" "a" "g" "c" "a"
    [4105] "g" "c" "a" "t" "c" "t" "g" "g" "g" "t" "g" "t" "g" "a" "g" "a" "g" "t"
    [4123] "g" "a" "a" "a" "c" "a" "a" "g" "c" "g" "t" "c" "t" "c" "t" "g" "a" "a"
    [4141] "g" "a" "c" "t" "g" "c" "t" "c" "a" "g" "g" "g" "c" "t" "a" "t" "c" "c"
    [4159] "t" "c" "t" "c" "a" "g" "a" "g" "t" "g" "a" "c" "a" "t" "t" "t" "t" "a"
    [4177] "a" "c" "c" "a" "c" "t" "c" "a" "g" "c" "a" "g" "a" "g" "g" "g" "a" "t"
    [4195] "a" "c" "c" "a" "t" "g" "c" "a" "a" "c" "a" "t" "a" "a" "c" "c" "t" "g"
    [4213] "a" "t" "a" "a" "a" "g" "c" "t" "c" "c" "a" "g" "c" "a" "g" "g" "a" "a"
    [4231] "a" "t" "g" "g" "c" "t" "g" "a" "a" "c" "t" "a" "g" "a" "a" "g" "c" "t"
    [4249] "g" "t" "g" "t" "t" "a" "g" "a" "a" "c" "a" "g" "c" "a" "t" "g" "g" "g"
    [4267] "a" "g" "c" "c" "a" "g" "c" "c" "t" "t" "c" "t" "a" "a" "c" "a" "g" "c"
    [4285] "t" "a" "c" "c" "c" "t" "t" "c" "c" "a" "t" "c" "a" "t" "a" "a" "g" "t"
    [4303] "g" "a" "c" "t" "c" "t" "t" "c" "t" "g" "c" "c" "c" "t" "t" "g" "a" "g"
    [4321] "g" "a" "c" "c" "t" "g" "c" "g" "a" "a" "a" "t" "c" "c" "a" "g" "a" "a"
    [4339] "c" "a" "a" "a" "g" "c" "a" "c" "a" "t" "c" "a" "g" "a" "a" "a" "a" "a"
    [4357] "g" "c" "a" "g" "t" "a" "t" "t" "a" "a" "c" "t" "t" "c" "a" "c" "a" "g"
    [4375] "a" "a" "a" "a" "g" "t" "a" "g" "t" "g" "a" "a" "t" "a" "c" "c" "c" "t"
    [4393] "a" "t" "a" "a" "g" "c" "c" "a" "g" "a" "a" "t" "c" "c" "a" "g" "a" "a"
    [4411] "g" "g" "c" "c" "t" "t" "t" "c" "t" "g" "c" "t" "g" "a" "c" "a" "a" "g"
    [4429] "t" "t" "t" "g" "a" "g" "g" "t" "g" "t" "c" "t" "g" "c" "a" "g" "a" "t"
    [4447] "a" "g" "t" "t" "c" "t" "a" "c" "c" "a" "g" "t" "a" "a" "a" "a" "a" "t"
    [4465] "a" "a" "a" "g" "a" "a" "c" "c" "a" "g" "g" "a" "g" "t" "g" "g" "a" "a"
    [4483] "a" "g" "g" "t" "c" "a" "t" "c" "c" "c" "c" "t" "t" "c" "t" "a" "a" "a"
    [4501] "t" "g" "c" "c" "c" "a" "t" "c" "a" "t" "t" "a" "g" "a" "t" "g" "a" "t"
    [4519] "a" "g" "g" "t" "g" "g" "t" "a" "c" "a" "t" "g" "c" "a" "c" "a" "g" "t"
    [4537] "t" "g" "c" "t" "c" "t" "g" "g" "g" "a" "g" "t" "c" "t" "t" "c" "a" "g"
    [4555] "a" "a" "t" "a" "g" "a" "a" "a" "c" "t" "a" "c" "c" "c" "a" "t" "c" "t"
    [4573] "c" "a" "a" "g" "a" "g" "g" "a" "g" "c" "t" "c" "a" "t" "t" "a" "a" "g"
    [4591] "g" "t" "t" "g" "t" "t" "g" "a" "t" "g" "t" "g" "g" "a" "g" "g" "a" "g"
    [4609] "c" "a" "a" "c" "a" "g" "c" "t" "g" "g" "a" "a" "g" "a" "g" "t" "c" "t"
    [4627] "g" "g" "g" "c" "c" "a" "c" "a" "c" "g" "a" "t" "t" "t" "g" "a" "c" "g"
    [4645] "g" "a" "a" "a" "c" "a" "t" "c" "t" "t" "a" "c" "t" "t" "g" "c" "c" "a"
    [4663] "a" "g" "g" "c" "a" "a" "g" "a" "t" "c" "t" "a" "g" "a" "g" "g" "g" "a"
    [4681] "a" "c" "c" "c" "c" "t" "t" "a" "c" "c" "t" "g" "g" "a" "a" "t" "c" "t"
    [4699] "g" "g" "a" "a" "t" "c" "a" "g" "c" "c" "t" "c" "t" "t" "c" "t" "c" "t"
    [4717] "g" "a" "t" "g" "a" "c" "c" "c" "t" "g" "a" "a" "t" "c" "t" "g" "a" "t"
    [4735] "c" "c" "t" "t" "c" "t" "g" "a" "a" "g" "a" "c" "a" "g" "a" "g" "c" "c"
    [4753] "c" "c" "a" "g" "a" "g" "t" "c" "a" "g" "c" "t" "c" "g" "t" "g" "t" "t"
    [4771] "g" "g" "c" "a" "a" "c" "a" "t" "a" "c" "c" "a" "t" "c" "t" "t" "c" "a"
    [4789] "a" "c" "c" "t" "c" "t" "g" "c" "a" "t" "t" "g" "a" "a" "a" "g" "t" "t"
    [4807] "c" "c" "c" "c" "a" "a" "t" "t" "g" "a" "a" "a" "g" "t" "t" "g" "c" "a"
    [4825] "g" "a" "a" "t" "c" "t" "g" "c" "c" "c" "a" "g" "a" "g" "t" "c" "c" "a"
    [4843] "g" "c" "t" "g" "c" "t" "g" "c" "t" "c" "a" "t" "a" "c" "t" "a" "c" "t"
    [4861] "g" "a" "t" "a" "c" "t" "g" "c" "t" "g" "g" "g" "t" "a" "t" "a" "a" "t"
    [4879] "g" "c" "a" "a" "t" "g" "g" "a" "a" "g" "a" "a" "a" "g" "t" "g" "t" "g"
    [4897] "a" "g" "c" "a" "g" "g" "g" "a" "g" "a" "a" "g" "c" "c" "a" "g" "a" "a"
    [4915] "t" "t" "g" "a" "c" "a" "g" "c" "t" "t" "c" "a" "a" "c" "a" "g" "a" "a"
    [4933] "a" "g" "g" "g" "t" "c" "a" "a" "c" "a" "a" "a" "a" "g" "a" "a" "t" "g"
    [4951] "t" "c" "c" "a" "t" "g" "g" "t" "g" "g" "t" "g" "t" "c" "t" "g" "g" "c"
    [4969] "c" "t" "g" "a" "c" "c" "c" "c" "a" "g" "a" "a" "g" "a" "a" "t" "t" "t"
    [4987] "a" "t" "g" "c" "t" "c" "g" "t" "g" "t" "a" "c" "a" "a" "g" "t" "t" "t"
    [5005] "g" "c" "c" "a" "g" "a" "a" "a" "a" "c" "a" "c" "c" "a" "c" "a" "t" "c"
    [5023] "a" "c" "t" "t" "t" "a" "a" "c" "t" "a" "a" "t" "c" "t" "a" "a" "t" "t"
    [5041] "a" "c" "t" "g" "a" "a" "g" "a" "g" "a" "c" "t" "a" "c" "t" "c" "a" "t"
    [5059] "g" "t" "t" "g" "t" "t" "a" "t" "g" "a" "a" "a" "a" "c" "a" "g" "a" "t"
    [5077] "g" "c" "t" "g" "a" "g" "t" "t" "t" "g" "t" "g" "t" "g" "t" "g" "a" "a"
    [5095] "c" "g" "g" "a" "c" "a" "c" "t" "g" "a" "a" "a" "t" "a" "t" "t" "t" "t"
    [5113] "c" "t" "a" "g" "g" "a" "a" "t" "t" "g" "c" "g" "g" "g" "a" "g" "g" "a"
    [5131] "a" "a" "a" "t" "g" "g" "g" "t" "a" "g" "t" "t" "a" "g" "c" "t" "a" "t"
    [5149] "t" "t" "c" "t" "g" "g" "g" "t" "g" "a" "c" "c" "c" "a" "g" "t" "c" "t"
    [5167] "a" "t" "t" "a" "a" "a" "g" "a" "a" "a" "g" "a" "a" "a" "a" "a" "t" "g"
    [5185] "c" "t" "g" "a" "a" "t" "g" "a" "g" "c" "a" "t" "g" "a" "t" "t" "t" "t"
    [5203] "g" "a" "a" "g" "t" "c" "a" "g" "a" "g" "g" "a" "g" "a" "t" "g" "t" "g"
    [5221] "g" "t" "c" "a" "a" "t" "g" "g" "a" "a" "g" "a" "a" "a" "c" "c" "a" "c"
    [5239] "c" "a" "a" "g" "g" "t" "c" "c" "a" "a" "a" "g" "c" "g" "a" "g" "c" "a"
    [5257] "a" "g" "a" "g" "a" "a" "t" "c" "c" "c" "a" "g" "g" "a" "c" "a" "g" "a"
    [5275] "a" "a" "g" "a" "t" "c" "t" "t" "c" "a" "g" "g" "g" "g" "g" "c" "t" "a"
    [5293] "g" "a" "a" "a" "t" "c" "t" "g" "t" "t" "g" "c" "t" "a" "t" "g" "g" "g"
    [5311] "c" "c" "c" "t" "t" "c" "a" "c" "c" "a" "a" "c" "a" "t" "g" "c" "c" "c"
    [5329] "a" "c" "a" "g" "a" "t" "c" "a" "a" "c" "t" "g" "g" "a" "a" "t" "g" "g"
    [5347] "a" "t" "g" "g" "t" "a" "c" "a" "g" "c" "t" "g" "t" "g" "t" "g" "g" "t"
    [5365] "g" "c" "t" "t" "c" "t" "g" "t" "g" "g" "t" "g" "a" "a" "g" "g" "a" "g"
    [5383] "c" "t" "t" "t" "c" "a" "t" "c" "a" "t" "t" "c" "a" "c" "c" "c" "t" "t"
    [5401] "g" "g" "c" "a" "c" "a" "g" "g" "t" "g" "t" "c" "c" "a" "c" "c" "c" "a"
    [5419] "a" "t" "t" "g" "t" "g" "g" "t" "t" "g" "t" "g" "c" "a" "g" "c" "c" "a"
    [5437] "g" "a" "t" "g" "c" "c" "t" "g" "g" "a" "c" "a" "g" "a" "g" "g" "a" "c"
    [5455] "a" "a" "t" "g" "g" "c" "t" "t" "c" "c" "a" "t" "g" "c" "a" "a" "t" "t"
    [5473] "g" "g" "g" "c" "a" "g" "a" "t" "g" "t" "g" "t" "g" "a" "g" "g" "c" "a"
    [5491] "c" "c" "t" "g" "t" "g" "g" "t" "g" "a" "c" "c" "c" "g" "a" "g" "a" "g"
    [5509] "t" "g" "g" "g" "t" "g" "t" "t" "g" "g" "a" "c" "a" "g" "t" "g" "t" "a"
    [5527] "g" "c" "a" "c" "t" "c" "t" "a" "c" "c" "a" "g" "t" "g" "c" "c" "a" "g"
    [5545] "g" "a" "g" "c" "t" "g" "g" "a" "c" "a" "c" "c" "t" "a" "c" "c" "t" "g"
    [5563] "a" "t" "a" "c" "c" "c" "c" "a" "g" "a" "t" "c" "c" "c" "c" "c" "a" "c"
    [5581] "a" "g" "c" "c" "a" "c" "t" "a" "c" "t" "g" "a"
    attr(,"name")
    [1] "DQ190453.BRCA1"
    attr(,"Annot")
    [1] ">DQ190453.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $DQ190454.BRCA1
       [1] "a" "t" "g" "g" "a" "t" "t" "t" "a" "t" "c" "t" "g" "c" "t" "c" "t" "t"
      [19] "c" "g" "c" "g" "t" "t" "g" "a" "a" "g" "a" "a" "g" "t" "a" "c" "a" "a"
      [37] "a" "a" "t" "g" "t" "c" "a" "t" "t" "a" "a" "t" "g" "c" "t" "a" "t" "g"
      [55] "c" "a" "g" "a" "a" "a" "a" "t" "c" "t" "t" "a" "g" "a" "g" "t" "g" "t"
      [73] "c" "c" "c" "a" "t" "c" "t" "g" "t" "c" "t" "g" "g" "a" "g" "t" "t" "g"
      [91] "a" "t" "c" "a" "a" "g" "g" "a" "a" "c" "c" "t" "g" "t" "c" "t" "c" "c"
     [109] "a" "c" "a" "a" "a" "g" "t" "g" "t" "g" "a" "c" "c" "a" "c" "a" "t" "a"
     [127] "t" "t" "t" "t" "g" "c" "a" "a" "a" "t" "t" "t" "t" "g" "c" "a" "t" "g"
     [145] "c" "t" "g" "a" "a" "a" "c" "t" "t" "c" "t" "c" "a" "a" "c" "c" "a" "g"
     [163] "a" "a" "g" "a" "a" "a" "g" "g" "g" "c" "c" "t" "t" "c" "a" "c" "a" "g"
     [181] "t" "g" "t" "c" "c" "t" "t" "t" "a" "t" "g" "t" "a" "a" "g" "a" "a" "t"
     [199] "g" "a" "t" "a" "t" "a" "a" "c" "c" "a" "a" "a" "a" "g" "g" "a" "g" "c"
     [217] "c" "t" "a" "c" "a" "a" "g" "a" "a" "a" "g" "t" "a" "c" "g" "a" "g" "a"
     [235] "t" "t" "t" "a" "g" "t" "c" "a" "a" "c" "t" "t" "g" "t" "t" "g" "a" "a"
     [253] "g" "a" "g" "c" "t" "a" "t" "t" "g" "a" "a" "a" "a" "t" "c" "a" "t" "t"
     [271] "t" "g" "t" "g" "c" "t" "t" "t" "t" "c" "a" "g" "c" "t" "t" "g" "a" "c"
     [289] "a" "c" "a" "g" "g" "t" "t" "t" "g" "g" "a" "g" "t" "a" "t" "g" "c" "a"
     [307] "a" "a" "c" "a" "g" "c" "t" "a" "t" "a" "a" "t" "t" "t" "t" "g" "c" "a"
     [325] "a" "a" "a" "a" "a" "g" "g" "a" "a" "a" "a" "t" "a" "a" "c" "t" "c" "t"
     [343] "c" "c" "t" "g" "a" "a" "c" "a" "t" "c" "t" "a" "a" "a" "a" "g" "a" "t"
     [361] "g" "a" "a" "g" "t" "t" "t" "c" "t" "a" "t" "c" "a" "t" "c" "c" "a" "a"
     [379] "a" "g" "t" "a" "t" "g" "g" "g" "c" "t" "a" "c" "a" "g" "a" "a" "a" "c"
     [397] "c" "g" "t" "g" "c" "c" "a" "a" "a" "a" "g" "a" "c" "t" "t" "c" "t" "a"
     [415] "c" "a" "g" "a" "g" "t" "g" "a" "a" "c" "c" "c" "g" "a" "a" "a" "a" "t"
     [433] "c" "c" "t" "t" "c" "c" "t" "t" "g" "c" "a" "g" "g" "a" "a" "a" "c" "c"
     [451] "a" "g" "t" "c" "t" "c" "a" "g" "t" "g" "t" "c" "c" "a" "a" "c" "t" "c"
     [469] "t" "c" "t" "a" "a" "c" "c" "t" "t" "g" "g" "a" "a" "c" "t" "g" "t" "g"
     [487] "a" "g" "a" "a" "c" "t" "c" "t" "g" "a" "g" "g" "a" "c" "a" "a" "a" "g"
     [505] "c" "a" "g" "c" "g" "g" "a" "t" "a" "c" "a" "a" "c" "c" "t" "c" "a" "a"
     [523] "a" "a" "g" "a" "c" "g" "t" "c" "t" "g" "t" "c" "t" "a" "c" "a" "t" "t"
     [541] "g" "a" "a" "t" "t" "g" "g" "g" "a" "t" "c" "t" "g" "a" "t" "t" "c" "t"
     [559] "t" "c" "t" "g" "a" "a" "g" "a" "t" "a" "c" "c" "g" "t" "t" "a" "a" "t"
     [577] "a" "a" "g" "g" "c" "a" "a" "c" "t" "t" "a" "t" "t" "g" "c" "a" "g" "t"
     [595] "g" "t" "g" "g" "g" "a" "g" "a" "t" "c" "a" "a" "g" "a" "a" "t" "t" "g"
     [613] "t" "t" "a" "c" "a" "a" "a" "t" "c" "a" "c" "c" "c" "c" "t" "c" "a" "a"
     [631] "g" "g" "a" "a" "c" "c" "a" "g" "g" "g" "a" "t" "g" "a" "a" "a" "t" "c"
     [649] "a" "g" "t" "t" "t" "g" "g" "a" "t" "t" "c" "t" "g" "c" "a" "a" "a" "a"
     [667] "a" "a" "g" "g" "c" "t" "g" "c" "t" "t" "g" "t" "g" "a" "a" "t" "t" "t"
     [685] "t" "c" "t" "g" "a" "g" "a" "c" "g" "g" "a" "t" "g" "t" "a" "a" "c" "a"
     [703] "a" "a" "t" "a" "c" "t" "g" "a" "a" "c" "a" "t" "c" "a" "t" "c" "a" "a"
     [721] "c" "c" "c" "a" "g" "t" "a" "a" "t" "a" "a" "t" "g" "a" "t" "t" "t" "g"
     [739] "a" "a" "c" "a" "c" "c" "a" "c" "t" "g" "a" "g" "a" "a" "g" "c" "g" "t"
     [757] "g" "c" "a" "g" "c" "t" "g" "a" "g" "a" "g" "g" "c" "a" "t" "c" "c" "a"
     [775] "g" "a" "a" "a" "a" "g" "t" "a" "t" "c" "a" "g" "g" "g" "t" "a" "g" "t"
     [793] "t" "c" "t" "g" "t" "t" "t" "c" "a" "a" "a" "c" "t" "t" "g" "c" "a" "t"
     [811] "g" "t" "g" "g" "a" "g" "c" "c" "a" "t" "g" "t" "g" "g" "c" "a" "c" "a"
     [829] "a" "a" "t" "a" "c" "t" "c" "a" "t" "g" "c" "c" "a" "g" "c" "t" "c" "a"
     [847] "t" "t" "a" "c" "a" "g" "c" "a" "t" "g" "a" "g" "a" "a" "c" "a" "g" "c"
     [865] "a" "g" "t" "t" "t" "a" "t" "t" "a" "c" "t" "c" "a" "c" "t" "a" "a" "a"
     [883] "g" "a" "c" "a" "g" "a" "a" "t" "g" "a" "a" "t" "g" "t" "a" "g" "a" "a"
     [901] "a" "a" "g" "g" "c" "t" "g" "a" "a" "t" "t" "c" "t" "g" "t" "a" "a" "t"
     [919] "a" "a" "a" "a" "g" "c" "a" "a" "a" "c" "a" "g" "c" "c" "t" "g" "g" "c"
     [937] "t" "t" "a" "g" "c" "a" "a" "g" "g" "a" "g" "c" "c" "a" "a" "c" "a" "t"
     [955] "a" "a" "c" "a" "g" "a" "t" "g" "g" "g" "c" "t" "g" "g" "a" "a" "g" "t"
     [973] "a" "a" "g" "g" "a" "a" "a" "c" "a" "t" "g" "t" "a" "a" "t" "g" "a" "t"
     [991] "a" "g" "g" "c" "g" "g" "a" "c" "t" "c" "c" "c" "a" "g" "c" "a" "c" "a"
    [1009] "g" "a" "a" "a" "a" "a" "a" "a" "g" "g" "t" "a" "g" "a" "t" "c" "t" "g"
    [1027] "a" "a" "t" "g" "c" "t" "g" "a" "t" "c" "c" "c" "c" "t" "g" "t" "g" "t"
    [1045] "g" "a" "g" "a" "g" "a" "a" "a" "a" "g" "a" "a" "t" "g" "g" "a" "a" "t"
    [1063] "a" "a" "g" "c" "a" "g" "a" "a" "a" "c" "t" "g" "c" "c" "a" "t" "g" "c"
    [1081] "t" "c" "a" "g" "a" "g" "a" "a" "t" "c" "c" "t" "a" "g" "a" "g" "a" "t"
    [1099] "a" "c" "t" "g" "a" "a" "g" "a" "t" "g" "t" "t" "c" "c" "t" "t" "g" "g"
    [1117] "a" "t" "a" "a" "c" "a" "c" "t" "a" "a" "a" "t" "a" "g" "c" "a" "g" "c"
    [1135] "a" "t" "t" "c" "a" "g" "a" "a" "a" "g" "t" "t" "a" "a" "t" "g" "a" "g"
    [1153] "t" "g" "g" "t" "t" "t" "t" "c" "c" "a" "g" "a" "a" "g" "t" "g" "a" "t"
    [1171] "g" "a" "a" "c" "t" "g" "t" "t" "a" "g" "g" "t" "t" "c" "t" "g" "a" "t"
    [1189] "g" "a" "c" "t" "c" "a" "c" "a" "t" "g" "a" "t" "g" "g" "g" "g" "a" "g"
    [1207] "t" "c" "t" "g" "a" "a" "t" "c" "a" "a" "a" "t" "g" "c" "c" "a" "a" "a"
    [1225] "g" "t" "a" "g" "c" "t" "g" "a" "t" "g" "t" "a" "t" "t" "g" "g" "a" "c"
    [1243] "g" "t" "t" "c" "t" "a" "a" "a" "t" "g" "a" "g" "g" "t" "a" "g" "a" "t"
    [1261] "g" "a" "a" "t" "a" "t" "t" "c" "t" "g" "g" "t" "t" "c" "t" "t" "c" "a"
    [1279] "g" "a" "g" "a" "a" "a" "a" "t" "a" "g" "a" "c" "t" "t" "a" "c" "t" "g"
    [1297] "g" "c" "c" "a" "g" "t" "g" "a" "t" "c" "c" "t" "c" "a" "t" "g" "a" "g"
    [1315] "g" "c" "t" "t" "t" "a" "a" "t" "a" "t" "g" "t" "a" "a" "a" "a" "g" "t"
    [1333] "g" "a" "a" "a" "g" "a" "g" "t" "t" "c" "a" "c" "t" "c" "c" "a" "a" "a"
    [1351] "t" "c" "a" "g" "t" "a" "g" "a" "g" "a" "g" "t" "a" "a" "t" "a" "t" "t"
    [1369] "g" "a" "a" "g" "a" "c" "a" "a" "a" "a" "t" "a" "t" "t" "t" "g" "g" "g"
    [1387] "a" "a" "a" "a" "c" "c" "t" "a" "t" "c" "g" "g" "a" "a" "g" "a" "a" "g"
    [1405] "g" "c" "a" "a" "g" "c" "c" "t" "c" "c" "c" "c" "a" "a" "c" "t" "t" "a"
    [1423] "a" "g" "c" "c" "a" "t" "g" "t" "a" "a" "c" "t" "g" "a" "a" "a" "a" "t"
    [1441] "c" "t" "a" "a" "t" "t" "a" "t" "a" "g" "g" "a" "g" "c" "a" "t" "t" "t"
    [1459] "g" "t" "t" "a" "c" "t" "g" "a" "g" "c" "c" "a" "c" "a" "g" "a" "t" "a"
    [1477] "a" "t" "a" "c" "a" "a" "g" "a" "g" "c" "g" "t" "c" "c" "c" "c" "t" "c"
    [1495] "a" "c" "a" "a" "a" "t" "a" "a" "a" "t" "t" "a" "a" "a" "g" "c" "g" "t"
    [1513] "a" "a" "a" "a" "g" "g" "a" "g" "a" "c" "c" "t" "a" "c" "a" "t" "c" "a"
    [1531] "g" "g" "c" "c" "t" "t" "c" "a" "t" "c" "c" "t" "g" "a" "g" "g" "a" "t"
    [1549] "t" "t" "t" "a" "t" "c" "a" "a" "g" "a" "a" "a" "g" "c" "a" "g" "a" "t"
    [1567] "t" "t" "g" "g" "c" "a" "g" "t" "t" "c" "a" "a" "a" "a" "g" "a" "c" "t"
    [1585] "c" "c" "t" "g" "a" "a" "a" "t" "g" "a" "t" "a" "a" "a" "t" "c" "a" "g"
    [1603] "g" "g" "a" "a" "c" "t" "a" "a" "c" "c" "a" "a" "a" "c" "g" "g" "a" "g"
    [1621] "c" "a" "g" "a" "a" "t" "g" "g" "t" "c" "a" "a" "g" "t" "g" "a" "t" "g"
    [1639] "a" "a" "t" "a" "t" "t" "a" "c" "t" "a" "a" "t" "a" "g" "t" "g" "g" "t"
    [1657] "c" "a" "t" "g" "a" "g" "a" "a" "t" "a" "a" "a" "a" "c" "a" "a" "a" "a"
    [1675] "g" "g" "t" "g" "a" "t" "t" "c" "t" "a" "t" "t" "c" "a" "g" "a" "a" "t"
    [1693] "g" "a" "g" "a" "a" "a" "a" "a" "t" "c" "c" "t" "a" "a" "c" "c" "c" "a"
    [1711] "a" "t" "a" "g" "a" "a" "t" "c" "a" "c" "t" "c" "g" "a" "a" "a" "a" "a"
    [1729] "g" "a" "a" "t" "c" "t" "g" "c" "t" "t" "t" "c" "a" "a" "a" "a" "c" "g"
    [1747] "a" "a" "a" "g" "c" "t" "g" "a" "a" "c" "c" "t" "a" "t" "a" "a" "g" "c"
    [1765] "a" "g" "c" "a" "g" "t" "a" "t" "a" "a" "g" "c" "a" "a" "t" "a" "t" "g"
    [1783] "g" "a" "a" "c" "t" "c" "g" "a" "a" "t" "t" "a" "a" "a" "t" "a" "t" "c"
    [1801] "c" "a" "c" "a" "a" "t" "t" "c" "a" "a" "a" "a" "g" "c" "a" "c" "c" "t"
    [1819] "a" "a" "a" "a" "a" "g" "a" "a" "t" "a" "g" "g" "c" "t" "g" "a" "g" "g"
    [1837] "a" "g" "g" "a" "a" "g" "t" "c" "t" "t" "c" "t" "a" "c" "c" "a" "g" "g"
    [1855] "c" "a" "t" "a" "t" "t" "c" "a" "t" "g" "c" "g" "c" "t" "t" "g" "a" "a"
    [1873] "c" "t" "a" "g" "t" "a" "g" "t" "c" "a" "g" "t" "a" "g" "a" "a" "a" "t"
    [1891] "c" "t" "a" "a" "g" "c" "c" "c" "a" "c" "c" "t" "a" "a" "t" "t" "g" "t"
    [1909] "a" "c" "t" "g" "a" "a" "t" "t" "g" "c" "a" "a" "a" "t" "t" "g" "a" "t"
    [1927] "a" "g" "t" "t" "g" "t" "t" "c" "t" "a" "g" "c" "a" "g" "t" "g" "a" "a"
    [1945] "g" "a" "g" "a" "t" "a" "a" "a" "g" "a" "a" "a" "a" "a" "a" "a" "a" "g"
    [1963] "t" "a" "c" "a" "a" "c" "c" "a" "a" "a" "t" "g" "c" "c" "a" "g" "t" "c"
    [1981] "a" "g" "g" "c" "a" "c" "a" "g" "c" "a" "g" "a" "a" "a" "c" "c" "t" "a"
    [1999] "c" "a" "a" "c" "t" "c" "a" "t" "g" "g" "a" "a" "g" "g" "t" "a" "a" "a"
    [2017] "g" "a" "a" "c" "c" "t" "g" "c" "a" "a" "c" "t" "g" "g" "a" "g" "c" "c"
    [2035] "a" "a" "g" "a" "a" "g" "a" "g" "t" "a" "a" "c" "a" "a" "g" "c" "c" "a"
    [2053] "a" "a" "t" "g" "a" "a" "c" "a" "g" "a" "c" "a" "a" "g" "t" "a" "a" "a"
    [2071] "a" "g" "a" "c" "a" "t" "g" "a" "c" "a" "g" "c" "g" "a" "t" "a" "c" "t"
    [2089] "t" "t" "c" "c" "c" "a" "g" "a" "g" "c" "t" "g" "a" "a" "g" "t" "t" "a"
    [2107] "a" "c" "a" "a" "a" "t" "g" "c" "a" "c" "c" "t" "g" "g" "t" "t" "c" "t"
    [2125] "t" "t" "t" "a" "c" "t" "a" "a" "g" "t" "g" "t" "t" "c" "a" "a" "a" "t"
    [2143] "a" "c" "c" "a" "g" "t" "g" "a" "a" "c" "t" "t" "a" "a" "a" "g" "a" "a"
    [2161] "t" "t" "t" "g" "t" "c" "a" "a" "t" "c" "c" "t" "a" "g" "c" "c" "t" "t"
    [2179] "c" "c" "a" "a" "g" "a" "g" "a" "a" "g" "a" "a" "a" "a" "a" "g" "a" "a"
    [2197] "g" "a" "g" "a" "a" "a" "c" "t" "a" "g" "a" "a" "a" "c" "a" "g" "t" "t"
    [2215] "a" "a" "a" "g" "t" "g" "t" "c" "t" "a" "a" "t" "a" "a" "t" "g" "c" "t"
    [2233] "g" "a" "a" "g" "a" "c" "c" "c" "c" "a" "a" "a" "g" "a" "t" "c" "t" "c"
    [2251] "a" "t" "g" "t" "t" "a" "a" "g" "t" "g" "g" "a" "g" "a" "a" "a" "g" "g"
    [2269] "g" "t" "t" "t" "t" "g" "c" "a" "a" "a" "c" "t" "g" "a" "a" "a" "g" "a"
    [2287] "t" "c" "t" "g" "t" "a" "g" "a" "g" "a" "g" "t" "a" "g" "c" "a" "g" "t"
    [2305] "a" "t" "t" "t" "c" "a" "t" "t" "g" "g" "t" "a" "c" "c" "t" "g" "g" "t"
    [2323] "a" "c" "t" "g" "a" "t" "t" "a" "t" "g" "g" "c" "a" "c" "t" "c" "a" "g"
    [2341] "g" "a" "a" "a" "g" "t" "a" "t" "c" "t" "c" "g" "t" "t" "a" "c" "t" "g"
    [2359] "g" "a" "a" "g" "t" "t" "a" "g" "c" "a" "c" "t" "c" "t" "a" "g" "g" "g"
    [2377] "a" "a" "g" "g" "c" "a" "a" "a" "a" "a" "c" "a" "g" "a" "a" "c" "c" "a"
    [2395] "a" "a" "t" "a" "a" "a" "t" "g" "t" "g" "t" "g" "a" "g" "t" "c" "a" "g"
    [2413] "t" "g" "t" "g" "c" "a" "g" "c" "a" "t" "t" "t" "g" "a" "a" "a" "a" "c"
    [2431] "c" "c" "c" "a" "a" "g" "g" "g" "a" "c" "t" "a" "a" "t" "t" "c" "a" "t"
    [2449] "g" "g" "t" "t" "g" "t" "t" "c" "c" "a" "a" "a" "g" "a" "t" "a" "a" "t"
    [2467] "a" "g" "a" "a" "a" "t" "g" "a" "c" "a" "c" "a" "g" "a" "a" "g" "g" "c"
    [2485] "t" "t" "t" "a" "a" "g" "t" "a" "t" "c" "c" "a" "t" "t" "g" "g" "g" "a"
    [2503] "c" "a" "t" "g" "a" "a" "g" "t" "t" "a" "a" "c" "c" "a" "c" "a" "g" "t"
    [2521] "c" "g" "g" "g" "a" "a" "a" "c" "a" "a" "g" "c" "a" "t" "a" "g" "a" "a"
    [2539] "a" "t" "g" "g" "a" "a" "g" "a" "a" "a" "g" "t" "g" "a" "a" "c" "t" "t"
    [2557] "g" "a" "t" "g" "c" "t" "c" "a" "g" "t" "a" "t" "t" "t" "g" "c" "a" "g"
    [2575] "a" "a" "t" "a" "c" "a" "t" "t" "c" "a" "a" "g" "g" "t" "t" "t" "c" "a"
    [2593] "a" "a" "g" "c" "g" "c" "c" "a" "g" "t" "c" "a" "t" "t" "t" "g" "c" "t"
    [2611] "c" "c" "g" "t" "t" "t" "t" "c" "a" "a" "a" "t" "c" "c" "a" "g" "g" "a"
    [2629] "a" "a" "t" "g" "c" "a" "g" "a" "a" "g" "a" "g" "g" "a" "a" "t" "g" "t"
    [2647] "g" "c" "a" "a" "c" "a" "t" "t" "c" "t" "c" "t" "g" "c" "c" "c" "a" "c"
    [2665] "t" "c" "t" "g" "g" "g" "t" "c" "c" "t" "t" "a" "a" "a" "g" "a" "a" "a"
    [2683] "c" "a" "a" "a" "g" "t" "c" "c" "a" "a" "a" "a" "g" "t" "c" "a" "c" "t"
    [2701] "t" "t" "t" "g" "a" "a" "t" "g" "t" "g" "a" "a" "c" "a" "a" "a" "a" "g"
    [2719] "g" "a" "a" "g" "a" "a" "a" "a" "t" "c" "a" "a" "g" "g" "a" "a" "a" "g"
    [2737] "a" "a" "t" "g" "a" "g" "t" "c" "t" "a" "a" "t" "a" "t" "c" "a" "a" "g"
    [2755] "c" "c" "t" "g" "t" "a" "c" "a" "g" "a" "c" "a" "g" "t" "t" "a" "a" "t"
    [2773] "a" "t" "c" "a" "c" "t" "g" "c" "a" "g" "g" "c" "t" "t" "t" "c" "c" "t"
    [2791] "g" "t" "g" "g" "t" "t" "g" "g" "t" "c" "a" "g" "a" "a" "a" "g" "a" "t"
    [2809] "a" "a" "g" "c" "c" "a" "g" "t" "t" "g" "a" "t" "a" "a" "t" "g" "c" "c"
    [2827] "a" "a" "a" "t" "g" "t" "a" "g" "t" "a" "t" "c" "a" "a" "a" "g" "g" "a"
    [2845] "g" "g" "c" "t" "c" "t" "a" "g" "g" "t" "t" "t" "t" "g" "t" "c" "t" "a"
    [2863] "t" "c" "a" "t" "c" "t" "c" "a" "g" "t" "t" "c" "a" "g" "a" "g" "g" "c"
    [2881] "a" "a" "c" "g" "a" "a" "a" "c" "t" "g" "g" "a" "c" "t" "c" "a" "t" "t"
    [2899] "a" "c" "t" "c" "c" "a" "a" "a" "t" "a" "a" "a" "c" "a" "t" "g" "g" "a"
    [2917] "c" "t" "t" "t" "t" "a" "c" "a" "a" "a" "a" "c" "c" "c" "a" "t" "a" "t"
    [2935] "c" "g" "t" "a" "t" "a" "c" "c" "a" "c" "c" "a" "c" "t" "t" "t" "t" "t"
    [2953] "c" "c" "c" "a" "t" "c" "a" "a" "g" "t" "c" "a" "t" "t" "t" "g" "t" "t"
    [2971] "a" "a" "a" "a" "c" "t" "a" "a" "a" "t" "g" "t" "a" "a" "g" "a" "a" "a"
    [2989] "a" "a" "t" "c" "t" "g" "c" "t" "a" "g" "a" "g" "g" "a" "a" "a" "a" "c"
    [3007] "t" "t" "t" "g" "a" "g" "g" "a" "a" "c" "a" "t" "t" "c" "a" "a" "t" "g"
    [3025] "t" "c" "a" "c" "c" "t" "g" "a" "a" "a" "g" "a" "g" "a" "a" "a" "t" "g"
    [3043] "g" "g" "a" "a" "a" "t" "g" "a" "g" "a" "a" "c" "a" "t" "t" "c" "c" "a"
    [3061] "a" "g" "t" "a" "c" "a" "g" "t" "g" "a" "g" "c" "a" "c" "a" "a" "t" "t"
    [3079] "a" "g" "c" "c" "g" "t" "a" "a" "t" "a" "a" "c" "a" "t" "t" "a" "g" "a"
    [3097] "g" "a" "a" "a" "a" "t" "g" "t" "t" "t" "t" "t" "a" "a" "a" "g" "a" "a"
    [3115] "g" "c" "c" "a" "g" "c" "t" "c" "a" "a" "g" "c" "a" "a" "t" "a" "t" "t"
    [3133] "a" "a" "t" "g" "a" "a" "g" "t" "a" "g" "g" "t" "t" "c" "c" "a" "g" "t"
    [3151] "a" "c" "t" "a" "a" "t" "g" "a" "a" "g" "t" "g" "g" "g" "c" "t" "c" "c"
    [3169] "a" "g" "t" "a" "t" "t" "a" "a" "t" "g" "a" "a" "a" "t" "a" "g" "g" "t"
    [3187] "t" "c" "c" "a" "g" "t" "g" "a" "t" "g" "a" "a" "a" "a" "c" "a" "t" "t"
    [3205] "c" "a" "a" "g" "c" "a" "g" "a" "a" "c" "t" "a" "g" "g" "t" "a" "g" "a"
    [3223] "a" "a" "c" "a" "g" "a" "g" "g" "g" "c" "c" "a" "a" "a" "a" "t" "t" "g"
    [3241] "a" "a" "t" "g" "c" "t" "a" "t" "g" "c" "t" "t" "a" "g" "a" "t" "t" "a"
    [3259] "g" "g" "g" "g" "t" "t" "t" "t" "g" "c" "a" "a" "c" "c" "t" "g" "a" "g"
    [3277] "g" "t" "c" "t" "a" "t" "a" "a" "a" "c" "a" "a" "a" "g" "t" "c" "t" "t"
    [3295] "c" "c" "t" "g" "g" "a" "a" "g" "t" "a" "a" "t" "t" "g" "t" "a" "a" "g"
    [3313] "c" "a" "t" "c" "c" "t" "g" "a" "a" "a" "t" "a" "a" "a" "a" "a" "a" "g"
    [3331] "c" "a" "a" "g" "a" "a" "t" "a" "t" "g" "a" "a" "g" "a" "a" "g" "t" "a"
    [3349] "g" "t" "t" "c" "a" "g" "a" "c" "t" "g" "t" "t" "a" "a" "t" "a" "c" "a"
    [3367] "g" "a" "t" "t" "t" "c" "t" "c" "t" "c" "c" "a" "t" "a" "t" "c" "t" "g"
    [3385] "a" "t" "t" "t" "c" "a" "g" "a" "t" "a" "a" "c" "t" "t" "a" "g" "a" "a"
    [3403] "c" "a" "g" "c" "c" "t" "a" "t" "g" "g" "g" "a" "a" "g" "t" "a" "g" "t"
    [3421] "c" "a" "t" "g" "c" "a" "t" "c" "t" "c" "a" "g" "g" "t" "t" "t" "g" "t"
    [3439] "t" "c" "t" "g" "a" "g" "a" "c" "a" "c" "c" "t" "g" "a" "t" "g" "a" "c"
    [3457] "c" "t" "g" "t" "t" "a" "g" "a" "t" "g" "a" "t" "g" "g" "t" "g" "a" "a"
    [3475] "a" "t" "a" "a" "a" "g" "g" "a" "a" "g" "a" "t" "a" "c" "t" "a" "g" "t"
    [3493] "t" "t" "t" "g" "c" "t" "g" "a" "a" "a" "a" "t" "g" "a" "c" "a" "t" "t"
    [3511] "a" "a" "g" "g" "a" "a" "a" "g" "t" "t" "c" "t" "g" "c" "t" "g" "t" "t"
    [3529] "t" "t" "t" "a" "g" "c" "a" "a" "a" "a" "g" "c" "g" "t" "c" "c" "a" "g"
    [3547] "a" "a" "a" "g" "g" "a" "g" "a" "g" "c" "t" "t" "a" "g" "c" "a" "g" "g"
    [3565] "a" "g" "t" "c" "c" "t" "a" "g" "c" "c" "c" "t" "t" "t" "c" "a" "c" "c"
    [3583] "c" "a" "t" "a" "c" "a" "c" "a" "t" "t" "t" "g" "g" "c" "t" "c" "a" "g"
    [3601] "g" "g" "t" "t" "a" "c" "c" "g" "a" "a" "g" "a" "g" "g" "g" "g" "c" "c"
    [3619] "a" "a" "g" "a" "a" "a" "t" "t" "a" "g" "a" "g" "t" "c" "c" "t" "c" "a"
    [3637] "g" "a" "a" "g" "a" "g" "a" "a" "c" "t" "t" "a" "t" "c" "t" "a" "g" "t"
    [3655] "g" "a" "g" "g" "a" "t" "g" "a" "a" "g" "a" "g" "c" "t" "t" "c" "c" "c"
    [3673] "t" "g" "c" "t" "t" "c" "c" "a" "a" "c" "a" "c" "t" "t" "g" "t" "t" "a"
    [3691] "t" "t" "t" "g" "g" "t" "a" "a" "a" "g" "t" "a" "a" "a" "c" "a" "a" "t"
    [3709] "a" "t" "a" "c" "c" "t" "t" "c" "t" "c" "a" "g" "t" "c" "t" "a" "c" "t"
    [3727] "a" "g" "g" "c" "a" "t" "a" "g" "c" "a" "c" "c" "g" "t" "t" "g" "c" "t"
    [3745] "a" "c" "c" "g" "a" "g" "t" "g" "t" "c" "t" "g" "t" "c" "t" "a" "a" "g"
    [3763] "a" "a" "c" "a" "c" "a" "g" "a" "g" "g" "a" "g" "a" "a" "t" "t" "t" "a"
    [3781] "t" "t" "a" "t" "c" "a" "t" "t" "g" "a" "a" "g" "a" "a" "t" "a" "g" "c"
    [3799] "t" "t" "a" "a" "a" "t" "g" "a" "c" "t" "g" "c" "a" "g" "t" "a" "a" "c"
    [3817] "c" "a" "g" "g" "t" "a" "a" "t" "a" "t" "t" "g" "g" "c" "a" "a" "a" "g"
    [3835] "g" "c" "a" "t" "c" "t" "c" "a" "g" "g" "a" "a" "c" "a" "t" "c" "a" "c"
    [3853] "c" "t" "t" "a" "g" "t" "g" "a" "g" "g" "a" "a" "a" "c" "a" "a" "a" "a"
    [3871] "t" "g" "t" "t" "c" "t" "g" "c" "t" "a" "g" "c" "t" "t" "g" "t" "t" "t"
    [3889] "t" "c" "t" "t" "c" "a" "c" "a" "g" "t" "g" "c" "a" "g" "t" "g" "a" "a"
    [3907] "t" "t" "g" "g" "a" "a" "g" "a" "c" "t" "t" "g" "a" "c" "t" "g" "c" "a"
    [3925] "a" "a" "t" "a" "c" "a" "a" "a" "c" "a" "c" "c" "c" "a" "g" "g" "a" "t"
    [3943] "c" "c" "t" "t" "t" "c" "t" "t" "g" "a" "t" "t" "g" "g" "t" "t" "c" "t"
    [3961] "t" "c" "c" "a" "a" "a" "c" "a" "a" "a" "t" "g" "a" "g" "g" "c" "a" "t"
    [3979] "c" "a" "g" "t" "c" "t" "g" "a" "a" "a" "g" "c" "c" "a" "g" "g" "g" "a"
    [3997] "g" "t" "t" "g" "g" "t" "c" "t" "g" "a" "g" "t" "g" "a" "c" "a" "a" "g"
    [4015] "g" "a" "a" "t" "t" "g" "g" "t" "t" "t" "c" "a" "g" "a" "t" "g" "a" "t"
    [4033] "g" "a" "a" "g" "a" "a" "a" "g" "a" "g" "g" "a" "a" "c" "g" "g" "g" "c"
    [4051] "t" "t" "g" "g" "a" "a" "g" "a" "a" "a" "a" "t" "a" "a" "t" "c" "a" "a"
    [4069] "g" "a" "a" "g" "a" "g" "c" "a" "a" "a" "g" "c" "a" "t" "g" "g" "a" "t"
    [4087] "t" "c" "a" "a" "a" "c" "t" "t" "a" "g" "g" "t" "g" "a" "a" "g" "c" "a"
    [4105] "g" "c" "a" "t" "c" "t" "g" "g" "g" "t" "g" "t" "g" "a" "g" "a" "g" "t"
    [4123] "g" "a" "a" "a" "c" "a" "a" "g" "c" "g" "t" "c" "t" "c" "t" "g" "a" "a"
    [4141] "g" "a" "c" "t" "g" "c" "t" "c" "a" "g" "g" "g" "c" "t" "a" "t" "c" "c"
    [4159] "t" "c" "t" "c" "a" "g" "a" "g" "t" "g" "a" "c" "a" "t" "t" "t" "t" "a"
    [4177] "a" "c" "c" "a" "c" "t" "c" "a" "g" "c" "a" "g" "a" "g" "g" "g" "a" "t"
    [4195] "a" "c" "c" "a" "t" "g" "c" "a" "a" "c" "a" "t" "a" "a" "c" "c" "t" "g"
    [4213] "a" "t" "a" "a" "a" "g" "c" "t" "c" "c" "a" "g" "c" "a" "g" "g" "a" "a"
    [4231] "a" "t" "g" "g" "c" "t" "g" "a" "a" "c" "t" "a" "g" "a" "a" "g" "c" "t"
    [4249] "g" "t" "g" "t" "t" "a" "g" "a" "a" "c" "a" "g" "c" "a" "t" "g" "g" "g"
    [4267] "a" "g" "c" "c" "a" "g" "c" "c" "t" "t" "c" "t" "a" "a" "c" "a" "g" "c"
    [4285] "t" "a" "c" "c" "c" "t" "t" "c" "c" "a" "t" "c" "a" "t" "a" "a" "g" "t"
    [4303] "g" "a" "c" "t" "c" "t" "t" "c" "t" "g" "c" "c" "c" "t" "t" "g" "a" "g"
    [4321] "g" "a" "c" "c" "t" "g" "c" "g" "a" "a" "a" "t" "c" "c" "a" "g" "a" "a"
    [4339] "c" "a" "a" "a" "g" "c" "a" "c" "a" "t" "c" "a" "g" "a" "a" "a" "a" "a"
    [4357] "g" "c" "a" "g" "t" "a" "t" "t" "a" "a" "c" "t" "t" "c" "a" "c" "a" "g"
    [4375] "a" "a" "a" "a" "g" "t" "a" "g" "t" "g" "a" "a" "t" "a" "c" "c" "c" "t"
    [4393] "a" "t" "a" "a" "g" "c" "c" "a" "g" "a" "a" "t" "c" "c" "a" "g" "a" "a"
    [4411] "g" "g" "c" "c" "t" "t" "t" "c" "t" "g" "c" "t" "g" "a" "c" "a" "a" "g"
    [4429] "t" "t" "t" "g" "a" "g" "g" "t" "g" "t" "c" "t" "g" "c" "a" "g" "a" "t"
    [4447] "a" "g" "t" "t" "c" "t" "a" "c" "c" "a" "g" "t" "a" "a" "a" "a" "a" "t"
    [4465] "a" "a" "a" "g" "a" "a" "c" "c" "a" "g" "g" "a" "g" "t" "g" "g" "a" "a"
    [4483] "a" "g" "g" "t" "c" "a" "t" "c" "c" "c" "c" "t" "t" "c" "t" "a" "a" "a"
    [4501] "t" "g" "c" "c" "c" "a" "t" "c" "a" "t" "t" "a" "g" "a" "t" "g" "a" "t"
    [4519] "a" "g" "g" "t" "g" "g" "t" "a" "c" "a" "t" "g" "c" "a" "c" "a" "g" "t"
    [4537] "t" "g" "c" "t" "c" "t" "g" "g" "g" "a" "g" "t" "c" "t" "t" "c" "a" "g"
    [4555] "a" "a" "t" "a" "g" "a" "a" "a" "c" "t" "a" "c" "c" "c" "a" "t" "c" "t"
    [4573] "c" "a" "a" "g" "a" "g" "g" "a" "g" "c" "t" "c" "a" "t" "t" "a" "a" "g"
    [4591] "g" "t" "t" "g" "t" "t" "g" "a" "t" "g" "t" "g" "g" "a" "g" "g" "a" "g"
    [4609] "c" "a" "a" "c" "a" "g" "c" "t" "g" "g" "a" "a" "g" "a" "g" "t" "c" "t"
    [4627] "g" "g" "g" "c" "c" "a" "c" "a" "c" "g" "a" "t" "t" "t" "g" "a" "c" "g"
    [4645] "g" "a" "a" "a" "c" "a" "t" "c" "t" "t" "a" "c" "t" "t" "g" "c" "c" "a"
    [4663] "a" "g" "g" "c" "a" "a" "g" "a" "t" "c" "t" "a" "g" "a" "g" "g" "g" "a"
    [4681] "a" "c" "c" "c" "c" "t" "t" "a" "c" "c" "t" "g" "g" "a" "a" "t" "c" "t"
    [4699] "g" "g" "a" "a" "t" "c" "a" "g" "c" "c" "t" "c" "t" "t" "c" "t" "c" "t"
    [4717] "g" "a" "t" "g" "a" "c" "c" "c" "t" "g" "a" "a" "t" "c" "t" "g" "a" "t"
    [4735] "c" "c" "t" "t" "c" "t" "g" "a" "a" "g" "a" "c" "a" "g" "a" "g" "c" "c"
    [4753] "c" "c" "a" "g" "a" "g" "t" "c" "a" "g" "c" "t" "c" "g" "t" "g" "t" "t"
    [4771] "g" "g" "c" "a" "a" "c" "a" "t" "a" "c" "c" "a" "t" "c" "t" "t" "c" "a"
    [4789] "a" "c" "c" "t" "c" "t" "g" "c" "a" "t" "t" "g" "a" "a" "a" "g" "t" "t"
    [4807] "c" "c" "c" "c" "a" "a" "t" "t" "g" "a" "a" "a" "g" "t" "t" "g" "c" "a"
    [4825] "g" "a" "a" "t" "c" "t" "g" "c" "c" "c" "a" "g" "a" "g" "t" "c" "c" "a"
    [4843] "g" "c" "t" "g" "c" "t" "g" "c" "t" "c" "a" "t" "a" "c" "t" "a" "c" "t"
    [4861] "g" "a" "t" "a" "c" "t" "g" "c" "t" "g" "g" "g" "t" "a" "t" "a" "a" "t"
    [4879] "g" "c" "a" "a" "t" "g" "g" "a" "a" "g" "a" "a" "a" "g" "t" "g" "t" "g"
    [4897] "a" "g" "c" "a" "g" "g" "g" "a" "g" "a" "a" "g" "c" "c" "a" "g" "a" "a"
    [4915] "t" "t" "g" "a" "c" "a" "g" "c" "t" "t" "c" "a" "a" "c" "a" "g" "a" "a"
    [4933] "a" "g" "g" "g" "t" "c" "a" "a" "c" "a" "a" "a" "a" "g" "a" "a" "t" "g"
    [4951] "t" "c" "c" "a" "t" "g" "g" "t" "g" "g" "t" "g" "t" "c" "t" "g" "g" "c"
    [4969] "c" "t" "g" "a" "c" "c" "c" "c" "a" "g" "a" "a" "g" "a" "a" "t" "t" "t"
    [4987] "a" "t" "g" "c" "t" "c" "g" "t" "g" "t" "a" "c" "a" "a" "g" "t" "t" "t"
    [5005] "g" "c" "c" "a" "g" "a" "a" "a" "a" "c" "a" "c" "c" "a" "c" "a" "t" "c"
    [5023] "a" "c" "t" "t" "t" "a" "a" "c" "t" "a" "a" "t" "c" "t" "a" "a" "t" "t"
    [5041] "a" "c" "t" "g" "a" "a" "g" "a" "g" "a" "c" "t" "a" "c" "t" "c" "a" "t"
    [5059] "g" "t" "t" "g" "t" "t" "a" "t" "g" "a" "a" "a" "a" "c" "a" "g" "a" "t"
    [5077] "g" "c" "t" "g" "a" "g" "t" "t" "t" "g" "t" "g" "t" "g" "t" "g" "a" "a"
    [5095] "c" "g" "g" "a" "c" "a" "c" "t" "g" "a" "a" "a" "t" "a" "t" "t" "t" "t"
    [5113] "c" "t" "a" "g" "g" "a" "a" "t" "t" "g" "c" "g" "g" "g" "a" "g" "g" "a"
    [5131] "a" "a" "a" "t" "g" "g" "g" "t" "a" "g" "t" "t" "a" "g" "c" "t" "a" "t"
    [5149] "t" "t" "c" "t" "g" "g" "g" "t" "g" "a" "c" "c" "c" "a" "g" "t" "c" "t"
    [5167] "a" "t" "t" "a" "a" "a" "g" "a" "a" "a" "g" "a" "a" "a" "a" "a" "t" "g"
    [5185] "c" "t" "g" "a" "a" "t" "g" "a" "g" "c" "a" "t" "g" "a" "t" "t" "t" "t"
    [5203] "g" "a" "a" "g" "t" "c" "a" "g" "a" "g" "g" "a" "g" "a" "t" "g" "t" "g"
    [5221] "g" "t" "c" "a" "a" "t" "g" "g" "a" "a" "g" "a" "a" "a" "c" "c" "a" "c"
    [5239] "c" "a" "a" "g" "g" "t" "c" "c" "a" "a" "a" "g" "c" "g" "a" "g" "c" "a"
    [5257] "a" "g" "a" "g" "a" "a" "t" "c" "c" "c" "a" "g" "g" "a" "c" "a" "g" "a"
    [5275] "a" "a" "g" "a" "t" "c" "t" "t" "c" "a" "g" "g" "g" "g" "g" "c" "t" "a"
    [5293] "g" "a" "a" "a" "t" "c" "t" "g" "t" "t" "g" "c" "t" "a" "t" "g" "g" "g"
    [5311] "c" "c" "c" "t" "t" "c" "a" "c" "c" "a" "a" "c" "a" "t" "g" "c" "c" "c"
    [5329] "a" "c" "a" "g" "a" "t" "c" "a" "a" "c" "t" "g" "g" "a" "a" "t" "g" "g"
    [5347] "a" "t" "g" "g" "t" "a" "c" "a" "g" "c" "t" "g" "t" "g" "t" "g" "g" "t"
    [5365] "g" "c" "t" "t" "c" "t" "g" "t" "g" "g" "t" "g" "a" "a" "g" "g" "a" "g"
    [5383] "c" "t" "t" "t" "c" "a" "t" "c" "a" "t" "t" "c" "a" "c" "c" "c" "t" "t"
    [5401] "g" "g" "c" "a" "c" "a" "g" "g" "t" "g" "t" "c" "c" "a" "c" "c" "c" "a"
    [5419] "a" "t" "t" "g" "t" "g" "g" "t" "t" "g" "t" "g" "c" "a" "g" "c" "c" "a"
    [5437] "g" "a" "t" "g" "c" "c" "t" "g" "g" "a" "c" "a" "g" "a" "g" "g" "a" "c"
    [5455] "a" "a" "t" "g" "g" "c" "t" "t" "c" "c" "a" "t" "g" "c" "a" "a" "t" "t"
    [5473] "g" "g" "g" "c" "a" "g" "a" "t" "g" "t" "g" "t" "g" "a" "g" "g" "c" "a"
    [5491] "c" "c" "t" "g" "t" "g" "g" "t" "g" "a" "c" "c" "c" "g" "a" "g" "a" "g"
    [5509] "t" "g" "g" "g" "t" "g" "t" "t" "g" "g" "a" "c" "a" "g" "t" "g" "t" "a"
    [5527] "g" "c" "a" "c" "t" "c" "t" "a" "c" "c" "a" "g" "t" "g" "c" "c" "a" "g"
    [5545] "g" "a" "g" "c" "t" "g" "g" "a" "c" "a" "c" "c" "t" "a" "c" "c" "t" "g"
    [5563] "a" "t" "a" "c" "c" "c" "c" "a" "g" "a" "t" "c" "c" "c" "c" "c" "a" "c"
    [5581] "a" "g" "c" "c" "a" "c" "t" "a" "c" "t" "g" "a"
    attr(,"name")
    [1] "DQ190454.BRCA1"
    attr(,"Annot")
    [1] ">DQ190454.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $DQ190455.BRCA1
       [1] "a" "t" "g" "g" "a" "t" "t" "t" "a" "t" "c" "t" "g" "c" "t" "c" "t" "t"
      [19] "c" "g" "c" "g" "t" "t" "g" "a" "a" "g" "a" "a" "g" "t" "a" "c" "a" "a"
      [37] "a" "a" "t" "g" "t" "c" "a" "t" "t" "a" "a" "t" "g" "c" "t" "a" "t" "g"
      [55] "c" "a" "g" "a" "a" "a" "a" "t" "c" "t" "t" "a" "g" "a" "g" "t" "g" "t"
      [73] "c" "c" "c" "a" "t" "c" "t" "g" "t" "c" "t" "g" "g" "a" "g" "t" "t" "g"
      [91] "a" "t" "c" "a" "a" "g" "g" "a" "a" "c" "c" "t" "g" "t" "c" "t" "c" "c"
     [109] "a" "c" "a" "a" "a" "g" "t" "g" "t" "g" "a" "c" "c" "a" "c" "a" "t" "a"
     [127] "t" "t" "t" "t" "g" "c" "a" "a" "a" "t" "t" "t" "t" "g" "c" "a" "t" "g"
     [145] "c" "t" "g" "a" "a" "a" "c" "t" "t" "c" "t" "c" "a" "a" "c" "c" "a" "g"
     [163] "a" "a" "g" "a" "a" "a" "g" "g" "g" "c" "c" "t" "t" "c" "a" "c" "a" "g"
     [181] "t" "g" "t" "c" "c" "t" "t" "t" "a" "t" "g" "t" "a" "a" "g" "a" "a" "t"
     [199] "g" "a" "t" "a" "t" "a" "a" "c" "c" "a" "a" "a" "a" "g" "g" "a" "g" "c"
     [217] "c" "t" "a" "c" "a" "a" "g" "a" "a" "a" "g" "t" "a" "c" "g" "a" "g" "a"
     [235] "t" "t" "t" "a" "g" "t" "c" "a" "a" "c" "t" "t" "g" "t" "t" "g" "a" "a"
     [253] "g" "a" "g" "c" "t" "a" "t" "t" "g" "a" "a" "a" "a" "t" "c" "a" "t" "t"
     [271] "t" "g" "t" "g" "c" "t" "t" "t" "t" "c" "a" "g" "c" "t" "t" "g" "a" "c"
     [289] "a" "c" "a" "g" "g" "t" "t" "t" "g" "g" "a" "g" "t" "a" "t" "g" "c" "a"
     [307] "a" "a" "c" "a" "g" "c" "t" "a" "t" "a" "a" "t" "t" "t" "t" "g" "c" "a"
     [325] "a" "a" "a" "a" "a" "g" "g" "a" "a" "a" "a" "t" "a" "a" "c" "t" "c" "t"
     [343] "c" "c" "t" "g" "a" "a" "c" "a" "t" "c" "t" "a" "a" "a" "a" "g" "a" "t"
     [361] "g" "a" "a" "g" "t" "t" "t" "c" "t" "a" "t" "c" "a" "t" "c" "c" "a" "a"
     [379] "a" "g" "t" "a" "t" "g" "g" "g" "c" "t" "a" "c" "a" "g" "a" "a" "a" "c"
     [397] "c" "g" "t" "g" "c" "c" "a" "a" "a" "a" "g" "a" "c" "t" "t" "c" "t" "a"
     [415] "c" "a" "g" "a" "g" "t" "g" "a" "a" "c" "c" "c" "g" "a" "a" "a" "a" "t"
     [433] "c" "c" "t" "t" "c" "c" "t" "t" "g" "c" "a" "g" "g" "a" "a" "a" "c" "c"
     [451] "a" "g" "t" "c" "t" "c" "a" "g" "t" "g" "t" "c" "c" "a" "a" "c" "t" "c"
     [469] "t" "c" "t" "a" "a" "c" "c" "t" "t" "g" "g" "a" "a" "c" "t" "g" "t" "g"
     [487] "a" "g" "a" "a" "c" "t" "c" "t" "g" "a" "g" "g" "a" "c" "a" "a" "a" "g"
     [505] "c" "a" "g" "c" "g" "g" "a" "t" "a" "c" "a" "a" "c" "c" "t" "c" "a" "a"
     [523] "a" "a" "g" "a" "c" "g" "t" "c" "t" "g" "t" "c" "t" "a" "c" "a" "t" "t"
     [541] "g" "a" "a" "t" "t" "g" "g" "g" "a" "t" "c" "t" "g" "a" "t" "t" "c" "t"
     [559] "t" "c" "t" "g" "a" "a" "g" "a" "t" "a" "c" "c" "g" "t" "t" "a" "a" "t"
     [577] "a" "a" "g" "g" "c" "a" "a" "c" "t" "t" "a" "t" "t" "g" "c" "a" "g" "t"
     [595] "g" "t" "g" "g" "g" "a" "g" "a" "t" "c" "a" "a" "g" "a" "a" "t" "t" "g"
     [613] "t" "t" "a" "c" "a" "a" "a" "t" "c" "a" "c" "c" "c" "c" "t" "c" "a" "a"
     [631] "g" "g" "a" "a" "c" "c" "a" "g" "g" "g" "a" "t" "g" "a" "a" "a" "t" "c"
     [649] "a" "g" "t" "t" "t" "g" "g" "a" "t" "t" "c" "t" "g" "c" "a" "a" "a" "a"
     [667] "a" "a" "g" "g" "c" "t" "g" "c" "t" "t" "g" "t" "g" "a" "a" "t" "t" "t"
     [685] "t" "c" "t" "g" "a" "g" "a" "c" "g" "g" "a" "t" "g" "t" "a" "a" "c" "a"
     [703] "a" "a" "t" "a" "c" "t" "g" "a" "a" "c" "a" "t" "c" "a" "t" "c" "a" "a"
     [721] "c" "c" "c" "a" "g" "t" "a" "a" "t" "a" "a" "t" "g" "a" "t" "t" "t" "g"
     [739] "a" "a" "c" "a" "c" "c" "a" "c" "t" "g" "a" "g" "a" "a" "g" "c" "g" "t"
     [757] "g" "c" "a" "g" "c" "t" "g" "a" "g" "a" "g" "g" "c" "a" "t" "c" "c" "a"
     [775] "g" "a" "a" "a" "a" "g" "t" "a" "t" "c" "a" "g" "g" "g" "t" "a" "g" "t"
     [793] "t" "c" "t" "g" "t" "t" "t" "c" "a" "a" "a" "c" "t" "t" "g" "c" "a" "t"
     [811] "g" "t" "g" "g" "a" "g" "c" "c" "a" "t" "g" "t" "g" "g" "c" "a" "c" "a"
     [829] "a" "a" "t" "a" "c" "t" "c" "a" "t" "g" "c" "c" "a" "g" "c" "t" "c" "a"
     [847] "t" "t" "a" "c" "a" "g" "c" "a" "t" "g" "a" "g" "a" "a" "c" "a" "g" "c"
     [865] "a" "g" "t" "t" "t" "a" "t" "t" "a" "c" "t" "c" "a" "c" "t" "a" "a" "a"
     [883] "g" "a" "c" "a" "g" "a" "a" "t" "g" "a" "a" "t" "g" "t" "a" "g" "a" "a"
     [901] "a" "a" "g" "g" "c" "t" "g" "a" "a" "t" "t" "c" "t" "g" "t" "a" "a" "t"
     [919] "a" "a" "a" "a" "g" "c" "a" "a" "a" "c" "a" "g" "c" "c" "t" "g" "g" "c"
     [937] "t" "t" "a" "g" "c" "a" "a" "g" "g" "a" "g" "c" "c" "a" "a" "c" "a" "t"
     [955] "a" "a" "c" "a" "g" "a" "t" "g" "g" "g" "c" "t" "g" "g" "a" "a" "g" "t"
     [973] "a" "a" "g" "g" "a" "a" "a" "c" "a" "t" "g" "t" "a" "a" "t" "g" "a" "t"
     [991] "a" "g" "g" "c" "g" "g" "a" "c" "t" "c" "c" "c" "a" "g" "c" "a" "c" "a"
    [1009] "g" "a" "a" "a" "a" "a" "a" "a" "g" "g" "t" "a" "g" "a" "t" "c" "t" "g"
    [1027] "a" "a" "t" "g" "c" "t" "g" "a" "t" "c" "c" "c" "c" "t" "g" "t" "g" "t"
    [1045] "g" "a" "g" "a" "g" "a" "a" "a" "a" "g" "a" "a" "t" "g" "g" "a" "a" "t"
    [1063] "a" "a" "g" "c" "a" "g" "a" "a" "a" "c" "t" "g" "c" "c" "a" "t" "g" "c"
    [1081] "t" "c" "a" "g" "a" "g" "a" "a" "t" "c" "c" "t" "a" "g" "a" "g" "a" "t"
    [1099] "a" "c" "t" "g" "a" "a" "g" "a" "t" "g" "t" "t" "c" "c" "t" "t" "g" "g"
    [1117] "a" "t" "a" "a" "c" "a" "c" "t" "a" "a" "a" "t" "a" "g" "c" "a" "g" "c"
    [1135] "a" "t" "t" "c" "a" "g" "a" "a" "a" "g" "t" "t" "a" "a" "t" "g" "a" "g"
    [1153] "t" "g" "g" "t" "t" "t" "t" "c" "c" "a" "g" "a" "a" "g" "t" "g" "a" "t"
    [1171] "g" "a" "a" "c" "t" "g" "t" "t" "a" "g" "g" "t" "t" "c" "t" "g" "a" "t"
    [1189] "g" "a" "c" "t" "c" "a" "c" "a" "t" "g" "a" "t" "g" "g" "g" "g" "a" "g"
    [1207] "t" "c" "t" "g" "a" "a" "t" "c" "a" "a" "a" "t" "g" "c" "c" "a" "a" "a"
    [1225] "g" "t" "a" "g" "c" "t" "g" "a" "t" "g" "t" "a" "t" "t" "g" "g" "a" "c"
    [1243] "g" "t" "t" "c" "t" "a" "a" "a" "t" "g" "a" "g" "g" "t" "a" "g" "a" "t"
    [1261] "g" "a" "a" "t" "a" "t" "t" "c" "t" "g" "g" "t" "t" "c" "t" "t" "c" "a"
    [1279] "g" "a" "g" "a" "a" "a" "a" "t" "a" "g" "a" "c" "t" "t" "a" "c" "t" "g"
    [1297] "g" "c" "c" "a" "g" "t" "g" "a" "t" "c" "c" "t" "c" "a" "t" "g" "a" "g"
    [1315] "g" "c" "t" "t" "t" "a" "a" "t" "a" "t" "g" "t" "a" "a" "a" "a" "g" "t"
    [1333] "g" "a" "a" "a" "g" "a" "g" "t" "t" "c" "a" "c" "t" "c" "c" "a" "a" "a"
    [1351] "t" "c" "a" "g" "t" "a" "g" "a" "g" "a" "g" "t" "a" "a" "t" "a" "t" "t"
    [1369] "g" "a" "a" "g" "a" "c" "a" "a" "a" "a" "t" "a" "t" "t" "t" "g" "g" "g"
    [1387] "a" "a" "a" "a" "c" "c" "t" "a" "t" "c" "g" "g" "a" "a" "g" "a" "a" "g"
    [1405] "g" "c" "a" "a" "g" "c" "c" "t" "c" "c" "c" "c" "a" "a" "c" "t" "t" "a"
    [1423] "a" "g" "c" "c" "a" "t" "g" "t" "a" "a" "c" "t" "g" "a" "a" "a" "a" "t"
    [1441] "c" "t" "a" "a" "t" "t" "a" "t" "a" "g" "g" "a" "g" "c" "a" "t" "t" "t"
    [1459] "g" "t" "t" "a" "c" "t" "g" "a" "g" "c" "c" "a" "c" "a" "g" "a" "t" "a"
    [1477] "a" "t" "a" "c" "a" "a" "g" "a" "g" "c" "g" "t" "c" "c" "c" "c" "t" "c"
    [1495] "a" "c" "a" "a" "a" "t" "a" "a" "a" "t" "t" "a" "a" "a" "g" "c" "g" "t"
    [1513] "a" "a" "a" "a" "g" "g" "a" "g" "a" "c" "c" "t" "a" "c" "a" "t" "c" "a"
    [1531] "g" "g" "c" "c" "t" "t" "c" "a" "t" "c" "c" "t" "g" "a" "g" "g" "a" "t"
    [1549] "t" "t" "t" "a" "t" "c" "a" "a" "g" "a" "a" "a" "g" "c" "a" "g" "a" "t"
    [1567] "t" "t" "g" "g" "c" "a" "g" "t" "t" "c" "a" "a" "a" "a" "g" "a" "c" "t"
    [1585] "c" "c" "t" "g" "a" "a" "a" "t" "g" "a" "t" "a" "a" "a" "t" "c" "a" "g"
    [1603] "g" "g" "a" "a" "c" "t" "a" "a" "c" "c" "a" "a" "a" "c" "g" "g" "a" "g"
    [1621] "c" "a" "g" "a" "a" "t" "g" "g" "t" "c" "a" "a" "g" "t" "g" "a" "t" "g"
    [1639] "a" "a" "t" "a" "t" "t" "a" "c" "t" "a" "a" "t" "a" "g" "t" "g" "g" "t"
    [1657] "c" "a" "t" "g" "a" "g" "a" "a" "t" "a" "a" "a" "a" "c" "a" "a" "a" "a"
    [1675] "g" "g" "t" "g" "a" "t" "t" "c" "t" "a" "t" "t" "c" "a" "g" "a" "a" "t"
    [1693] "g" "a" "g" "a" "a" "a" "a" "a" "t" "c" "c" "t" "a" "a" "c" "c" "c" "a"
    [1711] "a" "t" "a" "g" "a" "a" "t" "c" "a" "c" "t" "c" "g" "a" "a" "a" "a" "a"
    [1729] "g" "a" "a" "t" "c" "t" "g" "c" "t" "t" "t" "c" "a" "a" "a" "a" "c" "g"
    [1747] "a" "a" "a" "g" "c" "t" "g" "a" "a" "c" "c" "t" "a" "t" "a" "a" "g" "c"
    [1765] "a" "g" "c" "a" "g" "t" "a" "t" "a" "a" "g" "c" "a" "a" "t" "a" "t" "g"
    [1783] "g" "a" "a" "c" "t" "c" "g" "a" "a" "t" "t" "a" "a" "a" "t" "a" "t" "c"
    [1801] "c" "a" "c" "a" "a" "t" "t" "c" "a" "a" "a" "a" "g" "c" "a" "c" "c" "t"
    [1819] "a" "a" "a" "a" "a" "g" "a" "a" "t" "a" "g" "g" "c" "t" "g" "a" "g" "g"
    [1837] "a" "g" "g" "a" "a" "g" "t" "c" "t" "t" "c" "t" "a" "c" "c" "a" "g" "g"
    [1855] "c" "a" "t" "a" "t" "t" "c" "a" "t" "g" "c" "g" "c" "t" "t" "g" "a" "a"
    [1873] "c" "t" "a" "g" "t" "a" "g" "t" "c" "a" "g" "t" "a" "g" "a" "a" "a" "t"
    [1891] "c" "t" "a" "a" "g" "c" "c" "c" "a" "c" "c" "t" "a" "a" "t" "t" "g" "t"
    [1909] "a" "c" "t" "g" "a" "a" "t" "t" "g" "c" "a" "a" "a" "t" "t" "g" "a" "t"
    [1927] "a" "g" "t" "t" "g" "t" "t" "c" "t" "a" "g" "c" "a" "g" "t" "g" "a" "a"
    [1945] "g" "a" "g" "a" "t" "a" "a" "a" "g" "a" "a" "a" "a" "a" "a" "a" "a" "g"
    [1963] "t" "a" "c" "a" "a" "c" "c" "a" "a" "a" "t" "g" "c" "c" "a" "g" "t" "c"
    [1981] "a" "g" "g" "c" "a" "c" "a" "g" "c" "a" "g" "a" "a" "a" "c" "c" "t" "a"
    [1999] "c" "a" "a" "c" "t" "c" "a" "t" "g" "g" "a" "a" "g" "g" "t" "a" "a" "a"
    [2017] "g" "a" "a" "c" "c" "t" "g" "c" "a" "a" "c" "t" "g" "g" "a" "g" "c" "c"
    [2035] "a" "a" "g" "a" "a" "g" "a" "g" "t" "a" "a" "c" "a" "a" "g" "c" "c" "a"
    [2053] "a" "a" "t" "g" "a" "a" "c" "a" "g" "a" "c" "a" "a" "g" "t" "a" "a" "a"
    [2071] "a" "g" "a" "c" "a" "t" "g" "a" "c" "a" "g" "c" "g" "a" "t" "a" "c" "t"
    [2089] "t" "t" "c" "c" "c" "a" "g" "a" "g" "c" "t" "g" "a" "a" "g" "t" "t" "a"
    [2107] "a" "c" "a" "a" "a" "t" "g" "c" "a" "c" "c" "t" "g" "g" "t" "t" "c" "t"
    [2125] "t" "t" "t" "a" "c" "t" "a" "a" "g" "t" "g" "t" "t" "c" "a" "a" "a" "t"
    [2143] "a" "c" "c" "a" "g" "t" "g" "a" "a" "c" "t" "t" "a" "a" "a" "g" "a" "a"
    [2161] "t" "t" "t" "g" "t" "c" "a" "a" "t" "c" "c" "t" "a" "g" "c" "c" "t" "t"
    [2179] "c" "c" "a" "a" "g" "a" "g" "a" "a" "g" "a" "a" "a" "a" "a" "g" "a" "a"
    [2197] "g" "a" "g" "a" "a" "a" "c" "t" "a" "g" "a" "a" "a" "c" "a" "g" "t" "t"
    [2215] "a" "a" "a" "g" "t" "g" "t" "c" "t" "a" "a" "t" "a" "a" "t" "g" "c" "t"
    [2233] "g" "a" "a" "g" "a" "c" "c" "c" "c" "a" "a" "a" "g" "a" "t" "c" "t" "c"
    [2251] "a" "t" "g" "t" "t" "a" "a" "g" "t" "g" "g" "a" "g" "a" "a" "a" "g" "g"
    [2269] "g" "t" "t" "t" "t" "g" "c" "a" "a" "a" "c" "t" "g" "a" "a" "a" "g" "a"
    [2287] "t" "c" "t" "g" "t" "a" "g" "a" "g" "a" "g" "t" "a" "g" "c" "a" "g" "t"
    [2305] "a" "t" "t" "t" "c" "a" "t" "t" "g" "g" "t" "a" "c" "c" "t" "g" "g" "t"
    [2323] "a" "c" "t" "g" "a" "t" "t" "a" "t" "g" "g" "c" "a" "c" "t" "c" "a" "g"
    [2341] "g" "a" "a" "a" "g" "t" "a" "t" "c" "t" "c" "g" "t" "t" "a" "c" "t" "g"
    [2359] "g" "a" "a" "g" "t" "t" "a" "g" "c" "a" "c" "t" "c" "t" "a" "g" "g" "g"
    [2377] "a" "a" "g" "g" "c" "a" "a" "a" "a" "a" "c" "a" "g" "a" "a" "c" "c" "a"
    [2395] "a" "a" "t" "a" "a" "a" "t" "g" "t" "g" "t" "g" "a" "g" "t" "c" "a" "g"
    [2413] "t" "g" "t" "g" "c" "a" "g" "c" "a" "t" "t" "t" "g" "a" "a" "a" "a" "c"
    [2431] "c" "c" "c" "a" "a" "g" "g" "g" "a" "c" "t" "a" "a" "t" "t" "c" "a" "t"
    [2449] "g" "g" "t" "t" "g" "t" "t" "c" "c" "a" "a" "a" "g" "a" "t" "a" "a" "t"
    [2467] "a" "g" "a" "a" "a" "t" "g" "a" "c" "a" "c" "a" "g" "a" "a" "g" "g" "c"
    [2485] "t" "t" "t" "a" "a" "g" "t" "a" "t" "c" "c" "a" "t" "t" "g" "g" "g" "a"
    [2503] "c" "a" "t" "g" "a" "a" "g" "t" "t" "a" "a" "c" "c" "a" "c" "a" "g" "t"
    [2521] "c" "g" "g" "g" "a" "a" "a" "c" "a" "a" "g" "c" "a" "t" "a" "g" "a" "a"
    [2539] "a" "t" "g" "g" "a" "a" "g" "a" "a" "a" "g" "t" "g" "a" "a" "c" "t" "t"
    [2557] "g" "a" "t" "g" "c" "t" "c" "a" "g" "t" "a" "t" "t" "t" "g" "c" "a" "g"
    [2575] "a" "a" "t" "a" "c" "a" "t" "t" "c" "a" "a" "g" "g" "t" "t" "t" "c" "a"
    [2593] "a" "a" "g" "c" "g" "c" "c" "a" "g" "t" "c" "a" "t" "t" "t" "g" "c" "t"
    [2611] "c" "c" "g" "t" "t" "t" "t" "c" "a" "a" "a" "t" "c" "c" "a" "g" "g" "a"
    [2629] "a" "a" "t" "g" "c" "a" "g" "a" "a" "g" "a" "g" "g" "a" "a" "t" "g" "t"
    [2647] "g" "c" "a" "a" "c" "a" "t" "t" "c" "t" "c" "t" "g" "c" "c" "c" "a" "c"
    [2665] "t" "c" "t" "g" "g" "g" "t" "c" "c" "t" "t" "a" "a" "a" "g" "a" "a" "a"
    [2683] "c" "a" "a" "a" "g" "t" "c" "c" "a" "a" "a" "a" "g" "t" "c" "a" "c" "t"
    [2701] "t" "t" "t" "g" "a" "a" "t" "g" "t" "g" "a" "a" "c" "a" "a" "a" "a" "g"
    [2719] "g" "a" "a" "g" "a" "a" "a" "a" "t" "c" "a" "a" "g" "g" "a" "a" "a" "g"
    [2737] "a" "a" "t" "g" "a" "g" "t" "c" "t" "a" "a" "t" "a" "t" "c" "a" "a" "g"
    [2755] "c" "c" "t" "g" "t" "a" "c" "a" "g" "a" "c" "a" "g" "t" "t" "a" "a" "t"
    [2773] "a" "t" "c" "a" "c" "t" "g" "c" "a" "g" "g" "c" "t" "t" "t" "c" "c" "t"
    [2791] "g" "t" "g" "g" "t" "t" "g" "g" "t" "c" "a" "g" "a" "a" "a" "g" "a" "t"
    [2809] "a" "a" "g" "c" "c" "a" "g" "t" "t" "g" "a" "t" "a" "a" "t" "g" "c" "c"
    [2827] "a" "a" "a" "t" "g" "t" "a" "g" "t" "a" "t" "c" "a" "a" "a" "g" "g" "a"
    [2845] "g" "g" "c" "t" "c" "t" "a" "g" "g" "t" "t" "t" "t" "g" "t" "c" "t" "a"
    [2863] "t" "c" "a" "t" "c" "t" "c" "a" "g" "t" "t" "c" "a" "g" "a" "g" "g" "c"
    [2881] "a" "a" "c" "g" "a" "a" "a" "c" "t" "g" "g" "a" "c" "t" "c" "a" "t" "t"
    [2899] "a" "c" "t" "c" "c" "a" "a" "a" "t" "a" "a" "a" "c" "a" "t" "g" "g" "a"
    [2917] "c" "t" "t" "t" "t" "a" "c" "a" "a" "a" "a" "c" "c" "c" "a" "t" "a" "t"
    [2935] "c" "g" "t" "a" "t" "a" "c" "c" "a" "c" "c" "a" "c" "t" "t" "t" "t" "t"
    [2953] "c" "c" "c" "a" "t" "c" "a" "a" "g" "t" "c" "a" "t" "t" "t" "g" "t" "t"
    [2971] "a" "a" "a" "a" "c" "t" "a" "a" "a" "t" "g" "t" "a" "a" "g" "a" "a" "a"
    [2989] "a" "a" "t" "c" "t" "g" "c" "t" "a" "g" "a" "g" "g" "a" "a" "a" "a" "c"
    [3007] "t" "t" "t" "g" "a" "g" "g" "a" "a" "c" "a" "t" "t" "c" "a" "a" "t" "g"
    [3025] "t" "c" "a" "c" "c" "t" "g" "a" "a" "a" "g" "a" "g" "a" "a" "a" "t" "g"
    [3043] "g" "g" "a" "a" "a" "t" "g" "a" "g" "a" "a" "c" "a" "t" "t" "c" "c" "a"
    [3061] "a" "g" "t" "a" "c" "a" "g" "t" "g" "a" "g" "c" "a" "c" "a" "a" "t" "t"
    [3079] "a" "g" "c" "c" "g" "t" "a" "a" "t" "a" "a" "c" "a" "t" "t" "a" "g" "a"
    [3097] "g" "a" "a" "a" "a" "t" "g" "t" "t" "t" "t" "t" "a" "a" "a" "g" "a" "a"
    [3115] "g" "c" "c" "a" "g" "c" "t" "c" "a" "a" "g" "c" "a" "a" "t" "a" "t" "t"
    [3133] "a" "a" "t" "g" "a" "a" "g" "t" "a" "g" "g" "t" "t" "c" "c" "a" "g" "t"
    [3151] "a" "c" "t" "a" "a" "t" "g" "a" "a" "g" "t" "g" "g" "g" "c" "t" "c" "c"
    [3169] "a" "g" "t" "a" "t" "t" "a" "a" "t" "g" "a" "a" "a" "t" "a" "g" "g" "t"
    [3187] "t" "c" "c" "a" "g" "t" "g" "a" "t" "g" "a" "a" "a" "a" "c" "a" "t" "t"
    [3205] "c" "a" "a" "g" "c" "a" "g" "a" "a" "c" "t" "a" "g" "g" "t" "a" "g" "a"
    [3223] "a" "a" "c" "a" "g" "a" "g" "g" "g" "c" "c" "a" "a" "a" "a" "t" "t" "g"
    [3241] "a" "a" "t" "g" "c" "t" "a" "t" "g" "c" "t" "t" "a" "g" "a" "t" "t" "a"
    [3259] "g" "g" "g" "g" "t" "t" "t" "t" "g" "c" "a" "a" "c" "c" "t" "g" "a" "g"
    [3277] "g" "t" "c" "t" "a" "t" "a" "a" "a" "c" "a" "a" "a" "g" "t" "c" "t" "t"
    [3295] "c" "c" "t" "g" "g" "a" "a" "g" "t" "a" "a" "t" "t" "g" "t" "a" "a" "g"
    [3313] "c" "a" "t" "c" "c" "t" "g" "a" "a" "a" "t" "a" "a" "a" "a" "a" "a" "g"
    [3331] "c" "a" "a" "g" "a" "a" "t" "a" "t" "g" "a" "a" "g" "a" "a" "g" "t" "a"
    [3349] "g" "t" "t" "c" "a" "g" "a" "c" "t" "g" "t" "t" "a" "a" "t" "a" "c" "a"
    [3367] "g" "a" "t" "t" "t" "c" "t" "c" "t" "c" "c" "a" "t" "a" "t" "c" "t" "g"
    [3385] "a" "t" "t" "t" "c" "a" "g" "a" "t" "a" "a" "c" "t" "t" "a" "g" "a" "a"
    [3403] "c" "a" "g" "c" "c" "t" "a" "t" "g" "g" "g" "a" "a" "g" "t" "a" "g" "t"
    [3421] "c" "a" "t" "g" "c" "a" "t" "c" "t" "c" "a" "g" "g" "t" "t" "t" "g" "t"
    [3439] "t" "c" "t" "g" "a" "g" "a" "c" "a" "c" "c" "t" "g" "a" "t" "g" "a" "c"
    [3457] "c" "t" "g" "t" "t" "a" "g" "a" "t" "g" "a" "t" "g" "g" "t" "g" "a" "a"
    [3475] "a" "t" "a" "a" "a" "g" "g" "a" "a" "g" "a" "t" "a" "c" "t" "a" "g" "t"
    [3493] "t" "t" "t" "g" "c" "t" "g" "a" "a" "a" "a" "t" "g" "a" "c" "a" "t" "t"
    [3511] "a" "a" "g" "g" "a" "a" "a" "g" "t" "t" "c" "t" "g" "c" "t" "g" "t" "t"
    [3529] "t" "t" "t" "a" "g" "c" "a" "a" "a" "a" "g" "c" "g" "t" "c" "c" "a" "g"
    [3547] "a" "a" "a" "g" "g" "a" "g" "a" "g" "c" "t" "t" "a" "g" "c" "a" "g" "g"
    [3565] "a" "g" "t" "c" "c" "t" "a" "g" "c" "c" "c" "t" "t" "t" "c" "a" "c" "c"
    [3583] "c" "a" "t" "a" "c" "a" "c" "a" "t" "t" "t" "g" "g" "c" "t" "c" "a" "g"
    [3601] "g" "g" "t" "t" "a" "c" "c" "g" "a" "a" "g" "a" "g" "g" "g" "g" "c" "c"
    [3619] "a" "a" "g" "a" "a" "a" "t" "t" "a" "g" "a" "g" "t" "c" "c" "t" "c" "a"
    [3637] "g" "a" "a" "g" "a" "g" "a" "a" "c" "t" "t" "a" "t" "c" "t" "a" "g" "t"
    [3655] "g" "a" "g" "g" "a" "t" "g" "a" "a" "g" "a" "g" "c" "t" "t" "c" "c" "c"
    [3673] "t" "g" "c" "t" "t" "c" "c" "a" "a" "c" "a" "c" "t" "t" "g" "t" "t" "a"
    [3691] "t" "t" "t" "g" "g" "t" "a" "a" "a" "g" "t" "a" "a" "a" "c" "a" "a" "t"
    [3709] "a" "t" "a" "c" "c" "t" "t" "c" "t" "c" "a" "g" "t" "c" "t" "a" "c" "t"
    [3727] "a" "g" "g" "c" "a" "t" "a" "g" "c" "a" "c" "c" "g" "t" "t" "g" "c" "t"
    [3745] "a" "c" "c" "g" "a" "g" "t" "g" "t" "c" "t" "g" "t" "c" "t" "a" "a" "g"
    [3763] "a" "a" "c" "a" "c" "a" "g" "a" "g" "g" "a" "g" "a" "a" "t" "t" "t" "a"
    [3781] "t" "t" "a" "t" "c" "a" "t" "t" "g" "a" "a" "g" "a" "a" "t" "a" "g" "c"
    [3799] "t" "t" "a" "a" "a" "t" "g" "a" "c" "t" "g" "c" "a" "g" "t" "a" "a" "c"
    [3817] "c" "a" "g" "g" "t" "a" "a" "t" "a" "t" "t" "g" "g" "c" "a" "a" "a" "g"
    [3835] "g" "c" "a" "t" "c" "t" "c" "a" "g" "g" "a" "a" "c" "a" "t" "c" "a" "c"
    [3853] "c" "t" "t" "a" "g" "t" "g" "a" "g" "g" "a" "a" "a" "c" "a" "a" "a" "a"
    [3871] "t" "g" "t" "t" "c" "t" "g" "c" "t" "a" "g" "c" "t" "t" "g" "t" "t" "t"
    [3889] "t" "c" "t" "t" "c" "a" "c" "a" "g" "t" "g" "c" "a" "g" "t" "g" "a" "a"
    [3907] "t" "t" "g" "g" "a" "a" "g" "a" "c" "t" "t" "g" "a" "c" "t" "g" "c" "a"
    [3925] "a" "a" "t" "a" "c" "a" "a" "a" "c" "a" "c" "c" "c" "a" "g" "g" "a" "t"
    [3943] "c" "c" "t" "t" "t" "c" "t" "t" "g" "a" "t" "t" "g" "g" "t" "t" "c" "t"
    [3961] "t" "c" "c" "a" "a" "a" "c" "a" "a" "a" "t" "g" "a" "g" "g" "c" "a" "t"
    [3979] "c" "a" "g" "t" "c" "t" "g" "a" "a" "a" "g" "c" "c" "a" "g" "g" "g" "a"
    [3997] "g" "t" "t" "g" "g" "t" "c" "t" "g" "a" "g" "t" "g" "a" "c" "a" "a" "g"
    [4015] "g" "a" "a" "t" "t" "g" "g" "t" "t" "t" "c" "a" "g" "a" "t" "g" "a" "t"
    [4033] "g" "a" "a" "g" "a" "a" "a" "g" "a" "g" "g" "a" "a" "c" "g" "g" "g" "c"
    [4051] "t" "t" "g" "g" "a" "a" "g" "a" "a" "a" "a" "t" "a" "a" "t" "c" "a" "a"
    [4069] "g" "a" "a" "g" "a" "g" "c" "a" "a" "a" "g" "c" "a" "t" "g" "g" "a" "t"
    [4087] "t" "c" "a" "a" "a" "c" "t" "t" "a" "g" "g" "t" "g" "a" "a" "g" "c" "a"
    [4105] "g" "c" "a" "t" "c" "t" "g" "g" "g" "t" "g" "t" "g" "a" "g" "a" "g" "t"
    [4123] "g" "a" "a" "a" "c" "a" "a" "g" "c" "g" "t" "c" "t" "c" "t" "g" "a" "a"
    [4141] "g" "a" "c" "t" "g" "c" "t" "c" "a" "g" "g" "g" "c" "t" "a" "t" "c" "c"
    [4159] "t" "c" "t" "c" "a" "g" "a" "g" "t" "g" "a" "c" "a" "t" "t" "t" "t" "a"
    [4177] "a" "c" "c" "a" "c" "t" "c" "a" "g" "c" "a" "g" "a" "g" "g" "g" "a" "t"
    [4195] "a" "c" "c" "a" "t" "g" "c" "a" "a" "c" "a" "t" "a" "a" "c" "c" "t" "g"
    [4213] "a" "t" "a" "a" "a" "g" "c" "t" "c" "c" "a" "g" "c" "a" "g" "g" "a" "a"
    [4231] "a" "t" "g" "g" "c" "t" "g" "a" "a" "c" "t" "a" "g" "a" "a" "g" "c" "t"
    [4249] "g" "t" "g" "t" "t" "a" "g" "a" "a" "c" "a" "g" "c" "a" "t" "g" "g" "g"
    [4267] "a" "g" "c" "c" "a" "g" "c" "c" "t" "t" "c" "t" "a" "a" "c" "a" "g" "c"
    [4285] "t" "a" "c" "c" "c" "t" "t" "c" "c" "a" "t" "c" "a" "t" "a" "a" "g" "t"
    [4303] "g" "a" "c" "t" "c" "t" "t" "c" "t" "g" "c" "c" "c" "t" "t" "g" "a" "g"
    [4321] "g" "a" "c" "c" "t" "g" "c" "g" "a" "a" "a" "t" "c" "c" "a" "g" "a" "a"
    [4339] "c" "a" "a" "a" "g" "c" "a" "c" "a" "t" "c" "a" "g" "a" "a" "a" "a" "a"
    [4357] "g" "c" "a" "g" "t" "a" "t" "t" "a" "a" "c" "t" "t" "c" "a" "c" "a" "g"
    [4375] "a" "a" "a" "a" "g" "t" "a" "g" "t" "g" "a" "a" "t" "a" "c" "c" "c" "t"
    [4393] "a" "t" "a" "a" "g" "c" "c" "a" "g" "a" "a" "t" "c" "c" "a" "g" "a" "a"
    [4411] "g" "g" "c" "c" "t" "t" "t" "c" "t" "g" "c" "t" "g" "a" "c" "a" "a" "g"
    [4429] "t" "t" "t" "g" "a" "g" "g" "t" "g" "t" "c" "t" "g" "c" "a" "g" "a" "t"
    [4447] "a" "g" "t" "t" "c" "t" "a" "c" "c" "a" "g" "t" "a" "a" "a" "a" "a" "t"
    [4465] "a" "a" "a" "g" "a" "a" "c" "c" "a" "g" "g" "a" "g" "t" "g" "g" "a" "a"
    [4483] "a" "g" "g" "t" "c" "a" "t" "c" "c" "c" "c" "t" "t" "c" "t" "a" "a" "a"
    [4501] "t" "g" "c" "c" "c" "a" "t" "c" "a" "t" "t" "a" "g" "a" "t" "g" "a" "t"
    [4519] "a" "g" "g" "t" "g" "g" "t" "a" "c" "a" "t" "g" "c" "a" "c" "a" "g" "t"
    [4537] "t" "g" "c" "t" "c" "t" "g" "g" "g" "a" "g" "t" "c" "t" "t" "c" "a" "g"
    [4555] "a" "a" "t" "a" "g" "a" "a" "a" "c" "t" "a" "c" "c" "c" "a" "t" "c" "t"
    [4573] "c" "a" "a" "g" "a" "g" "g" "a" "g" "c" "t" "c" "a" "t" "t" "a" "a" "g"
    [4591] "g" "t" "t" "g" "t" "t" "g" "a" "t" "g" "t" "g" "g" "a" "g" "g" "a" "g"
    [4609] "c" "a" "a" "c" "a" "g" "c" "t" "g" "g" "a" "a" "g" "a" "g" "t" "c" "t"
    [4627] "g" "g" "g" "c" "c" "a" "c" "a" "c" "g" "a" "t" "t" "t" "g" "a" "c" "g"
    [4645] "g" "a" "a" "a" "c" "a" "t" "c" "t" "t" "a" "c" "t" "t" "g" "c" "c" "a"
    [4663] "a" "g" "g" "c" "a" "a" "g" "a" "t" "c" "t" "a" "g" "a" "g" "g" "g" "a"
    [4681] "a" "c" "c" "c" "c" "t" "t" "a" "c" "c" "t" "g" "g" "a" "a" "t" "c" "t"
    [4699] "g" "g" "a" "a" "t" "c" "a" "g" "c" "c" "t" "c" "t" "t" "c" "t" "c" "t"
    [4717] "g" "a" "t" "g" "a" "c" "c" "c" "t" "g" "a" "a" "t" "c" "t" "g" "a" "t"
    [4735] "c" "c" "t" "t" "c" "t" "g" "a" "a" "g" "a" "c" "a" "g" "a" "g" "c" "c"
    [4753] "c" "c" "a" "g" "a" "g" "t" "c" "a" "g" "c" "t" "c" "g" "t" "g" "t" "t"
    [4771] "g" "g" "c" "a" "a" "c" "a" "t" "a" "c" "c" "a" "t" "c" "t" "t" "c" "a"
    [4789] "a" "c" "c" "t" "c" "t" "g" "c" "a" "t" "t" "g" "a" "a" "a" "g" "t" "t"
    [4807] "c" "c" "c" "c" "a" "a" "t" "t" "g" "a" "a" "a" "g" "t" "t" "g" "c" "a"
    [4825] "g" "a" "a" "t" "c" "t" "g" "c" "c" "c" "a" "g" "a" "g" "t" "c" "c" "a"
    [4843] "g" "c" "t" "g" "c" "t" "g" "c" "t" "c" "a" "t" "a" "c" "t" "a" "c" "t"
    [4861] "g" "a" "t" "a" "c" "t" "g" "c" "t" "g" "g" "g" "t" "a" "t" "a" "a" "t"
    [4879] "g" "c" "a" "a" "t" "g" "g" "a" "a" "g" "a" "a" "a" "g" "t" "g" "t" "g"
    [4897] "a" "g" "c" "a" "g" "g" "g" "a" "g" "a" "a" "g" "c" "c" "a" "g" "a" "a"
    [4915] "t" "t" "g" "a" "c" "a" "g" "c" "t" "t" "c" "a" "a" "c" "a" "g" "a" "a"
    [4933] "a" "g" "g" "g" "t" "c" "a" "a" "c" "a" "a" "a" "a" "g" "a" "a" "t" "g"
    [4951] "t" "c" "c" "a" "t" "g" "g" "t" "g" "g" "t" "g" "t" "c" "t" "g" "g" "c"
    [4969] "c" "t" "g" "a" "c" "c" "c" "c" "a" "g" "a" "a" "g" "a" "a" "t" "t" "t"
    [4987] "a" "t" "g" "c" "t" "c" "g" "t" "g" "t" "a" "c" "a" "a" "g" "t" "t" "t"
    [5005] "g" "c" "c" "a" "g" "a" "a" "a" "a" "c" "a" "c" "c" "a" "c" "a" "t" "c"
    [5023] "a" "c" "t" "t" "t" "a" "a" "c" "t" "a" "a" "t" "c" "t" "a" "a" "t" "t"
    [5041] "a" "c" "t" "g" "a" "a" "g" "a" "g" "a" "c" "t" "a" "c" "t" "c" "a" "t"
    [5059] "g" "t" "t" "g" "t" "t" "a" "t" "g" "a" "a" "a" "a" "c" "a" "g" "a" "t"
    [5077] "g" "c" "t" "g" "a" "g" "t" "t" "t" "g" "t" "g" "t" "g" "t" "g" "a" "a"
    [5095] "c" "g" "g" "a" "c" "a" "c" "t" "g" "a" "a" "a" "t" "a" "t" "t" "t" "t"
    [5113] "c" "t" "a" "g" "g" "a" "a" "t" "t" "g" "c" "g" "g" "g" "a" "g" "g" "a"
    [5131] "a" "a" "a" "t" "g" "g" "g" "t" "a" "g" "t" "t" "a" "g" "c" "t" "a" "t"
    [5149] "t" "t" "c" "t" "g" "g" "g" "t" "g" "a" "c" "c" "c" "a" "g" "t" "c" "t"
    [5167] "a" "t" "t" "a" "a" "a" "g" "a" "a" "a" "g" "a" "a" "a" "a" "a" "t" "g"
    [5185] "c" "t" "g" "a" "a" "t" "g" "a" "g" "c" "a" "t" "g" "a" "t" "t" "t" "t"
    [5203] "g" "a" "a" "g" "t" "c" "a" "g" "a" "g" "g" "a" "g" "a" "t" "g" "t" "g"
    [5221] "g" "t" "c" "a" "a" "t" "g" "g" "a" "a" "g" "a" "a" "a" "c" "c" "a" "c"
    [5239] "c" "a" "a" "g" "g" "t" "c" "c" "a" "a" "a" "g" "c" "g" "a" "g" "c" "a"
    [5257] "a" "g" "a" "g" "a" "a" "t" "c" "c" "c" "a" "g" "g" "a" "c" "a" "g" "a"
    [5275] "a" "a" "g" "a" "t" "c" "t" "t" "c" "a" "g" "g" "g" "g" "g" "c" "t" "a"
    [5293] "g" "a" "a" "a" "t" "c" "t" "g" "t" "t" "g" "c" "t" "a" "t" "g" "g" "g"
    [5311] "c" "c" "c" "t" "t" "c" "a" "c" "c" "a" "a" "c" "a" "t" "g" "c" "c" "c"
    [5329] "a" "c" "a" "g" "a" "t" "c" "a" "a" "c" "t" "g" "g" "a" "a" "t" "g" "g"
    [5347] "a" "t" "g" "g" "t" "a" "c" "a" "g" "c" "t" "g" "t" "g" "t" "g" "g" "t"
    [5365] "g" "c" "t" "t" "c" "t" "g" "t" "g" "g" "t" "g" "a" "a" "g" "g" "a" "g"
    [5383] "c" "t" "t" "t" "c" "a" "t" "c" "a" "t" "t" "c" "a" "c" "c" "c" "t" "t"
    [5401] "g" "g" "c" "a" "c" "a" "g" "g" "t" "g" "t" "c" "c" "a" "c" "c" "c" "a"
    [5419] "a" "t" "t" "g" "t" "g" "g" "t" "t" "g" "t" "g" "c" "a" "g" "c" "c" "a"
    [5437] "g" "a" "t" "g" "c" "c" "t" "g" "g" "a" "c" "a" "g" "a" "g" "g" "a" "c"
    [5455] "a" "a" "t" "g" "g" "c" "t" "t" "c" "c" "a" "t" "g" "c" "a" "a" "t" "t"
    [5473] "g" "g" "g" "c" "a" "g" "a" "t" "g" "t" "g" "t" "g" "a" "g" "g" "c" "a"
    [5491] "c" "c" "t" "g" "t" "g" "g" "t" "g" "a" "c" "c" "c" "g" "a" "g" "a" "g"
    [5509] "t" "g" "g" "g" "t" "g" "t" "t" "g" "g" "a" "c" "a" "g" "t" "g" "t" "a"
    [5527] "g" "c" "a" "c" "t" "c" "t" "a" "c" "c" "a" "g" "t" "g" "c" "c" "a" "g"
    [5545] "g" "a" "g" "c" "t" "g" "g" "a" "c" "a" "c" "c" "t" "a" "c" "c" "t" "g"
    [5563] "a" "t" "a" "c" "c" "c" "c" "a" "g" "a" "t" "c" "c" "c" "c" "c" "a" "c"
    [5581] "a" "g" "c" "c" "a" "c" "t" "a" "c" "t" "g" "a"
    attr(,"name")
    [1] "DQ190455.BRCA1"
    attr(,"Annot")
    [1] ">DQ190455.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $DQ190456.BRCA1
       [1] "a" "t" "g" "g" "a" "t" "t" "t" "a" "t" "c" "t" "g" "c" "t" "c" "t" "t"
      [19] "c" "g" "c" "g" "t" "t" "g" "a" "a" "g" "a" "a" "g" "t" "a" "c" "a" "a"
      [37] "a" "a" "t" "g" "t" "c" "a" "t" "t" "a" "a" "t" "g" "c" "t" "a" "t" "g"
      [55] "c" "a" "g" "a" "a" "a" "a" "t" "c" "t" "t" "a" "g" "a" "g" "t" "g" "t"
      [73] "c" "c" "c" "a" "t" "c" "t" "g" "t" "c" "t" "g" "g" "a" "g" "t" "t" "g"
      [91] "a" "t" "c" "a" "a" "g" "g" "a" "a" "c" "c" "t" "g" "t" "c" "t" "c" "c"
     [109] "a" "c" "a" "a" "a" "g" "t" "g" "t" "g" "a" "c" "c" "a" "c" "a" "t" "a"
     [127] "t" "t" "t" "t" "g" "c" "a" "a" "a" "t" "t" "t" "t" "g" "c" "a" "t" "g"
     [145] "c" "t" "g" "a" "a" "a" "c" "t" "t" "c" "t" "c" "a" "a" "c" "c" "a" "g"
     [163] "a" "a" "g" "a" "a" "a" "g" "g" "g" "c" "c" "t" "t" "c" "a" "c" "a" "g"
     [181] "t" "g" "t" "c" "c" "t" "t" "t" "a" "t" "g" "t" "a" "a" "g" "a" "a" "t"
     [199] "g" "a" "t" "a" "t" "a" "a" "c" "c" "a" "a" "a" "a" "g" "g" "a" "g" "c"
     [217] "c" "t" "a" "c" "a" "a" "g" "a" "a" "a" "g" "t" "a" "c" "g" "a" "g" "a"
     [235] "t" "t" "t" "a" "g" "t" "c" "a" "a" "c" "t" "t" "g" "t" "t" "g" "a" "a"
     [253] "g" "a" "g" "c" "t" "a" "t" "t" "g" "a" "a" "a" "a" "t" "c" "a" "t" "t"
     [271] "t" "g" "t" "g" "c" "t" "t" "t" "t" "c" "a" "g" "c" "t" "t" "g" "a" "c"
     [289] "a" "c" "a" "g" "g" "t" "t" "t" "g" "g" "a" "g" "t" "a" "t" "g" "c" "a"
     [307] "a" "a" "c" "a" "g" "c" "t" "a" "t" "a" "a" "t" "t" "t" "t" "g" "c" "a"
     [325] "a" "a" "a" "a" "a" "g" "g" "a" "a" "a" "a" "t" "a" "a" "c" "t" "c" "t"
     [343] "c" "c" "t" "g" "a" "a" "c" "a" "t" "c" "t" "a" "a" "a" "a" "g" "a" "t"
     [361] "g" "a" "a" "g" "t" "t" "t" "c" "t" "a" "t" "c" "a" "t" "c" "c" "a" "a"
     [379] "a" "g" "t" "a" "t" "g" "g" "g" "c" "t" "a" "c" "a" "g" "a" "a" "a" "c"
     [397] "c" "g" "t" "g" "c" "c" "a" "a" "a" "a" "g" "a" "c" "t" "t" "c" "t" "a"
     [415] "c" "a" "g" "a" "g" "t" "g" "a" "a" "c" "c" "c" "g" "a" "a" "a" "a" "t"
     [433] "c" "c" "t" "t" "c" "c" "t" "t" "g" "c" "a" "g" "g" "a" "a" "a" "c" "c"
     [451] "a" "g" "t" "c" "t" "c" "a" "g" "t" "g" "t" "c" "c" "a" "a" "c" "t" "c"
     [469] "t" "c" "t" "a" "a" "c" "c" "t" "t" "g" "g" "a" "a" "c" "t" "g" "t" "g"
     [487] "a" "g" "a" "a" "c" "t" "c" "t" "g" "a" "g" "g" "a" "c" "a" "a" "a" "g"
     [505] "c" "a" "g" "c" "g" "g" "a" "t" "a" "c" "a" "a" "c" "c" "t" "c" "a" "a"
     [523] "a" "a" "g" "a" "c" "g" "t" "c" "t" "g" "t" "c" "t" "a" "c" "a" "t" "t"
     [541] "g" "a" "a" "t" "t" "g" "g" "g" "a" "t" "c" "t" "g" "a" "t" "t" "c" "t"
     [559] "t" "c" "t" "g" "a" "a" "g" "a" "t" "a" "c" "c" "g" "t" "t" "a" "a" "t"
     [577] "a" "a" "g" "g" "c" "a" "a" "c" "t" "t" "a" "t" "t" "g" "c" "a" "g" "t"
     [595] "g" "t" "g" "g" "g" "a" "g" "a" "t" "c" "a" "a" "g" "a" "a" "t" "t" "g"
     [613] "t" "t" "a" "c" "a" "a" "a" "t" "c" "a" "c" "c" "c" "c" "t" "c" "a" "a"
     [631] "g" "g" "a" "a" "c" "c" "a" "g" "g" "g" "a" "t" "g" "a" "a" "a" "t" "c"
     [649] "a" "g" "t" "t" "t" "g" "g" "a" "t" "t" "c" "t" "g" "c" "a" "a" "a" "a"
     [667] "a" "a" "g" "g" "c" "t" "g" "c" "t" "t" "g" "t" "g" "a" "a" "t" "t" "t"
     [685] "t" "c" "t" "g" "a" "g" "a" "c" "g" "g" "a" "t" "g" "t" "a" "a" "c" "a"
     [703] "a" "a" "t" "a" "c" "t" "g" "a" "a" "c" "a" "t" "c" "a" "t" "c" "a" "a"
     [721] "c" "c" "c" "a" "g" "t" "a" "a" "t" "a" "a" "t" "g" "a" "t" "t" "t" "g"
     [739] "a" "a" "c" "a" "c" "c" "a" "c" "t" "g" "a" "g" "a" "a" "g" "c" "g" "t"
     [757] "g" "c" "a" "g" "c" "t" "g" "a" "g" "a" "g" "g" "c" "a" "t" "c" "c" "a"
     [775] "g" "a" "a" "a" "a" "g" "t" "a" "t" "c" "a" "g" "g" "g" "t" "a" "g" "t"
     [793] "t" "c" "t" "g" "t" "t" "t" "c" "a" "a" "a" "c" "t" "t" "g" "c" "a" "t"
     [811] "g" "t" "g" "g" "a" "g" "c" "c" "a" "t" "g" "t" "g" "g" "c" "a" "c" "a"
     [829] "a" "a" "t" "a" "c" "t" "c" "a" "t" "g" "c" "c" "a" "g" "c" "t" "c" "a"
     [847] "t" "t" "a" "c" "a" "g" "c" "a" "t" "g" "a" "g" "a" "a" "c" "a" "g" "c"
     [865] "a" "g" "t" "t" "t" "a" "t" "t" "a" "c" "t" "c" "a" "c" "t" "a" "a" "a"
     [883] "g" "a" "c" "a" "g" "a" "a" "t" "g" "a" "a" "t" "g" "t" "a" "g" "a" "a"
     [901] "a" "a" "g" "g" "c" "t" "g" "a" "a" "t" "t" "c" "t" "g" "t" "a" "a" "t"
     [919] "a" "a" "a" "a" "g" "c" "a" "a" "a" "c" "a" "g" "c" "c" "t" "g" "g" "c"
     [937] "t" "t" "a" "g" "c" "a" "a" "g" "g" "a" "g" "c" "c" "a" "a" "c" "a" "t"
     [955] "a" "a" "c" "a" "g" "a" "t" "g" "g" "g" "c" "t" "g" "g" "a" "a" "g" "t"
     [973] "a" "a" "g" "g" "a" "a" "a" "c" "a" "t" "g" "t" "a" "a" "t" "g" "a" "t"
     [991] "a" "g" "g" "c" "g" "g" "a" "c" "t" "c" "c" "c" "a" "g" "c" "a" "c" "a"
    [1009] "g" "a" "a" "a" "a" "a" "a" "a" "g" "g" "t" "a" "g" "a" "t" "c" "t" "g"
    [1027] "a" "a" "t" "g" "c" "t" "g" "a" "t" "c" "c" "c" "c" "t" "g" "t" "g" "t"
    [1045] "g" "a" "g" "a" "g" "a" "a" "a" "a" "g" "a" "a" "t" "g" "g" "a" "a" "t"
    [1063] "a" "a" "g" "c" "a" "g" "a" "a" "a" "c" "t" "g" "c" "c" "a" "t" "g" "c"
    [1081] "t" "c" "a" "g" "a" "g" "a" "a" "t" "c" "c" "t" "a" "g" "a" "g" "a" "t"
    [1099] "a" "c" "t" "g" "a" "a" "g" "a" "t" "g" "t" "t" "c" "c" "t" "t" "g" "g"
    [1117] "a" "t" "a" "a" "c" "a" "c" "t" "a" "a" "a" "t" "a" "g" "c" "a" "g" "c"
    [1135] "a" "t" "t" "c" "a" "g" "a" "a" "a" "g" "t" "t" "a" "a" "t" "g" "a" "g"
    [1153] "t" "g" "g" "t" "t" "t" "t" "c" "c" "a" "g" "a" "a" "g" "t" "g" "a" "t"
    [1171] "g" "a" "a" "c" "t" "g" "t" "t" "a" "g" "g" "t" "t" "c" "t" "g" "a" "t"
    [1189] "g" "a" "c" "t" "c" "a" "c" "a" "t" "g" "a" "t" "g" "g" "g" "g" "a" "g"
    [1207] "t" "c" "t" "g" "a" "a" "t" "c" "a" "a" "a" "t" "g" "c" "c" "a" "a" "a"
    [1225] "g" "t" "a" "g" "c" "t" "g" "a" "t" "g" "t" "a" "t" "t" "g" "g" "a" "c"
    [1243] "g" "t" "t" "c" "t" "a" "a" "a" "t" "g" "a" "g" "g" "t" "a" "g" "a" "t"
    [1261] "g" "a" "a" "t" "a" "t" "t" "c" "t" "g" "g" "t" "t" "c" "t" "t" "c" "a"
    [1279] "g" "a" "g" "a" "a" "a" "a" "t" "a" "g" "a" "c" "t" "t" "a" "c" "t" "g"
    [1297] "g" "c" "c" "a" "g" "t" "g" "a" "t" "c" "c" "t" "c" "a" "t" "g" "a" "g"
    [1315] "g" "c" "t" "t" "t" "a" "a" "t" "a" "t" "g" "t" "a" "a" "a" "a" "g" "t"
    [1333] "g" "a" "a" "a" "g" "a" "g" "t" "t" "c" "a" "c" "t" "c" "c" "a" "a" "a"
    [1351] "t" "c" "a" "g" "t" "a" "g" "a" "g" "a" "g" "t" "a" "a" "t" "a" "t" "t"
    [1369] "g" "a" "a" "g" "a" "c" "a" "a" "a" "a" "t" "a" "t" "t" "t" "g" "g" "g"
    [1387] "a" "a" "a" "a" "c" "c" "t" "a" "t" "c" "g" "g" "a" "a" "g" "a" "a" "g"
    [1405] "g" "c" "a" "a" "g" "c" "c" "t" "c" "c" "c" "c" "a" "a" "c" "t" "t" "a"
    [1423] "a" "g" "c" "c" "a" "t" "g" "t" "a" "a" "c" "t" "g" "a" "a" "a" "a" "t"
    [1441] "c" "t" "a" "a" "t" "t" "a" "t" "a" "g" "g" "a" "g" "c" "a" "t" "t" "t"
    [1459] "g" "t" "t" "a" "c" "t" "g" "a" "g" "c" "c" "a" "c" "a" "g" "a" "t" "a"
    [1477] "a" "t" "a" "c" "a" "a" "g" "a" "g" "c" "g" "t" "c" "c" "c" "c" "t" "c"
    [1495] "a" "c" "a" "a" "a" "t" "a" "a" "a" "t" "t" "a" "a" "a" "g" "c" "g" "t"
    [1513] "a" "a" "a" "a" "g" "g" "a" "g" "a" "c" "c" "t" "a" "c" "a" "t" "c" "a"
    [1531] "g" "g" "c" "c" "t" "t" "c" "a" "t" "c" "c" "t" "g" "a" "g" "g" "a" "t"
    [1549] "t" "t" "t" "a" "t" "c" "a" "a" "g" "a" "a" "a" "g" "c" "a" "g" "a" "t"
    [1567] "t" "t" "g" "g" "c" "a" "g" "t" "t" "c" "a" "a" "a" "a" "g" "a" "c" "t"
    [1585] "c" "c" "t" "g" "a" "a" "a" "t" "g" "a" "t" "a" "a" "a" "t" "c" "a" "g"
    [1603] "g" "g" "a" "a" "c" "t" "a" "a" "c" "c" "a" "a" "a" "c" "g" "g" "a" "g"
    [1621] "c" "a" "g" "a" "a" "t" "g" "g" "t" "c" "a" "a" "g" "t" "g" "a" "t" "g"
    [1639] "a" "a" "t" "a" "t" "t" "a" "c" "t" "a" "a" "t" "a" "g" "t" "g" "g" "t"
    [1657] "c" "a" "t" "g" "a" "g" "a" "a" "t" "a" "a" "a" "a" "c" "a" "a" "a" "a"
    [1675] "g" "g" "t" "g" "a" "t" "t" "c" "t" "a" "t" "t" "c" "a" "g" "a" "a" "t"
    [1693] "g" "a" "g" "a" "a" "a" "a" "a" "t" "c" "c" "t" "a" "a" "c" "c" "c" "a"
    [1711] "a" "t" "a" "g" "a" "a" "t" "c" "a" "c" "t" "c" "g" "a" "a" "a" "a" "a"
    [1729] "g" "a" "a" "t" "c" "t" "g" "c" "t" "t" "t" "c" "a" "a" "a" "a" "c" "g"
    [1747] "a" "a" "a" "g" "c" "t" "g" "a" "a" "c" "c" "t" "a" "t" "a" "a" "g" "c"
    [1765] "a" "g" "c" "a" "g" "t" "a" "t" "a" "a" "g" "c" "a" "a" "t" "a" "t" "g"
    [1783] "g" "a" "a" "c" "t" "c" "g" "a" "a" "t" "t" "a" "a" "a" "t" "a" "t" "c"
    [1801] "c" "a" "c" "a" "a" "t" "t" "c" "a" "a" "a" "a" "g" "c" "a" "c" "c" "t"
    [1819] "a" "a" "a" "a" "a" "g" "a" "a" "t" "a" "g" "g" "c" "t" "g" "a" "g" "g"
    [1837] "a" "g" "g" "a" "a" "g" "t" "c" "t" "t" "c" "t" "a" "c" "c" "a" "g" "g"
    [1855] "c" "a" "t" "a" "t" "t" "c" "a" "t" "g" "c" "g" "c" "t" "t" "g" "a" "a"
    [1873] "c" "t" "a" "g" "t" "a" "g" "t" "c" "a" "g" "t" "a" "g" "a" "a" "a" "t"
    [1891] "c" "t" "a" "a" "g" "c" "c" "c" "a" "c" "c" "t" "a" "a" "t" "t" "g" "t"
    [1909] "a" "c" "t" "g" "a" "a" "t" "t" "g" "c" "a" "a" "a" "t" "t" "g" "a" "t"
    [1927] "a" "g" "t" "t" "g" "t" "t" "c" "t" "a" "g" "c" "a" "g" "t" "g" "a" "a"
    [1945] "g" "a" "g" "a" "t" "a" "a" "a" "g" "a" "a" "a" "a" "a" "a" "a" "a" "g"
    [1963] "t" "a" "c" "a" "a" "c" "c" "a" "a" "a" "t" "g" "c" "c" "a" "g" "t" "c"
    [1981] "a" "g" "g" "c" "a" "c" "a" "g" "c" "a" "g" "a" "a" "a" "c" "c" "t" "a"
    [1999] "c" "a" "a" "c" "t" "c" "a" "t" "g" "g" "a" "a" "g" "g" "t" "a" "a" "a"
    [2017] "g" "a" "a" "c" "c" "t" "g" "c" "a" "a" "c" "t" "g" "g" "a" "g" "c" "c"
    [2035] "a" "a" "g" "a" "a" "g" "a" "g" "t" "a" "a" "c" "a" "a" "g" "c" "c" "a"
    [2053] "a" "a" "t" "g" "a" "a" "c" "a" "g" "a" "c" "a" "a" "g" "t" "a" "a" "a"
    [2071] "a" "g" "a" "c" "a" "t" "g" "a" "c" "a" "g" "c" "g" "a" "t" "a" "c" "t"
    [2089] "t" "t" "c" "c" "c" "a" "g" "a" "g" "c" "t" "g" "a" "a" "g" "t" "t" "a"
    [2107] "a" "c" "a" "a" "a" "t" "g" "c" "a" "c" "c" "t" "g" "g" "t" "t" "c" "t"
    [2125] "t" "t" "t" "a" "c" "t" "a" "a" "g" "t" "g" "t" "t" "c" "a" "a" "a" "t"
    [2143] "a" "c" "c" "a" "g" "t" "g" "a" "a" "c" "t" "t" "a" "a" "a" "g" "a" "a"
    [2161] "t" "t" "t" "g" "t" "c" "a" "a" "t" "c" "c" "t" "a" "g" "c" "c" "t" "t"
    [2179] "c" "c" "a" "a" "g" "a" "g" "a" "a" "g" "a" "a" "a" "a" "a" "g" "a" "a"
    [2197] "g" "a" "g" "a" "a" "a" "c" "t" "a" "g" "a" "a" "a" "c" "a" "g" "t" "t"
    [2215] "a" "a" "a" "g" "t" "g" "t" "c" "t" "a" "a" "t" "a" "a" "t" "g" "c" "t"
    [2233] "g" "a" "a" "g" "a" "c" "c" "c" "c" "a" "a" "a" "g" "a" "t" "c" "t" "c"
    [2251] "a" "t" "g" "t" "t" "a" "a" "g" "t" "g" "g" "a" "g" "a" "a" "a" "g" "g"
    [2269] "g" "t" "t" "t" "t" "g" "c" "a" "a" "a" "c" "t" "g" "a" "a" "a" "g" "a"
    [2287] "t" "c" "t" "g" "t" "a" "g" "a" "g" "a" "g" "t" "a" "g" "c" "a" "g" "t"
    [2305] "a" "t" "t" "t" "c" "a" "t" "t" "g" "g" "t" "a" "c" "c" "t" "g" "g" "t"
    [2323] "a" "c" "t" "g" "a" "t" "t" "a" "t" "g" "g" "c" "a" "c" "t" "c" "a" "g"
    [2341] "g" "a" "a" "a" "g" "t" "a" "t" "c" "t" "c" "g" "t" "t" "a" "c" "t" "g"
    [2359] "g" "a" "a" "g" "t" "t" "a" "g" "c" "a" "c" "t" "c" "t" "a" "g" "g" "g"
    [2377] "a" "a" "g" "g" "c" "a" "a" "a" "a" "a" "c" "a" "g" "a" "a" "c" "c" "a"
    [2395] "a" "a" "t" "a" "a" "a" "t" "g" "t" "g" "t" "g" "a" "g" "t" "c" "a" "g"
    [2413] "t" "g" "t" "g" "c" "a" "g" "c" "a" "t" "t" "t" "g" "a" "a" "a" "a" "c"
    [2431] "c" "c" "c" "a" "a" "g" "g" "g" "a" "c" "t" "a" "a" "t" "t" "c" "a" "t"
    [2449] "g" "g" "t" "t" "g" "t" "t" "c" "c" "a" "a" "a" "g" "a" "t" "a" "a" "t"
    [2467] "a" "g" "a" "a" "a" "t" "g" "a" "c" "a" "c" "a" "g" "a" "a" "g" "g" "c"
    [2485] "t" "t" "t" "a" "a" "g" "t" "a" "t" "c" "c" "a" "t" "t" "g" "g" "g" "a"
    [2503] "c" "a" "t" "g" "a" "a" "g" "t" "t" "a" "a" "c" "c" "a" "c" "a" "g" "t"
    [2521] "c" "g" "g" "g" "a" "a" "a" "c" "a" "a" "g" "c" "a" "t" "a" "g" "a" "a"
    [2539] "a" "t" "g" "g" "a" "a" "g" "a" "a" "a" "g" "t" "g" "a" "a" "c" "t" "t"
    [2557] "g" "a" "t" "g" "c" "t" "c" "a" "g" "t" "a" "t" "t" "t" "g" "c" "a" "g"
    [2575] "a" "a" "t" "a" "c" "a" "t" "t" "c" "a" "a" "g" "g" "t" "t" "t" "c" "a"
    [2593] "a" "a" "g" "c" "g" "c" "c" "a" "g" "t" "c" "a" "t" "t" "t" "g" "c" "t"
    [2611] "c" "c" "g" "t" "t" "t" "t" "c" "a" "a" "a" "t" "c" "c" "a" "g" "g" "a"
    [2629] "a" "a" "t" "g" "c" "a" "g" "a" "a" "g" "a" "g" "g" "a" "a" "t" "g" "t"
    [2647] "g" "c" "a" "a" "c" "a" "t" "t" "c" "t" "c" "t" "g" "c" "c" "c" "a" "c"
    [2665] "t" "c" "t" "g" "g" "g" "t" "c" "c" "t" "t" "a" "a" "a" "g" "a" "a" "a"
    [2683] "c" "a" "a" "a" "g" "t" "c" "c" "a" "a" "a" "a" "g" "t" "c" "a" "c" "t"
    [2701] "t" "t" "t" "g" "a" "a" "t" "g" "t" "g" "a" "a" "c" "a" "a" "a" "a" "g"
    [2719] "g" "a" "a" "g" "a" "a" "a" "a" "t" "c" "a" "a" "g" "g" "a" "a" "a" "g"
    [2737] "a" "a" "t" "g" "a" "g" "t" "c" "t" "a" "a" "t" "a" "t" "c" "a" "a" "g"
    [2755] "c" "c" "t" "g" "t" "a" "c" "a" "g" "a" "c" "a" "g" "t" "t" "a" "a" "t"
    [2773] "a" "t" "c" "a" "c" "t" "g" "c" "a" "g" "g" "c" "t" "t" "t" "c" "c" "t"
    [2791] "g" "t" "g" "g" "t" "t" "g" "g" "t" "c" "a" "g" "a" "a" "a" "g" "a" "t"
    [2809] "a" "a" "g" "c" "c" "a" "g" "t" "t" "g" "a" "t" "a" "a" "t" "g" "c" "c"
    [2827] "a" "a" "a" "t" "g" "t" "a" "g" "t" "a" "t" "c" "a" "a" "a" "g" "g" "a"
    [2845] "g" "g" "c" "t" "c" "t" "a" "g" "g" "t" "t" "t" "t" "g" "t" "c" "t" "a"
    [2863] "t" "c" "a" "t" "c" "t" "c" "a" "g" "t" "t" "c" "a" "g" "a" "g" "g" "c"
    [2881] "a" "a" "c" "g" "a" "a" "a" "c" "t" "g" "g" "a" "c" "t" "c" "a" "t" "t"
    [2899] "a" "c" "t" "c" "c" "a" "a" "a" "t" "a" "a" "a" "c" "a" "t" "g" "g" "a"
    [2917] "c" "t" "t" "t" "t" "a" "c" "a" "a" "a" "a" "c" "c" "c" "a" "t" "a" "t"
    [2935] "c" "g" "t" "a" "t" "a" "c" "c" "a" "c" "c" "a" "c" "t" "t" "t" "t" "t"
    [2953] "c" "c" "c" "a" "t" "c" "a" "a" "g" "t" "c" "a" "t" "t" "t" "g" "t" "t"
    [2971] "a" "a" "a" "a" "c" "t" "a" "a" "a" "t" "g" "t" "a" "a" "g" "a" "a" "a"
    [2989] "a" "a" "t" "c" "t" "g" "c" "t" "a" "g" "a" "g" "g" "a" "a" "a" "a" "c"
    [3007] "t" "t" "t" "g" "a" "g" "g" "a" "a" "c" "a" "t" "t" "c" "a" "a" "t" "g"
    [3025] "t" "c" "a" "c" "c" "t" "g" "a" "a" "a" "g" "a" "g" "a" "a" "a" "t" "g"
    [3043] "g" "g" "a" "a" "a" "t" "g" "a" "g" "a" "a" "c" "a" "t" "t" "c" "c" "a"
    [3061] "a" "g" "t" "a" "c" "a" "g" "t" "g" "a" "g" "c" "a" "c" "a" "a" "t" "t"
    [3079] "a" "g" "c" "c" "g" "t" "a" "a" "t" "a" "a" "c" "a" "t" "t" "a" "g" "a"
    [3097] "g" "a" "a" "a" "a" "t" "g" "t" "t" "t" "t" "t" "a" "a" "a" "g" "a" "a"
    [3115] "g" "c" "c" "a" "g" "c" "t" "c" "a" "a" "g" "c" "a" "a" "t" "a" "t" "t"
    [3133] "a" "a" "t" "g" "a" "a" "g" "t" "a" "g" "g" "t" "t" "c" "c" "a" "g" "t"
    [3151] "a" "c" "t" "a" "a" "t" "g" "a" "a" "g" "t" "g" "g" "g" "c" "t" "c" "c"
    [3169] "a" "g" "t" "a" "t" "t" "a" "a" "t" "g" "a" "a" "a" "t" "a" "g" "g" "t"
    [3187] "t" "c" "c" "a" "g" "t" "g" "a" "t" "g" "a" "a" "a" "a" "c" "a" "t" "t"
    [3205] "c" "a" "a" "g" "c" "a" "g" "a" "a" "c" "t" "a" "g" "g" "t" "a" "g" "a"
    [3223] "a" "a" "c" "a" "g" "a" "g" "g" "g" "c" "c" "a" "a" "a" "a" "t" "t" "g"
    [3241] "a" "a" "t" "g" "c" "t" "a" "t" "g" "c" "t" "t" "a" "g" "a" "t" "t" "a"
    [3259] "g" "g" "g" "g" "t" "t" "t" "t" "g" "c" "a" "a" "c" "c" "t" "g" "a" "g"
    [3277] "g" "t" "c" "t" "a" "t" "a" "a" "a" "c" "a" "a" "a" "g" "t" "c" "t" "t"
    [3295] "c" "c" "t" "g" "g" "a" "a" "g" "t" "a" "a" "t" "t" "g" "t" "a" "a" "g"
    [3313] "c" "a" "t" "c" "c" "t" "g" "a" "a" "a" "t" "a" "a" "a" "a" "a" "a" "g"
    [3331] "c" "a" "a" "g" "a" "a" "t" "a" "t" "g" "a" "a" "g" "a" "a" "g" "t" "a"
    [3349] "g" "t" "t" "c" "a" "g" "a" "c" "t" "g" "t" "t" "a" "a" "t" "a" "c" "a"
    [3367] "g" "a" "t" "t" "t" "c" "t" "c" "t" "c" "c" "a" "t" "a" "t" "c" "t" "g"
    [3385] "a" "t" "t" "t" "c" "a" "g" "a" "t" "a" "a" "c" "t" "t" "a" "g" "a" "a"
    [3403] "c" "a" "g" "c" "c" "t" "a" "t" "g" "g" "g" "a" "a" "g" "t" "a" "g" "t"
    [3421] "c" "a" "t" "g" "c" "a" "t" "c" "t" "c" "a" "g" "g" "t" "t" "t" "g" "t"
    [3439] "t" "c" "t" "g" "a" "g" "a" "c" "a" "c" "c" "t" "g" "a" "t" "g" "a" "c"
    [3457] "c" "t" "g" "t" "t" "a" "g" "a" "t" "g" "a" "t" "g" "g" "t" "g" "a" "a"
    [3475] "a" "t" "a" "a" "a" "g" "g" "a" "a" "g" "a" "t" "a" "c" "t" "a" "g" "t"
    [3493] "t" "t" "t" "g" "c" "t" "g" "a" "a" "a" "a" "t" "g" "a" "c" "a" "t" "t"
    [3511] "a" "a" "g" "g" "a" "a" "a" "g" "t" "t" "c" "t" "g" "c" "t" "g" "t" "t"
    [3529] "t" "t" "t" "a" "g" "c" "a" "a" "a" "a" "g" "c" "g" "t" "c" "c" "a" "g"
    [3547] "a" "a" "a" "g" "g" "a" "g" "a" "g" "c" "t" "t" "a" "g" "c" "a" "g" "g"
    [3565] "a" "g" "t" "c" "c" "t" "a" "g" "c" "c" "c" "t" "t" "t" "c" "a" "c" "c"
    [3583] "c" "a" "t" "a" "c" "a" "c" "a" "t" "t" "t" "g" "g" "c" "t" "c" "a" "g"
    [3601] "g" "g" "t" "t" "a" "c" "c" "g" "a" "a" "g" "a" "g" "g" "g" "g" "c" "c"
    [3619] "a" "a" "g" "a" "a" "a" "t" "t" "a" "g" "a" "g" "t" "c" "c" "t" "c" "a"
    [3637] "g" "a" "a" "g" "a" "g" "a" "a" "c" "t" "t" "a" "t" "c" "t" "a" "g" "t"
    [3655] "g" "a" "g" "g" "a" "t" "g" "a" "a" "g" "a" "g" "c" "t" "t" "c" "c" "c"
    [3673] "t" "g" "c" "t" "t" "c" "c" "a" "a" "c" "a" "c" "t" "t" "g" "t" "t" "a"
    [3691] "t" "t" "t" "g" "g" "t" "a" "a" "a" "g" "t" "a" "a" "a" "c" "a" "a" "t"
    [3709] "a" "t" "a" "c" "c" "t" "t" "c" "t" "c" "a" "g" "t" "c" "t" "a" "c" "t"
    [3727] "a" "g" "g" "c" "a" "t" "a" "g" "c" "a" "c" "c" "g" "t" "t" "g" "c" "t"
    [3745] "a" "c" "c" "g" "a" "g" "t" "g" "t" "c" "t" "g" "t" "c" "t" "a" "a" "g"
    [3763] "a" "a" "c" "a" "c" "a" "g" "a" "g" "g" "a" "g" "a" "a" "t" "t" "t" "a"
    [3781] "t" "t" "a" "t" "c" "a" "t" "t" "g" "a" "a" "g" "a" "a" "t" "a" "g" "c"
    [3799] "t" "t" "a" "a" "a" "t" "g" "a" "c" "t" "g" "c" "a" "g" "t" "a" "a" "c"
    [3817] "c" "a" "g" "g" "t" "a" "a" "t" "a" "t" "t" "g" "g" "c" "a" "a" "a" "g"
    [3835] "g" "c" "a" "t" "c" "t" "c" "a" "g" "g" "a" "a" "c" "a" "t" "c" "a" "c"
    [3853] "c" "t" "t" "a" "g" "t" "g" "a" "g" "g" "a" "a" "a" "c" "a" "a" "a" "a"
    [3871] "t" "g" "t" "t" "c" "t" "g" "c" "t" "a" "g" "c" "t" "t" "g" "t" "t" "t"
    [3889] "t" "c" "t" "t" "c" "a" "c" "a" "g" "t" "g" "c" "a" "g" "t" "g" "a" "a"
    [3907] "t" "t" "g" "g" "a" "a" "g" "a" "c" "t" "t" "g" "a" "c" "t" "g" "c" "a"
    [3925] "a" "a" "t" "a" "c" "a" "a" "a" "c" "a" "c" "c" "c" "a" "g" "g" "a" "t"
    [3943] "c" "c" "t" "t" "t" "c" "t" "t" "g" "a" "t" "t" "g" "g" "t" "t" "c" "t"
    [3961] "t" "c" "c" "a" "a" "a" "c" "a" "a" "a" "t" "g" "a" "g" "g" "c" "a" "t"
    [3979] "c" "a" "g" "t" "c" "t" "g" "a" "a" "a" "g" "c" "c" "a" "g" "g" "g" "a"
    [3997] "g" "t" "t" "g" "g" "t" "c" "t" "g" "a" "g" "t" "g" "a" "c" "a" "a" "g"
    [4015] "g" "a" "a" "t" "t" "g" "g" "t" "t" "t" "c" "a" "g" "a" "t" "g" "a" "t"
    [4033] "g" "a" "a" "g" "a" "a" "a" "g" "a" "g" "g" "a" "a" "c" "g" "g" "g" "c"
    [4051] "t" "t" "g" "g" "a" "a" "g" "a" "a" "a" "a" "t" "a" "a" "t" "c" "a" "a"
    [4069] "g" "a" "a" "g" "a" "g" "c" "a" "a" "a" "g" "c" "a" "t" "g" "g" "a" "t"
    [4087] "t" "c" "a" "a" "a" "c" "t" "t" "a" "g" "g" "t" "g" "a" "a" "g" "c" "a"
    [4105] "g" "c" "a" "t" "c" "t" "g" "g" "g" "t" "g" "t" "g" "a" "g" "a" "g" "t"
    [4123] "g" "a" "a" "a" "c" "a" "a" "g" "c" "g" "t" "c" "t" "c" "t" "g" "a" "a"
    [4141] "g" "a" "c" "t" "g" "c" "t" "c" "a" "g" "g" "g" "c" "t" "a" "t" "c" "c"
    [4159] "t" "c" "t" "c" "a" "g" "a" "g" "t" "g" "a" "c" "a" "t" "t" "t" "t" "a"
    [4177] "a" "c" "c" "a" "c" "t" "c" "a" "g" "c" "a" "g" "a" "g" "g" "g" "a" "t"
    [4195] "a" "c" "c" "a" "t" "g" "c" "a" "a" "c" "a" "t" "a" "a" "c" "c" "t" "g"
    [4213] "a" "t" "a" "a" "a" "g" "c" "t" "c" "c" "a" "g" "c" "a" "g" "g" "a" "a"
    [4231] "a" "t" "g" "g" "c" "t" "g" "a" "a" "c" "t" "a" "g" "a" "a" "g" "c" "t"
    [4249] "g" "t" "g" "t" "t" "a" "g" "a" "a" "c" "a" "g" "c" "a" "t" "g" "g" "g"
    [4267] "a" "g" "c" "c" "a" "g" "c" "c" "t" "t" "c" "t" "a" "a" "c" "a" "g" "c"
    [4285] "t" "a" "c" "c" "c" "t" "t" "c" "c" "a" "t" "c" "a" "t" "a" "a" "g" "t"
    [4303] "g" "a" "c" "t" "c" "t" "t" "c" "t" "g" "c" "c" "c" "t" "t" "g" "a" "g"
    [4321] "g" "a" "c" "c" "t" "g" "c" "g" "a" "a" "a" "t" "c" "c" "a" "g" "a" "a"
    [4339] "c" "a" "a" "a" "g" "c" "a" "c" "a" "t" "c" "a" "g" "a" "a" "a" "a" "a"
    [4357] "g" "c" "a" "g" "t" "a" "t" "t" "a" "a" "c" "t" "t" "c" "a" "c" "a" "g"
    [4375] "a" "a" "a" "a" "g" "t" "a" "g" "t" "g" "a" "a" "t" "a" "c" "c" "c" "t"
    [4393] "a" "t" "a" "a" "g" "c" "c" "a" "g" "a" "a" "t" "c" "c" "a" "g" "a" "a"
    [4411] "g" "g" "c" "c" "t" "t" "t" "c" "t" "g" "c" "t" "g" "a" "c" "a" "a" "g"
    [4429] "t" "t" "t" "g" "a" "g" "g" "t" "g" "t" "c" "t" "g" "c" "a" "g" "a" "t"
    [4447] "a" "g" "t" "t" "c" "t" "a" "c" "c" "a" "g" "t" "a" "a" "a" "a" "a" "t"
    [4465] "a" "a" "a" "g" "a" "a" "c" "c" "a" "g" "g" "a" "g" "t" "g" "g" "a" "a"
    [4483] "a" "g" "g" "t" "c" "a" "t" "c" "c" "c" "c" "t" "t" "c" "t" "a" "a" "a"
    [4501] "t" "g" "c" "c" "c" "a" "t" "c" "a" "t" "t" "a" "g" "a" "t" "g" "a" "t"
    [4519] "a" "g" "g" "t" "g" "g" "t" "a" "c" "a" "t" "g" "c" "a" "c" "a" "g" "t"
    [4537] "t" "g" "c" "t" "c" "t" "g" "g" "g" "a" "g" "t" "c" "t" "t" "c" "a" "g"
    [4555] "a" "a" "t" "a" "g" "a" "a" "a" "c" "t" "a" "c" "c" "c" "a" "t" "c" "t"
    [4573] "c" "a" "a" "g" "a" "g" "g" "a" "g" "c" "t" "c" "a" "t" "t" "a" "a" "g"
    [4591] "g" "t" "t" "g" "t" "t" "g" "a" "t" "g" "t" "g" "g" "a" "g" "g" "a" "g"
    [4609] "c" "a" "a" "c" "a" "g" "c" "t" "g" "g" "a" "a" "g" "a" "g" "t" "c" "t"
    [4627] "g" "g" "g" "c" "c" "a" "c" "a" "c" "g" "a" "t" "t" "t" "g" "a" "c" "g"
    [4645] "g" "a" "a" "a" "c" "a" "t" "c" "t" "t" "a" "c" "t" "t" "g" "c" "c" "a"
    [4663] "a" "g" "g" "c" "a" "a" "g" "a" "t" "c" "t" "a" "g" "a" "g" "g" "g" "a"
    [4681] "a" "c" "c" "c" "c" "t" "t" "a" "c" "c" "t" "g" "g" "a" "a" "t" "c" "t"
    [4699] "g" "g" "a" "a" "t" "c" "a" "g" "c" "c" "t" "c" "t" "t" "c" "t" "c" "t"
    [4717] "g" "a" "t" "g" "a" "c" "c" "c" "t" "g" "a" "a" "t" "c" "t" "g" "a" "t"
    [4735] "c" "c" "t" "t" "c" "t" "g" "a" "a" "g" "a" "c" "a" "g" "a" "g" "c" "c"
    [4753] "c" "c" "a" "g" "a" "g" "t" "c" "a" "g" "c" "t" "c" "g" "t" "g" "t" "t"
    [4771] "g" "g" "c" "a" "a" "c" "a" "t" "a" "c" "c" "a" "t" "c" "t" "t" "c" "a"
    [4789] "a" "c" "c" "t" "c" "t" "g" "c" "a" "t" "t" "g" "a" "a" "a" "g" "t" "t"
    [4807] "c" "c" "c" "c" "a" "a" "t" "t" "g" "a" "a" "a" "g" "t" "t" "g" "c" "a"
    [4825] "g" "a" "a" "t" "c" "t" "g" "c" "c" "c" "a" "g" "a" "g" "t" "c" "c" "a"
    [4843] "g" "c" "t" "g" "c" "t" "g" "c" "t" "c" "a" "t" "a" "c" "t" "a" "c" "t"
    [4861] "g" "a" "t" "a" "c" "t" "g" "c" "t" "g" "g" "g" "t" "a" "t" "a" "a" "t"
    [4879] "g" "c" "a" "a" "t" "g" "g" "a" "a" "g" "a" "a" "a" "g" "t" "g" "t" "g"
    [4897] "a" "g" "c" "a" "g" "g" "g" "a" "g" "a" "a" "g" "c" "c" "a" "g" "a" "a"
    [4915] "t" "t" "g" "a" "c" "a" "g" "c" "t" "t" "c" "a" "a" "c" "a" "g" "a" "a"
    [4933] "a" "g" "g" "g" "t" "c" "a" "a" "c" "a" "a" "a" "a" "g" "a" "a" "t" "g"
    [4951] "t" "c" "c" "a" "t" "g" "g" "t" "g" "g" "t" "g" "t" "c" "t" "g" "g" "c"
    [4969] "c" "t" "g" "a" "c" "c" "c" "c" "a" "g" "a" "a" "g" "a" "a" "t" "t" "t"
    [4987] "a" "t" "g" "c" "t" "c" "g" "t" "g" "t" "a" "c" "a" "a" "g" "t" "t" "t"
    [5005] "g" "c" "c" "a" "g" "a" "a" "a" "a" "c" "a" "c" "c" "a" "c" "a" "t" "c"
    [5023] "a" "c" "t" "t" "t" "a" "a" "c" "t" "a" "a" "t" "c" "t" "a" "a" "t" "t"
    [5041] "a" "c" "t" "g" "a" "a" "g" "a" "g" "a" "c" "t" "a" "c" "t" "c" "a" "t"
    [5059] "g" "t" "t" "g" "t" "t" "a" "t" "g" "a" "a" "a" "a" "c" "a" "g" "a" "t"
    [5077] "g" "c" "t" "g" "a" "g" "t" "t" "t" "g" "t" "g" "t" "g" "t" "g" "a" "a"
    [5095] "c" "g" "g" "a" "c" "a" "c" "t" "g" "a" "a" "a" "t" "a" "t" "t" "t" "t"
    [5113] "c" "t" "a" "g" "g" "a" "a" "t" "t" "g" "c" "g" "g" "g" "a" "g" "g" "a"
    [5131] "a" "a" "a" "t" "g" "g" "g" "t" "a" "g" "t" "t" "a" "g" "c" "t" "a" "t"
    [5149] "t" "t" "c" "t" "g" "g" "g" "t" "g" "a" "c" "c" "c" "a" "g" "t" "c" "t"
    [5167] "a" "t" "t" "a" "a" "a" "g" "a" "a" "a" "g" "a" "a" "a" "a" "a" "t" "g"
    [5185] "c" "t" "g" "a" "a" "t" "g" "a" "g" "c" "a" "t" "g" "a" "t" "t" "t" "t"
    [5203] "g" "a" "a" "g" "t" "c" "a" "g" "a" "g" "g" "a" "g" "a" "t" "g" "t" "g"
    [5221] "g" "t" "c" "a" "a" "t" "g" "g" "a" "a" "g" "a" "a" "a" "c" "c" "a" "c"
    [5239] "c" "a" "a" "g" "g" "t" "c" "c" "a" "a" "a" "g" "c" "g" "a" "g" "c" "a"
    [5257] "a" "g" "a" "g" "a" "a" "t" "c" "c" "c" "a" "g" "g" "a" "c" "a" "g" "a"
    [5275] "a" "a" "g" "a" "t" "c" "t" "t" "c" "a" "g" "g" "g" "g" "g" "c" "t" "a"
    [5293] "g" "a" "a" "a" "t" "c" "t" "g" "t" "t" "g" "c" "t" "a" "t" "g" "g" "g"
    [5311] "c" "c" "c" "t" "t" "c" "a" "c" "c" "a" "a" "c" "a" "t" "g" "c" "c" "c"
    [5329] "a" "c" "a" "g" "a" "t" "c" "a" "a" "c" "t" "g" "g" "a" "a" "t" "g" "g"
    [5347] "a" "t" "g" "g" "t" "a" "c" "a" "g" "c" "t" "g" "t" "g" "t" "g" "g" "t"
    [5365] "g" "c" "t" "t" "c" "t" "g" "t" "g" "g" "t" "g" "a" "a" "g" "g" "a" "g"
    [5383] "c" "t" "t" "t" "c" "a" "t" "c" "a" "t" "t" "c" "a" "c" "c" "c" "t" "t"
    [5401] "g" "g" "c" "a" "c" "a" "g" "g" "t" "g" "t" "c" "c" "a" "c" "c" "c" "a"
    [5419] "a" "t" "t" "g" "t" "g" "g" "t" "t" "g" "t" "g" "c" "a" "g" "c" "c" "a"
    [5437] "g" "a" "t" "g" "c" "c" "t" "g" "g" "a" "c" "a" "g" "a" "g" "g" "a" "c"
    [5455] "a" "a" "t" "g" "g" "c" "t" "t" "c" "c" "a" "t" "g" "c" "a" "a" "t" "t"
    [5473] "g" "g" "g" "c" "a" "g" "a" "t" "g" "t" "g" "t" "g" "a" "g" "g" "c" "a"
    [5491] "c" "c" "t" "g" "t" "g" "g" "t" "g" "a" "c" "c" "c" "g" "a" "g" "a" "g"
    [5509] "t" "g" "g" "g" "t" "g" "t" "t" "g" "g" "a" "c" "a" "g" "t" "g" "t" "a"
    [5527] "g" "c" "a" "c" "t" "c" "t" "a" "c" "c" "a" "g" "t" "g" "c" "c" "a" "g"
    [5545] "g" "a" "g" "c" "t" "g" "g" "a" "c" "a" "c" "c" "t" "a" "c" "c" "t" "g"
    [5563] "a" "t" "a" "c" "c" "c" "c" "a" "g" "a" "t" "c" "c" "c" "c" "c" "a" "c"
    [5581] "a" "g" "c" "c" "a" "c" "t" "a" "c" "t" "g" "a"
    attr(,"name")
    [1] "DQ190456.BRCA1"
    attr(,"Annot")
    [1] ">DQ190456.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $DQ190457.BRCA1
       [1] "a" "t" "g" "g" "a" "t" "t" "t" "a" "t" "c" "t" "g" "c" "t" "c" "t" "t"
      [19] "c" "g" "c" "g" "t" "t" "g" "a" "a" "g" "a" "a" "g" "t" "a" "c" "a" "a"
      [37] "a" "a" "t" "g" "t" "c" "a" "t" "t" "a" "a" "t" "g" "c" "t" "a" "t" "g"
      [55] "c" "a" "g" "a" "a" "a" "a" "t" "c" "t" "t" "a" "g" "a" "g" "t" "g" "t"
      [73] "c" "c" "c" "a" "t" "c" "t" "g" "t" "c" "t" "g" "g" "a" "g" "t" "t" "g"
      [91] "a" "t" "c" "a" "a" "g" "g" "a" "a" "c" "c" "t" "g" "t" "c" "t" "c" "c"
     [109] "a" "c" "a" "a" "a" "g" "t" "g" "t" "g" "a" "c" "c" "a" "c" "a" "t" "a"
     [127] "t" "t" "t" "t" "g" "c" "a" "a" "a" "t" "t" "t" "t" "g" "c" "a" "t" "g"
     [145] "c" "t" "g" "a" "a" "a" "c" "t" "t" "c" "t" "c" "a" "a" "c" "c" "a" "g"
     [163] "a" "a" "g" "a" "a" "a" "g" "g" "g" "c" "c" "t" "t" "c" "a" "c" "a" "g"
     [181] "t" "g" "t" "c" "c" "t" "t" "t" "a" "t" "g" "t" "a" "a" "g" "a" "a" "t"
     [199] "g" "a" "t" "a" "t" "a" "a" "c" "c" "a" "a" "a" "a" "g" "g" "a" "g" "c"
     [217] "c" "t" "a" "c" "a" "a" "g" "a" "a" "a" "g" "t" "a" "c" "g" "a" "g" "a"
     [235] "t" "t" "t" "a" "g" "t" "c" "a" "a" "c" "t" "t" "g" "t" "t" "g" "a" "a"
     [253] "g" "a" "g" "c" "t" "a" "t" "t" "g" "a" "a" "a" "a" "t" "c" "a" "t" "t"
     [271] "t" "g" "t" "g" "c" "t" "t" "t" "t" "c" "a" "g" "c" "t" "t" "g" "a" "c"
     [289] "a" "c" "a" "g" "g" "t" "t" "t" "g" "g" "a" "g" "t" "a" "t" "g" "c" "a"
     [307] "a" "a" "c" "a" "g" "c" "t" "a" "t" "a" "a" "t" "t" "t" "t" "g" "c" "a"
     [325] "a" "a" "a" "a" "a" "g" "g" "a" "a" "a" "a" "t" "a" "a" "c" "t" "c" "t"
     [343] "c" "c" "t" "g" "a" "a" "c" "a" "t" "c" "t" "a" "a" "a" "a" "g" "a" "t"
     [361] "g" "a" "a" "g" "t" "t" "t" "c" "t" "a" "t" "c" "a" "t" "c" "c" "a" "a"
     [379] "a" "g" "t" "a" "t" "g" "g" "g" "c" "t" "a" "c" "a" "g" "a" "a" "a" "c"
     [397] "c" "g" "t" "g" "c" "c" "a" "a" "a" "a" "g" "a" "c" "t" "t" "c" "t" "a"
     [415] "c" "a" "g" "a" "g" "t" "g" "a" "a" "c" "c" "c" "g" "a" "a" "a" "a" "t"
     [433] "c" "c" "t" "t" "c" "c" "t" "t" "g" "c" "a" "g" "g" "a" "a" "a" "c" "c"
     [451] "a" "g" "t" "c" "t" "c" "a" "g" "t" "g" "t" "c" "c" "a" "a" "c" "t" "c"
     [469] "t" "c" "t" "a" "a" "c" "c" "t" "t" "g" "g" "a" "a" "c" "t" "g" "t" "g"
     [487] "a" "g" "a" "a" "c" "t" "c" "t" "g" "a" "g" "g" "a" "c" "a" "a" "a" "g"
     [505] "c" "a" "g" "c" "g" "g" "a" "t" "a" "c" "a" "a" "c" "c" "t" "c" "a" "a"
     [523] "a" "a" "g" "a" "c" "g" "t" "c" "t" "g" "t" "c" "t" "a" "c" "a" "t" "t"
     [541] "g" "a" "a" "t" "t" "g" "g" "g" "a" "t" "c" "t" "g" "a" "t" "t" "c" "t"
     [559] "t" "c" "t" "g" "a" "a" "g" "a" "t" "a" "c" "c" "g" "t" "t" "a" "a" "t"
     [577] "a" "a" "g" "g" "c" "a" "a" "c" "t" "t" "a" "t" "t" "g" "c" "a" "g" "t"
     [595] "g" "t" "g" "g" "g" "a" "g" "a" "t" "c" "a" "a" "g" "a" "a" "t" "t" "g"
     [613] "t" "t" "a" "c" "a" "a" "a" "t" "c" "a" "c" "c" "c" "c" "t" "c" "a" "a"
     [631] "g" "g" "a" "a" "c" "c" "a" "g" "g" "g" "a" "t" "g" "a" "a" "a" "t" "c"
     [649] "a" "g" "t" "t" "t" "g" "g" "a" "t" "t" "c" "t" "g" "c" "a" "a" "a" "a"
     [667] "a" "a" "g" "g" "c" "t" "g" "c" "t" "t" "g" "t" "g" "a" "a" "t" "t" "t"
     [685] "t" "c" "t" "g" "a" "g" "a" "c" "g" "g" "a" "t" "g" "t" "a" "a" "c" "a"
     [703] "a" "a" "t" "a" "c" "t" "g" "a" "a" "c" "a" "t" "c" "a" "t" "c" "a" "a"
     [721] "c" "c" "c" "a" "g" "t" "a" "a" "t" "a" "a" "t" "g" "a" "t" "t" "t" "g"
     [739] "a" "a" "c" "a" "c" "c" "a" "c" "t" "g" "a" "g" "a" "a" "g" "c" "g" "t"
     [757] "g" "c" "a" "g" "c" "t" "g" "a" "g" "a" "g" "g" "c" "a" "t" "c" "c" "a"
     [775] "g" "a" "a" "a" "a" "g" "t" "a" "t" "c" "a" "g" "g" "g" "t" "a" "g" "t"
     [793] "t" "c" "t" "g" "t" "t" "t" "c" "a" "a" "a" "c" "t" "t" "g" "c" "a" "t"
     [811] "g" "t" "g" "g" "a" "g" "c" "c" "a" "t" "g" "t" "g" "g" "c" "a" "c" "a"
     [829] "a" "a" "t" "a" "c" "t" "c" "a" "t" "g" "c" "c" "a" "g" "c" "t" "c" "a"
     [847] "t" "t" "a" "c" "a" "g" "c" "a" "t" "g" "a" "g" "a" "a" "c" "a" "g" "c"
     [865] "a" "g" "t" "t" "t" "a" "t" "t" "a" "c" "t" "c" "a" "c" "t" "a" "a" "a"
     [883] "g" "a" "c" "a" "g" "a" "a" "t" "g" "a" "a" "t" "g" "t" "a" "g" "a" "a"
     [901] "a" "a" "g" "g" "c" "t" "g" "a" "a" "t" "t" "c" "t" "g" "t" "a" "a" "t"
     [919] "a" "a" "a" "a" "g" "c" "a" "a" "a" "c" "a" "g" "c" "c" "t" "g" "g" "c"
     [937] "t" "t" "a" "g" "c" "a" "a" "g" "g" "a" "g" "c" "c" "a" "a" "c" "a" "t"
     [955] "a" "a" "c" "a" "g" "a" "t" "g" "g" "g" "c" "t" "g" "g" "a" "a" "g" "t"
     [973] "a" "a" "g" "g" "a" "a" "a" "c" "a" "t" "g" "t" "a" "a" "t" "g" "a" "t"
     [991] "a" "g" "g" "c" "g" "g" "a" "c" "t" "c" "c" "c" "a" "g" "c" "a" "c" "a"
    [1009] "g" "a" "a" "a" "a" "a" "a" "a" "g" "g" "t" "a" "g" "a" "t" "c" "t" "g"
    [1027] "a" "a" "t" "g" "c" "t" "g" "a" "t" "c" "c" "c" "c" "t" "g" "t" "g" "t"
    [1045] "g" "a" "g" "a" "g" "a" "a" "a" "a" "g" "a" "a" "t" "g" "g" "a" "a" "t"
    [1063] "a" "a" "g" "c" "a" "g" "a" "a" "a" "c" "t" "g" "c" "c" "a" "t" "g" "c"
    [1081] "t" "c" "a" "g" "a" "g" "a" "a" "t" "c" "c" "t" "a" "g" "a" "g" "a" "t"
    [1099] "a" "c" "t" "g" "a" "a" "g" "a" "t" "g" "t" "t" "c" "c" "t" "t" "g" "g"
    [1117] "a" "t" "a" "a" "c" "a" "c" "t" "a" "a" "a" "t" "a" "g" "c" "a" "g" "c"
    [1135] "a" "t" "t" "c" "a" "g" "a" "a" "a" "g" "t" "t" "a" "a" "t" "g" "a" "g"
    [1153] "t" "g" "g" "t" "t" "t" "t" "c" "c" "a" "g" "a" "a" "g" "t" "g" "a" "t"
    [1171] "g" "a" "a" "c" "t" "g" "t" "t" "a" "g" "g" "t" "t" "c" "t" "g" "a" "t"
    [1189] "g" "a" "c" "t" "c" "a" "c" "a" "t" "g" "a" "t" "g" "g" "g" "g" "a" "g"
    [1207] "t" "c" "t" "g" "a" "a" "t" "c" "a" "a" "a" "t" "g" "c" "c" "a" "a" "a"
    [1225] "g" "t" "a" "g" "c" "t" "g" "a" "t" "g" "t" "a" "t" "t" "g" "g" "a" "c"
    [1243] "g" "t" "t" "c" "t" "a" "a" "a" "t" "g" "a" "g" "g" "t" "a" "g" "a" "t"
    [1261] "g" "a" "a" "t" "a" "t" "t" "c" "t" "g" "g" "t" "t" "c" "t" "t" "c" "a"
    [1279] "g" "a" "g" "a" "a" "a" "a" "t" "a" "g" "a" "c" "t" "t" "a" "c" "t" "g"
    [1297] "g" "c" "c" "a" "g" "t" "g" "a" "t" "c" "c" "t" "c" "a" "t" "g" "a" "g"
    [1315] "g" "c" "t" "t" "t" "a" "a" "t" "a" "t" "g" "t" "a" "a" "a" "a" "g" "t"
    [1333] "g" "a" "a" "a" "g" "a" "g" "t" "t" "c" "a" "c" "t" "c" "c" "a" "a" "a"
    [1351] "t" "c" "a" "g" "t" "a" "g" "a" "g" "a" "g" "t" "a" "a" "t" "a" "t" "t"
    [1369] "g" "a" "a" "g" "a" "c" "a" "a" "a" "a" "t" "a" "t" "t" "t" "g" "g" "g"
    [1387] "a" "a" "a" "a" "c" "c" "t" "a" "t" "c" "g" "g" "a" "a" "g" "a" "a" "g"
    [1405] "g" "c" "a" "a" "g" "c" "c" "t" "c" "c" "c" "c" "a" "a" "c" "t" "t" "a"
    [1423] "a" "g" "c" "c" "a" "t" "g" "t" "a" "a" "c" "t" "g" "a" "a" "a" "a" "t"
    [1441] "c" "t" "a" "a" "t" "t" "a" "t" "a" "g" "g" "a" "g" "c" "a" "t" "t" "t"
    [1459] "g" "t" "t" "a" "c" "t" "g" "a" "g" "c" "c" "a" "c" "a" "g" "a" "t" "a"
    [1477] "a" "t" "a" "c" "a" "a" "g" "a" "g" "c" "g" "t" "c" "c" "c" "c" "t" "c"
    [1495] "a" "c" "a" "a" "a" "t" "a" "a" "a" "t" "t" "a" "a" "a" "g" "c" "g" "t"
    [1513] "a" "a" "a" "a" "g" "g" "a" "g" "a" "c" "c" "t" "a" "c" "a" "t" "c" "a"
    [1531] "g" "g" "c" "c" "t" "t" "c" "a" "t" "c" "c" "t" "g" "a" "g" "g" "a" "t"
    [1549] "t" "t" "t" "a" "t" "c" "a" "a" "g" "a" "a" "a" "g" "c" "a" "g" "a" "t"
    [1567] "t" "t" "g" "g" "c" "a" "g" "t" "t" "c" "a" "a" "a" "a" "g" "a" "c" "t"
    [1585] "c" "c" "t" "g" "a" "a" "a" "t" "g" "a" "t" "a" "a" "a" "t" "c" "a" "g"
    [1603] "g" "g" "a" "a" "c" "t" "a" "a" "c" "c" "a" "a" "a" "c" "g" "g" "a" "g"
    [1621] "c" "a" "g" "a" "a" "t" "g" "g" "t" "c" "a" "a" "g" "t" "g" "a" "t" "g"
    [1639] "a" "a" "t" "a" "t" "t" "a" "c" "t" "a" "a" "t" "a" "g" "t" "g" "g" "t"
    [1657] "c" "a" "t" "g" "a" "g" "a" "a" "t" "a" "a" "a" "a" "c" "a" "a" "a" "a"
    [1675] "g" "g" "t" "g" "a" "t" "t" "c" "t" "a" "t" "t" "c" "a" "g" "a" "a" "t"
    [1693] "g" "a" "g" "a" "a" "a" "a" "a" "t" "c" "c" "t" "a" "a" "c" "c" "c" "a"
    [1711] "a" "t" "a" "g" "a" "a" "t" "c" "a" "c" "t" "c" "g" "a" "a" "a" "a" "a"
    [1729] "g" "a" "a" "t" "c" "t" "g" "c" "t" "t" "t" "c" "a" "a" "a" "a" "c" "g"
    [1747] "a" "a" "a" "g" "c" "t" "g" "a" "a" "c" "c" "t" "a" "t" "a" "a" "g" "c"
    [1765] "a" "g" "c" "a" "g" "t" "a" "t" "a" "a" "g" "c" "a" "a" "t" "a" "t" "g"
    [1783] "g" "a" "a" "c" "t" "c" "g" "a" "a" "t" "t" "a" "a" "a" "t" "a" "t" "c"
    [1801] "c" "a" "c" "a" "a" "t" "t" "c" "a" "a" "a" "a" "g" "c" "a" "c" "c" "t"
    [1819] "a" "a" "a" "a" "a" "g" "a" "a" "t" "a" "g" "g" "c" "t" "g" "a" "g" "g"
    [1837] "a" "g" "g" "a" "a" "g" "t" "c" "t" "t" "c" "t" "a" "c" "c" "a" "g" "g"
    [1855] "c" "a" "t" "a" "t" "t" "c" "a" "t" "g" "c" "g" "c" "t" "t" "g" "a" "a"
    [1873] "c" "t" "a" "g" "t" "a" "g" "t" "c" "a" "g" "t" "a" "g" "a" "a" "a" "t"
    [1891] "c" "t" "a" "a" "g" "c" "c" "c" "a" "c" "c" "t" "a" "a" "t" "t" "g" "t"
    [1909] "a" "c" "t" "g" "a" "a" "t" "t" "g" "c" "a" "a" "a" "t" "t" "g" "a" "t"
    [1927] "a" "g" "t" "t" "g" "t" "t" "c" "t" "a" "g" "c" "a" "g" "t" "g" "a" "a"
    [1945] "g" "a" "g" "a" "t" "a" "a" "a" "g" "a" "a" "a" "a" "a" "a" "a" "a" "g"
    [1963] "t" "a" "c" "a" "a" "c" "c" "a" "a" "a" "t" "g" "c" "c" "a" "g" "t" "c"
    [1981] "a" "g" "g" "c" "a" "c" "a" "g" "c" "a" "g" "a" "a" "a" "c" "c" "t" "a"
    [1999] "c" "a" "a" "c" "t" "c" "a" "t" "g" "g" "a" "a" "g" "g" "t" "a" "a" "a"
    [2017] "g" "a" "a" "c" "c" "t" "g" "c" "a" "a" "c" "t" "g" "g" "a" "g" "c" "c"
    [2035] "a" "a" "g" "a" "a" "g" "a" "g" "t" "a" "a" "c" "a" "a" "g" "c" "c" "a"
    [2053] "a" "a" "t" "g" "a" "a" "c" "a" "g" "a" "c" "a" "a" "g" "t" "a" "a" "a"
    [2071] "a" "g" "a" "c" "a" "t" "g" "a" "c" "a" "g" "c" "g" "a" "t" "a" "c" "t"
    [2089] "t" "t" "c" "c" "c" "a" "g" "a" "g" "c" "t" "g" "a" "a" "g" "t" "t" "a"
    [2107] "a" "c" "a" "a" "a" "t" "g" "c" "a" "c" "c" "t" "g" "g" "t" "t" "c" "t"
    [2125] "t" "t" "t" "a" "c" "t" "a" "a" "g" "t" "g" "t" "t" "c" "a" "a" "a" "t"
    [2143] "a" "c" "c" "a" "g" "t" "g" "a" "a" "c" "t" "t" "a" "a" "a" "g" "a" "a"
    [2161] "t" "t" "t" "g" "t" "c" "a" "a" "t" "c" "c" "t" "a" "g" "c" "c" "t" "t"
    [2179] "c" "c" "a" "a" "g" "a" "g" "a" "a" "g" "a" "a" "a" "a" "a" "g" "a" "a"
    [2197] "g" "a" "g" "a" "a" "a" "c" "t" "a" "g" "a" "a" "a" "c" "a" "g" "t" "t"
    [2215] "a" "a" "a" "g" "t" "g" "t" "c" "t" "a" "a" "t" "a" "a" "t" "g" "c" "t"
    [2233] "g" "a" "a" "g" "a" "c" "c" "c" "c" "a" "a" "a" "g" "a" "t" "c" "t" "c"
    [2251] "a" "t" "g" "t" "t" "a" "a" "g" "t" "g" "g" "a" "g" "a" "a" "a" "g" "g"
    [2269] "g" "t" "t" "t" "t" "g" "c" "a" "a" "a" "c" "t" "g" "a" "a" "a" "g" "a"
    [2287] "t" "c" "t" "g" "t" "a" "g" "a" "g" "a" "g" "t" "a" "g" "c" "a" "g" "t"
    [2305] "a" "t" "t" "t" "c" "a" "t" "t" "g" "g" "t" "a" "c" "c" "t" "g" "g" "t"
    [2323] "a" "c" "t" "g" "a" "t" "t" "a" "t" "g" "g" "c" "a" "c" "t" "c" "a" "g"
    [2341] "g" "a" "a" "a" "g" "t" "a" "t" "c" "t" "c" "g" "t" "t" "a" "c" "t" "g"
    [2359] "g" "a" "a" "g" "t" "t" "a" "g" "c" "a" "c" "t" "c" "t" "a" "g" "g" "g"
    [2377] "a" "a" "g" "g" "c" "a" "a" "a" "a" "a" "c" "a" "g" "a" "a" "c" "c" "a"
    [2395] "a" "a" "t" "a" "a" "a" "t" "g" "t" "g" "t" "g" "a" "g" "t" "c" "a" "g"
    [2413] "t" "g" "t" "g" "c" "a" "g" "c" "a" "t" "t" "t" "g" "a" "a" "a" "a" "c"
    [2431] "c" "c" "c" "a" "a" "g" "g" "g" "a" "c" "t" "a" "a" "t" "t" "c" "a" "t"
    [2449] "g" "g" "t" "t" "g" "t" "t" "c" "c" "a" "a" "a" "g" "a" "t" "a" "a" "t"
    [2467] "a" "g" "a" "a" "a" "t" "g" "a" "c" "a" "c" "a" "g" "a" "a" "g" "g" "c"
    [2485] "t" "t" "t" "a" "a" "g" "t" "a" "t" "c" "c" "a" "t" "t" "g" "g" "g" "a"
    [2503] "c" "a" "t" "g" "a" "a" "g" "t" "t" "a" "a" "c" "c" "a" "c" "a" "g" "t"
    [2521] "c" "g" "g" "g" "a" "a" "a" "c" "a" "a" "g" "c" "a" "t" "a" "g" "a" "a"
    [2539] "a" "t" "g" "g" "a" "a" "g" "a" "a" "a" "g" "t" "g" "a" "a" "c" "t" "t"
    [2557] "g" "a" "t" "g" "c" "t" "c" "a" "g" "t" "a" "t" "t" "t" "g" "c" "a" "g"
    [2575] "a" "a" "t" "a" "c" "a" "t" "t" "c" "a" "a" "g" "g" "t" "t" "t" "c" "a"
    [2593] "a" "a" "g" "c" "g" "c" "c" "a" "g" "t" "c" "a" "t" "t" "t" "g" "c" "t"
    [2611] "c" "c" "g" "t" "t" "t" "t" "c" "a" "a" "a" "t" "c" "c" "a" "g" "g" "a"
    [2629] "a" "a" "t" "g" "c" "a" "g" "a" "a" "g" "a" "g" "g" "a" "a" "t" "g" "t"
    [2647] "g" "c" "a" "a" "c" "a" "t" "t" "c" "t" "c" "t" "g" "c" "c" "c" "a" "c"
    [2665] "t" "c" "t" "g" "g" "g" "t" "c" "c" "t" "t" "a" "a" "a" "g" "a" "a" "a"
    [2683] "c" "a" "a" "a" "g" "t" "c" "c" "a" "a" "a" "a" "g" "t" "c" "a" "c" "t"
    [2701] "t" "t" "t" "g" "a" "a" "t" "g" "t" "g" "a" "a" "c" "a" "a" "a" "a" "g"
    [2719] "g" "a" "a" "g" "a" "a" "a" "a" "t" "c" "a" "a" "g" "g" "a" "a" "a" "g"
    [2737] "a" "a" "t" "g" "a" "g" "t" "c" "t" "a" "a" "t" "a" "t" "c" "a" "a" "g"
    [2755] "c" "c" "t" "g" "t" "a" "c" "a" "g" "a" "c" "a" "g" "t" "t" "a" "a" "t"
    [2773] "a" "t" "c" "a" "c" "t" "g" "c" "a" "g" "g" "c" "t" "t" "t" "c" "c" "t"
    [2791] "g" "t" "g" "g" "t" "t" "g" "g" "t" "c" "a" "g" "a" "a" "a" "g" "a" "t"
    [2809] "a" "a" "g" "c" "c" "a" "g" "t" "t" "g" "a" "t" "a" "a" "t" "g" "c" "c"
    [2827] "a" "a" "a" "t" "g" "t" "a" "g" "t" "a" "t" "c" "a" "a" "a" "g" "g" "a"
    [2845] "g" "g" "c" "t" "c" "t" "a" "g" "g" "t" "t" "t" "t" "g" "t" "c" "t" "a"
    [2863] "t" "c" "a" "t" "c" "t" "c" "a" "g" "t" "t" "c" "a" "g" "a" "g" "g" "c"
    [2881] "a" "a" "c" "g" "a" "a" "a" "c" "t" "g" "g" "a" "c" "t" "c" "a" "t" "t"
    [2899] "a" "c" "t" "c" "c" "a" "a" "a" "t" "a" "a" "a" "c" "a" "t" "g" "g" "a"
    [2917] "c" "t" "t" "t" "t" "a" "c" "a" "a" "a" "a" "c" "c" "c" "a" "t" "a" "t"
    [2935] "c" "g" "t" "a" "t" "a" "c" "c" "a" "c" "c" "a" "c" "t" "t" "t" "t" "t"
    [2953] "c" "c" "c" "a" "t" "c" "a" "a" "g" "t" "c" "a" "t" "t" "t" "g" "t" "t"
    [2971] "a" "a" "a" "a" "c" "t" "a" "a" "a" "t" "g" "t" "a" "a" "g" "a" "a" "a"
    [2989] "a" "a" "t" "c" "t" "g" "c" "t" "a" "g" "a" "g" "g" "a" "a" "a" "a" "c"
    [3007] "t" "t" "t" "g" "a" "g" "g" "a" "a" "c" "a" "t" "t" "c" "a" "a" "t" "g"
    [3025] "t" "c" "a" "c" "c" "t" "g" "a" "a" "a" "g" "a" "g" "a" "a" "a" "t" "g"
    [3043] "g" "g" "a" "a" "a" "t" "g" "a" "g" "a" "a" "c" "a" "t" "t" "c" "c" "a"
    [3061] "a" "g" "t" "a" "c" "a" "g" "t" "g" "a" "g" "c" "a" "c" "a" "a" "t" "t"
    [3079] "a" "g" "c" "c" "g" "t" "a" "a" "t" "a" "a" "c" "a" "t" "t" "a" "g" "a"
    [3097] "g" "a" "a" "a" "a" "t" "g" "t" "t" "t" "t" "t" "a" "a" "a" "g" "a" "a"
    [3115] "g" "c" "c" "a" "g" "c" "t" "c" "a" "a" "g" "c" "a" "a" "t" "a" "t" "t"
    [3133] "a" "a" "t" "g" "a" "a" "g" "t" "a" "g" "g" "t" "t" "c" "c" "a" "g" "t"
    [3151] "a" "c" "t" "a" "a" "t" "g" "a" "a" "g" "t" "g" "g" "g" "c" "t" "c" "c"
    [3169] "a" "g" "t" "a" "t" "t" "a" "a" "t" "g" "a" "a" "a" "t" "a" "g" "g" "t"
    [3187] "t" "c" "c" "a" "g" "t" "g" "a" "t" "g" "a" "a" "a" "a" "c" "a" "t" "t"
    [3205] "c" "a" "a" "g" "c" "a" "g" "a" "a" "c" "t" "a" "g" "g" "t" "a" "g" "a"
    [3223] "a" "a" "c" "a" "g" "a" "g" "g" "g" "c" "c" "a" "a" "a" "a" "t" "t" "g"
    [3241] "a" "a" "t" "g" "c" "t" "a" "t" "g" "c" "t" "t" "a" "g" "a" "t" "t" "a"
    [3259] "g" "g" "g" "g" "t" "t" "t" "t" "g" "c" "a" "a" "c" "c" "t" "g" "a" "g"
    [3277] "g" "t" "c" "t" "a" "t" "a" "a" "a" "c" "a" "a" "a" "g" "t" "c" "t" "t"
    [3295] "c" "c" "t" "g" "g" "a" "a" "g" "t" "a" "a" "t" "t" "g" "t" "a" "a" "g"
    [3313] "c" "a" "t" "c" "c" "t" "g" "a" "a" "a" "t" "a" "a" "a" "a" "a" "a" "g"
    [3331] "c" "a" "a" "g" "a" "a" "t" "a" "t" "g" "a" "a" "g" "a" "a" "g" "t" "a"
    [3349] "g" "t" "t" "c" "a" "g" "a" "c" "t" "g" "t" "t" "a" "a" "t" "a" "c" "a"
    [3367] "g" "a" "t" "t" "t" "c" "t" "c" "t" "c" "c" "a" "t" "a" "t" "c" "t" "g"
    [3385] "a" "t" "t" "t" "c" "a" "g" "a" "t" "a" "a" "c" "t" "t" "a" "g" "a" "a"
    [3403] "c" "a" "g" "c" "c" "t" "a" "t" "g" "g" "g" "a" "a" "g" "t" "a" "g" "t"
    [3421] "c" "a" "t" "g" "c" "a" "t" "c" "t" "c" "a" "g" "g" "t" "t" "t" "g" "t"
    [3439] "t" "c" "t" "g" "a" "g" "a" "c" "a" "c" "c" "t" "g" "a" "t" "g" "a" "c"
    [3457] "c" "t" "g" "t" "t" "a" "g" "a" "t" "g" "a" "t" "g" "g" "t" "g" "a" "a"
    [3475] "a" "t" "a" "a" "a" "g" "g" "a" "a" "g" "a" "t" "a" "c" "t" "a" "g" "t"
    [3493] "t" "t" "t" "g" "c" "t" "g" "a" "a" "a" "a" "t" "g" "a" "c" "a" "t" "t"
    [3511] "a" "a" "g" "g" "a" "a" "a" "g" "t" "t" "c" "t" "g" "c" "t" "g" "t" "t"
    [3529] "t" "t" "t" "a" "g" "c" "a" "a" "a" "a" "g" "c" "g" "t" "c" "c" "a" "g"
    [3547] "a" "a" "a" "g" "g" "a" "g" "a" "g" "c" "t" "t" "a" "g" "c" "a" "g" "g"
    [3565] "a" "g" "t" "c" "c" "t" "a" "g" "c" "c" "c" "t" "t" "t" "c" "a" "c" "c"
    [3583] "c" "a" "t" "a" "c" "a" "c" "a" "t" "t" "t" "g" "g" "c" "t" "c" "a" "g"
    [3601] "g" "g" "t" "t" "a" "c" "c" "g" "a" "a" "g" "a" "g" "g" "g" "g" "c" "c"
    [3619] "a" "a" "g" "a" "a" "a" "t" "t" "a" "g" "a" "g" "t" "c" "c" "t" "c" "a"
    [3637] "g" "a" "a" "g" "a" "g" "a" "a" "c" "t" "t" "a" "t" "c" "t" "a" "g" "t"
    [3655] "g" "a" "g" "g" "a" "t" "g" "a" "a" "g" "a" "g" "c" "t" "t" "c" "c" "c"
    [3673] "t" "g" "c" "t" "t" "c" "c" "a" "a" "c" "a" "c" "t" "t" "g" "t" "t" "a"
    [3691] "t" "t" "t" "g" "g" "t" "a" "a" "a" "g" "t" "a" "a" "a" "c" "a" "a" "t"
    [3709] "a" "t" "a" "c" "c" "t" "t" "c" "t" "c" "a" "g" "t" "c" "t" "a" "c" "t"
    [3727] "a" "g" "g" "c" "a" "t" "a" "g" "c" "a" "c" "c" "g" "t" "t" "g" "c" "t"
    [3745] "a" "c" "c" "g" "a" "g" "t" "g" "t" "c" "t" "g" "t" "c" "t" "a" "a" "g"
    [3763] "a" "a" "c" "a" "c" "a" "g" "a" "g" "g" "a" "g" "a" "a" "t" "t" "t" "a"
    [3781] "t" "t" "a" "t" "c" "a" "t" "t" "g" "a" "a" "g" "a" "a" "t" "a" "g" "c"
    [3799] "t" "t" "a" "a" "a" "t" "g" "a" "c" "t" "g" "c" "a" "g" "t" "a" "a" "c"
    [3817] "c" "a" "g" "g" "t" "a" "a" "t" "a" "t" "t" "g" "g" "c" "a" "a" "a" "g"
    [3835] "g" "c" "a" "t" "c" "t" "c" "a" "g" "g" "a" "a" "c" "a" "t" "c" "a" "c"
    [3853] "c" "t" "t" "a" "g" "t" "g" "a" "g" "g" "a" "a" "a" "c" "a" "a" "a" "a"
    [3871] "t" "g" "t" "t" "c" "t" "g" "c" "t" "a" "g" "c" "t" "t" "g" "t" "t" "t"
    [3889] "t" "c" "t" "t" "c" "a" "c" "a" "g" "t" "g" "c" "a" "g" "t" "g" "a" "a"
    [3907] "t" "t" "g" "g" "a" "a" "g" "a" "c" "t" "t" "g" "a" "c" "t" "g" "c" "a"
    [3925] "a" "a" "t" "a" "c" "a" "a" "a" "c" "a" "c" "c" "c" "a" "g" "g" "a" "t"
    [3943] "c" "c" "t" "t" "t" "c" "t" "t" "g" "a" "t" "t" "g" "g" "t" "t" "c" "t"
    [3961] "t" "c" "c" "a" "a" "a" "c" "a" "a" "a" "t" "g" "a" "g" "g" "c" "a" "t"
    [3979] "c" "a" "g" "t" "c" "t" "g" "a" "a" "a" "g" "c" "c" "a" "g" "g" "g" "a"
    [3997] "g" "t" "t" "g" "g" "t" "c" "t" "g" "a" "g" "t" "g" "a" "c" "a" "a" "g"
    [4015] "g" "a" "a" "t" "t" "g" "g" "t" "t" "t" "c" "a" "g" "a" "t" "g" "a" "t"
    [4033] "g" "a" "a" "g" "a" "a" "a" "g" "a" "g" "g" "a" "a" "c" "g" "g" "g" "c"
    [4051] "t" "t" "g" "g" "a" "a" "g" "a" "a" "a" "a" "t" "a" "a" "t" "c" "a" "a"
    [4069] "g" "a" "a" "g" "a" "g" "c" "a" "a" "a" "g" "c" "a" "t" "g" "g" "a" "t"
    [4087] "t" "c" "a" "a" "a" "c" "t" "t" "a" "g" "g" "t" "g" "a" "a" "g" "c" "a"
    [4105] "g" "c" "a" "t" "c" "t" "g" "g" "g" "t" "g" "t" "g" "a" "g" "a" "g" "t"
    [4123] "g" "a" "a" "a" "c" "a" "a" "g" "c" "g" "t" "c" "t" "c" "t" "g" "a" "a"
    [4141] "g" "a" "c" "t" "g" "c" "t" "c" "a" "g" "g" "g" "c" "t" "a" "t" "c" "c"
    [4159] "t" "c" "t" "c" "a" "g" "a" "g" "t" "g" "a" "c" "a" "t" "t" "t" "t" "a"
    [4177] "a" "c" "c" "a" "c" "t" "c" "a" "g" "c" "a" "g" "a" "g" "g" "g" "a" "t"
    [4195] "a" "c" "c" "a" "t" "g" "c" "a" "a" "c" "a" "t" "a" "a" "c" "c" "t" "g"
    [4213] "a" "t" "a" "a" "a" "g" "c" "t" "c" "c" "a" "g" "c" "a" "g" "g" "a" "a"
    [4231] "a" "t" "g" "g" "c" "t" "g" "a" "a" "c" "t" "a" "g" "a" "a" "g" "c" "t"
    [4249] "g" "t" "g" "t" "t" "a" "g" "a" "a" "c" "a" "g" "c" "a" "t" "g" "g" "g"
    [4267] "a" "g" "c" "c" "a" "g" "c" "c" "t" "t" "c" "t" "a" "a" "c" "a" "g" "c"
    [4285] "t" "a" "c" "c" "c" "t" "t" "c" "c" "a" "t" "c" "a" "t" "a" "a" "g" "t"
    [4303] "g" "a" "c" "t" "c" "t" "t" "c" "t" "g" "c" "c" "c" "t" "t" "g" "a" "g"
    [4321] "g" "a" "c" "c" "t" "g" "c" "g" "a" "a" "a" "t" "c" "c" "a" "g" "a" "a"
    [4339] "c" "a" "a" "a" "g" "c" "a" "c" "a" "t" "c" "a" "g" "a" "a" "a" "a" "a"
    [4357] "g" "c" "a" "g" "t" "a" "t" "t" "a" "a" "c" "t" "t" "c" "a" "c" "a" "g"
    [4375] "a" "a" "a" "a" "g" "t" "a" "g" "t" "g" "a" "a" "t" "a" "c" "c" "c" "t"
    [4393] "a" "t" "a" "a" "g" "c" "c" "a" "g" "a" "a" "t" "c" "c" "a" "g" "a" "a"
    [4411] "g" "g" "c" "c" "t" "t" "t" "c" "t" "g" "c" "t" "g" "a" "c" "a" "a" "g"
    [4429] "t" "t" "t" "g" "a" "g" "g" "t" "g" "t" "c" "t" "g" "c" "a" "g" "a" "t"
    [4447] "a" "g" "t" "t" "c" "t" "a" "c" "c" "a" "g" "t" "a" "a" "a" "a" "a" "t"
    [4465] "a" "a" "a" "g" "a" "a" "c" "c" "a" "g" "g" "a" "g" "t" "g" "g" "a" "a"
    [4483] "a" "g" "g" "t" "c" "a" "t" "c" "c" "c" "c" "t" "t" "c" "t" "a" "a" "a"
    [4501] "t" "g" "c" "c" "c" "a" "t" "c" "a" "t" "t" "a" "g" "a" "t" "g" "a" "t"
    [4519] "a" "g" "g" "t" "g" "g" "t" "a" "c" "a" "t" "g" "c" "a" "c" "a" "g" "t"
    [4537] "t" "g" "c" "t" "c" "t" "g" "g" "g" "a" "g" "t" "c" "t" "t" "c" "a" "g"
    [4555] "a" "a" "t" "a" "g" "a" "a" "a" "c" "t" "a" "c" "c" "c" "a" "t" "c" "t"
    [4573] "c" "a" "a" "g" "a" "g" "g" "a" "g" "c" "t" "c" "a" "t" "t" "a" "a" "g"
    [4591] "g" "t" "t" "g" "t" "t" "g" "a" "t" "g" "t" "g" "g" "a" "g" "g" "a" "g"
    [4609] "c" "a" "a" "c" "a" "g" "c" "t" "g" "g" "a" "a" "g" "a" "g" "t" "c" "t"
    [4627] "g" "g" "g" "c" "c" "a" "c" "a" "c" "g" "a" "t" "t" "t" "g" "a" "c" "g"
    [4645] "g" "a" "a" "a" "c" "a" "t" "c" "t" "t" "a" "c" "t" "t" "g" "c" "c" "a"
    [4663] "a" "g" "g" "c" "a" "a" "g" "a" "t" "c" "t" "a" "g" "a" "g" "g" "g" "a"
    [4681] "a" "c" "c" "c" "c" "t" "t" "a" "c" "c" "t" "g" "g" "a" "a" "t" "c" "t"
    [4699] "g" "g" "a" "a" "t" "c" "a" "g" "c" "c" "t" "c" "t" "t" "c" "t" "c" "t"
    [4717] "g" "a" "t" "g" "a" "c" "c" "c" "t" "g" "a" "a" "t" "c" "t" "g" "a" "t"
    [4735] "c" "c" "t" "t" "c" "t" "g" "a" "a" "g" "a" "c" "a" "g" "a" "g" "c" "c"
    [4753] "c" "c" "a" "g" "a" "g" "t" "c" "a" "g" "c" "t" "c" "g" "t" "g" "t" "t"
    [4771] "g" "g" "c" "a" "a" "c" "a" "t" "a" "c" "c" "a" "t" "c" "t" "t" "c" "a"
    [4789] "a" "c" "c" "t" "c" "t" "g" "c" "a" "t" "t" "g" "a" "a" "a" "g" "t" "t"
    [4807] "c" "c" "c" "c" "a" "a" "t" "t" "g" "a" "a" "a" "g" "t" "t" "g" "c" "a"
    [4825] "g" "a" "a" "t" "c" "t" "g" "c" "c" "c" "a" "g" "a" "g" "t" "c" "c" "a"
    [4843] "g" "c" "t" "g" "c" "t" "g" "c" "t" "c" "a" "t" "a" "c" "t" "a" "c" "t"
    [4861] "g" "a" "t" "a" "c" "t" "g" "c" "t" "g" "g" "g" "t" "a" "t" "a" "a" "t"
    [4879] "g" "c" "a" "a" "t" "g" "g" "a" "a" "g" "a" "a" "a" "g" "t" "g" "t" "g"
    [4897] "a" "g" "c" "a" "g" "g" "g" "a" "g" "a" "a" "g" "c" "c" "a" "g" "a" "a"
    [4915] "t" "t" "g" "a" "c" "a" "g" "c" "t" "t" "c" "a" "a" "c" "a" "g" "a" "a"
    [4933] "a" "g" "g" "g" "t" "c" "a" "a" "c" "a" "a" "a" "a" "g" "a" "a" "t" "g"
    [4951] "t" "c" "c" "a" "t" "g" "g" "t" "g" "g" "t" "g" "t" "c" "t" "g" "g" "c"
    [4969] "c" "t" "g" "a" "c" "c" "c" "c" "a" "g" "a" "a" "g" "a" "a" "t" "t" "t"
    [4987] "a" "t" "g" "c" "t" "c" "g" "t" "g" "t" "a" "c" "a" "a" "g" "t" "t" "t"
    [5005] "g" "c" "c" "a" "g" "a" "a" "a" "a" "c" "a" "c" "c" "a" "c" "a" "t" "c"
    [5023] "a" "c" "t" "t" "t" "a" "a" "c" "t" "a" "a" "t" "c" "t" "a" "a" "t" "t"
    [5041] "a" "c" "t" "g" "a" "a" "g" "a" "g" "a" "c" "t" "a" "c" "t" "c" "a" "t"
    [5059] "g" "t" "t" "g" "t" "t" "a" "t" "g" "a" "a" "a" "a" "c" "a" "g" "a" "t"
    [5077] "g" "c" "t" "g" "a" "g" "t" "t" "t" "g" "t" "g" "t" "g" "t" "g" "a" "a"
    [5095] "c" "g" "g" "a" "c" "a" "c" "t" "g" "a" "a" "a" "t" "a" "t" "t" "t" "t"
    [5113] "c" "t" "a" "g" "g" "a" "a" "t" "t" "g" "c" "g" "g" "g" "a" "g" "g" "a"
    [5131] "a" "a" "a" "t" "g" "g" "g" "t" "a" "g" "t" "t" "a" "g" "c" "t" "a" "t"
    [5149] "t" "t" "c" "t" "g" "g" "g" "t" "g" "a" "c" "c" "c" "a" "g" "t" "c" "t"
    [5167] "a" "t" "t" "a" "a" "a" "g" "a" "a" "a" "g" "a" "a" "a" "a" "a" "t" "g"
    [5185] "c" "t" "g" "a" "a" "t" "g" "a" "g" "c" "a" "t" "g" "a" "t" "t" "t" "t"
    [5203] "g" "a" "a" "g" "t" "c" "a" "g" "a" "g" "g" "a" "g" "a" "t" "g" "t" "g"
    [5221] "g" "t" "c" "a" "a" "t" "g" "g" "a" "a" "g" "a" "a" "a" "c" "c" "a" "c"
    [5239] "c" "a" "a" "g" "g" "t" "c" "c" "a" "a" "a" "g" "c" "g" "a" "g" "c" "a"
    [5257] "a" "g" "a" "g" "a" "a" "t" "c" "c" "c" "a" "g" "g" "a" "c" "a" "g" "a"
    [5275] "a" "a" "g" "a" "t" "c" "t" "t" "c" "a" "g" "g" "g" "g" "g" "c" "t" "a"
    [5293] "g" "a" "a" "a" "t" "c" "t" "g" "t" "t" "g" "c" "t" "a" "t" "g" "g" "g"
    [5311] "c" "c" "c" "t" "t" "c" "a" "c" "c" "a" "a" "c" "a" "t" "g" "c" "c" "c"
    [5329] "a" "c" "a" "g" "a" "t" "c" "a" "a" "c" "t" "g" "g" "a" "a" "t" "g" "g"
    [5347] "a" "t" "g" "g" "t" "a" "c" "a" "g" "c" "t" "g" "t" "g" "t" "g" "g" "t"
    [5365] "g" "c" "t" "t" "c" "t" "g" "t" "g" "g" "t" "g" "a" "a" "g" "g" "a" "g"
    [5383] "c" "t" "t" "t" "c" "a" "t" "c" "a" "t" "t" "c" "a" "c" "c" "c" "t" "t"
    [5401] "g" "g" "c" "a" "c" "a" "g" "g" "t" "g" "t" "c" "c" "a" "c" "c" "c" "a"
    [5419] "a" "t" "t" "g" "t" "g" "g" "t" "t" "g" "t" "g" "c" "a" "g" "c" "c" "a"
    [5437] "g" "a" "t" "g" "c" "c" "t" "g" "g" "a" "c" "a" "g" "a" "g" "g" "a" "c"
    [5455] "a" "a" "t" "g" "g" "c" "t" "t" "c" "c" "a" "t" "g"
    attr(,"name")
    [1] "DQ190457.BRCA1"
    attr(,"Annot")
    [1] ">DQ190457.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $DQ299305
      [1] "a" "a" "a" "t" "g" "a" "c" "a" "t" "t" "a" "a" "g" "g" "a" "a" "a" "g"
     [19] "t" "t" "c" "t" "g" "c" "t" "g" "t" "t" "t" "t" "t" "a" "g" "c" "a" "a"
     [37] "a" "a" "g" "c" "g" "t" "c" "c" "a" "g" "a" "a" "a" "g" "g" "a" "g" "a"
     [55] "g" "c" "t" "t" "a" "g" "c" "a" "g" "g" "a" "g" "t" "c" "c" "t" "a" "g"
     [73] "c" "c" "c" "t" "t" "t" "c" "a" "c" "c" "c" "a" "t" "a" "c" "a" "c" "a"
     [91] "t" "t" "t" "g" "g" "c" "t" "c" "a" "g" "g" "g" "t" "t" "a" "c" "c" "g"
    [109] "a" "a" "g" "a" "g" "g" "g" "g" "c" "c" "a" "a"
    attr(,"name")
    [1] "DQ299305"
    attr(,"Annot")
    [1] ">DQ299305"
    attr(,"class")
    [1] "SeqFastadna"
    
    $DQ299306
      [1] "a" "a" "a" "t" "g" "a" "c" "a" "t" "t" "a" "a" "g" "g" "a" "a" "a" "g"
     [19] "t" "t" "c" "t" "g" "c" "t" "g" "t" "t" "t" "t" "t" "a" "g" "c" "a" "a"
     [37] "a" "a" "g" "c" "g" "t" "c" "c" "a" "g" "a" "a" "a" "g" "g" "a" "g" "a"
     [55] "g" "c" "t" "t" "a" "g" "c" "a" "g" "g" "a" "g" "t" "c" "c" "t" "a" "g"
     [73] "c" "c" "c" "t" "t" "t" "c" "a" "c" "c" "c" "a" "t" "a" "c" "a" "c" "a"
     [91] "t" "t" "t" "g" "g" "c" "t" "c" "a" "g" "g" "g" "t" "t" "a" "c" "c" "g"
    [109] "a" "a" "g" "a" "g" "g" "g" "g" "c" "c" "a" "a"
    attr(,"name")
    [1] "DQ299306"
    attr(,"Annot")
    [1] ">DQ299306"
    attr(,"class")
    [1] "SeqFastadna"
    
    $DQ299307
      [1] "a" "a" "a" "t" "g" "a" "c" "a" "t" "t" "a" "a" "g" "g" "a" "a" "a" "g"
     [19] "t" "t" "c" "t" "g" "c" "t" "g" "t" "t" "t" "t" "t" "a" "g" "c" "a" "a"
     [37] "a" "a" "g" "c" "g" "t" "c" "c" "a" "g" "a" "a" "a" "g" "g" "a" "g" "a"
     [55] "g" "c" "t" "t" "a" "g" "c" "a" "g" "g" "a" "g" "t" "c" "c" "t" "a" "g"
     [73] "c" "c" "c" "t" "t" "t" "c" "a" "c" "c" "c" "a" "t" "a" "c" "a" "c" "a"
     [91] "t" "t" "t" "g" "g" "c" "t" "c" "a" "g" "g" "g" "t" "t" "a" "c" "c" "g"
    [109] "a" "a" "g" "a" "g" "g" "g" "g" "c" "c" "a" "a"
    attr(,"name")
    [1] "DQ299307"
    attr(,"Annot")
    [1] ">DQ299307"
    attr(,"class")
    [1] "SeqFastadna"
    
    $DQ299308
      [1] "a" "a" "a" "t" "g" "a" "c" "a" "t" "t" "a" "a" "g" "g" "a" "a" "a" "g"
     [19] "t" "t" "c" "t" "g" "c" "t" "g" "t" "t" "t" "t" "t" "a" "g" "c" "a" "a"
     [37] "a" "a" "g" "c" "g" "t" "c" "c" "a" "g" "a" "a" "a" "g" "g" "a" "g" "a"
     [55] "g" "c" "t" "t" "a" "g" "c" "a" "g" "g" "a" "g" "t" "c" "c" "t" "a" "g"
     [73] "c" "c" "c" "t" "t" "t" "c" "a" "c" "c" "c" "a" "t" "a" "c" "a" "c" "a"
     [91] "t" "t" "t" "g" "g" "c" "t" "c" "a" "g" "g" "g" "t" "t" "a" "c" "c" "g"
    [109] "a" "a" "g" "a" "g" "g" "g" "g" "c" "c" "a" "a"
    attr(,"name")
    [1] "DQ299308"
    attr(,"Annot")
    [1] ">DQ299308"
    attr(,"class")
    [1] "SeqFastadna"
    
    $DQ299309
      [1] "a" "a" "a" "t" "g" "a" "c" "a" "t" "t" "a" "a" "g" "g" "a" "a" "a" "g"
     [19] "t" "t" "c" "t" "g" "c" "t" "g" "t" "t" "t" "t" "t" "a" "g" "c" "a" "a"
     [37] "a" "a" "g" "c" "g" "t" "c" "c" "a" "g" "a" "a" "a" "g" "g" "a" "g" "a"
     [55] "g" "c" "t" "t" "a" "g" "c" "a" "g" "g" "a" "g" "t" "c" "c" "t" "a" "g"
     [73] "c" "c" "c" "t" "t" "t" "c" "a" "c" "c" "c" "a" "t" "a" "c" "a" "c" "a"
     [91] "t" "t" "t" "g" "g" "c" "t" "c" "a" "g" "g" "g" "t" "t" "a" "c" "c" "g"
    [109] "a" "a" "g" "a" "g" "g" "g" "g" "c" "c" "a" "a"
    attr(,"name")
    [1] "DQ299309"
    attr(,"Annot")
    [1] ">DQ299309"
    attr(,"class")
    [1] "SeqFastadna"
    
    $DQ299310
      [1] "a" "a" "a" "t" "g" "a" "c" "a" "t" "t" "a" "a" "g" "g" "a" "a" "a" "g"
     [19] "t" "t" "c" "t" "g" "c" "t" "g" "t" "t" "t" "t" "t" "a" "g" "c" "a" "a"
     [37] "a" "a" "g" "c" "g" "t" "c" "c" "a" "g" "a" "a" "a" "g" "g" "a" "g" "a"
     [55] "g" "c" "t" "t" "a" "g" "c" "a" "g" "g" "a" "g" "t" "c" "c" "t" "a" "g"
     [73] "c" "c" "c" "t" "t" "t" "c" "a" "c" "c" "c" "a" "t" "a" "c" "a" "c" "a"
     [91] "t" "t" "t" "g" "g" "c" "t" "c" "a" "g" "g" "g" "t" "t" "a" "c" "c" "g"
    [109] "a" "a" "g" "a" "g" "g" "g" "g" "c" "c" "a" "a"
    attr(,"name")
    [1] "DQ299310"
    attr(,"Annot")
    [1] ">DQ299310"
    attr(,"class")
    [1] "SeqFastadna"
    
    $DQ299311
      [1] "a" "a" "a" "t" "g" "a" "c" "a" "t" "t" "a" "a" "g" "g" "a" "a" "a" "g"
     [19] "t" "t" "c" "t" "g" "c" "t" "g" "t" "t" "t" "t" "t" "a" "g" "c" "a" "a"
     [37] "a" "a" "g" "c" "g" "t" "c" "c" "a" "g" "a" "a" "a" "g" "g" "a" "g" "a"
     [55] "g" "c" "t" "t" "a" "g" "c" "a" "g" "g" "a" "g" "t" "c" "c" "t" "a" "g"
     [73] "c" "c" "c" "t" "t" "t" "c" "a" "c" "c" "c" "a" "t" "a" "c" "a" "c" "a"
     [91] "t" "t" "t" "g" "g" "c" "t" "c" "a" "g" "g" "g" "t" "t" "a" "c" "c" "g"
    [109] "a" "a" "g" "a" "g" "g" "g" "g" "c" "c" "a" "a"
    attr(,"name")
    [1] "DQ299311"
    attr(,"Annot")
    [1] ">DQ299311"
    attr(,"class")
    [1] "SeqFastadna"
    
    $DQ299312
      [1] "a" "a" "a" "t" "g" "a" "c" "a" "t" "t" "a" "a" "g" "g" "a" "a" "a" "g"
     [19] "t" "t" "c" "t" "g" "c" "t" "g" "t" "t" "t" "t" "t" "a" "g" "c" "a" "a"
     [37] "a" "a" "g" "c" "g" "t" "c" "c" "a" "g" "a" "a" "a" "g" "g" "a" "g" "a"
     [55] "g" "c" "t" "t" "a" "g" "c" "a" "g" "g" "a" "g" "t" "c" "c" "t" "a" "g"
     [73] "c" "c" "c" "t" "t" "t" "c" "a" "c" "c" "c" "a" "t" "a" "c" "a" "c" "a"
     [91] "t" "t" "t" "g" "g" "c" "t" "c" "a" "g" "g" "g" "t" "t" "a" "c" "c" "g"
    [109] "a" "a" "g" "a" "g" "g" "g" "g" "c" "c" "a" "a"
    attr(,"name")
    [1] "DQ299312"
    attr(,"Annot")
    [1] ">DQ299312"
    attr(,"class")
    [1] "SeqFastadna"
    
    $DQ299313
      [1] "a" "a" "a" "t" "g" "a" "c" "a" "t" "t" "a" "a" "g" "g" "a" "a" "a" "g"
     [19] "t" "t" "c" "t" "g" "c" "t" "g" "t" "t" "t" "t" "t" "a" "g" "c" "a" "a"
     [37] "a" "a" "g" "c" "g" "t" "c" "c" "a" "g" "a" "a" "a" "g" "g" "a" "g" "a"
     [55] "g" "c" "t" "t" "a" "g" "c" "a" "g" "g" "a" "g" "t" "c" "c" "t" "a" "g"
     [73] "c" "c" "c" "t" "t" "t" "c" "a" "c" "c" "c" "a" "t" "a" "c" "a" "c" "a"
     [91] "t" "t" "t" "g" "g" "c" "t" "c" "a" "g" "g" "g" "t" "t" "a" "c" "c" "g"
    [109] "a" "a" "g" "a" "g" "g" "g" "g" "c" "c" "a" "a"
    attr(,"name")
    [1] "DQ299313"
    attr(,"Annot")
    [1] ">DQ299313"
    attr(,"class")
    [1] "SeqFastadna"
    
    $DQ299314
      [1] "a" "a" "a" "t" "g" "a" "c" "a" "t" "t" "a" "a" "g" "g" "a" "a" "a" "g"
     [19] "t" "t" "c" "t" "g" "c" "t" "g" "t" "t" "t" "t" "t" "a" "g" "c" "a" "a"
     [37] "a" "a" "g" "c" "g" "t" "c" "c" "a" "g" "a" "a" "a" "g" "g" "a" "g" "a"
     [55] "g" "c" "t" "t" "a" "g" "c" "a" "g" "g" "a" "g" "t" "c" "c" "t" "a" "g"
     [73] "c" "c" "c" "t" "t" "t" "c" "a" "c" "c" "c" "a" "t" "a" "c" "a" "c" "a"
     [91] "t" "t" "t" "g" "g" "c" "t" "c" "a" "g" "g" "g" "t" "t" "a" "c" "c" "g"
    [109] "a" "a" "g" "a" "g" "g" "g" "g" "c" "c" "a" "a"
    attr(,"name")
    [1] "DQ299314"
    attr(,"Annot")
    [1] ">DQ299314"
    attr(,"class")
    [1] "SeqFastadna"
    
    $DQ299315
      [1] "a" "a" "a" "t" "g" "a" "c" "a" "t" "t" "a" "a" "g" "g" "a" "a" "a" "g"
     [19] "t" "t" "c" "t" "g" "c" "t" "g" "t" "t" "t" "t" "t" "a" "g" "c" "a" "a"
     [37] "a" "a" "g" "c" "g" "t" "c" "c" "a" "g" "a" "a" "a" "g" "g" "a" "g" "a"
     [55] "g" "c" "t" "t" "a" "g" "c" "a" "g" "g" "a" "g" "t" "c" "c" "t" "a" "g"
     [73] "c" "c" "c" "t" "t" "t" "c" "a" "c" "c" "c" "a" "t" "a" "c" "a" "c" "a"
     [91] "t" "t" "t" "g" "g" "c" "t" "c" "a" "g" "g" "g" "t" "t" "a" "c" "c" "g"
    [109] "a" "a" "g" "a" "g" "g" "g" "g" "c" "c" "a" "a"
    attr(,"name")
    [1] "DQ299315"
    attr(,"Annot")
    [1] ">DQ299315"
    attr(,"class")
    [1] "SeqFastadna"
    
    $DQ299316
      [1] "a" "g" "c" "a" "c" "c" "g" "t" "t" "g" "c" "t" "a" "c" "c" "g" "a" "g"
     [19] "t" "g" "t" "c" "t" "g" "t" "c" "t" "a" "a" "g" "a" "a" "c" "a" "c" "a"
     [37] "g" "a" "g" "g" "a" "g" "a" "a" "t" "t" "t" "a" "t" "t" "a" "t" "c" "a"
     [55] "t" "t" "g" "a" "a" "g" "a" "a" "t" "a" "g" "c" "t" "t" "a" "a" "a" "t"
     [73] "g" "a" "c" "t" "g" "c" "g" "g" "t" "a" "a" "c" "c" "a" "g" "g" "t" "a"
     [91] "a" "t" "a" "t" "t" "g" "g" "c" "a" "a" "a"
    attr(,"name")
    [1] "DQ299316"
    attr(,"Annot")
    [1] ">DQ299316"
    attr(,"class")
    [1] "SeqFastadna"
    
    $DQ299317
      [1] "a" "g" "c" "a" "c" "c" "g" "t" "t" "g" "c" "t" "a" "c" "c" "g" "a" "g"
     [19] "t" "g" "t" "c" "t" "g" "t" "c" "t" "a" "a" "g" "a" "a" "c" "a" "c" "a"
     [37] "g" "a" "g" "g" "a" "g" "a" "a" "t" "t" "t" "a" "t" "t" "a" "t" "c" "a"
     [55] "t" "t" "g" "a" "a" "g" "a" "a" "t" "a" "g" "c" "t" "t" "a" "a" "a" "t"
     [73] "g" "a" "c" "t" "g" "c" "g" "g" "t" "a" "a" "c" "c" "a" "g" "g" "t" "a"
     [91] "a" "t" "a" "t" "t" "g" "g" "c" "a" "a" "a"
    attr(,"name")
    [1] "DQ299317"
    attr(,"Annot")
    [1] ">DQ299317"
    attr(,"class")
    [1] "SeqFastadna"
    
    $DQ299318
      [1] "a" "g" "c" "a" "c" "c" "g" "t" "t" "g" "c" "t" "a" "c" "c" "g" "a" "g"
     [19] "t" "g" "t" "c" "t" "g" "t" "c" "t" "a" "a" "g" "a" "a" "c" "a" "c" "a"
     [37] "g" "a" "g" "g" "a" "g" "a" "a" "t" "t" "t" "a" "t" "t" "a" "t" "c" "a"
     [55] "t" "t" "g" "a" "a" "g" "a" "a" "t" "a" "g" "c" "t" "t" "a" "a" "a" "t"
     [73] "g" "a" "c" "t" "g" "c" "g" "g" "t" "a" "a" "c" "c" "a" "g" "g" "t" "a"
     [91] "a" "t" "a" "t" "t" "g" "g" "c" "a" "a" "a"
    attr(,"name")
    [1] "DQ299318"
    attr(,"Annot")
    [1] ">DQ299318"
    attr(,"class")
    [1] "SeqFastadna"
    
    $DQ299319
      [1] "a" "g" "c" "a" "c" "c" "g" "t" "t" "g" "c" "t" "a" "c" "c" "g" "a" "g"
     [19] "t" "g" "t" "c" "t" "g" "t" "c" "t" "a" "a" "g" "a" "a" "c" "a" "c" "a"
     [37] "g" "a" "g" "g" "a" "g" "a" "a" "t" "t" "t" "a" "t" "t" "a" "t" "c" "a"
     [55] "t" "t" "g" "a" "a" "g" "a" "a" "t" "a" "g" "c" "t" "t" "a" "a" "a" "t"
     [73] "g" "a" "c" "t" "g" "c" "g" "g" "t" "a" "a" "c" "c" "a" "g" "g" "t" "a"
     [91] "a" "t" "a" "t" "t" "g" "g" "c" "a" "a" "a"
    attr(,"name")
    [1] "DQ299319"
    attr(,"Annot")
    [1] ">DQ299319"
    attr(,"class")
    [1] "SeqFastadna"
    
    $DQ299320
      [1] "c" "a" "a" "t" "a" "t" "a" "c" "c" "t" "t" "c" "t" "c" "a" "g" "t" "c"
     [19] "t" "a" "c" "t" "a" "g" "g" "c" "a" "t" "a" "g" "c" "a" "c" "c" "g" "t"
     [37] "t" "g" "c" "t" "a" "c" "c" "g" "a" "g" "t" "g" "t" "c" "t" "g" "t" "c"
     [55] "t" "a" "a" "g" "a" "a" "c" "a" "c" "a" "g" "a" "g" "g" "a" "g" "a" "a"
     [73] "t" "t" "t" "a" "t" "t" "a" "t" "c" "a" "t" "t" "g" "a" "a" "g" "a" "a"
     [91] "t" "a" "g" "c" "t" "t" "a" "a" "a" "t" "g" "a" "g" "t" "g" "c" "g" "g"
    [109] "t" "a" "a" "c" "c" "a"
    attr(,"name")
    [1] "DQ299320"
    attr(,"Annot")
    [1] ">DQ299320"
    attr(,"class")
    [1] "SeqFastadna"
    
    $DQ299321
      [1] "t" "t" "a" "a" "a" "g" "t" "g" "t" "c" "t" "a" "a" "t" "a" "a" "t" "g"
     [19] "c" "t" "g" "a" "a" "g" "a" "c" "c" "c" "c" "a" "a" "a" "g" "a" "t" "c"
     [37] "t" "c" "a" "t" "g" "t" "t" "a" "a" "g" "t" "g" "g" "a" "g" "a" "a" "a"
     [55] "g" "g" "g" "t" "t" "t" "t" "g" "c" "a" "a" "a" "c" "t" "g" "a" "a" "a"
     [73] "g" "a" "t" "c" "t" "g" "t" "a" "g" "a" "g" "a" "g" "t" "a" "g" "c" "a"
     [91] "g" "t" "a" "t" "t" "t" "c" "a" "t" "t" "g" "g" "t" "a" "c" "c" "t" "g"
    [109] "g" "t" "a" "c" "t" "g" "a" "t" "t" "a" "t" "g"
    attr(,"name")
    [1] "DQ299321"
    attr(,"Annot")
    [1] ">DQ299321"
    attr(,"class")
    [1] "SeqFastadna"
    
    $DQ299322
      [1] "t" "t" "a" "a" "a" "g" "t" "g" "t" "c" "t" "a" "a" "t" "a" "a" "t" "g"
     [19] "c" "t" "g" "a" "a" "g" "a" "c" "c" "c" "c" "a" "a" "a" "g" "a" "t" "c"
     [37] "t" "c" "a" "t" "g" "t" "t" "a" "a" "g" "t" "g" "g" "a" "g" "a" "a" "a"
     [55] "g" "g" "g" "t" "t" "t" "t" "g" "c" "a" "a" "a" "c" "t" "g" "a" "a" "a"
     [73] "g" "a" "t" "c" "t" "g" "t" "a" "g" "a" "g" "a" "g" "t" "a" "g" "c" "a"
     [91] "g" "t" "a" "t" "t" "t" "c" "a" "t" "t" "g" "g" "t" "a" "c" "c" "t" "g"
    [109] "g" "t" "a" "c" "t" "g" "a" "t" "t" "a" "t" "g"
    attr(,"name")
    [1] "DQ299322"
    attr(,"Annot")
    [1] ">DQ299322"
    attr(,"class")
    [1] "SeqFastadna"
    
    $DQ299323
      [1] "t" "t" "a" "a" "a" "g" "t" "g" "t" "c" "t" "a" "a" "t" "a" "a" "t" "g"
     [19] "c" "t" "g" "a" "a" "g" "a" "c" "c" "c" "c" "a" "a" "a" "g" "a" "t" "c"
     [37] "t" "c" "a" "t" "g" "t" "t" "a" "a" "g" "t" "g" "g" "a" "g" "a" "a" "a"
     [55] "g" "g" "g" "t" "t" "t" "t" "g" "c" "a" "a" "a" "c" "t" "g" "a" "a" "a"
     [73] "g" "a" "t" "c" "t" "g" "t" "a" "g" "a" "g" "a" "g" "t" "a" "g" "c" "a"
     [91] "g" "t" "a" "t" "t" "t" "c" "a" "t" "t" "g" "g" "t" "a" "c" "c" "t" "g"
    [109] "g" "t" "a" "c" "t" "g" "a" "t" "t" "a" "t" "g"
    attr(,"name")
    [1] "DQ299323"
    attr(,"Annot")
    [1] ">DQ299323"
    attr(,"class")
    [1] "SeqFastadna"
    
    $DQ299324
      [1] "t" "t" "a" "a" "a" "g" "t" "g" "t" "c" "t" "a" "a" "t" "a" "a" "t" "g"
     [19] "c" "t" "g" "a" "a" "g" "a" "c" "c" "c" "c" "a" "a" "a" "g" "a" "t" "c"
     [37] "t" "c" "a" "t" "g" "t" "t" "a" "a" "g" "t" "g" "g" "a" "g" "a" "a" "a"
     [55] "g" "g" "g" "t" "t" "t" "t" "g" "c" "a" "a" "a" "c" "t" "g" "a" "a" "a"
     [73] "g" "a" "t" "c" "t" "g" "t" "a" "g" "a" "g" "a" "g" "t" "a" "g" "c" "a"
     [91] "g" "t" "a" "t" "t" "t" "c" "a" "t" "t" "g" "g" "t" "a" "c" "c" "t" "g"
    [109] "g" "t" "a" "c" "t" "g" "a" "t" "t" "a" "t" "g"
    attr(,"name")
    [1] "DQ299324"
    attr(,"Annot")
    [1] ">DQ299324"
    attr(,"class")
    [1] "SeqFastadna"
    
    $DQ299325
      [1] "t" "t" "a" "a" "a" "g" "t" "g" "t" "c" "t" "a" "a" "t" "a" "a" "t" "g"
     [19] "c" "t" "g" "a" "a" "g" "a" "c" "c" "c" "c" "a" "a" "a" "g" "a" "t" "c"
     [37] "t" "c" "a" "t" "g" "t" "t" "a" "a" "g" "t" "g" "g" "a" "g" "a" "a" "a"
     [55] "g" "g" "g" "t" "t" "t" "t" "g" "c" "a" "a" "a" "c" "t" "g" "a" "a" "a"
     [73] "g" "a" "t" "c" "t" "g" "t" "a" "g" "a" "g" "a" "g" "t" "a" "g" "c" "a"
     [91] "g" "t" "a" "t" "t" "t" "c" "a" "t" "t" "g" "g" "t" "a" "c" "c" "t" "g"
    [109] "g" "t" "a" "c" "t" "g" "a" "t" "t" "a" "t" "g"
    attr(,"name")
    [1] "DQ299325"
    attr(,"Annot")
    [1] ">DQ299325"
    attr(,"class")
    [1] "SeqFastadna"
    
    $DQ299326
      [1] "t" "t" "a" "a" "a" "g" "t" "g" "t" "c" "t" "a" "a" "t" "a" "a" "t" "g"
     [19] "c" "t" "g" "a" "a" "g" "a" "c" "c" "c" "c" "a" "a" "a" "g" "a" "t" "c"
     [37] "t" "c" "a" "t" "g" "t" "t" "a" "a" "g" "t" "g" "g" "a" "g" "a" "a" "a"
     [55] "g" "g" "g" "t" "t" "t" "t" "g" "c" "a" "a" "a" "c" "t" "g" "a" "a" "a"
     [73] "g" "a" "t" "c" "t" "g" "t" "a" "g" "a" "g" "a" "g" "t" "a" "g" "c" "a"
     [91] "g" "t" "a" "t" "t" "t" "c" "a" "t" "t" "g" "g" "t" "a" "c" "c" "t" "g"
    [109] "g" "t" "a" "c" "t" "g" "a" "t" "t" "a" "t" "g"
    attr(,"name")
    [1] "DQ299326"
    attr(,"Annot")
    [1] ">DQ299326"
    attr(,"class")
    [1] "SeqFastadna"
    
    $DQ299327
      [1] "t" "t" "a" "a" "a" "g" "t" "g" "t" "c" "t" "a" "a" "t" "a" "a" "t" "g"
     [19] "c" "t" "g" "a" "a" "g" "a" "c" "c" "c" "c" "a" "a" "a" "g" "a" "t" "c"
     [37] "t" "c" "a" "t" "g" "t" "t" "a" "a" "g" "t" "g" "g" "a" "g" "a" "a" "a"
     [55] "g" "g" "g" "t" "t" "t" "t" "g" "c" "a" "a" "a" "c" "t" "g" "a" "a" "a"
     [73] "g" "a" "t" "c" "t" "g" "t" "a" "g" "a" "g" "a" "g" "t" "a" "g" "c" "a"
     [91] "g" "t" "a" "t" "t" "t" "c" "a" "t" "t" "g" "g" "t" "a" "c" "c" "t" "g"
    [109] "g" "t" "a" "c" "t" "g" "a" "t" "t" "a" "t" "g"
    attr(,"name")
    [1] "DQ299327"
    attr(,"Annot")
    [1] ">DQ299327"
    attr(,"class")
    [1] "SeqFastadna"
    
    $DQ299328
      [1] "t" "t" "a" "a" "a" "g" "t" "g" "t" "c" "t" "a" "a" "t" "a" "a" "t" "g"
     [19] "c" "t" "g" "a" "a" "g" "a" "c" "c" "c" "c" "a" "a" "a" "g" "a" "t" "c"
     [37] "t" "c" "a" "t" "g" "t" "t" "a" "a" "g" "t" "g" "g" "a" "g" "a" "a" "a"
     [55] "g" "g" "g" "t" "t" "t" "t" "g" "c" "a" "a" "a" "c" "t" "g" "a" "a" "a"
     [73] "g" "a" "t" "c" "t" "g" "t" "a" "g" "a" "g" "a" "g" "t" "a" "g" "c" "a"
     [91] "g" "t" "a" "t" "t" "t" "c" "a" "t" "t" "g" "g" "t" "a" "c" "c" "t" "g"
    [109] "g" "t" "a" "c" "t" "g" "a" "t" "t" "a" "t" "g"
    attr(,"name")
    [1] "DQ299328"
    attr(,"Annot")
    [1] ">DQ299328"
    attr(,"class")
    [1] "SeqFastadna"
    
    $DQ299330
      [1] "t" "t" "a" "a" "a" "g" "t" "g" "t" "c" "t" "a" "a" "t" "a" "a" "t" "g"
     [19] "c" "t" "g" "a" "a" "g" "a" "c" "c" "c" "c" "a" "a" "a" "g" "a" "t" "c"
     [37] "t" "c" "a" "t" "g" "t" "t" "a" "a" "g" "t" "g" "g" "a" "g" "a" "a" "a"
     [55] "g" "g" "g" "t" "t" "t" "t" "g" "c" "a" "a" "a" "c" "t" "g" "a" "a" "a"
     [73] "g" "a" "t" "c" "t" "g" "t" "a" "g" "a" "g" "a" "g" "t" "a" "g" "c" "a"
     [91] "g" "t" "a" "t" "t" "t" "c" "a" "t" "t" "g" "g" "t" "a" "c" "c" "t" "g"
    [109] "g" "t" "a" "c" "t" "g" "a" "t" "t" "a" "t" "g"
    attr(,"name")
    [1] "DQ299330"
    attr(,"Annot")
    [1] ">DQ299330"
    attr(,"class")
    [1] "SeqFastadna"
    
    $DQ299331
      [1] "t" "t" "a" "a" "a" "g" "t" "g" "t" "c" "t" "a" "a" "t" "a" "a" "t" "g"
     [19] "c" "t" "g" "a" "a" "g" "a" "c" "c" "c" "c" "a" "a" "a" "g" "a" "t" "c"
     [37] "t" "c" "a" "t" "g" "t" "t" "a" "a" "g" "t" "g" "g" "a" "g" "a" "a" "a"
     [55] "g" "g" "g" "t" "t" "t" "t" "g" "c" "a" "a" "a" "c" "t" "g" "a" "a" "a"
     [73] "g" "a" "t" "c" "t" "g" "t" "a" "g" "a" "g" "a" "g" "t" "a" "g" "c" "a"
     [91] "g" "t" "a" "t" "t" "t" "c" "a" "t" "t" "g" "g" "t" "a" "c" "c" "t" "g"
    [109] "g" "t" "a" "c" "t" "g" "a" "t" "t" "a" "t" "g"
    attr(,"name")
    [1] "DQ299331"
    attr(,"Annot")
    [1] ">DQ299331"
    attr(,"class")
    [1] "SeqFastadna"
    
    $DQ363751.BRCA1
       [1] "a" "t" "g" "c" "a" "c" "a" "g" "t" "t" "g" "c" "t" "c" "t" "g" "g" "g"
      [19] "a" "g" "t" "c" "t" "t" "c" "a" "g" "a" "a" "t" "a" "g" "a" "a" "a" "c"
      [37] "t" "a" "c" "c" "c" "a" "t" "c" "t" "c" "a" "a" "g" "a" "g" "g" "a" "g"
      [55] "c" "t" "c" "a" "t" "t" "a" "a" "g" "g" "t" "t" "g" "t" "t" "g" "a" "t"
      [73] "g" "t" "g" "g" "a" "g" "g" "a" "g" "c" "a" "a" "c" "a" "g" "c" "t" "g"
      [91] "g" "a" "a" "g" "a" "g" "t" "c" "t" "g" "g" "g" "c" "c" "a" "c" "a" "c"
     [109] "g" "a" "t" "t" "t" "g" "a" "c" "g" "g" "a" "a" "a" "c" "a" "t" "c" "t"
     [127] "t" "a" "c" "t" "t" "g" "c" "c" "a" "a" "g" "g" "c" "a" "a" "g" "a" "t"
     [145] "c" "t" "a" "g" "a" "g" "g" "g" "a" "a" "c" "c" "c" "c" "t" "t" "a" "c"
     [163] "c" "t" "g" "g" "a" "a" "t" "c" "t" "g" "g" "a" "a" "t" "c" "a" "g" "c"
     [181] "c" "t" "c" "t" "t" "c" "t" "c" "t" "g" "a" "t" "g" "a" "c" "c" "c" "t"
     [199] "g" "a" "a" "t" "c" "t" "g" "a" "t" "c" "c" "t" "t" "c" "t" "g" "a" "a"
     [217] "g" "a" "c" "a" "g" "a" "g" "c" "c" "c" "c" "a" "g" "a" "g" "t" "c" "a"
     [235] "g" "c" "t" "c" "g" "t" "g" "t" "t" "g" "g" "c" "a" "a" "c" "a" "t" "a"
     [253] "c" "c" "a" "t" "c" "t" "t" "c" "a" "a" "c" "c" "t" "c" "t" "g" "c" "a"
     [271] "t" "t" "g" "a" "a" "g" "g" "t" "t" "c" "c" "c" "c" "a" "a" "t" "t" "g"
     [289] "a" "a" "a" "g" "t" "t" "g" "c" "a" "g" "a" "a" "t" "c" "t" "g" "c" "c"
     [307] "c" "a" "g" "g" "g" "t" "c" "c" "a" "g" "c" "t" "g" "c" "t" "g" "c" "t"
     [325] "c" "a" "c" "a" "c" "t" "a" "c" "t" "g" "a" "t" "a" "c" "t" "g" "c" "t"
     [343] "g" "g" "g" "t" "a" "t" "a" "a" "t" "g" "c" "a" "a" "t" "g" "g" "a" "a"
     [361] "g" "a" "a" "a" "g" "t" "g" "t" "g" "a" "g" "c" "a" "g" "g" "g" "a" "g"
     [379] "a" "a" "g" "c" "c" "a" "g" "a" "a" "t" "t" "g" "a" "c" "a" "g" "c" "t"
     [397] "t" "c" "a" "a" "c" "a" "g" "a" "a" "a" "g" "g" "g" "t" "c" "a" "a" "c"
     [415] "a" "a" "a" "a" "g" "a" "a" "t" "g" "t" "c" "c" "a" "t" "g" "g" "t" "g"
     [433] "g" "t" "g" "t" "c" "t" "g" "g" "c" "c" "t" "g" "a" "c" "c" "c" "c" "a"
     [451] "g" "a" "a" "g" "a" "a" "t" "t" "t" "a" "t" "g" "c" "t" "c" "g" "t" "g"
     [469] "t" "a" "c" "a" "a" "g" "t" "t" "t" "g" "c" "c" "a" "g" "a" "a" "a" "a"
     [487] "c" "a" "c" "c" "a" "c" "a" "t" "c" "a" "c" "t" "t" "t" "a" "a" "c" "t"
     [505] "a" "a" "t" "c" "t" "a" "a" "t" "t" "a" "c" "t" "g" "a" "a" "g" "a" "g"
     [523] "a" "c" "t" "a" "c" "t" "c" "a" "t" "g" "t" "t" "g" "t" "t" "a" "t" "g"
     [541] "a" "a" "a" "a" "c" "a" "g" "a" "t" "g" "c" "t" "g" "a" "g" "t" "t" "t"
     [559] "g" "t" "g" "t" "g" "t" "g" "a" "a" "c" "g" "g" "a" "c" "a" "c" "t" "g"
     [577] "a" "a" "a" "t" "a" "t" "t" "t" "t" "c" "t" "a" "g" "g" "a" "a" "t" "t"
     [595] "g" "c" "g" "g" "g" "a" "g" "g" "a" "a" "a" "a" "t" "g" "g" "g" "t" "a"
     [613] "g" "t" "t" "a" "g" "c" "t" "a" "t" "t" "t" "c" "t" "g" "g" "g" "t" "g"
     [631] "a" "c" "c" "c" "a" "g" "t" "c" "t" "a" "t" "t" "a" "a" "a" "g" "a" "a"
     [649] "a" "g" "a" "a" "a" "a" "a" "t" "g" "c" "t" "g" "a" "a" "t" "g" "a" "g"
     [667] "c" "a" "t" "g" "a" "t" "t" "t" "t" "g" "a" "a" "g" "t" "c" "a" "g" "a"
     [685] "g" "g" "a" "g" "a" "t" "g" "t" "g" "g" "t" "c" "a" "a" "t" "g" "g" "a"
     [703] "a" "g" "a" "a" "a" "c" "c" "a" "c" "c" "a" "a" "g" "g" "t" "c" "c" "a"
     [721] "a" "a" "g" "c" "g" "a" "g" "c" "a" "a" "g" "a" "g" "a" "a" "t" "c" "c"
     [739] "c" "a" "g" "g" "a" "c" "a" "g" "a" "a" "a" "g" "a" "t" "c" "t" "t" "c"
     [757] "a" "g" "g" "g" "g" "g" "c" "t" "a" "g" "a" "a" "a" "t" "c" "t" "g" "t"
     [775] "t" "g" "c" "t" "a" "t" "g" "g" "g" "c" "c" "c" "t" "t" "c" "a" "c" "c"
     [793] "a" "a" "c" "a" "t" "g" "c" "c" "c" "a" "c" "a" "g" "a" "t" "c" "a" "a"
     [811] "c" "t" "g" "g" "a" "a" "t" "g" "g" "a" "t" "g" "g" "t" "a" "c" "a" "g"
     [829] "c" "t" "g" "t" "g" "t" "g" "g" "t" "g" "c" "t" "t" "c" "t" "g" "t" "g"
     [847] "g" "t" "g" "a" "a" "g" "g" "a" "g" "c" "t" "t" "t" "c" "a" "t" "c" "a"
     [865] "t" "t" "c" "a" "c" "c" "c" "t" "t" "g" "g" "c" "a" "c" "a" "g" "g" "t"
     [883] "g" "t" "c" "c" "a" "c" "c" "c" "a" "a" "t" "t" "g" "t" "g" "g" "t" "t"
     [901] "g" "t" "g" "c" "a" "g" "c" "c" "a" "g" "a" "t" "g" "c" "c" "t" "g" "g"
     [919] "a" "c" "a" "g" "a" "g" "g" "a" "c" "a" "a" "t" "g" "g" "c" "t" "t" "c"
     [937] "c" "a" "t" "g" "c" "a" "a" "t" "t" "g" "g" "g" "c" "a" "g" "a" "t" "g"
     [955] "t" "g" "t" "g" "a" "g" "a" "c" "a" "c" "c" "t" "g" "t" "g" "g" "t" "g"
     [973] "a" "c" "c" "c" "g" "a" "g" "a" "g" "t" "g" "g" "g" "t" "g" "t" "t" "g"
     [991] "g" "a" "c" "a" "g" "t" "g" "t" "a" "g" "c" "a" "c" "t" "c" "t" "a" "c"
    [1009] "c" "a" "g" "t" "g" "c" "c" "a" "g" "g" "a" "g" "c" "t" "g" "g" "a" "c"
    [1027] "a" "c" "c" "t" "a" "c" "c" "t" "g" "a" "t" "a" "c" "c" "c" "c" "a" "g"
    [1045] "a" "t" "c" "c" "c" "c" "c" "a" "c" "a" "g" "t" "c" "a" "c" "t" "a" "c"
    [1063] "t" "g" "a"
    attr(,"name")
    [1] "DQ363751.BRCA1"
    attr(,"Annot")
    [1] ">DQ363751.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $DQ478408.BRCA1
       [1] "a" "t" "g" "g" "a" "t" "t" "t" "a" "t" "c" "t" "g" "c" "t" "c" "t" "t"
      [19] "c" "g" "c" "g" "t" "t" "g" "a" "a" "g" "a" "a" "g" "t" "a" "c" "a" "a"
      [37] "a" "a" "t" "g" "t" "c" "a" "t" "t" "a" "a" "t" "g" "c" "t" "a" "t" "g"
      [55] "c" "a" "g" "a" "a" "a" "a" "t" "c" "t" "t" "a" "g" "a" "g" "t" "g" "t"
      [73] "c" "c" "c" "a" "t" "c" "t" "g" "t" "c" "t" "g" "g" "a" "g" "t" "t" "g"
      [91] "a" "t" "c" "a" "a" "g" "g" "a" "a" "c" "c" "t" "g" "t" "c" "t" "c" "c"
     [109] "a" "c" "a" "a" "a" "g" "t" "g" "t" "g" "a" "c" "c" "a" "c" "a" "t" "a"
     [127] "t" "t" "t" "t" "g" "c" "a" "a" "a" "t" "t" "t" "t" "g" "c" "a" "t" "g"
     [145] "c" "t" "g" "a" "a" "a" "c" "t" "t" "c" "t" "c" "a" "a" "c" "c" "a" "g"
     [163] "a" "a" "g" "a" "a" "a" "g" "g" "g" "c" "c" "t" "t" "c" "a" "c" "a" "g"
     [181] "t" "g" "t" "c" "c" "t" "t" "t" "a" "t" "g" "t" "a" "a" "g" "a" "a" "t"
     [199] "g" "a" "t" "a" "t" "a" "a" "c" "c" "a" "a" "a" "a" "g" "g" "a" "g" "c"
     [217] "c" "t" "a" "c" "a" "a" "g" "a" "a" "a" "g" "t" "a" "c" "g" "a" "g" "a"
     [235] "t" "t" "t" "a" "g" "t" "c" "a" "a" "c" "t" "t" "g" "t" "t" "g" "a" "a"
     [253] "g" "a" "g" "c" "t" "a" "t" "t" "g" "a" "a" "a" "a" "t" "c" "a" "t" "t"
     [271] "t" "g" "t" "g" "c" "t" "t" "t" "t" "c" "a" "g" "c" "t" "t" "g" "a" "c"
     [289] "a" "c" "a" "g" "g" "t" "t" "t" "g" "g" "a" "g" "t" "a" "t" "g" "c" "a"
     [307] "a" "a" "c" "a" "g" "c" "t" "a" "t" "a" "a" "t" "t" "t" "t" "g" "c" "a"
     [325] "a" "a" "a" "a" "a" "g" "g" "a" "a" "a" "a" "t" "a" "a" "c" "t" "c" "t"
     [343] "c" "c" "t" "g" "a" "a" "c" "a" "t" "c" "t" "a" "a" "a" "a" "g" "a" "t"
     [361] "g" "a" "a" "g" "t" "t" "t" "c" "t" "a" "t" "c" "a" "t" "c" "c" "a" "a"
     [379] "a" "g" "t" "a" "t" "g" "g" "g" "c" "t" "a" "c" "a" "g" "a" "a" "a" "c"
     [397] "c" "g" "t" "g" "c" "c" "a" "a" "a" "a" "g" "a" "c" "t" "t" "c" "t" "a"
     [415] "c" "a" "g" "a" "g" "t" "g" "a" "a" "c" "c" "c" "g" "a" "a" "a" "a" "t"
     [433] "c" "c" "t" "t" "c" "c" "t" "t" "g" "c" "a" "g" "g" "a" "a" "a" "c" "c"
     [451] "a" "g" "t" "c" "t" "c" "a" "g" "t" "g" "t" "c" "c" "a" "a" "c" "t" "c"
     [469] "t" "c" "t" "a" "a" "c" "c" "t" "t" "g" "g" "a" "a" "c" "t" "g" "t" "g"
     [487] "a" "g" "a" "a" "c" "t" "c" "t" "g" "a" "g" "g" "a" "c" "a" "a" "a" "g"
     [505] "c" "a" "g" "c" "g" "g" "a" "t" "a" "c" "a" "a" "c" "c" "t" "c" "a" "a"
     [523] "a" "a" "g" "a" "c" "g" "t" "c" "t" "g" "t" "c" "t" "a" "c" "a" "t" "t"
     [541] "g" "a" "a" "t" "t" "g" "g" "g" "a" "t" "c" "t" "g" "a" "t" "t" "c" "t"
     [559] "t" "c" "t" "g" "a" "a" "g" "a" "t" "a" "c" "c" "g" "t" "t" "a" "a" "t"
     [577] "a" "a" "g" "g" "c" "a" "a" "c" "t" "t" "a" "t" "t" "g" "c" "a" "g" "t"
     [595] "g" "t" "g" "g" "g" "a" "g" "a" "t" "c" "a" "a" "g" "a" "a" "t" "t" "g"
     [613] "t" "t" "a" "c" "a" "a" "a" "t" "c" "a" "c" "c" "c" "c" "t" "c" "a" "a"
     [631] "g" "g" "a" "a" "c" "c" "a" "g" "g" "g" "a" "t" "g" "a" "a" "a" "t" "c"
     [649] "a" "g" "t" "t" "t" "g" "g" "a" "t" "t" "c" "t" "g" "c" "a" "a" "a" "a"
     [667] "a" "a" "g" "g" "c" "t" "g" "c" "t" "t" "g" "t" "g" "a" "a" "t" "t" "t"
     [685] "t" "c" "t" "g" "a" "g" "a" "c" "g" "g" "a" "t" "g" "t" "a" "a" "c" "a"
     [703] "a" "a" "t" "a" "c" "t" "g" "a" "a" "c" "a" "t" "c" "a" "t" "c" "a" "a"
     [721] "c" "c" "c" "a" "g" "t" "a" "a" "t" "a" "a" "t" "g" "a" "t" "t" "t" "g"
     [739] "a" "a" "c" "a" "c" "c" "a" "c" "t" "g" "a" "g" "a" "a" "g" "c" "g" "t"
     [757] "g" "c" "a" "g" "c" "t" "g" "a" "g" "a" "g" "g" "c" "a" "t" "c" "c" "a"
     [775] "g" "a" "a" "a" "a" "g" "t" "a" "t" "c" "a" "g" "g" "g" "t" "a" "g" "t"
     [793] "t" "c" "t" "g" "t" "t" "t" "c" "a" "a" "a" "c" "t" "t" "g" "c" "a" "t"
     [811] "g" "t" "g" "g" "a" "g" "c" "c" "a" "t" "g" "t" "g" "g" "c" "a" "c" "a"
     [829] "a" "a" "t" "a" "c" "t" "c" "a" "t" "g" "c" "c" "a" "g" "c" "t" "c" "a"
     [847] "t" "t" "a" "c" "a" "g" "c" "a" "t" "g" "a" "g" "a" "a" "c" "a" "g" "c"
     [865] "a" "g" "t" "t" "t" "a" "t" "t" "a" "c" "t" "c" "a" "c" "t" "a" "a" "a"
     [883] "g" "a" "c" "a" "g" "a" "a" "t" "g" "a" "a" "t" "g" "t" "a" "g" "a" "a"
     [901] "a" "a" "g" "g" "c" "t" "g" "a" "a" "t" "t" "c" "t" "g" "t" "a" "a" "t"
     [919] "a" "a" "a" "a" "g" "c" "a" "a" "a" "c" "a" "g" "c" "c" "t" "g" "g" "c"
     [937] "t" "t" "a" "g" "c" "a" "a" "g" "g" "a" "g" "c" "c" "a" "a" "c" "a" "t"
     [955] "a" "a" "c" "a" "g" "a" "t" "g" "g" "g" "c" "t" "g" "g" "a" "a" "g" "t"
     [973] "a" "a" "g" "g" "a" "a" "a" "c" "a" "t" "g" "t" "a" "a" "t" "g" "a" "t"
     [991] "a" "g" "g" "c" "g" "g" "a" "c" "t" "c" "c" "c" "a" "g" "c" "a" "c" "a"
    [1009] "g" "a" "a" "a" "a" "a" "a" "a" "g" "g" "t" "a" "g" "a" "t" "c" "t" "g"
    [1027] "a" "a" "t" "g" "c" "t" "g" "a" "t" "c" "c" "c" "c" "t" "g" "t" "g" "t"
    [1045] "g" "a" "g" "a" "g" "a" "a" "a" "a" "g" "a" "a" "t" "g" "g" "a" "a" "t"
    [1063] "a" "a" "g" "c" "a" "g" "a" "a" "a" "c" "t" "g" "c" "c" "a" "t" "g" "c"
    [1081] "t" "c" "a" "g" "a" "g" "a" "a" "t" "c" "c" "t" "a" "g" "a" "g" "a" "t"
    [1099] "a" "c" "t" "g" "a" "a" "g" "a" "t" "g" "t" "t" "c" "c" "t" "t" "g" "g"
    [1117] "a" "t" "a" "a" "c" "a" "c" "t" "a" "a" "a" "t" "a" "g" "c" "a" "g" "c"
    [1135] "a" "t" "t" "c" "a" "g" "a" "a" "a" "g" "t" "t" "a" "a" "t" "g" "a" "g"
    [1153] "t" "g" "g" "t" "t" "t" "t" "c" "c" "a" "g" "a" "a" "g" "t" "g" "a" "t"
    [1171] "g" "a" "a" "c" "t" "g" "t" "t" "a" "g" "g" "t" "t" "c" "t" "g" "a" "t"
    [1189] "g" "a" "c" "t" "c" "a" "c" "a" "t" "g" "a" "t" "g" "g" "g" "g" "a" "g"
    [1207] "t" "c" "t" "g" "a" "a" "t" "c" "a" "a" "a" "t" "g" "c" "c" "a" "a" "a"
    [1225] "g" "t" "a" "g" "c" "t" "g" "a" "t" "g" "t" "a" "t" "t" "g" "g" "a" "c"
    [1243] "g" "t" "t" "c" "t" "a" "a" "a" "t" "g" "a" "g" "g" "t" "a" "g" "a" "t"
    [1261] "g" "a" "a" "t" "a" "t" "t" "c" "t" "g" "g" "t" "t" "c" "t" "t" "c" "a"
    [1279] "g" "a" "g" "a" "a" "a" "a" "t" "a" "g" "a" "c" "t" "t" "a" "c" "t" "g"
    [1297] "g" "c" "c" "a" "g" "t" "g" "a" "t" "c" "c" "t" "c" "a" "t" "g" "a" "g"
    [1315] "g" "c" "t" "t" "t" "a" "a" "t" "a" "t" "g" "t" "a" "a" "a" "a" "g" "t"
    [1333] "g" "a" "a" "a" "g" "a" "g" "t" "t" "c" "a" "c" "t" "c" "c" "a" "a" "a"
    [1351] "t" "c" "a" "g" "t" "a" "g" "a" "g" "a" "g" "t" "a" "a" "t" "a" "t" "t"
    [1369] "g" "a" "a" "g" "a" "c" "a" "a" "a" "a" "t" "a" "t" "t" "t" "g" "g" "g"
    [1387] "a" "a" "a" "a" "c" "c" "t" "a" "t" "c" "g" "g" "a" "a" "g" "a" "a" "g"
    [1405] "g" "c" "a" "a" "g" "c" "c" "t" "c" "c" "c" "c" "a" "a" "c" "t" "t" "a"
    [1423] "a" "g" "c" "c" "a" "t" "g" "t" "a" "a" "c" "t" "g" "a" "a" "a" "a" "t"
    [1441] "c" "t" "a" "a" "t" "t" "a" "t" "a" "g" "g" "a" "g" "c" "a" "t" "t" "t"
    [1459] "g" "t" "t" "a" "c" "t" "g" "a" "g" "c" "c" "a" "c" "a" "g" "a" "t" "a"
    [1477] "a" "t" "a" "c" "a" "a" "g" "a" "g" "c" "g" "t" "c" "c" "c" "c" "t" "c"
    [1495] "a" "c" "a" "a" "a" "t" "a" "a" "a" "t" "t" "a" "a" "a" "g" "c" "g" "t"
    [1513] "a" "a" "a" "a" "g" "g" "a" "g" "a" "c" "c" "t" "a" "c" "a" "t" "c" "a"
    [1531] "g" "g" "c" "c" "t" "t" "c" "a" "t" "c" "c" "t" "g" "a" "g" "g" "a" "t"
    [1549] "t" "t" "t" "a" "t" "c" "a" "a" "g" "a" "a" "a" "g" "c" "a" "g" "a" "t"
    [1567] "t" "t" "g" "g" "c" "a" "g" "t" "t" "c" "a" "a" "a" "a" "g" "a" "c" "t"
    [1585] "c" "c" "t" "g" "a" "a" "a" "t" "g" "a" "t" "a" "a" "a" "t" "c" "a" "g"
    [1603] "g" "g" "a" "a" "c" "t" "a" "a" "c" "c" "a" "a" "a" "c" "g" "g" "a" "g"
    [1621] "c" "a" "g" "a" "a" "t" "g" "g" "t" "c" "a" "a" "g" "t" "g" "a" "t" "g"
    [1639] "a" "a" "t" "a" "t" "t" "a" "c" "t" "a" "a" "t" "a" "g" "t" "g" "g" "t"
    [1657] "c" "a" "t" "g" "a" "g" "a" "a" "t" "a" "a" "a" "a" "c" "a" "a" "a" "a"
    [1675] "g" "g" "t" "g" "a" "t" "t" "c" "t" "a" "t" "t" "c" "a" "g" "a" "a" "t"
    [1693] "g" "a" "g" "a" "a" "a" "a" "a" "t" "c" "c" "t" "a" "a" "c" "c" "c" "a"
    [1711] "a" "t" "a" "g" "a" "a" "t" "c" "a" "c" "t" "c" "g" "a" "a" "a" "a" "a"
    [1729] "g" "a" "a" "t" "c" "t" "g" "c" "t" "t" "t" "c" "a" "a" "a" "a" "c" "g"
    [1747] "a" "a" "a" "g" "c" "t" "g" "a" "a" "c" "c" "t" "a" "t" "a" "a" "g" "c"
    [1765] "a" "g" "c" "a" "g" "t" "a" "t" "a" "a" "g" "c" "a" "a" "t" "a" "t" "g"
    [1783] "g" "a" "a" "c" "t" "c" "g" "a" "a" "t" "t" "a" "a" "a" "t" "a" "t" "c"
    [1801] "c" "a" "c" "a" "a" "t" "t" "c" "a" "a" "a" "a" "g" "c" "a" "c" "c" "t"
    [1819] "a" "a" "a" "a" "a" "g" "a" "a" "t" "a" "g" "g" "c" "t" "g" "a" "g" "g"
    [1837] "a" "g" "g" "a" "a" "g" "t" "c" "t" "t" "c" "t" "a" "c" "c" "a" "g" "g"
    [1855] "c" "a" "t" "a" "t" "t" "c" "a" "t" "g" "c" "g" "c" "t" "t" "g" "a" "a"
    [1873] "c" "t" "a" "g" "t" "a" "g" "t" "c" "a" "g" "t" "a" "g" "a" "a" "a" "t"
    [1891] "c" "t" "a" "a" "g" "c" "c" "c" "a" "c" "c" "t" "a" "a" "t" "t" "g" "t"
    [1909] "a" "c" "t" "g" "a" "a" "t" "t" "g" "c" "a" "a" "a" "t" "t" "g" "a" "t"
    [1927] "a" "g" "t" "t" "g" "t" "t" "c" "t" "a" "g" "c" "a" "g" "t" "g" "a" "a"
    [1945] "g" "a" "g" "a" "t" "a" "a" "a" "g" "a" "a" "a" "a" "a" "a" "a" "a" "g"
    [1963] "t" "a" "c" "a" "a" "c" "c" "a" "a" "a" "t" "g" "c" "c" "a" "g" "t" "c"
    [1981] "a" "g" "g" "c" "a" "c" "a" "g" "c" "a" "g" "a" "a" "a" "c" "c" "t" "a"
    [1999] "c" "a" "a" "c" "t" "c" "a" "t" "g" "g" "a" "a" "g" "g" "t" "a" "a" "a"
    [2017] "g" "a" "a" "c" "c" "t" "g" "c" "a" "a" "c" "t" "g" "g" "a" "g" "c" "c"
    [2035] "a" "a" "g" "a" "a" "g" "a" "g" "t" "a" "a" "c" "a" "a" "g" "c" "c" "a"
    [2053] "a" "a" "t" "g" "a" "a" "c" "a" "g" "a" "c" "a" "a" "g" "t" "a" "a" "a"
    [2071] "a" "g" "a" "c" "a" "t" "g" "a" "c" "a" "g" "c" "g" "a" "t" "a" "c" "t"
    [2089] "t" "t" "c" "c" "c" "a" "g" "a" "g" "c" "t" "g" "a" "a" "g" "t" "t" "a"
    [2107] "a" "c" "a" "a" "a" "t" "g" "c" "a" "c" "c" "t" "g" "g" "t" "t" "c" "t"
    [2125] "t" "t" "t" "a" "c" "t" "a" "a" "g" "t" "g" "t" "t" "c" "a" "a" "a" "t"
    [2143] "a" "c" "c" "a" "g" "t" "g" "a" "a" "c" "t" "t" "a" "a" "a" "g" "a" "a"
    [2161] "t" "t" "t" "g" "t" "c" "a" "a" "t" "c" "c" "t" "a" "g" "c" "c" "t" "t"
    [2179] "c" "c" "a" "a" "g" "a" "g" "a" "a" "g" "a" "a" "a" "a" "a" "g" "a" "a"
    [2197] "g" "a" "g" "a" "a" "a" "c" "t" "a" "g" "a" "a" "a" "c" "a" "g" "t" "t"
    [2215] "a" "a" "a" "g" "t" "g" "t" "c" "t" "a" "a" "t" "a" "a" "t" "g" "c" "t"
    [2233] "g" "a" "a" "g" "a" "c" "c" "c" "c" "a" "a" "a" "g" "a" "t" "c" "t" "c"
    [2251] "a" "t" "g" "t" "t" "a" "a" "g" "t" "g" "g" "a" "g" "a" "a" "a" "g" "g"
    [2269] "g" "t" "t" "t" "t" "g" "c" "a" "a" "a" "c" "t" "g" "a" "a" "a" "g" "a"
    [2287] "t" "c" "t" "g" "t" "a" "g" "a" "g" "a" "g" "t" "a" "g" "c" "a" "g" "t"
    [2305] "a" "t" "t" "t" "c" "a" "t" "t" "g" "g" "t" "a" "c" "c" "t" "g" "g" "t"
    [2323] "a" "c" "t" "g" "a" "t" "t" "a" "t" "g" "g" "c" "a" "c" "t" "c" "a" "g"
    [2341] "g" "a" "a" "a" "g" "t" "a" "t" "c" "t" "c" "g" "t" "t" "a" "c" "t" "g"
    [2359] "g" "a" "a" "g" "t" "t" "a" "g" "c" "a" "c" "t" "c" "t" "a" "g" "g" "g"
    [2377] "a" "a" "g" "g" "c" "a" "a" "a" "a" "a" "c" "a" "g" "a" "a" "c" "c" "a"
    [2395] "a" "a" "t" "a" "a" "a" "t" "g" "t" "g" "t" "g" "a" "g" "t" "c" "a" "g"
    [2413] "t" "g" "t" "g" "c" "a" "g" "c" "a" "t" "t" "t" "g" "a" "a" "a" "a" "c"
    [2431] "c" "c" "c" "a" "a" "g" "g" "g" "a" "c" "t" "a" "a" "t" "t" "c" "a" "t"
    [2449] "g" "g" "t" "t" "g" "t" "t" "c" "c" "a" "a" "a" "g" "a" "t" "a" "a" "t"
    [2467] "a" "g" "a" "a" "a" "t" "g" "a" "c" "a" "c" "a" "g" "a" "a" "g" "g" "c"
    [2485] "t" "t" "t" "a" "a" "g" "t" "a" "t" "c" "c" "a" "t" "t" "g" "g" "g" "a"
    [2503] "c" "a" "t" "g" "a" "a" "g" "t" "t" "a" "a" "c" "c" "a" "c" "a" "g" "t"
    [2521] "c" "g" "g" "g" "a" "a" "a" "c" "a" "a" "g" "c" "a" "t" "a" "g" "a" "a"
    [2539] "a" "t" "g" "g" "a" "a" "g" "a" "a" "a" "g" "t" "g" "a" "a" "c" "t" "t"
    [2557] "g" "a" "t" "g" "c" "t" "c" "a" "g" "t" "a" "t" "t" "t" "g" "c" "a" "g"
    [2575] "a" "a" "t" "a" "c" "a" "t" "t" "c" "a" "a" "g" "g" "t" "t" "t" "c" "a"
    [2593] "a" "a" "g" "c" "g" "c" "c" "a" "g" "t" "c" "a" "t" "t" "t" "g" "c" "t"
    [2611] "c" "c" "g" "t" "t" "t" "t" "c" "a" "a" "a" "t" "c" "c" "a" "g" "g" "a"
    [2629] "a" "a" "t" "g" "c" "a" "g" "a" "a" "g" "a" "g" "g" "a" "a" "t" "g" "t"
    [2647] "g" "c" "a" "a" "c" "a" "t" "t" "c" "t" "c" "t" "g" "c" "c" "c" "a" "c"
    [2665] "t" "c" "t" "g" "g" "g" "t" "c" "c" "t" "t" "a" "a" "a" "g" "a" "a" "a"
    [2683] "c" "a" "a" "a" "g" "t" "c" "c" "a" "a" "a" "a" "g" "t" "c" "a" "c" "t"
    [2701] "t" "t" "t" "g" "a" "a" "t" "g" "t" "g" "a" "a" "c" "a" "a" "a" "a" "g"
    [2719] "g" "a" "a" "g" "a" "a" "a" "a" "t" "c" "a" "a" "g" "g" "a" "a" "a" "g"
    [2737] "a" "a" "t" "g" "a" "g" "t" "c" "t" "a" "a" "t" "a" "t" "c" "a" "a" "g"
    [2755] "c" "c" "t" "g" "t" "a" "c" "a" "g" "a" "c" "a" "g" "t" "t" "a" "a" "t"
    [2773] "a" "t" "c" "a" "c" "t" "g" "c" "a" "g" "g" "c" "t" "t" "t" "c" "c" "t"
    [2791] "g" "t" "g" "g" "t" "t" "g" "g" "t" "c" "a" "g" "a" "a" "a" "g" "a" "t"
    [2809] "a" "a" "g" "c" "c" "a" "g" "t" "t" "g" "a" "t" "a" "a" "t" "g" "c" "c"
    [2827] "a" "a" "a" "t" "g" "t" "a" "g" "t" "a" "t" "c" "a" "a" "a" "g" "g" "a"
    [2845] "g" "g" "c" "t" "c" "t" "a" "g" "g" "t" "t" "t" "t" "g" "t" "c" "t" "a"
    [2863] "t" "c" "a" "t" "c" "t" "c" "a" "g" "t" "t" "c" "a" "g" "a" "g" "g" "c"
    [2881] "a" "a" "c" "g" "a" "a" "a" "c" "t" "g" "g" "a" "c" "t" "c" "a" "t" "t"
    [2899] "a" "c" "t" "c" "c" "a" "a" "a" "t" "a" "a" "a" "c" "a" "t" "g" "g" "a"
    [2917] "c" "t" "t" "t" "t" "a" "c" "a" "a" "a" "a" "c" "c" "c" "a" "t" "a" "t"
    [2935] "c" "g" "t" "a" "t" "a" "c" "c" "a" "c" "c" "a" "c" "t" "t" "t" "t" "t"
    [2953] "c" "c" "c" "a" "t" "c" "a" "a" "g" "t" "c" "a" "t" "t" "t" "g" "t" "t"
    [2971] "a" "a" "a" "a" "c" "t" "a" "a" "a" "t" "g" "t" "a" "a" "g" "a" "a" "a"
    [2989] "a" "a" "t" "c" "t" "g" "c" "t" "a" "g" "a" "g" "g" "a" "a" "a" "a" "c"
    [3007] "t" "t" "t" "g" "a" "g" "g" "a" "a" "c" "a" "t" "t" "c" "a" "a" "t" "g"
    [3025] "t" "c" "a" "c" "c" "t" "g" "a" "a" "a" "g" "a" "g" "a" "a" "a" "t" "g"
    [3043] "g" "g" "a" "a" "a" "t" "g" "a" "g" "a" "a" "c" "a" "t" "t" "c" "c" "a"
    [3061] "a" "g" "t" "a" "c" "a" "g" "t" "g" "a" "g" "c" "a" "c" "a" "a" "t" "t"
    [3079] "a" "g" "c" "c" "g" "t" "a" "a" "t" "a" "a" "c" "a" "t" "t" "a" "g" "a"
    [3097] "g" "a" "a" "a" "a" "t" "g" "t" "t" "t" "t" "t" "a" "a" "a" "g" "a" "a"
    [3115] "g" "c" "c" "a" "g" "c" "t" "c" "a" "a" "g" "c" "a" "a" "t" "a" "t" "t"
    [3133] "a" "a" "t" "g" "a" "a" "g" "t" "a" "g" "g" "t" "t" "c" "c" "a" "g" "t"
    [3151] "a" "c" "t" "a" "a" "t" "g" "a" "a" "g" "t" "g" "g" "g" "c" "t" "c" "c"
    [3169] "a" "g" "t" "a" "t" "t" "a" "a" "t" "g" "a" "a" "a" "t" "a" "g" "g" "t"
    [3187] "t" "c" "c" "a" "g" "t" "g" "a" "t" "g" "a" "a" "a" "a" "c" "a" "t" "t"
    [3205] "c" "a" "a" "g" "c" "a" "g" "a" "a" "c" "t" "a" "g" "g" "t" "a" "g" "a"
    [3223] "a" "a" "c" "a" "g" "a" "g" "g" "g" "c" "c" "a" "a" "a" "a" "t" "t" "g"
    [3241] "a" "a" "t" "g" "c" "t" "a" "t" "g" "c" "t" "t" "a" "g" "a" "t" "t" "a"
    [3259] "g" "g" "g" "g" "t" "t" "t" "t" "g" "c" "a" "a" "c" "c" "t" "g" "a" "g"
    [3277] "g" "t" "c" "t" "a" "t" "a" "a" "a" "c" "a" "a" "a" "g" "t" "c" "t" "t"
    [3295] "c" "c" "t" "g" "g" "a" "a" "g" "t" "a" "a" "t" "t" "g" "t" "a" "a" "g"
    [3313] "c" "a" "t" "c" "c" "t" "g" "a" "a" "a" "t" "a" "a" "a" "a" "a" "a" "g"
    [3331] "c" "a" "a" "g" "a" "a" "t" "a" "t" "g" "a" "a" "g" "a" "a" "g" "t" "a"
    [3349] "g" "t" "t" "c" "a" "g" "a" "c" "t" "g" "t" "t" "a" "a" "t" "a" "c" "a"
    [3367] "g" "a" "t" "t" "t" "c" "t" "c" "t" "c" "c" "a" "t" "a" "t" "c" "t" "g"
    [3385] "a" "t" "t" "t" "c" "a" "g" "a" "t" "a" "a" "c" "t" "t" "a" "g" "a" "a"
    [3403] "c" "a" "g" "c" "c" "t" "a" "t" "g" "g" "g" "a" "a" "g" "t" "a" "g" "t"
    [3421] "c" "a" "t" "g" "c" "a" "t" "c" "t" "c" "a" "g" "g" "t" "t" "t" "g" "t"
    [3439] "t" "c" "t" "g" "a" "g" "a" "c" "a" "c" "c" "t" "g" "a" "t" "g" "a" "c"
    [3457] "c" "t" "g" "t" "t" "a" "g" "a" "t" "g" "a" "t" "g" "g" "t" "g" "a" "a"
    [3475] "a" "t" "a" "a" "a" "g" "g" "a" "a" "g" "a" "t" "a" "c" "t" "a" "g" "t"
    [3493] "t" "t" "t" "g" "c" "t" "g" "a" "a" "a" "a" "t" "g" "a" "c" "a" "t" "t"
    [3511] "a" "a" "g" "g" "a" "a" "a" "g" "t" "t" "c" "t" "g" "c" "t" "g" "t" "t"
    [3529] "t" "t" "t" "a" "g" "c" "a" "a" "a" "a" "g" "c" "g" "t" "c" "c" "a" "g"
    [3547] "a" "a" "a" "g" "g" "a" "g" "a" "g" "c" "t" "t" "a" "g" "c" "a" "g" "g"
    [3565] "a" "g" "t" "c" "c" "t" "a" "g" "c" "c" "c" "t" "t" "t" "c" "a" "c" "c"
    [3583] "c" "a" "t" "a" "c" "a" "c" "a" "t" "t" "t" "g" "g" "c" "t" "c" "a" "g"
    [3601] "g" "g" "t" "t" "a" "c" "c" "g" "a" "a" "g" "a" "g" "g" "g" "g" "c" "c"
    [3619] "a" "a" "g" "a" "a" "a" "t" "t" "a" "g" "a" "g" "t" "c" "c" "t" "c" "a"
    [3637] "g" "a" "a" "g" "a" "g" "a" "a" "c" "t" "t" "a" "t" "c" "t" "a" "g" "t"
    [3655] "g" "a" "g" "g" "a" "t" "g" "a" "a" "g" "a" "g" "c" "t" "t" "c" "c" "c"
    [3673] "t" "g" "c" "t" "t" "c" "c" "a" "a" "c" "a" "c" "t" "t" "g" "t" "t" "a"
    [3691] "t" "t" "t" "g" "g" "t" "a" "a" "a" "g" "t" "a" "a" "a" "c" "a" "a" "t"
    [3709] "a" "t" "a" "c" "c" "t" "t" "c" "t" "c" "a" "g" "t" "c" "t" "a" "c" "t"
    [3727] "a" "g" "g" "c" "a" "t" "a" "g" "c" "a" "c" "c" "g" "t" "t" "g" "c" "t"
    [3745] "a" "c" "c" "g" "a" "g" "t" "g" "t" "c" "t" "g" "t" "c" "t" "a" "a" "g"
    [3763] "a" "a" "c" "a" "c" "a" "g" "a" "g" "g" "a" "g" "a" "a" "t" "t" "t" "a"
    [3781] "t" "t" "a" "t" "c" "a" "t" "t" "g" "a" "a" "g" "a" "a" "t" "a" "g" "c"
    [3799] "t" "t" "a" "a" "a" "t" "g" "a" "c" "t" "g" "c" "a" "g" "t" "a" "a" "c"
    [3817] "c" "a" "g" "g" "t" "a" "a" "t" "a" "t" "t" "g" "g" "c" "a" "a" "a" "g"
    [3835] "g" "c" "a" "t" "c" "t" "c" "a" "g" "g" "a" "a" "c" "a" "t" "c" "a" "c"
    [3853] "c" "t" "t" "a" "g" "t" "g" "a" "g" "g" "a" "a" "a" "c" "a" "a" "a" "a"
    [3871] "t" "g" "t" "t" "c" "t" "g" "c" "t" "a" "g" "c" "t" "t" "g" "t" "t" "t"
    [3889] "t" "c" "t" "t" "c" "a" "c" "a" "g" "t" "g" "c" "a" "g" "t" "g" "a" "a"
    [3907] "t" "t" "g" "g" "a" "a" "g" "a" "c" "t" "t" "g" "a" "c" "t" "g" "c" "a"
    [3925] "a" "a" "t" "a" "c" "a" "a" "a" "c" "a" "c" "c" "c" "a" "g" "g" "a" "t"
    [3943] "c" "c" "t" "t" "t" "c" "t" "t" "g" "a" "t" "t" "g" "g" "t" "t" "c" "t"
    [3961] "t" "c" "c" "a" "a" "a" "c" "a" "a" "a" "t" "g" "a" "g" "g" "c" "a" "t"
    [3979] "c" "a" "g" "t" "c" "t" "g" "a" "a" "a" "g" "c" "c" "a" "g" "g" "g" "a"
    [3997] "g" "t" "t" "g" "g" "t" "c" "t" "g" "a" "g" "t" "g" "a" "c" "a" "a" "g"
    [4015] "g" "a" "a" "t" "t" "g" "g" "t" "t" "t" "c" "a" "g" "a" "t" "g" "a" "t"
    [4033] "g" "a" "a" "g" "a" "a" "a" "g" "a" "g" "g" "a" "a" "c" "g" "g" "g" "c"
    [4051] "t" "t" "g" "g" "a" "a" "g" "a" "a" "a" "a" "t" "a" "a" "t" "c" "a" "a"
    [4069] "g" "a" "a" "g" "a" "g" "c" "a" "a" "a" "g" "c" "a" "t" "g" "g" "a" "t"
    [4087] "t" "c" "a" "a" "a" "c" "t" "t" "a" "g" "g" "t" "g" "a" "a" "g" "c" "a"
    [4105] "g" "c" "a" "t" "c" "t" "g" "g" "g" "t" "g" "t" "g" "a" "g" "a" "g" "t"
    [4123] "g" "a" "a" "a" "c" "a" "a" "g" "c" "g" "t" "c" "t" "c" "t" "g" "a" "a"
    [4141] "g" "a" "c" "t" "g" "c" "t" "c" "a" "g" "g" "g" "c" "t" "a" "t" "c" "c"
    [4159] "t" "c" "t" "c" "a" "g" "a" "g" "t" "g" "a" "c" "a" "t" "t" "t" "t" "a"
    [4177] "a" "c" "c" "a" "c" "t" "c" "a" "g" "c" "a" "g" "a" "g" "g" "g" "a" "t"
    [4195] "a" "c" "c" "a" "t" "g" "c" "a" "a" "c" "a" "t" "a" "a" "c" "c" "t" "g"
    [4213] "a" "t" "a" "a" "a" "g" "c" "t" "c" "c" "a" "g" "c" "a" "g" "g" "a" "a"
    [4231] "a" "t" "g" "g" "c" "t" "g" "a" "a" "c" "t" "a" "g" "a" "a" "g" "c" "t"
    [4249] "g" "t" "g" "t" "t" "a" "g" "a" "a" "c" "a" "g" "c" "a" "t" "g" "g" "g"
    [4267] "a" "g" "c" "c" "a" "g" "c" "c" "t" "t" "c" "t" "a" "a" "c" "a" "g" "c"
    [4285] "t" "a" "c" "c" "c" "t" "t" "c" "c" "a" "t" "c" "a" "t" "a" "a" "g" "t"
    [4303] "g" "a" "c" "t" "c" "t" "t" "c" "t" "g" "c" "c" "c" "t" "t" "g" "a" "g"
    [4321] "g" "a" "c" "c" "t" "g" "c" "g" "a" "a" "a" "t" "c" "c" "a" "g" "a" "a"
    [4339] "c" "a" "a" "a" "g" "c" "a" "c" "a" "t" "c" "a" "g" "a" "a" "a" "a" "a"
    [4357] "g" "c" "a" "g" "t" "a" "t" "t" "a" "a" "c" "t" "t" "c" "a" "c" "a" "g"
    [4375] "a" "a" "a" "a" "g" "t" "a" "g" "t" "g" "a" "a" "t" "a" "c" "c" "c" "t"
    [4393] "a" "t" "a" "a" "g" "c" "c" "a" "g" "a" "a" "t" "c" "c" "a" "g" "a" "a"
    [4411] "g" "g" "c" "c" "t" "t" "t" "c" "t" "g" "c" "t" "g" "a" "c" "a" "a" "g"
    [4429] "t" "t" "t" "g" "a" "g" "g" "t" "g" "t" "c" "t" "g" "c" "a" "g" "a" "t"
    [4447] "a" "g" "t" "t" "c" "t" "a" "c" "c" "a" "g" "t" "a" "a" "a" "a" "a" "t"
    [4465] "a" "a" "a" "g" "a" "a" "c" "c" "a" "g" "g" "a" "g" "t" "g" "g" "a" "a"
    [4483] "a" "g" "g" "t" "c" "a" "t" "c" "c" "c" "c" "t" "t" "c" "t" "a" "a" "a"
    [4501] "t" "g" "c" "c" "c" "a" "t" "c" "a" "t" "t" "a" "g" "a" "t" "g" "a" "t"
    [4519] "a" "g" "g" "t" "g" "g" "t" "a" "c" "a" "t" "g" "c" "a" "c" "a" "g" "t"
    [4537] "t" "g" "c" "t" "c" "t" "g" "g" "g" "a" "g" "t" "c" "t" "t" "c" "a" "g"
    [4555] "a" "a" "t" "a" "g" "a" "a" "a" "c" "t" "a" "c" "c" "c" "a" "t" "c" "t"
    [4573] "c" "a" "a" "g" "a" "g" "g" "a" "g" "c" "t" "c" "a" "t" "t" "a" "a" "g"
    [4591] "g" "t" "t" "g" "t" "t" "g" "a" "t" "g" "t" "g" "g" "a" "g" "g" "a" "g"
    [4609] "c" "a" "a" "c" "a" "g" "c" "t" "g" "g" "a" "a" "g" "a" "g" "t" "c" "t"
    [4627] "g" "g" "g" "c" "c" "a" "c" "a" "c" "g" "a" "t" "t" "t" "g" "a" "c" "g"
    [4645] "g" "a" "a" "a" "c" "a" "t" "c" "t" "t" "a" "c" "t" "t" "g" "c" "c" "a"
    [4663] "a" "g" "g" "c" "a" "a" "g" "a" "t" "c" "t" "a" "g" "a" "g" "g" "g" "a"
    [4681] "a" "c" "c" "c" "c" "t" "t" "a" "c" "c" "t" "g" "g" "a" "a" "t" "c" "t"
    [4699] "g" "g" "a" "a" "t" "c" "a" "g" "c" "c" "t" "c" "t" "t" "c" "t" "c" "t"
    [4717] "g" "a" "t" "g" "a" "c" "c" "c" "t" "g" "a" "a" "t" "c" "t" "g" "a" "t"
    [4735] "c" "c" "t" "t" "c" "t" "g" "a" "a" "g" "a" "c" "a" "g" "a" "g" "c" "c"
    [4753] "c" "c" "a" "g" "a" "g" "t" "c" "a" "g" "c" "t" "c" "g" "t" "g" "t" "t"
    [4771] "g" "g" "c" "a" "a" "c" "a" "t" "a" "c" "c" "a" "t" "c" "t" "t" "c" "a"
    [4789] "a" "c" "c" "t" "c" "t" "g" "c" "a" "t" "t" "g" "a" "a" "a" "g" "t" "t"
    [4807] "c" "c" "c" "c" "a" "a" "t" "t" "g" "a" "a" "a" "g" "t" "t" "g" "c" "a"
    [4825] "g" "a" "a" "t" "c" "t" "g" "c" "c" "c" "a" "g" "a" "g" "t" "c" "c" "a"
    [4843] "g" "c" "t" "g" "c" "t" "g" "c" "t" "c" "a" "t" "a" "c" "t" "a" "c" "t"
    [4861] "g" "a" "t" "a" "c" "t" "g" "c" "t" "g" "g" "g" "t" "a" "t" "a" "a" "t"
    [4879] "g" "c" "a" "a" "t" "g" "g" "a" "a" "g" "a" "a" "a" "g" "t" "g" "t" "g"
    [4897] "a" "g" "c" "a" "g" "g" "g" "a" "g" "a" "a" "g" "c" "c" "a" "g" "a" "a"
    [4915] "t" "t" "g" "a" "c" "a" "g" "c" "t" "t" "c" "a" "a" "c" "a" "g" "a" "a"
    [4933] "a" "g" "g" "g" "t" "c" "a" "a" "c" "a" "a" "a" "a" "g" "a" "a" "t" "g"
    [4951] "t" "c" "c" "a" "t" "g" "g" "t" "g" "g" "t" "g" "t" "c" "t" "g" "g" "c"
    [4969] "c" "t" "g" "a" "c" "c" "c" "c" "a" "g" "a" "a" "g" "a" "a" "t" "t" "t"
    [4987] "a" "t" "g" "c" "t" "c" "g" "t" "g" "t" "a" "c" "a" "a" "g" "t" "t" "t"
    [5005] "g" "c" "c" "a" "g" "a" "a" "a" "a" "c" "a" "c" "c" "a" "c" "a" "t" "c"
    [5023] "a" "c" "t" "t" "t" "a" "a" "c" "t" "a" "a" "t" "c" "t" "a" "a" "t" "t"
    [5041] "a" "c" "t" "g" "a" "a" "g" "a" "g" "a" "c" "t" "a" "c" "t" "c" "a" "t"
    [5059] "g" "t" "t" "g" "t" "t" "a" "t" "g" "a" "a" "a" "a" "c" "a" "g" "a" "t"
    [5077] "g" "c" "t" "g" "a" "g" "t" "t" "t" "g" "t" "g" "t" "g" "t" "g" "a" "a"
    [5095] "c" "g" "g" "a" "c" "a" "c" "t" "g" "a" "a" "a" "t" "a" "t" "t" "t" "t"
    [5113] "c" "t" "a" "g" "g" "a" "a" "t" "t" "g" "c" "g" "g" "g" "a" "g" "g" "a"
    [5131] "a" "a" "a" "t" "g" "g" "g" "t" "a" "g" "t" "t" "a" "g" "c" "t" "a" "t"
    [5149] "t" "t" "c" "t" "g" "g" "g" "t" "g" "a" "c" "c" "c" "a" "g" "t" "c" "t"
    [5167] "a" "t" "t" "a" "a" "a" "g" "a" "a" "a" "g" "a" "a" "a" "a" "a" "t" "g"
    [5185] "c" "t" "g" "a" "a" "t" "g" "a" "g" "c" "a" "t" "g" "a" "t" "t" "t" "t"
    [5203] "g" "a" "a" "g" "t" "c" "a" "g" "a" "g" "g" "a" "g" "a" "t" "g" "t" "g"
    [5221] "g" "t" "c" "a" "a" "t" "g" "g" "a" "a" "g" "a" "a" "a" "c" "c" "a" "c"
    [5239] "c" "a" "a" "g" "g" "t" "c" "c" "a" "a" "a" "g" "c" "g" "a" "g" "c" "a"
    [5257] "a" "g" "a" "g" "a" "a" "t" "c" "c" "c" "a" "g" "g" "a" "c" "a" "g" "a"
    [5275] "a" "a" "g" "a" "t" "c" "t" "t" "c" "a" "g" "g" "g" "g" "g" "c" "t" "a"
    [5293] "g" "a" "a" "a" "t" "c" "t" "g" "t" "t" "g" "c" "t" "a" "t" "g" "g" "g"
    [5311] "c" "c" "c" "t" "t" "c" "a" "c" "c" "a" "a" "c" "a" "t" "g" "c" "c" "c"
    [5329] "a" "c" "a" "g" "a" "t" "c" "a" "a" "c" "t" "g" "g" "a" "a" "t" "g" "g"
    [5347] "a" "t" "g" "g" "t" "a" "c" "a" "g" "c" "t" "g" "t" "g" "t" "g" "g" "t"
    [5365] "g" "c" "t" "t" "c" "t" "g" "t" "g" "g" "t" "g" "a" "a" "g" "g" "a" "g"
    [5383] "c" "t" "t" "t" "c" "a" "t" "c" "a" "t" "t" "c" "a" "c" "c" "c" "t" "t"
    [5401] "g" "g" "c" "a" "c" "a" "g" "g" "t" "g" "t" "c" "c" "a" "c" "c" "c" "a"
    [5419] "a" "t" "t" "g" "t" "g" "g" "t" "t" "g" "t" "g" "c" "a" "g" "c" "c" "a"
    [5437] "g" "a" "t" "g" "c" "c" "t" "g" "g" "a" "c" "a" "g" "a" "g" "g" "a" "c"
    [5455] "a" "a" "t" "g" "g" "c" "t" "t" "c" "c" "a" "t" "g"
    attr(,"name")
    [1] "DQ478408.BRCA1"
    attr(,"Annot")
    [1] ">DQ478408.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $FJ940752.BRCA1
     [1] "a" "t" "g" "c" "t" "g" "a" "g" "t" "t" "t" "g" "t" "g" "t" "g" "t" "g" "a"
    [20] "a" "c" "g" "g" "a" "c" "a" "c" "t" "g" "a" "a" "a" "t" "a" "t" "t" "t" "t"
    [39] "c" "t" "a" "g" "g" "t" "g" "c" "g" "g" "g" "a" "g" "g" "a" "a" "a" "a" "t"
    [58] "g" "g" "g" "t" "a" "g" "t" "t" "a" "g" "c" "t" "a" "t" "t" "t" "c" "t"
    attr(,"name")
    [1] "FJ940752.BRCA1"
    attr(,"Annot")
    [1] ">FJ940752.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $HE600032
      [1] "t" "c" "t" "a" "c" "t" "a" "g" "g" "c" "a" "t" "a" "g" "c" "a" "c" "c"
     [19] "g" "t" "t" "g" "c" "t" "a" "c" "c" "g" "a" "g" "t" "g" "t" "c" "t" "g"
     [37] "t" "c" "t" "a" "a" "g" "a" "a" "c" "a" "c" "a" "g" "a" "g" "g" "a" "g"
     [55] "a" "a" "t" "t" "t" "a" "t" "t" "a" "t" "c" "a" "t" "t" "g" "a" "a" "g"
     [73] "a" "a" "t" "a" "g" "c" "t" "t" "a" "a" "a" "t" "g" "a" "c" "t" "g" "c"
     [91] "a" "g" "t" "a" "a" "c" "c" "a" "g" "g" "t" "a" "a" "t" "a" "t" "t" "g"
    [109] "g" "c" "a" "a" "a" "g" "g" "c" "a" "t" "c" "t" "c" "a" "g" "g" "a" "a"
    [127] "c" "a" "t" "c" "a" "c" "c" "t" "t" "a" "g"
    attr(,"name")
    [1] "HE600032"
    attr(,"Annot")
    [1] ">HE600032"
    attr(,"class")
    [1] "SeqFastadna"
    
    $HE600033.BRCA1
      [1] "t" "c" "t" "a" "c" "t" "a" "g" "g" "c" "a" "t" "a" "g" "c" "a" "c" "c"
     [19] "g" "t" "t" "g" "c" "t" "a" "c" "c" "g" "a" "g" "t" "g" "t" "c" "t" "g"
     [37] "t" "c" "t" "a" "a" "g" "a" "a" "c" "a" "c" "a" "g" "a" "g" "g" "a" "g"
     [55] "a" "a" "t" "t" "t" "a" "t" "t" "a" "t" "c" "a" "t" "t" "g" "a" "a" "g"
     [73] "a" "a" "t" "a" "g" "c" "t" "t" "a" "a" "a" "t" "g" "a" "c" "t" "g" "c"
     [91] "a" "g" "t" "a" "a" "c" "c" "a" "g" "g" "t" "a" "a" "t" "a" "t" "t" "g"
    [109] "g" "c" "a" "a" "a" "g" "g" "c" "a" "t" "c" "t" "c" "a" "g" "g" "a" "a"
    [127] "c" "a" "t" "c" "a" "c" "c" "t" "t" "a" "g" "t" "g" "a" "g" "g" "a" "a"
    [145] "a" "c" "a" "a" "a" "a"
    attr(,"name")
    [1] "HE600033.BRCA1"
    attr(,"Annot")
    [1] ">HE600033.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $HE600034
      [1] "t" "a" "c" "a" "g" "c" "a" "t" "g" "a" "g" "a" "a" "c" "a" "g" "c" "a"
     [19] "g" "t" "t" "t" "a" "t" "t" "a" "c" "t" "c" "a" "c" "t" "a" "a" "a" "g"
     [37] "a" "c" "a" "g" "a" "a" "t" "g" "a" "a" "t" "g" "t" "a" "g" "a" "a" "a"
     [55] "a" "g" "g" "c" "t" "g" "a" "a" "t" "t" "c" "t" "g" "t" "a" "a" "t" "a"
     [73] "a" "a" "a" "g" "c" "a" "a" "a" "c" "a" "g" "c" "c" "t" "g" "g" "c" "t"
     [91] "t" "a" "g" "c" "a" "a" "g" "g" "a" "g" "c" "c" "a" "a" "c" "a" "t" "a"
    [109] "a" "c" "a" "g" "a" "t" "g" "g" "g" "c" "t" "g" "g" "a" "a" "g" "t" "a"
    [127] "a" "g" "g" "a" "a" "a" "c" "a" "t" "g" "t" "a" "a"
    attr(,"name")
    [1] "HE600034"
    attr(,"Annot")
    [1] ">HE600034"
    attr(,"class")
    [1] "SeqFastadna"
    
    $HE600035
      [1] "t" "a" "c" "a" "g" "c" "a" "t" "g" "a" "g" "a" "a" "c" "a" "g" "c" "a"
     [19] "g" "t" "t" "t" "a" "t" "t" "a" "c" "t" "c" "a" "c" "t" "a" "a" "a" "g"
     [37] "a" "c" "a" "g" "a" "a" "t" "g" "a" "a" "t" "g" "t" "a" "g" "a" "a" "a"
     [55] "a" "g" "g" "c" "t" "g" "a" "a" "t" "t" "c" "t" "g" "t" "a" "a" "t" "a"
     [73] "a" "a" "a" "g" "c" "a" "a" "a" "c" "a" "g" "c" "c" "t" "g" "g" "c" "t"
     [91] "t" "a" "g" "c" "a" "a" "g" "g" "a" "g" "c" "c" "a" "a" "c" "a" "t" "a"
    [109] "a" "c" "a" "g" "a" "t" "g" "g" "g" "c" "t" "g" "g" "a" "a" "g" "t" "a"
    [127] "a" "g" "g" "a" "a" "a" "c" "a" "t" "g" "t" "a" "a"
    attr(,"name")
    [1] "HE600035"
    attr(,"Annot")
    [1] ">HE600035"
    attr(,"class")
    [1] "SeqFastadna"
    
    $HE600036
      [1] "t" "a" "c" "a" "g" "c" "a" "t" "g" "a" "g" "a" "a" "c" "a" "g" "c" "a"
     [19] "g" "t" "t" "t" "a" "t" "t" "a" "c" "t" "c" "a" "c" "t" "a" "a" "a" "g"
     [37] "a" "c" "a" "g" "a" "a" "t" "g" "a" "a" "t" "g" "t" "a" "g" "a" "a" "a"
     [55] "a" "g" "g" "c" "t" "g" "a" "a" "t" "t" "c" "t" "g" "t" "a" "a" "t" "a"
     [73] "a" "a" "a" "g" "c" "a" "a" "a" "c" "a" "g" "c" "c" "t" "g" "g" "c" "t"
     [91] "t" "a" "g" "c" "a" "a" "g" "g" "a" "g" "c" "c" "a" "a" "c" "a" "t" "a"
    [109] "a" "c" "a" "g" "a" "t" "g" "g" "g" "c" "t" "g" "g" "a" "a" "g" "t" "a"
    [127] "a" "g" "g" "a" "a" "a" "c" "a" "t" "g" "t" "a" "a"
    attr(,"name")
    [1] "HE600036"
    attr(,"Annot")
    [1] ">HE600036"
    attr(,"class")
    [1] "SeqFastadna"
    
    $HE600037
      [1] "a" "c" "a" "g" "c" "a" "t" "g" "a" "g" "a" "a" "c" "a" "g" "c" "a" "g"
     [19] "t" "t" "t" "a" "t" "t" "a" "c" "t" "c" "a" "c" "t" "a" "a" "a" "g" "a"
     [37] "c" "a" "g" "a" "a" "t" "g" "a" "a" "t" "g" "t" "a" "g" "a" "a" "a" "a"
     [55] "g" "g" "c" "t" "g" "a" "a" "t" "t" "c" "t" "g" "t" "a" "a" "t" "a" "a"
     [73] "a" "a" "g" "c" "a" "a" "a" "c" "a" "g" "c" "c" "t" "g" "g" "c" "t" "t"
     [91] "a" "g" "c" "a" "a" "g" "g" "a" "g" "c" "c" "a" "a" "c" "a" "t" "a" "a"
    [109] "c" "a" "g" "a" "t" "g" "g" "g" "c" "t" "g" "g" "a" "a" "g" "t" "a" "a"
    [127] "g" "g" "a" "a" "a" "c"
    attr(,"name")
    [1] "HE600037"
    attr(,"Annot")
    [1] ">HE600037"
    attr(,"class")
    [1] "SeqFastadna"
    
    $HE600038
      [1] "t" "t" "a" "c" "a" "g" "c" "a" "t" "g" "a" "g" "a" "a" "c" "a" "g" "c"
     [19] "a" "g" "t" "t" "t" "a" "t" "t" "a" "c" "t" "c" "a" "c" "t" "a" "a" "a"
     [37] "g" "a" "c" "a" "g" "a" "a" "t" "g" "a" "a" "t" "g" "t" "a" "g" "a" "a"
     [55] "a" "a" "g" "g" "c" "t" "g" "a" "a" "t" "t" "c" "t" "g" "t" "a" "a" "t"
     [73] "a" "a" "a" "a" "g" "c" "a" "a" "a" "c" "a" "g" "c" "c" "t" "g" "g" "c"
     [91] "t" "t" "a" "g" "c" "a" "a" "g" "g" "a" "g" "c" "c" "a" "a" "c" "a" "t"
    [109] "a" "a" "c" "a" "g" "a" "t" "g" "g" "g" "c" "t" "g" "g" "a" "a" "g" "t"
    [127] "a" "a" "g" "g" "a" "a" "a" "c" "a" "t" "g" "t" "a" "a"
    attr(,"name")
    [1] "HE600038"
    attr(,"Annot")
    [1] ">HE600038"
    attr(,"class")
    [1] "SeqFastadna"
    
    $HSU14680.BRCA1
       [1] "a" "t" "g" "g" "a" "t" "t" "t" "a" "t" "c" "t" "g" "c" "t" "c" "t" "t"
      [19] "c" "g" "c" "g" "t" "t" "g" "a" "a" "g" "a" "a" "g" "t" "a" "c" "a" "a"
      [37] "a" "a" "t" "g" "t" "c" "a" "t" "t" "a" "a" "t" "g" "c" "t" "a" "t" "g"
      [55] "c" "a" "g" "a" "a" "a" "a" "t" "c" "t" "t" "a" "g" "a" "g" "t" "g" "t"
      [73] "c" "c" "c" "a" "t" "c" "t" "g" "t" "c" "t" "g" "g" "a" "g" "t" "t" "g"
      [91] "a" "t" "c" "a" "a" "g" "g" "a" "a" "c" "c" "t" "g" "t" "c" "t" "c" "c"
     [109] "a" "c" "a" "a" "a" "g" "t" "g" "t" "g" "a" "c" "c" "a" "c" "a" "t" "a"
     [127] "t" "t" "t" "t" "g" "c" "a" "a" "a" "t" "t" "t" "t" "g" "c" "a" "t" "g"
     [145] "c" "t" "g" "a" "a" "a" "c" "t" "t" "c" "t" "c" "a" "a" "c" "c" "a" "g"
     [163] "a" "a" "g" "a" "a" "a" "g" "g" "g" "c" "c" "t" "t" "c" "a" "c" "a" "g"
     [181] "t" "g" "t" "c" "c" "t" "t" "t" "a" "t" "g" "t" "a" "a" "g" "a" "a" "t"
     [199] "g" "a" "t" "a" "t" "a" "a" "c" "c" "a" "a" "a" "a" "g" "g" "a" "g" "c"
     [217] "c" "t" "a" "c" "a" "a" "g" "a" "a" "a" "g" "t" "a" "c" "g" "a" "g" "a"
     [235] "t" "t" "t" "a" "g" "t" "c" "a" "a" "c" "t" "t" "g" "t" "t" "g" "a" "a"
     [253] "g" "a" "g" "c" "t" "a" "t" "t" "g" "a" "a" "a" "a" "t" "c" "a" "t" "t"
     [271] "t" "g" "t" "g" "c" "t" "t" "t" "t" "c" "a" "g" "c" "t" "t" "g" "a" "c"
     [289] "a" "c" "a" "g" "g" "t" "t" "t" "g" "g" "a" "g" "t" "a" "t" "g" "c" "a"
     [307] "a" "a" "c" "a" "g" "c" "t" "a" "t" "a" "a" "t" "t" "t" "t" "g" "c" "a"
     [325] "a" "a" "a" "a" "a" "g" "g" "a" "a" "a" "a" "t" "a" "a" "c" "t" "c" "t"
     [343] "c" "c" "t" "g" "a" "a" "c" "a" "t" "c" "t" "a" "a" "a" "a" "g" "a" "t"
     [361] "g" "a" "a" "g" "t" "t" "t" "c" "t" "a" "t" "c" "a" "t" "c" "c" "a" "a"
     [379] "a" "g" "t" "a" "t" "g" "g" "g" "c" "t" "a" "c" "a" "g" "a" "a" "a" "c"
     [397] "c" "g" "t" "g" "c" "c" "a" "a" "a" "a" "g" "a" "c" "t" "t" "c" "t" "a"
     [415] "c" "a" "g" "a" "g" "t" "g" "a" "a" "c" "c" "c" "g" "a" "a" "a" "a" "t"
     [433] "c" "c" "t" "t" "c" "c" "t" "t" "g" "c" "a" "g" "g" "a" "a" "a" "c" "c"
     [451] "a" "g" "t" "c" "t" "c" "a" "g" "t" "g" "t" "c" "c" "a" "a" "c" "t" "c"
     [469] "t" "c" "t" "a" "a" "c" "c" "t" "t" "g" "g" "a" "a" "c" "t" "g" "t" "g"
     [487] "a" "g" "a" "a" "c" "t" "c" "t" "g" "a" "g" "g" "a" "c" "a" "a" "a" "g"
     [505] "c" "a" "g" "c" "g" "g" "a" "t" "a" "c" "a" "a" "c" "c" "t" "c" "a" "a"
     [523] "a" "a" "g" "a" "c" "g" "t" "c" "t" "g" "t" "c" "t" "a" "c" "a" "t" "t"
     [541] "g" "a" "a" "t" "t" "g" "g" "g" "a" "t" "c" "t" "g" "a" "t" "t" "c" "t"
     [559] "t" "c" "t" "g" "a" "a" "g" "a" "t" "a" "c" "c" "g" "t" "t" "a" "a" "t"
     [577] "a" "a" "g" "g" "c" "a" "a" "c" "t" "t" "a" "t" "t" "g" "c" "a" "g" "t"
     [595] "g" "t" "g" "g" "g" "a" "g" "a" "t" "c" "a" "a" "g" "a" "a" "t" "t" "g"
     [613] "t" "t" "a" "c" "a" "a" "a" "t" "c" "a" "c" "c" "c" "c" "t" "c" "a" "a"
     [631] "g" "g" "a" "a" "c" "c" "a" "g" "g" "g" "a" "t" "g" "a" "a" "a" "t" "c"
     [649] "a" "g" "t" "t" "t" "g" "g" "a" "t" "t" "c" "t" "g" "c" "a" "a" "a" "a"
     [667] "a" "a" "g" "g" "c" "t" "g" "c" "t" "t" "g" "t" "g" "a" "a" "t" "t" "t"
     [685] "t" "c" "t" "g" "a" "g" "a" "c" "g" "g" "a" "t" "g" "t" "a" "a" "c" "a"
     [703] "a" "a" "t" "a" "c" "t" "g" "a" "a" "c" "a" "t" "c" "a" "t" "c" "a" "a"
     [721] "c" "c" "c" "a" "g" "t" "a" "a" "t" "a" "a" "t" "g" "a" "t" "t" "t" "g"
     [739] "a" "a" "c" "a" "c" "c" "a" "c" "t" "g" "a" "g" "a" "a" "g" "c" "g" "t"
     [757] "g" "c" "a" "g" "c" "t" "g" "a" "g" "a" "g" "g" "c" "a" "t" "c" "c" "a"
     [775] "g" "a" "a" "a" "a" "g" "t" "a" "t" "c" "a" "g" "g" "g" "t" "a" "g" "t"
     [793] "t" "c" "t" "g" "t" "t" "t" "c" "a" "a" "a" "c" "t" "t" "g" "c" "a" "t"
     [811] "g" "t" "g" "g" "a" "g" "c" "c" "a" "t" "g" "t" "g" "g" "c" "a" "c" "a"
     [829] "a" "a" "t" "a" "c" "t" "c" "a" "t" "g" "c" "c" "a" "g" "c" "t" "c" "a"
     [847] "t" "t" "a" "c" "a" "g" "c" "a" "t" "g" "a" "g" "a" "a" "c" "a" "g" "c"
     [865] "a" "g" "t" "t" "t" "a" "t" "t" "a" "c" "t" "c" "a" "c" "t" "a" "a" "a"
     [883] "g" "a" "c" "a" "g" "a" "a" "t" "g" "a" "a" "t" "g" "t" "a" "g" "a" "a"
     [901] "a" "a" "g" "g" "c" "t" "g" "a" "a" "t" "t" "c" "t" "g" "t" "a" "a" "t"
     [919] "a" "a" "a" "a" "g" "c" "a" "a" "a" "c" "a" "g" "c" "c" "t" "g" "g" "c"
     [937] "t" "t" "a" "g" "c" "a" "a" "g" "g" "a" "g" "c" "c" "a" "a" "c" "a" "t"
     [955] "a" "a" "c" "a" "g" "a" "t" "g" "g" "g" "c" "t" "g" "g" "a" "a" "g" "t"
     [973] "a" "a" "g" "g" "a" "a" "a" "c" "a" "t" "g" "t" "a" "a" "t" "g" "a" "t"
     [991] "a" "g" "g" "c" "g" "g" "a" "c" "t" "c" "c" "c" "a" "g" "c" "a" "c" "a"
    [1009] "g" "a" "a" "a" "a" "a" "a" "a" "g" "g" "t" "a" "g" "a" "t" "c" "t" "g"
    [1027] "a" "a" "t" "g" "c" "t" "g" "a" "t" "c" "c" "c" "c" "t" "g" "t" "g" "t"
    [1045] "g" "a" "g" "a" "g" "a" "a" "a" "a" "g" "a" "a" "t" "g" "g" "a" "a" "t"
    [1063] "a" "a" "g" "c" "a" "g" "a" "a" "a" "c" "t" "g" "c" "c" "a" "t" "g" "c"
    [1081] "t" "c" "a" "g" "a" "g" "a" "a" "t" "c" "c" "t" "a" "g" "a" "g" "a" "t"
    [1099] "a" "c" "t" "g" "a" "a" "g" "a" "t" "g" "t" "t" "c" "c" "t" "t" "g" "g"
    [1117] "a" "t" "a" "a" "c" "a" "c" "t" "a" "a" "a" "t" "a" "g" "c" "a" "g" "c"
    [1135] "a" "t" "t" "c" "a" "g" "a" "a" "a" "g" "t" "t" "a" "a" "t" "g" "a" "g"
    [1153] "t" "g" "g" "t" "t" "t" "t" "c" "c" "a" "g" "a" "a" "g" "t" "g" "a" "t"
    [1171] "g" "a" "a" "c" "t" "g" "t" "t" "a" "g" "g" "t" "t" "c" "t" "g" "a" "t"
    [1189] "g" "a" "c" "t" "c" "a" "c" "a" "t" "g" "a" "t" "g" "g" "g" "g" "a" "g"
    [1207] "t" "c" "t" "g" "a" "a" "t" "c" "a" "a" "a" "t" "g" "c" "c" "a" "a" "a"
    [1225] "g" "t" "a" "g" "c" "t" "g" "a" "t" "g" "t" "a" "t" "t" "g" "g" "a" "c"
    [1243] "g" "t" "t" "c" "t" "a" "a" "a" "t" "g" "a" "g" "g" "t" "a" "g" "a" "t"
    [1261] "g" "a" "a" "t" "a" "t" "t" "c" "t" "g" "g" "t" "t" "c" "t" "t" "c" "a"
    [1279] "g" "a" "g" "a" "a" "a" "a" "t" "a" "g" "a" "c" "t" "t" "a" "c" "t" "g"
    [1297] "g" "c" "c" "a" "g" "t" "g" "a" "t" "c" "c" "t" "c" "a" "t" "g" "a" "g"
    [1315] "g" "c" "t" "t" "t" "a" "a" "t" "a" "t" "g" "t" "a" "a" "a" "a" "g" "t"
    [1333] "g" "a" "a" "a" "g" "a" "g" "t" "t" "c" "a" "c" "t" "c" "c" "a" "a" "a"
    [1351] "t" "c" "a" "g" "t" "a" "g" "a" "g" "a" "g" "t" "a" "a" "t" "a" "t" "t"
    [1369] "g" "a" "a" "g" "a" "c" "a" "a" "a" "a" "t" "a" "t" "t" "t" "g" "g" "g"
    [1387] "a" "a" "a" "a" "c" "c" "t" "a" "t" "c" "g" "g" "a" "a" "g" "a" "a" "g"
    [1405] "g" "c" "a" "a" "g" "c" "c" "t" "c" "c" "c" "c" "a" "a" "c" "t" "t" "a"
    [1423] "a" "g" "c" "c" "a" "t" "g" "t" "a" "a" "c" "t" "g" "a" "a" "a" "a" "t"
    [1441] "c" "t" "a" "a" "t" "t" "a" "t" "a" "g" "g" "a" "g" "c" "a" "t" "t" "t"
    [1459] "g" "t" "t" "a" "c" "t" "g" "a" "g" "c" "c" "a" "c" "a" "g" "a" "t" "a"
    [1477] "a" "t" "a" "c" "a" "a" "g" "a" "g" "c" "g" "t" "c" "c" "c" "c" "t" "c"
    [1495] "a" "c" "a" "a" "a" "t" "a" "a" "a" "t" "t" "a" "a" "a" "g" "c" "g" "t"
    [1513] "a" "a" "a" "a" "g" "g" "a" "g" "a" "c" "c" "t" "a" "c" "a" "t" "c" "a"
    [1531] "g" "g" "c" "c" "t" "t" "c" "a" "t" "c" "c" "t" "g" "a" "g" "g" "a" "t"
    [1549] "t" "t" "t" "a" "t" "c" "a" "a" "g" "a" "a" "a" "g" "c" "a" "g" "a" "t"
    [1567] "t" "t" "g" "g" "c" "a" "g" "t" "t" "c" "a" "a" "a" "a" "g" "a" "c" "t"
    [1585] "c" "c" "t" "g" "a" "a" "a" "t" "g" "a" "t" "a" "a" "a" "t" "c" "a" "g"
    [1603] "g" "g" "a" "a" "c" "t" "a" "a" "c" "c" "a" "a" "a" "c" "g" "g" "a" "g"
    [1621] "c" "a" "g" "a" "a" "t" "g" "g" "t" "c" "a" "a" "g" "t" "g" "a" "t" "g"
    [1639] "a" "a" "t" "a" "t" "t" "a" "c" "t" "a" "a" "t" "a" "g" "t" "g" "g" "t"
    [1657] "c" "a" "t" "g" "a" "g" "a" "a" "t" "a" "a" "a" "a" "c" "a" "a" "a" "a"
    [1675] "g" "g" "t" "g" "a" "t" "t" "c" "t" "a" "t" "t" "c" "a" "g" "a" "a" "t"
    [1693] "g" "a" "g" "a" "a" "a" "a" "a" "t" "c" "c" "t" "a" "a" "c" "c" "c" "a"
    [1711] "a" "t" "a" "g" "a" "a" "t" "c" "a" "c" "t" "c" "g" "a" "a" "a" "a" "a"
    [1729] "g" "a" "a" "t" "c" "t" "g" "c" "t" "t" "t" "c" "a" "a" "a" "a" "c" "g"
    [1747] "a" "a" "a" "g" "c" "t" "g" "a" "a" "c" "c" "t" "a" "t" "a" "a" "g" "c"
    [1765] "a" "g" "c" "a" "g" "t" "a" "t" "a" "a" "g" "c" "a" "a" "t" "a" "t" "g"
    [1783] "g" "a" "a" "c" "t" "c" "g" "a" "a" "t" "t" "a" "a" "a" "t" "a" "t" "c"
    [1801] "c" "a" "c" "a" "a" "t" "t" "c" "a" "a" "a" "a" "g" "c" "a" "c" "c" "t"
    [1819] "a" "a" "a" "a" "a" "g" "a" "a" "t" "a" "g" "g" "c" "t" "g" "a" "g" "g"
    [1837] "a" "g" "g" "a" "a" "g" "t" "c" "t" "t" "c" "t" "a" "c" "c" "a" "g" "g"
    [1855] "c" "a" "t" "a" "t" "t" "c" "a" "t" "g" "c" "g" "c" "t" "t" "g" "a" "a"
    [1873] "c" "t" "a" "g" "t" "a" "g" "t" "c" "a" "g" "t" "a" "g" "a" "a" "a" "t"
    [1891] "c" "t" "a" "a" "g" "c" "c" "c" "a" "c" "c" "t" "a" "a" "t" "t" "g" "t"
    [1909] "a" "c" "t" "g" "a" "a" "t" "t" "g" "c" "a" "a" "a" "t" "t" "g" "a" "t"
    [1927] "a" "g" "t" "t" "g" "t" "t" "c" "t" "a" "g" "c" "a" "g" "t" "g" "a" "a"
    [1945] "g" "a" "g" "a" "t" "a" "a" "a" "g" "a" "a" "a" "a" "a" "a" "a" "a" "g"
    [1963] "t" "a" "c" "a" "a" "c" "c" "a" "a" "a" "t" "g" "c" "c" "a" "g" "t" "c"
    [1981] "a" "g" "g" "c" "a" "c" "a" "g" "c" "a" "g" "a" "a" "a" "c" "c" "t" "a"
    [1999] "c" "a" "a" "c" "t" "c" "a" "t" "g" "g" "a" "a" "g" "g" "t" "a" "a" "a"
    [2017] "g" "a" "a" "c" "c" "t" "g" "c" "a" "a" "c" "t" "g" "g" "a" "g" "c" "c"
    [2035] "a" "a" "g" "a" "a" "g" "a" "g" "t" "a" "a" "c" "a" "a" "g" "c" "c" "a"
    [2053] "a" "a" "t" "g" "a" "a" "c" "a" "g" "a" "c" "a" "a" "g" "t" "a" "a" "a"
    [2071] "a" "g" "a" "c" "a" "t" "g" "a" "c" "a" "g" "c" "g" "a" "t" "a" "c" "t"
    [2089] "t" "t" "c" "c" "c" "a" "g" "a" "g" "c" "t" "g" "a" "a" "g" "t" "t" "a"
    [2107] "a" "c" "a" "a" "a" "t" "g" "c" "a" "c" "c" "t" "g" "g" "t" "t" "c" "t"
    [2125] "t" "t" "t" "a" "c" "t" "a" "a" "g" "t" "g" "t" "t" "c" "a" "a" "a" "t"
    [2143] "a" "c" "c" "a" "g" "t" "g" "a" "a" "c" "t" "t" "a" "a" "a" "g" "a" "a"
    [2161] "t" "t" "t" "g" "t" "c" "a" "a" "t" "c" "c" "t" "a" "g" "c" "c" "t" "t"
    [2179] "c" "c" "a" "a" "g" "a" "g" "a" "a" "g" "a" "a" "a" "a" "a" "g" "a" "a"
    [2197] "g" "a" "g" "a" "a" "a" "c" "t" "a" "g" "a" "a" "a" "c" "a" "g" "t" "t"
    [2215] "a" "a" "a" "g" "t" "g" "t" "c" "t" "a" "a" "t" "a" "a" "t" "g" "c" "t"
    [2233] "g" "a" "a" "g" "a" "c" "c" "c" "c" "a" "a" "a" "g" "a" "t" "c" "t" "c"
    [2251] "a" "t" "g" "t" "t" "a" "a" "g" "t" "g" "g" "a" "g" "a" "a" "a" "g" "g"
    [2269] "g" "t" "t" "t" "t" "g" "c" "a" "a" "a" "c" "t" "g" "a" "a" "a" "g" "a"
    [2287] "t" "c" "t" "g" "t" "a" "g" "a" "g" "a" "g" "t" "a" "g" "c" "a" "g" "t"
    [2305] "a" "t" "t" "t" "c" "a" "t" "t" "g" "g" "t" "a" "c" "c" "t" "g" "g" "t"
    [2323] "a" "c" "t" "g" "a" "t" "t" "a" "t" "g" "g" "c" "a" "c" "t" "c" "a" "g"
    [2341] "g" "a" "a" "a" "g" "t" "a" "t" "c" "t" "c" "g" "t" "t" "a" "c" "t" "g"
    [2359] "g" "a" "a" "g" "t" "t" "a" "g" "c" "a" "c" "t" "c" "t" "a" "g" "g" "g"
    [2377] "a" "a" "g" "g" "c" "a" "a" "a" "a" "a" "c" "a" "g" "a" "a" "c" "c" "a"
    [2395] "a" "a" "t" "a" "a" "a" "t" "g" "t" "g" "t" "g" "a" "g" "t" "c" "a" "g"
    [2413] "t" "g" "t" "g" "c" "a" "g" "c" "a" "t" "t" "t" "g" "a" "a" "a" "a" "c"
    [2431] "c" "c" "c" "a" "a" "g" "g" "g" "a" "c" "t" "a" "a" "t" "t" "c" "a" "t"
    [2449] "g" "g" "t" "t" "g" "t" "t" "c" "c" "a" "a" "a" "g" "a" "t" "a" "a" "t"
    [2467] "a" "g" "a" "a" "a" "t" "g" "a" "c" "a" "c" "a" "g" "a" "a" "g" "g" "c"
    [2485] "t" "t" "t" "a" "a" "g" "t" "a" "t" "c" "c" "a" "t" "t" "g" "g" "g" "a"
    [2503] "c" "a" "t" "g" "a" "a" "g" "t" "t" "a" "a" "c" "c" "a" "c" "a" "g" "t"
    [2521] "c" "g" "g" "g" "a" "a" "a" "c" "a" "a" "g" "c" "a" "t" "a" "g" "a" "a"
    [2539] "a" "t" "g" "g" "a" "a" "g" "a" "a" "a" "g" "t" "g" "a" "a" "c" "t" "t"
    [2557] "g" "a" "t" "g" "c" "t" "c" "a" "g" "t" "a" "t" "t" "t" "g" "c" "a" "g"
    [2575] "a" "a" "t" "a" "c" "a" "t" "t" "c" "a" "a" "g" "g" "t" "t" "t" "c" "a"
    [2593] "a" "a" "g" "c" "g" "c" "c" "a" "g" "t" "c" "a" "t" "t" "t" "g" "c" "t"
    [2611] "c" "c" "g" "t" "t" "t" "t" "c" "a" "a" "a" "t" "c" "c" "a" "g" "g" "a"
    [2629] "a" "a" "t" "g" "c" "a" "g" "a" "a" "g" "a" "g" "g" "a" "a" "t" "g" "t"
    [2647] "g" "c" "a" "a" "c" "a" "t" "t" "c" "t" "c" "t" "g" "c" "c" "c" "a" "c"
    [2665] "t" "c" "t" "g" "g" "g" "t" "c" "c" "t" "t" "a" "a" "a" "g" "a" "a" "a"
    [2683] "c" "a" "a" "a" "g" "t" "c" "c" "a" "a" "a" "a" "g" "t" "c" "a" "c" "t"
    [2701] "t" "t" "t" "g" "a" "a" "t" "g" "t" "g" "a" "a" "c" "a" "a" "a" "a" "g"
    [2719] "g" "a" "a" "g" "a" "a" "a" "a" "t" "c" "a" "a" "g" "g" "a" "a" "a" "g"
    [2737] "a" "a" "t" "g" "a" "g" "t" "c" "t" "a" "a" "t" "a" "t" "c" "a" "a" "g"
    [2755] "c" "c" "t" "g" "t" "a" "c" "a" "g" "a" "c" "a" "g" "t" "t" "a" "a" "t"
    [2773] "a" "t" "c" "a" "c" "t" "g" "c" "a" "g" "g" "c" "t" "t" "t" "c" "c" "t"
    [2791] "g" "t" "g" "g" "t" "t" "g" "g" "t" "c" "a" "g" "a" "a" "a" "g" "a" "t"
    [2809] "a" "a" "g" "c" "c" "a" "g" "t" "t" "g" "a" "t" "a" "a" "t" "g" "c" "c"
    [2827] "a" "a" "a" "t" "g" "t" "a" "g" "t" "a" "t" "c" "a" "a" "a" "g" "g" "a"
    [2845] "g" "g" "c" "t" "c" "t" "a" "g" "g" "t" "t" "t" "t" "g" "t" "c" "t" "a"
    [2863] "t" "c" "a" "t" "c" "t" "c" "a" "g" "t" "t" "c" "a" "g" "a" "g" "g" "c"
    [2881] "a" "a" "c" "g" "a" "a" "a" "c" "t" "g" "g" "a" "c" "t" "c" "a" "t" "t"
    [2899] "a" "c" "t" "c" "c" "a" "a" "a" "t" "a" "a" "a" "c" "a" "t" "g" "g" "a"
    [2917] "c" "t" "t" "t" "t" "a" "c" "a" "a" "a" "a" "c" "c" "c" "a" "t" "a" "t"
    [2935] "c" "g" "t" "a" "t" "a" "c" "c" "a" "c" "c" "a" "c" "t" "t" "t" "t" "t"
    [2953] "c" "c" "c" "a" "t" "c" "a" "a" "g" "t" "c" "a" "t" "t" "t" "g" "t" "t"
    [2971] "a" "a" "a" "a" "c" "t" "a" "a" "a" "t" "g" "t" "a" "a" "g" "a" "a" "a"
    [2989] "a" "a" "t" "c" "t" "g" "c" "t" "a" "g" "a" "g" "g" "a" "a" "a" "a" "c"
    [3007] "t" "t" "t" "g" "a" "g" "g" "a" "a" "c" "a" "t" "t" "c" "a" "a" "t" "g"
    [3025] "t" "c" "a" "c" "c" "t" "g" "a" "a" "a" "g" "a" "g" "a" "a" "a" "t" "g"
    [3043] "g" "g" "a" "a" "a" "t" "g" "a" "g" "a" "a" "c" "a" "t" "t" "c" "c" "a"
    [3061] "a" "g" "t" "a" "c" "a" "g" "t" "g" "a" "g" "c" "a" "c" "a" "a" "t" "t"
    [3079] "a" "g" "c" "c" "g" "t" "a" "a" "t" "a" "a" "c" "a" "t" "t" "a" "g" "a"
    [3097] "g" "a" "a" "a" "a" "t" "g" "t" "t" "t" "t" "t" "a" "a" "a" "g" "a" "a"
    [3115] "g" "c" "c" "a" "g" "c" "t" "c" "a" "a" "g" "c" "a" "a" "t" "a" "t" "t"
    [3133] "a" "a" "t" "g" "a" "a" "g" "t" "a" "g" "g" "t" "t" "c" "c" "a" "g" "t"
    [3151] "a" "c" "t" "a" "a" "t" "g" "a" "a" "g" "t" "g" "g" "g" "c" "t" "c" "c"
    [3169] "a" "g" "t" "a" "t" "t" "a" "a" "t" "g" "a" "a" "a" "t" "a" "g" "g" "t"
    [3187] "t" "c" "c" "a" "g" "t" "g" "a" "t" "g" "a" "a" "a" "a" "c" "a" "t" "t"
    [3205] "c" "a" "a" "g" "c" "a" "g" "a" "a" "c" "t" "a" "g" "g" "t" "a" "g" "a"
    [3223] "a" "a" "c" "a" "g" "a" "g" "g" "g" "c" "c" "a" "a" "a" "a" "t" "t" "g"
    [3241] "a" "a" "t" "g" "c" "t" "a" "t" "g" "c" "t" "t" "a" "g" "a" "t" "t" "a"
    [3259] "g" "g" "g" "g" "t" "t" "t" "t" "g" "c" "a" "a" "c" "c" "t" "g" "a" "g"
    [3277] "g" "t" "c" "t" "a" "t" "a" "a" "a" "c" "a" "a" "a" "g" "t" "c" "t" "t"
    [3295] "c" "c" "t" "g" "g" "a" "a" "g" "t" "a" "a" "t" "t" "g" "t" "a" "a" "g"
    [3313] "c" "a" "t" "c" "c" "t" "g" "a" "a" "a" "t" "a" "a" "a" "a" "a" "a" "g"
    [3331] "c" "a" "a" "g" "a" "a" "t" "a" "t" "g" "a" "a" "g" "a" "a" "g" "t" "a"
    [3349] "g" "t" "t" "c" "a" "g" "a" "c" "t" "g" "t" "t" "a" "a" "t" "a" "c" "a"
    [3367] "g" "a" "t" "t" "t" "c" "t" "c" "t" "c" "c" "a" "t" "a" "t" "c" "t" "g"
    [3385] "a" "t" "t" "t" "c" "a" "g" "a" "t" "a" "a" "c" "t" "t" "a" "g" "a" "a"
    [3403] "c" "a" "g" "c" "c" "t" "a" "t" "g" "g" "g" "a" "a" "g" "t" "a" "g" "t"
    [3421] "c" "a" "t" "g" "c" "a" "t" "c" "t" "c" "a" "g" "g" "t" "t" "t" "g" "t"
    [3439] "t" "c" "t" "g" "a" "g" "a" "c" "a" "c" "c" "t" "g" "a" "t" "g" "a" "c"
    [3457] "c" "t" "g" "t" "t" "a" "g" "a" "t" "g" "a" "t" "g" "g" "t" "g" "a" "a"
    [3475] "a" "t" "a" "a" "a" "g" "g" "a" "a" "g" "a" "t" "a" "c" "t" "a" "g" "t"
    [3493] "t" "t" "t" "g" "c" "t" "g" "a" "a" "a" "a" "t" "g" "a" "c" "a" "t" "t"
    [3511] "a" "a" "g" "g" "a" "a" "a" "g" "t" "t" "c" "t" "g" "c" "t" "g" "t" "t"
    [3529] "t" "t" "t" "a" "g" "c" "a" "a" "a" "a" "g" "c" "g" "t" "c" "c" "a" "g"
    [3547] "a" "a" "a" "g" "g" "a" "g" "a" "g" "c" "t" "t" "a" "g" "c" "a" "g" "g"
    [3565] "a" "g" "t" "c" "c" "t" "a" "g" "c" "c" "c" "t" "t" "t" "c" "a" "c" "c"
    [3583] "c" "a" "t" "a" "c" "a" "c" "a" "t" "t" "t" "g" "g" "c" "t" "c" "a" "g"
    [3601] "g" "g" "t" "t" "a" "c" "c" "g" "a" "a" "g" "a" "g" "g" "g" "g" "c" "c"
    [3619] "a" "a" "g" "a" "a" "a" "t" "t" "a" "g" "a" "g" "t" "c" "c" "t" "c" "a"
    [3637] "g" "a" "a" "g" "a" "g" "a" "a" "c" "t" "t" "a" "t" "c" "t" "a" "g" "t"
    [3655] "g" "a" "g" "g" "a" "t" "g" "a" "a" "g" "a" "g" "c" "t" "t" "c" "c" "c"
    [3673] "t" "g" "c" "t" "t" "c" "c" "a" "a" "c" "a" "c" "t" "t" "g" "t" "t" "a"
    [3691] "t" "t" "t" "g" "g" "t" "a" "a" "a" "g" "t" "a" "a" "a" "c" "a" "a" "t"
    [3709] "a" "t" "a" "c" "c" "t" "t" "c" "t" "c" "a" "g" "t" "c" "t" "a" "c" "t"
    [3727] "a" "g" "g" "c" "a" "t" "a" "g" "c" "a" "c" "c" "g" "t" "t" "g" "c" "t"
    [3745] "a" "c" "c" "g" "a" "g" "t" "g" "t" "c" "t" "g" "t" "c" "t" "a" "a" "g"
    [3763] "a" "a" "c" "a" "c" "a" "g" "a" "g" "g" "a" "g" "a" "a" "t" "t" "t" "a"
    [3781] "t" "t" "a" "t" "c" "a" "t" "t" "g" "a" "a" "g" "a" "a" "t" "a" "g" "c"
    [3799] "t" "t" "a" "a" "a" "t" "g" "a" "c" "t" "g" "c" "a" "g" "t" "a" "a" "c"
    [3817] "c" "a" "g" "g" "t" "a" "a" "t" "a" "t" "t" "g" "g" "c" "a" "a" "a" "g"
    [3835] "g" "c" "a" "t" "c" "t" "c" "a" "g" "g" "a" "a" "c" "a" "t" "c" "a" "c"
    [3853] "c" "t" "t" "a" "g" "t" "g" "a" "g" "g" "a" "a" "a" "c" "a" "a" "a" "a"
    [3871] "t" "g" "t" "t" "c" "t" "g" "c" "t" "a" "g" "c" "t" "t" "g" "t" "t" "t"
    [3889] "t" "c" "t" "t" "c" "a" "c" "a" "g" "t" "g" "c" "a" "g" "t" "g" "a" "a"
    [3907] "t" "t" "g" "g" "a" "a" "g" "a" "c" "t" "t" "g" "a" "c" "t" "g" "c" "a"
    [3925] "a" "a" "t" "a" "c" "a" "a" "a" "c" "a" "c" "c" "c" "a" "g" "g" "a" "t"
    [3943] "c" "c" "t" "t" "t" "c" "t" "t" "g" "a" "t" "t" "g" "g" "t" "t" "c" "t"
    [3961] "t" "c" "c" "a" "a" "a" "c" "a" "a" "a" "t" "g" "a" "g" "g" "c" "a" "t"
    [3979] "c" "a" "g" "t" "c" "t" "g" "a" "a" "a" "g" "c" "c" "a" "g" "g" "g" "a"
    [3997] "g" "t" "t" "g" "g" "t" "c" "t" "g" "a" "g" "t" "g" "a" "c" "a" "a" "g"
    [4015] "g" "a" "a" "t" "t" "g" "g" "t" "t" "t" "c" "a" "g" "a" "t" "g" "a" "t"
    [4033] "g" "a" "a" "g" "a" "a" "a" "g" "a" "g" "g" "a" "a" "c" "g" "g" "g" "c"
    [4051] "t" "t" "g" "g" "a" "a" "g" "a" "a" "a" "a" "t" "a" "a" "t" "c" "a" "a"
    [4069] "g" "a" "a" "g" "a" "g" "c" "a" "a" "a" "g" "c" "a" "t" "g" "g" "a" "t"
    [4087] "t" "c" "a" "a" "a" "c" "t" "t" "a" "g" "g" "t" "g" "a" "a" "g" "c" "a"
    [4105] "g" "c" "a" "t" "c" "t" "g" "g" "g" "t" "g" "t" "g" "a" "g" "a" "g" "t"
    [4123] "g" "a" "a" "a" "c" "a" "a" "g" "c" "g" "t" "c" "t" "c" "t" "g" "a" "a"
    [4141] "g" "a" "c" "t" "g" "c" "t" "c" "a" "g" "g" "g" "c" "t" "a" "t" "c" "c"
    [4159] "t" "c" "t" "c" "a" "g" "a" "g" "t" "g" "a" "c" "a" "t" "t" "t" "t" "a"
    [4177] "a" "c" "c" "a" "c" "t" "c" "a" "g" "c" "a" "g" "a" "g" "g" "g" "a" "t"
    [4195] "a" "c" "c" "a" "t" "g" "c" "a" "a" "c" "a" "t" "a" "a" "c" "c" "t" "g"
    [4213] "a" "t" "a" "a" "a" "g" "c" "t" "c" "c" "a" "g" "c" "a" "g" "g" "a" "a"
    [4231] "a" "t" "g" "g" "c" "t" "g" "a" "a" "c" "t" "a" "g" "a" "a" "g" "c" "t"
    [4249] "g" "t" "g" "t" "t" "a" "g" "a" "a" "c" "a" "g" "c" "a" "t" "g" "g" "g"
    [4267] "a" "g" "c" "c" "a" "g" "c" "c" "t" "t" "c" "t" "a" "a" "c" "a" "g" "c"
    [4285] "t" "a" "c" "c" "c" "t" "t" "c" "c" "a" "t" "c" "a" "t" "a" "a" "g" "t"
    [4303] "g" "a" "c" "t" "c" "t" "t" "c" "t" "g" "c" "c" "c" "t" "t" "g" "a" "g"
    [4321] "g" "a" "c" "c" "t" "g" "c" "g" "a" "a" "a" "t" "c" "c" "a" "g" "a" "a"
    [4339] "c" "a" "a" "a" "g" "c" "a" "c" "a" "t" "c" "a" "g" "a" "a" "a" "a" "a"
    [4357] "g" "c" "a" "g" "t" "a" "t" "t" "a" "a" "c" "t" "t" "c" "a" "c" "a" "g"
    [4375] "a" "a" "a" "a" "g" "t" "a" "g" "t" "g" "a" "a" "t" "a" "c" "c" "c" "t"
    [4393] "a" "t" "a" "a" "g" "c" "c" "a" "g" "a" "a" "t" "c" "c" "a" "g" "a" "a"
    [4411] "g" "g" "c" "c" "t" "t" "t" "c" "t" "g" "c" "t" "g" "a" "c" "a" "a" "g"
    [4429] "t" "t" "t" "g" "a" "g" "g" "t" "g" "t" "c" "t" "g" "c" "a" "g" "a" "t"
    [4447] "a" "g" "t" "t" "c" "t" "a" "c" "c" "a" "g" "t" "a" "a" "a" "a" "a" "t"
    [4465] "a" "a" "a" "g" "a" "a" "c" "c" "a" "g" "g" "a" "g" "t" "g" "g" "a" "a"
    [4483] "a" "g" "g" "t" "c" "a" "t" "c" "c" "c" "c" "t" "t" "c" "t" "a" "a" "a"
    [4501] "t" "g" "c" "c" "c" "a" "t" "c" "a" "t" "t" "a" "g" "a" "t" "g" "a" "t"
    [4519] "a" "g" "g" "t" "g" "g" "t" "a" "c" "a" "t" "g" "c" "a" "c" "a" "g" "t"
    [4537] "t" "g" "c" "t" "c" "t" "g" "g" "g" "a" "g" "t" "c" "t" "t" "c" "a" "g"
    [4555] "a" "a" "t" "a" "g" "a" "a" "a" "c" "t" "a" "c" "c" "c" "a" "t" "c" "t"
    [4573] "c" "a" "a" "g" "a" "g" "g" "a" "g" "c" "t" "c" "a" "t" "t" "a" "a" "g"
    [4591] "g" "t" "t" "g" "t" "t" "g" "a" "t" "g" "t" "g" "g" "a" "g" "g" "a" "g"
    [4609] "c" "a" "a" "c" "a" "g" "c" "t" "g" "g" "a" "a" "g" "a" "g" "t" "c" "t"
    [4627] "g" "g" "g" "c" "c" "a" "c" "a" "c" "g" "a" "t" "t" "t" "g" "a" "c" "g"
    [4645] "g" "a" "a" "a" "c" "a" "t" "c" "t" "t" "a" "c" "t" "t" "g" "c" "c" "a"
    [4663] "a" "g" "g" "c" "a" "a" "g" "a" "t" "c" "t" "a" "g" "a" "g" "g" "g" "a"
    [4681] "a" "c" "c" "c" "c" "t" "t" "a" "c" "c" "t" "g" "g" "a" "a" "t" "c" "t"
    [4699] "g" "g" "a" "a" "t" "c" "a" "g" "c" "c" "t" "c" "t" "t" "c" "t" "c" "t"
    [4717] "g" "a" "t" "g" "a" "c" "c" "c" "t" "g" "a" "a" "t" "c" "t" "g" "a" "t"
    [4735] "c" "c" "t" "t" "c" "t" "g" "a" "a" "g" "a" "c" "a" "g" "a" "g" "c" "c"
    [4753] "c" "c" "a" "g" "a" "g" "t" "c" "a" "g" "c" "t" "c" "g" "t" "g" "t" "t"
    [4771] "g" "g" "c" "a" "a" "c" "a" "t" "a" "c" "c" "a" "t" "c" "t" "t" "c" "a"
    [4789] "a" "c" "c" "t" "c" "t" "g" "c" "a" "t" "t" "g" "a" "a" "a" "g" "t" "t"
    [4807] "c" "c" "c" "c" "a" "a" "t" "t" "g" "a" "a" "a" "g" "t" "t" "g" "c" "a"
    [4825] "g" "a" "a" "t" "c" "t" "g" "c" "c" "c" "a" "g" "a" "g" "t" "c" "c" "a"
    [4843] "g" "c" "t" "g" "c" "t" "g" "c" "t" "c" "a" "t" "a" "c" "t" "a" "c" "t"
    [4861] "g" "a" "t" "a" "c" "t" "g" "c" "t" "g" "g" "g" "t" "a" "t" "a" "a" "t"
    [4879] "g" "c" "a" "a" "t" "g" "g" "a" "a" "g" "a" "a" "a" "g" "t" "g" "t" "g"
    [4897] "a" "g" "c" "a" "g" "g" "g" "a" "g" "a" "a" "g" "c" "c" "a" "g" "a" "a"
    [4915] "t" "t" "g" "a" "c" "a" "g" "c" "t" "t" "c" "a" "a" "c" "a" "g" "a" "a"
    [4933] "a" "g" "g" "g" "t" "c" "a" "a" "c" "a" "a" "a" "a" "g" "a" "a" "t" "g"
    [4951] "t" "c" "c" "a" "t" "g" "g" "t" "g" "g" "t" "g" "t" "c" "t" "g" "g" "c"
    [4969] "c" "t" "g" "a" "c" "c" "c" "c" "a" "g" "a" "a" "g" "a" "a" "t" "t" "t"
    [4987] "a" "t" "g" "c" "t" "c" "g" "t" "g" "t" "a" "c" "a" "a" "g" "t" "t" "t"
    [5005] "g" "c" "c" "a" "g" "a" "a" "a" "a" "c" "a" "c" "c" "a" "c" "a" "t" "c"
    [5023] "a" "c" "t" "t" "t" "a" "a" "c" "t" "a" "a" "t" "c" "t" "a" "a" "t" "t"
    [5041] "a" "c" "t" "g" "a" "a" "g" "a" "g" "a" "c" "t" "a" "c" "t" "c" "a" "t"
    [5059] "g" "t" "t" "g" "t" "t" "a" "t" "g" "a" "a" "a" "a" "c" "a" "g" "a" "t"
    [5077] "g" "c" "t" "g" "a" "g" "t" "t" "t" "g" "t" "g" "t" "g" "t" "g" "a" "a"
    [5095] "c" "g" "g" "a" "c" "a" "c" "t" "g" "a" "a" "a" "t" "a" "t" "t" "t" "t"
    [5113] "c" "t" "a" "g" "g" "a" "a" "t" "t" "g" "c" "g" "g" "g" "a" "g" "g" "a"
    [5131] "a" "a" "a" "t" "g" "g" "g" "t" "a" "g" "t" "t" "a" "g" "c" "t" "a" "t"
    [5149] "t" "t" "c" "t" "g" "g" "g" "t" "g" "a" "c" "c" "c" "a" "g" "t" "c" "t"
    [5167] "a" "t" "t" "a" "a" "a" "g" "a" "a" "a" "g" "a" "a" "a" "a" "a" "t" "g"
    [5185] "c" "t" "g" "a" "a" "t" "g" "a" "g" "c" "a" "t" "g" "a" "t" "t" "t" "t"
    [5203] "g" "a" "a" "g" "t" "c" "a" "g" "a" "g" "g" "a" "g" "a" "t" "g" "t" "g"
    [5221] "g" "t" "c" "a" "a" "t" "g" "g" "a" "a" "g" "a" "a" "a" "c" "c" "a" "c"
    [5239] "c" "a" "a" "g" "g" "t" "c" "c" "a" "a" "a" "g" "c" "g" "a" "g" "c" "a"
    [5257] "a" "g" "a" "g" "a" "a" "t" "c" "c" "c" "a" "g" "g" "a" "c" "a" "g" "a"
    [5275] "a" "a" "g" "a" "t" "c" "t" "t" "c" "a" "g" "g" "g" "g" "g" "c" "t" "a"
    [5293] "g" "a" "a" "a" "t" "c" "t" "g" "t" "t" "g" "c" "t" "a" "t" "g" "g" "g"
    [5311] "c" "c" "c" "t" "t" "c" "a" "c" "c" "a" "a" "c" "a" "t" "g" "c" "c" "c"
    [5329] "a" "c" "a" "g" "a" "t" "c" "a" "a" "c" "t" "g" "g" "a" "a" "t" "g" "g"
    [5347] "a" "t" "g" "g" "t" "a" "c" "a" "g" "c" "t" "g" "t" "g" "t" "g" "g" "t"
    [5365] "g" "c" "t" "t" "c" "t" "g" "t" "g" "g" "t" "g" "a" "a" "g" "g" "a" "g"
    [5383] "c" "t" "t" "t" "c" "a" "t" "c" "a" "t" "t" "c" "a" "c" "c" "c" "t" "t"
    [5401] "g" "g" "c" "a" "c" "a" "g" "g" "t" "g" "t" "c" "c" "a" "c" "c" "c" "a"
    [5419] "a" "t" "t" "g" "t" "g" "g" "t" "t" "g" "t" "g" "c" "a" "g" "c" "c" "a"
    [5437] "g" "a" "t" "g" "c" "c" "t" "g" "g" "a" "c" "a" "g" "a" "g" "g" "a" "c"
    [5455] "a" "a" "t" "g" "g" "c" "t" "t" "c" "c" "a" "t" "g" "c" "a" "a" "t" "t"
    [5473] "g" "g" "g" "c" "a" "g" "a" "t" "g" "t" "g" "t" "g" "a" "g" "g" "c" "a"
    [5491] "c" "c" "t" "g" "t" "g" "g" "t" "g" "a" "c" "c" "c" "g" "a" "g" "a" "g"
    [5509] "t" "g" "g" "g" "t" "g" "t" "t" "g" "g" "a" "c" "a" "g" "t" "g" "t" "a"
    [5527] "g" "c" "a" "c" "t" "c" "t" "a" "c" "c" "a" "g" "t" "g" "c" "c" "a" "g"
    [5545] "g" "a" "g" "c" "t" "g" "g" "a" "c" "a" "c" "c" "t" "a" "c" "c" "t" "g"
    [5563] "a" "t" "a" "c" "c" "c" "c" "a" "g" "a" "t" "c" "c" "c" "c" "c" "a" "c"
    [5581] "a" "g" "c" "c" "a" "c" "t" "a" "c" "t" "g" "a"
    attr(,"name")
    [1] "HSU14680.BRCA1"
    attr(,"Annot")
    [1] ">HSU14680.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $HSU18009
       [1] "g" "a" "a" "t" "t" "c" "g" "g" "c" "a" "c" "g" "a" "g" "g" "c" "g" "g"
      [19] "g" "c" "t" "c" "a" "a" "c" "t" "t" "c" "g" "c" "a" "g" "a" "c" "c" "t"
      [37] "c" "a" "t" "g" "g" "c" "t" "a" "g" "g" "c" "a" "g" "g" "g" "g" "c" "t"
      [55] "g" "t" "a" "c" "g" "a" "c" "c" "g" "t" "c" "t" "c" "c" "c" "g" "c" "c"
      [73] "t" "c" "t" "g" "c" "c" "t" "g" "t" "c" "a" "c" "t" "c" "c" "g" "g" "g"
      [91] "c" "a" "t" "g" "g" "a" "g" "g" "g" "c" "g" "c" "g" "g" "g" "t" "g" "t"
     [109] "t" "g" "t" "g" "a" "t" "c" "g" "c" "a" "g" "t" "g" "g" "g" "c" "g" "a"
     [127] "g" "g" "g" "a" "g" "t" "c" "a" "g" "c" "g" "a" "c" "c" "g" "c" "a" "a"
     [145] "g" "g" "c" "a" "g" "g" "a" "g" "a" "c" "c" "g" "g" "g" "t" "g" "a" "t"
     [163] "g" "g" "t" "g" "t" "t" "g" "a" "a" "c" "c" "g" "g" "t" "c" "a" "g" "g"
     [181] "g" "a" "t" "g" "t" "g" "g" "c" "a" "g" "g" "a" "a" "g" "a" "g" "g" "t"
     [199] "g" "a" "c" "t" "g" "t" "g" "c" "c" "c" "t" "c" "g" "g" "t" "c" "c" "a"
     [217] "g" "a" "c" "c" "t" "t" "c" "c" "t" "g" "a" "t" "t" "c" "c" "t" "g" "a"
     [235] "g" "g" "c" "c" "a" "t" "g" "a" "c" "c" "t" "t" "t" "g" "a" "g" "g" "a"
     [253] "a" "g" "c" "t" "g" "c" "t" "g" "c" "c" "t" "t" "g" "c" "t" "c" "g" "t"
     [271] "c" "a" "a" "t" "t" "a" "c" "a" "t" "t" "a" "c" "a" "g" "c" "c" "t" "a"
     [289] "c" "a" "t" "g" "g" "t" "c" "c" "t" "c" "t" "t" "t" "g" "a" "c" "t" "t"
     [307] "c" "g" "g" "c" "a" "a" "c" "c" "t" "a" "c" "a" "g" "c" "c" "t" "g" "g"
     [325] "c" "c" "a" "c" "a" "g" "c" "g" "t" "c" "t" "t" "g" "g" "t" "a" "c" "a"
     [343] "c" "a" "t" "g" "g" "c" "t" "g" "c" "a" "g" "g" "g" "g" "g" "t" "g" "t"
     [361] "g" "g" "g" "t" "a" "t" "g" "g" "c" "t" "g" "c" "c" "g" "t" "g" "c" "a"
     [379] "c" "g" "t" "g" "t" "g" "c" "c" "g" "t" "a" "c" "a" "g" "t" "g" "g" "a"
     [397] "g" "a" "a" "t" "g" "t" "g" "a" "c" "a" "g" "t" "g" "t" "t" "c" "g" "g"
     [415] "a" "a" "c" "g" "g" "c" "c" "t" "c" "g" "g" "c" "c" "a" "g" "c" "a" "a"
     [433] "g" "c" "a" "c" "g" "a" "g" "g" "c" "a" "c" "t" "g" "a" "a" "g" "g" "a"
     [451] "g" "a" "a" "t" "g" "g" "g" "g" "t" "c" "a" "c" "a" "c" "a" "t" "c" "c"
     [469] "c" "a" "t" "c" "g" "a" "c" "t" "a" "t" "c" "a" "c" "a" "c" "g" "a" "c"
     [487] "t" "g" "a" "c" "t" "a" "c" "g" "t" "g" "g" "a" "t" "g" "a" "g" "a" "t"
     [505] "c" "a" "a" "g" "a" "a" "g" "a" "t" "t" "t" "c" "c" "c" "c" "t" "a" "a"
     [523] "a" "g" "g" "a" "g" "t" "g" "g" "a" "c" "a" "t" "t" "g" "t" "c" "a" "t"
     [541] "g" "g" "a" "c" "c" "c" "t" "c" "t" "g" "g" "g" "t" "g" "g" "g" "t" "c"
     [559] "a" "g" "a" "t" "a" "c" "t" "g" "c" "c" "a" "a" "g" "g" "g" "c" "t" "a"
     [577] "c" "a" "a" "c" "c" "t" "c" "c" "t" "g" "a" "a" "a" "c" "c" "c" "a" "t"
     [595] "g" "g" "g" "c" "a" "a" "a" "g" "t" "c" "g" "t" "c" "a" "c" "c" "t" "a"
     [613] "t" "g" "g" "a" "a" "t" "g" "g" "c" "c" "a" "a" "c" "c" "t" "g" "c" "t"
     [631] "g" "a" "c" "g" "g" "g" "c" "c" "c" "c" "a" "a" "a" "c" "g" "g" "a" "a"
     [649] "c" "c" "t" "g" "a" "t" "g" "g" "c" "c" "c" "t" "g" "g" "c" "c" "c" "g"
     [667] "g" "a" "c" "a" "t" "g" "g" "t" "g" "g" "a" "a" "t" "c" "a" "g" "t" "t"
     [685] "c" "a" "g" "c" "g" "t" "g" "a" "c" "a" "g" "c" "t" "c" "t" "g" "c" "a"
     [703] "g" "c" "t" "g" "c" "t" "g" "c" "a" "g" "g" "c" "c" "a" "a" "c" "c" "g"
     [721] "g" "g" "c" "t" "g" "t" "g" "t" "g" "t" "g" "g" "c" "t" "t" "c" "c" "a"
     [739] "c" "c" "t" "g" "g" "g" "c" "t" "a" "c" "c" "t" "g" "g" "a" "t" "g" "g"
     [757] "t" "g" "a" "g" "g" "t" "g" "g" "a" "g" "c" "t" "g" "g" "t" "c" "a" "g"
     [775] "t" "g" "g" "t" "g" "t" "g" "g" "t" "g" "g" "c" "c" "c" "g" "c" "c" "t"
     [793] "c" "c" "t" "g" "g" "c" "t" "c" "t" "g" "t" "a" "c" "a" "a" "c" "c" "a"
     [811] "g" "g" "g" "c" "c" "a" "c" "a" "t" "c" "a" "a" "g" "c" "c" "c" "c" "a"
     [829] "c" "a" "t" "t" "g" "a" "c" "t" "c" "a" "g" "t" "c" "t" "g" "g" "c" "c"
     [847] "c" "t" "t" "c" "g" "a" "g" "a" "a" "g" "g" "t" "g" "g" "c" "t" "g" "a"
     [865] "t" "g" "c" "c" "a" "t" "g" "a" "a" "a" "c" "a" "g" "a" "t" "g" "c" "a"
     [883] "g" "g" "a" "g" "a" "a" "g" "a" "a" "g" "a" "a" "t" "g" "t" "g" "g" "g"
     [901] "c" "a" "a" "g" "g" "t" "c" "c" "t" "c" "c" "t" "g" "g" "t" "t" "c" "c"
     [919] "a" "g" "g" "g" "c" "c" "a" "g" "a" "g" "a" "a" "g" "g" "a" "g" "a" "g"
     [937] "c" "t" "a" "g" "g" "g" "c" "a" "a" "g" "t" "g" "g" "c" "t" "g" "t" "g"
     [955] "a" "g" "a" "c" "c" "c" "t" "a" "g" "a" "g" "a" "c" "c" "a" "g" "c" "g"
     [973] "a" "a" "g" "g" "g" "a" "g" "a" "a" "g" "t" "t" "g" "g" "g" "a" "a" "g"
     [991] "c" "t" "a" "c" "g" "t" "t" "c" "t" "g" "t" "t" "g" "g" "c" "c" "a" "c"
    [1009] "c" "a" "g" "a" "c" "t" "t" "g" "c" "a" "t" "t" "t" "c" "a" "g" "c" "c"
    [1027] "t" "c" "t" "g" "t" "c" "a" "t" "a" "a" "t" "g" "c" "t" "c" "t" "g" "c"
    [1045] "c" "c" "t" "c" "c" "c" "t" "c" "c" "c" "c" "c" "g" "a" "a" "g" "g" "t"
    [1063] "c" "t" "c" "t" "g" "t" "g" "g" "t" "g" "a" "t" "g" "a" "c" "c" "g" "c"
    [1081] "t" "c" "t" "t" "c" "c" "c" "t" "g" "c" "c" "c" "c" "t" "c" "c" "c" "c"
    [1099] "g" "c" "t" "t" "c" "c" "t" "g" "a" "c" "c" "t" "c" "t" "g" "a" "a" "g"
    [1117] "a" "g" "g" "t" "t" "g" "g" "g" "a" "a" "g" "t" "g" "a" "c" "c" "a" "t"
    [1135] "t" "t" "g" "g" "a" "t" "g" "t" "c" "t" "g" "g" "g" "c" "c" "c" "t" "g"
    [1153] "c" "c" "a" "a" "g" "g" "c" "g" "a" "c" "a" "g" "g" "g" "a" "g" "g" "g"
    [1171] "t" "c" "a" "g" "a" "g" "g" "g" "a" "g" "g" "c" "c" "g" "g" "c" "t" "g"
    [1189] "c" "t" "t" "c" "c" "t" "g" "c" "c" "c" "c" "c" "a" "c" "c" "c" "t" "t"
    [1207] "t" "c" "c" "c" "c" "g" "g" "g" "c" "c" "t" "g" "c" "t" "g" "t" "g" "c"
    [1225] "t" "g" "c" "t" "t" "t" "t" "g" "t" "g" "c" "c" "a" "a" "g" "g" "t" "t"
    [1243] "a" "g" "c" "c" "a" "g" "t" "c" "c" "c" "c" "c" "c" "t" "g" "t" "t" "g"
    [1261] "t" "g" "t" "t" "c" "c" "a" "t" "g" "t" "g" "c" "t" "t" "t" "c" "a" "c"
    [1279] "c" "t" "c" "t" "g" "c" "c" "t" "c" "a" "t" "c" "t" "t" "t" "c" "c" "t"
    [1297] "c" "c" "c" "g" "t" "c" "c" "c" "t" "g" "c" "c" "c" "c" "g" "c" "c" "a"
    [1315] "c" "c" "t" "c" "c" "c" "c" "a" "a" "a" "g" "a" "a" "t" "t" "g" "a" "a"
    [1333] "g" "a" "c" "g" "t" "c" "a" "g" "c" "t" "c" "a" "g" "g" "a" "t" "a" "t"
    [1351] "g" "g" "g" "g" "c" "c" "a" "t" "c" "t" "c" "t" "g" "t" "g" "a" "g" "t"
    [1369] "c" "a" "g" "c" "a" "t" "g" "t" "a" "c" "c" "t" "g" "t" "c" "t" "c" "t"
    [1387] "c" "c" "t" "a" "g" "t" "g" "t" "c" "c" "t" "t" "c" "a" "g" "c" "c" "t"
    [1405] "g" "g" "g" "c" "t" "g" "a" "c" "c" "a" "g" "t" "g" "c" "c" "g" "c" "c"
    [1423] "t" "c" "t" "g" "g" "g" "c" "t" "t" "g" "a" "c" "c" "a" "g" "t" "t" "c"
    [1441] "c" "a" "a" "t" "c" "t" "g" "t" "c" "t" "g" "t" "c" "c" "a" "a" "c" "t"
    [1459] "t" "c" "t" "t" "a" "a" "g" "c" "a" "c" "a" "a" "t" "t" "g" "g" "g" "c"
    [1477] "t" "t" "c" "t" "t" "c" "a" "t" "c" "t" "c" "c" "a" "g" "g" "t" "t" "t"
    [1495] "t" "c" "t" "g" "c" "c" "a" "t" "t" "c" "t" "t" "a" "a" "c" "c" "a" "a"
    [1513] "g" "g" "c" "t" "g" "c" "c" "t" "c" "t" "t" "c" "c" "a" "a" "c" "a" "g"
    [1531] "g" "g" "c" "g" "g" "g" "a" "a" "t" "c" "a" "g" "a" "c" "c" "t" "a" "c"
    [1549] "t" "c" "c" "c" "c" "t" "a" "g" "g" "t" "c" "a" "c" "a" "a" "c" "t" "c"
    [1567] "t" "g" "g" "g" "a" "a" "g" "g" "a" "t" "a" "c" "a" "g" "a" "g" "c" "c"
    [1585] "c" "c" "c" "a" "c" "c" "c" "t" "t" "c" "a" "c" "t" "g" "a" "g" "t" "t"
    [1603] "c" "t" "c" "t" "g" "g" "a" "t" "t" "t" "g" "t" "t" "c" "t" "c" "a" "g"
    [1621] "t" "g" "c" "c" "t" "t" "a" "g" "c" "a" "a" "c" "g" "a" "a" "a" "a" "c"
    [1639] "c" "t" "g" "t" "g" "c" "t" "t" "g" "t" "g" "t" "g" "t" "g" "t" "g" "t"
    [1657] "g" "g" "c" "g" "g" "c" "g" "g" "g" "g" "a" "g" "g" "g" "a" "g" "g" "a"
    [1675] "t" "c" "c" "t" "g" "t" "t" "t" "c" "c" "c" "a" "c" "c" "t" "c" "c" "t"
    [1693] "t" "c" "t" "c" "c" "t" "c" "c" "c" "c" "t" "g" "t" "a" "c" "t" "c" "c"
    [1711] "c" "c" "a" "g" "t" "g" "c" "c" "t" "t" "c" "c" "t" "t" "g" "t" "t" "c"
    [1729] "t" "g" "g" "t" "g" "g" "a" "g" "c" "t" "g" "g" "g" "g" "t" "t" "t" "c"
    [1747] "t" "c" "t" "c" "c" "t" "c" "c" "c" "c" "a" "g" "t" "c" "c" "a" "c" "a"
    [1765] "a" "c" "a" "c" "t" "g" "c" "c" "a" "a" "a" "a" "t" "c" "t" "g" "t" "g"
    [1783] "t" "a" "t" "g" "t" "g" "g" "c" "c" "a" "t" "t" "g" "g" "g" "t" "g" "g"
    [1801] "g" "g" "c" "a" "g" "c" "c" "c" "c" "a" "a" "g" "c" "c" "t" "c" "c" "t"
    [1819] "g" "g" "g" "g" "a" "g" "g" "c" "a" "g" "g" "g" "c" "a" "a" "a" "a" "a"
    [1837] "c" "a" "g" "g" "t" "g" "c" "c" "c" "t" "c" "a" "t" "c" "g" "t" "g" "g"
    [1855] "t" "c" "t" "g" "t" "g" "c" "c" "a" "t" "g" "t" "c" "c" "c" "g" "t" "c"
    [1873] "t" "c" "t" "a" "t" "g" "g" "t" "g" "g" "t" "t" "g" "a" "g" "g" "a" "g"
    [1891] "a" "a" "a" "g" "g" "c" "g" "g" "g" "g" "a" "a" "g" "c" "t" "t" "c" "c"
    [1909] "t" "c" "a" "g" "c" "c" "t" "t" "g" "c" "a" "g" "a" "t" "a" "t" "g" "t"
    [1927] "g" "t" "g" "g" "c" "a" "t" "t" "t" "a" "c" "t" "a" "g" "c" "c" "a" "g"
    [1945] "a" "g" "c" "t" "c" "c" "g" "a" "a" "a" "g" "g" "c" "a" "g" "t" "g" "c"
    [1963] "t" "g" "t" "c" "t" "g" "t" "t" "t" "c" "t" "t" "g" "t" "a" "c" "t" "g"
    [1981] "g" "g" "a" "c" "c" "a" "a" "a" "g" "t" "a" "a" "a" "a" "a" "t" "c" "c"
    [1999] "a" "a" "g" "c" "a" "c" "a" "t" "c" "c" "c" "c" "t" "t" "g" "c" "a" "c"
    [2017] "t" "t" "a" "g" "g" "g" "g" "a" "g" "g" "c" "c" "c" "t" "a" "c" "t" "g"
    [2035] "c" "c" "t" "c" "t" "c" "t" "c" "a" "a" "a" "g" "c" "a" "g" "a" "g" "a"
    [2053] "g" "g" "c" "a" "g" "c" "t" "t" "a" "t" "c" "a" "a" "a" "c" "c" "t" "c"
    [2071] "a" "g" "c" "c" "c" "a" "a" "a" "a" "c" "t" "c" "t" "g" "t" "t" "t" "a"
    [2089] "c" "a" "t" "g" "g" "g" "t" "g" "g" "g" "g" "a" "g" "a" "t" "g" "g" "a"
    [2107] "g" "c" "a" "g" "g" "g" "a" "a" "c" "t" "a" "c" "a" "g" "a" "g" "t" "g"
    [2125] "g" "g" "a" "t" "g" "g" "t" "c" "a" "g" "g" "a" "c" "c" "t" "g" "g" "g"
    [2143] "c" "c" "a" "t" "t" "g" "c" "a" "a" "c" "c" "a" "a" "a" "a" "t" "g" "g"
    [2161] "g" "g" "a" "c" "t" "t" "c" "c" "t" "g" "g" "g" "t" "a" "g" "g" "g" "a"
    [2179] "g" "g" "t" "c" "a" "c" "t" "c" "c" "c" "t" "c" "t" "a" "c" "t" "c" "a"
    [2197] "c" "t" "a" "g" "c" "t" "a" "g" "g" "a" "t" "t" "a" "g" "g" "g" "a" "g"
    [2215] "g" "g" "t" "t" "a" "t" "t" "g" "c" "c" "c" "c" "a" "c" "c" "a" "t" "t"
    [2233] "g" "c" "a" "a" "t" "g" "g" "g" "a" "g" "g" "t" "g" "g" "a" "g" "g" "g"
    [2251] "a" "c" "a" "g" "g" "c" "t" "c" "a" "g" "c" "c" "t" "c" "c" "t" "c" "a"
    [2269] "t" "t" "g" "t" "c" "t" "a" "a" "a" "t" "g" "a" "g" "g" "c" "c" "t" "a"
    [2287] "a" "a" "t" "g" "t" "g" "t" "g" "a" "a" "g" "t" "g" "c" "g" "a" "t" "t"
    [2305] "t" "c" "t" "g" "c" "t" "t" "t" "t" "g" "t" "g" "t" "a" "c" "c" "c" "c"
    [2323] "a" "c" "c" "a" "c" "c" "c" "c" "a" "t" "t" "a" "c" "c" "a" "c" "a" "g"
    [2341] "c" "t" "g" "c" "c" "t" "t" "t" "g" "t" "g" "t" "g" "t" "t" "t" "g" "g"
    [2359] "g" "t" "c" "a" "a" "t" "a" "a" "a" "a" "a" "g" "c" "c" "a" "a" "a" "c"
    [2377] "c" "c" "t" "g" "a"
    attr(,"name")
    [1] "HSU18009"
    attr(,"Annot")
    [1] ">HSU18009"
    attr(,"class")
    [1] "SeqFastadna"
    
    $HSU18018
       [1] "a" "c" "a" "a" "c" "t" "g" "t" "c" "t" "g" "c" "t" "g" "c" "g" "c" "c"
      [19] "c" "g" "a" "a" "a" "a" "a" "c" "a" "a" "g" "t" "c" "g" "g" "t" "g" "c"
      [37] "g" "c" "t" "g" "g" "g" "g" "a" "c" "c" "c" "g" "g" "g" "g" "c" "c" "g"
      [55] "g" "g" "g" "c" "c" "g" "c" "c" "t" "t" "a" "c" "t" "c" "c" "g" "g" "c"
      [73] "c" "t" "a" "g" "c" "c" "c" "c" "g" "c" "g" "g" "c" "c" "c" "t" "c" "g"
      [91] "g" "t" "g" "c" "g" "g" "g" "c" "t" "c" "c" "a" "g" "g" "g" "c" "a" "t"
     [109] "g" "c" "t" "c" "g" "g" "t" "a" "c" "c" "c" "c" "c" "c" "g" "c" "g" "g"
     [127] "c" "t" "c" "c" "a" "g" "c" "c" "c" "a" "g" "a" "c" "g" "c" "c" "c" "c"
     [145] "g" "g" "c" "c" "t" "c" "a" "g" "g" "t" "c" "t" "c" "g" "g" "c" "c" "c"
     [163] "c" "c" "g" "c" "t" "t" "g" "g" "g" "g" "c" "c" "c" "c" "g" "g" "c" "c"
     [181] "g" "t" "g" "c" "g" "g" "c" "g" "c" "g" "a" "g" "g" "g" "a" "g" "c" "g"
     [199] "g" "c" "c" "g" "g" "a" "t" "g" "g" "a" "g" "c" "g" "g" "a" "g" "g" "a"
     [217] "t" "g" "a" "a" "a" "g" "c" "c" "g" "g" "a" "t" "a" "c" "t" "t" "g" "g"
     [235] "a" "c" "c" "a" "g" "c" "a" "a" "g" "t" "g" "c" "c" "c" "t" "a" "c" "a"
     [253] "c" "c" "t" "t" "c" "a" "g" "c" "a" "g" "c" "a" "a" "a" "t" "c" "g" "c"
     [271] "c" "c" "g" "g" "a" "a" "a" "t" "g" "g" "g" "a" "g" "c" "t" "t" "g" "c"
     [289] "g" "c" "g" "a" "a" "g" "c" "g" "c" "t" "g" "a" "t" "c" "g" "g" "c" "c"
     [307] "c" "g" "c" "t" "g" "g" "g" "g" "a" "a" "g" "c" "t" "c" "a" "t" "g" "g"
     [325] "a" "c" "c" "c" "g" "g" "g" "c" "t" "c" "c" "c" "t" "g" "c" "c" "g" "c"
     [343] "c" "c" "c" "t" "c" "g" "a" "c" "t" "c" "t" "g" "a" "a" "g" "a" "t" "c"
     [361] "t" "c" "t" "t" "c" "c" "a" "g" "g" "a" "t" "c" "t" "a" "a" "g" "t" "c"
     [379] "a" "c" "t" "t" "c" "c" "a" "g" "g" "a" "g" "a" "c" "g" "t" "g" "g" "c"
     [397] "t" "c" "g" "c" "t" "g" "a" "a" "g" "c" "t" "c" "a" "g" "g" "t" "a" "c"
     [415] "c" "a" "g" "a" "c" "a" "g" "t" "g" "a" "t" "g" "a" "g" "c" "a" "g" "t"
     [433] "t" "t" "g" "t" "t" "c" "c" "t" "g" "a" "t" "t" "t" "c" "c" "a" "t" "t"
     [451] "c" "a" "g" "a" "a" "a" "a" "c" "c" "t" "a" "g" "c" "t" "t" "t" "c" "c"
     [469] "a" "c" "a" "g" "c" "c" "c" "c" "a" "c" "c" "a" "c" "c" "a" "g" "g" "a"
     [487] "t" "c" "a" "a" "g" "a" "a" "g" "g" "a" "g" "c" "c" "c" "c" "a" "g" "a"
     [505] "g" "t" "c" "c" "c" "c" "g" "c" "a" "c" "a" "g" "a" "c" "c" "c" "g" "g"
     [523] "c" "c" "c" "t" "g" "t" "c" "c" "t" "g" "c" "a" "g" "c" "a" "g" "g" "a"
     [541] "a" "g" "c" "c" "g" "c" "c" "a" "c" "t" "c" "c" "c" "c" "t" "a" "c" "c"
     [559] "a" "c" "c" "a" "t" "g" "g" "c" "g" "a" "g" "c" "a" "g" "t" "g" "c" "c"
     [577] "t" "t" "t" "a" "c" "t" "c" "c" "a" "g" "t" "g" "c" "c" "t" "a" "t" "g"
     [595] "a" "c" "c" "c" "c" "c" "c" "c" "a" "g" "a" "c" "a" "a" "a" "t" "c" "g"
     [613] "c" "c" "a" "t" "c" "a" "a" "g" "t" "c" "c" "c" "c" "t" "g" "c" "c" "c"
     [631] "c" "t" "g" "g" "t" "g" "c" "c" "c" "t" "t" "g" "g" "a" "c" "a" "g" "t"
     [649] "c" "g" "c" "c" "c" "c" "t" "a" "c" "a" "g" "c" "c" "c" "t" "t" "t" "c"
     [667] "c" "c" "c" "g" "g" "g" "c" "a" "g" "a" "g" "c" "a" "a" "c" "g" "g" "a"
     [685] "a" "t" "t" "t" "c" "c" "t" "g" "a" "g" "a" "t" "c" "c" "t" "c" "t" "g"
     [703] "g" "c" "a" "c" "c" "t" "c" "c" "c" "a" "g" "c" "c" "c" "c" "a" "c" "c"
     [721] "c" "t" "g" "g" "c" "c" "a" "t" "g" "g" "g" "t" "a" "c" "c" "t" "c" "g"
     [739] "g" "g" "g" "a" "a" "c" "a" "t" "a" "g" "c" "t" "c" "c" "g" "t" "c" "t"
     [757] "t" "c" "c" "a" "g" "c" "a" "g" "c" "c" "c" "c" "t" "g" "g" "a" "c" "a"
     [775] "t" "t" "t" "g" "c" "c" "a" "c" "t" "c" "c" "t" "t" "c" "a" "c" "a" "t"
     [793] "c" "t" "c" "a" "g" "g" "g" "a" "g" "g" "g" "g" "g" "c" "c" "g" "g" "g"
     [811] "a" "a" "c" "c" "c" "c" "t" "c" "c" "c" "a" "g" "c" "c" "c" "c" "c" "t"
     [829] "a" "c" "c" "a" "a" "c" "a" "c" "c" "a" "g" "c" "t" "g" "t" "c" "g" "g"
     [847] "a" "g" "c" "c" "c" "t" "g" "c" "c" "c" "a" "c" "c" "c" "t" "a" "t" "c"
     [865] "c" "c" "c" "a" "g" "c" "a" "g" "a" "g" "c" "t" "t" "t" "a" "a" "g" "c"
     [883] "a" "a" "g" "a" "a" "t" "a" "c" "c" "a" "t" "g" "a" "t" "c" "c" "c" "c"
     [901] "t" "g" "t" "a" "t" "g" "a" "a" "c" "a" "g" "g" "c" "g" "g" "g" "c" "c"
     [919] "a" "g" "c" "c" "a" "g" "c" "c" "g" "t" "g" "g" "a" "c" "c" "a" "g" "g"
     [937] "g" "t" "g" "g" "g" "g" "t" "c" "a" "a" "t" "g" "g" "g" "c" "a" "c" "a"
     [955] "g" "g" "t" "a" "c" "c" "c" "a" "g" "g" "g" "g" "c" "g" "g" "g" "g" "g"
     [973] "t" "g" "g" "t" "g" "a" "t" "c" "a" "a" "a" "c" "a" "g" "g" "a" "a" "c"
     [991] "a" "g" "a" "c" "g" "g" "a" "c" "t" "t" "c" "g" "c" "c" "t" "a" "c" "g"
    [1009] "a" "c" "t" "c" "a" "g" "a" "t" "g" "t" "c" "a" "c" "c" "g" "g" "g" "t"
    [1027] "g" "c" "g" "c" "a" "t" "c" "a" "a" "t" "g" "t" "a" "c" "c" "t" "c" "c"
    [1045] "a" "c" "a" "c" "a" "g" "a" "g" "g" "g" "c" "t" "t" "c" "t" "c" "t" "g"
    [1063] "g" "g" "c" "c" "c" "t" "c" "t" "c" "c" "a" "g" "g" "t" "g" "a" "c" "g"
    [1081] "g" "g" "g" "c" "c" "a" "t" "g" "g" "g" "c" "t" "a" "t" "g" "g" "c" "t"
    [1099] "a" "t" "g" "a" "g" "a" "a" "a" "c" "c" "t" "c" "t" "g" "c" "g" "a" "c"
    [1117] "c" "a" "t" "t" "c" "c" "c" "a" "g" "a" "t" "g" "a" "t" "g" "t" "c" "t"
    [1135] "g" "c" "g" "t" "t" "g" "t" "c" "c" "c" "t" "g" "a" "g" "a" "a" "a" "t"
    [1153] "t" "t" "g" "a" "a" "g" "g" "a" "g" "a" "c" "a" "t" "c" "a" "a" "g" "c"
    [1171] "a" "g" "g" "a" "a" "g" "g" "g" "g" "t" "c" "g" "g" "t" "g" "c" "a" "t"
    [1189] "t" "t" "c" "g" "a" "g" "a" "g" "g" "g" "g" "c" "c" "g" "c" "c" "c" "t"
    [1207] "a" "c" "c" "a" "g" "c" "g" "c" "c" "g" "g" "g" "g" "t" "g" "c" "c" "c"
    [1225] "t" "g" "c" "a" "g" "c" "t" "g" "t" "g" "g" "c" "a" "a" "t" "t" "t" "c"
    [1243] "t" "g" "g" "t" "g" "g" "c" "c" "t" "t" "g" "c" "t" "g" "g" "a" "t" "g"
    [1261] "a" "c" "c" "c" "a" "a" "c" "a" "a" "a" "t" "g" "c" "c" "c" "a" "t" "t"
    [1279] "t" "c" "a" "t" "t" "g" "c" "c" "t" "g" "g" "a" "c" "g" "g" "g" "c" "c"
    [1297] "g" "g" "g" "g" "a" "a" "t" "g" "g" "a" "g" "t" "t" "c" "a" "a" "g" "c"
    [1315] "t" "c" "a" "t" "t" "g" "a" "g" "c" "c" "t" "g" "a" "g" "g" "a" "g" "g"
    [1333] "t" "c" "g" "c" "c" "a" "g" "g" "c" "t" "c" "t" "g" "g" "g" "g" "c" "a"
    [1351] "t" "c" "c" "a" "g" "a" "a" "g" "a" "a" "c" "c" "g" "g" "c" "c" "a" "g"
    [1369] "c" "c" "a" "t" "g" "a" "a" "t" "t" "a" "c" "g" "a" "c" "a" "a" "g" "c"
    [1387] "t" "g" "a" "g" "c" "c" "g" "c" "t" "c" "g" "c" "t" "c" "c" "g" "a" "t"
    [1405] "a" "c" "t" "a" "t" "t" "a" "t" "g" "a" "g" "a" "a" "a" "g" "g" "c" "a"
    [1423] "t" "c" "a" "t" "g" "c" "a" "g" "a" "a" "g" "g" "t" "g" "g" "c" "t" "g"
    [1441] "g" "t" "g" "a" "g" "c" "g" "t" "t" "a" "c" "g" "t" "g" "t" "a" "c" "a"
    [1459] "a" "g" "t" "t" "t" "g" "t" "g" "t" "g" "t" "g" "a" "g" "c" "c" "c" "g"
    [1477] "a" "g" "g" "c" "c" "c" "t" "c" "t" "t" "c" "t" "c" "t" "t" "t" "g" "g"
    [1495] "c" "c" "t" "t" "c" "c" "c" "g" "g" "a" "c" "a" "a" "t" "c" "a" "g" "c"
    [1513] "g" "t" "c" "c" "a" "g" "c" "t" "c" "t" "c" "a" "a" "g" "g" "c" "t" "g"
    [1531] "a" "g" "t" "t" "t" "g" "a" "c" "c" "g" "g" "c" "c" "t" "g" "t" "c" "a"
    [1549] "g" "t" "g" "a" "g" "g" "a" "g" "g" "a" "c" "a" "c" "a" "g" "t" "c" "c"
    [1567] "c" "t" "t" "t" "g" "t" "c" "c" "c" "a" "c" "t" "t" "g" "g" "a" "t" "g"
    [1585] "a" "g" "a" "g" "c" "c" "c" "c" "g" "c" "c" "t" "a" "c" "c" "t" "c" "c"
    [1603] "c" "a" "g" "a" "g" "c" "t" "g" "g" "c" "t" "g" "g" "c" "c" "c" "c" "g"
    [1621] "c" "c" "c" "a" "g" "c" "c" "a" "t" "t" "t" "g" "g" "c" "c" "c" "c" "a"
    [1639] "a" "g" "g" "g" "t" "g" "g" "c" "t" "a" "c" "t" "c" "t" "t" "a" "c" "t"
    [1657] "a" "g" "c" "c" "c" "c" "c" "a" "g" "c" "g" "g" "c" "t" "g" "t" "t" "c"
    [1675] "c" "c" "c" "c" "t" "g" "c" "c" "g" "c" "a" "g" "g" "t" "g" "g" "g" "t"
    [1693] "g" "c" "t" "g" "c" "c" "c" "t" "g" "t" "g" "t" "a" "c" "a" "t" "a" "t"
    [1711] "a" "a" "a" "t" "g" "a" "a" "t" "c" "t" "g" "g" "t" "g" "t" "t" "g" "g"
    [1729] "g" "g" "a" "a" "a" "c" "c" "t" "t" "c" "a" "t" "c" "t" "g" "a" "a" "a"
    [1747] "c" "c" "c" "a" "c" "a" "g" "a" "t" "g" "t" "c" "t" "c" "t" "g" "g" "g"
    [1765] "g" "c" "a" "g" "a" "t" "c" "c" "c" "c" "a" "c" "t" "g" "t" "c" "c" "t"
    [1783] "a" "c" "c" "a" "g" "t" "t" "g" "c" "c" "c" "t" "a" "g" "c" "c" "c" "a"
    [1801] "g" "a" "c" "t" "c" "t" "g" "a" "g" "c" "t" "g" "c" "t" "c" "a" "c" "c"
    [1819] "g" "g" "a" "g" "t" "c" "a" "t" "t" "g" "g" "g" "a" "a" "g" "g" "a" "a"
    [1837] "a" "a" "g" "t" "g" "g" "a" "g" "a" "a" "a" "t" "g" "g" "c" "a" "a" "g"
    [1855] "t" "c" "t" "a" "g" "a" "g" "t" "c" "t" "c" "a" "g" "a" "a" "a" "c" "t"
    [1873] "c" "c" "c" "c" "t" "g" "g" "g" "g" "g" "t" "t" "t" "c" "a" "c" "c" "t"
    [1891] "g" "g" "g" "c" "c" "c" "t" "g" "g" "a" "g" "g" "a" "a" "t" "t" "c" "a"
    [1909] "g" "c" "t" "c" "a" "g" "c" "t" "t" "c" "t" "t" "c" "c" "t" "a" "g" "g"
    [1927] "t" "c" "c" "a" "a" "g" "c" "c" "c" "c" "c" "c" "a" "c" "a" "c" "c" "t"
    [1945] "t" "t" "t" "c" "c" "c" "c" "a" "a" "c" "c" "a" "c" "a" "g" "a" "g" "a"
    [1963] "a" "c" "a" "a" "g" "a" "g" "t" "t" "t" "g" "t" "t" "c" "t" "g" "t" "t"
    [1981] "c" "t" "g" "g" "g" "g" "g" "a" "c" "a" "g" "a" "g" "a" "a" "g" "g" "c"
    [1999] "g" "c" "t" "t" "c" "c" "c" "a" "a" "c" "t" "t" "c" "a" "t" "a" "c" "t"
    [2017] "g" "g" "c" "a" "g" "g" "a" "g" "g" "g" "t" "g" "a" "g" "g" "a" "g" "g"
    [2035] "t" "t" "c" "a" "c" "t" "g" "a" "g" "c" "t" "c" "c" "c" "c" "a" "g" "a"
    [2053] "t" "c" "t" "c" "c" "c" "a" "c" "t" "g" "c" "g" "g" "g" "g" "a" "g" "a"
    [2071] "c" "a" "g" "a" "a" "g" "c" "c" "t" "g" "g" "a" "c" "t" "c" "t" "g" "c"
    [2089] "c" "c" "c" "a" "c" "g" "c" "t" "g" "t" "g" "g" "c" "c" "c" "t" "g" "g"
    [2107] "a" "g" "g" "g" "t" "c" "c" "c" "g" "g" "t" "t" "t" "g" "t" "c" "a" "g"
    [2125] "t" "t" "c" "t" "t" "g" "g" "t" "g" "c" "t" "c" "t" "g" "t" "g" "t" "t"
    [2143] "c" "c" "c" "a" "g" "a" "g" "g" "c" "a" "g" "g" "c" "g" "g" "a" "g" "g"
    [2161] "t" "t" "g" "a" "a" "g" "a" "a" "a" "g" "g" "a" "a" "c" "c" "t" "g" "g"
    [2179] "g" "a" "t" "g" "a" "g" "g" "g" "g" "t" "g" "c" "t" "g" "g" "g" "t" "a"
    [2197] "t" "a" "a" "g" "c" "a" "g" "a" "g" "a" "g" "g" "g" "a" "t" "g" "g" "g"
    [2215] "t" "t" "c" "c" "t" "g" "c" "t" "c" "c" "a" "a" "g" "g" "g" "a" "c" "c"
    [2233] "c" "t" "t" "t" "g" "c" "c" "t" "t" "t" "c" "t" "t" "c" "t" "g" "c" "c"
    [2251] "c" "t" "t" "t" "c" "c" "t" "a" "g" "g" "c" "c" "c" "a" "g" "g" "c" "c"
    [2269] "t" "g" "g" "g" "t" "t" "t" "g" "t" "a" "c" "t" "t" "c" "c" "a" "c" "c"
    [2287] "t" "c" "c" "a" "c" "c" "a" "c" "a" "t" "c" "t" "g" "c" "c" "a" "g" "a"
    [2305] "c" "c" "t" "t" "a" "a" "t" "a" "a" "a" "g" "g" "c" "c" "c" "c" "c" "a"
    [2323] "c" "t" "t" "c" "t" "c" "c" "c" "a" "t" "t"
    attr(,"name")
    [1] "HSU18018"
    attr(,"Annot")
    [1] ">HSU18018"
    attr(,"class")
    [1] "SeqFastadna"
    
    $HSU37574.BRCA1
     [1] "a" "t" "g" "g" "a" "t" "t" "t" "a" "t" "c" "t" "g" "c" "t" "c" "t" "t" "c"
    [20] "g" "c" "g" "t" "t" "g" "a" "a" "g" "a" "a" "g" "t" "a" "c" "a" "a" "a" "a"
    [39] "t" "g" "t" "c" "a" "t" "t" "a" "a" "t" "g" "c" "t" "a" "t" "g" "c" "a" "g"
    [58] "a" "a" "a" "a" "t" "c" "t" "t" "a" "g" "a" "g" "t" "g" "t" "c" "c" "c" "a"
    [77] "t" "c" "t" "g"
    attr(,"name")
    [1] "HSU37574.BRCA1"
    attr(,"Annot")
    [1] ">HSU37574.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $HSU61268.BRCA1
     [1] "a" "t" "g" "g" "a" "t" "t" "t" "a" "t" "c" "t" "g" "c" "t" "c" "t" "t" "c"
    [20] "g" "c" "g" "t" "t" "g" "a" "a" "g" "a" "a" "g" "t" "a" "c" "a" "a" "a" "a"
    [39] "t" "g" "t" "c" "a" "t" "t" "a" "a" "t" "g" "c" "t" "a" "t" "g" "c" "a" "g"
    [58] "a" "a" "a" "a" "t" "c" "t" "t" "a" "g" "a" "g" "t" "g" "t" "c" "c" "c" "a"
    [77] "t" "c" "t" "g"
    attr(,"name")
    [1] "HSU61268.BRCA1"
    attr(,"Annot")
    [1] ">HSU61268.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $HSU64805
       [1] "a" "t" "g" "g" "a" "t" "t" "t" "a" "t" "c" "t" "g" "c" "t" "c" "t" "t"
      [19] "c" "g" "c" "g" "t" "t" "g" "a" "a" "g" "a" "a" "g" "t" "a" "c" "a" "a"
      [37] "a" "a" "t" "g" "t" "c" "a" "t" "t" "a" "a" "t" "g" "c" "t" "a" "t" "g"
      [55] "c" "a" "g" "a" "a" "a" "a" "t" "c" "t" "t" "a" "g" "a" "g" "t" "g" "t"
      [73] "c" "c" "c" "a" "t" "c" "t" "g" "t" "c" "t" "g" "g" "a" "g" "t" "t" "g"
      [91] "a" "t" "c" "a" "a" "g" "g" "a" "a" "c" "c" "t" "g" "t" "c" "t" "c" "c"
     [109] "a" "c" "a" "a" "a" "g" "t" "g" "t" "g" "a" "c" "c" "a" "c" "a" "t" "a"
     [127] "t" "t" "t" "t" "g" "c" "a" "a" "a" "t" "t" "t" "t" "g" "c" "a" "t" "g"
     [145] "c" "t" "g" "a" "a" "a" "c" "t" "t" "c" "t" "c" "a" "a" "c" "c" "a" "g"
     [163] "a" "a" "g" "a" "a" "a" "g" "g" "g" "c" "c" "t" "t" "c" "a" "c" "a" "g"
     [181] "t" "g" "t" "c" "c" "t" "t" "t" "a" "t" "g" "t" "a" "a" "g" "a" "a" "t"
     [199] "g" "a" "t" "a" "t" "a" "a" "c" "c" "a" "a" "a" "a" "g" "g" "a" "g" "c"
     [217] "c" "t" "a" "c" "a" "a" "g" "a" "a" "a" "g" "t" "a" "c" "g" "a" "g" "a"
     [235] "t" "t" "t" "a" "g" "t" "c" "a" "a" "c" "t" "t" "g" "t" "t" "g" "a" "a"
     [253] "g" "a" "g" "c" "t" "a" "t" "t" "g" "a" "a" "a" "a" "t" "c" "a" "t" "t"
     [271] "t" "g" "t" "g" "c" "t" "t" "t" "t" "c" "a" "g" "c" "t" "t" "g" "a" "c"
     [289] "a" "c" "a" "g" "g" "t" "t" "t" "g" "g" "a" "g" "t" "a" "t" "g" "c" "a"
     [307] "a" "a" "c" "a" "g" "c" "t" "a" "t" "a" "a" "t" "t" "t" "t" "g" "c" "a"
     [325] "a" "a" "a" "a" "a" "g" "g" "a" "a" "a" "a" "t" "a" "a" "c" "t" "c" "t"
     [343] "c" "c" "t" "g" "a" "a" "c" "a" "t" "c" "t" "a" "a" "a" "a" "g" "a" "t"
     [361] "g" "a" "a" "g" "t" "t" "t" "c" "t" "a" "t" "c" "a" "t" "c" "c" "a" "a"
     [379] "a" "g" "t" "a" "t" "g" "g" "g" "c" "t" "a" "c" "a" "g" "a" "a" "a" "c"
     [397] "c" "g" "t" "g" "c" "c" "a" "a" "a" "a" "g" "a" "c" "t" "t" "c" "t" "a"
     [415] "c" "a" "g" "a" "g" "t" "g" "a" "a" "c" "c" "c" "g" "a" "a" "a" "a" "t"
     [433] "c" "c" "t" "t" "c" "c" "t" "t" "g" "c" "a" "g" "g" "a" "a" "a" "c" "c"
     [451] "a" "g" "t" "c" "t" "c" "a" "g" "t" "g" "t" "c" "c" "a" "a" "c" "t" "c"
     [469] "t" "c" "t" "a" "a" "c" "c" "t" "t" "g" "g" "a" "a" "c" "t" "g" "t" "g"
     [487] "a" "g" "a" "a" "c" "t" "c" "t" "g" "a" "g" "g" "a" "c" "a" "a" "a" "g"
     [505] "c" "a" "g" "c" "g" "g" "a" "t" "a" "c" "a" "a" "c" "c" "t" "c" "a" "a"
     [523] "a" "a" "g" "a" "c" "g" "t" "c" "t" "g" "t" "c" "t" "a" "c" "a" "t" "t"
     [541] "g" "a" "a" "t" "t" "g" "g" "g" "a" "t" "c" "t" "g" "a" "t" "t" "c" "t"
     [559] "t" "c" "t" "g" "a" "a" "g" "a" "t" "a" "c" "c" "g" "t" "t" "a" "a" "t"
     [577] "a" "a" "g" "g" "c" "a" "a" "c" "t" "t" "a" "t" "t" "g" "c" "a" "g" "t"
     [595] "g" "t" "g" "g" "g" "a" "g" "a" "t" "c" "a" "a" "g" "a" "a" "t" "t" "g"
     [613] "t" "t" "a" "c" "a" "a" "a" "t" "c" "a" "c" "c" "c" "c" "t" "c" "a" "a"
     [631] "g" "g" "a" "a" "c" "c" "a" "g" "g" "g" "a" "t" "g" "a" "a" "a" "t" "c"
     [649] "a" "g" "t" "t" "t" "g" "g" "a" "t" "t" "c" "t" "g" "c" "a" "a" "a" "a"
     [667] "a" "a" "g" "g" "c" "t" "g" "c" "t" "t" "g" "t" "g" "a" "a" "t" "t" "t"
     [685] "t" "c" "t" "g" "a" "g" "a" "c" "g" "g" "a" "t" "g" "t" "a" "a" "c" "a"
     [703] "a" "a" "t" "a" "c" "t" "g" "a" "a" "c" "a" "t" "c" "g" "t" "c" "a" "a"
     [721] "c" "c" "c" "a" "g" "t" "a" "a" "t" "a" "a" "t" "g" "a" "t" "t" "t" "g"
     [739] "a" "a" "c" "a" "c" "c" "a" "c" "t" "g" "a" "g" "a" "a" "g" "c" "g" "t"
     [757] "g" "t" "a" "g" "c" "t" "g" "a" "g" "a" "g" "g" "c" "a" "t" "c" "c" "a"
     [775] "g" "a" "a" "a" "a" "g" "t" "a" "t" "c" "a" "g" "g" "g" "t" "g" "a" "a"
     [793] "g" "c" "a" "g" "c" "a" "t" "c" "t" "g" "g" "g" "t" "g" "t" "g" "a" "g"
     [811] "a" "g" "t" "g" "a" "a" "a" "c" "a" "a" "g" "c" "g" "t" "c" "t" "c" "t"
     [829] "g" "a" "a" "g" "a" "c" "t" "g" "c" "t" "c" "a" "g" "g" "g" "c" "t" "a"
     [847] "t" "c" "c" "t" "c" "t" "c" "a" "g" "a" "g" "t" "g" "a" "c" "a" "t" "t"
     [865] "t" "t" "a" "a" "c" "c" "a" "c" "t" "c" "a" "g" "c" "a" "g" "a" "g" "g"
     [883] "g" "a" "t" "a" "c" "c" "a" "t" "g" "c" "a" "a" "c" "a" "t" "a" "a" "c"
     [901] "c" "t" "g" "a" "t" "a" "a" "a" "g" "c" "t" "c" "c" "a" "g" "c" "a" "g"
     [919] "g" "a" "a" "a" "t" "g" "g" "c" "t" "g" "a" "a" "c" "t" "a" "g" "a" "a"
     [937] "g" "c" "t" "g" "t" "g" "t" "t" "a" "g" "a" "a" "c" "a" "g" "c" "a" "t"
     [955] "g" "g" "g" "a" "g" "c" "c" "a" "g" "c" "c" "t" "c" "c" "t" "a" "a" "c"
     [973] "a" "g" "c" "t" "a" "c" "c" "c" "t" "t" "c" "c" "a" "t" "c" "a" "t" "a"
     [991] "a" "g" "t" "g" "a" "c" "t" "c" "c" "t" "c" "t" "g" "c" "c" "c" "t" "t"
    [1009] "g" "a" "g" "g" "a" "c" "c" "t" "g" "c" "g" "a" "a" "a" "t" "c" "c" "a"
    [1027] "g" "a" "a" "c" "a" "a" "a" "g" "c" "a" "c" "a" "t" "c" "a" "g" "a" "a"
    [1045] "a" "a" "a" "g" "t" "a" "t" "t" "a" "a" "c" "t" "t" "c" "a" "c" "a" "g"
    [1063] "a" "a" "a" "a" "g" "t" "a" "g" "t" "g" "a" "a" "t" "a" "c" "c" "c" "t"
    [1081] "a" "t" "a" "a" "g" "c" "c" "a" "g" "a" "a" "t" "c" "c" "a" "g" "a" "a"
    [1099] "g" "g" "c" "c" "t" "t" "t" "c" "t" "g" "c" "t" "g" "a" "c" "a" "a" "g"
    [1117] "t" "t" "t" "g" "a" "g" "g" "t" "g" "t" "c" "t" "g" "c" "a" "g" "a" "t"
    [1135] "a" "g" "t" "t" "c" "t" "a" "c" "c" "a" "g" "t" "a" "a" "a" "a" "a" "t"
    [1153] "a" "a" "a" "g" "a" "a" "c" "c" "a" "g" "g" "a" "g" "t" "g" "g" "a" "a"
    [1171] "a" "g" "g" "t" "c" "a" "t" "c" "c" "c" "c" "t" "t" "c" "t" "a" "a" "a"
    [1189] "t" "g" "c" "c" "c" "a" "t" "c" "a" "t" "t" "a" "g" "a" "t" "g" "a" "t"
    [1207] "a" "g" "g" "t" "g" "g" "t" "a" "c" "a" "t" "g" "c" "a" "c" "a" "g" "t"
    [1225] "t" "g" "c" "t" "c" "t" "g" "g" "g" "a" "g" "t" "c" "t" "t" "c" "a" "g"
    [1243] "a" "a" "t" "a" "g" "a" "a" "a" "c" "t" "a" "c" "c" "c" "a" "t" "c" "t"
    [1261] "c" "a" "a" "g" "a" "g" "g" "a" "g" "c" "t" "c" "a" "t" "t" "a" "a" "g"
    [1279] "g" "t" "t" "g" "t" "t" "g" "a" "t" "g" "t" "g" "g" "a" "g" "g" "a" "g"
    [1297] "c" "a" "a" "c" "a" "g" "c" "t" "g" "g" "a" "a" "g" "a" "g" "t" "c" "t"
    [1315] "g" "g" "g" "c" "c" "a" "c" "a" "c" "g" "a" "t" "t" "t" "g" "a" "c" "g"
    [1333] "g" "a" "a" "a" "c" "a" "t" "c" "t" "t" "a" "c" "t" "t" "g" "c" "c" "a"
    [1351] "a" "g" "g" "c" "a" "a" "g" "a" "t" "c" "t" "a" "g" "a" "g" "g" "g" "a"
    [1369] "a" "c" "c" "c" "c" "t" "t" "a" "c" "c" "t" "g" "g" "a" "a" "t" "c" "t"
    [1387] "g" "g" "a" "a" "t" "c" "a" "g" "c" "c" "t" "c" "t" "t" "c" "t" "c" "t"
    [1405] "g" "a" "t" "g" "a" "c" "c" "c" "t" "g" "a" "a" "t" "c" "t" "g" "a" "t"
    [1423] "c" "c" "t" "t" "c" "t" "g" "a" "a" "g" "a" "c" "a" "g" "a" "g" "c" "c"
    [1441] "c" "c" "a" "g" "a" "g" "t" "c" "a" "g" "c" "t" "c" "g" "t" "g" "t" "t"
    [1459] "g" "g" "c" "a" "a" "c" "a" "t" "a" "c" "c" "a" "t" "c" "t" "t" "c" "a"
    [1477] "a" "c" "c" "t" "c" "t" "g" "c" "a" "t" "t" "g" "a" "a" "a" "g" "t" "t"
    [1495] "c" "c" "c" "c" "a" "a" "t" "t" "g" "a" "a" "a" "g" "t" "t" "g" "c" "a"
    [1513] "g" "a" "a" "t" "c" "t" "g" "c" "c" "c" "a" "g" "g" "g" "t" "c" "c" "a"
    [1531] "g" "c" "t" "g" "c" "t" "g" "c" "t" "c" "a" "t" "a" "c" "t" "a" "c" "t"
    [1549] "g" "a" "t" "a" "c" "t" "g" "c" "t" "g" "g" "g" "t" "a" "t" "a" "a" "t"
    [1567] "g" "c" "a" "a" "t" "g" "g" "a" "a" "g" "a" "a" "a" "g" "t" "g" "t" "g"
    [1585] "a" "g" "c" "a" "g" "g" "g" "a" "g" "a" "a" "g" "c" "c" "a" "g" "a" "a"
    [1603] "t" "t" "g" "a" "c" "a" "g" "c" "t" "t" "c" "a" "a" "c" "a" "g" "a" "a"
    [1621] "a" "g" "g" "g" "t" "c" "a" "a" "c" "a" "a" "a" "a" "g" "a" "a" "t" "g"
    [1639] "t" "c" "c" "a" "t" "g" "g" "t" "g" "g" "t" "g" "t" "c" "t" "g" "g" "c"
    [1657] "c" "t" "g" "a" "c" "c" "c" "c" "a" "g" "a" "a" "g" "a" "a" "t" "t" "t"
    [1675] "a" "t" "g" "c" "t" "c" "g" "t" "g" "t" "a" "c" "a" "a" "g" "t" "t" "t"
    [1693] "g" "c" "c" "a" "g" "a" "a" "a" "a" "c" "a" "c" "c" "a" "c" "a" "t" "c"
    [1711] "a" "c" "t" "t" "t" "a" "a" "c" "t" "a" "a" "t" "c" "t" "a" "a" "t" "t"
    [1729] "a" "c" "t" "g" "a" "a" "g" "a" "g" "a" "c" "t" "a" "c" "t" "c" "a" "t"
    [1747] "g" "t" "t" "g" "t" "t" "a" "t" "g" "a" "a" "a" "a" "c" "a" "g" "a" "t"
    [1765] "g" "c" "t" "g" "a" "g" "t" "t" "t" "g" "t" "g" "t" "g" "t" "g" "a" "a"
    [1783] "c" "g" "g" "a" "c" "a" "c" "t" "g" "a" "a" "a" "t" "a" "t" "t" "t" "t"
    [1801] "c" "t" "a" "g" "g" "a" "a" "t" "t" "g" "c" "g" "g" "g" "a" "g" "g" "a"
    [1819] "a" "a" "a" "t" "g" "g" "g" "t" "a" "g" "t" "t" "a" "g" "c" "t" "a" "t"
    [1837] "t" "t" "c" "t" "g" "g" "g" "t" "g" "a" "c" "c" "c" "a" "g" "t" "c" "t"
    [1855] "a" "t" "t" "a" "a" "a" "g" "a" "a" "a" "g" "a" "a" "a" "a" "a" "t" "g"
    [1873] "c" "t" "g" "a" "a" "t" "g" "a" "g" "c" "a" "t" "g" "a" "t" "t" "t" "t"
    [1891] "g" "a" "a" "g" "t" "c" "a" "g" "a" "g" "g" "a" "g" "a" "t" "g" "t" "g"
    [1909] "g" "t" "c" "a" "a" "t" "g" "g" "a" "a" "g" "a" "a" "a" "c" "c" "a" "c"
    [1927] "c" "a" "a" "g" "g" "t" "c" "c" "a" "a" "a" "g" "c" "g" "a" "g" "c" "a"
    [1945] "a" "g" "a" "g" "a" "a" "t" "c" "c" "c" "a" "g" "g" "a" "c" "a" "g" "a"
    [1963] "a" "a" "g" "a" "t" "c" "t" "t" "c" "a" "g" "g" "g" "g" "g" "c" "t" "a"
    [1981] "g" "a" "a" "a" "t" "c" "t" "g" "t" "t" "g" "c" "t" "a" "t" "g" "g" "g"
    [1999] "c" "c" "c" "t" "t" "c" "a" "c" "c" "a" "a" "c" "a" "t" "g" "c" "c" "c"
    [2017] "a" "c" "a" "g" "a" "t" "c" "a" "a" "c" "t" "g" "g" "a" "a" "t" "g" "g"
    [2035] "a" "t" "g" "g" "t" "a" "c" "a" "g" "c" "t" "g" "t" "g" "t" "g" "g" "t"
    [2053] "g" "c" "t" "t" "c" "t" "g" "t" "g" "g" "t" "g" "a" "a" "g" "g" "a" "g"
    [2071] "c" "t" "t" "t" "c" "a" "t" "c" "a" "t" "t" "c" "a" "c" "c" "c" "t" "t"
    [2089] "g" "g" "c" "a" "c" "a" "g" "g" "t" "g" "t" "c" "c" "a" "c" "c" "c" "a"
    [2107] "a" "t" "t" "g" "t" "g" "g" "t" "t" "g" "t" "g" "c" "a" "g" "c" "c" "a"
    [2125] "g" "a" "t" "g" "c" "c" "t" "g" "g" "a" "c" "a" "g" "a" "g" "g" "a" "c"
    [2143] "a" "a" "t" "g" "g" "c" "t" "t" "c" "c" "a" "t" "g" "c" "a" "a" "t" "t"
    [2161] "g" "g" "g" "c" "a" "g" "a" "t" "g" "t" "g" "t" "g" "a" "g" "g" "c" "a"
    [2179] "c" "c" "t" "g" "t" "g" "g" "t" "g" "a" "c" "c" "c" "g" "a" "g" "a" "g"
    [2197] "t" "g" "g" "g" "t" "g" "t" "t" "g" "g" "a" "c" "a" "g" "t" "g" "t" "a"
    [2215] "g" "c" "a" "c" "t" "c" "t" "a" "c" "c" "a" "g" "t" "g" "c" "c" "a" "g"
    [2233] "g" "a" "g" "c" "t" "g" "g" "a" "c" "a" "c" "c" "t" "a" "c" "c" "t" "g"
    [2251] "a" "t" "a" "c" "c" "c" "c" "a" "g" "a" "t" "c" "c" "c" "c" "c" "a" "c"
    [2269] "a" "g" "c" "c" "a" "c" "t" "a" "c" "t" "g" "a"
    attr(,"name")
    [1] "HSU64805"
    attr(,"Annot")
    [1] ">HSU64805"
    attr(,"class")
    [1] "SeqFastadna"
    
    $HSU68041.BRCA1
      [1] "a" "g" "t" "g" "t" "g" "a" "g" "c" "a" "g" "g" "g" "a" "g" "a" "a" "g"
     [19] "c" "c" "a" "g" "a" "a" "t" "t" "g" "a" "c" "a" "g" "c" "t" "t" "c" "a"
     [37] "a" "c" "a" "g" "a" "a" "a" "g" "g" "g" "t" "c" "a" "a" "c" "a" "a" "a"
     [55] "a" "g" "a" "a" "t" "g" "t" "c" "c" "a" "t" "g" "g" "t" "g" "g" "t" "g"
     [73] "t" "c" "t" "g" "g" "c" "c" "t" "g" "a" "c" "c" "c" "c" "a" "g" "a" "a"
     [91] "g" "a" "a" "t" "t" "t" "a" "t" "g" "c" "t" "c" "g" "t" "g" "t" "a" "c"
    [109] "a" "a" "g" "t" "t" "t" "g" "c" "c" "a" "g" "a" "a" "a" "a" "c" "a" "c"
    [127] "c" "a" "c" "a" "t" "c" "a" "c" "t" "t" "t" "a" "a" "c" "t" "a" "a" "t"
    [145] "c" "t" "a" "a" "t" "t" "a" "c" "t" "g" "a" "a" "g" "a" "g" "a" "c" "t"
    [163] "a" "c" "t" "c" "a" "t" "g" "t" "t" "g" "t" "t" "a" "t" "g" "a" "a" "a"
    [181] "a" "c" "a" "g" "a" "t" "g" "c" "t" "g" "a" "g" "t" "t" "t" "g" "t" "g"
    [199] "t" "g" "t" "g" "a" "a" "c" "g" "g" "a" "c" "a" "c" "t" "g" "a" "a" "a"
    [217] "t" "a" "t" "t" "t" "t" "c" "t" "a" "g" "g" "a" "a" "t" "t" "g" "c" "g"
    [235] "g" "g" "a" "g" "g" "a" "a" "a" "a" "t" "g" "g" "g" "t" "a" "g" "t" "t"
    [253] "a" "g" "c" "t" "a" "t" "t" "t" "c" "t" "g" "g" "g" "t" "g" "a" "c" "c"
    [271] "c" "a" "g" "t" "c" "t" "a" "t" "t" "a" "a" "a" "g" "a" "a" "a" "g" "a"
    [289] "a" "a" "a" "a" "t" "g" "c" "t" "g" "a" "a" "t" "g" "a" "g" "c" "a" "t"
    [307] "g" "a" "t" "t" "t" "t" "g" "a" "a" "g" "t" "c" "a" "g" "a" "g" "g" "a"
    [325] "g" "a" "t" "g" "t" "g" "g" "t" "c" "a" "a" "t" "g" "g" "a" "a" "g" "a"
    [343] "a" "a" "c" "c" "a" "c" "c" "a" "a" "g" "g" "t" "c" "c" "a" "a" "a" "g"
    [361] "c" "g" "a" "g" "c" "a" "a" "g" "a" "g" "a" "a" "t" "c" "c" "c" "a" "g"
    [379] "g" "a" "c" "a" "g" "a" "a" "a" "g" "a" "t" "c" "t" "t" "c" "a" "g" "g"
    [397] "g" "g" "g" "c" "t" "a" "g" "a" "a" "a" "t" "c" "t" "g" "t" "t" "g" "c"
    [415] "t" "a" "t" "g" "g" "g" "c" "c" "c" "t" "t" "c" "a" "c" "c" "a" "a" "c"
    [433] "a" "t" "g" "c" "c" "c" "a" "c" "a" "g" "a" "t" "c" "a" "a" "c" "t" "g"
    [451] "g" "a" "a" "t" "g" "g" "a" "t" "g" "g" "t" "a" "c" "a" "g" "c" "t" "g"
    [469] "t" "g" "t" "g" "g" "t" "g" "c" "t" "t" "c" "t" "g" "t" "g" "g" "t" "g"
    [487] "a" "a" "g" "g" "a" "g" "c" "t" "t" "t" "c" "a" "t" "c" "a" "t" "t" "c"
    [505] "a" "c" "c" "c" "t" "t" "g" "g" "c" "a" "c" "a" "g" "g" "t" "g" "t" "c"
    [523] "c" "a" "c" "c" "c" "a" "a" "t" "t" "g" "t" "g" "g" "t" "t" "g" "t" "g"
    [541] "c" "a" "g" "c" "c" "a" "g" "a" "t" "g" "c" "c" "t" "g" "g" "a" "c" "a"
    [559] "g" "a" "g" "g" "a" "c" "a" "a" "t" "g" "g" "c" "t" "t" "c" "c" "a" "t"
    [577] "g" "c" "a" "a" "t" "t" "g" "g" "g" "c" "a" "g" "a" "t" "g" "t" "g" "t"
    [595] "g" "a" "g" "g" "c" "a" "c" "c" "t" "g" "t" "g" "g" "t" "g" "a" "c" "c"
    [613] "c" "g" "a" "g" "a" "g" "t" "g" "g" "g" "t" "g" "t" "t" "g" "g" "a" "c"
    [631] "a" "g" "t" "g" "t" "a" "g" "c" "a" "c" "t" "c" "t" "a" "c" "c" "a" "g"
    [649] "t" "g" "c" "c" "a" "g" "g" "a" "g" "c" "t" "g" "g" "a" "c" "a" "c" "c"
    [667] "t" "a" "c" "c" "t" "c" "a" "t" "a" "c" "c" "c" "c" "a" "g" "a" "t" "c"
    [685] "c" "c" "c" "c" "a" "c" "a" "g" "c" "c" "a" "c" "t" "a" "c" "t" "g" "a"
    attr(,"name")
    [1] "HSU68041.BRCA1"
    attr(,"Annot")
    [1] ">HSU68041.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $JN384124.BRCA1
      [1] "c" "t" "g" "a" "t" "g" "a" "c" "t" "c" "a" "c" "a" "t" "g" "a" "t" "g"
     [19] "c" "g" "c" "c" "c" "c" "c" "c" "c" "c" "c" "a" "t" "c" "a" "a" "a" "t"
     [37] "g" "c" "c" "a" "a" "a" "g" "t" "a" "g" "c" "t" "g" "a" "t" "g" "t" "a"
     [55] "t" "t" "g" "g" "a" "c" "g" "t" "t" "c" "t" "a" "a" "a" "t" "g" "a" "g"
     [73] "g" "t" "a" "g" "a" "t" "g" "a" "a" "t" "a" "t" "t" "c" "t" "g" "g" "t"
     [91] "t" "c" "t" "t" "c" "a" "g" "a" "g" "a" "a" "a" "a" "t" "a" "g" "a" "c"
    [109] "t" "t" "a" "c" "t" "g" "g" "c" "c" "a" "g" "t" "g" "a" "t" "c" "c" "t"
    [127] "c" "a" "t" "g" "a" "g" "g" "c" "t" "t" "t" "a" "a" "t" "a" "t" "g" "t"
    [145] "a" "a" "a" "a" "g" "t" "g" "a" "a" "a" "g" "a" "g" "t" "t" "c" "a" "c"
    [163] "t" "c" "c" "a" "a" "a" "t" "c" "a" "g" "t" "a" "g" "a" "g" "a" "g" "t"
    [181] "a" "a" "t" "a" "t" "t" "g" "a" "a" "g" "a" "c" "a" "a" "a" "a" "t" "g"
    [199] "a" "t" "t" "t" "g" "g" "g" "a" "a" "a" "a" "c" "c" "t" "a" "t" "c" "g"
    [217] "g" "a" "a" "g" "a" "a" "g" "g" "c" "a" "a" "g" "c" "c" "t" "c" "c" "c"
    [235] "c" "a" "a" "c" "t" "a" "a" "g" "c" "c" "a" "t" "t" "t" "a" "t" "c" "t"
    [253] "g" "a" "a" "a" "a" "t" "t"
    attr(,"name")
    [1] "JN384124.BRCA1"
    attr(,"Annot")
    [1] ">JN384124.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $JN686490
       [1] "a" "t" "g" "g" "a" "t" "t" "t" "a" "t" "c" "t" "g" "c" "t" "c" "t" "t"
      [19] "c" "g" "c" "g" "t" "t" "g" "a" "a" "g" "a" "a" "g" "t" "a" "c" "a" "a"
      [37] "a" "a" "t" "g" "t" "c" "a" "t" "t" "a" "a" "t" "g" "c" "t" "a" "t" "g"
      [55] "c" "a" "g" "a" "a" "a" "a" "t" "c" "t" "t" "a" "g" "a" "g" "t" "g" "t"
      [73] "c" "c" "c" "a" "t" "c" "t" "g" "t" "c" "t" "g" "g" "a" "g" "t" "t" "g"
      [91] "a" "t" "c" "a" "a" "g" "g" "a" "a" "c" "c" "t" "g" "t" "c" "t" "c" "c"
     [109] "a" "c" "a" "a" "a" "g" "t" "g" "t" "g" "a" "c" "c" "a" "c" "a" "t" "a"
     [127] "t" "t" "t" "t" "g" "c" "a" "a" "a" "t" "t" "t" "t" "g" "c" "a" "t" "g"
     [145] "c" "t" "g" "a" "a" "a" "c" "t" "t" "c" "t" "c" "a" "a" "c" "c" "a" "g"
     [163] "a" "a" "g" "a" "a" "a" "g" "g" "g" "c" "c" "t" "t" "c" "a" "c" "a" "g"
     [181] "t" "g" "t" "c" "c" "t" "t" "t" "a" "t" "g" "t" "a" "a" "g" "a" "a" "t"
     [199] "g" "a" "t" "a" "t" "a" "a" "c" "c" "a" "a" "a" "a" "g" "g" "a" "g" "c"
     [217] "c" "t" "a" "c" "a" "a" "g" "a" "a" "a" "g" "t" "a" "c" "g" "a" "g" "a"
     [235] "t" "t" "t" "a" "g" "t" "c" "a" "a" "c" "t" "t" "g" "t" "t" "g" "a" "a"
     [253] "g" "a" "g" "c" "t" "a" "t" "t" "g" "a" "a" "a" "a" "t" "c" "a" "t" "t"
     [271] "t" "g" "t" "g" "c" "t" "t" "t" "t" "c" "a" "g" "c" "t" "t" "g" "a" "c"
     [289] "a" "c" "a" "g" "g" "t" "t" "t" "g" "g" "a" "g" "t" "a" "t" "g" "c" "a"
     [307] "a" "a" "c" "a" "g" "c" "t" "a" "t" "a" "a" "t" "t" "t" "t" "g" "c" "a"
     [325] "a" "a" "a" "a" "a" "g" "g" "a" "a" "a" "a" "t" "a" "a" "c" "t" "c" "t"
     [343] "c" "c" "t" "g" "a" "a" "c" "a" "t" "c" "t" "a" "a" "a" "a" "g" "a" "t"
     [361] "g" "a" "a" "g" "t" "t" "t" "c" "t" "a" "t" "c" "a" "t" "c" "c" "a" "a"
     [379] "a" "g" "t" "a" "t" "g" "g" "g" "c" "t" "a" "c" "a" "g" "a" "a" "a" "c"
     [397] "c" "g" "t" "g" "c" "c" "a" "a" "a" "a" "g" "a" "c" "t" "t" "c" "t" "a"
     [415] "c" "a" "g" "a" "g" "t" "g" "a" "a" "c" "c" "c" "g" "a" "a" "a" "a" "t"
     [433] "c" "c" "t" "t" "c" "c" "t" "t" "g" "c" "a" "g" "g" "a" "a" "a" "c" "c"
     [451] "a" "g" "t" "c" "t" "c" "a" "g" "t" "g" "t" "c" "c" "a" "a" "c" "t" "c"
     [469] "t" "c" "t" "a" "a" "c" "c" "t" "t" "g" "g" "a" "a" "c" "t" "g" "t" "g"
     [487] "a" "g" "a" "a" "c" "t" "c" "t" "g" "a" "g" "g" "a" "c" "a" "a" "a" "g"
     [505] "c" "a" "g" "c" "g" "g" "a" "t" "a" "c" "a" "a" "c" "c" "t" "c" "a" "a"
     [523] "a" "a" "g" "a" "c" "g" "t" "c" "t" "g" "t" "c" "t" "a" "c" "a" "t" "t"
     [541] "g" "a" "a" "t" "t" "g" "g" "g" "a" "t" "c" "t" "g" "a" "t" "t" "c" "t"
     [559] "t" "c" "t" "g" "a" "a" "g" "a" "t" "a" "c" "c" "g" "t" "t" "a" "a" "t"
     [577] "a" "a" "g" "g" "c" "a" "a" "c" "t" "t" "a" "t" "t" "g" "c" "a" "g" "t"
     [595] "g" "t" "g" "g" "g" "a" "g" "a" "t" "c" "a" "a" "g" "a" "a" "t" "t" "g"
     [613] "t" "t" "a" "c" "a" "a" "a" "t" "c" "a" "c" "c" "c" "c" "t" "c" "a" "a"
     [631] "g" "g" "a" "a" "c" "c" "a" "g" "g" "g" "a" "t" "g" "a" "a" "a" "t" "c"
     [649] "a" "g" "t" "t" "t" "g" "g" "a" "t" "t" "c" "t" "g" "c" "a" "a" "a" "a"
     [667] "a" "a" "g" "g" "c" "t" "g" "c" "t" "t" "g" "t" "g" "a" "a" "t" "t" "t"
     [685] "t" "c" "t" "g" "a" "g" "a" "c" "g" "g" "a" "t" "g" "t" "a" "a" "c" "a"
     [703] "a" "a" "t" "a" "c" "t" "g" "a" "a" "c" "a" "t" "c" "a" "t" "c" "a" "a"
     [721] "c" "c" "c" "a" "g" "t" "a" "a" "t" "a" "a" "t" "g" "a" "t" "t" "t" "g"
     [739] "a" "a" "c" "a" "c" "c" "a" "c" "t" "g" "a" "g" "a" "a" "g" "c" "g" "t"
     [757] "g" "c" "a" "g" "c" "t" "g" "a" "g" "a" "g" "g" "c" "a" "t" "c" "c" "a"
     [775] "g" "a" "a" "a" "a" "g" "t" "a" "t" "c" "a" "g" "g" "g" "t" "a" "g" "t"
     [793] "t" "c" "t" "g" "t" "t" "t" "c" "a" "a" "a" "c" "t" "t" "g" "c" "a" "t"
     [811] "g" "t" "g" "g" "a" "g" "c" "c" "a" "t" "g" "t" "g" "g" "c" "a" "c" "a"
     [829] "a" "a" "t" "a" "c" "t" "c" "a" "t" "g" "c" "c" "a" "g" "c" "t" "c" "a"
     [847] "t" "t" "a" "c" "a" "g" "c" "a" "t" "g" "a" "g" "a" "a" "c" "a" "g" "c"
     [865] "a" "g" "t" "t" "t" "a" "t" "t" "a" "c" "t" "c" "a" "c" "t" "a" "a" "a"
     [883] "g" "a" "c" "a" "g" "a" "a" "t" "g" "a" "a" "t" "g" "t" "a" "g" "a" "a"
     [901] "a" "a" "g" "g" "c" "t" "g" "a" "a" "t" "t" "c" "t" "g" "t" "a" "a" "t"
     [919] "a" "a" "a" "a" "g" "c" "a" "a" "a" "c" "a" "g" "c" "c" "t" "g" "g" "c"
     [937] "t" "t" "a" "g" "c" "a" "a" "g" "g" "a" "g" "c" "c" "a" "a" "c" "a" "t"
     [955] "a" "a" "c" "a" "g" "a" "t" "g" "g" "g" "c" "t" "g" "g" "a" "a" "g" "t"
     [973] "a" "a" "g" "g" "a" "a" "a" "c" "a" "t" "g" "t" "a" "a" "t" "g" "a" "t"
     [991] "a" "g" "g" "c" "g" "g" "a" "c" "t" "c" "c" "c" "a" "g" "c" "a" "c" "a"
    [1009] "g" "a" "a" "a" "a" "a" "a" "a" "g" "g" "t" "a" "g" "a" "t" "c" "t" "g"
    [1027] "a" "a" "t" "g" "c" "t" "g" "a" "t" "c" "c" "c" "c" "t" "g" "t" "g" "t"
    [1045] "g" "a" "g" "a" "g" "a" "a" "a" "a" "g" "a" "a" "t" "g" "g" "a" "a" "t"
    [1063] "a" "a" "g" "c" "a" "g" "a" "a" "a" "c" "t" "g" "c" "c" "a" "t" "g" "c"
    [1081] "t" "c" "a" "g" "a" "g" "a" "a" "t" "c" "c" "t" "a" "g" "a" "g" "a" "t"
    [1099] "a" "c" "t" "g" "a" "a" "g" "a" "t" "g" "t" "t" "c" "c" "t" "t" "g" "g"
    [1117] "a" "t" "a" "a" "c" "a" "c" "t" "a" "a" "a" "t" "a" "g" "c" "a" "g" "c"
    [1135] "a" "t" "t" "c" "a" "g" "a" "a" "a" "g" "t" "t" "a" "a" "t" "g" "a" "g"
    [1153] "t" "g" "g" "t" "t" "t" "t" "c" "c" "a" "g" "a" "a" "g" "t" "g" "a" "t"
    [1171] "g" "a" "a" "c" "t" "g" "t" "t" "a" "g" "g" "t" "t" "c" "t" "g" "a" "t"
    [1189] "g" "a" "c" "t" "c" "a" "c" "a" "t" "g" "a" "t" "g" "g" "g" "g" "a" "g"
    [1207] "t" "c" "t" "g" "a" "a" "t" "c" "a" "a" "a" "t" "g" "c" "c" "a" "a" "a"
    [1225] "g" "t" "a" "g" "c" "t" "g" "a" "t" "g" "t" "a" "t" "t" "g" "g" "a" "c"
    [1243] "g" "t" "t" "c" "t" "a" "a" "a" "t" "g" "a" "g" "g" "t" "a" "g" "a" "t"
    [1261] "g" "a" "a" "t" "a" "t" "t" "c" "t" "g" "g" "t" "t" "c" "t" "t" "c" "a"
    [1279] "g" "a" "g" "a" "a" "a" "a" "t" "a" "g" "a" "c" "t" "t" "a" "c" "t" "g"
    [1297] "g" "c" "c" "a" "g" "t" "g" "a" "t" "c" "c" "t" "c" "a" "t" "g" "a" "g"
    [1315] "g" "c" "t" "t" "t" "a" "a" "t" "a" "t" "g" "t" "a" "a" "a" "a" "g" "t"
    [1333] "g" "a" "a" "a" "g" "a" "g" "t" "t" "c" "a" "c" "t" "c" "c" "a" "a" "a"
    [1351] "t" "c" "a" "g" "t" "a" "g" "a" "g" "a" "g" "t" "a" "a" "t" "a" "t" "t"
    [1369] "g" "a" "a" "g" "a" "c" "a" "a" "a" "a" "t" "a" "t" "t" "t" "g" "g" "g"
    [1387] "a" "a" "a" "a" "c" "c" "t" "a" "t" "c" "g" "g" "a" "a" "g" "a" "a" "g"
    [1405] "g" "c" "a" "a" "g" "c" "c" "t" "c" "c" "c" "c" "a" "a" "c" "t" "t" "a"
    [1423] "a" "g" "c" "c" "a" "t" "g" "t" "a" "a" "c" "t" "g" "a" "a" "a" "a" "t"
    [1441] "c" "t" "a" "a" "t" "t" "a" "t" "a" "g" "g" "a" "g" "c" "a" "t" "t" "t"
    [1459] "g" "t" "t" "a" "c" "t" "g" "a" "g" "c" "c" "a" "c" "a" "g" "a" "t" "a"
    [1477] "a" "t" "a" "c" "a" "a" "g" "a" "g" "c" "g" "t" "c" "c" "c" "c" "t" "c"
    [1495] "a" "c" "a" "a" "a" "t" "a" "a" "a" "t" "t" "a" "a" "a" "g" "c" "g" "t"
    [1513] "a" "a" "a" "a" "g" "g" "a" "g" "a" "c" "c" "t" "a" "c" "a" "t" "c" "a"
    [1531] "g" "g" "c" "c" "t" "t" "c" "a" "t" "c" "c" "t" "g" "a" "g" "g" "a" "t"
    [1549] "t" "t" "t" "a" "t" "c" "a" "a" "g" "a" "a" "a" "g" "c" "a" "g" "a" "t"
    [1567] "t" "t" "g" "g" "c" "a" "g" "t" "t" "c" "a" "a" "a" "a" "g" "a" "c" "t"
    [1585] "c" "c" "t" "g" "a" "a" "a" "t" "g" "a" "t" "a" "a" "a" "t" "c" "a" "g"
    [1603] "g" "g" "a" "a" "c" "t" "a" "a" "c" "c" "a" "a" "a" "c" "g" "g" "a" "g"
    [1621] "c" "a" "g" "a" "a" "t" "g" "g" "t" "c" "a" "a" "g" "t" "g" "a" "t" "g"
    [1639] "a" "a" "t" "a" "t" "t" "a" "c" "t" "a" "a" "t" "a" "g" "t" "g" "g" "t"
    [1657] "c" "a" "t" "g" "a" "g" "a" "a" "t" "a" "a" "a" "a" "c" "a" "a" "a" "a"
    [1675] "g" "g" "t" "g" "a" "t" "t" "c" "t" "a" "t" "t" "c" "a" "g" "a" "a" "t"
    [1693] "g" "a" "g" "a" "a" "a" "a" "a" "t" "c" "c" "t" "a" "a" "c" "c" "c" "a"
    [1711] "a" "t" "a" "g" "a" "a" "t" "c" "a" "c" "t" "c" "g" "a" "a" "a" "a" "a"
    [1729] "g" "a" "a" "t" "c" "t" "g" "c" "t" "t" "t" "c" "a" "a" "a" "a" "c" "g"
    [1747] "a" "a" "a" "g" "c" "t" "g" "a" "a" "c" "c" "t" "a" "t" "a" "a" "g" "c"
    [1765] "a" "g" "c" "a" "g" "t" "a" "t" "a" "a" "g" "c" "a" "a" "t" "a" "t" "g"
    [1783] "g" "a" "a" "c" "t" "c" "g" "a" "a" "t" "t" "a" "a" "a" "t" "a" "t" "c"
    [1801] "c" "a" "c" "a" "a" "t" "t" "c" "a" "a" "a" "a" "g" "c" "a" "c" "c" "t"
    [1819] "a" "a" "a" "a" "a" "g" "a" "a" "t" "a" "g" "g" "c" "t" "g" "a" "g" "g"
    [1837] "a" "g" "g" "a" "a" "g" "t" "c" "t" "t" "c" "t" "a" "c" "c" "a" "g" "g"
    [1855] "c" "a" "t" "a" "t" "t" "c" "a" "t" "g" "c" "g" "c" "t" "t" "g" "a" "a"
    [1873] "c" "t" "a" "g" "t" "a" "g" "t" "c" "a" "g" "t" "a" "g" "a" "a" "a" "t"
    [1891] "c" "t" "a" "a" "g" "c" "c" "c" "a" "c" "c" "t" "a" "a" "t" "t" "g" "t"
    [1909] "a" "c" "t" "g" "a" "a" "t" "t" "g" "c" "a" "a" "a" "t" "t" "g" "a" "t"
    [1927] "a" "g" "t" "t" "g" "t" "t" "c" "t" "a" "g" "c" "a" "g" "t" "g" "a" "a"
    [1945] "g" "a" "g" "a" "t" "a" "a" "a" "g" "a" "a" "a" "a" "a" "a" "a" "a" "g"
    [1963] "t" "a" "c" "a" "a" "c" "c" "a" "a" "a" "t" "g" "c" "c" "a" "g" "t" "c"
    [1981] "a" "g" "g" "c" "a" "c" "a" "g" "c" "a" "g" "a" "a" "a" "c" "c" "t" "a"
    [1999] "c" "a" "a" "c" "t" "c" "a" "t" "g" "g" "a" "a" "g" "g" "t" "a" "a" "a"
    [2017] "g" "a" "a" "c" "c" "t" "g" "c" "a" "a" "c" "t" "g" "g" "a" "g" "c" "c"
    [2035] "a" "a" "g" "a" "a" "g" "a" "g" "t" "a" "a" "c" "a" "a" "g" "c" "c" "a"
    [2053] "a" "a" "t" "g" "a" "a" "c" "a" "g" "a" "c" "a" "a" "g" "t" "a" "a" "a"
    [2071] "a" "g" "a" "c" "a" "t" "g" "a" "c" "a" "g" "c" "g" "a" "t" "a" "c" "t"
    [2089] "t" "t" "c" "c" "c" "a" "g" "a" "g" "c" "t" "g" "a" "a" "g" "t" "t" "a"
    [2107] "a" "c" "a" "a" "a" "t" "g" "c" "a" "c" "c" "t" "g" "g" "t" "t" "c" "t"
    [2125] "t" "t" "t" "a" "c" "t" "a" "a" "g" "t" "g" "t" "t" "c" "a" "a" "a" "t"
    [2143] "a" "c" "c" "a" "g" "t" "g" "a" "a" "c" "t" "t" "a" "a" "a" "g" "a" "a"
    [2161] "t" "t" "t" "g" "t" "c" "a" "a" "t" "c" "c" "t" "a" "g" "c" "c" "t" "t"
    [2179] "c" "c" "a" "a" "g" "a" "g" "a" "a" "g" "a" "a" "a" "a" "a" "g" "a" "a"
    [2197] "g" "a" "g" "a" "a" "a" "c" "t" "a" "g" "a" "a" "a" "c" "a" "g" "t" "t"
    [2215] "a" "a" "a" "g" "t" "g" "t" "c" "t" "a" "a" "t" "a" "a" "t" "g" "c" "t"
    [2233] "g" "a" "a" "g" "a" "c" "c" "c" "c" "a" "a" "a" "g" "a" "t" "c" "t" "c"
    [2251] "a" "t" "g" "t" "t" "a" "a" "g" "t" "g" "g" "a" "g" "a" "a" "a" "g" "g"
    [2269] "g" "t" "t" "t" "t" "g" "c" "a" "a" "a" "c" "t" "g" "a" "a" "a" "g" "t"
    [2287] "t" "c" "t" "g" "t" "a" "g" "a" "g" "a" "g" "t" "a" "g" "c" "a" "g" "t"
    [2305] "a" "t" "t" "t" "c" "a" "t" "t" "g" "g" "t" "a" "c" "c" "t" "g" "g" "t"
    [2323] "a" "c" "t" "g" "a" "t" "t" "a" "t" "g" "g" "c" "a" "c" "t" "c" "a" "g"
    [2341] "g" "a" "a" "a" "g" "t" "a" "t" "c" "t" "c" "g" "t" "t" "a" "c" "t" "g"
    [2359] "g" "a" "a" "g" "t" "t" "a" "g" "c" "a" "c" "t" "c" "t" "a" "g" "g" "g"
    [2377] "a" "a" "g" "g" "c" "a" "a" "a" "a" "a" "c" "a" "g" "a" "a" "c" "c" "a"
    [2395] "a" "a" "t" "a" "a" "a" "t" "g" "t" "g" "t" "g" "a" "g" "t" "c" "a" "g"
    [2413] "t" "g" "t" "g" "c" "a" "g" "c" "a" "t" "t" "t" "g" "a" "a" "a" "a" "c"
    [2431] "c" "c" "c" "a" "a" "g" "g" "g" "a" "c" "t" "a" "a" "t" "t" "c" "a" "t"
    [2449] "g" "g" "t" "t" "g" "t" "t" "c" "c" "a" "a" "a" "g" "a" "t" "a" "a" "t"
    [2467] "a" "g" "a" "a" "a" "t" "g" "a" "c" "a" "c" "a" "g" "a" "a" "g" "g" "c"
    [2485] "t" "t" "t" "a" "a" "g" "t" "a" "t" "c" "c" "a" "t" "t" "g" "g" "g" "a"
    [2503] "c" "a" "t" "g" "a" "a" "g" "t" "t" "a" "a" "c" "c" "a" "c" "a" "g" "t"
    [2521] "c" "g" "g" "g" "a" "a" "a" "c" "a" "a" "g" "c" "a" "t" "a" "g" "a" "a"
    [2539] "a" "t" "g" "g" "a" "a" "g" "a" "a" "a" "g" "t" "g" "a" "a" "c" "t" "t"
    [2557] "g" "a" "t" "g" "c" "t" "c" "a" "g" "t" "a" "t" "t" "t" "g" "c" "a" "g"
    [2575] "a" "a" "t" "a" "c" "a" "t" "t" "c" "a" "a" "g" "g" "t" "t" "t" "c" "a"
    [2593] "a" "a" "g" "c" "g" "c" "c" "a" "g" "t" "c" "a" "t" "t" "t" "g" "c" "t"
    [2611] "c" "c" "g" "t" "t" "t" "t" "c" "a" "a" "a" "t" "c" "c" "a" "g" "g" "a"
    [2629] "a" "a" "t" "g" "c" "a" "g" "a" "a" "g" "a" "g" "g" "a" "a" "t" "g" "t"
    [2647] "g" "c" "a" "a" "c" "a" "t" "t" "c" "t" "c" "t" "g" "c" "c" "c" "a" "c"
    [2665] "t" "c" "t" "g" "g" "g" "t" "c" "c" "t" "t" "a" "a" "a" "g" "a" "a" "a"
    [2683] "c" "a" "a" "a" "g" "t" "c" "c" "a" "a" "a" "a" "g" "t" "c" "a" "c" "t"
    [2701] "t" "t" "t" "g" "a" "a" "t" "g" "t" "g" "a" "a" "c" "a" "a" "a" "a" "g"
    [2719] "g" "a" "a" "g" "a" "a" "a" "a" "t" "c" "a" "a" "g" "g" "a" "a" "a" "g"
    [2737] "a" "a" "t" "g" "a" "g" "t" "c" "t" "a" "a" "t" "a" "t" "c" "a" "a" "g"
    [2755] "c" "c" "t" "g" "t" "a" "c" "a" "g" "a" "c" "a" "g" "t" "t" "a" "a" "t"
    [2773] "a" "t" "c" "a" "c" "t" "g" "c" "a" "g" "g" "c" "t" "t" "t" "c" "c" "t"
    [2791] "g" "t" "g" "g" "t" "t" "g" "g" "t" "c" "a" "g" "a" "a" "a" "g" "a" "t"
    [2809] "a" "a" "g" "c" "c" "a" "g" "t" "t" "g" "a" "t" "a" "a" "t" "g" "c" "c"
    [2827] "a" "a" "a" "t" "g" "t" "a" "g" "t" "a" "t" "c" "a" "a" "a" "g" "g" "a"
    [2845] "g" "g" "c" "t" "c" "t" "a" "g" "g" "t" "t" "t" "t" "g" "t" "c" "t" "a"
    [2863] "t" "c" "a" "t" "c" "t" "c" "a" "g" "t" "t" "c" "a" "g" "a" "g" "g" "c"
    [2881] "a" "a" "c" "g" "a" "a" "a" "c" "t" "g" "g" "a" "c" "t" "c" "a" "t" "t"
    [2899] "a" "c" "t" "c" "c" "a" "a" "a" "t" "a" "a" "a" "c" "a" "t" "g" "g" "a"
    [2917] "c" "t" "t" "t" "t" "a" "c" "a" "a" "a" "a" "c" "c" "c" "a" "t" "a" "t"
    [2935] "c" "g" "t" "a" "t" "a" "c" "c" "a" "c" "c" "a" "c" "t" "t" "t" "t" "t"
    [2953] "c" "c" "c" "a" "t" "c" "a" "a" "g" "t" "c" "a" "t" "t" "t" "g" "t" "t"
    [2971] "a" "a" "a" "a" "c" "t" "a" "a" "a" "t" "g" "t" "a" "a" "g" "a" "a" "a"
    [2989] "a" "a" "t" "c" "t" "g" "c" "t" "a" "g" "a" "g" "g" "a" "a" "a" "a" "c"
    [3007] "t" "t" "t" "g" "a" "g" "g" "a" "a" "c" "a" "t" "t" "c" "a" "a" "t" "g"
    [3025] "t" "c" "a" "c" "c" "t" "g" "a" "a" "a" "g" "a" "g" "a" "a" "a" "t" "g"
    [3043] "g" "g" "a" "a" "a" "t" "g" "a" "g" "a" "a" "c" "a" "t" "t" "c" "c" "a"
    [3061] "a" "g" "t" "a" "c" "a" "g" "t" "g" "a" "g" "c" "a" "c" "a" "a" "t" "t"
    [3079] "a" "g" "c" "c" "g" "t" "a" "a" "t" "a" "a" "c" "a" "t" "t" "a" "g" "a"
    [3097] "g" "a" "a" "a" "a" "t" "g" "t" "t" "t" "t" "t" "a" "a" "a" "g" "a" "a"
    [3115] "g" "c" "c" "a" "g" "c" "t" "c" "a" "a" "g" "c" "a" "a" "t" "a" "t" "t"
    [3133] "a" "a" "t" "g" "a" "a" "g" "t" "a" "g" "g" "t" "t" "c" "c" "a" "g" "t"
    [3151] "a" "c" "t" "a" "a" "t" "g" "a" "a" "g" "t" "g" "g" "g" "c" "t" "c" "c"
    [3169] "a" "g" "t" "a" "t" "t" "a" "a" "t" "g" "a" "a" "a" "t" "a" "g" "g" "t"
    [3187] "t" "c" "c" "a" "g" "t" "g" "a" "t" "g" "a" "a" "a" "a" "c" "a" "t" "t"
    [3205] "c" "a" "a" "g" "c" "a" "g" "a" "a" "c" "t" "a" "g" "g" "t" "a" "g" "a"
    [3223] "a" "a" "c" "a" "g" "a" "g" "g" "g" "c" "c" "a" "a" "a" "a" "t" "t" "g"
    [3241] "a" "a" "t" "g" "c" "t" "a" "t" "g" "c" "t" "t" "a" "g" "a" "t" "t" "a"
    [3259] "g" "g" "g" "g" "t" "t" "t" "t" "g" "c" "a" "a" "c" "c" "t" "g" "a" "g"
    [3277] "g" "t" "c" "t" "a" "t" "a" "a" "a" "c" "a" "a" "a" "g" "t" "c" "t" "t"
    [3295] "c" "c" "t" "g" "g" "a" "a" "g" "t" "a" "a" "t" "t" "g" "t" "a" "a" "g"
    [3313] "c" "a" "t" "c" "c" "t" "g" "a" "a" "a" "t" "a" "a" "a" "a" "a" "a" "g"
    [3331] "c" "a" "a" "g" "a" "a" "t" "a" "t" "g" "a" "a" "g" "a" "a" "g" "t" "a"
    [3349] "g" "t" "t" "c" "a" "g" "a" "c" "t" "g" "t" "t" "a" "a" "t" "a" "c" "a"
    [3367] "g" "a" "t" "t" "t" "c" "t" "c" "t" "c" "c" "a" "t" "a" "t" "c" "t" "g"
    [3385] "a" "t" "t" "t" "c" "a" "g" "a" "t" "a" "a" "c" "t" "t" "a" "g" "a" "a"
    [3403] "c" "a" "g" "c" "c" "t" "a" "t" "g" "g" "g" "a" "a" "g" "t" "a" "g" "t"
    [3421] "c" "a" "t" "g" "c" "a" "t" "c" "t" "c" "a" "g" "g" "t" "t" "t" "g" "t"
    [3439] "t" "c" "t" "g" "a" "g" "a" "c" "a" "c" "c" "t" "g" "a" "t" "g" "a" "c"
    [3457] "c" "t" "g" "t" "t" "a" "g" "a" "t" "g" "a" "t" "g" "g" "t" "g" "a" "a"
    [3475] "a" "t" "a" "a" "a" "g" "g" "a" "a" "g" "a" "t" "a" "c" "t" "a" "g" "t"
    [3493] "t" "t" "t" "g" "c" "t" "g" "a" "a" "a" "a" "t" "g" "a" "c" "a" "t" "t"
    [3511] "a" "a" "g" "g" "a" "a" "a" "g" "t" "t" "c" "t" "g" "c" "t" "g" "t" "t"
    [3529] "t" "t" "t" "a" "g" "c" "a" "a" "a" "a" "g" "c" "g" "t" "c" "c" "a" "g"
    [3547] "a" "a" "a" "g" "g" "a" "g" "a" "g" "c" "t" "t" "a" "g" "c" "a" "g" "g"
    [3565] "a" "g" "t" "c" "c" "t" "a" "g" "c" "c" "c" "t" "t" "t" "c" "a" "c" "c"
    [3583] "c" "a" "t" "a" "c" "a" "c" "a" "t" "t" "t" "g" "g" "c" "t" "c" "a" "g"
    [3601] "g" "g" "t" "t" "a" "c" "c" "g" "a" "a" "g" "a" "g" "g" "g" "g" "c" "c"
    [3619] "a" "a" "g" "a" "a" "a" "t" "t" "a" "g" "a" "g" "t" "c" "c" "t" "c" "a"
    [3637] "g" "a" "a" "g" "a" "g" "a" "a" "c" "t" "t" "a" "t" "c" "t" "a" "g" "t"
    [3655] "g" "a" "g" "g" "a" "t" "g" "a" "a" "g" "a" "g" "c" "t" "t" "c" "c" "c"
    [3673] "t" "g" "c" "t" "t" "c" "c" "a" "a" "c" "a" "c" "t" "t" "g" "t" "t" "a"
    [3691] "t" "t" "t" "g" "g" "t" "a" "a" "a" "g" "t" "a" "a" "a" "c" "a" "a" "t"
    [3709] "a" "t" "a" "c" "c" "t" "t" "c" "t" "c" "a" "g" "t" "c" "t" "a" "c" "t"
    [3727] "a" "g" "g" "c" "a" "t" "a" "g" "c" "a" "c" "c" "g" "t" "t" "g" "c" "t"
    [3745] "a" "c" "c" "g" "a" "g" "t" "g" "t" "c" "t" "g" "t" "c" "t" "a" "a" "g"
    [3763] "a" "a" "c" "a" "c" "a" "g" "a" "g" "g" "a" "g" "a" "a" "t" "t" "t" "a"
    [3781] "t" "t" "a" "t" "c" "a" "t" "t" "g" "a" "a" "g" "a" "a" "t" "a" "g" "c"
    [3799] "t" "t" "a" "a" "a" "t" "g" "a" "c" "t" "g" "c" "a" "g" "t" "a" "a" "c"
    [3817] "c" "a" "g" "g" "t" "a" "a" "t" "a" "t" "t" "g" "g" "c" "a" "a" "a" "g"
    [3835] "g" "c" "a" "t" "c" "t" "c" "a" "g" "g" "a" "a" "c" "a" "t" "c" "a" "c"
    [3853] "c" "t" "t" "a" "g" "t" "g" "a" "g" "g" "a" "a" "a" "c" "a" "a" "a" "a"
    [3871] "t" "g" "t" "t" "c" "t" "g" "c" "t" "a" "g" "c" "t" "t" "g" "t" "t" "t"
    [3889] "t" "c" "t" "t" "c" "a" "c" "a" "g" "t" "g" "c" "a" "g" "t" "g" "a" "a"
    [3907] "t" "t" "g" "g" "a" "a" "g" "a" "c" "t" "t" "g" "a" "c" "t" "g" "c" "a"
    [3925] "a" "a" "t" "a" "c" "a" "a" "a" "c" "a" "c" "c" "c" "a" "g" "g" "a" "t"
    [3943] "c" "c" "t" "t" "t" "c" "t" "t" "g" "a" "t" "t" "g" "g" "t" "t" "c" "t"
    [3961] "t" "c" "c" "a" "a" "a" "c" "a" "a" "a" "t" "g" "a" "g" "g" "c" "a" "t"
    [3979] "c" "a" "g" "t" "c" "t" "g" "a" "a" "a" "g" "c" "c" "a" "g" "g" "g" "a"
    [3997] "g" "t" "t" "g" "g" "t" "c" "t" "g" "a" "g" "t" "g" "a" "c" "a" "a" "g"
    [4015] "g" "a" "a" "t" "t" "g" "g" "t" "t" "t" "c" "a" "g" "a" "t" "g" "a" "t"
    [4033] "g" "a" "a" "g" "a" "a" "a" "g" "a" "g" "g" "a" "a" "c" "g" "g" "g" "c"
    [4051] "t" "t" "g" "g" "a" "a" "g" "a" "a" "a" "a" "t" "a" "a" "t" "c" "a" "a"
    [4069] "g" "a" "a" "g" "a" "g" "c" "a" "a" "a" "g" "c" "a" "t" "g" "g" "a" "t"
    [4087] "t" "c" "a" "a" "a" "c" "t" "t" "a" "g" "g" "t" "g" "a" "a" "g" "c" "a"
    [4105] "g" "c" "a" "t" "c" "t" "g" "g" "g" "t" "g" "t" "g" "a" "g" "a" "g" "t"
    [4123] "g" "a" "a" "a" "c" "a" "a" "g" "c" "g" "t" "c" "t" "c" "t" "g" "a" "a"
    [4141] "g" "a" "c" "t" "g" "c" "t" "c" "a" "g" "g" "g" "c" "t" "a" "t" "c" "c"
    [4159] "t" "c" "t" "c" "a" "g" "a" "g" "t" "g" "a" "c" "a" "t" "t" "t" "t" "a"
    [4177] "a" "c" "c" "a" "c" "t" "c" "a" "g" "c" "a" "g" "a" "g" "g" "g" "a" "t"
    [4195] "a" "c" "c" "a" "t" "g" "c" "a" "a" "c" "a" "t" "a" "a" "c" "c" "t" "g"
    [4213] "a" "t" "a" "a" "a" "g" "c" "t" "c" "c" "a" "g" "c" "a" "g" "g" "a" "a"
    [4231] "a" "t" "g" "g" "c" "t" "g" "a" "a" "c" "t" "a" "g" "a" "a" "g" "c" "t"
    [4249] "g" "t" "g" "t" "t" "a" "g" "a" "a" "c" "a" "g" "c" "a" "t" "g" "g" "g"
    [4267] "a" "g" "c" "c" "a" "g" "c" "c" "t" "t" "c" "t" "a" "a" "c" "a" "g" "c"
    [4285] "t" "a" "c" "c" "c" "t" "t" "c" "c" "a" "t" "c" "a" "t" "a" "a" "g" "t"
    [4303] "g" "a" "c" "t" "c" "t" "t" "c" "t" "g" "c" "c" "c" "t" "t" "g" "a" "g"
    [4321] "g" "a" "c" "c" "t" "g" "c" "g" "a" "a" "a" "t" "c" "c" "a" "g" "a" "a"
    [4339] "c" "a" "a" "a" "g" "c" "a" "c" "a" "t" "c" "a" "g" "a" "a" "a" "a" "a"
    [4357] "g" "c" "a" "g" "t" "a" "t" "t" "a" "a" "c" "t" "t" "c" "a" "c" "a" "g"
    [4375] "a" "a" "a" "a" "g" "t" "a" "g" "t" "g" "a" "a" "t" "a" "c" "c" "c" "t"
    [4393] "a" "t" "a" "a" "g" "c" "c" "a" "g" "a" "a" "t" "c" "c" "a" "g" "a" "a"
    [4411] "g" "g" "c" "c" "t" "t" "t" "c" "t" "g" "c" "t" "g" "a" "c" "a" "a" "g"
    [4429] "t" "t" "t" "g" "a" "g" "g" "t" "g" "t" "c" "t" "g" "c" "a" "g" "a" "t"
    [4447] "a" "g" "t" "t" "c" "t" "a" "c" "c" "a" "g" "t" "a" "a" "a" "a" "a" "t"
    [4465] "a" "a" "a" "g" "a" "a" "c" "c" "a" "g" "g" "a" "g" "t" "g" "g" "a" "a"
    [4483] "a" "g" "g" "t" "c" "a" "t" "c" "c" "c" "c" "t" "t" "c" "t" "a" "a" "a"
    [4501] "t" "g" "c" "c" "c" "a" "t" "c" "a" "t" "t" "a" "g" "a" "t" "g" "a" "t"
    [4519] "a" "g" "g" "t" "g" "g" "t" "a" "c" "a" "t" "g" "c" "a" "c" "a" "g" "t"
    [4537] "t" "g" "c" "t" "c" "t" "g" "g" "g" "a" "g" "t" "c" "t" "t" "c" "a" "g"
    [4555] "a" "a" "t" "a" "g" "a" "a" "a" "c" "t" "a" "c" "c" "c" "a" "t" "c" "t"
    [4573] "c" "a" "a" "g" "a" "g" "g" "a" "g" "c" "t" "c" "a" "t" "t" "a" "a" "g"
    [4591] "g" "t" "t" "g" "t" "t" "g" "a" "t" "g" "t" "g" "g" "a" "g" "g" "a" "g"
    [4609] "c" "a" "a" "c" "a" "g" "c" "t" "g" "g" "a" "a" "g" "a" "g" "t" "c" "t"
    [4627] "g" "g" "g" "c" "c" "a" "c" "a" "c" "g" "a" "t" "t" "t" "g" "a" "c" "g"
    [4645] "g" "a" "a" "a" "c" "a" "t" "c" "t" "t" "a" "c" "t" "t" "g" "c" "c" "a"
    [4663] "a" "g" "g" "c" "a" "a" "g" "a" "t" "c" "t" "a" "g" "a" "g" "g" "g" "a"
    [4681] "a" "c" "c" "c" "c" "t" "t" "a" "c" "c" "t" "g" "g" "a" "a" "t" "c" "t"
    [4699] "g" "g" "a" "a" "t" "c" "a" "g" "c" "c" "t" "c" "t" "t" "c" "t" "c" "t"
    [4717] "g" "a" "t" "g" "a" "c" "c" "c" "t" "g" "a" "a" "t" "c" "t" "g" "a" "t"
    [4735] "c" "c" "t" "t" "c" "t" "g" "a" "a" "g" "a" "c" "a" "g" "a" "g" "c" "c"
    [4753] "c" "c" "a" "g" "a" "g" "t" "c" "a" "g" "c" "t" "c" "g" "t" "g" "t" "t"
    [4771] "g" "g" "c" "a" "a" "c" "a" "t" "a" "c" "c" "a" "t" "c" "t" "t" "c" "a"
    [4789] "a" "c" "c" "t" "c" "t" "g" "c" "a" "t" "t" "g" "a" "a" "a" "g" "t" "t"
    [4807] "c" "c" "c" "c" "a" "a" "t" "t" "g" "a" "a" "a" "g" "t" "t" "g" "c" "a"
    [4825] "g" "a" "a" "t" "c" "t" "g" "c" "c" "c" "a" "g" "a" "g" "t" "c" "c" "a"
    [4843] "g" "c" "t" "g" "c" "t" "g" "c" "t" "c" "a" "t" "a" "c" "t" "a" "c" "t"
    [4861] "g" "a" "t" "a" "c" "t" "g" "c" "t" "g" "g" "g" "t" "a" "t" "a" "a" "t"
    [4879] "g" "c" "a" "a" "t" "g" "g" "a" "a" "g" "a" "a" "a" "g" "t" "g" "t" "g"
    [4897] "a" "g" "c" "a" "g" "g" "g" "a" "g" "a" "a" "g" "c" "c" "a" "g" "a" "a"
    [4915] "t" "t" "g" "a" "c" "a" "g" "c" "t" "t" "c" "a" "a" "c" "a" "g" "a" "a"
    [4933] "a" "g" "g" "g" "t" "c" "a" "a" "c" "a" "a" "a" "a" "g" "a" "a" "t" "g"
    [4951] "t" "c" "c" "a" "t" "g" "g" "t" "g" "g" "t" "g" "t" "c" "t" "g" "g" "c"
    [4969] "c" "t" "g" "a" "c" "c" "c" "c" "a" "g" "a" "a" "g" "a" "a" "t" "t" "t"
    [4987] "a" "t" "g" "c" "t" "c" "g" "t" "g" "t" "a" "c" "a" "a" "g" "t" "t" "t"
    [5005] "g" "c" "c" "a" "g" "a" "a" "a" "a" "c" "a" "c" "c" "a" "c" "a" "t" "c"
    [5023] "a" "c" "t" "t" "t" "a" "a" "c" "t" "a" "a" "t" "c" "t" "a" "a" "t" "t"
    [5041] "a" "c" "t" "g" "a" "a" "g" "a" "g" "a" "c" "t" "a" "c" "t" "c" "a" "t"
    [5059] "g" "t" "t" "g" "t" "t" "a" "t" "g" "a" "a" "a" "a" "c" "a" "g" "a" "t"
    [5077] "g" "c" "t" "g" "a" "g" "t" "t" "t" "g" "t" "g" "t" "g" "t" "g" "a" "a"
    [5095] "c" "g" "g" "a" "c" "a" "c" "t" "g" "a" "a" "a" "t" "a" "t" "t" "t" "t"
    [5113] "c" "t" "a" "g" "g" "a" "a" "t" "t" "g" "c" "g" "g" "g" "a" "g" "g" "a"
    [5131] "a" "a" "a" "t" "g" "g" "g" "t" "a" "g" "t" "t" "a" "g" "c" "t" "a" "t"
    [5149] "t" "t" "c" "t" "g" "g" "g" "t" "g" "a" "c" "c" "c" "a" "g" "t" "c" "t"
    [5167] "a" "t" "t" "a" "a" "a" "g" "a" "a" "a" "g" "a" "a" "a" "a" "a" "t" "g"
    [5185] "c" "t" "g" "a" "a" "t" "g" "a" "g" "c" "a" "t" "g" "a" "t" "t" "t" "t"
    [5203] "g" "a" "a" "g" "t" "c" "a" "g" "a" "g" "g" "a" "g" "a" "t" "g" "t" "g"
    [5221] "g" "t" "c" "a" "a" "t" "g" "g" "a" "a" "g" "a" "a" "a" "c" "c" "a" "c"
    [5239] "c" "a" "a" "g" "g" "t" "c" "c" "a" "a" "a" "g" "c" "g" "a" "g" "c" "a"
    [5257] "a" "g" "a" "g" "a" "a" "t" "c" "c" "c" "a" "g" "g" "a" "c" "a" "g" "a"
    [5275] "a" "a" "g" "a" "t" "c" "t" "t" "c" "a" "g" "g" "g" "g" "g" "c" "t" "a"
    [5293] "g" "a" "a" "a" "t" "c" "t" "g" "t" "t" "g" "c" "t" "a" "t" "g" "g" "g"
    [5311] "c" "c" "c" "t" "t" "c" "a" "c" "c" "a" "a" "c" "a" "t" "g" "c" "c" "c"
    [5329] "a" "c" "a" "g" "a" "t" "c" "a" "a" "c" "t" "g" "g" "a" "a" "t" "g" "g"
    [5347] "a" "t" "g" "g" "t" "a" "c" "a" "g" "c" "t" "g" "t" "g" "t" "g" "g" "t"
    [5365] "g" "c" "t" "t" "c" "t" "g" "t" "g" "g" "t" "g" "a" "a" "g" "g" "a" "g"
    [5383] "c" "t" "t" "t" "c" "a" "t" "c" "a" "t" "t" "c" "a" "c" "c" "c" "t" "t"
    [5401] "g" "g" "c" "a" "c" "a" "g" "g" "t" "g" "t" "c" "c" "a" "c" "c" "c" "a"
    [5419] "a" "t" "t" "g" "t" "g" "g" "t" "t" "g" "t" "g" "c" "a" "g" "c" "c" "a"
    [5437] "g" "a" "t" "g" "c" "c" "t" "g" "g" "a" "c" "a" "g" "a" "g" "g" "a" "c"
    [5455] "a" "a" "t" "g" "g" "c" "t" "t" "c" "c" "a" "t" "g" "c" "a" "a" "t" "t"
    [5473] "g" "g" "g" "c" "a" "g" "a" "t" "g" "t" "g" "t" "g" "a" "g" "g" "c" "a"
    [5491] "c" "c" "t" "g" "t" "g" "g" "t" "g" "a" "c" "c" "c" "g" "a" "g" "a" "g"
    [5509] "t" "g" "g" "g" "t" "g" "t" "t" "g" "g" "a" "c" "a" "g" "t" "g" "t" "a"
    [5527] "g" "c" "a" "c" "t" "c" "t" "a" "c" "c" "a" "g" "t" "g" "c" "c" "a" "g"
    [5545] "g" "a" "g" "c" "t" "g" "g" "a" "c" "a" "c" "c" "t" "a" "c" "c" "t" "g"
    [5563] "a" "t" "a" "c" "c" "c" "c" "a" "g" "a" "t" "c" "c" "c" "c" "c" "a" "c"
    [5581] "a" "g" "c" "c" "a" "c" "t" "a" "c" "t" "g" "a"
    attr(,"name")
    [1] "JN686490"
    attr(,"Annot")
    [1] ">JN686490"
    attr(,"class")
    [1] "SeqFastadna"
    
    $JX480460.BRCA1
     [1] "c" "a" "a" "c" "c" "t" "g" "a" "g" "g" "t" "c" "t" "a" "t" "a" "a" "a" "c"
    [20] "a" "a" "a" "g" "t" "c" "t" "t" "c" "c" "t" "g" "g" "a" "a" "g" "t" "a" "a"
    [39] "t" "t" "g" "t" "a" "a" "g" "c" "a" "t" "c" "c" "t" "g" "a" "a" "a" "t" "a"
    [58] "a" "a" "a" "a" "a" "g" "c" "a" "a" "g" "a" "a" "t" "a" "t" "g" "a" "g" "t"
    [77] "a" "g"
    attr(,"name")
    [1] "JX480460.BRCA1"
    attr(,"Annot")
    [1] ">JX480460.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $JX480461
      [1] "c" "c" "t" "t" "c" "c" "t" "t" "g" "g" "a" "a" "a" "c" "c" "a" "g" "t"
     [19] "c" "t" "c" "a" "g" "t" "g" "t" "c" "c" "a" "a" "c" "t" "c" "t" "c" "t"
     [37] "a" "a" "c" "c" "t" "t" "g" "g" "a" "a" "c" "t" "g" "t" "g" "a" "g" "a"
     [55] "a" "c" "t" "c" "t" "g" "a" "g" "g" "a" "c" "a" "a" "a" "g" "c" "a" "g"
     [73] "c" "g" "g" "a" "t" "a" "c" "a" "a" "c" "c" "t" "c" "a" "a" "a" "a" "g"
     [91] "a" "c" "g" "t" "c" "t" "g" "t" "c" "t" "a" "c" "a" "t" "t" "g" "a" "a"
    [109] "t" "t" "g" "g"
    attr(,"name")
    [1] "JX480461"
    attr(,"Annot")
    [1] ">JX480461"
    attr(,"class")
    [1] "SeqFastadna"
    
    $JX480462.BRCA1
      [1] "g" "g" "a" "a" "c" "c" "c" "c" "t" "t" "a" "c" "c" "t" "g" "g" "a" "a"
     [19] "t" "c" "t" "g" "g" "a" "a" "t" "c" "a" "g" "c" "c" "t" "c" "t" "t" "c"
     [37] "t" "c" "t" "g" "a" "t" "g" "a" "c" "c" "c" "t" "g" "a" "a" "t" "c" "t"
     [55] "g" "a" "t" "c" "c" "t" "t" "c" "t" "g" "a" "a" "g" "a" "c" "a" "g" "a"
     [73] "g" "c" "c" "c" "c" "a" "g" "a" "g" "t" "c" "a" "g" "c" "t" "c" "g" "t"
     [91] "g" "t" "t" "g" "g" "c" "a" "a" "c" "a" "t" "a" "c" "c" "a" "t" "c" "t"
    [109] "t" "c" "a" "a" "c" "c" "t" "c" "t" "g" "c" "a" "t" "t" "g" "a" "a" "a"
    [127] "g" "t" "t" "c" "c" "c" "c" "a" "a" "t" "t" "g" "a" "a" "a" "g" "t" "t"
    [145] "g" "c" "a" "g" "a" "a" "t" "c" "t" "g" "c" "c" "c" "a" "g" "a" "g" "t"
    [163] "c" "c" "a" "g" "c" "t" "g" "c" "t" "g" "c" "t" "c" "a" "t" "a" "c" "t"
    [181] "a" "c" "t" "g" "a" "t" "a" "c" "t" "g" "c" "t" "g" "g" "g" "t" "a" "t"
    [199] "a" "a" "t" "g" "c" "a" "a" "t" "g" "g" "a" "a" "g" "a" "a" "a" "g" "t"
    [217] "g" "t" "g" "a" "g" "c" "a" "g" "g" "t" "a" "g"
    attr(,"name")
    [1] "JX480462.BRCA1"
    attr(,"Annot")
    [1] ">JX480462.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $JX480463.BRCA1
     [1] "c" "t" "g" "a" "g" "t" "g" "a" "c" "a" "a" "g" "g" "a" "a" "t" "t" "g" "g"
    [20] "t" "t" "t" "c" "a" "g" "a" "t" "g" "a" "t" "g" "a" "a" "g" "a" "a" "a" "g"
    [39] "a" "g" "g" "a" "a" "c" "g" "g" "g" "g" "c" "t" "t" "g" "g" "a" "a" "g" "a"
    [58] "a" "a" "a" "t" "a" "a"
    attr(,"name")
    [1] "JX480463.BRCA1"
    attr(,"Annot")
    [1] ">JX480463.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $JX480464.BRCA1
      [1] "t" "g" "g" "t" "t" "t" "t" "c" "t" "c" "c" "t" "t" "c" "c" "a" "t" "t"
     [19] "t" "a" "t" "c" "t" "t" "t" "c" "t" "a" "g" "g" "t" "c" "a" "t" "c" "c"
     [37] "c" "c" "t" "t" "c" "t" "a" "a" "a" "t" "g" "c" "c" "c" "a" "t" "c" "a"
     [55] "t" "t" "a" "g" "a" "t" "g" "a" "t" "a" "g" "g" "t" "g" "g" "t" "a" "c"
     [73] "a" "t" "g" "c" "a" "c" "a" "g" "t" "t" "g" "c" "t" "c" "t" "g" "g" "g"
     [91] "a" "g" "t" "c" "t" "t" "c" "a" "g" "a" "a" "t" "a" "g" "a" "a" "a" "c"
    [109] "t" "a" "c" "c" "c" "a" "t" "c" "t" "c" "a" "a" "g" "a" "g" "g" "a" "g"
    [127] "c" "t" "c" "a" "t" "t" "a" "a" "g" "g" "t" "t" "g" "t" "t" "g" "a" "t"
    [145] "g" "t" "g" "g" "a" "g" "g" "a" "g" "c" "a" "a" "c" "a" "g" "c" "t" "g"
    [163] "g" "a" "a" "g" "a" "g" "t" "c" "t" "g" "g" "g" "c" "c" "a" "c" "a" "c"
    [181] "g" "a" "t" "t" "t" "g" "a" "c" "g" "g" "a" "a" "a" "c" "a" "t" "c" "t"
    [199] "t" "a" "g"
    attr(,"name")
    [1] "JX480464.BRCA1"
    attr(,"Annot")
    [1] ">JX480464.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $JX480465.BRCA1
     [1] "g" "g" "a" "a" "c" "c" "c" "c" "t" "t" "a" "c" "c" "t" "g" "g" "a" "a" "a"
    [20] "t" "c" "t" "g" "g" "a" "a" "t" "c" "a" "g" "c" "c" "t" "c" "t" "t" "c" "t"
    [39] "c" "t" "g" "a"
    attr(,"name")
    [1] "JX480465.BRCA1"
    attr(,"Annot")
    [1] ">JX480465.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $JX480466.BRCA1
      [1] "g" "a" "a" "a" "t" "g" "g" "g" "a" "a" "a" "t" "g" "a" "g" "a" "a" "c"
     [19] "a" "t" "t" "c" "c" "a" "a" "g" "t" "a" "c" "a" "g" "t" "g" "a" "g" "c"
     [37] "a" "c" "a" "a" "t" "t" "a" "g" "c" "c" "g" "t" "a" "a" "t" "a" "a" "c"
     [55] "a" "t" "t" "a" "g" "a" "g" "a" "a" "a" "a" "t" "g" "t" "t" "t" "t" "t"
     [73] "a" "a" "a" "g" "a" "a" "g" "c" "c" "a" "g" "c" "t" "c" "a" "a" "g" "c"
     [91] "a" "a" "t" "a" "t" "t" "a" "a" "t" "g" "a" "a" "g" "t" "a" "g" "g" "t"
    [109] "t" "c" "c" "a" "g" "t" "a" "c" "t" "a" "a" "t" "g" "a" "a" "g" "t" "g"
    [127] "g" "g" "c" "t" "c" "c" "a" "g" "t" "a" "t" "t" "a" "a" "t" "g" "a" "a"
    [145] "a" "t" "a" "g" "g" "t" "t" "c" "c" "a" "g" "t" "g" "a" "t" "g" "a" "a"
    [163] "a" "a" "c" "a" "t" "t" "c" "a" "a" "g" "c" "a" "g" "a" "a" "t" "a" "g"
    attr(,"name")
    [1] "JX480466.BRCA1"
    attr(,"Annot")
    [1] ">JX480466.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $JX480467
     [1] "a" "t" "g" "t" "c" "c" "a" "t" "g" "g" "t" "g" "g" "t" "g" "t" "c" "t" "g"
    [20] "g" "c" "c" "t" "g" "a" "c" "c" "c" "c" "a" "g" "a" "a" "g" "a" "a" "t" "t"
    [39] "t" "a" "t" "c" "t" "t" "c" "a" "g" "g" "g" "g" "g" "c" "t" "a" "g" "a" "a"
    [58] "a" "t" "c" "t" "g" "t" "t" "g" "c" "t" "a" "t" "g" "g" "g" "c" "c" "c" "t"
    [77] "t" "c" "a" "c" "c" "a" "a" "c" "a" "t" "g" "c" "c" "c"
    attr(,"name")
    [1] "JX480467"
    attr(,"Annot")
    [1] ">JX480467"
    attr(,"class")
    [1] "SeqFastadna"
    
    $KJ625149.BRCA1
     [1] "a" "t" "g" "g" "a" "t" "t" "t" "a" "t" "c" "t" "g" "c" "t" "c" "t" "t" "c"
    [20] "g" "c" "g" "t" "t" "g" "a" "a" "g" "a" "a" "g" "t" "a" "t" "a" "a"
    attr(,"name")
    [1] "KJ625149.BRCA1"
    attr(,"Annot")
    [1] ">KJ625149.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $KJ625150.BRCA1
     [1] "a" "t" "g" "g" "a" "t" "t" "t" "a" "t" "c" "t" "g" "c" "t" "c" "t" "t" "c"
    [20] "g" "c" "g" "t" "t" "g" "a" "a" "g" "a" "a" "g" "t" "a" "c" "a" "a" "a" "a"
    [39] "t" "g" "t" "c" "a" "t" "t" "a" "a" "t" "g" "c" "t" "a" "t" "g" "c" "a" "g"
    [58] "a" "a" "a" "a" "t" "c" "t" "t" "a" "a" "g" "a" "g" "t" "g" "t" "c" "c" "c"
    [77] "a" "t" "c" "t" "g"
    attr(,"name")
    [1] "KJ625150.BRCA1"
    attr(,"Annot")
    [1] ">KJ625150.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $KJ625151.BRCA1
     [1] "g" "a" "g" "c" "c" "t" "a" "c" "a" "a" "g" "a" "a" "a" "g" "t" "a" "c" "g"
    [20] "a" "g" "a" "t" "t" "t" "a" "g" "t" "c" "a" "a" "c" "t" "t" "g" "t" "t" "g"
    [39] "a" "a" "g" "a" "g" "c" "t" "a" "t" "t" "g" "a" "a" "a" "a" "t" "c" "a" "t"
    [58] "t" "t" "g" "t" "g" "c" "t" "t" "t" "t" "t" "a" "g"
    attr(,"name")
    [1] "KJ625151.BRCA1"
    attr(,"Annot")
    [1] ">KJ625151.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $KJ625152.BRCA1
     [1] "c" "a" "g" "g" "a" "a" "a" "c" "c" "a" "g" "t" "c" "t" "c" "a" "g" "t" "g"
    [20] "t" "c" "c" "a" "a" "c" "t" "c" "t" "a" "a"
    attr(,"name")
    [1] "KJ625152.BRCA1"
    attr(,"Annot")
    [1] ">KJ625152.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $KJ625153.BRCA1
     [1] "c" "a" "g" "g" "a" "a" "a" "c" "c" "a" "g" "t" "c" "t" "c" "a" "g" "t" "g"
    [20] "t" "c" "c" "a" "a" "c" "t" "c" "t" "c" "t" "a" "a" "c" "c" "t" "t" "g" "g"
    [39] "a" "a" "c" "t" "g" "t" "g" "a" "g" "a" "a" "c" "t" "c" "t" "g" "a" "g" "g"
    [58] "a" "c" "a" "t" "a" "g"
    attr(,"name")
    [1] "KJ625153.BRCA1"
    attr(,"Annot")
    [1] ">KJ625153.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $KJ625154.BRCA1
     [1] "c" "a" "g" "g" "a" "a" "a" "c" "c" "a" "g" "t" "c" "t" "c" "a" "g" "t" "g"
    [20] "t" "c" "c" "a" "a" "c" "t" "c" "t" "c" "t" "a" "a" "c" "c" "t" "t" "g" "g"
    [39] "a" "a" "c" "t" "g" "t" "g" "a" "g" "a" "a" "c" "t" "c" "t" "g" "a" "g" "g"
    [58] "a" "c" "a" "a" "a" "g" "c" "a" "g" "c" "g" "g" "a" "t" "a" "t" "a" "a"
    attr(,"name")
    [1] "KJ625154.BRCA1"
    attr(,"Annot")
    [1] ">KJ625154.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $KJ625155.BRCA1
      [1] "c" "t" "g" "c" "t" "t" "g" "t" "g" "a" "a" "t" "t" "t" "t" "c" "t" "g"
     [19] "a" "g" "a" "c" "g" "g" "a" "t" "g" "t" "a" "a" "c" "a" "a" "a" "t" "a"
     [37] "c" "t" "g" "a" "a" "c" "a" "t" "c" "a" "t" "c" "a" "a" "c" "c" "c" "a"
     [55] "g" "t" "a" "a" "t" "a" "a" "t" "g" "a" "t" "t" "t" "g" "a" "a" "c" "a"
     [73] "c" "c" "a" "c" "t" "g" "a" "g" "a" "a" "g" "c" "g" "t" "g" "c" "a" "g"
     [91] "c" "t" "g" "a" "g" "a" "g" "g" "c" "a" "t" "c" "c" "a" "g" "a" "a" "a"
    [109] "a" "g" "t" "a" "t" "c" "a" "g" "g" "g" "t" "a" "g" "t" "t" "c" "t" "g"
    [127] "t" "t" "t" "c" "a" "a" "a" "c" "t" "t" "g" "c" "a" "t" "g" "t" "g" "g"
    [145] "a" "g" "c" "c" "a" "t" "g" "t" "g" "g" "c" "a" "c" "a" "a" "a" "t" "a"
    [163] "c" "t" "c" "a" "t" "g" "c" "c" "a" "g" "c" "t" "c" "a" "t" "t" "a" "c"
    [181] "a" "g" "c" "a" "t" "g" "a" "g" "a" "a" "c" "a" "g" "c" "a" "g" "t" "t"
    [199] "t" "a" "t" "t" "a" "c" "t" "c" "a" "c" "t" "a" "a" "a" "g" "a" "c" "a"
    [217] "g" "a" "a" "t" "g" "a" "a" "t" "g" "t" "a" "g" "a" "a" "a" "a" "g" "g"
    [235] "c" "t" "g" "a" "a" "t" "t" "c" "t" "g" "t" "a" "a" "t" "a" "a" "a" "a"
    [253] "g" "c" "a" "a" "a" "c" "g" "g" "c" "t" "t" "a" "g" "c" "a" "a" "g" "g"
    [271] "a" "g" "c" "c" "a" "a" "c" "a" "t" "a" "a"
    attr(,"name")
    [1] "KJ625155.BRCA1"
    attr(,"Annot")
    [1] ">KJ625155.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $KJ625156.BRCA1
      [1] "t" "t" "a" "c" "t" "c" "a" "c" "t" "a" "a" "a" "g" "a" "c" "a" "g" "a"
     [19] "a" "t" "g" "a" "a" "t" "g" "t" "a" "g" "a" "a" "a" "a" "g" "g" "c" "t"
     [37] "g" "a" "a" "t" "t" "c" "t" "g" "t" "a" "a" "t" "a" "a" "a" "a" "g" "c"
     [55] "a" "a" "a" "c" "a" "g" "c" "c" "t" "g" "g" "c" "t" "t" "a" "g" "c" "a"
     [73] "a" "g" "g" "a" "g" "c" "c" "a" "a" "c" "a" "t" "a" "a" "c" "a" "g" "a"
     [91] "t" "g" "g" "c" "t" "g" "g" "a" "a" "g" "t" "a" "a" "g" "g" "a" "a" "a"
    [109] "c" "a" "t" "g" "t" "a" "a" "t" "g" "a" "t" "a" "g" "g" "c" "g" "g" "a"
    [127] "c" "t" "c" "c" "c" "a" "g" "c" "a" "c" "a" "g" "a" "a" "a" "a" "a" "a"
    [145] "a" "g" "g" "t" "a" "g"
    attr(,"name")
    [1] "KJ625156.BRCA1"
    attr(,"Annot")
    [1] ">KJ625156.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $KJ625157.BRCA1
      [1] "t" "t" "a" "c" "t" "c" "a" "c" "t" "a" "a" "a" "g" "a" "c" "a" "g" "a"
     [19] "a" "t" "g" "a" "a" "t" "g" "t" "a" "g" "a" "a" "a" "a" "g" "g" "c" "t"
     [37] "g" "a" "a" "t" "t" "c" "t" "g" "t" "a" "a" "t" "a" "a" "a" "a" "g" "c"
     [55] "a" "a" "a" "c" "a" "g" "c" "c" "t" "g" "g" "c" "t" "t" "a" "g" "c" "a"
     [73] "a" "g" "g" "a" "g" "c" "c" "a" "a" "c" "a" "t" "a" "a" "c" "a" "g" "a"
     [91] "t" "g" "g" "g" "c" "t" "g" "g" "a" "a" "g" "t" "a" "a" "g" "g" "a" "a"
    [109] "a" "c" "g" "t" "a" "a"
    attr(,"name")
    [1] "KJ625157.BRCA1"
    attr(,"Annot")
    [1] ">KJ625157.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $KJ625158.BRCA1
      [1] "t" "t" "a" "c" "t" "c" "a" "c" "t" "a" "a" "a" "g" "a" "c" "a" "g" "a"
     [19] "a" "t" "g" "a" "a" "t" "g" "t" "a" "g" "a" "a" "a" "a" "g" "g" "c" "t"
     [37] "g" "a" "a" "t" "t" "c" "t" "g" "t" "a" "a" "t" "a" "a" "a" "a" "g" "c"
     [55] "a" "a" "a" "c" "a" "g" "c" "c" "t" "g" "g" "c" "t" "t" "a" "g" "c" "a"
     [73] "a" "g" "g" "a" "g" "c" "c" "a" "a" "c" "a" "t" "a" "a" "c" "a" "g" "a"
     [91] "t" "g" "g" "g" "c" "t" "g" "g" "a" "a" "g" "t" "a" "a" "g" "g" "a" "a"
    [109] "a" "c" "a" "t" "g" "t" "a" "a" "t" "g" "a" "t" "a" "g" "g" "c" "g" "g"
    [127] "a" "c" "t" "c" "c" "c" "a" "g" "c" "a" "c" "a" "g" "a" "a" "a" "a" "a"
    [145] "a" "g" "g" "t" "a" "g"
    attr(,"name")
    [1] "KJ625158.BRCA1"
    attr(,"Annot")
    [1] ">KJ625158.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $KJ625159.BRCA1
      [1] "g" "t" "a" "t" "a" "a" "g" "c" "a" "a" "t" "a" "t" "g" "g" "a" "a" "c"
     [19] "t" "c" "g" "a" "a" "t" "t" "a" "a" "a" "t" "a" "t" "c" "c" "a" "c" "a"
     [37] "a" "t" "t" "c" "a" "a" "a" "a" "g" "c" "a" "c" "c" "t" "a" "a" "a" "a"
     [55] "a" "g" "a" "a" "t" "a" "g" "g" "c" "t" "g" "a" "g" "g" "a" "g" "g" "a"
     [73] "a" "g" "t" "c" "t" "t" "c" "t" "a" "c" "c" "a" "g" "g" "c" "a" "t" "a"
     [91] "t" "t" "c" "a" "t" "g" "c" "g" "c" "t" "t" "g" "a" "a" "c" "t" "a" "g"
    [109] "t" "a" "g" "t" "c" "a" "g" "t" "a" "g" "a" "a" "a" "t" "c" "t" "a" "a"
    [127] "g" "c" "c" "c" "a" "c" "c" "t" "a" "a" "t" "t" "g" "t" "a" "c" "t" "g"
    [145] "a" "a" "t" "t" "g" "c" "a" "a" "a" "t" "t" "g" "a" "t" "a" "g" "t" "t"
    [163] "g" "t" "t" "c" "t" "a" "g" "c" "a" "g" "t" "g" "a" "a" "g" "a" "g" "a"
    [181] "t" "a" "a" "a" "g" "a" "a" "a" "a" "a" "a" "a" "a" "a" "g" "t" "a" "c"
    [199] "a" "a" "c" "c" "a" "a" "a" "t" "g" "c" "c" "a" "g" "t" "c" "a" "g" "g"
    [217] "c" "a" "c" "a" "g" "c" "a" "g" "a" "a" "a" "c" "c" "t" "a" "c" "a" "a"
    [235] "c" "t" "c" "a" "t" "g" "g" "a" "a" "g" "g" "t" "a" "a"
    attr(,"name")
    [1] "KJ625159.BRCA1"
    attr(,"Annot")
    [1] ">KJ625159.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $KJ625160.BRCA1
      [1] "g" "g" "a" "a" "g" "t" "t" "a" "g" "c" "a" "c" "t" "c" "t" "a" "g" "g"
     [19] "g" "a" "a" "g" "g" "c" "a" "a" "a" "a" "a" "c" "a" "g" "a" "a" "c" "c"
     [37] "a" "a" "a" "t" "a" "a" "a" "t" "g" "t" "g" "t" "g" "a" "g" "t" "c" "a"
     [55] "g" "t" "g" "t" "g" "c" "a" "g" "c" "a" "t" "t" "t" "g" "a" "a" "a" "a"
     [73] "c" "c" "c" "c" "a" "a" "g" "g" "g" "a" "c" "t" "a" "a" "t" "t" "c" "a"
     [91] "t" "g" "g" "t" "t" "g" "t" "t" "c" "c" "a" "a" "a" "g" "a" "t" "a" "a"
    [109] "t" "a" "g" "a" "a" "a" "t" "g" "a" "c" "a" "c" "a" "g" "a" "a" "g" "g"
    [127] "c" "t" "t" "t" "a" "a" "g" "t" "a" "t" "c" "c" "a" "t" "t" "g" "g" "g"
    [145] "a" "c" "a" "t" "g" "a" "a" "g" "t" "t" "a" "a" "c" "c" "a" "c" "a" "g"
    [163] "t" "c" "g" "g" "g" "a" "a" "a" "c" "a" "a" "g" "c" "a" "t" "a" "g" "a"
    [181] "a" "a" "t" "g" "g" "a" "a" "g" "a" "a" "a" "g" "t" "g" "a" "a" "c" "t"
    [199] "t" "g" "a" "t" "g" "c" "t" "c" "a" "g" "t" "a" "t" "t" "t" "g" "c" "a"
    [217] "g" "a" "a" "t" "t" "t" "a" "c" "a" "t" "t" "c" "a" "a" "g" "g" "t" "t"
    [235] "t" "c" "a" "a" "a" "g" "c" "g" "c" "c" "a" "g" "t" "c" "a" "t" "t" "t"
    [253] "g" "c" "t" "c" "c" "g" "t" "t" "t" "t" "c" "a" "a" "a" "t" "c" "c" "a"
    [271] "g" "g" "a" "a" "a" "t" "g" "c" "a" "g" "a" "a" "g" "a" "g" "g" "a" "a"
    [289] "t" "g" "t" "g" "c" "a" "a" "c" "a" "t" "t" "c" "t" "c" "t" "g" "c" "c"
    [307] "c" "a" "c" "t" "c" "t" "g" "g" "g" "t" "c" "c" "t" "t" "a" "a"
    attr(,"name")
    [1] "KJ625160.BRCA1"
    attr(,"Annot")
    [1] ">KJ625160.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $KJ625161.BRCA1
      [1] "g" "g" "a" "a" "g" "t" "t" "a" "g" "c" "a" "c" "t" "c" "t" "a" "g" "g"
     [19] "g" "a" "a" "g" "g" "c" "a" "a" "a" "a" "a" "c" "a" "g" "a" "a" "c" "c"
     [37] "a" "a" "a" "t" "a" "a" "a" "t" "g" "t" "g" "t" "g" "a" "g" "t" "c" "a"
     [55] "g" "t" "g" "t" "g" "c" "a" "g" "c" "a" "t" "t" "t" "g" "a" "a" "a" "a"
     [73] "c" "c" "c" "c" "a" "a" "g" "g" "g" "a" "c" "t" "a" "a" "t" "t" "c" "a"
     [91] "t" "g" "g" "t" "t" "g" "t" "t" "c" "c" "a" "a" "a" "g" "a" "t" "a" "a"
    [109] "t" "a" "g" "a" "a" "a" "t" "g" "a" "c" "a" "c" "a" "g" "a" "a" "g" "g"
    [127] "c" "t" "t" "t" "a" "a" "g" "t" "a" "t" "c" "c" "a" "t" "t" "g" "g" "g"
    [145] "a" "c" "a" "t" "g" "a" "a" "g" "t" "t" "a" "a" "c" "c" "a" "c" "a" "g"
    [163] "t" "c" "g" "g" "g" "a" "a" "a" "c" "a" "a" "g" "c" "a" "t" "a" "g" "a"
    [181] "a" "a" "t" "g" "g" "a" "a" "g" "a" "a" "a" "g" "t" "g" "a" "a" "c" "t"
    [199] "t" "g" "a" "t" "g" "c" "t" "c" "a" "g" "t" "a" "t" "t" "t" "g" "c" "a"
    [217] "g" "a" "a" "t" "a" "c" "a" "t" "t" "c" "a" "a" "g" "g" "t" "t" "t" "c"
    [235] "a" "a" "a" "g" "c" "g" "c" "c" "a" "g" "t" "c" "a" "t" "t" "t" "g" "c"
    [253] "t" "c" "c" "g" "t" "t" "t" "t" "c" "a" "a" "a" "t" "c" "c" "a" "g" "g"
    [271] "a" "a" "a" "t" "g" "c" "a" "t" "a" "a"
    attr(,"name")
    [1] "KJ625161.BRCA1"
    attr(,"Annot")
    [1] ">KJ625161.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $KJ625162
      [1] "g" "g" "a" "a" "g" "t" "t" "a" "g" "c" "a" "c" "t" "c" "t" "a" "g" "g"
     [19] "g" "a" "a" "g" "g" "c" "a" "a" "a" "a" "a" "c" "a" "g" "a" "a" "c" "c"
     [37] "a" "a" "a" "t" "a" "a" "a" "t" "g" "t" "g" "t" "g" "a" "g" "t" "c" "a"
     [55] "g" "t" "g" "t" "g" "c" "a" "g" "c" "a" "t" "t" "t" "g" "a" "a" "a" "a"
     [73] "c" "c" "c" "c" "a" "a" "g" "g" "g" "a" "c" "t" "a" "a" "t" "t" "c" "a"
     [91] "t" "g" "g" "t" "t" "g" "t" "t" "c" "c" "a" "a" "a" "g" "a" "t" "a" "a"
    [109] "t" "a" "g" "a" "a" "a" "t" "g" "a" "c" "a" "c" "a" "g" "a" "a" "g" "g"
    [127] "c" "t" "t" "t" "a" "a" "g" "t" "a" "t" "c" "c" "a" "t" "t" "g" "g" "g"
    [145] "a" "c" "a" "t" "g" "a" "a" "g" "t" "t" "a" "a" "c" "c" "a" "c" "a" "g"
    [163] "t" "c" "g" "g" "g" "a" "a" "a" "c" "a" "a" "g" "c" "a" "t" "a" "g" "a"
    [181] "a" "a" "t" "g" "g" "a" "a" "g" "a" "a" "a" "g" "t" "g" "a" "a" "c" "t"
    [199] "t" "g" "a" "t" "g" "c" "t" "c" "a" "g" "t" "a" "t" "t" "t" "g" "c" "a"
    [217] "g" "a" "a" "t" "a" "c" "a" "t" "t" "c" "a" "a" "g" "g" "t" "t" "t" "c"
    [235] "a" "a" "a" "g" "c" "g" "c" "c" "a" "g" "t" "c" "a" "t" "t" "t" "g" "c"
    [253] "t" "c" "c" "g" "t" "t" "t" "t" "c" "a" "a" "a" "t" "c" "c" "a" "g" "g"
    [271] "a" "a" "a" "t" "g" "c" "a" "g" "a" "a" "g" "a" "g" "g" "a" "a" "t" "g"
    [289] "t" "g" "c" "a" "a" "c" "a" "t" "t" "c" "t" "c" "t" "g" "c" "c" "c" "a"
    [307] "c" "t" "c" "t" "g" "g" "g" "t" "c" "c" "t" "t" "a" "a" "a" "g" "a" "a"
    [325] "a" "c" "a" "a" "a" "g" "t" "c" "c" "a" "a" "a" "a" "t" "c" "a" "c" "t"
    [343] "t" "t" "t" "g" "a" "a" "t" "g" "t" "g" "a" "a" "c" "a" "a" "a" "a" "g"
    [361] "g" "a" "a" "g" "a" "a" "a" "a" "t" "c" "a" "a" "g" "g" "a" "a" "a" "g"
    [379] "a" "a" "t" "g" "a" "g" "t" "c" "t" "a" "a" "t" "a" "t" "c" "a" "a" "g"
    [397] "c" "c" "t" "g" "t" "a" "c" "a" "g" "a" "c" "a" "g" "t" "t" "a" "a" "t"
    [415] "a" "t" "c" "a" "c" "t" "g" "c"
    attr(,"name")
    [1] "KJ625162"
    attr(,"Annot")
    [1] ">KJ625162"
    attr(,"class")
    [1] "SeqFastadna"
    
    $KJ625163.BRCA1
      [1] "c" "c" "a" "c" "t" "t" "t" "t" "t" "c" "c" "c" "a" "t" "c" "a" "a" "g"
     [19] "t" "c" "a" "t" "t" "t" "g" "t" "t" "a" "a" "a" "a" "c" "t" "a" "a" "a"
     [37] "t" "g" "t" "a" "a" "g" "a" "a" "a" "a" "a" "t" "c" "t" "g" "c" "t" "a"
     [55] "g" "a" "g" "g" "a" "a" "a" "a" "c" "t" "t" "t" "g" "a" "g" "g" "a" "a"
     [73] "c" "a" "t" "t" "c" "a" "a" "t" "g" "t" "c" "a" "c" "c" "t" "g" "a" "a"
     [91] "a" "g" "a" "g" "a" "a" "a" "t" "g" "g" "g" "a" "a" "a" "t" "g" "a" "g"
    [109] "a" "a" "c" "a" "t" "t" "c" "c" "a" "a" "g" "t" "a" "c" "a" "g" "t" "g"
    [127] "a" "g" "c" "a" "c" "a" "a" "t" "t" "a" "g" "c" "c" "g" "t" "a" "a" "t"
    [145] "a" "a" "c" "a" "t" "t" "a" "g" "a" "g" "a" "a" "a" "a" "t" "g" "t" "t"
    [163] "t" "t" "t" "a" "a" "a" "g" "a" "a" "g" "c" "c" "a" "g" "c" "t" "c" "a"
    [181] "a" "g" "c" "a" "a" "t" "a" "t" "t" "a" "a" "t" "g" "a" "a" "g" "t" "a"
    [199] "g" "g" "t" "t" "c" "c" "a" "g" "t" "a" "c" "t" "a" "a" "t" "g" "a" "a"
    [217] "g" "t" "g" "g" "g" "c" "t" "c" "c" "a" "g" "t" "a" "t" "t" "a" "a" "t"
    [235] "g" "a" "a" "a" "t" "a" "g" "g" "t" "t" "c" "c" "a" "g" "t" "g" "a" "t"
    [253] "g" "a" "a" "a" "a" "c" "a" "t" "t" "c" "a" "a" "g" "c" "a" "g" "a" "a"
    [271] "c" "t" "a" "g" "g" "t" "a" "g" "a" "a" "a" "c" "a" "g" "a" "g" "g" "g"
    [289] "c" "c" "a" "a" "a" "a" "t" "t" "g" "a" "a" "t" "g" "c" "t" "a" "t" "g"
    [307] "c" "t" "t" "a" "g" "a" "t" "t" "a" "g" "g" "g" "g" "t" "t" "t" "t" "g"
    [325] "c" "a" "a" "c" "c" "t" "g" "a" "g" "g" "t" "c" "t" "a" "t" "a" "a" "a"
    [343] "t" "a" "a"
    attr(,"name")
    [1] "KJ625163.BRCA1"
    attr(,"Annot")
    [1] ">KJ625163.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $KJ625164.BRCA1
      [1] "c" "a" "a" "a" "a" "t" "t" "g" "a" "a" "t" "g" "c" "t" "a" "t" "g" "c"
     [19] "t" "t" "a" "g" "a" "t" "t" "a" "g" "g" "g" "g" "t" "t" "t" "t" "g" "c"
     [37] "a" "a" "c" "c" "t" "g" "a" "g" "g" "t" "c" "t" "a" "t" "a" "a" "a" "c"
     [55] "a" "a" "a" "g" "t" "c" "t" "t" "c" "c" "t" "g" "g" "a" "a" "g" "t" "a"
     [73] "a" "t" "t" "g" "t" "a" "a" "g" "c" "a" "t" "c" "c" "t" "g" "a" "a" "a"
     [91] "t" "a" "a" "a" "a" "a" "a" "g" "c" "a" "g" "a" "a" "t" "a" "t" "g" "a"
    [109] "a" "g" "a" "a" "g" "t" "a" "g"
    attr(,"name")
    [1] "KJ625164.BRCA1"
    attr(,"Annot")
    [1] ">KJ625164.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $KJ625165.BRCA1
      [1] "g" "c" "a" "a" "a" "a" "g" "c" "g" "t" "c" "c" "a" "g" "a" "a" "a" "g"
     [19] "g" "a" "g" "a" "g" "c" "t" "t" "a" "g" "c" "a" "g" "g" "a" "g" "t" "c"
     [37] "c" "t" "a" "g" "c" "c" "c" "t" "t" "t" "c" "a" "c" "c" "c" "a" "t" "a"
     [55] "c" "a" "c" "a" "t" "t" "t" "g" "g" "c" "t" "c" "a" "g" "g" "g" "t" "t"
     [73] "a" "c" "c" "g" "a" "a" "g" "a" "g" "g" "g" "g" "c" "c" "a" "a" "g" "a"
     [91] "a" "a" "t" "t" "a" "g" "a" "g" "t" "c" "c" "t" "c" "a" "g" "a" "a" "g"
    [109] "a" "g" "a" "a" "c" "t" "t" "a" "t" "c" "t" "a" "g" "t" "g" "a" "g" "g"
    [127] "a" "t" "g" "a" "a" "g" "a" "g" "c" "t" "t" "c" "c" "c" "t" "g" "c" "t"
    [145] "t" "c" "c" "a" "a" "c" "a" "c" "t" "t" "g" "t" "t" "a" "t" "t" "t" "g"
    [163] "g" "t" "a" "a" "a" "g" "t" "a" "a" "a" "c" "a" "a" "t" "a" "t" "a" "c"
    [181] "c" "t" "t" "c" "t" "c" "a" "g" "t" "c" "t" "a" "c" "t" "a" "g" "g" "c"
    [199] "a" "t" "a" "g" "c" "a" "c" "c" "g" "t" "t" "g" "c" "t" "a" "c" "c" "g"
    [217] "a" "g" "t" "g" "t" "c" "t" "a" "a" "g" "a" "a" "c" "a" "c" "a" "g" "a"
    [235] "g" "g" "a" "g" "a" "a" "t" "t" "t" "a" "t" "t" "a" "t" "c" "a" "t" "t"
    [253] "g" "a"
    attr(,"name")
    [1] "KJ625165.BRCA1"
    attr(,"Annot")
    [1] ">KJ625165.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $KJ625166.BRCA1
      [1] "g" "c" "a" "a" "a" "a" "g" "c" "g" "t" "c" "c" "a" "g" "a" "a" "a" "g"
     [19] "g" "a" "g" "a" "g" "c" "t" "t" "a" "g" "c" "a" "g" "g" "a" "g" "t" "c"
     [37] "c" "t" "a" "g" "c" "c" "c" "t" "t" "t" "c" "a" "c" "c" "c" "a" "t" "a"
     [55] "c" "a" "c" "a" "t" "t" "t" "g" "g" "c" "t" "c" "a" "g" "g" "g" "t" "t"
     [73] "a" "c" "c" "g" "a" "a" "g" "a" "g" "g" "g" "g" "c" "c" "a" "a" "g" "a"
     [91] "a" "a" "t" "t" "a" "g" "a" "g" "t" "c" "c" "t" "c" "a" "g" "a" "a" "g"
    [109] "a" "g" "a" "a" "c" "t" "t" "a" "t" "c" "t" "a" "g" "t" "g" "a" "g" "g"
    [127] "a" "t" "g" "a" "a" "g" "a" "g" "c" "t" "t" "c" "c" "c" "t" "g" "c" "t"
    [145] "t" "c" "c" "a" "a" "c" "a" "c" "t" "t" "g" "t" "t" "a" "t" "t" "t" "g"
    [163] "g" "t" "a" "a" "a" "g" "t" "a" "a" "a" "c" "a" "a" "t" "a" "t" "a" "c"
    [181] "c" "t" "t" "c" "t" "c" "a" "g" "t" "c" "t" "a" "c" "t" "a" "g" "g" "c"
    [199] "a" "t" "a" "g" "c" "a" "c" "c" "g" "t" "t" "g" "c" "t" "a" "c" "c" "g"
    [217] "a" "g" "t" "g" "t" "c" "t" "g" "t" "c" "t" "a" "a" "g" "a" "a" "c" "a"
    [235] "g" "a" "g" "g" "a" "g" "a" "a" "t" "t" "t" "a" "t" "t" "a" "t" "c" "a"
    [253] "t" "t" "g" "a" "a" "g" "a" "a" "t" "a" "g"
    attr(,"name")
    [1] "KJ625166.BRCA1"
    attr(,"Annot")
    [1] ">KJ625166.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $KJ625167.BRCA1
      [1] "g" "c" "a" "a" "a" "a" "g" "c" "g" "t" "c" "c" "a" "g" "a" "a" "a" "g"
     [19] "g" "a" "g" "a" "g" "c" "t" "t" "a" "g" "c" "a" "g" "g" "a" "g" "t" "c"
     [37] "c" "t" "a" "g" "c" "c" "c" "t" "t" "t" "c" "a" "c" "c" "c" "a" "t" "a"
     [55] "c" "a" "c" "a" "t" "t" "t" "g" "g" "c" "t" "c" "a" "g" "g" "g" "t" "t"
     [73] "a" "c" "c" "g" "a" "a" "g" "a" "g" "g" "g" "g" "c" "c" "a" "a" "g" "a"
     [91] "a" "a" "t" "t" "a" "g" "a" "g" "t" "c" "c" "t" "c" "a" "g" "a" "a" "g"
    [109] "a" "g" "a" "a" "c" "t" "t" "a" "t" "c" "t" "a" "g" "t" "g" "a" "g" "g"
    [127] "a" "t" "g" "a" "a" "g" "a" "g" "c" "t" "t" "c" "c" "c" "t" "g" "c" "t"
    [145] "t" "c" "c" "a" "a" "c" "a" "c" "t" "t" "g" "t" "t" "a" "t" "t" "t" "g"
    [163] "g" "t" "a" "a" "a" "g" "t" "a" "a" "a" "c" "a" "a" "t" "a" "t" "a" "c"
    [181] "c" "t" "t" "c" "t" "c" "a" "g" "t" "c" "t" "a" "c" "t" "a" "g" "g" "c"
    [199] "a" "t" "a" "g" "c" "a" "c" "c" "g" "t" "t" "g" "c" "t" "a" "c" "c" "g"
    [217] "a" "g" "t" "g" "t" "c" "t" "g" "t" "c" "t" "a" "a" "g" "a" "a" "c" "a"
    [235] "c" "a" "g" "a" "g" "g" "a" "g" "a" "a" "t" "t" "t" "a" "t" "t" "a" "t"
    [253] "c" "a" "t" "t" "g" "a" "a" "g" "a" "a" "t" "a" "g" "c" "t" "t" "a" "a"
    [271] "a" "t" "g" "a" "c" "t" "g" "c" "a" "g" "t" "a" "a" "c" "c" "a" "g" "g"
    [289] "t" "a" "a" "t" "a" "t" "t" "g" "g" "c" "a" "a" "a" "g" "g" "c" "a" "t"
    [307] "c" "t" "c" "a" "g" "g" "a" "a" "c" "a" "t" "c" "a" "c" "c" "t" "t" "a"
    [325] "g" "g" "a" "a" "a" "c" "a" "a" "a" "a" "t" "g" "t" "t" "c" "t" "g" "c"
    [343] "t" "a" "g" "c" "t" "t" "g" "t" "t" "t" "t" "c" "t" "t" "c" "a" "c" "a"
    [361] "g" "t" "g" "c" "a" "g" "t" "g" "a" "a" "t" "t" "g" "g" "a" "a" "g" "a"
    [379] "c" "t" "t" "g" "a"
    attr(,"name")
    [1] "KJ625167.BRCA1"
    attr(,"Annot")
    [1] ">KJ625167.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $KJ625168.BRCA1
      [1] "g" "c" "a" "a" "a" "a" "g" "c" "g" "t" "c" "c" "a" "g" "a" "a" "a" "g"
     [19] "g" "a" "g" "a" "g" "c" "t" "t" "a" "g" "c" "a" "g" "g" "a" "g" "t" "c"
     [37] "c" "t" "a" "g" "c" "c" "c" "t" "t" "t" "c" "a" "c" "c" "c" "a" "t" "a"
     [55] "c" "a" "c" "a" "t" "t" "t" "g" "g" "c" "t" "c" "a" "g" "g" "g" "t" "t"
     [73] "a" "c" "c" "g" "a" "a" "g" "a" "g" "g" "g" "g" "c" "c" "a" "a" "g" "a"
     [91] "a" "a" "t" "t" "a" "g" "a" "g" "t" "c" "c" "t" "c" "a" "g" "a" "a" "g"
    [109] "a" "g" "a" "a" "c" "t" "t" "a" "t" "c" "t" "a" "g" "t" "g" "a" "g" "g"
    [127] "a" "t" "g" "a" "a" "g" "a" "g" "c" "t" "t" "c" "c" "c" "t" "g" "c" "t"
    [145] "t" "c" "c" "a" "a" "c" "a" "c" "t" "t" "g" "t" "t" "a" "t" "t" "t" "g"
    [163] "g" "t" "a" "a" "a" "g" "t" "a" "a" "a" "c" "a" "a" "t" "a" "t" "a" "c"
    [181] "c" "t" "t" "c" "t" "c" "a" "g" "t" "c" "t" "a" "c" "t" "a" "g" "g" "c"
    [199] "a" "t" "a" "g" "c" "a" "c" "c" "g" "t" "t" "g" "c" "t" "a" "c" "c" "g"
    [217] "a" "g" "t" "g" "t" "c" "t" "g" "t" "c" "t" "a" "a" "g" "a" "a" "c" "a"
    [235] "c" "a" "g" "a" "g" "g" "a" "g" "a" "a" "t" "t" "t" "a" "t" "t" "a" "t"
    [253] "c" "a" "t" "t" "g" "a" "a" "g" "a" "a" "t" "a" "g" "c" "t" "t" "a" "a"
    [271] "a" "t" "g" "a" "c" "t" "g" "c" "a" "g" "t" "a" "a" "c" "c" "a" "g" "g"
    [289] "t" "a" "a" "t" "a" "t" "t" "g" "g" "c" "a" "a" "a" "g" "g" "c" "a" "t"
    [307] "c" "t" "c" "a" "g" "g" "a" "a" "c" "a" "t" "c" "a" "c" "c" "t" "t" "a"
    [325] "g" "t" "g" "a" "g" "g" "a" "a" "a" "c" "a" "a" "t" "g" "t" "t" "c" "t"
    [343] "g" "c" "t" "a" "g"
    attr(,"name")
    [1] "KJ625168.BRCA1"
    attr(,"Annot")
    [1] ">KJ625168.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $KJ625169.BRCA1
      [1] "g" "t" "a" "a" "t" "a" "t" "t" "g" "g" "c" "a" "a" "a" "g" "g" "c" "a"
     [19] "t" "c" "t" "c" "a" "g" "g" "a" "a" "c" "a" "t" "c" "a" "c" "c" "t" "t"
     [37] "a" "g" "t" "g" "a" "g" "g" "a" "a" "a" "c" "a" "a" "a" "a" "t" "g" "t"
     [55] "t" "c" "t" "g" "c" "t" "a" "g" "c" "t" "t" "g" "t" "t" "t" "t" "c" "t"
     [73] "t" "c" "a" "c" "a" "g" "t" "g" "c" "a" "g" "t" "g" "a" "a" "t" "t" "g"
     [91] "g" "a" "a" "g" "a" "c" "t" "t" "g" "a" "c" "t" "g" "c" "a" "a" "a" "t"
    [109] "a" "c" "a" "a" "a" "c" "a" "c" "c" "c" "a" "g" "g" "a" "t" "c" "c" "t"
    [127] "t" "t" "c" "t" "t" "g" "a" "t" "t" "g" "g" "t" "t" "c" "t" "t" "c" "c"
    [145] "a" "a" "a" "c" "a" "a" "a" "t" "g" "a" "g" "g" "c" "a" "t" "c" "a" "g"
    [163] "t" "c" "t" "g" "a" "a" "a" "g" "c" "c" "a" "g" "g" "g" "a" "g" "t" "t"
    [181] "g" "g" "t" "c" "t" "g" "a" "g" "t" "g" "a" "c" "a" "a" "g" "g" "a" "a"
    [199] "t" "t" "g" "g" "t" "t" "t" "c" "a" "g" "a" "t" "g" "a" "t" "g" "a" "a"
    [217] "g" "a" "a" "a" "g" "a" "g" "g" "a" "a" "c" "g" "g" "g" "c" "t" "t" "g"
    [235] "g" "a" "a" "g" "a" "a" "a" "a" "t" "a" "a" "g" "a" "a" "g" "a" "g" "c"
    [253] "a" "a" "a" "g" "c" "a" "t" "g" "g" "a" "t" "t" "c" "a" "a" "a" "c" "t"
    [271] "t" "a" "g"
    attr(,"name")
    [1] "KJ625169.BRCA1"
    attr(,"Annot")
    [1] ">KJ625169.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $KJ625170.BRCA1
     [1] "g" "t" "g" "a" "a" "g" "c" "a" "g" "c" "a" "t" "c" "t" "g" "g" "g" "t" "g"
    [20] "t" "g" "a" "g" "a" "g" "t" "g" "a" "a" "a" "c" "a" "a" "g" "c" "g" "t" "c"
    [39] "t" "c" "t" "g" "a" "a" "g" "a" "c" "t" "g" "c" "t" "g" "a"
    attr(,"name")
    [1] "KJ625170.BRCA1"
    attr(,"Annot")
    [1] ">KJ625170.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $KJ625171.BRCA1
     [1] "c" "a" "g" "t" "a" "t" "t" "a" "a" "c" "t" "t" "c" "a" "t" "a" "g"
    attr(,"name")
    [1] "KJ625171.BRCA1"
    attr(,"Annot")
    [1] ">KJ625171.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $KJ625172.BRCA1
      [1] "a" "g" "g" "g" "a" "a" "c" "c" "c" "c" "t" "t" "a" "c" "c" "t" "g" "g"
     [19] "a" "a" "t" "c" "t" "g" "g" "a" "a" "t" "c" "a" "g" "c" "c" "t" "c" "t"
     [37] "t" "c" "t" "c" "t" "g" "a" "t" "g" "a" "c" "c" "c" "t" "g" "a" "a" "t"
     [55] "c" "t" "g" "a" "t" "c" "c" "t" "t" "c" "t" "g" "a" "a" "g" "a" "c" "a"
     [73] "g" "a" "g" "c" "c" "c" "c" "a" "g" "a" "g" "t" "c" "a" "g" "c" "t" "c"
     [91] "g" "t" "g" "t" "t" "g" "g" "c" "a" "a" "c" "a" "t" "a" "t" "g" "c" "a"
    [109] "t" "t" "g" "a" "a" "a" "g" "t" "t" "c" "c" "c" "c" "a" "a" "t" "t" "g"
    [127] "a" "a" "a" "g" "t" "t" "g" "c" "a" "g" "a" "a" "t" "c" "t" "g" "c" "c"
    [145] "c" "a" "g" "a" "g" "t" "c" "c" "a" "g" "c" "t" "g" "c" "t" "g" "c" "t"
    [163] "c" "a" "t" "a" "c" "t" "a" "c" "t" "g" "a"
    attr(,"name")
    [1] "KJ625172.BRCA1"
    attr(,"Annot")
    [1] ">KJ625172.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $KJ625173.BRCA1
     [1] "a" "t" "c" "a" "a" "c" "t" "g" "g" "a" "a" "t" "g" "g" "a" "t" "g" "g" "t"
    [20] "a" "t" "a" "g"
    attr(,"name")
    [1] "KJ625173.BRCA1"
    attr(,"Annot")
    [1] ">KJ625173.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $KJ625174.BRCA1
     [1] "c" "a" "a" "t" "t" "g" "g" "g" "c" "a" "g" "a" "t" "g" "t" "g" "t" "g" "a"
    [20] "g" "g" "c" "a" "c" "c" "t" "g" "t" "g" "g" "t" "g" "a" "c" "c" "t" "g" "a"
    attr(,"name")
    [1] "KJ625174.BRCA1"
    attr(,"Annot")
    [1] ">KJ625174.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $KJ625175.BRCA1
      [1] "t" "c" "t" "g" "g" "a" "g" "t" "t" "g" "a" "t" "c" "a" "a" "g" "g" "a"
     [19] "a" "c" "c" "t" "g" "t" "c" "t" "c" "c" "a" "c" "a" "a" "a" "g" "t" "g"
     [37] "t" "g" "a" "c" "c" "a" "c" "a" "t" "a" "t" "t" "t" "t" "g" "c" "a" "a"
     [55] "a" "t" "t" "t" "t" "g" "c" "a" "t" "g" "c" "t" "g" "a" "a" "a" "c" "t"
     [73] "t" "c" "t" "c" "a" "a" "c" "c" "a" "g" "a" "a" "g" "a" "a" "a" "g" "g"
     [91] "g" "c" "c" "t" "t" "c" "a" "c" "a" "g" "t" "g" "t" "c" "c" "t" "t" "t"
    [109] "a" "t" "g" "a"
    attr(,"name")
    [1] "KJ625175.BRCA1"
    attr(,"Annot")
    [1] ">KJ625175.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $KJ625176.BRCA1
      [1] "t" "c" "t" "g" "g" "a" "g" "t" "t" "g" "a" "t" "c" "a" "a" "g" "g" "a"
     [19] "a" "c" "c" "t" "g" "t" "c" "t" "c" "c" "a" "c" "a" "a" "a" "g" "t" "g"
     [37] "t" "g" "a" "c" "c" "a" "c" "a" "t" "a" "t" "t" "t" "t" "g" "c" "a" "a"
     [55] "g" "a" "g" "c" "c" "t" "a" "c" "a" "a" "g" "a" "a" "a" "g" "t" "a" "c"
     [73] "g" "a" "g" "a" "t" "t" "t" "a" "g" "t" "c" "a" "a" "c" "t" "t" "g" "t"
     [91] "t" "g" "a" "a" "g" "a" "g" "c" "t" "a" "t" "t" "g" "a" "a" "a" "a" "t"
    [109] "c" "a" "t" "t" "t" "g" "t" "g" "c" "t" "t" "t" "t" "c" "a" "g" "c" "t"
    [127] "t" "g" "a" "c" "a" "c" "a" "g" "g" "t" "t" "t" "g" "g" "a" "g" "t"
    attr(,"name")
    [1] "KJ625176.BRCA1"
    attr(,"Annot")
    [1] ">KJ625176.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $KJ625176.PE2
      [1] "t" "c" "t" "g" "g" "a" "g" "t" "t" "g" "a" "t" "c" "a" "a" "g" "g" "a"
     [19] "a" "c" "c" "t" "g" "t" "c" "t" "c" "c" "a" "c" "a" "a" "a" "g" "t" "g"
     [37] "t" "g" "a" "c" "c" "a" "c" "a" "t" "a" "t" "t" "t" "t" "g" "c" "a" "a"
     [55] "a" "t" "t" "t" "t" "g" "c" "a" "t" "g" "c" "t" "g" "a" "a" "a" "c" "t"
     [73] "t" "c" "t" "c" "a" "a" "c" "c" "a" "g" "a" "a" "g" "a" "a" "a" "g" "g"
     [91] "g" "c" "c" "t" "t" "c" "a" "c" "a" "g" "t" "g" "t" "c" "c" "t" "t" "t"
    [109] "a" "t" "g" "a"
    attr(,"name")
    [1] "KJ625176.PE2"
    attr(,"Annot")
    [1] ">KJ625176.PE2"
    attr(,"class")
    [1] "SeqFastadna"
    
    $KJ625177.BRCA1
     [1] "a" "t" "g" "c" "t" "g" "a" "g" "t" "t" "t" "g" "t" "g" "t" "g" "t" "g" "a"
    [20] "a" "c" "g" "g" "a" "c" "a" "c" "t" "g" "a" "a" "a" "t" "a" "t" "t" "t" "t"
    [39] "c" "t" "a" "g" "g" "a" "a" "t" "t" "g" "c" "g" "g" "g" "a" "g" "g" "a" "a"
    [58] "a" "a" "t" "g" "g" "g" "t" "a" "g" "t" "t" "a" "g" "c" "t" "a" "t" "t" "t"
    [77] "c" "t" "c" "a" "t" "g" "a"
    attr(,"name")
    [1] "KJ625177.BRCA1"
    attr(,"Annot")
    [1] ">KJ625177.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $KJ625178.BRCA1
      [1] "a" "t" "c" "t" "t" "c" "a" "g" "g" "g" "g" "g" "c" "t" "a" "g" "a" "a"
     [19] "a" "t" "c" "t" "g" "t" "t" "g" "c" "t" "a" "t" "g" "g" "g" "c" "c" "c"
     [37] "t" "t" "c" "a" "c" "c" "a" "a" "c" "a" "t" "g" "c" "c" "c" "a" "c" "a"
     [55] "g" "g" "g" "t" "g" "t" "c" "c" "a" "c" "c" "c" "a" "a" "t" "t" "g" "t"
     [73] "g" "g" "t" "t" "g" "t" "g" "c" "a" "g" "c" "c" "a" "g" "a" "t" "g" "c"
     [91] "c" "t" "g" "g" "a" "c" "a" "g" "a" "g" "g" "a" "c" "a" "a" "t" "g" "g"
    [109] "c" "t" "t" "c" "c" "a" "t" "g"
    attr(,"name")
    [1] "KJ625178.BRCA1"
    attr(,"Annot")
    [1] ">KJ625178.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $KJ625179.BRCA1
      [1] "a" "t" "c" "t" "t" "c" "a" "g" "g" "g" "g" "g" "c" "t" "a" "g" "a" "a"
     [19] "a" "t" "c" "t" "g" "t" "t" "g" "c" "t" "a" "t" "g" "g" "g" "c" "c" "c"
     [37] "t" "t" "c" "a" "c" "c" "a" "a" "c" "a" "t" "g" "c" "c" "c" "a" "c" "a"
     [55] "g" "g" "g" "t" "g" "t" "c" "c" "a" "c" "c" "c" "a" "a" "t" "t" "g" "t"
     [73] "g" "g" "t" "t" "g" "t" "g" "c" "a" "g" "c" "c" "a" "g" "a" "t" "g" "c"
     [91] "c" "t" "g" "g" "a" "c" "a" "g" "a" "g" "g" "a" "c" "a" "a" "t" "g" "g"
    [109] "c" "t" "t" "c" "c" "a" "t" "g"
    attr(,"name")
    [1] "KJ625179.BRCA1"
    attr(,"Annot")
    [1] ">KJ625179.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $KM434065
      [1] "g" "a" "t" "c" "t" "g" "a" "a" "t" "g" "c" "t" "g" "a" "t" "c" "c" "c"
     [19] "c" "t" "g" "t" "g" "t" "g" "a" "g" "a" "g" "a" "a" "a" "a" "g" "a" "a"
     [37] "t" "g" "g" "a" "a" "t" "a" "a" "g" "c" "a" "g" "a" "a" "a" "c" "t" "g"
     [55] "c" "c" "a" "t" "g" "c" "t" "c" "a" "g" "a" "g" "a" "a" "t" "c" "c" "t"
     [73] "a" "g" "a" "g" "a" "t" "a" "c" "t" "g" "a" "a" "g" "a" "t" "g" "t" "t"
     [91] "c" "c" "t" "t" "g" "g" "a" "t" "a" "a" "c" "a" "c" "t" "a" "a" "a" "t"
    [109] "a" "g" "c" "a" "g" "c" "a" "t" "t" "c" "a" "g" "a" "a" "a" "g" "t" "t"
    [127] "a" "a" "t" "g" "a" "g" "t" "g" "g" "t" "t" "t" "t" "c" "c" "a" "g" "a"
    [145] "a" "g" "t" "g" "a" "t" "g" "a" "a" "c" "t" "g" "t" "t" "a" "g" "g" "t"
    [163] "t" "c" "t" "g" "a" "t" "g" "a" "c" "t" "c" "a" "c" "a" "t" "g" "a" "t"
    [181] "g" "g" "g" "g" "a" "g" "t" "c" "t" "g" "a" "a" "t" "c" "a" "a" "a" "t"
    [199] "g" "c" "c" "a" "a" "a" "g" "t" "a" "g" "c" "t" "g" "a" "t" "g" "t" "a"
    [217] "t" "t" "g" "g" "a" "c" "g" "t" "t" "c" "t" "a" "a" "a" "t" "g" "a" "g"
    [235] "g" "t" "a" "g" "a" "t" "g" "a" "a" "t" "a" "t" "t" "c" "t" "g" "g" "t"
    [253] "t" "c" "t" "t" "c" "a" "g" "a" "g" "a" "a" "a" "a" "t" "a" "g" "a" "c"
    [271] "t" "t" "a" "c" "t" "g" "g" "c" "c" "a" "g" "t" "g" "a" "t" "c" "c" "t"
    [289] "c" "a" "t" "g" "a" "g" "g" "c" "t" "t" "t" "a" "a" "t" "a" "t" "g" "t"
    [307] "a" "a" "a" "a" "g" "t" "g" "a" "a" "a" "g" "a" "g" "t" "t" "c" "a" "c"
    [325] "t" "c" "c" "a" "a" "a" "t" "c" "a" "g" "t" "a" "g" "a" "g" "a" "g" "t"
    [343] "a" "a" "t" "a" "t" "t" "g" "a" "a" "g" "a" "c" "a" "a" "a" "a" "t" "a"
    [361] "t" "t" "t" "g" "g" "g" "a" "a" "a" "a" "c" "c" "t" "a" "t" "c" "g" "g"
    [379] "a" "a" "g" "a" "a" "g" "g" "c" "a" "a" "g" "c" "c" "t" "c" "c" "c" "c"
    [397] "a" "a" "c" "t" "t" "a" "a" "g" "c" "c" "a" "t" "g" "t" "a" "a" "c" "t"
    [415] "g" "a" "a" "a" "a" "t" "c" "t" "a" "a" "t" "t" "a" "t" "a" "g" "g" "a"
    [433] "g" "c" "a" "t" "t" "t" "g" "t" "t" "a" "c" "t" "g" "a" "g" "c" "c" "a"
    [451] "c" "a" "g" "a" "t" "a" "a" "t" "a" "c" "a" "a" "g" "a" "g" "c" "g" "t"
    [469] "c" "c" "c" "c" "t" "c" "a" "c" "a" "a" "a" "t" "a" "a" "a" "t" "t" "a"
    [487] "a" "a" "g" "c" "g" "t" "a" "a" "a" "a" "g" "g" "a" "g" "a" "c" "c" "t"
    [505] "a" "c" "a" "t" "c" "a" "g" "g" "c" "c" "t" "t" "c" "a" "t" "c" "c" "t"
    [523] "g" "a" "g" "g" "a" "t" "t" "t" "t" "a" "t" "c" "a" "a" "g" "a" "a" "a"
    [541] "g" "c" "a" "g" "a" "t" "t" "g" "g" "g" "a" "a" "g" "t" "c" "c" "a" "a"
    [559] "a" "g" "a" "a" "t" "t" "c" "t" "g" "g" "a" "a" "t" "g" "g" "a" "t" "a"
    [577] "a"
    attr(,"name")
    [1] "KM434065"
    attr(,"Annot")
    [1] ">KM434065"
    attr(,"class")
    [1] "SeqFastadna"
    
    $KP255396.BRCA1
      [1] "t" "t" "a" "c" "t" "c" "a" "c" "t" "a" "a" "a" "g" "a" "c" "a" "g" "a"
     [19] "a" "t" "g" "a" "a" "t" "g" "t" "a" "g" "a" "a" "a" "a" "g" "g" "c" "t"
     [37] "g" "a" "a" "t" "t" "c" "t" "g" "t" "a" "a" "t" "a" "a" "a" "a" "g" "c"
     [55] "a" "a" "a" "c" "a" "g" "c" "c" "t" "g" "g" "c" "t" "t" "a" "g" "c" "a"
     [73] "a" "g" "g" "a" "g" "c" "c" "a" "a" "c" "a" "t" "a" "a" "c" "a" "t" "g"
     [91] "g" "g" "c" "t" "g" "g" "a" "a" "g" "t" "a" "a"
    attr(,"name")
    [1] "KP255396.BRCA1"
    attr(,"Annot")
    [1] ">KP255396.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $KP255397.BRCA1
      [1] "g" "t" "a" "t" "a" "a" "g" "c" "a" "a" "t" "a" "t" "g" "g" "a" "a" "c"
     [19] "t" "c" "g" "a" "a" "t" "t" "a" "a" "a" "t" "a" "t" "c" "c" "a" "c" "a"
     [37] "a" "t" "t" "c" "a" "a" "a" "a" "g" "c" "a" "c" "c" "t" "a" "a" "a" "a"
     [55] "a" "g" "a" "a" "t" "a" "g" "g" "c" "t" "g" "a" "g" "g" "a" "g" "g" "a"
     [73] "a" "g" "t" "c" "t" "t" "c" "t" "a" "c" "c" "a" "g" "g" "c" "a" "t" "a"
     [91] "t" "t" "c" "a" "t" "g" "c" "g" "c" "t" "t" "g" "a" "a" "c" "t" "a" "g"
    [109] "t" "t" "a" "g" "t" "a" "g" "t" "c" "a" "g" "t" "a" "g"
    attr(,"name")
    [1] "KP255397.BRCA1"
    attr(,"Annot")
    [1] ">KP255397.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $KP255398.BRCA1
      [1] "a" "c" "a" "g" "c" "a" "g" "a" "a" "a" "c" "c" "t" "a" "c" "a" "a" "c"
     [19] "t" "c" "a" "t" "g" "g" "a" "a" "g" "g" "t" "a" "a" "a" "g" "a" "a" "c"
     [37] "c" "t" "g" "c" "a" "a" "c" "t" "g" "g" "a" "g" "c" "c" "a" "a" "g" "a"
     [55] "a" "g" "a" "g" "t" "a" "a" "c" "a" "a" "g" "c" "c" "a" "a" "a" "t" "g"
     [73] "a" "a" "c" "a" "g" "a" "c" "a" "a" "g" "t" "a" "a" "a" "a" "g" "a" "c"
     [91] "a" "t" "g" "a" "c" "a" "g" "c" "g" "a" "t" "a" "c" "t" "t" "t" "c" "c"
    [109] "c" "a" "g" "a" "g" "c" "t" "g" "a" "a" "g" "t" "t" "a" "a" "c" "a" "a"
    [127] "a" "t" "g" "c" "a" "c" "c" "t" "g" "g" "t" "t" "c" "t" "t" "t" "t" "a"
    [145] "c" "t" "a" "a" "g" "t" "g" "t" "t" "c" "a" "a" "a" "t" "a" "c" "c" "a"
    [163] "g" "t" "g" "a" "a" "c" "t" "t" "a" "a" "a" "g" "a" "a" "t" "t" "t" "g"
    [181] "t" "a" "a" "t" "c" "c" "t" "a" "g" "c" "c" "t" "t" "c" "c" "a" "a" "g"
    [199] "a" "g" "a" "a" "g" "a" "a" "a" "a" "a" "g" "a" "a" "g" "a" "g" "a" "a"
    [217] "a" "c" "t" "a" "g"
    attr(,"name")
    [1] "KP255398.BRCA1"
    attr(,"Annot")
    [1] ">KP255398.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $KP255399.BRCA1
      [1] "a" "c" "a" "g" "c" "a" "g" "a" "a" "a" "c" "c" "t" "a" "c" "a" "a" "c"
     [19] "t" "c" "a" "t" "g" "g" "a" "a" "g" "g" "t" "a" "a" "a" "g" "a" "a" "c"
     [37] "c" "t" "g" "c" "a" "a" "c" "t" "g" "g" "a" "g" "c" "c" "a" "a" "g" "a"
     [55] "a" "g" "a" "g" "t" "a" "a" "c" "a" "a" "g" "c" "c" "a" "a" "a" "t" "g"
     [73] "a" "a" "c" "a" "g" "a" "c" "a" "a" "g" "t" "a" "a" "a" "a" "g" "a" "c"
     [91] "a" "t" "g" "a" "c" "a" "g" "c" "g" "a" "t" "a" "c" "t" "t" "t" "c" "c"
    [109] "c" "a" "g" "a" "g" "c" "t" "g" "a" "a" "g" "t" "t" "a" "a" "c" "a" "a"
    [127] "a" "t" "g" "c" "a" "c" "c" "t" "g" "g" "t" "t" "c" "t" "t" "t" "t" "a"
    [145] "c" "t" "a" "a" "g" "t" "g" "t" "t" "c" "a" "a" "a" "t" "a" "c" "c" "a"
    [163] "g" "t" "g" "a" "a" "c" "t" "t" "a" "a" "a" "g" "a" "a" "t" "t" "t" "g"
    [181] "t" "c" "a" "a" "t" "c" "c" "t" "a" "g" "c" "c" "t" "t" "c" "c" "a" "a"
    [199] "g" "a" "g" "a" "a" "g" "a" "a" "a" "a" "a" "g" "a" "a" "g" "a" "g" "a"
    [217] "a" "a" "c" "t" "a" "g" "a" "a" "a" "c" "a" "g" "t" "t" "a" "a" "a" "g"
    [235] "t" "g" "t" "c" "t" "a" "a" "t" "a" "a" "t" "g" "c" "t" "g" "a" "a" "g"
    [253] "a" "c" "c" "c" "c" "a" "a" "a" "g" "a" "t" "c" "t" "c" "a" "t" "g" "t"
    [271] "t" "a" "a" "g" "t" "g" "g" "a" "g" "a" "a" "a" "g" "g" "g" "t" "t" "t"
    [289] "t" "g" "c" "a" "a" "a" "c" "t" "g" "a" "a" "a" "g" "a" "t" "c" "t" "g"
    [307] "t" "a" "g" "a" "g" "a" "g" "t" "a" "g" "c" "a" "g" "t" "a" "t" "t" "t"
    [325] "c" "a" "t" "t" "g" "g" "t" "a" "c" "c" "t" "g" "g" "t" "a" "c" "t" "g"
    [343] "a" "t" "t" "a" "t" "g" "g" "c" "a" "c" "t" "t" "a" "g"
    attr(,"name")
    [1] "KP255399.BRCA1"
    attr(,"Annot")
    [1] ">KP255399.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $KP255400.BRCA1
     [1] "g" "t" "c" "a" "t" "c" "c" "c" "c" "t" "t" "c" "t" "a" "a" "a" "t" "g" "c"
    [20] "c" "c" "a" "t" "c" "a" "t" "t" "a" "g" "a" "t" "g" "a" "t" "a" "g" "g" "t"
    [39] "g" "g" "t" "a" "c" "a" "t" "g" "c" "a" "c" "a" "g" "t" "t" "g" "c" "t" "c"
    [58] "t" "g" "g" "g" "a" "g" "t" "c" "t" "t" "c" "a" "g" "a" "a" "t" "a" "g" "a"
    [77] "a" "a" "c" "t" "a" "c" "c" "c" "a" "t" "c" "t" "t" "a" "a"
    attr(,"name")
    [1] "KP255400.BRCA1"
    attr(,"Annot")
    [1] ">KP255400.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $KP255401.BRCA1
      [1] "a" "g" "g" "g" "a" "a" "c" "c" "c" "c" "t" "t" "a" "c" "c" "t" "g" "g"
     [19] "a" "a" "t" "c" "t" "g" "g" "a" "a" "t" "c" "a" "g" "c" "c" "t" "c" "t"
     [37] "t" "c" "t" "c" "t" "g" "a" "t" "g" "a" "c" "c" "c" "t" "g" "a" "a" "t"
     [55] "c" "t" "g" "a" "t" "c" "c" "t" "t" "c" "t" "g" "a" "a" "g" "a" "c" "a"
     [73] "g" "a" "g" "c" "c" "c" "c" "a" "g" "a" "g" "t" "c" "a" "g" "c" "c" "g"
     [91] "t" "g" "t" "t" "g" "g" "c" "a" "a" "c" "a" "t" "a" "c" "c" "a" "t" "c"
    [109] "t" "t" "c" "a" "a" "c" "c" "t" "c" "t" "g" "c" "a" "t" "t" "g" "a"
    attr(,"name")
    [1] "KP255401.BRCA1"
    attr(,"Annot")
    [1] ">KP255401.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $KP255402.BRCA1
      [1] "a" "c" "a" "g" "c" "a" "g" "a" "a" "a" "c" "c" "t" "a" "c" "a" "a" "c"
     [19] "t" "c" "a" "t" "g" "g" "a" "a" "g" "g" "t" "a" "a" "a" "g" "a" "a" "c"
     [37] "c" "t" "g" "c" "a" "a" "c" "t" "g" "g" "a" "g" "c" "c" "a" "a" "g" "a"
     [55] "a" "g" "a" "g" "t" "a" "a" "c" "a" "a" "g" "c" "c" "a" "a" "a" "t" "g"
     [73] "a" "a" "c" "a" "g" "a" "c" "a" "a" "g" "t" "a" "a" "a" "a" "g" "a" "c"
     [91] "a" "t" "g" "a" "c" "a" "g" "c" "g" "a" "t" "a" "c" "t" "t" "t" "c" "c"
    [109] "c" "a" "g" "a" "g" "c" "t" "g" "a" "a" "g" "t" "t" "a" "a" "c" "a" "a"
    [127] "a" "t" "g" "c" "a" "c" "c" "t" "g" "g" "t" "t" "c" "t" "t" "t" "t" "a"
    [145] "c" "t" "a" "a" "g" "t" "g" "t" "t" "c" "a" "a" "a" "t" "a" "c" "c" "a"
    [163] "g" "t" "g" "a" "a" "c" "t" "t" "a" "a" "a" "g" "a" "a" "t" "t" "t" "g"
    [181] "t" "c" "a" "a" "t" "c" "c" "t" "a" "g" "c" "c" "t" "t" "c" "c" "a" "a"
    [199] "g" "a" "g" "a" "a" "g" "a" "a" "a" "a" "a" "g" "a" "a" "g" "a" "g" "a"
    [217] "a" "a" "c" "t" "a" "g" "a" "a" "a" "c" "a" "g" "t" "t" "a" "a" "a" "g"
    [235] "t" "g" "t" "c" "t" "a" "a" "t" "a" "a" "t" "g" "c" "t" "g" "a" "a" "g"
    [253] "a" "c" "c" "c" "c" "a" "a" "a" "g" "a" "t" "c" "t" "c" "a" "t" "t" "a"
    [271] "a" "g" "t" "g" "g" "a" "g" "a" "a" "a" "g" "g" "g" "t" "t" "t" "t" "g"
    [289] "c" "a" "a" "a" "c" "t" "g" "a"
    attr(,"name")
    [1] "KP255402.BRCA1"
    attr(,"Annot")
    [1] ">KP255402.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $KP255403.BRCA1
      [1] "t" "t" "g" "t" "t" "c" "t" "g" "a" "g" "a" "c" "a" "c" "c" "t" "g" "a"
     [19] "t" "g" "a" "c" "c" "t" "g" "t" "t" "a" "g" "a" "t" "g" "a" "t" "g" "g"
     [37] "t" "g" "a" "a" "a" "t" "a" "a" "a" "g" "g" "a" "a" "g" "a" "t" "a" "c"
     [55] "t" "a" "g" "t" "t" "t" "t" "g" "c" "t" "g" "a" "a" "a" "a" "t" "g" "a"
     [73] "c" "a" "t" "t" "a" "a" "g" "g" "a" "a" "a" "g" "t" "t" "c" "t" "g" "c"
     [91] "t" "g" "t" "t" "t" "t" "t" "a" "g" "c" "a" "a" "a" "a" "g" "c" "g" "t"
    [109] "c" "c" "a" "g" "a" "a" "a" "g" "g" "a" "g" "a" "g" "c" "t" "t" "a" "g"
    [127] "c" "a" "g" "g" "a" "g" "t" "c" "c" "t" "a" "g" "c" "c" "c" "t" "t" "t"
    [145] "c" "a" "c" "c" "c" "a" "t" "a" "c" "a" "c" "a" "t" "t" "t" "g" "g" "c"
    [163] "t" "c" "a" "g" "g" "g" "t" "t" "a" "c" "t" "g" "a"
    attr(,"name")
    [1] "KP255403.BRCA1"
    attr(,"Annot")
    [1] ">KP255403.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $KP272102.BRCA1
      [1] "t" "c" "t" "g" "g" "a" "g" "t" "t" "g" "a" "t" "c" "a" "a" "g" "g" "a"
     [19] "a" "c" "c" "t" "g" "t" "c" "t" "c" "c" "a" "c" "a" "a" "a" "g" "t" "g"
     [37] "t" "g" "a" "c" "c" "a" "c" "a" "t" "a" "t" "t" "t" "t" "g" "c" "a" "a"
     [55] "g" "a" "g" "c" "c" "t" "a" "c" "a" "a" "g" "a" "a" "a" "g" "t" "a" "c"
     [73] "g" "a" "g" "a" "t" "t" "t" "a" "g" "t" "c" "a" "a" "c" "t" "t" "g" "t"
     [91] "t" "g" "a" "a" "g" "a" "g" "c" "t" "a" "t" "t" "g" "a" "a" "a" "a" "t"
    [109] "c" "a" "t" "t" "t" "g" "t" "g" "c" "t" "t" "t" "t" "c" "a" "g" "c" "t"
    [127] "t" "g" "a" "c" "a" "c" "a" "g" "g" "t" "t" "t" "g" "g" "a" "g" "t"
    attr(,"name")
    [1] "KP272102.BRCA1"
    attr(,"Annot")
    [1] ">KP272102.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $KP272103.BRCA1
     [1] "a" "t" "t" "t" "t" "g" "c" "a" "t" "g" "c" "t" "g" "a" "a" "a" "c" "t" "t"
    [20] "c" "t" "c" "a" "a" "c" "c" "a" "g" "a" "a" "g" "a" "a" "a" "g" "g" "g" "c"
    [39] "c" "t" "t" "c" "a" "c" "a" "g" "t" "g" "t" "c" "c" "t" "t" "t" "a" "t" "g"
    [58] "t" "a" "a" "g" "a" "a" "t" "g" "a" "t" "a" "t" "a" "a" "c" "c" "a" "a" "a"
    [77] "a" "g" "g" "a" "g" "c" "c" "t" "a" "t" "a" "a"
    attr(,"name")
    [1] "KP272103.BRCA1"
    attr(,"Annot")
    [1] ">KP272103.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $KP272104.BRCA1
      [1] "c" "a" "g" "t" "a" "t" "t" "a" "a" "c" "t" "t" "c" "a" "c" "a" "g" "a"
     [19] "a" "a" "a" "g" "t" "a" "g" "t" "g" "a" "a" "t" "a" "c" "c" "c" "t" "a"
     [37] "t" "a" "a" "g" "c" "c" "a" "g" "a" "a" "t" "c" "c" "a" "g" "a" "a" "g"
     [55] "g" "c" "c" "t" "t" "t" "c" "t" "g" "c" "t" "g" "a" "c" "a" "a" "g" "t"
     [73] "t" "t" "g" "a" "g" "g" "t" "g" "t" "c" "t" "g" "c" "a" "g" "a" "t" "a"
     [91] "g" "t" "t" "c" "t" "a" "c" "c" "a" "g" "t" "a" "a" "a" "a" "a" "t" "a"
    [109] "a" "a" "g" "a" "a" "c" "c" "a" "g" "g" "a" "g" "t" "g" "g" "a" "a" "a"
    [127] "g" "a" "g" "g" "g" "a" "a" "c" "c" "c" "c" "t" "t" "a" "c" "c" "t" "g"
    [145] "g" "a" "a" "t" "c" "t" "g" "g" "a" "a" "t" "c" "a" "g" "c" "c" "t" "c"
    [163] "t" "t" "c" "t" "c" "t" "g" "a"
    attr(,"name")
    [1] "KP272104.BRCA1"
    attr(,"Annot")
    [1] ">KP272104.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $KP272105.BRCA1
      [1] "t" "c" "t" "g" "g" "a" "g" "t" "t" "g" "a" "t" "c" "a" "a" "g" "g" "a"
     [19] "a" "c" "c" "t" "g" "t" "c" "t" "c" "c" "a" "c" "a" "a" "a" "g" "t" "g"
     [37] "t" "g" "a" "c" "c" "a" "c" "a" "t" "a" "t" "t" "t" "t" "g" "c" "a" "a"
     [55] "a" "t" "t" "t" "t" "g" "c" "a" "t" "g" "c" "t" "g" "a" "a" "a" "c" "t"
     [73] "t" "c" "t" "c" "a" "a" "c" "c" "a" "g" "a" "a" "g" "a" "a" "a" "g" "g"
     [91] "g" "c" "c" "t" "t" "c" "a" "c" "a" "g" "t" "g" "t" "c" "c" "t" "t" "t"
    [109] "a" "t" "g" "t" "a" "a" "g" "a" "a" "t" "g" "a" "t" "a" "t" "a" "a" "c"
    [127] "c" "a" "a" "a" "a" "g" "g" "a" "g" "c" "c" "t" "a" "c" "a" "a" "g" "a"
    [145] "a" "a" "g" "t" "a" "c" "g" "a" "g" "a" "t" "t" "t" "a" "g" "t" "c" "a"
    [163] "a" "c" "t" "t" "g" "t" "t" "g" "a" "a" "g" "a" "g" "c" "t" "a" "t" "t"
    [181] "g" "a" "a" "a" "a" "t" "c" "a" "t" "t" "t" "g" "t" "g" "c" "t" "t" "t"
    [199] "t" "c" "a" "g" "c" "t" "t" "g" "a" "c" "a" "c" "a" "g" "g" "t" "t" "t"
    [217] "g" "g" "a" "g" "t" "a" "t" "g" "c" "a" "a" "a" "c" "a" "g" "c" "t" "a"
    [235] "t" "a" "a" "t" "t" "t" "t" "g" "c" "a" "a" "a" "a" "a" "a" "g" "g" "a"
    [253] "a" "a" "a" "t" "a" "a" "c" "t" "c" "t" "c" "c" "t" "g" "a" "a" "c" "a"
    [271] "t" "c" "t" "a" "a" "a" "a" "g" "a" "t" "g" "a" "a" "g" "t" "t" "t" "c"
    [289] "t" "a" "t" "c" "a" "t" "c" "c" "a" "a" "a" "g" "t" "a" "t" "g" "g" "g"
    [307] "c" "t" "a" "c" "a" "g" "a" "a" "a" "c" "c" "g" "t" "g" "c" "c" "a" "a"
    [325] "a" "a" "g" "a" "c" "t" "t" "c" "t" "a" "c" "a" "g" "a" "g" "t" "g" "a"
    [343] "a" "c" "c" "c" "g" "a" "a" "a" "a" "t" "c" "c" "t" "t" "c" "c" "t" "t"
    [361] "g" "a" "t" "t" "t" "t" "g" "c" "a" "t" "g" "c" "t" "g" "a" "a" "a" "c"
    [379] "t" "t" "c" "t" "c" "a" "a" "c" "c" "a" "g" "a" "a" "g" "a" "a" "a" "g"
    [397] "g" "g" "c" "c" "t" "t" "c" "a" "c" "a" "g" "t" "g" "t" "c" "c" "t" "t"
    [415] "t" "a" "t" "g" "t" "a" "a"
    attr(,"name")
    [1] "KP272105.BRCA1"
    attr(,"Annot")
    [1] ">KP272105.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $KP272106.BRCA1
      [1] "g" "g" "g" "t" "g" "a" "c" "c" "c" "a" "g" "t" "c" "t" "a" "t" "t" "a"
     [19] "a" "a" "g" "a" "a" "a" "g" "a" "a" "a" "a" "a" "t" "g" "c" "t" "g" "a"
     [37] "a" "t" "g" "a" "g" "g" "g" "t" "g" "t" "c" "c" "a" "c" "c" "c" "a" "a"
     [55] "t" "t" "g" "t" "g" "g" "t" "t" "g" "t" "g" "c" "a" "g" "c" "c" "a" "g"
     [73] "a" "t" "g" "c" "c" "t" "g" "g" "a" "c" "a" "g" "a" "g" "g" "a" "c" "a"
     [91] "a" "t" "g" "g" "c" "t" "t" "c" "c" "a" "t" "g" "c" "a" "a" "t" "t" "g"
    [109] "g" "g" "c" "a" "g" "a" "t" "g" "t" "g" "t" "g" "a" "g" "g" "c" "a" "c"
    [127] "c" "t" "g" "t" "g" "g" "t" "g" "a" "c" "c" "c" "g" "a" "g" "a" "g" "t"
    [145] "g" "g" "g" "t" "g" "t" "t" "g" "g" "a" "c" "a" "g" "t" "g" "t" "a" "g"
    [163] "c" "a" "c" "t" "c" "t" "a" "c" "c" "a" "g" "t" "g" "c" "c" "a" "g" "g"
    [181] "a" "g" "c" "t" "g" "g" "a" "c" "a" "c" "c" "t" "a" "c" "c" "t" "g" "a"
    [199] "t" "a" "c" "c" "c" "c" "a" "g" "a" "t" "c" "c" "c" "c" "c" "a" "c" "a"
    [217] "g" "c" "c" "a" "c" "t" "a" "c" "t" "g" "a"
    attr(,"name")
    [1] "KP272106.BRCA1"
    attr(,"Annot")
    [1] ">KP272106.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $KP404097
      [1] "g" "a" "a" "c" "a" "g" "a" "c" "a" "a" "g" "t" "a" "a" "a" "a" "g" "a"
     [19] "c" "a" "t" "g" "a" "c" "a" "g" "t" "g" "a" "t" "a" "c" "t" "t" "t" "c"
     [37] "c" "c" "a" "g" "a" "g" "c" "t" "g" "a" "a" "g" "t" "t" "a" "a" "c" "a"
     [55] "a" "a" "t" "g" "c" "a" "c" "c" "t" "g" "g" "t" "t" "c" "t" "t" "t" "t"
     [73] "a" "c" "t" "a" "a" "g" "t" "g" "t" "t" "c" "a" "a" "a" "t" "a" "c" "c"
     [91] "a" "g" "t" "g" "a" "a" "c" "t" "t" "a" "a" "a" "g" "a" "a" "t" "t" "t"
    [109] "g" "t" "c" "a" "a" "t" "c" "c" "t" "a" "g" "c" "c" "t" "t" "c" "c" "a"
    [127] "a" "g" "a" "g" "a" "a" "g" "a" "a" "a" "a" "a" "g" "a" "a" "g" "a" "g"
    [145] "a" "a" "a" "c" "t" "a" "g" "a" "a" "a" "c" "a" "g" "t" "t" "a" "a" "a"
    [163] "g" "t" "g" "t" "c" "t" "a" "a" "t" "a" "a" "t" "g" "c" "t" "g" "a" "a"
    [181] "g" "a" "c" "c" "c" "c" "a" "a" "a" "g" "a" "t" "c" "t" "c" "a" "t" "g"
    [199] "t" "t" "a" "a" "g" "t" "g" "g" "a" "g" "a" "a" "a" "g" "g" "g" "t" "t"
    [217] "t" "t" "g" "c" "a" "a" "a" "c" "t" "g" "a" "a" "a" "g" "a" "t" "c" "t"
    [235] "g" "t" "a" "g" "a" "g" "a" "g" "t" "a" "g" "c" "a" "g" "t" "a" "t" "t"
    [253] "t" "c" "a" "t" "t" "g" "g" "t" "a" "c" "c" "t" "g" "g" "t" "a" "c" "t"
    [271] "g" "a" "t" "t" "a" "t" "g" "g" "c" "a" "c" "t" "c" "a" "g" "g" "a" "a"
    [289] "a" "g" "t" "a" "t" "c" "t" "c" "g" "t" "t" "a" "c" "t" "g" "g" "a" "a"
    [307] "g" "t" "t" "a" "g" "c" "a" "c" "t" "c" "t" "a" "g" "g" "g" "a" "a" "g"
    [325] "g" "c" "a" "a" "a" "a" "a" "c" "a" "g" "a" "a" "c" "c" "a" "a" "a" "t"
    [343] "a" "a" "a" "t" "g" "t" "g" "t" "g" "a" "g" "t" "c" "a" "g" "t" "g" "t"
    [361] "g" "c" "a" "g" "c" "a" "t" "t" "t" "g" "a" "a" "a" "a" "c" "c" "c" "c"
    [379] "a" "a" "g" "g" "g" "a" "c" "t" "a" "a" "t" "t" "c" "a" "t" "g" "g" "t"
    [397] "t" "g" "t" "t" "c" "c" "a" "a" "a" "g" "a" "t" "a" "a" "t" "a" "g" "a"
    [415] "a" "a" "t" "g" "a" "c" "a" "c" "a" "g" "a" "a" "g" "g" "c" "t" "t" "t"
    [433] "a" "a" "g" "t" "a" "t" "c" "c" "a" "t" "t" "g" "g" "g" "a" "c" "a" "t"
    [451] "g" "a" "a" "g" "t" "t" "a" "a" "c" "c" "a" "c" "a" "g" "t" "c" "g" "g"
    [469] "g" "a" "a" "a" "c" "a" "a" "g" "c" "a" "t" "a" "g" "a" "a" "a" "t" "g"
    [487] "g" "a" "a" "g" "a" "a" "a" "g" "t" "g" "a" "a" "c" "t" "t" "g" "a" "t"
    [505] "g" "c" "t" "c" "a" "g" "t" "a" "t" "t" "t" "g" "c" "a" "g" "a" "a" "t"
    [523] "a" "c" "a" "t" "t" "c" "a" "a" "g" "g" "t" "t" "t" "c" "a" "a" "a" "g"
    [541] "c" "g" "c" "c" "a" "g" "t" "c" "a" "t" "t" "t" "g" "c" "t" "c" "c" "g"
    [559] "t" "t" "t" "t" "c" "a" "a" "a" "t" "c" "c" "a" "g" "g" "a" "a" "a" "t"
    [577] "g" "c" "a" "g" "a" "a" "g" "a" "g" "g" "a" "a" "t" "g" "t" "g" "c" "a"
    [595] "a" "c" "a" "t" "t" "c" "t" "c" "t" "g" "c" "c" "c" "a" "c" "t" "c" "t"
    [613] "g" "g" "g" "t" "c" "c" "t" "t" "a" "a" "a" "g" "a" "a" "a" "c" "a" "a"
    attr(,"name")
    [1] "KP404097"
    attr(,"Annot")
    [1] ">KP404097"
    attr(,"class")
    [1] "SeqFastadna"
    
    $KP455327.BRCA1
      [1] "c" "a" "g" "a" "g" "g" "g" "a" "t" "a" "c" "c" "a" "t" "g" "c" "a" "a"
     [19] "c" "a" "t" "a" "a" "c" "c" "t" "g" "a" "t" "a" "a" "a" "g" "c" "t" "c"
     [37] "c" "a" "g" "c" "a" "g" "g" "a" "a" "a" "t" "g" "g" "c" "t" "g" "a" "a"
     [55] "c" "t" "a" "g" "a" "a" "g" "c" "t" "g" "t" "g" "t" "t" "a" "g" "a" "a"
     [73] "c" "a" "g" "c" "a" "t" "g" "g" "g" "a" "g" "c" "c" "a" "g" "c" "c" "t"
     [91] "t" "c" "t" "a" "a" "c" "a" "g" "c" "t" "a" "c" "c" "c" "t" "t" "c" "c"
    [109] "a" "t" "c" "a" "t" "a" "a" "g" "t" "g" "a" "c" "t" "c" "c" "t" "c" "t"
    [127] "g" "c" "c" "c" "t" "t" "g" "a" "g" "g" "a" "c" "c" "t" "g" "c" "g" "a"
    [145] "a" "a" "t" "c" "c" "a" "g" "a" "a" "c" "a" "a" "a" "g" "c" "a" "c" "a"
    [163] "t" "c" "a" "g" "a" "a" "a" "a" "a" "g"
    attr(,"name")
    [1] "KP455327.BRCA1"
    attr(,"Annot")
    [1] ">KP455327.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $KP701015
      [1] "a" "a" "a" "a" "c" "a" "g" "a" "a" "c" "c" "a" "a" "a" "t" "a" "a" "a"
     [19] "t" "g" "t" "g" "t" "g" "a" "g" "t" "c" "a" "g" "t" "g" "t" "g" "c" "a"
     [37] "g" "c" "a" "t" "t" "t" "g" "a" "a" "a" "a" "c" "c" "c" "c" "a" "a" "g"
     [55] "g" "g" "a" "c" "t" "a" "a" "t" "t" "c" "a" "t" "g" "g" "t" "t" "g" "t"
     [73] "t" "c" "c" "a" "a" "a" "g" "a" "t" "a" "a" "t" "a" "g" "a" "a" "a" "t"
     [91] "g" "a" "c" "a" "c" "a" "g" "a" "a" "g" "g" "c" "t" "t" "t" "a" "a" "g"
    [109] "t" "a" "t" "c" "c" "a" "t" "t" "g" "g" "g" "a" "c" "a" "t" "g" "a" "a"
    [127] "g" "t" "t" "a" "a" "c" "c" "a" "c" "a" "g" "t" "c" "g" "g" "g" "a" "a"
    [145] "a" "c" "a" "a" "g" "c" "a" "t" "a" "g" "a" "a" "a" "t" "g" "g" "a" "a"
    [163] "g" "a" "a" "a" "g" "t" "g" "a" "a" "c" "t" "t" "g" "a" "t" "g" "c" "t"
    [181] "c" "a" "g" "t" "a" "t" "t" "t" "g" "c" "a" "g" "a" "a" "t" "a" "c" "a"
    [199] "t" "t" "c" "a" "a" "g" "g" "t" "t" "t" "c" "a" "a" "a" "g" "c" "g" "c"
    [217] "c" "a" "g" "t" "c" "a" "t" "t" "t" "g" "c" "t" "c" "c" "g" "t" "t" "t"
    [235] "t" "c" "a" "a" "a" "t" "c" "c" "a" "g" "g" "a" "a" "a" "t" "g" "c" "a"
    [253] "g" "a" "a" "g" "a" "g" "g" "a" "a" "t" "g" "t" "g" "c" "a" "g" "g" "c"
    [271] "a" "a" "c" "a" "t" "t" "c" "t" "c" "t" "g" "c" "c" "c" "a" "c" "t" "c"
    [289] "t" "g" "g" "g" "t" "c" "c" "t" "t" "a" "a" "a" "g" "a" "a" "a" "c" "a"
    [307] "a" "a" "g" "t" "c" "c" "a" "a" "a" "a" "g" "t" "c" "a" "c" "t" "t" "t"
    [325] "t" "g" "a"
    attr(,"name")
    [1] "KP701015"
    attr(,"Annot")
    [1] ">KP701015"
    attr(,"class")
    [1] "SeqFastadna"
    
    $KP701016.BRCA1
      [1] "t" "c" "a" "t" "c" "c" "c" "c" "t" "t" "c" "t" "a" "a" "a" "t" "g" "c"
     [19] "c" "c" "a" "t" "c" "a" "t" "t" "a" "g" "a" "t" "g" "a" "t" "a" "g" "g"
     [37] "t" "g" "g" "t" "a" "c" "a" "t" "g" "c" "a" "c" "a" "g" "t" "t" "g" "c"
     [55] "t" "c" "t" "g" "g" "g" "a" "g" "t" "c" "t" "t" "c" "a" "g" "a" "a" "t"
     [73] "a" "g" "a" "a" "a" "c" "t" "a" "c" "c" "c" "a" "t" "c" "t" "c" "a" "a"
     [91] "g" "a" "g" "g" "a" "g" "c" "t" "c" "a" "t" "t" "a" "a" "g" "g" "t" "t"
    [109] "g" "t" "t" "g" "a" "t" "g" "t" "g" "g" "a" "g" "g" "a" "g" "t" "a" "a"
    attr(,"name")
    [1] "KP701016.BRCA1"
    attr(,"Annot")
    [1] ">KP701016.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $KP729136
      [1] "g" "a" "g" "a" "a" "c" "a" "t" "t" "c" "c" "a" "a" "g" "t" "a" "c" "a"
     [19] "g" "t" "g" "a" "g" "c" "a" "c" "a" "a" "t" "t" "a" "g" "c" "c" "g" "t"
     [37] "a" "a" "t" "a" "a" "c" "a" "t" "t" "a" "g" "a" "g" "a" "a" "a" "a" "t"
     [55] "g" "t" "t" "t" "t" "t" "a" "a" "a" "g" "a" "a" "g" "c" "c" "a" "g" "c"
     [73] "t" "c" "a" "a" "g" "c" "a" "a" "t" "a" "t" "t" "a" "a" "t" "g" "a" "a"
     [91] "g" "t" "a" "g" "g" "t" "t" "c" "c" "a" "g" "t" "a" "c" "t" "a" "a" "t"
    [109] "g" "a" "a" "g" "t" "g" "g" "g" "c" "t" "c" "c" "a" "g" "t" "a" "t" "t"
    [127] "a" "a" "t" "g" "a" "a" "a" "t" "a" "g" "g" "t" "t" "c" "c" "a" "g" "t"
    [145] "g" "a" "t" "g" "a" "a" "a" "a" "c" "a" "t" "t" "c" "a" "a" "g" "c" "a"
    [163] "g" "a" "a" "c" "t" "a" "g" "g" "t" "a" "g" "a" "a" "a" "c" "a" "g" "a"
    [181] "g" "g" "g" "c" "c" "a" "a" "a" "a" "t" "t" "g" "a" "a" "t" "g" "c" "t"
    [199] "a" "t" "g" "c" "t" "t" "a" "g" "a" "t" "t" "a" "g" "g" "g" "g" "t" "t"
    [217] "t" "t" "g" "c" "a" "a" "c" "c" "t" "g" "a" "g" "g" "t" "c" "t" "a" "t"
    [235] "a" "a" "a" "c" "a" "a" "a" "g" "t" "c" "t" "t" "c" "c" "t" "g" "g" "a"
    [253] "a" "g" "t" "a" "a" "t" "t" "g" "t" "a" "a" "g" "c" "a" "t" "c" "c" "t"
    [271] "g" "a" "a" "a" "t" "a" "a" "a" "a" "a" "a" "g" "c" "a" "a" "g" "a" "a"
    [289] "t" "a" "t" "g" "a" "a" "g" "a" "a" "g" "t" "a" "g" "t" "t" "c" "a" "g"
    [307] "a" "c" "t" "g" "t" "t" "a" "a" "t" "a" "c" "a" "g" "a" "t" "t" "t" "c"
    [325] "t" "c" "t" "c" "c" "a" "t" "a" "t" "c" "t" "g" "a" "t" "t" "t" "c" "a"
    [343] "g" "a" "t" "a" "a" "c" "t" "t" "a" "g" "a" "a" "c" "a" "g" "c" "c" "t"
    [361] "a" "t" "g" "g" "g" "a" "a" "g" "t" "a" "g" "t" "c" "a" "t" "g" "c" "a"
    [379] "t" "c" "t" "c" "a" "g" "g" "t" "t" "t" "g" "t" "t" "c" "t" "g" "a" "g"
    [397] "a" "c" "a" "c" "c" "t" "g" "a" "t" "g" "a" "c" "c" "t" "g" "t" "t" "a"
    [415] "g" "a" "t" "g" "a" "t" "g" "g" "t" "g" "a" "a" "a" "t" "a" "a" "a" "g"
    [433] "g" "a" "a" "g" "a" "t" "a" "c" "t" "a" "g" "t" "t" "t" "t" "g" "c" "t"
    [451] "g" "a" "a" "a" "a" "t" "g" "a" "c" "a" "t" "t" "a" "a" "g" "g" "a" "a"
    [469] "a" "g" "t" "t" "c" "t" "g" "c" "t" "g" "t" "t" "t" "t" "t" "a" "g" "c"
    [487] "a" "a" "a" "a" "g" "c" "g" "t" "c" "c" "a" "g" "a" "g" "a" "g" "g" "a"
    [505] "g" "a" "g" "c" "t" "t" "a" "g" "c" "a" "g" "g" "a" "g" "t" "c" "c" "t"
    [523] "a" "g" "c" "c" "c" "t" "t" "t" "c" "a" "c" "c" "c" "a" "t" "a" "c" "a"
    [541] "c" "a" "t" "t" "t" "g" "g" "c" "t" "c" "a" "g" "g" "g" "t" "t" "a" "c"
    [559] "c" "g" "a" "a" "g" "a" "g" "g" "g" "g" "c" "c" "a" "a" "g" "a" "a" "a"
    [577] "t" "t" "a" "g" "a" "g" "t" "c" "c" "t" "c" "a" "g" "a" "a" "g" "a" "g"
    [595] "a" "a" "c" "t" "t" "a" "t" "c" "t" "a" "g" "t" "g" "a" "g" "g" "a" "t"
    [613] "g" "a" "a" "g" "a" "g" "c" "t" "t" "c" "c" "c" "t" "g" "c" "t" "t" "c"
    [631] "c" "a" "a" "c" "a" "c" "t" "t" "g" "t" "t" "a" "t" "t" "t" "g" "g" "t"
    [649] "a" "a" "a" "g" "t" "a" "a" "a" "c" "a" "a" "t" "a" "t" "a" "c" "c" "t"
    [667] "t" "c" "t" "c" "a" "g" "t" "c" "t" "a" "c" "t" "a" "g" "g" "c" "a" "t"
    [685] "a" "g" "c" "a" "c" "c" "g" "t" "t" "g" "c" "t" "a" "c" "c" "g" "a" "g"
    [703] "t" "g" "t" "c" "t" "g" "t" "c" "t" "a" "a" "g" "a" "a" "c"
    attr(,"name")
    [1] "KP729136"
    attr(,"Annot")
    [1] ">KP729136"
    attr(,"class")
    [1] "SeqFastadna"
    
    $KP729137
      [1] "a" "g" "g" "g" "a" "a" "c" "c" "c" "c" "t" "t" "a" "c" "c" "t" "g" "g"
     [19] "a" "a" "t" "c" "t" "g" "g" "a" "a" "t" "c" "a" "g" "c" "c" "t" "c" "t"
     [37] "t" "c" "t" "c" "t" "g" "a" "t" "g" "a" "c" "c" "c" "t" "g" "a" "a" "t"
     [55] "c" "t" "g" "a" "t" "c" "c" "t" "t" "c" "t" "g" "a" "a" "g" "a" "c" "a"
     [73] "g" "a" "g" "c" "c" "c" "c" "a" "g" "a" "g" "t" "c" "a" "g" "c" "t" "c"
     [91] "g" "t" "g" "t" "t" "g" "g" "c" "a" "a" "c" "a" "t" "a" "c" "c" "a" "t"
    [109] "c" "t" "t" "c" "a" "a" "c" "c" "t" "c" "t" "g" "c" "a" "t" "t" "g" "a"
    [127] "a" "a" "g" "t" "t" "c" "c" "c" "c" "a" "a" "t" "t" "g" "a" "a" "a" "g"
    [145] "t" "t" "g" "c" "a" "g" "a" "a" "t" "c" "t" "g" "c" "c" "c" "a" "g" "g"
    [163] "g" "t" "c" "c" "a" "g" "c" "t" "g" "c" "t" "g" "c" "t" "c" "a" "t" "a"
    [181] "c" "t" "a" "c" "t" "g" "a" "t" "a" "c" "t" "g" "c" "t" "g" "g" "g" "t"
    [199] "a" "t" "a" "a" "t" "g" "c" "a" "a" "t" "g" "g" "a" "a" "g" "a" "a" "a"
    [217] "g" "t" "g" "t" "g" "a" "g" "c" "a" "g" "g" "g" "a" "g" "a" "a" "g" "c"
    [235] "c" "a" "g" "a" "a" "t" "t" "g" "a" "c" "a" "g" "c" "t" "t" "c" "a" "a"
    [253] "c" "a" "g" "a" "a" "a" "g" "g" "g" "t" "c" "a" "a" "c" "a" "a" "a" "a"
    [271] "g" "a" "a" "t" "g" "t" "c" "c" "a" "t" "g" "g" "t" "g" "g" "t" "g" "t"
    [289] "c" "t" "g" "g" "c" "c" "t" "g" "a" "c" "c" "c" "c" "a" "g" "a" "a" "g"
    [307] "a" "a" "t" "t" "t"
    attr(,"name")
    [1] "KP729137"
    attr(,"Annot")
    [1] ">KP729137"
    attr(,"class")
    [1] "SeqFastadna"
    
    $KP744861
      [1] "g" "g" "t" "t" "g" "t" "t" "c" "c" "a" "a" "a" "g" "a" "t" "a" "a" "t"
     [19] "a" "g" "a" "a" "a" "t" "g" "a" "c" "a" "c" "a" "g" "a" "a" "g" "g" "c"
     [37] "t" "t" "t" "a" "a" "g" "t" "a" "t" "c" "c" "a" "t" "t" "g" "g" "g" "a"
     [55] "c" "a" "t" "g" "a" "a" "g" "t" "t" "a" "a" "c" "c" "a" "c" "a" "g" "t"
     [73] "c" "g" "g" "g" "a" "a" "a" "c" "a" "a" "g" "c" "a" "t" "a" "g" "a" "a"
     [91] "a" "t" "g" "g" "a" "a" "g" "a" "a" "a" "g" "t" "g" "a" "a" "c" "t" "t"
    [109] "g" "a" "t" "g" "c" "t" "c" "a" "g" "t" "a" "t" "t" "t" "g" "c" "a" "g"
    [127] "a" "a" "t" "a" "c" "a" "t" "t" "c" "a" "a" "g" "g" "t" "t" "t" "c" "a"
    [145] "a" "a" "g" "c" "g" "c" "c" "a" "g" "t" "c" "a" "t" "t" "t" "g" "c" "t"
    [163] "c" "c" "g" "t" "t" "t" "t" "c" "a" "a" "a" "t" "c" "c" "a" "g" "g" "a"
    [181] "a" "a" "t" "g" "c" "a" "g" "a" "a" "g" "a" "g" "g" "a" "a" "t" "g" "t"
    [199] "g" "c" "a" "a" "c" "a" "t" "t" "c" "t" "c" "t" "g" "c" "c" "c" "a" "c"
    [217] "t" "c" "t" "g" "g" "g" "t" "c" "c" "t" "t" "a" "a" "a" "g" "a" "a" "a"
    [235] "c" "a" "a" "a" "g" "t" "c" "c" "a" "a" "a" "a" "g" "t" "c" "a" "c" "t"
    [253] "t" "t" "t" "g" "a" "a" "t" "g" "t" "g" "a" "a" "c" "a" "a" "a" "a" "g"
    [271] "g" "a" "a" "g" "a" "a" "a" "a" "t" "c" "a" "a" "g" "g" "a" "a" "a" "g"
    [289] "a" "a" "t" "g" "a" "g" "t" "c" "t" "a" "a" "t" "a" "t" "c" "a" "a" "g"
    [307] "c" "c" "t" "g" "t" "a" "c" "a" "g" "a" "c" "a" "g" "t" "t" "a" "a" "t"
    [325] "a" "t" "c" "a" "c" "t" "g" "c" "a" "g" "g" "c" "t" "t" "t" "c" "c" "t"
    [343] "g" "t" "g" "g" "t" "t" "g" "g" "t" "c" "a" "g" "a" "a" "a" "g" "a" "t"
    [361] "a" "a" "g" "c" "c" "a" "g" "t" "t" "g" "a" "t" "a" "a" "t" "g" "c" "c"
    [379] "a" "a" "a" "t" "g" "t" "a" "g" "t" "a" "t" "c" "a" "a" "a" "g" "g" "a"
    [397] "g" "g" "c" "t" "c" "t" "a" "g" "g" "t" "t" "t" "t" "g" "t" "c" "t" "a"
    [415] "t" "c" "a" "t" "c" "t" "c" "a" "g" "t" "t" "c" "a" "g" "a" "g" "g" "c"
    [433] "a" "a" "c" "g" "a" "a" "a" "c" "t" "g" "g" "a" "c" "t" "c" "a" "t" "t"
    [451] "a" "c" "t" "c" "c" "a" "a" "a" "t" "a" "a" "a" "c" "a" "t" "g" "g" "a"
    [469] "c" "t" "t" "t" "t" "a" "c" "a" "a" "a" "a" "c" "c" "c" "a" "t" "a" "t"
    [487] "c" "g" "t" "a" "t" "a" "c" "c" "a" "c" "c" "a" "c" "t" "t" "t" "t" "t"
    [505] "c" "c" "c" "a" "t" "c" "a" "a" "g" "t" "c" "a" "t" "t" "t" "g" "t" "t"
    [523] "a" "a" "a" "a" "c" "t" "a" "a" "a" "t" "g" "t" "a" "a" "g" "a" "a" "a"
    [541] "a" "a" "t" "c" "t" "g" "c" "t" "a" "g" "a" "g" "g" "a" "a" "a" "a" "c"
    [559] "t" "t" "t" "g" "a" "g" "g" "a" "a" "c" "a" "t" "t" "c" "a" "a" "t" "g"
    [577] "t" "c" "a" "c" "c" "t" "g" "a" "a" "a" "g" "a" "g" "a" "a" "a" "t" "g"
    [595] "g" "g" "a" "a" "a" "t" "g" "a" "g" "a" "a" "c" "a" "t" "t" "c" "c" "a"
    [613] "a" "g" "t" "a" "c" "a" "g" "t" "g" "a" "g" "c" "a" "c" "a" "a" "t" "t"
    [631] "a" "g" "c" "c" "g" "t" "a" "a" "t" "a" "a" "c" "a" "t" "t" "a" "g" "a"
    [649] "g" "a" "a" "a" "a" "t" "g" "t" "t" "t" "t" "t" "a" "a" "a" "g" "g" "a"
    [667] "g" "c" "c" "a" "g" "c" "t" "c" "a" "a" "g" "c" "a" "a" "t" "a" "t" "t"
    [685] "a" "a" "t" "g" "a" "a" "g" "t" "a" "g" "g" "t" "t" "c" "c" "a" "g" "t"
    [703] "a" "c" "t" "a" "a" "t" "g" "a" "a" "g" "t" "g" "g" "g" "c" "t" "c" "c"
    [721] "a" "g" "t" "a" "t" "t" "a" "a" "t" "g" "a" "a" "a" "t" "a" "g" "g" "t"
    [739] "t" "c" "c" "a" "g" "t" "g" "a" "t" "g" "a" "a" "a" "a" "c" "a" "t" "t"
    [757] "c" "a" "a" "g" "c" "a" "g" "a" "a" "c" "t" "a" "g" "g" "t" "a" "g" "a"
    [775] "a" "a" "c"
    attr(,"name")
    [1] "KP744861"
    attr(,"Annot")
    [1] ">KP744861"
    attr(,"class")
    [1] "SeqFastadna"
    
    $KP753383.BRCA1
     [1] "g" "t" "g" "a" "a" "g" "c" "a" "g" "c" "a" "t" "c" "t" "g" "g" "g" "t" "g"
    [20] "t" "g" "a" "g" "a" "g" "t" "g" "a" "a" "a" "c" "a" "a" "g" "c" "g" "t" "c"
    [39] "t" "c" "t" "g" "a" "a" "g" "a" "c" "t" "g" "c" "t" "c" "a" "g" "g" "g" "c"
    [58] "t" "a" "t" "c" "c" "t" "c" "t" "c" "a" "g" "a" "g" "t" "g" "a" "c" "a" "t"
    [77] "t" "t" "t" "a" "a" "c" "c" "a" "c" "t" "c" "a" "g"
    attr(,"name")
    [1] "KP753383.BRCA1"
    attr(,"Annot")
    [1] ">KP753383.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $KT120061.BRCA1
     [1] "a" "t" "g" "g" "a" "t" "t" "t" "a" "t" "c" "t" "g" "c" "t" "c" "t" "t" "c"
    [20] "g" "c" "g" "t" "t" "g" "a" "a" "g" "a" "a" "g" "t" "a" "c" "a" "a" "a" "a"
    [39] "t" "g" "t" "c" "a" "t" "t" "a" "a" "t" "g" "c" "t" "a" "t" "g" "c" "a" "g"
    [58] "a" "a" "a" "a" "t" "c" "t" "t" "a" "g" "a" "g" "t" "g" "t" "c" "c" "c" "a"
    [77] "t" "c" "t" "c"
    attr(,"name")
    [1] "KT120061.BRCA1"
    attr(,"Annot")
    [1] ">KT120061.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $KT152888
      [1] "a" "a" "t" "t" "g" "t" "a" "a" "g" "c" "a" "t" "c" "c" "t" "g" "a" "a"
     [19] "a" "t" "a" "a" "a" "a" "a" "a" "g" "c" "a" "a" "g" "a" "a" "t" "a" "t"
     [37] "g" "a" "a" "g" "a" "a" "g" "t" "a" "g" "t" "t" "c" "a" "g" "a" "c" "t"
     [55] "g" "t" "t" "a" "a" "t" "a" "c" "a" "g" "a" "t" "t" "t" "c" "t" "c" "t"
     [73] "c" "c" "a" "t" "a" "t" "c" "t" "g" "a" "t" "t" "t" "c" "a" "g" "a" "t"
     [91] "a" "a" "c" "t" "t" "a" "g" "a" "a" "c" "a" "g" "c" "c" "t" "a" "t" "g"
    [109] "g" "g" "a" "a" "g" "t" "a" "g" "t" "c" "a" "t" "g" "c" "a" "t" "c" "t"
    [127] "c" "a" "g" "g" "t" "t" "t" "g" "t" "t" "c" "t" "g" "a" "g" "a" "c" "a"
    [145] "c" "c" "t" "g" "a" "t" "g" "a" "c" "c" "t" "g" "t" "t" "a" "g" "a" "t"
    [163] "g" "a" "t" "g" "g" "t" "g" "a" "a" "a" "t" "a" "a" "a" "g" "g" "a" "a"
    [181] "g" "a" "t" "a" "c" "t" "a" "g" "t" "t" "t" "t" "g" "c" "t" "g" "a" "a"
    [199] "a" "a" "t" "g" "a" "c" "a" "t" "t" "a" "a" "g" "g" "a" "a" "a" "g" "t"
    [217] "t" "c" "t" "g" "c" "t" "g" "t" "t" "t" "t" "t" "a" "g" "c" "a" "a" "a"
    [235] "a" "g" "c" "g" "t" "c" "c" "g" "a" "a" "a" "g" "g" "a" "g" "a" "g" "c"
    [253] "t" "t" "a" "g" "c" "a" "g" "g" "a" "g" "t" "c" "c" "t" "a" "g" "c" "c"
    [271] "c" "t" "t" "t" "c" "a" "c" "c" "c" "a" "t" "a" "c" "a" "c" "a" "t" "t"
    [289] "t" "g" "g" "c" "t" "c" "a" "g" "g" "g" "t" "t" "a" "c" "c" "g" "a" "a"
    [307] "g" "a" "g" "g" "g" "g" "c" "c" "a" "a" "g" "a" "a" "a" "t" "t" "a" "g"
    attr(,"name")
    [1] "KT152888"
    attr(,"Annot")
    [1] ">KT152888"
    attr(,"class")
    [1] "SeqFastadna"
    
    $KT152889.BRCA1
      [1] "a" "t" "g" "c" "t" "g" "a" "g" "t" "t" "t" "g" "t" "g" "t" "g" "t" "g"
     [19] "a" "a" "c" "g" "g" "a" "c" "a" "c" "t" "g" "a" "a" "a" "t" "a" "t" "t"
     [37] "t" "t" "c" "t" "a" "g" "g" "a" "a" "t" "t" "g" "c" "g" "g" "g" "a" "g"
     [55] "g" "a" "a" "a" "a" "t" "g" "g" "g" "t" "a" "g" "t" "t" "a" "g" "c" "t"
     [73] "a" "t" "t" "t" "c" "t" "g" "g" "g" "t" "g" "a" "c" "c" "c" "a" "g" "t"
     [91] "c" "t" "a" "t" "t" "a" "a" "a" "g" "a" "a" "a" "g" "a" "a" "a" "a" "a"
    [109] "t" "g" "c" "t" "g" "a" "a" "t" "g" "a" "g" "c" "a" "t" "g" "a" "t" "t"
    [127] "t" "t" "g" "a" "a" "g" "t" "c" "a" "g" "a" "g" "g" "a" "g" "a" "t" "g"
    [145] "t" "g" "g" "t" "c" "a" "a" "t" "g" "g" "a" "a" "g" "a" "a" "a" "c" "c"
    [163] "a" "c" "c" "a" "a" "g" "g" "t" "c" "c" "a" "t" "a" "a"
    attr(,"name")
    [1] "KT152889.BRCA1"
    attr(,"Annot")
    [1] ">KT152889.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $KT152890.BRCA1
       [1] "t" "t" "c" "t" "c" "t" "c" "c" "a" "t" "a" "t" "c" "t" "g" "a" "t" "t"
      [19] "t" "c" "a" "g" "a" "t" "a" "a" "c" "t" "t" "a" "g" "a" "a" "c" "a" "g"
      [37] "c" "c" "t" "a" "t" "g" "g" "g" "a" "a" "g" "t" "a" "g" "t" "c" "a" "t"
      [55] "g" "c" "a" "t" "c" "t" "c" "a" "g" "g" "t" "t" "t" "g" "t" "t" "c" "t"
      [73] "g" "a" "g" "a" "c" "a" "c" "c" "t" "g" "a" "t" "g" "a" "c" "c" "t" "g"
      [91] "t" "t" "a" "g" "a" "t" "g" "a" "t" "g" "g" "t" "g" "a" "a" "a" "t" "a"
     [109] "a" "a" "g" "g" "a" "a" "g" "a" "t" "a" "c" "t" "a" "g" "t" "t" "t" "t"
     [127] "g" "c" "t" "g" "a" "a" "a" "a" "t" "g" "a" "c" "a" "t" "t" "a" "a" "g"
     [145] "g" "a" "a" "a" "g" "t" "t" "c" "t" "g" "c" "t" "g" "t" "t" "t" "t" "t"
     [163] "a" "g" "c" "a" "a" "a" "a" "g" "c" "g" "t" "c" "c" "a" "g" "a" "g" "a"
     [181] "g" "g" "a" "g" "a" "g" "c" "t" "t" "a" "g" "c" "a" "g" "g" "a" "g" "t"
     [199] "c" "c" "t" "a" "g" "c" "c" "c" "t" "t" "t" "c" "a" "c" "c" "c" "a" "t"
     [217] "a" "c" "a" "c" "a" "t" "t" "t" "g" "g" "c" "t" "c" "a" "g" "g" "g" "t"
     [235] "t" "a" "c" "c" "g" "a" "a" "g" "a" "g" "g" "g" "g" "c" "c" "a" "a" "g"
     [253] "a" "a" "a" "t" "t" "a" "g" "a" "g" "t" "c" "c" "t" "c" "a" "g" "a" "a"
     [271] "g" "a" "g" "a" "a" "c" "t" "t" "a" "t" "c" "t" "a" "g" "t" "g" "a" "g"
     [289] "g" "a" "t" "g" "a" "a" "g" "a" "g" "c" "t" "t" "c" "c" "c" "t" "g" "c"
     [307] "t" "t" "c" "c" "a" "a" "c" "a" "c" "t" "t" "g" "t" "t" "a" "t" "t" "t"
     [325] "g" "g" "t" "a" "a" "a" "g" "t" "a" "a" "a" "c" "a" "a" "t" "a" "t" "a"
     [343] "c" "c" "t" "t" "c" "t" "c" "a" "g" "t" "c" "t" "a" "c" "t" "a" "g" "g"
     [361] "c" "a" "t" "a" "g" "c" "a" "c" "c" "g" "t" "t" "g" "c" "t" "a" "c" "c"
     [379] "g" "a" "g" "t" "g" "t" "c" "t" "g" "t" "c" "t" "a" "a" "g" "a" "a" "c"
     [397] "a" "c" "a" "g" "a" "g" "g" "a" "g" "a" "a" "t" "t" "t" "a" "t" "t" "a"
     [415] "t" "c" "a" "t" "t" "g" "a" "a" "g" "a" "a" "t" "a" "g" "c" "t" "t" "a"
     [433] "a" "a" "t" "g" "a" "c" "t" "g" "c" "a" "g" "t" "a" "a" "c" "c" "a" "g"
     [451] "g" "t" "a" "a" "t" "a" "t" "t" "g" "g" "c" "a" "a" "a" "g" "g" "c" "a"
     [469] "t" "c" "t" "c" "a" "g" "g" "a" "a" "c" "a" "t" "c" "a" "c" "c" "t" "t"
     [487] "a" "g" "t" "g" "a" "g" "g" "a" "a" "a" "c" "a" "a" "a" "a" "t" "g" "t"
     [505] "t" "c" "t" "g" "c" "t" "a" "g" "c" "t" "t" "g" "t" "t" "t" "t" "c" "t"
     [523] "t" "c" "a" "c" "a" "g" "t" "g" "c" "a" "g" "t" "g" "a" "a" "t" "t" "g"
     [541] "g" "a" "a" "g" "a" "c" "t" "t" "g" "a" "c" "t" "g" "c" "a" "a" "a" "t"
     [559] "a" "c" "a" "a" "a" "c" "a" "c" "c" "c" "a" "g" "g" "a" "t" "c" "c" "t"
     [577] "t" "t" "c" "t" "t" "g" "a" "t" "t" "g" "g" "t" "t" "c" "t" "t" "c" "c"
     [595] "a" "a" "a" "c" "a" "a" "a" "t" "g" "a" "g" "g" "c" "a" "t" "c" "a" "g"
     [613] "t" "c" "t" "g" "a" "a" "a" "g" "c" "c" "a" "g" "g" "g" "a" "g" "t" "t"
     [631] "g" "g" "t" "c" "t" "g" "a" "g" "t" "g" "a" "c" "a" "a" "g" "g" "a" "a"
     [649] "t" "t" "g" "g" "t" "t" "t" "c" "a" "g" "a" "t" "g" "a" "t" "g" "a" "a"
     [667] "g" "a" "a" "a" "g" "a" "g" "g" "a" "a" "c" "g" "g" "g" "c" "t" "t" "g"
     [685] "g" "a" "a" "g" "a" "a" "a" "a" "t" "a" "a" "t" "c" "a" "a" "g" "a" "a"
     [703] "g" "a" "g" "c" "a" "a" "a" "g" "c" "a" "t" "g" "g" "a" "t" "t" "c" "a"
     [721] "a" "a" "c" "t" "t" "a" "g" "g" "t" "g" "a" "a" "g" "c" "a" "g" "c" "a"
     [739] "t" "c" "t" "g" "g" "g" "t" "g" "t" "g" "a" "g" "a" "g" "t" "g" "a" "a"
     [757] "a" "c" "a" "a" "g" "c" "g" "t" "c" "t" "c" "t" "g" "a" "a" "g" "a" "c"
     [775] "t" "g" "c" "t" "c" "a" "g" "g" "g" "c" "t" "a" "t" "c" "c" "t" "c" "t"
     [793] "c" "a" "g" "a" "g" "t" "g" "a" "c" "a" "t" "t" "t" "t" "a" "a" "c" "c"
     [811] "a" "c" "t" "c" "a" "g" "c" "a" "g" "a" "g" "g" "g" "a" "t" "a" "c" "c"
     [829] "a" "t" "g" "c" "a" "a" "c" "a" "t" "a" "a" "c" "c" "t" "g" "a" "t" "a"
     [847] "a" "a" "g" "c" "t" "c" "c" "a" "g" "c" "a" "g" "g" "a" "a" "a" "t" "g"
     [865] "g" "c" "t" "g" "a" "a" "c" "t" "a" "g" "a" "a" "g" "c" "t" "g" "t" "g"
     [883] "t" "t" "a" "g" "a" "a" "c" "a" "g" "c" "a" "t" "g" "g" "g" "a" "g" "c"
     [901] "c" "a" "g" "c" "c" "t" "t" "c" "t" "a" "a" "c" "a" "g" "c" "t" "a" "c"
     [919] "c" "c" "t" "t" "c" "c" "a" "t" "c" "a" "t" "a" "a" "g" "t" "g" "a" "c"
     [937] "t" "c" "t" "t" "c" "t" "g" "c" "c" "c" "t" "t" "g" "a" "g" "g" "a" "c"
     [955] "c" "t" "g" "c" "g" "a" "a" "a" "t" "c" "c" "a" "g" "a" "a" "c" "a" "a"
     [973] "a" "g" "c" "a" "c" "a" "t" "c" "a" "g" "a" "a" "a" "a" "a" "g" "c" "a"
     [991] "g" "t" "a" "t" "t" "a" "a" "c" "t" "t" "c" "a" "c" "a" "g" "a" "a" "a"
    [1009] "a" "g" "t" "a" "g" "t" "g" "a" "a" "t" "a" "c" "c" "c" "t" "a" "t" "a"
    [1027] "a" "g" "c" "c" "a" "g" "a" "a" "t" "c" "c" "a" "g" "a" "a" "g" "g" "c"
    [1045] "c" "t" "t" "t" "c" "t" "g" "c" "t" "g" "a" "c" "a" "a" "g" "t" "t" "t"
    [1063] "g" "a" "g" "g" "t" "g" "t" "c" "t" "g" "c" "a" "g" "a" "t" "a" "g" "t"
    [1081] "t" "c" "t" "a" "c" "c" "a" "g" "t" "a" "a" "a" "a" "a" "t" "a" "a" "a"
    [1099] "g" "a" "a" "c" "c" "a" "g" "g" "a" "g" "t" "g" "g" "a" "a" "a" "g" "g"
    [1117] "t" "c" "a" "t" "c" "c" "c" "c" "t" "t" "c" "t" "a" "a" "a" "t" "g" "c"
    [1135] "c" "c" "a" "t" "c" "a" "t" "t" "a" "g" "a" "t" "g" "a" "t" "a" "g" "g"
    [1153] "t" "g" "g" "t" "a" "c" "a" "t" "g" "c" "a" "c" "a" "g" "t" "t" "g" "c"
    [1171] "t" "c" "t" "g" "g" "g" "a" "g" "t" "c" "t" "t" "c" "a" "g" "a" "a" "t"
    [1189] "a" "g" "a" "a" "a" "c" "t" "a" "c" "c" "c" "a" "t" "c" "t" "c" "a" "a"
    [1207] "g" "a" "g" "g" "a" "g" "c" "t" "c" "a" "t" "t" "a" "a" "g" "g" "t" "t"
    [1225] "g" "t" "t" "g" "a" "t" "g" "t" "g" "g" "a" "g" "g" "a" "g" "c" "a" "a"
    [1243] "c" "a" "g" "c" "t" "g" "g" "a" "a" "g" "a" "g" "t" "c" "t" "g" "g" "g"
    [1261] "c" "c" "a" "c" "a" "c" "g" "a" "t" "t" "t" "g" "a" "c" "g" "g" "a" "a"
    [1279] "a" "c" "a" "t" "c" "t" "t" "a" "c" "t" "t" "g" "c" "c" "a" "a" "g" "g"
    [1297] "c" "a" "a" "g" "a" "t" "c" "t" "a" "g" "a" "g" "g" "g" "a" "a" "c" "c"
    [1315] "c" "c" "t" "t" "a" "c" "c" "t" "g" "g" "a" "a" "t" "c" "t" "g" "g" "a"
    [1333] "a" "t" "c" "a" "g" "c" "c" "t" "c" "t" "t" "c" "t" "c" "t" "g" "a" "t"
    [1351] "g" "a" "c" "c" "c" "t" "g" "a" "a" "t" "c" "t" "g" "a" "t" "c" "c" "t"
    [1369] "t" "c" "t" "g" "a" "a" "g" "a" "c" "a" "g" "a" "g" "c" "c" "c" "c" "a"
    [1387] "g" "a" "g" "t" "c" "a" "g" "c" "t" "c" "g" "t" "g" "t" "t" "g" "g" "c"
    [1405] "a" "a" "c" "a" "t" "a" "c" "c" "a" "t" "c" "t" "t" "c" "a" "a" "c" "c"
    [1423] "t" "c" "t" "g" "c" "a" "t" "t" "g" "a" "a" "a" "g" "t" "t" "c" "c" "c"
    [1441] "c" "a" "a" "t" "t" "g" "a" "a" "a" "g" "t" "t" "g" "c" "a" "g" "a" "a"
    [1459] "t" "c" "t" "g" "c" "c" "c" "a" "g" "a" "g" "t" "c" "c" "a" "g" "c" "t"
    [1477] "g" "c" "t" "g" "c" "t" "c" "a" "t" "a" "c" "t" "a" "c" "t" "g" "a" "t"
    [1495] "a" "c" "t" "g" "c" "t" "g" "g" "g" "t" "a" "t" "a" "a" "t" "g" "c" "a"
    [1513] "a" "t" "g" "g" "a" "a" "g" "a" "a" "a" "g" "t" "g" "t" "g" "a" "g" "c"
    [1531] "a" "g" "g" "g" "a" "g" "a" "a" "g" "c" "c" "a" "g" "a" "a" "t" "t" "g"
    [1549] "a" "c" "a" "g" "c" "t" "t" "c" "a" "a" "c" "a" "g" "a" "a" "a" "g" "g"
    [1567] "g" "t" "c" "a" "a" "c" "a" "a" "a" "a" "g" "a" "a" "t" "g" "t" "c" "c"
    [1585] "a" "t" "g" "g" "t" "g" "g" "t" "g" "t" "c" "t" "g" "g" "c" "c" "t" "g"
    [1603] "a" "c" "c" "c" "c" "a" "g" "a" "a" "g" "a" "a" "t" "t" "t" "a" "t" "g"
    [1621] "c" "t" "c" "g" "t" "g" "t" "a" "c" "a" "a" "g" "t" "t" "t" "g" "c" "c"
    [1639] "a" "g" "a" "a" "a" "a" "c" "a" "c" "c" "a" "c" "a" "t" "c" "a" "c" "t"
    [1657] "t" "t" "a" "a" "c" "t" "a" "a" "t" "c" "t" "a" "a" "t" "t" "a" "c" "t"
    [1675] "g" "a" "a" "g" "a" "g" "a" "c" "t" "a" "c" "t" "c" "a" "t" "g" "t" "t"
    [1693] "g" "t" "t" "a" "t" "g" "a" "a" "a" "a" "c" "a" "g" "a" "t" "g" "c" "t"
    [1711] "g" "a" "g" "t" "t" "t" "g" "t" "g" "t" "g" "t" "g" "a" "a" "c" "g" "g"
    [1729] "a" "c" "a" "c" "t" "g" "a" "a" "a" "t" "a" "t" "t" "t" "t" "c" "t" "a"
    [1747] "g" "g" "a" "a" "t" "t" "g" "c" "g" "g" "g" "a" "g" "g" "a" "a" "a" "a"
    [1765] "t" "g" "g" "g" "t" "a" "g" "t" "t" "a" "g" "c" "t" "a" "t" "t" "t" "c"
    [1783] "t" "g" "g" "g" "t" "g" "a" "c" "c" "c" "a" "g" "t" "c" "t" "a" "t" "t"
    [1801] "a" "a" "a" "g" "a" "a" "a" "g" "a" "a" "a" "a" "a" "t" "g" "c" "t" "g"
    [1819] "a" "a" "t" "g" "a" "g" "c" "a" "t" "g" "a" "t" "t" "t" "t" "g" "a" "a"
    [1837] "g" "t" "c" "a" "g" "a" "g" "g" "a" "g" "a" "t" "g" "t" "g" "g" "t" "c"
    [1855] "a" "a" "t" "g" "g" "a" "a" "g" "a" "a" "a" "c" "c" "a" "c" "c" "a" "a"
    [1873] "g" "g" "t" "c" "c" "a" "a" "a" "g" "c" "g" "a" "g" "c" "a" "a" "g" "a"
    [1891] "g" "a" "a" "t" "c" "c" "c" "a" "g" "g" "a" "c" "a" "g" "a" "a" "a" "g"
    [1909] "a" "t" "c" "t" "t" "c" "a" "g" "g" "g" "g" "g" "c" "t" "a" "g" "a" "a"
    [1927] "a" "t" "c" "t" "g" "t" "t" "g" "c" "t" "a" "t" "g" "g" "g" "c" "c" "c"
    [1945] "t" "t" "c" "a" "c" "c" "a" "a" "c" "a" "t" "g" "c" "c" "c" "a" "c" "a"
    [1963] "g" "a" "t" "c" "a" "a" "c" "t" "g" "g" "a" "a" "t" "g" "g" "a" "t" "g"
    [1981] "g" "t" "a" "c" "a" "g" "c" "t" "g" "t" "g" "t" "g" "g" "t" "g" "c" "t"
    [1999] "t" "c" "t" "g" "t" "g" "g" "t" "g" "a" "a" "g" "g" "a" "g" "c" "t" "t"
    [2017] "t" "c" "a" "t" "c" "a" "t" "t" "c" "a" "c" "c" "c" "t" "t" "g" "g" "c"
    [2035] "a" "c" "a" "g" "g" "t" "g" "t" "c" "c" "a" "c" "c" "c" "a" "a" "t" "t"
    [2053] "g" "t" "g" "g" "t" "t" "g" "t" "g" "c" "a" "g" "c" "c" "a" "g" "a" "t"
    [2071] "g" "c" "c" "t" "g" "g" "a" "c" "a" "g" "a" "g" "g" "a" "c" "a" "a" "t"
    [2089] "g" "g" "c" "t" "t" "c" "c" "a" "t" "g" "c" "a" "a" "t" "t" "g" "g" "g"
    [2107] "c" "a" "g" "a" "t" "g" "t" "g" "t" "g" "a" "g" "g" "c" "a" "c" "c" "t"
    [2125] "g" "t" "g" "g" "t" "g" "a" "c" "c" "c" "g" "a" "g" "a" "g" "t" "g" "g"
    [2143] "g" "t" "g" "t" "t" "g" "g" "a" "c" "a" "g" "t" "g" "t" "a" "g" "c" "a"
    [2161] "c" "t" "c" "t" "a" "c" "c" "a" "g" "t" "g" "c" "c" "a" "g" "g" "a" "g"
    [2179] "c" "t" "g" "g" "a" "c" "a" "c" "c" "t" "a" "c" "c" "t" "g" "a" "t" "a"
    [2197] "c" "c" "c" "c" "a" "g" "a" "t" "c" "c" "c" "c" "c" "a" "c" "a" "g" "c"
    [2215] "c" "a" "c" "t" "a" "c" "t" "g" "a"
    attr(,"name")
    [1] "KT152890.BRCA1"
    attr(,"Annot")
    [1] ">KT152890.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $KT844468.BRCA1
     [1] "a" "t" "g" "g" "a" "t" "t" "t" "a" "t" "c" "t" "g" "c" "t" "c" "t" "t" "c"
    [20] "g" "c" "g" "t" "t" "g" "a" "a" "g" "a" "a" "g" "t" "a" "c" "a" "a" "a" "a"
    [39] "t" "g" "t" "c" "a" "t" "t" "a" "a" "t" "g" "c" "t" "a" "g" "c" "a" "g" "a"
    [58] "a" "a" "a" "t" "c" "t" "t" "a" "g"
    attr(,"name")
    [1] "KT844468.BRCA1"
    attr(,"Annot")
    [1] ">KT844468.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $KT844469.BRCA1
     [1] "g" "a" "g" "c" "c" "t" "a" "c" "a" "a" "g" "a" "a" "a" "g" "t" "a" "c" "g"
    [20] "a" "g" "a" "t" "t" "t" "a" "g" "t" "c" "a" "a" "c" "t" "t" "g" "t" "t" "g"
    [39] "a" "a" "g" "a" "g" "c" "t" "a" "t" "t" "g" "a" "a" "a" "a" "t" "c" "a" "t"
    [58] "t" "t" "g" "t" "g" "c" "t" "t" "t" "t" "c" "a" "g" "c" "t" "t" "g" "a" "c"
    [77] "a" "c" "a" "g" "g" "t" "t" "t" "c" "t" "c" "a" "a"
    attr(,"name")
    [1] "KT844469.BRCA1"
    attr(,"Annot")
    [1] ">KT844469.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $KU359055
     [1] "g" "a" "g" "c" "a" "t" "g" "a" "t" "t" "t" "t" "g" "a" "a" "g" "a" "c" "a"
    [20] "g" "a" "g" "g" "a" "g" "a" "t" "g" "t" "g" "g" "t" "c" "a" "a" "t" "g" "g"
    [39] "a" "a" "g" "a" "a" "a" "c" "c" "a" "c" "c" "a" "a" "g" "g" "t" "c" "c" "a"
    [58] "a" "a" "g" "c" "g" "a" "g" "c" "a" "a" "g" "a" "g" "a" "a" "t" "c" "c" "c"
    [77] "a" "g" "g" "a" "c" "a" "g" "a" "a" "a" "g" "a"
    attr(,"name")
    [1] "KU359055"
    attr(,"Annot")
    [1] ">KU359055"
    attr(,"class")
    [1] "SeqFastadna"
    
    $KU359056
     [1] "g" "a" "g" "c" "a" "t" "g" "a" "t" "t" "t" "t" "g" "a" "a" "g" "a" "c" "a"
    [20] "g" "a" "g" "g" "a" "g" "a" "t" "g" "t" "g" "g" "t" "c" "a" "a" "t" "g" "g"
    [39] "a" "a" "g" "a" "a" "a" "c" "c" "a" "c" "c" "a" "a" "g" "g" "t" "c" "c" "a"
    [58] "a" "a" "g" "c" "g" "a" "g" "c" "a" "a" "g" "a" "g" "a" "a" "t" "c" "c" "c"
    [77] "a" "g" "g" "a" "c" "a" "g" "a" "a" "a" "g" "a"
    attr(,"name")
    [1] "KU359056"
    attr(,"Annot")
    [1] ">KU359056"
    attr(,"class")
    [1] "SeqFastadna"
    
    $KU359057
     [1] "g" "a" "g" "c" "a" "t" "g" "a" "t" "t" "t" "t" "g" "a" "a" "g" "a" "c" "a"
    [20] "g" "a" "g" "g" "a" "g" "a" "t" "g" "t" "g" "g" "t" "c" "a" "a" "t" "g" "g"
    [39] "a" "a" "g" "a" "a" "a" "c" "c" "a" "c" "c" "a" "a" "g" "g" "t" "c" "c" "a"
    [58] "a" "a" "g" "c" "g" "a" "g" "c" "a" "a" "g" "a" "g" "a" "a" "t" "c" "c" "c"
    [77] "a" "g" "g" "a" "c" "a" "g" "a" "a" "a" "g" "a"
    attr(,"name")
    [1] "KU359057"
    attr(,"Annot")
    [1] ">KU359057"
    attr(,"class")
    [1] "SeqFastadna"
    
    $KU359058
     [1] "g" "a" "g" "c" "a" "t" "g" "a" "t" "t" "t" "t" "g" "a" "a" "g" "t" "c" "a"
    [20] "g" "a" "g" "g" "a" "g" "a" "t" "g" "t" "g" "g" "t" "c" "a" "a" "t" "g" "g"
    [39] "a" "a" "g" "a" "a" "a" "c" "c" "a" "c" "c" "a" "a" "g" "g" "t" "c" "c" "a"
    [58] "a" "a" "g" "c" "g" "a" "g" "c" "a" "a" "g" "a" "g" "a" "a" "t" "c" "c" "c"
    [77] "a" "g" "g" "a" "c" "a" "g" "a" "a" "a" "g" "a"
    attr(,"name")
    [1] "KU359058"
    attr(,"Annot")
    [1] ">KU359058"
    attr(,"class")
    [1] "SeqFastadna"
    
    $KU359059
     [1] "g" "a" "g" "c" "a" "t" "g" "a" "t" "t" "t" "t" "g" "a" "a" "g" "t" "c" "a"
    [20] "g" "a" "g" "g" "a" "g" "a" "t" "g" "t" "g" "g" "t" "c" "a" "a" "t" "g" "g"
    [39] "a" "a" "g" "a" "a" "a" "c" "c" "a" "c" "c" "a" "a" "g" "g" "t" "c" "c" "a"
    [58] "a" "a" "g" "c" "g" "a" "g" "c" "a" "a" "g" "a" "g" "a" "a" "t" "c" "c" "c"
    [77] "a" "g" "g" "a" "c" "a" "g" "a" "a" "a" "g" "a"
    attr(,"name")
    [1] "KU359059"
    attr(,"Annot")
    [1] ">KU359059"
    attr(,"class")
    [1] "SeqFastadna"
    
    $KU359060
     [1] "g" "a" "g" "c" "a" "t" "g" "a" "t" "t" "t" "t" "g" "a" "a" "g" "t" "c" "a"
    [20] "g" "a" "g" "g" "a" "g" "a" "t" "g" "t" "g" "g" "t" "c" "a" "a" "t" "g" "g"
    [39] "a" "a" "g" "a" "a" "a" "c" "c" "a" "c" "c" "a" "a" "g" "g" "t" "c" "c" "a"
    [58] "a" "a" "g" "c" "g" "a" "g" "c" "a" "a" "g" "a" "g" "a" "a" "t" "c" "c" "c"
    [77] "a" "g" "g" "a" "c" "a" "g" "a" "a" "a" "g" "a"
    attr(,"name")
    [1] "KU359060"
    attr(,"Annot")
    [1] ">KU359060"
    attr(,"class")
    [1] "SeqFastadna"
    
    $KU359061
     [1] "g" "a" "g" "c" "a" "t" "g" "a" "t" "t" "t" "t" "g" "a" "a" "g" "t" "c" "a"
    [20] "g" "a" "g" "g" "a" "g" "a" "t" "g" "t" "g" "g" "t" "c" "a" "a" "t" "g" "g"
    [39] "a" "a" "g" "a" "a" "a" "c" "c" "a" "c" "c" "a" "a" "g" "g" "t" "c" "c" "a"
    [58] "a" "a" "g" "c" "g" "a" "g" "c" "a" "a" "g" "a" "g" "a" "a" "t" "c" "c" "c"
    [77] "a" "g" "g" "a" "c" "a" "g" "a" "a" "a" "g" "a"
    attr(,"name")
    [1] "KU359061"
    attr(,"Annot")
    [1] ">KU359061"
    attr(,"class")
    [1] "SeqFastadna"
    
    $KU359062
     [1] "g" "a" "g" "c" "a" "t" "g" "a" "t" "t" "t" "t" "g" "a" "a" "g" "t" "c" "a"
    [20] "g" "a" "g" "g" "a" "g" "a" "t" "g" "t" "g" "g" "t" "c" "a" "a" "t" "g" "g"
    [39] "a" "a" "g" "a" "a" "a" "c" "c" "a" "c" "c" "a" "a" "g" "g" "t" "c" "c" "a"
    [58] "a" "a" "g" "c" "g" "a" "g" "c" "a" "a" "g" "a" "g" "a" "a" "t" "c" "c" "c"
    [77] "a" "g" "g" "a" "c" "a" "g" "a" "a" "a" "g" "a"
    attr(,"name")
    [1] "KU359062"
    attr(,"Annot")
    [1] ">KU359062"
    attr(,"class")
    [1] "SeqFastadna"
    
    $KU359063
     [1] "g" "a" "g" "c" "a" "t" "g" "a" "t" "t" "t" "t" "g" "a" "a" "g" "t" "c" "a"
    [20] "g" "a" "g" "g" "a" "g" "a" "t" "g" "t" "g" "g" "t" "c" "a" "a" "t" "g" "g"
    [39] "a" "a" "g" "a" "a" "a" "c" "c" "a" "c" "c" "a" "a" "g" "g" "t" "c" "c" "a"
    [58] "a" "a" "g" "c" "g" "a" "g" "c" "a" "a" "g" "a" "g" "a" "a" "t" "c" "c" "c"
    [77] "a" "g" "g" "a" "c" "a" "g" "a" "a" "a" "g" "a"
    attr(,"name")
    [1] "KU359063"
    attr(,"Annot")
    [1] ">KU359063"
    attr(,"class")
    [1] "SeqFastadna"
    
    $KU359064
     [1] "g" "a" "g" "c" "a" "t" "g" "a" "t" "t" "t" "t" "g" "a" "a" "g" "t" "c" "a"
    [20] "g" "a" "g" "g" "a" "g" "a" "t" "g" "t" "g" "g" "t" "c" "a" "a" "t" "g" "g"
    [39] "a" "a" "g" "a" "a" "a" "c" "c" "a" "c" "c" "a" "a" "g" "g" "t" "c" "c" "a"
    [58] "a" "a" "g" "c" "g" "a" "g" "c" "a" "a" "g" "a" "g" "a" "a" "t" "c" "c" "c"
    [77] "a" "g" "g" "a" "c" "a" "g" "a" "a" "a" "g" "a"
    attr(,"name")
    [1] "KU359064"
    attr(,"Annot")
    [1] ">KU359064"
    attr(,"class")
    [1] "SeqFastadna"
    
    $KX580312.BRCA1
      [1] "g" "t" "c" "a" "t" "c" "c" "c" "c" "t" "t" "c" "t" "a" "a" "a" "t" "g"
     [19] "c" "c" "c" "a" "t" "c" "a" "t" "t" "a" "g" "a" "t" "g" "a" "t" "a" "g"
     [37] "g" "t" "g" "g" "t" "a" "c" "a" "t" "g" "c" "a" "c" "a" "g" "t" "t" "g"
     [55] "c" "t" "c" "t" "g" "g" "g" "a" "g" "t" "c" "t" "t" "c" "a" "g" "a" "a"
     [73] "t" "a" "g" "a" "a" "a" "c" "t" "a" "c" "c" "c" "a" "t" "c" "t" "c" "a"
     [91] "a" "g" "a" "g" "g" "a" "g" "c" "t" "c" "a" "t" "t" "a" "a" "g" "g" "t"
    [109] "t" "g" "t" "t" "g" "a" "t" "g" "t" "g" "g" "a" "g" "g" "a" "g" "t" "a"
    [127] "a"
    attr(,"name")
    [1] "KX580312.BRCA1"
    attr(,"Annot")
    [1] ">KX580312.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $KX944478.BRCA1
       [1] "c" "t" "g" "c" "t" "t" "g" "t" "g" "a" "a" "t" "t" "t" "t" "c" "t" "g"
      [19] "a" "g" "a" "c" "g" "g" "a" "t" "g" "t" "a" "a" "c" "a" "a" "a" "t" "a"
      [37] "c" "t" "g" "a" "a" "c" "a" "t" "c" "a" "t" "c" "a" "a" "c" "c" "c" "a"
      [55] "g" "t" "a" "a" "t" "a" "a" "t" "g" "a" "t" "t" "t" "g" "a" "a" "c" "a"
      [73] "c" "c" "a" "c" "t" "g" "a" "g" "a" "a" "g" "c" "g" "t" "g" "c" "a" "g"
      [91] "c" "t" "g" "a" "g" "a" "g" "g" "c" "a" "t" "c" "c" "a" "g" "a" "a" "a"
     [109] "a" "g" "t" "a" "t" "c" "a" "g" "g" "g" "t" "a" "g" "t" "t" "c" "t" "g"
     [127] "t" "t" "t" "c" "a" "a" "a" "c" "t" "t" "g" "c" "a" "t" "g" "t" "g" "g"
     [145] "a" "g" "c" "c" "a" "t" "g" "t" "g" "g" "c" "a" "c" "a" "a" "a" "t" "a"
     [163] "c" "t" "c" "a" "t" "g" "c" "c" "a" "g" "c" "t" "c" "a" "t" "t" "a" "c"
     [181] "a" "g" "c" "a" "t" "g" "a" "g" "a" "a" "c" "a" "g" "c" "a" "g" "t" "t"
     [199] "t" "a" "t" "t" "a" "c" "t" "c" "a" "c" "t" "a" "a" "a" "g" "a" "c" "a"
     [217] "g" "a" "a" "t" "g" "a" "a" "t" "g" "t" "a" "g" "a" "a" "a" "a" "g" "g"
     [235] "c" "t" "g" "a" "a" "t" "t" "c" "t" "g" "t" "a" "a" "t" "a" "a" "a" "a"
     [253] "g" "c" "a" "a" "a" "c" "a" "g" "c" "c" "t" "g" "g" "c" "t" "t" "a" "g"
     [271] "c" "a" "a" "g" "g" "a" "g" "c" "c" "a" "a" "c" "a" "t" "a" "a" "c" "a"
     [289] "g" "a" "t" "g" "g" "g" "c" "t" "g" "g" "a" "a" "g" "t" "a" "a" "g" "g"
     [307] "a" "a" "a" "c" "a" "t" "g" "t" "a" "a" "t" "g" "a" "t" "a" "g" "g" "c"
     [325] "g" "g" "a" "c" "t" "c" "c" "c" "a" "g" "c" "a" "c" "a" "g" "a" "a" "a"
     [343] "a" "a" "a" "a" "g" "g" "t" "a" "g" "a" "t" "c" "t" "g" "a" "a" "t" "g"
     [361] "c" "t" "g" "a" "t" "c" "c" "c" "c" "t" "g" "t" "g" "t" "g" "a" "g" "a"
     [379] "g" "a" "a" "a" "a" "g" "a" "a" "t" "g" "g" "a" "a" "t" "a" "a" "g" "c"
     [397] "a" "g" "a" "a" "a" "c" "t" "g" "c" "c" "a" "t" "g" "c" "t" "c" "a" "g"
     [415] "a" "g" "a" "a" "t" "c" "c" "t" "a" "g" "a" "g" "a" "t" "a" "c" "t" "g"
     [433] "a" "a" "g" "a" "t" "g" "t" "t" "c" "c" "t" "t" "g" "g" "a" "t" "a" "a"
     [451] "c" "a" "c" "t" "a" "a" "a" "t" "a" "g" "c" "a" "g" "c" "a" "t" "t" "c"
     [469] "a" "g" "a" "a" "a" "g" "t" "t" "a" "a" "t" "g" "a" "g" "t" "g" "g" "t"
     [487] "t" "t" "t" "c" "c" "a" "g" "a" "a" "g" "t" "g" "a" "t" "g" "a" "a" "c"
     [505] "t" "g" "t" "t" "a" "g" "g" "t" "t" "c" "t" "g" "a" "t" "g" "a" "c" "t"
     [523] "c" "a" "c" "a" "t" "g" "a" "t" "g" "g" "g" "g" "a" "g" "t" "c" "t" "g"
     [541] "a" "a" "t" "c" "a" "a" "a" "t" "g" "c" "c" "a" "a" "a" "g" "t" "a" "g"
     [559] "c" "t" "g" "a" "t" "g" "t" "a" "t" "t" "g" "g" "a" "c" "g" "t" "t" "c"
     [577] "t" "a" "a" "a" "t" "g" "a" "g" "g" "t" "a" "g" "a" "t" "g" "a" "a" "t"
     [595] "a" "t" "t" "c" "t" "g" "g" "t" "t" "c" "t" "t" "c" "a" "g" "a" "g" "a"
     [613] "a" "a" "a" "t" "a" "g" "a" "c" "t" "t" "a" "c" "t" "g" "g" "c" "c" "a"
     [631] "g" "t" "g" "a" "t" "c" "c" "t" "c" "a" "t" "g" "a" "g" "g" "c" "t" "t"
     [649] "t" "a" "a" "t" "a" "t" "g" "t" "a" "a" "a" "a" "g" "t" "g" "a" "a" "a"
     [667] "g" "a" "g" "t" "t" "c" "a" "c" "t" "c" "c" "a" "a" "a" "t" "c" "a" "g"
     [685] "t" "a" "g" "a" "g" "a" "g" "t" "a" "a" "t" "a" "t" "t" "g" "a" "a" "g"
     [703] "a" "c" "a" "a" "a" "a" "t" "a" "t" "t" "t" "g" "g" "g" "a" "a" "a" "a"
     [721] "c" "c" "t" "a" "t" "c" "g" "g" "a" "a" "g" "a" "a" "g" "g" "c" "a" "a"
     [739] "g" "c" "c" "t" "c" "c" "c" "c" "a" "a" "c" "t" "t" "a" "a" "g" "c" "c"
     [757] "a" "t" "g" "t" "a" "a" "c" "t" "g" "a" "a" "a" "a" "t" "c" "t" "a" "a"
     [775] "t" "t" "a" "t" "a" "g" "g" "a" "g" "c" "a" "t" "t" "t" "g" "t" "t" "a"
     [793] "c" "t" "g" "a" "g" "c" "c" "a" "c" "a" "g" "a" "t" "a" "a" "t" "a" "c"
     [811] "a" "a" "g" "a" "g" "c" "g" "t" "c" "c" "c" "c" "t" "c" "a" "c" "a" "a"
     [829] "a" "t" "a" "a" "a" "t" "t" "a" "a" "a" "g" "c" "g" "t" "a" "a" "a" "a"
     [847] "g" "g" "a" "g" "a" "c" "c" "t" "a" "c" "a" "t" "c" "a" "g" "g" "c" "c"
     [865] "t" "t" "c" "a" "t" "c" "c" "t" "g" "a" "g" "g" "a" "t" "t" "t" "t" "a"
     [883] "t" "c" "a" "a" "g" "a" "a" "a" "g" "c" "a" "g" "a" "t" "t" "t" "g" "g"
     [901] "c" "a" "g" "t" "t" "c" "a" "a" "a" "a" "g" "a" "c" "t" "c" "c" "t" "g"
     [919] "a" "a" "a" "t" "g" "a" "t" "a" "a" "a" "t" "c" "a" "g" "g" "g" "a" "a"
     [937] "c" "t" "a" "a" "c" "c" "a" "a" "a" "c" "g" "g" "a" "g" "c" "a" "g" "a"
     [955] "a" "t" "g" "g" "t" "c" "a" "a" "g" "t" "g" "a" "t" "g" "a" "a" "t" "a"
     [973] "t" "t" "a" "c" "t" "a" "a" "t" "a" "g" "t" "g" "g" "t" "c" "a" "t" "g"
     [991] "a" "g" "a" "a" "t" "a" "a" "a" "a" "c" "a" "a" "a" "a" "g" "g" "t" "g"
    [1009] "a" "t" "t" "c" "t" "a" "a" "t" "t" "c" "a" "g" "a" "a" "t" "g" "a"
    attr(,"name")
    [1] "KX944478.BRCA1"
    attr(,"Annot")
    [1] ">KX944478.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $L78833.BRCA1
       [1] "a" "t" "g" "g" "a" "t" "t" "t" "a" "t" "c" "t" "g" "c" "t" "c" "t" "t"
      [19] "c" "g" "c" "g" "t" "t" "g" "a" "a" "g" "a" "a" "g" "t" "a" "c" "a" "a"
      [37] "a" "a" "t" "g" "t" "c" "a" "t" "t" "a" "a" "t" "g" "c" "t" "a" "t" "g"
      [55] "c" "a" "g" "a" "a" "a" "a" "t" "c" "t" "t" "a" "g" "a" "g" "t" "g" "t"
      [73] "c" "c" "c" "a" "t" "c" "t" "g" "t" "c" "t" "g" "g" "a" "g" "t" "t" "g"
      [91] "a" "t" "c" "a" "a" "g" "g" "a" "a" "c" "c" "t" "g" "t" "c" "t" "c" "c"
     [109] "a" "c" "a" "a" "a" "g" "t" "g" "t" "g" "a" "c" "c" "a" "c" "a" "t" "a"
     [127] "t" "t" "t" "t" "g" "c" "a" "a" "a" "t" "t" "t" "t" "g" "c" "a" "t" "g"
     [145] "c" "t" "g" "a" "a" "a" "c" "t" "t" "c" "t" "c" "a" "a" "c" "c" "a" "g"
     [163] "a" "a" "g" "a" "a" "a" "g" "g" "g" "c" "c" "t" "t" "c" "a" "c" "a" "g"
     [181] "t" "g" "t" "c" "c" "t" "t" "t" "a" "t" "g" "t" "a" "a" "g" "a" "a" "t"
     [199] "g" "a" "t" "a" "t" "a" "a" "c" "c" "a" "a" "a" "a" "g" "g" "a" "g" "c"
     [217] "c" "t" "a" "c" "a" "a" "g" "a" "a" "a" "g" "t" "a" "c" "g" "a" "g" "a"
     [235] "t" "t" "t" "a" "g" "t" "c" "a" "a" "c" "t" "t" "g" "t" "t" "g" "a" "a"
     [253] "g" "a" "g" "c" "t" "a" "t" "t" "g" "a" "a" "a" "a" "t" "c" "a" "t" "t"
     [271] "t" "g" "t" "g" "c" "t" "t" "t" "t" "c" "a" "g" "c" "t" "t" "g" "a" "c"
     [289] "a" "c" "a" "g" "g" "t" "t" "t" "g" "g" "a" "g" "t" "a" "t" "g" "c" "a"
     [307] "a" "a" "c" "a" "g" "c" "t" "a" "t" "a" "a" "t" "t" "t" "t" "g" "c" "a"
     [325] "a" "a" "a" "a" "a" "g" "g" "a" "a" "a" "a" "t" "a" "a" "c" "t" "c" "t"
     [343] "c" "c" "t" "g" "a" "a" "c" "a" "t" "c" "t" "a" "a" "a" "a" "g" "a" "t"
     [361] "g" "a" "a" "g" "t" "t" "t" "c" "t" "a" "t" "c" "a" "t" "c" "c" "a" "a"
     [379] "a" "g" "t" "a" "t" "g" "g" "g" "c" "t" "a" "c" "a" "g" "a" "a" "a" "c"
     [397] "c" "g" "t" "g" "c" "c" "a" "a" "a" "a" "g" "a" "c" "t" "t" "c" "t" "a"
     [415] "c" "a" "g" "a" "g" "t" "g" "a" "a" "c" "c" "c" "g" "a" "a" "a" "a" "t"
     [433] "c" "c" "t" "t" "c" "c" "t" "t" "g" "c" "a" "g" "g" "a" "a" "a" "c" "c"
     [451] "a" "g" "t" "c" "t" "c" "a" "g" "t" "g" "t" "c" "c" "a" "a" "c" "t" "c"
     [469] "t" "c" "t" "a" "a" "c" "c" "t" "t" "g" "g" "a" "a" "c" "t" "g" "t" "g"
     [487] "a" "g" "a" "a" "c" "t" "c" "t" "g" "a" "g" "g" "a" "c" "a" "a" "a" "g"
     [505] "c" "a" "g" "c" "g" "g" "a" "t" "a" "c" "a" "a" "c" "c" "t" "c" "a" "a"
     [523] "a" "a" "g" "a" "c" "g" "t" "c" "t" "g" "t" "c" "t" "a" "c" "a" "t" "t"
     [541] "g" "a" "a" "t" "t" "g" "g" "g" "a" "t" "c" "t" "g" "a" "t" "t" "c" "t"
     [559] "t" "c" "t" "g" "a" "a" "g" "a" "t" "a" "c" "c" "g" "t" "t" "a" "a" "t"
     [577] "a" "a" "g" "g" "c" "a" "a" "c" "t" "t" "a" "t" "t" "g" "c" "a" "g" "t"
     [595] "g" "t" "g" "g" "g" "a" "g" "a" "t" "c" "a" "a" "g" "a" "a" "t" "t" "g"
     [613] "t" "t" "a" "c" "a" "a" "a" "t" "c" "a" "c" "c" "c" "c" "t" "c" "a" "a"
     [631] "g" "g" "a" "a" "c" "c" "a" "g" "g" "g" "a" "t" "g" "a" "a" "a" "t" "c"
     [649] "a" "g" "t" "t" "t" "g" "g" "a" "t" "t" "c" "t" "g" "c" "a" "a" "a" "a"
     [667] "a" "a" "g" "g" "c" "t" "g" "c" "t" "t" "g" "t" "g" "a" "a" "t" "t" "t"
     [685] "t" "c" "t" "g" "a" "g" "a" "c" "g" "g" "a" "t" "g" "t" "a" "a" "c" "a"
     [703] "a" "a" "t" "a" "c" "t" "g" "a" "a" "c" "a" "t" "c" "a" "t" "c" "a" "a"
     [721] "c" "c" "c" "a" "g" "t" "a" "a" "t" "a" "a" "t" "g" "a" "t" "t" "t" "g"
     [739] "a" "a" "c" "a" "c" "c" "a" "c" "t" "g" "a" "g" "a" "a" "g" "c" "g" "t"
     [757] "g" "c" "a" "g" "c" "t" "g" "a" "g" "a" "g" "g" "c" "a" "t" "c" "c" "a"
     [775] "g" "a" "a" "a" "a" "g" "t" "a" "t" "c" "a" "g" "g" "g" "t" "a" "g" "t"
     [793] "t" "c" "t" "g" "t" "t" "t" "c" "a" "a" "a" "c" "t" "t" "g" "c" "a" "t"
     [811] "g" "t" "g" "g" "a" "g" "c" "c" "a" "t" "g" "t" "g" "g" "c" "a" "c" "a"
     [829] "a" "a" "t" "a" "c" "t" "c" "a" "t" "g" "c" "c" "a" "g" "c" "t" "c" "a"
     [847] "t" "t" "a" "c" "a" "g" "c" "a" "t" "g" "a" "g" "a" "a" "c" "a" "g" "c"
     [865] "a" "g" "t" "t" "t" "a" "t" "t" "a" "c" "t" "c" "a" "c" "t" "a" "a" "a"
     [883] "g" "a" "c" "a" "g" "a" "a" "t" "g" "a" "a" "t" "g" "t" "a" "g" "a" "a"
     [901] "a" "a" "g" "g" "c" "t" "g" "a" "a" "t" "t" "c" "t" "g" "t" "a" "a" "t"
     [919] "a" "a" "a" "a" "g" "c" "a" "a" "a" "c" "a" "g" "c" "c" "t" "g" "g" "c"
     [937] "t" "t" "a" "g" "c" "a" "a" "g" "g" "a" "g" "c" "c" "a" "a" "c" "a" "t"
     [955] "a" "a" "c" "a" "g" "a" "t" "g" "g" "g" "c" "t" "g" "g" "a" "a" "g" "t"
     [973] "a" "a" "g" "g" "a" "a" "a" "c" "a" "t" "g" "t" "a" "a" "t" "g" "a" "t"
     [991] "a" "g" "g" "c" "g" "g" "a" "c" "t" "c" "c" "c" "a" "g" "c" "a" "c" "a"
    [1009] "g" "a" "a" "a" "a" "a" "a" "a" "g" "g" "t" "a" "g" "a" "t" "c" "t" "g"
    [1027] "a" "a" "t" "g" "c" "t" "g" "a" "t" "c" "c" "c" "c" "t" "g" "t" "g" "t"
    [1045] "g" "a" "g" "a" "g" "a" "a" "a" "a" "g" "a" "a" "t" "g" "g" "a" "a" "t"
    [1063] "a" "a" "g" "c" "a" "g" "a" "a" "a" "c" "t" "g" "c" "c" "a" "t" "g" "c"
    [1081] "t" "c" "a" "g" "a" "g" "a" "a" "t" "c" "c" "t" "a" "g" "a" "g" "a" "t"
    [1099] "a" "c" "t" "g" "a" "a" "g" "a" "t" "g" "t" "t" "c" "c" "t" "t" "g" "g"
    [1117] "a" "t" "a" "a" "c" "a" "c" "t" "a" "a" "a" "t" "a" "g" "c" "a" "g" "c"
    [1135] "a" "t" "t" "c" "a" "g" "a" "a" "a" "g" "t" "t" "a" "a" "t" "g" "a" "g"
    [1153] "t" "g" "g" "t" "t" "t" "t" "c" "c" "a" "g" "a" "a" "g" "t" "g" "a" "t"
    [1171] "g" "a" "a" "c" "t" "g" "t" "t" "a" "g" "g" "t" "t" "c" "t" "g" "a" "t"
    [1189] "g" "a" "c" "t" "c" "a" "c" "a" "t" "g" "a" "t" "g" "g" "g" "g" "a" "g"
    [1207] "t" "c" "t" "g" "a" "a" "t" "c" "a" "a" "a" "t" "g" "c" "c" "a" "a" "a"
    [1225] "g" "t" "a" "g" "c" "t" "g" "a" "t" "g" "t" "a" "t" "t" "g" "g" "a" "c"
    [1243] "g" "t" "t" "c" "t" "a" "a" "a" "t" "g" "a" "g" "g" "t" "a" "g" "a" "t"
    [1261] "g" "a" "a" "t" "a" "t" "t" "c" "t" "g" "g" "t" "t" "c" "t" "t" "c" "a"
    [1279] "g" "a" "g" "a" "a" "a" "a" "t" "a" "g" "a" "c" "t" "t" "a" "c" "t" "g"
    [1297] "g" "c" "c" "a" "g" "t" "g" "a" "t" "c" "c" "t" "c" "a" "t" "g" "a" "g"
    [1315] "g" "c" "t" "t" "t" "a" "a" "t" "a" "t" "g" "t" "a" "a" "a" "a" "g" "t"
    [1333] "g" "a" "a" "a" "g" "a" "g" "t" "t" "c" "a" "c" "t" "c" "c" "a" "a" "a"
    [1351] "t" "c" "a" "g" "t" "a" "g" "a" "g" "a" "g" "t" "a" "a" "t" "a" "t" "t"
    [1369] "g" "a" "a" "g" "a" "c" "a" "a" "a" "a" "t" "a" "t" "t" "t" "g" "g" "g"
    [1387] "a" "a" "a" "a" "c" "c" "t" "a" "t" "c" "g" "g" "a" "a" "g" "a" "a" "g"
    [1405] "g" "c" "a" "a" "g" "c" "c" "t" "c" "c" "c" "c" "a" "a" "c" "t" "t" "a"
    [1423] "a" "g" "c" "c" "a" "t" "g" "t" "a" "a" "c" "t" "g" "a" "a" "a" "a" "t"
    [1441] "c" "t" "a" "a" "t" "t" "a" "t" "a" "g" "g" "a" "g" "c" "a" "t" "t" "t"
    [1459] "g" "t" "t" "a" "c" "t" "g" "a" "g" "c" "c" "a" "c" "a" "g" "a" "t" "a"
    [1477] "a" "t" "a" "c" "a" "a" "g" "a" "g" "c" "g" "t" "c" "c" "c" "c" "t" "c"
    [1495] "a" "c" "a" "a" "a" "t" "a" "a" "a" "t" "t" "a" "a" "a" "g" "c" "g" "t"
    [1513] "a" "a" "a" "a" "g" "g" "a" "g" "a" "c" "c" "t" "a" "c" "a" "t" "c" "a"
    [1531] "g" "g" "c" "c" "t" "t" "c" "a" "t" "c" "c" "t" "g" "a" "g" "g" "a" "t"
    [1549] "t" "t" "t" "a" "t" "c" "a" "a" "g" "a" "a" "a" "g" "c" "a" "g" "a" "t"
    [1567] "t" "t" "g" "g" "c" "a" "g" "t" "t" "c" "a" "a" "a" "a" "g" "a" "c" "t"
    [1585] "c" "c" "t" "g" "a" "a" "a" "t" "g" "a" "t" "a" "a" "a" "t" "c" "a" "g"
    [1603] "g" "g" "a" "a" "c" "t" "a" "a" "c" "c" "a" "a" "a" "c" "g" "g" "a" "g"
    [1621] "c" "a" "g" "a" "a" "t" "g" "g" "t" "c" "a" "a" "g" "t" "g" "a" "t" "g"
    [1639] "a" "a" "t" "a" "t" "t" "a" "c" "t" "a" "a" "t" "a" "g" "t" "g" "g" "t"
    [1657] "c" "a" "t" "g" "a" "g" "a" "a" "t" "a" "a" "a" "a" "c" "a" "a" "a" "a"
    [1675] "g" "g" "t" "g" "a" "t" "t" "c" "t" "a" "t" "t" "c" "a" "g" "a" "a" "t"
    [1693] "g" "a" "g" "a" "a" "a" "a" "a" "t" "c" "c" "t" "a" "a" "c" "c" "c" "a"
    [1711] "a" "t" "a" "g" "a" "a" "t" "c" "a" "c" "t" "c" "g" "a" "a" "a" "a" "a"
    [1729] "g" "a" "a" "t" "c" "t" "g" "c" "t" "t" "t" "c" "a" "a" "a" "a" "c" "g"
    [1747] "a" "a" "a" "g" "c" "t" "g" "a" "a" "c" "c" "t" "a" "t" "a" "a" "g" "c"
    [1765] "a" "g" "c" "a" "g" "t" "a" "t" "a" "a" "g" "c" "a" "a" "t" "a" "t" "g"
    [1783] "g" "a" "a" "c" "t" "c" "g" "a" "a" "t" "t" "a" "a" "a" "t" "a" "t" "c"
    [1801] "c" "a" "c" "a" "a" "t" "t" "c" "a" "a" "a" "a" "g" "c" "a" "c" "c" "t"
    [1819] "a" "a" "a" "a" "a" "g" "a" "a" "t" "a" "g" "g" "c" "t" "g" "a" "g" "g"
    [1837] "a" "g" "g" "a" "a" "g" "t" "c" "t" "t" "c" "t" "a" "c" "c" "a" "g" "g"
    [1855] "c" "a" "t" "a" "t" "t" "c" "a" "t" "g" "c" "g" "c" "t" "t" "g" "a" "a"
    [1873] "c" "t" "a" "g" "t" "a" "g" "t" "c" "a" "g" "t" "a" "g" "a" "a" "a" "t"
    [1891] "c" "t" "a" "a" "g" "c" "c" "c" "a" "c" "c" "t" "a" "a" "t" "t" "g" "t"
    [1909] "a" "c" "t" "g" "a" "a" "t" "t" "g" "c" "a" "a" "a" "t" "t" "g" "a" "t"
    [1927] "a" "g" "t" "t" "g" "t" "t" "c" "t" "a" "g" "c" "a" "g" "t" "g" "a" "a"
    [1945] "g" "a" "g" "a" "t" "a" "a" "a" "g" "a" "a" "a" "a" "a" "a" "a" "a" "g"
    [1963] "t" "a" "c" "a" "a" "c" "c" "a" "a" "a" "t" "g" "c" "c" "a" "g" "t" "c"
    [1981] "a" "g" "g" "c" "a" "c" "a" "g" "c" "a" "g" "a" "a" "a" "c" "c" "t" "a"
    [1999] "c" "a" "a" "c" "t" "c" "a" "t" "g" "g" "a" "a" "g" "g" "t" "a" "a" "a"
    [2017] "g" "a" "a" "c" "c" "t" "g" "c" "a" "a" "c" "t" "g" "g" "a" "g" "c" "c"
    [2035] "a" "a" "g" "a" "a" "g" "a" "g" "t" "a" "a" "c" "a" "a" "g" "c" "c" "a"
    [2053] "a" "a" "t" "g" "a" "a" "c" "a" "g" "a" "c" "a" "a" "g" "t" "a" "a" "a"
    [2071] "a" "g" "a" "c" "a" "t" "g" "a" "c" "a" "g" "c" "g" "a" "t" "a" "c" "t"
    [2089] "t" "t" "c" "c" "c" "a" "g" "a" "g" "c" "t" "g" "a" "a" "g" "t" "t" "a"
    [2107] "a" "c" "a" "a" "a" "t" "g" "c" "a" "c" "c" "t" "g" "g" "t" "t" "c" "t"
    [2125] "t" "t" "t" "a" "c" "t" "a" "a" "g" "t" "g" "t" "t" "c" "a" "a" "a" "t"
    [2143] "a" "c" "c" "a" "g" "t" "g" "a" "a" "c" "t" "t" "a" "a" "a" "g" "a" "a"
    [2161] "t" "t" "t" "g" "t" "c" "a" "a" "t" "c" "c" "t" "a" "g" "c" "c" "t" "t"
    [2179] "c" "c" "a" "a" "g" "a" "g" "a" "a" "g" "a" "a" "a" "a" "a" "g" "a" "a"
    [2197] "g" "a" "g" "a" "a" "a" "c" "t" "a" "g" "a" "a" "a" "c" "a" "g" "t" "t"
    [2215] "a" "a" "a" "g" "t" "g" "t" "c" "t" "a" "a" "t" "a" "a" "t" "g" "c" "t"
    [2233] "g" "a" "a" "g" "a" "c" "c" "c" "c" "a" "a" "a" "g" "a" "t" "c" "t" "c"
    [2251] "a" "t" "g" "t" "t" "a" "a" "g" "t" "g" "g" "a" "g" "a" "a" "a" "g" "g"
    [2269] "g" "t" "t" "t" "t" "g" "c" "a" "a" "a" "c" "t" "g" "a" "a" "a" "g" "a"
    [2287] "t" "c" "t" "g" "t" "a" "g" "a" "g" "a" "g" "t" "a" "g" "c" "a" "g" "t"
    [2305] "a" "t" "t" "t" "c" "a" "t" "t" "g" "g" "t" "a" "c" "c" "t" "g" "g" "t"
    [2323] "a" "c" "t" "g" "a" "t" "t" "a" "t" "g" "g" "c" "a" "c" "t" "c" "a" "g"
    [2341] "g" "a" "a" "a" "g" "t" "a" "t" "c" "t" "c" "g" "t" "t" "a" "c" "t" "g"
    [2359] "g" "a" "a" "g" "t" "t" "a" "g" "c" "a" "c" "t" "c" "t" "a" "g" "g" "g"
    [2377] "a" "a" "g" "g" "c" "a" "a" "a" "a" "a" "c" "a" "g" "a" "a" "c" "c" "a"
    [2395] "a" "a" "t" "a" "a" "a" "t" "g" "t" "g" "t" "g" "a" "g" "t" "c" "a" "g"
    [2413] "t" "g" "t" "g" "c" "a" "g" "c" "a" "t" "t" "t" "g" "a" "a" "a" "a" "c"
    [2431] "c" "c" "c" "a" "a" "g" "g" "g" "a" "c" "t" "a" "a" "t" "t" "c" "a" "t"
    [2449] "g" "g" "t" "t" "g" "t" "t" "c" "c" "a" "a" "a" "g" "a" "t" "a" "a" "t"
    [2467] "a" "g" "a" "a" "a" "t" "g" "a" "c" "a" "c" "a" "g" "a" "a" "g" "g" "c"
    [2485] "t" "t" "t" "a" "a" "g" "t" "a" "t" "c" "c" "a" "t" "t" "g" "g" "g" "a"
    [2503] "c" "a" "t" "g" "a" "a" "g" "t" "t" "a" "a" "c" "c" "a" "c" "a" "g" "t"
    [2521] "c" "g" "g" "g" "a" "a" "a" "c" "a" "a" "g" "c" "a" "t" "a" "g" "a" "a"
    [2539] "a" "t" "g" "g" "a" "a" "g" "a" "a" "a" "g" "t" "g" "a" "a" "c" "t" "t"
    [2557] "g" "a" "t" "g" "c" "t" "c" "a" "g" "t" "a" "t" "t" "t" "g" "c" "a" "g"
    [2575] "a" "a" "t" "a" "c" "a" "t" "t" "c" "a" "a" "g" "g" "t" "t" "t" "c" "a"
    [2593] "a" "a" "g" "c" "g" "c" "c" "a" "g" "t" "c" "a" "t" "t" "t" "g" "c" "t"
    [2611] "c" "c" "g" "t" "t" "t" "t" "c" "a" "a" "a" "t" "c" "c" "a" "g" "g" "a"
    [2629] "a" "a" "t" "g" "c" "a" "g" "a" "a" "g" "a" "g" "g" "a" "a" "t" "g" "t"
    [2647] "g" "c" "a" "a" "c" "a" "t" "t" "c" "t" "c" "t" "g" "c" "c" "c" "a" "c"
    [2665] "t" "c" "t" "g" "g" "g" "t" "c" "c" "t" "t" "a" "a" "a" "g" "a" "a" "a"
    [2683] "c" "a" "a" "a" "g" "t" "c" "c" "a" "a" "a" "a" "g" "t" "c" "a" "c" "t"
    [2701] "t" "t" "t" "g" "a" "a" "t" "g" "t" "g" "a" "a" "c" "a" "a" "a" "a" "g"
    [2719] "g" "a" "a" "g" "a" "a" "a" "a" "t" "c" "a" "a" "g" "g" "a" "a" "a" "g"
    [2737] "a" "a" "t" "g" "a" "g" "t" "c" "t" "a" "a" "t" "a" "t" "c" "a" "a" "g"
    [2755] "c" "c" "t" "g" "t" "a" "c" "a" "g" "a" "c" "a" "g" "t" "t" "a" "a" "t"
    [2773] "a" "t" "c" "a" "c" "t" "g" "c" "a" "g" "g" "c" "t" "t" "t" "c" "c" "t"
    [2791] "g" "t" "g" "g" "t" "t" "g" "g" "t" "c" "a" "g" "a" "a" "a" "g" "a" "t"
    [2809] "a" "a" "g" "c" "c" "a" "g" "t" "t" "g" "a" "t" "a" "a" "t" "g" "c" "c"
    [2827] "a" "a" "a" "t" "g" "t" "a" "g" "t" "a" "t" "c" "a" "a" "a" "g" "g" "a"
    [2845] "g" "g" "c" "t" "c" "t" "a" "g" "g" "t" "t" "t" "t" "g" "t" "c" "t" "a"
    [2863] "t" "c" "a" "t" "c" "t" "c" "a" "g" "t" "t" "c" "a" "g" "a" "g" "g" "c"
    [2881] "a" "a" "c" "g" "a" "a" "a" "c" "t" "g" "g" "a" "c" "t" "c" "a" "t" "t"
    [2899] "a" "c" "t" "c" "c" "a" "a" "a" "t" "a" "a" "a" "c" "a" "t" "g" "g" "a"
    [2917] "c" "t" "t" "t" "t" "a" "c" "a" "a" "a" "a" "c" "c" "c" "a" "t" "a" "t"
    [2935] "c" "g" "t" "a" "t" "a" "c" "c" "a" "c" "c" "a" "c" "t" "t" "t" "t" "t"
    [2953] "c" "c" "c" "a" "t" "c" "a" "a" "g" "t" "c" "a" "t" "t" "t" "g" "t" "t"
    [2971] "a" "a" "a" "a" "c" "t" "a" "a" "a" "t" "g" "t" "a" "a" "g" "a" "a" "a"
    [2989] "a" "a" "t" "c" "t" "g" "c" "t" "a" "g" "a" "g" "g" "a" "a" "a" "a" "c"
    [3007] "t" "t" "t" "g" "a" "g" "g" "a" "a" "c" "a" "t" "t" "c" "a" "a" "t" "g"
    [3025] "t" "c" "a" "c" "c" "t" "g" "a" "a" "a" "g" "a" "g" "a" "a" "a" "t" "g"
    [3043] "g" "g" "a" "a" "a" "t" "g" "a" "g" "a" "a" "c" "a" "t" "t" "c" "c" "a"
    [3061] "a" "g" "t" "a" "c" "a" "g" "t" "g" "a" "g" "c" "a" "c" "a" "a" "t" "t"
    [3079] "a" "g" "c" "c" "g" "t" "a" "a" "t" "a" "a" "c" "a" "t" "t" "a" "g" "a"
    [3097] "g" "a" "a" "a" "a" "t" "g" "t" "t" "t" "t" "t" "a" "a" "a" "g" "a" "a"
    [3115] "g" "c" "c" "a" "g" "c" "t" "c" "a" "a" "g" "c" "a" "a" "t" "a" "t" "t"
    [3133] "a" "a" "t" "g" "a" "a" "g" "t" "a" "g" "g" "t" "t" "c" "c" "a" "g" "t"
    [3151] "a" "c" "t" "a" "a" "t" "g" "a" "a" "g" "t" "g" "g" "g" "c" "t" "c" "c"
    [3169] "a" "g" "t" "a" "t" "t" "a" "a" "t" "g" "a" "a" "a" "t" "a" "g" "g" "t"
    [3187] "t" "c" "c" "a" "g" "t" "g" "a" "t" "g" "a" "a" "a" "a" "c" "a" "t" "t"
    [3205] "c" "a" "a" "g" "c" "a" "g" "a" "a" "c" "t" "a" "g" "g" "t" "a" "g" "a"
    [3223] "a" "a" "c" "a" "g" "a" "g" "g" "g" "c" "c" "a" "a" "a" "a" "t" "t" "g"
    [3241] "a" "a" "t" "g" "c" "t" "a" "t" "g" "c" "t" "t" "a" "g" "a" "t" "t" "a"
    [3259] "g" "g" "g" "g" "t" "t" "t" "t" "g" "c" "a" "a" "c" "c" "t" "g" "a" "g"
    [3277] "g" "t" "c" "t" "a" "t" "a" "a" "a" "c" "a" "a" "a" "g" "t" "c" "t" "t"
    [3295] "c" "c" "t" "g" "g" "a" "a" "g" "t" "a" "a" "t" "t" "g" "t" "a" "a" "g"
    [3313] "c" "a" "t" "c" "c" "t" "g" "a" "a" "a" "t" "a" "a" "a" "a" "a" "a" "g"
    [3331] "c" "a" "a" "g" "a" "a" "t" "a" "t" "g" "a" "a" "g" "a" "a" "g" "t" "a"
    [3349] "g" "t" "t" "c" "a" "g" "a" "c" "t" "g" "t" "t" "a" "a" "t" "a" "c" "a"
    [3367] "g" "a" "t" "t" "t" "c" "t" "c" "t" "c" "c" "a" "t" "a" "t" "c" "t" "g"
    [3385] "a" "t" "t" "t" "c" "a" "g" "a" "t" "a" "a" "c" "t" "t" "a" "g" "a" "a"
    [3403] "c" "a" "g" "c" "c" "t" "a" "t" "g" "g" "g" "a" "a" "g" "t" "a" "g" "t"
    [3421] "c" "a" "t" "g" "c" "a" "t" "c" "t" "c" "a" "g" "g" "t" "t" "t" "g" "t"
    [3439] "t" "c" "t" "g" "a" "g" "a" "c" "a" "c" "c" "t" "g" "a" "t" "g" "a" "c"
    [3457] "c" "t" "g" "t" "t" "a" "g" "a" "t" "g" "a" "t" "g" "g" "t" "g" "a" "a"
    [3475] "a" "t" "a" "a" "a" "g" "g" "a" "a" "g" "a" "t" "a" "c" "t" "a" "g" "t"
    [3493] "t" "t" "t" "g" "c" "t" "g" "a" "a" "a" "a" "t" "g" "a" "c" "a" "t" "t"
    [3511] "a" "a" "g" "g" "a" "a" "a" "g" "t" "t" "c" "t" "g" "c" "t" "g" "t" "t"
    [3529] "t" "t" "t" "a" "g" "c" "a" "a" "a" "a" "g" "c" "g" "t" "c" "c" "a" "g"
    [3547] "a" "a" "a" "g" "g" "a" "g" "a" "g" "c" "t" "t" "a" "g" "c" "a" "g" "g"
    [3565] "a" "g" "t" "c" "c" "t" "a" "g" "c" "c" "c" "t" "t" "t" "c" "a" "c" "c"
    [3583] "c" "a" "t" "a" "c" "a" "c" "a" "t" "t" "t" "g" "g" "c" "t" "c" "a" "g"
    [3601] "g" "g" "t" "t" "a" "c" "c" "g" "a" "a" "g" "a" "g" "g" "g" "g" "c" "c"
    [3619] "a" "a" "g" "a" "a" "a" "t" "t" "a" "g" "a" "g" "t" "c" "c" "t" "c" "a"
    [3637] "g" "a" "a" "g" "a" "g" "a" "a" "c" "t" "t" "a" "t" "c" "t" "a" "g" "t"
    [3655] "g" "a" "g" "g" "a" "t" "g" "a" "a" "g" "a" "g" "c" "t" "t" "c" "c" "c"
    [3673] "t" "g" "c" "t" "t" "c" "c" "a" "a" "c" "a" "c" "t" "t" "g" "t" "t" "a"
    [3691] "t" "t" "t" "g" "g" "t" "a" "a" "a" "g" "t" "a" "a" "a" "c" "a" "a" "t"
    [3709] "a" "t" "a" "c" "c" "t" "t" "c" "t" "c" "a" "g" "t" "c" "t" "a" "c" "t"
    [3727] "a" "g" "g" "c" "a" "t" "a" "g" "c" "a" "c" "c" "g" "t" "t" "g" "c" "t"
    [3745] "a" "c" "c" "g" "a" "g" "t" "g" "t" "c" "t" "g" "t" "c" "t" "a" "a" "g"
    [3763] "a" "a" "c" "a" "c" "a" "g" "a" "g" "g" "a" "g" "a" "a" "t" "t" "t" "a"
    [3781] "t" "t" "a" "t" "c" "a" "t" "t" "g" "a" "a" "g" "a" "a" "t" "a" "g" "c"
    [3799] "t" "t" "a" "a" "a" "t" "g" "a" "c" "t" "g" "c" "a" "g" "t" "a" "a" "c"
    [3817] "c" "a" "g" "g" "t" "a" "a" "t" "a" "t" "t" "g" "g" "c" "a" "a" "a" "g"
    [3835] "g" "c" "a" "t" "c" "t" "c" "a" "g" "g" "a" "a" "c" "a" "t" "c" "a" "c"
    [3853] "c" "t" "t" "a" "g" "t" "g" "a" "g" "g" "a" "a" "a" "c" "a" "a" "a" "a"
    [3871] "t" "g" "t" "t" "c" "t" "g" "c" "t" "a" "g" "c" "t" "t" "g" "t" "t" "t"
    [3889] "t" "c" "t" "t" "c" "a" "c" "a" "g" "t" "g" "c" "a" "g" "t" "g" "a" "a"
    [3907] "t" "t" "g" "g" "a" "a" "g" "a" "c" "t" "t" "g" "a" "c" "t" "g" "c" "a"
    [3925] "a" "a" "t" "a" "c" "a" "a" "a" "c" "a" "c" "c" "c" "a" "g" "g" "a" "t"
    [3943] "c" "c" "t" "t" "t" "c" "t" "t" "g" "a" "t" "t" "g" "g" "t" "t" "c" "t"
    [3961] "t" "c" "c" "a" "a" "a" "c" "a" "a" "a" "t" "g" "a" "g" "g" "c" "a" "t"
    [3979] "c" "a" "g" "t" "c" "t" "g" "a" "a" "a" "g" "c" "c" "a" "g" "g" "g" "a"
    [3997] "g" "t" "t" "g" "g" "t" "c" "t" "g" "a" "g" "t" "g" "a" "c" "a" "a" "g"
    [4015] "g" "a" "a" "t" "t" "g" "g" "t" "t" "t" "c" "a" "g" "a" "t" "g" "a" "t"
    [4033] "g" "a" "a" "g" "a" "a" "a" "g" "a" "g" "g" "a" "a" "c" "g" "g" "g" "c"
    [4051] "t" "t" "g" "g" "a" "a" "g" "a" "a" "a" "a" "t" "a" "a" "t" "c" "a" "a"
    [4069] "g" "a" "a" "g" "a" "g" "c" "a" "a" "a" "g" "c" "a" "t" "g" "g" "a" "t"
    [4087] "t" "c" "a" "a" "a" "c" "t" "t" "a" "g" "g" "t" "g" "a" "a" "g" "c" "a"
    [4105] "g" "c" "a" "t" "c" "t" "g" "g" "g" "t" "g" "t" "g" "a" "g" "a" "g" "t"
    [4123] "g" "a" "a" "a" "c" "a" "a" "g" "c" "g" "t" "c" "t" "c" "t" "g" "a" "a"
    [4141] "g" "a" "c" "t" "g" "c" "t" "c" "a" "g" "g" "g" "c" "t" "a" "t" "c" "c"
    [4159] "t" "c" "t" "c" "a" "g" "a" "g" "t" "g" "a" "c" "a" "t" "t" "t" "t" "a"
    [4177] "a" "c" "c" "a" "c" "t" "c" "a" "g" "c" "a" "g" "a" "g" "g" "g" "a" "t"
    [4195] "a" "c" "c" "a" "t" "g" "c" "a" "a" "c" "a" "t" "a" "a" "c" "c" "t" "g"
    [4213] "a" "t" "a" "a" "a" "g" "c" "t" "c" "c" "a" "g" "c" "a" "g" "g" "a" "a"
    [4231] "a" "t" "g" "g" "c" "t" "g" "a" "a" "c" "t" "a" "g" "a" "a" "g" "c" "t"
    [4249] "g" "t" "g" "t" "t" "a" "g" "a" "a" "c" "a" "g" "c" "a" "t" "g" "g" "g"
    [4267] "a" "g" "c" "c" "a" "g" "c" "c" "t" "t" "c" "t" "a" "a" "c" "a" "g" "c"
    [4285] "t" "a" "c" "c" "c" "t" "t" "c" "c" "a" "t" "c" "a" "t" "a" "a" "g" "t"
    [4303] "g" "a" "c" "t" "c" "t" "t" "c" "t" "g" "c" "c" "c" "t" "t" "g" "a" "g"
    [4321] "g" "a" "c" "c" "t" "g" "c" "g" "a" "a" "a" "t" "c" "c" "a" "g" "a" "a"
    [4339] "c" "a" "a" "a" "g" "c" "a" "c" "a" "t" "c" "a" "g" "a" "a" "a" "a" "a"
    [4357] "g" "c" "a" "g" "t" "a" "t" "t" "a" "a" "c" "t" "t" "c" "a" "c" "a" "g"
    [4375] "a" "a" "a" "a" "g" "t" "a" "g" "t" "g" "a" "a" "t" "a" "c" "c" "c" "t"
    [4393] "a" "t" "a" "a" "g" "c" "c" "a" "g" "a" "a" "t" "c" "c" "a" "g" "a" "a"
    [4411] "g" "g" "c" "c" "t" "t" "t" "c" "t" "g" "c" "t" "g" "a" "c" "a" "a" "g"
    [4429] "t" "t" "t" "g" "a" "g" "g" "t" "g" "t" "c" "t" "g" "c" "a" "g" "a" "t"
    [4447] "a" "g" "t" "t" "c" "t" "a" "c" "c" "a" "g" "t" "a" "a" "a" "a" "a" "t"
    [4465] "a" "a" "a" "g" "a" "a" "c" "c" "a" "g" "g" "a" "g" "t" "g" "g" "a" "a"
    [4483] "a" "g" "g" "t" "c" "a" "t" "c" "c" "c" "c" "t" "t" "c" "t" "a" "a" "a"
    [4501] "t" "g" "c" "c" "c" "a" "t" "c" "a" "t" "t" "a" "g" "a" "t" "g" "a" "t"
    [4519] "a" "g" "g" "t" "g" "g" "t" "a" "c" "a" "t" "g" "c" "a" "c" "a" "g" "t"
    [4537] "t" "g" "c" "t" "c" "t" "g" "g" "g" "a" "g" "t" "c" "t" "t" "c" "a" "g"
    [4555] "a" "a" "t" "a" "g" "a" "a" "a" "c" "t" "a" "c" "c" "c" "a" "t" "c" "t"
    [4573] "c" "a" "a" "g" "a" "g" "g" "a" "g" "c" "t" "c" "a" "t" "t" "a" "a" "g"
    [4591] "g" "t" "t" "g" "t" "t" "g" "a" "t" "g" "t" "g" "g" "a" "g" "g" "a" "g"
    [4609] "c" "a" "a" "c" "a" "g" "c" "t" "g" "g" "a" "a" "g" "a" "g" "t" "c" "t"
    [4627] "g" "g" "g" "c" "c" "a" "c" "a" "c" "g" "a" "t" "t" "t" "g" "a" "c" "g"
    [4645] "g" "a" "a" "a" "c" "a" "t" "c" "t" "t" "a" "c" "t" "t" "g" "c" "c" "a"
    [4663] "a" "g" "g" "c" "a" "a" "g" "a" "t" "c" "t" "a" "g" "a" "g" "g" "g" "a"
    [4681] "a" "c" "c" "c" "c" "t" "t" "a" "c" "c" "t" "g" "g" "a" "a" "t" "c" "t"
    [4699] "g" "g" "a" "a" "t" "c" "a" "g" "c" "c" "t" "c" "t" "t" "c" "t" "c" "t"
    [4717] "g" "a" "t" "g" "a" "c" "c" "c" "t" "g" "a" "a" "t" "c" "t" "g" "a" "t"
    [4735] "c" "c" "t" "t" "c" "t" "g" "a" "a" "g" "a" "c" "a" "g" "a" "g" "c" "c"
    [4753] "c" "c" "a" "g" "a" "g" "t" "c" "a" "g" "c" "t" "c" "g" "t" "g" "t" "t"
    [4771] "g" "g" "c" "a" "a" "c" "a" "t" "a" "c" "c" "a" "t" "c" "t" "t" "c" "a"
    [4789] "a" "c" "c" "t" "c" "t" "g" "c" "a" "t" "t" "g" "a" "a" "a" "g" "t" "t"
    [4807] "c" "c" "c" "c" "a" "a" "t" "t" "g" "a" "a" "a" "g" "t" "t" "g" "c" "a"
    [4825] "g" "a" "a" "t" "c" "t" "g" "c" "c" "c" "a" "g" "a" "g" "t" "c" "c" "a"
    [4843] "g" "c" "t" "g" "c" "t" "g" "c" "t" "c" "a" "t" "a" "c" "t" "a" "c" "t"
    [4861] "g" "a" "t" "a" "c" "t" "g" "c" "t" "g" "g" "g" "t" "a" "t" "a" "a" "t"
    [4879] "g" "c" "a" "a" "t" "g" "g" "a" "a" "g" "a" "a" "a" "g" "t" "g" "t" "g"
    [4897] "a" "g" "c" "a" "g" "g" "g" "a" "g" "a" "a" "g" "c" "c" "a" "g" "a" "a"
    [4915] "t" "t" "g" "a" "c" "a" "g" "c" "t" "t" "c" "a" "a" "c" "a" "g" "a" "a"
    [4933] "a" "g" "g" "g" "t" "c" "a" "a" "c" "a" "a" "a" "a" "g" "a" "a" "t" "g"
    [4951] "t" "c" "c" "a" "t" "g" "g" "t" "g" "g" "t" "g" "t" "c" "t" "g" "g" "c"
    [4969] "c" "t" "g" "a" "c" "c" "c" "c" "a" "g" "a" "a" "g" "a" "a" "t" "t" "t"
    [4987] "a" "t" "g" "c" "t" "c" "g" "t" "g" "t" "a" "c" "a" "a" "g" "t" "t" "t"
    [5005] "g" "c" "c" "a" "g" "a" "a" "a" "a" "c" "a" "c" "c" "a" "c" "a" "t" "c"
    [5023] "a" "c" "t" "t" "t" "a" "a" "c" "t" "a" "a" "t" "c" "t" "a" "a" "t" "t"
    [5041] "a" "c" "t" "g" "a" "a" "g" "a" "g" "a" "c" "t" "a" "c" "t" "c" "a" "t"
    [5059] "g" "t" "t" "g" "t" "t" "a" "t" "g" "a" "a" "a" "a" "c" "a" "g" "a" "t"
    [5077] "g" "c" "t" "g" "a" "g" "t" "t" "t" "g" "t" "g" "t" "g" "t" "g" "a" "a"
    [5095] "c" "g" "g" "a" "c" "a" "c" "t" "g" "a" "a" "a" "t" "a" "t" "t" "t" "t"
    [5113] "c" "t" "a" "g" "g" "a" "a" "t" "t" "g" "c" "g" "g" "g" "a" "g" "g" "a"
    [5131] "a" "a" "a" "t" "g" "g" "g" "t" "a" "g" "t" "t" "a" "g" "c" "t" "a" "t"
    [5149] "t" "t" "c" "t" "g" "g" "g" "t" "g" "a" "c" "c" "c" "a" "g" "t" "c" "t"
    [5167] "a" "t" "t" "a" "a" "a" "g" "a" "a" "a" "g" "a" "a" "a" "a" "a" "t" "g"
    [5185] "c" "t" "g" "a" "a" "t" "g" "a" "g" "c" "a" "t" "g" "a" "t" "t" "t" "t"
    [5203] "g" "a" "a" "g" "t" "c" "a" "g" "a" "g" "g" "a" "g" "a" "t" "g" "t" "g"
    [5221] "g" "t" "c" "a" "a" "t" "g" "g" "a" "a" "g" "a" "a" "a" "c" "c" "a" "c"
    [5239] "c" "a" "a" "g" "g" "t" "c" "c" "a" "a" "a" "g" "c" "g" "a" "g" "c" "a"
    [5257] "a" "g" "a" "g" "a" "a" "t" "c" "c" "c" "a" "g" "g" "a" "c" "a" "g" "a"
    [5275] "a" "a" "g" "a" "t" "c" "t" "t" "c" "a" "g" "g" "g" "g" "g" "c" "t" "a"
    [5293] "g" "a" "a" "a" "t" "c" "t" "g" "t" "t" "g" "c" "t" "a" "t" "g" "g" "g"
    [5311] "c" "c" "c" "t" "t" "c" "a" "c" "c" "a" "a" "c" "a" "t" "g" "c" "c" "c"
    [5329] "a" "c" "a" "g" "a" "t" "c" "a" "a" "c" "t" "g" "g" "a" "a" "t" "g" "g"
    [5347] "a" "t" "g" "g" "t" "a" "c" "a" "g" "c" "t" "g" "t" "g" "t" "g" "g" "t"
    [5365] "g" "c" "t" "t" "c" "t" "g" "t" "g" "g" "t" "g" "a" "a" "g" "g" "a" "g"
    [5383] "c" "t" "t" "t" "c" "a" "t" "c" "a" "t" "t" "c" "a" "c" "c" "c" "t" "t"
    [5401] "g" "g" "c" "a" "c" "a" "g" "g" "t" "g" "t" "c" "c" "a" "c" "c" "c" "a"
    [5419] "a" "t" "t" "g" "t" "g" "g" "t" "t" "g" "t" "g" "c" "a" "g" "c" "c" "a"
    [5437] "g" "a" "t" "g" "c" "c" "t" "g" "g" "a" "c" "a" "g" "a" "g" "g" "a" "c"
    [5455] "a" "a" "t" "g" "g" "c" "t" "t" "c" "c" "a" "t" "g" "c" "a" "a" "t" "t"
    [5473] "g" "g" "g" "c" "a" "g" "a" "t" "g" "t" "g" "t" "g" "a" "g" "g" "c" "a"
    [5491] "c" "c" "t" "g" "t" "g" "g" "t" "g" "a" "c" "c" "c" "g" "a" "g" "a" "g"
    [5509] "t" "g" "g" "g" "t" "g" "t" "t" "g" "g" "a" "c" "a" "g" "t" "g" "t" "a"
    [5527] "g" "c" "a" "c" "t" "c" "t" "a" "c" "c" "a" "g" "t" "g" "c" "c" "a" "g"
    [5545] "g" "a" "g" "c" "t" "g" "g" "a" "c" "a" "c" "c" "t" "a" "c" "c" "t" "g"
    [5563] "a" "t" "a" "c" "c" "c" "c" "a" "g" "a" "t" "c" "c" "c" "c" "c" "a" "c"
    [5581] "a" "g" "c" "c" "a" "c" "t" "a" "c" "t" "g" "a"
    attr(,"name")
    [1] "L78833.BRCA1"
    attr(,"Annot")
    [1] ">L78833.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $MG969402
      [1] "t" "g" "g" "a" "a" "g" "t" "a" "a" "g" "g" "a" "a" "a" "c" "a" "t" "g"
     [19] "t" "a" "a" "t" "g" "a" "t" "a" "g" "g" "c" "g" "g" "a" "c" "t" "c" "c"
     [37] "c" "a" "g" "c" "a" "c" "a" "g" "a" "a" "a" "a" "a" "a" "a" "g" "g" "t"
     [55] "a" "g" "a" "t" "c" "t" "g" "a" "a" "t" "g" "c" "t" "g" "a" "t" "c" "c"
     [73] "c" "c" "t" "g" "t" "g" "t" "g" "a" "g" "a" "g" "a" "a" "a" "a" "g" "a"
     [91] "a" "t" "g" "g" "a" "a" "t" "a" "a" "g" "c" "a" "g" "a" "a" "a" "c" "t"
    [109] "g" "c" "c" "a" "t" "g" "c" "t" "c" "a" "g" "a" "g" "a" "a" "t" "c" "c"
    [127] "t" "a" "g" "a" "g" "a" "t" "a" "c" "t" "g" "a" "a" "g" "a" "t" "g" "t"
    [145] "t" "c" "c" "t" "t" "g" "g" "a" "t" "a" "a" "c" "a" "c" "t" "a" "a" "a"
    [163] "t" "a" "g" "c" "a" "g" "c" "a" "t" "t" "c" "a" "g" "a" "a" "a" "g" "t"
    [181] "t" "a" "a" "t" "g" "a" "g" "t" "g" "g" "t" "t" "t" "t" "c" "c" "a" "g"
    [199] "a" "a" "g" "t" "g" "a" "t" "g" "a" "a" "c" "t" "g" "t" "t" "a" "g" "g"
    [217] "t" "t" "c" "t" "g" "a" "t" "g" "a" "c" "t" "c" "a" "c" "a" "t" "g" "a"
    [235] "t" "g" "g" "g" "g" "a" "g" "t" "c" "t" "g" "a" "a" "t" "c" "a" "a" "a"
    [253] "t" "g" "c" "c" "a" "a" "a" "g" "t" "a" "g" "c" "t" "g" "a" "t" "g" "t"
    [271] "a" "t"
    attr(,"name")
    [1] "MG969402"
    attr(,"Annot")
    [1] ">MG969402"
    attr(,"class")
    [1] "SeqFastadna"
    
    $S78558.BRCA1
     [1] "c" "c" "a" "g" "a" "t" "c" "c" "t" "g" "g" "a" "c" "a" "g" "a" "g" "g" "a"
    [20] "c" "a" "a" "t" "g" "g" "c" "t" "t" "c" "c" "a" "t" "g" "c" "a" "a" "t" "t"
    [39] "g" "g" "g" "c" "a" "g" "a" "t" "g" "t" "g" "t" "g" "a" "g" "g" "c" "a" "c"
    [58] "c" "t" "g" "t" "g" "g" "t" "g" "a"
    attr(,"name")
    [1] "S78558.BRCA1"
    attr(,"Annot")
    [1] ">S78558.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    
    $Y08757.BRCA1
     [1] "a" "t" "g" "c" "t" "c" "g" "t" "g" "t" "a" "c" "a" "a" "g" "t" "t" "t" "g"
    [20] "c" "c" "a" "g" "a" "a" "a" "a" "c" "a" "c" "c" "a" "c" "a" "t" "c" "a" "c"
    [39] "t" "t" "t" "a" "a" "c" "t" "a" "a" "t" "c" "t" "a" "a" "t" "t" "a" "c" "t"
    [58] "g" "a" "a" "g" "a" "g" "a" "c" "t" "a" "c" "t" "c" "a" "t" "g" "t" "t" "g"
    [77] "t" "t" "a" "t" "g" "a" "a" "a" "a" "c" "a" "g"
    attr(,"name")
    [1] "Y08757.BRCA1"
    attr(,"Annot")
    [1] ">Y08757.BRCA1"
    attr(,"class")
    [1] "SeqFastadna"
    



```R
closebank()
```


```R

```
