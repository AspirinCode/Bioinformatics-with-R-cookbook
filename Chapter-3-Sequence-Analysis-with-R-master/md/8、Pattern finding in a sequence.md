
# Pattern finding in a sequence
We know about fetching a sequence, aligning, and matching two or more sequences. And checking for some patterns in the gene apart from the general nucleotide content analysis is also important. This can help us to detect various interesting subsets of characters in nucleotide or protein sequences, for example, the start and stop codons in nucleic acids. Start and stop codons mark the site at which the translation into protein sequences begins and the site at which the translation ends. The common example of a start codon in mRNAs is AUG and for a stop codon, some examples are UAG,
UGA, and UAA. Even in the case of protein sequences, we can search conserved motifs. We will try in finding patterns in the sequence that can be interesting from the point of view of functional and molecular biology, such as gene finding in a sequence or ORF finding in a sequence. ORFs are the frames in a genomic sequence that are not interrupted by
stop codons.<br>
`Biostrings` library is necessary in thia part.


```R
library(Biostrings)
```


```R
# create the sequence to be analyzed(or we can use the sequences fetched from GenBank as well):
mynucleotide <- DNAString("aacataatgcagtagaacccatgagccc")
mynucleotide
```


      28-letter "DNAString" instance
    seq: AACATAATGCAGTAGAACCCATGAGCCC



```R
# Look for a pattern of the sequence, such as a start codon ATG:
matchPattern(DNAString("ATG"), mynucleotide)
```


      Views on a 28-letter DNAString subject
    subject: AACATAATGCAGTAGAACCCATGAGCCC
    views:
        start end width
    [1]     7   9     3 [ATG]
    [2]    21  23     3 [ATG]



```R
# Similarly, look for the pattern for the stop codons, such as TAA or other stop codons:
matchPattern("TAA", mynucleotide)
```


      Views on a 28-letter DNAString subject
    subject: AACATAATGCAGTAGAACCCATGAGCCC
    views:
        start end width
    [1]     5   7     3 [TAA]



```R
# combine these two aspects into a single function to return the overall results for all the codons:
myCodonFinder <- function(sequence){
    startCodon = DNAString("ATG") # Assign start codons 
    stopCodons = list("TAA", "TAG", "TGA") # Assign stop codons
    codonPosition = list() #initialize the output to be returned as a list
    codonPosition$Start = matchPattern(startCodon, sequence) # search start codons
    x=list()
    for(i in 1:3){ # iterate over all stop codons
        x[[i]]= matchPattern(DNAString (stopCodons[[i]]), sequence)
        codonPosition$Stop=x
    }
    return(codonPosition) # returns results
}
```


```R
# paste the code for the previous function into the R session and run it with your sequence object mynucleotide :
StartStops <- myCodonFinder(mynucleotide)
StartStops
# Alternatively, save the source code as a file, `myCodonFinder.R`, and source it into the R session:`source("myCodonFinder.R")`
```


    $Start
      Views on a 28-letter DNAString subject
    subject: AACATAATGCAGTAGAACCCATGAGCCC
    views:
        start end width
    [1]     7   9     3 [ATG]
    [2]    21  23     3 [ATG]
    
    $Stop
    $Stop[[1]]
      Views on a 28-letter DNAString subject
    subject: AACATAATGCAGTAGAACCCATGAGCCC
    views:
        start end width
    [1]     5   7     3 [TAA]
    
    $Stop[[2]]
      Views on a 28-letter DNAString subject
    subject: AACATAATGCAGTAGAACCCATGAGCCC
    views:
        start end width
    [1]    13  15     3 [TAG]
    
    $Stop[[3]]
      Views on a 28-letter DNAString subject
    subject: AACATAATGCAGTAGAACCCATGAGCCC
    views:
        start end width
    [1]    22  24     3 [TGA]
    
    


We can use this to find *open reading frames (ORFs)*.<br>
The `matchPattern` function of `Biostrings` is an implementation to identify the occurrences of a particular pattern or motif in a sequence. It requires a string as an input(not a vector of characters) that is created by the `DNAString` function. If you intend to use a vector of characters from the *Retrieving a sequence , you must convert it into a
string via the use of the `c2s` and `DNAString` functions (`AAString` for a protein sequence and `RNAString` for RNAs). The `matchPattern` function returns a table with columns that represent the start, end, and width of the hit or match (the width is obviously the same as the length of the pattern). All the hits found in the sequence are arranged in the rows of the returned object.<br>
In our function, to find all the start and stop codons, we combine two pattern searches; first, we combine the search for the start codon, and second, we look iteratively for all the three stop codons. We can modify this function to find an ORF. This is left for the reader to try. The `Biostrings` package has many other interesting functions to solve more complicated issues regarding patterns in sequences. Another similar example is the `matchLRPatterns` function that finds the paired matches in a sequence. <br>
These two packages (`Biostrings` and `seqinr`) are important to learn more about such methods.


```R
mytarget <- DNAString("AAATTAACCCTT")
matchLRPatterns("AA", "TT", 5, mytarget)
```


      Views on a 12-letter DNAString subject
    subject: AAATTAACCCTT
    views:
        start end width
    [1]     1   5     5 [AAATT]
    [2]     2   5     4 [AATT]
    [3]     6  12     7 [AACCCTT]



```R

```
