
# Phylogenetic analysis and tree plotting
Phylogenetic analysis is about finding the evolutionary relationship among species(organisms) based on sequence data. Once we have a set of sequences from different sources, we can use the plots to understand how close or distant they are in terms of molecular evolution. As a result of the mutations during evolution, there emerge differences at the sequence level. These differences can be represented in terms of distance measures (check the See also section at the end of this recipe). These measures can then be used to estimate the evolutionary relations among the species, often represented as phylogenetic trees. This relation is very often depicted in terms of a phylogenetic tree. Now we can look into the various aspects of sequence retrieval, alignment, and analysis. <br>
`ape` package is important here.


```R
install.packages("ape", repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
```

    Installing package into 'C:/Users/Master1/Documents/R/win-library/3.4'
    (as 'lib' is unspecified)
    also installing the dependencies 'nlme', 'lattice', 'Rcpp'
    
    

    package 'nlme' successfully unpacked and MD5 sums checked
    package 'lattice' successfully unpacked and MD5 sums checked
    package 'Rcpp' successfully unpacked and MD5 sums checked
    package 'ape' successfully unpacked and MD5 sums checked
    
    The downloaded binary packages are in
    	C:\Users\Master1\AppData\Local\Temp\Rtmpa2HoEf\downloaded_packages
    


```R
library(ape)
```


```R
# Define the sequences of our interest as a set in terms of IDs:
myset <- c("U15717", "U15718", "U15719", "U15720", "U15721", "U15722")
myset
```


<ol class=list-inline>
	<li>'U15717'</li>
	<li>'U15718'</li>
	<li>'U15719'</li>
	<li>'U15720'</li>
	<li>'U15721'</li>
	<li>'U15722'</li>
</ol>




```R
# Fetch the sequences that we have defined, and check 'read.GenBnak'(connects to the GenBank database, and reads nucleotide sequences using accession numbers given as arguments):
myseqs <- read.GenBank(myset)
myseqs
?read.GenBank
```


    6 DNA sequences in binary format stored in a list.
    
    All sequences of same length: 1045 
    
    Labels:
    U15717
    U15718
    U15719
    U15720
    U15721
    U15722
    
    Base composition:
        a     c     g     t 
    0.267 0.351 0.134 0.248 



```R
# Compute the distance matrix for the sequences with the 'dist.dna' function, just check 'dist.dna'
mydist <- dist.dna(myseqs)
mydist
?dist.dna
```


                 U15717       U15718       U15719       U15720       U15721
    U15718 0.0963968720                                                    
    U15719 0.0519601191 0.0821667377                                       
    U15720 0.0155113115 0.0932071967 0.0489261666                          
    U15721 0.0624466459 0.0797534702 0.0498979851 0.0551246867             
    U15722 0.0164965334 0.0920721224 0.0478790458 0.0009578547 0.0540676038



```R
# Use the following 'triangMtd' function to reconstruct the tree to get the 'phylo' object for the phylogenetic trees:
myphylo <- triangMtd(mydist)
myphylo
```


    
    Phylogenetic tree with 6 tips and 4 internal nodes.
    
    Tip labels:
    [1] "U15717" "U15718" "U15719" "U15720" "U15721" "U15722"
    
    Unrooted; includes branch lengths.



```R
# there are 4 types of phylogenetic trees for nanlysis
plot(myphylo, type="phylogram", edge.color="red", cex=1, edge.width=1,main="(A) Phylogram")
```


[plot](https://github.com/Chengshu21/Bioinformatics-with-R-cookbook/blob/master/Chapter-3-Sequence-Analysis-with-R-master/md/output_7-Phylogenetic%20analysis%20and%20tree%20plotting.png)



```R
plot(myphylo, type="cladogram", edge.color="red", cex=1, edge.width=1, main="(B) Cladogram")
```


[plot](https://github.com/Chengshu21/Bioinformatics-with-R-cookbook/blob/master/Chapter-3-Sequence-Analysis-with-R-master/md/output_8-Phylogenetic%20analysis%20and%20tree%20plotting.png)



```R
plot(myphylo, type="fan", edge.color="red", cex=1, edge.width=1, main="(C) Fan")
```


[plot](https://github.com/Chengshu21/Bioinformatics-with-R-cookbook/blob/master/Chapter-3-Sequence-Analysis-with-R-master/md/output_9-Phylogenetic%20analysis%20and%20tree%20plotting.png)



```R
plot(myphylo, type="unrooted", edge.color="red", cex=1, edge.width=1, main="(D) Unrooted")
```


[plot](https://github.com/Chengshu21/Bioinformatics-with-R-cookbook/blob/master/Chapter-3-Sequence-Analysis-with-R-master/md/output_10-Phylogenetic%20analysis%20and%20tree%20plotting.png)



```R
?boot.phylo
```


```R

```
