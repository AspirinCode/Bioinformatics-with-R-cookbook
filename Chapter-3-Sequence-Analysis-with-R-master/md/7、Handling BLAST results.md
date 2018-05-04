
# Handling BLAST results
Basic Local Alignment Search Tool (BLAST) is a basic tool to look for local similarities across sequences against a database of relevant sequences. It takes an input sequence (nucleotide or protein) and a reference database (genomic data or EST) and performs the local alignment of the former against the latter. We can understand the function and comparison of new sequences if we undertsand the relation among them. We can BLAST any of the sequences that we retrieved in the preceding sections on the [NCBI BLAST](http://blast.ncbi.nlm.nih.gov/Blast.cgi) web page.
BLAST results are usually an HTML page but can be downloaded as a table. <br>
At first, we should have a BLAST result. Input the sequence we retrieved(The input can either be made in the textbox by pasting the sequence or path to a `FASTA` file or a sequence ID). Then, we can run the `BLAST` program of our retrieved sequence to get the results. The result table can then be downloaded as a `table`. <br>
If the results are downloaded from the NCBI standalone BLAST, we can easily use the `RFLPtools` package to do DNA fragment analysis and BLAST report analysis.


```R
source("http://www.bioconductor.org/biocLite.R")
```

    Bioconductor version 3.6 (BiocInstaller 1.28.0), ?biocLite for help
    


```R
biocLite("RFLPtools"）
```

    BioC_mirror: https://bioconductor.org
    Using Bioconductor 3.6 (BiocInstaller 1.28.0), R 3.4.2 (2017-09-28).
    Installing package(s) 'RFLPtools'
    

    package 'RFLPtools' successfully unpacked and MD5 sums checked
    
    The downloaded binary packages are in
    	C:\Users\Master1\AppData\Local\Temp\Rtmpwx21tx\downloaded_packages
    

    Old packages: 'GenomicFeatures', 'GenomicRanges', 'matrixStats', 'pbdZMQ',
      'tibble', 'tidyr', 'digest', 'stringi'
    


```R
library(RFLPtools)
```


```R
# Use the built-in blast result from the package for a demo :
DIR <- system.file("extdata", package = "RFLPtools")
MyFile <- file.path(DIR, "BLASTexample.txt")
MyFile
```


'C:/Users/Master1/Documents/R/win-library/3.4/RFLPtools/extdata/BLASTexample.txt'



```R
# read the file,use 'read.blast':
MyBlast<- read.blast(file = MyFile)
head(MyBlast)
```


<table>
<thead><tr><th scope=col>query.id</th><th scope=col>subject.id</th><th scope=col>identity</th><th scope=col>alignment.length</th><th scope=col>mismatches</th><th scope=col>gap.opens</th><th scope=col>q.start</th><th scope=col>q.end</th><th scope=col>s.start</th><th scope=col>s.end</th><th scope=col>evalue</th><th scope=col>bit.score</th></tr></thead>
<tbody>
	<tr><td>agrFF002</td><td>agrFF002</td><td>100.00  </td><td>544     </td><td> 0      </td><td>0       </td><td>  1     </td><td>544     </td><td>  1     </td><td>544     </td><td> 0.0e+00</td><td>944.0   </td></tr>
	<tr><td>agrFF002</td><td>agrFF148</td><td> 93.42  </td><td>243     </td><td>14      </td><td>2       </td><td>199     </td><td>439     </td><td>671     </td><td>913     </td><td>6.0e-102</td><td>360.0   </td></tr>
	<tr><td>agrFF002</td><td>agrFF148</td><td>100.00  </td><td> 11     </td><td> 0      </td><td>0       </td><td>462     </td><td>472     </td><td>785     </td><td>795     </td><td> 6.7e+00</td><td> 21.1   </td></tr>
	<tr><td>agrFF002</td><td>agrFF176</td><td> 91.37  </td><td>255     </td><td>20      </td><td>2       </td><td>187     </td><td>439     </td><td>123     </td><td>377     </td><td>2.0e-100</td><td>354.0   </td></tr>
	<tr><td>agrFF002</td><td>agrFF176</td><td>100.00  </td><td> 11     </td><td> 0      </td><td>0       </td><td>462     </td><td>472     </td><td>250     </td><td>260     </td><td> 6.7e+00</td><td> 21.1   </td></tr>
	<tr><td>agrFF002</td><td>agrFF040</td><td> 91.37  </td><td>255     </td><td>20      </td><td>2       </td><td>187     </td><td>439     </td><td>121     </td><td>375     </td><td>2.0e-100</td><td>354.0   </td></tr>
</tbody>
</table>



If the BLAST results are from a standalone BLAST, use the RFLP library's ‘read.blast’ function, to read the MyblastRes.txt BLAST file into a data frame:
* MyBLAST <- read.blast(file="MyblastRes.txt")<br>

In case you did not use the BLAST program but used the web-based tools to get the ‘myAlign.txt’ file, you can still read the file as a data frame：<br>
* MyBLAST2 <- read.csv(file="myAlign.txt", head=TRUE, sep=",")


```R
# The 'RFLPtool' package allows the computing/extracting of similar matrices out of the data frame:
mySimMat <- simMatrix(MyBlast)
head(mySimMat)
```


<table>
<thead><tr><th></th><th scope=col>agrFF002</th><th scope=col>agrFF003</th><th scope=col>agrFF005</th><th scope=col>agrFF023</th><th scope=col>agrFF036</th><th scope=col>agrFF040</th><th scope=col>agrFF042</th><th scope=col>agrFF043</th><th scope=col>agrFF044</th><th scope=col>agrFF049</th><th scope=col>...</th><th scope=col>agrFF185</th><th scope=col>agrFF190</th><th scope=col>agrFF191</th><th scope=col>agrFF192</th><th scope=col>agrFF200</th><th scope=col>agrFF202</th><th scope=col>agrFF203</th><th scope=col>agrFF204</th><th scope=col>agrFF206</th><th scope=col>agrFF212</th></tr></thead>
<tbody>
	<tr><th scope=row>agrFF002</th><td>1.0000000</td><td>0.5142857</td><td>0.3645418</td><td>0.6389497</td><td>0.6409692</td><td>0.5065217</td><td>0.5154867</td><td>0.4458078</td><td>0.5125348</td><td>0.3995680</td><td>...      </td><td>0.5087336</td><td>0.3682093</td><td>0.6198704</td><td>0.5200000</td><td>0.3549696</td><td>0.3809524</td><td>0.3825364</td><td>0.2904412</td><td>0.4478528</td><td>0.3731343</td></tr>
	<tr><th scope=row>agrFF003</th><td>0.5142857</td><td>1.0000000</td><td>0.3846154</td><td>0.9912088</td><td>0.9801762</td><td>0.9736264</td><td>0.9800885</td><td>0.7274725</td><td>0.6935933</td><td>0.4351648</td><td>...      </td><td>0.9846154</td><td>0.3846154</td><td>0.9780220</td><td>0.9866667</td><td>0.3956044</td><td>0.4329670</td><td>0.4329670</td><td>0.3252747</td><td>0.7076923</td><td>0.4175824</td></tr>
	<tr><th scope=row>agrFF005</th><td>0.3645418</td><td>0.3846154</td><td>1.0000000</td><td>0.3829322</td><td>0.3876652</td><td>0.3804348</td><td>0.3871681</td><td>0.3680982</td><td>0.5041783</td><td>0.7386609</td><td>...      </td><td>0.3820961</td><td>0.9979879</td><td>0.3779698</td><td>0.3888889</td><td>0.4036511</td><td>0.7039337</td><td>0.7068607</td><td>0.2988048</td><td>0.3701431</td><td>0.4243070</td></tr>
	<tr><th scope=row>agrFF023</th><td>0.6389497</td><td>0.9912088</td><td>0.3829322</td><td>1.0000000</td><td>0.9911894</td><td>0.9759300</td><td>0.9889381</td><td>0.7264770</td><td>0.8495822</td><td>0.4332604</td><td>...      </td><td>0.9868709</td><td>0.3829322</td><td>0.9737418</td><td>0.9977778</td><td>0.3938731</td><td>0.4310722</td><td>0.4310722</td><td>0.3238512</td><td>0.7067834</td><td>0.4157549</td></tr>
	<tr><th scope=row>agrFF036</th><td>0.6409692</td><td>0.9801762</td><td>0.3876652</td><td>0.9911894</td><td>1.0000000</td><td>0.9691630</td><td>0.9845133</td><td>0.7290749</td><td>0.8467967</td><td>0.4339207</td><td>...      </td><td>0.9801762</td><td>0.3876652</td><td>0.9669604</td><td>0.9888889</td><td>0.3986784</td><td>0.4317181</td><td>0.4317181</td><td>0.3259912</td><td>0.7092511</td><td>0.4207048</td></tr>
	<tr><th scope=row>agrFF040</th><td>0.5065217</td><td>0.9736264</td><td>0.3804348</td><td>0.9759300</td><td>0.9691630</td><td>1.0000000</td><td>0.9823009</td><td>0.7195652</td><td>0.8551532</td><td>0.4282609</td><td>...      </td><td>0.9825328</td><td>0.3804348</td><td>0.9695652</td><td>0.9800000</td><td>0.3782609</td><td>0.4260870</td><td>0.4260870</td><td>0.3217391</td><td>0.7000000</td><td>0.4108696</td></tr>
</tbody>
</table>



Then we can name the columns of the data per your convenience and use it for further analysis like the `colnames`、`dim.names` function.

The web-based BLAST algorithm primarily outputs an `HTML` file together with some other formats. However, the alignment is available for download in a tabular format. The recipe just explains a way to import these files as a data frame. The `read.blast` functions do the same with the file generated by the standalone BLAST. For instance, such data is available as an example in the `RFLPtools` package. 


```R
data(BLASTdata)
head(BLASTdata)
```


<table>
<thead><tr><th scope=col>query.id</th><th scope=col>subject.id</th><th scope=col>identity</th><th scope=col>alignment.length</th><th scope=col>mismatches</th><th scope=col>gap.opens</th><th scope=col>q.start</th><th scope=col>q.end</th><th scope=col>s.start</th><th scope=col>s.end</th><th scope=col>evalue</th><th scope=col>bit.score</th></tr></thead>
<tbody>
	<tr><td>agrFF002</td><td>agrFF002</td><td>100.00  </td><td>544     </td><td> 0      </td><td>0       </td><td>  1     </td><td>544     </td><td>  1     </td><td>544     </td><td> 0.0e+00</td><td>944.0   </td></tr>
	<tr><td>agrFF002</td><td>agrFF148</td><td> 93.42  </td><td>243     </td><td>14      </td><td>2       </td><td>199     </td><td>439     </td><td>671     </td><td>913     </td><td>6.0e-102</td><td>360.0   </td></tr>
	<tr><td>agrFF002</td><td>agrFF148</td><td>100.00  </td><td> 11     </td><td> 0      </td><td>0       </td><td>462     </td><td>472     </td><td>785     </td><td>795     </td><td> 6.7e+00</td><td> 21.1   </td></tr>
	<tr><td>agrFF002</td><td>agrFF176</td><td> 91.37  </td><td>255     </td><td>20      </td><td>2       </td><td>187     </td><td>439     </td><td>123     </td><td>377     </td><td>2.0e-100</td><td>354.0   </td></tr>
	<tr><td>agrFF002</td><td>agrFF176</td><td>100.00  </td><td> 11     </td><td> 0      </td><td>0       </td><td>462     </td><td>472     </td><td>250     </td><td>260     </td><td> 6.7e+00</td><td> 21.1   </td></tr>
	<tr><td>agrFF002</td><td>agrFF040</td><td> 91.37  </td><td>255     </td><td>20      </td><td>2       </td><td>187     </td><td>439     </td><td>121     </td><td>375     </td><td>2.0e-100</td><td>354.0   </td></tr>
</tbody>
</table>



This gives us an idea about the contents of the BLAST result data. The data consists of various elements, such as the name, identity, sequence length, alignment length, E-values, and some other attributes of the BLAST search. 



```R
colnames(BLASTdata)
```


<ol class=list-inline>
	<li>'query.id'</li>
	<li>'subject.id'</li>
	<li>'identity'</li>
	<li>'alignment.length'</li>
	<li>'mismatches'</li>
	<li>'gap.opens'</li>
	<li>'q.start'</li>
	<li>'q.end'</li>
	<li>'s.start'</li>
	<li>'s.end'</li>
	<li>'evalue'</li>
	<li>'bit.score'</li>
</ol>




```R

```
