
# 11、Analyzing ChipSeq data

Chromatin immuneprecipitation sequencing (ChipSeq) is a method to identify genome-wide DNA binding sites for a protein of interest. It is often used to determine the binding sites for transcription factors, DNA-binding enzymes, histones, chaperones, or nucleosomes. The workflow to produce the ChipSeq data starts from the cross-links bound proteins and chromatin. The chromatin is fragmented, and the DNA fragments bound to one protein are captured using an antibody specific to it. The ends of the captured fragments are sequenced using NGS. The computational mapping of the sequenced DNA leads to the identification of the genomic locations of these fragments, illuminating their role in DNA
protein interactions and epigenetics research. The `ChipSeq` data consists of short reads in a `FASTQ` file format. There is a short read after every fifth line of the file. This part will deal with the analysis of the ChipSeq data in R.<br>

Here requires the `chipseq` package and the built-in data in the package for demonstration purposes. And it also needs the mouse genome data for the biological mapping of results.<br>

染色质免疫沉淀测序（ChipSeq）是鉴定感兴趣蛋白质的全基因组DNA结合位点的方法。它通常用于确定转录因子，DNA结合酶，组蛋白，分子伴侣或核小体的结合位点。产生ChipSeq数据的工作流程从交联蛋白和染色质开始。染色质被片段化，并且使用对其特异性的抗体捕获与一种蛋白质结合的DNA片段。捕获片段的末端使用NGS进行测序。测序DNA的计算作图导致鉴定这些片段的基因组位置，阐明它们在DNA蛋白质相互作用和表观遗传学研究中的作用。`ChipSeq`数据由`FASTQ`文件格式的短读取组成。文件的每五行之后都有一个简短的读取。这部分内容是在R中分析ChipSeq数据。

这里需要`chipseq`包和包中的内置数据。而且它还需要将小鼠基因组数据用于结果的生物映射。

1、Load `chipseq` and `TxDb.Mmusculus.UCSC.mm9.knownGene` libraries.


```R
source("http://bioconductor.org/biocLite.R")
biocLite(c("chipseq", "TxDb.Mmusculus.UCSC.mm9.knownGene"))
```

    Bioconductor version 3.6 (BiocInstaller 1.28.0), ?biocLite for help
    BioC_mirror: https://bioconductor.org
    Using Bioconductor 3.6 (BiocInstaller 1.28.0), R 3.4.4 (2018-03-15).
    Installing package(s) 'chipseq', 'TxDb.Mmusculus.UCSC.mm9.knownGene'
    

    package 'chipseq' successfully unpacked and MD5 sums checked
    
    The downloaded binary packages are in
    	C:\Users\Master1\AppData\Local\Temp\RtmpyajxuD\downloaded_packages
    

    installing the source package 'TxDb.Mmusculus.UCSC.mm9.knownGene'
    
    Old packages: 'foreign', 'futile.options', 'htmlwidgets', 'httpuv', 'lambda.r',
      'party', 'R.oo', 'readxl', 'robustbase', 'XML'
    


```R
library(TxDb.Mmusculus.UCSC.mm9.knownGene)   # load the library
library(chipseq)    #load the library
data(cstest)        #use the built-in package dataset named `cstest`
cstest              #take a look at the data 
```


    GRangesList object of length 2:
    $ctcf 
    GRanges object with 450096 ranges and 0 metadata columns:
               seqnames                 ranges strand
                  <Rle>              <IRanges>  <Rle>
           [1]    chr10     [3012936, 3012959]      +
           [2]    chr10     [3012941, 3012964]      +
           [3]    chr10     [3012944, 3012967]      +
           [4]    chr10     [3012955, 3012978]      +
           [5]    chr10     [3012963, 3012986]      +
           ...      ...                    ...    ...
      [450092]    chr12 [121239376, 121239399]      -
      [450093]    chr12 [121245849, 121245872]      -
      [450094]    chr12 [121245895, 121245918]      -
      [450095]    chr12 [121246344, 121246367]      -
      [450096]    chr12 [121253499, 121253522]      -
    
    ...
    <1 more element>
    -------
    seqinfo: 35 sequences from an unspecified genome


2、Estimate the length of the fragments in the data.评估数据中片段的长度。


```R
estimate.mean.fraglen(cstest$ctcf)
```


<dl class=dl-horizontal>
	<dt>chr10</dt>
		<dd>179.679378486558</dd>
	<dt>chr11</dt>
		<dd>172.488350300294</dd>
	<dt>chr12</dt>
		<dd>181.673228507249</dd>
</dl>



3、Extend the fragments to cover the binding sites in the sequences. Use an extension length inspired by the fragment lengths in the last step (length =200 > estimate.mean.fraglen(cstest$ctcf)).<br>
扩展片段以覆盖序列中的结合位点。


```R
ctcf_ext <- resize(cstest$ctcf, width = 200)
gfp_ext <- resize(cstest$gfp, width = 200)
# A useful summary of this information is the coverage, that is, how many times each
# base in the genome was covered by one of these intervals and can be computed
cov_ctcf <- coverage(ctcf_ext)
cov_gfp <- coverage(gfp_ext)
cov_ctcf
```


    RleList of length 35
    $chr1
    integer-Rle of length 197195432 with 1 run
      Lengths: 197195432
      Values :         0
    
    $chr2
    integer-Rle of length 181748087 with 1 run
      Lengths: 181748087
      Values :         0
    
    $chr3
    integer-Rle of length 159599783 with 1 run
      Lengths: 159599783
      Values :         0
    
    $chr4
    integer-Rle of length 155630120 with 1 run
      Lengths: 155630120
      Values :         0
    
    $chr5
    integer-Rle of length 152537259 with 1 run
      Lengths: 152537259
      Values :         0
    
    ...
    <30 more elements>


4、Create a plot called islands for the regions of interest. They are contiguous segments.
对感兴趣的区域创建名为岛的图。他们是连续的部分。


```R
library(lattice)
par(mfrow = c(2, 1))
islandDepthPlot(cov_ctcf)
islandDepthPlot(cov_gfp)
# In the following plot, the x axis shows the depth, whereas the y axis shows the corresponding log counts 
# overall, the plot shows the coverage for the sample data
```






[output](https://github.com/Chengshu21/Chapter-8-Analyzing-NGS-Data/blob/master/MD/pic/output_2_Analyzing%20ChipSeq%20data.png)



[output](https://github.com/Chengshu21/Chapter-8-Analyzing-NGS-Data/blob/master/MD/pic/output_3_Analyzing%20ChipSeq%20data.png)


5、Compute the peak's cut off for a desired `fdr` (here, 0.01).<br>
计算所需fdr的峰值的阈值（此处为0.01）。


```R
peakCutoff(cov_ctcf, fdr = 0.01)
peakCutoff(cov_gfp, fdr = 0.01)
```


5.09190915706007



6.97927298645879


6、With this value, decide a cut-off value (7) to be used to get the peaks with high coverage in the segments of the data for both lanes, `ctcf` and `gfp`(note that the chosen cut off is based on computation above 5.09 and 6.97).<br>
用这个值，决定一个阈值（7），用来获得数据段中包含高覆盖率的峰值，ctcf和gfp（注意，选择的截止值是基于5.09以上的计算和6.97）。


```R
peaks_ctcf <- slice(cov_ctcf, lower = 7)
peaks_gfp <- slice(cov_gfp, lower = 7)
```

7、Compute the differential peaks with `diffPeakSummary` function to determine which peaks are different in the two samples. And visualize the results in terms of an XY plot.<br>
使用diffPeakSummary函数计算差分峰，以确定两个样本中哪些峰不同。并根据XY图显示结果。


```R
peakSummary <- diffPeakSummary(peaks_gfp, peaks_ctcf)
head(data.frame(peakSummary))
xyplot(asinh(sums2) ~ asinh(sums1) | space, data = as.data.frame(peakSummary))
# In the following plot, the summary of the peaks for three chromosomes
```

    Warning message in as.data.frame(x, row.names = NULL, optional = optional, ...):
    "'optional' and arguments in '...' ignored"


<table>
<thead><tr><th scope=col>space</th><th scope=col>start</th><th scope=col>end</th><th scope=col>width</th><th scope=col>comb.max</th><th scope=col>sums1</th><th scope=col>sums2</th><th scope=col>maxs1</th><th scope=col>maxs2</th></tr></thead>
<tbody>
	<tr><td>chr10  </td><td>3012944</td><td>3013140</td><td>197    </td><td>11     </td><td>  0    </td><td>1911   </td><td>0      </td><td>11     </td></tr>
	<tr><td>chr10  </td><td>3135027</td><td>3135029</td><td>  3    </td><td> 7     </td><td>  0    </td><td>  21   </td><td>0      </td><td> 7     </td></tr>
	<tr><td>chr10  </td><td>3234798</td><td>3234896</td><td> 99    </td><td>10     </td><td>  0    </td><td> 910   </td><td>0      </td><td>10     </td></tr>
	<tr><td>chr10  </td><td>3234924</td><td>3234933</td><td> 10    </td><td> 7     </td><td>  0    </td><td>  70   </td><td>0      </td><td> 7     </td></tr>
	<tr><td>chr10  </td><td>3270010</td><td>3270301</td><td>292    </td><td>20     </td><td>164    </td><td>4072   </td><td>1      </td><td>19     </td></tr>
	<tr><td>chr10  </td><td>3277660</td><td>3277861</td><td>202    </td><td>13     </td><td>  0    </td><td>1897   </td><td>0      </td><td>13     </td></tr>
</tbody>
</table>






[output](https://github.com/Chengshu21/Chapter-8-Analyzing-NGS-Data/blob/master/MD/pic/output_4_Analyzing%20ChipSeq%20data.png)


8、Here are the peaks. Find if they are in the region of interest (promoter region). And take a look at the peaks in the promoter region by looking at the object created in the previous step.<br>
这里是高峰。找到他们是否在所需要的区域（启动子区域）。通过查看上一步创建的对象，查看启动子区域中的峰。


```R
gregions <- transcripts(TxDb.Mmusculus.UCSC.mm9.knownGene)
promoters <- flank(gregions, 1000, both = TRUE)
peakSummary$inPromoter <- peakSummary %over% promoters
which(peakSummary$inPromoter)
```


<ol class=list-inline>
	<li>2</li>
	<li>14</li>
	<li>41</li>
	<li>51</li>
	<li>100</li>
	<li>101</li>
	<li>110</li>
	<li>124</li>
	<li>126</li>
	<li>130</li>
	<li>133</li>
	<li>141</li>
	<li>155</li>
	<li>156</li>
	<li>157</li>
	<li>158</li>
	<li>175</li>
	<li>176</li>
	<li>177</li>
	<li>180</li>
	<li>184</li>
	<li>185</li>
	<li>211</li>
	<li>212</li>
	<li>213</li>
	<li>217</li>
	<li>218</li>
	<li>257</li>
	<li>273</li>
	<li>276</li>
	<li>299</li>
	<li>306</li>
	<li>307</li>
	<li>315</li>
	<li>318</li>
	<li>319</li>
	<li>325</li>
	<li>327</li>
	<li>328</li>
	<li>358</li>
	<li>361</li>
	<li>362</li>
	<li>370</li>
	<li>371</li>
	<li>375</li>
	<li>381</li>
	<li>386</li>
	<li>432</li>
	<li>456</li>
	<li>460</li>
	<li>461</li>
	<li>493</li>
	<li>499</li>
	<li>502</li>
	<li>513</li>
	<li>527</li>
	<li>555</li>
	<li>556</li>
	<li>566</li>
	<li>595</li>
	<li>596</li>
	<li>604</li>
	<li>606</li>
	<li>607</li>
	<li>621</li>
	<li>622</li>
	<li>623</li>
	<li>628</li>
	<li>652</li>
	<li>657</li>
	<li>663</li>
	<li>665</li>
	<li>676</li>
	<li>696</li>
	<li>697</li>
	<li>701</li>
	<li>710</li>
	<li>715</li>
	<li>725</li>
	<li>824</li>
	<li>825</li>
	<li>877</li>
	<li>900</li>
	<li>914</li>
	<li>915</li>
	<li>918</li>
	<li>919</li>
	<li>925</li>
	<li>938</li>
	<li>942</li>
	<li>963</li>
	<li>976</li>
	<li>977</li>
	<li>993</li>
	<li>996</li>
	<li>998</li>
	<li>1008</li>
	<li>1027</li>
	<li>1033</li>
	<li>1034</li>
	<li>1036</li>
	<li>1037</li>
	<li>1040</li>
	<li>1051</li>
	<li>1052</li>
	<li>1055</li>
	<li>1056</li>
	<li>1057</li>
	<li>1059</li>
	<li>1067</li>
	<li>1068</li>
	<li>1071</li>
	<li>1072</li>
	<li>1076</li>
	<li>1079</li>
	<li>1080</li>
	<li>1087</li>
	<li>1091</li>
	<li>1108</li>
	<li>1109</li>
	<li>1110</li>
	<li>1112</li>
	<li>1114</li>
	<li>1115</li>
	<li>1124</li>
	<li>1127</li>
	<li>1131</li>
	<li>1144</li>
	<li>1159</li>
	<li>1166</li>
	<li>1175</li>
	<li>1182</li>
	<li>1195</li>
	<li>1211</li>
	<li>1226</li>
	<li>1263</li>
	<li>1288</li>
	<li>1289</li>
	<li>1290</li>
	<li>1292</li>
	<li>1308</li>
	<li>1339</li>
	<li>1346</li>
	<li>1370</li>
	<li>1386</li>
	<li>1412</li>
	<li>1413</li>
	<li>1414</li>
	<li>1425</li>
	<li>1436</li>
	<li>1463</li>
	<li>1491</li>
	<li>1495</li>
	<li>1496</li>
	<li>1497</li>
	<li>1498</li>
	<li>1549</li>
	<li>1556</li>
	<li>1557</li>
	<li>1573</li>
	<li>1580</li>
	<li>1590</li>
	<li>1600</li>
	<li>1608</li>
	<li>1609</li>
	<li>1620</li>
	<li>1647</li>
	<li>1651</li>
	<li>1685</li>
	<li>1686</li>
	<li>1705</li>
	<li>1715</li>
	<li>1716</li>
	<li>1721</li>
	<li>1740</li>
	<li>1741</li>
	<li>1776</li>
	<li>1783</li>
	<li>1795</li>
	<li>1796</li>
	<li>1797</li>
	<li>1802</li>
	<li>1803</li>
	<li>1804</li>
	<li>1831</li>
	<li>1834</li>
	<li>1836</li>
	<li>1838</li>
	<li>1845</li>
	<li>1847</li>
	<li>1848</li>
	<li>1850</li>
	<li>1855</li>
	<li>1857</li>
	<li>1865</li>
	<li>1866</li>
	<li>1869</li>
	<li>1870</li>
	<li>1875</li>
	<li>1878</li>
	<li>1881</li>
	<li>1882</li>
	<li>1883</li>
	<li>1937</li>
	<li>1938</li>
	<li>1939</li>
	<li>1940</li>
	<li>1941</li>
	<li>1942</li>
	<li>1943</li>
	<li>1944</li>
	<li>1966</li>
	<li>1967</li>
	<li>1968</li>
	<li>1969</li>
	<li>2034</li>
	<li>2035</li>
	<li>2036</li>
	<li>2037</li>
	<li>2038</li>
	<li>2039</li>
	<li>2040</li>
	<li>2051</li>
	<li>2053</li>
	<li>2058</li>
	<li>2063</li>
	<li>2066</li>
	<li>2081</li>
	<li>2093</li>
	<li>2094</li>
	<li>2095</li>
	<li>2098</li>
	<li>2102</li>
	<li>2127</li>
	<li>2134</li>
	<li>2137</li>
	<li>2150</li>
	<li>2153</li>
	<li>2154</li>
	<li>2159</li>
	<li>2164</li>
	<li>2172</li>
	<li>2182</li>
	<li>2183</li>
	<li>2202</li>
	<li>2218</li>
	<li>2235</li>
	<li>2236</li>
	<li>2237</li>
	<li>2238</li>
	<li>2240</li>
	<li>2248</li>
	<li>2268</li>
	<li>2300</li>
	<li>2302</li>
	<li>2330</li>
	<li>2334</li>
	<li>2336</li>
	<li>2338</li>
	<li>2339</li>
	<li>2365</li>
	<li>2367</li>
	<li>2372</li>
	<li>2375</li>
	<li>2378</li>
	<li>2379</li>
	<li>2380</li>
	<li>2412</li>
	<li>2413</li>
	<li>2419</li>
	<li>2452</li>
	<li>2503</li>
	<li>2560</li>
	<li>2586</li>
	<li>2591</li>
	<li>2592</li>
	<li>2637</li>
	<li>2656</li>
	<li>2659</li>
	<li>2666</li>
	<li>2667</li>
	<li>2675</li>
	<li>2697</li>
	<li>2703</li>
	<li>2711</li>
	<li>2712</li>
	<li>2714</li>
	<li>2723</li>
	<li>2763</li>
	<li>2770</li>
	<li>2771</li>
	<li>2785</li>
	<li>2804</li>
	<li>2832</li>
	<li>2846</li>
	<li>2893</li>
	<li>2899</li>
	<li>2900</li>
	<li>2911</li>
	<li>2918</li>
	<li>2919</li>
	<li>2929</li>
	<li>2932</li>
	<li>2946</li>
	<li>2948</li>
	<li>2954</li>
	<li>2970</li>
	<li>2974</li>
	<li>2988</li>
	<li>2996</li>
	<li>3006</li>
	<li>3022</li>
	<li>3029</li>
	<li>3034</li>
	<li>3036</li>
	<li>3047</li>
	<li>3048</li>
	<li>3129</li>
	<li>3132</li>
	<li>3138</li>
	<li>3149</li>
	<li>3150</li>
	<li>3158</li>
	<li>3162</li>
	<li>3164</li>
	<li>3172</li>
	<li>3176</li>
	<li>3190</li>
	<li>3193</li>
	<li>3210</li>
	<li>3211</li>
	<li>3215</li>
	<li>3216</li>
	<li>3218</li>
	<li>3226</li>
	<li>3229</li>
	<li>3230</li>
	<li>3233</li>
	<li>3236</li>
	<li>3242</li>
	<li>3245</li>
	<li>3255</li>
	<li>3256</li>
	<li>3260</li>
	<li>3264</li>
	<li>3277</li>
	<li>3280</li>
	<li>3281</li>
	<li>3290</li>
	<li>3310</li>
	<li>3316</li>
	<li>3317</li>
	<li>3318</li>
	<li>3326</li>
	<li>3329</li>
	<li>3330</li>
	<li>3331</li>
	<li>3334</li>
	<li>3353</li>
	<li>3356</li>
	<li>3357</li>
	<li>3389</li>
	<li>3391</li>
	<li>3392</li>
	<li>3393</li>
	<li>3394</li>
	<li>3402</li>
	<li>3423</li>
	<li>3427</li>
	<li>3432</li>
	<li>3442</li>
	<li>3452</li>
	<li>3479</li>
	<li>3488</li>
	<li>3497</li>
	<li>3504</li>
	<li>3505</li>
	<li>3526</li>
	<li>3560</li>
	<li>3564</li>
	<li>3570</li>
	<li>3573</li>
	<li>3579</li>
	<li>3586</li>
	<li>3587</li>
	<li>3588</li>
	<li>3598</li>
	<li>3599</li>
	<li>3631</li>
	<li>3632</li>
	<li>3633</li>
	<li>3635</li>
	<li>3638</li>
	<li>3648</li>
	<li>3680</li>
	<li>3688</li>
	<li>3693</li>
	<li>3719</li>
	<li>3730</li>
	<li>3751</li>
	<li>3754</li>
	<li>3755</li>
	<li>3766</li>
	<li>3768</li>
	<li>3771</li>
	<li>3778</li>
	<li>3792</li>
	<li>3793</li>
	<li>3801</li>
	<li>3810</li>
	<li>3847</li>
	<li>3862</li>
	<li>3864</li>
	<li>3865</li>
	<li>3866</li>
	<li>3870</li>
	<li>3871</li>
	<li>3884</li>
	<li>3890</li>
	<li>3895</li>
	<li>3900</li>
	<li>3902</li>
	<li>3903</li>
	<li>3919</li>
	<li>3959</li>
	<li>3961</li>
	<li>3975</li>
	<li>3996</li>
	<li>3997</li>
	<li>4000</li>
	<li>4001</li>
	<li>4006</li>
	<li>4007</li>
	<li>4008</li>
	<li>4016</li>
	<li>4017</li>
	<li>4043</li>
	<li>4050</li>
	<li>4051</li>
	<li>4055</li>
	<li>4061</li>
	<li>4062</li>
	<li>4067</li>
	<li>4068</li>
	<li>4076</li>
	<li>4078</li>
	<li>4094</li>
	<li>4101</li>
	<li>4102</li>
	<li>4130</li>
	<li>4133</li>
	<li>4182</li>
	<li>4184</li>
	<li>4186</li>
	<li>4191</li>
	<li>4192</li>
	<li>4199</li>
	<li>4201</li>
	<li>4217</li>
	<li>4221</li>
	<li>4226</li>
	<li>4230</li>
	<li>4236</li>
	<li>4238</li>
	<li>4244</li>
	<li>4263</li>
	<li>4265</li>
	<li>4266</li>
	<li>4270</li>
	<li>4280</li>
	<li>4281</li>
	<li>4290</li>
	<li>4295</li>
	<li>4301</li>
	<li>4321</li>
	<li>4322</li>
	<li>4329</li>
	<li>4331</li>
	<li>4332</li>
	<li>4333</li>
	<li>4340</li>
	<li>4341</li>
	<li>4356</li>
	<li>4368</li>
	<li>4371</li>
	<li>4380</li>
	<li>4381</li>
	<li>4389</li>
	<li>4402</li>
	<li>4411</li>
	<li>4414</li>
	<li>4415</li>
	<li>4419</li>
	<li>4442</li>
	<li>4445</li>
	<li>4477</li>
	<li>4493</li>
	<li>4495</li>
	<li>4499</li>
	<li>4536</li>
	<li>4558</li>
	<li>4559</li>
	<li>4562</li>
	<li>4602</li>
	<li>4603</li>
	<li>4605</li>
	<li>4622</li>
	<li>4628</li>
	<li>4630</li>
	<li>4631</li>
	<li>4632</li>
	<li>4633</li>
	<li>4634</li>
	<li>4640</li>
	<li>4645</li>
	<li>4655</li>
	<li>4659</li>
	<li>4661</li>
	<li>4668</li>
	<li>4672</li>
	<li>4675</li>
	<li>4676</li>
	<li>4678</li>
	<li>4679</li>
	<li>4695</li>
	<li>4696</li>
	<li>4699</li>
	<li>4706</li>
	<li>4707</li>
	<li>4713</li>
	<li>4733</li>
	<li>4755</li>
	<li>4756</li>
	<li>4759</li>
	<li>4761</li>
	<li>4762</li>
	<li>4765</li>
	<li>4766</li>
	<li>4768</li>
	<li>4769</li>
	<li>4787</li>
	<li>4792</li>
	<li>4816</li>
	<li>4820</li>
	<li>4825</li>
	<li>4836</li>
	<li>4846</li>
	<li>4850</li>
	<li>4861</li>
	<li>4862</li>
	<li>4874</li>
	<li>4875</li>
	<li>4892</li>
	<li>4894</li>
	<li>4895</li>
	<li>4897</li>
	<li>4912</li>
	<li>4914</li>
	<li>4919</li>
	<li>4920</li>
	<li>4938</li>
	<li>4939</li>
	<li>4946</li>
	<li>4950</li>
	<li>4952</li>
	<li>4966</li>
	<li>4979</li>
	<li>4980</li>
	<li>4982</li>
	<li>4989</li>
	<li>4994</li>
	<li>5011</li>
	<li>5012</li>
	<li>5017</li>
	<li>5019</li>
	<li>5027</li>
	<li>5037</li>
	<li>5050</li>
	<li>5066</li>
	<li>5079</li>
	<li>5100</li>
	<li>5103</li>
	<li>5124</li>
	<li>5132</li>
	<li>5133</li>
	<li>5139</li>
	<li>5140</li>
	<li>5153</li>
	<li>5165</li>
	<li>5166</li>
	<li>5180</li>
	<li>5181</li>
	<li>5229</li>
	<li>5246</li>
	<li>5247</li>
	<li>5260</li>
	<li>5273</li>
	<li>5308</li>
	<li>5312</li>
	<li>5335</li>
	<li>5362</li>
	<li>5395</li>
	<li>5405</li>
	<li>5410</li>
	<li>5411</li>
	<li>5412</li>
	<li>5413</li>
	<li>5438</li>
	<li>5441</li>
	<li>5455</li>
	<li>5456</li>
	<li>5467</li>
	<li>5479</li>
	<li>5480</li>
	<li>5489</li>
	<li>5491</li>
	<li>5505</li>
	<li>5506</li>
	<li>5525</li>
	<li>5526</li>
	<li>5535</li>
	<li>5567</li>
	<li>5603</li>
	<li>5604</li>
	<li>5627</li>
	<li>5645</li>
	<li>5647</li>
	<li>5651</li>
	<li>5656</li>
	<li>5657</li>
	<li>5660</li>
	<li>5661</li>
	<li>5674</li>
	<li>5683</li>
	<li>5684</li>
	<li>5697</li>
	<li>5698</li>
	<li>5699</li>
	<li>5701</li>
	<li>5706</li>
	<li>5728</li>
	<li>5731</li>
	<li>5732</li>
	<li>5739</li>
	<li>5740</li>
	<li>5758</li>
	<li>5759</li>
	<li>5763</li>
	<li>5779</li>
	<li>5792</li>
	<li>5794</li>
	<li>5815</li>
	<li>5825</li>
	<li>5909</li>
	<li>5947</li>
	<li>5948</li>
	<li>5962</li>
	<li>5979</li>
	<li>5996</li>
	<li>6022</li>
	<li>6047</li>
	<li>6068</li>
	<li>6130</li>
	<li>6131</li>
	<li>6148</li>
	<li>6161</li>
	<li>6182</li>
	<li>6224</li>
	<li>6259</li>
	<li>6267</li>
	<li>6268</li>
	<li>6309</li>
	<li>6313</li>
	<li>6318</li>
	<li>6322</li>
	<li>6327</li>
	<li>6328</li>
	<li>6349</li>
	<li>6361</li>
	<li>6362</li>
</ol>



In this part, the input data was from three chromosomes (10, 11, and 12) and in two lanes representing `ctcf` and `gfp` pull-down in mouse. The aim of the analysis is to find which peaks are different in two lanes (samples). The method explained in this part finds the coverage vectors for the data and computes statistical scores for the peaks therein. We start with finding the coverage for both the lanes by extending the reads to cover the binding sites. Thereafter, we define our peaks based on a cut-off value according to the desired `fdr` (0.01, as seen in step 5). The `diffPeakSummary` function then combines a set of peaks for the two lanes and summarizes them. Finally, we map the peaks onto the reference genome in step 7, and select the peaks that lie in the promoter region, which are of
interest to find the binding sites (step 8).<br>
在这部分中，输入数据来自三条染色体（10,11和12），并且在两条泳道中代表小鼠中的ctcf和gfp pull-down。分析的目的是找出两条泳道（样本）中哪些峰是不同的。本部分介绍的方法找到数据的覆盖矢量并计算其中的峰值的统计分数。我们首先通过扩展阅读来覆盖绑定站点，从而找到两条通道的覆盖范围。此后，我们根据所需的fdr（0.01，如步骤5所示）基于截止值定义我们的峰值。 diffPeakSummary函数然后为两条通道组合一组峰并汇总它们。最后，我们在步骤7中将峰映射到参考基因组上，并选择位于启动子区域中的峰，这些峰对于找到结合位点是有意义的（步骤8）。


```R

```
