
# 1、The SNP association analysis

GWAS is used to identify genetic variants associated with a disease phenotype. It's performed by statistical analysis of the GWAS data. To perform the SNP association analysis, we take the data built in package `SNPsocc` as an example.<br>

GWAS用于鉴定与疾病表型相关的遗传变异。它通过GWAS数据的统计分析来执行。为了实现SNP关联分析，我们以“SNPassoc”包中的数据为例。

1、Install a package `SNPssoc` for association studies, load the built-in data for this package, and take a look at the data objects.<br>
安装`SNPassoc`包并载入内置数据，查看数据中的对象（前六行）。


```R
source("http://www.bioconductor.org/biocLite.R")
biocLite("SNPassoc")
library(SNPassoc)
data(SNPs)               # load the data
table1 <- head(SNPs)               # take a look at the objects
table1
table2 <- head(SNPs.info.pos)
table2 
```


<table>
<thead><tr><th scope=col>id</th><th scope=col>casco</th><th scope=col>sex</th><th scope=col>blood.pre</th><th scope=col>protein</th><th scope=col>snp10001</th><th scope=col>snp10002</th><th scope=col>snp10003</th><th scope=col>snp10004</th><th scope=col>snp10005</th><th scope=col>...</th><th scope=col>snp100026</th><th scope=col>snp100027</th><th scope=col>snp100028</th><th scope=col>snp100029</th><th scope=col>snp100030</th><th scope=col>snp100031</th><th scope=col>snp100032</th><th scope=col>snp100033</th><th scope=col>snp100034</th><th scope=col>snp100035</th></tr></thead>
<tbody>
	<tr><td>1       </td><td>1       </td><td>Female  </td><td>13.7    </td><td>75640.52</td><td>TT      </td><td>CC      </td><td>GG      </td><td>GG      </td><td>GG      </td><td>...     </td><td>GG      </td><td>CC      </td><td>CC      </td><td>GG      </td><td>AA      </td><td>TT      </td><td>AA      </td><td>AA      </td><td>TT      </td><td>TT      </td></tr>
	<tr><td>2       </td><td>1       </td><td>Female  </td><td>12.7    </td><td>28688.22</td><td>TT      </td><td>AC      </td><td>GG      </td><td>GG      </td><td>AG      </td><td>...     </td><td>GG      </td><td>CG      </td><td>CT      </td><td>GG      </td><td>AA      </td><td>TT      </td><td>AG      </td><td>AG      </td><td>TT      </td><td>TT      </td></tr>
	<tr><td>3       </td><td>1       </td><td>Female  </td><td>12.9    </td><td>17279.59</td><td>TT      </td><td>CC      </td><td>GG      </td><td>GG      </td><td>GG      </td><td>...     </td><td>GG      </td><td>CC      </td><td>CC      </td><td>GG      </td><td>AA      </td><td>TT      </td><td>AA      </td><td>AA      </td><td>TT      </td><td>TT      </td></tr>
	<tr><td>4       </td><td>1       </td><td>Male    </td><td>14.6    </td><td>27253.99</td><td>CT      </td><td>CC      </td><td>GG      </td><td>GG      </td><td>GG      </td><td>...     </td><td>GG      </td><td>CC      </td><td>CT      </td><td>AG      </td><td>AA      </td><td>TT      </td><td>AG      </td><td>AG      </td><td>CT      </td><td>TT      </td></tr>
	<tr><td>5       </td><td>1       </td><td>Female  </td><td>13.4    </td><td>38066.57</td><td>TT      </td><td>AC      </td><td>GG      </td><td>GG      </td><td>GG      </td><td>...     </td><td>GG      </td><td>CG      </td><td>CT      </td><td>GG      </td><td>AA      </td><td>TT      </td><td>AG      </td><td>AG      </td><td>TT      </td><td>TT      </td></tr>
	<tr><td>6       </td><td>1       </td><td>Female  </td><td>11.3    </td><td> 9872.46</td><td>TT      </td><td>CC      </td><td>GG      </td><td>GG      </td><td>GG      </td><td>...     </td><td>GG      </td><td>CC      </td><td>CC      </td><td>GG      </td><td>AA      </td><td>TT      </td><td>AA      </td><td>AA      </td><td>TT      </td><td>NA      </td></tr>
</tbody>
</table>




<table>
<thead><tr><th scope=col>snp</th><th scope=col>chr</th><th scope=col>pos</th></tr></thead>
<tbody>
	<tr><td>snp10001</td><td>Chr1    </td><td>2987398 </td></tr>
	<tr><td>snp10002</td><td>Chr1    </td><td>1913558 </td></tr>
	<tr><td>snp10003</td><td>Chr1    </td><td>1982067 </td></tr>
	<tr><td>snp10004</td><td>Chr1    </td><td> 447403 </td></tr>
	<tr><td>snp10005</td><td>Chr1    </td><td>2212031 </td></tr>
	<tr><td>snp10006</td><td>Chr1    </td><td>2515720 </td></tr>
</tbody>
</table>



The SNP data is a data.frame object that consists of 157 rows and 40 columns. Each row of the data represents a sample that is either a case or control. The first five columns represent the attributes for each sample and the columns from 6:40 are the SNP information. That is why we used these columns while creating the SNP class object. The first five columns are for sample identifiers, for case or control identification (1 is for case and 0 is for control, which means casco). This is followed by information about the blood pressure and protein level of each gender. The table1 shows the first six entries of the SNP data.<br>
Table2 shows information about the SNP name, chromosome name, and genomic position.<br>

SNP数据是一个data.frame对象，由157行和40列组成。数据的每一行都代表一个案例或对照。前五列代表每个样本的属性，而6:40的列是SNP信息。这就是我们在创建SNP类对象时使用这些列的原因。 前五列用于样本标识符，用于案例或控制标识（1代表案例，0代表控制，代表casco）。随后是关于每个性别的血压和蛋白质水平的信息。表1显示了SNP数据的前六个条目。<br>
表2显示了关于SNP名称，染色体名称和基因组位置的信息。

2、Create an `SNP` object using the `setupSNP` function with the `SNP` information in columns 6 to 40.<br>
使用“setupSNP”函数和第6至40列中的“SNP”信息创建一个SNP对象。


```R
mySNP <- setupSNP(SNPs, 6:40, sep="")
head(mySNP)
```


<table>
<thead><tr><th scope=col>id</th><th scope=col>casco</th><th scope=col>sex</th><th scope=col>blood.pre</th><th scope=col>protein</th><th scope=col>snp10001</th><th scope=col>snp10002</th><th scope=col>snp10003</th><th scope=col>snp10004</th><th scope=col>snp10005</th><th scope=col>...</th><th scope=col>snp100026</th><th scope=col>snp100027</th><th scope=col>snp100028</th><th scope=col>snp100029</th><th scope=col>snp100030</th><th scope=col>snp100031</th><th scope=col>snp100032</th><th scope=col>snp100033</th><th scope=col>snp100034</th><th scope=col>snp100035</th></tr></thead>
<tbody>
	<tr><td>1       </td><td>1       </td><td>Female  </td><td>13.7    </td><td>75640.52</td><td>T/T     </td><td>C/C     </td><td>G/G     </td><td>G/G     </td><td>G/G     </td><td>...     </td><td>G/G     </td><td>C/C     </td><td>C/C     </td><td>G/G     </td><td>A/A     </td><td>T/T     </td><td>A/A     </td><td>A/A     </td><td>T/T     </td><td>T/T     </td></tr>
	<tr><td>2       </td><td>1       </td><td>Female  </td><td>12.7    </td><td>28688.22</td><td>T/T     </td><td>A/C     </td><td>G/G     </td><td>G/G     </td><td>A/G     </td><td>...     </td><td>G/G     </td><td>C/G     </td><td>C/T     </td><td>G/G     </td><td>A/A     </td><td>T/T     </td><td>A/G     </td><td>A/G     </td><td>T/T     </td><td>T/T     </td></tr>
	<tr><td>3       </td><td>1       </td><td>Female  </td><td>12.9    </td><td>17279.59</td><td>T/T     </td><td>C/C     </td><td>G/G     </td><td>G/G     </td><td>G/G     </td><td>...     </td><td>G/G     </td><td>C/C     </td><td>C/C     </td><td>G/G     </td><td>A/A     </td><td>T/T     </td><td>A/A     </td><td>A/A     </td><td>T/T     </td><td>T/T     </td></tr>
	<tr><td>4       </td><td>1       </td><td>Male    </td><td>14.6    </td><td>27253.99</td><td>C/T     </td><td>C/C     </td><td>G/G     </td><td>G/G     </td><td>G/G     </td><td>...     </td><td>G/G     </td><td>C/C     </td><td>C/T     </td><td>A/G     </td><td>A/A     </td><td>T/T     </td><td>A/G     </td><td>A/G     </td><td>C/T     </td><td>T/T     </td></tr>
	<tr><td>5       </td><td>1       </td><td>Female  </td><td>13.4    </td><td>38066.57</td><td>T/T     </td><td>A/C     </td><td>G/G     </td><td>G/G     </td><td>G/G     </td><td>...     </td><td>G/G     </td><td>C/G     </td><td>C/T     </td><td>G/G     </td><td>A/A     </td><td>T/T     </td><td>A/G     </td><td>A/G     </td><td>T/T     </td><td>T/T     </td></tr>
	<tr><td>6       </td><td>1       </td><td>Female  </td><td>11.3    </td><td> 9872.46</td><td>T/T     </td><td>C/C     </td><td>G/G     </td><td>G/G     </td><td>G/G     </td><td>...     </td><td>G/G     </td><td>C/C     </td><td>C/C     </td><td>G/G     </td><td>A/A     </td><td>T/T     </td><td>A/A     </td><td>A/A     </td><td>T/T     </td><td>NA      </td></tr>
</tbody>
</table>



3、To do an association test for the trait and variables of your interest with `association` function.<br>
用`association`函数对你感兴趣的特征和变体进行关联检验。


```R
myres <- association(casco~sex+snp10001+blood.pre, data = mySNP,
                     model.interaction = c("dominant","codominant"))
myres   # Take a look at the results
```


    
    SNP: snp10001  adjusted by: sex blood.pre 
                  0    %   1    %   OR lower upper p-value   AIC
    Codominant                                                  
    T/T          24 51.1  68 61.8 1.00             0.15410 195.8
    C/T          21 44.7  32 29.1 0.55  0.26  1.14              
    C/C           2  4.3  10  9.1 1.74  0.35  8.63              
    Dominant                                                    
    T/T          24 51.1  68 61.8 1.00             0.22859 196.1
    C/T-C/C      23 48.9  42 38.2 0.65  0.32  1.31              
    Recessive                                                   
    T/T-C/T      45 95.7 100 90.9 1.00             0.28494 196.4
    C/C           2  4.3  10  9.1 2.22  0.46 10.70              
    Overdominant                                                
    T/T-C/C      26 55.3  78 70.9 1.00             0.07188 194.3
    C/T          21 44.7  32 29.1 0.52  0.25  1.06              
    log-Additive                                                
    0,1,2        47 29.9 110 70.1 0.87  0.51  1.49 0.60861 197.3


The `association` function perfroms an association analysis between a single SNP and a dependent variable via model fitting and computing statistics. The arguments include the model to be fitted (formula object). This needs the train of interest and the SNP under investigation separated by a `~` operator. It is possible to add other variables to the formula using the `+` operator. Adding more SNPs is possible, but the analysis is done only for the first SNP and adjusted by the others. The `association` function can model the dependence based on the codominant, dominant, recessive, overdominant, and log-additive genetic models of inheritance. This is supplied to the function under the `model.interaction` argument.<br>

In this recipe, we tried to find the dependence of case and control on the first SNP—snp10001—the sex, and the blood pressure. We used two genetic models: dominant and codominant. The result looks like a matrix. The columns contain information about the sample size and percentages for each genotype, the odds ratio and its 95 percent confidence interval (which takes the most frequent homozygous genotype as the reference), and the p-value corresponding to the likelihood ratio test obtained from a comparison with the null model. Besides, the matrix also has the `Akaike Information Criterion`(AIC) of each genetic model. Thus, the `myres` table reflects the dependence of phenotypes on SNP
(together with other factors), which shows the result of the association analysis test for `snp10001`.<br>
There are other packages for performing similar analyses. Some examples are `snp.plotter` and `GenABEL`.<br>

`association`函数通过模型拟合和计算统计来进行单个SNP与因变量之间的关联分析。参数包括要拟合的模型（公式对象）。这需要利用`〜`运算符分隔正在调查的SNP和目标列表。使用`+`运算符将其他变量添加到公式中。可以添加更多的SNP，但分析仅对第一个SNP进行，并由其他人进行调整。 `association`函数可以基于遗传的共显性，显性，隐性，超显性和对数加性遗传模型对依赖性进行建模。 `model.interaction`参数下的函数支持建模。<br>

我们试图找到病例和对照中对第一个SNP-snp10001- 性别和血压的依赖性。这里使用了两种遗传模型：显性和共显性。结果看起来像一个矩阵。这些列包含关于每个基因型的样本量和百分比的信息，比值比和其95％置信区间（以最常见的纯合基因型作为参考），以及对应于似然比检验的p值 与空模型进行比较。此外，矩阵还具有每个遗传模型的“Akaike信息标准”（AIC）。因此，myres表反映了表型对SNP的依赖性（连同其他因素）并且显示`snp10001`的关联分析测试结果。<br>



```R

```
