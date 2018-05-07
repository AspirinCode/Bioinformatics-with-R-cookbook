
# Generating artificial expression data
Generally, to make it economical and available, analysts always use simulated data (artificial) sets to get well-characterized, synthetic datasets that allow thorough testing of learning algorithms in a fast and reproducible manner.<br>
一般来说，为了使其经济可用，分析师总是使用模拟数据（人工）集合来获得充分表征的合成数据集，从而以快速和可重现的方式对学习算法进行全面测试。<br>

Before we start generating the artificial data, we must know how do we want our data to look in terms of the bounds (upper and lower), fraction of DE genes, amount of data, and associated statistical parameters, which will be explained in the How it works...<br>
在我们开始生成人工数据之前，我们必须知道我们希望如何根据边界（上限和下限），DE基因的比例，数据量和相关统计参数来查看我们的数据...<br>

For example, we will generate a dataset of 35,000 genes with 1 percent of DE genes in this part.<br>
(DE means differentially expressed)<br>
举个例子，生成一个含有1％DE基因的35,000个基因的数据集。

In this part, `madsim` labrary is important. More details by typing `?madsim`.

1、Install and load the `madsim` library from the CRAN repository


```R
install.packages("madsim", repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
```

    Installing package into 'C:/Users/Master1/Documents/R/win-library/3.4'
    (as 'lib' is unspecified)
    

    package 'madsim' successfully unpacked and MD5 sums checked
    
    The downloaded binary packages are in
    	C:\Users\Master1\AppData\Local\Temp\RtmpYh7ETq\downloaded_packages
    


```R
library(madsim)
```

2、Define the first set of parameters for the simulation process for the overall statistical parameters and the distribution along arrays.<br>
(定义整个统计参数和阵列分布的仿真过程的第一组参数。)


```R
fparams <- data.frame(m1 = 7, m2 = 7, shape2 = 4, lb = 4, ub = 14, pde = 0.02, sym = 0.5)
fparams
```


<table>
<thead><tr><th scope=col>m1</th><th scope=col>m2</th><th scope=col>shape2</th><th scope=col>lb</th><th scope=col>ub</th><th scope=col>pde</th><th scope=col>sym</th></tr></thead>
<tbody>
	<tr><td>7   </td><td>7   </td><td>4   </td><td>4   </td><td>14  </td><td>0.02</td><td>0.5 </td></tr>
</tbody>
</table>



3、Define the second set of parameters that consists of the statistical parameters that define the level of expression in the genes.<br>
(对由定义基因中表达水平的统计参数组成的第二组参数进行定义)


```R
dparams <- data.frame(lambda1 = 0.13, lambda2 = 2, muminde = 1, sdde = 0.5)
dparams
```


<table>
<thead><tr><th scope=col>lambda1</th><th scope=col>lambda2</th><th scope=col>muminde</th><th scope=col>sdde</th></tr></thead>
<tbody>
	<tr><td>0.13</td><td>2   </td><td>1   </td><td>0.5 </td></tr>
</tbody>
</table>




```R
sdn <- 0.4
rseed <- 50
```

4、Define the number of genes required in the expression data.<br>
(定义报答数据中的基因数)


```R
n <- 35000
```

5、generate the synthetic data


```R
atidata <- madsim(mdata=NULL, n=35000, ratio=0, fparams, dparams, sdn, rseed)
```

6、look at the structure of the object created by `madsim`, we can find three components.


```R
str(atidata)
```

    List of 3
     $ xdata: num [1:35000, 1:14] 10.1 5.99 8.12 9.94 4.1 ...
      ..- attr(*, "dimnames")=List of 2
      .. ..$ : NULL
      .. ..$ : chr [1:14] "cont1" "cont2" "cont3" "cont4" ...
     $ xid  : num [1:35000, 1] 0 0 0 0 0 0 0 0 0 0 ...
     $ xsd  : num 2.5
    

7、Create an MA plot for any sample,（say, #sample 1）to visualize the data


```R
library(limma)
```


```R
plotMA(atidata[[1]], 1)
```


[MAplot](https://github.com/Chengshu21/Chapter-5-Analyzing-Microarray-Data-with-R/blob/master/md/pic/MAplot-Generating%20artificial%20expression%20data.png)



```R
boxplot(atidata[[1]], 1)
```


[boxplot](https://github.com/Chengshu21/Chapter-5-Analyzing-Microarray-Data-with-R/blob/master/md/pic/boxplot-Generating%20artificial%20expression%20data.png)

