
# Building the ExpressionSet object
The `ExpressionSet` class combines several different sources of information into one data structure including 
the intensities, phenotype data, and experiment information and annotation information. When we read a set of CEL 
files using the `ReadAffy` or `read.affyBatch` function, an `AffyBatch` object is created that extends the `ExpressionSet` structure. The `AffyBatch` object is probe-level data, whereas `ExpressionSet` is probeset-level data,
which is extended to a probe level by `AffyBatch`. The intensity values are always in the form of a table, matrix, or
data frame together with phenotype data, experiment details and annotations as separate objects (or files). We must create an `ExpressionSet` object from these individual files from scratch to facilitate the analysis work.<br>

ExpressionSet类将几种不同来源的信息组合成一个数据结构，包括强度，表型数据，实验信息和注释信息。 当我们使用ReadAffy或read.affyBatch函数读取一组`CEL`文件时，会创建一个扩展'ExpressionSet`结构的`AffyBatch`对象。 'AffyBatch`对象是探测针级别的数据，而`ExpressionSet`是探测针集合级别的数据，通过`AffyBatch`扩展到探测级别。 强度值始终以表格，矩阵或数据框架的形式与表型数据，实验细节和注释一起作为单独的对象（或文件）。我们必须从头开始从这些单独的文件创建一个`ExpressionSet`对象，以促进分析工作。<br>

In this part, we will use the package named `Biobase`. And the files contain different types of information, such as assay data, phenotypic metadata, feature annotations and metadata, and a description of the experiment. Take `myData` object as an example.<br>


在这一部分中，我们将使用名为`Biobase`的软件包。 这些文件包含不同类型的信息，如化验数据，表型元数据，特征注释和元数据以及实验描述。在这里以读取CEL文件中的对象`myData`为例。

1、install and load the `Biobase` library(安装并载入Biobase库).                         


```R
source("http://bioconductor.org/biocLite.R")
```

    Bioconductor version 3.6 (BiocInstaller 1.28.0), ?biocLite for help
    


```R
biocLite("Biobase")
```

    BioC_mirror: https://bioconductor.org
    Using Bioconductor 3.6 (BiocInstaller 1.28.0), R 3.4.3 (2017-11-30).
    Installing package(s) 'Biobase'
    Warning message:
    "package 'Biobase' is in use and will not be installed"Old packages: 'GenomicRanges'
    


```R
library(Biobase)
```

2、Load the built-in data for `Biobase` library as an demo expression file and the phenotypic data (pData) file.<br>
(`Biobase`库中的内置数据作为演示表达式文件和表型数据（pData）文件)


```R
DIR <- system.file("extdata", package="Biobase")
exprsLoc <- file.path(DIR, "exprsData.txt") 
pDataLoc <- file.path(DIR, "pData.txt")
DIR
exprsLoc
pDataLoc
```


'C:/Users/Master1/Documents/R/win-library/3.4/Biobase/extdata'



'C:/Users/Master1/Documents/R/win-library/3.4/Biobase/extdata/exprsData.txt'



'C:/Users/Master1/Documents/R/win-library/3.4/Biobase/extdata/pData.txt'


3、Read the table from the text file that contains the expression values using the usual `read.table` or `read.csv` function.<br>(使用通常的read.table或read.csv函数从包含表达式值的文本文件中读取表)


```R
exprs <- as.matrix(read.csv(exprsLoc, header = TRUE, sep = "\t", row.names = 1, as.is = TRUE))
exprs
```


<table>
<thead><tr><th></th><th scope=col>A</th><th scope=col>B</th><th scope=col>C</th><th scope=col>D</th><th scope=col>E</th><th scope=col>F</th><th scope=col>G</th><th scope=col>H</th><th scope=col>I</th><th scope=col>J</th><th scope=col>...</th><th scope=col>Q</th><th scope=col>R</th><th scope=col>S</th><th scope=col>T</th><th scope=col>U</th><th scope=col>V</th><th scope=col>W</th><th scope=col>X</th><th scope=col>Y</th><th scope=col>Z</th></tr></thead>
<tbody>
	<tr><th scope=row>AFFX-MurIL2_at</th><td> 192.7420  </td><td>  85.75330 </td><td> 176.7570  </td><td> 135.57500 </td><td>  64.49390 </td><td>  76.35690 </td><td>  160.5050 </td><td>  65.96310 </td><td>  56.90390 </td><td> 135.608000</td><td>...        </td><td> 179.84500 </td><td> 152.46700 </td><td> 180.83400 </td><td>  85.41460 </td><td> 157.98900 </td><td>  146.8000 </td><td>  93.8829  </td><td> 103.85500 </td><td>  64.4340  </td><td> 175.61500 </td></tr>
	<tr><th scope=row>AFFX-MurIL10_at</th><td>  97.1370  </td><td> 126.19600 </td><td>  77.9216  </td><td>  93.37130 </td><td>  24.39860 </td><td>  85.50880 </td><td>   98.9086 </td><td>  81.69320 </td><td>  97.80150 </td><td>  90.483800</td><td>...        </td><td>  87.68060 </td><td> 108.03200 </td><td> 134.26300 </td><td>  91.40310 </td><td>  -8.68811 </td><td>   85.0212 </td><td>  79.2998  </td><td>  71.65520 </td><td>  64.2369  </td><td>  78.70680 </td></tr>
	<tr><th scope=row>AFFX-MurIL4_at</th><td>  45.8192  </td><td>   8.83135 </td><td>  33.0632  </td><td>  28.70720 </td><td>   5.94492 </td><td>  28.29250 </td><td>   30.9694 </td><td>  14.79230 </td><td>  14.23990 </td><td>  34.487400</td><td>...        </td><td>  32.79110 </td><td>  33.52920 </td><td>  19.81720 </td><td>  20.41900 </td><td>  26.87200 </td><td>   31.1488 </td><td>  22.3420  </td><td>  19.01350 </td><td>  12.1686  </td><td>  17.37800 </td></tr>
	<tr><th scope=row>AFFX-MurFAS_at</th><td>  22.5445  </td><td>   3.60093 </td><td>  14.6883  </td><td>  12.33970 </td><td>  36.86630 </td><td>  11.25680 </td><td>   23.0034 </td><td>  16.21340 </td><td>  12.03750 </td><td>   4.549780</td><td>...        </td><td>  15.94880 </td><td>  14.67530 </td><td>  -7.91911 </td><td>  12.88750 </td><td>  11.91860 </td><td>   12.8324 </td><td>  11.1390  </td><td>   7.55564 </td><td>  19.9849  </td><td>   8.96849 </td></tr>
	<tr><th scope=row>AFFX-BioB-5_at</th><td>  96.7875  </td><td>  30.43800 </td><td>  46.1271  </td><td>  70.93190 </td><td>  56.17440 </td><td>  42.67560 </td><td>   86.5156 </td><td>  30.79270 </td><td>  19.71830 </td><td>  46.352000</td><td>...        </td><td>  58.62390 </td><td> 114.06200 </td><td>  93.44020 </td><td>  22.51680 </td><td>  48.64620 </td><td>   90.2215 </td><td>  42.0053  </td><td>  57.57380 </td><td>  44.8216  </td><td>  61.70440 </td></tr>
	<tr><th scope=row>AFFX-BioB-M_at</th><td>  89.0730  </td><td>  25.84610 </td><td>  57.2033  </td><td>  69.97660 </td><td>  49.58220 </td><td>  26.12620 </td><td>   75.0083 </td><td>  42.33520 </td><td>  41.12070 </td><td>  91.530700</td><td>...        </td><td>  58.13310 </td><td> 104.12200 </td><td> 115.83100 </td><td>  58.12240 </td><td>  73.42210 </td><td>   64.6066 </td><td>  40.3068  </td><td>  41.82090 </td><td>  46.1087  </td><td>  49.41220 </td></tr>
	<tr><th scope=row>AFFX-BioB-3_at</th><td> 265.9640  </td><td> 181.08000 </td><td> 164.9260  </td><td> 161.46900 </td><td> 236.97600 </td><td> 156.80300 </td><td>  211.2570 </td><td> 235.99400 </td><td> 175.64000 </td><td> 229.671000</td><td>...        </td><td> 192.22100 </td><td> 305.56700 </td><td> 300.68900 </td><td> 146.08100 </td><td> 142.91300 </td><td>  187.1320 </td><td> 170.5830  </td><td> 133.27900 </td><td> 187.4070  </td><td> 144.78400 </td></tr>
	<tr><th scope=row>AFFX-BioC-5_at</th><td> 110.1360  </td><td>  57.28890 </td><td>  67.3980  </td><td>  77.22070 </td><td>  41.34880 </td><td>  37.97800 </td><td>  110.5510 </td><td>  47.76900 </td><td>  24.78750 </td><td>  66.730200</td><td>...        </td><td>  53.27110 </td><td> 107.23700 </td><td> 119.66600 </td><td>  24.06540 </td><td>  98.84250 </td><td>   92.0846 </td><td>  53.3866  </td><td>  52.01640 </td><td>  65.9154  </td><td>  75.00430 </td></tr>
	<tr><th scope=row>AFFX-BioC-3_at</th><td>  43.0794  </td><td>  16.80060 </td><td>  37.6002  </td><td>  46.52720 </td><td>  22.24750 </td><td>  61.64010 </td><td>   33.6623 </td><td>  31.44230 </td><td>  23.10080 </td><td>  39.741900</td><td>...        </td><td>  57.50780 </td><td>  41.13370 </td><td>  79.98290 </td><td>  23.49530 </td><td>  51.56090 </td><td>   48.1247 </td><td>  31.8358  </td><td>  29.92640 </td><td>  37.8611  </td><td>  60.47720 </td></tr>
	<tr><th scope=row>AFFX-BioDn-5_at</th><td>  10.9187  </td><td>  16.17890 </td><td>  10.1495  </td><td>   9.73639 </td><td>  16.90280 </td><td>   5.33328 </td><td>   25.1182 </td><td>  38.75760 </td><td>  31.40410 </td><td>   0.398779</td><td>...        </td><td>  21.50910 </td><td>   3.10536 </td><td>   5.95347 </td><td>   5.66012 </td><td>  52.93380 </td><td>   15.7267 </td><td>  15.2116  </td><td>  -5.35282 </td><td>  13.1884  </td><td>  10.03850 </td></tr>
	<tr><th scope=row>AFFX-BioDn-3_at</th><td> 751.2270  </td><td> 515.00400 </td><td> 622.9010  </td><td> 669.85900 </td><td> 414.16500 </td><td> 654.07800 </td><td>  704.7810 </td><td> 472.08700 </td><td> 456.49600 </td><td> 601.335000</td><td>...        </td><td> 401.43000 </td><td> 757.49500 </td><td> 595.90800 </td><td> 381.23000 </td><td> 501.74400 </td><td>  659.6130 </td><td> 590.1560  </td><td> 461.23000 </td><td> 367.4330  </td><td> 790.94300 </td></tr>
	<tr><th scope=row>AFFX-CreX-5_at</th><td>  76.9437  </td><td>  40.90700 </td><td>  62.0314  </td><td>  54.42180 </td><td>  29.07040 </td><td>  19.52710 </td><td>   56.3164 </td><td>  36.20440 </td><td>  34.41180 </td><td>  54.076500</td><td>...        </td><td>  57.84270 </td><td>  83.19140 </td><td>  66.67830 </td><td>  24.88520 </td><td>  61.95480 </td><td>   58.0325 </td><td>  28.7707  </td><td>  48.53460 </td><td>  31.1489  </td><td>  72.45670 </td></tr>
	<tr><th scope=row>AFFX-CreX-3_at</th><td> 105.3780  </td><td>  97.49320 </td><td>  74.0299  </td><td>  54.52770 </td><td>  54.98490 </td><td>  58.08770 </td><td>   96.6320 </td><td>  52.73100 </td><td>  35.45880 </td><td>  60.264200</td><td>...        </td><td>  53.48370 </td><td> 108.54500 </td><td> 136.04400 </td><td>  43.86190 </td><td>  49.82890 </td><td>   88.3787 </td><td>  91.5539  </td><td>  40.22980 </td><td>  36.2151  </td><td>  73.53090 </td></tr>
	<tr><th scope=row>AFFX-BioB-5_st</th><td>  40.4826  </td><td>   7.45801 </td><td>  19.4069  </td><td>  20.62460 </td><td>  25.04960 </td><td>  12.48040 </td><td>   21.9102 </td><td>  23.77200 </td><td>  24.18400 </td><td>  29.703200</td><td>...        </td><td>   6.53565 </td><td>  42.46960 </td><td>  41.46690 </td><td>  21.65480 </td><td>  34.61080 </td><td>   36.9725 </td><td>  16.9959  </td><td>   8.39896 </td><td>  41.0271  </td><td>  21.08800 </td></tr>
	<tr><th scope=row>AFFX-BioB-M_st</th><td>  58.1706  </td><td>  15.79260 </td><td>  25.1962  </td><td>  46.50570 </td><td>  15.31570 </td><td>  16.68330 </td><td>   93.1759 </td><td>  -2.28600 </td><td>   9.00485 </td><td>  13.125300</td><td>...        </td><td>  26.92140 </td><td>  40.38730 </td><td>  17.68820 </td><td>  38.31000 </td><td>  18.02970 </td><td>   28.0120 </td><td>  24.8743  </td><td>  17.74280 </td><td>  34.7923  </td><td>  40.22130 </td></tr>
	<tr><th scope=row>AFFX-BioB-3_st</th><td> 257.6190  </td><td> 113.69000 </td><td> 187.7960  </td><td> 210.58000 </td><td> 137.39000 </td><td> 104.15900 </td><td>  296.2870 </td><td> 110.53600 </td><td> 123.76700 </td><td> 165.210000</td><td>...        </td><td> 204.75000 </td><td> 265.77100 </td><td> 317.31400 </td><td>  88.07730 </td><td> 156.17900 </td><td>  249.8720 </td><td> 141.8030  </td><td> 117.47300 </td><td> 126.1140  </td><td> 210.56100 </td></tr>
	<tr><th scope=row>AFFX-BioC-5_st</th><td> 129.0560  </td><td>  74.60950 </td><td>  82.8271  </td><td> 101.53400 </td><td>  83.49860 </td><td>  73.19860 </td><td>  110.6310 </td><td> 116.74200 </td><td> 149.32900 </td><td> 113.737000</td><td>...        </td><td> 103.28900 </td><td> 140.76000 </td><td> 177.44100 </td><td>  75.88880 </td><td>  87.36990 </td><td>  130.5090 </td><td>  91.7096  </td><td>  77.92770 </td><td> 129.6270  </td><td>  91.78030 </td></tr>
	<tr><th scope=row>AFFX-BioC-3_st</th><td>  61.7251  </td><td>  50.23720 </td><td>  61.6710  </td><td>  93.22350 </td><td>  38.11300 </td><td>  51.08690 </td><td>   69.0242 </td><td>  51.73520 </td><td>  48.49430 </td><td>  66.332400</td><td>...        </td><td>  74.19100 </td><td>  74.71060 </td><td> 112.96400 </td><td>  63.23490 </td><td>  37.48920 </td><td>   69.9946 </td><td>  42.8123  </td><td>  63.80590 </td><td>  50.1246  </td><td> 134.58800 </td></tr>
	<tr><th scope=row>AFFX-BioDn-5_st</th><td> -40.9349  </td><td> -83.93020 </td><td> -28.7050  </td><td> -27.99790 </td><td> -29.90970 </td><td> -26.90040 </td><td>  -45.6312 </td><td> -62.94740 </td><td> -31.43590 </td><td> -26.325300</td><td>...        </td><td> -30.48800 </td><td> -54.56340 </td><td> -49.08790 </td><td> -22.59160 </td><td> -21.23940 </td><td>  -38.7340 </td><td> -44.2136  </td><td> -23.69930 </td><td> -38.6532  </td><td> -40.72510 </td></tr>
	<tr><th scope=row>AFFX-BioDn-3_st</th><td> 284.4070  </td><td> 208.09900 </td><td> 239.0390  </td><td> 236.42800 </td><td> 152.32700 </td><td> 159.50500 </td><td>  316.9310 </td><td> 152.18800 </td><td> 182.80300 </td><td> 275.020000</td><td>...        </td><td> 144.42100 </td><td> 316.03700 </td><td> 269.48500 </td><td> 148.11400 </td><td> 193.20500 </td><td>  294.2380 </td><td> 197.7450  </td><td> 203.53100 </td><td> 181.6230  </td><td> 234.95300 </td></tr>
	<tr><th scope=row>AFFX-CreX-5_st</th><td> 178.7450  </td><td> 101.30000 </td><td> 118.6990  </td><td> 131.83400 </td><td> 109.35500 </td><td>  98.17990 </td><td>  177.5330 </td><td> 124.79500 </td><td>  86.07680 </td><td> 143.596000</td><td>...        </td><td> 107.08800 </td><td> 188.79200 </td><td> 240.35000 </td><td>  94.97540 </td><td> 104.14700 </td><td>  168.9920 </td><td> 109.5860  </td><td> 106.75700 </td><td>  90.1958  </td><td> 116.06900 </td></tr>
	<tr><th scope=row>AFFX-CreX-3_st</th><td>  79.7368  </td><td>  55.56320 </td><td>  68.5976  </td><td>  55.68810 </td><td>  56.39600 </td><td>  31.30030 </td><td>   84.8437 </td><td>  33.42830 </td><td>  42.31720 </td><td> 124.882000</td><td>...        </td><td>  36.64480 </td><td>  94.36010 </td><td>  71.47270 </td><td>  24.36270 </td><td>  48.48060 </td><td>   79.3437 </td><td>  42.6647  </td><td> -12.64890 </td><td>  37.4389  </td><td>  62.07060 </td></tr>
	<tr><th scope=row>AFFX-hum_alu_at</th><td>9903.1900  </td><td>8501.62000 </td><td>9453.0000  </td><td>8595.65000 </td><td>9198.53000 </td><td>8729.83000 </td><td>10085.3000 </td><td>5398.15000 </td><td>7851.25000 </td><td>9906.750000</td><td>...        </td><td>6945.46000 </td><td>9186.23000 </td><td>9889.05000 </td><td>8872.19000 </td><td>9682.71000 </td><td>10396.1000 </td><td>8365.8000  </td><td>9345.94000 </td><td>8252.6500  </td><td>4652.85000 </td></tr>
	<tr><th scope=row>AFFX-DapX-5_at</th><td>  61.2671  </td><td>  37.47400 </td><td>  44.7525  </td><td>  43.90200 </td><td>  40.56370 </td><td>  28.58190 </td><td>   49.2893 </td><td>   7.59488 </td><td>  23.62900 </td><td>  57.042900</td><td>...        </td><td>  26.35510 </td><td>  70.37660 </td><td>  72.77700 </td><td>  16.69440 </td><td>  36.38470 </td><td>   58.1014 </td><td>  40.3792  </td><td>  33.92190 </td><td>  23.7046  </td><td>  45.14850 </td></tr>
	<tr><th scope=row>AFFX-DapX-M_at</th><td> 120.5440  </td><td>  75.98540 </td><td> 126.3740  </td><td>  90.40210 </td><td>  99.62140 </td><td>  59.88540 </td><td>  129.4190 </td><td>  52.93500 </td><td>  64.18610 </td><td> 101.061000</td><td>...        </td><td>  72.23230 </td><td> 131.60200 </td><td> 136.20800 </td><td>  25.96240 </td><td>  76.89740 </td><td>  104.4280 </td><td>  68.4470  </td><td>  84.32500 </td><td>  46.6156  </td><td>  74.53150 </td></tr>
	<tr><th scope=row>AFFX-DapX-3_at</th><td>  50.0962  </td><td>  27.95320 </td><td>  29.2610  </td><td>  38.64360 </td><td>  34.48540 </td><td>  15.93390 </td><td>   55.8445 </td><td>  20.09040 </td><td>  24.73830 </td><td>  38.872500</td><td>...        </td><td>  34.97870 </td><td>  56.82460 </td><td>  46.60260 </td><td>  14.22910 </td><td>  24.45630 </td><td>   59.6125 </td><td>  22.6792  </td><td>  26.54100 </td><td>  30.6529  </td><td>  29.44770 </td></tr>
	<tr><th scope=row>AFFX-LysX-5_at</th><td>  42.5285  </td><td>  33.71860 </td><td>  35.8420  </td><td>  43.81730 </td><td>  21.10380 </td><td>  26.00260 </td><td>   47.7015 </td><td>  23.80350 </td><td>  11.37370 </td><td>  40.582200</td><td>...        </td><td>  31.92710 </td><td>  40.08590 </td><td>  49.11440 </td><td>  26.78740 </td><td>  50.96660 </td><td>   27.7013 </td><td>  30.9739  </td><td>  20.44470 </td><td>  14.6951  </td><td>  33.85440 </td></tr>
	<tr><th scope=row>AFFX-LysX-M_at</th><td>  36.8936  </td><td>  35.16970 </td><td>  36.5703  </td><td>  37.02740 </td><td>  24.05680 </td><td>  19.86490 </td><td>   57.4157 </td><td>  18.31420 </td><td>  13.96590 </td><td>  37.389700</td><td>...        </td><td>  39.36620 </td><td>  48.33230 </td><td>  56.42690 </td><td>  18.70090 </td><td>  51.72850 </td><td>   60.7026 </td><td>  43.9691  </td><td>  12.95760 </td><td>  16.9243  </td><td>  35.29280 </td></tr>
	<tr><th scope=row>AFFX-LysX-3_at</th><td> 234.6980  </td><td> 102.46700 </td><td>  97.9010  </td><td> 146.23900 </td><td> 127.06800 </td><td>  65.07980 </td><td>  262.5790 </td><td>  67.98070 </td><td>  67.95660 </td><td>  85.489600</td><td>...        </td><td>  88.35100 </td><td> 187.83100 </td><td> 247.96600 </td><td>  67.10260 </td><td>  73.02550 </td><td>  189.1860 </td><td> 104.8180  </td><td> 129.61900 </td><td>  91.4007  </td><td> 102.62100 </td></tr>
	<tr><th scope=row>AFFX-PheX-5_at</th><td>  26.9561  </td><td>  -2.37297 </td><td>  34.4333  </td><td>  17.59470 </td><td>  30.90680 </td><td>  19.45640 </td><td>   18.5628 </td><td>   5.85058 </td><td>  26.22780 </td><td> -16.781100</td><td>...        </td><td>  20.20630 </td><td>  40.03670 </td><td>  30.28970 </td><td>  22.47410 </td><td>  24.21060 </td><td>   37.1239 </td><td>  23.3193  </td><td>  39.30200 </td><td>  25.4601  </td><td>  28.48380 </td></tr>
	<tr><th scope=row>...</th><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>   </td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td></tr>
	<tr><th scope=row>31710_at</th><td> 545.37700  </td><td> 508.60800  </td><td> 634.472000 </td><td> 548.60300  </td><td> 440.82200  </td><td> 442.06700  </td><td> 667.5870   </td><td> 396.013000 </td><td> 354.5690000</td><td> 589.90000  </td><td>...         </td><td> 433.1090   </td><td> 567.512000 </td><td> 486.907000 </td><td> 464.39400  </td><td> 613.63700  </td><td> 590.86500  </td><td> 390.3610000</td><td> 382.37200  </td><td> 376.90900  </td><td> 611.53400  </td></tr>
	<tr><th scope=row>31711_at</th><td>  66.11170  </td><td>  93.93000  </td><td>  64.739200 </td><td>  31.17310  </td><td> 101.20600  </td><td>  55.94350  </td><td>  45.1672   </td><td>  56.418700 </td><td>  97.9919000</td><td>  78.50150  </td><td>...         </td><td>  63.4572   </td><td>  93.040500 </td><td>  64.478300 </td><td> 110.78000  </td><td>  52.57460  </td><td>  83.57380  </td><td>  66.1281000</td><td>  50.58970  </td><td>  92.05640  </td><td> 133.23500  </td></tr>
	<tr><th scope=row>31712_at</th><td>  82.63310  </td><td>  31.55360  </td><td>  62.054300 </td><td>  80.84610  </td><td>  71.89850  </td><td>  39.88630  </td><td>  90.5216   </td><td>  31.348600 </td><td>  39.7226000</td><td>  39.03820  </td><td>...         </td><td>  54.1816   </td><td>  88.320000 </td><td>  66.735400 </td><td>  19.73900  </td><td>  39.55710  </td><td>  96.21790  </td><td>  66.8283000</td><td>  52.39640  </td><td>  59.05830  </td><td>  53.93870  </td></tr>
	<tr><th scope=row>31713_s_at</th><td>  54.94340  </td><td>  31.03550  </td><td>  37.418400 </td><td>  33.72460  </td><td>  27.30930  </td><td>  38.20140  </td><td>  47.5207   </td><td>  48.081400 </td><td>  -4.7357900</td><td>  42.76900  </td><td>...         </td><td>  33.1410   </td><td>  58.962000 </td><td>  46.408600 </td><td>  21.77740  </td><td>  49.83840  </td><td>  65.20520  </td><td>  41.6442000</td><td>  32.07010  </td><td>  31.76400  </td><td>  41.82590  </td></tr>
	<tr><th scope=row>31714_at</th><td>  11.35110  </td><td>   3.83873  </td><td>  17.200200 </td><td>  10.11530  </td><td>   2.20914  </td><td>   6.79520  </td><td>  17.3392   </td><td>  13.669000 </td><td>  71.3370000</td><td>  29.06290  </td><td>...         </td><td>  18.6417   </td><td>  30.366000 </td><td>  33.819600 </td><td>  -1.61634  </td><td> -12.18010  </td><td>  28.10820  </td><td>  17.8877000</td><td>   5.56768  </td><td>  16.22110  </td><td>   6.69247  </td></tr>
	<tr><th scope=row>31715_at</th><td> 259.18300  </td><td> 230.84900  </td><td> 260.407000 </td><td> 268.56100  </td><td> 193.89100  </td><td> 184.69400  </td><td> 314.5470   </td><td> 168.408000 </td><td> 140.9250000</td><td> 276.21100  </td><td>...         </td><td> 195.8320   </td><td> 271.145000 </td><td> 233.154000 </td><td> 189.76000  </td><td> 195.25500  </td><td> 274.23900  </td><td> 207.3590000</td><td> 212.71200  </td><td> 215.90800  </td><td> 225.88800  </td></tr>
	<tr><th scope=row>31716_at</th><td> 153.75900  </td><td> 242.94700  </td><td> 185.728000 </td><td> 155.21200  </td><td> 125.92500  </td><td> 172.74700  </td><td> 112.6270   </td><td> 159.795000 </td><td> 253.2670000</td><td> 135.46500  </td><td>...         </td><td> 131.3290   </td><td> 149.116000 </td><td> 136.060000 </td><td> 178.32000  </td><td> 178.57700  </td><td> 144.53200  </td><td> 185.9380000</td><td> 106.36600  </td><td> 142.98400  </td><td> 150.87600  </td></tr>
	<tr><th scope=row>31717_at</th><td>  16.88730  </td><td>  13.91510  </td><td>  16.128700 </td><td>   5.48352  </td><td>  20.23880  </td><td>  11.69760  </td><td>  11.9066   </td><td>  29.286000 </td><td>  39.1377000</td><td>   3.88573  </td><td>...         </td><td>  26.9925   </td><td>  27.749700 </td><td>  47.679700 </td><td>  10.17730  </td><td>  11.73180  </td><td>  18.84590  </td><td>  15.4422000</td><td>   9.38988  </td><td>  16.40430  </td><td>  19.33940  </td></tr>
	<tr><th scope=row>31718_at</th><td>  30.25690  </td><td>   3.29854  </td><td>  21.148000 </td><td>  16.94180  </td><td>  16.03220  </td><td>  -5.41524  </td><td>  15.2477   </td><td>   0.992278 </td><td>  -0.0626833</td><td>  12.61000  </td><td>...         </td><td>  15.6661   </td><td>  28.916000 </td><td>   6.527990 </td><td>  29.46810  </td><td>  12.52550  </td><td>  15.91590  </td><td>   0.0721809</td><td> -15.67780  </td><td>   7.04334  </td><td>  17.73640  </td></tr>
	<tr><th scope=row>31719_at</th><td>1505.73000  </td><td>1556.68000  </td><td>1761.110000 </td><td> 361.75000  </td><td>1773.40000  </td><td>1350.98000  </td><td>1232.7500   </td><td>1970.630000 </td><td>2094.6100000</td><td> 912.20800  </td><td>...         </td><td>2709.4000   </td><td>1564.670000 </td><td> 340.142000 </td><td>2330.43000  </td><td>1648.03000  </td><td> 568.05000  </td><td>2390.0200000</td><td>2116.01000  </td><td>1406.45000  </td><td> 497.67000  </td></tr>
	<tr><th scope=row>31720_s_at</th><td>1937.31000  </td><td>1911.55000  </td><td>1776.920000 </td><td> 179.56700  </td><td>1622.36000  </td><td>1425.73000  </td><td>2182.1300   </td><td>2206.930000 </td><td>2349.0600000</td><td> 663.22800  </td><td>...         </td><td>1378.6300   </td><td>1988.770000 </td><td> 141.886000 </td><td>1832.11000  </td><td>1427.10000  </td><td> 826.80300  </td><td>2842.8000000</td><td> 678.94900  </td><td>1583.27000  </td><td> 345.19400  </td></tr>
	<tr><th scope=row>31721_at</th><td> 233.21400  </td><td> 146.58000  </td><td> 257.100000 </td><td> 217.95200  </td><td> 167.25200  </td><td> 166.82300  </td><td> 212.5340   </td><td>  93.718500 </td><td> 100.4540000</td><td> 288.02200  </td><td>...         </td><td> 237.1510   </td><td> 229.076000 </td><td> 266.149000 </td><td> 163.47000  </td><td> 243.45800  </td><td> 177.08200  </td><td> 188.2110000</td><td> 264.59100  </td><td> 131.58100  </td><td> 157.83300  </td></tr>
	<tr><th scope=row>31722_at</th><td>3250.14000  </td><td>4013.99000  </td><td>3288.880000 </td><td>4233.25000  </td><td>4208.54000  </td><td>3207.93000  </td><td>2586.2800   </td><td>4695.350000 </td><td>5027.3400000</td><td>3034.40000  </td><td>...         </td><td>4236.6700   </td><td>4052.490000 </td><td>5408.500000 </td><td>4295.36000  </td><td>2091.88000  </td><td>2839.70000  </td><td>3583.6700000</td><td>3558.35000  </td><td>3462.20000  </td><td>3569.07000  </td></tr>
	<tr><th scope=row>31723_at</th><td> -22.24480  </td><td> -13.42290  </td><td>  -8.507690 </td><td> -13.52920  </td><td>  -1.83895  </td><td>  -5.10412  </td><td> -28.1626   </td><td>  -8.870770 </td><td>  -5.0921000</td><td> 189.17900  </td><td>...         </td><td> -12.2979   </td><td> -20.084900 </td><td>  -0.330254 </td><td>   1.24901  </td><td>  -6.21422  </td><td> -25.39210  </td><td> -13.5936000</td><td> -10.31500  </td><td>   7.78682  </td><td> -15.66590  </td></tr>
	<tr><th scope=row>31724_at</th><td> 269.91200  </td><td> 211.98000  </td><td> 338.949000 </td><td> 212.17300  </td><td> 194.46700  </td><td> 250.14100  </td><td> 273.3520   </td><td> 220.340000 </td><td> 197.6480000</td><td> 286.52400  </td><td>...         </td><td> 277.6690   </td><td> 304.424000 </td><td> 249.207000 </td><td> 169.03600  </td><td> 272.61100  </td><td> 214.29600  </td><td> 213.8010000</td><td> 182.81300  </td><td> 172.98700  </td><td> 259.54900  </td></tr>
	<tr><th scope=row>31725_s_at</th><td>  84.60640  </td><td>  72.72370  </td><td>  66.667400 </td><td>  60.01190  </td><td>  52.42020  </td><td>  35.72430  </td><td>  88.6858   </td><td>  81.275500 </td><td>  33.0178000</td><td>  62.35560  </td><td>...         </td><td>  60.6225   </td><td>  94.669200 </td><td>  91.413600 </td><td>  51.29220  </td><td>  59.90220  </td><td>  88.95820  </td><td>  72.5324000</td><td>  51.27280  </td><td>  45.30990  </td><td>  77.81270  </td></tr>
	<tr><th scope=row>31726_at</th><td> 233.00300  </td><td> 186.69200  </td><td> 336.641000 </td><td> 338.30700  </td><td> 155.88900  </td><td> 188.31600  </td><td> 235.7160   </td><td> 132.113000 </td><td> 154.7630000</td><td> 379.11400  </td><td>...         </td><td> 364.2820   </td><td> 234.158000 </td><td> 217.254000 </td><td> 162.70900  </td><td> 424.44300  </td><td> 262.80500  </td><td> 184.5900000</td><td> 221.99700  </td><td> 220.59800  </td><td> 390.22000  </td></tr>
	<tr><th scope=row>31727_at</th><td> 248.21600  </td><td> 165.05600  </td><td> 298.400000 </td><td> 315.37600  </td><td> 170.92900  </td><td> 202.63700  </td><td> 288.5600   </td><td> 159.900000 </td><td> 215.4460000</td><td> 292.43400  </td><td>...         </td><td> 250.0750   </td><td> 280.213000 </td><td> 223.499000 </td><td> 210.07800  </td><td> 250.47600  </td><td> 197.05300  </td><td> 223.8150000</td><td> 214.56100  </td><td> 198.22100  </td><td> 296.98200  </td></tr>
	<tr><th scope=row>31728_at</th><td> 150.12700  </td><td> 214.07300  </td><td> 195.258000 </td><td> 177.60300  </td><td> 146.26800  </td><td> 156.02200  </td><td> 147.0610   </td><td>  87.984000 </td><td>  74.4776000</td><td> 197.34700  </td><td>...         </td><td> 168.2120   </td><td> 163.260000 </td><td> 125.535000 </td><td> 101.10100  </td><td> 211.56300  </td><td> 159.42500  </td><td> 144.3630000</td><td> 231.79000  </td><td> 165.06800  </td><td> 277.79700  </td></tr>
	<tr><th scope=row>31729_at</th><td>   4.59592  </td><td>   9.80107  </td><td>   0.368562 </td><td>   1.57145  </td><td>  -3.25423  </td><td> -13.13720  </td><td> -52.4786   </td><td>  -9.886800 </td><td> -30.5011000</td><td>  -4.49651  </td><td>...         </td><td> -24.2041   </td><td>  -0.523823 </td><td> -26.370100 </td><td> -44.46980  </td><td> -42.47630  </td><td>  -1.24586  </td><td>  -0.4919960</td><td>  -3.55962  </td><td> -19.36460  </td><td> -26.67930  </td></tr>
	<tr><th scope=row>31730_at</th><td> 129.86700  </td><td>  84.41120  </td><td> 116.449000 </td><td> 166.47800  </td><td>  92.33060  </td><td>  78.08650  </td><td> 182.9700   </td><td>  88.292300 </td><td>  71.1762000</td><td> 120.19600  </td><td>...         </td><td>  93.1122   </td><td>  95.161000 </td><td> 124.581000 </td><td>  60.07000  </td><td>  97.41350  </td><td> 152.39800  </td><td>  95.4303000</td><td>  98.05440  </td><td> 100.47300  </td><td> 124.43500  </td></tr>
	<tr><th scope=row>31731_at</th><td>  19.75050  </td><td>  89.00200  </td><td>  16.639000 </td><td>  29.04000  </td><td>  28.28210  </td><td>  51.56130  </td><td>  36.2152   </td><td>  38.653500 </td><td>  78.4340000</td><td>  27.30890  </td><td>...         </td><td>  16.1807   </td><td>  -5.098980 </td><td>  23.090200 </td><td>  31.33940  </td><td>  21.34330  </td><td>  20.85870  </td><td>  53.5745000</td><td>  21.02060  </td><td>  51.60150  </td><td>  23.65370  </td></tr>
	<tr><th scope=row>31732_at</th><td>  24.81900  </td><td>  26.97430  </td><td>  39.559300 </td><td>  46.95140  </td><td>  36.10260  </td><td>  33.05720  </td><td>  33.1716   </td><td>   7.802740 </td><td>   9.6353500</td><td>  59.29320  </td><td>...         </td><td>  38.0609   </td><td>  35.584800 </td><td>  27.387700 </td><td>  30.05440  </td><td>  34.67080  </td><td>  32.89640  </td><td>  23.8470000</td><td>  52.56220  </td><td>  20.08900  </td><td>  59.35990  </td></tr>
	<tr><th scope=row>31733_at</th><td>  63.57600  </td><td>  11.67840  </td><td>  55.677300 </td><td>  42.20100  </td><td>  52.90310  </td><td>  43.56170  </td><td>  58.5177   </td><td>   8.338340 </td><td>  38.4290000</td><td>  68.85230  </td><td>...         </td><td>  49.8754   </td><td>  51.419100 </td><td>  63.523500 </td><td>  26.29020  </td><td>  59.00110  </td><td>  87.76500  </td><td>  43.1828000</td><td>  46.73890  </td><td>  28.98340  </td><td>  63.83320  </td></tr>
	<tr><th scope=row>31734_at</th><td> 190.53300  </td><td> 169.96200  </td><td> 156.710000 </td><td> 211.62400  </td><td> 105.93900  </td><td> 171.99200  </td><td> 181.2790   </td><td> 164.635000 </td><td>  95.5045000</td><td> 178.68300  </td><td>...         </td><td> 118.7170   </td><td> 218.051000 </td><td>  91.589900 </td><td> 162.19000  </td><td> 191.44300  </td><td> 159.64000  </td><td> 122.6620000</td><td> 144.45000  </td><td> 102.52100  </td><td> 114.49000  </td></tr>
	<tr><th scope=row>31735_at</th><td>  26.70160  </td><td>  33.15780  </td><td>  31.711300 </td><td>  36.62170  </td><td>  26.82840  </td><td>  38.99100  </td><td>  70.0673   </td><td>  37.418100 </td><td>  35.9276000</td><td>  37.80290  </td><td>...         </td><td>  32.9751   </td><td>  43.318400 </td><td>  30.248500 </td><td>  23.65490  </td><td>  28.99500  </td><td>  58.85010  </td><td>  50.3786000</td><td>  30.42690  </td><td>  35.53090  </td><td>  39.68020  </td></tr>
	<tr><th scope=row>31736_at</th><td> 446.51200  </td><td> 271.49400  </td><td> 304.809000 </td><td> 340.97800  </td><td> 356.12700  </td><td> 279.01500  </td><td> 393.7300   </td><td> 173.935000 </td><td> 350.8060000</td><td> 285.40900  </td><td>...         </td><td> 239.8520   </td><td> 501.531000 </td><td> 414.683000 </td><td> 193.00400  </td><td> 282.50700  </td><td> 504.18900  </td><td> 318.6810000</td><td> 304.96000  </td><td> 318.73600  </td><td> 291.70900  </td></tr>
	<tr><th scope=row>31737_at</th><td>  22.46410  </td><td>  23.45890  </td><td>  39.417000 </td><td>  25.74520  </td><td>  14.40320  </td><td>  27.16600  </td><td>  35.3598   </td><td>  21.224800 </td><td>  41.0651000</td><td>  25.25010  </td><td>...         </td><td>  20.6530   </td><td>  27.493000 </td><td>  32.080500 </td><td>  15.05670  </td><td>  21.90040  </td><td>  31.65340  </td><td>  22.9250000</td><td>  14.90170  </td><td>  22.73410  </td><td>  24.27460  </td></tr>
	<tr><th scope=row>31738_at</th><td> 299.43400  </td><td> 233.13800  </td><td> 355.204000 </td><td> 314.81800  </td><td> 238.68400  </td><td> 205.69700  </td><td> 400.9550   </td><td> 218.935000 </td><td> 155.0070000</td><td> 320.90400  </td><td>...         </td><td> 237.0430   </td><td> 287.333000 </td><td> 292.537000 </td><td> 284.32800  </td><td> 223.67900  </td><td> 273.78000  </td><td> 283.9620000</td><td> 230.87800  </td><td> 216.65200  </td><td> 297.74000  </td></tr>
	<tr><th scope=row>31739_at</th><td> 253.69200  </td><td> 183.30600  </td><td> 291.385000 </td><td> 270.71900  </td><td> 212.02500  </td><td> 225.35700  </td><td> 267.0190   </td><td> 213.479000 </td><td> 147.5640000</td><td> 258.65800  </td><td>...         </td><td> 310.9600   </td><td> 262.567000 </td><td> 337.140000 </td><td> 304.22000  </td><td> 247.84300  </td><td> 202.28400  </td><td> 180.5560000</td><td> 247.59000  </td><td> 137.95900  </td><td> 287.74900  </td></tr>
</tbody>
</table>



4、check the object created previously.<br>(检查确认刚才创建的对象)


```R
class(exprs)
```


'matrix'



```R
dim(exprs)
```


<ol class=list-inline>
	<li>500</li>
	<li>26</li>
</ol>



5、read the phenotype information file in a similar way using the read.csv function.<br>
（使用read.csv函数以类似的方式读取表型信息文件）


```R
pData <- read.table(pDataLoc, row.names = 1, header = TRUE, sep = "\t")
pData
```


<table>
<thead><tr><th></th><th scope=col>gender</th><th scope=col>type</th><th scope=col>score</th></tr></thead>
<tbody>
	<tr><th scope=row>A</th><td>Female </td><td>Control</td><td>0.75   </td></tr>
	<tr><th scope=row>B</th><td>Male   </td><td>Case   </td><td>0.40   </td></tr>
	<tr><th scope=row>C</th><td>Male   </td><td>Control</td><td>0.73   </td></tr>
	<tr><th scope=row>D</th><td>Male   </td><td>Case   </td><td>0.42   </td></tr>
	<tr><th scope=row>E</th><td>Female </td><td>Case   </td><td>0.93   </td></tr>
	<tr><th scope=row>F</th><td>Male   </td><td>Control</td><td>0.22   </td></tr>
	<tr><th scope=row>G</th><td>Male   </td><td>Case   </td><td>0.96   </td></tr>
	<tr><th scope=row>H</th><td>Male   </td><td>Case   </td><td>0.79   </td></tr>
	<tr><th scope=row>I</th><td>Female </td><td>Case   </td><td>0.37   </td></tr>
	<tr><th scope=row>J</th><td>Male   </td><td>Control</td><td>0.63   </td></tr>
	<tr><th scope=row>K</th><td>Male   </td><td>Case   </td><td>0.26   </td></tr>
	<tr><th scope=row>L</th><td>Female </td><td>Control</td><td>0.36   </td></tr>
	<tr><th scope=row>M</th><td>Male   </td><td>Case   </td><td>0.41   </td></tr>
	<tr><th scope=row>N</th><td>Male   </td><td>Case   </td><td>0.80   </td></tr>
	<tr><th scope=row>O</th><td>Female </td><td>Case   </td><td>0.10   </td></tr>
	<tr><th scope=row>P</th><td>Female </td><td>Control</td><td>0.41   </td></tr>
	<tr><th scope=row>Q</th><td>Female </td><td>Case   </td><td>0.16   </td></tr>
	<tr><th scope=row>R</th><td>Male   </td><td>Control</td><td>0.72   </td></tr>
	<tr><th scope=row>S</th><td>Male   </td><td>Case   </td><td>0.17   </td></tr>
	<tr><th scope=row>T</th><td>Female </td><td>Case   </td><td>0.74   </td></tr>
	<tr><th scope=row>U</th><td>Male   </td><td>Control</td><td>0.35   </td></tr>
	<tr><th scope=row>V</th><td>Female </td><td>Control</td><td>0.77   </td></tr>
	<tr><th scope=row>W</th><td>Male   </td><td>Control</td><td>0.27   </td></tr>
	<tr><th scope=row>X</th><td>Male   </td><td>Control</td><td>0.98   </td></tr>
	<tr><th scope=row>Y</th><td>Female </td><td>Case   </td><td>0.94   </td></tr>
	<tr><th scope=row>Z</th><td>Female </td><td>Case   </td><td>0.32   </td></tr>
</tbody>
</table>




```R
pData <- new("AnnotatedDataFrame", data = pData)
pData
```


    An object of class 'AnnotatedDataFrame'
      rowNames: A B ... Z (26 total)
      varLabels: gender type score
      varMetadata: labelDescription


6、Compile the experiment information, an object of the MIAME class with slots for investigator name, lab name, and so on.<br>
(编辑实验信息，这是MIAME类的一个对象，带有调查员姓名，实验室名称等)<br>
微阵列实验最小信息集MIAME( Minimum Information About a Microarray Experiment)


```R
exData <- new("MIAME", name="ABCabc", lab="XYZ Lab", contact="abc@xyz", title="", abstract="", url="www.xyz")
exData
```


    Experiment data
      Experimenter name: ABCabc 
      Laboratory: XYZ Lab 
      Contact information: abc@xyz 
      Title:  
      URL: www.xyz 
      PMIDs:  
      No abstract available.


7、Because it's important to know the chip annotation as it is a part of the `ExpressionSet` object for this data.
Take the hgu95av2 chip annotation as an example. Create a new `ExpressionSet` object using the information compiled 
in the previous steps.<br>
(因为了解芯片注释非常重要，因为它是此数据的`ExpressionSet`对象的一部分。以hgu95av2芯片注释为例。使用前面步骤中编译的信息创建一个新的`ExpressionSet`对象。)


```R
exampleSet <- new("ExpressionSet", exprs = exprs, phenoData = pData, experimentData = exData, annotation = "hgu133a2")
exampleSet
```


    ExpressionSet (storageMode: lockedEnvironment)
    assayData: 500 features, 26 samples 
      element names: exprs 
    protocolData: none
    phenoData
      sampleNames: A B ... Z (26 total)
      varLabels: gender type score
      varMetadata: labelDescription
    featureData: none
    experimentData: use 'experimentData(object)'
    Annotation: hgu133a2 


8、Check the object


```R
 str (exampleSet)
```

    Formal class 'ExpressionSet' [package "Biobase"] with 7 slots
      ..@ experimentData   :Formal class 'MIAME' [package "Biobase"] with 13 slots
      .. .. ..@ name             : chr "ABCabc"
      .. .. ..@ lab              : chr "XYZ Lab"
      .. .. ..@ contact          : chr "abc@xyz"
      .. .. ..@ title            : chr ""
      .. .. ..@ abstract         : chr ""
      .. .. ..@ url              : chr "www.xyz"
      .. .. ..@ pubMedIds        : chr ""
      .. .. ..@ samples          : list()
      .. .. ..@ hybridizations   : list()
      .. .. ..@ normControls     : list()
      .. .. ..@ preprocessing    : list()
      .. .. ..@ other            : list()
      .. .. ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slot
      .. .. .. .. ..@ .Data:List of 2
      .. .. .. .. .. ..$ : int [1:3] 1 0 0
      .. .. .. .. .. ..$ : int [1:3] 1 1 0
      ..@ assayData        :<environment: 0x000000000bbbbe00> 
      ..@ phenoData        :Formal class 'AnnotatedDataFrame' [package "Biobase"] with 4 slots
      .. .. ..@ varMetadata      :'data.frame':	3 obs. of  1 variable:
      .. .. .. ..$ labelDescription: chr [1:3] NA NA NA
      .. .. ..@ data             :'data.frame':	26 obs. of  3 variables:
      .. .. .. ..$ gender: Factor w/ 2 levels "Female","Male": 1 2 2 2 1 2 2 2 1 2 ...
      .. .. .. ..$ type  : Factor w/ 2 levels "Case","Control": 2 1 2 1 1 2 1 1 1 2 ...
      .. .. .. ..$ score : num [1:26] 0.75 0.4 0.73 0.42 0.93 0.22 0.96 0.79 0.37 0.63 ...
      .. .. ..@ dimLabels        : chr [1:2] "sampleNames" "sampleColumns"
      .. .. ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slot
      .. .. .. .. ..@ .Data:List of 1
      .. .. .. .. .. ..$ : int [1:3] 1 1 0
      ..@ featureData      :Formal class 'AnnotatedDataFrame' [package "Biobase"] with 4 slots
      .. .. ..@ varMetadata      :'data.frame':	0 obs. of  1 variable:
      .. .. .. ..$ labelDescription: chr(0) 
      .. .. ..@ data             :'data.frame':	500 obs. of  0 variables
      .. .. ..@ dimLabels        : chr [1:2] "featureNames" "featureColumns"
      .. .. ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slot
      .. .. .. .. ..@ .Data:List of 1
      .. .. .. .. .. ..$ : int [1:3] 1 1 0
      ..@ annotation       : chr "hgu133a2"
      ..@ protocolData     :Formal class 'AnnotatedDataFrame' [package "Biobase"] with 4 slots
      .. .. ..@ varMetadata      :'data.frame':	0 obs. of  1 variable:
      .. .. .. ..$ labelDescription: chr(0) 
      .. .. ..@ data             :'data.frame':	26 obs. of  0 variables
      .. .. ..@ dimLabels        : chr [1:2] "sampleNames" "sampleColumns"
      .. .. ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slot
      .. .. .. .. ..@ .Data:List of 1
      .. .. .. .. .. ..$ : int [1:3] 1 1 0
      ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slot
      .. .. ..@ .Data:List of 4
      .. .. .. ..$ : int [1:3] 3 4 2
      .. .. .. ..$ : int [1:3] 2 38 0
      .. .. .. ..$ : int [1:3] 1 3 0
      .. .. .. ..$ : int [1:3] 1 0 0
    

9、Test the validity of the object created before continuing with the analysis.<br>
(在继续分析之前测试创建的对象的有效性)


```R
validObject(exampleSet)
```


TRUE


In this chapter, we read different information files individually using the conventional `read.csv` function in a matrix or data frame. The expression data is a matrix that contains the intensities measured, whereas the phenotypic data carries information about the conditions (for example, control or disease) of the data and samples. The experimental data simply has certain formal information, and it is not obligatory to fill it in. As the order is very important for the final eSet, we check the validity of the created object. The annotation chip used is because the built-in data for the package actually comes from the hgu133a2Affymetrix chip. For example, if the sample names in the expression data and phenotypic data are different, the function will return the object as invalid. These individual objects are then assembled into ExpressionSet by creating a new object. Each component of ExpressionSet has its own role. The exprs object is the expression data, the phenotypic data summarizes information about the samples (for example, the sex, age, and treatment status—referred to as covariates), and the annotated package provides basic data manipulation tools for the metadata packages. This can be done with any platform, be it Affymetrix or Illumina.<br>
	在本章中，我们使用传统的`read.csv`函数在矩阵或数据框中分别读取不同的信息文件。表达数据一种所测量的强度值的矩阵，而表型数据携带关于数据和样本的条件（例如，对照或疾病）的信息。实验数据只是具有一定的形式信息，并非必须填写。由于顺序对于最终的eSet对象非常重要，我们需要检查所创建对象的有效性。因为所采用的内置数据实际上来自hgu133a2Affymetrix芯片，所以采用示例中使用的注释芯片。例如，如果表达式数据和表型数据中的样本名称不同，则该函数将返回该对象为无效。然后通过创建一个新对象将这些单独的对象组装到ExpressionSet中。ExpressionSet的每个组件都有其自己的角色。exprs对象是表达数据，表型数据总结了关于样本的信息（例如，性别，年龄和治疗状态 - 称为协变量），并且注释的包为元数据包提供了基本的数据处理工具。这可以通过任何平台完成，无论是Affymetrix还是Illumina。
