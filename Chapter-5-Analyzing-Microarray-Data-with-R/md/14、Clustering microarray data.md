
# Clustering microarray data
   Clustering is about aggregating similar genes together in a group (called cluster) and away from other such groups. When genes get clustered together (falling in the same group/cluster), it means they follow a similar pattern based on the expression data under the given conditions. This section presents the widely used concept of hierarchical clustering in gene
expression analysis.<br>
聚类是关于将一组（称为聚类）中的相似基因聚集在一起并远离其他这样的聚类。当基因聚集在一起时（落在同一组/簇中），这意味着它们在给定条件下基于表达数据遵循类似的模式。本节内容介绍了基因表达分析中广泛使用的等级聚类的概念。<br>
Requires:<br>
1、The normalized breast cancer data from the earlier part. We use only part of it, the top 1500 genes—for a faster computation.<br>
2、`EMA` package.

1、Create a dataset for clustering purposes from the leukemia data again. Use only the first 100 data instances for demonstration purposes. 


```R
library(affy)
library(affydata)
library(hgu133a2cdf)
mydata <- ReadAffy(celfile.path= "D:/Try-practice/Chapter 5/GSE24460_RAW/")
# mydata <- ReadAffy(filenames="D:/Try-practice/Chapter 5/GSE24460_RAW/GSM602658_MCF71.CEL")
library(Biobase)
DIR <- system.file("extdata", package="Biobase")
exprsLoc <- file.path(DIR, "exprsData.txt") 
pDataLoc <- file.path(DIR, "pData.txt")
DIR
exprsLoc
pDataLoc
exprs <- as.matrix(read.csv(exprsLoc, header = TRUE, sep = "\t", row.names = 1, as.is = TRUE))
pData <- read.table(pDataLoc, row.names = 1, header = TRUE, sep = "\t")
pData <- new("AnnotatedDataFrame", data = pData)
exData <- new("MIAME", name="ABCabc", lab="XYZ Lab", contact="abc@xyz", title="", abstract="", url="www.xyz")
eset <- new("ExpressionSet", exprs = exprs, phenoData = pData, experimentData = exData, annotation = "hgu133a2")
```


'C:/Users/Master1/Documents/R/win-library/3.4/Biobase/extdata'



'C:/Users/Master1/Documents/R/win-library/3.4/Biobase/extdata/exprsData.txt'



'C:/Users/Master1/Documents/R/win-library/3.4/Biobase/extdata/pData.txt'



```R
c.data <- exprs(eset[1:100,])
c.data
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
	<tr><th scope=row>31310_at</th><td>  17.15290  </td><td>  27.17220  </td><td>  -0.338279 </td><td> -14.10610  </td><td>  23.02630  </td><td>  17.929400 </td><td>  35.720900 </td><td>  -3.43857  </td><td>   7.00862  </td><td>  14.68320  </td><td>...         </td><td>  -2.7933000</td><td>  -0.905975 </td><td>  12.07360  </td><td>  16.0306000</td><td>   9.011890 </td><td>   9.60877  </td><td>   7.265800 </td><td> -12.0555000</td><td>   3.83654  </td><td>   8.301990 </td></tr>
	<tr><th scope=row>31311_at</th><td> 216.93200  </td><td> 162.56300  </td><td> 166.482000 </td><td> 248.73500  </td><td> 183.23700  </td><td> 129.778000 </td><td> 270.832000 </td><td> 126.72000  </td><td> 159.53900  </td><td> 173.03300  </td><td>...         </td><td> 143.7310000</td><td> 233.632000 </td><td> 249.66700  </td><td> 175.3850000</td><td> 195.923000 </td><td> 242.32700  </td><td> 132.225000 </td><td> 204.7950000</td><td> 188.77500  </td><td> 201.313000 </td></tr>
	<tr><th scope=row>31312_at</th><td> 190.17800  </td><td> 120.16100  </td><td> 169.273000 </td><td> 173.70200  </td><td> 102.38500  </td><td>  93.438600 </td><td> 183.729000 </td><td> 113.82100  </td><td>  77.68790  </td><td> 173.09200  </td><td>...         </td><td> 151.5810000</td><td> 204.266000 </td><td> 152.51100  </td><td> 123.6800000</td><td> 148.046000 </td><td> 149.53400  </td><td> 113.424000 </td><td> 136.2080000</td><td>  87.39970  </td><td> 154.514000 </td></tr>
	<tr><th scope=row>31313_at</th><td> 113.47600  </td><td>  84.91850  </td><td> 177.832000 </td><td> 117.49300  </td><td>  88.02680  </td><td>  90.219300 </td><td> 117.556000 </td><td> 124.49000  </td><td>  71.38880  </td><td> 236.38800  </td><td>...         </td><td> 122.3920000</td><td> 172.993000 </td><td> 106.06400  </td><td>  85.9635000</td><td> 137.080000 </td><td> 153.57600  </td><td>  90.781200 </td><td>  67.8371000</td><td> 136.29700  </td><td> 120.003000 </td></tr>
	<tr><th scope=row>31314_at</th><td>  50.79450  </td><td>  53.81800  </td><td>  58.241200 </td><td>  59.73410  </td><td>  75.85510  </td><td>  65.724700 </td><td>  93.811700 </td><td>  48.99770  </td><td>  44.65160  </td><td>  95.00140  </td><td>...         </td><td>  64.1095000</td><td> 103.938000 </td><td>  91.26490  </td><td>  69.9593000</td><td>  88.721300 </td><td> 132.63300  </td><td>  49.423600 </td><td>  47.1633000</td><td> 100.92300  </td><td> 128.271000 </td></tr>
	<tr><th scope=row>31315_at</th><td>1139.56000  </td><td>5154.31000  </td><td> 600.494000 </td><td> 658.42100  </td><td>3378.38000  </td><td> 265.693000 </td><td>2347.970000 </td><td>1474.59000  </td><td> 966.33900  </td><td>1603.25000  </td><td>...         </td><td> 241.9610000</td><td>1614.460000 </td><td> 163.92900  </td><td> 838.7970000</td><td>1731.100000 </td><td>3006.34000  </td><td> 194.233000 </td><td> 130.2770000</td><td>3595.46000  </td><td>3949.660000 </td></tr>
	<tr><th scope=row>31316_at</th><td> -18.28700  </td><td>   3.74538  </td><td>  -3.035990 </td><td>  -5.98093  </td><td>   3.59555  </td><td>   2.956180 </td><td>  -6.845870 </td><td>   6.59068  </td><td>  18.00900  </td><td>  -4.79112  </td><td>...         </td><td>   9.3319800</td><td>  20.499700 </td><td> -15.27830  </td><td>  42.3079000</td><td>   5.709700 </td><td>   2.60968  </td><td>  -2.227830 </td><td>  -3.0579800</td><td>   7.79922  </td><td>  -3.908320 </td></tr>
	<tr><th scope=row>31317_r_at</th><td>  17.27610  </td><td>  36.00950  </td><td>   0.533814 </td><td>  11.40780  </td><td>  65.02420  </td><td>   1.607010 </td><td>  -3.516010 </td><td>  26.29300  </td><td>  14.35940  </td><td>  14.39310  </td><td>...         </td><td>   2.8432100</td><td>  12.068400 </td><td>  22.20970  </td><td>   9.2864500</td><td>  65.745700 </td><td> 109.94300  </td><td>  12.381200 </td><td>  15.5240000</td><td> 187.65600  </td><td>   1.375560 </td></tr>
	<tr><th scope=row>31318_at</th><td>  17.70390  </td><td>  15.23570  </td><td>  10.275300 </td><td>  12.96480  </td><td>  12.96000  </td><td>  -0.761071 </td><td>  16.788300 </td><td>  15.17440  </td><td>  16.26390  </td><td>  12.91680  </td><td>...         </td><td>  -0.0813966</td><td>  15.249900 </td><td>  32.28400  </td><td>   7.8017900</td><td>  11.692400 </td><td>  25.94010  </td><td>   2.208690 </td><td>   0.3940580</td><td>  16.36970  </td><td>   0.115654 </td></tr>
	<tr><th scope=row>31319_at</th><td>  15.86550  </td><td>  23.91830  </td><td>  22.797800 </td><td>  24.81950  </td><td>  99.21890  </td><td>  25.180500 </td><td>  -6.772620 </td><td>  31.09880  </td><td>  32.86500  </td><td>  45.52500  </td><td>...         </td><td>  13.0512000</td><td>  37.118200 </td><td>  19.45690  </td><td>  41.7727000</td><td> 164.186000 </td><td>  70.51820  </td><td>  11.569500 </td><td>  -0.5751810</td><td> 141.18300  </td><td> 186.293000 </td></tr>
	<tr><th scope=row>31320_at</th><td> 123.62100  </td><td>  95.35860  </td><td>  94.446800 </td><td> 124.79400  </td><td>  89.24970  </td><td> 113.248000 </td><td> 143.695000 </td><td> 164.83900  </td><td> 135.09800  </td><td> 123.52100  </td><td>...         </td><td> 137.2740000</td><td> 133.572000 </td><td>  90.36660  </td><td> 147.4480000</td><td>  62.757500 </td><td> 118.33600  </td><td>  99.764700 </td><td> 116.9370000</td><td>  79.45960  </td><td> 120.317000 </td></tr>
	<tr><th scope=row>31321_at</th><td> 174.66100  </td><td> 103.08000  </td><td> 177.465000 </td><td> 174.55900  </td><td> 121.01000  </td><td> 128.501000 </td><td> 173.604000 </td><td> 105.42700  </td><td> 124.68800  </td><td> 160.28400  </td><td>...         </td><td> 152.3120000</td><td> 145.011000 </td><td> 154.36300  </td><td> 124.6520000</td><td> 170.711000 </td><td> 173.69400  </td><td> 128.898000 </td><td> 137.4110000</td><td> 112.34400  </td><td> 131.966000 </td></tr>
	<tr><th scope=row>31322_at</th><td>  25.41690  </td><td>  42.35770  </td><td>  28.938600 </td><td>  19.36230  </td><td>  22.20280  </td><td>  13.287200 </td><td>  25.716500 </td><td>  15.35150  </td><td>  10.43240  </td><td>  63.72350  </td><td>...         </td><td>  17.9966000</td><td>  29.006900 </td><td>  -2.09740  </td><td>  14.7371000</td><td>  27.103000 </td><td>  17.66400  </td><td>  19.561500 </td><td>  15.8249000</td><td>  20.76490  </td><td>  33.655400 </td></tr>
	<tr><th scope=row>31323_r_at</th><td>  17.03390  </td><td>  23.47410  </td><td>  18.933200 </td><td>  26.61890  </td><td>  26.81390  </td><td>  33.094200 </td><td>  22.559300 </td><td>   7.64309  </td><td>-158.62400  </td><td>  52.24130  </td><td>...         </td><td>  24.0255000</td><td>  20.752600 </td><td>  28.54180  </td><td>  31.7192000</td><td>  19.026500 </td><td>   5.00616  </td><td>  35.292000 </td><td>  27.7821000</td><td>  16.97820  </td><td>   4.972290 </td></tr>
	<tr><th scope=row>31324_at</th><td> 102.59100  </td><td>  92.34860  </td><td> 120.461000 </td><td> 113.72300  </td><td>  90.03440  </td><td> 143.403000 </td><td>  95.818900 </td><td>  75.70230  </td><td>  77.41820  </td><td> 115.72500  </td><td>...         </td><td>  56.6119000</td><td>  72.851600 </td><td>  56.72670  </td><td>  37.8179000</td><td>  77.292200 </td><td>  93.25110  </td><td>  94.533200 </td><td>  56.4698000</td><td>  84.32640  </td><td>  30.532700 </td></tr>
	<tr><th scope=row>31325_at</th><td>  67.86060  </td><td>  24.43680  </td><td>  84.328500 </td><td>  84.25960  </td><td>  74.55310  </td><td>  62.906100 </td><td> 100.242000 </td><td>  58.94120  </td><td>  52.93510  </td><td>  76.43690  </td><td>...         </td><td>  19.8654000</td><td> 102.421000 </td><td>  66.85600  </td><td>  42.2796000</td><td>  43.557400 </td><td>  73.81480  </td><td>  67.459100 </td><td>  72.6469000</td><td>  44.89180  </td><td>  38.467700 </td></tr>
	<tr><th scope=row>31326_at</th><td> 754.36100  </td><td> 777.10500  </td><td>1208.680000 </td><td> 730.31400  </td><td> 679.41000  </td><td> 472.386000 </td><td> 875.514000 </td><td> 426.06900  </td><td> 541.81200  </td><td> 647.94600  </td><td>...         </td><td> 602.8920000</td><td> 734.952000 </td><td> 807.89000  </td><td> 488.8850000</td><td> 689.458000 </td><td> 848.67200  </td><td> 554.534000 </td><td> 793.0190000</td><td> 629.05500  </td><td> 695.874000 </td></tr>
	<tr><th scope=row>31327_at</th><td>   1.86526  </td><td>   5.41691  </td><td>  11.763100 </td><td>   8.62522  </td><td>  13.40060  </td><td>   9.466650 </td><td>   5.644840 </td><td>  22.78630  </td><td>   9.68571  </td><td>  10.82870  </td><td>...         </td><td>  18.9863000</td><td>  10.860800 </td><td>  16.85840  </td><td> -46.2152000</td><td>  17.202200 </td><td>  14.98240  </td><td>   7.661780 </td><td>  -3.2077700</td><td>  22.39660  </td><td>  12.178000 </td></tr>
	<tr><th scope=row>31328_at</th><td> 121.89200  </td><td>  92.03100  </td><td> 201.403000 </td><td> 134.99300  </td><td> 117.13400  </td><td> 146.648000 </td><td> 159.059000 </td><td>  93.86750  </td><td> 110.69400  </td><td> 196.50400  </td><td>...         </td><td> 128.2560000</td><td> 148.424000 </td><td> 107.98300  </td><td> 132.0580000</td><td> 158.953000 </td><td> 104.30300  </td><td> 125.216000 </td><td> 166.7850000</td><td>  99.98540  </td><td> 130.613000 </td></tr>
	<tr><th scope=row>31329_at</th><td>  14.55860  </td><td>  34.33590  </td><td>  16.304600 </td><td>  18.44850  </td><td>  13.41280  </td><td>  20.175100 </td><td>  12.877700 </td><td>  36.69980  </td><td>  21.35950  </td><td>  30.84720  </td><td>...         </td><td>   5.1899400</td><td>  14.288700 </td><td>  20.47140  </td><td>   5.2728700</td><td>  19.079600 </td><td>  18.50450  </td><td>  14.458100 </td><td>  20.0154000</td><td>  26.02740  </td><td>  15.966400 </td></tr>
	<tr><th scope=row>31330_at</th><td>3175.57000  </td><td>3548.02000  </td><td>1820.510000 </td><td>2612.13000  </td><td>3201.35000  </td><td>2055.840000 </td><td>1445.570000 </td><td>3752.44000  </td><td>3766.64000  </td><td>2082.26000  </td><td>...         </td><td>3225.7800000</td><td>1699.090000 </td><td>3460.35000  </td><td>4440.3300000</td><td>1570.760000 </td><td>2350.16000  </td><td>2937.070000 </td><td>3257.5700000</td><td>2532.94000  </td><td>5821.410000 </td></tr>
	<tr><th scope=row>31331_at</th><td>  40.80680  </td><td>   3.23610  </td><td>  18.700900 </td><td>  28.57110  </td><td>  23.36410  </td><td>  10.412800 </td><td>  56.507500 </td><td>  15.53750  </td><td>  11.79480  </td><td>  24.98310  </td><td>...         </td><td>  20.0576000</td><td>  40.093800 </td><td>  47.61550  </td><td>  -0.0718337</td><td>  11.124300 </td><td>  52.17060  </td><td>  16.209700 </td><td>  27.9982000</td><td>  18.32180  </td><td>  16.724800 </td></tr>
	<tr><th scope=row>31332_at</th><td>  18.12730  </td><td>   9.65582  </td><td>   8.060630 </td><td>   5.85333  </td><td>  -5.30554  </td><td>   1.514980 </td><td>   6.677060 </td><td> -37.26260  </td><td> -41.51810  </td><td>   7.45908  </td><td>...         </td><td> -22.8684000</td><td>   4.555100 </td><td>  -5.06083  </td><td> -57.6240000</td><td> -47.842500 </td><td>   4.98633  </td><td> -13.329700 </td><td>  -0.0363981</td><td> -23.32180  </td><td>   6.459900 </td></tr>
	<tr><th scope=row>31333_at</th><td>  22.59120  </td><td>   9.13811  </td><td>  17.792700 </td><td>  19.19990  </td><td>  20.12530  </td><td>   7.614780 </td><td>  20.367200 </td><td>  40.97200  </td><td>  32.13780  </td><td>  29.29040  </td><td>...         </td><td>  24.0201000</td><td>  22.204200 </td><td>  36.69950  </td><td>   4.7169500</td><td>  10.652400 </td><td>  24.29240  </td><td>   8.896130 </td><td>  25.4436000</td><td>  33.75870  </td><td>  24.294900 </td></tr>
	<tr><th scope=row>31334_at</th><td>  14.40760  </td><td>  -1.38469  </td><td>   4.914040 </td><td>  11.13250  </td><td>   9.54839  </td><td>  10.127300 </td><td>  16.617500 </td><td>   2.83778  </td><td>  18.39910  </td><td>   3.96332  </td><td>...         </td><td>  23.2550000</td><td>  20.286200 </td><td>  25.77820  </td><td>   5.4642600</td><td>  -0.824531 </td><td>  15.92690  </td><td>   0.268704 </td><td>  10.4671000</td><td>  11.06420  </td><td>   9.854540 </td></tr>
	<tr><th scope=row>31335_at</th><td>  64.55550  </td><td>  56.56760  </td><td>  64.335100 </td><td>  56.77370  </td><td>  38.77180  </td><td>  50.104600 </td><td>  71.037800 </td><td>  29.80970  </td><td>  44.21400  </td><td>  63.39640  </td><td>...         </td><td>  52.3609000</td><td>  66.578500 </td><td>  67.98570  </td><td>  39.1079000</td><td>  67.961300 </td><td>  63.37910  </td><td>  53.331200 </td><td>  38.5609000</td><td>  56.98480  </td><td>  57.845800 </td></tr>
	<tr><th scope=row>31336_at</th><td>  17.84680  </td><td>   1.12837  </td><td>   9.332220 </td><td>  -1.59376  </td><td>  -4.62199  </td><td>  17.191900 </td><td>   0.121618 </td><td>   7.85899  </td><td>  25.68680  </td><td>   5.05391  </td><td>...         </td><td>   9.7868100</td><td>   8.129380 </td><td>  -7.19370  </td><td>  12.6737000</td><td>  18.069700 </td><td>   4.59506  </td><td>  -3.185410 </td><td>   3.5801800</td><td>  14.23930  </td><td>   5.815030 </td></tr>
	<tr><th scope=row>31337_at</th><td> 412.28300  </td><td> 440.35000  </td><td> 495.373000 </td><td> 501.23600  </td><td> 348.54000  </td><td> 416.904000 </td><td> 770.457000 </td><td> 383.39300  </td><td> 332.74700  </td><td> 523.43400  </td><td>...         </td><td> 378.5320000</td><td> 438.409000 </td><td> 510.06600  </td><td> 321.9940000</td><td> 456.784000 </td><td> 336.98200  </td><td> 400.768000 </td><td> 404.8340000</td><td> 293.74800  </td><td> 386.676000 </td></tr>
	<tr><th scope=row>31338_at</th><td>   8.93417  </td><td>   4.28296  </td><td>  11.343400 </td><td>  13.89760  </td><td>  13.88210  </td><td>  16.398600 </td><td>  16.053600 </td><td>   8.59187  </td><td>  20.88130  </td><td>   3.81431  </td><td>...         </td><td>   7.8640900</td><td>   9.160170 </td><td>   8.34297  </td><td>   9.9734400</td><td>   4.193710 </td><td>  17.33300  </td><td>   7.918640 </td><td>  -4.7036100</td><td>  15.96070  </td><td>   1.522470 </td></tr>
	<tr><th scope=row>31339_at</th><td> -28.99850  </td><td> -30.05320  </td><td> -26.972700 </td><td> -23.00420  </td><td> -18.31410  </td><td> 316.922000 </td><td> -21.241000 </td><td> -14.67470  </td><td> -22.34640  </td><td> -26.95820  </td><td>...         </td><td> -10.9199000</td><td> -26.021700 </td><td> -22.64670  </td><td> -17.4564000</td><td> -26.200100 </td><td> -22.03100  </td><td> -22.027600 </td><td> -11.9295000</td><td> -14.14710  </td><td> -33.954000 </td></tr>
</tbody>
</table>



2、Use `EMA` package to do an array clustering. (使用`EMA`包做一个数组集群.)


```R
library(EMA)
```

    Warning message:
    "package 'EMA' was built under R version 3.4.4"
    
    ################################################################################
    
    Easy Microarray Analysis
    
    EMA stable version
    
    Current release : v1.4.4 - march 2014
    
    ################################################################################
    
    

3、Simply use the `c.data` object with clustering from the `EMA` library to perform the clustering of arrays.<br>
仅仅使用`c.data`对象和来自`EMA`库的聚类来执行数组的聚类.


```R
c.array <- clustering(data=c.data, metric="pearson", method="ward")
c.array
```


    Call:	 agnes(x = DIS, diss = TRUE, method = method, keep.diss = TRUE) 
    Agglomerative coefficient:  0.9255212 
    Order of objects:
     [1] A C M D J G O R U V Q T S X B E P Y F W I K H L N Z
    Height (summary):
        Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
    0.001988 0.004644 0.007296 0.023240 0.020449 0.149199 
    
    Available components:
    [1] "order"     "height"    "ac"        "merge"     "diss"      "call"     
    [7] "method"    "order.lab"


4、Create the dendrogram plot for the cluster by plotting the clusters.通过绘制群集来为群集创建树状图.


```R
plot(c.array)
```


[plot1](https://github.com/Chengshu21/Chapter-5-Analyzing-Microarray-Data-with-R/blob/master/md/pic/output-1Clustering%20microarray%20data.png)



[plot2](https://github.com/Chengshu21/Chapter-5-Analyzing-Microarray-Data-with-R/blob/master/md/pic/output-2Clustering%20microarray%20data.png)


5、Simply transpose the data matrix and use it as input for the data argument in the clustering function, and define the similarity metric and clustering method to cluster the gene.<br>
简单地转置数据矩阵并将其用作聚类函数中数据参数的输入，并定义相似性度量和聚类方法以对基因进行聚类.


```R
c.gene <- clustering(data=t(c.data), metric="pearsonabs", method="ward")
c.gene
```


    Call:	 agnes(x = DIS, diss = TRUE, method = method, keep.diss = TRUE) 
    Agglomerative coefficient:  0.8726019 
    Order of objects:
      [1] AFFX-MurIL2_at              AFFX-YEL024w/RIP1_at       
      [3] AFFX-CreX-5_at              31312_at                   
      [5] AFFX-HSAC07/X00351_5_st     AFFX-YEL021w/URA3_at       
      [7] AFFX-LysX-5_at              AFFX-LysX-M_at             
      [9] AFFX-MurIL4_at              31321_at                   
     [11] 31328_at                    AFFX-CreX-3_st             
     [13] 31313_at                    AFFX-TrpnX-3_at            
     [15] 31335_at                    AFFX-PheX-M_at             
     [17] AFFX-HSAC07/X00351_3_at     31314_at                   
     [19] AFFX-BioB-5_at              AFFX-DapX-3_at             
     [21] AFFX-LysX-3_at              31331_at                   
     [23] AFFX-ThrX-M_at              31307_at                   
     [25] 31311_at                    AFFX-BioB-M_at             
     [27] AFFX-BioC-5_at              AFFX-BioB-3_st             
     [29] AFFX-CreX-5_st              AFFX-BioDn-3_st            
     [31] AFFX-DapX-5_at              AFFX-DapX-M_at             
     [33] AFFX-TrpnX-5_at             AFFX-BioB-3_at             
     [35] AFFX-BioB-5_st              AFFX-BioC-5_st             
     [37] 31318_at                    AFFX-ThrX-3_at             
     [39] AFFX-HUMRGE/M10098_3_at     31327_at                   
     [41] 31333_at                    AFFX-MurIL10_at            
     [43] AFFX-BioDn-5_st             AFFX-CreX-3_at             
     [45] AFFX-HUMRGE/M10098_5_at     31326_at                   
     [47] AFFX-BioDn-5_at             AFFX-HUMRGE/M10098_M_at    
     [49] 31323_r_at                  31336_at                   
     [51] AFFX-BioDn-3_at             AFFX-M27830_3_at           
     [53] AFFX-ThrX-5_at              AFFX-YEL018w/_at           
     [55] AFFX-HSAC07/X00351_M_st     31325_at                   
     [57] AFFX-HUMGAPDH/M33197_5_st   AFFX-HUMGAPDH/M33197_3_st  
     [59] 31316_at                    31332_at                   
     [61] 31337_at                    AFFX-HSAC07/X00351_3_st    
     [63] 31308_at                    31324_at                   
     [65] 31339_at                    AFFX-MurFAS_at             
     [67] AFFX-HUMTFRR/M11507_5_at    AFFX-HUMTFRR/M11507_M_at   
     [69] AFFX-HUMTFRR/M11507_3_at    AFFX-BioC-3_at             
     [71] AFFX-BioC-3_st              AFFX-M27830_5_at           
     [73] AFFX-M27830_M_at            AFFX-YEL002c/WBP1_at       
     [75] 31309_r_at                  AFFX-HSAC07/X00351_5_at    
     [77] AFFX-HSAC07/X00351_M_at     AFFX-BioB-M_st             
     [79] 31329_at                    AFFX-PheX-5_at             
     [81] AFFX-TrpnX-M_at             31310_at                   
     [83] AFFX-hum_alu_at             AFFX-PheX-3_at             
     [85] AFFX-HUMGAPDH/M33197_5_at   AFFX-HUMGAPDH/M33197_M_at  
     [87] AFFX-HUMGAPDH/M33197_3_at   31330_at                   
     [89] AFFX-HUMGAPDH/M33197_M_st   31338_at                   
     [91] AFFX-HUMISGF3A/M97935_5_at  AFFX-HUMISGF3A/M97935_MA_at
     [93] AFFX-HUMISGF3A/M97935_MB_at AFFX-HUMISGF3A/M97935_3_at 
     [95] 31315_at                    31322_at                   
     [97] 31334_at                    31317_r_at                 
     [99] 31319_at                    31320_at                   
    Height (summary):
       Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    0.03771 0.26474 0.46398 0.58815 0.75655 2.52614 
    
    Available components:
    [1] "order"     "height"    "ac"        "merge"     "diss"      "call"     
    [7] "method"    "order.lab"


6、Plot the results (note that for readability issues, the following screenshot shows the results for only 100 genes).<br>
根据数据结果绘制图形.


```R
plot(c.gene)
```


[plot3](https://github.com/Chengshu21/Chapter-5-Analyzing-Microarray-Data-with-R/blob/master/md/pic/output-3Clustering%20microarray%20data.png)



[plot4(https://github.com/Chengshu21/Chapter-5-Analyzing-Microarray-Data-with-R/blob/master/md/pic/output-4Clustering%20microarray%20data.png)

