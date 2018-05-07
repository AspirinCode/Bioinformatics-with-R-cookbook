
# More visualizations for gene expression data
   In this part, there will be some interesting and useful visualizations for expression data. The visualizations or plots shown are heatmaps, Venn diagrams, and volcano plots.<br>
   Followed by the `Clustering microarray data` section.


```R
library(affy)
library(affydata)
library(hgu133a2cdf)
library(EMA)
# mydata <- ReadAffy(celfile.path= "D:/Try-practice/Chapter 5/GSE24460_RAW/")
mydata <- ReadAffy(filenames="D:/Try-practice/Chapter 5/GSE24460_RAW/GSM602658_MCF71.CEL")
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
c.data <- exprs(eset[1:50,])
c.data
c.array <- clustering(data=c.data, metric="pearson", method="ward")
c.array
c.gene <- clustering(data=t(c.data), metric="pearsonabs", method="ward")
c.gene
```


'C:/Users/Master1/Documents/R/win-library/3.4/Biobase/extdata'



'C:/Users/Master1/Documents/R/win-library/3.4/Biobase/extdata/exprsData.txt'



'C:/Users/Master1/Documents/R/win-library/3.4/Biobase/extdata/pData.txt'



<table>
<thead><tr><th></th><th scope=col>A</th><th scope=col>B</th><th scope=col>C</th><th scope=col>D</th><th scope=col>E</th><th scope=col>F</th><th scope=col>G</th><th scope=col>H</th><th scope=col>I</th><th scope=col>J</th><th scope=col>...</th><th scope=col>Q</th><th scope=col>R</th><th scope=col>S</th><th scope=col>T</th><th scope=col>U</th><th scope=col>V</th><th scope=col>W</th><th scope=col>X</th><th scope=col>Y</th><th scope=col>Z</th></tr></thead>
<tbody>
	<tr><th scope=row>AFFX-MurIL2_at</th><td> 192.7420  </td><td>  85.75330 </td><td> 176.7570  </td><td> 135.57500 </td><td>  64.49390 </td><td>  76.35690 </td><td>  160.50500</td><td>  65.96310 </td><td>  56.90390 </td><td> 135.608000</td><td>...        </td><td> 179.84500 </td><td> 152.46700 </td><td> 180.83400 </td><td>  85.41460 </td><td> 157.98900 </td><td>  146.8000 </td><td>  93.8829  </td><td> 103.85500 </td><td>  64.43400 </td><td> 175.61500 </td></tr>
	<tr><th scope=row>AFFX-MurIL10_at</th><td>  97.1370  </td><td> 126.19600 </td><td>  77.9216  </td><td>  93.37130 </td><td>  24.39860 </td><td>  85.50880 </td><td>   98.90860</td><td>  81.69320 </td><td>  97.80150 </td><td>  90.483800</td><td>...        </td><td>  87.68060 </td><td> 108.03200 </td><td> 134.26300 </td><td>  91.40310 </td><td>  -8.68811 </td><td>   85.0212 </td><td>  79.2998  </td><td>  71.65520 </td><td>  64.23690 </td><td>  78.70680 </td></tr>
	<tr><th scope=row>AFFX-MurIL4_at</th><td>  45.8192  </td><td>   8.83135 </td><td>  33.0632  </td><td>  28.70720 </td><td>   5.94492 </td><td>  28.29250 </td><td>   30.96940</td><td>  14.79230 </td><td>  14.23990 </td><td>  34.487400</td><td>...        </td><td>  32.79110 </td><td>  33.52920 </td><td>  19.81720 </td><td>  20.41900 </td><td>  26.87200 </td><td>   31.1488 </td><td>  22.3420  </td><td>  19.01350 </td><td>  12.16860 </td><td>  17.37800 </td></tr>
	<tr><th scope=row>AFFX-MurFAS_at</th><td>  22.5445  </td><td>   3.60093 </td><td>  14.6883  </td><td>  12.33970 </td><td>  36.86630 </td><td>  11.25680 </td><td>   23.00340</td><td>  16.21340 </td><td>  12.03750 </td><td>   4.549780</td><td>...        </td><td>  15.94880 </td><td>  14.67530 </td><td>  -7.91911 </td><td>  12.88750 </td><td>  11.91860 </td><td>   12.8324 </td><td>  11.1390  </td><td>   7.55564 </td><td>  19.98490 </td><td>   8.96849 </td></tr>
	<tr><th scope=row>AFFX-BioB-5_at</th><td>  96.7875  </td><td>  30.43800 </td><td>  46.1271  </td><td>  70.93190 </td><td>  56.17440 </td><td>  42.67560 </td><td>   86.51560</td><td>  30.79270 </td><td>  19.71830 </td><td>  46.352000</td><td>...        </td><td>  58.62390 </td><td> 114.06200 </td><td>  93.44020 </td><td>  22.51680 </td><td>  48.64620 </td><td>   90.2215 </td><td>  42.0053  </td><td>  57.57380 </td><td>  44.82160 </td><td>  61.70440 </td></tr>
	<tr><th scope=row>AFFX-BioB-M_at</th><td>  89.0730  </td><td>  25.84610 </td><td>  57.2033  </td><td>  69.97660 </td><td>  49.58220 </td><td>  26.12620 </td><td>   75.00830</td><td>  42.33520 </td><td>  41.12070 </td><td>  91.530700</td><td>...        </td><td>  58.13310 </td><td> 104.12200 </td><td> 115.83100 </td><td>  58.12240 </td><td>  73.42210 </td><td>   64.6066 </td><td>  40.3068  </td><td>  41.82090 </td><td>  46.10870 </td><td>  49.41220 </td></tr>
	<tr><th scope=row>AFFX-BioB-3_at</th><td> 265.9640  </td><td> 181.08000 </td><td> 164.9260  </td><td> 161.46900 </td><td> 236.97600 </td><td> 156.80300 </td><td>  211.25700</td><td> 235.99400 </td><td> 175.64000 </td><td> 229.671000</td><td>...        </td><td> 192.22100 </td><td> 305.56700 </td><td> 300.68900 </td><td> 146.08100 </td><td> 142.91300 </td><td>  187.1320 </td><td> 170.5830  </td><td> 133.27900 </td><td> 187.40700 </td><td> 144.78400 </td></tr>
	<tr><th scope=row>AFFX-BioC-5_at</th><td> 110.1360  </td><td>  57.28890 </td><td>  67.3980  </td><td>  77.22070 </td><td>  41.34880 </td><td>  37.97800 </td><td>  110.55100</td><td>  47.76900 </td><td>  24.78750 </td><td>  66.730200</td><td>...        </td><td>  53.27110 </td><td> 107.23700 </td><td> 119.66600 </td><td>  24.06540 </td><td>  98.84250 </td><td>   92.0846 </td><td>  53.3866  </td><td>  52.01640 </td><td>  65.91540 </td><td>  75.00430 </td></tr>
	<tr><th scope=row>AFFX-BioC-3_at</th><td>  43.0794  </td><td>  16.80060 </td><td>  37.6002  </td><td>  46.52720 </td><td>  22.24750 </td><td>  61.64010 </td><td>   33.66230</td><td>  31.44230 </td><td>  23.10080 </td><td>  39.741900</td><td>...        </td><td>  57.50780 </td><td>  41.13370 </td><td>  79.98290 </td><td>  23.49530 </td><td>  51.56090 </td><td>   48.1247 </td><td>  31.8358  </td><td>  29.92640 </td><td>  37.86110 </td><td>  60.47720 </td></tr>
	<tr><th scope=row>AFFX-BioDn-5_at</th><td>  10.9187  </td><td>  16.17890 </td><td>  10.1495  </td><td>   9.73639 </td><td>  16.90280 </td><td>   5.33328 </td><td>   25.11820</td><td>  38.75760 </td><td>  31.40410 </td><td>   0.398779</td><td>...        </td><td>  21.50910 </td><td>   3.10536 </td><td>   5.95347 </td><td>   5.66012 </td><td>  52.93380 </td><td>   15.7267 </td><td>  15.2116  </td><td>  -5.35282 </td><td>  13.18840 </td><td>  10.03850 </td></tr>
	<tr><th scope=row>AFFX-BioDn-3_at</th><td> 751.2270  </td><td> 515.00400 </td><td> 622.9010  </td><td> 669.85900 </td><td> 414.16500 </td><td> 654.07800 </td><td>  704.78100</td><td> 472.08700 </td><td> 456.49600 </td><td> 601.335000</td><td>...        </td><td> 401.43000 </td><td> 757.49500 </td><td> 595.90800 </td><td> 381.23000 </td><td> 501.74400 </td><td>  659.6130 </td><td> 590.1560  </td><td> 461.23000 </td><td> 367.43300 </td><td> 790.94300 </td></tr>
	<tr><th scope=row>AFFX-CreX-5_at</th><td>  76.9437  </td><td>  40.90700 </td><td>  62.0314  </td><td>  54.42180 </td><td>  29.07040 </td><td>  19.52710 </td><td>   56.31640</td><td>  36.20440 </td><td>  34.41180 </td><td>  54.076500</td><td>...        </td><td>  57.84270 </td><td>  83.19140 </td><td>  66.67830 </td><td>  24.88520 </td><td>  61.95480 </td><td>   58.0325 </td><td>  28.7707  </td><td>  48.53460 </td><td>  31.14890 </td><td>  72.45670 </td></tr>
	<tr><th scope=row>AFFX-CreX-3_at</th><td> 105.3780  </td><td>  97.49320 </td><td>  74.0299  </td><td>  54.52770 </td><td>  54.98490 </td><td>  58.08770 </td><td>   96.63200</td><td>  52.73100 </td><td>  35.45880 </td><td>  60.264200</td><td>...        </td><td>  53.48370 </td><td> 108.54500 </td><td> 136.04400 </td><td>  43.86190 </td><td>  49.82890 </td><td>   88.3787 </td><td>  91.5539  </td><td>  40.22980 </td><td>  36.21510 </td><td>  73.53090 </td></tr>
	<tr><th scope=row>AFFX-BioB-5_st</th><td>  40.4826  </td><td>   7.45801 </td><td>  19.4069  </td><td>  20.62460 </td><td>  25.04960 </td><td>  12.48040 </td><td>   21.91020</td><td>  23.77200 </td><td>  24.18400 </td><td>  29.703200</td><td>...        </td><td>   6.53565 </td><td>  42.46960 </td><td>  41.46690 </td><td>  21.65480 </td><td>  34.61080 </td><td>   36.9725 </td><td>  16.9959  </td><td>   8.39896 </td><td>  41.02710 </td><td>  21.08800 </td></tr>
	<tr><th scope=row>AFFX-BioB-M_st</th><td>  58.1706  </td><td>  15.79260 </td><td>  25.1962  </td><td>  46.50570 </td><td>  15.31570 </td><td>  16.68330 </td><td>   93.17590</td><td>  -2.28600 </td><td>   9.00485 </td><td>  13.125300</td><td>...        </td><td>  26.92140 </td><td>  40.38730 </td><td>  17.68820 </td><td>  38.31000 </td><td>  18.02970 </td><td>   28.0120 </td><td>  24.8743  </td><td>  17.74280 </td><td>  34.79230 </td><td>  40.22130 </td></tr>
	<tr><th scope=row>AFFX-BioB-3_st</th><td> 257.6190  </td><td> 113.69000 </td><td> 187.7960  </td><td> 210.58000 </td><td> 137.39000 </td><td> 104.15900 </td><td>  296.28700</td><td> 110.53600 </td><td> 123.76700 </td><td> 165.210000</td><td>...        </td><td> 204.75000 </td><td> 265.77100 </td><td> 317.31400 </td><td>  88.07730 </td><td> 156.17900 </td><td>  249.8720 </td><td> 141.8030  </td><td> 117.47300 </td><td> 126.11400 </td><td> 210.56100 </td></tr>
	<tr><th scope=row>AFFX-BioC-5_st</th><td> 129.0560  </td><td>  74.60950 </td><td>  82.8271  </td><td> 101.53400 </td><td>  83.49860 </td><td>  73.19860 </td><td>  110.63100</td><td> 116.74200 </td><td> 149.32900 </td><td> 113.737000</td><td>...        </td><td> 103.28900 </td><td> 140.76000 </td><td> 177.44100 </td><td>  75.88880 </td><td>  87.36990 </td><td>  130.5090 </td><td>  91.7096  </td><td>  77.92770 </td><td> 129.62700 </td><td>  91.78030 </td></tr>
	<tr><th scope=row>AFFX-BioC-3_st</th><td>  61.7251  </td><td>  50.23720 </td><td>  61.6710  </td><td>  93.22350 </td><td>  38.11300 </td><td>  51.08690 </td><td>   69.02420</td><td>  51.73520 </td><td>  48.49430 </td><td>  66.332400</td><td>...        </td><td>  74.19100 </td><td>  74.71060 </td><td> 112.96400 </td><td>  63.23490 </td><td>  37.48920 </td><td>   69.9946 </td><td>  42.8123  </td><td>  63.80590 </td><td>  50.12460 </td><td> 134.58800 </td></tr>
	<tr><th scope=row>AFFX-BioDn-5_st</th><td> -40.9349  </td><td> -83.93020 </td><td> -28.7050  </td><td> -27.99790 </td><td> -29.90970 </td><td> -26.90040 </td><td>  -45.63120</td><td> -62.94740 </td><td> -31.43590 </td><td> -26.325300</td><td>...        </td><td> -30.48800 </td><td> -54.56340 </td><td> -49.08790 </td><td> -22.59160 </td><td> -21.23940 </td><td>  -38.7340 </td><td> -44.2136  </td><td> -23.69930 </td><td> -38.65320 </td><td> -40.72510 </td></tr>
	<tr><th scope=row>AFFX-BioDn-3_st</th><td> 284.4070  </td><td> 208.09900 </td><td> 239.0390  </td><td> 236.42800 </td><td> 152.32700 </td><td> 159.50500 </td><td>  316.93100</td><td> 152.18800 </td><td> 182.80300 </td><td> 275.020000</td><td>...        </td><td> 144.42100 </td><td> 316.03700 </td><td> 269.48500 </td><td> 148.11400 </td><td> 193.20500 </td><td>  294.2380 </td><td> 197.7450  </td><td> 203.53100 </td><td> 181.62300 </td><td> 234.95300 </td></tr>
	<tr><th scope=row>AFFX-CreX-5_st</th><td> 178.7450  </td><td> 101.30000 </td><td> 118.6990  </td><td> 131.83400 </td><td> 109.35500 </td><td>  98.17990 </td><td>  177.53300</td><td> 124.79500 </td><td>  86.07680 </td><td> 143.596000</td><td>...        </td><td> 107.08800 </td><td> 188.79200 </td><td> 240.35000 </td><td>  94.97540 </td><td> 104.14700 </td><td>  168.9920 </td><td> 109.5860  </td><td> 106.75700 </td><td>  90.19580 </td><td> 116.06900 </td></tr>
	<tr><th scope=row>AFFX-CreX-3_st</th><td>  79.7368  </td><td>  55.56320 </td><td>  68.5976  </td><td>  55.68810 </td><td>  56.39600 </td><td>  31.30030 </td><td>   84.84370</td><td>  33.42830 </td><td>  42.31720 </td><td> 124.882000</td><td>...        </td><td>  36.64480 </td><td>  94.36010 </td><td>  71.47270 </td><td>  24.36270 </td><td>  48.48060 </td><td>   79.3437 </td><td>  42.6647  </td><td> -12.64890 </td><td>  37.43890 </td><td>  62.07060 </td></tr>
	<tr><th scope=row>AFFX-hum_alu_at</th><td>9903.1900  </td><td>8501.62000 </td><td>9453.0000  </td><td>8595.65000 </td><td>9198.53000 </td><td>8729.83000 </td><td>10085.30000</td><td>5398.15000 </td><td>7851.25000 </td><td>9906.750000</td><td>...        </td><td>6945.46000 </td><td>9186.23000 </td><td>9889.05000 </td><td>8872.19000 </td><td>9682.71000 </td><td>10396.1000 </td><td>8365.8000  </td><td>9345.94000 </td><td>8252.65000 </td><td>4652.85000 </td></tr>
	<tr><th scope=row>AFFX-DapX-5_at</th><td>  61.2671  </td><td>  37.47400 </td><td>  44.7525  </td><td>  43.90200 </td><td>  40.56370 </td><td>  28.58190 </td><td>   49.28930</td><td>   7.59488 </td><td>  23.62900 </td><td>  57.042900</td><td>...        </td><td>  26.35510 </td><td>  70.37660 </td><td>  72.77700 </td><td>  16.69440 </td><td>  36.38470 </td><td>   58.1014 </td><td>  40.3792  </td><td>  33.92190 </td><td>  23.70460 </td><td>  45.14850 </td></tr>
	<tr><th scope=row>AFFX-DapX-M_at</th><td> 120.5440  </td><td>  75.98540 </td><td> 126.3740  </td><td>  90.40210 </td><td>  99.62140 </td><td>  59.88540 </td><td>  129.41900</td><td>  52.93500 </td><td>  64.18610 </td><td> 101.061000</td><td>...        </td><td>  72.23230 </td><td> 131.60200 </td><td> 136.20800 </td><td>  25.96240 </td><td>  76.89740 </td><td>  104.4280 </td><td>  68.4470  </td><td>  84.32500 </td><td>  46.61560 </td><td>  74.53150 </td></tr>
	<tr><th scope=row>AFFX-DapX-3_at</th><td>  50.0962  </td><td>  27.95320 </td><td>  29.2610  </td><td>  38.64360 </td><td>  34.48540 </td><td>  15.93390 </td><td>   55.84450</td><td>  20.09040 </td><td>  24.73830 </td><td>  38.872500</td><td>...        </td><td>  34.97870 </td><td>  56.82460 </td><td>  46.60260 </td><td>  14.22910 </td><td>  24.45630 </td><td>   59.6125 </td><td>  22.6792  </td><td>  26.54100 </td><td>  30.65290 </td><td>  29.44770 </td></tr>
	<tr><th scope=row>AFFX-LysX-5_at</th><td>  42.5285  </td><td>  33.71860 </td><td>  35.8420  </td><td>  43.81730 </td><td>  21.10380 </td><td>  26.00260 </td><td>   47.70150</td><td>  23.80350 </td><td>  11.37370 </td><td>  40.582200</td><td>...        </td><td>  31.92710 </td><td>  40.08590 </td><td>  49.11440 </td><td>  26.78740 </td><td>  50.96660 </td><td>   27.7013 </td><td>  30.9739  </td><td>  20.44470 </td><td>  14.69510 </td><td>  33.85440 </td></tr>
	<tr><th scope=row>AFFX-LysX-M_at</th><td>  36.8936  </td><td>  35.16970 </td><td>  36.5703  </td><td>  37.02740 </td><td>  24.05680 </td><td>  19.86490 </td><td>   57.41570</td><td>  18.31420 </td><td>  13.96590 </td><td>  37.389700</td><td>...        </td><td>  39.36620 </td><td>  48.33230 </td><td>  56.42690 </td><td>  18.70090 </td><td>  51.72850 </td><td>   60.7026 </td><td>  43.9691  </td><td>  12.95760 </td><td>  16.92430 </td><td>  35.29280 </td></tr>
	<tr><th scope=row>AFFX-LysX-3_at</th><td> 234.6980  </td><td> 102.46700 </td><td>  97.9010  </td><td> 146.23900 </td><td> 127.06800 </td><td>  65.07980 </td><td>  262.57900</td><td>  67.98070 </td><td>  67.95660 </td><td>  85.489600</td><td>...        </td><td>  88.35100 </td><td> 187.83100 </td><td> 247.96600 </td><td>  67.10260 </td><td>  73.02550 </td><td>  189.1860 </td><td> 104.8180  </td><td> 129.61900 </td><td>  91.40070 </td><td> 102.62100 </td></tr>
	<tr><th scope=row>AFFX-PheX-5_at</th><td>  26.9561  </td><td>  -2.37297 </td><td>  34.4333  </td><td>  17.59470 </td><td>  30.90680 </td><td>  19.45640 </td><td>   18.56280</td><td>   5.85058 </td><td>  26.22780 </td><td> -16.781100</td><td>...        </td><td>  20.20630 </td><td>  40.03670 </td><td>  30.28970 </td><td>  22.47410 </td><td>  24.21060 </td><td>   37.1239 </td><td>  23.3193  </td><td>  39.30200 </td><td>  25.46010 </td><td>  28.48380 </td></tr>
	<tr><th scope=row>AFFX-PheX-M_at</th><td>  58.1240  </td><td>  -1.69785 </td><td>  52.6747  </td><td>  55.70560 </td><td>  36.44900 </td><td>  27.13700 </td><td>   59.64770</td><td>  26.22140 </td><td>   8.82339 </td><td>  58.485100</td><td>...        </td><td>  47.75100 </td><td>  63.25770 </td><td>  15.99760 </td><td>   9.27033 </td><td>  81.02650 </td><td>   87.2832 </td><td>  30.0672  </td><td>  32.43330 </td><td>  20.76030 </td><td>  89.58840 </td></tr>
	<tr><th scope=row>AFFX-PheX-3_at</th><td>  40.6160  </td><td>  34.71300 </td><td>  82.4409  </td><td>  42.45960 </td><td>  50.35630 </td><td>  34.74800 </td><td>   26.40460</td><td>  35.91000 </td><td>  25.12950 </td><td>  49.387800</td><td>...        </td><td>  39.38820 </td><td>  23.54240 </td><td>  43.31840 </td><td>  47.50630 </td><td>  47.44000 </td><td>   32.3454 </td><td>  36.1116  </td><td>  36.69020 </td><td>  30.39780 </td><td> -53.83520 </td></tr>
	<tr><th scope=row>AFFX-ThrX-5_at</th><td> 125.0630  </td><td>  98.03690 </td><td>  77.1769  </td><td> 107.86100 </td><td>  72.35610 </td><td> 104.56000 </td><td>  103.89800</td><td>  61.46780 </td><td>  95.49600 </td><td>  74.320000</td><td>...        </td><td>  93.42100 </td><td> 133.26700 </td><td>  98.58730 </td><td>  63.13810 </td><td>  82.93620 </td><td>  129.2340 </td><td>  96.6750  </td><td>  94.55720 </td><td>  59.07050 </td><td>  71.21600 </td></tr>
	<tr><th scope=row>AFFX-ThrX-M_at</th><td>  49.9943  </td><td>  19.87220 </td><td>  32.2058  </td><td>  37.01370 </td><td>  38.76890 </td><td>  19.41190 </td><td>   42.14600</td><td>  26.79500 </td><td>  20.64170 </td><td>  31.576600</td><td>...        </td><td>  38.42290 </td><td>  40.16760 </td><td>  52.82150 </td><td>  15.62870 </td><td>  22.74920 </td><td>   49.4815 </td><td>  24.9362  </td><td>  27.08470 </td><td>  24.71060 </td><td>  31.35450 </td></tr>
	<tr><th scope=row>AFFX-ThrX-3_at</th><td>  33.1246  </td><td>   4.91484 </td><td>  24.6518  </td><td>  30.32560 </td><td>  23.33830 </td><td>  20.58270 </td><td>   42.63310</td><td>  20.73300 </td><td>  25.99310 </td><td>  36.588400</td><td>...        </td><td>  24.44730 </td><td>  26.71210 </td><td>  29.89500 </td><td>  16.48740 </td><td>  30.17520 </td><td>   39.3013 </td><td>  18.2112  </td><td>  15.24790 </td><td>  16.14300 </td><td>  23.46620 </td></tr>
	<tr><th scope=row>AFFX-TrpnX-5_at</th><td> 148.4940  </td><td>  70.92190 </td><td> 118.6320  </td><td> 118.69100 </td><td>  82.87940 </td><td>  69.07690 </td><td>  147.79700</td><td>  53.94830 </td><td>  55.52770 </td><td> 127.344000</td><td>...        </td><td>  97.06140 </td><td> 157.43000 </td><td> 104.56600 </td><td>  48.36610 </td><td>  86.59140 </td><td>  144.8190 </td><td>  79.4140  </td><td>  64.40910 </td><td>  72.01060 </td><td>  68.81200 </td></tr>
	<tr><th scope=row>AFFX-TrpnX-M_at</th><td>  66.6936  </td><td>   5.60854 </td><td>  60.3028  </td><td>  63.52760 </td><td>  33.52340 </td><td>  38.04180 </td><td>   60.98280</td><td>  40.46610 </td><td>  39.50320 </td><td>  29.885700</td><td>...        </td><td>  48.28420 </td><td>  93.52280 </td><td>  51.28240 </td><td>  29.35140 </td><td>  36.82390 </td><td>   89.5470 </td><td>  73.7415  </td><td>  52.74970 </td><td>  30.11150 </td><td>  38.83300 </td></tr>
	<tr><th scope=row>AFFX-TrpnX-3_at</th><td>  19.1364  </td><td>  20.80990 </td><td>  10.7195  </td><td>  13.97640 </td><td>  10.78760 </td><td>   3.04198 </td><td>    9.54646</td><td>   5.08160 </td><td>   6.61344 </td><td>  16.906400</td><td>...        </td><td>  20.41140 </td><td>  27.28890 </td><td>  14.16100 </td><td>   5.19944 </td><td>  14.76460 </td><td>   19.1634 </td><td>  11.4693  </td><td>  -2.53228 </td><td>   1.63297 </td><td>  16.62070 </td></tr>
	<tr><th scope=row>AFFX-HUMISGF3A/M97935_5_at</th><td> 165.7500  </td><td> 300.02400 </td><td> 152.4240  </td><td> 118.69300 </td><td>  95.32720 </td><td>  81.51550 </td><td>  152.29800</td><td> 122.23000 </td><td>  91.95470 </td><td> 129.385000</td><td>...        </td><td> 149.48000 </td><td> 195.41200 </td><td> 105.06100 </td><td>  80.66990 </td><td> 162.67300 </td><td>  192.3140 </td><td> 125.1290  </td><td>  77.66400 </td><td> 197.79700 </td><td> 188.11600 </td></tr>
	<tr><th scope=row>AFFX-HUMISGF3A/M97935_MA_at</th><td> 179.9890  </td><td> 790.94300 </td><td> 150.2490  </td><td>  82.63720 </td><td> 163.14000 </td><td> 147.75700 </td><td>  169.07800</td><td> 255.64600 </td><td> 135.59200 </td><td> 134.197000</td><td>...        </td><td>  91.99740 </td><td> 267.39500 </td><td> 121.30800 </td><td> 130.44600 </td><td> 253.17600 </td><td>  233.3040 </td><td> 234.0520  </td><td>  85.43020 </td><td> 499.75600 </td><td> 312.65300 </td></tr>
	<tr><th scope=row>AFFX-HUMISGF3A/M97935_MB_at</th><td> 151.7200  </td><td> 546.34300 </td><td> 173.6240  </td><td> 135.16300 </td><td> 193.15900 </td><td> 163.45800 </td><td>  145.28700</td><td> 144.06700 </td><td> 138.77400 </td><td> 161.823000</td><td>...        </td><td>  99.00990 </td><td> 209.60500 </td><td>  98.40660 </td><td>  92.67890 </td><td> 234.59500 </td><td>  232.4160 </td><td> 180.5480  </td><td>  81.73550 </td><td> 352.70800 </td><td> 309.97300 </td></tr>
	<tr><th scope=row>AFFX-HUMISGF3A/M97935_3_at</th><td> 553.8700  </td><td>1758.42000 </td><td> 599.8570  </td><td> 426.56900 </td><td> 859.04500 </td><td> 552.00600 </td><td>  499.94300</td><td> 701.33900 </td><td> 662.18300 </td><td> 603.143000</td><td>...        </td><td> 489.61600 </td><td> 582.48800 </td><td> 449.08100 </td><td> 991.45500 </td><td> 516.14200 </td><td>  691.1340 </td><td> 780.2940  </td><td> 981.64600 </td><td>1662.26000 </td><td> 976.52800 </td></tr>
	<tr><th scope=row>AFFX-HUMRGE/M10098_5_at</th><td>  72.2579  </td><td> 159.48400 </td><td>  42.4559  </td><td>   3.26054 </td><td>  11.48790 </td><td>  13.40200 </td><td>   20.54810</td><td> -18.06810 </td><td> -71.79690 </td><td>  14.806600</td><td>...        </td><td> -47.58540 </td><td>  25.00430 </td><td>  90.94600 </td><td> -47.55680 </td><td>   3.69034 </td><td>   15.2424 </td><td>  20.5804  </td><td>  38.51420 </td><td>   6.94581 </td><td>  19.52670 </td></tr>
	<tr><th scope=row>AFFX-HUMRGE/M10098_M_at</th><td> -30.1595  </td><td> -64.76580 </td><td> -28.2104  </td><td> -37.34640 </td><td> -45.81950 </td><td> -30.63210 </td><td>  -56.74850</td><td> -50.58990 </td><td> -43.36540 </td><td> -24.196500</td><td>...        </td><td> -50.12080 </td><td> -29.40050 </td><td>  -3.12303 </td><td> -27.88090 </td><td> -48.93650 </td><td>  -29.8916 </td><td> -27.1269  </td><td> -17.58330 </td><td> -38.14360 </td><td> -29.07920 </td></tr>
	<tr><th scope=row>AFFX-HUMRGE/M10098_3_at</th><td>  65.1004  </td><td>  19.71630 </td><td>  39.1968  </td><td>  24.39870 </td><td>  38.90760 </td><td>  22.41970 </td><td>   67.83080</td><td>  46.57300 </td><td>  53.04830 </td><td>  49.625200</td><td>...        </td><td>  40.99940 </td><td>  41.42280 </td><td>  52.17800 </td><td>  34.87810 </td><td>  65.82440 </td><td>   37.2982 </td><td>  27.4847  </td><td>  27.78670 </td><td>  53.72290 </td><td>  51.88410 </td></tr>
	<tr><th scope=row>AFFX-HUMGAPDH/M33197_5_at</th><td>1781.9500  </td><td>2370.97000 </td><td>1693.0000  </td><td> 931.98100 </td><td>2813.41000 </td><td>2773.80000 </td><td> 1331.06000</td><td>3409.56000 </td><td>2500.59000 </td><td>1397.650000</td><td>...        </td><td>2056.61000 </td><td>1071.56000 </td><td> 787.85500 </td><td>1764.84000 </td><td>1617.96000 </td><td> 1198.8000 </td><td>2573.8700  </td><td>1552.46000 </td><td>2323.13000 </td><td>5396.98000 </td></tr>
	<tr><th scope=row>AFFX-HUMGAPDH/M33197_M_at</th><td>3311.1800  </td><td>3270.14000 </td><td>2670.9400  </td><td>1916.00000 </td><td>3973.08000 </td><td>3533.69000 </td><td> 3001.52000</td><td>3670.05000 </td><td>3411.11000 </td><td>2582.350000</td><td>...        </td><td>3297.70000 </td><td>2141.77000 </td><td>2671.19000 </td><td>3354.22000 </td><td>2640.46000 </td><td> 2378.5800 </td><td>3477.1700  </td><td>3626.69000 </td><td>3276.53000 </td><td>6268.57000 </td></tr>
	<tr><th scope=row>AFFX-HUMGAPDH/M33197_3_at</th><td>4478.9900  </td><td>3937.24000 </td><td>3822.9400  </td><td>2995.59000 </td><td>4775.69000 </td><td>4276.28000 </td><td> 3922.80000</td><td>4113.84000 </td><td>3853.04000 </td><td>4196.470000</td><td>...        </td><td>5375.36000 </td><td>2913.78000 </td><td>4881.61000 </td><td>5156.98000 </td><td>3657.65000 </td><td> 2860.8700 </td><td>4471.0900  </td><td>5900.96000 </td><td>3914.43000 </td><td>8007.31000 </td></tr>
	<tr><th scope=row>AFFX-HSAC07/X00351_5_at</th><td>3835.3100  </td><td>5529.02000 </td><td>2961.3900  </td><td>1712.31000 </td><td>3090.42000 </td><td>4859.32000 </td><td> 4656.49000</td><td>4652.41000 </td><td>4628.45000 </td><td>1871.550000</td><td>...        </td><td>1481.99000 </td><td>4131.25000 </td><td>1094.09000 </td><td>2924.77000 </td><td>3819.17000 </td><td> 3879.7600 </td><td>4331.7700  </td><td> 965.54000 </td><td>3999.14000 </td><td>2593.41000 </td></tr>
	<tr><th scope=row>AFFX-HSAC07/X00351_M_at</th><td>4252.9800  </td><td>5758.12000 </td><td>3739.8200  </td><td>2266.59000 </td><td>4237.75000 </td><td>5339.43000 </td><td> 5809.61000</td><td>5529.77000 </td><td>5465.07000 </td><td>2697.630000</td><td>...        </td><td>2634.74000 </td><td>4546.89000 </td><td>1992.91000 </td><td>4220.09000 </td><td>4090.65000 </td><td> 4507.6200 </td><td>5276.0400  </td><td>2335.83000 </td><td>4771.61000 </td><td>2792.50000 </td></tr>
</tbody>
</table>




    Call:	 agnes(x = DIS, diss = TRUE, method = method, keep.diss = TRUE) 
    Agglomerative coefficient:  0.9367669 
    Order of objects:
     [1] A C U M E T G O V R D J S X Q B P Y F K W I H L N Z
    Height (summary):
         Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
    0.0004355 0.0016153 0.0047194 0.0204141 0.0160734 0.1565105 
    
    Available components:
    [1] "order"     "height"    "ac"        "merge"     "diss"      "call"     
    [7] "method"    "order.lab"



    Call:	 agnes(x = DIS, diss = TRUE, method = method, keep.diss = TRUE) 
    Agglomerative coefficient:  0.8619217 
    Order of objects:
     [1] AFFX-MurIL2_at              AFFX-CreX-5_at             
     [3] AFFX-PheX-M_at              AFFX-MurIL4_at             
     [5] AFFX-BioDn-3_at             AFFX-CreX-3_st             
     [7] AFFX-LysX-5_at              AFFX-LysX-M_at             
     [9] AFFX-BioB-M_st              AFFX-BioB-5_at             
    [11] AFFX-DapX-3_at              AFFX-ThrX-M_at             
    [13] AFFX-BioC-5_at              AFFX-BioB-3_st             
    [15] AFFX-CreX-5_st              AFFX-LysX-3_at             
    [17] AFFX-BioB-M_at              AFFX-BioDn-3_st            
    [19] AFFX-DapX-5_at              AFFX-DapX-M_at             
    [21] AFFX-TrpnX-5_at             AFFX-BioB-3_at             
    [23] AFFX-PheX-5_at              AFFX-ThrX-5_at             
    [25] AFFX-TrpnX-M_at             AFFX-TrpnX-3_at            
    [27] AFFX-BioB-5_st              AFFX-BioC-5_st             
    [29] AFFX-ThrX-3_at              AFFX-HUMRGE/M10098_3_at    
    [31] AFFX-hum_alu_at             AFFX-PheX-3_at             
    [33] AFFX-HUMGAPDH/M33197_5_at   AFFX-HUMGAPDH/M33197_M_at  
    [35] AFFX-HUMGAPDH/M33197_3_at   AFFX-MurIL10_at            
    [37] AFFX-BioDn-5_st             AFFX-CreX-3_at             
    [39] AFFX-HUMRGE/M10098_5_at     AFFX-HUMISGF3A/M97935_5_at 
    [41] AFFX-HUMISGF3A/M97935_MA_at AFFX-HUMISGF3A/M97935_MB_at
    [43] AFFX-HUMISGF3A/M97935_3_at  AFFX-MurFAS_at             
    [45] AFFX-BioC-3_at              AFFX-BioC-3_st             
    [47] AFFX-HSAC07/X00351_5_at     AFFX-HSAC07/X00351_M_at    
    [49] AFFX-BioDn-5_at             AFFX-HUMRGE/M10098_M_at    
    Height (summary):
       Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    0.03771 0.21669 0.41126 0.52716 0.60976 2.17552 
    
    Available components:
    [1] "order"     "height"    "ac"        "merge"     "diss"      "call"     
    [7] "method"    "order.lab"


1、To do a heatmap, use the `clustering.plot` function of the `EMA` library and pass the previously mentioned clusters as an argument.


```R
clustering.plot(tree=c.array, tree.sup=c.gene, data=c.data)
# This generates the required heatmap.
```


[plot1](https://github.com/Chengshu21/Chapter-5-Analyzing-Microarray-Data-with-R/blob/master/md/pic/output-1More%20visualizations%20for%20gene%20expression%20data.png)



```R
# But `heatmap` can do the same
heatmap(c.data)
# The following figure shows a heatmap for selected genes from the leukemiaexpression data, 
# which features normal cells and the different types of leukemia.
```


[plot2](https://github.com/Chengshu21/Chapter-5-Analyzing-Microarray-Data-with-R/blob/master/md/pic/output-2More%20visualizations%20for%20gene%20expression%20data.png)


2、The followning is a Venn diagram. It can be used to show the common and unique contents of two variables. This will require the following VennDiagram library.


```R
source("http://www.bioconductor.org/biocLite.R")
biocLite("VennDiagram")
library(VennDiagram)
```

    Bioconductor version 3.6 (BiocInstaller 1.28.0), ?biocLite for help
    BioC_mirror: https://bioconductor.org
    Using Bioconductor 3.6 (BiocInstaller 1.28.0), R 3.4.4 (2018-03-15).
    Installing package(s) 'VennDiagram'
    

    package 'VennDiagram' successfully unpacked and MD5 sums checked
    
    The downloaded binary packages are in
    	C:\Users\Master1\AppData\Local\Temp\RtmpeuY1NK\downloaded_packages
    

    Old packages: 'bindr', 'GenomicRanges', 'gsubfn', 'hms', 'igraph', 'Rcpp',
      'stringi', 'withr', 'xts', 'yaml'
    Loading required package: grid
    Loading required package: futile.logger
    

3、Create artificial data that consists of five variables (as sets) as a named list object.


```R
set <- list()
for(i in 1:5){
    set[[i]]=sample(LETTERS[1:20], replace=TRUE, prob=rep(0.05,20))
}
names(set)=c(paste("S", 1:5, sep=""))
```

4、To plot the Venn diagram, use the five sets created in the list set and create the `gList` object.


```R
venn.plot <- venn.diagram(x = set, filename = NULL, cat.cex = 1.5, alpha = 0.50, col = "black",
                          fill = c("dodgerblue","goldenrod1", "darkorange1", "seagreen3", "orchid3"),
                          cex= c(1.5, 1.5, 1.5, 1.5, 1.5, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8,
                                 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 1, 1, 1, 1, 1.5),
                          cat.col = c("dodgerblue", "goldenrod1", "darkorange1","seagreen3", "orchid3"),
                          cat.fontface = "bold", margin = 0.05)
```

5、Create the plot onto a PDF file using the `grid.draw` function. 


```R
pdf("D:/Try-practice/venn.pdf")
grid.draw(venn.plot)
dev.off()
```
[Venn plot](https://github.com/Chengshu21/Chapter-5-Analyzing-Microarray-Data-with-R/blob/master/md/pic/output-3More%20visualizations%20for%20gene%20expression%20data.png)

<strong>png:</strong> 2


6、The volcano plots are used to plot the fold changes and p-values against each other. Here to deal with a more intuitive representation of the volcano plot using the `ggplot2` library.<br>
And then select the top-ranked genes from the `limma` analysis as a data frame.


```R
library(ggplot2)
library(leukemiasEset)
data(leukemiasEset)
pheno <- pData(leukemiasEset)
mydata <- leukemiasEset[, sampleNames(leukemiasEset)[c(1:3, 13:15, 25:27, 49:51)]]
design <- model.matrix(~0 + factor(pData(mydata)$LeukemiaType))
colnames(design) <- unique(as.character(pData(mydata)$LeukemiaType))
library(limma)
fit <- lmFit(mydata, design)
contrast.matrix <- makeContrasts(NoL- ALL, NoL- AML, NoL- CLL, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
tested2 <- topTable(fit2,adjust="fdr",sort.by="B",number=Inf, coef=1)
DE2 <- tested2[tested2$adj.P.Val < 0.01,]
```


```R
head(tested2)
```


<table>
<thead><tr><th></th><th scope=col>logFC</th><th scope=col>AveExpr</th><th scope=col>t</th><th scope=col>P.Value</th><th scope=col>adj.P.Val</th><th scope=col>B</th></tr></thead>
<tbody>
	<tr><th scope=row>ENSG00000152078</th><td> 4.510507   </td><td>4.856523    </td><td> 28.13988   </td><td>4.463747e-11</td><td>9.004270e-07</td><td>14.01472    </td></tr>
	<tr><th scope=row>ENSG00000117519</th><td>-4.185175   </td><td>4.791585    </td><td>-22.73888   </td><td>3.878292e-10</td><td>3.911645e-06</td><td>12.69738    </td></tr>
	<tr><th scope=row>ENSG00000145850</th><td> 4.142236   </td><td>4.507655    </td><td> 17.38636   </td><td>5.759942e-09</td><td>2.925048e-05</td><td>10.72782    </td></tr>
	<tr><th scope=row>ENSG00000170180</th><td> 5.681327   </td><td>5.734169    </td><td> 17.37423   </td><td>5.800214e-09</td><td>2.925048e-05</td><td>10.72231    </td></tr>
	<tr><th scope=row>ENSG00000087586</th><td> 3.952183   </td><td>5.720789    </td><td> 16.45393   </td><td>9.977396e-09</td><td>3.111188e-05</td><td>10.28705    </td></tr>
	<tr><th scope=row>ENSG00000047597</th><td> 5.362419   </td><td>5.108415    </td><td> 16.32474   </td><td>1.079114e-08</td><td>3.111188e-05</td><td>10.22315    </td></tr>
</tbody>
</table>



7、Use the threshold the p-values and log fold change to define colors in the plot.


```R
threshold <- as.factor(tested2$logFC>1.5 &tested2$P.Value< 0.1)
```

8、Use the data and the `threshold` object to create a `g` object.


```R
g <- ggplot(data = tested2, aes (x= logFC, y = -log10(P.Value), colour =threshold)) + geom_point() + xlab("log2 fold change") + ylab("-log10 p-value")
# the codes "opts(legend.position = "right")" that I have deleted made it run smoothly
```

9、Plot the object by simply typing in the object name or save it as an object.


```R
g
```




[plot4](https://github.com/Chengshu21/Chapter-5-Analyzing-Microarray-Data-with-R/blob/master/md/pic/output-4More%20visualizations%20for%20gene%20expression%20data.png)


The heatmap uses the clustering results from the corresponding section. It uses two-way clustering for the arrays and genes. More arguments can be used to add information to the plot, such as labels and titles, and even the color palette can be changed. For more information, you can seek help in R using `?clustering.plot`.<br>

The Venn diagram recipe uses a named list object that contains the sets that the Venn diagram has to be drawn for. The code explained in this section is for five sets. In order to use the function for a different number of sets, parameters such as `colors` and `cex` must be modified. An example with a different number of sets has been explained in the help file for the function. Type in `?venn.diagram` for a detailed explanation. These plots can be used to show different attributes of results from data analysis. One possible scenario is that while comparing multiple treatments with a control, we can visualize how many DE genes/probes are shared across the experiments. Another possible situation can be the sharing of overexpressed GO terms in two or more analyses.<br>

The volcano plot is rather simple and uses the `ggplot` function instead of the unsual `plot` function. It segregates the genes from other genes in terms of colors, showing the genes with p-values less than 0.05 and the absolute log fold change greater than 1.5.<br>

热图使用相应部分的聚类结果.它使用阵列和基因的双向聚类,可以添加更多参数信息到情节，如标签和标题，甚至可以更改调色板。有关更多信息，可以使用`？clustering.plot`在R中寻求帮助.<br>
venn图使用一个包含绘制venn图必须需要的集合，作为列表对象。在这个部分中代码阐述了5个数据集。为了将函数用于不同数量的集合，必须修改`colors`和`cex`等参数。在函数的帮助文件中已经解释了具有不同数量集合的示例。请输入`？venn.diagram`查看详细解释。这些图可以用来显示数据分析结果的不同属性。一种可能的情况是，在将多种处理组与对照组进行比较时，我们可以想象在实验中共用了多少差异表达基因/探针。另一种可能的情况可能是在两个或多个分析中共用过表达的GO术语.<br>
火山图相当简单，使用`ggplot`函数而不是``plot``函数。它从颜色方面分离其他基因的基因，显示p值小于0.05且绝对对数倍数变化大于1.5的基因。





