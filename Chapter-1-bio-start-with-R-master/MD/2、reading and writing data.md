
# Reading and Writing Data

Before we start with analyzing any data, we must load it into our R workspace. Not only an external R object (typical file extensions are .rda or .RData, but also an internal R object for a package or a TXT, CSV, or Excel file can be loaded into the R console. At first, I learn how to find and set the working directory, which is convenient for me to find the files.

### 1、Get the current working directory:


```R
getwd()
```


'C:/Users/Administrator'


### 2、Change to a desired directory, by the command `setwd("Path/to/desired directory")`, I change the directory to F:/Git/Chapter 1/Reading and Writing Data, it is important that there is a folder existed:


```R
setwd("F:/Git/Chapter 1/Reading and Writing Data")
```

### 3、Reading and Writing Data
I will learn this content form the following steps: 
* [Creating data frames](https://github.com/Chengshu21/bio-start-with-R/blob/master/Chapter%201/reading%20and%20writing%20data.md#3.1-creating-data-frames)
* [Writing data)](https://github.com/Chengshu21/bio-start-with-R/blob/master/Chapter%201/reading%20and%20writing%20data.md#3.2-writing-data)
* [Reading data](https://github.com/Chengshu21/bio-start-with-R/blob/master/Chapter%201/reading%20and%20writing%20data.md#3.3-reading-data)

If there is a dataset existed in the package, and we can use `data()` to get an access to the dataset.
#### 3.1 Creating data frame


```R
d <- data.frame(obs = c(1, 2, 3, 4, 5, 6), treat = c("A", "B", "A", "A", "O", "B"), weight = c(2.3, NA, 9, 8, 4, 7))
d
```


<table>
<thead><tr><th scope=col>obs</th><th scope=col>treat</th><th scope=col>weight</th></tr></thead>
<tbody>
	<tr><td>1  </td><td>A  </td><td>2.3</td></tr>
	<tr><td>2  </td><td>B  </td><td> NA</td></tr>
	<tr><td>3  </td><td>A  </td><td>9.0</td></tr>
	<tr><td>4  </td><td>A  </td><td>8.0</td></tr>
	<tr><td>5  </td><td>O  </td><td>4.0</td></tr>
	<tr><td>6  </td><td>B  </td><td>7.0</td></tr>
</tbody>
</table>



#### 3.2 Writing data
More detailes in [R](http://127.0.0.1:29909/library/utils/html/write.table.html).


```R
write.table(d, file = "F:/Git/Chapter 1/Reading and Writing Data/data1.txt", row.names = F, quote = F)
# you can check the data1.text in the path
write.table(d, file = "F:/Git/Chapter 1/Reading and Writing Data/data2.txt", row.names = F, quote = F, sep = "\t")
write.table(d, file = "F:/Git/Chapter 1/Reading and Writing Data/data3.txt", row.names = F, quote = T, sep = "\t")
```

`write.table` is used to write data frames or table objects into a table file or other files. Due to the three `.txt` files, we can fine the parameters' difference.<br>
Of course, there are `write.csv` command to write `.csv` files.<br>
But if you want to save as a `.r` file, the command should be changed into `save()`.<br>
If forgetting some parameters and the function, it's a good way to type `??write`.


```R
write.csv(d, file = "F:/Git/Chapter 1/Reading and Writing Data/data4.csv", row.names = F,quote = F)
save(d, file = "F:/Git/Chapter 1/Reading and Writing Data/data5.RData")
```

#### 3.3 Reading data
To read the tabular data in the form of a `.csv` file with `read.csv` or `read.table`, there are some differences in the two commands, running the following commands and check the functions by `?read.csv()` and `?read.table`.<br>
There are other ways to read data, like `read.fwf( )`, `read.delim` or use packages `xlsx`.


```R
mydata <- read.table("F:/Git/Chapter 1/Reading and Writing Data/data1.txt")
mydata
```


<table>
<thead><tr><th scope=col>V1</th><th scope=col>V2</th><th scope=col>V3</th></tr></thead>
<tbody>
	<tr><td>obs   </td><td>treat </td><td>weight</td></tr>
	<tr><td>1     </td><td>A     </td><td>2.3   </td></tr>
	<tr><td>2     </td><td>B     </td><td>NA    </td></tr>
	<tr><td>3     </td><td>A     </td><td>9     </td></tr>
	<tr><td>4     </td><td>A     </td><td>8     </td></tr>
	<tr><td>5     </td><td>O     </td><td>4     </td></tr>
	<tr><td>6     </td><td>B     </td><td>7     </td></tr>
</tbody>
</table>




```R
mydata <- read.table("F:/Git/Chapter 1/Reading and Writing Data/data1.txt", header = TRUE)
# you can try add different parameters to check the `read.table`
mydata
```


<table>
<thead><tr><th scope=col>obs</th><th scope=col>treat</th><th scope=col>weight</th></tr></thead>
<tbody>
	<tr><td>1  </td><td>A  </td><td>2.3</td></tr>
	<tr><td>2  </td><td>B  </td><td> NA</td></tr>
	<tr><td>3  </td><td>A  </td><td>9.0</td></tr>
	<tr><td>4  </td><td>A  </td><td>8.0</td></tr>
	<tr><td>5  </td><td>O  </td><td>4.0</td></tr>
	<tr><td>6  </td><td>B  </td><td>7.0</td></tr>
</tbody>
</table>




```R
mydata <- read.csv("F:/Git/Chapter 1/Reading and Writing Data/data1.txt")
mydata
```


<table>
<thead><tr><th scope=col>obs.treat.weight</th></tr></thead>
<tbody>
	<tr><td>1 A 2.3</td></tr>
	<tr><td>2 B NA </td></tr>
	<tr><td>3 A 9  </td></tr>
	<tr><td>4 A 8  </td></tr>
	<tr><td>5 O 4  </td></tr>
	<tr><td>6 B 7  </td></tr>
</tbody>
</table>




```R
?read.csv
```

It is also possible to read an Excel file in R. That needs the package like `xlsx` and `gdata`. The `xlsx` package requires Java settings, while `gdata` is relatively simple. However, the `xlsx` package offers more functionalities, such as read permissions for different sheets in a workbook and the newer versions of Excel files. I will take the `xlsx` package. Use the `read.xlsx` function to read an Excel file as follows:


```R
options(repos=structure(c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")))    
install.packages("xlsx", dependencies=TRUE)
```

    package 'xlsx' successfully unpacked and MD5 sums checked
    
    The downloaded binary packages are in
    	C:\Users\Administrator\AppData\Local\Temp\RtmpUzKBeW\downloaded_packages
    

```R
library(xlsx)
# Before you check and load 'xlsx' packages, you should set up the'Java' environment. If not, there will be warning messages to inform you istall and create the Java environment to load "xlsx" packages.
```


```R
mydata6 <- read.xls("F:/Git/Chapter 1/Reading and Writing Data/data6.xls")
mydata6
```


