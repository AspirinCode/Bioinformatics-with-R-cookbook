
# Filtering and Subsetting data
We always get enormous data after experiments or other ways, and they are stored in R as data frames. Data frames are
the primary structures of tabular data in R. We set up row-column format of these data frames. The columns of data can be of various types, such as numeric or factor. Through learning how to filter and subset data, we can easily extract the parts of these data frames we need, and add a new chunk, or filter a part that satisfies certain conditions.
### Filtering data
Just take the built-in data`iris` as an example.


```R
data(iris)
# load the `iris` data
```


```R
str(iris)
# To check the `iris` data's structure of its columns and rows
```

    'data.frame':	150 obs. of  5 variables:
     $ Sepal.Length: num  5.1 4.9 4.7 4.6 5 5.4 4.6 5 4.4 4.9 ...
     $ Sepal.Width : num  3.5 3 3.2 3.1 3.6 3.9 3.4 3.4 2.9 3.1 ...
     $ Petal.Length: num  1.4 1.4 1.3 1.5 1.4 1.7 1.4 1.5 1.4 1.5 ...
     $ Petal.Width : num  0.2 0.2 0.2 0.2 0.2 0.4 0.3 0.2 0.2 0.1 ...
     $ Species     : Factor w/ 3 levels "setosa","versicolor",..: 1 1 1 1 1 1 1 1 1 1 ...
    


```R
exiris <- head(iris)
exiris
# because the iris data has 150 rows, so I just take the first six rows as an example
```


<table>
<thead><tr><th scope=col>Sepal.Length</th><th scope=col>Sepal.Width</th><th scope=col>Petal.Length</th><th scope=col>Petal.Width</th><th scope=col>Species</th></tr></thead>
<tbody>
	<tr><td>5.1   </td><td>3.5   </td><td>1.4   </td><td>0.2   </td><td>setosa</td></tr>
	<tr><td>4.9   </td><td>3.0   </td><td>1.4   </td><td>0.2   </td><td>setosa</td></tr>
	<tr><td>4.7   </td><td>3.2   </td><td>1.3   </td><td>0.2   </td><td>setosa</td></tr>
	<tr><td>4.6   </td><td>3.1   </td><td>1.5   </td><td>0.2   </td><td>setosa</td></tr>
	<tr><td>5.0   </td><td>3.6   </td><td>1.4   </td><td>0.2   </td><td>setosa</td></tr>
	<tr><td>5.4   </td><td>3.9   </td><td>1.7   </td><td>0.4   </td><td>setosa</td></tr>
</tbody>
</table>




```R
myiris = data.frame(Sepal.Length = exiris$Sepal.Length, Sepal.Width = exiris$Sepal.Width, Species = exiris$Species)
myiris
# data.frame() can help to creates a data frame with the defined columns
# or the command: `myiris <- exiris[, c(1,2,5)]
# or the command: `myiris <- exiris[, -c[3, 4]]
```


<table>
<thead><tr><th scope=col>Sepal.Length</th><th scope=col>Sepal.Width</th><th scope=col>Species</th></tr></thead>
<tbody>
	<tr><td>5.1   </td><td>3.5   </td><td>setosa</td></tr>
	<tr><td>4.9   </td><td>3.0   </td><td>setosa</td></tr>
	<tr><td>4.7   </td><td>3.2   </td><td>setosa</td></tr>
	<tr><td>4.6   </td><td>3.1   </td><td>setosa</td></tr>
	<tr><td>5.0   </td><td>3.6   </td><td>setosa</td></tr>
	<tr><td>5.4   </td><td>3.9   </td><td>setosa</td></tr>
</tbody>
</table>




```R
myiris <- exiris[, c(1, 2, 5)]
myiris
```


<table>
<thead><tr><th scope=col>Sepal.Length</th><th scope=col>Sepal.Width</th><th scope=col>Species</th></tr></thead>
<tbody>
	<tr><td>5.1   </td><td>3.5   </td><td>setosa</td></tr>
	<tr><td>4.9   </td><td>3.0   </td><td>setosa</td></tr>
	<tr><td>4.7   </td><td>3.2   </td><td>setosa</td></tr>
	<tr><td>4.6   </td><td>3.1   </td><td>setosa</td></tr>
	<tr><td>5.0   </td><td>3.6   </td><td>setosa</td></tr>
	<tr><td>5.4   </td><td>3.9   </td><td>setosa</td></tr>
</tbody>
</table>




```R
Stalk.Length <- c(rnorm(6, 1, 0.1))
Stalk.Length
newdata <- data.frame(Sepal.Length = 4.8, Sepal.Width = 3.3, Petal.Length = 1.6, Petal.Width = 0.2, 
                      Species = "myspecies")
newdata
# the `rnorm` function generates a random sample from a normal distribution and more details in `?rnorm`
```


<ol class=list-inline>
	<li>1.05203835866992</li>
	<li>0.831959573312066</li>
	<li>0.919135216666077</li>
	<li>1.08858589787697</li>
	<li>0.95939020572669</li>
	<li>1.14613849357116</li>
</ol>




<table>
<thead><tr><th scope=col>Sepal.Length</th><th scope=col>Sepal.Width</th><th scope=col>Petal.Length</th><th scope=col>Petal.Width</th><th scope=col>Species</th></tr></thead>
<tbody>
	<tr><td>4.8      </td><td>3.3      </td><td>1.6      </td><td>0.2      </td><td>myspecies</td></tr>
</tbody>
</table>




```R
myiris <- cbind(exiris, Stalk.Length)
myiris
# `cbind` function adds a new column
```


<table>
<thead><tr><th scope=col>Sepal.Length</th><th scope=col>Sepal.Width</th><th scope=col>Petal.Length</th><th scope=col>Petal.Width</th><th scope=col>Species</th><th scope=col>Stalk.Length</th></tr></thead>
<tbody>
	<tr><td>5.1      </td><td>3.5      </td><td>1.4      </td><td>0.2      </td><td>setosa   </td><td>1.0520384</td></tr>
	<tr><td>4.9      </td><td>3.0      </td><td>1.4      </td><td>0.2      </td><td>setosa   </td><td>0.8319596</td></tr>
	<tr><td>4.7      </td><td>3.2      </td><td>1.3      </td><td>0.2      </td><td>setosa   </td><td>0.9191352</td></tr>
	<tr><td>4.6      </td><td>3.1      </td><td>1.5      </td><td>0.2      </td><td>setosa   </td><td>1.0885859</td></tr>
	<tr><td>5.0      </td><td>3.6      </td><td>1.4      </td><td>0.2      </td><td>setosa   </td><td>0.9593902</td></tr>
	<tr><td>5.4      </td><td>3.9      </td><td>1.7      </td><td>0.4      </td><td>setosa   </td><td>1.1461385</td></tr>
</tbody>
</table>




```R
myiris <- rbind(exiris, newdata)
myiris
# `rbind` function adds a new row 
```


<table>
<thead><tr><th scope=col>Sepal.Length</th><th scope=col>Sepal.Width</th><th scope=col>Petal.Length</th><th scope=col>Petal.Width</th><th scope=col>Species</th></tr></thead>
<tbody>
	<tr><td>5.1      </td><td>3.5      </td><td>1.4      </td><td>0.2      </td><td>setosa   </td></tr>
	<tr><td>4.9      </td><td>3.0      </td><td>1.4      </td><td>0.2      </td><td>setosa   </td></tr>
	<tr><td>4.7      </td><td>3.2      </td><td>1.3      </td><td>0.2      </td><td>setosa   </td></tr>
	<tr><td>4.6      </td><td>3.1      </td><td>1.5      </td><td>0.2      </td><td>setosa   </td></tr>
	<tr><td>5.0      </td><td>3.6      </td><td>1.4      </td><td>0.2      </td><td>setosa   </td></tr>
	<tr><td>5.4      </td><td>3.9      </td><td>1.7      </td><td>0.4      </td><td>setosa   </td></tr>
	<tr><td>4.8      </td><td>3.3      </td><td>1.6      </td><td>0.2      </td><td>myspecies</td></tr>
</tbody>
</table>



There is another way of typing commands to add a column or row as follows:


```R
myiris$Stalk.Length = c(rnorm(7,1,0.1))
myiris
# at this step, 'myiris' data has changed into a data frame with seven rows and five columns..
# so the rnorm's 'n' parameter should be 7.
# The $ sign placed after the data followed by the column name specifies the data in that column.
```


<table>
<thead><tr><th scope=col>Sepal.Length</th><th scope=col>Sepal.Width</th><th scope=col>Petal.Length</th><th scope=col>Petal.Width</th><th scope=col>Species</th><th scope=col>Stalk.Length</th></tr></thead>
<tbody>
	<tr><td>5.1      </td><td>3.5      </td><td>1.4      </td><td>0.2      </td><td>setosa   </td><td>1.1665048</td></tr>
	<tr><td>4.9      </td><td>3.0      </td><td>1.4      </td><td>0.2      </td><td>setosa   </td><td>0.9672891</td></tr>
	<tr><td>4.7      </td><td>3.2      </td><td>1.3      </td><td>0.2      </td><td>setosa   </td><td>1.1543450</td></tr>
	<tr><td>4.6      </td><td>3.1      </td><td>1.5      </td><td>0.2      </td><td>setosa   </td><td>1.0114124</td></tr>
	<tr><td>5.0      </td><td>3.6      </td><td>1.4      </td><td>0.2      </td><td>setosa   </td><td>1.1273507</td></tr>
	<tr><td>5.4      </td><td>3.9      </td><td>1.7      </td><td>0.4      </td><td>setosa   </td><td>1.0666074</td></tr>
	<tr><td>4.8      </td><td>3.3      </td><td>1.6      </td><td>0.2      </td><td>myspecies</td><td>0.9292590</td></tr>
</tbody>
</table>




```R
dim(myiris)
```


<ol class=list-inline>
	<li>7</li>
	<li>6</li>
</ol>




```R
colnames(myiris)
```


<ol class=list-inline>
	<li>'Sepal.Length'</li>
	<li>'Sepal.Width'</li>
	<li>'Petal.Length'</li>
	<li>'Petal.Width'</li>
	<li>'Species'</li>
	<li>'Stalk.Length'</li>
</ol>




```R
myiris[7, ]
# extract the seventh row
```


<table>
<thead><tr><th></th><th scope=col>Sepal.Length</th><th scope=col>Sepal.Width</th><th scope=col>Petal.Length</th><th scope=col>Petal.Width</th><th scope=col>Species</th><th scope=col>Stalk.Length</th></tr></thead>
<tbody>
	<tr><th scope=row>7</th><td>4.8      </td><td>3.3      </td><td>1.6      </td><td>0.2      </td><td>myspecies</td><td>1.040109 </td></tr>
</tbody>
</table>




```R
myiris[, 3]
# extract the third column, it just shows the numeric data from the table...
```


<ol class=list-inline>
	<li>1.4</li>
	<li>1.4</li>
	<li>1.3</li>
	<li>1.5</li>
	<li>1.4</li>
	<li>1.7</li>
	<li>1.6</li>
</ol>



### Subseting data
Via using `subset()`, we can extract a part from the data frame, which meets certain conditions.<br>
Take `iris` data as an example.


```R
myiris1 <- subset(iris, Sepal.Length == 5.1)
myiris1
# `subset()`
```


<table>
<thead><tr><th></th><th scope=col>Sepal.Length</th><th scope=col>Sepal.Width</th><th scope=col>Petal.Length</th><th scope=col>Petal.Width</th><th scope=col>Species</th></tr></thead>
<tbody>
	<tr><th scope=row>1</th><td>5.1       </td><td>3.5       </td><td>1.4       </td><td>0.2       </td><td>setosa    </td></tr>
	<tr><th scope=row>18</th><td>5.1       </td><td>3.5       </td><td>1.4       </td><td>0.3       </td><td>setosa    </td></tr>
	<tr><th scope=row>20</th><td>5.1       </td><td>3.8       </td><td>1.5       </td><td>0.3       </td><td>setosa    </td></tr>
	<tr><th scope=row>22</th><td>5.1       </td><td>3.7       </td><td>1.5       </td><td>0.4       </td><td>setosa    </td></tr>
	<tr><th scope=row>24</th><td>5.1       </td><td>3.3       </td><td>1.7       </td><td>0.5       </td><td>setosa    </td></tr>
	<tr><th scope=row>40</th><td>5.1       </td><td>3.4       </td><td>1.5       </td><td>0.2       </td><td>setosa    </td></tr>
	<tr><th scope=row>45</th><td>5.1       </td><td>3.8       </td><td>1.9       </td><td>0.4       </td><td>setosa    </td></tr>
	<tr><th scope=row>47</th><td>5.1       </td><td>3.8       </td><td>1.6       </td><td>0.2       </td><td>setosa    </td></tr>
	<tr><th scope=row>99</th><td>5.1       </td><td>2.5       </td><td>3.0       </td><td>1.1       </td><td>versicolor</td></tr>
</tbody>
</table>




```R
myiris2 <- iris[iris$Sepal.Length == 5.1, ]
myiris2
```


<table>
<thead><tr><th></th><th scope=col>Sepal.Length</th><th scope=col>Sepal.Width</th><th scope=col>Petal.Length</th><th scope=col>Petal.Width</th><th scope=col>Species</th></tr></thead>
<tbody>
	<tr><th scope=row>1</th><td>5.1       </td><td>3.5       </td><td>1.4       </td><td>0.2       </td><td>setosa    </td></tr>
	<tr><th scope=row>18</th><td>5.1       </td><td>3.5       </td><td>1.4       </td><td>0.3       </td><td>setosa    </td></tr>
	<tr><th scope=row>20</th><td>5.1       </td><td>3.8       </td><td>1.5       </td><td>0.3       </td><td>setosa    </td></tr>
	<tr><th scope=row>22</th><td>5.1       </td><td>3.7       </td><td>1.5       </td><td>0.4       </td><td>setosa    </td></tr>
	<tr><th scope=row>24</th><td>5.1       </td><td>3.3       </td><td>1.7       </td><td>0.5       </td><td>setosa    </td></tr>
	<tr><th scope=row>40</th><td>5.1       </td><td>3.4       </td><td>1.5       </td><td>0.2       </td><td>setosa    </td></tr>
	<tr><th scope=row>45</th><td>5.1       </td><td>3.8       </td><td>1.9       </td><td>0.4       </td><td>setosa    </td></tr>
	<tr><th scope=row>47</th><td>5.1       </td><td>3.8       </td><td>1.6       </td><td>0.2       </td><td>setosa    </td></tr>
	<tr><th scope=row>99</th><td>5.1       </td><td>2.5       </td><td>3.0       </td><td>1.1       </td><td>versicolor</td></tr>
</tbody>
</table>




```R
myiris3 <- subset(myiris, Species == "setosa")
myiris3
```


<table>
<thead><tr><th scope=col>Sepal.Length</th><th scope=col>Sepal.Width</th><th scope=col>Petal.Length</th><th scope=col>Petal.Width</th><th scope=col>Species</th><th scope=col>Stalk.Length</th></tr></thead>
<tbody>
	<tr><td>5.1      </td><td>3.5      </td><td>1.4      </td><td>0.2      </td><td>setosa   </td><td>0.8962507</td></tr>
	<tr><td>4.9      </td><td>3.0      </td><td>1.4      </td><td>0.2      </td><td>setosa   </td><td>0.9685493</td></tr>
	<tr><td>4.7      </td><td>3.2      </td><td>1.3      </td><td>0.2      </td><td>setosa   </td><td>0.8629482</td></tr>
	<tr><td>4.6      </td><td>3.1      </td><td>1.5      </td><td>0.2      </td><td>setosa   </td><td>0.8499977</td></tr>
	<tr><td>5.0      </td><td>3.6      </td><td>1.4      </td><td>0.2      </td><td>setosa   </td><td>0.9015226</td></tr>
	<tr><td>5.4      </td><td>3.9      </td><td>1.7      </td><td>0.4      </td><td>setosa   </td><td>0.7995041</td></tr>
</tbody>
</table>




```R
myiris3[1, ]
```


<table>
<thead><tr><th scope=col>Sepal.Length</th><th scope=col>Sepal.Width</th><th scope=col>Petal.Length</th><th scope=col>Petal.Width</th><th scope=col>Species</th><th scope=col>Stalk.Length</th></tr></thead>
<tbody>
	<tr><td>5.1      </td><td>3.5      </td><td>1.4      </td><td>0.2      </td><td>setosa   </td><td>0.8962507</td></tr>
</tbody>
</table>



###### Another way to select part of the data is using %in% operators with the data frame, as follows:


```R
mylength <- c(5, 6, 7.2)
mynew.iris <- iris[iris[, 1] %in% mylength, ]
mynew.iris
# This selects all the rows from the data that meet the defined condition. 
# The condition here means that the value in column 1 of iris is the same as (matching) any value in the mylength vector. 
# The extracted rows are then assigned to a new object, 'mynew.iris'.
```


<table>
<thead><tr><th></th><th scope=col>Sepal.Length</th><th scope=col>Sepal.Width</th><th scope=col>Petal.Length</th><th scope=col>Petal.Width</th><th scope=col>Species</th></tr></thead>
<tbody>
	<tr><th scope=row>5</th><td>5.0       </td><td>3.6       </td><td>1.4       </td><td>0.2       </td><td>setosa    </td></tr>
	<tr><th scope=row>8</th><td>5.0       </td><td>3.4       </td><td>1.5       </td><td>0.2       </td><td>setosa    </td></tr>
	<tr><th scope=row>26</th><td>5.0       </td><td>3.0       </td><td>1.6       </td><td>0.2       </td><td>setosa    </td></tr>
	<tr><th scope=row>27</th><td>5.0       </td><td>3.4       </td><td>1.6       </td><td>0.4       </td><td>setosa    </td></tr>
	<tr><th scope=row>36</th><td>5.0       </td><td>3.2       </td><td>1.2       </td><td>0.2       </td><td>setosa    </td></tr>
	<tr><th scope=row>41</th><td>5.0       </td><td>3.5       </td><td>1.3       </td><td>0.3       </td><td>setosa    </td></tr>
	<tr><th scope=row>44</th><td>5.0       </td><td>3.5       </td><td>1.6       </td><td>0.6       </td><td>setosa    </td></tr>
	<tr><th scope=row>50</th><td>5.0       </td><td>3.3       </td><td>1.4       </td><td>0.2       </td><td>setosa    </td></tr>
	<tr><th scope=row>61</th><td>5.0       </td><td>2.0       </td><td>3.5       </td><td>1.0       </td><td>versicolor</td></tr>
	<tr><th scope=row>63</th><td>6.0       </td><td>2.2       </td><td>4.0       </td><td>1.0       </td><td>versicolor</td></tr>
	<tr><th scope=row>79</th><td>6.0       </td><td>2.9       </td><td>4.5       </td><td>1.5       </td><td>versicolor</td></tr>
	<tr><th scope=row>84</th><td>6.0       </td><td>2.7       </td><td>5.1       </td><td>1.6       </td><td>versicolor</td></tr>
	<tr><th scope=row>86</th><td>6.0       </td><td>3.4       </td><td>4.5       </td><td>1.6       </td><td>versicolor</td></tr>
	<tr><th scope=row>94</th><td>5.0       </td><td>2.3       </td><td>3.3       </td><td>1.0       </td><td>versicolor</td></tr>
	<tr><th scope=row>110</th><td>7.2       </td><td>3.6       </td><td>6.1       </td><td>2.5       </td><td>virginica </td></tr>
	<tr><th scope=row>120</th><td>6.0       </td><td>2.2       </td><td>5.0       </td><td>1.5       </td><td>virginica </td></tr>
	<tr><th scope=row>126</th><td>7.2       </td><td>3.2       </td><td>6.0       </td><td>1.8       </td><td>virginica </td></tr>
	<tr><th scope=row>130</th><td>7.2       </td><td>3.0       </td><td>5.8       </td><td>1.6       </td><td>virginica </td></tr>
	<tr><th scope=row>139</th><td>6.0       </td><td>3.0       </td><td>4.8       </td><td>1.8       </td><td>virginica </td></tr>
</tbody>
</table>




