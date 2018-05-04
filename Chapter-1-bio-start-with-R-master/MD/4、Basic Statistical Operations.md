
# Basic Statistical Operations
R being a statistical programming environment has a number of built-in functionalities to perform statistics on data. Nevertheless, some specific functionalities are either available in packages or can easily be written. 


```R
summary(iris)
# 'summary()' is a generic function used to produce result summaries of the results of various model fitting functions.
# The function invokes particular methods which depend on the class of the first argument. 
```


      Sepal.Length    Sepal.Width     Petal.Length    Petal.Width   
     Min.   :4.300   Min.   :2.000   Min.   :1.000   Min.   :0.100  
     1st Qu.:5.100   1st Qu.:2.800   1st Qu.:1.600   1st Qu.:0.300  
     Median :5.800   Median :3.000   Median :4.350   Median :1.300  
     Mean   :5.843   Mean   :3.057   Mean   :3.758   Mean   :1.199  
     3rd Qu.:6.400   3rd Qu.:3.300   3rd Qu.:5.100   3rd Qu.:1.800  
     Max.   :7.900   Max.   :4.400   Max.   :6.900   Max.   :2.500  
           Species  
     setosa    :50  
     versicolor:50  
     virginica :50  
                    
                    
                    



```R
mean(iris[, 1])
# 'mean()' means arithmetic mean or type '?mean()' to know its meanning.
```


5.84333333333333



```R
sd(iris[, 1])
# This function computes the standard deviation of the values in the object, 'x'. 
```


0.828066127977863



```R
cor(iris[, 1], iris[, 2])
cor(iris[, 1], iris[, 3])
# 'cor' computes the correlation of x and y if these are vectors.  
```


-0.117569784133002



0.871753775886583



```R
Cov.mat <- cov(iris[, 1:4])
Cov.mat
# 'cov()' computes the covariance of x and y if these are data matrix. 
```


<table>
<thead><tr><th></th><th scope=col>Sepal.Length</th><th scope=col>Sepal.Width</th><th scope=col>Petal.Length</th><th scope=col>Petal.Width</th></tr></thead>
<tbody>
	<tr><th scope=row>Sepal.Length</th><td> 0.6856935</td><td>-0.0424340</td><td> 1.2743154</td><td> 0.5162707</td></tr>
	<tr><th scope=row>Sepal.Width</th><td>-0.0424340</td><td> 0.1899794</td><td>-0.3296564</td><td>-0.1216394</td></tr>
	<tr><th scope=row>Petal.Length</th><td> 1.2743154</td><td>-0.3296564</td><td> 3.1162779</td><td> 1.2956094</td></tr>
	<tr><th scope=row>Petal.Width</th><td> 0.5162707</td><td>-0.1216394</td><td> 1.2956094</td><td> 0.5810063</td></tr>
</tbody>
</table>



* It is important to know how to handle the missing data. If we have missing values (called NA in R) in our data, we can set the `na.rm` argument to `TRUE`, and the computation will be done only based on non-NA values. 


```R
a <- c(1:4, NA, 6)
```


```R
mean(a)
```


&lt;NA&gt;



```R
mean(a, na.rm = TRUE)
```


3.2



```R

```
