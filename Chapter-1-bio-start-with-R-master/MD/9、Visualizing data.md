
# Visualizing data 

To get a direct understanding of our various data, these plotting functions manipulated in R can help us make these boring and numerous<br>
data to be a vivid plot or graph. Take the `iris` data as an example.<br>
More details in R help files.


```R
sl <- iris[,1]
pl <- iris[,4]
plot(x=pl, y=sl, xlab="Petal length", ylab="Sepal length",col="black", main="Varition of sepal length with petal length")
# 'plot(with(iris, plot(x = Sepal.Length, y=Petal.Length))' can be another way to get the same results. 
# The argument for 'xlab' and 'ylab' respectively is used to change the axis labels.
# The main argument gives the plot a title..
```


[png2](https://github.com/Chengshu21/Bioinformatics-with-R-cookbook/blob/master/Chapter-1-bio-start-with-R-master/MD/9、Visualizing%20out_put_2_0.png)



```R
boxplot(Sepal.Length~Species, data=iris, ylab="sepal length", xlab="Species", main="Sepal length for different species")
# create a boxplot for the data
```


[png3](https://github.com/Chengshu21/Bioinformatics-with-R-cookbook/blob/master/Chapter-1-bio-start-with-R-master/MD/9、Visualizing%20output_3_0.png)



```R
genex <- c(rnorm(100, 1, 0.1), rnorm(100, 2, 0.1), rnorm(50, 3, 0.1))
plot(x=genex, xlim=c(1,5), type='l', main="line diagram")
# 'type' argument is used to decide the type of the plot.
```


[png4](https://github.com/Chengshu21/Bioinformatics-with-R-cookbook/blob/master/Chapter-1-bio-start-with-R-master/MD/9、Visualizing%20output_4_0.png)



```R
x <- rnorm(1000, 3, 0.02)
hist(x)
# Histograms can used to visualize the density of the data and the frequency of every category.
```


[png5](https://github.com/Chengshu21/Bioinformatics-with-R-cookbook/blob/master/Chapter-1-bio-start-with-R-master/MD/9、Visualizing%20output_5_0.png)

