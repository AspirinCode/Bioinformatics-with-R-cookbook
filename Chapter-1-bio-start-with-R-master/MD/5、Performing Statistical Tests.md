
# Performing Statistical Tests
When we handle statistical data, we want to know the significance of results in research or application. And by using R built-in functions, we can easily perform statistical tests to assess the significance and make quantitative decisions. The idea is to determine whether there is enough evidence to reject a conjecture about the results. In-built functions in R allow several such tests on data. The choice of test depends on the data and the question being asked. When we need to compare a group against a hypothetical value and our measurements follow the Gaussian distribution, we can use a one-sample t-test. However, if we have two paired groups (both measurements that follow the Gaussian distribution) being compared, it is better to use a paired t-test. R has built-in functions to carry out such tests, and in this recipe, we will try out some of these.


```R
data(sleep)
t.test(sleep[, 1]~sleep[, 2])
```


    
    	Welch Two Sample t-test
    
    data:  sleep[, 1] by sleep[, 2]
    t = -1.8608, df = 17.776, p-value = 0.07939
    alternative hypothesis: true difference in means is not equal to 0
    95 percent confidence interval:
     -3.3654832  0.2054832
    sample estimates:
    mean in group 1 mean in group 2 
               0.75            2.33 
    



```R
cont <- matrix(c(14, 33, 7, 3), ncol = 2)
cont
```


<table>
<tbody>
	<tr><td>14</td><td>7 </td></tr>
	<tr><td>33</td><td>3 </td></tr>
</tbody>
</table>




```R
colnames(cont) <- c("Sedan", "Convertible")
rownames(cont) <- c("Male", "Female")
cont
```


<table>
<thead><tr><th></th><th scope=col>Sedan</th><th scope=col>Convertible</th></tr></thead>
<tbody>
	<tr><th scope=row>Male</th><td>14</td><td>7 </td></tr>
	<tr><th scope=row>Female</th><td>33</td><td>3 </td></tr>
</tbody>
</table>




```R
chisq.test(as.table(cont))
```

    Warning message in chisq.test(as.table(cont)):
    "Chi-squared approximation may be incorrect"


    
    	Pearson's Chi-squared test with Yates' continuity correction
    
    data:  as.table(cont)
    X-squared = 4.1324, df = 1, p-value = 0.04207
    



```R
test <- chisq.test(as.table(cont))
test
```

    Warning message in chisq.test(as.table(cont)):
    "Chi-squared approximation may be incorrect"


    
    	Pearson's Chi-squared test with Yates' continuity correction
    
    data:  as.table(cont)
    X-squared = 4.1324, df = 1, p-value = 0.04207
    



```R
x <- c(1.83, 0.50, 1.62, 2.48, 1.68, 1.88, 1.55, 3.06, 1.30)
y <- c(0.878, 0.647, 0.598, 2.05, 1.06, 1.29, 1.06, 3.14, 1.29)
test <- wilcox.test(x, y, paried = T, alternative = "greater")
test
```

    Warning message in wilcox.test.default(x, y, paried = T, alternative = "greater"):
    "cannot compute exact p-value with ties"


    
    	Wilcoxon rank sum test with continuity correction
    
    data:  x and y
    W = 58, p-value = 0.06646
    alternative hypothesis: true location shift is greater than 0
    



```R
str(test)
# by str() function, you can see various parameters of the data.
```

    List of 7
     $ statistic  : Named num 58
      ..- attr(*, "names")= chr "W"
     $ parameter  : NULL
     $ p.value    : num 0.0665
     $ null.value : Named num 0
      ..- attr(*, "names")= chr "location shift"
     $ alternative: chr "greater"
     $ method     : chr "Wilcoxon rank sum test with continuity correction"
     $ data.name  : chr "x and y"
     - attr(*, "class")= chr "htest"
    


```R
test$p.value
```


0.0664597290926594




