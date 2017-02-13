---
layout: post
title: "Missing data mechanism and its impact on model estimates"
excerpt: "Simulation study from fitted linear model"
date: "2017-02-13 10:06:13"
category: statistical-methods
share: true
published: true
comments: true
tags: [simulation, missing data, dplyr, MCAR, multiple regression, bivariate, distribution, linear model, MNAR, MAR, bias, randomization] 
---

Today we will asses the impact of missing data on estimates from fitting linear regression model. 

First we will generate to covariate **X1** and **X2** from bivariate normal distribution with specific correlation **r**. And, then generate **Y** with known value of $$ \beta_{0} $$  and $$ \beta_{1} $$   and some errors.

Then in each dataset, we will impute some missing data with known missing data mechanism **(MCAR, MAR and MNAR)** and will see whether **simulation** (generating data 1000 times) will enable us to estimate true value of parameters or not ( means biased or unbiased) in data with missing values.




We will need following packages,


{% highlight r %}
require(dplyr)
require(ggplot2)
require(broom)
require(knitr)
require(cowplot)
require(pander)
require(mvtnorm)
{% endhighlight %}


Let first write function which will generate **x1** and **x2** with specified correlation **r**.


{% highlight r %}
bivariate_normal2 <- function(r,n) {

Cov_matrix <- matrix(c(1,r,r,1), nrow = 2)
create_data <-rmvnorm(mean = c(0,0), sig = Cov_matrix, n = n)
x1 <- create_data[,1]
x2 <-create_data[,2]
create_data <- data.frame(x1,x2)
return(create_data)

}
{% endhighlight %}

From above function we are generating **n=100** observations of **x1** and **x2** with correlation **0.5**.


{% highlight r %}
set.seed(12)
df <- bivariate_normal2(0.5,100)
cor(df$x1,df$x2)
{% endhighlight %}



{% highlight text %}
## [1] 0.4898011
{% endhighlight %}

Now generating y with $$ \beta_{0}=2 $$,$$ \beta_{1}=-1 $$,$$ \beta_{2}=0.5 $$ and error $$ e_{i} $$ with mean 0 and standard deviation 1.


{% highlight r %}
set.seed(1252)
error <- rnorm(100,0,1)
y <- 2 + (-1)*(df$x1) + 0.5*(df$x2) + error
df <- mutate(df,y=y)
{% endhighlight %}
Now we are ready with the data


{% highlight r %}
kable(summary(df))
{% endhighlight %}



|   |      x1         |      x2         |      y        |
|:--|:----------------|:----------------|:--------------|
|   |Min.   :-2.00005 |Min.   :-1.91294 |Min.   :-1.351 |
|   |1st Qu.:-0.71895 |1st Qu.:-0.63666 |1st Qu.: 1.339 |
|   |Median :-0.06823 |Median :-0.03420 |Median : 2.185 |
|   |Mean   :-0.05472 |Mean   : 0.02842 |Mean   : 2.106 |
|   |3rd Qu.: 0.57552 |3rd Qu.: 0.70633 |3rd Qu.: 3.047 |
|   |Max.   : 1.99234 |Max.   : 2.17587 |Max.   : 4.368 |


### **Scenarios:** 
#### **1. No missing data**

#### what if data did not contain any missing value? what will be the our estimates from fitted regression model?


{% highlight r %}
Model1 <- lm(y ~ x1 +x2, data=df)
kable(tidy(Model1, conf.int = T))
{% endhighlight %}



|term        |   estimate| std.error| statistic| p.value|  conf.low|  conf.high|
|:-----------|----------:|---------:|---------:|-------:|---------:|----------:|
|(Intercept) |  2.0442770| 0.1044346| 19.574713|   0e+00|  1.837003|  2.2515507|
|x1          | -0.7768768| 0.1315464| -5.905723|   1e-07| -1.037960| -0.5157935|
|x2          |  0.6613673| 0.1265442|  5.226373|   1e-06|  0.410212|  0.9125226|


Model provides estimate of $$ \beta_{0}=2.04 $$,$$ \beta_{1}=-0.77 $$,$$ \beta_{2}=0.66 $$. Keep note of standard errors and confidence intervals.We will require it later stage when we asses impact of missing data.  

What we see here estimates are not exactly the same from which we simulated the data. It is random variation.

**Let see what will be the estimates if we simulate above dataset 1000 times and fit model to each dataset.**


{% highlight r %}
Models <- list()
beta_0 <- vector()
beta_1 <- vector()
beta_2 <- vector()


for ( i in 1:1000) {
  
      dfs <- bivariate_normal2(0.5,100)
      error <- rnorm(100,0,1)
      dfs <- mutate(dfs,y=2 -1*dfs$x1 +0.5*dfs$x2 + error)
      Models[[i]] <- lm(y~x1+x2,data = dfs)
      beta_0[i] <- Models[[i]]$coef[1]
      beta_1[i] <- Models[[i]]$coef[2]
      beta_2[i] <- Models[[i]]$coef[3]
  
}
{% endhighlight %}



{% highlight r %}
results <-matrix(c(mean(beta_0),mean(beta_1),mean(beta_2),quantile(beta_0,probs = 0.025),quantile(beta_1,probs = 0.025),quantile(beta_2,probs = 0.025),quantile(beta_0,probs = 0.975),quantile(beta_1,probs = 0.975),quantile(beta_2,probs = 0.975)),nrow = 3,ncol=3)
dimnames(results) <-list(c("beta_0","beta_1","beta_2"),c("mean","CI-L","CI-U"))
{% endhighlight %}

{% highlight r %}
kable(results)
{% endhighlight %}



|       |      mean|       CI-L|       CI-U|
|:------|---------:|----------:|----------:|
|beta_0 |  1.999483|  1.8028743|  2.2064731|
|beta_1 | -1.003165| -1.2252990| -0.7644327|
|beta_2 |  0.502960|  0.2724443|  0.7253348|


As we see, estimates are equivalent to $$ \beta_{0}=2 $$,$$ \beta_{1}=-1 $$,$$ \beta_{2}=0.5 $$ through which we had fitted model.

#### **2. MCAR (Missing Completely At Random) dataset**

#### **what if data contain missing values as per MCAR mechnism?**

In this scenario, we first Make each observation of **x1** missing completely at random (MCAR) with **probability 0.5** and then fit the linear model to this reduced data set.


{% highlight r %}
MCAR_df <- data.frame(df)
set.seed(873)
mcar_x1 <- rbinom(100,1,prob = 0.5)
MCAR_df$x1[mcar_x1 == 1] <- NA
kable(summary(MCAR_df))
{% endhighlight %}



|   |      x1         |      x2         |      y        |
|:--|:----------------|:----------------|:--------------|
|   |Min.   :-2.00005 |Min.   :-1.91294 |Min.   :-1.351 |
|   |1st Qu.:-0.74991 |1st Qu.:-0.63666 |1st Qu.: 1.339 |
|   |Median : 0.14131 |Median :-0.03420 |Median : 2.185 |
|   |Mean   : 0.03105 |Mean   : 0.02842 |Mean   : 2.106 |
|   |3rd Qu.: 0.86724 |3rd Qu.: 0.70633 |3rd Qu.: 3.047 |
|   |Max.   : 1.99234 |Max.   : 2.17587 |Max.   : 4.368 |
|   |NA's   :43       |NA               |NA             |


Missingness mechanism in this dataset is missing completely at random(MCAR) because the missingness of observations **X1** is not depends on other observed (and unobserved) variables data. MCAR data have **less power** because estimation based on reduced sample size.


{% highlight r %}
Model2 <- lm(y ~ x1 +x2, data=MCAR_df)
kable(tidy(Model2, conf.int=T))
{% endhighlight %}



|term        |   estimate| std.error| statistic|   p.value|   conf.low|  conf.high|
|:-----------|----------:|---------:|---------:|---------:|----------:|----------:|
|(Intercept) |  1.9901913| 0.1420728| 14.008245| 0.0000000|  1.7053524|  2.2750302|
|x1          | -0.8147650| 0.1630754| -4.996245| 0.0000065| -1.1417116| -0.4878184|
|x2          |  0.6754888| 0.1659180|  4.071221| 0.0001538|  0.3428433|  1.0081343|


"lm" command in R automatically remove entire row with missing values (mean only includes complete cases).You can see, compare to previous scenario, now **standard errors are high** and **high uncertainty in estimates(wider confidence Interval)**.But still there is no evidence of estimates beings biased.

**Let see what will be the estimates if we simulate above dataset 1000 times and fit model to each dataset.**


{% highlight r %}
Models <- list()
beta_0 <- vector()
beta_1 <- vector()
beta_2 <- vector()


for ( i in 1:1000) {
  
      dfs <- bivariate_normal2(0.5,100)
      error <- rnorm(100,0,1)
      dfs <- mutate(dfs,y=2 -1*dfs$x1 +0.5*dfs$x2 + error)
      dfs <- mutate(dfs,MCAR = rbinom(100,1,prob = 0.5))
      
      dfs <- mutate(dfs,x1 = ifelse(MCAR==1,x1,NA))
      
      Models[[i]] <- lm(y~x1+x2,data = dfs)
      beta_0[i] <- Models[[i]]$coef[1]
      beta_1[i] <- Models[[i]]$coef[2]
      beta_2[i] <- Models[[i]]$coef[3]
  
}
{% endhighlight %}



{% highlight r %}
results <-matrix(c(mean(beta_0),mean(beta_1),mean(beta_2),quantile(beta_0,probs = 0.025),quantile(beta_1,probs = 0.025),quantile(beta_2,probs = 0.025),quantile(beta_0,probs = 0.975),quantile(beta_1,probs = 0.975),quantile(beta_2,probs = 0.975)),nrow = 3,ncol=3)
dimnames(results) <-list(c("beta_0","beta_1","beta_2"),c("mean","CI-L","CI-U"))
{% endhighlight %}

{% highlight r %}
kable(results)
{% endhighlight %}



|       |       mean|       CI-L|       CI-U|
|:------|----------:|----------:|----------:|
|beta_0 |  2.0081838|  1.7283638|  2.2844832|
|beta_1 | -1.0069078| -1.3304568| -0.6780697|
|beta_2 |  0.4956468|  0.1549021|  0.8118329|


Simulating this MCAR dataset 1000 times and fitting linear model to each dataset, shows that estimates are similler from which we have simulated data even there is presence  MCAR missing data. So simulation of MCAR dataset many times generally do not give biased estimates.



#### **3. MAR  (Missing At Random) dataset**

#### **what if data contain missing values as per MAR mechanism?**

In this scenario, we first Make each observation of **x1** missing dependent on the value of **x2**. To do this make x1 missing with **probability 0.5** if **x2 < 0**.


{% highlight r %}
MAR_df <- data.frame(df)
MAR_x1_x2 <- vector()
set.seed(4578)
for ( i in 1:100) {
  if (MAR_df$x2[i] < 0) {MAR_x1_x2[i] <- rbinom(1,1,prob=0.5)}
  if (MAR_df$x2[i] > 0) {MAR_x1_x2[i] <-  rbinom(1,1,prob = 1)}
}

MAR_df$x1[MAR_x1_x2 == 1] <- NA
kable(summary(MAR_df))
{% endhighlight %}



|   |      x1        |      x2         |      y        |
|:--|:---------------|:----------------|:--------------|
|   |Min.   :-1.7996 |Min.   :-1.91294 |Min.   :-1.351 |
|   |1st Qu.:-0.8946 |1st Qu.:-0.63666 |1st Qu.: 1.339 |
|   |Median :-0.4672 |Median :-0.03420 |Median : 2.185 |
|   |Mean   :-0.4767 |Mean   : 0.02842 |Mean   : 2.106 |
|   |3rd Qu.: 0.1088 |3rd Qu.: 0.70633 |3rd Qu.: 3.047 |
|   |Max.   : 0.9099 |Max.   : 2.17587 |Max.   : 4.368 |
|   |NA's   :71      |NA               |NA             |


missing values in **x1** which are dependent on value of **X2**. Missingness mechanism in this dataset is missing at random(MAR) because of probability of observation is missing in **X1** is depends on **observed values-(in X2)** in datasets. 


{% highlight r %}
Model3 <- lm(y ~ x1 +x2, data=MAR_df)
kable(tidy(Model3, conf.int=T))
{% endhighlight %}



|term        |   estimate| std.error| statistic|   p.value|   conf.low|  conf.high|
|:-----------|----------:|---------:|---------:|---------:|----------:|----------:|
|(Intercept) |  2.4716109| 0.3663254|  6.747037| 0.0000004|  1.7186183|  3.2246035|
|x1          | -0.7267735| 0.2824330| -2.573260| 0.0161296| -1.3073228| -0.1462241|
|x2          |  0.8901526| 0.3865287|  2.302940| 0.0295373|  0.0956314|  1.6846738|


You can see, compare to MCAR dataset, now standard errors are almost **3 times larger** and **very high uncertainty in estimates(wider confidence Interval)**.But still there is no evidence of estimates beings biased.

**Let see what will be the estimates if we simulate above dataset 1000 times and fit model to each dataset.**


{% highlight r %}
Models <- list()
beta_0 <- vector()
beta_1 <- vector()
beta_2 <- vector()


for ( i in 1:1000) {
  
      dfs <- bivariate_normal2(0.5,100)
      error <- rnorm(100,0,1)
      dfs <- mutate(dfs,y=2 -1*dfs$x1 +0.5*dfs$x2 + error)
      dfs <- mutate(dfs,MAR = ifelse(x2<0,rbinom(1,1,0.5),rbinom(1,1,1)))
      
      dfs <- mutate(dfs,x1 = ifelse(MAR==1,x1,NA))
      
      Models[[i]] <- lm(y~x1+x2,data = dfs)
      beta_0[i] <- Models[[i]]$coef[1]
      beta_1[i] <- Models[[i]]$coef[2]
      beta_2[i] <- Models[[i]]$coef[3]
  
}
{% endhighlight %}


{% highlight r %}
results <-matrix(c(mean(beta_0),mean(beta_1),mean(beta_2),quantile(beta_0,probs = 0.025),quantile(beta_1,probs = 0.025),quantile(beta_2,probs = 0.025),quantile(beta_0,probs = 0.975),quantile(beta_1,probs = 0.975),quantile(beta_2,probs = 0.975)),nrow = 3,ncol=3)
dimnames(results) <-list(c("beta_0","beta_1","beta_2"),c("mean","CI-L","CI-U"))
{% endhighlight %}

{% highlight r %}
kable(results)
{% endhighlight %}



|       |       mean|       CI-L|       CI-U|
|:------|----------:|----------:|----------:|
|beta_0 |  2.0023796|  1.6050618|  2.4275536|
|beta_1 | -0.9951489| -1.2844987| -0.7192845|
|beta_2 |  0.4976232|  0.0540765|  0.8982080|


#### **4. MNAR  (Missing Not At Random) dataset**

#### **what if data contain missing values as per MNAR mechanism?

In this scenario, we will make the **response variable**, **y**,
missing with **probability 0.75** if **y** is **less than 2**.


{% highlight r %}
MNAR_df <- df
set.seed(3659)
MNAR_y <- vector()
for ( i in 1:100) {
  if (MNAR_df$y[i] < 2){ MNAR_y[i] <- rbinom(1,1,0.75)}
  if (MNAR_df$y[i] > 2){ MNAR_y[i] <- rbinom(1,1,1)}
}
MNAR_df$y[MNAR_y == 1] <- NA
summary(MNAR_df)
{% endhighlight %}



{% highlight text %}
##        x1                 x2                 y         
##  Min.   :-2.00005   Min.   :-1.91294   Min.   :0.3776  
##  1st Qu.:-0.71895   1st Qu.:-0.63666   1st Qu.:1.3661  
##  Median :-0.06823   Median :-0.03420   Median :1.4215  
##  Mean   :-0.05472   Mean   : 0.02842   Mean   :1.4552  
##  3rd Qu.: 0.57552   3rd Qu.: 0.70633   3rd Qu.:1.8363  
##  Max.   : 1.99234   Max.   : 2.17587   Max.   :1.9826  
##                                        NA's   :93
{% endhighlight %}

 Missingness mechanism in dataset is missing not at random(MNAR). As word suggest missingness is **neither MCAR nor MAR** (as missingness MNAR depends on unobserved, unrecorded observations)


{% highlight r %}
Model4 <- lm(y ~ x1 +x2, data=MNAR_df)
kable(tidy(Model4, conf.int=T))
{% endhighlight %}



|term        |   estimate| std.error|  statistic|   p.value|   conf.low| conf.high|
|:-----------|----------:|---------:|----------:|---------:|----------:|---------:|
|(Intercept) |  1.5470059| 0.2787311|  5.5501733| 0.0051562|  0.7731243|  2.320887|
|x1          | -0.3078378| 0.6097679| -0.5048443| 0.6402122| -2.0008249|  1.385149|
|x2          |  0.0592018| 0.5198745|  0.1138771| 0.9148221| -1.3842013|  1.502605|


**Let see what will be the estimates if we simulate above dataset 1000 times and fit model to each dataset.**


{% highlight r %}
Models <- list()
beta_0 <- vector()
beta_1 <- vector()
beta_2 <- vector()


for ( i in 1:1000) {
  
      dfs <- bivariate_normal2(0.5,100)
      error <- rnorm(100,0,1)
      dfs <- mutate(dfs,y=2 -1*dfs$x1 +0.5*dfs$x2 + error)
      dfs <- mutate(dfs,MNAR = ifelse(y<2,rbinom(1,1,0.75),rbinom(1,1,1)))
      
      dfs <- mutate(dfs,y = ifelse(MNAR==1,y,NA))
      
      Models[[i]] <- lm(y~x1+x2,data = dfs)
      beta_0[i] <- Models[[i]]$coef[1]
      beta_1[i] <- Models[[i]]$coef[2]
      beta_2[i] <- Models[[i]]$coef[3]
  
}
{% endhighlight %}


{% highlight r %}
results <-matrix(c(mean(beta_0),mean(beta_1),mean(beta_2),quantile(beta_0,probs = 0.025),quantile(beta_1,probs = 0.025),quantile(beta_2,probs = 0.025),quantile(beta_0,probs = 0.975),quantile(beta_1,probs = 0.975),quantile(beta_2,probs = 0.975)),nrow = 3,ncol=3)
dimnames(results) <-list(c("beta_0","beta_1","beta_2"),c("mean","CI-L","CI-U"))
{% endhighlight %}


{% highlight r %}
kable(results)
{% endhighlight %}



|       |       mean|       CI-L|       CI-U|
|:------|----------:|----------:|----------:|
|beta_0 |  2.2257260|  1.8210183|  2.9458047|
|beta_1 | -0.8640819| -1.2026538| -0.3032271|
|beta_2 |  0.4290430|  0.0826673|  0.7119431|


As you see now, even simulation did not give original estimates- **did not remove bias**. That has relation with **randomisation**. If sample are **not random** ( random when each individual in population has equal chance of selection in sample), there will always possibility of bias which even **cant be corrected by simulation**. Lesson learned !
