# Non-linear regression model for DNA melting curves.

Melting curve analysis on double-stranded DNA is an assessment of the dissociation-characteristics during heating of the DNA. As the temperature is raised, the double strand begins to dissociate or melt, leading to a rise in the absorbance intensity and a lowering of the fluorescence intensity (Rn). The temperature at which 50% of DNA is denatured is known as the melting point (Tm). I am wondering, would it be possible to reproduce the melting curves as produced by the real-time PCR machine based on the raw data? 




```r
# libraries
library(ggplot2)
```


First, read the raw data from a .csv file

```r
smelt <- read.csv("d2.csv", header = TRUE, sep = ";")
summary(smelt)
```

```
##       temp            D        
##  Min.   :74.5   Min.   : 4968  
##  1st Qu.:79.1   1st Qu.: 5316  
##  Median :83.8   Median : 5726  
##  Mean   :83.8   Mean   :13963  
##  3rd Qu.:88.4   3rd Qu.:25173  
##  Max.   :93.0   Max.   :40212
```

Now determine the nonlinear (weighted) least-squares estimates of the parameters of a nonlinear model.

```r
pop.ss <- nls(D ~ SSlogis(temp, phi1, phi2, phi3), data = smelt)
summary(pop.ss)
```

```
## 
## Formula: D ~ SSlogis(temp, phi1, phi2, phi3)
## 
## Parameters:
##       Estimate Std. Error t value Pr(>|t|)    
## phi1 40003.841   1410.562   28.36  < 2e-16 ***
## phi2    79.330      0.131  606.68  < 2e-16 ***
## phi3    -0.640      0.111   -5.74  2.1e-07 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 4570 on 72 degrees of freedom
## 
## Number of iterations to convergence: 0 
## Achieved convergence tolerance: 5.57e-06
```



```r
psi1 <- coef(pop.ss)[[1]]
psi2 <- coef(pop.ss)[[2]]
psi3 <- coef(pop.ss)[[3]]
```

Add a column where a model is placed. Then calculate the new Rn values: place the estimated paramethers into a function of a non-lineair (sigmoide) regression model. The formula is: 
http://www.theparasitologist.com/wp-content/uploads/2013/11/Screen-Shot-2013-11-07-at-21.09.54-.png
with http://www.theparasitologist.com/wp-content/uploads/2013/11/Screen-Shot-2013-11-07-at-21.10.02-.png
From ‘Nonlinear Regression and Least Squares in R. John Fox & Sanford Weisberg. 2010′.

```r
smelt["modD"] <- NA
smelt$modD <- psi1/(1 + exp(-((smelt$temp - psi2)/psi3)))
```


To calculate the melting point -Tm-, an inverted derivative is needed from the sigmoid curve. The derivative of will look like a like a 'Gaussian curve'. The Tm will be the top of the gaussian curve.  


```r
# making derivative
fo <- function(x) (psi1/(1 + exp(-(x - psi2)/psi3)))
Dfo <- fo
body(Dfo <- deriv(body(fo), "x"))
```

```
## NULL
```

```r
Dfo
```

```
## expression({
##     .expr4 <- exp(-(x - psi2)/psi3)
##     .expr5 <- 1 + .expr4
##     .value <- psi1/.expr5
##     .grad <- array(0, c(length(.value), 1L), list(NULL, c("x")))
##     .grad[, "x"] <- psi1 * (.expr4 * (1/psi3))/.expr5^2
##     attr(.value, "gradient") <- .grad
##     .value
## })
```

Now that we know what the derivative looks like, a new column is added for the derivative values. Then the derivative values are calculated. 

```r
smelt["dmodD"] <- NA

expr4 <- exp(-(smelt$temp - psi2)/psi3)
expr5 <- 1 + expr4
value <- psi1/expr5
grad <- -1 * (psi1 * (expr4 * (1/psi3))/expr5^2)
smelt["dmodD"] <- grad
```


Let's see how the model (D_mod) and the derivative of the model (dmodD) fits next to the original (D) curve.

```r
# plotting the stuff
g <- ggplot(smelt, aes(temp)) + geom_line(aes(y = D, colour = "D")) + geom_line(aes(y = modD, 
    colour = "modD")) + geom_line(aes(y = dmodD, colour = "dmodD")) + # theme_bw(base_family = 'Avenir', base_size = 10) +
coord_cartesian() + labs(x = "ramp temperature") + labs(y = expression(paste(Delta(Rn)))) + 
    labs(title = "melting curves")
g
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5.png) 

Now calculate the second derivative and set Rn = 0. The corresponing temperature is the Tm.

```r

fo2 <- function(x) -1 * (psi1 * ((exp(-(x - psi2)/psi3)) * (1/psi3))/(1 + (exp(-(x - 
    psi2)/psi3)))^2)
dDfo <- fo2
body(dDfo <- deriv(body(fo2), "x"))
```

```
## NULL
```

```r
dDfo
```

```
## expression({
##     .expr5 <- exp(-(x - psi2)/psi3)
##     .expr6 <- 1/psi3
##     .expr7 <- .expr5 * .expr6
##     .expr8 <- psi1 * .expr7
##     .expr9 <- 1 + .expr5
##     .expr10 <- .expr9^2
##     .value <- -1 * (.expr8/.expr10)
##     .grad <- array(0, c(length(.value), 1L), list(NULL, c("x")))
##     .grad[, "x"] <- psi1 * (.expr7 * .expr6)/.expr10 - .expr8 * 
##         (2 * (.expr7 * .expr9))/.expr10^2
##     attr(.value, "gradient") <- .grad
##     .value
## })
```

```r

smelt["ddmodD"] <- NA  # make column for 2nd derivative

expr6 <- exp(-(smelt$temp - psi2)/psi3)
expr7 <- 1/psi3
expr8 <- expr6 * expr7
expr9 <- psi1 * expr8
expr10 <- 1 + expr6
expr11 <- expr10^2

gradd <- psi1 * (expr8 * expr7)/expr11 - expr9 * (2 * (expr8 * expr10))/expr11^2
smelt["ddmodD"] <- gradd
```

In this plot the first (dmodD) and the second derivative (ddmodD) are visible. The point where ddmodD crosses the temperuture axis, is the Tm.

```r
# plotting the 1st and 2nd derivatives
g2 <- ggplot(smelt, aes(temp)) + # geom_line(aes(y = D, colour = 'D')) + geom_line(aes(y = modD, colour =
# 'modD')) +
geom_line(aes(y = dmodD, colour = "dmodD")) + geom_line(aes(y = ddmodD, colour = "ddmodD")) + 
    # theme_bw(base_family = 'Avenir', base_size = 10) +
coord_cartesian() + labs(x = "ramp temperature") + labs(y = expression(paste(Delta(Rn)))) + 
    labs(title = "melting curves")
g2
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7.png) 

The exact Tm can be calculated by setting Rn'(temp) to zero. 

