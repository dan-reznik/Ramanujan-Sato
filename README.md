Computing Pi w/ Ramanujan-Sato
================

One of
[Ramanujan-Sato](https://en.wikipedia.org/wiki/Ramanujan%E2%80%93Sato_series)
Formulas for \(pi\) (1917):

<img src="pics/ramanujan-sato.png" width="50%" style="display: block; margin: auto;" />

To correctly calculate this series we will need floating point precision
far superior to R’s defaults. We will be using the `Rmpfr` package.

-----

Load libraries

``` r
library(tidyverse)
library(Rmpfr) # use this for arbitrary-precision floats
source("util.R")
```

We will be using 120 bits of representation precision.

``` r
bits <- 120
```

One term of the R-S series and its front coefficient

``` r
sqrt2 <- sqrt(mpfr(2L, bits))
rs_coeff <- 9801L/(2L*sqrt2)

rs_term <- function(k) {
  (factorial(4L*k)*(26390L*k+1103L))/
  ((factorial(k)*(396L^k))^4L)
}
```

Calculates \(pi\) from 1 to kmax+1 iterations:

``` r
rs_cumsum <- function(kmax) {
  terms <- rs_term(mpfr(0:kmax,bits))
  (rs_coeff/cumsum(terms))
}
```

Four iterations places the number pretty close to \(pi\)

``` r
rs_cumsum(4)
```

    ## 5 'mpfr' numbers of precision  120   bits 
    ## [1] 3.141592730013305660313996189025215515
    ## [2] 3.141592653589793877998905826306013092
    ## [3] 3.141592653589793238462649065702758895
    ## [4]  3.14159265358979323846264338327955527
    ## [5] 3.141592653589793238462643383279502882

Difference from \(pi\) computed internally by Rmpfr:

``` r
rs_cumsum(4)-Const("pi",bits)
```

    ## 5 'mpfr' numbers of precision  120   bits 
    ## [1]   7.642351242185135280574571263059095384e-8
    ## [2]   6.39536262443026510207448352094098835e-16
    ## [3]  5.682423256010174379301934938494086872e-24
    ## [4]  5.238529448733281520312260003831002306e-32
    ## [5] -3.009265538105056020399965535288948935e-36

Ramanujan-Sato convergence is astounding\! 😄

-----

Dan S. Reznik, May 2019
