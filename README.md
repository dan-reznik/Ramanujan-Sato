Computing π w/ Ramanujan-Sato
================

One of
[Ramanujan-Sato](https://en.wikipedia.org/wiki/Ramanujan%E2%80%93Sato_series)
Formulas for π (1917):

<img src="pics/ramanujan-sato.png" width="50%" style="display: block; margin: auto;" />

Consider the \(396^(4k)\) term in the series. Suppose k=4. R 53-bit
doubles can represent up to 22 significant digits:

``` r
print(396^(4*4),digits=22)
#> [1] 3.6569832807775449e+41
```

Arbitrary-precision numbers are available with the excellent
[Rmpfr](https://cran.r-project.org/web/packages/Rmpfr/vignettes/Rmpfr-pkg.pdf)
package. With 120-bits, all the 42 digits of the above number are shown:

``` r
Rmpfr::mpfr(396, 120)^(4*4)
#> 1 'mpfr' number of precision  120   bits 
#> [1] 365698328077754498546241794891999342493696
```

-----

Load the arbitrary-precision library

``` r
library(Rmpfr) # use this for arbitrary-precision floats
```

We will be using 120 bits of representation precision:

``` r
bits <- 120
```

One term of the R-S series and its front coefficient:

``` r
sqrt2 <- sqrt(mpfr(2L, bits))
rs_coeff <- 9801L/(2L*sqrt2)

rs_term <- function(k) {
  num <- factorial(4L*k)*(26390L*k+1103L)
  den <- (factorial(k)*(396L^k))^4L
  num/den
}
```

Calculates π w/ kmax iterations of the R-S formula:

``` r
rs_series <- function(kmax) {
  terms <- rs_term(mpfr(0:(kmax-1),bits))
  (rs_coeff/cumsum(terms))
}
```

Only 5 iterations and we’re pretty close:

``` r
rs_series(5)
#> 5 'mpfr' numbers of precision  120   bits 
#> [1] 3.141592730013305660313996189025215515
#> [2] 3.141592653589793877998905826306013092
#> [3] 3.141592653589793238462649065702758895
#> [4]  3.14159265358979323846264338327955527
#> [5] 3.141592653589793238462643383279502882
```

Error from π computed internally by Rmpfr:

``` r
rs_series(5)-Const("pi",bits)
#> 5 'mpfr' numbers of precision  120   bits 
#> [1]   7.642351242185135280574571263059095384e-8
#> [2]   6.39536262443026510207448352094098835e-16
#> [3]  5.682423256010174379301934938494086872e-24
#> [4]  5.238529448733281520312260003831002306e-32
#> [5] -3.009265538105056020399965535288948935e-36
```

Mind-blowing convergence of π digits through the Ramanujan-Sato series\!
😄

-----

### References:

  - N. Baruah, B. Berndt, H. Chan, “Ramanujan’s Series for 1/π: A
    Survey”, Mathematical Association of America, Aug-Sept 2009.
    [pdf](https://faculty.math.illinois.edu/~berndt/articles/monthly567-587.pdf)
  - Wikipedia, “[Srinivasa
    Ramanujan](https://en.wikipedia.org/wiki/Srinivasa_Ramanujan)”
  - M. Mächler, “Arbitrarily Accurate Computation with R: The Rmpfr
    Package”, ETH Zurich, 2019.
    [pdf](https://cran.r-project.org/web/packages/Rmpfr/vignettes/Rmpfr-pkg.pdf)

-----

© 2019 Dan S. Reznik
