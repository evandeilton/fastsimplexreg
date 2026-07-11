# Print a Simplex Regression Fit

Compactly prints a fitted `"simplex_fast"` object: the formula, links,
number of observations, fit statistics and the mean and dispersion
coefficients.

## Usage

``` r
# S3 method for class 'simplex_fast'
print(x, digits = max(3L, getOption("digits") - 3L), ...)
```

## Arguments

- x:

  An object of class `"simplex_fast"`.

- digits:

  Integer; the number of significant digits to display.

- ...:

  Additional arguments, currently ignored.

## Value

The object `x`, invisibly.

## See also

[`fastsimplexreg()`](https://evandeilton.github.io/fastsimplexreg/reference/fastsimplexreg.md),
[`summary.simplex_fast()`](https://evandeilton.github.io/fastsimplexreg/reference/summary.simplex_fast.md)

## Examples

``` r
set.seed(3)
dat <- data.frame(x1 = rnorm(200))
dat$y <- rsimplex(200, simplex_linkinv(0.3 + 0.5 * dat$x1, "logit"), 1)
fit <- fastsimplexreg(y ~ x1, data = dat, n_threads = 1L)
print(fit)
#> 
#> Fast simplex regression with variable dispersion
#> Formula: y ~ x1
#> <environment: 0x56543f1f6ff0>
#> Mean link: logit | Dispersion link: log 
#> Observations: 200 
#> Log-likelihood: 190.9 
#> AIC: -375.9 
#> BIC:  -366 
#> Convergence: 0 - Converged: relative objective tolerance satisfied. 
#> 
#> Mean coefficients [logit link]:
#> (Intercept)          x1 
#>      0.3648      0.5153 
#> 
#> Dispersion coefficients [log link]:
#> (Intercept) 
#>      -0.156 
```
