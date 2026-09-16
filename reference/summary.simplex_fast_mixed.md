# Summarise a Simplex Mixed-Model Fit

Produces a summary of a fitted `"simplex_fast_mixed"` object, including
coefficient tables with standard errors, Wald z-statistics and p-values
for the fixed-effect mean and dispersion submodels, together with the
estimated random-effect covariance.

## Usage

``` r
# S3 method for class 'simplex_fast_mixed'
summary(object, ...)

# S3 method for class 'summary.simplex_fast_mixed'
print(x, digits = max(3L, getOption("digits") - 3L), ...)
```

## Arguments

- object:

  A fitted `"simplex_fast_mixed"` object.

- ...:

  Additional arguments, currently ignored.

- x:

  A `"summary.simplex_fast_mixed"` object.

- digits:

  Number of significant digits.

## Value

An object of class `"summary.simplex_fast_mixed"`.

## See also

[`fastsimplexregmixed()`](https://evandeilton.github.io/fastsimplexreg/reference/fastsimplexregmixed.md)

## Examples

``` r
set.seed(1)
J <- 40; nj <- 8; n <- J * nj
dat <- data.frame(g = factor(rep(seq_len(J), each = nj)), x1 = rnorm(n))
b <- rnorm(J, 0, 0.7)[dat$g]
dat$y <- rsimplex(n, simplex_linkinv(0.3 - 0.6 * dat$x1 + b, "logit"), 1)
fit <- fastsimplexregmixed(y ~ x1, random = ~ 1 | g, data = dat,
                           nAGQ = 7, n_threads = 1)
summary(fit)
#> 
#> Call:
#> fastsimplexregmixed(formula = y ~ x1, data = dat, random = ~1 | 
#>     g, nAGQ = 7, n_threads = 1)
#> 
#> Quantile residuals:
#>      Min       1Q   Median       3Q      Max 
#> -2.64067 -0.61871  0.04348  0.60741  3.07652 
#> 
#> Coefficients (mean model with logit link):
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  0.34265    0.11449   2.993  0.00276 ** 
#> x1          -0.63305    0.02601 -24.343  < 2e-16 ***
#> 
#> Coefficients (dispersion model with log link):
#>             Estimate Std. Error z value Pr(>|z|)
#> (Intercept)  0.12005    0.08462   1.419    0.156
#> 
#> Random effects:
#> Random effects covariance (group: g)
#>             Variance Std.Dev.
#> (Intercept)   0.5001   0.7072
#> 
#> Log-likelihood: 265.6 | AIC: -523.1 | BIC: -508.1 
#> Observations: 320 | Groups: 40 | nAGQ: 7 | Iterations: 14 
#> Convergence: 0 - Converged: relative objective tolerance satisfied. 
```
