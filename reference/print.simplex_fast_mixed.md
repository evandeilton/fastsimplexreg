# Print a Simplex Mixed-Model Fit

Print a Simplex Mixed-Model Fit

## Usage

``` r
# S3 method for class 'simplex_fast_mixed'
print(x, digits = max(3L, getOption("digits") - 3L), ...)
```

## Arguments

- x:

  A fitted `"simplex_fast_mixed"` object.

- digits:

  Number of significant digits.

- ...:

  Additional arguments, currently ignored.

## Value

The object `x`, invisibly.

## See also

[`fastsimplexregmixed()`](https://evandeilton.github.io/fastsimplexreg/reference/fastsimplexregmixed.md),
[`summary.simplex_fast_mixed()`](https://evandeilton.github.io/fastsimplexreg/reference/summary.simplex_fast_mixed.md)

## Examples

``` r
set.seed(1)
J <- 40; nj <- 8; n <- J * nj
dat <- data.frame(g = factor(rep(seq_len(J), each = nj)), x1 = rnorm(n))
b <- rnorm(J, 0, 0.7)[dat$g]
dat$y <- rsimplex(n, simplex_linkinv(0.3 - 0.6 * dat$x1 + b, "logit"), 1)
fit <- fastsimplexregmixed(y ~ x1, random = ~ 1 | g, data = dat,
                           nAGQ = 7, n_threads = 1)
print(fit)
#> 
#> Fast simplex mixed model with variable dispersion
#> Formula: y ~ x1
#> <environment: 0x55692e6b3fe0>
#> Random:  ~1 | g
#> <environment: 0x55692e6b3fe0>
#> Mean link: logit | Dispersion link: log 
#> Observations: 320 | Groups: 40 | nAGQ: 7 
#> Log-likelihood: 265.6 | AIC: -523.1 | BIC: -508.1 
#> 
#> Mean coefficients [logit link]:
#> (Intercept)          x1 
#>      0.3427     -0.6330 
#> 
#> Dispersion coefficients [log link]:
#> (Intercept) 
#>      0.1201 
#> 
#> Random-effect covariance (group: g):
#>             (Intercept)
#> (Intercept)      0.5001
```
