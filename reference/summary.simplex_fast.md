# Summarise a Simplex Regression Fit

Produces a summary of a fitted `"simplex_fast"` object, including
coefficient tables with standard errors, Wald z-statistics and p-values
for both the mean and dispersion submodels.

## Usage

``` r
# S3 method for class 'simplex_fast'
summary(object, ...)

# S3 method for class 'summary.simplex_fast'
print(x, digits = max(3L, getOption("digits") - 3L), ...)
```

## Arguments

- object:

  An object of class `"simplex_fast"`.

- ...:

  Additional arguments, currently ignored.

- x:

  An object of class `"summary.simplex_fast"`.

- digits:

  Integer; the number of significant digits to display.

## Value

An object of class `"summary.simplex_fast"`, a list whose main component
`coefficients` is itself a list with the `mean` and `dispersion`
coefficient tables (each with columns `Estimate`, `Std. Error`,
`z value` and `Pr(>|z|)`), together with the Pearson residuals, the
links, fit statistics (log-likelihood, AIC, BIC, deviance) and optimiser
diagnostics. The `print` method returns its argument invisibly.

## See also

[`fastsimplexreg()`](https://evandeilton.github.io/fastsimplexreg/reference/fastsimplexreg.md),
[`print.simplex_fast()`](https://evandeilton.github.io/fastsimplexreg/reference/print.simplex_fast.md)

## Examples

``` r
set.seed(4)
dat <- data.frame(x1 = rnorm(300), z1 = rnorm(300))
mu <- simplex_linkinv(0.2 + 0.6 * dat$x1, "logit")
dat$y <- rsimplex(300, mu, exp(-0.4 + 0.3 * dat$z1))
fit <- fastsimplexreg(y ~ x1 | z1, data = dat, n_threads = 1L)
summary(fit)
#> 
#> Call:
#> fastsimplexreg(formula = y ~ x1 | z1, data = dat, n_threads = 1L)
#> 
#> Pearson residuals:
#>      Min       1Q   Median       3Q      Max 
#> -2.69348 -0.62167 -0.00168  0.60695  2.53326 
#> 
#> Coefficients (mean model with logit link):
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  0.18499    0.01987   9.309   <2e-16 ***
#> x1           0.60285    0.01927  31.280   <2e-16 ***
#> 
#> Coefficients (dispersion model with log link):
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept) -0.54654    0.08174  -6.686 2.29e-11 ***
#> z1           0.17699    0.08605   2.057   0.0397 *  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Log-likelihood: 335.4 | AIC: -662.9 | BIC:  -648 
#> Deviance:   300 | Observations: 300 | Iterations: 12 
#> Convergence: 0 - Converged: relative objective tolerance satisfied. 
```
