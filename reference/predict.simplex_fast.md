# Predictions from a Simplex Regression Fit

Computes predictions from a fitted `"simplex_fast"` model, either on the
data used for fitting or on new data. When `newdata` is supplied, the
stored `terms`, `xlevels` and `contrasts` are reused so that the design
matrices are built consistently with the fit.

## Usage

``` r
# S3 method for class 'simplex_fast'
predict(
  object,
  newdata = NULL,
  type = c("response", "mean", "dispersion", "link", "both"),
  ...
)
```

## Arguments

- object:

  An object of class `"simplex_fast"`.

- newdata:

  Optional `data.frame` of new observations. When `NULL`, the fitted
  values are returned.

- type:

  Type of prediction: `"response"` or `"mean"` (fitted mean \\\mu\\),
  `"dispersion"` (fitted \\\phi\\), `"link"` (a list with the linear
  predictors `mean` and `dispersion`) or `"both"` (a `data.frame` with
  columns `mu` and `phi`).

- ...:

  Additional arguments, currently ignored.

## Value

A numeric vector, a list or a `data.frame`, depending on `type`.

## See also

[`fastsimplexreg()`](https://evandeilton.github.io/fastsimplexreg/reference/fastsimplexreg.md)

## Examples

``` r
set.seed(2)
n <- 300
dat <- data.frame(x1 = rnorm(n), z1 = rnorm(n))
mu <- simplex_linkinv(0.1 + 0.6 * dat$x1, link = "logit")
dat$y <- rsimplex(n, mu, exp(-0.5 + 0.3 * dat$z1))
fit <- fastsimplexreg(y ~ x1 | z1, data = dat, n_threads = 1L)
head(predict(fit, type = "response"))
#> [1] 0.3841861 0.5525140 0.7495474 0.3499761 0.5108834 0.5443204
head(predict(fit, newdata = dat[1:5, ], type = "both"))
#>          mu       phi
#> 1 0.3841861 0.4715160
#> 2 0.5525140 0.4719501
#> 3 0.7495474 0.7181893
#> 4 0.3499761 0.2724608
#> 5 0.5108834 0.6809525
```
