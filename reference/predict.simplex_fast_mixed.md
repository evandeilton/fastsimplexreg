# Predictions from a Simplex Mixed-Model Fit

Predictions from a Simplex Mixed-Model Fit

## Usage

``` r
# S3 method for class 'simplex_fast_mixed'
predict(
  object,
  newdata = NULL,
  type = c("response", "mean", "dispersion", "link", "both"),
  re.form = NULL,
  ...
)
```

## Arguments

- object:

  A fitted `"simplex_fast_mixed"` object.

- newdata:

  Optional new data. When `NULL`, in-sample predictions are returned.

- type:

  Type of prediction: `"response"`/`"mean"`, `"dispersion"`, `"link"` or
  `"both"`.

- re.form:

  Controls the random effects. `NULL` (default) includes the estimated
  random effects for groups seen in the fit; `NA` (or `~0`) gives
  population-level predictions (random effects set to zero).

- ...:

  Additional arguments, currently ignored.

## Value

A numeric vector, list or `data.frame`, depending on `type`.

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
head(predict(fit))
#>         1         2         3         4         5         6 
#> 0.8358421 0.7530194 0.8532132 0.5550628 0.7354481 0.8520072 
head(predict(fit, re.form = NA))
#> [1] 0.6768264 0.5563583 0.7050874 0.3391156 0.5334640 0.7030878
```
