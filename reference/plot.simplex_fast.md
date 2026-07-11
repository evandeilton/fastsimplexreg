# Diagnostic Plots for a Simplex Regression Fit

Produces model-diagnostic plots for a fitted `"simplex_fast"` object
using ggplot2. Up to four panels are available:

1.  Residuals against the fitted mean \\\hat\mu\\.

2.  Normal quantile-quantile plot of the residuals, with a reference
    line through the first and third quartiles.

3.  Scale-location plot (\\\sqrt{\|\text{residual}\|}\\ against
    \\\hat\mu\\).

4.  Observed response against the fitted mean, with the identity line.

## Usage

``` r
# S3 method for class 'simplex_fast'
plot(
  x,
  which = 1:4,
  type = c("deviance", "pearson", "response"),
  smooth = TRUE,
  ...
)
```

## Arguments

- x:

  An object of class `"simplex_fast"`.

- which:

  Integer vector selecting the panels to draw, a subset of `1:4`.

- type:

  Type of residual used in panels 1-3: one of `"deviance"`, `"pearson"`
  or `"response"`. See
  [`residuals.simplex_fast()`](https://evandeilton.github.io/fastsimplexreg/reference/simplex_fast-methods.md).

- smooth:

  Logical; if `TRUE` (the default) a LOESS smoother is added to the
  residual and scale-location panels when there are enough observations.

- ...:

  Additional arguments, currently ignored.

## Value

Invisibly, a single `ggplot` object when one panel is requested, a
combined patchwork object when several panels are requested and
patchwork is available, or a named list of `ggplot` objects otherwise.

## Details

Each panel is a self-contained `ggplot` object. When more than one panel
is requested the panels are combined into a single figure with patchwork
if that package is installed; otherwise they are drawn one after
another. The (list of) `ggplot` object(s) is returned invisibly, so the
plots can be further customised or re-arranged by the caller.

## See also

[`fastsimplexreg()`](https://evandeilton.github.io/fastsimplexreg/reference/fastsimplexreg.md),
[`residuals.simplex_fast()`](https://evandeilton.github.io/fastsimplexreg/reference/simplex_fast-methods.md)

## Examples

``` r
set.seed(6)
n <- 400
dat <- data.frame(x1 = rnorm(n), z1 = rnorm(n))
mu <- simplex_linkinv(0.2 + 0.7 * dat$x1, link = "logit")
dat$y <- rsimplex(n, mu, exp(-0.5 + 0.4 * dat$z1))
fit <- fastsimplexreg(y ~ x1 | z1, data = dat, n_threads = 1L)
p <- plot(fit, which = 1:4)

# 'p' is a ggplot/patchwork object that can be further customised.
```
