# Diagnostic Plots for a Simplex Mixed-Model Fit

Diagnostic plots built with ggplot2, sharing the panels of
[`plot.simplex_fast()`](https://evandeilton.github.io/fastsimplexreg/reference/plot.simplex_fast.md)
but based on the mixed-model fit (residuals are conditional on the
empirical-Bayes random effects).

## Usage

``` r
# S3 method for class 'simplex_fast_mixed'
plot(
  x,
  which = 1:4,
  type = c("quantile", "deviance", "pearson", "response"),
  smooth = TRUE,
  ...
)
```

## Arguments

- x:

  A fitted `"simplex_fast_mixed"` object.

- which:

  Integer subset of `1:4` selecting panels.

- type:

  Type of residual used in panels 1-3: `"quantile"` (the default;
  randomised quantile residuals), `"deviance"`, `"pearson"` or
  `"response"`. See
  [`residuals.simplex_fast_mixed()`](https://evandeilton.github.io/fastsimplexreg/reference/simplex_fast_mixed-methods.md).

- smooth:

  Logical; add a LOESS smoother.

- ...:

  Additional arguments, currently ignored.

## Value

Invisibly, a `ggplot`/patchwork object or a list of `ggplot`s.

## See also

[`fastsimplexregmixed()`](https://evandeilton.github.io/fastsimplexreg/reference/fastsimplexregmixed.md),
[`plot.simplex_fast()`](https://evandeilton.github.io/fastsimplexreg/reference/plot.simplex_fast.md)

## Examples

``` r
set.seed(1)
J <- 40; nj <- 8; n <- J * nj
dat <- data.frame(g = factor(rep(seq_len(J), each = nj)), x1 = rnorm(n))
b <- rnorm(J, 0, 0.7)[dat$g]
dat$y <- rsimplex(n, simplex_linkinv(0.3 - 0.6 * dat$x1 + b, "logit"), 1)
fit <- fastsimplexregmixed(y ~ x1, random = ~ 1 | g, data = dat,
                           nAGQ = 7, n_threads = 1)
if (requireNamespace("ggplot2", quietly = TRUE)) p <- plot(fit, which = 1:2)
```
