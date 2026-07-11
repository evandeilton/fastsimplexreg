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
  type = c("deviance", "pearson", "response"),
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

  Residual type used in the panels.

- smooth:

  Logical; add a LOESS smoother.

- ...:

  Additional arguments, currently ignored.

## Value

Invisibly, a `ggplot`/patchwork object or a list of `ggplot`s.

## See also

[`fastsimplexregmixed()`](https://evandeilton.github.io/fastsimplexreg/reference/fastsimplexregmixed.md),
[`plot.simplex_fast()`](https://evandeilton.github.io/fastsimplexreg/reference/plot.simplex_fast.md)
