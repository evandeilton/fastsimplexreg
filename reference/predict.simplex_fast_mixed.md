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
