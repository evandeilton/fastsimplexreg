# Extractor Methods for Simplex Mixed-Model Fits

Standard extractor methods for objects of class `"simplex_fast_mixed"`
produced by
[`fastsimplexregmixed()`](https://evandeilton.github.io/fastsimplexreg/reference/fastsimplexregmixed.md).

- `coef`:

  Fixed-effect coefficients (`model = "all"`, `"mean"` or
  `"dispersion"`). Random effects are obtained with
  [`ranef()`](https://rdrr.io/pkg/nlme/man/random.effects.html).

- `vcov`:

  Covariance matrix of the estimated parameters `c(beta, gamma, omega)`.

- `logLik`:

  Maximised marginal log-likelihood, with `df = p + r + q(q+1)/2` and
  `nobs`.

- `nobs`:

  Number of observations.

- `ngrps`:

  Number of groups.

- `fitted`:

  Fitted means (conditional on the empirical-Bayes random effects) or
  fitted dispersions.

- `residuals`:

  Response, Pearson or deviance residuals, conditional on the
  empirical-Bayes random effects.

- `ranef`:

  Empirical-Bayes random-effect modes (a groups-by-`q` matrix); with
  `postVar = TRUE`, the posterior covariances are attached as the
  `"postVar"` attribute.

- `VarCorr`:

  The estimated random-effect covariance matrix \\\Sigma\\, with
  standard deviations and correlations.

## Usage

``` r
# S3 method for class 'simplex_fast_mixed'
coef(object, model = c("all", "mean", "dispersion"), ...)

# S3 method for class 'simplex_fast_mixed'
vcov(object, ...)

# S3 method for class 'simplex_fast_mixed'
logLik(object, ...)

# S3 method for class 'simplex_fast_mixed'
nobs(object, ...)

# S3 method for class 'simplex_fast_mixed'
fitted(object, model = c("mean", "dispersion"), ...)

# S3 method for class 'simplex_fast_mixed'
residuals(object, type = c("response", "pearson", "deviance"), ...)

# S3 method for class 'simplex_fast_mixed'
ranef(object, postVar = FALSE, ...)

# S3 method for class 'simplex_fast_mixed'
VarCorr(x, sigma = 1, ...)

# S3 method for class 'VarCorr.simplex_fast_mixed'
print(x, digits = max(3L, getOption("digits") - 3L), ...)
```

## Arguments

- object, x:

  A fitted `"simplex_fast_mixed"` object.

- model:

  For `coef`, one of `"all"`, `"mean"` or `"dispersion"`; for `fitted`,
  one of `"mean"` or `"dispersion"`.

- ...:

  Additional arguments, currently ignored.

- type:

  For `residuals`, one of `"response"`, `"pearson"` or `"deviance"`.

- postVar:

  For `ranef`, logical; attach posterior covariances.

- sigma:

  For `VarCorr`, an optional scale multiplier (kept for compatibility
  with the generic; defaults to 1).

- digits:

  For the `VarCorr` print method, the number of significant digits to
  display.

## Value

`coef`, `fitted` and `residuals` return numeric vectors; `vcov` returns
a matrix; `ranef` returns a matrix; `VarCorr` returns the covariance
matrix with `stddev`/`correlation` attributes; `logLik` returns a
`"logLik"` object.

## See also

[`fastsimplexregmixed()`](https://evandeilton.github.io/fastsimplexreg/reference/fastsimplexregmixed.md)
