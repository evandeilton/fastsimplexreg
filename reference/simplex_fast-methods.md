# Extractor Methods for Simplex Regression Fits

Standard extractor methods for objects of class `"simplex_fast"`
produced by
[`fastsimplexreg()`](https://evandeilton.github.io/fastsimplexreg/reference/fastsimplexreg.md).
They mirror the conventions of the corresponding generics for linear and
generalised linear models.

- `coef`:

  Returns the estimated coefficients. `model = "all"` returns the full
  vector, `"mean"` the mean submodel and `"dispersion"` the dispersion
  submodel.

- `vcov`:

  Returns the variance-covariance matrix of the estimates. Errors if the
  model was fitted with `inference = FALSE`.

- `logLik`:

  Returns the maximised log-likelihood, with attributes `df` (number of
  estimated parameters) and `nobs`, and class `"logLik"`.

- `nobs`:

  Returns the number of observations used in the fit.

- `fitted`:

  Returns the fitted means (`model = "mean"`) or fitted dispersions
  (`model = "dispersion"`).

- `residuals`:

  Returns residuals. `type = "response"` gives \\y - \hat\mu\\;
  `type = "pearson"` gives \\(y - \hat\mu) / \sqrt{\hat\phi\\
  V(\hat\mu)}\\ with the simplex unit variance function \\V(\mu) =
  \\\mu(1-\mu)\\^3\\, i.e. the first-order dispersion-model
  approximation \\\mathrm{Var}(Y) \approx \phi\\V(\mu)\\;
  `type = "deviance"` gives the signed deviance residual
  \\\mathrm{sign}(y-\hat\mu)\sqrt{d(y;\hat\mu)/\hat\phi}\\ based on the
  simplex unit deviance \\d(y;\mu) = (y-\mu)^2 /
  \\y(1-y)\mu^2(1-\mu)^2\\\\.

- `confint`:

  Returns Wald confidence intervals \\\hat\theta \pm
  z\_{1-\alpha/2}\\\mathrm{SE}(\hat\theta)\\, using the covariance
  matrix from `vcov`.

- `deviance`:

  Returns the unscaled deviance \\\sum_i d(y_i;\hat\mu_i)\\.
  `type = "scaled"` gives \\\sum_i d(y_i;\hat\mu_i)/\hat\phi_i\\, the
  sum of squared deviance residuals – but note that this is identically
  `nobs(object)` whenever the dispersion submodel contains an intercept,
  since that is what the score equation for \\\gamma\\ forces, so it
  cannot distinguish two models.

- `model.matrix`:

  Returns the mean (`model = "mean"`) or dispersion
  (`model = "dispersion"`) design matrix. Requires the fit to have
  stored the design (`x = TRUE`) or the model frame (`model = TRUE`).

- `terms`:

  Returns the `"terms"` object for the mean or dispersion submodel.

- `formula`:

  Returns the (multi-part) model formula.

- `model.frame`:

  Returns the stored model frame (requires `model = TRUE` at fitting
  time).

- `update`:

  Refits the model with a modified formula and/or arguments. The
  multi-part (`|`) structure is updated correctly through the Formula
  package.

## Usage

``` r
# S3 method for class 'simplex_fast'
coef(object, model = c("all", "mean", "dispersion"), ...)

# S3 method for class 'simplex_fast'
vcov(object, ...)

# S3 method for class 'simplex_fast'
logLik(object, ...)

# S3 method for class 'simplex_fast'
nobs(object, ...)

# S3 method for class 'simplex_fast'
fitted(object, model = c("mean", "dispersion"), ...)

# S3 method for class 'simplex_fast'
residuals(object, type = c("quantile", "response", "pearson", "deviance"), ...)

# S3 method for class 'simplex_fast'
deviance(object, type = c("unscaled", "scaled"), ...)

# S3 method for class 'simplex_fast'
model.matrix(object, model = c("mean", "dispersion"), ...)

# S3 method for class 'simplex_fast'
terms(x, model = c("mean", "dispersion"), ...)

# S3 method for class 'simplex_fast'
formula(x, ...)

# S3 method for class 'simplex_fast'
model.frame(formula, ...)

# S3 method for class 'simplex_fast'
update(object, formula., ..., evaluate = TRUE)

# S3 method for class 'simplex_fast'
confint(object, parm, level = 0.95, ...)
```

## Arguments

- object, x, formula:

  An object of class `"simplex_fast"` (the argument is named `object`,
  `x` or `formula` to match the corresponding generic).

- model:

  For `coef`, one of `"all"`, `"mean"` or `"dispersion"`; for `fitted`,
  `model.matrix` and `terms`, one of `"mean"` or `"dispersion"`.

- ...:

  Additional arguments, currently ignored.

- type:

  For `residuals`, one of `"quantile"` (the default; randomized quantile
  residuals in the sense of Dunn and Smyth, 1996, which are exactly
  standard normal under a correct model), `"response"`, `"pearson"` or
  `"deviance"`. For `deviance`, `"unscaled"` (default) or `"scaled"`.

- formula.:

  For `update`, a change to the model formula, following the conventions
  of
  [`stats::update.formula()`](https://rdrr.io/r/stats/update.formula.html);
  the two-part mean/dispersion structure is handled via the Formula
  package.

- evaluate:

  For `update`, logical; if `TRUE` (default) the updated call is
  evaluated and the refitted model returned, otherwise the updated call
  is returned unevaluated.

- parm:

  For `confint`, a specification of which parameters to report, either a
  vector of numeric indices or of names. Defaults to all.

- level:

  For `confint`, the confidence level.

## Value

`coef`, `fitted`, `residuals` and `nobs` return numeric vectors; `vcov`,
`confint` and `model.matrix` return matrices; `logLik` returns a
`"logLik"` object; `deviance` returns a single number; `terms`,
`formula` and `model.frame` return the corresponding model-description
objects.

## See also

[`fastsimplexreg()`](https://evandeilton.github.io/fastsimplexreg/reference/fastsimplexreg.md),
[`summary.simplex_fast()`](https://evandeilton.github.io/fastsimplexreg/reference/summary.simplex_fast.md),
[`predict.simplex_fast()`](https://evandeilton.github.io/fastsimplexreg/reference/predict.simplex_fast.md),
[`plot.simplex_fast()`](https://evandeilton.github.io/fastsimplexreg/reference/plot.simplex_fast.md),
[`simulate.simplex_fast()`](https://evandeilton.github.io/fastsimplexreg/reference/simulate.simplex_fast.md)

## Examples

``` r
set.seed(1)
n <- 300
dat <- data.frame(x1 = rnorm(n), z1 = rnorm(n))
mu <- simplex_linkinv(0.2 + 0.7 * dat$x1, link = "logit")
dat$y <- rsimplex(n, mu, exp(-0.5 + 0.4 * dat$z1))
fit <- fastsimplexreg(y ~ x1 | z1, data = dat, n_threads = 1L, x = TRUE)
coef(fit)
#>       (Intercept)                x1 (phi)_(Intercept)          (phi)_z1 
#>         0.2143248         0.6415763        -0.5119182         0.3018420 
coef(fit, model = "mean")
#> (Intercept)          x1 
#>   0.2143248   0.6415763 
vcov(fit)
#>                     (Intercept)            x1 (phi)_(Intercept)      (phi)_z1
#> (Intercept)        3.967455e-04 -3.906049e-05     -2.114781e-06 -1.979147e-04
#> x1                -3.906049e-05  3.371292e-04      4.858320e-07  4.546636e-05
#> (phi)_(Intercept) -2.114781e-06  4.858320e-07      6.667492e-03  7.725782e-05
#> (phi)_z1          -1.979147e-04  4.546636e-05      7.725782e-05  7.230704e-03
logLik(fit)
#> 'log Lik.' 338.8326 (df=4)
AIC(fit)
#> [1] -669.6653
nobs(fit)
#> [1] 300
deviance(fit)
#> [1] 186.9003
head(fitted(fit))
#>         1         2         3         4         5         6 
#> 0.4532388 0.5822815 0.4202384 0.7751842 0.6048534 0.4226100 
head(residuals(fit, type = "deviance"))
#>          1          2          3          4          5          6 
#>  0.5268737 -0.1993219  0.1995191  0.6009331  0.7250419 -0.3582407 
confint(fit)
#>                        2.5 %     97.5 %
#> (Intercept)        0.1752854  0.2533643
#> x1                 0.6055893  0.6775634
#> (phi)_(Intercept) -0.6719585 -0.3518779
#> (phi)_z1           0.1351793  0.4685047
head(model.matrix(fit, model = "mean"))
#>   (Intercept)         x1
#> 1           1 -0.6264538
#> 2           1  0.1836433
#> 3           1 -0.8356286
#> 4           1  1.5952808
#> 5           1  0.3295078
#> 6           1 -0.8204684
formula(fit)
#> y ~ x1 | z1
#> <environment: 0x55df04b51528>
```
