# fastsimplexreg

A high-performance implementation of simplex regression for continuous
proportions in the open interval $`(0, 1)`$, with separate submodels for
the mean and the dispersion, and both fixed- and mixed-effects fitting.
The numerical core (log-likelihood, analytic score, native BFGS,
adaptive Gauss-Hermite quadrature, density, random generation,
prediction) is written in C++ with `RcppArmadillo`, `BLAS`/`LAPACK` and
optional `OpenMP`.

`fastsimplexreg` provides high-performance maximum-likelihood estimation
of **simplex regression** models for continuous proportions in the open
interval $`(0, 1)`$, following Barndorff-Nielsen and Jorgensen (1991).
It fits separate submodels for the mean and the dispersion,

``` math
Y_i \sim \mathrm{Simplex}(\mu_i, \phi_i), \qquad
g(\mu_i) = \mathbf{x}_i^\top \boldsymbol\beta, \qquad
\log(\phi_i) = \mathbf{z}_i^\top \boldsymbol\gamma,
```

with the density

``` math
f(y; \mu, \phi) = \left[2\pi\phi\,(y(1-y))^3\right]^{-1/2}
\exp\!\left\lbrace -\frac{1}{2\phi}\,
\frac{(y-\mu)^2}{y(1-y)\,\mu^2(1-\mu)^2}\right\rbrace, \qquad 0 < y < 1.
```

For nested or clustered proportion data, the package also fits **simplex
mixed-effects models** with cluster-specific random effects
$`\mathbf{b}_j \sim N_q(\mathbf{0}, \boldsymbol\Sigma)`$ in the mean
submodel:

``` math
Y_{ij} \mid \mathbf{b}_j \sim \mathrm{Simplex}(\mu_{ij}, \phi_{ij}), \qquad
g(\mu_{ij}) = \mathbf{x}_{ij}^\top \boldsymbol\beta + \mathbf{z}_{ij}^\top \mathbf{b}_j, \qquad
\log(\phi_{ij}) = \mathbf{w}_{ij}^\top \boldsymbol\gamma.
```

Since the random effects are unobserved, estimation is performed by
maximizing the marginal likelihood:

``` math
L(\boldsymbol\beta, \boldsymbol\gamma, \boldsymbol\Sigma) = \prod_j \int_{\mathbb{R}^q} \left[ \prod_i f(y_{ij} \mid \mathbf{b}_j; \mu_{ij}, \phi_{ij}) \right] f(\mathbf{b}_j; \boldsymbol\Sigma) \, \mathrm{d}\mathbf{b}_j,
```

where $`f(\mathbf{b}_j; \boldsymbol\Sigma)`$ is the multivariate normal
density. This high-dimensional integral is approximated using
high-performance **adaptive Gauss-Hermite quadrature** (AGHQ).

The entire numerical hot path – log-likelihood, analytic score, a native
BFGS optimiser, density, random generation, prediction and link inverses
– is implemented in C++ with `RcppArmadillo`, BLAS/LAPACK and optional
OpenMP parallelism, so that models scale to large data sets.

## Installation

You can install the development version from
[GitHub](https://github.com/evandeilton/fastsimplexreg):

``` r

# install.packages("remotes")
remotes::install_github("evandeilton/fastsimplexreg")
```

## The multi-part formula interface

The API uses the `Formula` package and separates the two submodels with
the `|` operator:

``` r

fit <- fastsimplexreg(y ~ x1 + x2 | z1 + z2, data = dat, link = "logit")
```

The first right-hand side component models the mean $`\mu`$; the second
models the dispersion $`\phi`$. When the second component is omitted, as
in `y ~ x1 + x2`, the dispersion is constant (equivalent to `| 1`).

## Mean links

The mean supports four links; the dispersion always uses a log link.

| Link      | $`g(\mu)`$               | $`g^{-1}(\eta)`$          |
|-----------|--------------------------|---------------------------|
| `logit`   | $`\log\{\mu/(1-\mu)\}`$  | $`1/(1 + e^{-\eta})`$     |
| `probit`  | $`\Phi^{-1}(\mu)`$       | $`\Phi(\eta)`$            |
| `cloglog` | $`\log\{-\log(1-\mu)\}`$ | $`1 - \exp(-\exp(\eta))`$ |
| `neglog`  | $`-\log\{-\log(\mu)\}`$  | $`\exp(-\exp(-\eta))`$    |

The `neglog` definition follows Zhang et al. (2016).

## A worked example

We use the `ReadingSkills` data from the **betareg** package: the
response `accuracy` is a reading-accuracy score in $`(0, 1)`$ for
children with and without `dyslexia`, together with a standardised `iq`.
We let the mean depend on both and the dispersion depend on dyslexia
status.

``` r

library(fastsimplexreg)
data("ReadingSkills", package = "betareg")

fit <- fastsimplexreg(accuracy ~ dyslexia + iq | dyslexia,
                      data = ReadingSkills, link = "logit")
summary(fit)
#> 
#> Call:
#> fastsimplexreg(formula = accuracy ~ dyslexia + iq | dyslexia, 
#>     data = ReadingSkills, link = "logit")
#> 
#> Pearson residuals:
#>      Min       1Q   Median       3Q      Max 
#> -2.39081 -0.62295  0.24243  0.43805  1.48447 
#> 
#> Coefficients (mean model with logit link):
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  1.37697    0.15352   8.969  < 2e-16 ***
#> dyslexia    -0.97657    0.15485  -6.307 2.85e-10 ***
#> iq          -0.04369    0.07130  -0.613     0.54    
#> 
#> Coefficients (dispersion model with log link):
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)   1.4242     0.2152   6.617 3.66e-11 ***
#> dyslexia     -2.6917     0.2162 -12.450  < 2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Log-likelihood: 68.01 | AIC:  -126 | BIC: -117.1 
#> Deviance:    44 | Observations: 44 | Iterations: 15 
#> Convergence: 0 - Converged: relative objective tolerance satisfied.
```

Coefficients, the variance-covariance matrix, the log-likelihood and
predictions are available through the usual extractor methods:

``` r

coef(fit)
#> (Intercept)    dyslexia          iq (Intercept)    dyslexia 
#>  1.37696520 -0.97656510 -0.04369204  1.42419092 -2.69174257
confint(fit)
#>                  2.5 %      97.5 %
#> (Intercept)  1.0760642  1.67786626
#> dyslexia    -1.2800609 -0.67306929
#> iq          -0.1834433  0.09605926
#> (Intercept)  1.0023482  1.84603360
#> dyslexia    -3.1155107 -2.26797445
logLik(fit)
#> 'log Lik.' 68.00509 (df=5)
head(predict(fit, type = "both"))
#>          mu      phi
#> 1 0.9103076 61.30942
#> 2 0.9111495 61.30942
#> 3 0.9115695 61.30942
#> 4 0.9091703 61.30942
#> 5 0.9155269 61.30942
#> 6 0.9159281 61.30942
```

## Mixed-effects models

For nested or clustered proportion data,
[`fastsimplexregmixed()`](https://evandeilton.github.io/fastsimplexreg/reference/fastsimplexregmixed.md)
adds a cluster-specific random effect, estimated by adaptive
Gauss-Hermite quadrature. The `GasolineYield` data (also from
**betareg**) records the proportion `yield` of crude oil converted to
gasoline across several experimental conditions, in 10 crude-oil
`batch`es. We fit a random intercept per batch:

``` r

data("GasolineYield", package = "betareg")

mfit <- fastsimplexregmixed(yield ~ temp, random = ~ 1 | batch,
                            data = GasolineYield, link = "logit", nAGQ = 15)
summary(mfit)
#> 
#> Call:
#> fastsimplexregmixed(formula = yield ~ temp, data = GasolineYield, 
#>     random = ~1 | batch, link = "logit", nAGQ = 15)
#> 
#> Pearson residuals:
#>      Min       1Q   Median       3Q      Max 
#> -2.00837 -0.50439  0.09397  0.48606  1.32804 
#> 
#> Coefficients (mean model with logit link):
#>               Estimate Std. Error z value Pr(>|z|)    
#> (Intercept) -5.5767677  0.2709674  -20.58   <2e-16 ***
#> temp         0.0119791  0.0006187   19.36   <2e-16 ***
#> 
#> Coefficients (dispersion model with log link):
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  -1.1049     0.3047  -3.626 0.000288 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Random effects:
#> Random effects covariance (group: batch)
#>             Variance Std.Dev.
#> (Intercept)   0.3406   0.5836
#> 
#> Log-likelihood: 53.17 | AIC: -98.34 | BIC: -92.48 
#> Observations: 32 | Groups: 10 | nAGQ: 15 | Iterations: 18 
#> Convergence: 0 - Converged: relative objective tolerance satisfied.
VarCorr(mfit)
#> Random effects covariance (group: batch)
#>             Variance Std.Dev.
#> (Intercept)   0.3406   0.5836
```

## Prediction

The available prediction types are `"response"`/`"mean"`,
`"dispersion"`, `"link"` and `"both"`. Prediction on new data reuses the
`terms`, `xlevels` and `contrasts` stored at fitting time:

``` r

predict(fit, newdata = new_dat, type = "both")
```

## Performance notes

- **Analytic score** for every link, avoiding numerical differentiation
  during optimisation.
- **Native BFGS in C++**, so the optimiser does not cross the R/C++
  boundary at each objective evaluation.
- **Armadillo/BLAS** linear algebra for the linear predictors.
- **OpenMP** parallelism over observations; `n_threads = 0` requests all
  available threads.
- **A single `model.frame`** for both submodels, guaranteeing consistent
  handling of `subset`, `NA`, factors and levels.
- **Optional inference**: pass `inference = FALSE` to skip the Hessian
  for exploratory fits on massive data.

In head-to-head benchmarks, `fastsimplexreg` returns **numerically
identical maximum-likelihood estimates** to the CRAN package
`simplexreg` (fitted-mean correlation = 1) while running roughly **5-13x
faster** for sample sizes up to 1e5, and is ~2-4x faster than `betareg`
while fitting a richer model. See the [benchmark
article](https://evandeilton.github.io/fastsimplexreg/articles/benchmark.html)
(the full reproducible study ships in `inst/benchmark/`).

## References

Barndorff-Nielsen, O. E. and Jorgensen, B. (1991). Some parametric
models on the simplex. *Journal of Multivariate Analysis*, **39**(1),
106–116.

Zhang, P., Qiu, Z. and Shi, C. (2016). simplexreg: An R Package for
Regression Analysis of Proportional Data Using the Simplex Distribution.
*Journal of Statistical Software*, **71**(11), 1–21.

## License

MIT © José Evandeilton Lopes.
