# fastsimplexreg

A high-performance implementation of simplex regression for continuous
proportions in the open interval $`(0, 1)`$, with separate submodels for
the mean and the dispersion, and both fixed- and mixed-effects fitting.
The numerical core (log-likelihood, analytic score, native BFGS,
adaptive Gauss-Hermite quadrature, density, random generation,
prediction) is written in C++ with `RcppArmadillo`, `BLAS`/`LAPACK` and
optional `OpenMP`. Documentation website:
<https://evandeilton.github.io/fastsimplexreg/>.

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
\exp\!\left\{-\frac{1}{2\phi}\,
\frac{(y-\mu)^2}{y(1-y)\,\mu^2(1-\mu)^2}\right\}, \qquad 0 < y < 1.
```

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

``` r

library(fastsimplexreg)

set.seed(20260710)
n <- 2000
dat <- data.frame(
  x1 = rnorm(n),
  x2 = rbinom(n, size = 1, prob = 0.35),
  z1 = rnorm(n)
)

mu_true  <- simplex_linkinv(-0.30 + 0.90 * dat$x1 - 0.55 * dat$x2, link = "logit")
phi_true <- exp(-1.10 + 0.65 * dat$z1)
dat$y <- rsimplex(n, mu = mu_true, phi = phi_true)

fit <- fastsimplexreg(y ~ x1 + x2 | z1, data = dat, link = "logit", n_threads = 1L)
summary(fit)
#> 
#> Call:
#> fastsimplexreg(formula = y ~ x1 + x2 | z1, data = dat, link = "logit", 
#>     n_threads = 1L)
#> 
#> Pearson residuals:
#>      Min       1Q   Median       3Q      Max 
#> -3.45204 -0.68367 -0.04273  0.64817  3.69627 
#> 
#> Coefficients (mean model with logit link):
#>              Estimate Std. Error z value Pr(>|z|)    
#> (Intercept) -0.295723   0.006053  -48.85   <2e-16 ***
#> x1           0.906999   0.003760  241.24   <2e-16 ***
#> x2          -0.539969   0.009846  -54.84   <2e-16 ***
#> 
#> Coefficients (dispersion model with log link):
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept) -1.13841    0.03162  -36.00   <2e-16 ***
#> z1           0.67458    0.03161   21.34   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Log-likelihood:  3298 | AIC: -6585 | BIC: -6557 
#> Deviance:  2000 | Observations: 2000 | Iterations: 19 
#> Convergence: 0 - Converged: relative objective tolerance satisfied.
```

Coefficients, the variance-covariance matrix, the log-likelihood and
predictions are available through the usual extractor methods:

``` r

coef(fit)
#> (Intercept)          x1          x2 (Intercept)          z1 
#>  -0.2957225   0.9069989  -0.5399693  -1.1384071   0.6745761
confint(fit)
#>                  2.5 %     97.5 %
#> (Intercept) -0.3075866 -0.2838585
#> x1           0.8996299  0.9143678
#> x2          -0.5592677 -0.5206709
#> (Intercept) -1.2003874 -1.0764269
#> z1           0.6126137  0.7365385
logLik(fit)
#> 'log Lik.' 3297.736 (df=5)
head(predict(fit, type = "both"))
#>          mu       phi
#> 1 0.2530652 0.5034738
#> 2 0.3016079 0.3234770
#> 3 0.4111718 0.4302684
#> 4 0.2548136 0.2273266
#> 5 0.2216869 0.1832636
#> 6 0.2087973 0.1607353
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

## References

Barndorff-Nielsen, O. E. and Jorgensen, B. (1991). Some parametric
models on the simplex. *Journal of Multivariate Analysis*, **39**(1),
106–116.

Zhang, P., Qiu, Z. and Shi, C. (2016). simplexreg: An R Package for
Regression Analysis of Proportional Data Using the Simplex Distribution.
*Journal of Statistical Software*, **71**(11), 1–21.

## License

GPL-3.
