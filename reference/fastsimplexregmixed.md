# Fast Simplex Mixed-Effects Regression with Variable Dispersion

Fits a two-level (nested) simplex mixed model for continuous proportions
in the open interval \\(0, 1)\\ by maximum marginal likelihood, using
adaptive Gauss-Hermite quadrature (AGHQ). Conditional on a cluster
random effect \\b_j \sim N_q(0, \Sigma)\\, \$\$y\_{ij} \mid b_j \sim
\mathrm{Simplex}(\mu\_{ij}, \phi\_{ij}), \quad g(\mu\_{ij}) =
x\_{ij}^\top\beta + z\_{ij}^\top b_j, \quad \log\phi\_{ij} =
w\_{ij}^\top\gamma.\$\$ The fixed-effects mean and dispersion submodels
use the same multi-part Formula interface as
[`fastsimplexreg()`](https://evandeilton.github.io/fastsimplexreg/reference/fastsimplexreg.md);
the random effects and the grouping factor are given through `random`.

## Usage

``` r
fastsimplexregmixed(
  formula,
  data,
  random,
  link = c("logit", "probit", "cloglog", "neglog"),
  nAGQ = 11L,
  start = NULL,
  maxit = 300L,
  rel_tol = 1e-09,
  grad_tol = 1e-06,
  n_threads = 1L,
  inference = TRUE,
  hessian_rel_step = 1e-05,
  inner_maxit = 50L,
  inner_tol = 1e-08,
  trace = FALSE,
  subset = NULL,
  na.action = stats::na.omit,
  model = TRUE,
  x = FALSE,
  y = TRUE
)
```

## Arguments

- formula:

  A multi-part formula `y ~ mean_terms | dispersion_terms`. When the
  dispersion part is omitted the dispersion is constant.

- data:

  A `data.frame` containing the response, covariates and grouping
  factor.

- random:

  A one-sided formula giving the random-effects design and the grouping
  factor, `~ z1 + z2 | group` (lme4-style bar). `~ 1 | group` is a
  random intercept; `~ 1 + x1 | group` a random intercept and slope.

- link:

  Mean link: one of `"logit"`, `"probit"`, `"cloglog"` or `"neglog"`.
  The dispersion uses a log link.

- nAGQ:

  Number of adaptive Gauss-Hermite quadrature points per random- effect
  dimension. `nAGQ = 1` is the Laplace approximation.

- start:

  Optional starting vector `c(beta, gamma, omega)`. When `NULL`, fast
  link-specific values are used.

- maxit:

  Maximum number of BFGS iterations.

- rel_tol:

  Relative objective tolerance.

- grad_tol:

  Infinity-norm gradient tolerance.

- n_threads:

  Number of OpenMP threads (over clusters). Zero uses all available
  threads.

- inference:

  Logical; compute the Hessian, covariance matrix and standard errors.

- hessian_rel_step:

  Relative step for the finite-difference Hessian.

- inner_maxit:

  Maximum iterations of the per-cluster inner solver.

- inner_tol:

  Convergence tolerance of the inner solver.

- trace:

  Logical; print optimiser progress.

- subset:

  Optional index vector selecting observations.

- na.action:

  Missing-data handler.

- model, x, y:

  Logical; store the model frame, the design matrices, and the response
  in the fitted object.

## Value

An object of class `"simplex_fast_mixed"`.

## Details

The marginal likelihood integrates the cluster random effects out with
AGHQ (`nAGQ` points per dimension; `nAGQ = 1` gives the Laplace
approximation). The per-cluster inner mode-finding, the quadrature and
the analytic score are implemented in C++ (RcppArmadillo, BLAS) and
parallelised over clusters with OpenMP, so the fit scales to large
nested data sets. The random-effect covariance \\\Sigma = D D^\top\\ is
estimated on an unconstrained log-Cholesky scale, guaranteeing a
positive-definite estimate.

This version supports a single grouping factor (two-level nesting),
Gaussian random effects in the mean submodel, and fixed-effect
(variable) dispersion.

## References

Barndorff-Nielsen, O. E. and Jorgensen, B. (1991). Some parametric
models on the simplex. *Journal of Multivariate Analysis*, **39**(1),
106-116.

Pinheiro, J. C. and Bates, D. M. (1995). Approximations to the
log-likelihood function in the nonlinear mixed-effects model. *Journal
of Computational and Graphical Statistics*, **4**(1), 12-35.

## See also

[`fastsimplexreg()`](https://evandeilton.github.io/fastsimplexreg/reference/fastsimplexreg.md),
[`ranef()`](https://rdrr.io/pkg/nlme/man/random.effects.html),
[`VarCorr()`](https://rdrr.io/pkg/nlme/man/VarCorr.html)

## Examples

``` r
set.seed(1)
J <- 60; nj <- 8; n <- J * nj
dat <- data.frame(
  g  = factor(rep(seq_len(J), each = nj)),
  x1 = rnorm(n),
  z1 = rnorm(n)
)
b <- rnorm(J, 0, 0.7)[dat$g]
mu <- simplex_linkinv(0.3 - 0.6 * dat$x1 + b, "logit")
dat$y <- rsimplex(n, mu, exp(-0.4 + 0.3 * dat$z1))
fit <- fastsimplexregmixed(y ~ x1 | z1, random = ~ 1 | g, data = dat,
                           nAGQ = 7, n_threads = 1)
summary(fit)
#> 
#> Call:
#> fastsimplexregmixed(formula = y ~ x1 | z1, data = dat, random = ~1 | 
#>     g, nAGQ = 7, n_threads = 1)
#> 
#> Pearson residuals:
#>      Min       1Q   Median       3Q      Max 
#> -2.44761 -0.63577  0.02067  0.59711  2.73018 
#> 
#> Coefficients (mean model with logit link):
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  0.16879    0.09652   1.749   0.0803 .  
#> x1          -0.59902    0.01612 -37.162   <2e-16 ***
#> 
#> Coefficients (dispersion model with log link):
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept) -0.36203    0.06904  -5.244 1.57e-07 ***
#> z1           0.29533    0.06615   4.465 8.02e-06 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Random effects:
#> Random effects covariance (group: g)
#>             Variance Std.Dev.
#> (Intercept)   0.5437   0.7374
#> 
#> Log-likelihood: 485.8 | AIC: -961.7 | BIC: -940.8 
#> Observations: 480 | Groups: 60 | nAGQ: 7 | Iterations: 17 
#> Convergence: 0 - Converged: relative objective tolerance satisfied. 
VarCorr(fit)
#> Random effects covariance (group: g)
#>             Variance Std.Dev.
#> (Intercept)   0.5437   0.7374
```
