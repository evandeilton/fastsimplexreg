# Simplex mixed models for nested proportion data

## 1. Nested proportion data

Continuous proportions in $`(0, 1)`$ are frequently **nested**: repeated
measurements within a subject, students within schools, plots within
sites. The observations within a cluster are correlated, and treating
them as independent understates uncertainty and can bias both the mean
and the dispersion estimates.
[`fastsimplexregmixed()`](https://evandeilton.github.io/fastsimplexreg/reference/fastsimplexregmixed.md)
extends the simplex regression of
[`fastsimplexreg()`](https://evandeilton.github.io/fastsimplexreg/articles/fastsimplexreg.md)
with a cluster-specific random effect, while keeping the
variable-dispersion structure.

## 2. The two-level simplex mixed model

For observation $`i`$ in cluster $`j`$, conditionally on a cluster
random effect $`b_j \sim N_q(0, \Sigma)`$,
``` math
y_{ij} \mid b_j \sim \mathrm{Simplex}(\mu_{ij}, \phi_{ij}), \qquad
g(\mu_{ij}) = x_{ij}^\top\beta + z_{ij}^\top b_j, \qquad
\log\phi_{ij} = w_{ij}^\top\gamma .
```
The fixed effects $`\beta`$ (mean) and $`\gamma`$ (dispersion) use the
same multi-part `Formula` as the fixed-effects fit; the
$`q`$-dimensional random effect $`b_j`$ enters the mean submodel and
induces the within-cluster correlation. The mean link $`g`$ is one of
`logit`, `probit`, `cloglog`, `neglog`; the dispersion uses a log link
and carries fixed effects only.

Because the random effects are unobserved, the likelihood is the
**marginal** likelihood obtained by integrating them out,
``` math
\ell(\theta) = \sum_j \log \int_{\mathbb{R}^q}
  \exp\!\Big\{\textstyle\sum_i \log f(y_{ij}\mid b)\Big\}\,
  N(b; 0, \Sigma)\, \mathrm{d}b ,
```
which has no closed form.

### Adaptive Gauss-Hermite quadrature

The integral is approximated by **adaptive Gauss-Hermite quadrature**
(AGHQ). For each cluster the integrand is centred at its posterior mode
and rescaled by the curvature there, after which a Gauss-Hermite rule
with `nAGQ` points per dimension is applied. `nAGQ = 1` recovers the
**Laplace approximation**; larger `nAGQ` refines the integral, at the
cost of `nAGQ^q` evaluations per cluster. The random-effect covariance
$`\Sigma = D D^\top`$ is estimated on an unconstrained log-Cholesky
scale, so the estimate is always positive definite. The per-cluster
inner mode-finding, the quadrature and the analytic score run in C++ and
are parallelised over clusters with OpenMP.

### The `random` interface

Random effects and the grouping factor are supplied through `random`,
using the lme4-style bar:

- `random = ~ 1 | subject` — a random intercept;
- `random = ~ 1 + t | subject` — a random intercept and a random slope
  in `t`;
- `random = ~ 0 + t | subject` — a random slope only.

## 3. A worked example: simulated intraocular gas decay

We simulate longitudinal data loosely inspired by intraocular gas decay
after retinal surgery: the response `Gas` is the fraction of gas
remaining in the eye, measured **repeatedly over time** for each of 31
patients (variable `ID`), with `LogT`/`LogT2` the (log) time since
surgery and `Time` the raw time. Successive measurements on the same eye
are correlated — a textbook nested structure.

``` r

set.seed(1)
J <- 31; nj <- 6; n <- J * nj
ID <- rep(seq_len(J), each = nj); Time <- rep(seq_len(nj) * 10, J)
LogT <- log(Time); LogT2 <- LogT^2
b <- rnorm(J, 0, 0.8)[ID]
mu <- simplex_linkinv(2.5 - 1.2 * LogT + 0.1 * LogT2 + b, "logit")
retinal <- data.frame(Gas = rsimplex(n, mu, exp(-0.5)),
                      Time = Time, LogT = LogT, LogT2 = LogT2, ID = ID)
retinal$ID <- factor(retinal$ID)
str(retinal)
#> 'data.frame':    186 obs. of  5 variables:
#>  $ Gas  : num  0.516 0.288 0.336 0.23 0.286 ...
#>  $ Time : num  10 20 30 40 50 60 10 20 30 40 ...
#>  $ LogT : num  2.3 3 3.4 3.69 3.91 ...
#>  $ LogT2: num  5.3 8.97 11.57 13.61 15.3 ...
#>  $ ID   : Factor w/ 31 levels "1","2","3","4",..: 1 1 1 1 1 1 2 2 2 2 ...
```

The gas fraction decays with time; we model the mean decay through
`LogT` and `LogT2`, let the dispersion depend on `Time`, and add a
**random intercept per eye** to capture between-subject heterogeneity:

``` r

fit <- fastsimplexregmixed(
  Gas ~ LogT + LogT2 | Time,
  random = ~ 1 | ID,
  data = retinal,
  link = "logit",
  nAGQ = 15,
  n_threads = 1
)
summary(fit)
#> 
#> Call:
#> fastsimplexregmixed(formula = Gas ~ LogT + LogT2 | Time, data = retinal, 
#>     random = ~1 | ID, link = "logit", nAGQ = 15, n_threads = 1)
#> 
#> Pearson residuals:
#>      Min       1Q   Median       3Q      Max 
#> -1.95864 -0.63377 -0.03089  0.51770  2.69270 
#> 
#> Coefficients (mean model with logit link):
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)   3.6650     0.7319   5.008 5.51e-07 ***
#> LogT         -1.9227     0.4672  -4.115 3.87e-05 ***
#> LogT2         0.2158     0.0734   2.940  0.00328 ** 
#> 
#> Coefficients (dispersion model with log link):
#>              Estimate Std. Error z value Pr(>|z|)   
#> (Intercept) -0.897638   0.288661  -3.110  0.00187 **
#> Time         0.007232   0.007514   0.963  0.33579   
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Random effects:
#> Random effects covariance (group: ID)
#>             Variance Std.Dev.
#> (Intercept)   0.5242    0.724
#> 
#> Log-likelihood: 183.6 | AIC: -355.1 | BIC: -335.8 
#> Observations: 186 | Groups: 31 | nAGQ: 15 | Iterations: 21 
#> Convergence: 0 - Converged: relative objective tolerance satisfied.
```

The estimated between-subject variance and the group-level random
effects:

``` r

VarCorr(fit)
#> Random effects covariance (group: ID)
#>             Variance Std.Dev.
#> (Intercept)   0.5242    0.724
head(ranef(fit))
#>   (Intercept)
#> 1 -0.44088405
#> 2  0.07652446
#> 3 -0.72924923
#> 4  1.02256040
#> 5  0.04428093
#> 6 -0.70221986
```

A non-negligible random-intercept variance confirms that eyes differ
systematically in their baseline gas retention beyond what the time
covariates explain. The number of groups and the marginal fit
statistics:

``` r

c(groups = ngrps(fit), nobs = nobs(fit))
#> groups   nobs 
#>     31    186
c(logLik = as.numeric(logLik(fit)), AIC = AIC(fit), BIC = BIC(fit))
#>    logLik       AIC       BIC 
#>  183.5747 -355.1495 -335.7950
```

### Conditional vs population predictions

By default [`predict()`](https://rdrr.io/r/stats/predict.html) is
**conditional** on the estimated random effects (a subject-specific
curve); `re.form = NA` gives the **population-level** (marginal-mode)
prediction with $`b = 0`$:

``` r

newd <- data.frame(LogT = log(c(5, 20, 60)), LogT2 = log(c(5, 20, 60))^2,
                   Time = c(5, 20, 60), ID = factor(levels(retinal$ID)[1],
                                                    levels = levels(retinal$ID)))
data.frame(
  Time        = c(5, 20, 60),
  conditional = predict(fit, newdata = newd, type = "response"),
  population  = predict(fit, newdata = newd, type = "response", re.form = NA)
)
#>   Time conditional population
#> 1    5   0.6656749  0.7557598
#> 2   20   0.3546044  0.4605877
#> 3   60   0.2630857  0.3568396
```

### Diagnostics

The [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method
mirrors the fixed-effects diagnostics, with residuals computed
**conditional on the empirical-Bayes random effects**:

``` r

plot(fit, which = 1:4)
```

![](fastsimplexregmixed_files/figure-html/diagnostics-1.png)

## 4. Choosing `nAGQ` and performance

`nAGQ` trades accuracy for cost. `nAGQ = 1` (Laplace) is fastest and is
often adequate when clusters are large; `nAGQ = 9`-`15` gives high
accuracy for small clusters. Since the tensor grid has
$`\mathtt{nAGQ}^q`$ nodes, keep the number of random effects small
($`q \le 3`$) and lower `nAGQ` as $`q`$ grows (good defaults: `11` for
$`q = 1`$, `9` for $`q = 2`$, `7` for $`q = 3`$). Because clusters are
conditionally independent, the work is parallelised over clusters — set
`n_threads = 0` to use all available cores.

## 5. Scope and extensions

This version supports a single grouping factor (two-level nesting),
Gaussian random effects in the mean submodel, and fixed-effect
(variable) dispersion. Crossed or three-level random effects, and random
effects in the dispersion submodel, are natural extensions for future
releases.

## References

Barndorff-Nielsen, O. E. and Jørgensen, B. (1991). Some parametric
models on the simplex. *Journal of Multivariate Analysis*, **39**(1),
106-116.

Pinheiro, J. C. and Bates, D. M. (1995). Approximations to the
log-likelihood function in the nonlinear mixed-effects model. *Journal
of Computational and Graphical Statistics*, **4**(1), 12-35.

Song, P. X.-K. and Tan, M. (2000). Marginal models for longitudinal
continuous proportional data. *Biometrics*, **56**(2), 496-502.
