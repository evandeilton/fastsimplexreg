# Simulate Responses from a Simplex Regression Fit

Simulates new response vectors from a fitted `"simplex_fast"` model by
drawing from the simplex distribution at the fitted means \\\hat\mu_i\\
and dispersions \\\hat\phi_i\\, using
[`rsimplex()`](https://evandeilton.github.io/fastsimplexreg/reference/rsimplex.md).

## Usage

``` r
# S3 method for class 'simplex_fast'
simulate(object, nsim = 1, seed = NULL, ...)
```

## Arguments

- object:

  An object of class `"simplex_fast"`.

- nsim:

  Number of response vectors to simulate.

- seed:

  Optional seed for the random number generator. Handled following the
  convention of
  [`stats::simulate()`](https://rdrr.io/r/stats/simulate.html): when
  supplied, the current RNG state is restored on exit and the seed is
  recorded in the `"seed"` attribute of the result.

- ...:

  Additional arguments, currently ignored.

## Value

A `data.frame` with `nsim` columns (`sim_1`, `sim_2`, ...), each a
simulated response of length `nobs(object)`.

## See also

[`rsimplex()`](https://evandeilton.github.io/fastsimplexreg/reference/rsimplex.md),
[`fastsimplexreg()`](https://evandeilton.github.io/fastsimplexreg/reference/fastsimplexreg.md)

## Examples

``` r
set.seed(5)
dat <- data.frame(x1 = rnorm(200))
dat$y <- rsimplex(200, simplex_linkinv(0.3 + 0.5 * dat$x1, "logit"), 1)
fit <- fastsimplexreg(y ~ x1, data = dat, n_threads = 1L)
sims <- simulate(fit, nsim = 3, seed = 42)
str(sims)
#> 'data.frame':    200 obs. of  3 variables:
#>  $ sim_1: num  0.667 0.894 0.503 0.599 0.788 ...
#>  $ sim_2: num  0.36 0.764 0.365 0.593 0.837 ...
#>  $ sim_3: num  0.287 0.768 0.569 0.66 0.713 ...
#>  - attr(*, "seed")= num 42
#>   ..- attr(*, "kind")=List of 3
#>   .. ..$ : chr "Mersenne-Twister"
#>   .. ..$ : chr "Inversion"
#>   .. ..$ : chr "Rejection"
```
