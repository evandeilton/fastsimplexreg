# Number of Groups in a Mixed-Model Fit

Generic and method returning the number of groups (clusters) of the
single grouping factor in a fitted `"simplex_fast_mixed"` model.

## Usage

``` r
ngrps(object, ...)

# S3 method for class 'simplex_fast_mixed'
ngrps(object, ...)
```

## Arguments

- object:

  A fitted model object.

- ...:

  Additional arguments, currently ignored.

## Value

An integer, the number of groups.

## Examples

``` r
set.seed(1)
J <- 40; nj <- 8; n <- J * nj
dat <- data.frame(g = factor(rep(seq_len(J), each = nj)), x1 = rnorm(n))
b <- rnorm(J, 0, 0.7)[dat$g]
dat$y <- rsimplex(n, simplex_linkinv(0.3 - 0.6 * dat$x1 + b, "logit"), 1)
fit <- fastsimplexregmixed(y ~ x1, random = ~ 1 | g, data = dat,
                           nAGQ = 7, n_threads = 1)
ngrps(fit)
#> [1] 40
```
