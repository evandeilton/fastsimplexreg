# The Simplex Distribution

Density, distribution function, quantile function and random generation
for the simplex distribution of Barndorff-Nielsen and Jorgensen (1991),
with mean `mu` and dispersion `phi` (the parameter often written
\\\sigma^2\\). The density is \$\$f(x; \mu, \phi) =
\[2\pi\phi\\(x(1-x))^3\]^{-1/2} \exp\\\left\\-\frac{1}{2\phi}\\
\frac{(x-\mu)^2}{x(1-x)\\\mu^2(1-\mu)^2}\right\\, \qquad 0 \< x \<
1.\$\$

## Usage

``` r
dsimplex(x, mu, phi, log = FALSE, n_threads = 1L)

psimplex(q, mu, phi, lower.tail = TRUE, log.p = FALSE, n_threads = 1L)

qsimplex(p, mu, phi, lower.tail = TRUE, log.p = FALSE, n_threads = 1L)

rsimplex(n, mu, phi)
```

## Arguments

- x, q:

  Numeric vector of quantiles.

- mu:

  Numeric vector of means in \\(0, 1)\\.

- phi:

  Numeric vector of positive dispersion values.

- log, log.p:

  Logical; if `TRUE`, probabilities/densities are given as \\\log(p)\\.

- n_threads:

  Integer number of OpenMP threads. Use `0` to request all threads
  available to the backend. Defaults to `1L` (serial).

- lower.tail:

  Logical; if `TRUE` (default), probabilities are \\P(X \le x)\\,
  otherwise \\P(X \> x)\\.

- p:

  Numeric vector of probabilities.

- n:

  Number of observations to generate. If `length(n) > 1`, the length is
  taken to be the number required (the base-R convention).

## Value

`dsimplex()` gives the density, `psimplex()` the distribution function,
`qsimplex()` the quantile function, and `rsimplex()` generates random
deviates. The length of the result of `rsimplex()` is `n`; for the other
functions it is the maximum of the lengths of the numeric arguments.

## Details

These functions follow the conventions of base R's distribution family:
`x`/`q`/`p`, `mu` and `phi` are recycled to their common length; `NA`
propagates as `NA` and `NaN` as `NaN`; and a parameter outside its
domain (`mu` outside \\(0,1)\\, `phi` not positive, or either
non-finite) produces `NaN` with a warning, rather than an error or a
silent zero. Values of `x` outside the open support \\(0, 1)\\ have
density `0` (`-Inf` on the log scale), which is a genuine density value
and is therefore not a warning.

The distribution function is available in **closed form**. Mapping to
the odds scale \\X = Y/(1-Y)\\ turns the simplex density into a mixture
of an inverse Gaussian and its size-biased version, both of which
integrate exactly, giving \$\$F(y; \mu, \phi) = \Phi(a) + (1 -
2\mu)\\e^{k}\\\Phi(b),\$\$ with \$\$a = \frac{y -
\mu}{\mu(1-\mu)\sqrt{\phi\\ y (1-y)}}, \qquad b = \frac{-(y + \mu - 2 y
\mu)}{\mu(1-\mu)\sqrt{\phi\\ y (1-y)}}, \qquad k = \frac{2}{\phi\\ \mu
(1-\mu)}.\$\$ The product \\e^{k}\Phi(b)\\ is formed on the log scale
and never evaluated directly, so it neither overflows nor loses `log.p`:
`psimplex(0.15, mu = 0.5, phi = 0.01, log.p = TRUE)` returns about
`-773.2` rather than `-Inf`. The upper tail uses the exact reflection
\\P(Y \> y \mid \mu) = F(1 - y \mid 1 - \mu)\\, which follows from the
unit deviance satisfying \\d(y; \mu) = d(1-y; 1-\mu)\\ and keeps full
relative accuracy in both tails.

`qsimplex()` inverts `psimplex()` by safeguarded Newton-bisection. Both
are more expensive than `dsimplex()` and accept `n_threads` for that
reason.

`rsimplex()` uses the exact inverse-Gaussian-mixture representation:
with \\\epsilon = \mu/(1-\mu)\\ and \\\tau = \phi (1-\mu)^2\\, a variate
\\x\\ is built from an inverse-Gaussian draw plus, with probability
\\\mu\\, a chi-squared(1) term, and mapped back to \\(0,1)\\ through
\\x/(1+x)\\.

## References

Barndorff-Nielsen, O. E. and Jorgensen, B. (1991). Some parametric
models on the simplex. *Journal of Multivariate Analysis*, **39**(1),
106–116.
[doi:10.1016/0047-259X(91)90008-P](https://doi.org/10.1016/0047-259X%2891%2990008-P)

## See also

[`fastsimplexreg()`](https://evandeilton.github.io/fastsimplexreg/reference/fastsimplexreg.md)

## Examples

``` r
dsimplex(c(0.2, 0.5, 0.8), mu = 0.5, phi = 1)
#> [1] 0.06924763 3.19153824 0.06924763
dsimplex(c(0.2, 0.5, 0.8), mu = 0.5, phi = 1, log = TRUE)
#> [1] -2.670066  1.160503 -2.670066

# Integrates to one over the support.
psimplex(1, mu = 0.4, phi = 2)
#> [1] 1

# q is the inverse of p.
psimplex(qsimplex(c(0.1, 0.5, 0.9), mu = 0.4, phi = 2), mu = 0.4, phi = 2)
#> [1] 0.1 0.5 0.9

set.seed(123)
y <- rsimplex(1000, mu = 0.35, phi = 0.8)
summary(y)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  0.1308  0.2800  0.3436  0.3469  0.4030  0.6443 
```
