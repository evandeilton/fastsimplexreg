# Simplex Distribution Density

Evaluates the probability density function of the simplex distribution
of Barndorff-Nielsen and Jorgensen (1991) for a mean `mu` and a
dispersion `phi` (the parameter often written \\\sigma^2\\). The density
is \$\$f(x; \mu, \phi) = \[2\pi\phi\\(x(1-x))^3\]^{-1/2}
\exp\\\left\\-\frac{1}{2\phi}\\
\frac{(x-\mu)^2}{x(1-x)\\\mu^2(1-\mu)^2}\right\\, \qquad 0 \< x \<
1.\$\$ Values of `x` outside the open interval \\(0, 1)\\ return `0` (or
`-Inf` on the log scale). The arguments `mu` and `phi` are recycled
against `x`. The computation is carried out in C++ and may use OpenMP
threads.

## Usage

``` r
dsimplex(x, mu, phi, log = FALSE, n_threads = 1L)
```

## Arguments

- x:

  Numeric vector of observations. Values must lie strictly inside \\(0,
  1)\\ to receive positive density.

- mu:

  Numeric vector of means in \\(0, 1)\\, of length one or `length(x)`.

- phi:

  Numeric vector of positive dispersion values, of length one or
  `length(x)`.

- log:

  Logical; if `TRUE`, log-densities are returned.

- n_threads:

  Integer number of OpenMP threads. Use `0` to request all threads
  available to the backend. Defaults to `1L` (serial).

## Value

A numeric vector of densities (or log-densities when `log = TRUE`) with
the recycled length of `x`, `mu` and `phi`.

## References

Barndorff-Nielsen, O. E. and Jorgensen, B. (1991). Some parametric
models on the simplex. *Journal of Multivariate Analysis*, **39**(1),
106–116.

## See also

[`rsimplex()`](https://evandeilton.github.io/fastsimplexreg/reference/rsimplex.md),
[`fastsimplexreg()`](https://evandeilton.github.io/fastsimplexreg/reference/fastsimplexreg.md)

## Examples

``` r
dsimplex(c(0.2, 0.5, 0.8), mu = 0.5, phi = 1)
#> [1] 0.06924763 3.19153824 0.06924763
dsimplex(c(0.2, 0.5, 0.8), mu = 0.5, phi = 1, log = TRUE)
#> [1] -2.670066  1.160503 -2.670066

# Integrates to one over the support.
integrate(function(u) dsimplex(u, mu = 0.4, phi = 2), 0, 1)$value
#> [1] 1
```
