# Simplex Distribution Random Generation

Generates random deviates from the simplex distribution. The sampler is
implemented in C++ using an exact transformation based on an
inverse-Gaussian mixture representation. The arguments `mu` and `phi`
are recycled to length `n`.

## Usage

``` r
rsimplex(n, mu, phi)
```

## Arguments

- n:

  Integer number of observations to generate.

- mu:

  Numeric vector of means in \\(0, 1)\\, of length one or `n`.

- phi:

  Numeric vector of positive dispersion values, of length one or `n`.

## Value

A numeric vector of length `n` with values in \\(0, 1)\\.

## References

Barndorff-Nielsen, O. E. and Jorgensen, B. (1991). Some parametric
models on the simplex. *Journal of Multivariate Analysis*, **39**(1),
106–116.

## See also

[`dsimplex()`](https://evandeilton.github.io/fastsimplexreg/reference/dsimplex.md),
[`fastsimplexreg()`](https://evandeilton.github.io/fastsimplexreg/reference/fastsimplexreg.md)

## Examples

``` r
set.seed(123)
y <- rsimplex(1000, mu = 0.35, phi = 0.8)
summary(y)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  0.1308  0.2800  0.3436  0.3469  0.4030  0.6443 
```
