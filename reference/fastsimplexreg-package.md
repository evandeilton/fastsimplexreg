# fastsimplexreg: Fast Simplex Regression with Variable Dispersion

High-performance maximum-likelihood estimation of simplex regression
models for continuous proportions in the open interval \\(0, 1)\\. The
package supports separate submodels for the mean and the dispersion
through a multi-part
[Formula::Formula](https://rdrr.io/pkg/Formula/man/Formula.html)
interface (`y ~ x1 + x2 | z1 + z2`), four mean links (`logit`, `probit`,
`cloglog`, `neglog`) and a log link for the dispersion. The entire
numerical hot path – log-likelihood, analytic score, native BFGS
optimiser, density, random generation, prediction and link inverses – is
implemented in C++ with RcppArmadillo, BLAS/LAPACK and optional OpenMP
parallelism, so that models scale to large data sets.

## References

Barndorff-Nielsen, O. E. and Jorgensen, B. (1991). Some parametric
models on the simplex. *Journal of Multivariate Analysis*, **39**(1),
106–116.
[doi:10.1016/0047-259X(91)90008-P](https://doi.org/10.1016/0047-259X%2891%2990008-P)

Jorgensen, B. (1997). *The Theory of Dispersion Models*. Chapman & Hall,
London.

Zhang, P., Qiu, Z. and Shi, C. (2016). simplexreg: An R Package for
Regression Analysis of Proportional Data Using the Simplex Distribution.
*Journal of Statistical Software*, **71**(11), 1–21.
[doi:10.18637/jss.v071.i11](https://doi.org/10.18637/jss.v071.i11)

## See also

[`fastsimplexreg()`](https://evandeilton.github.io/fastsimplexreg/reference/fastsimplexreg.md),
[`dsimplex()`](https://evandeilton.github.io/fastsimplexreg/reference/dsimplex.md),
[`rsimplex()`](https://evandeilton.github.io/fastsimplexreg/reference/rsimplex.md),
[`simplex_linkinv()`](https://evandeilton.github.io/fastsimplexreg/reference/simplex_linkinv.md)

## Author

**Maintainer**: José Evandeilton Lopes <evandeilton@gmail.com>
