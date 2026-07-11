# Inverse of the Simplex Mean Link

Computes the mean \\\mu = g^{-1}(\eta)\\ from a linear predictor `eta`,
using the same C++ backend employed when fitting a model with
[`fastsimplexreg()`](https://evandeilton.github.io/fastsimplexreg/reference/fastsimplexreg.md).
Four links are supported:

- `logit`:

  \\g^{-1}(\eta) = 1 / (1 + e^{-\eta})\\.

- `probit`:

  \\g^{-1}(\eta) = \Phi(\eta)\\.

- `cloglog`:

  \\g^{-1}(\eta) = 1 - \exp(-\exp(\eta))\\.

- `neglog`:

  \\g^{-1}(\eta) = \exp(-\exp(-\eta))\\, following the definition in
  Zhang et al. (2016).

## Usage

``` r
simplex_linkinv(eta, link = c("logit", "probit", "cloglog", "neglog"))
```

## Arguments

- eta:

  Numeric vector; the linear predictor for the mean.

- link:

  Character string selecting the mean link. One of `"logit"`,
  `"probit"`, `"cloglog"` or `"neglog"`.

## Value

A numeric vector of means in \\(0, 1)\\, of the same length as `eta`.

## References

Zhang, P., Qiu, Z. and Shi, C. (2016). simplexreg: An R Package for
Regression Analysis of Proportional Data Using the Simplex Distribution.
*Journal of Statistical Software*, **71**(11), 1–21.

## See also

[`fastsimplexreg()`](https://evandeilton.github.io/fastsimplexreg/reference/fastsimplexreg.md)

## Examples

``` r
eta <- seq(-3, 3, length.out = 7)
simplex_linkinv(eta, link = "logit")
#> [1] 0.04742587 0.11920292 0.26894142 0.50000000 0.73105858 0.88079708 0.95257413
simplex_linkinv(eta, link = "probit")
#> [1] 0.001349898 0.022750132 0.158655254 0.500000000 0.841344746 0.977249868
#> [7] 0.998650102
simplex_linkinv(eta, link = "cloglog")
#> [1] 0.04856801 0.12657698 0.30779937 0.63212056 0.93401196 0.99938202 1.00000000
simplex_linkinv(eta, link = "neglog")
#> [1] 1.892179e-09 6.179790e-04 6.598804e-02 3.678794e-01 6.922006e-01
#> [6] 8.734230e-01 9.514320e-01
```
