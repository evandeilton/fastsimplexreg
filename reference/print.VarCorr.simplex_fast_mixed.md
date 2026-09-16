# Print a Simplex Mixed-Model Variance-Covariance Matrix

Prints the random-effect covariance matrix returned by
[`VarCorr.simplex_fast_mixed()`](https://evandeilton.github.io/fastsimplexreg/reference/simplex_fast_mixed-methods.md),
together with the standard deviations and, for more than one random
effect, the correlation matrix.

## Usage

``` r
# S3 method for class 'VarCorr.simplex_fast_mixed'
print(x, digits = max(3L, getOption("digits") - 3L), ...)
```

## Arguments

- x:

  An object of class `"VarCorr.simplex_fast_mixed"`, as returned by
  [`VarCorr.simplex_fast_mixed()`](https://evandeilton.github.io/fastsimplexreg/reference/simplex_fast_mixed-methods.md)
  – not a fitted model, but the covariance matrix itself, with its
  `stddev`/`correlation`/`group` attributes.

- digits:

  Number of significant digits to display.

- ...:

  Additional arguments, currently ignored.

## Value

The object `x`, invisibly.

## See also

[`VarCorr.simplex_fast_mixed()`](https://evandeilton.github.io/fastsimplexreg/reference/simplex_fast_mixed-methods.md)
