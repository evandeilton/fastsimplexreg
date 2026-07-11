# Package index

## Model fitting

Fit fixed- and mixed-effects simplex regression models.

- [`fastsimplexreg()`](https://evandeilton.github.io/fastsimplexreg/reference/fastsimplexreg.md)
  : Fit a Fast Simplex Regression with Variable Dispersion
- [`fastsimplexregmixed()`](https://evandeilton.github.io/fastsimplexreg/reference/fastsimplexregmixed.md)
  : Fast Simplex Mixed-Effects Regression with Variable Dispersion

## The simplex distribution

Density, random generation and the mean link inverse.

- [`dsimplex()`](https://evandeilton.github.io/fastsimplexreg/reference/dsimplex.md)
  : Simplex Distribution Density
- [`rsimplex()`](https://evandeilton.github.io/fastsimplexreg/reference/rsimplex.md)
  : Simplex Distribution Random Generation
- [`simplex_linkinv()`](https://evandeilton.github.io/fastsimplexreg/reference/simplex_linkinv.md)
  : Inverse of the Simplex Mean Link

## S3 methods: fixed-effects fits

Methods for objects of class simplex_fast.

- [`coef(`*`<simplex_fast>`*`)`](https://evandeilton.github.io/fastsimplexreg/reference/simplex_fast-methods.md)
  [`vcov(`*`<simplex_fast>`*`)`](https://evandeilton.github.io/fastsimplexreg/reference/simplex_fast-methods.md)
  [`logLik(`*`<simplex_fast>`*`)`](https://evandeilton.github.io/fastsimplexreg/reference/simplex_fast-methods.md)
  [`nobs(`*`<simplex_fast>`*`)`](https://evandeilton.github.io/fastsimplexreg/reference/simplex_fast-methods.md)
  [`fitted(`*`<simplex_fast>`*`)`](https://evandeilton.github.io/fastsimplexreg/reference/simplex_fast-methods.md)
  [`residuals(`*`<simplex_fast>`*`)`](https://evandeilton.github.io/fastsimplexreg/reference/simplex_fast-methods.md)
  [`deviance(`*`<simplex_fast>`*`)`](https://evandeilton.github.io/fastsimplexreg/reference/simplex_fast-methods.md)
  [`model.matrix(`*`<simplex_fast>`*`)`](https://evandeilton.github.io/fastsimplexreg/reference/simplex_fast-methods.md)
  [`terms(`*`<simplex_fast>`*`)`](https://evandeilton.github.io/fastsimplexreg/reference/simplex_fast-methods.md)
  [`formula(`*`<simplex_fast>`*`)`](https://evandeilton.github.io/fastsimplexreg/reference/simplex_fast-methods.md)
  [`model.frame(`*`<simplex_fast>`*`)`](https://evandeilton.github.io/fastsimplexreg/reference/simplex_fast-methods.md)
  [`update(`*`<simplex_fast>`*`)`](https://evandeilton.github.io/fastsimplexreg/reference/simplex_fast-methods.md)
  [`confint(`*`<simplex_fast>`*`)`](https://evandeilton.github.io/fastsimplexreg/reference/simplex_fast-methods.md)
  : Extractor Methods for Simplex Regression Fits
- [`summary(`*`<simplex_fast>`*`)`](https://evandeilton.github.io/fastsimplexreg/reference/summary.simplex_fast.md)
  [`print(`*`<summary.simplex_fast>`*`)`](https://evandeilton.github.io/fastsimplexreg/reference/summary.simplex_fast.md)
  : Summarise a Simplex Regression Fit
- [`print(`*`<simplex_fast>`*`)`](https://evandeilton.github.io/fastsimplexreg/reference/print.simplex_fast.md)
  : Print a Simplex Regression Fit
- [`predict(`*`<simplex_fast>`*`)`](https://evandeilton.github.io/fastsimplexreg/reference/predict.simplex_fast.md)
  : Predictions from a Simplex Regression Fit
- [`plot(`*`<simplex_fast>`*`)`](https://evandeilton.github.io/fastsimplexreg/reference/plot.simplex_fast.md)
  : Diagnostic Plots for a Simplex Regression Fit
- [`simulate(`*`<simplex_fast>`*`)`](https://evandeilton.github.io/fastsimplexreg/reference/simulate.simplex_fast.md)
  : Simulate Responses from a Simplex Regression Fit

## S3 methods: mixed-effects fits

Methods for objects of class simplex_fast_mixed.

- [`coef(`*`<simplex_fast_mixed>`*`)`](https://evandeilton.github.io/fastsimplexreg/reference/simplex_fast_mixed-methods.md)
  [`vcov(`*`<simplex_fast_mixed>`*`)`](https://evandeilton.github.io/fastsimplexreg/reference/simplex_fast_mixed-methods.md)
  [`logLik(`*`<simplex_fast_mixed>`*`)`](https://evandeilton.github.io/fastsimplexreg/reference/simplex_fast_mixed-methods.md)
  [`nobs(`*`<simplex_fast_mixed>`*`)`](https://evandeilton.github.io/fastsimplexreg/reference/simplex_fast_mixed-methods.md)
  [`fitted(`*`<simplex_fast_mixed>`*`)`](https://evandeilton.github.io/fastsimplexreg/reference/simplex_fast_mixed-methods.md)
  [`residuals(`*`<simplex_fast_mixed>`*`)`](https://evandeilton.github.io/fastsimplexreg/reference/simplex_fast_mixed-methods.md)
  [`ranef(`*`<simplex_fast_mixed>`*`)`](https://evandeilton.github.io/fastsimplexreg/reference/simplex_fast_mixed-methods.md)
  [`VarCorr(`*`<simplex_fast_mixed>`*`)`](https://evandeilton.github.io/fastsimplexreg/reference/simplex_fast_mixed-methods.md)
  [`print(`*`<VarCorr.simplex_fast_mixed>`*`)`](https://evandeilton.github.io/fastsimplexreg/reference/simplex_fast_mixed-methods.md)
  : Extractor Methods for Simplex Mixed-Model Fits
- [`summary(`*`<simplex_fast_mixed>`*`)`](https://evandeilton.github.io/fastsimplexreg/reference/summary.simplex_fast_mixed.md)
  [`print(`*`<summary.simplex_fast_mixed>`*`)`](https://evandeilton.github.io/fastsimplexreg/reference/summary.simplex_fast_mixed.md)
  : Summarise a Simplex Mixed-Model Fit
- [`print(`*`<simplex_fast_mixed>`*`)`](https://evandeilton.github.io/fastsimplexreg/reference/print.simplex_fast_mixed.md)
  : Print a Simplex Mixed-Model Fit
- [`predict(`*`<simplex_fast_mixed>`*`)`](https://evandeilton.github.io/fastsimplexreg/reference/predict.simplex_fast_mixed.md)
  : Predictions from a Simplex Mixed-Model Fit
- [`plot(`*`<simplex_fast_mixed>`*`)`](https://evandeilton.github.io/fastsimplexreg/reference/plot.simplex_fast_mixed.md)
  : Diagnostic Plots for a Simplex Mixed-Model Fit
- [`ngrps()`](https://evandeilton.github.io/fastsimplexreg/reference/ngrps.md)
  : Number of Groups in a Mixed-Model Fit
- [`reexports`](https://evandeilton.github.io/fastsimplexreg/reference/reexports.md)
  [`ranef`](https://evandeilton.github.io/fastsimplexreg/reference/reexports.md)
  [`VarCorr`](https://evandeilton.github.io/fastsimplexreg/reference/reexports.md)
  : Objects exported from other packages

## Package

- [`fastsimplexreg-package`](https://evandeilton.github.io/fastsimplexreg/reference/fastsimplexreg-package.md)
  : fastsimplexreg: Fast Simplex Regression with Variable Dispersion
