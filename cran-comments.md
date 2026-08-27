## Test environments

* local Linux, R 4.6.1
* win-builder (devel and release)
* R-hub (Windows, macOS, Linux)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.

## Notes

* The package compiles C++ source via 'Rcpp' and 'RcppArmadillo' and uses
  'OpenMP' when available; every OpenMP use is guarded with `#ifdef _OPENMP`.
  R's random-number generator is only invoked from serial code.
* Examples, tests and the vignette run with `n_threads = 1` and modest sample
  sizes so that check times remain short.
* The installed size (~7 Mb, almost all of it `libs/`) is the usual footprint of
  a compiled 'RcppArmadillo' package; R >= 4.5 reports it as INFO rather than a
  NOTE.

## Downstream dependencies

There are currently no downstream dependencies.
