// simplex_fast.cpp
// High-performance simplex regression with variable dispersion.
// Native C++ backend for the fastsimplexreg package. These functions are
// internal package routines exported to R via // [[Rcpp::export]]; the C++
// standard (C++17) and OpenMP flags are supplied by src/Makevars, not by
// Rcpp plugins.

// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include "simplex_common.h"
#include <cmath>
#include <limits>
#include <algorithm>
#include <string>
#include <utility>
#include <vector>
#ifdef _OPENMP
  #include <omp.h>
#endif

using arma::mat;
using arma::vec;
using arma::uword;
using Rcpp::List;
using Rcpp::Named;

namespace simplex_fast {

// The numeric core (link constants and functions, the per-observation kernel,
// the EvalResult bundle and the native BFGS driver) lives in simplex_common.h,
// shared with the mixed-effects backend. This file keeps only the fixed-effects
// evaluator, the RNG-dependent generator, and the Rcpp-exported wrappers.

// Core evaluator of the negative log-likelihood and (optionally) its analytic
// gradient for the simplex regression model with variable dispersion.
//
// For observation i with mean mu_i and dispersion phi_i the simplex
// log-density is
//   ld_i = -0.5*(log(2*pi) + log(phi_i))
//          -1.5*(log(y_i) + log(1 - y_i))
//          -0.5 * dev_i / phi_i,
// where the unit deviance is
//   dev_i = (y_i - mu_i)^2 / [ y_i (1 - y_i) mu_i^2 (1 - mu_i)^2 ].
// The submodels are g(mu_i) = x_i' beta and log(phi_i) = z_i' gamma.
//
// The gradient of the negative log-likelihood is accumulated per observation.
// For the mean submodel the chain rule gives
//   d log L / d beta_j = (d log L / d mu) * (d mu / d eta) * x_ij,
// and for the dispersion submodel, because log(phi) = z' gamma,
//   d log L / d gamma_j = (d log L / d eta_phi) * z_ij,
//   with d log L / d eta_phi = -1/2 + dev/(2*phi).
//
// The observation loop is optionally parallelized with OpenMP. Only pure C++
// arithmetic runs inside the parallel region (no R API calls); each thread
// accumulates into its own private buffers which are reduced afterwards.
// Precompute log(y) + log(1-y), the part of the simplex log-density that
// depends only on the data. Entries for y outside the open support are left at
// zero; evaluate_impl() rejects those observations before ever reading them.
inline vec make_log_yv(const vec& y) {
  vec out(y.n_elem, arma::fill::zeros);
  for (uword i = 0; i < y.n_elem; ++i) {
    const double yi = y[i];
    if (yi > 0.0 && yi < 1.0) out[i] = std::log(yi) + std::log(1.0 - yi);
  }
  return out;
}

EvalResult evaluate_impl(
    const vec& theta,
    const vec& y,
    const mat& X,
    const mat& Z,
    const vec& off_mu,
    const vec& off_phi,
    const vec& log_yv,
    const int mean_link,
    const int n_threads,
    const bool need_grad = true) {

  const uword n = y.n_elem;
  const uword p = X.n_cols;
  const uword q = Z.n_cols;
  const uword d = p + q;

  // An offset of length 0 means "no offset"; anything else must match n.
  if (theta.n_elem != d || X.n_rows != n || Z.n_rows != n ||
      (off_mu.n_elem != 0 && off_mu.n_elem != n) ||
      (off_phi.n_elem != 0 && off_phi.n_elem != n)) {
    return {std::numeric_limits<double>::infinity(), vec(d, arma::fill::zeros), false};
  }

  const vec beta = theta.head(p);
  const vec gamma = theta.tail(q);

  // The two linear predictors are computed INSIDE the parallel region below,
  // one row at a time, rather than by two BLAS gemv calls here. An offset
  // shifts eta by a known constant; it carries no parameter, so the score is
  // unchanged in form -- d log L / d beta_j is still (d log L / d eta) * x_ij.
  //
  // Why not BLAS: these gemvs are memory-bound, and a threaded BLAS makes them
  // slower, not faster. Measured on a 1e6 x 2 matvec, OpenBLAS took 14.9 ms
  // across 24 threads against 5.4 ms pinned to one -- 2.8x slower FOR being
  // parallelized. Worse, the gemvs ran serially with respect to our own OpenMP
  // region, so they became an Amdahl ceiling: at n = 2e6 they were 24% of an
  // evaluation on one thread but 71% on sixteen, capping total speedup near 2x
  // while the observation loop itself scaled 7.5x. Folding them into the loop
  // removes both problems and leaves the package independent of how the user's
  // BLAS is configured.
  const bool use_off_mu = (off_mu.n_elem == n);
  const bool use_off_phi = (off_phi.n_elem == n);

  int threads = 1;
#ifdef _OPENMP
  threads = (n_threads > 0) ? n_threads : omp_get_max_threads();
  threads = std::max(1, threads);
#else
  (void)n_threads;
#endif

  // Per-thread accumulators to avoid data races; reduced after the region.
  // The totals are blocked and compensated sums (see two_sum_acc and BLOCK
  // below), not plain running sums. A plain sum of n = 1e6 terms reaches 8.7e5
  // and each addition rounds at ~1e-10, so the objective was noisy in theta by
  // ~1e-8..1e-7. At one thread the line search then could not see the ~4e-10
  // decrease it was after: the fit stalled at a gradient of 2.7e-2 and took 608
  // evaluations (19.2 s) against 96 (0.9 s) on eight threads, whose partial
  // sums are 8x shorter. With this change both take 96 evaluations, and the
  // objective is identical at 1 and 8 threads.
  std::vector<double> nll_local(static_cast<std::size_t>(threads), 0.0);
  std::vector<double> nll_comp(static_cast<std::size_t>(threads), 0.0);
  std::vector<int> sat_local(static_cast<std::size_t>(threads), 0);
  std::vector<vec> grad_local, grad_comp, grad_blk;
  if (need_grad) {
    grad_local.reserve(static_cast<std::size_t>(threads));
    grad_comp.reserve(static_cast<std::size_t>(threads));
    grad_blk.reserve(static_cast<std::size_t>(threads));
    for (int t = 0; t < threads; ++t) {
      grad_local.emplace_back(d, arma::fill::zeros);
      grad_comp.emplace_back(d, arma::fill::zeros);
      grad_blk.emplace_back(d, arma::fill::zeros);
    }
  }

  int invalid = 0;

#ifdef _OPENMP
  #pragma omp parallel num_threads(threads) reduction(|:invalid)
#endif
  {
    int tid = 0;
#ifdef _OPENMP
    tid = omp_get_thread_num();
#endif
    double local_nll = 0.0, local_nll_c = 0.0;
    int local_sat = 0;
    double* gs = need_grad ? grad_local[static_cast<std::size_t>(tid)].memptr() : nullptr;
    double* gc = need_grad ? grad_comp[static_cast<std::size_t>(tid)].memptr() : nullptr;
    double* gb = need_grad ? grad_blk[static_cast<std::size_t>(tid)].memptr() : nullptr;

    // Terms go into a short plain block first (<= BLOCK of them, so its partial
    // sums stay small and its rounding is negligible); only the block totals
    // enter the compensated thread totals. A TwoSum on every term and component
    // gave the same accuracy but cost +16% per evaluation at n = 1e6; this costs
    // nothing measurable (37 ms either way).
    constexpr uword BLOCK = 256;
    double blk_nll = 0.0;
    uword in_blk = 0;
    auto flush_block = [&]() {
      two_sum_acc(local_nll, local_nll_c, blk_nll);
      blk_nll = 0.0;
      if (need_grad) {
        for (uword j = 0; j < d; ++j) {
          two_sum_acc(gs[j], gc[j], gb[j]);
          gb[j] = 0.0;
        }
      }
      in_blk = 0;
    };

#ifdef _OPENMP
    #pragma omp for schedule(static)
#endif
    for (uword i = 0; i < n; ++i) {
      const double yi = y[i];
      if (!(yi > 0.0 && yi < 1.0) || !std::isfinite(yi)) {
        invalid = 1;
        continue;
      }

      // eta_i = x_i' beta and eta_phi_i = z_i' gamma, accumulated in the same
      // column order a column-major gemv uses, so the result is bit-identical
      // to the BLAS call this replaces (verified across n, p and link).
      double eta_mu_i = 0.0;
      for (uword j = 0; j < p; ++j) eta_mu_i += X(i, j) * beta[j];
      if (use_off_mu) eta_mu_i += off_mu[i];
      double eta_phi_i = 0.0;
      for (uword j = 0; j < q; ++j) eta_phi_i += Z(i, j) * gamma[j];
      if (use_off_phi) eta_phi_i += off_phi[i];

      double mu = 0.0;
      double dmu_deta = 0.0;
      bool sat = false;
      if (!mean_from_eta(eta_mu_i, mean_link, mu, dmu_deta, &sat)) {
        invalid = 1;
        continue;
      }
      if (sat) ++local_sat;

      // Dispersion link is log, so phi = exp(eta_phi) is strictly positive.
      const double phi = safe_exp(eta_phi_i);
      if (!(phi > 0.0) || !std::isfinite(phi)) {
        invalid = 1;
        continue;
      }

      const double one_y = 1.0 - yi;
      const double one_mu = 1.0 - mu;
      const double qmu = mu * one_mu;
      const double diff = yi - mu;
      const double inv_yvar = 1.0 / (yi * one_y);
      const double qmu2 = qmu * qmu;
      // Unit deviance dev = (y-mu)^2 / [y(1-y) (mu(1-mu))^2].
      const double dev = diff * diff * inv_yvar / qmu2;

      // log(y_i) + log(1-y_i) depends only on the DATA, so it is computed once
      // per call site rather than on every one of the ~76 objective evaluations
      // a fit makes. A profile attributed ~20% of a single-threaded fixed-effects
      // fit to libm's log, and two of the three log calls per observation were
      // this constant. The arithmetic is unchanged term for term, so the result
      // is bit-identical -- only the transcendental calls move.
      const double loglik_i = -0.5 * (LOG_2PI + std::log(phi))
                            -1.5 * log_yv[i]
                            -0.5 * dev / phi;

      if (!std::isfinite(loglik_i)) {
        invalid = 1;
        continue;
      }
      // Accumulate the negative log-likelihood.
      blk_nll -= loglik_i;

      if (need_grad) {
        // Analytic score for the mean submodel via the chain rule.
        // d log L / d mu = (y-mu)(mu^2 - 2*mu*y + y)
        //                    / { phi * y * (1-y) * [mu*(1-mu)]^3 }.
        const double score_mu_raw = diff * (mu * mu - 2.0 * mu * yi + yi)
                                  * inv_yvar / (phi * qmu2 * qmu);
        // Multiply by dmu/deta to obtain the score with respect to eta_mu.
        const double score_mu = score_mu_raw * dmu_deta;

        // Score for the dispersion submodel under the log link:
        // d log L / d eta_phi = -1/2 + dev/(2*phi).
        const double score_phi = -0.5 + 0.5 * dev / phi;

        // Gradient of the NEGATIVE log-likelihood accumulates -score * design.
        for (uword j = 0; j < p; ++j) {
          gb[j] -= X(i, j) * score_mu;
        }
        for (uword j = 0; j < q; ++j) {
          gb[p + j] -= Z(i, j) * score_phi;
        }
      }
      if (++in_blk == BLOCK) flush_block();
    }
    flush_block();

    nll_local[static_cast<std::size_t>(tid)] = local_nll;
    nll_comp[static_cast<std::size_t>(tid)] = local_nll_c;
    sat_local[static_cast<std::size_t>(tid)] = local_sat;
  }

  if (invalid != 0) {
    return {std::numeric_limits<double>::infinity(), vec(d, arma::fill::zeros), false};
  }

  // Cross-thread reduction, compensated as well: thread partials are summed
  // with two_sum_acc and every thread's own error term is added at the end.
  double nll = 0.0, nll_c = 0.0;
  for (std::size_t t = 0; t < nll_local.size(); ++t) {
    two_sum_acc(nll, nll_c, nll_local[t]);
    nll_c += nll_comp[t];
  }
  nll += nll_c;

  int n_saturated = 0;
  for (const int value : sat_local) n_saturated += value;

  vec grad(d, arma::fill::zeros);
  if (need_grad) {
    vec grad_c(d, arma::fill::zeros);
    for (std::size_t t = 0; t < grad_local.size(); ++t) {
      for (uword j = 0; j < d; ++j) {
        two_sum_acc(grad[j], grad_c[j], grad_local[t][j]);
        grad_c[j] += grad_comp[t][j];
      }
    }
    grad += grad_c;
  }

  // Only loglik_i is checked per observation, but the score can still overflow
  // (e.g. 0 * Inf where the mean is clamped and phi is tiny, which gave a NaN
  // gradient with valid = TRUE). A non-finite term never becomes finite by
  // addition, so one check on the totals rejects the point like any other
  // invalid one and the line search simply backs off.
  if (!std::isfinite(nll) || (need_grad && !grad.is_finite())) {
    return {std::numeric_limits<double>::infinity(), vec(d, arma::fill::zeros), false};
  }

  return {nll, std::move(grad), true, n_saturated};
}

// Draw a single inverse-Gaussian variate by the Michael-Schucany-Haas
// algorithm. A chi-squared(1) draw z produces a candidate root x; with the
// acceptance probability mean/(mean + x) the smaller root is kept, otherwise
// the reflected root mean^2/x is returned. This calls R's RNG and therefore
// must run in serial code only.
inline double inv_gaussian_one(const double mean, const double tau) {
  // Parameterization inherited from the simplex-regression mixture generator.
  const double z = R::rchisq(1.0);
  const double root = std::sqrt(4.0 * mean * z / tau + (mean * z) * (mean * z));
  double x = mean + 0.5 * mean * mean * tau * z - 0.5 * mean * tau * root;
  x = std::max(x, std::numeric_limits<double>::min());
  if (R::runif(0.0, 1.0) > mean / (mean + x)) {
    x = mean * mean / x;
  }
  return x;
}

constexpr double HALF_LOG_2PI = 0.5 * LOG_2PI;
constexpr double LOG_2 = 0.693147180559945309417232121458;

// log(1 - exp(x)) for x < 0, switching form at -log(2) so that neither a tiny
// nor an O(1) exp(x) loses digits (Maechler, 2012).
inline double log1mexp(const double x) noexcept {
  return (x > -LOG_2) ? std::log(-std::expm1(x)) : std::log1p(-std::exp(x));
}

// log of the Mills ratio m(x) = Phi(-x) / phi(x) for x >= 0, with m(Inf) -> 1/x.
//
// Below x = 3 the direct form log Phi(-x) + x^2/2 + log sqrt(2 pi) is accurate
// to ~4e-16. Above it the rounding of log Phi(-x) itself (|log Phi| ~ x^2/2)
// grows to ~2e-15 at x = 6 and worse beyond, so the Laplace continued fraction
//   m(x) = 1 / (x + 1/(x + 2/(x + 3/(x + ...)))) = 1 / t_0,
//   t_k = x + (k+1)/t_{k+1},
// is evaluated backward instead. The tail is started at the fixed point of
// t = x + (n+1)/t rather than at x, which saves ~30% of the terms; the depth n
// per band is the smallest tested one plus 2, and reproduces a 20000-term
// evaluation to <= 1.5 ulp on a dense grid over [3, 1e8]. For x > 1e8,
// m(x) = 1/x (1 - 1/x^2 + ...) is 1/x to 1e-16 (and Inf gives -Inf).
inline double log_mills(const double x) noexcept {
  if (x < 3.0) return R::pnorm(-x, 0.0, 1.0, 1, 1) + 0.5 * x * x + HALF_LOG_2PI;
  if (x > 1e8) return -std::log(x);
  const int n = (x < 3.5) ? 38 : (x < 4.0) ? 31 : (x < 5.0) ? 26 : (x < 6.0) ? 20
              : (x < 8.0) ? 17 : (x < 10.0) ? 13 : (x < 15.0) ? 11 : (x < 20.0) ? 9
              : (x < 30.0) ? 8 : (x < 50.0) ? 6 : 5;
  double t = 0.5 * (x + std::sqrt(x * x + 4.0 * (n + 1)));
  for (int k = n; k >= 1; --k) t = x + k / t;
  return -std::log(t);
}

// log of one tail of the simplex distribution, in CLOSED FORM.
//
// Writing X = Y/(1-Y) maps the simplex onto the odds scale, where its density
// is the mixture (1-mu) f_IG(x) + mu * (size-biased f_IG)(x) with an inverse
// Gaussian of mean mu/(1-mu) and shape 1/(phi (1-mu)^2). Both pieces integrate
// in closed form, and after collecting terms
//
//   F(y; mu, phi) = Phi(a) + (1 - 2 mu) exp(k) Phi(b),
//     a = (y - mu) / s,    b = -(y + mu - 2 y mu) / s,
//     s = mu(1-mu) sqrt(phi y (1-y)),    k = 2 / (phi mu (1-mu)).
//
// Forming exp(k) overflows already at k = 800, and k + log Phi(b) subtracts two
// numbers of size k, so for tiny mu (k ~ 1e8) it carries an absolute error of
// ~1e-8 or worse. Neither is necessary. With A = y(1-mu) and B = mu(1-y) one has
// a = (A-B)/s, b = -(A+B)/s and (A+B)^2 - (A-B)^2 = 4AB = 4 mu(1-mu) y(1-y), so
//
//   k - b^2/2 = -a^2/2    EXACTLY,
//
// and with beta = -b > 0 and the Mills ratio m(x) = Phi(-x)/phi(x)
//
//   exp(k) Phi(b) = exp(-a^2/2) m(beta) / sqrt(2 pi) = phi(a) m(beta),
//   F = Phi(a) + (1 - 2 mu) phi(a) m(beta).
//
// Nothing of size k is ever formed. For a < 0 the identity Phi(a) = phi(a) m(-a)
// also makes the ratio of the two terms a pure difference of log Mills ratios,
// d = log|c| + log m(beta) - log m(-a), in which the a^2/2 terms cancel exactly.
//
// The arguments are the quantities that enter, supplied by the caller so that no
// "1 - x" is ever formed from a rounded complement: u = mu(1-mu), v = y(1-y),
// A, B as above, diff = y - mu, c = 1 - 2 mu. The upper tail is the exact
// reflection S(y; mu) = F(1-y; 1-mu) (the unit deviance satisfies
// d(y; mu) = d(1-y; 1-mu)), which merely swaps A and B and flips the signs of
// diff and c -- so a tiny mu is never converted into mu' = 1 - mu ~ 1, where the
// complement 1 - mu' would keep only ~16 - log10(1/mu) digits.
//
// Remaining limit: for mu > 1/2 and a < 0 the result is log Phi(a) + log(1 - r)
// with r = |c| m(beta)/m(-a) < 1, and 1 - r is a genuine cancellation when mu is
// within ~1e-8 of 1 (relative error ~ 1e-16 / (1 - mu)); the value stays finite
// and in [0, 1], and is -Inf only if rounding makes r >= 1.
//
// Measured: F + S = 1 to 2.2e-16 for mu from 1e-12 to 1 - 1e-12 (the previous
// code: 1 at mu = 1e-8), and both tails match a quadrature of the density to
// 2e-15 down to mu = 1e-10.
inline double simplex_logtail_core(const double u, const double v,
                                   const double A, const double B,
                                   const double diff, const double c,
                                   const double phi) noexcept {
  constexpr double NINF = -std::numeric_limits<double>::infinity();
  const double s = u * std::sqrt(phi * v);
  // Degenerate (s underflows or overflows): the distribution is a point mass.
  if (!(s > 0.0) || !std::isfinite(s)) return (diff >= 0.0) ? 0.0 : NINF;
  const double a = diff / s;
  if (!std::isfinite(a)) return (a > 0.0) ? 0.0 : NINF;
  const double t1 = R::pnorm(a, 0.0, 1.0, 1, 1);          // log Phi(a)
  if (c == 0.0) return t1;                                // mu = 1/2: the term vanishes

  const double lc = std::log(std::fabs(c));
  const double lmb = log_mills((A + B) / s);              // log m(beta); -Inf if beta = Inf
  if (a >= 0.0) {
    // Phi(a) >= 1/2 dominates, so work with the log of the second term directly.
    const double lt2 = lc - 0.5 * a * a - HALF_LOG_2PI + lmb;
    if (c > 0.0) {                                        // mu < 1/2: log-sum-exp
      const double m = std::max(t1, lt2);
      return m + std::log(std::exp(t1 - m) + std::exp(lt2 - m));
    }
    // mu > 1/2: log-difference-exp. t1 dominates in exact arithmetic; if rounding
    // says otherwise the true value has underflowed to zero.
    if (lt2 >= t1) return NINF;
    return t1 + log1mexp(lt2 - t1);
  }
  // a < 0: factor Phi(a) out; the a^2/2 terms of the two pieces cancel exactly.
  const double d = lc + lmb - log_mills(-a);
  if (c > 0.0) return t1 + std::log1p(std::exp(d));       // mu < 1/2
  if (d >= 0.0) return NINF;                              // mu > 1/2, rounding
  return t1 + log1mexp(d);
}

// log of ONE tail at a single point: lower = TRUE gives log P(Y <= q) and
// FALSE gives log P(Y > q). Both are evaluated directly, never as log(1 - .),
// so each keeps full relative accuracy in its own far tail.
inline double simplex_logtail_one(const double q, const double mu, const double phi,
                                  const bool lower) noexcept {
  constexpr double NINF = -std::numeric_limits<double>::infinity();
  if (q <= 0.0) return lower ? NINF : 0.0;
  if (q >= 1.0) return lower ? 0.0 : NINF;
  const double qbar = 1.0 - q, mubar = 1.0 - mu;
  const double u = mu * mubar, v = q * qbar;
  const double A = q * mubar, B = mu * qbar;
  const double diff = q - mu, c = 1.0 - 2.0 * mu;
  // Upper tail: reflection (y, mu) -> (1-y, 1-mu), i.e. A <-> B, diff -> -diff, c -> -c.
  const double lt = lower ? simplex_logtail_core(u, v, A, B, diff, c, phi)
                          : simplex_logtail_core(u, v, B, A, -diff, -c, phi);
  return (lt > 0.0) ? 0.0 : lt;                 // a log probability is <= 0; trim rounding
}

// Simplex quantile from the log of a tail probability, by safeguarded
// Newton-bisection (Press et al., "rtsafe"): a Newton step is taken only when
// it stays inside the current bracket and halves the step, otherwise the method
// bisects. Bisection alone would converge; the Newton steps cut the number of
// (expensive) CDF evaluations to roughly a dozen.
//
// The equation solved is  log T(y) - lp = 0,  T the requested tail, so that
// probabilities down to 1e-300 (any lp <= 0, which is how log.p = TRUE arrives)
// stay exact, and the upper tail is solved in y itself -- never as 1 - q(1-mu) --
// so a tiny quantile keeps its relative precision. The derivative of log T is
// +-f/T. A probability above 1/2 is solved on the opposite tail with
// log(1 - p), which is exact for p near 1 (p = 1 - 1e-12 would otherwise be
// pinned to 1e-4 relative accuracy by log F = -1e-12).
//
// The step test is RELATIVE (10 eps * y): an absolute 1e-15 is about one
// standard deviation when mu ~ 1e-10, where sd = sqrt(phi mu^3) ~ 1e-15.
inline double simplex_quantile_one(double lp, const double mu, const double phi,
                                   bool lower) noexcept {
  if (!(lp > -std::numeric_limits<double>::infinity())) return lower ? 0.0 : 1.0;
  if (!(lp < 0.0)) return lower ? 1.0 : 0.0;
  if (lp > -LOG_2) { lp = log1mexp(lp); lower = !lower; }

  constexpr double TOL = 10.0 * std::numeric_limits<double>::epsilon();
  double lo = 0.0, hi = 1.0;
  double y = mu;
  double dy_old = 1.0, dy = 1.0;

  for (int it = 0; it < 200; ++it) {
    const double lt = simplex_logtail_one(y, mu, phi, lower);
    const double g = lt - lp;                    // zero at the root
    if (std::isnan(g)) return std::numeric_limits<double>::quiet_NaN();
    if (g == 0.0) return y;
    // g > 0: the lower tail is too large (y too big) or the upper one is (y too small).
    if ((g > 0.0) == lower) hi = y; else lo = y;

    // d g / d y = +-f/T, the hazard-type ratio of density to tail probability.
    const double r = std::exp(simplex_logpdf(y, mu, phi) - lt);
    const double dg = lower ? r : -r;

    const bool newton_ok =
        (r > 0.0) && std::isfinite(r) &&
        (((y - hi) * dg - g) * ((y - lo) * dg - g) <= 0.0) &&
        (std::abs(2.0 * g) <= std::abs(dy_old * dg));

    dy_old = dy;
    if (newton_ok) {
      dy = g / dg;
      y -= dy;
    } else {
      dy = 0.5 * (hi - lo);
      y = lo + dy;
    }
    if (std::abs(dy) <= TOL * y || hi - lo <= TOL * y) break;
  }
  return y;
}

} // namespace simplex_fast


// Fast simplex density in C++.
// Follows the base-R d*() contract exactly:
//   * NA / NaN in any argument propagates (NA wins over NaN, via the standard
//     `x + mu + phi` idiom used throughout R's own d*() sources);
//   * an out-of-domain mu or phi yields NaN and is counted, so the R wrapper can
//     raise the canonical "NaNs produced" warning;
//   * x outside the open support (0, 1) yields 0 (or -Inf on the log scale),
//     which is a genuine density value, not a missing one;
//   * x, mu and phi are recycled to their common maximum length.
// The per-observation loop is optionally parallelized with OpenMP; it touches no
// R API state and is therefore thread-safe.
// [[Rcpp::export]]
Rcpp::NumericVector dsimplex_cpp(
    const Rcpp::NumericVector& y,
    const Rcpp::NumericVector& mu,
    const Rcpp::NumericVector& phi,
    const bool log = false,
    const int n_threads = 1) {

  const R_xlen_t ny = y.size(), nm = mu.size(), np = phi.size();
  if (ny == 0 || nm == 0 || np == 0) return Rcpp::NumericVector(0);
  const R_xlen_t n = std::max(ny, std::max(nm, np));

  Rcpp::NumericVector out(n);
  int threads = 1;
#ifdef _OPENMP
  threads = (n_threads > 0) ? n_threads : omp_get_max_threads();
  threads = std::max(1, threads);
#else
  (void)n_threads;
#endif

  int n_bad = 0;
  // Chunked so that an interrupt check can sit BETWEEN parallel regions --
  // Rcpp::checkUserInterrupt() throws, so it must never be called from inside
  // one. 65536 elements is ~0.10 s of work for the slowest of these kernels
  // (qsimplex, measured at 1.58 us/element single-threaded), which bounds the
  // interrupt latency, while libgomp's region-entry cost of a few microseconds
  // is under 0.01% of a chunk. n_bad becomes a sum of partial reductions, which
  // is an integer count, so nothing numerical changes.
  constexpr R_xlen_t CHUNK = 1 << 16;
  for (R_xlen_t base = 0; base < n; base += CHUNK) {
    Rcpp::checkUserInterrupt();
    const R_xlen_t hi = std::min(n, base + CHUNK);
#ifdef _OPENMP
  #pragma omp parallel for num_threads(threads) schedule(static) reduction(+:n_bad)
#endif
    for (R_xlen_t i = base; i < hi; ++i) {
    const double yi = y[i % ny];
    const double mui = mu[i % nm];
    const double phii = phi[i % np];

    if (ISNAN(yi) || ISNAN(mui) || ISNAN(phii)) {
      out[i] = yi + mui + phii;        // propagates NA, else NaN
      continue;
    }
    if (!(mui > 0.0 && mui < 1.0) || !(phii > 0.0) ||
        !R_FINITE(mui) || !R_FINITE(phii)) {
      out[i] = R_NaN;                  // invalid parameter, R convention
      ++n_bad;
      continue;
    }
    if (!(yi > 0.0 && yi < 1.0)) {
      out[i] = log ? R_NegInf : 0.0;   // outside the support: density is zero
      continue;
    }
    const double ld = simplex_fast::simplex_logpdf(yi, mui, phii);
    out[i] = log ? ld : std::exp(ld);
    }
  }

  out.attr("n_invalid_par") = n_bad;
  return out;
}


// Simplex distribution function in C++.
// Evaluated from the closed form in simplex_logtail_core(); log_p is computed on
// the log scale throughout rather than as log() of a linear value. Argument
// conventions follow base R's p*(): `lower_tail` and `log_p` are honored,
// NA/NaN propagate, and invalid parameters give NaN.
// [[Rcpp::export]]
Rcpp::NumericVector psimplex_cpp(
    const Rcpp::NumericVector& q,
    const Rcpp::NumericVector& mu,
    const Rcpp::NumericVector& phi,
    const bool lower_tail = true,
    const bool log_p = false,
    const int n_threads = 1) {

  const R_xlen_t nq = q.size(), nm = mu.size(), np = phi.size();
  if (nq == 0 || nm == 0 || np == 0) return Rcpp::NumericVector(0);
  const R_xlen_t n = std::max(nq, std::max(nm, np));

  Rcpp::NumericVector out(n);
  int threads = 1;
#ifdef _OPENMP
  threads = (n_threads > 0) ? n_threads : omp_get_max_threads();
  threads = std::max(1, threads);
#else
  (void)n_threads;
#endif

  int n_bad = 0;
  // Chunked so that an interrupt check can sit BETWEEN parallel regions --
  // Rcpp::checkUserInterrupt() throws, so it must never be called from inside
  // one. 65536 elements is ~0.10 s of work for the slowest of these kernels
  // (qsimplex, measured at 1.58 us/element single-threaded), which bounds the
  // interrupt latency, while libgomp's region-entry cost of a few microseconds
  // is under 0.01% of a chunk. n_bad becomes a sum of partial reductions, which
  // is an integer count, so nothing numerical changes.
  constexpr R_xlen_t CHUNK = 1 << 16;
  for (R_xlen_t base = 0; base < n; base += CHUNK) {
    Rcpp::checkUserInterrupt();
    const R_xlen_t hi = std::min(n, base + CHUNK);
#ifdef _OPENMP
  #pragma omp parallel for num_threads(threads) schedule(static) reduction(+:n_bad)
#endif
    for (R_xlen_t i = base; i < hi; ++i) {
    const double qi = q[i % nq];
    const double mui = mu[i % nm];
    const double phii = phi[i % np];

    if (ISNAN(qi) || ISNAN(mui) || ISNAN(phii)) {
      out[i] = qi + mui + phii;
      continue;
    }
    if (!(mui > 0.0 && mui < 1.0) || !(phii > 0.0) ||
        !R_FINITE(mui) || !R_FINITE(phii)) {
      out[i] = R_NaN;
      ++n_bad;
      continue;
    }
    // The log scale is computed DIRECTLY, never as log() of a linear value
    // that has already underflowed: psimplex(0.15, 0.5, 0.01, log.p = TRUE)
    // used to return -Inf where the true value is -773.216.
    const double lv = simplex_fast::simplex_logtail_one(qi, mui, phii, lower_tail);
    out[i] = log_p ? lv : std::min(1.0, std::max(0.0, std::exp(lv)));
    }
  }

  out.attr("n_invalid_par") = n_bad;
  return out;
}


// Simplex quantile function in C++.
// Inverts psimplex_cpp() by safeguarded Newton-bisection. Argument conventions
// follow base R's q*(): `lower_tail` and `log_p` are honored, NA/NaN
// propagate, and a probability outside [0, 1] (or an invalid parameter) gives
// NaN.
// [[Rcpp::export]]
Rcpp::NumericVector qsimplex_cpp(
    const Rcpp::NumericVector& p,
    const Rcpp::NumericVector& mu,
    const Rcpp::NumericVector& phi,
    const bool lower_tail = true,
    const bool log_p = false,
    const int n_threads = 1) {

  const R_xlen_t npr = p.size(), nm = mu.size(), np = phi.size();
  if (npr == 0 || nm == 0 || np == 0) return Rcpp::NumericVector(0);
  const R_xlen_t n = std::max(npr, std::max(nm, np));

  Rcpp::NumericVector out(n);
  int threads = 1;
#ifdef _OPENMP
  threads = (n_threads > 0) ? n_threads : omp_get_max_threads();
  threads = std::max(1, threads);
#else
  (void)n_threads;
#endif

  int n_bad = 0;
  // Chunked so that an interrupt check can sit BETWEEN parallel regions --
  // Rcpp::checkUserInterrupt() throws, so it must never be called from inside
  // one. 65536 elements is ~0.10 s of work for the slowest of these kernels
  // (qsimplex, measured at 1.58 us/element single-threaded), which bounds the
  // interrupt latency, while libgomp's region-entry cost of a few microseconds
  // is under 0.01% of a chunk. n_bad becomes a sum of partial reductions, which
  // is an integer count, so nothing numerical changes.
  constexpr R_xlen_t CHUNK = 1 << 16;
  for (R_xlen_t base = 0; base < n; base += CHUNK) {
    Rcpp::checkUserInterrupt();
    const R_xlen_t hi = std::min(n, base + CHUNK);
#ifdef _OPENMP
  #pragma omp parallel for num_threads(threads) schedule(static) reduction(+:n_bad)
#endif
    for (R_xlen_t i = base; i < hi; ++i) {
    const double pi_ = p[i % npr];
    const double mui = mu[i % nm];
    const double phii = phi[i % np];

    if (ISNAN(pi_) || ISNAN(mui) || ISNAN(phii)) {
      out[i] = pi_ + mui + phii;
      continue;
    }
    if (!(mui > 0.0 && mui < 1.0) || !(phii > 0.0) ||
        !R_FINITE(mui) || !R_FINITE(phii)) {
      out[i] = R_NaN;
      ++n_bad;
      continue;
    }
    // The solver works on log p of the REQUESTED tail, so neither 1 - p nor
    // exp(log p) is ever formed: p = 1e-300 and any log p <= 0 are exact.
    if (log_p ? (pi_ > 0.0) : (pi_ < 0.0 || pi_ > 1.0)) {
      out[i] = R_NaN;
      ++n_bad;
      continue;
    }
    const double lp = log_p ? pi_ : std::log(pi_);

    out[i] = simplex_fast::simplex_quantile_one(lp, mui, phii, lower_tail);
    }
  }

  out.attr("n_invalid_par") = n_bad;
  return out;
}


// Fast random generation from the simplex distribution in C++.
// Uses the exact inverse-Gaussian-mixture transformation: with
// epsilon = mu/(1-mu) and tau = phi (1-mu)^2, a variate x is built from an
// inverse-Gaussian draw plus, with probability mu, a chi-squared(1) term; the
// result is mapped back to (0,1) via x/(1+x). Because it calls R's RNG, the
// loop is kept strictly serial (never parallelize R API calls).
// Follows the base-R r*() contract: mu and phi are recycled to length n, and an
// invalid or missing parameter yields NaN for that draw (counted, so the R
// wrapper can raise the canonical "NAs produced" warning) rather than aborting
// the whole call.
// [[Rcpp::export]]
Rcpp::NumericVector rsimplex_cpp(
    const R_xlen_t n,
    const Rcpp::NumericVector& mu,
    const Rcpp::NumericVector& phi) {

  if (n < 0) Rcpp::stop("'n' must be non-negative.");
  const R_xlen_t nm = mu.size(), np = phi.size();
  if (n > 0 && (nm == 0 || np == 0)) {
    Rcpp::stop("'mu' and 'phi' must have positive length.");
  }

  Rcpp::RNGScope scope;
  Rcpp::NumericVector out(n);
  int n_bad = 0;

  for (R_xlen_t i = 0; i < n; ++i) {
    const double mui = mu[i % nm];
    const double phii = phi[i % np];
    if (ISNAN(mui) || ISNAN(phii) ||
        !(mui > 0.0 && mui < 1.0) || !(phii > 0.0) ||
        !R_FINITE(mui) || !R_FINITE(phii)) {
      out[i] = R_NaN;
      ++n_bad;
      continue;
    }

    const double epsilon = mui / (1.0 - mui);
    const double tau = phii * (1.0 - mui) * (1.0 - mui);

    const double x1 = simplex_fast::inv_gaussian_one(epsilon, tau);
    const double x3 = R::rchisq(1.0) * tau * epsilon * epsilon;
    const double x = (R::runif(0.0, 1.0) < mui) ? (x1 + x3) : x1;
    out[i] = x / (1.0 + x);
  }

  out.attr("n_invalid_par") = n_bad;
  return out;
}


// Evaluate the negative log-likelihood and its analytic gradient at a given
// parameter vector theta = c(beta, gamma). Returns a list with elements
// "value" (the NLL), "gradient" (its gradient), and "valid" (false when any
// observation falls outside the model support).
// [[Rcpp::export]]
Rcpp::List simplex_eval_cpp(
    const arma::vec& theta,
    const arma::vec& y,
    const arma::mat& X,
    const arma::mat& Z,
    const int mean_link = 1,
    const int n_threads = 1,
    Rcpp::Nullable<Rcpp::NumericVector> off_mu_ = R_NilValue,
    Rcpp::Nullable<Rcpp::NumericVector> off_phi_ = R_NilValue) {

  const arma::vec off_mu = off_mu_.isNotNull()
    ? Rcpp::as<arma::vec>(off_mu_.get()) : arma::vec();
  const arma::vec off_phi = off_phi_.isNotNull()
    ? Rcpp::as<arma::vec>(off_phi_.get()) : arma::vec();

  const arma::vec log_yv = simplex_fast::make_log_yv(y);
  const auto res = simplex_fast::evaluate_impl(theta, y, X, Z, off_mu, off_phi,
                                               log_yv, mean_link, n_threads, true);
  return List::create(
    Named("value") = res.nll,
    Named("gradient") = res.grad,
    Named("valid") = res.valid,
    Named("n_saturated") = res.n_saturated
  );
}


// Fit simplex regression by a native BFGS optimizer with an Armijo
// backtracking line search.
//
// At each iteration the search direction is -H * grad, where H is the current
// BFGS approximation to the inverse Hessian. If that direction is not a descent
// direction the method resets H to the identity and falls back to steepest
// descent. The Armijo line search shrinks the step by a factor of 0.5 until the
// sufficient-decrease condition nll(theta + step*dir) <= nll + c1*step*slope
// holds. The inverse Hessian is updated by the BFGS formula
//   H <- (I - rho s y') H (I - rho y s') + rho s s',  rho = 1/(y's),
// with the curvature safeguard y's > 0, and symmetrized each step. Convergence
// is declared on the infinity-norm of the gradient or on a small relative change
// in the objective combined with a modest gradient norm.
// [[Rcpp::export]]
Rcpp::List simplex_bfgs_cpp(
    const arma::vec& start,
    const arma::vec& y,
    const arma::mat& X,
    const arma::mat& Z,
    const int mean_link = 1,
    const int maxit = 300,
    const double rel_tol = 1e-9,
    const double grad_tol = 1e-6,
    const int n_threads = 1,
    const bool trace = false,
    Rcpp::Nullable<Rcpp::NumericVector> off_mu_ = R_NilValue,
    Rcpp::Nullable<Rcpp::NumericVector> off_phi_ = R_NilValue) {

  const arma::vec off_mu = off_mu_.isNotNull()
    ? Rcpp::as<arma::vec>(off_mu_.get()) : arma::vec();
  const arma::vec off_phi = off_phi_.isNotNull()
    ? Rcpp::as<arma::vec>(off_phi_.get()) : arma::vec();

  // Delegate to the shared native BFGS driver, wrapping the fixed-effects
  // evaluator as the objective. The optimizer logic is identical to before;
  // it now lives once in simplex_common.h and is reused by the mixed backend.
  const arma::vec log_yv = simplex_fast::make_log_yv(y);
  auto objective = [&](const arma::vec& th) {
    return simplex_fast::evaluate_impl(th, y, X, Z, off_mu, off_phi, log_yv,
                                       mean_link, n_threads, true);
  };
  return simplex_fast::bfgs_minimize(start, objective, maxit, rel_tol, grad_tol, trace);
}


// Finite-difference Hessian of the negative log-likelihood.
// Forms the Hessian by central differencing of the analytic gradient: column j
// is (grad(theta + h e_j) - grad(theta - h e_j)) / (2h) with a relative step h.
// If either perturbed evaluation is non-finite the step is shrunk adaptively.
// The result is symmetrized. Intended for post-fit inference only.
// [[Rcpp::export]]
arma::mat simplex_hessian_fd_cpp(
    const arma::vec& theta,
    const arma::vec& y,
    const arma::mat& X,
    const arma::mat& Z,
    const int mean_link = 1,
    const double rel_step = 1e-5,
    const int n_threads = 1,
    Rcpp::Nullable<Rcpp::NumericVector> off_mu_ = R_NilValue,
    Rcpp::Nullable<Rcpp::NumericVector> off_phi_ = R_NilValue) {

  const arma::vec off_mu = off_mu_.isNotNull()
    ? Rcpp::as<arma::vec>(off_mu_.get()) : arma::vec();
  const arma::vec off_phi = off_phi_.isNotNull()
    ? Rcpp::as<arma::vec>(off_phi_.get()) : arma::vec();

  const uword d = theta.n_elem;
  mat H(d, d, arma::fill::zeros);
  const arma::vec log_yv = simplex_fast::make_log_yv(y);

  for (uword j = 0; j < d; ++j) {
    // Serial loop of 2d objective evaluations; a meaningful share of a long
    // fit, and previously unreachable by Ctrl-C.
    Rcpp::checkUserInterrupt();
    double h = rel_step * std::max(1.0, std::abs(theta[j]));
    bool success = false;

    for (int attempt = 0; attempt < 12; ++attempt) {
      vec plus = theta;
      vec minus = theta;
      plus[j] += h;
      minus[j] -= h;

      const auto gp = simplex_fast::evaluate_impl(plus, y, X, Z, off_mu, off_phi,
                                                 log_yv, mean_link, n_threads, true);
      const auto gm = simplex_fast::evaluate_impl(minus, y, X, Z, off_mu, off_phi,
                                                 log_yv, mean_link, n_threads, true);

      if (gp.valid && gm.valid && std::isfinite(gp.nll) && std::isfinite(gm.nll)) {
        H.col(j) = (gp.grad - gm.grad) / (2.0 * h);
        success = true;
        break;
      }
      h *= 0.25;
    }

    if (!success) {
      Rcpp::stop("Non-finite evaluation while computing the Hessian after adaptive step reduction.");
    }
  }

  return 0.5 * (H + H.t());
}


// Compute fitted mean and dispersion vectors from a parameter vector.
// Applies the mean link inverse to eta_mu = X * beta and the log link inverse
// (exp) to eta_phi = Z * gamma. Returns a list with "mu", "phi", and both
// linear predictors "eta_mu" and "eta_phi".
// [[Rcpp::export]]
Rcpp::List simplex_predict_cpp(
    const arma::vec& theta,
    const arma::mat& X,
    const arma::mat& Z,
    const int mean_link = 1,
    Rcpp::Nullable<Rcpp::NumericVector> off_mu_ = R_NilValue,
    Rcpp::Nullable<Rcpp::NumericVector> off_phi_ = R_NilValue) {

  const arma::vec off_mu = off_mu_.isNotNull()
    ? Rcpp::as<arma::vec>(off_mu_.get()) : arma::vec();
  const arma::vec off_phi = off_phi_.isNotNull()
    ? Rcpp::as<arma::vec>(off_phi_.get()) : arma::vec();

  const uword p = X.n_cols;
  const uword q = Z.n_cols;
  const uword n = X.n_rows;
  if (theta.n_elem != p + q || X.n_rows != Z.n_rows) {
    Rcpp::stop("Non-conformable parameter vector and design matrices.");
  }
  if ((off_mu.n_elem != 0 && off_mu.n_elem != n) ||
      (off_phi.n_elem != 0 && off_phi.n_elem != n)) {
    Rcpp::stop("Offset length does not match the number of rows.");
  }

  vec eta_mu = X * theta.head(p);
  vec eta_phi = Z * theta.tail(q);
  if (off_mu.n_elem == n) eta_mu += off_mu;
  if (off_phi.n_elem == n) eta_phi += off_phi;
  vec mu(eta_mu.n_elem);
  vec phi(eta_phi.n_elem);

  for (uword i = 0; i < eta_mu.n_elem; ++i) {
    // Unclamped: fitted means keep their full dynamic range and agree with base
    // R (plogis/pnorm/...). The epsilon floor belongs to the likelihood path.
    if (!simplex_fast::mean_from_eta_exact(eta_mu[i], mean_link, mu[i])) {
      Rcpp::stop("The linear predictor is outside the valid domain of the selected mean link.");
    }
    phi[i] = simplex_fast::safe_exp(eta_phi[i]);
  }

  return List::create(
    Named("mu") = mu,
    Named("phi") = phi,
    Named("eta_mu") = eta_mu,
    Named("eta_phi") = eta_phi
  );
}


// Inverse mean-link transformation in C++.
// Maps a linear predictor eta to the mean in (0,1) using the mean link
// identified by mean_link (1 logit, 2 probit, 3 cloglog, 4 neglog).
// [[Rcpp::export]]
Rcpp::NumericVector simplex_linkinv_cpp(
    const Rcpp::NumericVector& eta,
    const int mean_link = 1) {

  const R_xlen_t n = eta.size();
  Rcpp::NumericVector out(n);
  for (R_xlen_t i = 0; i < n; ++i) {
    double mu = 0.0;
    if (Rcpp::NumericVector::is_na(eta[i])) { out[i] = eta[i]; continue; }
    // Unclamped, so simplex_linkinv() reproduces stats::plogis()/pnorm()/...
    // exactly instead of flooring at 1e-12.
    if (!simplex_fast::mean_from_eta_exact(eta[i], mean_link, mu)) {
      Rcpp::stop("Linear predictor outside the valid domain of the selected link.");
    }
    out[i] = mu;
  }
  return out;
}
