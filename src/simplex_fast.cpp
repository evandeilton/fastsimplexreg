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
EvalResult evaluate_impl(
    const vec& theta,
    const vec& y,
    const mat& X,
    const mat& Z,
    const int mean_link,
    const int n_threads,
    const bool need_grad = true) {

  const uword n = y.n_elem;
  const uword p = X.n_cols;
  const uword q = Z.n_cols;
  const uword d = p + q;

  if (theta.n_elem != d || X.n_rows != n || Z.n_rows != n) {
    return {std::numeric_limits<double>::infinity(), vec(d, arma::fill::zeros), false};
  }

  const vec beta = theta.head(p);
  const vec gamma = theta.tail(q);

  // BLAS-backed matrix-vector products for the two linear predictors.
  const vec eta_mu = X * beta;
  const vec eta_phi = Z * gamma;

  int threads = 1;
#ifdef _OPENMP
  threads = (n_threads > 0) ? n_threads : omp_get_max_threads();
  threads = std::max(1, threads);
#else
  (void)n_threads;
#endif

  // Per-thread accumulators to avoid data races; reduced after the region.
  std::vector<double> nll_local(static_cast<std::size_t>(threads), 0.0);
  std::vector<int> sat_local(static_cast<std::size_t>(threads), 0);
  std::vector<vec> grad_local;
  if (need_grad) {
    grad_local.reserve(static_cast<std::size_t>(threads));
    for (int t = 0; t < threads; ++t) {
      grad_local.emplace_back(d, arma::fill::zeros);
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
    double local_nll = 0.0;
    int local_sat = 0;
    vec* local_grad = need_grad ? &grad_local[static_cast<std::size_t>(tid)] : nullptr;

#ifdef _OPENMP
    #pragma omp for schedule(static)
#endif
    for (uword i = 0; i < n; ++i) {
      const double yi = y[i];
      if (!(yi > 0.0 && yi < 1.0) || !std::isfinite(yi)) {
        invalid = 1;
        continue;
      }

      double mu = 0.0;
      double dmu_deta = 0.0;
      bool sat = false;
      if (!mean_from_eta(eta_mu[i], mean_link, mu, dmu_deta, &sat)) {
        invalid = 1;
        continue;
      }
      if (sat) ++local_sat;

      // Dispersion link is log, so phi = exp(eta_phi) is strictly positive.
      const double phi = safe_exp(eta_phi[i]);
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

      const double loglik_i = -0.5 * (LOG_2PI + std::log(phi))
                            -1.5 * (std::log(yi) + std::log(one_y))
                            -0.5 * dev / phi;

      if (!std::isfinite(loglik_i)) {
        invalid = 1;
        continue;
      }
      // Accumulate the negative log-likelihood.
      local_nll -= loglik_i;

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
          (*local_grad)[j] -= X(i, j) * score_mu;
        }
        for (uword j = 0; j < q; ++j) {
          (*local_grad)[p + j] -= Z(i, j) * score_phi;
        }
      }
    }

    nll_local[static_cast<std::size_t>(tid)] = local_nll;
    sat_local[static_cast<std::size_t>(tid)] = local_sat;
  }

  if (invalid != 0) {
    return {std::numeric_limits<double>::infinity(), vec(d, arma::fill::zeros), false};
  }

  double nll = 0.0;
  for (const double value : nll_local) nll += value;

  int n_saturated = 0;
  for (const int value : sat_local) n_saturated += value;

  vec grad(d, arma::fill::zeros);
  if (need_grad) {
    for (const auto& g : grad_local) grad += g;
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

// Integrate the simplex density over [a, b] using an adaptive Gauss-Legendre
// rule seeded with knots placed at multiples of the first-order standard
// deviation around the mean. Seeding matters: for small phi the density is a
// narrow spike, and a naive adaptive rule started on the whole interval can
// bisect into two panels that both miss the spike, converge on a near-zero
// estimate and stop. The knots guarantee the peak is always resolved.
inline double simplex_integrate(const double a, const double b,
                                const double mu, const double phi,
                                const arma::vec& gln, const arma::vec& glw) {
  if (!(b > a)) return 0.0;
  const double u = mu * (1.0 - mu);
  double sd = std::sqrt(phi * u * u * u);          // Var(Y) ~ phi * V(mu)
  if (!(sd > 0.0) || !std::isfinite(sd)) sd = 0.1;

  static const double mult[] = {-12.0, -8.0, -6.0, -4.0, -3.0, -2.0, -1.5, -1.0,
                                -0.5, -0.25, 0.0, 0.25, 0.5, 1.0, 1.5, 2.0, 3.0,
                                4.0, 6.0, 8.0, 12.0};
  std::vector<double> knots;
  knots.reserve(24);
  knots.push_back(a);
  for (const double m : mult) {
    const double k = mu + m * sd;
    if (k > a && k < b) knots.push_back(k);
  }
  knots.push_back(b);
  std::sort(knots.begin(), knots.end());
  knots.erase(std::unique(knots.begin(), knots.end()), knots.end());

  const auto integrand = [&](const double t) { 
    return std::exp(simplex_logpdf(t, mu, phi));
  };

  double total = 0.0;
  for (std::size_t k = 0; k + 1 < knots.size(); ++k) {
    total += integrate_adaptive(integrand, knots[k], knots[k + 1], gln, glw,
                                1e-15, 1e-12, 40);
  }
  return total;
}

// Simplex CDF at a single point. The smaller tail is always the one integrated,
// so both tails keep full relative accuracy instead of being formed as
// 1 - (something close to 1).
inline void simplex_cdf_one(const double q, const double mu, const double phi,
                            const arma::vec& gln, const arma::vec& glw,
                            double& lower, double& upper) {
  if (q <= 0.0) { lower = 0.0; upper = 1.0; return; }
  if (q >= 1.0) { lower = 1.0; upper = 0.0; return; }
  if (q <= mu) {
    lower = simplex_integrate(0.0, q, mu, phi, gln, glw);
    lower = std::min(1.0, std::max(0.0, lower));
    upper = 1.0 - lower;
  } else {
    upper = simplex_integrate(q, 1.0, mu, phi, gln, glw);
    upper = std::min(1.0, std::max(0.0, upper));
    lower = 1.0 - upper;
  }
}

// Simplex quantile at a single lower-tail probability, by safeguarded
// Newton-bisection (Press et al., "rtsafe"): a Newton step is taken only when
// it stays inside the current bracket and halves the step, otherwise the method
// bisects. Bisection alone would converge; the Newton steps cut the number of
// (expensive) CDF evaluations to roughly a dozen.
inline double simplex_quantile_one(const double p, const double mu, const double phi,
                                   const arma::vec& gln, const arma::vec& glw) {
  if (!(p > 0.0)) return 0.0;
  if (!(p < 1.0)) return 1.0;

  double lo = 0.0, hi = 1.0;
  double y = mu;
  double dy_old = 1.0, dy = 1.0;

  for (int it = 0; it < 200; ++it) {
    double lower, upper;
    simplex_cdf_one(y, mu, phi, gln, glw, lower, upper);
    const double err = lower - p;
    const double dens = std::exp(simplex_logpdf(y, mu, phi));

    if (err > 0.0) hi = y; else lo = y;
    if (err == 0.0) return y;

    const bool newton_ok =
        (dens > 0.0) && std::isfinite(dens) &&
        (((y - hi) * dens - err) * ((y - lo) * dens - err) <= 0.0) &&
        (std::abs(2.0 * err) <= std::abs(dy_old * dens));

    dy_old = dy;
    if (newton_ok) {
      dy = err / dens;
      y -= dy;
    } else {
      dy = 0.5 * (hi - lo);
      y = lo + dy;
    }
    if (std::abs(dy) < 1e-15 || hi - lo < 1e-15) break;
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
#ifdef _OPENMP
  #pragma omp parallel for num_threads(threads) schedule(static) reduction(+:n_bad)
#endif
  for (R_xlen_t i = 0; i < n; ++i) {
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

  out.attr("n_invalid_par") = n_bad;
  return out;
}


// Simplex distribution function in C++.
// There is no closed form, so the density is integrated numerically with an
// adaptive Gauss-Legendre rule whose panels are seeded around the mean (see
// simplex_integrate). Argument conventions follow base R's p*(): `lower_tail`
// and `log_p` are honoured, NA/NaN propagate, and invalid parameters give NaN.
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

  arma::vec gln, glw;
  simplex_fast::gauss_legendre(15, gln, glw);

  Rcpp::NumericVector out(n);
  int threads = 1;
#ifdef _OPENMP
  threads = (n_threads > 0) ? n_threads : omp_get_max_threads();
  threads = std::max(1, threads);
#else
  (void)n_threads;
#endif

  int n_bad = 0;
#ifdef _OPENMP
  #pragma omp parallel for num_threads(threads) schedule(dynamic, 1) reduction(+:n_bad)
#endif
  for (R_xlen_t i = 0; i < n; ++i) {
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
    double lower = 0.0, upper = 1.0;
    simplex_fast::simplex_cdf_one(qi, mui, phii, gln, glw, lower, upper);
    const double val = lower_tail ? lower : upper;
    out[i] = log_p ? std::log(val) : val;
  }

  out.attr("n_invalid_par") = n_bad;
  return out;
}


// Simplex quantile function in C++.
// Inverts psimplex_cpp() by safeguarded Newton-bisection. Argument conventions
// follow base R's q*(): `lower_tail` and `log_p` are honoured, NA/NaN
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

  arma::vec gln, glw;
  simplex_fast::gauss_legendre(15, gln, glw);

  Rcpp::NumericVector out(n);
  int threads = 1;
#ifdef _OPENMP
  threads = (n_threads > 0) ? n_threads : omp_get_max_threads();
  threads = std::max(1, threads);
#else
  (void)n_threads;
#endif

  int n_bad = 0;
#ifdef _OPENMP
  #pragma omp parallel for num_threads(threads) schedule(dynamic, 1) reduction(+:n_bad)
#endif
  for (R_xlen_t i = 0; i < n; ++i) {
    double pi_ = p[i % npr];
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
    if (log_p) {
      if (pi_ > 0.0) { out[i] = R_NaN; ++n_bad; continue; }
      pi_ = std::exp(pi_);
    }
    if (!lower_tail) pi_ = 0.5 - pi_ + 0.5;   // 1 - p, guarding cancellation
    if (pi_ < 0.0 || pi_ > 1.0) { out[i] = R_NaN; ++n_bad; continue; }

    out[i] = simplex_fast::simplex_quantile_one(pi_, mui, phii, gln, glw);
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
    const int n_threads = 1) {

  const auto res = simplex_fast::evaluate_impl(theta, y, X, Z, mean_link, n_threads, true);
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
    const bool trace = false) {

  // Delegate to the shared native BFGS driver, wrapping the fixed-effects
  // evaluator as the objective. The optimizer logic is identical to before;
  // it now lives once in simplex_common.h and is reused by the mixed backend.
  auto objective = [&](const arma::vec& th) {
    return simplex_fast::evaluate_impl(th, y, X, Z, mean_link, n_threads, true);
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
    const int n_threads = 1) {

  const uword d = theta.n_elem;
  mat H(d, d, arma::fill::zeros);

  for (uword j = 0; j < d; ++j) {
    double h = rel_step * std::max(1.0, std::abs(theta[j]));
    bool success = false;

    for (int attempt = 0; attempt < 12; ++attempt) {
      vec plus = theta;
      vec minus = theta;
      plus[j] += h;
      minus[j] -= h;

      const auto gp = simplex_fast::evaluate_impl(plus, y, X, Z, mean_link, n_threads, true);
      const auto gm = simplex_fast::evaluate_impl(minus, y, X, Z, mean_link, n_threads, true);

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
    const int mean_link = 1) {

  const uword p = X.n_cols;
  const uword q = Z.n_cols;
  if (theta.n_elem != p + q || X.n_rows != Z.n_rows) {
    Rcpp::stop("Non-conformable parameter vector and design matrices.");
  }

  const vec eta_mu = X * theta.head(p);
  const vec eta_phi = Z * theta.tail(q);
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
