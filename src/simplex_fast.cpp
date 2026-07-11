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
      if (!mean_from_eta(eta_mu[i], mean_link, mu, dmu_deta)) {
        invalid = 1;
        continue;
      }

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
  }

  if (invalid != 0) {
    return {std::numeric_limits<double>::infinity(), vec(d, arma::fill::zeros), false};
  }

  double nll = 0.0;
  for (const double value : nll_local) nll += value;

  vec grad(d, arma::fill::zeros);
  if (need_grad) {
    for (const auto& g : grad_local) grad += g;
  }

  return {nll, std::move(grad), true};
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

} // namespace simplex_fast


// Fast simplex density in C++.
// Evaluates the simplex density (or log-density) at each y[i] given mu and phi,
// which are recycled when supplied with length one. Values outside the support
// (y or mu outside (0,1), phi <= 0, or non-finite) map to 0 (or -Inf on the log
// scale). The per-observation loop is optionally parallelized with OpenMP; it
// touches no R API state and is therefore thread-safe.
// [[Rcpp::export]]
Rcpp::NumericVector dsimplex_cpp(
    const Rcpp::NumericVector& y,
    const Rcpp::NumericVector& mu,
    const Rcpp::NumericVector& phi,
    const bool log = false,
    const int n_threads = 1) {

  const R_xlen_t n = y.size();
  if (!((mu.size() == 1 || mu.size() == n) && (phi.size() == 1 || phi.size() == n))) {
    Rcpp::stop("'mu' and 'phi' must have length 1 or length(y).");
  }

  Rcpp::NumericVector out(n);
  int threads = 1;
#ifdef _OPENMP
  threads = (n_threads > 0) ? n_threads : omp_get_max_threads();
#else
  (void)n_threads;
#endif

#ifdef _OPENMP
  #pragma omp parallel for num_threads(threads) schedule(static)
#endif
  for (R_xlen_t i = 0; i < n; ++i) {
    const double yi = y[i];
    const double mui = mu[(mu.size() == 1) ? 0 : i];
    const double phii = phi[(phi.size() == 1) ? 0 : i];

    if (!(yi > 0.0 && yi < 1.0 && mui > 0.0 && mui < 1.0 && phii > 0.0) ||
        !std::isfinite(yi) || !std::isfinite(mui) || !std::isfinite(phii)) {
      out[i] = log ? R_NegInf : 0.0;
      continue;
    }

    const double one_y = 1.0 - yi;
    const double qmu = mui * (1.0 - mui);
    const double diff = yi - mui;
    // Unit deviance dev = (y-mu)^2 / [y(1-y)(mu(1-mu))^2].
    const double dev = (diff * diff) / (yi * one_y * qmu * qmu);
    const double ld = -0.5 * (simplex_fast::LOG_2PI + std::log(phii))
                    -1.5 * (std::log(yi) + std::log(one_y))
                    -0.5 * dev / phii;
    out[i] = log ? ld : std::exp(ld);
  }

  return out;
}


// Fast random generation from the simplex distribution in C++.
// Uses the exact inverse-Gaussian-mixture transformation: with
// epsilon = mu/(1-mu) and tau = phi (1-mu)^2, a variate x is built from an
// inverse-Gaussian draw plus, with probability mu, a chi-squared(1) term; the
// result is mapped back to (0,1) via x/(1+x). Because it calls R's RNG, the
// loop is kept strictly serial (never parallelize R API calls).
// [[Rcpp::export]]
Rcpp::NumericVector rsimplex_cpp(
    const int n,
    const Rcpp::NumericVector& mu,
    const Rcpp::NumericVector& phi) {

  if (n < 0) Rcpp::stop("'n' must be non-negative.");
  if (!((mu.size() == 1 || mu.size() == n) && (phi.size() == 1 || phi.size() == n))) {
    Rcpp::stop("'mu' and 'phi' must have length 1 or n.");
  }

  Rcpp::RNGScope scope;
  Rcpp::NumericVector out(n);

  for (int i = 0; i < n; ++i) {
    const double mui = mu[(mu.size() == 1) ? 0 : i];
    const double phii = phi[(phi.size() == 1) ? 0 : i];
    if (!(mui > 0.0 && mui < 1.0 && phii > 0.0) ||
        !std::isfinite(mui) || !std::isfinite(phii)) {
      Rcpp::stop("All 'mu' values must lie in (0, 1) and all 'phi' values must be positive.");
    }

    const double epsilon = mui / (1.0 - mui);
    const double tau = phii * (1.0 - mui) * (1.0 - mui);

    const double x1 = simplex_fast::inv_gaussian_one(epsilon, tau);
    const double x3 = R::rchisq(1.0) * tau * epsilon * epsilon;
    const double x = (R::runif(0.0, 1.0) < mui) ? (x1 + x3) : x1;
    out[i] = x / (1.0 + x);
  }

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
    Named("valid") = res.valid
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
    double dmu_deta = 0.0;
    if (!simplex_fast::mean_from_eta(eta_mu[i], mean_link, mu[i], dmu_deta)) {
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
    double dmu_deta = 0.0;
    if (!simplex_fast::mean_from_eta(eta[i], mean_link, mu, dmu_deta)) {
      Rcpp::stop("Linear predictor outside the valid domain of the selected link.");
    }
    out[i] = mu;
  }
  return out;
}
