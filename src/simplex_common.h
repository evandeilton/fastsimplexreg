// simplex_common.h
// Shared, header-only numeric core for the fastsimplexreg C++ backend.
//
// Everything here is pure numeric (no R API state, no RNG) and is either
// `inline` (functions), `constexpr` (constants) or a template, so the header
// may be included by several translation units without violating the ODR.
// `compileAttributes()` never scans headers, so nothing here becomes an R
// export; the `// [[Rcpp::export]]` wrappers live in the .cpp files.
//
// The fixed-effects backend (simplex_fast.cpp) and the mixed-effects backend
// (simplex_mixed.cpp) both include this file and share the link functions, the
// per-observation simplex kernel, the native BFGS driver and the Gauss-Hermite
// rule generator.

#ifndef FASTSIMPLEXREG_SIMPLEX_COMMON_H
#define FASTSIMPLEXREG_SIMPLEX_COMMON_H

#include <RcppArmadillo.h>
#include <cmath>
#include <limits>
#include <algorithm>
#include <string>
#include <utility>
#include <vector>

namespace simplex_fast {

// log(2*pi); constant term of the Gaussian-like normalizing factor.
constexpr double LOG_2PI = 1.837877066409345483560659472811;
// Numerical floor used to keep probabilities strictly inside (0, 1).
constexpr double DEFAULT_EPS = 1e-12;
// 1 / sqrt(2*pi), the standard normal density peak (used by the probit link).
constexpr double INV_SQRT_2PI = 0.398942280401432677939946059934;
// 1 / sqrt(2), used with erfc to form the standard normal CDF.
constexpr double INV_SQRT_2 = 0.707106781186547524400844362105;
// sqrt(pi): total mass of the physicists' Gauss-Hermite weight (avoids the
// non-portable M_PI macro).
constexpr double SQRT_PI = 1.772453850905516027298167483341;

// Integer codes for the mean link functions. These MUST match the codes used
// on the R side: logit = 1, probit = 2, cloglog = 3, neglog = 4.
enum MeanLink : int {
  LOGIT = 1,
  PROBIT = 2,
  CLOGLOG = 3,
  NEGLOG = 4
};

// Clamp a probability to the closed interval [eps, 1 - eps] so that
// downstream logs and reciprocals stay finite.
inline double clamp_prob(const double x, const double eps = DEFAULT_EPS) noexcept {
  return std::min(1.0 - eps, std::max(eps, x));
}

// Clamp mu to [eps, 1 - eps] and, when clamping actually occurs, ZERO the
// supplied derivatives. At saturation the mean lies on a frozen (constant) part
// of the inverse-link surface, so d mu / d eta is genuinely zero there. Zeroing
// it keeps the analytic score and Hessian consistent with the (clamped)
// objective and removes the 1 / (mu (1-mu))^3 blow-up of the score that would
// otherwise make the analytic gradient disagree with the objective by many
// orders of magnitude, breaking the line search near the boundary.
inline void clamp_prob_deriv(double& mu, double& dmu) noexcept {
  if (mu <= DEFAULT_EPS) { mu = DEFAULT_EPS; dmu = 0.0; }
  else if (mu >= 1.0 - DEFAULT_EPS) { mu = 1.0 - DEFAULT_EPS; dmu = 0.0; }
}

inline void clamp_prob_deriv2(double& mu, double& dmu, double& d2mu) noexcept {
  if (mu <= DEFAULT_EPS) { mu = DEFAULT_EPS; dmu = 0.0; d2mu = 0.0; }
  else if (mu >= 1.0 - DEFAULT_EPS) { mu = 1.0 - DEFAULT_EPS; dmu = 0.0; d2mu = 0.0; }
}

// Numerically stable logistic (inverse-logit) function. The two branches
// avoid overflow of exp() for large-magnitude arguments of either sign.
inline double logistic_stable(const double x) noexcept {
  if (x >= 0.0) {
    const double z = std::exp(-x);
    return 1.0 / (1.0 + z);
  }
  const double z = std::exp(x);
  return z / (1.0 + z);
}

// exp() with the argument clamped to a safe range to avoid Inf/underflow.
// Used for the dispersion link (phi = exp(eta_phi)) and the cloglog/neglog
// mean links.
inline double safe_exp(const double x) noexcept {
  constexpr double LO = -700.0;
  constexpr double HI = 700.0;
  return std::exp(std::min(HI, std::max(LO, x)));
}

// Map a linear predictor eta to the mean mu via the chosen link and also
// return the derivative dmu/deta, needed by the chain rule in the score.
// Returns false when the transformation produces a non-finite value.
inline bool mean_from_eta(
    const double eta,
    const int mean_link,
    double& mu,
    double& dmu_deta) noexcept {

  switch (mean_link) {
    case LOGIT:
      // g(mu) = log(mu/(1-mu)); mu = 1/(1+exp(-eta)); dmu/deta = mu(1-mu).
      mu = logistic_stable(eta);
      dmu_deta = mu * (1.0 - mu);
      clamp_prob_deriv(mu, dmu_deta);
      return true;

    case PROBIT:
      // g(mu) = Phi^{-1}(mu); mu = Phi(eta) = 0.5*erfc(-eta/sqrt(2));
      // dmu/deta = phi(eta) = exp(-eta^2/2)/sqrt(2*pi).
      mu = 0.5 * std::erfc(-eta * INV_SQRT_2);
      dmu_deta = INV_SQRT_2PI * std::exp(-0.5 * eta * eta);
      clamp_prob_deriv(mu, dmu_deta);
      return std::isfinite(mu) && std::isfinite(dmu_deta);

    case CLOGLOG: {
      // g(mu) = log(-log(1-mu)); mu = 1 - exp(-exp(eta));
      // dmu/deta = exp(eta) * exp(-exp(eta)).
      const double exp_eta = safe_exp(eta);
      const double survival = std::exp(-exp_eta);
      mu = -std::expm1(-exp_eta);
      dmu_deta = exp_eta * survival;
      clamp_prob_deriv(mu, dmu_deta);
      return std::isfinite(mu) && std::isfinite(dmu_deta);
    }

    case NEGLOG: {
      // neglog link per Zhang et al. (2016):
      // g(mu) = -log(-log(mu)); mu = exp(-exp(-eta));
      // dmu/deta = mu * exp(-eta).
      const double exp_minus_eta = safe_exp(-eta);
      mu = std::exp(-exp_minus_eta);
      dmu_deta = mu * exp_minus_eta;
      clamp_prob_deriv(mu, dmu_deta);
      return std::isfinite(mu) && std::isfinite(dmu_deta);
    }

    default:
      return false;
  }
}

// Second-derivative-capable link map. In addition to mu and dmu/deta it returns
// d2mu/deta2, required by the observed second derivative of the log-likelihood
// with respect to eta_mu (needed by the mixed-model inner solver). The extra
// argument makes this a superset of mean_from_eta; the latter is kept unchanged
// so its existing callers are unaffected.
//
// Per-link second derivatives (verified against finite differences):
//   logit:   mu'' = mu' (1 - 2 mu)
//   probit:  mu'' = -eta mu'
//   cloglog: mu'' = mu' (1 - exp(eta))
//   neglog:  mu'' = mu' (exp(-eta) - 1)
inline bool mean_deriv2_from_eta(
    const double eta,
    const int mean_link,
    double& mu,
    double& dmu,
    double& d2mu) noexcept {

  switch (mean_link) {
    case LOGIT: {
      mu = logistic_stable(eta);
      dmu = mu * (1.0 - mu);
      d2mu = dmu * (1.0 - 2.0 * mu);
      clamp_prob_deriv2(mu, dmu, d2mu);
      return true;
    }
    case PROBIT: {
      mu = 0.5 * std::erfc(-eta * INV_SQRT_2);
      dmu = INV_SQRT_2PI * std::exp(-0.5 * eta * eta);
      d2mu = -eta * dmu;
      clamp_prob_deriv2(mu, dmu, d2mu);
      return std::isfinite(mu) && std::isfinite(dmu) && std::isfinite(d2mu);
    }
    case CLOGLOG: {
      const double exp_eta = safe_exp(eta);
      const double survival = std::exp(-exp_eta);
      mu = -std::expm1(-exp_eta);
      dmu = exp_eta * survival;
      d2mu = dmu * (1.0 - exp_eta);
      clamp_prob_deriv2(mu, dmu, d2mu);
      return std::isfinite(mu) && std::isfinite(dmu) && std::isfinite(d2mu);
    }
    case NEGLOG: {
      const double exp_minus_eta = safe_exp(-eta);
      mu = std::exp(-exp_minus_eta);
      dmu = mu * exp_minus_eta;
      d2mu = dmu * (exp_minus_eta - 1.0);
      clamp_prob_deriv2(mu, dmu, d2mu);
      return std::isfinite(mu) && std::isfinite(dmu) && std::isfinite(d2mu);
    }
    default:
      return false;
  }
}

// Per-observation simplex log-density and its derivatives with respect to the
// linear predictors, at a given (mu, dmu, d2mu, phi). All quantities are on the
// eta scale where relevant. The second derivative uses the division-free
// expanded form of d^2(dev)/d mu^2 so that it stays finite even when y == mu.
//
//   dev          = (y-mu)^2 / [ y(1-y) (mu(1-mu))^2 ]
//   dl/dmu       = (y-mu)(mu^2 - 2 mu y + y) / [ phi y(1-y) (mu(1-mu))^3 ]
//   dl/deta_mu   = dl/dmu * dmu
//   d2l/dmu^2    = -d''(mu) / (2 phi),   d''(mu) = d^2 dev / d mu^2
//   d2l/deta_mu2 = d2l/dmu^2 * dmu^2 + dl/dmu * d2mu     (observed)
//   dl/deta_phi  = -1/2 + dev/(2 phi)
struct ObsKernel {
  double logf;          // simplex log-density
  double dl_deta_mu;    // score wrt eta_mu
  double d2l_deta_mu2;  // observed second derivative wrt eta_mu
  double dl_deta_phi;   // score wrt eta_phi
  double dev;           // unit deviance
  bool ok;              // false when any quantity is non-finite / out of support
};

inline ObsKernel simplex_obs_kernel(
    const double y,
    const double mu,
    const double dmu,
    const double d2mu,
    const double phi) noexcept {

  ObsKernel k;
  k.ok = false;
  k.logf = k.dl_deta_mu = k.d2l_deta_mu2 = k.dl_deta_phi = k.dev = 0.0;

  if (!(y > 0.0 && y < 1.0) || !(phi > 0.0) || !std::isfinite(phi)) {
    return k;
  }

  const double one_y = 1.0 - y;
  const double one_mu = 1.0 - mu;
  const double v = y * one_y;          // y(1-y)
  const double u = mu * one_mu;        // mu(1-mu)
  const double u2 = u * u;
  const double u3 = u2 * u;
  const double u4 = u2 * u2;
  const double diff = y - mu;
  const double inv_v = 1.0 / v;

  k.dev = diff * diff * inv_v / u2;
  k.logf = -0.5 * (LOG_2PI + std::log(phi))
           - 1.5 * (std::log(y) + std::log(one_y))
           - 0.5 * k.dev / phi;

  // First derivative wrt mu (and chain rule to eta_mu). P = (y-mu)(mu^2-2 mu y+y).
  const double P = diff * (mu * mu - 2.0 * mu * y + y);
  const double dl_dmu = P * inv_v / (phi * u3);
  k.dl_deta_mu = dl_dmu * dmu;

  // Second derivative wrt mu via the expanded (division-by-diff-free) form of
  // d''(mu) = d^2 dev / d mu^2 = (-2/v) u^{-4} [ P'(mu) u - 3 P(mu) u' ],
  // with P'(mu) = -3 mu^2 + 6 mu y - 2 y^2 - y and u' = 1 - 2 mu.
  const double Pprime = -3.0 * mu * mu + 6.0 * mu * y - 2.0 * y * y - y;
  const double uprime = 1.0 - 2.0 * mu;
  const double ddev2 = (-2.0 * inv_v) * (Pprime * u - 3.0 * P * uprime) / u4;
  const double d2l_dmu2 = -0.5 * ddev2 / phi;
  k.d2l_deta_mu2 = d2l_dmu2 * dmu * dmu + dl_dmu * d2mu;

  // Dispersion score (log link): d log f / d eta_phi = -1/2 + dev/(2 phi).
  k.dl_deta_phi = -0.5 + 0.5 * k.dev / phi;

  k.ok = std::isfinite(k.logf) && std::isfinite(k.dl_deta_mu) &&
         std::isfinite(k.d2l_deta_mu2) && std::isfinite(k.dl_deta_phi);
  return k;
}

// Expected (Fisher) information for eta_mu, a strictly positive small-dispersion
// surrogate I_{eta_mu} = (dmu)^2 / (phi mu^3 (1-mu)^3). Used by the mixed-model
// inner solver as a guaranteed-SPD fallback for the observed curvature.
inline double simplex_fisher_eta_mu(const double mu, const double dmu, const double phi) noexcept {
  const double u = mu * (1.0 - mu);
  const double u3 = u * u * u;
  return (dmu * dmu) / (phi * u3);
}

// Bundle returned by an objective evaluation: the negative log-likelihood, its
// analytic gradient with respect to the full parameter vector, and a validity
// flag that is false when any observation falls outside the model support.
struct EvalResult {
  double nll;
  arma::vec grad;
  bool valid;
};

// Native BFGS minimizer with an Armijo backtracking line search, shared by both
// backends. `eval` is any callable taking an arma::vec theta and returning an
// EvalResult (nll + analytic gradient + validity). The algorithm is identical
// to the one previously inlined in simplex_bfgs_cpp: quasi-Newton direction
// -H*grad with an identity reset on a non-descent direction, Armijo line
// search, curvature-safeguarded inverse-Hessian update, and the same tolerances
// and messages. It returns the same Rcpp::List shape.
template <class Objective>
inline Rcpp::List bfgs_minimize(
    const arma::vec& start,
    Objective&& eval,
    const int maxit,
    const double rel_tol,
    const double grad_tol,
    const bool trace) {

  using arma::vec;
  using arma::mat;
  using arma::uword;

  vec theta = start;
  const uword d = theta.n_elem;
  mat H = arma::eye<mat>(d, d);

  EvalResult current = eval(theta);
  if (!current.valid || !std::isfinite(current.nll)) {
    Rcpp::stop("Initial parameter vector produces a non-finite objective.");
  }

  int convergence = 1;
  int iter_done = 0;
  int fn_evals = 1;
  int grad_evals = 1;
  double last_rel_change = std::numeric_limits<double>::infinity();
  std::string message = "Maximum number of iterations reached.";

  for (int iter = 0; iter < maxit; ++iter) {
    iter_done = iter + 1;
    if (arma::abs(current.grad).max() <= grad_tol) {
      convergence = 0;
      message = "Converged: gradient tolerance satisfied.";
      break;
    }

    vec direction = -H * current.grad;
    double slope = arma::dot(current.grad, direction);
    if (!std::isfinite(slope) || slope >= -1e-14) {
      H.eye();
      direction = -current.grad;
      slope = -arma::dot(current.grad, current.grad);
    }

    constexpr double c1 = 1e-4;
    constexpr double shrink = 0.5;
    constexpr int max_ls = 40;
    double step = 1.0;
    EvalResult candidate;
    bool accepted = false;

    for (int ls = 0; ls < max_ls; ++ls) {
      candidate = eval(theta + step * direction);
      ++fn_evals;
      ++grad_evals;
      if (candidate.valid && std::isfinite(candidate.nll) &&
          candidate.nll <= current.nll + c1 * step * slope) {
        accepted = true;
        break;
      }
      step *= shrink;
    }

    if (!accepted) {
      // Soft convergence: if the objective was already numerically stationary
      // on the previous accepted step, the line search cannot improve because we
      // have reached the objective's floor (e.g. the adaptive-quadrature noise
      // floor of the marginal likelihood, where the gradient tolerance is
      // unreachable). Report this as convergence rather than a hard failure.
      if (iter > 0 && last_rel_change <= rel_tol) {
        convergence = 0;
        message = "Converged: objective stationary (line search reached its floor).";
      } else {
        convergence = 2;
        message = "Line search failed to find a sufficient decrease.";
      }
      break;
    }

    const vec theta_new = theta + step * direction;
    const vec s = theta_new - theta;
    const vec yk = candidate.grad - current.grad;
    const double ys = arma::dot(yk, s);

    if (ys > 1e-12 * arma::norm(s, 2) * arma::norm(yk, 2)) {
      const double rho = 1.0 / ys;
      const mat I = arma::eye<mat>(d, d);
      const mat V = I - rho * s * yk.t();
      H = V * H * V.t() + rho * s * s.t();
      H = 0.5 * (H + H.t());
    } else {
      H.eye();
    }

    const double rel_change = std::abs(current.nll - candidate.nll) /
                              (1.0 + std::abs(current.nll));
    last_rel_change = rel_change;

    theta = theta_new;
    current = std::move(candidate);

    if (trace) {
      Rcpp::Rcout << "iter=" << iter_done
                  << " nll=" << current.nll
                  << " grad_inf=" << arma::abs(current.grad).max()
                  << " step=" << step << '\n';
    }

    if (rel_change <= rel_tol && arma::abs(current.grad).max() <= std::sqrt(grad_tol)) {
      convergence = 0;
      message = "Converged: relative objective tolerance satisfied.";
      break;
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("par") = theta,
    Rcpp::Named("value") = current.nll,
    Rcpp::Named("gradient") = current.grad,
    Rcpp::Named("inverse_hessian_bfgs") = H,
    Rcpp::Named("convergence") = convergence,
    Rcpp::Named("message") = message,
    Rcpp::Named("iterations") = iter_done,
    Rcpp::Named("function_evaluations") = fn_evals,
    Rcpp::Named("gradient_evaluations") = grad_evals
  );
}

// Gauss-Hermite (physicists' weight e^{-t^2}) nodes and weights via the
// Golub-Welsch algorithm: the nodes are the eigenvalues of the symmetric
// tridiagonal Jacobi matrix (zero diagonal, off-diagonal sqrt(m/2)), and the
// weights are sqrt(pi) times the squared first component of each normalized
// eigenvector. Exact for polynomials of degree <= 2*nAGQ - 1; sum of weights
// equals sqrt(pi). This avoids any dependency on tabulated rules.
inline void gauss_hermite(const int nAGQ, arma::vec& nodes, arma::vec& weights) {
  const int M = std::max(1, nAGQ);
  if (M == 1) {
    nodes = arma::vec(1, arma::fill::zeros);
    weights = arma::vec(1);
    weights[0] = SQRT_PI;
    return;
  }
  arma::mat J(M, M, arma::fill::zeros);
  for (int m = 1; m < M; ++m) {
    const double b = std::sqrt(m / 2.0);
    J(m - 1, m) = b;
    J(m, m - 1) = b;
  }
  arma::vec eval;
  arma::mat evec;
  arma::eig_sym(eval, evec, J);
  nodes = eval;                                   // ascending order
  weights = SQRT_PI * arma::square(evec.row(0)).t();
}

} // namespace simplex_fast

#endif // FASTSIMPLEXREG_SIMPLEX_COMMON_H
