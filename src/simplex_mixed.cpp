// simplex_mixed.cpp
// Native C++ backend for the two-level simplex mixed model with variable
// dispersion, estimated by adaptive Gauss-Hermite quadrature (AGHQ).
//
// Design conventions (see src/simplex_common.h for the shared numeric core):
//   X = mean fixed-effects design (N x p)
//   Z = random-effects design     (N x q),  q in 1..3
//   W = dispersion design         (N x r)
// (Note: this differs from simplex_fast.cpp, where Z denotes the dispersion
// design. Here we follow the mixed-model convention X=mean, Z=random,
// W=dispersion.)
//
// Model, conditional on the cluster random effect b_j ~ N_q(0, Sigma):
//   g(mu_ij)   = x_ij' beta + z_ij' b_j        (mean; link = mean_link)
//   log phi_ij = w_ij' gamma                    (dispersion; log link)
// Parameter vector theta = c(beta, gamma, omega), where omega is the
// unconstrained log-Cholesky packing of the lower-triangular factor D of
// Sigma = D D'. Rows are supplied group-contiguous with CSR offsets `starts`.

// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include "simplex_common.h"
#include <cmath>
#include <limits>
#include <vector>
#ifdef _OPENMP
  #include <omp.h>
#endif

using arma::mat;
using arma::vec;
using arma::uvec;
using arma::uword;
using Rcpp::List;
using Rcpp::Named;

namespace simplex_fast {

// ---------------------------------------------------------------------------
// Unconstrained <-> Cholesky-factor packing.
// omega packs the lower-triangular D column by column: for each column c the
// log of the (positive) diagonal entry first, then the free sub-diagonal
// entries. Sigma = D D' is SPD for every omega.
// ---------------------------------------------------------------------------
inline mat build_D(const vec& omega, const int q) {
  mat D(q, q, arma::fill::zeros);
  int idx = 0;
  for (int c = 0; c < q; ++c) {
    D(c, c) = std::exp(omega[idx++]);
    for (int rr = c + 1; rr < q; ++rr) {
      D(rr, c) = omega[idx++];
    }
  }
  return D;
}

inline vec pack_omega(const mat& D, const int q) {
  const int m = q * (q + 1) / 2;
  vec omega(m);
  int idx = 0;
  for (int c = 0; c < q; ++c) {
    omega[idx++] = std::log(D(c, c));
    for (int rr = c + 1; rr < q; ++rr) {
      omega[idx++] = D(rr, c);
    }
  }
  return omega;
}

// Gradient of log N(b; 0, Sigma) with respect to omega, evaluated at b.
// Uses M = 2 G D = -(Sigma^{-1} - Sigma^{-1} b b' Sigma^{-1}) D and reads off
// the packed entries: diagonal omega gets D(c,c) * M(c,c), off-diagonal gets
// M(r,c). Writes into `out` (length m).
inline void omega_grad(const vec& b, const mat& Sigma_inv, const mat& D,
                       const int q, vec& out) {
  const vec Sib = Sigma_inv * b;                       // Sigma^{-1} b
  const mat Gmat = -0.5 * (Sigma_inv - Sib * Sib.t()); // dlogN/dSigma
  const mat Mmat = 2.0 * Gmat * D;
  int idx = 0;
  for (int c = 0; c < q; ++c) {
    out[idx++] = D(c, c) * Mmat(c, c);
    for (int rr = c + 1; rr < q; ++rr) {
      out[idx++] = Mmat(rr, c);
    }
  }
}

// Build the tensor-product Gauss-Hermite rule for dimension q: node matrix
// T (K x q), tensor log-weights logW (K) and squared norms t2 (K), K = M^q.
inline void build_tensor(const int nAGQ, const int q,
                         mat& T, vec& logW, vec& t2) {
  vec nodes, wts;
  gauss_hermite(nAGQ, nodes, wts);
  const int M = nodes.n_elem;
  const vec logw = arma::log(wts);

  uword K = 1;
  for (int d = 0; d < q; ++d) K *= static_cast<uword>(M);

  T.set_size(K, q);
  logW.set_size(K);
  t2.set_size(K);

  std::vector<int> mi(q, 0);              // multi-index odometer
  for (uword k = 0; k < K; ++k) {
    double lw = 0.0, s2 = 0.0;
    for (int d = 0; d < q; ++d) {
      const double t = nodes[mi[d]];
      T(k, d) = t;
      lw += logw[mi[d]];
      s2 += t * t;
    }
    logW[k] = lw;
    t2[k] = s2;
    for (int d = 0; d < q; ++d) {          // increment odometer
      if (++mi[d] < M) break;
      mi[d] = 0;
    }
  }
}

// Log-sum-exp of a vector.
inline double log_sum_exp(const vec& a) {
  const double amax = a.max();
  if (!std::isfinite(amax)) return amax;
  return amax + std::log(arma::sum(arma::exp(a - amax)));
}

// ---------------------------------------------------------------------------
// Core AGHQ evaluator: marginal negative log-likelihood and analytic gradient.
// Updates `Bhat` (J x q) with the per-cluster empirical-Bayes modes (warm
// start in / out). Parallel over clusters.
// ---------------------------------------------------------------------------
EvalResult mixed_core(
    const vec& theta, const vec& y, const mat& X, const mat& Z, const mat& W,
    const uvec& starts, const int q, const int mean_link,
    const mat& T, const vec& logW, const vec& t2,
    const int n_threads, const int inner_maxit, const double inner_tol,
    const bool need_grad, mat& Bhat) {

  const uword p = X.n_cols;
  const uword r = W.n_cols;
  const uword m = static_cast<uword>(q) * (q + 1) / 2;
  const uword dim = p + r + m;
  const uword J = starts.n_elem - 1;
  const uword K = T.n_rows;

  auto fail = [&]() {
    return EvalResult{std::numeric_limits<double>::infinity(),
                      vec(dim, arma::fill::zeros), false};
  };

  if (theta.n_elem != dim) return fail();

  const vec beta = theta.head(p);
  const vec gamma = theta.subvec(p, p + r - 1);
  const vec omega = theta.tail(m);

  // Sigma = D D'; keep Sigma_inv and log|Sigma| from the Cholesky factor D.
  const mat D = build_D(omega, q);
  for (int c = 0; c < q; ++c) {
    if (!(D(c, c) > 0.0) || !std::isfinite(D(c, c))) return fail();
  }
  const mat Dinv = arma::inv(arma::trimatl(D));
  const mat Sigma_inv = Dinv.t() * Dinv;
  double logdetSigma = 0.0;
  for (int c = 0; c < q; ++c) logdetSigma += 2.0 * std::log(D(c, c));

  const double sqrt2 = std::sqrt(2.0);
  const double half_q_log2 = 0.5 * q * std::log(2.0);
  const double half_q_log2pi = 0.5 * q * LOG_2PI;

  int threads = 1;
#ifdef _OPENMP
  threads = (n_threads > 0) ? n_threads : omp_get_max_threads();
  threads = std::max(1, threads);
#else
  (void)n_threads;
#endif

  std::vector<double> nll_local(static_cast<std::size_t>(threads), 0.0);
  std::vector<vec> grad_local;
  if (need_grad) {
    grad_local.reserve(static_cast<std::size_t>(threads));
    for (int t = 0; t < threads; ++t) grad_local.emplace_back(dim, arma::fill::zeros);
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
    vec* lg = need_grad ? &grad_local[static_cast<std::size_t>(tid)] : nullptr;

#ifdef _OPENMP
    #pragma omp for schedule(dynamic, 8)
#endif
    for (uword j = 0; j < J; ++j) {
      if (invalid) continue;
      const uword a = starts[j];
      const uword bb = starts[j + 1];
      const uword nj = bb - a;

      const mat Xj = X.rows(a, bb - 1);
      const mat Zj = Z.rows(a, bb - 1);
      const mat Wj = W.rows(a, bb - 1);
      const vec yj = y.subvec(a, bb - 1);

      const vec eta_mu_fixed = Xj * beta;
      const vec eta_phi = Wj * gamma;
      vec phi(nj);
      for (uword i = 0; i < nj; ++i) phi[i] = safe_exp(eta_phi[i]);

      // ---- inner Newton (Fisher scoring) for the posterior mode ----
      vec b = Bhat.row(j).t();

      auto hval = [&](const vec& bv, bool& ok) -> double {
        const vec eta = eta_mu_fixed + Zj * bv;
        double s = 0.0;
        ok = true;
        for (uword i = 0; i < nj; ++i) {
          double mu, dmu;
          if (!mean_from_eta(eta[i], mean_link, mu, dmu)) { ok = false; return 0.0; }
          const double one_y = 1.0 - yj[i];
          const double u = mu * (1.0 - mu);
          const double diff = yj[i] - mu;
          const double dev = diff * diff / (yj[i] * one_y * u * u);
          const double lf = -0.5 * (LOG_2PI + std::log(phi[i]))
                            - 1.5 * (std::log(yj[i]) + std::log(one_y))
                            - 0.5 * dev / phi[i];
          if (!std::isfinite(lf)) { ok = false; return 0.0; }
          s += lf;
        }
        s -= 0.5 * arma::dot(bv, Sigma_inv * bv);
        return s;
      };

      bool ok_cluster = true;
      for (int it = 0; it < inner_maxit; ++it) {
        const vec eta = eta_mu_fixed + Zj * b;
        vec s_mu(nj);
        vec Iinfo(nj);
        bool ok = true;
        for (uword i = 0; i < nj; ++i) {
          double mu, dmu;
          if (!mean_from_eta(eta[i], mean_link, mu, dmu)) { ok = false; break; }
          const double u = mu * (1.0 - mu);
          const double u3 = u * u * u;
          const double one_y = 1.0 - yj[i];
          const double diff = yj[i] - mu;
          const double Pterm = diff * (mu * mu - 2.0 * mu * yj[i] + yj[i]);
          s_mu[i] = (Pterm / (yj[i] * one_y)) / (phi[i] * u3) * dmu;   // dl/deta_mu
          Iinfo[i] = (dmu * dmu) / (phi[i] * u3);                      // Fisher info
        }
        if (!ok) { ok_cluster = false; break; }

        const vec g = Zj.t() * s_mu - Sigma_inv * b;
        if (arma::abs(g).max() < inner_tol) break;

        mat Qf = Sigma_inv;
        Qf += Zj.t() * (Zj.each_col() % Iinfo);   // Z' diag(I) Z + Sigma_inv (SPD)
        vec delta;
        if (!arma::solve(delta, Qf, g, arma::solve_opts::likely_sympd)) { ok_cluster = false; break; }

        // step-halving line search on h_j
        double alpha = 1.0;
        bool okb = false;
        double h0 = hval(b, okb);
        vec bnew;
        int hs = 0;
        for (; hs < 30; ++hs) {
          bnew = b + alpha * delta;
          bool okn = false;
          const double h1 = hval(bnew, okn);
          if (okn && h1 >= h0 - 1e-12) break;
          alpha *= 0.5;
        }
        b = bnew;
      }
      if (!ok_cluster) { invalid = 1; continue; }
      Bhat.row(j) = b.t();

      // ---- observed curvature Q_j at the mode (fall back to Fisher) ----
      mat Q(q, q, arma::fill::zeros);
      {
        const vec eta = eta_mu_fixed + Zj * b;
        vec w2(nj);       // observed d2l/deta_mu2
        vec Iinfo(nj);
        bool ok = true;
        for (uword i = 0; i < nj; ++i) {
          double mu, dmu, d2mu;
          if (!mean_deriv2_from_eta(eta[i], mean_link, mu, dmu, d2mu)) { ok = false; break; }
          const ObsKernel kk = simplex_obs_kernel(yj[i], mu, dmu, d2mu, phi[i]);
          if (!kk.ok) { ok = false; break; }
          w2[i] = kk.d2l_deta_mu2;
          Iinfo[i] = (dmu * dmu) / (phi[i] * (mu * (1.0 - mu)) * (mu * (1.0 - mu)) * (mu * (1.0 - mu)));
        }
        if (!ok) { invalid = 1; continue; }
        Q = Sigma_inv - Zj.t() * (Zj.each_col() % w2);   // -Hessian (observed)
        mat Rchk;
        if (!arma::chol(Rchk, Q)) {
          Q = Sigma_inv + Zj.t() * (Zj.each_col() % Iinfo);  // Fisher fallback (SPD)
        }
      }

      mat R;
      if (!arma::chol(R, Q)) { invalid = 1; continue; }   // Q = R' R (upper R)
      double logdetQ = 0.0;
      for (int d = 0; d < q; ++d) logdetQ += 2.0 * std::log(R(d, d));
      const mat C = arma::inv(arma::trimatu(R));           // C C' = Q^{-1}, |C| = |Q|^{-1/2}

      // ---- AGHQ sweep: unnormalized log-weights a_k ----
      vec avec(K);
      mat Bnodes(q, K);
      bool ok_nodes = true;
      for (uword k = 0; k < K; ++k) {
        const vec tk = T.row(k).t();
        const vec bk = b + sqrt2 * (C * tk);
        Bnodes.col(k) = bk;
        const vec eta = eta_mu_fixed + Zj * bk;
        double sumlogf = 0.0;
        for (uword i = 0; i < nj; ++i) {
          double mu, dmu;
          if (!mean_from_eta(eta[i], mean_link, mu, dmu)) { ok_nodes = false; break; }
          const double one_y = 1.0 - yj[i];
          const double u = mu * (1.0 - mu);
          const double diff = yj[i] - mu;
          const double dev = diff * diff / (yj[i] * one_y * u * u);
          const double lf = -0.5 * (LOG_2PI + std::log(phi[i]))
                            - 1.5 * (std::log(yj[i]) + std::log(one_y))
                            - 0.5 * dev / phi[i];
          sumlogf += lf;
        }
        if (!ok_nodes) break;
        const double quad = 0.5 * arma::dot(bk, Sigma_inv * bk);
        avec[k] = logW[k] + t2[k] + sumlogf - quad;
      }
      if (!ok_nodes) { invalid = 1; continue; }

      const double lse = log_sum_exp(avec);
      const double loglik_j = half_q_log2 - 0.5 * logdetQ
                              - half_q_log2pi - 0.5 * logdetSigma + lse;
      if (!std::isfinite(loglik_j)) { invalid = 1; continue; }
      local_nll -= loglik_j;

      // ---- posterior-weighted analytic gradient ----
      if (need_grad) {
        const vec Wjk = arma::exp(avec - lse);   // normalized posterior weights
        vec gbeta(p, arma::fill::zeros);
        vec ggamma(r, arma::fill::zeros);
        vec gomega(m, arma::fill::zeros);
        vec og(m);
        for (uword k = 0; k < K; ++k) {
          const double wk = Wjk[k];
          if (wk <= 0.0) continue;
          const vec bk = Bnodes.col(k);
          const vec eta = eta_mu_fixed + Zj * bk;
          vec s_mu(nj), s_phi(nj);
          for (uword i = 0; i < nj; ++i) {
            double mu, dmu;
            mean_from_eta(eta[i], mean_link, mu, dmu);
            const double u = mu * (1.0 - mu);
            const double u3 = u * u * u;
            const double one_y = 1.0 - yj[i];
            const double diff = yj[i] - mu;
            const double dev = diff * diff / (yj[i] * one_y * u * u);
            const double Pterm = diff * (mu * mu - 2.0 * mu * yj[i] + yj[i]);
            s_mu[i] = (Pterm / (yj[i] * one_y)) / (phi[i] * u3) * dmu;
            s_phi[i] = -0.5 + 0.5 * dev / phi[i];
          }
          gbeta += wk * (Xj.t() * s_mu);
          ggamma += wk * (Wj.t() * s_phi);
          omega_grad(bk, Sigma_inv, D, q, og);
          gomega += wk * og;
        }
        // gradient of the NEGATIVE log-likelihood
        lg->subvec(0, p - 1) -= gbeta;
        lg->subvec(p, p + r - 1) -= ggamma;
        lg->subvec(p + r, dim - 1) -= gomega;
      }
    }

    nll_local[static_cast<std::size_t>(tid)] = local_nll;
  }

  if (invalid != 0) return fail();

  double nll = 0.0;
  for (const double v : nll_local) nll += v;
  vec grad(dim, arma::fill::zeros);
  if (need_grad) for (const auto& g : grad_local) grad += g;

  return EvalResult{nll, std::move(grad), true};
}

} // namespace simplex_fast


// Marginal negative log-likelihood and analytic gradient at theta (AGHQ).
// [[Rcpp::export]]
Rcpp::List simplex_mixed_eval_cpp(
    const arma::vec& theta, const arma::vec& y,
    const arma::mat& X, const arma::mat& Z, const arma::mat& W,
    const arma::uvec& starts, const int q,
    const int mean_link = 1, const int nAGQ = 11,
    const int n_threads = 1, const int inner_maxit = 50,
    const double inner_tol = 1e-8) {

  mat T; vec logW, t2;
  simplex_fast::build_tensor(nAGQ, q, T, logW, t2);
  mat Bhat(starts.n_elem - 1, q, arma::fill::zeros);
  const auto res = simplex_fast::mixed_core(theta, y, X, Z, W, starts, q, mean_link,
                                            T, logW, t2, n_threads, inner_maxit,
                                            inner_tol, true, Bhat);
  return List::create(Named("value") = res.nll,
                      Named("gradient") = res.grad,
                      Named("valid") = res.valid);
}


// Native BFGS on the marginal NLL (reuses the shared bfgs_minimize driver, with
// warm-started per-cluster modes across evaluations).
// [[Rcpp::export]]
Rcpp::List simplex_mixed_bfgs_cpp(
    const arma::vec& start, const arma::vec& y,
    const arma::mat& X, const arma::mat& Z, const arma::mat& W,
    const arma::uvec& starts, const int q,
    const int mean_link = 1, const int nAGQ = 11,
    const int maxit = 300, const double rel_tol = 1e-9,
    const double grad_tol = 1e-6, const int n_threads = 1,
    const int inner_maxit = 50, const double inner_tol = 1e-8,
    const bool trace = false) {

  mat T; vec logW, t2;
  simplex_fast::build_tensor(nAGQ, q, T, logW, t2);
  mat Bhat(starts.n_elem - 1, q, arma::fill::zeros);

  auto objective = [&](const arma::vec& th) {
    return simplex_fast::mixed_core(th, y, X, Z, W, starts, q, mean_link,
                                    T, logW, t2, n_threads, inner_maxit,
                                    inner_tol, true, Bhat);
  };
  return simplex_fast::bfgs_minimize(start, objective, maxit, rel_tol, grad_tol, trace);
}


// Finite-difference Hessian of the analytic marginal gradient (central diff).
// [[Rcpp::export]]
arma::mat simplex_mixed_hessian_fd_cpp(
    const arma::vec& theta, const arma::vec& y,
    const arma::mat& X, const arma::mat& Z, const arma::mat& W,
    const arma::uvec& starts, const int q,
    const int mean_link = 1, const int nAGQ = 11,
    const double rel_step = 1e-5, const int n_threads = 1,
    const int inner_maxit = 50, const double inner_tol = 1e-8) {

  mat T; vec logW, t2;
  simplex_fast::build_tensor(nAGQ, q, T, logW, t2);
  const uword d = theta.n_elem;
  mat Bhat(starts.n_elem - 1, q, arma::fill::zeros);

  auto grad_at = [&](const arma::vec& th) {
    return simplex_fast::mixed_core(th, y, X, Z, W, starts, q, mean_link,
                                    T, logW, t2, n_threads, inner_maxit,
                                    inner_tol, true, Bhat).grad;
  };

  mat H(d, d, arma::fill::zeros);
  for (uword jcol = 0; jcol < d; ++jcol) {
    double h = rel_step * std::max(1.0, std::abs(theta[jcol]));
    bool success = false;
    for (int attempt = 0; attempt < 12; ++attempt) {
      vec plus = theta, minus = theta;
      plus[jcol] += h; minus[jcol] -= h;
      const vec gp = grad_at(plus);
      const vec gm = grad_at(minus);
      if (gp.is_finite() && gm.is_finite()) {
        H.col(jcol) = (gp - gm) / (2.0 * h);
        success = true;
        break;
      }
      h *= 0.25;
    }
    if (!success) Rcpp::stop("Non-finite gradient while forming the mixed-model Hessian.");
  }
  return 0.5 * (H + H.t());
}


// Empirical-Bayes modes b_hat_j (J x q) and posterior covariances Q_j^{-1}
// (q x q x J) for ranef()/predict().
// [[Rcpp::export]]
Rcpp::List simplex_mixed_ranef_cpp(
    const arma::vec& theta, const arma::vec& y,
    const arma::mat& X, const arma::mat& Z, const arma::mat& W,
    const arma::uvec& starts, const int q,
    const int mean_link = 1, const int n_threads = 1,
    const int inner_maxit = 50, const double inner_tol = 1e-8) {

  // One AGHQ node reproduces the mode-finding path; then read the modes back.
  mat T; vec logW, t2;
  simplex_fast::build_tensor(1, q, T, logW, t2);
  const uword J = starts.n_elem - 1;
  mat Bhat(J, q, arma::fill::zeros);
  const auto res = simplex_fast::mixed_core(theta, y, X, Z, W, starts, q, mean_link,
                                            T, logW, t2, n_threads, inner_maxit,
                                            inner_tol, false, Bhat);
  if (!res.valid) Rcpp::stop("Random-effects prediction produced a non-finite value.");

  // Posterior covariances: recompute Q_j^{-1} at the modes.
  const uword p = X.n_cols, r = W.n_cols, m = static_cast<uword>(q) * (q + 1) / 2;
  const vec beta = theta.head(p);
  const vec gamma = theta.subvec(p, p + r - 1);
  const vec omega = theta.tail(m);
  const mat D = simplex_fast::build_D(omega, q);
  const mat Dinv = arma::inv(arma::trimatl(D));
  const mat Sigma_inv = Dinv.t() * Dinv;

  arma::cube postvar(q, q, J);
  for (uword j = 0; j < J; ++j) {
    const uword a = starts[j], bb = starts[j + 1], nj = bb - a;
    const mat Xj = X.rows(a, bb - 1);
    const mat Zj = Z.rows(a, bb - 1);
    const mat Wj = W.rows(a, bb - 1);
    const vec yj = y.subvec(a, bb - 1);
    const vec eta_phi = Wj * gamma;
    const vec bmode = Bhat.row(j).t();
    const vec eta = Xj * beta + Zj * bmode;
    vec w2(nj), Iinfo(nj);
    bool okobs = true;
    for (uword i = 0; i < nj; ++i) {
      const double phii = simplex_fast::safe_exp(eta_phi[i]);
      double mu, dmu, d2mu;
      simplex_fast::mean_deriv2_from_eta(eta[i], mean_link, mu, dmu, d2mu);
      const simplex_fast::ObsKernel kk = simplex_fast::simplex_obs_kernel(yj[i], mu, dmu, d2mu, phii);
      if (!kk.ok) { okobs = false; break; }
      w2[i] = kk.d2l_deta_mu2;
      const double u = mu * (1.0 - mu);
      Iinfo[i] = (dmu * dmu) / (phii * u * u * u);
    }
    mat Q = okobs ? (mat)(Sigma_inv - Zj.t() * (Zj.each_col() % w2)) : Sigma_inv;
    mat Rchk;
    if (!okobs || !arma::chol(Rchk, Q)) {
      Q = Sigma_inv + Zj.t() * (Zj.each_col() % Iinfo);
    }
    postvar.slice(j) = arma::inv_sympd(Q);
  }

  return List::create(Named("b") = Bhat, Named("postvar") = postvar);
}


// Fitted mu/phi and linear predictors, optionally including the random effects.
// `b` is the J x q matrix of modes aligned to clusters; when include_re is
// false the random-effect contribution is dropped (population level).
// [[Rcpp::export]]
Rcpp::List simplex_mixed_predict_cpp(
    const arma::vec& theta,
    const arma::mat& X, const arma::mat& Z, const arma::mat& W,
    const arma::uvec& starts, const int q,
    const arma::mat& b,
    const int mean_link = 1, const bool include_re = true) {

  const uword p = X.n_cols, r = W.n_cols, N = X.n_rows;
  const uword J = starts.n_elem - 1;
  const vec beta = theta.head(p);
  const vec gamma = theta.subvec(p, p + r - 1);

  vec eta_mu = X * beta;
  if (include_re) {
    for (uword j = 0; j < J; ++j) {
      const uword a = starts[j], bb = starts[j + 1];
      if (bb > a) eta_mu.subvec(a, bb - 1) += Z.rows(a, bb - 1) * b.row(j).t();
    }
  }
  const vec eta_phi = W * gamma;

  vec mu(N), phi(N);
  for (uword i = 0; i < N; ++i) {
    double m_, dmu;
    if (!simplex_fast::mean_from_eta(eta_mu[i], mean_link, m_, dmu)) {
      Rcpp::stop("Linear predictor outside the valid domain of the selected mean link.");
    }
    mu[i] = m_;
    phi[i] = simplex_fast::safe_exp(eta_phi[i]);
  }
  return List::create(Named("mu") = mu, Named("phi") = phi,
                      Named("eta_mu") = eta_mu, Named("eta_phi") = eta_phi);
}


// Build the Cholesky factor D from the unconstrained omega packing.
// [[Rcpp::export]]
arma::mat simplex_mixed_D_from_omega_cpp(const arma::vec& omega, const int q) {
  return simplex_fast::build_D(omega, q);
}

// Recover the unconstrained omega packing from a lower-triangular D.
// [[Rcpp::export]]
arma::vec simplex_mixed_omega_from_D_cpp(const arma::mat& D) {
  return simplex_fast::pack_omega(D, D.n_rows);
}
