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

  // NO node pruning. The previous version dropped nodes whose product
  // Gauss-Hermite weight logW sat below a relative floor, on the reasoning that
  // the corners of the tensor grid "contribute nothing". That reasoning is
  // wrong for ADAPTIVE Gauss-Hermite and was measurably harmful.
  //
  // The adaptive transformation undoes the e^{-t^2} factor, so what a node k
  // actually multiplies into the log-sum-exp below is logW[k] + t2[k], not
  // logW[k]. Outer nodes carry a tiny weight and an enormous e^{t^2} and the
  // two very nearly cancel: across the WHOLE tensor grid the effective
  // multiplier spans a factor of only 1.6 (nAGQ = 5), 2.4 (11), 3.5 (21) and
  // 5.2 (41). There is no negligible tail to drop -- every node matters to
  // within a factor of five.
  //
  // Pruning on logW alone therefore discarded nodes carrying 5.7% of the
  // effective quadrature mass at nAGQ = 11, 48.7% at nAGQ = 21 and 66.3% at
  // nAGQ = 31, which made the AGHQ sequence stop converging: raising nAGQ
  // moved the marginal log-likelihood monotonically AWAY from the unpruned
  // value instead of towards it. The node-count guard that pruning was meant
  // to relieve lives in R (fastsimplexregmixed() rejects nAGQ^q > 1e5).
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
// Validate the CSR cluster offsets. A malformed 'starts' (empty, not starting
// at 0, not ending at N, or not strictly increasing) would otherwise index out
// of bounds; inside the OpenMP region a thrown exception aborts the process, so
// we reject it here, in serial code, as a clean R error.
inline void check_starts(const uvec& starts, const uword N) {
  if (starts.n_elem < 2 || starts[0] != 0 || starts[starts.n_elem - 1] != N) {
    Rcpp::stop("Invalid cluster offsets 'starts': must run 0, ..., nrow(data).");
  }
  for (uword j = 0; j + 1 < starts.n_elem; ++j) {
    if (starts[j + 1] <= starts[j]) {
      Rcpp::stop("Invalid cluster offsets 'starts': must be strictly increasing (no empty clusters).");
    }
  }
}

EvalResult mixed_core(
    const vec& theta, const vec& y, const mat& X, const mat& Z, const mat& W,
    const uvec& starts, const int q, const int mean_link,
    const mat& T, const vec& logW, const vec& t2,
    const int n_threads, const int inner_maxit, const double inner_tol,
    const bool need_grad, mat& Bhat,
    const vec& off_mu = vec(), const vec& off_phi = vec()) {

  check_starts(starts, y.n_elem);
  // Length 0 means "no offset"; anything else must match N.
  if ((off_mu.n_elem != 0 && off_mu.n_elem != y.n_elem) ||
      (off_phi.n_elem != 0 && off_phi.n_elem != y.n_elem)) {
    Rcpp::stop("Offset length does not match the number of observations.");
  }
  const bool use_off_mu = (off_mu.n_elem == y.n_elem);
  const bool use_off_phi = (off_phi.n_elem == y.n_elem);
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
  std::vector<int> sat_local(static_cast<std::size_t>(threads), 0);
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
    int local_sat = 0;
    vec* lg = need_grad ? &grad_local[static_cast<std::size_t>(tid)] : nullptr;

#ifdef _OPENMP
    #pragma omp for schedule(dynamic, 8)
#endif
    for (uword j = 0; j < J; ++j) {
      // `invalid` is a reduction variable, so this reads THIS thread's private
      // copy: it short-circuits the rest of this thread's clusters, not the
      // other threads'. That is intentional and the final reduction is still
      // correct; making it genuinely shared would need an atomic read and buy
      // nothing.
      if (invalid) continue;
      // Exception containment. The body below performs on the order of twenty
      // Armadillo allocations per cluster, every one of which can throw
      // std::bad_alloc under memory pressure -- and a throw crossing an OpenMP
      // structured block terminates the process rather than unwinding. Turning
      // any such failure into the existing `invalid` path converts a SIGABRT
      // into the clean "non-finite objective" R error. Table-based unwinding
      // costs nothing on the non-throwing path.
      try {
      const uword a = starts[j];
      const uword bb = starts[j + 1];
      const uword nj = bb - a;

      // Zero-copy views over the cluster's contiguous rows (no per-cluster
      // allocation/memcpy). Zj is kept as a materialised matrix because it is
      // consumed with .each_col().
      const auto Xj = X.rows(a, bb - 1);
      const mat  Zj = Z.rows(a, bb - 1);
      const auto Wj = W.rows(a, bb - 1);
      const auto yj = y.subvec(a, bb - 1);

      vec eta_mu_fixed = Xj * beta;
      vec eta_phi = Wj * gamma;
      if (use_off_mu)  eta_mu_fixed += off_mu.subvec(a, bb - 1);
      if (use_off_phi) eta_phi      += off_phi.subvec(a, bb - 1);
      vec phi(nj);
      for (uword i = 0; i < nj; ++i) phi[i] = safe_exp(eta_phi[i]);

      // Constants of the inner objective, hoisted out of hval().
      //
      // log(phi_i), log(y_i) and log(1-y_i) do not depend on b, yet hval()
      // recomputed all three on EVERY call -- and hval() is called twice for the
      // warm-start guard and then once plus up to thirty step-halvings per inner
      // Newton iteration, so a cluster paid hundreds of redundant log() calls.
      // A profile attributed 27% of a single-threaded mixed fit to libm's
      // log/exp, the single largest category. The node loop below already
      // hoisted exactly these into `cst`; the inner solver did not, which was an
      // inconsistency rather than a decision.
      //
      // The per-observation arithmetic is unchanged term for term, so the
      // objective is bit-identical -- only the transcendental calls move.
      vec log_phi(nj), log_yv(nj);
      for (uword i = 0; i < nj; ++i) {
        log_phi[i] = std::log(phi[i]);
        log_yv[i] = std::log(yj[i]) + std::log(1.0 - yj[i]);
      }

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
          const double lf = -0.5 * (LOG_2PI + log_phi[i])
                            - 1.5 * log_yv[i]
                            - 0.5 * dev / phi[i];
          if (!std::isfinite(lf)) { ok = false; return 0.0; }
          s += lf;
        }
        s -= 0.5 * arma::dot(bv, Sigma_inv * bv);
        return s;
      };

      // Warm-start guard. Bhat persists across objective evaluations, so the
      // mode stored here may have been produced at a theta that the OUTER line
      // search subsequently rejected, and can sit far from the mode at the
      // current theta. Keep it only when it actually beats the always-admissible
      // cold start b = 0. Without this the objective inherits a path dependence
      // on the rejected trial points, and a single wild trial can leave every
      // later evaluation starting from a useless mode -- enough to stall the
      // outer line search completely.
      {
        bool ok_warm = false;
        const double h_warm = hval(b, ok_warm);
        vec b_cold(static_cast<uword>(q), arma::fill::zeros);
        bool ok_cold = false;
        const double h_cold = hval(b_cold, ok_cold);
        if (!ok_warm || (ok_cold && h_cold > h_warm)) b = b_cold;
      }

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
          Iinfo[i] = simplex_fisher_eta_mu(mu, dmu, phi[i]);           // Fisher info
        }
        if (!ok) { ok_cluster = false; break; }

        const vec g = Zj.t() * s_mu - Sigma_inv * b;
        if (arma::abs(g).max() < inner_tol) break;

        mat Qf = Sigma_inv;
        Qf += Zj.t() * (Zj.each_col() % Iinfo);   // Z' diag(I) Z + Sigma_inv (SPD)
        vec delta;
        if (!arma::solve(delta, Qf, g, arma::solve_opts::likely_sympd)) { ok_cluster = false; break; }

        // Step-halving line search on h_j. Two safeguards matter here:
        // (1) the objective at the CURRENT b must itself be valid, otherwise
        //     h0 would silently be the sentinel 0.0 and the acceptance test
        //     would compare against a meaningless level;
        // (2) b is updated only when a step was actually accepted -- committing
        //     the last rejected trial would move the mode to a point the line
        //     search just refused.
        double alpha = 1.0;
        bool okb = false;
        const double h0 = hval(b, okb);
        if (!okb) { ok_cluster = false; break; }

        vec bnew;
        bool accepted = false;
        for (int hs = 0; hs < 30; ++hs) {
          const vec btry = b + alpha * delta;
          bool okn = false;
          const double h1 = hval(btry, okn);
          if (okn && h1 >= h0 - 1e-12) { bnew = btry; accepted = true; break; }
          alpha *= 0.5;
        }
        if (!accepted) break;   // already at the mode within this direction
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
          bool sat = false;
          if (!mean_deriv2_from_eta(eta[i], mean_link, mu, dmu, d2mu, &sat)) { ok = false; break; }
          if (sat) ++local_sat;
          const ObsKernel kk = simplex_obs_kernel(yj[i], mu, dmu, d2mu, phi[i]);
          if (!kk.ok) { ok = false; break; }
          w2[i] = kk.d2l_deta_mu2;
          Iinfo[i] = simplex_fisher_eta_mu(mu, dmu, phi[i]);
        }
        if (!ok) { invalid = 1; continue; }
        // symmatu() forces exact symmetry: Q is symmetric in exact arithmetic
        // but floating-point rounding leaves it slightly asymmetric, which makes
        // arma::chol warn and spuriously fail, triggering a needless fallback.
        Q = arma::symmatu(Sigma_inv - Zj.t() * (Zj.each_col() % w2));  // -Hessian (observed)
        // If the observed curvature is not PD, RIDGE it towards the (always
        // SPD) Fisher information instead of SWITCHING to it outright.
        //
        // Switching made the objective DISCONTINUOUS in theta: evaluated at the
        // same theta, the two choices of Q differed by 0.04 to 2.0 nats, and a
        // cluster could flip between them as the outer line search moved. A
        // jump like that violates what the line search assumes, and with the
        // Fisher Q the quadrature had still not converged at nAGQ = 21 in the
        // large-sigma cases where the observed Q converged by nAGQ = 5.
        //
        // The ridge is continuous in theta: lambda grows from zero only as far
        // as it must, so a cluster whose curvature is merely borderline gets a
        // Q that is arbitrarily close to the observed one rather than a
        // different matrix altogether. In a sweep of 1800 clusters (sigma up to
        // 4, mu down to 4e-11) the observed curvature was PD every time, so
        // this path is rare -- which is precisely why it must not distort the
        // objective when it does fire.
        mat Rchk;
        if (!arma::chol(Rchk, Q)) {
          const mat Qf = arma::symmatu(Sigma_inv + Zj.t() * (Zj.each_col() % Iinfo));
          double lambda = 1e-8;
          bool fixed = false;
          for (int t = 0; t < 40; ++t) {
            const mat Qr = arma::symmatu((1.0 - lambda) * Q + lambda * Qf);
            if (lambda < 1.0 && arma::chol(Rchk, Qr)) { Q = Qr; fixed = true; break; }
            lambda *= 4.0;
          }
          if (!fixed) Q = Qf;   // fully Fisher: the last resort, not the default
        }
      }

      mat R;
      if (!arma::chol(R, Q)) { invalid = 1; continue; }   // Q = R' R (upper R)
      double logdetQ = 0.0;
      for (int d = 0; d < q; ++d) logdetQ += 2.0 * std::log(R(d, d));
      // Non-throwing form. arma::inv()'s value-returning overload calls
      // arma_stop_runtime_error on failure, i.e. it THROWS -- and a throw
      // crossing an OpenMP structured block is undefined behaviour that
      // terminates the process on this toolchain (verified: SIGABRT, the outer
      // catch never runs). Every other decomposition in this region already
      // uses the bool form; this one did not. R came from a successful chol,
      // so failure is unlikely, but a denormal pivot that passes dpotrf can
      // still make dtrtri report a zero diagonal, and low probability times
      // process abort is not an acceptable trade.
      mat C;
      if (!arma::inv(C, arma::trimatu(R))) { invalid = 1; continue; }   // C C' = Q^{-1}

      // ---- Node-invariant precomputations, hoisted out of the node loop:
      //      the constant part of the log-density, 1/phi, 1/(y(1-y)), the fixed
      //      linear predictor eta_b = X beta + Z bhat, and the scaled node basis
      //      Zj * (sqrt2 C), so a node's eta is eta_b + (Zj sqrt2 C) t_k. ----
      const mat Cs = sqrt2 * C;                 // q x q
      const mat ZjCs = Zj * Cs;                 // nj x q
      const vec eta_b = eta_mu_fixed + Zj * b;   // nj
      vec inv_phi(nj), inv_yv(nj);
      double cst = 0.0;
      for (uword i = 0; i < nj; ++i) {
        inv_phi[i] = 1.0 / phi[i];
        const double one_y = 1.0 - yj[i];
        inv_yv[i] = 1.0 / (yj[i] * one_y);
        cst += -0.5 * (LOG_2PI + log_phi[i]) - 1.5 * log_yv[i];
      }

      // ---- Single AGHQ pass: build the unnormalized log-weights a_k and, when
      //      needed, cache the per-observation scores for a BLAS-based gradient
      //      (the kernel is thus evaluated once per node, not twice). ----
      vec avec(K);
      mat Bnodes(q, K);
      mat Smu, Sphi;
      if (need_grad) { Smu.set_size(nj, K); Sphi.set_size(nj, K); }
      bool ok_nodes = true;
      for (uword k = 0; k < K; ++k) {
        const vec tk = T.row(k).t();
        const vec bk = b + Cs * tk;
        Bnodes.col(k) = bk;
        const vec eta = eta_b + ZjCs * tk;
        double devsum = 0.0;
        for (uword i = 0; i < nj; ++i) {
          double mu, dmu;
          if (!mean_from_eta(eta[i], mean_link, mu, dmu)) { ok_nodes = false; break; }
          const double u = mu * (1.0 - mu);
          const double diff = yj[i] - mu;
          const double dev = diff * diff * inv_yv[i] / (u * u);
          devsum += dev * inv_phi[i];
          if (need_grad) {
            const double u3 = u * u * u;
            const double Pterm = diff * (mu * mu - 2.0 * mu * yj[i] + yj[i]);
            Smu(i, k) = Pterm * inv_yv[i] * inv_phi[i] / u3 * dmu;  // dl/deta_mu
            Sphi(i, k) = -0.5 + 0.5 * dev * inv_phi[i];            // dl/deta_phi
          }
        }
        if (!ok_nodes) break;
        const double quad = 0.5 * arma::dot(bk, Sigma_inv * bk);
        avec[k] = logW[k] + t2[k] + cst - 0.5 * devsum - quad;
      }
      if (!ok_nodes) { invalid = 1; continue; }

      const double lse = log_sum_exp(avec);
      const double loglik_j = half_q_log2 - 0.5 * logdetQ
                              - half_q_log2pi - 0.5 * logdetSigma + lse;
      if (!std::isfinite(loglik_j)) { invalid = 1; continue; }
      local_nll -= loglik_j;

      // ---- posterior-weighted analytic gradient: contract the cached scores
      //      with the normalized weights via BLAS. ----
      if (need_grad) {
        const vec Wjk = arma::exp(avec - lse);        // normalized weights, sum 1
        const vec gbeta = Xj.t() * (Smu * Wjk);
        const vec ggamma = Wj.t() * (Sphi * Wjk);
        vec gomega(m, arma::fill::zeros), og(m);
        for (uword k = 0; k < K; ++k) {
          const double wk = Wjk[k];
          if (wk <= 0.0) continue;
          omega_grad(Bnodes.col(k), Sigma_inv, D, q, og);
          gomega += wk * og;
        }
        // gradient of the NEGATIVE log-likelihood
        lg->subvec(0, p - 1) -= gbeta;
        lg->subvec(p, p + r - 1) -= ggamma;
        lg->subvec(p + r, dim - 1) -= gomega;
      }
      } catch (...) {
        invalid = 1;
      }
    }

    nll_local[static_cast<std::size_t>(tid)] = local_nll;
    sat_local[static_cast<std::size_t>(tid)] = local_sat;
  }

  if (invalid != 0) return fail();

  double nll = 0.0;
  for (const double v : nll_local) nll += v;
  int n_saturated = 0;
  for (const int v : sat_local) n_saturated += v;
  vec grad(dim, arma::fill::zeros);
  if (need_grad) for (const auto& g : grad_local) grad += g;

  return EvalResult{nll, std::move(grad), true, n_saturated};
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
    const double inner_tol = 1e-8,
    Rcpp::Nullable<Rcpp::NumericVector> off_mu_ = R_NilValue,
    Rcpp::Nullable<Rcpp::NumericVector> off_phi_ = R_NilValue) {

  const arma::vec off_mu = off_mu_.isNotNull()
    ? Rcpp::as<arma::vec>(off_mu_.get()) : arma::vec();
  const arma::vec off_phi = off_phi_.isNotNull()
    ? Rcpp::as<arma::vec>(off_phi_.get()) : arma::vec();

  mat T; vec logW, t2;
  simplex_fast::check_starts(starts, y.n_elem);
  simplex_fast::build_tensor(nAGQ, q, T, logW, t2);
  mat Bhat(starts.n_elem - 1, q, arma::fill::zeros);
  const auto res = simplex_fast::mixed_core(theta, y, X, Z, W, starts, q, mean_link,
                                            T, logW, t2, n_threads, inner_maxit,
                                            inner_tol, true, Bhat, off_mu, off_phi);
  return List::create(Named("value") = res.nll,
                      Named("gradient") = res.grad,
                      Named("valid") = res.valid,
                      Named("n_saturated") = res.n_saturated);
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
    const bool trace = false,
    Rcpp::Nullable<Rcpp::NumericVector> off_mu_ = R_NilValue,
    Rcpp::Nullable<Rcpp::NumericVector> off_phi_ = R_NilValue) {

  const arma::vec off_mu = off_mu_.isNotNull()
    ? Rcpp::as<arma::vec>(off_mu_.get()) : arma::vec();
  const arma::vec off_phi = off_phi_.isNotNull()
    ? Rcpp::as<arma::vec>(off_phi_.get()) : arma::vec();

  mat T; vec logW, t2;
  simplex_fast::check_starts(starts, y.n_elem);
  simplex_fast::build_tensor(nAGQ, q, T, logW, t2);
  mat Bhat(starts.n_elem - 1, q, arma::fill::zeros);

  auto objective = [&](const arma::vec& th) {
    return simplex_fast::mixed_core(th, y, X, Z, W, starts, q, mean_link,
                                    T, logW, t2, n_threads, inner_maxit,
                                    inner_tol, true, Bhat, off_mu, off_phi);
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
    const int inner_maxit = 50, const double inner_tol = 1e-8,
    Rcpp::Nullable<Rcpp::NumericVector> off_mu_ = R_NilValue,
    Rcpp::Nullable<Rcpp::NumericVector> off_phi_ = R_NilValue) {

  const arma::vec off_mu = off_mu_.isNotNull()
    ? Rcpp::as<arma::vec>(off_mu_.get()) : arma::vec();
  const arma::vec off_phi = off_phi_.isNotNull()
    ? Rcpp::as<arma::vec>(off_phi_.get()) : arma::vec();

  mat T; vec logW, t2;
  simplex_fast::check_starts(starts, y.n_elem);
  simplex_fast::build_tensor(nAGQ, q, T, logW, t2);
  const uword d = theta.n_elem;
  mat Bhat(starts.n_elem - 1, q, arma::fill::zeros);

  auto grad_at = [&](const arma::vec& th, bool& ok) {
    const auto res = simplex_fast::mixed_core(th, y, X, Z, W, starts, q, mean_link,
                                              T, logW, t2, n_threads, inner_maxit,
                                              inner_tol, true, Bhat, off_mu, off_phi);
    ok = res.valid;
    return res.grad;
  };

  mat H(d, d, arma::fill::zeros);
  for (uword jcol = 0; jcol < d; ++jcol) {
    // Serial loop of 2d objective evaluations; a meaningful share of a long
    // fit, and previously unreachable by Ctrl-C.
    Rcpp::checkUserInterrupt();
    double h = rel_step * std::max(1.0, std::abs(theta[jcol]));
    bool success = false;
    for (int attempt = 0; attempt < 12; ++attempt) {
      vec plus = theta, minus = theta;
      plus[jcol] += h; minus[jcol] -= h;
      bool okp = false, okm = false;
      const vec gp = grad_at(plus, okp);
      const vec gm = grad_at(minus, okm);
      // Guard on VALIDITY, not just finiteness: an invalid evaluation returns an
      // all-zero (finite) gradient, which would silently produce a zero Hessian
      // column and hence wrong standard errors.
      if (okp && okm && gp.is_finite() && gm.is_finite()) {
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
    const int inner_maxit = 50, const double inner_tol = 1e-8,
    Rcpp::Nullable<Rcpp::NumericVector> off_mu_ = R_NilValue,
    Rcpp::Nullable<Rcpp::NumericVector> off_phi_ = R_NilValue) {

  const arma::vec off_mu = off_mu_.isNotNull()
    ? Rcpp::as<arma::vec>(off_mu_.get()) : arma::vec();
  const arma::vec off_phi = off_phi_.isNotNull()
    ? Rcpp::as<arma::vec>(off_phi_.get()) : arma::vec();

  simplex_fast::check_starts(starts, y.n_elem);
  // One AGHQ node reproduces the mode-finding path; then read the modes back.
  mat T; vec logW, t2;
  simplex_fast::build_tensor(1, q, T, logW, t2);
  const uword J = starts.n_elem - 1;
  mat Bhat(J, q, arma::fill::zeros);
  const auto res = simplex_fast::mixed_core(theta, y, X, Z, W, starts, q, mean_link,
                                            T, logW, t2, n_threads, inner_maxit,
                                            inner_tol, false, Bhat, off_mu, off_phi);
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
    vec eta_phi = Wj * gamma;
    const vec bmode = Bhat.row(j).t();
    vec eta = Xj * beta + Zj * bmode;
    if (off_mu.n_elem == y.n_elem)  eta     += off_mu.subvec(a, bb - 1);
    if (off_phi.n_elem == y.n_elem) eta_phi += off_phi.subvec(a, bb - 1);
    // Both buffers are zero-filled: arma::vec(n) leaves its memory
    // uninitialised, and the Fisher fallback below consumes Iinfo even on the
    // path where the observed loop breaks early.
    vec w2(nj, arma::fill::zeros), Iinfo(nj, arma::fill::zeros);
    // The Fisher information depends only on (mu, dmu, phi) and is always
    // well defined, so it is built in its own pass and is never left partial.
    bool oklink = true;
    for (uword i = 0; i < nj; ++i) {
      const double phii = simplex_fast::safe_exp(eta_phi[i]);
      double mu, dmu, d2mu;
      if (!simplex_fast::mean_deriv2_from_eta(eta[i], mean_link, mu, dmu, d2mu)) {
        oklink = false;
        break;
      }
      Iinfo[i] = simplex_fast::simplex_fisher_eta_mu(mu, dmu, phii);
    }
    if (!oklink) {
      Rcpp::stop("The linear predictor is outside the valid domain of the selected mean link.");
    }
    bool okobs = true;
    for (uword i = 0; i < nj; ++i) {
      const double phii = simplex_fast::safe_exp(eta_phi[i]);
      double mu, dmu, d2mu;
      if (!simplex_fast::mean_deriv2_from_eta(eta[i], mean_link, mu, dmu, d2mu)) {
        okobs = false; break;
      }
      const simplex_fast::ObsKernel kk = simplex_fast::simplex_obs_kernel(yj[i], mu, dmu, d2mu, phii);
      if (!kk.ok) { okobs = false; break; }
      w2[i] = kk.d2l_deta_mu2;
    }
    mat Q = okobs ? arma::symmatu(Sigma_inv - Zj.t() * (Zj.each_col() % w2)) : Sigma_inv;
    const mat Qf = arma::symmatu(Sigma_inv + Zj.t() * (Zj.each_col() % Iinfo));
    mat Rchk;
    if (!okobs) {
      Q = Qf;
    } else if (!arma::chol(Rchk, Q)) {
      // Ridge towards Fisher rather than switching, matching mixed_core().
      double lambda = 1e-8;
      bool fixed = false;
      for (int t = 0; t < 40; ++t) {
        const mat Qr = arma::symmatu((1.0 - lambda) * Q + lambda * Qf);
        if (lambda < 1.0 && arma::chol(Rchk, Qr)) { Q = Qr; fixed = true; break; }
        lambda *= 4.0;
      }
      if (!fixed) Q = Qf;
    }
    mat Qinv;
    if (!arma::inv_sympd(Qinv, Q)) {
      Rcpp::stop("Posterior covariance of the random effects is not positive definite.");
    }
    postvar.slice(j) = Qinv;
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
    const int mean_link = 1, const bool include_re = true,
    Rcpp::Nullable<Rcpp::NumericVector> off_mu_ = R_NilValue,
    Rcpp::Nullable<Rcpp::NumericVector> off_phi_ = R_NilValue) {

  const arma::vec off_mu = off_mu_.isNotNull()
    ? Rcpp::as<arma::vec>(off_mu_.get()) : arma::vec();
  const arma::vec off_phi = off_phi_.isNotNull()
    ? Rcpp::as<arma::vec>(off_phi_.get()) : arma::vec();

  const uword p = X.n_cols, r = W.n_cols, N = X.n_rows;
  simplex_fast::check_starts(starts, N);
  if ((off_mu.n_elem != 0 && off_mu.n_elem != N) ||
      (off_phi.n_elem != 0 && off_phi.n_elem != N)) {
    Rcpp::stop("Offset length does not match the number of rows.");
  }
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
  vec eta_phi = W * gamma;
  if (off_mu.n_elem == N)  eta_mu  += off_mu;
  if (off_phi.n_elem == N) eta_phi += off_phi;

  vec mu(N), phi(N);
  for (uword i = 0; i < N; ++i) {
    double m_;
    // Unclamped (reporting path), matching simplex_predict_cpp().
    if (!simplex_fast::mean_from_eta_exact(eta_mu[i], mean_link, m_)) {
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
