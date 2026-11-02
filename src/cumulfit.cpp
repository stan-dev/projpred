#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

//[[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// Cumulative-link model fits for the augmented-data projection, with the
// objective of MASS::polr(): for observation i with linear predictor eta_i and
// category weights w_ic,
//   L(beta, zeta) = -sum_i sum_c w_ic log(F(z1_ic) - F(z2_ic)),
//   z1_ic = min(100, zeta_c - eta_i), z2_ic = max(-100, zeta_{c-1} - eta_i),
// with zeta_0 = -Inf and zeta_C = Inf. The objective, with its bounds of +-100,
// and the latent distributions transcribe the negative log-likelihood and the
// `pfun` / `dfun` of MASS::polr() (R package MASS, GPL-2 | GPL-3). The fit is
// new: polr() minimizes with optim() (BFGS), and this code uses Newton steps
// with the analytic gradient and Hessian, and a backtracking line search that
// keeps the thresholds increasing.

// Link codes (as polr()'s `method`):
enum { CL_LOGISTIC = 1, CL_PROBIT = 2, CL_LOGLOG = 3, CL_CLOGLOG = 4,
       CL_CAUCHIT = 5 };

// Distribution function F, survival function S = 1 - F (computed directly, so
// that upper-tail probabilities do not cancel), density f and derivative of the
// density f' of the latent variable, as in polr()'s `pfun` and `dfun`:
static inline void cl_funs(int link, double z, double& F, double& S,
                           double& f, double& fd) {
  switch (link) {
  case CL_LOGISTIC: {
    // One exp() for all four values (e <= 1, so nothing overflows):
    double e = std::exp(-std::abs(z));
    double d = 1.0 / (1.0 + e);
    double tail = e * d;
    F = (z >= 0.0) ? d : tail;
    S = (z >= 0.0) ? tail : d;
    f = tail * d;
    fd = f * (S - F);
    break;
  }
  case CL_PROBIT: {
    F = R::pnorm(z, 0.0, 1.0, 1, 0);
    S = R::pnorm(z, 0.0, 1.0, 0, 0);
    f = R::dnorm(z, 0.0, 1.0, 0);
    fd = -z * f;
    break;
  }
  case CL_LOGLOG: {
    // MASS::pgumbel(): exp(-exp(-z))
    double e = std::exp(-z);
    F = std::exp(-e);
    S = -std::expm1(-e);
    f = std::exp(-z - e);
    fd = f * (e - 1.0);
    break;
  }
  case CL_CLOGLOG: {
    // MASS:::pGumbel(): 1 - exp(-exp(z))
    double e = std::exp(z);
    F = -std::expm1(-e);
    S = std::exp(-e);
    f = std::exp(z - e);
    fd = f * (1.0 - e);
    break;
  }
  case CL_CAUCHIT: {
    F = R::pcauchy(z, 0.0, 1.0, 1, 0);
    S = R::pcauchy(z, 0.0, 1.0, 0, 0);
    f = R::dcauchy(z, 0.0, 1.0, 0);
    fd = -2.0 * z * M_PI * f * f;
    break;
  }
  default:
    Rcpp::stop("Unknown link code.");
  }
}

static inline double cl_quantile(int link, double p) {
  switch (link) {
  case CL_LOGISTIC: return R::qlogis(p, 0.0, 1.0, 1, 0);
  case CL_PROBIT: return R::qnorm(p, 0.0, 1.0, 1, 0);
  case CL_LOGLOG: return -std::log(-std::log(p));
  case CL_CLOGLOG: return std::log(-std::log1p(-p));
  case CL_CAUCHIT: return R::qcauchy(p, 0.0, 1.0, 1, 0);
  default: Rcpp::stop("Unknown link code.");
  }
  return 0.0;
}

// Objective (and optionally gradient and Hessian) at (beta, zeta). Returns
// Inf if a category probability with positive weight is not positive.
static double cl_eval(const mat& x, const mat& w, const vec& offset,
                      const vec& beta, const vec& zeta, int link,
                      bool deriv, vec& grad, mat& hess) {
  const uword n = w.n_rows;
  const uword ncat = w.n_cols;
  const uword q = ncat - 1;
  const uword p = x.n_cols;
  vec eta = offset;
  if (p > 0) eta += x * beta;
  // std::min() and std::max() would turn NaN arguments into the bounds +-100:
  if (!eta.is_finite()) return std::numeric_limits<double>::infinity();
  vec g_eta, h_etaeta;
  mat h_etazeta;
  vec g_zeta;
  mat h_zetazeta;
  if (deriv) {
    g_eta.zeros(n);
    h_etaeta.zeros(n);
    h_etazeta.zeros(n, q);
    g_zeta.zeros(q);
    h_zetazeta.zeros(q, q);
  }
  double loss = 0.0;
  // The functions at zeta_j - eta_i are computed once per threshold and
  // observation (the cache key is z, so values clamped at +-100 are computed
  // again).
  double F_hi, S_hi, f_hi, fd_hi, F_lo, S_lo, f_lo, fd_lo;
  cl_funs(link, 100.0, F_hi, S_hi, f_hi, fd_hi);
  cl_funs(link, -100.0, F_lo, S_lo, f_lo, fd_lo);
  std::vector<double> t_z(q), t_F(q), t_S(q), t_f(q), t_fd(q);
  std::vector<char> t_done(q);
  auto thr_funs = [&](uword j, double z, double& F, double& S, double& f,
                      double& fd) {
    if (!(t_done[j] && t_z[j] == z)) {
      cl_funs(link, z, t_F[j], t_S[j], t_f[j], t_fd[j]);
      t_z[j] = z;
      t_done[j] = 1;
    }
    F = t_F[j];
    S = t_S[j];
    f = t_f[j];
    fd = t_fd[j];
  };
  for (uword i = 0; i < n; ++i) {
    std::fill(t_done.begin(), t_done.end(), 0);
    for (uword c = 0; c < ncat; ++c) {
      double wic = w(i, c);
      if (wic == 0.0) continue;
      // Upper and lower threshold of category c (0-based): zeta_c and
      // zeta_{c-1}:
      double z1 = (c < q) ? std::min(100.0, zeta(c) - eta(i)) : 100.0;
      double z2 = (c > 0) ? std::max(-100.0, zeta(c - 1) - eta(i)) : -100.0;
      double F1, S1, f1, fd1, F2, S2, f2, fd2;
      if (c < q) {
        thr_funs(c, z1, F1, S1, f1, fd1);
      } else {
        F1 = F_hi;
        S1 = S_hi;
        f1 = f_hi;
        fd1 = fd_hi;
      }
      if (c > 0) {
        thr_funs(c - 1, z2, F2, S2, f2, fd2);
      } else {
        F2 = F_lo;
        S2 = S_lo;
        f2 = f_lo;
        fd2 = fd_lo;
      }
      if (c == q) {
        f1 = 0.0;
        fd1 = 0.0;
      }
      if (c == 0) {
        f2 = 0.0;
        fd2 = 0.0;
      }
      // F(z1) - F(z2), from the upper tails if both arguments are positive:
      double pr = (z2 > 0.0) ? (S2 - S1) : (F1 - F2);
      if (!(pr > 0.0)) return std::numeric_limits<double>::infinity();
      loss -= wic * std::log(pr);
      if (!deriv) continue;
      double a = f1 / pr, b = f2 / pr;
      double ad = fd1 / pr, bd = fd2 / pr;
      // Derivatives of -w * log(pr):
      g_eta(i) -= wic * (-(a - b));
      h_etaeta(i) -= wic * ((ad - bd) - (a - b) * (a - b));
      if (c < q) {
        g_zeta(c) -= wic * a;
        h_zetazeta(c, c) -= wic * (ad - a * a);
        h_etazeta(i, c) -= wic * (-ad + a * (a - b));
      }
      if (c > 0) {
        g_zeta(c - 1) -= wic * (-b);
        h_zetazeta(c - 1, c - 1) -= wic * (-bd - b * b);
        h_etazeta(i, c - 1) -= wic * (bd - b * (a - b));
      }
      if (c > 0 && c < q) {
        h_zetazeta(c, c - 1) -= wic * (a * b);
        h_zetazeta(c - 1, c) -= wic * (a * b);
      }
    }
  }
  if (deriv) {
    grad.set_size(p + q);
    hess.set_size(p + q, p + q);
    if (p > 0) {
      grad.head(p) = x.t() * g_eta;
      hess.submat(0, 0, p - 1, p - 1) = x.t() * (x.each_col() % h_etaeta);
      mat hbz = x.t() * h_etazeta;
      hess.submat(0, p, p - 1, p + q - 1) = hbz;
      hess.submat(p, 0, p + q - 1, p - 1) = hbz.t();
    }
    grad.tail(q) = g_zeta;
    hess.submat(p, p, p + q - 1, p + q - 1) = h_zetazeta;
  }
  return loss;
}

// Fits one weight matrix `w` (n x C). Writes the estimates into `beta` and
// `zeta`; returns the objective. `conv` is 0 on convergence.
static double cl_fit_one(const mat& x, const mat& w, const vec& offset,
                         int link, double thresh, int maxit, vec& beta,
                         vec& zeta, int& conv, int& niter) {
  const uword p = x.n_cols;
  const uword ncat = w.n_cols;
  const uword q = ncat - 1;
  // Start: beta = 0 and the thresholds of the model without predictors:
  beta.zeros(p);
  zeta.set_size(q);
  rowvec wc = sum(w, 0);
  double wtot = accu(wc);
  double cum = 0.0;
  for (uword c = 0; c < q; ++c) {
    cum += wc(c);
    double pc = std::min(std::max(cum / wtot, 1e-8), 1.0 - 1e-8);
    // Within +-50, so that each category has a positive probability at the
    // start (heavy-tailed links such as cauchit have quantiles beyond the
    // bounds +-100 for small weight shares):
    zeta(c) = std::min(50.0, std::max(-50.0, cl_quantile(link, pc)));
    if (c > 0 && zeta(c) <= zeta(c - 1)) zeta(c) = zeta(c - 1) + 1e-6;
  }
  vec grad;
  mat hess;
  double loss = cl_eval(x, w, offset, beta, zeta, link, true, grad, hess);
  conv = 1;
  niter = 0;
  // Infinite objective at the start values (no derivatives): not converged.
  if (!std::isfinite(loss)) return loss;
  vec theta(p + q), dir, theta_new;
  for (int it = 0; it < maxit; ++it) {
    niter = it + 1;
    // Newton direction (damped if the Hessian is not positive definite):
    mat hreg = hess;
    double damp = 0.0;
    bool ok = false;
    double dscale = std::max(1e-10, max(abs(hess.diag())));
    for (int k = 0; k < 30 && !ok; ++k) {
      mat r;
      if (damp > 0.0) hreg = hess + damp * eye(p + q, p + q);
      ok = chol(r, hreg);
      if (ok) {
        dir = -solve(trimatu(r), solve(trimatl(r.t()), grad));
      } else {
        damp = (damp == 0.0) ? 1e-8 * dscale : 10.0 * damp;
      }
    }
    if (!ok) break;
    double decrement = -dot(grad, dir);
    if (decrement < thresh * std::max(1.0, std::abs(loss))) {
      conv = 0;
      break;
    }
    if (p > 0) theta.head(p) = beta;
    theta.tail(q) = zeta;
    // Backtracking line search (Armijo), keeping the thresholds increasing.
    // The full-step trial (t = 1, usually accepted) also computes the
    // gradient and the Hessian, so that an accepted full step needs no second
    // evaluation at the same point.
    double t = 1.0;
    bool accepted = false;
    bool accepted_deriv = false;
    double loss_accepted = 0.0;
    vec grad_new;
    mat hess_new;
    for (int ls = 0; ls < 60; ++ls) {
      theta_new = theta + t * dir;
      vec zeta_new = theta_new.tail(q);
      bool ordered = true;
      for (uword c = 1; c < q; ++c) {
        if (!(zeta_new(c) > zeta_new(c - 1))) {
          ordered = false;
          break;
        }
      }
      if (ordered) {
        vec beta_new = (p > 0) ? vec(theta_new.head(p)) : vec();
        bool with_deriv = (ls == 0);
        double loss_new = cl_eval(x, w, offset, beta_new, zeta_new, link,
                                  with_deriv, grad_new, hess_new);
        if (std::isfinite(loss_new) &&
            loss_new <= loss - 1e-4 * t * decrement) {
          accepted = true;
          accepted_deriv = with_deriv;
          loss_accepted = loss_new;
          break;
        }
      }
      t *= 0.5;
    }
    if (!accepted) {
      // No decrease possible along the direction: at the optimum up to
      // rounding if the decrement is tiny, otherwise not converged.
      if (decrement < 1e3 * thresh * std::max(1.0, std::abs(loss))) conv = 0;
      break;
    }
    if (p > 0) beta = theta_new.head(p);
    zeta = theta_new.tail(q);
    if (accepted_deriv) {
      loss = loss_accepted;
      grad.swap(grad_new);
      hess.swap(hess_new);
    } else {
      loss = cl_eval(x, w, offset, beta, zeta, link, true, grad, hess);
    }
  }
  return loss;
}

// Fits the cumulative-link model for S weight matrices, given as the n x C
// column blocks of `w_all` (n x (C * S)). The model matrix `x` (n x p, without
// intercept) and the offsets are shared.
// [[Rcpp::export]]
List cumul_fit_multi_c(arma::mat x, arma::mat w_all, int ncat,
                       arma::vec offset, int link, double thresh,
                       int maxit) {
  const uword n = x.n_rows;
  const uword p = x.n_cols;
  const uword q = ncat - 1;
  const uword S = w_all.n_cols / ncat;
  if (w_all.n_rows != n || S * ncat != w_all.n_cols) {
    Rcpp::stop("Dimensions of `w_all` do not match.");
  }
  mat beta_out(p, S);
  mat zeta_out(q, S);
  vec loss_out(S);
  IntegerVector conv_out(S), niter_out(S);
  if (!x.is_finite() || !offset.is_finite()) {
    // Not fitted (the linear predictor would not be finite):
    beta_out.zeros();
    zeta_out.fill(NA_REAL);
    loss_out.fill(std::numeric_limits<double>::infinity());
    std::fill(conv_out.begin(), conv_out.end(), 1);
    return List::create(Named("beta") = beta_out, Named("zeta") = zeta_out,
                        Named("loss") = loss_out,
                        Named("convergence") = conv_out,
                        Named("niter") = niter_out);
  }
  for (uword s = 0; s < S; ++s) {
    mat w = w_all.cols(s * ncat, (s + 1) * ncat - 1);
    vec beta, zeta;
    int conv, niter;
    loss_out(s) = cl_fit_one(x, w, offset, link, thresh, maxit, beta, zeta,
                             conv, niter);
    if (p > 0) beta_out.col(s) = beta;
    zeta_out.col(s) = zeta;
    conv_out[s] = conv;
    niter_out[s] = niter;
  }
  return List::create(Named("beta") = beta_out, Named("zeta") = zeta_out,
                      Named("loss") = loss_out, Named("convergence") = conv_out,
                      Named("niter") = niter_out);
}
