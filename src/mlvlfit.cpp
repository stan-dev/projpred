#include <RcppArmadillo.h>
#include <cmath>

//[[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// Laplace approximation of binomial (logit) and Poisson (log) multilevel
// models with random intercepts, as in lme4::glmer() with `nAGQ = 1`:
//   eta = offset + X beta + Z Lambda u,  Lambda = diag(theta[factor]),
//   u ~ N(0, I),
//   D(theta, beta) = sum(dev.resids(y, mu, w)) + |u_hat|^2 + log det(H),
//   H = A' W A + I,  A = Z Lambda,  W = diag(w dmu/deta),
// where u_hat minimizes the first two terms (penalized IRLS).  Both links are
// canonical, so dD/deta through the deviance is -2 w (y - mu).  The gradient
// of D with respect to theta and beta is exact: u_hat is optimal, and the
// log-determinant term uses the implicit derivative of u_hat.
//
// Z has one indicator column per level of each grouping factor: `grp` is the
// n x k matrix of the 0-based level of observation i in factor j, `nlev` the
// number of levels per factor.  Level l of factor j is column off_j + l of Z.
// Cost per IRLS step O(n k^2 + q^3), q = sum(nlev).

// Family codes:
enum { ML_BINOMIAL = 0, ML_POISSON = 1 };

static vec ml_linkinv(const vec& eta, int fam) {
  if (fam == ML_POISSON) return exp(eta);
  vec mu(eta.n_elem);
  for (uword i = 0; i < eta.n_elem; i++) {
    double e = eta[i];
    mu[i] = (e >= 0) ? 1.0 / (1.0 + std::exp(-e)) :
      std::exp(e) / (1.0 + std::exp(e));
  }
  return mu;
}

// w dmu/deta:
static vec ml_wvar(const vec& mu, const vec& w, int fam) {
  if (fam == ML_POISSON) return w % mu;
  return w % mu % (1 - mu);
}

// Terms of the deviance residuals without eta, per observation (see ml_dev()):
static vec ml_dev_y(const vec& y, int fam) {
  vec c(y.n_elem);
  for (uword i = 0; i < y.n_elem; i++) {
    double yi = y[i];
    double t = (yi > 0) ? yi * std::log(yi) : 0.0;
    if (fam == ML_POISSON) {
      t -= yi;
    } else if (yi < 1) {
      t += (1 - yi) * std::log1p(-yi);
    }
    c[i] = t;
  }
  return c;
}

// The sum of the deviance residuals is computed from eta, so that it stays
// finite and accurate where mu rounds to 0 or 1:
//   binomial: y log(y) + (1 - y) log(1 - y) + softplus(eta) - y eta,
//   Poisson:  y log(y) - y + exp(eta) - y eta.
// `dev_y` holds the terms without eta (ml_dev_y()).  Each observation adds
// them before the sum, so that rounding stays relative to the deviance near
// the optimum.  softplus(eta) - y eta is (1 - y) eta + log1p(exp(-eta)) for
// eta > 0.  A product with y = 0 or 1 - y = 0 is skipped, so that an infinite
// eta gives the limit 0.
static double ml_dev(const vec& y, const vec& eta, const vec& w, int fam,
                     const vec& dev_y) {
  double d = 0;
  for (uword i = 0; i < y.n_elem; i++) {
    if (w[i] == 0) continue;
    double yi = y[i], e = eta[i];
    double t;
    if (fam == ML_POISSON) {
      t = std::exp(e) - ((yi > 0) ? yi * e : 0.0);
    } else if (e > 0) {
      t = ((yi < 1) ? (1 - yi) * e : 0.0) + std::log1p(std::exp(-e));
    } else {
      t = std::log1p(std::exp(e)) - ((yi > 0) ? yi * e : 0.0);
    }
    d += 2 * w[i] * (dev_y[i] + t);
  }
  return d;
}

// Objective D and (for `want_grad`) its gradient with respect to (theta,
// beta), at the conditional modes u_hat found from the start `u`.
// [[Rcpp::export]]
List mlvl_laplace_c(const arma::vec& theta, const arma::vec& beta,
                    const arma::mat& x, const IntegerMatrix& grp,
                    const IntegerVector& nlev, const arma::vec& y,
                    const arma::vec& w, const arma::vec& offset, arma::vec u,
                    int fam, bool want_grad, double tol = 1e-10,
                    int maxit = 50) {
  const uword n = x.n_rows, k = theta.n_elem;
  uvec offs(k);
  uword q = 0;
  for (uword j = 0; j < k; j++) {
    offs[j] = q;
    q += nlev[j];
  }
  vec lam(q);
  for (uword j = 0; j < k; j++) {
    lam.subvec(offs[j], offs[j] + nlev[j] - 1).fill(theta[j]);
  }
  // Column of Z for observation i in factor j:
  umat col(n, k);
  for (uword j = 0; j < k; j++) {
    for (uword i = 0; i < n; i++) col(i, j) = offs[j] + grp(i, j);
  }
  vec eta0 = offset + x * beta;
  // A u:
  auto a_u = [&](const vec& uu) -> vec {
    vec out(n, fill::zeros);
    for (uword j = 0; j < k; j++) {
      for (uword i = 0; i < n; i++) out[i] += lam[col(i, j)] * uu[col(i, j)];
    }
    return out;
  };
  // A' v:
  auto at_v = [&](const vec& v) -> vec {
    vec out(q, fill::zeros);
    for (uword j = 0; j < k; j++) {
      for (uword i = 0; i < n; i++) out[col(i, j)] += v[i];
    }
    return out % lam;
  };
  // Z' W Z:
  auto zwz = [&](const vec& wd) -> mat {
    mat s(q, q, fill::zeros);
    for (uword i = 0; i < n; i++) {
      for (uword j = 0; j < k; j++) {
        for (uword jj = 0; jj < k; jj++) s(col(i, j), col(i, jj)) += wd[i];
      }
    }
    return s;
  };
  const vec dev_y = ml_dev_y(y, fam);
  auto obj = [&](const vec& uu) -> double {
    return ml_dev(y, eta0 + a_u(uu), w, fam, dev_y) + dot(uu, uu);
  };
  double f = obj(u);
  if (!std::isfinite(f)) {
    u.zeros();
    f = obj(u);
  }
  if (!std::isfinite(f)) {
    return List::create(Named("D") = R_PosInf, Named("u") = u);
  }
  int n_irls = 0;
  mat h(q, q);
  for (int it = 0; it < maxit; it++) {
    n_irls++;
    vec mu = ml_linkinv(eta0 + a_u(u), fam);
    vec wd = ml_wvar(mu, w, fam);
    vec r = w % (y - mu);
    h = zwz(wd) % (lam * lam.t());
    h.diag() += 1;
    vec step;
    if (!solve(step, h, at_v(r) - u, solve_opts::likely_sympd)) {
      return List::create(Named("D") = R_PosInf, Named("u") = u);
    }
    // Step halving; without a step that does not increase the objective, keep
    // u:
    double s = 1, f_new;
    vec u_new;
    bool accepted = false;
    while (s >= 1e-10) {
      u_new = u + s * step;
      f_new = obj(u_new);
      if (std::isfinite(f_new) && f_new <= f + 1e-12 * std::abs(f)) {
        accepted = true;
        break;
      }
      s /= 2;
    }
    if (!accepted) break;
    bool conv = accu(abs(u_new - u)) < tol * (1 + accu(abs(u)));
    u = u_new;
    f = f_new;
    if (conv) break;
  }
  vec eta = eta0 + a_u(u);
  vec mu = ml_linkinv(eta, fam);
  vec wd = ml_wvar(mu, w, fam);
  vec r = w % (y - mu);
  mat s_zwz = zwz(wd);
  h = s_zwz % (lam * lam.t());
  h.diag() += 1;
  mat rch;
  if (!chol(rch, h)) {
    return List::create(Named("D") = R_PosInf, Named("u") = u);
  }
  double d_val = f + 2 * accu(log(rch.diag()));
  if (!want_grad) {
    return List::create(Named("D") = d_val, Named("u") = u,
                        Named("n_irls") = n_irls);
  }
  mat h_inv = inv_sympd(h);
  // dW/deta:
  vec w1 = (fam == ML_POISSON) ? wd : vec(wd % (1 - 2 * mu));
  // diag(A H^-1 A'):
  vec cdiag(n, fill::zeros);
  for (uword i = 0; i < n; i++) {
    for (uword j = 0; j < k; j++) {
      for (uword jj = 0; jj < k; jj++) {
        cdiag[i] += lam[col(i, j)] * lam[col(i, jj)] *
          h_inv(col(i, j), col(i, jj));
      }
    }
  }
  vec v = cdiag % w1;
  // beta:
  vec s_v = h_inv * at_v(v);
  vec g_beta = -2 * x.t() * r + x.t() * v - x.t() * (wd % a_u(s_v));
  // theta: diagonal of H^-1 Lambda Z'WZ:
  mat ls = s_zwz.each_col() % lam;
  vec d_hm = sum(h_inv % ls.t(), 1);
  vec g_theta(k);
  for (uword j = 0; j < k; j++) {
    uword a = offs[j], b = offs[j] + nlev[j] - 1;
    vec zu(n);
    for (uword i = 0; i < n; i++) zu[i] = u[col(i, j)];
    double part1 = -2 * dot(r, zu);
    double part2 = 2 * accu(d_hm.subvec(a, b));
    // Derivative of u_hat with respect to theta_j:
    vec rhs = -at_v(wd % zu);
    vec zr(q, fill::zeros);
    for (uword i = 0; i < n; i++) zr[col(i, j)] += r[i];
    rhs.subvec(a, b) += zr.subvec(a, b);
    vec du = h_inv * rhs;
    vec deta = zu + a_u(du);
    g_theta[j] = part1 + part2 + dot(v, deta);
  }
  return List::create(Named("D") = d_val, Named("u") = u,
                      Named("grad") = join_cols(g_theta, g_beta),
                      Named("n_irls") = n_irls);
}
