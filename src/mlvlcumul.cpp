#include <RcppArmadillo.h>
#include <cmath>
#include <vector>

//[[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// Laplace fit of a cumulative-link submodel with one group-level intercept
// term `(1 | g)` (augmented-data projection). Model, in the parameterization
// of ordinal::clmm():
//   P(y <= j) = F(alpha_j - x' beta - tau u_g),  u_g ~ N(0, 1),
// with alpha_0 = -Inf and alpha_K = Inf. For row i with category y_i and
// weight w_i,
//   pi_i = F(e1_i - s_i) - F(e2_i - s_i),
//   e1_i = alpha_{y_i} - x_i' beta, e2_i = alpha_{y_i - 1} - x_i' beta,
//   s_i = tau u_{g(i)}.
// Objective (minus the Laplace approximation of the log-likelihood):
//   L(alpha, beta, tau) = -sum_i w_i log pi_i + |u_hat|^2 / 2 +
//                         sum_g log(D_g) / 2,
//   D_g = 1 + tau^2 sum_{i in g} h_i,
// where h_i is the second derivative of -w_i log pi_i with respect to s_i, and
// u_hat minimizes the first two terms. There is one grouping variable, so the
// Hessian with respect to u is diagonal, and the Newton iteration for u_hat
// works per group. The gradient with respect to (alpha, beta, tau) is exact:
// it includes the derivative of u_hat in the log-determinant term.

// Link codes (as in `cumulfit.cpp`):
enum { MC_LOGISTIC = 1, MC_PROBIT = 2, MC_LOGLOG = 3, MC_CLOGLOG = 4,
       MC_CAUCHIT = 5 };

// Distribution function F, survival function S = 1 - F and the first three
// derivatives of F (f, g, k) of the latent variable:
static inline void mc_funs(int link, double z, double& F, double& S, double& f,
                           double& g, double& k) {
  if (z == R_PosInf) {
    F = 1; S = 0; f = 0; g = 0; k = 0;
    return;
  }
  if (z == R_NegInf) {
    F = 0; S = 1; f = 0; g = 0; k = 0;
    return;
  }
  switch (link) {
  case MC_LOGISTIC: {
    F = R::plogis(z, 0.0, 1.0, 1, 0);
    S = R::plogis(z, 0.0, 1.0, 0, 0);
    f = F * S;
    g = f * (S - F);
    k = f * (1 - 6 * f);
    break;
  }
  case MC_PROBIT: {
    F = R::pnorm(z, 0.0, 1.0, 1, 0);
    S = R::pnorm(z, 0.0, 1.0, 0, 0);
    f = R::dnorm(z, 0.0, 1.0, 0);
    g = -z * f;
    k = (z * z - 1) * f;
    break;
  }
  case MC_LOGLOG: {
    double e = std::exp(-z);
    F = std::exp(-e);
    S = -std::expm1(-e);
    f = std::exp(-z - e);
    // In the far tail, f is 0 and e is infinite:
    g = (f > 0) ? f * (e - 1) : 0.0;
    k = (f > 0) ? f * ((e - 1) * (e - 1) - e) : 0.0;
    break;
  }
  case MC_CLOGLOG: {
    double e = std::exp(z);
    F = -std::expm1(-e);
    S = std::exp(-e);
    f = std::exp(z - e);
    g = (f > 0) ? f * (1 - e) : 0.0;
    k = (f > 0) ? f * ((1 - e) * (1 - e) - e) : 0.0;
    break;
  }
  case MC_CAUCHIT: {
    double d = 1 + z * z;
    F = R::pcauchy(z, 0.0, 1.0, 1, 0);
    S = R::pcauchy(z, 0.0, 1.0, 0, 0);
    f = 1 / (M_PI * d);
    g = -2 * z / (M_PI * d * d);
    k = (6 * z * z - 2) / (M_PI * d * d * d);
    break;
  }
  default:
    Rcpp::stop("Unknown link code.");
  }
}

// -w log pi of one row and its derivatives: m, h with respect to the shift s;
// l1, l2 with respect to e1, e2; m1, m2, h1, h2: the derivatives of m and h
// with respect to e1, e2.
struct McParts {
  double l, m, h, l1, l2, m1, m2, h1, h2;
  bool finite() const {
    return std::isfinite(l + m + h + l1 + l2 + m1 + m2 + h1 + h2);
  }
};

// F(z1) - F(z2), from the upper tails if both arguments are positive:
static inline double mc_prob(double F1, double S1, double F2, double S2,
                             double z2) {
  return (z2 > 0) ? (S2 - S1) : (F1 - F2);
}

static inline McParts mc_parts(int link, double z1, double z2, double w) {
  McParts r;
  double F1, S1, f1, g1, k1, F2, S2, f2, g2, k2;
  mc_funs(link, z1, F1, S1, f1, g1, k1);
  mc_funs(link, z2, F2, S2, f2, g2, k2);
  double ip = 1 / mc_prob(F1, S1, F2, S2, z2);
  double ip2 = ip * ip, ip3 = ip2 * ip;
  double a = f1 - f2, c = g1 - g2;
  r.l = w * std::log(ip);
  r.m = w * a * ip;
  r.h = w * (a * a * ip2 - c * ip);
  r.l1 = -w * f1 * ip;
  r.l2 = w * f2 * ip;
  r.m1 = w * (g1 * ip - a * f1 * ip2);
  r.m2 = w * (-g2 * ip + a * f2 * ip2);
  r.h1 = w * (2 * a * g1 * ip2 - 2 * a * a * f1 * ip3 - k1 * ip + c * f1 * ip2);
  r.h2 = w * (-2 * a * g2 * ip2 + 2 * a * a * f2 * ip3 + k2 * ip - c * f2 * ip2);
  return r;
}

// Objective L and (for `want_grad`) its gradient with respect to (alpha, beta,
// tau), at the conditional modes u_hat found from the start `u`. `grp` holds
// 0-based level codes, `y` 1-based category codes. Rows with weight 0 are
// ignored. L is `Inf` if the thresholds do not increase, or if a category
// probability with positive weight is not positive or so small that the
// derivatives are not finite. `gmax` is the largest absolute element of the
// gradient with respect to u at the returned `u` (0 at the conditional
// modes).
// [[Rcpp::export]]
List mlvl_cumul_c(const arma::vec& alpha, const arma::vec& beta, double tau,
                  const arma::mat& x, const IntegerVector& grp, int nlev,
                  const IntegerVector& y, const arma::vec& w, arma::vec u,
                  int link, bool want_grad, double tol = 1e-10,
                  int maxit = 50) {
  const uword n = x.n_rows;
  const int ncat = alpha.n_elem + 1;
  List out_inf = List::create(Named("L") = R_PosInf, Named("u") = u);
  for (uword j = 1; j < alpha.n_elem; j++) {
    if (!(alpha[j] > alpha[j - 1])) return out_inf;
  }
  vec xb = (beta.n_elem > 0) ? vec(x * beta) : vec(n, fill::zeros);
  if (!xb.is_finite() || !std::isfinite(tau)) return out_inf;
  vec e1(n), e2(n);
  for (uword i = 0; i < n; i++) {
    int yi = y[i];
    e1[i] = ((yi < ncat) ? alpha[yi - 1] : R_PosInf) - xb[i];
    e2[i] = ((yi > 1) ? alpha[yi - 2] : R_NegInf) - xb[i];
  }
  // The first two terms of L at `uu`:
  auto fval = [&](const vec& uu) -> double {
    double f = 0.5 * dot(uu, uu);
    for (uword i = 0; i < n; i++) {
      if (w[i] == 0) continue;
      double s = tau * uu[grp[i]];
      double F1, S1, f1, g1, k1, F2, S2, f2, g2, k2;
      mc_funs(link, e1[i] - s, F1, S1, f1, g1, k1);
      mc_funs(link, e2[i] - s, F2, S2, f2, g2, k2);
      double pr = mc_prob(F1, S1, F2, S2, e2[i] - s);
      if (!(pr > 0)) return R_PosInf;
      f -= w[i] * std::log(pr);
    }
    return f;
  };
  double f = fval(u);
  if (!std::isfinite(f)) {
    u.zeros();
    f = fval(u);
    if (!std::isfinite(f)) return out_inf;
  }
  vec Mg(nlev), Hg(nlev);
  for (int it = 0; it < maxit; it++) {
    Mg.zeros();
    Hg.zeros();
    for (uword i = 0; i < n; i++) {
      if (w[i] == 0) continue;
      double s = tau * u[grp[i]];
      McParts r = mc_parts(link, e1[i] - s, e2[i] - s, w[i]);
      Mg[grp[i]] += r.m;
      Hg[grp[i]] += r.h;
    }
    if (!Mg.is_finite() || !Hg.is_finite()) return out_inf;
    vec grad_u = tau * Mg + u;
    if (max(abs(grad_u)) < tol) break;
    // Newton step per group; a gradient step where the Hessian is not
    // positive:
    vec dg = 1 + tau * tau * Hg;
    dg.elem(find(dg <= 0)).ones();
    vec step = grad_u / dg;
    double st = 1, f_new = f;
    vec u_new = u;
    bool accepted = false;
    while (st >= 1e-10) {
      u_new = u - st * step;
      f_new = fval(u_new);
      if (std::isfinite(f_new) && f_new <= f + 1e-12 * std::abs(f)) {
        accepted = true;
        break;
      }
      st /= 2;
    }
    if (!accepted) break;
    u = u_new;
    f = f_new;
  }
  // The parts at u_hat:
  vec HSg(nlev, fill::zeros);
  Mg.zeros();
  Hg.zeros();
  std::vector<McParts> parts(n);
  double lsum = 0;
  for (uword i = 0; i < n; i++) {
    if (w[i] == 0) continue;
    double s = tau * u[grp[i]];
    parts[i] = mc_parts(link, e1[i] - s, e2[i] - s, w[i]);
    if (!parts[i].finite()) return out_inf;
    lsum += parts[i].l;
    Mg[grp[i]] += parts[i].m;
    Hg[grp[i]] += parts[i].h;
    HSg[grp[i]] -= parts[i].h1 + parts[i].h2;
  }
  vec Dg = 1 + tau * tau * Hg;
  if (!Mg.is_finite() || !Dg.is_finite() || !(Dg.min() > 0)) return out_inf;
  double L = lsum + 0.5 * dot(u, u) + 0.5 * accu(log(Dg));
  double gmax = max(abs(tau * Mg + u));
  if (!want_grad) {
    return List::create(Named("L") = L, Named("u") = u,
                        Named("gmax") = gmax);
  }
  vec cg = 0.5 / Dg;
  vec eg = cg * std::pow(tau, 3) % HSg / Dg;
  vec g_alpha(alpha.n_elem, fill::zeros);
  vec qsum(n, fill::zeros);
  double g_tau = accu(2 * tau * cg % Hg - eg % Mg);
  for (uword i = 0; i < n; i++) {
    if (w[i] == 0) continue;
    const McParts& r = parts[i];
    int gi = grp[i];
    double q1 = r.l1 + cg[gi] * tau * tau * r.h1 - eg[gi] * tau * r.m1;
    double q2 = r.l2 + cg[gi] * tau * tau * r.h2 - eg[gi] * tau * r.m2;
    if (y[i] < ncat) g_alpha[y[i] - 1] += q1;
    if (y[i] > 1) g_alpha[y[i] - 2] += q2;
    qsum[i] = q1 + q2;
    g_tau -= u[gi] * qsum[i];
  }
  vec grad(alpha.n_elem + beta.n_elem + 1);
  grad.head(alpha.n_elem) = g_alpha;
  if (beta.n_elem > 0) {
    grad.subvec(alpha.n_elem, alpha.n_elem + beta.n_elem - 1) = -x.t() * qsum;
  }
  grad[grad.n_elem - 1] = g_tau;
  if (!std::isfinite(L) || !grad.is_finite()) return out_inf;
  return List::create(Named("L") = L, Named("u") = u, Named("grad") = grad,
                      Named("gmax") = gmax);
}
