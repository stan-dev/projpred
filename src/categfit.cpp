#include <RcppArmadillo.h>
#include <cmath>
#include <limits>

//[[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// Multinomial-logit fits for the augmented-data projection, with the objective
// of nnet::multinom() (without weight decay): for observation i with category
// weights w_ic, c = 0, ..., C - 1, and linear predictors eta_i0 = 0 and
// eta_ik = x_i' b_k (k = 1, ..., C - 1),
//   L(b) = -sum_i sum_c w_ic log(pi_ic),  pi_ic = exp(eta_ic) / sum_c exp(eta_ic).
// The objective is convex. nnet::multinom() minimizes it with the BFGS method
// of nnet::nnet(); this code uses Newton steps with the analytic gradient and
// Hessian and a backtracking line search.

// Objective at the coefficients `b` (p x K, K = C - 1; stored column-wise in
// the parameter vector). Also writes the probabilities pi_ik (k >= 1) into
// `pr` (n x K), from which ml_deriv() computes the derivatives.
static double ml_loss(const mat& x, const mat& w, const mat& b, mat& pr) {
  const uword n = x.n_rows;
  const uword K = b.n_cols;
  mat eta = x * b;  // n x K
  pr.set_size(n, K);
  double loss = 0.0;
  for (uword i = 0; i < n; ++i) {
    double m = 0.0;
    for (uword k = 0; k < K; ++k) m = std::max(m, eta(i, k));
    double den = std::exp(-m);
    for (uword k = 0; k < K; ++k) den += std::exp(eta(i, k) - m);
    double lden = m + std::log(den);
    // log(pi_i0) = -lden, log(pi_ik) = eta_ik - lden:
    if (w(i, 0) != 0.0) loss -= w(i, 0) * (-lden);
    for (uword k = 0; k < K; ++k) {
      double lp = eta(i, k) - lden;
      pr(i, k) = std::exp(lp);
      if (w(i, k + 1) != 0.0) loss -= w(i, k + 1) * lp;
    }
  }
  return loss;
}

// Gradient and Hessian of the objective, from the probabilities `pr` of
// ml_loss():
static void ml_deriv(const mat& x, const mat& w, const vec& nw, const mat& pr,
                     vec& grad, mat& hess) {
  const uword p = x.n_cols;
  const uword K = pr.n_cols;
  grad.set_size(p * K);
  hess.set_size(p * K, p * K);
  for (uword k = 0; k < K; ++k) {
    vec r = nw % pr.col(k) - w.col(k + 1);
    grad.subvec(k * p, (k + 1) * p - 1) = x.t() * r;
    for (uword l = k; l < K; ++l) {
      vec d = -nw % pr.col(k) % pr.col(l);
      if (l == k) d += nw % pr.col(k);
      mat h = x.t() * (x.each_col() % d);
      hess.submat(k * p, l * p, (k + 1) * p - 1, (l + 1) * p - 1) = h;
      if (l != k) {
        hess.submat(l * p, k * p, (l + 1) * p - 1, (k + 1) * p - 1) = h.t();
      }
    }
  }
}

// Fits one weight matrix `w` (n x C). The model matrix `x` contains the
// intercept column (if any) at position `icpt` (-1 if none). Writes the
// estimates into `b` (p x K); returns the objective. `conv` is 0 on
// convergence.
static double ml_fit_one(const mat& x, const mat& w, int icpt, double thresh,
                         int maxit, mat& b, int& conv, int& niter) {
  const uword p = x.n_cols;
  const uword K = w.n_cols - 1;
  vec nw = sum(w, 1);
  // Start: the intercepts of the model without predictors, otherwise zero:
  b.zeros(p, K);
  if (icpt >= 0) {
    rowvec wc = sum(w, 0);
    for (uword k = 0; k < K; ++k) {
      if (wc(0) > 0.0 && wc(k + 1) > 0.0) {
        b(icpt, k) = std::log(wc(k + 1) / wc(0));
      }
    }
  }
  vec grad;
  mat hess;
  mat pr;
  double loss = ml_loss(x, w, b, pr);
  ml_deriv(x, w, nw, pr, grad, hess);
  conv = 1;
  niter = 0;
  // Non-finite objective at the start values: not converged.
  if (!std::isfinite(loss)) return loss;
  const uword np = p * K;
  for (int it = 0; it < maxit; ++it) {
    niter = it + 1;
    vec dir;
    mat hreg = hess;
    double damp = 0.0;
    bool ok = false;
    double dscale = std::max(1e-10, max(abs(hess.diag())));
    for (int k = 0; k < 30 && !ok; ++k) {
      mat r;
      if (damp > 0.0) hreg = hess + damp * eye(np, np);
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
    mat bdir = reshape(dir, p, K);
    // Backtracking line search. The trials compute only the objective; the
    // derivatives at the accepted point come from its stored probabilities,
    // so that the point needs no second evaluation.
    double t = 1.0;
    bool accepted = false;
    double loss_accepted = 0.0;
    mat b_new, pr_new;
    for (int ls = 0; ls < 60; ++ls) {
      b_new = b + t * bdir;
      double loss_new = ml_loss(x, w, b_new, pr_new);
      if (std::isfinite(loss_new) && loss_new <= loss - 1e-4 * t * decrement) {
        accepted = true;
        loss_accepted = loss_new;
        break;
      }
      t *= 0.5;
    }
    if (!accepted) {
      if (decrement < 1e3 * thresh * std::max(1.0, std::abs(loss))) conv = 0;
      break;
    }
    b = b_new;
    loss = loss_accepted;
    pr.swap(pr_new);
    ml_deriv(x, w, nw, pr, grad, hess);
  }
  return loss;
}

// Fits the multinomial-logit model for S weight matrices, given as the n x C
// column blocks of `w_all` (n x (C * S)). The model matrix `x` (n x p) is
// shared; `icpt` is the 0-based index of its intercept column (-1 if none).
// Returns the coefficients as a p x K x S array.
// [[Rcpp::export]]
List categ_fit_multi_c(arma::mat x, arma::mat w_all, int ncat, int icpt,
                       double thresh, int maxit) {
  const uword n = x.n_rows;
  const uword p = x.n_cols;
  const uword K = ncat - 1;
  const uword S = w_all.n_cols / ncat;
  if (w_all.n_rows != n || S * ncat != w_all.n_cols) {
    Rcpp::stop("Dimensions of `w_all` do not match.");
  }
  cube b_out(p, K, S);
  vec loss_out(S);
  IntegerVector conv_out(S), niter_out(S);
  const bool x_finite = x.is_finite();
  for (uword s = 0; s < S; ++s) {
    mat w = w_all.cols(s * ncat, (s + 1) * ncat - 1);
    if (!x_finite || !w.is_finite()) {
      // Not fitted (a product with non-finite values need not give NaN in
      // every BLAS, so ml_fit_one() could not detect this case):
      b_out.slice(s).zeros();
      loss_out(s) = std::numeric_limits<double>::infinity();
      conv_out[s] = 1;
      niter_out[s] = 0;
      continue;
    }
    mat b;
    int conv, niter;
    loss_out(s) = ml_fit_one(x, w, icpt, thresh, maxit, b, conv, niter);
    b_out.slice(s) = b;
    conv_out[s] = conv;
    niter_out[s] = niter;
  }
  return List::create(Named("coef") = b_out, Named("loss") = loss_out,
                      Named("convergence") = conv_out,
                      Named("niter") = niter_out);
}
