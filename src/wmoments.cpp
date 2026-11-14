#include <RcppArmadillo.h>

//[[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// Weighted moments of scaled columns of `mu` (observations x draws) with the
// observation weights `wobs` (sum `sum_w`), without temporary observations x
// draws matrices. For k = 1, ..., length(idx), the column is
// `scale[k] * mu[, idx[k]]` (`idx` 0-based). The function returns the weighted
// mean of each column (`mean`) and the weighted sum of squares about that mean
// (`ss`, second pass). A non-finite value in a column gives a non-finite `ss`
// for it, also with a weight of 0 (perf_eval_gauss() relies on this).
// [[Rcpp::export]]
List wcol_moments_c(const arma::mat& mu, const IntegerVector& idx,
                    const arma::vec& scale, const arma::vec& wobs,
                    double sum_w) {
  const uword n_cols = idx.size();
  const uword n_obs = mu.n_rows;
  if (scale.n_elem != n_cols || wobs.n_elem != n_obs) {
    stop("wcol_moments_c(): the dimensions of the inputs do not match");
  }
  for (uword k = 0; k < n_cols; k++) {
    if (idx[k] < 0 || static_cast<uword>(idx[k]) >= mu.n_cols) {
      stop("wcol_moments_c(): a column index is out of range");
    }
  }
  NumericVector mean_mu(n_cols), ss(n_cols);
  // Plain loops: dot() can call BLAS once per column, which costs more than
  // the sums for short columns.
  const double* w = wobs.memptr();
  for (uword k = 0; k < n_cols; k++) {
    const double* col = mu.colptr(idx[k]);
    const double sc = scale(k);
    double acc_mu = 0;
    for (uword i = 0; i < n_obs; i++) {
      acc_mu += w[i] * (sc * col[i]);
    }
    const double m = acc_mu / sum_w;
    double acc_ss = 0;
    for (uword i = 0; i < n_obs; i++) {
      const double d = sc * col[i] - m;
      acc_ss += w[i] * d * d;
    }
    mean_mu[k] = m;
    ss[k] = acc_ss;
  }
  return List::create(Named("mean") = mean_mu, Named("ss") = ss);
}
