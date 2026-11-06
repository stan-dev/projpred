#include <RcppArmadillo.h>

//[[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// Weighted moments of the draws (columns) within each cluster, in one pass over
// the draws for the means and, if `need_var` is `true`, one for the variances.
// `cl` holds the 0-based cluster of each draw (-1: the draw is in no cluster)
// and `wn` the weights of the draws, normalized to sum to 1 within each
// cluster. For each cluster, the function returns the weighted means of `mu`,
// `mu_offs`, `dis` and `dis^2`, and the weighted variance of `mu_offs` (about
// its weighted mean). If `need_var` is `false`, the mean of `dis^2` and the
// variance of `mu_offs` are `NA`.
// [[Rcpp::export]]
List clust_moments_c(const arma::mat& mu, const arma::mat& mu_offs,
                     const arma::vec& dis, const arma::vec& wn,
                     const IntegerVector& cl, int nclusters,
                     bool need_var) {
  const uword n_draws = mu.n_cols;
  if (mu_offs.n_rows != mu.n_rows || mu_offs.n_cols != n_draws ||
      dis.n_elem != n_draws || wn.n_elem != n_draws ||
      static_cast<uword>(cl.size()) != n_draws) {
    stop("clust_moments_c(): the dimensions of the inputs do not match");
  }
  for (uword s = 0; s < n_draws; s++) {
    if (cl[s] >= nclusters || (cl[s] < 0 && cl[s] != -1)) {
      stop("clust_moments_c(): a cluster index is out of range");
    }
  }
  mat mean_mu(mu.n_rows, nclusters, fill::zeros);
  mat mean_offs(mu_offs.n_rows, nclusters, fill::zeros);
  const double var_init = need_var ? 0.0 : NA_REAL;
  mat var_offs(mu_offs.n_rows, nclusters);
  var_offs.fill(var_init);
  vec mean_dis(nclusters, fill::zeros);
  vec mean_dis2(nclusters);
  mean_dis2.fill(var_init);
  for (uword s = 0; s < n_draws; s++) {
    const int k = cl[s];
    if (k < 0) {
      continue;
    }
    const double w = wn(s);
    mean_mu.col(k) += w * mu.col(s);
    mean_offs.col(k) += w * mu_offs.col(s);
    mean_dis(k) += w * dis(s);
    if (need_var) {
      mean_dis2(k) += w * dis(s) * dis(s);
    }
  }
  if (need_var) {
    for (uword s = 0; s < n_draws; s++) {
      const int k = cl[s];
      if (k < 0) {
        continue;
      }
      var_offs.col(k) += wn(s) * square(mu_offs.col(s) - mean_offs.col(k));
    }
  }
  return List::create(Named("mu") = mean_mu,
                      Named("mu_offs") = mean_offs,
                      Named("var_offs") = var_offs,
                      Named("dis") = NumericVector(mean_dis.begin(),
                                                   mean_dis.end()),
                      Named("dis2") = NumericVector(mean_dis2.begin(),
                                                    mean_dis2.end()));
}
