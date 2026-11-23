// The aggregated reference probability mass function of the negative binomial
// (NB) projection (see R/negbinomial.R): for observations i and draws s with
// weights w[i] * v[s], the sum of the NB(mu(i, s), shape(i, s)) probabilities
// for y = 0, ..., y0 - 1, and for y >= y0 the sums of the mass, of (y - a_j)
// and of (y - a_j)^2 in the bins [a_j, a_{j+1} - 1] with a_0 = y0 and
// a_{j+1} = a_j + max(1, floor(c_bin * a_j^(3/4))).
//
// method 0 ("exact"): the probabilities of all y, from the ratio
// p(y + 1) / p(y) = (y + r) / (y + 1) * mu / (mu + r) (Poisson: mu / (y + 1)),
// from the mode in both directions. Each direction stops where the remaining
// mass is below `tail`.
// method 1 ("quadrature"): the same for y < y0; a bin uses a quadrature of
// the continuous extension of the pmf (add_bin_quad()).
#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <algorithm>

namespace {

const double kLogSqrt2Pi = 0.918938533204672742;

// lgamma(), digamma() and trigamma() by the Stirling series, with log(x)
// given. All calls have x >= 99 (the bins start at y0 >= 100), where the
// error is below 1e-20.
inline double lgamma_big(double x, double lx) {
  const double z = 1.0 / x;
  const double z2 = z * z;
  return (x - 0.5) * lx - x + kLogSqrt2Pi +
    z * (1.0 / 12 - z2 * (1.0 / 360 - z2 * (1.0 / 1260 - z2 / 1680)));
}
inline double digamma_big(double x, double lx) {
  const double z = 1.0 / x;
  const double z2 = z * z;
  return lx - 0.5 * z - z2 * (1.0 / 12 - z2 * (1.0 / 120 - z2 / 252));
}
// The series part of lgamma_big():
inline double lgamma_series(double x) {
  const double z = 1.0 / x;
  const double z2 = z * z;
  return z * (1.0 / 12 - z2 * (1.0 / 360 - z2 * (1.0 / 1260 - z2 / 1680)));
}
inline double trigamma_big(double x) {
  const double z = 1.0 / x;
  const double z2 = z * z;
  return z + 0.5 * z2 + z * z2 * (1.0 / 6 - z2 * (1.0 / 30 - z2 / 42));
}

struct Accum {
  int y0;
  double c_bin;
  std::vector<double> p0, edges, sm, s1, s2;
  Accum(int y0_, double c_bin_) : y0(y0_), c_bin(c_bin_), p0(y0_, 0.0),
    edges(1, static_cast<double>(y0_)) {}
  void extend_edges(std::size_t n) {
    while (edges.size() < n) {
      const double a = edges.back();
      edges.push_back(a + std::max(1.0, std::floor(c_bin *
                                                   std::pow(a, 0.75))));
    }
  }
  void size_sums(std::size_t j) {
    if (sm.size() <= j) {
      sm.resize(j + 1, 0.0);
      s1.resize(j + 1, 0.0);
      s2.resize(j + 1, 0.0);
    }
  }
  // Index of the bin of y >= y0, from the start index j:
  std::size_t bin_of(double y, std::size_t j) {
    while (true) {
      extend_edges(j + 2);
      if (y < edges[j]) {
        --j;
      } else if (y >= edges[j + 1]) {
        ++j;
      } else {
        break;
      }
    }
    size_sums(j);
    return j;
  }
  // Index of the bin of y >= y0 by bisection:
  std::size_t bin_index(double y) {
    while (edges.back() <= y) {
      extend_edges(edges.size() + 1);
    }
    const std::size_t j = static_cast<std::size_t>(
      std::upper_bound(edges.begin(), edges.end(), y) - edges.begin()) - 1;
    return bin_of(y, j);
  }
  void add_value(double y, double wp, std::size_t& j) {
    if (y < y0) {
      p0[static_cast<std::size_t>(y)] += wp;
    } else {
      j = bin_of(y, j);
      const double d = y - edges[j];
      sm[j] += wp;
      s1[j] += wp * d;
      s2[j] += wp * d * d;
    }
  }
};

// Checks for a user interrupt every 2^20 steps:
inline void check_interrupt(std::size_t& steps) {
  if ((++steps & 0xFFFFF) == 0) {
    Rcpp::checkUserInterrupt();
  }
}

// The NB(m, r) distribution (r = Inf: Poisson). For r >= 99 (`big_r`),
// lgamma(t + r) - lgamma(r) is one Stirling difference
// (t + r - 1/2) log1p(t / r) + t log(r) - t + series(t + r) - series(r), so
// that the two large terms r log(r) do not cancel.
struct NbDist {
  double m, r, q, lq, lc, mode, lr0, s0;
  bool pois, big_r;
  NbDist(double m_, double r_) : m(m_), r(r_), pois(!R_FINITE(r_)),
    big_r(R_FINITE(r_) && r_ >= 99) {
    q = pois ? 0.0 : m / (m + r);
    lq = pois ? std::log(m) : -std::log1p(r / m);
    lc = pois ? -m : -r * std::log1p(m / r) - (big_r ? 0.0 : R::lgammafn(r));
    mode = pois ? std::floor(m) : (r > 1 ? std::floor((r - 1) * m / r) : 0.0);
    lr0 = big_r ? std::log(r) : 0.0;
    s0 = big_r ? lgamma_series(r) : 0.0;
  }
  // lgamma(t + r), minus lgamma(r) if big_r (log(t + r) = lr is used only
  // otherwise):
  double lgamma_tr(double t, double lr) const {
    if (big_r) {
      return (t + r - 0.5) * std::log1p(t / r) + t * lr0 - t +
        lgamma_series(t + r) - s0;
    }
    return lgamma_big(t + r, lr);
  }
  double ratio(double y) const {
    return pois ? m / (y + 1.0) : (y + r) / (y + 1.0) * q;
  }
  // For r >= 99, R::dnbinom_mu() loses accuracy (about 1e-9 relative for
  // r = 1e8), so the pmf comes from the Stirling difference:
  double pmf(double y) const {
    if (pois) {
      return R::dpois(y, m, 0);
    }
    if (big_r) {
      return std::exp(lc + lgamma_tr(y, 0.0) - R::lgammafn(y + 1) + y * lq);
    }
    return R::dnbinom_mu(y, r, m, 0);
  }
  // The continuous extension for t >= 99 (Stirling series): the log pmf and
  // its first derivative in `lf` and `dlf`.
  void lf_dlf(double t, double& lf, double& dlf) const {
    const double x1 = t + 1.0;
    const double l1 = std::log(x1);
    if (pois) {
      lf = lc + t * lq - lgamma_big(x1, l1);
      dlf = lq - digamma_big(x1, l1);
    } else {
      const double xr = t + r;
      const double lr = std::log(xr);
      lf = lc + lgamma_tr(t, lr) - lgamma_big(x1, l1) + t * lq;
      dlf = digamma_big(xr, lr) - digamma_big(x1, l1) + lq;
    }
  }
  double lf(double t) const {
    const double x1 = t + 1.0;
    if (pois) {
      return lc + t * lq - lgamma_big(x1, std::log(x1));
    }
    const double xr = t + r;
    return lc + lgamma_tr(t, big_r ? 0.0 : std::log(xr)) -
      lgamma_big(x1, std::log(x1)) + t * lq;
  }
  double d2lf(double t) const {
    return pois ? -trigamma_big(t + 1.0) :
      trigamma_big(t + r) - trigamma_big(t + 1.0);
  }
};

// The pmf and its log derivative at a bin edge t (a half-integer), with a
// cache of the last two edges (neighboring bins share an edge, in both
// directions of the bin loops).
struct EdgeCache {
  double t[2] = {-1.0, -1.0}, f[2] = {0.0, 0.0}, d[2] = {0.0, 0.0};
  int next = 0;
  void get(const NbDist& nb, double t_new, double& f_out, double& d_out) {
    for (int k = 0; k < 2; ++k) {
      if (t[k] == t_new) {
        f_out = f[k];
        d_out = d[k];
        return;
      }
    }
    double lf;
    nb.lf_dlf(t_new, lf, d[next]);
    f[next] = std::exp(lf);
    t[next] = t_new;
    f_out = f[next];
    d_out = d[next];
    next = 1 - next;
  }
};

// Adds the NB(m, r) probabilities of all y (method 0).
void add_pair(Accum& acc, const NbDist& nb, double wt, double tail,
              std::size_t& steps) {
  const double p_mode = nb.pmf(nb.mode);
  std::size_t j = 0;
  std::size_t j_mode = 0;
  // From the mode upwards:
  double p = p_mode;
  for (double y = nb.mode;; y += 1.0) {
    acc.add_value(y, wt * p, j);
    if (y == nb.mode) {
      j_mode = j;
    }
    const double ratio = nb.ratio(y);
    p *= ratio;
    // Beyond the mode, the later ratios are at most `rmax` (the current ratio
    // for r >= 1 and Poisson, where the ratio decreases; q for r < 1, where it
    // increases towards q), so the remaining mass is at most p / (1 - rmax):
    const double rmax = (!nb.pois && nb.r < 1) ? nb.q : ratio;
    if ((rmax < 1 && p <= tail * (1 - rmax)) || p == 0) {
      break;
    }
    check_interrupt(steps);
  }
  // From the mode downwards (below the mode, the probabilities decrease at
  // least as fast as a geometric series with the ratio of the current step):
  p = p_mode;
  j = j_mode;
  for (double y = nb.mode - 1.0; y >= 0; y -= 1.0) {
    const double ratio = nb.ratio(y);
    p /= ratio;
    acc.add_value(y, wt * p, j);
    const double back = 1.0 / ratio;
    if ((back < 1 && p <= tail * (1 - back)) || p == 0) {
      break;
    }
    check_interrupt(steps);
  }
}

// Adds bin j by the quadrature. The integral of the continuous extension f
// over [a - 1/2, b + 1/2] with the Euler-Maclaurin correction
// -(g'(b + 1/2) - g'(a - 1/2)) / 24 gives the midpoint sum, for g = f,
// (t - a) f and (t - a)^2 f. The integral uses the 4-point Gauss-Legendre
// rule on k panels of width at most a quarter of the local SD
// |d2lf|^(-1/2) at the center. The exact sum is used for a bin with fewer
// than 16 values, if the log pmf is steep at an edge (|dlf| > 0.5) or if the
// panels need at least as many points as the bin has values. (At |dlf| = 0.5,
// the next Euler-Maclaurin term, 7/5760 times the third derivative of g, is
// about 1.5e-4 times f at the edge.) Returns the pmf at the edges in `fl` and
// `fu` (for the tail bounds).
void add_bin_quad(Accum& acc, const NbDist& nb, std::size_t j, double wt,
                  EdgeCache& cache, double& fl, double& fu) {
  const double a = acc.edges[j];
  const double b = acc.edges[j + 1] - 1.0;
  const double tl = a - 0.5;
  const double tu = b + 0.5;
  double dl, du;
  cache.get(nb, tl, fl, dl);
  cache.get(nb, tu, fu, du);
  const double len = tu - tl;
  double sm = 0, s1 = 0, s2 = 0;
  double k_panels = 0;
  bool exact = len < 16 || std::fabs(dl) > 0.5 || std::fabs(du) > 0.5;
  if (!exact) {
    const double curv = std::fabs(nb.d2lf((tl + tu) / 2));
    k_panels = std::max(1.0, std::ceil(4.0 * len * std::sqrt(curv)));
    exact = 4.0 * k_panels >= len;
  }
  if (exact) {
    double p = std::exp(nb.lf(a));
    for (double y = a; y <= b; y += 1.0) {
      const double dd = y - a;
      sm += p;
      s1 += p * dd;
      s2 += p * dd * dd;
      p *= nb.ratio(y);
    }
  } else {
    static const double xg[4] = {-0.8611363115940526, -0.3399810435848563,
                                 0.3399810435848563, 0.8611363115940526};
    static const double wg[4] = {0.3478548451374538, 0.6521451548625461,
                                 0.6521451548625461, 0.3478548451374538};
    const int np = static_cast<int>(k_panels);
    const double h = len / (2.0 * np);
    for (int ip = 0; ip < np; ++ip) {
      const double c = tl + (2 * ip + 1) * h;
      for (int k = 0; k < 4; ++k) {
        const double t = c + h * xg[k];
        const double f = h * wg[k] * std::exp(nb.lf(t));
        const double dd = t - a;
        sm += f;
        s1 += f * dd;
        s2 += f * dd * dd;
      }
    }
    const double el = tl - a;
    const double eu = tu - a;
    sm -= (fu * du - fl * dl) / 24.0;
    s1 -= ((fu + eu * fu * du) - (fl + el * fl * dl)) / 24.0;
    s2 -= ((2 * eu * fu + eu * eu * fu * du) -
           (2 * el * fl + el * el * fl * dl)) / 24.0;
  }
  acc.sm[j] += wt * sm;
  acc.s1[j] += wt * s1;
  acc.s2[j] += wt * s2;
}

// Adds the NB(m, r) probabilities with the quadrature for the bins (method 1).
void add_pair_quad(Accum& acc, const NbDist& nb, double wt, double tail,
                   std::size_t& steps) {
  const double y0 = acc.y0;
  const double mode = nb.mode;
  // The bound of the mass above y for y beyond the mode, from p(y + 1) (or
  // from an upper bound of it):
  auto up_bound = [&](double p_next, double y) {
    const double rt = nb.ratio(y + 1.0);
    const double rmax = (!nb.pois && nb.r < 1) ? nb.q : rt;
    return rmax < 1 ? p_next / (1 - rmax) : R_PosInf;
  };
  bool need_bins = true;
  std::size_t jdummy = 0;
  if (mode < y0) {
    // y = mode, ..., y0 - 1 upwards, then y < mode downwards:
    const double p_mode = nb.pmf(mode);
    double p = p_mode;
    for (double y = mode; y < y0; y += 1.0) {
      acc.add_value(y, wt * p, jdummy);
      p *= nb.ratio(y);
      if (up_bound(p, y) <= tail || p == 0) {
        need_bins = false;
        break;
      }
    }
    p = p_mode;
    for (double y = mode - 1.0; y >= 0; y -= 1.0) {
      const double ratio = nb.ratio(y);
      p /= ratio;
      acc.add_value(y, wt * p, jdummy);
      const double back = 1.0 / ratio;
      if ((back < 1 && p <= tail * (1 - back)) || p == 0) {
        break;
      }
    }
  } else {
    // y = y0 - 1, ..., 0 downwards (all below the mode):
    double p = std::exp(nb.lf(y0 - 1.0));
    for (double y = y0 - 1.0; y >= 0; y -= 1.0) {
      if (y < y0 - 1.0) {
        p /= nb.ratio(y);
      }
      acc.add_value(y, wt * p, jdummy);
      const double back = y >= 1 ? 1.0 / nb.ratio(y - 1.0) : 0.0;
      if ((back < 1 && p <= tail * (1 - back)) || p == 0) {
        break;
      }
    }
  }
  if (!need_bins) {
    return;
  }
  EdgeCache cache;
  double fl, fu;
  const std::size_t j_start = acc.bin_index(std::max(mode, y0));
  // Upwards from the bin of the mode. Beyond the mode, the continuous
  // extension decreases, so p(b + 1) <= f(b + 1/2) = fu:
  for (std::size_t j = j_start;; ++j) {
    acc.extend_edges(j + 2);
    acc.size_sums(j);
    add_bin_quad(acc, nb, j, wt, cache, fl, fu);
    const double b = acc.edges[j + 1] - 1.0;
    if (b + 1.0 > mode && (up_bound(fu, b) <= tail || fu == 0)) {
      break;
    }
    check_interrupt(steps);
  }
  // Downwards to bin 0. Below the mode, the continuous extension increases,
  // so the mass below the lower edge a of the last bin is at most
  // f(a - 1/2) / (1 - back):
  double d_unused;
  cache.get(nb, acc.edges[j_start] - 0.5, fl, d_unused);
  for (std::size_t j = j_start; j > 0; --j) {
    const double a_prev = acc.edges[j];
    const double back = 1.0 / nb.ratio(a_prev - 2.0);
    if ((back < 1 && fl <= tail * (1 - back)) || fl == 0) {
      break;
    }
    add_bin_quad(acc, nb, j - 1, wt, cache, fl, fu);
    check_interrupt(steps);
  }
}

}  // namespace

// [[Rcpp::export]]
Rcpp::List nb_gsum_c(const Rcpp::NumericMatrix& mu,
                     const Rcpp::NumericMatrix& shape,
                     const Rcpp::NumericVector& w,
                     const Rcpp::NumericVector& v,
                     int y0, double c_bin, double tail, int method) {
  const int n = mu.nrow();
  const int n_draws = mu.ncol();
  if (shape.nrow() != n || shape.ncol() != n_draws || w.size() != n ||
      v.size() != n_draws) {
    Rcpp::stop("nb_gsum_c(): the dimensions of the arguments do not match.");
  }
  if (method != 0 && method != 1) {
    Rcpp::stop("nb_gsum_c(): unknown method.");
  }
  // The Stirling series of the quadrature needs y0 >= 100:
  if (y0 < 1 || (method == 1 && y0 < 100) || !(c_bin > 0)) {
    Rcpp::stop("nb_gsum_c(): invalid `y0` or `c_bin`.");
  }
  Accum acc(y0, c_bin);
  std::size_t steps = 0;
  for (int s = 0; s < n_draws; ++s) {
    for (int i = 0; i < n; ++i) {
      const double wt = w[i] * v[s];
      // A pair with weight 0 (for example the left-out observation of a LOO
      // fold) does not enter, and its mean and shape are not checked:
      if (!(wt > 0)) {
        continue;
      }
      const double m = mu(i, s);
      const double r = shape(i, s);
      if (!R_FINITE(m) || m < 0 || ISNAN(r) || r <= 0) {
        Rcpp::stop("nb_gsum_c(): invalid mean or shape.");
      }
      if (m == 0) {
        acc.p0[0] += wt;
        continue;
      }
      const NbDist nb(m, r);
      if (method == 1) {
        add_pair_quad(acc, nb, wt, tail, steps);
      } else {
        add_pair(acc, nb, wt, tail, steps);
      }
    }
  }
  const std::size_t nb = acc.sm.size();
  Rcpp::NumericMatrix sums(nb, 3);
  for (std::size_t j = 0; j < nb; ++j) {
    sums(j, 0) = acc.sm[j];
    sums(j, 1) = acc.s1[j];
    sums(j, 2) = acc.s2[j];
  }
  acc.extend_edges(nb + 1);
  return Rcpp::List::create(
    Rcpp::Named("p0") = Rcpp::NumericVector(acc.p0.begin(), acc.p0.end()),
    Rcpp::Named("sums") = sums,
    Rcpp::Named("edges") = Rcpp::NumericVector(acc.edges.begin(),
                                               acc.edges.begin() + nb + 1));
}
