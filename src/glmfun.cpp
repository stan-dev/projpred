#include <iostream>
#include <RcppArmadillo.h>
#include <cmath>
#include <limits>
#include <string>
#include <sstream>
#include <vector>

//[[Rcpp::depends(RcppArmadillo)]]


using namespace Rcpp;
using namespace arma;



/**
 * Returns the value of the penalty term in the elastic net regularization.
 */
double elnet_penalty(vec beta, // coefficients
                     double lambda, // regularization parameter
                     double alpha, // elastic net mixing parameter
                     vec penalty) // relative penalties for the variables 
{
  double value;
  uvec fin = find_finite(penalty);
  value = lambda*sum(penalty.elem(fin) % (0.5*(1-alpha)*square(beta.elem(fin))
                                            + alpha*(abs(beta.elem(fin))) ) );
  return(value);
}


/** Returns the value of the regularized quadratic approximation to the loss function
that is to be minimized iteratively:
L = 0.5*sum_i{ w_i*(z_i - f_i)^2 } + lambda*{ 0.5*(1-alpha)*||beta||_2^2 + alpha*||beta||_1 }
*/
double loss_approx(vec& beta,    // coefficients
                   vec& f,       // latent values
                   vec& z,       // locations for pseudo obsevations
                   vec& w,       // weights of the pseudo observations (inverse-variances)
                   double lambda, // regularization parameter
                   double alpha, // elastic net mixing parameter
                   vec& penalty) // relative penalties for the variables  
{
  double loss;
  uvec fin = find_finite(penalty);
  loss = 0.5*sum(w % square(z-f)) + elnet_penalty(beta,lambda,alpha,penalty);
  return loss;
}




/** Updates the regression coefficients and the intercept (unless excluded) based on the
current quadratic approximation to the loss function. This is done via the 'soft-thresholding'
as described by Friedman et. al (2009). Performs either one pass through the specified set
of varibles or iterates until convergence.
*/
void coord_descent(	vec& beta, // regression coefficients
                    double& beta0, // intercept
                    vec& f, // latent values
                    mat& x, // input matrix
                    vec& z, // locations for pseudo obsevations
                    vec& w, // weights of the pseudo observations (inverse-variances)
                    double& lambda, // regularization parameter
                    double& alpha, // elastic net mixing parameter
                    vec& penalty, // relative penalties for the variables
                    bool intercept, // whether to use intercept
                    std::set<size_t>& varind, // which coefficients are updated
                    std::set<size_t>& active_set, // active set, may change if some variables enter or leave
                    bool until_convergence, // true = until convergence, false = one pass through varind
                    int& npasses, // counts total passes through the variables
                    double tol, // stop when change in the loss is smaller than this
                    int maxiter = 1000) // maximum number of iterations (passes) through varind
{
  
  int iter = 0;
  double loss,loss_old;
  size_t j;
  double h;
  
  // initial loss
  loss_old = loss_approx(beta,f,z,w,lambda,alpha,penalty);
  
  // auxiliary that will be used later on
  // double lam_alpha = lambda*alpha;
  // double lam_oneminus_alpha = lambda*(1-alpha);
  
  while (iter < maxiter) {
    
    // update the intercept
    if (intercept) {
      f = f - beta0;
      beta0 = sum(w % (z - f)) / sum(w);
      f = f + beta0;
    }
    
    active_set.clear();
    
    for (std::set<size_t>::iterator it=varind.begin(); it!=varind.end(); ++it) {
      
      // update the regression coefficients via 'soft-thresholding'
      
      // varible index
      j = *it;
      
      f = f - beta(j)*x.col(j);
      h = sum( w % x.col(j) % (z - f) ); // auxiliary variable
      
      if (fabs(h) <= penalty(j)*lambda*alpha) {
        beta(j) = 0.0;
      } else if (h > 0) {
        beta(j) = (h - penalty(j)*lambda*alpha) / ( sum(w % square(x.col(j))) + penalty(j)*lambda*(1-alpha) );
        active_set.insert(j);
      } else {
        beta(j) = (h + penalty(j)*lambda*alpha) / ( sum(w % square(x.col(j))) + penalty(j)*lambda*(1-alpha) );
        active_set.insert(j);
      }
      f = f + beta(j)*x.col(j);
    }
    
    ++iter;
    ++npasses;
    loss = loss_approx(beta,f,z,w,lambda,alpha,penalty);
    
    if (until_convergence) {
      if (loss_old-loss < tol) {
        break;
      } else {
        // continue iterating
        loss_old = loss;
      }
    } else {
      break;
    }
  }
  
  if (iter == maxiter)
    Rcpp::Rcout << "Warning: maximum number of iterations reached in coordinate descent. Results can be inaccurate!\n";
  
}



/** Computes the whole elastic-net regularization path given the grid of values to lambda.
Assumes that the lambda grid is selected carefully and utilizes the function pseudo_obs
that returns the pseudo-observations corresponding to the quadratic approximation to the
loss function for a given vector of latent values (see elnetfun.R).
*/
// [[Rcpp::export]]
List glm_elnet_c(arma::mat x, // input matrix
                 Function pseudo_obs, // R-function returning the pseudo-data based on the quadratic approximation
                 arma::vec lambda, // grid for the regularization parameter
                 double alpha, // elastic net mixing parameter
                 bool intercept, // whether to use intercept
                 arma::vec penalty, // relative penalties for the variables
                 double thresh, // threshold for determining the convergence
                 int qa_updates_max, // maximum for the total number of quadratic approximation updates
                 int pmax, // stop computation when the active set size is equal or greater than this
                 bool pmax_strict, // if true, then the active set size of the last beta is always at most pmax
                 arma::vec beta, // initial value for the regression coefficients
                 double beta0, // initial value for the intercept
                 arma::vec w0, // initial guess for the weights of the pseudo-gaussian observations (needed for Student-t model)
                 int as_updates_max = 50) // maximum number of active set updates for one quadratic approximation
{
  
  // for gaussian pseudo data
  List obs; 
  vec z; // observations
  vec w; // weights (inverse variances)
  
  size_t D = x.n_cols; // number of inputs
  size_t pmaxu = (size_t) pmax; // converting pmax to unsigned int (avoids some compiler warnings)
  int nlam = lambda.size();
  double lam; // temporary varible for fixed lambda
  int k; // lambda index
  int qau; // number of quadratic approximation updates
  int asu; // number of active set updates (for a given quadratic approximation)
  
  
  // for storing the whole solution path
  rowvec beta0_path(nlam);
  mat beta_path(D,nlam);
  mat w_path(x.n_rows,nlam);
  beta0_path.zeros();
  beta_path.zeros();
  int npasses = 0; // counts how many times the coefficient vector is looped through
  urowvec qa_updates(nlam);
  qa_updates.zeros();
  urowvec as_updates(nlam);
  as_updates.zeros();
  
  
  // initialization
  if (!intercept)	beta0 = 0; // ensure intercept is zero when it is not used
  vec f = x*beta + beta0;
  std::set<size_t> active_set; 
  std::set<size_t> active_set_old;
  std::set<size_t> varind_all; // a constant set containing indices of all the variables
  for (size_t j=0; j<D; j++)
    varind_all.insert(j);
  
  
  obs = pseudo_obs(f,w0);
  z = as<vec>(obs["z"]);
  w = as<vec>(obs["wobs"]);
  double loss_initial = loss_approx(beta, f, z, w, lambda(0), alpha, penalty); // initial loss
  double loss_old = loss_initial; // will be updated iteratively
  double loss; // will be updated iteratively
  double tol = thresh*fabs(loss_initial); // convergence criterion for coordinate descent
  
  // loop over lambda values
  for (k=0; k<nlam; ++k) {
    
    lam = lambda(k);
    
    qau = 0;
    while (qau < qa_updates_max) {
      
      // update the quadratic likelihood approximation
      obs = pseudo_obs(f,w);
      z = as<vec>(obs["z"]);
      w = as<vec>(obs["wobs"]);
      ++qau;
      
      // current value of the (approximate) loss function
      loss_old = loss_approx(beta, f, z, w, lam, alpha, penalty);
      // loss_old = ((double) obs["loss"]) + elnet_penalty(beta, lam, alpha, penalty);
      
      // run the coordinate descent until convergence for the current
      // quadratic approximation
      asu = 0;
      while (asu < as_updates_max) {
      	
        // iterate within the current active set until convergence (this might update 
        // active_set_old, if some variable goes to zero)
        coord_descent(beta, beta0, f, x, z, w, lam, alpha, penalty, intercept, active_set, active_set_old, true, npasses, tol);
      	
        // perfom one pass over all the variables and check if the active set changes 
        // (this might update active_set)
        coord_descent(beta, beta0, f, x, z, w, lam, alpha, penalty, intercept, varind_all, active_set, false, npasses, tol);
        
        ++asu;
        
        if (active_set==active_set_old) {
          // active set did not change so convergence reached
          // (for the current quadratic approximation to the loss function)
          break;
        }
      }
      as_updates(k) = as_updates(k) + asu;
      
      // the loss after updating the coefficients
      loss = loss_approx(beta, f, z, w, lam, alpha, penalty);
      // obs = pseudo_obs(f,w);
      // loss = ((double) obs["loss"]) + elnet_penalty(beta, lam, alpha, penalty);

      // check if converged
      if (fabs(loss_old-loss) < tol) {
      // if (loss_old-loss < tol) {
        // convergence reached; proceed to the next lambda value
        break;
      }
    }
    // store the current solution
    beta0_path(k) = beta0;
    beta_path.col(k) = beta;
    w_path.col(k) = w;
    qa_updates(k) = qau;
    
    if (qau == qa_updates_max && qa_updates_max > 1)
      Rcpp::Rcout << "glm_elnet warning: maximum number of quadratic approximation updates reached. Results can be inaccurate.\n";
    
    if ((alpha > 0.0) && (active_set.size() >= pmaxu)) {
      // obtained solution with at least pmax variables and penalty is not ridge, so terminate
      if (pmax_strict) {
        // return solutions only up to the previous lambda value
        beta0_path = beta0_path.head(k);
        beta_path = beta_path.head_cols(k);
      } else {
        // return solutions up to the current lambda value
        beta0_path = beta0_path.head(k+1);
        beta_path = beta_path.head_cols(k+1);
      }
      break;
    }
  }
  
  return List::create(beta_path, beta0_path, w_path, npasses, qa_updates, as_updates);
}





/** Bounded least squares: minimizes 0.5*b'Qb - q'b subject to b(cons) >= 0,
 * with a primal active-set method (Lawson and Hanson, 1974). On input, `b` is
 * the start (its constrained elements are set to max(b, 0)); on output, the
 * solution. Returns false if the maximum number of iterations is reached.
 */
bool bls_active_set(vec& b,
                    const mat& Q,
                    const vec& q,
                    const uvec& cons,
                    double tol = 1e-10,
                    double tol_grad_rel = 1e-8,
                    int iter_max = 500)
{
  int D = q.n_elem;
  int n_cons = cons.n_elem;
  int j, k;
  // act(j) = 1: element j is held at its bound 0
  uvec act(D); act.zeros();
  for (k=0; k<n_cons; ++k) {
    j = cons(k);
    if (b(j) <= 0) {
      b(j) = 0;
      act(j) = 1;
    }
  }
  double tol_grad = tol_grad_rel*std::max(1.0, max(abs(q)));
  iter_max = std::max(iter_max, 3*n_cons);
  int j_rel = -1; // the bound released in the previous iteration, if any
  vec b_new(D);
  vec grad(D);
  for (int iter=0; iter<iter_max; ++iter) {
    uvec free_idx = find(act == 0);
    b_new.zeros();
    if (free_idx.n_elem > 0)
      b_new.elem(free_idx) = solve(Q.submat(free_idx, free_idx), q.elem(free_idx));
    // if the released element does not become positive (possible only by
    // rounding errors), the current `b` is the solution (this avoids cycling)
    if (j_rel >= 0 && b_new(j_rel) <= tol) {
      act(j_rel) = 1;
      return true;
    }
    j_rel = -1;
    // the step size to the first constrained element that becomes negative
    double alpha = 1.0;
    bool blocked = false;
    for (k=0; k<n_cons; ++k) {
      j = cons(k);
      if (act(j) == 0 && b_new(j) < -tol) {
        alpha = std::min(alpha, b(j)/(b(j) - b_new(j)));
        blocked = true;
      }
    }
    if (!blocked) {
      // values in [-tol, 0) are rounding errors; set them to the bound so that
      // `b` is feasible
      b = b_new;
      for (k=0; k<n_cons; ++k) {
        j = cons(k);
        if (b(j) < 0)
          b(j) = 0;
      }
      // release the bound with the most negative gradient, if any
      grad = Q*b - q;
      double grad_min = -tol_grad;
      for (k=0; k<n_cons; ++k) {
        j = cons(k);
        if (act(j) == 1 && grad(j) < grad_min) {
          grad_min = grad(j);
          j_rel = j;
        }
      }
      if (j_rel < 0)
        return true;
      act(j_rel) = 0;
    } else {
      b = b + alpha*(b_new - b);
      for (k=0; k<n_cons; ++k) {
        j = cons(k);
        if (b(j) <= tol) {
          b(j) = 0;
          act(j) = 1;
        }
      }
    }
  }
  return false;
}

/**
 * Wrapper for bls_active_set that can be called from R: minimizes
 * 0.5*b'Qb - q'b subject to b(cons) >= 0 (`cons`: 1-based indices), starting
 * from `b_init`.
 */
// [[Rcpp::export]]
List bls_c(arma::mat Q, arma::vec q, Rcpp::IntegerVector cons, arma::vec b_init)
{
  int D = q.n_elem;
  if ((int) Q.n_rows != D || (int) Q.n_cols != D || (int) b_init.n_elem != D)
    Rcpp::stop("Q, q and b_init must have matching dimensions.");
  uvec cons0(cons.size());
  for (int k=0; k<cons.size(); ++k) {
    if (cons[k] < 1 || cons[k] > D)
      Rcpp::stop("cons must contain indices of q.");
    cons0(k) = cons[k] - 1;
  }
  vec b = b_init;
  bool ok = bls_active_set(b, Q, q, cons0);
  return List::create(Named("b") = b, Named("converged") = ok);
}

/**
 * Pseudo-observations of the quadratic approximation to the loss function at
 * the latent values f (see pseudo_data() in R): locations z, weights w, the
 * pointwise gradient of the loss w.r.t. f, and the loss.
 */
struct PseudoObs {
  virtual ~PseudoObs() {}
  virtual void eval(const vec& f, const vec& wprev, vec& z, vec& w, vec& grad,
                    double& loss) = 0;
};

// Pseudo-observations from an R function (any family supported by
// pseudo_data()).
struct PseudoObsR : public PseudoObs {
  Function fn;
  PseudoObsR(Function fn_) : fn(fn_) {}
  void eval(const vec& f, const vec& wprev, vec& z, vec& w, vec& grad,
            double& loss) {
    List obs = fn(f, wprev);
    z = as<vec>(obs["z"]);
    w = as<vec>(obs["wobs"]);
    grad = as<vec>(obs["grad"]);
    loss = (double) obs["loss"];
  }
};

// Family and link codes of PseudoObsNative (see native_pseudo_codes() in R).
enum { FAM_GAUSSIAN = 1, FAM_BINOMIAL = 2, FAM_POISSON = 3 };
enum { LINK_IDENTITY = 1, LINK_LOG = 2, LINK_LOGIT = 3, LINK_PROBIT = 4,
       LINK_CLOGLOG = 5, LINK_CAUCHIT = 6 };

// pmax() and pmin() of R for a scalar bound (NaN propagates):
static inline double r_pmax(double x, double bound) {
  return (std::isnan(x) || x >= bound) ? x : bound;
}
static inline double r_pmin(double x, double bound) {
  return (std::isnan(x) || x <= bound) ? x : bound;
}

// Pseudo-observations computed in C++. linkinv() and mu.eta() transcribe the
// link functions of R's stats package (make.link(), and for the logit link
// logit_linkinv() and logit_mu_eta() in src/library/stats/src/family.c of R,
// GPL-2 | GPL-3). The variances and deviances are those of the stats families
// and of extend_family():
struct PseudoObsNative : public PseudoObs {
  vec y, weights, offset;
  int fam, link;
  PseudoObsNative(const vec& y_, const vec& weights_, const vec& offset_,
                  int fam_, int link_) :
    y(y_), weights(weights_), offset(offset_), fam(fam_), link(link_) {}

  // linkinv() and mu.eta() of the link:
  void link_funs(double eta, double& mu, double& dmu) const {
    const double eps = std::numeric_limits<double>::epsilon();
    switch (link) {
    case LINK_IDENTITY:
      mu = eta;
      dmu = 1.0;
      break;
    case LINK_LOG:
      mu = r_pmax(std::exp(eta), eps);
      dmu = mu;
      break;
    case LINK_LOGIT: {
      // As logit_linkinv() and logit_mu_eta() in R's stats/src/family.c:
      const double thr = 30.0;
      double tmp = (eta < -thr) ? eps : ((eta > thr) ? 1.0 / eps : std::exp(eta));
      mu = tmp / (1.0 + tmp);
      double opexp = 1.0 + std::exp(eta);
      dmu = (eta > thr || eta < -thr) ? eps : std::exp(eta) / (opexp * opexp);
      break;
    }
    case LINK_PROBIT: {
      const double thr = -R::qnorm(eps, 0.0, 1.0, 1, 0);
      double eta_c = r_pmin(r_pmax(eta, -thr), thr);
      mu = R::pnorm(eta_c, 0.0, 1.0, 1, 0);
      dmu = r_pmax(R::dnorm(eta, 0.0, 1.0, 0), eps);
      break;
    }
    case LINK_CLOGLOG: {
      mu = r_pmax(r_pmin(-std::expm1(-std::exp(eta)), 1.0 - eps), eps);
      double eta_c = r_pmin(eta, 700.0);
      dmu = r_pmax(std::exp(eta_c) * std::exp(-std::exp(eta_c)), eps);
      break;
    }
    case LINK_CAUCHIT: {
      const double thr = -R::qcauchy(eps, 0.0, 1.0, 1, 0);
      double eta_c = r_pmin(r_pmax(eta, -thr), thr);
      mu = R::pcauchy(eta_c, 0.0, 1.0, 1, 0);
      dmu = r_pmax(R::dcauchy(eta, 0.0, 1.0, 0), eps);
      break;
    }
    default:
      Rcpp::stop("Unknown link code.");
    }
  }

  void eval(const vec& f, const vec& wprev, vec& z, vec& w, vec& grad,
            double& loss) {
    size_t n = f.n_elem;
    z.set_size(n);
    w.set_size(n);
    grad.set_size(n);
    // R's sum() accumulates in long double:
    long double dev_sum = 0.0;
    double mu, dmu, var, dev;
    for (size_t i = 0; i < n; ++i) {
      link_funs(f(i) + offset(i), mu, dmu);
      z(i) = f(i) + (y(i) - mu) / dmu;
      switch (fam) {
      case FAM_GAUSSIAN:
        var = 1.0;
        dev = -2.0 * weights(i) * (-0.5 * ((y(i) - mu) * (y(i) - mu)) - 0.0);
        break;
      case FAM_BINOMIAL:
        var = mu * (1.0 - mu);
        dev = -2.0 * (weights(i) * (y(i) * std::log(mu) +
                                    (1.0 - y(i)) * std::log(1.0 - mu)));
        break;
      case FAM_POISSON:
        var = mu;
        dev = -2.0 * (weights(i) * (y(i) * std::log(mu) - mu));
        break;
      default:
        Rcpp::stop("Unknown family code.");
      }
      w(i) = (weights(i) * (dmu * dmu)) / var;
      grad(i) = (-w(i)) * (z(i) - f(i));
      dev_sum += dev;
    }
    loss = 0.5 * (double) dev_sum;
  }
};


/** Internal function that gives the output in c++ way (writes into allocated memory).
 * See glm_ridge_c for a wrapper that is being called from R.
 * The optional constraints sign(k)*beta(cons(k)) >= 0 (indices including the
 * intercept, if any) make each least-squares step a bounded least-squares step.
 */
void glm_ridge( vec& beta,      // output: regression coefficients (contains intercept)
                double& loss,   // output: value of the loss function
                vec& w, 				// output: weights of the pseudo-gaussian observations at the optimum (needed for Student-t model)
                int& qau,       // output: number of quadratic approximation updates 
                std::vector<std::string>& msgs, // output: warning messages
                arma::mat x,
                PseudoObs& pseudo_obs,
                double lambda,
                bool intercept,
                arma::vec penalty, // relative penalties for the variables
                double thresh,
                int qa_updates_max,
                int ls_iter_max=50,
                bool debug=false,
                uvec cons=uvec(),
                vec signs=vec())
{
  
  if (intercept) {
    // add a vector of ones to x and set zero penalty for the intercept
    x = join_horiz(ones<vec>(x.n_rows), x);
    penalty = join_vert(zeros<vec>(1), penalty); 
  }
  
  int n = x.n_rows;
  int D = x.n_cols;
  int n_cons = cons.n_elem;
  bool bls_ok = true;
  // with constraints, flip the signs of the constrained columns (and
  // coefficients) so that all bounds are beta(cons(k)) >= 0, and start from a
  // feasible point (the penalty does not depend on the signs)
  for (int k=0; k<n_cons; ++k) {
    x.col(cons(k)) *= signs(k);
    beta(cons(k)) = std::max(0.0, signs(k)*beta(cons(k)));
  }
  int ls_iter; // counts linesearch iterations
  int j;
  double t; // step size in line search 
  double a = 0.1; // backtracking line search parameters a and b (see Boyd and Vandenberghe, 2004)
  double b = 0.5; 
  
  // Variables with an infinite penalty (zero variance) keep the coefficient 0
  // and stay out of the linear system, which would otherwise contain Inf.
  uvec act = find_finite(penalty);
  uvec inact = find_nonfinite(penalty);
  vec penalty_grad = penalty;
  penalty_grad.elem(inact).zeros();
  beta.elem(inact).zeros();
  // positions of the constrained variables among the active ones (an inactive
  // variable keeps 0, which satisfies its sign constraint); used only if some
  // variables are inactive
  uvec cons_act;
  if (n_cons > 0 && inact.n_elem > 0) {
    ivec pos_act(D); pos_act.fill(-1);
    for (uword i=0; i<act.n_elem; ++i)
      pos_act(act(i)) = i;
    std::vector<uword> cons_tmp;
    for (int k=0; k<n_cons; ++k)
      if (pos_act(cons(k)) >= 0)
        cons_tmp.push_back(pos_act(cons(k)));
    cons_act = conv_to<uvec>::from(cons_tmp);
  }
  
  // initialization
  vec beta_new(D); beta_new.zeros();
  vec dbeta(D); dbeta.zeros();
  vec grad(D); grad.zeros(); // gradient of the negative log likelihood w.r.t. the regression coefficients
  vec grad_f(n); grad_f.zeros(); // pointwise gradient of the negative log likelihood w.r.t. the latent values f
  vec f = x*beta;
  
  mat xw(n,D); // this will be the weighted x
  mat regmat = lambda*diagmat(penalty.elem(act)); // regularization matrix (active variables)
  
  // initial quadratic approximation
  vec z(n);
  double loss_obs;
  pseudo_obs.eval(f, w, z, w, grad_f, loss_obs);
  vec z_ls(n), w_ls(n), grad_ls(n); // pseudo-observations in the line search
  double loss_initial = loss_obs + elnet_penalty(beta, lambda, 0, penalty);
  double loss_old = loss_initial; // will be updated iteratively
  loss = loss_initial; // will be updated iteratively
  double tol = thresh*fabs(loss_initial); // threshold for convergence
  double decrement = 0; // newton decrement, used to monitor convergence
  
  qau = 0;
  while (qau < qa_updates_max) {
    
    // weight the observations
    for (j=0; j<D; ++j)
      xw.col(j) = x.col(j) % sqrt(w);
    
    // weighted least squares solution
    // (with constraints: a bounded weighted least squares solution, starting
    // from the current (feasible) solution; the line search below then stays
    // feasible)
    if (inact.n_elem == 0) {
      if (n_cons == 0) {
        beta_new = solve( xw.t()*xw + regmat, xw.t()*(sqrt(w)%z) );
      } else {
        beta_new = beta;
        bls_ok = bls_active_set(beta_new, xw.t()*xw + regmat,
                                xw.t()*(sqrt(w)%z), cons) && bls_ok;
      }
    } else {
      beta_new.zeros();
      if (act.n_elem > 0) {
        mat xw_act = xw.cols(act);
        if (cons_act.n_elem == 0) {
          beta_new.elem(act) = solve( xw_act.t()*xw_act + regmat,
                                      xw_act.t()*(sqrt(w)%z) );
        } else {
          vec beta_act = beta.elem(act);
          bls_ok = bls_active_set(beta_act, xw_act.t()*xw_act + regmat,
                                  xw_act.t()*(sqrt(w)%z), cons_act) && bls_ok;
          beta_new.elem(act) = beta_act;
        }
      }
    }
    dbeta = beta_new - beta;
    grad = x.t()*grad_f + lambda*penalty_grad%beta; // gradient of negative log likelihood + gradient of penalty
    decrement = -sum(grad%dbeta); // newton decrement
    
    // check for convergence (with tol = 0: stop only at the solution, up to
    // rounding)
    if (decrement <= tol)
    	break;
    
    // backtracking line search
    t = 1.0/b;
    
    ls_iter = 0;
    
    while (ls_iter < ls_iter_max) {
      
      t = b*t;
      f = x*(beta+t*dbeta);
      pseudo_obs.eval(f, w, z_ls, w_ls, grad_ls, loss_obs);
      loss = loss_obs + elnet_penalty(beta+t*dbeta, lambda, 0, penalty);
      ++ls_iter;
      
      if (std::isnan(loss))
        continue;
      
      if (decrement > 0) {
      	if (loss < loss_old - a*t*decrement )
      		break;
      } else {
      	std::ostringstream msg;
      	msg << "The search direction is not a descent direction ";
      	msg << "(newton decrement = " << decrement << ", should be positive), ";
      	msg << ", this is likely a bug. Please report to the developers.";
      	msgs.push_back(msg.str());
      }
    }
    
    if (ls_iter == ls_iter_max && ls_iter_max > 1) {
      // beta.print("beta = ");
      // dbeta.print("dbeta = ");
      // grad.print("grad = ");
      // Rcpp::Rcout << "loss = " << loss << "\n";
      // Rcpp::Rcout << "loss_initial = " << loss_initial << "\n";
      // Rcpp::Rcout << "tol = " << tol << "\n";
      // Rcpp::Rcout << "decrement = " << decrement << "\n";
      // Rcpp::Rcout << "\n\n";
      msgs.push_back("glm_ridge warning: maximum number of line search iterations reached. The optimization can be ill-behaved.");
      break;
    }
    
    // update the solution
    beta = beta + t*dbeta;
    z = z_ls;
    w = w_ls;
    grad_f = grad_ls;
    loss_old = loss;
    ++qau;
  }
  
  if (qau == qa_updates_max && qa_updates_max > 1) {
  	if (decrement/fabs(loss_initial) > 100*tol) {
			// warn the user if the max number of iterations is reached and we are relatively far
			// (two orders of magnitude) from the given convergence threshold 
  		std::ostringstream msg;
  		msg << "glm_ridge warning: maximum number of quadratic approximation updates reached, within ";
  		msg << decrement << " from optimum (tolerance = " << thresh << ").";
  		msgs.push_back(msg.str());
  	}
  }
  
  if (!bls_ok)
    msgs.push_back("glm_ridge warning: maximum number of iterations of the bounded least-squares step reached. Results can be inaccurate.");
  
  // undo the sign flips of the constrained coefficients
  for (int k=0; k<n_cons; ++k)
    beta(cons(k)) *= signs(k);
}



/**
 * Wrapper for glm_ridge that can be called from R.
 */
// [[Rcpp::export]]
List glm_ridge_c( arma::mat x,
                  Function pseudo_obs,
                  double lambda,
                  bool intercept,
                  arma::vec penalty, // relative penalties for the variables
                  arma::vec beta_init, // initial value for the coefficients (containing the intercept as the first element)
                  arma::vec w_init, // initial guess for the weights of the pseudo-gaussian observations (needed for Student-t model)
                  double thresh,
                  int qa_updates_max,
                  int ls_iter_max=100,
                  bool debug=false,
                  Rcpp::IntegerVector cons_cols=Rcpp::IntegerVector::create(), // 1-based columns of x with a sign constraint
                  Rcpp::NumericVector cons_signs=Rcpp::NumericVector::create()) // their signs (1 or -1)
{
  int D = x.n_cols;
  if (cons_cols.size() != cons_signs.size())
    Rcpp::stop("cons_cols and cons_signs must have the same length.");
  int n_cons = cons_cols.size();
  uvec cons(n_cons);
  vec signs(n_cons);
  uvec seen(D); seen.zeros();
  for (int k=0; k<n_cons; ++k) {
    if (cons_cols[k] < 1 || cons_cols[k] > D || seen(cons_cols[k] - 1) == 1)
      Rcpp::stop("cons_cols must contain unique column indices of x.");
    seen(cons_cols[k] - 1) = 1;
    if (cons_signs[k] != 1.0 && cons_signs[k] != -1.0)
      Rcpp::stop("cons_signs must be 1 or -1.");
    // indices into the coefficient vector (which starts with the intercept, if any)
    cons(k) = cons_cols[k] - 1 + (intercept ? 1 : 0);
    signs(k) = cons_signs[k];
  }
  if (intercept)
    D++;
  
  vec beta = beta_init;
  vec w = w_init;
  int qau;
  double loss;
  std::vector<std::string> msgs;
  PseudoObsR po(pseudo_obs);
  glm_ridge(beta, loss, w, qau, msgs, x, po, lambda, intercept, penalty, thresh, qa_updates_max, ls_iter_max, debug, cons, signs);
    
  if (intercept) 
    return List::create(vec(beta.tail(D-1)), beta(0), w, loss, qau, wrap(msgs));
    // return List::create(vec(beta.tail(D-1)), beta(0), w, qau);
  else 
    return List::create(beta, 0.0, w, loss, qau, wrap(msgs));
    // return List::create(beta, 0.0, w, qau);
}


/**
 * As glm_ridge_c, but with the pseudo-observations computed in C++ for the
 * family and link given by their codes (see native_pseudo_codes() in R).
 */
// [[Rcpp::export]]
List glm_ridge_native_c( arma::mat x,
                         arma::vec y,
                         arma::vec weights,
                         arma::vec offset,
                         int family_code,
                         int link_code,
                         double lambda,
                         bool intercept,
                         arma::vec penalty,
                         arma::vec beta_init,
                         arma::vec w_init,
                         double thresh,
                         int qa_updates_max,
                         int ls_iter_max=100,
                         bool debug=false)
{
  int D = x.n_cols;
  if (intercept)
    D++;
  
  vec beta = beta_init;
  vec w = w_init;
  int qau;
  double loss;
  std::vector<std::string> msgs;
  PseudoObsNative po(y, weights, offset, family_code, link_code);
  glm_ridge(beta, loss, w, qau, msgs, x, po, lambda, intercept, penalty, thresh, qa_updates_max, ls_iter_max, debug);
  
  if (intercept) 
    return List::create(vec(beta.tail(D-1)), beta(0), w, loss, qau, wrap(msgs));
  else 
    return List::create(beta, 0.0, w, loss, qau, wrap(msgs));
}


/**
 * glm_ridge_native_c for several responses (the columns of y) with the same
 * model matrix, observation weights, offsets and penalty, as for the draw-wise
 * submodels of one projection. Each column is fitted separately, from the same
 * initial values.
 */
// [[Rcpp::export]]
List glm_ridge_native_multi_c( arma::mat x,
                               arma::mat y,
                               arma::vec weights,
                               arma::vec offset,
                               int family_code,
                               int link_code,
                               double lambda,
                               bool intercept,
                               arma::vec penalty,
                               arma::vec beta_init,
                               arma::vec w_init,
                               double thresh,
                               int qa_updates_max,
                               int ls_iter_max=100,
                               bool debug=false)
{
  int D = x.n_cols;
  if (intercept)
    D++;
  int S = y.n_cols;
  int n = y.n_rows;
  
  mat beta_out(intercept ? D-1 : D, S);
  vec beta0_out(S);
  mat w_out(n, S);
  vec loss_out(S);
  IntegerVector qau_out(S);
  List msgs_out(S);
  for (int s = 0; s < S; ++s) {
    vec beta = beta_init;
    vec w = w_init;
    int qau;
    double loss;
    std::vector<std::string> msgs;
    PseudoObsNative po(y.col(s), weights, offset, family_code, link_code);
    glm_ridge(beta, loss, w, qau, msgs, x, po, lambda, intercept, penalty, thresh, qa_updates_max, ls_iter_max, debug);
    if (intercept) {
      beta_out.col(s) = beta.tail(D-1);
      beta0_out(s) = beta(0);
    } else {
      beta_out.col(s) = beta;
      beta0_out(s) = 0.0;
    }
    w_out.col(s) = w;
    loss_out(s) = loss;
    qau_out[s] = qau;
    msgs_out[s] = wrap(msgs);
  }
  return List::create(Named("beta") = beta_out, Named("beta0") = beta0_out,
                      Named("w") = w_out, Named("loss") = loss_out,
                      Named("qa_updates") = qau_out, Named("msgs") = msgs_out);
}
