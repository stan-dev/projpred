#include <iostream>
#include <RcppArmadillo.h>
#include <cmath>
#include <limits>
#include <string>

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





/** Internal function that gives the output in c++ way (writes into allocated memory).
 * See glm_ridge_c for a wrapper that is being called from R.
 */
void glm_ridge( vec& beta,      // output: regression coefficients (contains intercept)
                double& loss,   // output: value of the loss function
                vec& w, 				// output: weights of the pseudo-gaussian observations at the optimum (needed for Student-t model)
                int& qau,       // output: number of quadratic approximation updates 
                arma::mat x,
                Function pseudo_obs,
                double lambda,
                bool intercept,
                arma::vec penalty, // relative penalties for the variables
                double thresh,
                int qa_updates_max,
                int ls_iter_max=50,
                bool debug=false)
{
  
  if (intercept) {
    // add a vector of ones to x and set zero penalty for the intercept
    x = join_horiz(ones<vec>(x.n_rows), x);
    penalty = join_vert(zeros<vec>(1), penalty); 
  }
  
  int n = x.n_rows;
  int D = x.n_cols;
  int ls_iter; // counts linesearch iterations
  int j;
  double t; // step size in line search 
  double a = 0.1; // backtracking line search parameters a and b (see Boyd and Vandenberghe, 2004)
  double b = 0.5; 
  
  // initialization
  vec beta_new(D); beta_new.zeros();
  vec dbeta(D); dbeta.zeros();
  vec grad(D); grad.zeros(); // gradient of the negative log likelihood w.r.t. the regression coefficients
  vec grad_f(n); grad_f.zeros(); // pointwise gradient of the negative log likelihood w.r.t. the latent values f
  vec f = x*beta;
  
  mat xw(n,D); // this will be the weighted x
  mat regmat = lambda*diagmat(penalty);//eye(D,D); // regularization matrix
  
  // initial quadratic approximation
  List obs = pseudo_obs(f,w);
  vec z = as<vec>(obs["z"]);
  w = as<vec>(obs["wobs"]);
  grad_f = as<vec>(obs["grad"]);
  double loss_initial = ((double) obs["loss"]) + elnet_penalty(beta, lambda, 0, penalty);
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
    beta_new = solve( xw.t()*xw + regmat, xw.t()*(sqrt(w)%z) ); 
    dbeta = beta_new - beta;
    grad = x.t()*grad_f + lambda*penalty%beta; // gradient of negative log likelihood + gradient of penalty
    decrement = -sum(grad%dbeta); // newton decrement
    
    // check for convergence
    if (decrement < tol)
    	break;
    
    // backtracking line search
    t = 1.0/b;
    
    ls_iter = 0;
    
    while (ls_iter < ls_iter_max) {
      
      t = b*t;
      f = x*(beta+t*dbeta);
      obs = pseudo_obs(f,w);
      loss = ((double) obs["loss"]) + elnet_penalty(beta+t*dbeta, lambda, 0, penalty);
      ++ls_iter;
      
      if (std::isnan(loss))
        continue;
      
      if (decrement > 0) {
      	if (loss < loss_old - a*t*decrement )
      		break;
      } else {
      	Rcpp::Rcout << "The search direction is not a descent direction ";
      	Rcpp::Rcout << "(newton decrement = " << decrement << ", should be positive), ";
      	Rcpp::Rcout << ", this is likely a bug. Please report to the developers." << '\n';
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
      Rcpp::Rcout << "glm_ridge warning: maximum number of line search iterations reached. The optimization can be ill-behaved.\n";
      break;
    }
    
    // update the solution
    beta = beta + t*dbeta;
    z = as<vec>(obs["z"]);
    w = as<vec>(obs["wobs"]);
    grad_f = as<vec>(obs["grad"]);
    loss_old = loss;
    ++qau;
  }
  
  if (qau == qa_updates_max && qa_updates_max > 1) {
  	if (decrement/fabs(loss_initial) > 100*tol) {
			// warn the user if the max number of iterations is reached and we are relatively far
			// (two orders of magnitude) from the given convergence threshold 
  		Rcpp::Rcout << "glm_ridge warning: maximum number of quadratic approximation updates reached, within ";
  		Rcpp::Rcout << decrement << " from optimum (tolerance = " << thresh << ").\n";
  	}
  }
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
                  bool debug=false)
{
  int D = x.n_cols;
  if (intercept)
    D++;
  
  vec beta = beta_init;
  vec w = w_init;
  int qau;
  double loss;
  glm_ridge(beta, loss, w, qau, x, pseudo_obs, lambda, intercept, penalty, thresh, qa_updates_max, ls_iter_max, debug);
    
  if (intercept) 
    return List::create(vec(beta.tail(D-1)), beta(0), w, loss, qau);
    // return List::create(vec(beta.tail(D-1)), beta(0), w, qau);
  else 
    return List::create(beta, 0.0, w, loss, qau);
    // return List::create(beta, 0.0, w, qau);
}


