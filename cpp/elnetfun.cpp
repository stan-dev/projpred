#include<iostream>
#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]


using namespace Rcpp;
using namespace arma;


double loss_approx(arma::vec beta,
                   arma::vec f,
                   arma::vec z,
                   arma::vec w,
                   double lambda,
                   double alpha)
{
	double loss;
	loss = 0.5*sum(w % square(z-f)) + lambda*( 0.5*(1-alpha)*sum(square(beta)) + alpha*(sum(abs(beta))) );
	return loss;
}



void coord_descent(	arma::vec& beta,
					double& beta0,
					arma::vec& f,
					arma::mat& x, 
					arma::vec& z,
					arma::vec& w,
					double& lambda,
					double& alpha,
					std::set<int>& varind,
					std::set<int>& active_set,
					bool until_convergence,
					int& npasses,
					double thresh_abs,
					int maxiter = 1000)
{
	
	int iter = 0;
	double loss,loss_old;
	int k,j;
	double h;
	
	// initial loss
	loss_old = loss_approx(beta,f,z,w,lambda,alpha);
	
	// auxiliary that will be used later on
	double lam_alpha = lambda*alpha;
	double lam_oneminus_alpha = lambda*(1-alpha);
	
	while (iter < maxiter) {
		
		// update the intercept
		f = f - beta0;
		beta0 = sum(w % (z - f)) / sum(w);
		f = f + beta0;
		
		active_set.clear();
		
		for (std::set<int>::iterator it=varind.begin(); it!=varind.end(); ++it) {
			
			// update the regression coefficients via 'soft-thresholding'
			
			// varible index
			j = *it;
			
			f = f - beta(j)*x.col(j);
			h = sum( w % x.col(j) % (z - f) ); // auxiliary variable
			
			if (fabs(h) <= lam_alpha) {
				beta(j) = 0.0;
			} else if (h > 0) {
				beta(j) = (h - lam_alpha) / ( sum(w % square(x.col(j))) + lam_oneminus_alpha );
		        active_set.insert(j);
			} else {
				beta(j) = (h + lam_alpha) / ( sum(w % square(x.col(j))) + lam_oneminus_alpha );
		        active_set.insert(j);
			}
			f = f + beta(j)*x.col(j);
		}
		
		++iter;
		++npasses;
		loss = loss_approx(beta,f,z,w,lambda,alpha);
		
		if (until_convergence) {
			if (fabs(loss-loss_old) < thresh_abs) {
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
		std::cout << "Warning: maximum number of iterations reached in coordinate descent. Results can be inaccurate!\n";
	
}




// [[Rcpp::export]]
List glm_elnet_c(mat x,
               Function pseudo_obs,
               vec lambda,
               double alpha,
               double thresh,
               int qa_updates_max,
               int pmax,
               bool pmax_strict,
               int as_updates_max = 50)
{
	
	
    // for gaussian pseudo data
    List obs; 
    vec z; // observations
    vec w; // weights (inverse variances)
    
    int n = x.n_rows;
    int D = x.n_cols;
    int nlam = lambda.size();
    double lam; // temporary varible for fixed lambda
    int k; // lambda index
    int qau; // number of quadratic approximation updates
    int asu; // number of active set updates (for a given quadratic approximation)
    
    
    // for storing the whole solution path
    rowvec beta0_path(nlam);
    mat beta_path(D,nlam);
    beta0_path.zeros();
    beta_path.zeros();
    int npasses = 0; // counts how many times the coefficient vector is looped through
    urowvec qa_updates(nlam);
    qa_updates.zeros();
    urowvec as_updates(nlam);
    as_updates.zeros();
    
    
    // initialization
    vec beta(D);
    beta.zeros(D);
    double beta0 = 0.0;
    vec f = x*beta + beta0;
    std::set<int> active_set; 
    std::set<int> active_set_old;
    std::set<int> varind_all; // a constant set containing indices of all the variables
    for (int j=0; j<D; j++)
        varind_all.insert(j);
    
    
    obs = pseudo_obs(f);
    z = as<arma::vec>(obs["z"]);
    w = as<arma::vec>(obs["w"]);
    double loss_initial = loss_approx(beta, f, z, w, lambda(0), alpha); // initial loss
    double loss_old = loss_initial; // will be updated iteratively
    double loss; // will be updated iteratively
    
    
    // loop over lambda values
    for (k=0; k<nlam; ++k) {
        
        lam = lambda(k);
    	
        // for (qau=1; qau<=qa_updates_max; ++qau) {
        qau = 0;
        while (qau < qa_updates_max) {
            
            /** update the quadratic likelihood approximation (would be needed only
             for likelihoods other than gaussian) */
            obs = pseudo_obs(f);
            z = as<vec>(obs["z"]);
            w = as<vec>(obs["w"]);
            ++qau;
            
            // run the coordinate descent until convergence for the current
            // quadratic approximation
            asu = 0;
            while (asu < as_updates_max) {

                // iterate within the current active set until convergence (this might update the variable active_set_old,
                // if some variable goes to zero)
                coord_descent(beta, beta0, f, x, z, w, lam, alpha, active_set, active_set_old, true, npasses, 0.01*thresh*loss_initial);
                
                // perfom one pass over all the variables and check if the active set changes (this might update the
                // variable active_set)
                coord_descent(beta, beta0, f, x, z, w, lam, alpha, varind_all, active_set, false, npasses, 0.01*thresh*loss_initial);
                
                ++asu;

				if (active_set==active_set_old) {
                    // active set did not change so convergence reached
                    // (for the current quadratic approximation to the likelihood)
                    break;
                }
            }
            as_updates(k) = as_updates(k) + asu;
            
            loss = loss_approx(beta, f, z, w, lam, alpha);
            
            // check if converged
            if (fabs(loss-loss_old) < thresh * fabs(loss_initial)) {
                // convergence reached; proceed to the next lambda value
                break;
            } else {
                // continue iterating
                loss_old = loss;
            }
        }
        // store the current solution
        beta0_path(k) = beta0;
        beta_path.col(k) = beta;
        qa_updates(k) = qau;
        
        if (qau == qa_updates_max && qa_updates_max > 1)
        	std::cout << "glm_elnet warning: maximum number of quadratic approximation updates reached. Results can be inaccurate!\n";
        
        if (active_set.size() >= pmax+1 || active_set.size() == D) {
			// obtained solution with more than pmax variables (or the number of columns in x), so terminate
			if (pmax_strict) {
			    // return solutions only up to the previous lambda value
			    beta0_path = beta0_path.head(k);
			    beta_path = beta_path.head_cols(k);
			    break;
			} else {
			    // return solutions up to the current lambda value
			    beta0_path = beta0_path.head(k+1);
			    beta_path = beta_path.head_cols(k+1);
			    break;
			}
        }
    }
    
    return List::create(beta_path, beta0_path, npasses, qa_updates, as_updates);
}









// [[Rcpp::export]]
List glm_ridge_c(arma::mat x,
               Function pseudo_obs,
               double lambda,
               double thresh,
               int qa_updates_max)
{
	// for gaussian pseudo data
	List obs;
	vec z; // observations
	vec w; // weights (inverse variances)

	int n = x.n_rows;
	int D = x.n_cols;
	int alpha = 0;
	int qau;
	int j;

	// initialization
	vec beta(D+1);
	beta.zeros();

	// add a vector of ones to x
	mat xm = join_horiz(ones<vec>(n), x);
	vec f = xm*beta;

	// this will be the weighted x
	mat xmw;

	obs = pseudo_obs(f);
	z = as<vec>(obs["z"]);
	w = as<vec>(obs["w"]);
	double loss_initial = loss_approx(beta, f, z, w, lambda, alpha); // initial loss
	double loss_old = loss_initial; // will be updated iteratively
	double loss; // will be updated iteratively

	for (qau=1; qau<=qa_updates_max; ++qau) {

		/** update the quadratic likelihood approximation (would be needed only
		for likelihoods other than gaussian) */
		obs = pseudo_obs(f);
		z = as<vec>(obs["z"]);
		w = as<vec>(obs["w"]);

		// weight the observations
		xmw = xm;
		for (j=0; j<D+1; ++j)
			xmw.col(j) = xmw.col(j) % sqrt(w);
		
		// weighted least squares
		beta = solve(xmw.t()*xmw, xmw.t()*(sqrt(w)%z) );
		f = xm*beta;

		loss = loss_approx(beta, f, z, w, lambda, alpha);

		// check if converged
		if (fabs(loss-loss_old) < thresh * fabs(loss_initial)) {
			// convergence reached
			break;
		} else {
			// continue iterating
			loss_old = loss;
		}
	}
	if (qau-1 == qa_updates_max && qa_updates_max > 1)
		std::cout << "Warning: maximum number of quadratic approximation updates reached. Results can be inaccurate!\n";

	// separate the intercept and the other coefficients
	return List::create(vec(beta.tail(D)), beta(0), qau);
}













