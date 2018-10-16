# Model-specific helper functions.
#
# \code{kl_helpers(fam)} returns a family object augmented with auxiliary functions that
# are needed for computing KL divergence, log predictive density, projecting dispersion etc.
#
# Missing: Quasi-families not implemented. If dis_gamma is the correct shape
# parameter for projected Gamma regression, everything should be OK for gamma.

kl_helpers <- function(fam) {

  # define the functions for all families but
  # return only the ones that are needed.
	
	if (.has.fam.extras(fam))
		# if the object already was created using this function, then return
		return(fam)

  # kl-divergences
  # for binomial and poisson it is the mean of the dev.resids divided by 2
	# NOTE: we should get rid off these, they are not much of a help..
  kl_dev <- function(pref, data, psub) {
    if(NCOL(pref$mu)>1) {
      w <- rep(data$weights, NCOL(pref$mu))
      colMeans(fam$dev.resids(pref$mu, psub$mu, w))/2
    } else {
      mean(fam$dev.resids(pref$mu, psub$mu, data$weights))/2
    }
  }
  kl_gauss <- function(pref, data, psub) colSums(data$weights*(psub$mu-pref$mu)^2) # not the actual kl but reasonable surrogate..
  kl_student_t <- function(pref, data, psub) log(psub$dis) #- 0.5*log(pref$var) # FIX THIS, NOT CORRECT
  kl_gamma <- function(pref, data, psub) {
  	stop('KL-divergence for gamma not implemented yet.')
    # mean(data$weights*(
    #   p_sub$dis*(log(pref$dis)-log(p_sub$dis)+log(psub$mu)-log(pref$mu)) +
    #     digamma(pref$dis)*(pref$dis - p_sub$dis) - lgamma(pref$dis) +
    #     lgamma(p_sub$dis) + pref$mu*p_sub$dis/p_sub$mu - pref$dis))
  }

  # dispersion parameters in draw-by-draw or clustered projection.
  # for gaussian and student-t dispersion is the noise scale, and for gamma it is the shape parameter.
  # in both cases pref is a list with field mu and var giving the mean and predictive
  # variance for each draw/cluster (columns) and each observation (rows).
  # psub is a list containing mu (analogous to pref$mu) and w, which give the weights
  # of thee pseudo-observations at optimal coefficients (needed for student-t projection).
  # wobs denote the observation weights. 
  dis_na <- function(pref, psub, wobs) rep(0, ncol(pref$mu)) 
  dis_gauss <- function(pref, psub, wobs) {
  	sqrt(colSums(wobs/sum(wobs)*(pref$var + (pref$mu-psub$mu)^2)))
  }
  dis_student_t <- function(pref, psub, wobs) { 
  	s2 <- colSums( psub$w/sum(wobs)*(pref$var+(pref$mu-psub$mu)^2) ) # CHECK THIS
  	sqrt(s2)
  	# stop('Projection of dispersion not yet implemented for student-t')
  }
  dis_gamma <- function(pref, psub, wobs) {
      # TODO, IMPLEMENT THIS
      stop('Projection of dispersion parameter not yet implemented for family Gamma.')
      #mean(data$weights*((pref$mu - p_sub$mu)/
      #                      fam$mu.eta(fam$linkfun(p_sub$mu))^2))
  }
  
  
  # functions for computing the predictive variance (taking into account 
  # the uncertainty in mu)
  predvar_na <- function(mu, dis, wsample=1) { 0 }
  predvar_gauss <- function(mu, dis, wsample=1) { 
  	wsample <- wsample/sum(wsample)
  	mu_mean <- mu %*% wsample
  	mu_var <- mu^2 %*% wsample - mu_mean^2
  	as.vector( sum(wsample*dis^2) + mu_var )
  }
  predvar_student_t <- function(mu, dis, wsample=1) { 
  	wsample <- wsample/sum(wsample)
  	mu_mean <- mu %*% wsample
  	mu_var <- mu^2 %*% wsample - mu_mean^2
  	as.vector( fam$nu/(fam$nu-2)*sum(wsample*dis^2) + mu_var )
  }
  predvar_gamma <- function(mu, dis, wsample) { stop('Family Gamma not implemented yet.')}

  # log likelihoods
  ll_binom <- function(mu, dis, y, weights=1) dbinom(y*weights, weights, mu, log=T)
  ll_poiss <- function(mu, dis, y, weights=1) weights*dpois(y, mu, log=T)
  ll_gauss <- function(mu, dis, y, weights=1) {
    dis <- matrix(rep(dis, each=length(y)), ncol=NCOL(mu))
    weights*dnorm(y, mu, dis, log=T)
  }
  ll_student_t <- function(mu, dis, y, weights=1) {
    dis <- matrix(rep(dis, each=length(y)), ncol=NCOL(mu))
    if (NCOL(y) < NCOL(mu))
    	y <- matrix(y, nrow=length(y), ncol=NCOL(mu))
    weights*(dt((y-mu)/dis, fam$nu, log=T) - log(dis))
  }
  ll_gamma <- function(mu, dis, y, weights=1) {
    dis <- matrix(rep(dis, each=length(y)), ncol=NCOL(mu))
    weights*dgamma(y, dis, dis/matrix(mu), log=T)
  }
  
  # loss functions for projection. these are defined to be -2*log-likelihood, ignoring any additional constants.
  # we need to define these separately from the log-likelihoods because some of the  log-likelihoods
  # or deviance functions do not work when given the fit of the reference model (float) in place of y (integer),
  # for instance binomial and poisson models.
  dev_binom <- function(mu, y, weights=1, dis=NULL) {
  	if (NCOL(y) < NCOL(mu))
  		y <- matrix(y, nrow=length(y), ncol=NCOL(mu))
  	-2*weights*(y*log(mu) + (1-y)*log(1-mu))
  }
  dev_poiss <- function(mu, y, weights=1, dis=NULL) {
  	if (NCOL(y) < NCOL(mu))
  		y <- matrix(y, nrow=length(y), ncol=NCOL(mu))
  	-2*weights*(y*log(mu) - mu)
  } 
  dev_gauss <- function(mu, y, weights=1, dis=NULL) {
  	if (is.null(dis))
  		dis <- 1
  	else
  		dis <- matrix(rep(dis, each=length(y)), ncol=NCOL(mu))
  	if (NCOL(y) < NCOL(mu))
  		y <- matrix(y, nrow=length(y), ncol=NCOL(mu))
  	-2*weights*(-0.5/dis*(y-mu)^2 - log(dis))
  }
  dev_student_t <- function(mu, y, weights=1, dis=NULL) {
  	if (is.null(dis))
  		dis <- 1
  	else
  		dis <- matrix(rep(dis, each=length(y)), ncol=NCOL(mu))
  	if (NCOL(y) < NCOL(mu))
  		y <- matrix(y, nrow=length(y), ncol=NCOL(mu))
  	-2*weights*(-0.5*(fam$nu+1)*log(1 + 1/fam$nu*((y-mu)/dis)^2) - log(dis))
  }
  dev_gamma <- function(mu, dis, y, weights=1) {
  	# dis <- matrix(rep(dis, each=length(y)), ncol=NCOL(mu))
  	# weights*dgamma(y, dis, dis/matrix(mu), log=T)
  	stop('Loss function not implemented for Gamma-family yet.')
  }

  # functions to sample from posterior predictive distribution
  ppd_gauss <- function(mu, dis, weights = 1) rnorm(length(mu), mu, dis)
  ppd_binom <- function(mu, dis, weights = 1) rbinom(length(mu), weights, mu)
  ppd_poiss <- function(mu, dis, weights = 1) rpois(length(mu), mu)
  ppd_student_t <- function(mu, dis, weights = 1) rt(length(mu), fam$nu)*dis + mu
  ppd_gamma <- function(mu, dis, weights = 1) rgamma(length(mu), dis, dis/mu)


  # function for computing mu = E(y)
  mu_fun <- function(x, alpha, beta, offset) {
    if (!is.matrix(x)) stop('x must be a matrix.')
    if (!is.matrix(beta)) stop('beta must be a matrix')
    fam$linkinv(cbind(1, x) %*% rbind(alpha, beta) + offset)
  }

  # return the family object with the correct function handles
  c(switch(fam$family,
           'binomial' = list(kl = kl_dev, ll_fun = ll_binom, deviance = dev_binom, dis_fun = dis_na,
                             predvar = predvar_na, ppd_fun = ppd_binom),
           'poisson' = list(kl = kl_dev, ll_fun = ll_poiss, deviance = dev_poiss, dis_fun = dis_na, 
                            predvar = predvar_na, ppd_fun = ppd_poiss),
           'gaussian' = list(kl = kl_gauss, ll_fun = ll_gauss, deviance = dev_gauss, dis_fun = dis_gauss,
                             predvar = predvar_gauss, ppd_fun = ppd_gauss),
           'Gamma' = list(kl = kl_gamma, ll_fun = ll_gamma, deviance = dev_gamma, dis_fun = dis_gamma, 
                          predvar_gamma, ppd_fun = ppd_gamma),
  				 'Student_t' = list(kl = kl_student_t, ll_fun = ll_student_t, deviance = dev_student_t, dis_fun = dis_student_t,
  				 									predvar = predvar_student_t, ppd_fun = ppd_student_t)
  				 ),
    list(mu_fun = mu_fun), fam)

}




.has.dispersion <- function(fam) {
	# a function for checking whether the family has a dispersion parameter
	fam$family %in% c('gaussian','Student_t','Gamma')
}

.has.fam.extras <- function(fam) {
  # check whether the family object has the extra functions, that is, whether it was
  # created by kl_helpers
  !is.null(fam$deviance)
}


















