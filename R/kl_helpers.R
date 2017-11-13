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

  # kl-divergences
  # for binomial and poisson it is the mean of the dev.resids divided by 2
  kl_dev <- function(p_full, data, p_sub) {
    if(NCOL(p_full$mu)>1) {
      w <- rep(data$weights, NCOL(p_full$mu))
      colMeans(fam$dev.resids(p_full$mu, p_sub$mu, w))/2
    } else {
      mean(fam$dev.resids(p_full$mu, p_sub$mu, data$weights))/2
    }
  }
  kl_gauss <- function(p_full, data, p_sub) log(p_sub$dis) - log(p_full$dis)
  kl_gamma <- function(p_full, data, p_sub) {
    mean(data$weights*(
      p_sub$dis*(log(p_full$dis)-log(p_sub$dis)+log(p_sub$mu)-log(p_full$mu)) +
        digamma(p_full$dis)*(p_full$dis - p_sub$dis) - lgamma(p_full$dis) +
        lgamma(p_sub$dis) + p_full$mu*p_sub$dis/p_sub$mu - p_full$dis))
  }

  # dispersion parameters in one-to-one projection.
  # for gaussian dispersion is sigma and for gamma it is the shape parameter.
  # NOTE: these could be functions of p_full = list(mu, var), p_sub = list(mu), wobs
  dis_na <- function(mu,v,musub,wobs) rep(0, ncol(mu)) 
  dis_gauss <- function(mu,v,musub,wobs) {
  	wobs <- wobs / sum(wobs)
  	sqrt(colSums(wobs*(v + (mu-musub)^2)))
  }
  dis_student_t <- function(mu,v,musub,wobs) { stop('Projection of dispersion not yet implemented for student-t') }
  dis_gamma <- function(mu,v,musub,wobs) {
      # TODO, IMPLEMENT THIS
      stop('Projection of dispersion parameter not yet implemented for family Gamma.')
      #mean(data$weights*((p_full$mu - p_sub$mu)/
      #                      fam$mu.eta(fam$linkfun(p_sub$mu))^2))
  }
  
  # # dispersion parameters for a given cluster in the sample clustering
  # discl_na <- function(mu, dis, wobs, wsample) { 1 }
  # discl_gauss <- function(mu, dis, wobs, wsample) {
  #     mu_mean <- mu %*% wsample
  #     mu_var <- mu^2 %*% wsample - mu_mean^2
  #     sqrt( sum(wsample*dis^2) + mean(wobs*mu_var) )
  # }
  # discl_gamma <- function(mu, dis, wobs, wsample) {
  #     stop('Projection of dispersion parameter not yet implemented for family Gamma.')
  # }
  
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
  	stop('not implemented for student_t yet.')
  	# mu_mean <- mu %*% wsample
  	# mu_var <- mu^2 %*% wsample - mu_mean^2
  	# sum(wsample*dis^2) + mu_var
  }
  predvar_gamma <- function(mu, dis, wsample) { stop('Family Gamma not implemented yet.')}

  # log likelihoods
  ll_binom <- function(mu, dis, y, weights=1) dbinom(weights*y, weights, mu, log=T)
  ll_poiss <- function(mu, dis, y, weights=1) weights*dpois(y, mu, log=T)
  ll_gauss <- function(mu, dis, y, weights=1) {
    dis <- matrix(rep(dis, each=length(y)), ncol=NCOL(mu))
    weights*dnorm(y, mu, dis, log=T)
  }
  ll_student_t <- function(mu, dis, y, weights=1) {
    dis <- matrix(rep(dis, each=length(y)), ncol=NCOL(mu))
    weights*(dt((y-mu)/dis, fam$nu, log=T) - log(dis))
    # weights*dnorm(y, mu, dis, log=T)
  }
  ll_gamma <- function(mu, dis, y, weights=1) {
    dis <- matrix(rep(dis, each=length(y)), ncol=NCOL(mu))
    weights*dgamma(y, dis, dis/matrix(mu), log=T)
  }

  # functions to sample from posterior predictive distribution
  ppd_gauss <- function(mu, dis, weights = 1) rnorm(length(mu), mu, dis)
  ppd_binom <- function(mu, dis, weights = 1) rbinom(length(mu), weights, mu)
  ppd_poiss <- function(mu, dis, weights = 1) rpois(length(mu), mu)
  ppd_gamma <- function(mu, dis, weights = 1) rgamma(length(mu), dis, dis/mu)


  # function for computing mu = E(y)
  mu_fun <- function(x, alpha, beta, offset) {
    if (!is.matrix(x)) stop('x must be a matrix.')
    if (!is.matrix(beta)) stop('beta must be a matrix')
    fam$linkinv(cbind(1, x) %*% rbind(alpha, beta) + offset)
  }

  # return the family object with the correct function handles
  c(switch(fam$family,
           'binomial' = list(kl = kl_dev, ll_fun = ll_binom, dis_fun = dis_na, #discl_fun = discl_na,
                             predvar = predvar_na, ppd_fun = ppd_binom),
           'poisson' = list(kl = kl_dev, ll_fun = ll_poiss, dis_fun = dis_na, #discl_fun = discl_na,
                            predvar = predvar_na, ppd_fun = ppd_poiss),
           'gaussian' = list(kl = kl_gauss, ll_fun = ll_gauss, dis_fun = dis_gauss, #discl_fun = discl_gauss,
                             predvar = predvar_gauss, ppd_fun = ppd_gauss),
           'Gamma' = list(kl = kl_gamma, ll_fun = ll_gamma, dis_fun = dis_gamma, #discl_fun = discl_gamma,
                          predvar_gamma, ppd_fun = ppd_gamma),
  				 'Student_t' = list(kl = NULL, ll_fun = ll_student_t, dis_fun = dis_student_t, #discl_fun = discl_gauss,
  				 									predvar = predvar_student_t, ppd_fun = NULL)
  				 ),
    list(mu_fun = mu_fun), fam)

}






# define a student-t family object. Dispersion is defined to be the scale parameter
# of the distribution
Student_t <- function(link='identity', nu=1) {
	
	if (link != 'identity')
		stop('Only identity link supported currently.')
	if (!is.character(link))
		stop('Link must be a string.')
	
	# fetch the link statistics
	stats <- make.link(link)
	
	# variance function # CHECK THIS!!!! 
	varfun <- function(mu) {
		if (nu > 2)
			rep(nu/(nu-2), length(mu))
		else
			rep(Inf, length(mu))
	}
	
	# create the object and append the relevant fields
	fam <- list(
		family = 'Student_t',
		nu = nu,
		link = link,
		linkfun = stats$linkfun,
		linkinv = stats$linkinv,
		variance = varfun, 
		dev.resids = function(y, mu, wt, dis=1) (nu+1) * log(1 + 1/nu*((y-mu)/dis)^2), 
		aic = function(y, n, mu, wt, dev) stop('aic not implemented yet.'),
		mu.eta = stats$mu.eta,
		initialize = expression({ stop('initialization not implemented yet.')	}),
		validmu = function(mu) TRUE,
		valideta = stats$valideta
	)
	
	structure(fam, class = 'family')
}
















