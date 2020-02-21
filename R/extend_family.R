# Model-specific helper functions.
#
# \code{extend_family(fam)} returns a family object augmented with auxiliary functions that
# are needed for computing KL divergence, log predictive density, projecting dispersion etc.
#
# Missing: Quasi-families not implemented. If dis_gamma is the correct shape
# parameter for projected Gamma regression, everything should be OK for gamma.

#' @export
extend_family <- function(fam) {
  if (.has_fam_extras(fam))
    ## if the object already was created using this function, then return
    return(fam)

  extend_family_specific <- match.fun(paste0("extend_family", fam))
  extend_family_specific(fam)
}

extend_family_binomial <- function(fam) {
  kl_dev <- function(pref, data, psub) {
    if(NCOL(pref$mu) > 1) {
      w <- rep(data$weights, NCOL(pref$mu))
      colMeans(fam$dev.resids(pref$mu, psub$mu, w))/2
    } else {
      mean(fam$dev.resids(pref$mu, psub$mu, data$weights))/2
    }
  }
  dis_na <- function(pref, psub, wobs=1) rep(0, ncol(pref$mu))
  predvar_na <- function(mu, dis, wsample=1) { 0 }
  ll_binom <- function(mu, dis, y, weights=1) dbinom(y*weights, weights, mu, log=T)
  dev_binom <- function(mu, y, weights=1, dis=NULL) {
  	if (NCOL(y) < NCOL(mu))
  		y <- matrix(y, nrow=length(y), ncol=NCOL(mu))
  	-2*weights*(y*log(mu) + (1-y)*log(1-mu))
  }
  ppd_binom <- function(mu, dis, weights = 1) rbinom(length(mu), weights, mu)

  fam$kl <- kl_dev
  fam$dis_fun <- dis_na
  fam$predvar <- predvar_na
  fam$ll_fun <- ll_binom
  fam$deviance <- dev_binom
  fam$ppd <- ppd_binom

  return(fam)
}

extend_family_poisson <- function(fam) {
  kl_dev <- function(pref, data, psub) {
    if(NCOL(pref$mu) > 1) {
      w <- rep(data$weights, NCOL(pref$mu))
      colMeans(fam$dev.resids(pref$mu, psub$mu, w))/2
    } else {
      mean(fam$dev.resids(pref$mu, psub$mu, data$weights))/2
    }
  }
  dis_na <- function(pref, psub, wobs=1) rep(0, ncol(pref$mu))
  predvar_na <- function(mu, dis, wsample=1) { 0 }
  ll_poiss <- function(mu, dis, y, weights=1) weights*dpois(y, mu, log=T)
  dev_poiss <- function(mu, y, weights=1, dis=NULL) {
  	if (NCOL(y) < NCOL(mu))
  		y <- matrix(y, nrow=length(y), ncol=NCOL(mu))
  	-2*weights*(y*log(y / mu) - (y - mu))
  }
  ppd_poiss <- function(mu, dis, weights = 1) rpois(length(mu), mu)

  fam$kl <- kl_dev
  fam$dis_fun <- dis_na
  fam$predvar <- predvar_na
  fam$ll_fun <- ll_poiss
  fam$deviance <- dev_poiss
  fam$ppd <- ppd_poiss

  return(fam)
}

extend_family_gaussian <- function(fam) {
  kl_gauss <- function(pref, data, psub) colSums(data$weights * (psub$mu-pref$mu)^2) # not the actual kl but reasonable surrogate..
  dis_gauss <- function(pref, psub, wobs=1) {
  	sqrt(colSums(wobs/sum(wobs)*(pref$var + (pref$mu-psub$mu)^2)))
  }
  predvar_gauss <- function(mu, dis, wsample=1) {
  	wsample <- wsample/sum(wsample)
  	mu_mean <- mu %*% wsample
  	mu_var <- mu^2 %*% wsample - mu_mean^2
  	as.vector( sum(wsample*dis^2) + mu_var )
  }
  ll_gauss <- function(mu, dis, y, weights=1) {
    dis <- matrix(rep(dis, each=length(y)), ncol=NCOL(mu))
    weights*dnorm(y, mu, dis, log=T)
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
  ppd_gauss <- function(mu, dis, weights = 1) rnorm(length(mu), mu, dis)

  fam$kl <- kl_gauss
  fam$dis_fun <- dis_gauss
  fam$predvar <- predvar_gauss
  fam$ll_fun <- ll_gauss
  fam$deviance <- dev_gauss
  fam$ppd <- ppd_gauss

  return(fam)
}

extend_family_gamma <- function(fam) {
  kl_gamma <- function(pref, data, psub) {
    stop('KL-divergence for gamma not implemented yet.')
    ## mean(data$weights*(
    ##   p_sub$dis*(log(pref$dis)-log(p_sub$dis)+log(psub$mu)-log(pref$mu)) +
    ##     digamma(pref$dis)*(pref$dis - p_sub$dis) - lgamma(pref$dis) +
    ##     lgamma(p_sub$dis) + pref$mu*p_sub$dis/p_sub$mu - pref$dis))
  }
  dis_gamma <- function(pref, psub, wobs=1) {
    ## TODO, IMPLEMENT THIS
    stop('Projection of dispersion parameter not yet implemented for family Gamma.')
    ## mean(data$weights*((pref$mu - p_sub$mu)/
    ##                      fam$mu.eta(fam$linkfun(p_sub$mu))^2))
  }
  predvar_gamma <- function(mu, dis, wsample=1) { stop('Family Gamma not implemented yet.')}
  ll_gamma <- function(mu, dis, y, weights=1) {
    dis <- matrix(rep(dis, each=length(y)), ncol=NCOL(mu))
    weights*dgamma(y, dis, dis/matrix(mu), log=T)
  }
  dev_gamma <- function(mu, dis, y, weights=1) {
    ## dis <- matrix(rep(dis, each=length(y)), ncol=NCOL(mu))
    ## weights*dgamma(y, dis, dis/matrix(mu), log=T)
    stop('Loss function not implemented for Gamma-family yet.')
  }
  ppd_gamma <- function(mu, dis, weights = 1) rgamma(length(mu), dis, dis/mu)

  fam$kl <- kl_gamma
  fam$dis_fun <- dis_gamma
  fam$predvar <- predvar_gamma
  fam$ll_fun <- ll_gamma
  fam$deviance <- dev_gamma
  fam$ppd <- ppd_gamma

  return(fam)
}

extend_family_student_t <- function(fam) {
  kl_student_t <- function(pref, data, psub) log(psub$dis) #- 0.5*log(pref$var) # FIX THIS, NOT CORRECT
  dis_student_t <- function(pref, psub, wobs=1) {
  	s2 <- colSums( psub$w/sum(wobs)*(pref$var+(pref$mu-psub$mu)^2) ) # CHECK THIS
  	sqrt(s2)
    ## stop('Projection of dispersion not yet implemented for student-t')
  }
  predvar_student_t <- function(mu, dis, wsample=1) {
  	wsample <- wsample/sum(wsample)
  	mu_mean <- mu %*% wsample
  	mu_var <- mu^2 %*% wsample - mu_mean^2
  	as.vector( fam$nu/(fam$nu-2)*sum(wsample*dis^2) + mu_var )
  }
  ll_student_t <- function(mu, dis, y, weights=1) {
    dis <- matrix(rep(dis, each=length(y)), ncol=NCOL(mu))
    if (NCOL(y) < NCOL(mu))
    	y <- matrix(y, nrow=length(y), ncol=NCOL(mu))
    weights*(dt((y-mu)/dis, fam$nu, log=T) - log(dis))
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
  ppd_student_t <- function(mu, dis, weights = 1) rt(length(mu), fam$nu)*dis + mu

  fam$kl <- kl_student_t
  fam$dis_fun <- dis_student_t
  fam$predvar <- predvar_student_t
  fam$ll_fun <- ll_student_t
  fam$deviance <- dev_student_t
  fam$ppd <- ppd_student_t

  return(fam)
}

.has.dispersion <- function(fam) {
	# a function for checking whether the family has a dispersion parameter
	fam$family %in% c('gaussian','Student_t','Gamma')
}

.has_fam_extras <- function(fam) {
  # check whether the family object has the extra functions, that is, whether it was
  # created by extend_family
  !is.null(fam$deviance)
}
