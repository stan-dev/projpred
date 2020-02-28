# Model-specific helper functions.
#
# \code{extend_family(family)} returns a family object augmented with auxiliary functions that
# are needed for computing KL divergence, log predictive density, projecting dispersion etc.
#
# Missing: Quasi-familyilies not implemented. If dis_gamma is the correct shape
# parameter for projected Gamma regression, everything should be OK for gamma.

#' Add extra fields to the family object.
#' @param family Family object.
#' @return Extended family object.
#' @export
extend_family <- function(family) {
  if (.has_family_extras(family)) {
    ## if the object already was created using this function, then return
    return(family)
  }
  extend_family_specific <- paste0("extend_family_", family$family)
  extend_family_specific <- get(extend_family_specific, mode = "function")
  extend_family_specific(family)
}

extend_family_binomial <- function(family) {
  kl_dev <- function(pref, data, psub) {
    if(NCOL(pref$mu) > 1) {
      w <- rep(data$weights, NCOL(pref$mu))
      colMeans(family$dev.resids(pref$mu, psub$mu, w))/2
    } else {
      mean(family$dev.resids(pref$mu, psub$mu, data$weights))/2
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

  family$kl <- kl_dev
  family$dis_fun <- dis_na
  family$predvar <- predvar_na
  family$ll_fun <- ll_binom
  family$deviance <- dev_binom
  family$ppd <- ppd_binom

  return(family)
}

extend_family_poisson <- function(family) {
  kl_dev <- function(pref, data, psub) {
    if(NCOL(pref$mu) > 1) {
      w <- rep(data$weights, NCOL(pref$mu))
      colMeans(family$dev.resids(pref$mu, psub$mu, w))/2
    } else {
      mean(family$dev.resids(pref$mu, psub$mu, data$weights))/2
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

  family$kl <- kl_dev
  family$dis_fun <- dis_na
  family$predvar <- predvar_na
  family$ll_fun <- ll_poiss
  family$deviance <- dev_poiss
  family$ppd <- ppd_poiss

  return(family)
}

extend_family_gaussian <- function(family) {
  kl_gauss <- function(pref, data, psub)
    colSums(data$weights * (psub$mu-pref$mu)^2) # not the actual kl but reasonable surrogate..
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

  family$kl <- kl_gauss
  family$dis_fun <- dis_gauss
  family$predvar <- predvar_gauss
  family$ll_fun <- ll_gauss
  family$deviance <- dev_gauss
  family$ppd <- ppd_gauss

  return(family)
}

extend_family_gamma <- function(family) {
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
    ##                      family$mu.eta(family$linkfun(p_sub$mu))^2))
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

  family$kl <- kl_gamma
  family$dis_fun <- dis_gamma
  family$predvar <- predvar_gamma
  family$ll_fun <- ll_gamma
  family$deviance <- dev_gamma
  family$ppd <- ppd_gamma

  return(family)
}

extend_family_student_t <- function(family) {
  kl_student_t <- function(pref, data, psub)
    log(psub$dis) #- 0.5*log(pref$var) # FIX THIS, NOT CORRECT
  dis_student_t <- function(pref, psub, wobs=1) {
  	s2 <- colSums( psub$w/sum(wobs)*(pref$var+(pref$mu-psub$mu)^2) ) # CHECK THIS
  	sqrt(s2)
    ## stop('Projection of dispersion not yet implemented for student-t')
  }
  predvar_student_t <- function(mu, dis, wsample=1) {
  	wsample <- wsample/sum(wsample)
  	mu_mean <- mu %*% wsample
  	mu_var <- mu^2 %*% wsample - mu_mean^2
  	as.vector( family$nu/(family$nu-2)*sum(wsample*dis^2) + mu_var )
  }
  ll_student_t <- function(mu, dis, y, weights=1) {
    dis <- matrix(rep(dis, each=length(y)), ncol=NCOL(mu))
    if (NCOL(y) < NCOL(mu))
    	y <- matrix(y, nrow=length(y), ncol=NCOL(mu))
    weights*(dt((y-mu)/dis, family$nu, log=T) - log(dis))
  }
  dev_student_t <- function(mu, y, weights=1, dis=NULL) {
  	if (is.null(dis))
  		dis <- 1
  	else
  		dis <- matrix(rep(dis, each=length(y)), ncol=NCOL(mu))
  	if (NCOL(y) < NCOL(mu))
  		y <- matrix(y, nrow=length(y), ncol=NCOL(mu))
  	-2*weights*(-0.5*(family$nu+1)*log(1 + 1/family$nu*((y-mu)/dis)^2) - log(dis))
  }
  ppd_student_t <- function(mu, dis, weights = 1) rt(length(mu), family$nu)*dis + mu

  family$kl <- kl_student_t
  family$dis_fun <- dis_student_t
  family$predvar <- predvar_student_t
  family$ll_fun <- ll_student_t
  family$deviance <- dev_student_t
  family$ppd <- ppd_student_t

  return(family)
}

.has_dispersion <- function(family) {
	# a function for checking whether the family has a dispersion parameter
	family$family %in% c('gaussian','Student_t','Gamma')
}

.has_family_extras <- function(family) {
  # check whether the family object has the extra functions, that is, whether it was
  # created by extend_family
  !is.null(family$deviance)
}
