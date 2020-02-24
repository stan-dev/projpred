## Function handles for the projection
##

project_submodel_poc <- function(vind, p_ref, refmodel, family_kl, intercept, regul = 1e-12) {

  mu <- p_ref$mu
  dis <- p_ref$dis

  wobs <- rep(1.0, NROW(mu))
  wsample <- rep(1.0, NCOL(mu))

  wobs <- wobs / sum(wobs)
  wsample <- wsample / sum(wsample)

  if (intercept) {
    vind <- c("1", vind)
  }

  pobs <- pseudo_data(0, mu, family_kl, weights=wobs)
  form <- refmodel$formula
  subset <- subset_formula_and_data(unique(unlist(vind)), refmodel$fetch_data(), y=pobs$z)
  capture.output(proj_refit <- refmodel$mle(flatten_formula(subset$formula), subset$data, regul=regul),
                 type = "message")
  musub <- family_kl$mu_fun(proj_refit)

  if (family_kl$family == "gaussian")
    ref <- list(mu=pobs$z, var=p_ref$var, w=pobs$w)
  else {
    ref <- p_ref
    ref$w <- rep(0, NROW(mu))
  }

  dis_sub <- family_kl$dis_fun(ref, list(mu=musub), ref$w)
  kl <- family_kl$kl(ref, list(weights=wobs), list(mu=musub, dis=dis_sub))
  submodel <- list(kl = kl, dis = dis_sub, weights = wsample)

  if(length(vind) == 1 && intercept) {
    submodel$vind <- integer(length=0)
  } else {
    submodel$vind <- vind[(1+intercept*1):length(vind)]
  }
  submodel$intercept <- intercept
  submodel$sub_fit <- proj_refit
  return(submodel)
}


## function handle for the projection over samples
.get_proj_handle_poc <- function(family_kl, regul=1e-9) {
  return(function(vind, p_ref, refmodel, intercept)
    project_submodel_poc(vind, p_ref, refmodel, family_kl, intercept, regul=regul))
}

.get_submodels_poc <- function(searchpath, nv, family_kl, p_ref,
                               refmodel, intercept, regul, cv_search=FALSE) {
  ##
  ##
  ## Project onto given model sizes nv. Returns a list of submodels. If cv_search=FALSE,
  ## submodels parameters will be as they were computed during the search, so there is
  ## no need to project anything anymore, and this function simply fetches the information
  ## from the searchpath list, which contains the parameter values.
  ##

  varorder <- searchpath$vind
  p_sel <- searchpath$p_sel

  if (!cv_search) {
    ## simply fetch the already computed quantities for each submodel size
    fetch_submodel <- function(nv) {
      submodel <- list()
      vind <- utils::head(varorder, nv)

      mu_ref <- p_sel$mu
      wobs <- rep(1.0, NROW(mu_ref))
      wsample <- rep(1.0, NCOL(mu_ref))

      wobs <- wobs / sum(wobs)
      wsample <- wsample / sum(wsample)

      pobs <- pseudo_data(0, mu_ref, family_kl, weights=wobs)

      if (nv == 0) {
        vind <- c("1")
        form <- refmodel$formula
        subset <- subset_formula_and_data(unique(unlist(vind)), refmodel$fetch_data(), y=pobs$z)
        sub_refit <- refmodel$mle(flatten_formula(subset$formula), subset$data, regul=regul)
      } else {
        ## reuse sub_fit as projected during search
        sub_refit <- searchpath$sub_fits[[nv]]
      }

      ## split b to alpha and beta, add it to submodel and return the result
      if (family_kl$family == "gaussian")
        ref <- list(mu=pobs$z, var=p_sel$var, w=pobs$w)
      else {
        ref <- p_sel
        ref$w <- rep(0, NROW(mu_ref))
      }

      mu <- family_kl$mu_fun(sub_refit)
      submodel$dis <- family_kl$dis_fun(ref, list(mu=mu), ref$w)
      submodel$kl <- family_kl$kl(ref, list(weights=wobs), list(mu=mu, dis=submodel$dis))
      submodel$weights <- wsample
      submodel$vind <- vind
      submodel$intercept <- intercept
      submodel$sub_fit <- sub_refit
      return(submodel)
    }
  } else {
    ## need to project again for each submodel size
    projfun <- .get_proj_handle_poc(family_kl, regul)
    fetch_submodel <- function(nv) {
      if (nv == 0)
        ## empty
        vind <- integer(length=0)
      else
        vind <- varorder[1:nv]
      return(projfun(vind, p_ref, refmodel, intercept))
    }
  }
  submodels <- lapply(nv, fetch_submodel)
  return(submodels)
}
