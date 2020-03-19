## Function handles for the projection
##

project_submodel <- function(vind, p_ref, refmodel, family, intercept, regul = 1e-12) {
  mu <- p_ref$mu
  dis <- p_ref$dis

  if (is.null(refmodel$wobs)) {
    wobs <- rep(1.0, NROW(mu))
  } else {
    wobs <- refmodel$wobs
  }

  if (is.null(p_ref$weights)) {
    wsample <- rep(1.0, NCOL(mu))
  } else {
    wsample <- p_ref$weights
  }

  wobs <- wobs / sum(wobs)
  wsample <- wsample / sum(wsample)

  form <- refmodel$formula
  pobs <- pseudo_data(0, mu, family, offset = refmodel$offset, weights = wobs)

  link <- function(f, wprev = NULL) {
    pseudo_data(f, mu, family, offset = refmodel$offset, wprev = wprev)
  }
  mle <- function(formula, data, weights) {
    refmodel$mle(formula, data, weights = weights)
  }
  linear_predict <- function(fit) {
    refmodel$proj_predfun(fit)
  }
  replace_response <- get_replace_response(form, vind)

  subset <- subset_formula_and_data(form, unique(unlist(vind)),
    refmodel$fetch_data(),
    y = pobs$z
  )
  ## capture.output(proj_refit <- refmodel$mle(flatten_formula(subset$formula),
  ##                                           subset$data),
  ##                type = "message")

  ## capture.output(proj_refit <- iterative_weighted_least_squares(
  ##   flatten_formula(subset$formula), refmodel$fetch_data(), 3, link,
  ##   replace_response, wprev = wobs, mle = mle),
  ##   type = "message")
  proj_refit <- iterative_weighted_least_squares(
    flatten_formula(subset$formula), refmodel$fetch_data(), 3, link,
    replace_response,
    wprev = wobs, mle = mle, linear_predict = linear_predict
  )
  musub <- family$mu_fun(proj_refit, offset = refmodel$offset, weights = wobs)
  if (family$family == "gaussian") {
    ref <- list(mu = pobs$z, var = p_ref$var, w = pobs$w)
  } else {
    ref <- p_ref
    ref$w <- rep(0, NROW(mu))
  }

  dis_sub <- family$dis_fun(ref, list(mu = musub), ref$w)
  kl <- family$kl(ref, list(weights = wobs), list(mu = musub, dis = dis_sub))
  submodel <- list(kl = kl, dis = dis_sub, weights = wsample)

  submodel$vind <- vind
  submodel$sub_fit <- proj_refit
  return(submodel)
}

iterative_weighted_least_squares <- function(formula, data, iters, link,
                                             replace_response, wprev = NULL,
                                             mle = lm, linear_predict) {
  pobs <- link(0, wprev)
  wprev <- pobs$w
  data <- replace_response(pobs$z, data)
  old_fit <- NULL
  for (i in seq_len(iters)) {
    fit <- mle(formula, cbind(data, weights = wprev), weights = wprev)
    pobs <- link(linear_predict(fit), wprev)
    if (any(is.na(pobs$z))) {
      break
    }
    old_fit <- fit
    data <- replace_response(pobs$z, data)
    wprev <- pobs$w[seq_len(NROW(data))]
  }
  if (is.null(old_fit)) {
    return(fit)
  }
  return(old_fit)
}

## function handle for the projection over samples
.get_proj_handle <- function(family, regul = 1e-9) {
  return(function(vind, p_ref, refmodel, intercept) {
    project_submodel(vind, p_ref, refmodel, family, intercept, regul = regul)
  })
}

.get_submodels <- function(searchpath, nv, family, p_ref,
                           refmodel, intercept, regul, cv_search = FALSE) {
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

      if (is.null(refmodel$wobs)) {
        wobs <- rep(1.0, NROW(mu_ref))
      } else {
        wobs <- refmodel$wobs
      }

      if (is.null(p_sel$weights)) {
        wsample <- rep(1.0, NCOL(mu_ref))
      } else {
        wsample <- p_sel$weights
      }

      wobs <- wobs / sum(wobs)
      wsample <- wsample / sum(wsample)

      pobs <- pseudo_data(0, mu_ref, family, weights = wobs)

      ## reuse sub_fit as projected during search
      sub_refit <- searchpath$sub_fits[[nv + 1]]

      ## split b to alpha and beta, add it to submodel and return the result
      if (family$family == "gaussian") {
        ref <- list(mu = pobs$z, var = p_sel$var, w = pobs$w)
      } else {
        ref <- p_sel
        ref$w <- rep(0, NROW(mu_ref))
      }

      mu <- family$mu_fun(sub_refit, offset = refmodel$offset, weights = wobs)
      submodel$dis <- family$dis_fun(ref, list(mu = mu), ref$w)
      submodel$kl <- family$kl(
        ref, list(weights = wobs),
        list(mu = mu, dis = submodel$dis)
      )
      submodel$weights <- wsample
      submodel$vind <- vind
      submodel$sub_fit <- sub_refit
      return(submodel)
    }
  } else {
    ## need to project again for each submodel size
    projfun <- .get_proj_handle(family, regul)
    fetch_submodel <- function(nv) {
      if (nv == 0) {
        ## empty
        vind <- c("1")
      } else {
        vind <- varorder[1:nv]
      }
      return(projfun(vind, p_ref, refmodel, intercept))
    }
  }
  submodels <- lapply(nv, fetch_submodel)
  return(submodels)
}
