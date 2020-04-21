## Function handles for the projection
##

project_submodel <- function(solution_terms, p_ref, refmodel, family, intercept,
                             regul = 1e-4) {
  mu <- p_ref$mu
  dis <- p_ref$dis

  validparams <- .validate_wobs_wsample(refmodel$wobs, p_ref$weights, mu)
  wobs <- validparams$wobs
  wsample <- validparams$wsample

  form <- refmodel$formula
  pobs <- pseudo_data(
    f = 0, y = mu, family = family, offset = refmodel$offset, weights = wobs
  )

  link <- function(f, wprev = NULL) {
    pseudo_data(
      f = f, y = mu, family = family, offset = refmodel$offset, wprev = wprev
    )
  }
  div_minimizer <- function(formula, data, weights) {
    refmodel$div_minimizer(formula, data, weights = weights, regul = regul)
  }
  linear_predict <- function(fit) {
    refmodel$proj_predfun(fit)
  }
  replace_response <- get_replace_response(form, solution_terms)

  subset <- subset_formula_and_data(
    formula = form, terms_ = unique(unlist(solution_terms)),
    data = refmodel$fetch_data(), y = pobs$z
  )
  ## capture.output(sub_fit <- iterative_weighted_least_squares(
  ##   flatten_formula(subset$formula), refmodel$fetch_data(), 3, link,
  ##   replace_response, wprev = wobs, div_minimizer = div_minimizer),
  ##   type = "message")
  sub_fit <- iterative_weighted_least_squares(
    flatten_formula(subset$formula), refmodel$fetch_data(), 3, link,
    replace_response,
    wprev = wobs, div_minimizer = div_minimizer,
    linear_predict = linear_predict
  )

  return(.init_submodel(
    sub_fit = sub_fit, p_ref = p_ref, refmodel = refmodel,
    family = family, solution_terms = solution_terms, ref_mu = mu,
    weights = wobs, wsample = wsample
  ))
}

iterative_weighted_least_squares <- function(formula, data, iters, link,
                                             replace_response, wprev = NULL,
                                             div_minimizer = lm,
                                             linear_predict) {
  pobs <- link(0, wprev)
  wprev <- pobs$w
  data <- replace_response(pobs$z, data)
  old_fit <- NULL
  for (i in seq_len(iters)) {
    fit <- div_minimizer(formula, cbind(data, weights = wprev), weights = wprev)
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
  return(function(solution_terms, p_ref, refmodel, intercept) {
    project_submodel(
      solution_terms = solution_terms, p_ref = p_ref, refmodel = refmodel,
      family = family, intercept = intercept, regul = regul
    )
  })
}

.get_submodels <- function(search_path, nterms, family, p_ref,
                           refmodel, intercept, regul, cv_search = FALSE) {
  ##
  ##
  ## Project onto given model sizes nterms. Returns a list of submodels. If
  ## cv_search=FALSE, submodels parameters will be as they were computed during
  ## the search, so there is no need to project anything anymore, and this
  ## function simply fetches the information from the search_path list, which
  ## contains the parameter values.

  varorder <- search_path$solution_terms
  p_sel <- search_path$p_sel

  if (!cv_search) {
    ## simply fetch the already computed quantities for each submodel size
    fetch_submodel <- function(nterms) {
      solution_terms <- utils::head(varorder, nterms)

      ref_mu <- p_sel$mu

      validparams <- .validate_wobs_wsample(
        refmodel$wobs, p_sel$weights, ref_mu
      )
      wobs <- validparams$wobs
      wsample <- validparams$wsample

      ## reuse sub_fit as projected during search
      sub_refit <- search_path$sub_fits[[nterms + 1]]

      return(.init_submodel(
        sub_fit = sub_refit, p_ref = p_sel, refmodel = refmodel,
        family = family, solution_terms = solution_terms, ref_mu = ref_mu,
        weights = wobs, wsample = wsample
      ))
    }
  } else {
    ## need to project again for each submodel size
    projfun <- .get_proj_handle(family, regul)
    fetch_submodel <- function(nterms) {
      if (nterms == 0) {
        ## empty
        solution_terms <- c("1")
      } else {
        solution_terms <- varorder[seq_len(nterms)]
      }
      return(projfun(solution_terms, p_ref, refmodel, intercept))
    }
  }
  submodels <- lapply(nterms, fetch_submodel)
  return(submodels)
}

.validate_wobs_wsample <- function(ref_wobs, ref_wsample, ref_mu) {
  if (is.null(ref_wobs)) {
    wobs <- rep(1.0, NROW(ref_mu))
  } else {
    wobs <- ref_wobs
  }

  if (is.null(ref_wsample)) {
    wsample <- rep(1.0, NCOL(ref_mu))
  } else {
    wsample <- ref_wsample
  }

  wobs <- wobs / sum(wobs)
  wsample <- wsample / sum(wsample)
  return(nlist(wobs, wsample))
}

.init_submodel <- function(sub_fit, p_ref, refmodel, family, solution_terms,
                           ref_mu, weights, wsample) {
  pobs <- pseudo_data(
    f = 0, y = ref_mu, family = family, weights = weights,
    offset = refmodel$offset
  )

  ## split b to alpha and beta, add it to submodel and return the result
  if (family$family == "gaussian") {
    ref <- list(mu = pobs$z, var = p_ref$var, w = pobs$w)
  } else {
    ref <- p_ref
    ref$w <- rep(0, NROW(ref_mu))
  }

  mu <- family$mu_fun(sub_fit, offset = refmodel$offset, weights = weights)
  dis <- family$dis_fun(ref, list(mu = mu), ref$w)
  kl <- family$kl(
    ref, list(weights = weights),
    nlist(mu, dis)
  )
  weights <- wsample
  solution_terms <- solution_terms
  submodel <- nlist(dis, kl, weights = wsample, solution_terms, sub_fit)
  return(submodel)
}
