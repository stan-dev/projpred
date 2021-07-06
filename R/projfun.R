## Function handles for the projection
##

project_submodel <- function(solution_terms, p_ref, refmodel, family, intercept,
                             regul = 1e-4) {
  validparams <- .validate_wobs_wsample(refmodel$wobs, p_ref$weights, p_ref$mu)
  wobs <- validparams$wobs
  wsample <- validparams$wsample

  subset <- subset_formula_and_data(
    formula = refmodel$formula, terms_ = unique(unlist(solution_terms)),
    data = refmodel$fetch_data(), y = p_ref$mu
  )

  sub_fit <- refmodel$div_minimizer(
    formula = flatten_formula(subset$formula),
    data = subset$data,
    weights = refmodel$wobs,
    family = family,
    regul = regul,
    var = p_ref$var
  )

  return(.init_submodel(
    sub_fit = sub_fit, p_ref = p_ref, refmodel = refmodel,
    family = family, solution_terms = solution_terms, wobs = wobs,
    wsample = wsample
  ))
}

## function handle for the projection over samples
.get_proj_handle <- function(refmodel, p_ref, family, regul = 1e-9,
                             intercept = TRUE) {
  return(function(solution_terms) {
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

  if (!cv_search) {
    ## simply fetch the already computed quantities for each submodel size
    fetch_submodel <- function(nterms) {
      solution_terms <- utils::head(varorder, nterms)

      validparams <- .validate_wobs_wsample(
        refmodel$wobs, search_path$p_sel$weights, search_path$p_sel$mu
      )
      wobs <- validparams$wobs
      wsample <- validparams$wsample

      ## reuse sub_fit as projected during search
      sub_refit <- search_path$sub_fits[[nterms + 1]]

      return(.init_submodel(
        sub_fit = sub_refit, p_ref = search_path$p_sel, refmodel = refmodel,
        family = family, solution_terms = solution_terms,
        wobs = wobs, wsample = wsample
      ))
    }
  } else {
    ## need to project again for each submodel size
    projfun <- .get_proj_handle(refmodel, p_ref, family, regul, intercept)
    fetch_submodel <- function(nterms) {
      solution_terms <- varorder[seq_len(nterms)]
      return(projfun(solution_terms))
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
                           wobs, wsample) {
  ## split b to alpha and beta, add it to submodel and return the result
  if (family$family == "gaussian") {
    pobs <- pseudo_data(
      f = 0, y = p_ref$mu, family = family, weights = wobs,
      offset = refmodel$offset
    )
    ref <- list(mu = pobs$z, var = p_ref$var, wobs = pobs$wobs)
  } else {
    ref <- p_ref
  }

  mu <- family$mu_fun(sub_fit, offset = refmodel$offset, weights = 1)
  dis <- family$dis_fun(ref, nlist(mu), ref$wobs)
  kl <- weighted.mean(family$kl(
    ref, nlist(weights = wobs),
    nlist(mu, dis)
  ), wsample)
  weights <- wsample
  submodel <- nlist(dis, kl, weights, solution_terms, sub_fit)
  return(submodel)
}
