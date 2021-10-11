# Function to project the reference model onto a single submodel with predictor
# terms given in `solution_terms`. Note that "single submodel" does not refer to
# a single fit (there are as many fits for this single submodel as there are
# projected draws).
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
    family = family,
    weights = refmodel$wobs,
    projpred_var = p_ref$var,
    projpred_regul = regul
  )

  return(.init_submodel(
    sub_fit = sub_fit, p_ref = p_ref, refmodel = refmodel,
    family = family, solution_terms = solution_terms, wobs = wobs,
    wsample = wsample
  ))
}

# Function to project the reference model onto the submodels of given model
# sizes `nterms`. Returns a list of submodels (each processed by
# .init_submodel()).
.get_submodels <- function(search_path, nterms, family, p_ref,
                           refmodel, intercept, regul, cv_search = FALSE) {
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
    fetch_submodel <- function(nterms) {
      project_submodel(
        solution_terms = utils::head(varorder, nterms), p_ref = p_ref,
        refmodel = refmodel, family = family, intercept = intercept,
        regul = regul
      )
    }
  }
  return(lapply(nterms, fetch_submodel))
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

  mu <- family$mu_fun(sub_fit, offset = refmodel$offset)
  dis <- family$dis_fun(ref, nlist(mu), ref$wobs)
  kl <- weighted.mean(family$kl(
    ref, nlist(weights = wobs),
    nlist(mu, dis)
  ), wsample)
  weights <- wsample
  return(nlist(dis, kl, weights, solution_terms, sub_fit))
}
