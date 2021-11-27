# Function to project the reference model onto a single submodel with predictor
# terms given in `solution_terms`. Note that "single submodel" does not refer to
# a single fit (there are as many fits for this single submodel as there are
# projected draws).
project_submodel <- function(solution_terms, p_ref, refmodel, regul = 1e-4) {
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
    family = refmodel$family,
    weights = refmodel$wobs,
    projpred_var = p_ref$var,
    projpred_regul = regul
  )

  if (isTRUE(getOption("projpred.check_conv", FALSE))) {
    check_conv(sub_fit)
  }

  return(.init_submodel(
    sub_fit = sub_fit, p_ref = p_ref, refmodel = refmodel,
    family = refmodel$family, solution_terms = solution_terms, wobs = wobs,
    wsample = wsample
  ))
}

# Function to project the reference model onto the submodels of given model
# sizes `nterms`. Returns a list of submodels (each processed by
# .init_submodel()).
.get_submodels <- function(search_path, nterms, family, p_ref, refmodel, regul,
                           cv_search = FALSE) {
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
        refmodel = refmodel, regul = regul
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
  p_ref$mu <- family$linkinv(family$linkfun(p_ref$mu) + refmodel$offset)
  if (!(all(is.na(p_ref$var)) ||
        family$family %in% c("gaussian", "Student_t"))) {
    stop("For family `", family$family, "()`, .init_submodel() might have to ",
         "be adapted, depending on whether family$predvar() is invariant with ",
         "respect to offsets (this would be OK and does not need an ",
         "adaptation) or not (this would need an adaptation).")
  }
  if (family$family == "Student_t") {
    stop("For the `Student_t()` family, .init_submodel() is not finished yet.")
    ### TODO (`Student_t()` family): Check if this is needed (perhaps with some
    ### modifications) or if something completely different is needed (there
    ### used to be no special handling of the `Student_t()` family here at all):
    # pobs <- pseudo_data(
    #   f = 0, y = p_ref$mu, family = family, weights = wobs,
    #   offset = refmodel$offset
    # )
    # ### TODO: Add `dis` and perhaps other elements here?:
    # p_ref <- list(mu = pobs$z, var = p_ref$var)
    # ###
    # p_ref$mu <- family$linkinv(family$linkfun(p_ref$mu) + refmodel$offset)
    # wobs <- pobs$wobs
    ###
  }

  mu <- family$mu_fun(sub_fit, offset = refmodel$offset)
  dis <- family$dis_fun(p_ref, nlist(mu), wobs)
  kl <- weighted.mean(
    family$kl(p_ref,
              nlist(weights = wobs),
              nlist(mu, dis)),
    wsample
  )
  return(nlist(dis, kl, weights = wsample, solution_terms, sub_fit))
}
