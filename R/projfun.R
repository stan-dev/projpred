# Function to project the reference model onto a single submodel with predictor
# terms given in `solution_terms`. Note that "single submodel" does not refer to
# a single fit (there are as many fits for this single submodel as there are
# projected draws).
project_submodel <- function(solution_terms, p_ref, refmodel, regul = 1e-4,
                             ...) {
  validparams <- .validate_wobs_wsample(refmodel$wobs, p_ref$weights, p_ref$mu)
  wobs <- validparams$wobs
  wsample <- validparams$wsample

  y_unqs_aug <- refmodel$family$cats
  if (refmodel$family$for_latent && !is.null(y_unqs_aug)) {
    y_unqs_aug <- NULL
  }
  subset <- subset_formula_and_data(
    formula = refmodel$formula, terms_ = unique(unlist(solution_terms)),
    data = refmodel$fetch_data(), y = p_ref$mu, y_unqs = y_unqs_aug
  )

  submodl <- refmodel$div_minimizer(
    formula = flatten_formula(subset$formula),
    data = subset$data,
    family = refmodel$family,
    weights = refmodel$wobs,
    projpred_var = p_ref$var,
    projpred_regul = regul,
    projpred_ws_aug = p_ref$mu,
    ...
  )

  if (isTRUE(getOption("projpred.check_conv", FALSE))) {
    check_conv(submodl)
  }

  return(.init_submodel(
    submodl = submodl, p_ref = p_ref, refmodel = refmodel,
    solution_terms = solution_terms, wobs = wobs, wsample = wsample
  ))
}

# Function to project the reference model onto the submodels of given model
# sizes `nterms`. Returns a list of submodels (each processed by
# .init_submodel(), so of class `initsubmodl`).
.get_submodels <- function(search_path, nterms, p_ref, refmodel, regul,
                           refit_prj = FALSE, ...) {
  if (!refit_prj) {
    # In this case, simply fetch the already computed projections, so don't
    # project again.
    fetch_submodel <- function(nterms, ...) {
      validparams <- .validate_wobs_wsample(
        refmodel$wobs, search_path$p_sel$weights, search_path$p_sel$mu
      )
      wobs <- validparams$wobs
      wsample <- validparams$wsample
      return(.init_submodel(
        # Re-use the submodel fits from the search:
        submodl = search_path$submodls[[nterms + 1]],
        p_ref = search_path$p_sel,
        refmodel = refmodel,
        solution_terms = utils::head(search_path$solution_terms, nterms),
        wobs = wobs,
        wsample = wsample
      ))
    }
  } else {
    # In this case, project again.
    fetch_submodel <- function(nterms, ...) {
      return(project_submodel(
        solution_terms = utils::head(search_path$solution_terms, nterms),
        p_ref = p_ref, refmodel = refmodel, regul = regul, ...
      ))
    }
  }
  return(lapply(nterms, fetch_submodel, ...))
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

.init_submodel <- function(submodl, p_ref, refmodel, solution_terms, wobs,
                           wsample) {
  # Take offsets into account (the `if ()` condition is added for efficiency):
  if (!all(refmodel$offset == 0)) {
    p_ref$mu <- refmodel$family$linkinv(
      refmodel$family$linkfun(p_ref$mu) + refmodel$offset
    )
  }
  if (!(all(is.na(p_ref$var)) ||
        refmodel$family$family %in% c("gaussian", "Student_t"))) {
    stop("For family `", refmodel$family$family, "()`, .init_submodel() might ",
         "have to be adapted, depending on whether family$predvar() is ",
         "invariant with respect to offsets (this would be OK and does not ",
         "need an adaptation) or not (this would need an adaptation).")
  }
  if (refmodel$family$family == "Student_t") {
    stop("For the `Student_t()` family, .init_submodel() is not finished yet.")
    ### TODO (`Student_t()` family): Check if this is needed (perhaps with some
    ### modifications) or if something completely different is needed (there
    ### used to be no special handling of the `Student_t()` family here at all):
    # pobs <- pseudo_data(
    #   f = 0, y = p_ref$mu, family = refmodel$family, weights = wobs,
    #   offset = refmodel$offset
    # )
    # ### TODO: Add `dis` and perhaps other elements here?:
    # p_ref <- list(mu = pobs$z, var = p_ref$var)
    # ###
    # if (!all(refmodel$offset == 0)) {
    #   p_ref$mu <- refmodel$family$linkinv(
    #     refmodel$family$linkfun(p_ref$mu) + refmodel$offset
    #   )
    # }
    # wobs <- pobs$wobs
    ###
  }

  mu <- refmodel$family$mu_fun(submodl, offset = refmodel$offset)
  dis <- refmodel$family$dis_fun(p_ref, nlist(mu), wobs)
  kl <- weighted.mean(
    refmodel$family$kl(p_ref,
                       nlist(weights = wobs),
                       nlist(mu, dis)),
    wsample
  )
  return(structure(
    nlist(dis, kl, weights = wsample, solution_terms, submodl, cl = p_ref$cl,
          wsample_orig = p_ref$wsample_orig),
    class = "initsubmodl"
  ))
}
