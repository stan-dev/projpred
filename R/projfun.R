# Function to project the reference model onto a single submodel with predictor
# terms given in `solution_terms`. Note that "single submodel" does not refer to
# a single fit (there are as many fits for this single submodel as there are
# projected draws). At the end, init_submodl() is called, so the output is of
# class `submodl`.
get_submodl_prj <- function(solution_terms, p_ref, refmodel, regul = 1e-4,
                            ...) {
  validparams <- validate_wobs_wdraws(refmodel$wobs, p_ref$wdraws_prj,
                                      p_ref$mu)
  wobs <- validparams$wobs
  wdraws_prj <- validparams$wdraws_prj

  y_unqs_aug <- refmodel$family$cats
  if (refmodel$family$for_latent && !is.null(y_unqs_aug)) {
    y_unqs_aug <- NULL
  }
  subset <- subset_formula_and_data(
    formula = refmodel$formula, terms_ = unique(unlist(solution_terms)),
    data = refmodel$fetch_data(), y = p_ref$mu, y_unqs = y_unqs_aug
  )
  fml_divmin <- flatten_formula(subset$formula)

  if (getOption("projpred.extra_verbose", FALSE)) {
    rhs_chr <- as.character(fml_divmin)
    if (length(rhs_chr) != 3) {
      rhs_chr <- paste("<EXCEPTION: Unexpected length of the character-coerced",
                       "formula passed to the divergence minimizer.>")
    }
    verb_out("  Projecting onto ", utils::tail(rhs_chr, 1))
  }

  outdmin <- refmodel$div_minimizer(
    formula = fml_divmin,
    data = subset$data,
    family = refmodel$family,
    weights = refmodel$wobs,
    projpred_var = p_ref$var,
    projpred_regul = regul,
    projpred_ws_aug = p_ref$mu,
    ...
  )

  if (isTRUE(getOption("projpred.check_conv", FALSE))) {
    check_conv(outdmin)
  }

  return(init_submodl(
    outdmin = outdmin, p_ref = p_ref, refmodel = refmodel,
    solution_terms = solution_terms, wobs = wobs, wdraws_prj = wdraws_prj
  ))
}

# Function to fetch init_submodl() output (of class `submodl`) for each of given
# model sizes `nterms`, so this gives a list of objects, each of class
# `submodl`.
get_submodls <- function(search_path, nterms, p_ref, refmodel, regul,
                         refit_prj = FALSE, ...) {
  if (!refit_prj) {
    # In this case, simply fetch the already computed projections, so don't
    # project again.
    fetch_submodl <- function(nterms, ...) {
      validparams <- validate_wobs_wdraws(refmodel$wobs,
                                          search_path$p_sel$wdraws_prj,
                                          search_path$p_sel$mu)
      wobs <- validparams$wobs
      wdraws_prj <- validparams$wdraws_prj
      return(init_submodl(
        # Re-use the submodel fits from the search:
        outdmin = search_path$outdmins[[nterms + 1]],
        p_ref = search_path$p_sel,
        refmodel = refmodel,
        solution_terms = utils::head(search_path$solution_terms, nterms),
        wobs = wobs,
        wdraws_prj = wdraws_prj
      ))
    }
  } else {
    # In this case, project again.
    fetch_submodl <- function(nterms, ...) {
      return(get_submodl_prj(
        solution_terms = utils::head(search_path$solution_terms, nterms),
        p_ref = p_ref, refmodel = refmodel, regul = regul, ...
      ))
    }
  }
  return(lapply(nterms, fetch_submodl, ...))
}

validate_wobs_wdraws <- function(ref_wobs, ref_wdraws, ref_mu) {
  if (is.null(ref_wobs)) {
    wobs <- rep(1.0, NROW(ref_mu))
  } else {
    wobs <- ref_wobs
  }

  if (is.null(ref_wdraws)) {
    wdraws_prj <- rep(1.0, NCOL(ref_mu))
  } else {
    wdraws_prj <- ref_wdraws
  }

  wdraws_prj <- wdraws_prj / sum(wdraws_prj)
  return(nlist(wobs, wdraws_prj))
}

# Process the output of the `divergence_minimizer` function (see
# init_refmodel()) to create an object of class `submodl`.
init_submodl <- function(outdmin, p_ref, refmodel, solution_terms, wobs,
                         wdraws_prj) {
  p_ref$mu <- p_ref$mu_offs
  if (!(all(is.na(p_ref$var)) ||
        refmodel$family$family %in% c("gaussian", "Student_t"))) {
    stop("For family `", refmodel$family$family, "()`, init_submodl() might ",
         "have to be adapted, depending on whether family$predvar() is ",
         "invariant with respect to offsets (this would be OK and does not ",
         "need an adaptation) or not (this would need an adaptation).")
  }
  if (refmodel$family$family == "Student_t") {
    stop("For the `Student_t()` family, init_submodl() is not finished yet.")
    ### TODO (Student_t()): Check if this is needed (perhaps with some
    ### modifications) or if something completely different is needed (there
    ### used to be no special handling of the `Student_t()` family here at all):
    # pobs <- pseudo_data(
    #   f = 0, y = p_ref$mu, family = refmodel$family, weights = wobs,
    #   offset = refmodel$offset
    # )
    # ### TODO (Student_t()): Add `dis` and perhaps other elements here?:
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

  mu <- refmodel$family$mu_fun(outdmin, offset = refmodel$offset)
  dis <- refmodel$family$dis_fun(p_ref, nlist(mu), wobs)
  ce <- weighted.mean(
    refmodel$family$ce(p_ref,
                       nlist(weights = wobs),
                       nlist(mu, dis)),
    wdraws_prj
  )
  return(structure(
    nlist(dis, ce, wdraws_prj = wdraws_prj, solution_terms, outdmin,
          cl_ref = p_ref$cl, wdraws_ref = p_ref$wdraws_orig),
    class = "submodl"
  ))
}
