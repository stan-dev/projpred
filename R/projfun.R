# Function to project the reference model onto a single submodel with predictor
# terms given in `predictor_terms`. Note that "single submodel" does not refer
# to a single fit (there are as many fits for this single submodel as there are
# projected draws). The case `is.null(search_control)` occurs in two situations:
# (i) when called from search_forward() with `...` as the intended control
# arguments and (ii) when called from perf_eval(). At the end, init_submodl() is
# called, so the output is of class `submodl`.
proj_to_submodl <- function(predictor_terms, p_ref, refmodel,
                            search_control = NULL, ...) {
  y_unqs_aug <- refmodel$family$cats
  if (refmodel$family$for_latent && !is.null(y_unqs_aug)) {
    y_unqs_aug <- NULL
  }
  subset <- subset_formula_and_data(
    formula = refmodel$formula, terms_ = unique(unlist(predictor_terms)),
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

  args_divmin <- list(formula = fml_divmin,
                      data = subset$data,
                      family = refmodel$family,
                      weights = refmodel$wobs,
                      projpred_var = p_ref$var,
                      projpred_ws_aug = p_ref$mu)
  if (!is.null(search_control)) {
    args_divmin <- c(args_divmin, search_control)
  } else {
    args_divmin <- c(args_divmin, list(...))
  }
  outdmin <- do.call(refmodel$div_minimizer, args_divmin)

  return(init_submodl(
    outdmin = outdmin, p_ref = p_ref, refmodel = refmodel,
    predictor_terms = predictor_terms, wobs = refmodel$wobs
  ))
}

# Function to prepare the performance evaluation. For each submodel size along
# the predictor ranking, it first fetches the init_submodl() output (of class
# `submodl`) and then calculates the submodel "summary" (precursor quantities)
# for the actual performance evaluation performed by summary.vsel() later.
perf_eval <- function(search_path,
                      nterms = c(0, seq_along(search_path$predictor_ranking)),
                      refmodel, refit_prj = FALSE, ndraws, nclusters,
                      reweighting_args = NULL, return_submodls = FALSE,
                      return_preds = FALSE, return_p_ref = FALSE,
                      refmodel_fulldata = refmodel,
                      indices_test, newdata_test = NULL,
                      offset_test = refmodel_fulldata$offset[indices_test],
                      wobs_test = refmodel_fulldata$wobs[indices_test],
                      y_test = refmodel_fulldata$y[indices_test],
                      y_oscale_test = refmodel_fulldata$y_oscale[indices_test],
                      ...) {
  if (!refit_prj) {
    p_ref <- search_path$p_sel
    # In this case, simply fetch the already computed projections, so don't
    # project again.
    fetch_submodl <- function(size_j, ...) {
      return(init_submodl(
        # Re-use the submodel fits from the search:
        outdmin = search_path$outdmins[[size_j + 1]],
        p_ref = p_ref,
        refmodel = refmodel,
        predictor_terms = utils::head(search_path$predictor_ranking, size_j),
        wobs = refmodel$wobs
      ))
    }
  } else {
    # In this case, project again.
    if (is.null(reweighting_args)) {
      p_ref <- get_refdist(refmodel, ndraws = ndraws, nclusters = nclusters)
    } else {
      # Reweight the clusters (or thinned draws) according to the PSIS weights:
      p_ref <- get_p_clust(
        family = refmodel$family, eta = refmodel$eta, mu = refmodel$mu,
        mu_offs = refmodel$mu_offs, dis = refmodel$dis,
        wdraws = reweighting_args$wdraws_ref, cl = reweighting_args$cl_ref
      )
    }
    fetch_submodl <- function(size_j, ...) {
      return(proj_to_submodl(
        predictor_terms = utils::head(search_path$predictor_ranking, size_j),
        p_ref = p_ref, refmodel = refmodel, ...
      ))
    }
  }
  out_by_size <- lapply(nterms, function(size_j) {
    # Fetch the init_submodl() output (of class `submodl`) for the submodel at
    # position `size_j + 1` of the predictor ranking:
    submodl <- fetch_submodl(size_j, ...)
    if (return_submodls) {
      # Currently only called in project().
      return(submodl)
    }
    if (return_preds) {
      # Currently only called in loo_varsel()'s `validate_search = FALSE` case.
      mu_j <- refmodel$mu_fun(submodl$outdmin, obs = indices_test,
                              offset = refmodel$offset[indices_test])
      lppd_j <- t(refmodel$family$ll_fun(
        mu_j, submodl$dis, refmodel$y[indices_test], refmodel$wobs[indices_test]
      ))
      out_j <- nlist(mu_j, lppd_j)
    } else {
      # Calculate precursor quantities for predictive performance statistic(s)
      # of the submodel at position `size_j + 1` of the predictor ranking:
      sub_summary <- weighted_summary_means(
        y_wobs_test = data.frame(y = y_test, y_oscale = y_oscale_test,
                                 wobs = wobs_test),
        family = refmodel_fulldata$family,
        wdraws = submodl$wdraws_prj,
        mu = refmodel_fulldata$mu_fun(submodl$outdmin,
                                      obs = indices_test,
                                      newdata = newdata_test,
                                      offset = offset_test),
        dis = submodl$dis,
        cl_ref = submodl$cl_ref,
        wdraws_ref = submodl$wdraws_ref
      )
      out_j <- nlist(sub_summary)
    }
    return(c(out_j, list(ce = submodl[["ce"]])))
  })
  if (return_submodls) {
    return(out_by_size)
  }
  if (return_preds) {
    out <- list(mu_by_size = lapply(out_by_size, "[[", "mu_j"),
                lppd_by_size = lapply(out_by_size, "[[", "lppd_j"))
  } else {
    out <- list(sub_summaries = lapply(out_by_size, "[[", "sub_summary"))
  }
  out <- c(out, list(ce = sapply(out_by_size, "[[", "ce"),
                     clust_used = p_ref$clust_used,
                     nprjdraws = p_ref$nprjdraws))
  if (return_p_ref) {
    # Currently only called in loo_varsel()'s `validate_search = FALSE` case.
    out <- c(out, nlist(p_ref))
  }
  return(out)
}

# Process the output of the `div_minimizer` function (see init_refmodel()) to
# create an object of class `submodl`.
init_submodl <- function(outdmin, p_ref, refmodel, predictor_terms, wobs) {
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

  mu <- refmodel$mu_fun(outdmin, offset = refmodel$offset)
  dis <- refmodel$family$dis_fun(p_ref, nlist(mu), wobs)
  ce <- weighted.mean(
    refmodel$family$ce(p_ref,
                       nlist(weights = wobs),
                       nlist(mu, dis)),
    p_ref$wdraws_prj
  )
  return(structure(
    nlist(dis, ce, wdraws_prj = p_ref$wdraws_prj,
          const_wdraws_prj = length(unique(p_ref$wdraws_prj)) == 1,
          predictor_terms, outdmin, cl_ref = p_ref$cl,
          wdraws_ref = p_ref$wdraws_orig, clust_used = p_ref$clust_used,
          nprjdraws = p_ref$nprjdraws),
    class = "submodl"
  ))
}
