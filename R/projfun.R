# Function to project the reference model onto a single submodel with predictor
# terms given in `predictor_terms`. Note that "single submodel" does not refer
# to a single fit (there are as many fits for this single submodel as there are
# projected draws). The case `is.null(search_control)` occurs in two situations:
# (i) when called from search_forward() with `...` as the intended control
# arguments and (ii) when called from perf_eval(). At the end, init_submodl() is
# called, so the output is of class `submodl`.
proj_to_submodl <- function(predictor_terms, p_ref, refmodel,
                            search_control = NULL, verbose, ...) {
  y_unqs_aug <- refmodel$family$cats
  if (refmodel$family$for_latent && !is.null(y_unqs_aug)) {
    y_unqs_aug <- NULL
  }
  terms_ <- unique(unlist(predictor_terms))
  if (!is.null(search_control)) {
    dot_args_divmin <- search_control
  } else {
    dot_args_divmin <- list(...)
  }
  # Fast path for the multi-draw fitters of divmin() (see divmin_multi()): the
  # submodel formula, the draw-wise formulas and the shared model matrix depend
  # only on the reference model and the predictor terms, so they are cached
  # (see proj_cache_store()). The submodel fits are the same as from divmin().
  use_fast <- is.null(y_unqs_aug) &&
    identical(refmodel$div_minimizer, divmin) &&
    proj_fast_ok(dot_args_divmin) &&
    !anyNA(p_ref$mu)
  store <- if (use_fast) proj_cache_store(refmodel) else NULL
  ent <- NULL
  if (!is.null(store)) {
    cache_key <- paste(c(terms_, NCOL(p_ref$mu)), collapse = "\r")
    ent <- store[[cache_key]]
  }
  subset <- NULL
  if (is.null(ent)) {
    subset <- subset_formula_and_data(
      formula = refmodel$formula, terms_ = terms_,
      data = refmodel$fetch_data(), y = p_ref$mu, y_unqs = y_unqs_aug
    )
    fml_divmin <- flatten_formula(subset$formula)
    prep <- if (use_fast) divmin_prepare(fml_divmin, subset$data) else NULL
    if (!is.null(store)) {
      assign(cache_key, list(fml_divmin = fml_divmin, prep = prep),
             envir = store)
    }
  } else {
    fml_divmin <- ent$fml_divmin
    prep <- ent$prep
  }

  if (verbose >= 2L) {
    rhs_chr <- as.character(fml_divmin)
    if (length(rhs_chr) != 3) {
      rhs_chr <- paste("<EXCEPTION: Unexpected length of the character-coerced",
                       "formula passed to the divergence minimizer.>")
    }
    verb_out("  Projecting onto `~ ", utils::tail(rhs_chr, 1), "`")
  }

  outdmin <- NULL
  if (use_fast && !is.null(prep)) {
    outdmin <- divmin_multi(
      prep$formulas, shared_x = prep$shared_x,
      dot_args = c(list(family = refmodel$family, weights = refmodel$wobs),
                   dot_args_divmin),
      throw_warn_sdivmin = getOption("projpred.warn_proj_drawwise", TRUE),
      do_check_conv = getOption("projpred.check_convergence", TRUE),
      y = unclass(p_ref$mu)
    )
  }
  if (is.null(outdmin)) {
    if (is.null(subset)) {
      subset <- subset_formula_and_data(
        formula = refmodel$formula, terms_ = terms_,
        data = refmodel$fetch_data(), y = p_ref$mu, y_unqs = y_unqs_aug
      )
    }
    projpred_ws_aug <- p_ref$mu
    if (refmodel$family$for_augdat && any(refmodel$wobs != 1)) {
      # The augmented-data submodel fitters use the probabilities of the
      # response categories as weights and ignore `weights`, so apply the
      # observation weights to these probabilities. The rows of the
      # augmented-rows matrix are the observations within the response
      # categories:
      projpred_ws_aug <- projpred_ws_aug *
        rep(refmodel$wobs,
            times = nrow(projpred_ws_aug) / length(refmodel$wobs))
    }
    args_divmin <- c(list(formula = fml_divmin,
                          data = subset$data,
                          family = refmodel$family,
                          weights = refmodel$wobs,
                          projpred_var = p_ref$var,
                          projpred_ws_aug = projpred_ws_aug,
                          verbose_divmin = verbose >= 3L),
                     dot_args_divmin)
    outdmin <- do.call(refmodel$div_minimizer, args_divmin)
  }

  return(init_submodl(
    outdmin = outdmin, p_ref = p_ref, refmodel = refmodel,
    predictor_terms = predictor_terms, wobs = refmodel$wobs
  ))
}

# Whether proj_to_submodl() may call divmin_multi() directly instead of divmin()
# (with the arguments `dot_args` for the divergence minimizer): divmin() has to
# use fit_glm_ridge_callback() and to take its warning and convergence settings
# from the current global options, and `dot_args` must not contain arguments of
# divmin() itself.
proj_fast_ok <- function(dot_args) {
  if (!identical(getOption("projpred.glm_fitter", "fit_glm_ridge_callback"),
                 "fit_glm_ridge_callback")) {
    return(FALSE)
  }
  if ((is.null(getOption("projpred.warn_proj_drawwise")) &&
       !is.null(getOption("projpred.warn_prj_drawwise"))) ||
      (is.null(getOption("projpred.check_convergence")) &&
       !is.null(getOption("projpred.check_conv")))) {
    return(FALSE)
  }
  return(!any(names(dot_args) %in% c("throw_warn_sdivmin", "do_check_conv",
                                     "formula", "data", "family", "weights",
                                     "projpred_var", "projpred_ws_aug",
                                     "verbose_divmin")))
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
                      return_dis_sub = FALSE, refmodel_fulldata = refmodel,
                      indices_test, newdata_test = NULL,
                      offset_test = refmodel_fulldata$offset[indices_test],
                      wobs_test = refmodel_fulldata$wobs[indices_test],
                      y_test = refmodel_fulldata$y[indices_test],
                      y_oscale_test = refmodel_fulldata$y_oscale[indices_test],
                      verbose = FALSE, verbose_line_length = 5,
                      verbose_txt_add = "", ...) {
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
        wdraws = reweighting_args$wdraws_ref, cl = reweighting_args$cl_ref,
        clust_used_forced = clust_info(
          ndraws = ndraws,
          nclusters = nclusters,
          S = length(refmodel$wdraws_ref)
        )[["clust_used"]]
      )
    }
    fetch_submodl <- function(size_j, ...) {
      return(proj_to_submodl(
        predictor_terms = utils::head(search_path$predictor_ranking, size_j),
        p_ref = p_ref, refmodel = refmodel, verbose = verbose, ...
      ))
    }
  }
  verb_out(rep("-", verbose_line_length), "\nRunning the performance evaluation ",
           verbose_txt_add, "with ",
           txt_clust_draws(p_ref[["clust_used"]], p_ref[["nprjdraws"]]),
           " (`refit_prj = ", refit_prj, "`) ...",
           verbose = verbose && !return_submodls)
  out_by_size <- NULL
  if (refit_prj && !is.null(reweighting_args) && !return_submodls &&
      !return_preds && getOption("projpred.search_gauss_crossprod", TRUE)) {
    # In the LOO CV with validated search, compute the evaluation of Gaussian
    # submodels from cross-products (`NULL` if this case does not apply):
    out_by_size <- perf_eval_gauss(
      search_path = search_path, nterms = nterms, refmodel = refmodel,
      p_ref = p_ref, refmodel_fulldata = refmodel_fulldata,
      indices_test = indices_test, newdata_test = newdata_test,
      offset_test = offset_test, wobs_test = wobs_test, y_test = y_test,
      y_oscale_test = y_oscale_test, divmin_args = list(...)
    )
  }
  if (is.null(out_by_size)) out_by_size <- lapply(nterms, function(size_j) {
    # Fetch the init_submodl() output (of class `submodl`) for the submodel at
    # position `size_j + 1` of the predictor ranking:
    submodl <- fetch_submodl(size_j, ...)
    if (return_submodls) {
      # Currently only the case if called from project().
      return(submodl)
    }
    if (return_preds) {
      # Currently only the case if called in loo_varsel()'s `validate_search =
      # FALSE` case.
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
        data_aux_test = newdata_test %||%
          refmodel_fulldata$fetch_data(obs = indices_test),
        family = refmodel$family,
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
    out_j <- c(out_j, list(ce = submodl[["ce"]]))
    if (return_dis_sub) {
      # Currently only called in loo_varsel()'s `validate_search = FALSE` case.
      out_j <- c(out_j, list(dis_sub = submodl[["dis"]]))
    }
    return(out_j)
  })
  verb_out(rep("-", verbose_line_length), verbose = verbose && !return_submodls)
  if (return_submodls) {
    # Currently only called in project().
    return(out_by_size)
  }
  if (return_preds) {
    # Currently only called in loo_varsel()'s `validate_search = FALSE` case.
    out <- list(mu_by_size = lapply(out_by_size, "[[", "mu_j"),
                lppd_by_size = lapply(out_by_size, "[[", "lppd_j"))
  } else {
    out <- list(sub_summaries = lapply(out_by_size, "[[", "sub_summary"))
  }
  out <- c(out, list(ce = sapply(out_by_size, "[[", "ce")),
           p_ref[c("clust_used", "nprjdraws")])
  if (return_p_ref) {
    # Currently only called in loo_varsel()'s `validate_search = FALSE` case.
    out <- c(out, nlist(p_ref))
  }
  if (return_dis_sub) {
    # Currently only called in loo_varsel()'s `validate_search = FALSE` case.
    out <- c(out, list(dis_sub = lapply(out_by_size, "[[", "dis_sub")))
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
          const_wdraws_prj = p_ref$const_wdraws_prj, predictor_terms, outdmin,
          cl_ref = p_ref$cl, wdraws_ref = p_ref$wdraws_orig,
          clust_used = p_ref$clust_used, nprjdraws = p_ref$nprjdraws),
    class = "submodl"
  ))
}
