search_forward <- function(p_ref, refmodel, nterms_max, verbose = TRUE,
                           search_terms, est_runtime = TRUE,
                           search_terms_was_null, ...) {
  nterms_max_with_icpt <- nterms_max + 1L
  iq <- ceiling(quantile(seq_len(nterms_max_with_icpt), 1:10 / 10))
  if (is.null(search_terms)) {
    stop("Did not expect `search_terms` to be `NULL`. Please report this.")
  }

  chosen <- character()
  outdmins <- c()

  for (size in seq_len(nterms_max_with_icpt)) {
    # Determine candidate predictors for the current size:
    cands <- select_possible_terms_size(chosen, search_terms, size = size)
    if (is.null(cands))
      next
    full_cands <- lapply(cands, function(cand) c(chosen, cand))

    # Perform the projections (for the intercept-only model, we measure the
    # runtime to estimate the total runtime for the search):
    if (size == 1 && est_runtime && getOption("projpred.mssg_time", TRUE)) {
      time_bef <- Sys.time()
    }
    submodls <- lapply(full_cands, proj_to_submodl, p_ref = p_ref,
                       refmodel = refmodel, ...)
    if (size == 1 && est_runtime && getOption("projpred.mssg_time", TRUE)) {
      time_aft <- Sys.time()
      dtime <- difftime(time_aft, time_bef, units = "secs")
      # Scale up to the remaining forward search where we have p * (p + 1) / 2
      # projections (with `p = nterms_max`) if `search_terms` is at its default:
      time_est_min <- dtime * nterms_max * (nterms_max + 1) / 2
      # Scale from the intercept-only submodel to GLM submodels (see PR #459 for
      # an empirical derivation of the factor):
      time_est_min <- time_est_min * 1.3
      # Initialize upper bound:
      time_est_max <- time_est_min
      # Adjust for multilevel and/or additive terms if necessary (see PR #459
      # for an empirical derivation of the factors):
      has_mlvl <- any(grepl("\\|", search_terms))
      has_smooth <- length(parse_additive_terms(search_terms)) > 0
      if (has_mlvl && !has_smooth) {
        time_est_max <- time_est_max * 26.9
      } else if (!has_mlvl && has_smooth) {
        time_est_max <- time_est_max * 9.8
      } else if (has_mlvl && has_smooth) {
        time_est_max <- time_est_max * 57.6
      }
      if (time_est_max > 3 * 60) {
        mssg_time_start <- "The remaining forward search is estimated to take "
        if (time_est_max == time_est_min) {
          mssg_time_est <- paste0(
            "ca. ", round(time_est_min / 60, 1), " minutes "
          )
        } else {
          mssg_time_est <- paste0(
            "between ca. ", round(time_est_min / 60, 1), " and ",
            round(time_est_max / 60, 1), " minutes "
          )
        }
        mssg_time_start_curr <- paste0(
          "(current time: ", format(time_aft, usetz = TRUE),
          "; estimated end time: "
        )
        if (time_est_max == time_est_min) {
          mssg_time_est_curr <- paste0(
            format(time_aft + time_est_min, usetz = TRUE), "). Note that ",
            "this is only a rough estimate."
          )
        } else {
          mssg_time_est_curr <- paste0(
            "between ", format(time_aft + time_est_min, usetz = TRUE), " and ",
            format(time_aft + time_est_max, usetz = TRUE), "). Note that ",
            "these are not guaranteed bounds but the bounds of a rough ",
            "interval estimate."
          )
        }
        mssg_time <- paste0(mssg_time_start, mssg_time_est,
                            mssg_time_start_curr, mssg_time_est_curr)
        if (has_mlvl) {
          # The empirically derived factor used above assumes that *all*
          # projections involve one multilevel term (and apart from that, the
          # empirical derivation used only group-level intercepts, not
          # group-level slopes):
          mssg_time <- paste0(
            mssg_time, " Furthermore, since there are multilevel predictor ",
            "terms, this estimate may be even more unreliable."
          )
        }
        if (has_smooth) {
          # The empirically derived factor used above assumes that *all*
          # projections involve one multilevel and one smooth term (and apart
          # from that, the empirical derivation used only group-level intercepts
          # and an s() smooth term, not more complex terms):
          mssg_time <- paste0(
            mssg_time, " Furthermore, since there are additive (\"smooth\") ",
            "predictor terms, this estimate may be even more unreliable."
          )
        }
        if (!search_terms_was_null) {
          # In case of custom `search_terms`, some model sizes may be skipped:
          mssg_time <- paste0(
            mssg_time, " Since argument `search_terms` differs from its ",
            "default, this estimate may be an overestimate."
          )
        }
        message(mssg_time)
      }
    }

    # Select best candidate:
    imin <- which.min(sapply(submodls, "[[", "ce"))
    chosen <- c(chosen, cands[imin])

    # Store `outdmin` (i.e., the object containing the projection results)
    # corresponding to the best candidate:
    outdmins <- c(outdmins, list(submodls[[imin]]$outdmin))

    # Verbose mode output:
    ct_chosen <- count_terms_chosen(chosen)
    if (verbose && ct_chosen %in% iq) {
      vtxt <- paste(names(iq)[max(which(ct_chosen == iq))], "of terms selected")
      if (getOption("projpred.extra_verbose", FALSE)) {
        vtxt <- paste0(vtxt, ": ", paste(chosen, collapse = " + "))
      }
      verb_out(vtxt)
    }

    # Free up some memory:
    rm(submodls)
    if (getOption("projpred.run_gc", FALSE)) {
      gc()
    }
  }

  # For `predictor_ranking`, `reduce_models(chosen)` used to be used instead of
  # `chosen`. However, `reduce_models(chosen)` and `chosen` should be identical
  # at this place because select_possible_terms_size() already avoids redundant
  # models. Thus, use `chosen` here because it matches `outdmins` (this
  # matching is necessary because later in perf_eval()'s `!refit_prj` case,
  # `outdmins` is indexed with integers which are based on `predictor_ranking`):
  return(nlist(predictor_ranking = setdiff(chosen, "1"), outdmins))
}

#' Force search terms
#'
#' A helper function to construct the input for argument `search_terms` of
#' [varsel()] or [cv_varsel()] if certain predictor terms should be forced to be
#' selected first whereas other predictor terms are optional (i.e., they are
#' subject to the variable selection, but only after the inclusion of the
#' "forced" terms).
#'
#' @param forced_terms A character vector of predictor terms that should be
#'   selected first.
#' @param optional_terms A character vector of predictor terms that should be
#'   subject to the variable selection after the inclusion of the "forced"
#'   terms.
#'
#' @return A character vector that may be used as input for argument
#'   `search_terms` of [varsel()] or [cv_varsel()].
#'
#' @seealso [varsel()], [cv_varsel()]
#'
#' @examplesIf requireNamespace("rstanarm", quietly = TRUE)
#' # Data:
#' dat_gauss <- data.frame(y = df_gaussian$y, df_gaussian$x)
#'
#' # The `stanreg` fit which will be used as the reference model (with small
#' # values for `chains` and `iter`, but only for technical reasons in this
#' # example; this is not recommended in general):
#' fit <- rstanarm::stan_glm(
#'   y ~ X1 + X2 + X3 + X4 + X5, family = gaussian(), data = dat_gauss,
#'   QR = TRUE, chains = 2, iter = 500, refresh = 0, seed = 9876
#' )
#'
#' # We will force X1 and X2 to be selected first:
#' search_terms_forced <- force_search_terms(
#'   forced_terms = paste0("X", 1:2),
#'   optional_terms = paste0("X", 3:5)
#' )
#'
#' # Run varsel() (here without cross-validation and with small values for
#' # `nterms_max`, `nclusters`, and `nclusters_pred`, but only for the sake of
#' # speed in this example; this is not recommended in general):
#' vs <- varsel(fit, nclusters = 5, nclusters_pred = 10,
#'              search_terms = search_terms_forced, seed = 5555)
#' # Now see, for example, `?print.vsel`, `?plot.vsel`, `?suggest_size.vsel`,
#' # and `?ranking` for possible post-processing functions.
#'
#' @export
force_search_terms <- function(forced_terms, optional_terms) {
  stopifnot(length(forced_terms) > 0)
  stopifnot(length(optional_terms) > 0)
  forced_terms <- paste(forced_terms, collapse = " + ")
  return(c(forced_terms, paste0(forced_terms, " + ", optional_terms)))
}

search_L1_surrogate <- function(p_ref, d_train, family, intercept, nterms_max,
                                penalty, opt) {

  ## predictive mean and variance of the reference model (with parameters
  ## integrated out)
  mu <- p_ref$mu
  v <- p_ref$var

  if (NCOL(mu) > 1 || NCOL(v) > 1) {
    stop("Internal error: search_L1 received multiple draws. ",
         "Please report to the developers.")
  }

  ## L1-penalized projection (projection path).
  ## (Notice: here we use pmax = nterms_max+1 so that the computation gets
  ## carried until all the way down to the least regularization also for model
  ## size nterms_max)
  search <- glm_elnet(d_train$x, mu, family,
                      lambda_min_ratio = opt$lambda_min_ratio,
                      nlambda = opt$nlambda,
                      pmax = nterms_max + 1, pmax_strict = FALSE,
                      weights = d_train$weights,
                      intercept = intercept, obsvar = v, penalty = penalty,
                      thresh = opt$thresh)

  ## sort the variables according to the order in which they enter the model in
  ## the L1-path
  entering_indices <- apply(search$beta != 0, 1, function(num) {
    which(num)[1] # na for those that did not enter
  })
  ## variables that entered at some point
  entered_variables <- c(seq_len(NCOL(d_train$x)))[!is.na(entering_indices)]
  ## variables that did not enter at any point
  notentered_variables <- c(seq_len(NCOL(d_train$x)))[is.na(entering_indices)]
  order_of_entered <- sort(entering_indices, index.return = TRUE)$ix
  order <- c(entered_variables[order_of_entered], notentered_variables)

  ## fetch the coefficients corresponding to those points at the search_path
  ## where new variable enters
  nvar <- length(order)
  n <- nrow(p_ref$mu)
  out <- list(
    alpha = rep(NA, nterms_max + 1),
    beta = matrix(0, nrow = nterms_max, ncol = nterms_max + 1),
    lambda = rep(NA, nterms_max + 1),
    w = matrix(NA, nrow = n, ncol = nterms_max + 1)
  )
  for (k in 0:nterms_max) {
    if (k == 0) {
      out$alpha[1] <- search$beta0[1]
      out$lambda[1] <- search$lambda[1]
      out$w[, 1] <- search$w[, 1]
    } else {
      ## find those points in the L1-path where only the k most relevant
      ## features can have nonzero coefficient, and then fetch their
      ## coefficients with least regularization
      ivar <- utils::tail(order, nvar - k)
      steps_k_var <- which(colSums(search$beta[ivar, , drop = FALSE] != 0) == 0)
      if (length(steps_k_var) > 0) {
        j <- utils::tail(steps_k_var, 1)
      } else {
        ## no steps where all the variables in set ivar would have zero
        ## coefficient (could be due to one or more of these variables having
        ## penalty = 0 so they are always in the model) so set the coefficients
        ## to be equal to the starting value
        j <- 1
      }
      out$alpha[k + 1] <- search$beta0[j]
      out$beta[1:k, k + 1] <- search$beta[order[1:k], j]
      out$lambda[k + 1] <- search$lambda[j]
      out$w[, k + 1] <- search$w[, j]
    }
  }

  out$predictor_ranking <- order[seq_len(nterms_max)]
  if (any(is.na(out$predictor_ranking)) &&
      length(entered_variables) < nterms_max) {
    if (length(setdiff(notentered_variables,
                       which(penalty == Inf))) > 0) {
      warning("Less than nterms_max variables entered L1-path. ",
              "Try reducing lambda_min_ratio. ")
    }
  }

  return(out)
}

search_L1 <- function(p_ref, refmodel, nterms_max, penalty, opt) {
  if (nterms_max == 0) {
    stop("L1 search cannot be used for an empty (i.e. intercept-only) ",
         "full-model formula or `nterms_max = 0`.")
  }
  # Preparations:
  fr <- model.frame(refmodel$formula, data = refmodel$fetch_data(),
                    drop.unused.levels = TRUE)
  da_classes <- attr(attr(fr, "terms"), "dataClasses")
  nms_chr_fac <- names(da_classes)[da_classes %in% c("character", "factor")]
  resp_nm <- all.vars(attr(fr, "terms"))[attr(attr(fr, "terms"), "response")]
  nms_chr_fac <- setdiff(nms_chr_fac, resp_nm)
  if (length(nms_chr_fac) > 0) {
    xlvls <- lapply(setNames(nm = nms_chr_fac), function(nm_chr_fac) {
      levels(as.factor(fr[[nm_chr_fac]]))
    })
  } else {
    xlvls <- NULL
  }
  # TODO: In the following model.matrix() call, allow user-specified contrasts
  # to be passed to argument `contrasts.arg`. The `contrasts.arg` default
  # (`NULL`) uses `options("contrasts")` internally, but it might be more
  # convenient to let users specify contrasts directly. At that occasion,
  # contrasts should also be tested thoroughly (not done until now).
  x <- model.matrix(refmodel$formula, data = fr)
  x <- x[, colnames(x) != "(Intercept)", drop = FALSE]
  ## it's important to keep the original order because that's the order
  ## in which lasso will estimate the parameters
  tt <- terms(refmodel$formula)
  terms_ <- attr(tt, "term.labels")
  search_path <- search_L1_surrogate(
    p_ref, nlist(x, weights = refmodel$wobs), refmodel$family,
    intercept = TRUE, ncol(x), penalty, opt
  )

  predictor_ranking_orig <- collapse_ranked_predictors(
    path = colnames(x)[search_path$predictor_ranking],
    formula = refmodel$formula, data = fr
  )
  predictor_ranking <- utils::head(predictor_ranking_orig, nterms_max)
  # Place lower-order interaction terms before higher-order interaction terms,
  # but otherwise preserve the ranking:
  ia_orders <- sapply(gregexpr(":", predictor_ranking), function(greg_colon) {
    sum(greg_colon != -1)
  })
  ia_order_max <- max(ia_orders)
  for (ia_order in rev(seq_len(ia_order_max))) {
    ias <- predictor_ranking[ia_orders == ia_order]
    stopifnot(!any(duplicated(ias))) # safety measure for which.max()
    for (ia in ias) {
      ia_idx <- which.max(predictor_ranking == ia)
      if (ia_idx > nterms_max) break
      main_terms_ia <- strsplit(ia, ":")[[1]]
      ias_lower_split <- utils::combn(main_terms_ia, m = ia_order,
                                      simplify = FALSE)
      ias_lower <- lapply(ias_lower_split, all_ia_perms, is_split = TRUE)
      ias_lower <- unlist(ias_lower)
      ias_lower <- intersect(ias_lower, predictor_ranking_orig)
      prev_terms <- utils::head(predictor_ranking, ia_idx - 1L)
      has_lower_after <- !all(ias_lower %in% prev_terms)
      if (has_lower_after) {
        if (getOption("projpred.warn_L1_interactions", TRUE)) {
          warning("Interaction term `", ia, "` was selected before all ",
                  "corresponding lower-order interaction terms have been ",
                  "selected. This is a known deficiency of L1 search. Use ",
                  "forward search to avoid this. Now ranking the lower-order ",
                  "interaction terms before this interaction term.")
        }
        ias_lower <- setdiff(ias_lower, prev_terms)
        ias_lower <- ias_lower[order(match(ias_lower, predictor_ranking_orig))]
        new_head <- c(prev_terms, ias_lower, ia)
        predictor_ranking <- c(new_head, setdiff(predictor_ranking, new_head))
        predictor_ranking <- utils::head(predictor_ranking, nterms_max)
      }
    }
  }

  outdmins <- lapply(0:length(predictor_ranking), function(nterms) {
    if (nterms == 0) {
      formula <- make_formula(c("1"))
      beta <- NULL
      x <- x[, numeric(), drop = FALSE]
    } else {
      formula <- make_formula(predictor_ranking[seq_len(nterms)])
      variables <- unlist(lapply(
        predictor_ranking[seq_len(nterms)],
        function(term) {
          # TODO: In the following model.matrix() call, allow user-specified
          # contrasts to be passed to argument `contrasts.arg`. The
          # `contrasts.arg` default (`NULL`) uses `options("contrasts")`
          # internally, but it might be more convenient to let users specify
          # contrasts directly. At that occasion, contrasts should also be
          # tested thoroughly (not done until now).
          mm <- model.matrix(as.formula(paste("~ 1 +", term)), data = fr)
          return(setdiff(colnames(mm), "(Intercept)"))
        }
      ))
      indices <- match(variables, colnames(x)[search_path$predictor_ranking])
      indices <- indices[!is.na(indices)]
      beta <- search_path$beta[indices, max(indices) + 1, drop = FALSE]
      # Also reduce `x` (important for coef.subfit(), for example); note that
      # `x <- x[, variables, drop = FALSE]` should also be possible, but the
      # re-use of `colnames(x)` should provide another sanity check:
      x <- x[
        , colnames(x)[search_path$predictor_ranking[indices]], drop = FALSE
      ]
      # For consistency with fit_glm_ridge_callback():
      rownames(beta) <- colnames(x)
    }
    # Avoid model.frame.default()'s warning "variable '<...>' is not a factor"
    # when calling predict.subfit() later:
    xlvls <- xlvls[intersect(names(xlvls), all.vars(formula))]
    if (length(xlvls) == 0) {
      xlvls <- NULL
    }

    sub <- nlist(alpha = search_path$alpha[nterms + 1], beta,
                 w = search_path$w[, nterms + 1], formula, x, xlvls)
    class(sub) <- "subfit"
    return(list(sub))
  })

  return(nlist(predictor_ranking, outdmins))
}
