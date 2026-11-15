# Latent screens (argument `search_screen` = m) for models without a score
# screen. As in
# search_forward_screen(), a cheap criterion ranks the candidates of a step,
# only the m best get an exact projection, and the exact projections choose
# the term. Returns the result of the first screen whose conditions hold, or
# `NULL`; then search_forward() continues with its own search.
search_forward_latent <- function(...) {
  out <- search_forward_latent_cum(...)
  if (is.null(out)) {
    out <- search_forward_latent_mlvl(...)
  }
  return(out)
}

# Verbose mode output of the latent screens, as in search_forward(). `iq` is
# the named vector of the sizes at which a line is printed.
latent_verb_progress <- function(chosen, iq, verbose) {
  ct_chosen <- count_terms_chosen(chosen)
  if (verbose && ct_chosen %in% iq) {
    vtxt <- paste(names(iq)[max(which(ct_chosen == iq))], "of terms selected")
    if (verbose >= 2L) {
      vtxt <- paste0(vtxt, ": `~ ", paste(chosen, collapse = " + "), "`")
    }
    verb_out(vtxt)
  }
  return(invisible(NULL))
}

# One step of the latent screens: projects the submodels `chosen` plus each
# candidate `cands[fit_idx]` and returns the index `k` (into `cands`) and the
# `outdmin` of the submodel with the smallest `ce`. Returns `NULL` if no `ce`
# is available (all `NA`).
latent_project_best <- function(chosen, cands, fit_idx, p_ref, refmodel,
                                verbose, ...) {
  submodls <- lapply(fit_idx, function(k) {
    proj_to_submodl(c(chosen, cands[k]), p_ref = p_ref, refmodel = refmodel,
                    verbose = verbose, ...)
  })
  ibest <- which.min(vapply(submodls, `[[`, numeric(1), "ce"))
  if (length(ibest) == 0) {
    return(NULL)
  }
  return(list(k = fit_idx[ibest], outdmin = submodls[[ibest]]$outdmin))
}

# Latent-screened forward search for cumulative models (family "cumulative",
# logit or probit link) with the augmented-data projection and the default
# divergence minimizer. For projected draw s, let p_c be the reference
# probabilities of the categories and z_c = F^-1(P(y <= c)) with the link F.
# The latent predictor is eta_ref = -mean_c z_c (the mean over the C - 1
# thresholds; equal to eta up to a constant), and the Fisher information of
# the cumulative model for eta is
#   V = w * sum_c (f(z_c) - f(z_{c-1}))^2 / p_c
# with f the density of F, f(z_0) = f(z_C) = 0 and the observation weights w.
# The approximate gain of a candidate with the new column x_J at the submodel
# with the columns X_A is the reduction of the weighted residual sum of squares
# when eta_ref is regressed on [1, X_A, x_J] with the weights V, averaged over
# the draws with `wdraws_prj`.
# Limits (the function returns `NULL`): offsets; group-level or smooth terms;
# a search term with more than one model-matrix column (for example, a factor
# with more than two levels).
search_forward_latent_cum <- function(p_ref, refmodel, nterms_max, verbose,
                                      search_terms, screen_m, ...) {
  fam <- refmodel$family
  if (!isTRUE(fam$for_augdat) || !identical(fam$family, "cumulative") ||
      !isTRUE(fam$link %in% c("logit", "probit")) ||
      !identical(refmodel$div_minimizer, divmin_augdat) ||
      any(refmodel$offset != 0) || anyNA(p_ref$mu) ||
      !isTRUE(screen_m >= 1L)) {
    return(NULL)
  }
  des <- gauss_crossprod_design(refmodel, search_terms, regul = 1)
  if (is.null(des) || any(lengths(des$cols_by_term) != 1L)) {
    return(NULL)
  }
  x_c <- des$x_c
  n_obs <- nrow(x_c)
  wdraws <- p_ref$wdraws_prj
  n_draws <- length(wdraws)
  n_cat <- attr(p_ref$mu, "ndiscrete")
  if (is.null(n_cat) || length(p_ref$mu) != n_obs * n_cat * n_draws ||
      length(des$wobs) != n_obs) {
    return(NULL)
  }
  lat <- latent_cum_response(
    array(unclass(p_ref$mu), dim = c(n_obs, n_cat, n_draws)),
    wobs = des$wobs, link = fam$link
  )
  if (!all(is.finite(lat$eta)) || !all(is.finite(lat$v))) {
    return(NULL)
  }
  step_cands <- fast_step_cands(search_terms, des$cols_by_term)
  verb_out("Screening the candidates of each step with the latent ",
           "approximation of the cumulative model.",
           verbose = verbose)

  iq <- ceiling(quantile(seq_len(nterms_max + 1L), 1:10 / 10))
  chosen <- character()
  labels_chosen <- character()
  cols_chosen <- integer(0)
  outdmins <- list()
  for (size in seq_len(nterms_max + 1L)) {
    st <- step_cands(size, chosen, labels_chosen)
    if (is.null(st)) {
      return(NULL)
    }
    cands <- st$cands
    if (is.null(cands)) {
      next
    }
    gain <- NULL
    if (size > 1L && length(cands) > screen_m) {
      gain <- latent_cum_gains(x_c = x_c, eta_ref = lat$eta, v_ref = lat$v,
                               wdraws = wdraws, cols_cands = st$cols,
                               cols_chosen = cols_chosen)
    }
    fit_idx <- if (is.null(gain)) {
      seq_along(cands)
    } else {
      utils::head(order(gain, decreasing = TRUE), screen_m)
    }
    best <- latent_project_best(chosen, cands, fit_idx, p_ref = p_ref,
                                refmodel = refmodel, verbose = verbose, ...)
    if (is.null(best)) {
      return(NULL)
    }
    k <- best$k
    outdmins <- c(outdmins, list(best$outdmin))
    chosen <- c(chosen, cands[k])
    labels_chosen <- c(labels_chosen, st$labels_new[[k]])
    cols_chosen <- c(cols_chosen, st$cols[[k]])
    latent_verb_progress(chosen, iq, verbose)
  }
  return(nlist(predictor_ranking = setdiff(chosen, "1"), outdmins))
}

# The latent predictor `eta` and the Fisher information `v` (N x S matrices) of
# search_forward_latent_cum() from the N x C x S array `p_arr` of the reference
# probabilities, the observation weights `wobs` and the link ("logit" or
# "probit").
latent_cum_response <- function(p_arr, wobs, link) {
  n_obs <- dim(p_arr)[1]
  n_cat <- dim(p_arr)[2]
  n_draws <- dim(p_arr)[3]
  q_fun <- if (link == "probit") stats::qnorm else stats::qlogis
  d_fun <- if (link == "probit") stats::dnorm else stats::dlogis
  eta <- matrix(NA_real_, n_obs, n_draws)
  v <- matrix(NA_real_, n_obs, n_draws)
  for (s in seq_len(n_draws)) {
    p <- pmax(matrix(p_arr[, , s], n_obs, n_cat), 1e-12)
    p <- p / rowSums(p)
    p_cum <- matrix(t(apply(p, 1, cumsum)), n_obs, n_cat)
    p_cum <- p_cum[, -n_cat, drop = FALSE]
    z <- q_fun(pmin(pmax(p_cum, 1e-12), 1 - 1e-12))
    dens <- cbind(0, d_fun(z), 0)
    dp <- dens[, -1, drop = FALSE] - dens[, -(n_cat + 1L), drop = FALSE]
    eta[, s] <- -rowMeans(z)
    v[, s] <- wobs * rowSums(dp^2 / p)
  }
  return(nlist(eta, v))
}

# The gains of search_forward_latent_cum() for the candidates with the new
# columns `cols_cands` (a list with one index into the columns of `x_c` per
# candidate) at the submodel with the columns `cols_chosen` and an intercept.
# Returns `NULL` if X_A' V X_A is not positive definite for a draw (then all
# candidates of the step are projected).
latent_cum_gains <- function(x_c, eta_ref, v_ref, wdraws, cols_cands,
                             cols_chosen) {
  x_a <- cbind(1, x_c[, cols_chosen, drop = FALSE])
  x_cand <- x_c[, unlist(cols_cands), drop = FALSE]
  x_cand_sq <- x_cand^2
  gain <- numeric(ncol(x_cand))
  for (s in seq_along(wdraws)) {
    v <- v_ref[, s]
    xav <- x_a * v
    r_a <- tryCatch(chol(crossprod(x_a, xav)), error = function(e) NULL)
    if (is.null(r_a)) {
      return(NULL)
    }
    b_a <- backsolve(r_a, backsolve(r_a, crossprod(xav, eta_ref[, s]),
                                    transpose = TRUE))
    resid <- eta_ref[, s] - drop(x_a %*% b_a)
    g <- drop(crossprod(x_cand, v * resid))
    u <- backsolve(r_a, crossprod(xav, x_cand), transpose = TRUE)
    info_tot <- drop(crossprod(v, x_cand_sq))
    info <- info_tot - colSums(u^2)
    # A column (almost) in the span of the chosen columns gets gain 0:
    gain <- gain + wdraws[s] * ifelse(info > 1e-10 * info_tot, g^2 / info, 0)
  }
  return(gain)
}

# Latent-screened forward search for the latent projection with group-level
# terms (and no smooth terms) and the default divergence minimizer. The
# submodels are Gaussian multilevel models. For the screen, the projected
# draws are combined to one response (the mean of the latent predictor over
# the projected draws with `wdraws_prj`). The criterion of a candidate is the
# weighted residual sum of squares (see latent_mlvl_rss()) of one fit of this
# response on the terms of the submodel plus the candidate. A candidate whose
# screen fit fails gets an exact projection in addition to the m best.
# Limits (the function returns `NULL`): offsets; smooth terms; a reference
# model without group-level terms.
search_forward_latent_mlvl <- function(p_ref, refmodel, nterms_max, verbose,
                                       search_terms, screen_m, ...) {
  fam <- refmodel$family
  if (!isTRUE(fam$for_latent) || isTRUE(fam$for_augdat) ||
      !identical(refmodel$div_minimizer, divmin) ||
      any(refmodel$offset != 0) || anyNA(p_ref$mu) ||
      !isTRUE(screen_m >= 1L) ||
      !formula_contains_group_terms(refmodel$formula) ||
      formula_contains_additive_terms(refmodel$formula)) {
    return(NULL)
  }
  dat <- refmodel$fetch_data()
  wobs <- refmodel$wobs
  mu_ref <- unclass(p_ref$mu)
  if (!is.matrix(mu_ref) || nrow(mu_ref) != nrow(dat) ||
      length(wobs) != nrow(dat) || any(wobs < 0)) {
    return(NULL)
  }
  wdraws <- p_ref$wdraws_prj
  dat$.projpred_eta <- drop(mu_ref %*% (wdraws / sum(wdraws)))
  dat$.projpred_w <- wobs
  # The screen fits use the observations with a positive weight (a fold of
  # LOO CV sets the weight of the left-out observation to zero):
  dat <- droplevels(dat[wobs > 0, , drop = FALSE])
  # The contrasts of the reference model, as in proj_to_submodl() and
  # drop_lost_levels(): a factor with a dropped level keeps the default
  # contrasts (see set_contrasts()). The projection of a submodel without
  # group-level terms keeps the zero-weight rows and their levels, so there a
  # reduced contrast matrix can give another column space than the screen fit:
  dat <- set_contrasts(dat, refmodel$contrasts)
  fml_base <- stats::as.formula(".projpred_eta ~ 1",
                                env = environment(refmodel$formula))
  verb_out("Screening the candidates of each step with one multilevel fit ",
           "to the mean latent predictor.", verbose = verbose)

  iq <- ceiling(quantile(seq_len(nterms_max + 1L), 1:10 / 10))
  chosen <- character()
  outdmins <- list()
  for (size in seq_len(nterms_max + 1L)) {
    cands <- select_possible_terms_size(chosen, search_terms, size = size,
                                        use_term_info = TRUE)
    if (is.null(cands)) {
      next
    }
    fit_idx <- seq_along(cands)
    if (size > 1L && length(cands) > screen_m) {
      rss <- vapply(cands, function(cand) {
        latent_mlvl_rss(c(chosen, cand), dat = dat, fml_base = fml_base)
      }, numeric(1))
      fit_idx <- sort(union(which(is.na(rss)),
                            utils::head(order(rss, na.last = NA), screen_m)))
      if (anyNA(rss)) {
        verb_out("Screen fit failed for ", sum(is.na(rss)), " candidate(s); ",
                 "projecting them, too.", verbose = verbose)
      }
    }
    best <- latent_project_best(chosen, cands, fit_idx, p_ref = p_ref,
                                refmodel = refmodel, verbose = verbose, ...)
    if (is.null(best)) {
      return(NULL)
    }
    outdmins <- c(outdmins, list(best$outdmin))
    chosen <- c(chosen, cands[best$k])
    latent_verb_progress(chosen, iq, verbose)
  }
  return(nlist(predictor_ranking = setdiff(chosen, "1"), outdmins))
}

# The formula of a screen fit of search_forward_latent_mlvl(): the submodel
# with the terms `terms_` and the response of `fml_base`, with the group-level
# terms of each grouping factor combined as for the projection.
latent_mlvl_formula <- function(terms_, fml_base) {
  return(flatten_formula(make_formula(terms_, formula = fml_base)))
}

# The criterion of search_forward_latent_mlvl() for the submodel with the terms
# `terms_`: the weighted residual sum of squares of the response
# `.projpred_eta` with the weights `.projpred_w` (columns of `dat`). With
# group-level terms, the fit is one lme4::lmer() fit (maximum likelihood) and
# the residuals use the conditional fitted values; `NA` if lme4::lmer() fails.
# Without group-level terms, the fit is weighted least squares. A variable
# that is not in `dat` is an error.
latent_mlvl_rss <- function(terms_, dat, fml_base) {
  fml <- latent_mlvl_formula(terms_, fml_base)
  vars_miss <- setdiff(all.vars(fml), names(dat))
  if (length(vars_miss) > 0) {
    stop("Unexpected variable(s) in a screen fit: ",
         paste(vars_miss, collapse = ", "), ". Please notify the package ",
         "maintainer.")
  }
  if (!formula_contains_group_terms(fml)) {
    mm <- stats::model.matrix(stats::delete.response(stats::terms(fml)), dat)
    fit <- stats::lm.wfit(mm, dat$.projpred_eta, dat$.projpred_w)
    return(sum(dat$.projpred_w * fit$residuals^2))
  }
  ctrl <- lme4::lmerControl(calc.derivs = FALSE,
                            check.conv.singular = "ignore")
  # `data` and `weights` are passed as names (lme4::lmer() evaluates `weights`
  # in `data`):
  fit <- tryCatch(
    suppressMessages(suppressWarnings(do.call(lme4::lmer, list(
      formula = fml, data = quote(dat), weights = quote(.projpred_w),
      REML = FALSE, control = ctrl
    )))),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(NA_real_)
  }
  return(sum(dat$.projpred_w * (dat$.projpred_eta - stats::fitted(fit))^2))
}
