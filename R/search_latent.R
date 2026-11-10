# Latent-screened forward search for cumulative models with the augmented-data
# projection (experimental; global option `projpred.search_screen_latent`,
# with argument `search_screen` = m). For projected draw s, let p_c be the
# reference probabilities of the categories and z_c = F^-1(P(y <= c)) with the
# link F. The latent predictor is eta_ref = -mean_c z_c (up to a constant), and
# the Fisher information of the cumulative model for eta is
#   V = w * sum_c (f(z_c) - f(z_{c-1}))^2 / p_c
# with f the density of F, f(z_0) = f(z_C) = 0 and the observation weights w.
# The approximate gain of a candidate with the new column x_J at the submodel
# with the columns X_A is the reduction of the weighted residual sum of squares
# when eta_ref is regressed on [1, X_A, x_J] with the weights V, averaged over
# the draws with `wdraws_prj`. The gains only select the candidates to project;
# the exact projections choose the term. Returns `NULL` if a condition does not
# hold; then search_forward() continues with its own search.
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
  wobs <- des$wobs
  n_obs <- nrow(x_c)
  wdraws <- p_ref$wdraws_prj
  n_draws <- length(wdraws)
  n_cat <- attr(p_ref$mu, "ndiscrete")
  if (is.null(n_cat) || length(p_ref$mu) != n_obs * n_cat * n_draws ||
      length(wobs) != n_obs) {
    return(NULL)
  }
  p_arr <- array(unclass(p_ref$mu), dim = c(n_obs, n_cat, n_draws))
  q_fun <- if (fam$link == "probit") stats::qnorm else stats::qlogis
  d_fun <- if (fam$link == "probit") stats::dnorm else stats::dlogis
  eta_ref <- matrix(NA_real_, n_obs, n_draws)
  v_ref <- matrix(NA_real_, n_obs, n_draws)
  for (s in seq_len(n_draws)) {
    p <- pmax(p_arr[, , s], 1e-12)
    p <- p / rowSums(p)
    p_cum <- t(apply(p, 1, cumsum))[, -n_cat, drop = FALSE]
    z <- q_fun(pmin(pmax(p_cum, 1e-12), 1 - 1e-12))
    dens <- cbind(0, d_fun(z), 0)
    dp <- dens[, -1, drop = FALSE] - dens[, -(n_cat + 1L), drop = FALSE]
    eta_ref[, s] <- -rowMeans(z)
    v_ref[, s] <- wobs * rowSums(dp^2 / p)
  }
  if (!all(is.finite(eta_ref)) || !all(is.finite(v_ref))) {
    return(NULL)
  }
  step_cands <- fast_step_cands(search_terms, des$cols_by_term)

  # The gains of the candidates with the new columns `cols_cands` at the
  # submodel with the columns `cols_chosen`; `NULL` if X_A' V X_A is not
  # positive definite for a draw:
  latent_gains <- function(cols_cands, cols_chosen) {
    x_a <- cbind(1, x_c[, cols_chosen, drop = FALSE])
    x_cand <- x_c[, unlist(cols_cands), drop = FALSE]
    x_cand_sq <- x_cand^2
    gain <- numeric(ncol(x_cand))
    for (s in seq_len(n_draws)) {
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
      gain <- gain + wdraws[s] * ifelse(info > 1e-10 * info_tot, g^2 / info, 0)
    }
    return(gain)
  }

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
      gain <- latent_gains(st$cols, cols_chosen)
    }
    fit_idx <- if (is.null(gain)) {
      seq_along(cands)
    } else {
      utils::head(order(gain, decreasing = TRUE), screen_m)
    }
    submodls <- lapply(fit_idx, function(k) {
      proj_to_submodl(c(chosen, cands[k]), p_ref = p_ref, refmodel = refmodel,
                      verbose = verbose, ...)
    })
    ibest <- which.min(vapply(submodls, `[[`, numeric(1), "ce"))
    if (length(ibest) == 0) {
      return(NULL)
    }
    k <- fit_idx[ibest]
    outdmins <- c(outdmins, list(submodls[[ibest]]$outdmin))
    rm(submodls)
    chosen <- c(chosen, cands[k])
    labels_chosen <- c(labels_chosen, st$labels_new[[k]])
    cols_chosen <- c(cols_chosen, st$cols[[k]])
  }
  return(nlist(predictor_ranking = setdiff(chosen, "1"), outdmins))
}

# Latent screens (experimental; global option `projpred.search_screen_latent`):
# the screen of the first function with a result, or `NULL`.
search_forward_latent <- function(...) {
  out <- search_forward_latent_cum(...)
  if (is.null(out)) {
    out <- search_forward_latent_mlvl(...)
  }
  return(out)
}

# Latent-screened forward search for the latent projection with group-level
# terms (experimental). The submodels are Gaussian multilevel models. For the
# screen, the projected draws are combined to one response (the mean of the
# latent predictor over the projected draws with `wdraws_prj`). The criterion
# of a candidate is the weighted residual sum of squares (conditional fitted
# values) of one lme4::lmer() fit (maximum likelihood) of this response on the
# terms of the submodel plus the candidate; without group-level terms, it is a
# weighted least-squares fit. The screen_m candidates with the smallest
# criterion get an exact projection, and the exact projections choose the
# term. Returns `NULL` if a condition does not hold; then search_forward()
# continues with its own search.
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
  dat <- dat[wobs > 0, , drop = FALSE]
  wobs <- wobs[wobs > 0]
  env_fml <- environment(refmodel$formula)

  # The right-hand side of a submodel. A term "(x | g)" includes the
  # group-level intercept, so "(1 | g)" of the same factor is dropped:
  rhs_of <- function(terms) {
    tt <- unique(trimws(unlist(strsplit(setdiff(terms, "1"), " + ",
                                        fixed = TRUE))))
    is_gl <- grepl("|", tt, fixed = TRUE)
    grp <- sub(".*\\|\\s*(.*)\\)$", "\\1", tt[is_gl])
    lhs <- trimws(sub("^\\((.*)\\|.*$", "\\1", tt[is_gl]))
    keep <- !(lhs == "1" & grp %in% grp[lhs != "1"])
    list(rhs = paste(c("1", tt[!is_gl], tt[is_gl][keep]), collapse = " + "),
         has_gl = any(is_gl))
  }
  # The criterion of a submodel (`NA` if the fit fails):
  rss_of <- function(terms) {
    sp <- rhs_of(terms)
    tryCatch({
      if (!sp$has_gl) {
        mm <- model.matrix(as.formula(paste("~", sp$rhs), env = env_fml), dat)
        fit <- lm.wfit(mm, dat$.projpred_eta, wobs)
        sum(wobs * fit$residuals^2)
      } else {
        fit <- suppressMessages(suppressWarnings(lme4::lmer(
          as.formula(paste(".projpred_eta ~", sp$rhs), env = env_fml),
          data = dat, weights = .projpred_w, REML = FALSE,
          control = lme4::lmerControl(calc.derivs = FALSE,
                                      check.conv.singular = "ignore")
        )))
        sum(wobs * (dat$.projpred_eta - stats::fitted(fit))^2)
      }
    }, error = function(e) NA_real_)
  }

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
      rss <- vapply(cands, function(cand) rss_of(c(chosen, cand)), numeric(1))
      # A failed screen fit counts as the best, so that the candidate gets
      # an exact projection:
      rss[is.na(rss)] <- -Inf
      fit_idx <- utils::head(order(rss), screen_m)
    }
    submodls <- lapply(fit_idx, function(k) {
      proj_to_submodl(c(chosen, cands[k]), p_ref = p_ref, refmodel = refmodel,
                      verbose = verbose, ...)
    })
    ibest <- which.min(vapply(submodls, `[[`, numeric(1), "ce"))
    if (length(ibest) == 0) {
      return(NULL)
    }
    outdmins <- c(outdmins, list(submodls[[ibest]]$outdmin))
    rm(submodls)
    chosen <- c(chosen, cands[fit_idx[ibest]])
  }
  return(nlist(predictor_ranking = setdiff(chosen, "1"), outdmins))
}
