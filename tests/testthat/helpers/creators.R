imperfect_alternation <- function(unq_vals, n_tail = 6L, ...) {
  x <- rep_len(unq_vals, ...)
  tail_idxs_x <- tail(seq_along(x), n = n_tail)
  x[tail_idxs_x] <- x[c(tail_idxs_x[-1],
                        tail_idxs_x[1])]
  return(x)
}

# A Bernoulli or Poisson reference model with simulated draws whose linear
# predictor reaches |eta| of about 100 to 300 for some draws (no model fit is
# necessary; the draws are not a posterior).
refm_extreme <- function(family, n = 40, ndraws = 60, seed = 5523,
                         offset = FALSE, n_extreme = 10) {
  withr::local_seed(seed)
  dat <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  dat$offs <- if (offset) runif(n, -1, 1) else 0
  eta_true <- 1.2 * dat$x1 - 0.8 * dat$x2 + dat$offs
  dat$y <- if (family$family == "binomial") {
    rbinom(n, 1, plogis(eta_true))
  } else {
    rpois(n, exp(eta_true))
  }
  # `n_extreme` of the draws have a coefficient of x1 between 50 and 100:
  beta <- rbind(rnorm(ndraws, 0, 0.1), rnorm(ndraws, 1.2, 0.2),
                rnorm(ndraws, -0.8, 0.2))
  beta[2, seq_len(n_extreme)] <- sample(c(-1, 1), n_extreme, replace = TRUE) *
    runif(n_extreme, 50, 100)
  fml <- y ~ x1 + x2
  ref_predfun <- function(fit, newdata = NULL) {
    nd <- if (is.null(newdata)) fit$data else newdata
    model.matrix(delete.response(terms(fml)), nd) %*% fit$beta
  }
  extract_model_data <- function(object, newdata, wrhs = NULL, orhs = NULL,
                                 extract_y = TRUE) {
    y_wobs_offs(newdata = newdata, wrhs = wrhs,
                orhs = orhs %||% if (offset) ~ offs else NULL,
                resp_form = if (extract_y) ~ y else NULL)
  }
  # For K-fold CV, each fold gets the full-data draws:
  fit <- list(beta = beta, data = dat)
  make_refm <- function(fit_k) {
    init_refmodel(fit_k, data = dat, formula = fml, family = family,
                  ref_predfun = ref_predfun,
                  extract_model_data = extract_model_data,
                  cvfun = function(folds) rep(list(fit), max(folds)),
                  cvrefbuilder = make_refm)
  }
  make_refm(fit)
}
