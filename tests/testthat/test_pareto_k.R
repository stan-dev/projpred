test_that(paste(
  "n_high_pareto_k() counts a k value that could not be estimated as `Inf`,",
  "as loo::psis() does"
), {
  k <- c(NA, NaN, -Inf, Inf, 0.2, 0.69, 0.71)
  expect_identical(n_high_pareto_k(k, 0.7), 5L)
  expect_identical(n_high_pareto_k(c(0.1, 0.5), 0.7), 0L)
  expect_identical(n_high_pareto_k(numeric(0), 0.7), 0L)
})

test_that(paste(
  ".ps_khat_threshold() gives the threshold of posterior, capped at 0.7"
), {
  for (S in c(40, 100, 1000, 2200, 4000, 1e5)) {
    expect_equal(.ps_khat_threshold(S), min(1 - 1 / log10(S), 0.7),
                 info = S)
  }
})

test_that(paste(
  "warn_pareto() counts and prints the k values with the threshold for `S`",
  "draws, and counts k values that could not be estimated"
), {
  txt <- "Some (%d / %d) k's are > %s."
  # The threshold for 100 draws is 0.5, for 4000 draws 0.7:
  k <- c(0.55, 0.65, 0.2, NA)
  expect_warning(warn_pareto(k, S = 100, n = 4, warn_txt = txt),
                 "Some (3 / 4) k's are > 0.5.", fixed = TRUE)
  expect_warning(warn_pareto(k, S = 4000, n = 4, warn_txt = txt),
                 "Some (1 / 4) k's are > 0.7.", fixed = TRUE)
  expect_warning(warn_pareto(c(0.2, -Inf), S = 4000, n = 2, warn_txt = txt),
                 "Some (1 / 2) k's are > 0.7.", fixed = TRUE)
  expect_silent(warn_pareto(c(0.2, 0.3), S = 4000, n = 2, warn_txt = txt))
  withr::local_options(projpred.warn_psis = FALSE)
  expect_silent(warn_pareto(k, S = 100, n = 4, warn_txt = txt))
})

test_that(paste(
  "cv_varsel() with LOO and the latent projection counts the",
  "expectation-specific Pareto k values that could not be estimated"
), {
  skip_on_cran()
  withr::local_seed(5629)
  n <- 30
  S <- 40
  x <- matrix(rnorm(n * 2), n, 2, dimnames = list(NULL, c("x1", "x2")))
  dat <- data.frame(x)
  dat$y <- rpois(n, exp(0.5 + 0.5 * x[, 1]))
  beta <- rbind(rnorm(S, 0.5, 0.1), rnorm(S, 0.5, 0.1), rnorm(S, 0, 0.1))
  # Draws 1 to 4 predict a mean far from y_1. If observation 1 is left out,
  # the PSIS weights of all other draws underflow to 0, so that the
  # expectation-specific Pareto k of observation 1 cannot be estimated:
  ref_predfun <- function(fit, newdata = NULL) {
    nd <- if (is.null(newdata)) dat else newdata
    eta <- cbind(1, as.matrix(nd[, colnames(x)])) %*% fit$beta
    if (is.null(newdata)) {
      eta[1, 1:4] <- 20
    }
    eta
  }
  ref <- suppressWarnings(suppressMessages(init_refmodel(
    object = list(beta = beta), data = dat, formula = y ~ x1 + x2,
    family = poisson(), ref_predfun = ref_predfun, dis = rep(1, S),
    latent = TRUE
  )))
  # The k values of loo::E_loo():
  k_eloo <- new.env()
  ns_loo <- asNamespace("loo")
  withr::defer(suppressMessages(untrace("E_loo.matrix", where = ns_loo)))
  suppressMessages(trace(
    "E_loo.matrix", where = ns_loo, print = FALSE,
    exit = bquote(assign("k", returnValue()$pareto_k, envir = .(k_eloo)))
  ))
  warns_lat <- character()
  cvvs <- withCallingHandlers(
    cv_varsel(ref, cv_method = "LOO", validate_search = FALSE,
              nclusters_pred = 5, nterms_max = 1, seed = 1, verbose = FALSE),
    warning = function(w) {
      if (grepl("recalculation of the latent response",
                conditionMessage(w), fixed = TRUE)) {
        warns_lat <<- c(warns_lat, conditionMessage(w))
      }
      # The other warnings (Pareto k of the reference model) are expected:
      invokeRestart("muffleWarning")
    }
  )
  expect_s3_class(cvvs, "vsel")
  k <- k_eloo$k
  expect_length(k, n)
  thres <- .ps_khat_threshold(S)
  # The k value of observation 1 could not be estimated (or is above the
  # threshold):
  expect_true(is.na(k[1]) || k[1] > thres)
  k_inf <- replace(k, is.na(k) | k == -Inf, Inf)
  expect_length(warns_lat, 1)
  expect_match(warns_lat, sprintf("some (%d / %d)", sum(k_inf > thres), n),
               fixed = TRUE)
})
