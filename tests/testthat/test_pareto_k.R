test_that(paste(
  ".ps_khat_threshold() gives the threshold of posterior, capped at 0.7"
), {
  for (S in c(40, 100, 1000, 2200, 4000, 1e5)) {
    expect_equal(.ps_khat_threshold(S),
                 min(posterior::ps_khat_threshold(S), 0.7), info = S)
    # The values of loo's internal function:
    expect_equal(.ps_khat_threshold(S), min(1 - 1 / log10(S), 0.7),
                 info = S)
  }
})

test_that(paste(
  "warn_pareto() counts and prints the k values with the threshold for `S`",
  "draws"
), {
  # setup.R turns off these warnings:
  withr::local_options(projpred.warn_psis = TRUE)
  txt <- "Some (%d / %d) k's are > %s."
  # The threshold for 100 draws is 0.5, for 4000 draws 0.7:
  k <- c(0.55, 0.65, 0.2, Inf)
  expect_warning(warn_pareto(k, S = 100, n = 4, warn_txt = txt),
                 "Some (3 / 4) k's are > 0.5.", fixed = TRUE)
  expect_warning(warn_pareto(k, S = 4000, n = 4, warn_txt = txt),
                 "Some (1 / 4) k's are > 0.7.", fixed = TRUE)
  expect_silent(warn_pareto(c(0.2, 0.3), S = 4000, n = 2, warn_txt = txt))
  withr::local_options(projpred.warn_psis = FALSE)
  expect_silent(warn_pareto(k, S = 100, n = 4, warn_txt = txt))
})

test_that(paste(
  "cv_varsel() with LOO and the latent projection counts an",
  "expectation-specific Pareto k with the Pareto k of loo::psis() if",
  "loo::E_loo() cannot estimate it"
), {
  skip_on_cran()
  # setup.R turns off the Pareto k warnings:
  withr::local_options(projpred.warn_psis = TRUE)
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
  # The k values of loo::E_loo() and of the loo::psis() object that it gets:
  k_eloo <- new.env()
  ns_loo <- asNamespace("loo")
  skip_if_not(exists("E_loo.matrix", ns_loo, inherits = FALSE))
  withr::defer(suppressMessages(untrace("E_loo.matrix", where = ns_loo)))
  suppressMessages(trace(
    "E_loo.matrix", where = ns_loo, print = FALSE,
    exit = bquote({
      assign("k", returnValue()$pareto_k, envir = .(k_eloo))
      assign("k_psis", loo::pareto_k_values(psis_object), envir = .(k_eloo))
    })
  ))
  warns_lat <- character()
  # With 5 projected draws, cv_varsel() gives a message about SIS:
  cvvs <- withCallingHandlers(
    suppressMessages(cv_varsel(
      ref, cv_method = "LOO", validate_search = FALSE, nclusters_pred = 5,
      nterms_max = 1, seed = 1, verbose = FALSE
    )),
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
  k_psis <- k_eloo$k_psis
  expect_length(k, n)
  expect_length(k_psis, n)
  thres <- .ps_khat_threshold(S)
  # For observation 1, loo::psis() gives `Inf`. loo::E_loo() gives `Inf`, or
  # `NA` or `-Inf` in loo 2.10.1:
  expect_identical(k_psis[1], Inf)
  expect_true(is.na(k[1]) || is.infinite(k[1]))
  expect_false(anyNA(k_psis))
  # An observation counts if one of its two k values is above the threshold:
  n_high <- sum(k_psis > thres | (!is.na(k) & k > thres))
  expect_length(warns_lat, 1)
  expect_match(warns_lat, sprintf("some (%d / %d)", n_high, n), fixed = TRUE)
})

test_that(paste(
  "cv_varsel() with LOO counts the Pareto k values of the submodel evaluation",
  "with the threshold for the number of projected draws"
), {
  skip_on_cran()
  # setup.R turns off the Pareto k warnings:
  withr::local_options(projpred.warn_psis = TRUE)
  withr::local_seed(8135)
  n <- 30
  S <- 4000
  S_prj <- 30
  dat <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  dat$y <- 0.5 * dat$x1 + rnorm(n)
  beta <- rbind(rnorm(S, 0, 0.2), rnorm(S, 0.5, 0.2), rnorm(S, 0, 0.2))
  # init_refmodel() warns about the default `cvrefbuilder`:
  ref <- suppressWarnings(init_refmodel(
    object = list(beta = beta), data = dat, formula = y ~ x1 + x2,
    family = gaussian(), dis = rep(1, S),
    ref_predfun = function(fit, newdata = NULL) {
      nd <- if (is.null(newdata)) dat else newdata
      cbind(1, as.matrix(nd[, c("x1", "x2")])) %*% fit$beta
    }
  ))
  # The number of draws and the k values of each loo::psis() call:
  psis_calls <- new.env()
  psis_calls$out <- list()
  ns_loo <- asNamespace("loo")
  skip_if_not(exists("psis.matrix", ns_loo, inherits = FALSE))
  withr::defer(suppressMessages(untrace("psis.matrix", where = ns_loo)))
  suppressMessages(trace(
    "psis.matrix", where = ns_loo, print = FALSE,
    exit = bquote(assign("out", c(
      get("out", envir = .(psis_calls)),
      list(list(S = dim(returnValue())[1],
                k = loo::pareto_k_values(returnValue())))
    ), envir = .(psis_calls)))
  ))
  warns_sub <- character()
  cvvs <- withCallingHandlers(
    cv_varsel(ref, cv_method = "LOO", validate_search = FALSE,
              nclusters_pred = S_prj, nterms_max = 1, seed = 1,
              verbose = FALSE),
    warning = function(w) {
      if (grepl("PSIS-LOO weights given", conditionMessage(w),
                fixed = TRUE)) {
        warns_sub <<- c(warns_sub, conditionMessage(w))
      }
      invokeRestart("muffleWarning")
    }
  )
  expect_s3_class(cvvs, "vsel")
  out_sub <- Filter(function(o) o$S == S_prj, psis_calls$out)
  expect_length(out_sub, 1)
  k_sub <- out_sub[[1]]$k
  thres_prj <- .ps_khat_threshold(S_prj)
  thres_ref <- .ps_khat_threshold(S)
  # The fixture has k values between the two thresholds, so that the count
  # depends on the threshold:
  expect_true(any(k_sub > thres_prj & k_sub <= thres_ref))
  expect_length(warns_sub, 1)
  expect_match(warns_sub, sprintf("Some (%d / %d) ", sum(k_sub > thres_prj), n),
               fixed = TRUE)
  expect_match(warns_sub, sprintf("are > %s.", round(thres_prj, 2)),
               fixed = TRUE)
})
