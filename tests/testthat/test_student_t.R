test_that("Student_t() is deprecated and its projection stops", {
  # In edition 3, expect_warning() does not hide other warnings:
  local_edition(3)
  withr::local_options(projpred.warn_cvrefbuilder_NULL = FALSE)
  expect_warning(fam <- Student_t(nu = 4), "Student_t\\(\\) is deprecated")
  expect_identical(fam$family, "Student_t")

  # A reference model with simulated draws (no model fit is necessary):
  withr::local_seed(4731)
  n <- 20
  ndraws <- 10
  dat <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  dat$y <- dat$x1 + rt(n, df = 4)
  fit <- list(beta = rbind(rnorm(ndraws, 0, 0.1), rnorm(ndraws, 1, 0.1),
                           rnorm(ndraws, 0, 0.1)))
  ref_predfun <- function(fit, newdata = NULL) {
    nd <- if (is.null(newdata)) dat else newdata
    cbind(1, nd$x1, nd$x2) %*% fit$beta
  }
  extract_model_data <- function(object, newdata = NULL, wrhs = NULL,
                                 orhs = NULL, extract_y = TRUE) {
    nd <- if (is.null(newdata)) dat else newdata
    list(y = if (extract_y) nd$y, weights = rep(1, nrow(nd)),
         offset = rep(0, nrow(nd)))
  }
  expect_warning(
    refm <- init_refmodel(fit, data = dat, formula = y ~ x1 + x2,
                          family = fam, ref_predfun = ref_predfun,
                          extract_model_data = extract_model_data,
                          dis = rep(1, ndraws)),
    "`Student_t\\(\\)` family is deprecated"
  )
  expect_error(project(refm, predictor_terms = "x1", verbose = 0),
               "use the latent projection")
  expect_error(varsel(refm, method = "forward", nterms_max = 1, verbose = 0),
               "use the latent projection")
  expect_error(varsel(refm, method = "L1", nterms_max = 1, verbose = 0),
               "use the latent projection")
})

test_that("the latent projection of brms::student() uses family$linkinv", {
  skip_if_not_installed("brms")
  mssgs <- character()
  withCallingHandlers(
    extend_family(brms::student(), latent = TRUE),
    message = function(m) {
      mssgs <<- c(mssgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  expect_false(any(grepl("Defining `latent_ilink`", mssgs)))
})
