context("project")


# tests for project

if (require(rstanarm)) {
  seed <- 1235
  set.seed(seed)
  n <- 40
  nterms <- 5
  ndraws <- 1
  ndraws_pred <- 5
  x <- matrix(rnorm(n * nterms, 0, 1), n, nterms)
  b <- runif(nterms) - 0.5
  dis <- runif(1, 1, 2)
  weights <- sample(1:4, n, replace = TRUE)
  offset <- rnorm(n)
  chains <- 2
  iter <- 500
  source(testthat::test_path("helpers", "SW.R"))

  f_gauss <- gaussian()
  df_gauss <- data.frame(y = rnorm(n, f_gauss$linkinv(x %*% b), dis), x = x)
  f_binom <- binomial()
  df_binom <- data.frame(y = rbinom(n, weights, f_binom$linkinv(x %*% b)),
                         x = x, weights = weights)
  f_poiss <- poisson()
  df_poiss <- data.frame(y = rpois(n, f_poiss$linkinv(x %*% b)), x = x)

  SW({
    fit_gauss <- stan_glm(y ~ x.1 + x.2 + x.3 + x.4 + x.5,
                          family = f_gauss, data = df_gauss, QR = TRUE,
                          weights = weights, offset = offset,
                          chains = chains, seed = seed, iter = iter)
    fit_binom <- stan_glm(cbind(y, weights - y) ~ x.1 + x.2 + x.3 + x.4 + x.5,
                          family = f_binom, weights = weights,
                          data = df_binom, chains = chains, seed = seed,
                          iter = iter)
    fit_poiss <- stan_glm(y ~ x.1 + x.2 + x.3 + x.4 + x.5,
                          family = f_poiss, data = df_poiss,
                          chains = chains, seed = seed, iter = iter)
  })

  fam_nms <- setNames(nm = c("gauss", "binom", "poiss"))
  fit_list <- lapply(fam_nms, function(fam_nm) {
    get(paste0("fit_", fam_nm))
  })

  solterms_tst <- c("x.2", "x.4")
  ndraws_default <- 400L # Adopt this if the default is changed.
  nclusters_tst <- 2L
  nclusters_pred_tst <- 3L
  seed_tst <- 866028

  projection_nms <- c(
    "dis", "kl", "weights", "solution_terms", "sub_fit", "family",
    "p_type", "intercept", "extract_model_data", "refmodel"
  )
  sub_fit_nms <- c("alpha", "beta", "w", "formula", "x", "y")

  # For the binomial family with > 1 trials, we currently expect the warning
  # "Using formula(x) is deprecated when x is a character vector of length > 1"
  # (see GitHub issue #136), so temporarily wrap the following call in SW():
  SW(vs_list <- lapply(fit_list, varsel,
                       nterms_max = nterms,
                       verbose = FALSE))

  test_that("\"vsel\" object as input leads to correct output structure", {
    for (i in fam_nms) {
      p <- project(vs_list[[i]], nterms = 0:nterms)
      expect_type(p, "list")
      expect_length(p, nterms + 1)
      expect_true(.is_proj_list(p), info = i)

      prjdraw_weights <- p[[1]]$weights
      for (j in seq_along(p)) {
        expect_s3_class(p[[!!j]], "projection")
        # Check the names using `ignore.order = FALSE` because an incorrect
        # order would mean that the documentation of project()'s return value
        # would have to be updated:
        expect_named(p[[!!j]], projection_nms, info = i)
        # Number of projected draws should be equal to the default of `ndraws`
        # (note that more extensive tests for as.matrix.projection() may be
        # found in "test_as_matrix.R"):
        expect_length(p[[!!j]]$sub_fit, ndraws_default)
        expect_length(p[[!!j]]$weights, ndraws_default)
        expect_length(p[[!!j]]$dis, ndraws_default)
        expect_identical(NROW(as.matrix(p[[!!j]])), ndraws_default, info = i)
        # The j-th element should have j solution terms (usually excluding the
        # intercept, but counting it for `j == 1`):
        expect_length(p[[!!j]]$solution_terms, max(j - 1, 1))
        # Same check, but using count_terms_chosen():
        expect_equal(count_terms_chosen(p[[!!j]]$solution_terms), !!j, info = i)
        expect_identical(p[[!!j]]$family, vs_list[[!!i]]$family, info = i)
        # All submodels should use the same clustering:
        expect_identical(p[[!!j]]$weights, prjdraw_weights, info = i)
      }
      # kl should be non-increasing on training data
      klseq <- sapply(p, function(x) sum(x$kl))
      # Remove intercept from the comparison:
      klseq <- klseq[-1]
      expect_identical(klseq, cummin(klseq), info = i)
      ### Check with tolerance:
      # expect_true(all(diff(klseq) - 1e-1 < 0), info = i)
      ###
    }
  })

  test_that(paste(
    "an error is thrown if object is not of class \"vsel\" and",
    "`solution_terms` is not specified"
  ), {
    expect_error(project(fit_gauss), "is not an object of class \"vsel\"")
  })

  test_that("specifying `nterms` incorrectly leads to an error", {
    expect_error(project(vs_list[[1]], nterms = 1000),
                 "Cannot perform the projection with 1000 variables")
    expect_error(project(vs_list[[1]], nterms = -1),
                 "must contain non-negative values")
    expect_error(project(vs_list[[1]], nterms = "a"),
                 "must contain non-negative values")
    expect_error(project(vs_list[[1]], nterms = df_gauss),
                 "must contain non-negative values")
  })

  test_that("specifying `nterms` correctly leads to correct output structure", {
    for (i in fam_nms) {
      for (nterms_tst in list(NULL, 0, 3, c(1, 3))) {
        p <- project(vs_list[[i]], nclusters = nclusters_tst,
                     nterms = nterms_tst)
        out_size <- if (is.null(nterms_tst)) {
          suggest_size(vs_list[[i]])
        } else {
          nterms_tst
        }
        if (length(out_size) == 1) {
          expect_s3_class(p, "projection")
          expect_named(p, projection_nms, info = i)
          expect_length(p$sub_fit, nclusters_tst)
          expect_length(p$weights, nclusters_tst)
          expect_length(p$dis, nclusters_tst)
          SW(nprjdraws <- NROW(as.matrix(p)))
          expect_identical(nprjdraws, nclusters_tst, info = i)
          expect_length(p$solution_terms, max(out_size, 1))
          expect_equal(count_terms_chosen(p$solution_terms) - 1, out_size,
                       info = i)
        } else {
          expect_type(p, "list")
          expect_length(p, length(out_size))
          expect_true(.is_proj_list(p), info = i)

          prjdraw_weights <- p[[1]]$weights
          for (j in seq_along(p)) {
            expect_s3_class(p[[!!j]], "projection")
            expect_named(p[[!!j]], projection_nms, info = i)
            expect_length(p[[!!j]]$sub_fit, nclusters_tst)
            expect_length(p[[!!j]]$weights, nclusters_tst)
            expect_length(p[[!!j]]$dis, nclusters_tst)
            SW(nprjdraws <- NROW(as.matrix(p[[!!j]])))
            expect_identical(nprjdraws, nclusters_tst, info = i)
            expect_length(p[[!!j]]$solution_terms, max(out_size[j], 1))
            expect_equal(count_terms_chosen(p[[!!j]]$solution_terms) - 1,
                         out_size[!!j], info = i)
            expect_identical(p[[!!j]]$family, vs_list[[!!i]]$family, info = i)
            expect_identical(p[[!!j]]$weights, prjdraw_weights, info = i)
          }
        }
      }
    }
  })

  test_that(paste(
    "specifying `solution_terms` incorrectly leads to a warning or an error"
  ), {
    expect_error(project(fit_list[[1]], nclusters = nclusters_tst,
                         solution_terms = NULL),
                 "is not an object of class \"vsel\"")
    for (solterms_tsttmp in list(2, 1:3, "1", list(c("x.3", "x.5"),
                                                   c("x.2", "x.4")))) {
      expect_warning(
        p <- project(fit_list[[1]], nclusters = nclusters_tst,
                     solution_terms = solterms_tsttmp),
        paste("At least one element of `solution_terms` could not be found",
              "among the terms in the reference model"),
        info = as.character(solterms_tsttmp)
      )
      expect_s3_class(p, "projection")
      expect_named(p, projection_nms, info = solterms_tsttmp)
      expect_length(p$sub_fit, nclusters_tst)
      expect_length(p$weights, nclusters_tst)
      expect_length(p$dis, nclusters_tst)
      SW(nprjdraws <- NROW(as.matrix(p)))
      expect_identical(nprjdraws, nclusters_tst, info = solterms_tsttmp)
      expect_identical(p$solution_terms, "1")
    }
  })

  test_that(paste(
    "specifying `solution_terms` correctly leads to correct output structure"
  ), {
    for (i in fam_nms) {
      if (i == "binom") {
        # For the binomial family with > 1 trials, we expect a warning (see
        # GitHub issue #136):
        warn_prj_expect <- paste("Using formula\\(x\\) is deprecated when x",
                                 "is a character vector of length > 1")
      } else {
        warn_prj_expect <- NA
      }
      for (solterms_tsttmp in list(character(), "x.3", c("x.2", "x.4"))) {
        expect_warning(
          p <- project(fit_list[[i]], nclusters = nclusters_tst,
                       solution_terms = solterms_tsttmp),
          warn_prj_expect,
          info = c(i, as.character(solterms_tsttmp))
        )
        expect_s3_class(p, "projection")
        expect_named(p, projection_nms, info = i)
        expect_length(p$sub_fit, nclusters_tst)
        expect_length(p$weights, nclusters_tst)
        expect_length(p$dis, nclusters_tst)
        SW(nprjdraws <- NROW(as.matrix(p)))
        expect_identical(nprjdraws, nclusters_tst, info = i)
        solterms_out <- if (length(solterms_tsttmp) == 0) {
          "1"
        } else {
          solterms_tsttmp
        }
        expect_identical(p$solution_terms, solterms_out)
      }
    }
  })

  test_that("specifying `ndraws` incorrectly leads to an error", {
    i <- "gauss"
    expect_error(project(fit_list[[i]],
                         ndraws = NULL,
                         solution_terms = solterms_tst),
                 "^!is\\.null\\(ndraws\\) is not TRUE$", info = i)
  })

  test_that(paste(
    "specifying `ndraws` and/or `nclusters` too big causes them to be cut off",
    "at the number of posterior draws in the reference model"
  ), {
    i <- "gauss"
    S <- nrow(as.matrix(fit_list[[i]]))
    for (ndraws_tsttmp in list(S + 1L)) {
      for (nclusters_tsttmp in list(NULL, S + 1L)) {
        tstsetup <- unlist(nlist(i, ndraws_tsttmp, nclusters_tsttmp))
        p <- project(fit_list[[i]],
                     ndraws = ndraws_tsttmp,
                     nclusters = nclusters_tsttmp,
                     solution_terms = solterms_tst)
        expect_s3_class(p, "projection")
        expect_named(p, projection_nms, info = tstsetup)
        nprjdraws_out <- S
        nprjdraws_sub_fit <- if (nprjdraws_out == 1) {
          length(sub_fit_nms)
        } else {
          nprjdraws_out
        }
        expect_length(p$sub_fit, nprjdraws_sub_fit)
        expect_length(p$weights, nprjdraws_out)
        expect_length(p$dis, nprjdraws_out)
        SW(nprjdraws <- NROW(as.matrix(p)))
        expect_identical(nprjdraws, nprjdraws_out, info = tstsetup)
        solterms_out <- if (length(solterms_tst) == 0) "1" else solterms_tst
        expect_identical(p$solution_terms, solterms_out)
        if (nprjdraws_out == 1) {
          expect_identical(p$weights, 1, info = tstsetup)
        }
      }
    }
  })

  test_that(paste(
    "specifying `ndraws` and/or `nclusters` correctly leads to correct output",
    "structure"
  ), {
    for (i in fam_nms) {
      if (i == "binom") {
        # For the binomial family with > 1 trials, we expect a warning (see
        # GitHub issue #136):
        warn_prj_expect <- paste("Using formula\\(x\\) is deprecated when x",
                                 "is a character vector of length > 1")
      } else {
        warn_prj_expect <- NA
      }
      for (ndraws_tsttmp in list(1L, 20L, 21L)) {
        for (nclusters_tsttmp in list(NULL, 1L, 2L, 3L)) {
          tstsetup <- unlist(nlist(i, ndraws_tsttmp, nclusters_tsttmp))
          expect_warning(
            p <- project(fit_list[[i]],
                         ndraws = ndraws_tsttmp,
                         nclusters = nclusters_tsttmp,
                         solution_terms = solterms_tst),
            warn_prj_expect,
            info = tstsetup
          )
          expect_s3_class(p, "projection")
          expect_named(p, projection_nms, info = tstsetup)
          nprjdraws_out <- if (!is.null(nclusters_tsttmp)) {
            nclusters_tsttmp
          } else {
            ndraws_tsttmp
          }
          nprjdraws_sub_fit <- if (nprjdraws_out == 1) {
            length(sub_fit_nms)
          } else {
            nprjdraws_out
          }
          expect_length(p$sub_fit, nprjdraws_sub_fit)
          expect_length(p$weights, nprjdraws_out)
          expect_length(p$dis, nprjdraws_out)
          SW(nprjdraws <- NROW(as.matrix(p)))
          expect_identical(nprjdraws, nprjdraws_out, info = tstsetup)
          solterms_out <- if (length(solterms_tst) == 0) "1" else solterms_tst
          expect_identical(p$solution_terms, solterms_out)
          if (nprjdraws_out == 1) {
            expect_identical(p$weights, 1, info = tstsetup)
          }
        }
      }
    }
  })

  test_that("specifying `seed` correctly leads to reproducible results", {
    i <- "gauss"
    p1 <- project(fit_list[[i]],
                  nclusters = nclusters_tst,
                  solution_terms = solterms_tst,
                  seed = seed_tst)
    p2 <- project(fit_list[[i]],
                  nclusters = nclusters_tst,
                  solution_terms = solterms_tst,
                  seed = seed_tst + 1L)
    p3 <- project(fit_list[[i]],
                  nclusters = nclusters_tst,
                  solution_terms = solterms_tst,
                  seed = seed_tst)
    p4 <- project(fit_list[[i]],
                  nclusters = nclusters_tst,
                  solution_terms = solterms_tst)

    # Expected equality:
    expect_true(isTRUE(all.equal(p1, p3)), info = i)
    # The resulting objects are even identical when ignoring the environments of
    # functions:
    expect_identical(p1, p3, info = i, ignore.environment = TRUE)

    # Expected inequality:
    expect_false(isTRUE(all.equal(p1, p2)), info = i)
    expect_false(isTRUE(all.equal(p1, p4)), info = i)
    expect_false(isTRUE(all.equal(p2, p3)), info = i)
    expect_false(isTRUE(all.equal(p2, p4)), info = i)
    expect_false(isTRUE(all.equal(p3, p4)), info = i)
  })

  test_that(paste(
    "projecting the reference model onto the full model (i.e.,",
    "itself) does not change results on average (even though this is not",
    "guaranteed; see the comments)"
  ), {
    # NOTE: Projecting the reference model onto the full model (i.e., itself)
    # does not necessarily have to give results close to the reference model's
    # since in contrast to the reference model, the projection is "fitting to
    # the fit" of the reference model, not to the observed response. The fact
    # that the tolerance for the Gaussian reference model needs to be
    # increased here (see below) might be an indicator for this inequality.
    tol <- setNames(rep(1e-3, length(fit_list)), fam_nms)
    tol["gauss"] <- 0.25

    for (i in fam_nms) {
      fit <- fit_list[[i]]
      draws <- as.data.frame(fit)
      alpha_ref <- draws[, "(Intercept)"]
      beta_ref <- draws[, 1 + seq_len(nterms), drop = FALSE]
      S <- nrow(draws)
      SW(vs <- varsel(fit))
      proj <- project(vs,
                      solution_terms = vs$solution_terms[1:nterms],
                      seed = seed, ndraws = S)

      # test alpha and beta
      coefs <- as.matrix(proj)
      dalpha <- abs(mean(coefs[, 1]) - mean(alpha_ref))
      order <- match(colnames(fit_list[[i]]$data), proj$solution_terms)
      order <- order[!is.na(order)]
      dbeta <- max(abs(colMeans(coefs[, -1, drop = FALSE][, order]) -
                         colMeans(beta_ref)))
      expect_lt(dalpha, tol[!!i])
      expect_lt(dbeta, tol[!!i])
    }
  })

  test_that("project() works as expected from a vsel object", {
    SW({
      cvs <- cv_varsel(fit_binom,
                       nterms_max = 3, verbose = FALSE, ndraws = ndraws,
                       ndraws_pred = ndraws_pred)
      p <- project(cvs, nterms = 3)
    })
    expect_length(p$solution_terms, 3)
  })
}
