context("datafit")

# Setup -------------------------------------------------------------------

if (!requireNamespace("glmnet", quietly = TRUE)) {
  stop("Package \"glmnet\" needed for this test to work. Please install it.",
       call. = FALSE)
}

.extrmoddat_datafit <- function(object, newdata = NULL, wrhs = NULL,
                                orhs = NULL, resp_form = NULL) {
  if (is.null(newdata)) {
    newdata <- object$data
  }

  if (inherits(wrhs, "formula")) {
    weights <- eval_rhs(wrhs, newdata)
  } else if (is.null(wrhs)) {
    weights <- newdata$wobs_col
  } else {
    weights <- wrhs
  }

  if (inherits(orhs, "formula")) {
    offset <- eval_rhs(orhs, newdata)
  } else if (is.null(orhs)) {
    offset <- newdata$offs_col
  } else {
    offset <- orhs
  }

  if (inherits(resp_form, "formula")) {
    y <- eval_rhs(resp_form, newdata)
  } else {
    y <- NULL
  }

  return(nlist(y, weights, offset))
}

## Reference model --------------------------------------------------------
## (actually "datafit"s)

# For the binomial family with > 1 trials, we currently expect the warning
# "Using formula(x) is deprecated when x is a character vector of length > 1"
# (see GitHub issue #136), so temporarily wrap the following call in SW():
SW(datafits <- lapply(mod_nms, function(mod_nm) {
  lapply(fam_nms, function(fam_nm) {
    formul_crr <- fits[[mod_nm]][[fam_nm]]$formula
    extrmoddat <- function(object, newdata = NULL, wrhs = NULL, orhs = NULL,
                           extract_y = TRUE) {
      resp_form <- if (!extract_y) NULL else lhs(formul_crr)
      if (is.null(newdata)) {
        newdata <- dat
      }
      args <- nlist(object, newdata, wrhs, orhs, resp_form)
      return(do.call(.extrmoddat_datafit, args))
    }
    return(init_refmodel(
      object = NULL,
      data = dat,
      formula = formul_crr,
      family = get(paste0("f_", fam_nm)),
      extract_model_data = extrmoddat
    ))
  })
}))

## Variable selection -----------------------------------------------------

### varsel() --------------------------------------------------------------

args_vs_datafit <- args_vs
args_vs_datafit <- lapply(args_vs_datafit, function(args_vs_i) {
  return(args_vs_i[setdiff(names(args_vs_i),
                           c("ndraws", "nclusters",
                             "ndraws_pred", "nclusters_pred"))])
})

if (run_vs) {
  vss_datafit <- lapply(args_vs_datafit, function(args_vs_i) {
    do.call(varsel, c(
      list(object = datafits[[args_vs_i$mod_nm]][[args_vs_i$fam_nm]]),
      args_vs_i[setdiff(names(args_vs_i), c("mod_nm", "fam_nm"))]
    ))
  })
}

### cv_varsel() -----------------------------------------------------------

# (PSIS-)LOO CV is not possible for `"datafit"`s, so only use K-fold CV:
args_cvvs_datafit <- args_cvvs[
  grep("kfold", names(args_cvvs), value = TRUE, invert = TRUE)
]
args_cvvs_datafit <- lapply(args_cvvs_datafit, "c",
                            list(cv_method = "kfold", K = 2))
names(args_cvvs_datafit) <- gsub("default_cvmeth", "kfold",
                                 names(args_cvvs_datafit))
args_cvvs_datafit <- lapply(args_cvvs_datafit, function(args_cvvs_i) {
  return(args_cvvs_i[setdiff(names(args_cvvs_i),
                             c("ndraws", "nclusters",
                               "ndraws_pred", "nclusters_pred"))])
})

if (run_cvvs) {
  cvvss_datafit <- lapply(args_cvvs_datafit, function(args_cvvs_i) {
    do.call(cv_varsel, c(
      list(object = datafits[[args_cvvs_i$mod_nm]][[args_cvvs_i$fam_nm]]),
      args_cvvs_i[setdiff(names(args_cvvs_i), c("mod_nm", "fam_nm"))]
    ))
  })
}

## Projection -------------------------------------------------------------

### From varsel() ---------------------------------------------------------

args_prj_vs_datafit <- args_prj_vs
args_prj_vs_datafit <- lapply(args_prj_vs_datafit, function(args_prj_vs_i) {
  return(args_prj_vs_i[setdiff(names(args_prj_vs_i), c("ndraws", "nclusters"))])
})

if (run_vs) {
  prjs_vs_datafit <- lapply(args_prj_vs_datafit, function(args_prj_vs_i) {
    do.call(project, c(
      list(object = vss_datafit[[args_prj_vs_i$tstsetup_vsel]]),
      args_prj_vs_i[setdiff(names(args_prj_vs_i), c("tstsetup_vsel"))]
    ))
  })
}

## Prediction -------------------------------------------------------------

### From "proj_list" ------------------------------------------------------

pls_vs_datafit <- lapply(prjs_vs_datafit, proj_linpred)
pps_vs_datafit <- lapply(prjs_vs_datafit, proj_predict, .seed = seed2_tst)

# Tests (projpred only) ---------------------------------------------------

## Reference model --------------------------------------------------------

test_that("predict.refmodel(): error if `object` is of class \"datafit\"", {
  for (mod_nm in mod_nms) {
    for (fam_nm in fam_nms) {
      expect_error(predict(datafits[[mod_nm]][[fam_nm]], newdata = dat),
                   "^Cannot make predictions with data reference only\\.$")
    }
  }
})

## Variable selection -----------------------------------------------------

test_that(paste(
  "varsel(): `object` of class \"datafit\", `method`, and `nterms_max` work"
), {
  skip_if_not(run_vs)
  for (tstsetup in names(vss_datafit)) {
    mod_crr <- args_vs_datafit[[tstsetup]]$mod
    fam_crr <- args_vs_datafit[[tstsetup]]$fam
    meth_exp_crr <- args_vs_datafit[[tstsetup]]$method
    if (is.null(meth_exp_crr)) {
      meth_exp_crr <- ifelse(mod_crr == "glm", "L1", "forward")
    }
    vsel_tester(
      vss_datafit[[tstsetup]],
      from_datafit = TRUE,
      refmod_expected = datafits[[mod_crr]][[fam_crr]],
      solterms_len_expected = args_vs_datafit[[tstsetup]]$nterms_max,
      method_expected = meth_exp_crr,
      nclusters_expected = 1L,
      nclusters_pred_expected = 1L,
      info_str = tstsetup
    )
  }
})

test_that(paste(
  "cv_varsel(): `object` of class \"datafit\", `method`, `cv_method`, and",
  "`nterms_max` work"
), {
  skip_if_not(run_cvvs)
  for (tstsetup in names(cvvss_datafit)) {
    mod_crr <- args_cvvs_datafit[[tstsetup]]$mod
    fam_crr <- args_cvvs_datafit[[tstsetup]]$fam
    meth_exp_crr <- args_cvvs_datafit[[tstsetup]]$method
    if (is.null(meth_exp_crr)) {
      meth_exp_crr <- ifelse(mod_crr == "glm", "L1", "forward")
    }
    vsel_tester(
      cvvss_datafit[[tstsetup]],
      with_cv = TRUE,
      from_datafit = TRUE,
      refmod_expected = datafits[[mod_crr]][[fam_crr]],
      solterms_len_expected = args_cvvs_datafit[[tstsetup]]$nterms_max,
      method_expected = meth_exp_crr,
      cv_method_expected = "kfold",
      valsearch_expected = args_cvvs_datafit[[tstsetup]]$validate_search,
      nclusters_expected = 1L,
      nclusters_pred_expected = 1L,
      info_str = tstsetup
    )
  }
})

test_that("summary.vsel(): error if `baseline = \"ref\"` and `deltas = TRUE`", {
  for (tstsetup in names(vss_datafit)[1]) {
    expect_error(
      summary(vss_datafit[[tstsetup]], baseline = "ref", deltas = TRUE),
      paste("^Cannot use deltas = TRUE and baseline = 'ref' when there is no",
            "reference model\\.$"),
      info = tstsetup
    )
  }
})

## Projection -------------------------------------------------------------

test_that("project(): error if `object` is of class \"datafit\"", {
  tstsetups <- grep("\\.solterms_x.*\\.clust$", names(args_prj), value = TRUE)
  for (tstsetup in tstsetups) {
    args_prj_i <- args_prj[[tstsetup]]
    expect_error(
      do.call(project, c(
        list(object = datafits[[args_prj_i$mod_nm]][[args_prj_i$fam_nm]]),
        args_prj_i[setdiff(names(args_prj_i), c("mod_nm", "fam_nm"))]
      )),
      paste("^no applicable method for 'predict' applied to an object of class",
            "\"NULL\"$"),
      info = tstsetup
    )
  }
})

test_that(paste(
  "project(): `object` of class \"vsel\" (created by varsel() applied to an",
  "`object` of class \"datafit\"), `nclusters`, and `nterms` work"
), {
  skip_if_not(run_vs)
  for (tstsetup in names(prjs_vs_datafit)) {
    tstsetup_vs <- args_prj_vs_datafit[[tstsetup]]$tstsetup_vsel
    stopifnot(length(tstsetup_vs) > 0)
    mod_crr <- args_vs_datafit[[tstsetup_vs]]$mod_nm
    fam_crr <- args_vs_datafit[[tstsetup_vs]]$fam_nm
    nterms_crr <- args_prj_vs_datafit[[tstsetup]]$nterms
    if (is.null(nterms_crr)) {
      nterms_crr <- vss_datafit[[tstsetup_vs]]$suggested_size
    }
    if (length(nterms_crr) == 1) {
      solterms_expected_crr <- vss_datafit[[tstsetup_vs]]$solution_terms[
        seq_len(nterms_crr)
      ]
      projection_tester(
        prjs_vs_datafit[[tstsetup]],
        solterms_expected = solterms_expected_crr,
        nprjdraws_expected = 1L,
        p_type_expected = TRUE,
        from_datafit = TRUE,
        info_str = tstsetup
      )
    } else {
      proj_list_tester(
        prjs_vs_datafit[[tstsetup]],
        len_expected = length(nterms_crr),
        is_seq = all(diff(nterms_crr) == 1),
        from_datafit = TRUE,
        info_str = tstsetup,
        nprjdraws_expected = 1L,
        p_type_expected = TRUE,
        fam_expected = vss_datafit[[tstsetup_vs]]$family,
        prjdraw_weights_expected = prjs_vs_datafit[[tstsetup]][[1]]$weights
      )
    }
  }
})

test_that(paste(
  "proj_linpred(): `object` of (informal) class \"proj_list\" (based on",
  "varsel()) works"
), {
  skip_if_not(run_vs)
  for (tstsetup in names(prjs_vs_datafit)) {
    tstsetup_vs <- args_prj_vs_datafit[[tstsetup]]$tstsetup_vsel
    nterms_crr <- args_prj_vs_datafit[[tstsetup]]$nterms
    if (is.null(nterms_crr)) {
      nterms_crr <- vss_datafit[[tstsetup_vs]]$suggested_size
    }
    pl_tester(pls_vs_datafit[[tstsetup]],
              len_expected = length(nterms_crr),
              nprjdraws_expected = 1L,
              info_str = tstsetup)
    pl_with_args <- proj_linpred(prjs_vs_datafit[[tstsetup]],
                                 newdata = head(dat, tail(nobsv_tst, 1)),
                                 weightsnew = ~ wobs_col,
                                 offsetnew = ~ offs_col,
                                 filter_nterms = nterms_crr[1])
    pl_tester(pl_with_args,
              len_expected = 1L,
              nprjdraws_expected = 1L,
              nobsv_expected = tail(nobsv_tst, 1),
              info_str = paste("with_args", tstsetup, sep = "__"))
  }
})

test_that(paste(
  "proj_predict(): `object` of (informal) class \"proj_list\" (based on",
  "varsel()) works"
), {
  skip_if_not(run_vs)
  for (tstsetup in names(prjs_vs_datafit)) {
    tstsetup_vs <- args_prj_vs_datafit[[tstsetup]]$tstsetup_vsel
    nterms_crr <- args_prj_vs_datafit[[tstsetup]]$nterms
    if (is.null(nterms_crr)) {
      nterms_crr <- vss_datafit[[tstsetup_vs]]$suggested_size
    }
    pp_tester(pps_vs_datafit[[tstsetup]],
              len_expected = length(nterms_crr),
              info_str = tstsetup)
    pp_with_args <- proj_predict(
      prjs_vs_datafit[[tstsetup]],
      newdata = head(dat, tail(nobsv_tst, 1)),
      weightsnew = ~ wobs_col,
      offsetnew = ~ offs_col,
      filter_nterms = nterms_crr[1],
      nresample_clusters = tail(nresample_clusters_tst, 1),
      .seed = seed2_tst
    )
    pp_tester(pp_with_args,
              len_expected = 1L,
              nprjdraws_out_expected = tail(nresample_clusters_tst, 1),
              nobsv_expected = tail(nobsv_tst, 1),
              info_str = paste("with_args", tstsetup, sep = "__"))
  }
})

# Comparison with glmnet --------------------------------------------------

# below are some tests that check Lasso solution computed with varsel is the
# same as that of glmnet. (notice that glm_ridge and glm_elnet are already
# tested separately, so these would only check that the results do not change
# due to varsel/cv_varsel etc.)

set.seed(1235)
n <- 100
nterms <- 10
x <- matrix(rnorm(n * nterms, 0, 1), n, nterms)
b <- seq(0, 1, length.out = nterms)
dis <- runif(1, 0.3, 0.5)
weights <- sample(1:4, n, replace = TRUE)
offset <- 0.1 * rnorm(n)

fams <- list(gaussian(), binomial(), poisson())
x_list <- lapply(fams, function(fam) x)
y_list <- lapply(fams, function(fam) {
  if (fam$family == "gaussian") {
    y <- rnorm(n, x %*% b, 0.5)
    weights <- NULL
    y_glmnet <- y
  } else if (fam$family == "binomial") {
    y <- rbinom(n, weights, fam$linkinv(x %*% b))
    ## y <- y / weights
    ## different way of specifying binomial y for glmnet
    y_glmnet <- cbind(1 - y / weights, y / weights)
    weights <- weights
  } else if (fam$family == "poisson") {
    y <- rpois(n, fam$linkinv(x %*% b))
    y_glmnet <- y
    weights <- NULL
  }
  nlist(y, y_glmnet, weights)
})

median_lasso_preds <- list(
  c(0.2774068, 0.2857059, 0.2878935, 0.2813947, 0.2237729,
    0.2895152, 0.3225808, 0.3799348),
  c(0.009607217, 0.015400719, -0.017591445, -0.009711566,
    -0.023867036, -0.038964983, -0.036081074, -0.045065655),
  c(1.8846845, 1.8830678, 1.8731548, 1.4232035, 0.9960167,
    0.9452660, 0.6216253, 0.5856283)
)

solution_terms_lasso <- list(
  c(10, 9, 6, 8, 7, 5, 4, 3, 1, 2),
  c(10, 9, 8, 6, 7, 5, 3, 4, 2, 1),
  c(9, 10, 6, 7, 3, 5, 2, 4, 3, 1)
)

test_that(paste(
  "L1-projection with data reference gives the same results as",
  "Lasso from glmnet."
), {

  extract_model_data <- function(object, newdata = NULL, wrhs = NULL,
                                 orhs = NULL, extract_y = FALSE) {
    if (!is.null(object)) {
      formula <- formula(object)
      tt <- extract_terms_response(formula)
      response_name <- tt$response
    } else {
      response_name <- NULL
    }

    if (is.null(newdata)) {
      newdata <- object$data
    }

    resp_form <- NULL
    if (is.null(object)) {
      if ("weights" %in% colnames(newdata)) {
        wrhs <- ~ weights
      }
      if ("offset" %in% colnames(newdata)) {
        orhs <- ~ offset
      }
      if ("y" %in% colnames(newdata)) {
        resp_form <- ~ y
      }
    }

    args <- nlist(object, newdata, wrhs, orhs, resp_form)
    return(do_call(.extract_model_data, args))
  }

  for (i in seq_along(fams)) {
    x <- x_list[[i]]
    y <- y_list[[i]]$y
    y_glmnet <- y_list[[i]]$y_glmnet
    fam <- fams[[i]]
    weights <- y_list[[i]]$weights
    if (is.null(weights)) {
      weights <- rep(1, NROW(y))
    }

    lambda_min_ratio <- 1e-7
    nlambda <- 1500

    df <- data.frame(y = y, x = x, weights = weights)
    formula <- y ~ x.1 + x.2 + x.3 + x.4 + x.5 + x.6 + x.7 + x.8 + x.9 + x.10
    # Lasso solution with projpred
    ref <- init_refmodel(
      object = NULL, data = df, formula = formula,
      family = fam, extract_model_data = extract_model_data
    )
    SW({
      vs <- varsel(ref,
                   method = "l1", lambda_min_ratio = lambda_min_ratio,
                   nlambda = nlambda, thresh = 1e-12)
    })
    pred1 <- proj_linpred(vs,
                          newdata = data.frame(x = x, offset = offset,
                                               weights = weights),
                          nterms = 0:nterms, transform = FALSE,
                          offsetnew = ~offset)

    # compute the results for the Lasso
    lasso <- glmnet::glmnet(x, y_glmnet,
                            family = fam$family, weights = weights,
                            offset = offset,
                            lambda.min.ratio = lambda_min_ratio,
                            nlambda = nlambda, thresh = 1e-12)
    solution_terms <- predict(lasso, type = "nonzero", s = lasso$lambda)
    nselected <- sapply(solution_terms, function(e) length(e))
    lambdainds <- sapply(unique(nselected), function(nterms) {
      max(which(nselected == nterms))
    })
    lambdaval <- lasso$lambda[lambdainds]
    ## pred2 <- predict(lasso,
    ##   newx = x, type = "link", s = lambdaval,
    ##   newoffset = offset
    ## )

    # check that the predictions agree (up to nterms-2 only, because glmnet
    # terminates the coefficient path computation too early for some reason)
    for (j in 1:(nterms - 2)) {
      expect_true(median(pred1[[j]]$pred) - median_lasso_preds[[i]][j] < 3e-1)
    }

    # check that the coefficients are similar
    ind <- match(vs$solution_terms, setdiff(split_formula(formula), "1"))
    if (Sys.getenv("NOT_CRAN") == "true") {
      betas <- sapply(vs$search_path$sub_fits, function(x) x$beta %||% 0)
      delta <- sapply(seq_len(nterms), function(i) {
        abs(t(betas[[i + 1]]) - lasso$beta[ind[1:i], lambdainds[i + 1]])
      })
      expect_true(median(unlist(delta)) < 6e-2)
      expect_true(median(abs(sapply(vs$search_path$sub_fits, function(x) {
        x$alpha
      }) - lasso$a0[lambdainds])) < 1.5e-1)
    } else {
      expect_true(sum(ind == solution_terms_lasso[[i]]) >= nterms / 2)
    }
  }
})
