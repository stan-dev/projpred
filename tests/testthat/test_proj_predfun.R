context("proj_predfun")

test_that("repair_re() works for GLMMs", {
  # Setup -------------------------------------------------------------------

  if (exists(".Random.seed", envir = .GlobalEnv)) {
    rng_old <- get(".Random.seed", envir = .GlobalEnv)
  }
  set.seed(seed2_tst)

  # Data --------------------------------------------------------------------

  dat <- matrix(rnorm(41 * 2, sd = 0.4), ncol = 2,
                dimnames = list(NULL, paste0("X", 1:2)))
  dat <- data.frame(dat)

  nsbj <- 15L
  dat$sbj <- gl(n = nsbj,
                k = floor(nrow(dat) / nsbj),
                length = nrow(dat),
                labels = paste0("subj", seq_len(nsbj)))
  dat$sbj <- as.integer(dat$sbj)

  ngrp <- 8L
  dat$grp <- gl(n = ngrp,
                k = floor(nrow(dat) / ngrp),
                length = nrow(dat),
                labels = paste0("gr", seq_len(ngrp)))
  dat$grp <- sample(dat$grp)

  nXf <- 3L
  dat$Xf <- gl(n = nXf,
               k = floor(nrow(dat) / nXf),
               length = nrow(dat),
               labels = paste0("l", seq_len(nXf)))
  dat$Xf <- as.character(dat$Xf)
  dat$Xf <- sample(dat$Xf)

  sbj_icpts_truth <- rnorm(nsbj, sd = 6)
  grp_icpts_truth <- rnorm(ngrp, sd = 6)
  icpt <- -4.2
  dat$y <- icpt + sbj_icpts_truth[dat$sbj] + grp_icpts_truth[dat$grp]
  dat$y <- rnorm(nrow(dat), mean = dat$y, sd = 4)

  ## New --------------------------------------------------------------------

  dat_new <- matrix(rnorm(3 * 2, sd = 0.4), ncol = 2,
                    dimnames = list(NULL, paste0("X", 1:2)))
  dat_new <- data.frame(dat_new)
  dat_new$sbj <- c(nsbj, rep(nsbj + 1L, 2))
  dat_new$grp <- c(ngrp, rep(ngrp + 1L, 2))
  # as.factor() could be applied here to test the case of a group variable which
  # is a `factor`:
  dat_new$grp <- paste0("gr", dat_new$grp)
  # as.factor() could be applied here to test the case of a categorical
  # predictor variable which is a `factor`:
  dat_new$Xf <- paste0("l", nXf)

  # Fit with lme4 -----------------------------------------------------------

  lmm_fit <- lme4::lmer(y ~ X1 + X2 + Xf + (1 | sbj) + (1 + X1 + Xf | grp),
                        data = dat)

  # Extract estimated parameters --------------------------------------------

  lmm_fixef <- lme4::fixef(lmm_fit)
  lmm_ranef <- lme4::ranef(lmm_fit, condVar = FALSE)
  lmm_VarCorr <- lme4::VarCorr(lmm_fit)

  # Predict -----------------------------------------------------------------

  lpreds_orig <- predict(lmm_fit, newdata = dat_new, allow.new.levels = TRUE)

  fixnms_b <- c(paste0("X", 1:2), paste0("Xfl", 2:3))
  dat_new_ch <- within(dat_new, {
    Xfl2 <- as.integer(Xf == "l2")
    Xfl3 <- as.integer(Xf == "l3")
  })
  ranslopes_orig <- list(
    X1 = c(lmm_ranef$grp[as.character(dat_new_ch$grp[1]), "X1"], rep(0, 2)),
    X2 = rep(0, 3),
    Xfl2 = c(lmm_ranef$grp[as.character(dat_new_ch$grp[1]), "Xfl2"], rep(0, 2)),
    Xfl3 = c(lmm_ranef$grp[as.character(dat_new_ch$grp[1]), "Xfl3"], rep(0, 2))
  )
  lpreds_orig_man <- lmm_fixef["(Intercept)"] +
    c(lmm_ranef$sbj[as.character(dat_new_ch$sbj[1]), "(Intercept)"], rep(0, 2)) +
    c(lmm_ranef$grp[as.character(dat_new_ch$grp[1]), "(Intercept)"], rep(0, 2)) +
    sapply(seq_len(nrow(dat_new_ch)), function(i) {
      drop(as.matrix(dat_new_ch[i, fixnms_b, drop = FALSE]) %*%
             (lmm_fixef[fixnms_b] + sapply(ranslopes_orig[fixnms_b], "[", i)))
    })
  expect_equal(lpreds_orig, lpreds_orig_man, tolerance = 1e-15)

  ## With repair_re() -------------------------------------------------------

  set.seed(seed3_tst)
  lpreds <- subprd(list(lmm_fit), newdata = dat_new)
  expect_identical(dim(lpreds), c(nrow(dat_new), 1L))
  lpreds <- lpreds[, 1]

  set.seed(seed3_tst)
  ranmulti_sbj <- mvtnorm::rmvnorm(
    n = 1, sigma = lmm_VarCorr$sbj[, , drop = FALSE], checkSymmetry = FALSE
  )
  ranmulti_grp <- mvtnorm::rmvnorm(
    n = 1, sigma = lmm_VarCorr$grp[, , drop = FALSE], checkSymmetry = FALSE
  )
  ranslopes <- list(
    X1 = c(lmm_ranef$grp[as.character(dat_new_ch$grp[1]), "X1"],
           rep(ranmulti_grp[1, 2], 2)),
    X2 = rep(0, 3),
    Xfl2 = c(lmm_ranef$grp[as.character(dat_new_ch$grp[1]), "Xfl2"],
             rep(ranmulti_grp[1, 3], 2)),
    Xfl3 = c(lmm_ranef$grp[as.character(dat_new_ch$grp[1]), "Xfl3"],
             rep(ranmulti_grp[1, 4], 2))
  )
  lpreds_man <- lmm_fixef["(Intercept)"] +
    c(lmm_ranef$sbj[as.character(dat_new_ch$sbj[1]), "(Intercept)"],
      rep(ranmulti_sbj[1, 1], 2)) +
    c(lmm_ranef$grp[as.character(dat_new_ch$grp[1]), "(Intercept)"],
      rep(ranmulti_grp[1, 1], 2)) +
    sapply(seq_len(nrow(dat_new_ch)), function(i) {
      drop(as.matrix(dat_new_ch[i, fixnms_b, drop = FALSE]) %*%
             (lmm_fixef[fixnms_b] + sapply(ranslopes[fixnms_b], "[", i)))
    })
  expect_equal(lpreds, lpreds_man, tolerance = 1e-15)

  # Teardown ----------------------------------------------------------------

  if (exists("rng_old")) assign(".Random.seed", rng_old, envir = .GlobalEnv)
})

test_that("repair_re() works for multilevel cumulative() models", {
  skip_if_not(requireNamespace("brms", quietly = TRUE))

  # Setup -------------------------------------------------------------------

  if (exists(".Random.seed", envir = .GlobalEnv)) {
    rng_old <- get(".Random.seed", envir = .GlobalEnv)
  }
  set.seed(seed2_tst)

  # Data --------------------------------------------------------------------

  data(inhaler, package = "brms")
  inhaler$rating <- as.factor(paste0("rate", inhaler$rating))
  # inhaler$subject <- paste0("subj", inhaler$subject)
  names(inhaler)[names(inhaler) == "subject"] <- "sbj"
  inhaler$period <- paste0("pd", inhaler$period + 0.5)
  ngrp <- 8L
  inhaler$grp <- gl(n = ngrp,
                    k = floor(nrow(inhaler) / ngrp),
                    length = nrow(inhaler),
                    labels = paste0("gr", seq_len(ngrp)))
  inhaler$grp <- sample(inhaler$grp)

  ## New --------------------------------------------------------------------

  inhaler_new <- head(inhaler, 3)
  inhaler_new$carry <- -1
  inhaler_new$sbj <- c(length(unique(inhaler$sbj)),
                       rep(length(unique(inhaler$sbj)) + 1L, 2))
  inhaler_new$grp <- c(ngrp, rep(ngrp + 1L, 2))
  # as.factor() could be applied here to test the case of a group variable which
  # is a `factor`:
  inhaler_new$grp <- paste0("gr", inhaler_new$grp)

  # Fit with ordinal::clmm() ------------------------------------------------

  expect_warning(
    ofit <- ordinal::clmm(
      rating ~ period + carry + treat + (1 | sbj) + (1 + period + carry | grp),
      data = inhaler,
      Hess = FALSE,
      model = FALSE
    ),
    paste("^Using formula\\(x\\) is deprecated when x is a character vector of",
          "length > 1")
  )
  # Needed for the ordinal:::predict.clm() workaround (the value `"negative"` is
  # the default, see `?ordinal::clm.control`):
  ofit$control$sign.location <- "negative"

  # Extract estimated parameters --------------------------------------------

  oparams <- coef(ofit)
  idx_thres <- grep("\\|", names(oparams))
  othres <- oparams[idx_thres]
  ofixef <- oparams[-idx_thres]
  oranef <- ordinal::ranef(ofit)
  oVarCorr <- ordinal::VarCorr(ofit)

  # Predict -----------------------------------------------------------------

  prbs <- predict(structure(ofit, class = c(oldClass(ofit), "clm")),
                  # Need to remove the response to make ordinal:::predict.clm()
                  # return predictions for *all* response categories:
                  newdata = within(inhaler_new, rating <- NULL),
                  type = "prob")$fit
  link_nm <- ofit$link
  prbs <- do.call(rbind, apply(prbs, 1, cumsum, simplify = FALSE))
  prbs <- prbs[, -ncol(prbs), drop = FALSE]
  if (link_nm %in% c("logistic", "logit")) {
    linkfun_raw <- function(x) qlogis(x)
  } else if (link_nm == "probit") {
    linkfun_raw <- function(x) qnorm(x)
  } else if (link_nm == "cloglog") {
    linkfun_raw <- function(x) log(-log1p(-x))
  } else if (link_nm == "cauchit") {
    linkfun_raw <- function(x) qcauchy(x)
  } else {
    stop("Unknown `link_nm`.")
  }
  lpreds_orig <- linkfun_raw(prbs)

  fixnms_b <- c("periodpd1", "carry", "treat")
  inhaler_new_ch <- within(inhaler_new, {
    periodpd1 <- as.integer(factor(period, levels = paste0("pd", 0:1))) - 1L
  })
  lpreds_orig_man <- matrix(othres, nrow = nrow(inhaler_new_ch),
                            ncol = length(othres), byrow = TRUE) -
    drop(as.matrix(inhaler_new_ch[, fixnms_b, drop = FALSE]) %*%
           ofixef[fixnms_b])
  expect_equal(unname(lpreds_orig), lpreds_orig_man, tolerance = 1e-14)

  ## With repair_re() -------------------------------------------------------

  set.seed(seed3_tst)
  lpreds <- subprd_augdat(list(ofit), newdata = inhaler_new)
  expect_identical(dim(lpreds),
                   c(nrow(inhaler_new), nlevels(inhaler$rating) - 1L, 1L))
  lpreds <- lpreds[, , 1]

  set.seed(seed3_tst)
  ranmulti_sbj <- mvtnorm::rmvnorm(
    n = 1, sigma = oVarCorr$sbj[, , drop = FALSE], checkSymmetry = FALSE
  )
  ranmulti_grp <- mvtnorm::rmvnorm(
    n = 1, sigma = oVarCorr$grp[, , drop = FALSE], checkSymmetry = FALSE
  )
  ranslopes <- list(
    periodpd1 = c(oranef$grp[as.character(inhaler_new_ch$grp[1]), "periodpd1"],
                  rep(ranmulti_grp[1, 2], 2)),
    carry = c(oranef$grp[as.character(inhaler_new_ch$grp[1]), "carry"],
              rep(ranmulti_grp[1, 3], 2)),
    treat = rep(0, 3)
  )
  lpreds_man <- matrix(othres, nrow = nrow(inhaler_new_ch),
                       ncol = length(othres), byrow = TRUE) -
    (
      c(oranef$sbj[as.character(inhaler_new_ch$sbj[1]), "(Intercept)"],
        rep(ranmulti_sbj[1, 1], 2)) +
        c(oranef$grp[as.character(inhaler_new_ch$grp[1]), "(Intercept)"],
          rep(ranmulti_grp[1, 1], 2)) +
        sapply(seq_len(nrow(inhaler_new_ch)), function(i) {
          drop(as.matrix(inhaler_new_ch[i, fixnms_b, drop = FALSE]) %*%
                 (ofixef[fixnms_b] + sapply(ranslopes[fixnms_b], "[", i)))
        })
    )
  expect_equal(unname(lpreds), lpreds_man, tolerance = 1e-14)

  # Teardown ----------------------------------------------------------------

  if (exists("rng_old")) assign(".Random.seed", rng_old, envir = .GlobalEnv)
})

# Multilevel categorical models -------------------------------------------

## Setup ------------------------------------------------------------------

# Needed to clean up the workspace afterwards (i.e, after this test file):
ls_bu <- ls()

# Set seed:
if (exists(".Random.seed", envir = .GlobalEnv)) {
  rng_old <- get(".Random.seed", envir = .GlobalEnv)
}
set.seed(seed2_tst)

## Data -------------------------------------------------------------------

data(VA, package = "MASS")
levels(VA$cell) <- paste0("lvl", levels(VA$cell))
VA$cell <- factor(VA$cell,
                  levels = paste0("lvl", seq_len(nlevels(VA$cell))),
                  labels = c("low", "medium", "mid-high", "high"))

nsbj <- 3L
VA$sbj <- gl(n = nsbj,
             k = floor(nrow(VA) / nsbj),
             length = nrow(VA),
             labels = paste0("subj", seq_len(nsbj)))
VA$sbj <- as.integer(VA$sbj)

ngrp <- 8L
VA$grp <- gl(n = ngrp,
             k = floor(nrow(VA) / ngrp),
             length = nrow(VA),
             labels = paste0("gr", seq_len(ngrp)))
VA$grp <- sample(VA$grp)

### New -------------------------------------------------------------------

VA_new <- head(VA, 3)
VA_new$treat[] <- "2"
VA_new$sbj <- c(nsbj, rep(nsbj + 1L, 2))
# VA_new$sbj <- paste0("subj", VA_new$sbj)
VA_new$grp <- c(ngrp, rep(ngrp + 1L, 2))
# as.factor() could be applied here to test the case of a group variable which
# is a `factor`:
VA_new$grp <- paste0("gr", VA_new$grp)

# Replace new group levels in `VA_new` by existing ones (here the first one) to
# be able to run predict.mmblogit() with `conditional = TRUE`:
VA_new_dummy <- VA_new
VA_new_dummy[!VA_new_dummy$grp %in% levels(VA$grp), "grp"] <-
  head(levels(VA$grp), 1)
VA_new_dummy[!VA_new_dummy$sbj %in% unique(VA$sbj), "sbj"] <-
  head(sort(unique(VA$sbj)), 1)

# More objects needed for the following tests:
fixnms_b <- c("treat2", "age", "Karn", "prior10")
VA_new_ch <- within(VA_new, {
  treat2 <- as.integer(factor(treat, levels = paste0("", 1:2))) - 1L
  prior10 <- as.integer(factor(prior, levels = paste0("", c(0, 10)))) - 1L
})
lvl_idx_new1_sbj <- match(VA_new_ch$sbj[1], sort(unique(VA$sbj)))
lvl_idx_new1 <- match(VA_new_ch$grp[1], levels(VA$grp))
lvl_idx1 <- 1L
nlats <- nlevels(VA$cell) - 1L

## Single group variable --------------------------------------------------

test_that(paste(
  "repair_re() works for multilevel categorical() models with only a single",
  "group variable"
), {
  ### Fit with mclogit::mblogit() -------------------------------------------

  warn_expected <- if (packageVersion("mclogit") <= "0.8.7.3") {
    "variable 'prior' is absent, its contrast will be ignored"
  } else {
    NA
  }
  expect_warning(
    out_capt <- capture.output({
      mfit <- mclogit::mblogit(
        formula = cell ~ treat + age + Karn + prior,
        data = VA,
        random = ~ 1 + treat | grp,
        model = FALSE,
        y = FALSE
      )
    }),
    warn_expected
  )
  expect_identical(tail(out_capt, 1), "converged")

  ### Extract estimated parameters ------------------------------------------

  mfixef <- mfit$coefmat
  mranef <- setNames(mfit$random.effects, names(mfit$groups))
  mVarCorr <- mfit$VarCov
  if (packageVersion("mclogit") < "0.9") {
    mVarCorr <- setNames(mVarCorr, names(mfit$groups))
  }

  ### Predict ---------------------------------------------------------------

  # predict.mmblogit() with `conditional = TRUE` (requires `VA_dummy_new`, the
  # dataset where new group levels were replaced by existing ones):
  if (packageVersion("mclogit") <= "0.8.7.3" ||
      packageVersion("mclogit") >= "0.9.4.1") {
    # Not run for mclogit versions between 0.8.7.3 and 0.9.4.1 because of
    # mclogit's GitHub issue #23.
    warn_expected <- if (packageVersion("mclogit") <= "0.8.7.3") {
      "variable 'prior' is absent, its contrast will be ignored"
    } else {
      NA
    }
    expect_warning(
      lpreds_orig_dummy <- predict(mfit, newdata = VA_new_dummy),
      warn_expected
    )
  }

  ranslopes_orig_dummy <- list(
    treat2 = rbind(
      mranef$grp[
        seq_len(nlats) + nlats + (lvl_idx_new1 - 1L) * 2L * nlats,
        1
      ],
      mranef$grp[
        seq_len(nlats) + nlats + (lvl_idx1 - 1L) * 2L * nlats,
        1
      ],
      mranef$grp[
        seq_len(nlats) + nlats + (lvl_idx1 - 1L) * 2L * nlats,
        1
      ]
    ),
    age = matrix(0, nrow = nrow(VA_new_ch), ncol = nlats),
    Karn = matrix(0, nrow = nrow(VA_new_ch), ncol = nlats),
    prior10 = matrix(0, nrow = nrow(VA_new_ch), ncol = nlats)
  )
  lpreds_orig_dummy_man <- matrix(mfixef[, "(Intercept)"],
                                  nrow = nrow(VA_new_ch),
                                  ncol = nlats, byrow = TRUE) +
    rbind(
      mranef$grp[
        seq_len(nlats) + (lvl_idx_new1 - 1L) * 2L * nlats,
        1
      ],
      mranef$grp[
        seq_len(nlats) + (lvl_idx1 - 1L) * 2L * nlats,
        1
      ],
      mranef$grp[
        seq_len(nlats) + (lvl_idx1 - 1L) * 2L * nlats,
        1
      ]
    ) +
    do.call(rbind, lapply(seq_len(nrow(VA_new_ch)), function(i) {
      as.matrix(VA_new_ch[i, fixnms_b, drop = FALSE]) %*%
        (t(mfixef[, fixnms_b]) +
           do.call(rbind, lapply(ranslopes_orig_dummy[fixnms_b],
                                 function(x) x[i, ])))
    }))
  if (packageVersion("mclogit") <= "0.8.7.3" ||
      packageVersion("mclogit") >= "0.9.4.1") {
    # Not run for mclogit versions between 0.8.7.3 and 0.9.4.1 because of
    # mclogit's GitHub issue #23.
    expect_equal(unname(lpreds_orig_dummy), unname(lpreds_orig_dummy_man),
                 tolerance = 1e-14)
  }

  # Corrected version of predict.mmblogit() with `conditional = TRUE`:
  mfit_dummy <- mfit
  # Set those random effects which correspond to the "dummy" level (used as
  # replacement for the new group level) to zero:
  mfit_dummy$random.effects[names(mfit_dummy$groups) == "grp"][[1]][
    seq_len(lvl_idx1 * 2L * nlats),
    1
  ] <- 0
  if (packageVersion("mclogit") <= "0.8.7.3" ||
      packageVersion("mclogit") >= "0.9.4.1") {
    # Not run for mclogit versions between 0.8.7.3 and 0.9.4.1 because of
    # mclogit's GitHub issue #23.
    warn_expected <- if (packageVersion("mclogit") <= "0.8.7.3") {
      "variable 'prior' is absent, its contrast will be ignored"
    } else {
      NA
    }
    expect_warning(
      lpreds_orig <- predict(mfit_dummy, newdata = VA_new_dummy),
      warn_expected
    )
  }

  ranslopes_orig <- list(
    treat2 = rbind(
      mranef$grp[
        seq_len(nlats) + nlats + (lvl_idx_new1 - 1L) * 2L * nlats,
        1
      ],
      matrix(0, nrow = nrow(VA_new_ch) - 1L, ncol = nlats)
    ),
    age = matrix(0, nrow = nrow(VA_new_ch), ncol = nlats),
    Karn = matrix(0, nrow = nrow(VA_new_ch), ncol = nlats),
    prior10 = matrix(0, nrow = nrow(VA_new_ch), ncol = nlats)
  )
  lpreds_orig_man <- matrix(mfixef[, "(Intercept)"],
                            nrow = nrow(VA_new_ch),
                            ncol = nlats, byrow = TRUE) +
    rbind(
      mranef$grp[
        seq_len(nlats) + (lvl_idx_new1 - 1L) * 2L * nlats,
        1
      ],
      matrix(0, nrow = nrow(VA_new_ch) - 1L, ncol = nlats)
    ) +
    do.call(rbind, lapply(seq_len(nrow(VA_new_ch)), function(i) {
      as.matrix(VA_new_ch[i, fixnms_b, drop = FALSE]) %*%
        (t(mfixef[, fixnms_b]) +
           do.call(rbind, lapply(ranslopes_orig[fixnms_b],
                                 function(x) x[i, ])))
    }))
  if (packageVersion("mclogit") <= "0.8.7.3" ||
      packageVersion("mclogit") >= "0.9.4.1") {
    # Not run for mclogit versions between 0.8.7.3 and 0.9.4.1 because of
    # mclogit's GitHub issue #23.
    expect_equal(unname(lpreds_orig), unname(lpreds_orig_man),
                 tolerance = 1e-14)
  }

  # With `conditional = FALSE` (this is the way used by projpred):
  lpreds_orig2 <- predict(mfit_dummy, newdata = VA_new, conditional = FALSE)

  lpreds_orig2_man <- matrix(mfixef[, "(Intercept)"],
                             nrow = nrow(VA_new_ch),
                             ncol = nlats, byrow = TRUE) +
    as.matrix(VA_new_ch[, fixnms_b, drop = FALSE]) %*%
    t(mfixef[, fixnms_b])
  expect_equal(unname(lpreds_orig2), unname(lpreds_orig2_man),
               tolerance = 1e-14)

  #### With repair_re() -----------------------------------------------------

  set.seed(seed3_tst)
  lpreds <- subprd_augdat(list(mfit), newdata = VA_new)
  expect_identical(dim(lpreds), c(nrow(VA_new), nlats, 1L))
  lpreds <- lpreds[, , 1]

  set.seed(seed3_tst)
  ranmulti_grp <- mvtnorm::rmvnorm(
    n = 1, sigma = mVarCorr$grp, checkSymmetry = FALSE
  )
  ranslopes <- list(
    treat2 = rbind(
      mranef$grp[
        seq_len(nlats) + nlats + (lvl_idx_new1 - 1L) * 2L * nlats,
        1
      ],
      matrix(ranmulti_grp[1, grep("~treat2$", colnames(mVarCorr$grp))],
             nrow = nrow(VA_new) - 1L, ncol = nlats, byrow = TRUE)
    ),
    age = matrix(0, nrow = nrow(VA_new), ncol = nlats),
    Karn = matrix(0, nrow = nrow(VA_new), ncol = nlats),
    prior10 = matrix(0, nrow = nrow(VA_new), ncol = nlats)
  )
  lpreds_man <- matrix(mfixef[, "(Intercept)"],
                       nrow = nrow(VA_new_ch),
                       ncol = nlats, byrow = TRUE) +
    rbind(
      mranef$grp[
        seq_len(nlats) + (lvl_idx_new1 - 1L) * 2L * nlats,
        1
      ],
      matrix(ranmulti_grp[1, grep("~1$", colnames(mVarCorr$grp))],
             nrow = nrow(VA_new) - 1L, ncol = nlats, byrow = TRUE)
    ) +
    do.call(rbind, lapply(seq_len(nrow(VA_new_ch)), function(i) {
      as.matrix(VA_new_ch[i, fixnms_b, drop = FALSE]) %*%
        (t(mfixef[, fixnms_b]) +
           do.call(rbind, lapply(ranslopes[fixnms_b],
                                 function(x) x[i, ])))
    }))
  dmnms <- list(NULL, dimnames(lpreds_man)[[2]])
  expect_equal(lpreds, structure(lpreds_man, dimnames = dmnms),
               tolerance = 1e-14)
})

## Multiple group variables -----------------------------------------------

test_that(paste(
  "repair_re() works for multilevel categorical() models with multiple group",
  "variables"
), {
  skip_if_not(packageVersion("mclogit") >= "0.9")

  ### Fit with mclogit::mblogit() -------------------------------------------

  warn_expected <- if (packageVersion("mclogit") <= "0.8.7.3") {
    "variable 'prior' is absent, its contrast will be ignored"
  } else {
    NA
  }
  expect_warning(
    out_capt <- capture.output({
      mfit <- mclogit::mblogit(
        formula = cell ~ treat + age + Karn + prior,
        data = VA,
        random = list(~ 1 | sbj,
                      ~ 1 + treat | grp),
        model = FALSE,
        y = FALSE
      )
    }),
    warn_expected
  )
  expect_identical(tail(out_capt, 1), "converged")

  ### Extract estimated parameters ------------------------------------------

  mfixef <- mfit$coefmat
  mranef <- setNames(mfit$random.effects, names(mfit$groups))
  mVarCorr <- mfit$VarCov
  if (packageVersion("mclogit") < "0.9") {
    mVarCorr <- setNames(mVarCorr, names(mfit$groups))
  }

  ### Predict ---------------------------------------------------------------

  # predict.mmblogit() with `conditional = TRUE` (requires `VA_dummy_new`, the
  # dataset where new group levels were replaced by existing ones):
  if (packageVersion("mclogit") <= "0.8.7.3" ||
      packageVersion("mclogit") >= "0.9.4.1") {
    # Not run for mclogit versions between 0.8.7.3 and 0.9.4.1 because of
    # mclogit's GitHub issue #23.
    warn_expected <- if (packageVersion("mclogit") <= "0.8.7.3") {
      "variable 'prior' is absent, its contrast will be ignored"
    } else {
      NA
    }
    expect_warning(
      lpreds_orig_dummy <- predict(mfit, newdata = VA_new_dummy),
      warn_expected
    )
  }

  ranslopes_orig_dummy <- list(
    treat2 = rbind(
      mranef$grp[
        seq_len(nlats) + nlats + (lvl_idx_new1 - 1L) * 2L * nlats,
        1
      ],
      mranef$grp[
        seq_len(nlats) + nlats + (lvl_idx1 - 1L) * 2L * nlats,
        1
      ],
      mranef$grp[
        seq_len(nlats) + nlats + (lvl_idx1 - 1L) * 2L * nlats,
        1
      ]
    ),
    age = matrix(0, nrow = nrow(VA_new_ch), ncol = nlats),
    Karn = matrix(0, nrow = nrow(VA_new_ch), ncol = nlats),
    prior10 = matrix(0, nrow = nrow(VA_new_ch), ncol = nlats)
  )
  lpreds_orig_dummy_man <- matrix(mfixef[, "(Intercept)"],
                                  nrow = nrow(VA_new_ch),
                                  ncol = nlats, byrow = TRUE) +
    rbind(
      mranef$sbj[
        seq_len(nlats) + (lvl_idx_new1_sbj - 1L) * nlats,
        1
      ],
      mranef$sbj[
        seq_len(nlats) + (lvl_idx1 - 1L) * nlats,
        1
      ],
      mranef$sbj[
        seq_len(nlats) + (lvl_idx1 - 1L) * nlats,
        1
      ]
    ) +
    rbind(
      mranef$grp[
        seq_len(nlats) + (lvl_idx_new1 - 1L) * 2L * nlats,
        1
      ],
      mranef$grp[
        seq_len(nlats) + (lvl_idx1 - 1L) * 2L * nlats,
        1
      ],
      mranef$grp[
        seq_len(nlats) + (lvl_idx1 - 1L) * 2L * nlats,
        1
      ]
    ) +
    do.call(rbind, lapply(seq_len(nrow(VA_new_ch)), function(i) {
      as.matrix(VA_new_ch[i, fixnms_b, drop = FALSE]) %*%
        (t(mfixef[, fixnms_b]) +
           do.call(rbind, lapply(ranslopes_orig_dummy[fixnms_b],
                                 function(x) x[i, ])))
    }))
  if ((packageVersion("mclogit") <= "0.8.7.3" ||
       packageVersion("mclogit") >= "0.9.4.1") &&
      packageVersion("mclogit") >= "0.9.4.2") {
    # Not run for mclogit versions between 0.8.7.3 and 0.9.4.1 because of
    # mclogit's GitHub issue #23.
    # Not run for mclogit versions < 0.9.4.2 because of mclogit's GitHub issue
    # #24.
    expect_equal(unname(lpreds_orig_dummy), unname(lpreds_orig_dummy_man),
                 tolerance = 1e-14)
  }

  # Corrected version of predict.mmblogit() with `conditional = TRUE`:
  mfit_dummy <- mfit
  # Set those random effects which correspond to the "dummy" level (used as
  # replacement for the new group level) to zero:
  mfit_dummy$random.effects[names(mfit_dummy$groups) == "sbj"][[1]][
    seq_len(lvl_idx1 * nlats),
    1
  ] <- 0
  mfit_dummy$random.effects[names(mfit_dummy$groups) == "grp"][[1]][
    seq_len(lvl_idx1 * 2L * nlats),
    1
  ] <- 0
  if (packageVersion("mclogit") <= "0.8.7.3" ||
      packageVersion("mclogit") >= "0.9.4.1") {
    # Not run for mclogit versions between 0.8.7.3 and 0.9.4.1 because of
    # mclogit's GitHub issue #23.
    warn_expected <- if (packageVersion("mclogit") <= "0.8.7.3") {
      "variable 'prior' is absent, its contrast will be ignored"
    } else {
      NA
    }
    expect_warning(
      lpreds_orig <- predict(mfit_dummy, newdata = VA_new_dummy),
      warn_expected
    )
  }

  ranslopes_orig <- list(
    treat2 = rbind(
      mranef$grp[
        seq_len(nlats) + nlats + (lvl_idx_new1 - 1L) * 2L * nlats,
        1
      ],
      matrix(0, nrow = nrow(VA_new_ch) - 1L, ncol = nlats)
    ),
    age = matrix(0, nrow = nrow(VA_new_ch), ncol = nlats),
    Karn = matrix(0, nrow = nrow(VA_new_ch), ncol = nlats),
    prior10 = matrix(0, nrow = nrow(VA_new_ch), ncol = nlats)
  )
  lpreds_orig_man <- matrix(mfixef[, "(Intercept)"],
                            nrow = nrow(VA_new_ch),
                            ncol = nlats, byrow = TRUE) +
    rbind(
      mranef$sbj[
        seq_len(nlats) + (lvl_idx_new1_sbj - 1L) * nlats,
        1
      ],
      matrix(0, nrow = nrow(VA_new_ch) - 1L, ncol = nlats)
    ) +
    rbind(
      mranef$grp[
        seq_len(nlats) + (lvl_idx_new1 - 1L) * 2L * nlats,
        1
      ],
      matrix(0, nrow = nrow(VA_new_ch) - 1L, ncol = nlats)
    ) +
    do.call(rbind, lapply(seq_len(nrow(VA_new_ch)), function(i) {
      as.matrix(VA_new_ch[i, fixnms_b, drop = FALSE]) %*%
        (t(mfixef[, fixnms_b]) +
           do.call(rbind, lapply(ranslopes_orig[fixnms_b],
                                 function(x) x[i, ])))
    }))
  if ((packageVersion("mclogit") <= "0.8.7.3" ||
       packageVersion("mclogit") >= "0.9.4.1") &&
      packageVersion("mclogit") >= "0.9.4.2") {
    # Not run for mclogit versions between 0.8.7.3 and 0.9.4.1 because of
    # mclogit's GitHub issue #23.
    # Not run for mclogit versions < 0.9.4.2 because of mclogit's GitHub issue
    # #24.
    expect_equal(unname(lpreds_orig), unname(lpreds_orig_man),
                 tolerance = 1e-14)
  }

  # With `conditional = FALSE` (this is the way used by projpred):
  lpreds_orig2 <- predict(mfit_dummy, newdata = VA_new, conditional = FALSE)

  lpreds_orig2_man <- matrix(mfixef[, "(Intercept)"],
                             nrow = nrow(VA_new_ch),
                             ncol = nlats, byrow = TRUE) +
    as.matrix(VA_new_ch[, fixnms_b, drop = FALSE]) %*%
    t(mfixef[, fixnms_b])
  expect_equal(unname(lpreds_orig2), unname(lpreds_orig2_man),
               tolerance = 1e-14)

  #### With repair_re() -----------------------------------------------------

  set.seed(seed3_tst)
  lpreds <- subprd_augdat(list(mfit), newdata = VA_new)
  expect_identical(dim(lpreds), c(nrow(VA_new), nlats, 1L))
  lpreds <- lpreds[, , 1]

  set.seed(seed3_tst)
  ranmulti_sbj <- mvtnorm::rmvnorm(
    n = 1, sigma = mVarCorr$sbj, checkSymmetry = FALSE
  )
  ranmulti_grp <- mvtnorm::rmvnorm(
    n = 1, sigma = mVarCorr$grp, checkSymmetry = FALSE
  )
  ranslopes <- list(
    treat2 = rbind(
      mranef$grp[
        seq_len(nlats) + nlats + (lvl_idx_new1 - 1L) * 2L * nlats,
        1
      ],
      matrix(ranmulti_grp[1, grep("~treat2$", colnames(mVarCorr$grp))],
             nrow = nrow(VA_new) - 1L, ncol = nlats, byrow = TRUE)
    ),
    age = matrix(0, nrow = nrow(VA_new), ncol = nlats),
    Karn = matrix(0, nrow = nrow(VA_new), ncol = nlats),
    prior10 = matrix(0, nrow = nrow(VA_new), ncol = nlats)
  )
  lpreds_man <- matrix(mfixef[, "(Intercept)"],
                       nrow = nrow(VA_new_ch),
                       ncol = nlats, byrow = TRUE) +
    rbind(
      mranef$sbj[
        seq_len(nlats) + (lvl_idx_new1_sbj - 1L) * nlats,
        1
      ],
      matrix(ranmulti_sbj[1, grep("~1$", colnames(mVarCorr$sbj))],
             nrow = nrow(VA_new) - 1L, ncol = nlats, byrow = TRUE)
    ) +
    rbind(
      mranef$grp[
        seq_len(nlats) + (lvl_idx_new1 - 1L) * 2L * nlats,
        1
      ],
      matrix(ranmulti_grp[1, grep("~1$", colnames(mVarCorr$grp))],
             nrow = nrow(VA_new) - 1L, ncol = nlats, byrow = TRUE)
    ) +
    do.call(rbind, lapply(seq_len(nrow(VA_new_ch)), function(i) {
      as.matrix(VA_new_ch[i, fixnms_b, drop = FALSE]) %*%
        (t(mfixef[, fixnms_b]) +
           do.call(rbind, lapply(ranslopes[fixnms_b],
                                 function(x) x[i, ])))
    }))
  dmnms <- list(NULL, dimnames(lpreds_man)[[2]])
  expect_equal(lpreds, structure(lpreds_man, dimnames = dmnms),
               tolerance = 1e-14)
})

# Teardown ----------------------------------------------------------------

if (exists("rng_old")) assign(".Random.seed", rng_old, envir = .GlobalEnv)

# Clean up the workspace:
rm(list = setdiff(ls(), ls_bu))
