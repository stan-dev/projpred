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
