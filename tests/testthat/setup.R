#_______________________________________________________________________________
# Setup for the unit tests
#_______________________________________________________________________________

# General setup -----------------------------------------------------------

options(warn = 1)

# These switches may be set to `FALSE` to save time (e.g., when debugging
# interactively):
# Run more tests, at the downside of increased runtime?:
run_more <- FALSE
# Run project()?:
run_prj <- identical(Sys.getenv("NOT_CRAN"), "true")
# Run varsel()?:
run_vs <- identical(Sys.getenv("NOT_CRAN"), "true")
# Run cv_varsel()?:
run_cvvs <- run_vs
# Run cv_varsel() with `validate_search = TRUE` always (`TRUE`) or just for L1
# search (`FALSE`)?:
run_valsearch_always <- FALSE
# Run cv_varsel() with `validate_search = TRUE` also for the latent and the
# augmented-data projection (and the corresponding traditional projection
# settings which are used for comparison)? (Only relevant if
# `run_valsearch_always = FALSE`.):
run_valsearch_aug_lat <- FALSE
# Run the `cvfits` test for all possible test setups (`TRUE`) or just for the
# first one among the GLMMs (`FALSE`; note that if there is no GLMM available in
# that test, the first test setup among those for K-fold CV is used)?:
run_cvfits_all <- run_more
# Run tests for "brmsfit"s?:
run_brms <- identical(Sys.getenv("NOT_CRAN"), "true")
# Run snapshot tests?:
# Notes:
#   * For general information about snapshot tests, see, e.g.,
#   `?expect_snapshot` and `vignette("snapshotting", package = "testthat")`.
#   * The snapshot tests are at least OS-dependent (perhaps even
#   machine-dependent), so they only make sense locally. Therefore, we don't run
#   the snapshot tests on CRAN or continuous integration (CI) systems. The
#   detection of a CI system by the help of environment variable `CI` needs
#   special care, see <https://github.com/r-lib/testthat/issues/825> and the
#   source code of `testthat:::on_ci()`.
#   * The last of the following conditions avoids that the snapshot tests are
#   run by a local `R CMD check` (at least in RStudio). The reason for avoiding
#   this is that in `R CMD check`, the previous snapshots are not available (at
#   least as long as they are listed in the `.Rbuildignore` file), so they would
#   be re-created, which would throw a lot of test warnings (which could obscure
#   potentially important warnings).
run_snaps <- Sys.getenv("RUN_SNAPS")
if (identical(run_snaps, "")) {
  run_snaps <- identical(Sys.getenv("NOT_CRAN"), "true") &&
    !identical(toupper(Sys.getenv("CI")), "TRUE") &&
    identical(Sys.getenv("_R_CHECK_FORCE_SUGGESTS_"), "")
} else {
  run_snaps <- as.logical(toupper(run_snaps))
  stopifnot(isTRUE(run_snaps) || isFALSE(run_snaps))
}
if (run_snaps && !requireNamespace("rlang", quietly = TRUE)) {
  warning("Package 'rlang' is needed for snapshot testing, but could not be ",
          "found. Deactivating snapshot testing now.")
  run_snaps <- FALSE
}
if (run_snaps) {
  testthat_ed_max2 <- edition_get() <= 2
}
# Run tests for the parallelization of the projection?:
# Notes:
#   * Throughout the tests, the terms "parallelization" and "parallel" refer to
#   the parallelization of the projection ("projection parallelization"), not
#   the parallelization of the CV ("CV parallelization").
#   * We don't run the parallel tests on CRAN or continuous integration (CI)
#   systems because parallelization might require special care there.
#   * Currently, parallelization on Windows takes longer than running
#   sequentially. This makes parallelization impractical on Windows, so we
#   don't run the tests on Windows by default.
#   * Currently, parallelization only works reliably for GLMs (because of
#   memory issues for more complex models like GLMMs, GAMs and GAMMs).
#   Therefore, we will only test GLMs here.
run_prll <- identical(Sys.getenv("NOT_CRAN"), "true") &&
  !identical(toupper(Sys.getenv("CI")), "TRUE") &&
  !identical(.Platform$OS.type, "windows")
if (run_prll) {
  ncores <- parallel::detectCores(logical = FALSE)
  if (ncores == 1) {
    warning("Deactivating the parallel tests because only a single worker ",
            "could be detected.")
    run_prll <- FALSE
  }
  # Do not run on more than 2 cores if requested so:
  if (identical(Sys.getenv("_R_CHECK_LIMIT_CORES_"), "TRUE")) {
    ncores <- min(ncores, 2L)
  }
  # Use the 'doParallel' package on all platforms except Windows. For Windows,
  # the 'doFuture' package provides a faster alternative via the 'future.callr'
  # package (which is still slower than a sequential run, though):
  if (!identical(.Platform$OS.type, "windows")) {
    if (!requireNamespace("doParallel", quietly = TRUE)) {
      warning("Package 'doParallel' is needed for the parallel tests, but ",
              "could not be found. Deactivating the parallel tests now.")
      run_prll <- FALSE
    } else {
      dopar_backend <- "doParallel"
    }
  } else {
    # This case (which should not be possible by default) is only included
    # here to demonstrate how parallelization should be used on Windows (but
    # currently, this makes no sense, as explained above).
    if (!requireNamespace("doFuture", quietly = TRUE)) {
      warning("Package 'doFuture' is needed for the parallel tests, but ",
              "could not be found. Deactivating the parallel tests now.")
      run_prll <- FALSE
    } else {
      dopar_backend <- "doFuture"
    }
    if (!identical(.Platform$OS.type, "windows")) {
      # This case (which should not be possible by default) is only included
      # here to demonstrate how other systems should be used with the 'doFuture'
      # package.
      future_plan <- "multicore"
    } else {
      ### Not used in this case because the 'future.callr' package provides a
      ### faster alternative on Windows (which is still slower than a sequential
      ### run, though):
      # future_plan <- "multisession"
      ###
      if (!requireNamespace("future.callr", quietly = TRUE)) {
        warning("Package 'future.callr' is needed for the parallel tests, but ",
                "could not be found. Deactivating the parallel tests now.")
        run_prll <- FALSE
      } else {
        future_plan <- "callr"
      }
    }
  }
}
# Run all test scripts (following this setup script) in a completely random RNG
# state? (The tests should still pass then, because in all situations where RNG
# is used, a specific seed is supposed to be set.):
run_randRNG <- identical(Sys.getenv("NOT_CRAN"), "true")
# Run tests for additive models (GAMs and GAMMs)?:
run_additive <- TRUE

# Use a factor or an integer response for ordinal and categorical families?:
use_fac <- TRUE

# Use polym() instead of poly()?:
use_polym <- FALSE

source(testthat::test_path("helpers", "unlist_cust.R"), local = TRUE)
source(testthat::test_path("helpers", "testers.R"), local = TRUE)
source(testthat::test_path("helpers", "args.R"), local = TRUE)
source(testthat::test_path("helpers", "getters.R"), local = TRUE)
source(testthat::test_path("helpers", "formul_handlers.R"), local = TRUE)
source(testthat::test_path("helpers", "predictor_handlers.R"), local = TRUE)
source(testthat::test_path("helpers", "dummies.R"), local = TRUE)

# Note: The following `mod_nms` refer to *generalized* (linear/additive,
# multilevel) models. This is due to history (when these tests were written,
# only such *generalized* models were supported by projpred). Now that more
# models are supported (even non-generalized ones), these model names are not
# really correct anymore. However, we keep them for simplicity.
mod_nms <- c("glm", "glmm", "gam", "gamm")
if (run_additive) {
  # Suppress the warning for additive models (GAMs and GAMMs) stating that their
  # implementation is currently only experimental:
  options(projpred.warn_additive_experimental = FALSE)
} else {
  mod_nms <- setdiff(mod_nms, c("gam", "gamm"))
}
mod_nms <- setNames(nm = mod_nms)

fam_nms_trad <- c("gauss", "brnll", "binom", "poiss")
fam_nms_ordin <- c("cumul", "srtio", "crtio", "adcat")
fam_nms_categ <- "categ"
fam_nms_aug <- c(fam_nms_ordin, fam_nms_categ)
fam_nms <- c(fam_nms_trad, fam_nms_aug)
fam_nms_unsupp <- setdiff(fam_nms_ordin, "cumul")
fam_nms_brms_only <- setdiff(fam_nms_aug, "cumul")
if (!run_brms) {
  fam_nms <- setdiff(fam_nms, fam_nms_brms_only)
}
fam_nms_trad <- setNames(nm = fam_nms_trad)
fam_nms_ordin <- setNames(nm = fam_nms_ordin)
fam_nms_categ <- setNames(nm = fam_nms_categ)
fam_nms_aug <- setNames(nm = fam_nms_aug)
fam_nms <- setNames(nm = fam_nms)
fam_nms_unsupp <- setNames(nm = fam_nms_unsupp)
fam_nms_brms_only <- setNames(nm = fam_nms_brms_only)
# Long names:
fam_nms_aug_long <- c(sapply(fam_nms_aug, get_fam_long),
                      cumul = "cumulative_rstanarm")
fam_nms_ordin_long <- c(sapply(fam_nms_ordin, get_fam_long),
                        cumul = "cumulative_rstanarm")
fam_nms_long <- c(sapply(fam_nms, get_fam_long_full),
                  cumul = "cumulative_rstanarm")
# Regular expressions:
fam_nms_aug_regex <- paste0("\\.(", paste(fam_nms_aug, collapse = "|"), ")\\.")
fam_nms_unsupp_regex <- paste0("\\.(", paste(fam_nms_unsupp, collapse = "|"),
                               ")\\.")

# Needed for package mclogit (providing the submodel fitter for multilevel
# brms::categorical() models):
warn_mclogit <- if (packageVersion("mclogit") >= "0.9.6") {
  "Inner iterations did not coverge"
} else {
  paste0("^step size truncated due to possible divergence$|",
         "^Algorithm stopped due to false convergence$")
}

# Data --------------------------------------------------------------------

## Setup ------------------------------------------------------------------

# Number of observations:
nobsv <- 41L
# Values for testing:
nobsv_tst <- c(1L, nobsv %/% 3L)

# For ordinal models (but also used for categorical models):
nthres <- 2L
ncat <- nthres + 1L
yunq_num <- seq_len(ncat)
yunq_chr <- paste0("y", yunq_num)
# The intercepts at centered predictors, also known as thresholds:
thres <- qlogis(seq_len(nthres) / ncat)
link_str <- "logit"

# Seed:
seed_dat <- 8541351
set.seed(seed_dat)

## GLMs --------------------------------------------------------------------
## Add population-level effects to the intercept-(and-offset-)only model

nterms_cont <- 3L
x_cont <- matrix(rnorm(nobsv * nterms_cont), nobsv, nterms_cont)
b_cont <- runif(nterms_cont, min = -0.5, max = 0.5)

nlvl_fix <- c(3L, 2L)
nlvl_fix <- setNames(nlvl_fix, seq_along(nlvl_fix))
if (length(nlvl_fix) <= 1) {
  names(nlvl_fix) <- paste0("z.", names(nlvl_fix))
}
nterms_cate <- length(nlvl_fix)
x_cate_list <- lapply(nlvl_fix, function(nlvl_fix_i) {
  x_cate <- gl(n = nlvl_fix_i, k = floor(nobsv / nlvl_fix_i), length = nobsv,
               labels = paste0("lvl", seq_len(nlvl_fix_i)))
  b_cate <- runif(nlvl_fix_i, min = -0.5, max = 0.5)
  ### Using a model.matrix() approach:
  # x_cate_mat <- model.matrix(~ 0 + x_cate)
  # eta_cate <- x_cate_mat %*% b_cate
  ###
  ### Using an indexing approach:
  eta_cate <- b_cate[x_cate]
  ###
  return(nlist(x_cate, eta_cate, b_cate))
})

# Intercept and offsets (offsets are added later):
icpt <- -0.42
offs_tst <- rnorm(nobsv)

eta_glm <- icpt +
  x_cont %*% b_cont +
  do.call("+", lapply(x_cate_list, "[[", "eta_cate"))

nterms_glm <- nterms_cont + nterms_cate

## GLMMs ------------------------------------------------------------------
## Add group-level effects to the GLMs, yielding GLMMs

nlvl_ran <- c(6L)
nlvl_ran <- setNames(nlvl_ran, seq_along(nlvl_ran))
if (length(nlvl_ran) <= 1) {
  names(nlvl_ran) <- paste0("z.", names(nlvl_ran))
}
# Multiply by 2 because of the random intercept as well as the random slope (the
# latter for `x_cont[, 1]`):
nterms_z <- length(nlvl_ran) * 2L
z_list <- lapply(nlvl_ran, function(nlvl_ran_i) {
  z <- gl(n = nlvl_ran_i, k = floor(nobsv / nlvl_ran_i), length = nobsv,
          labels = paste0("lvl", seq_len(nlvl_ran_i)))
  r_icpts <- rnorm(nlvl_ran_i, sd = 0.4)
  r_xco1 <- rnorm(nlvl_ran_i, sd = 0.4)
  eta_z <- r_icpts[z] + r_xco1[z] * x_cont[, 1]
  return(nlist(z, eta_z, r_icpts, r_xco1))
})
eta_glmm <- eta_glm + do.call("+", lapply(z_list, "[[", "eta_z"))

nterms_glmm <- nterms_glm + nterms_z

## GAMs -------------------------------------------------------------------
## Add nonlinear (smoothed) effects to the GLMs

# For simplicity, always use the same nonlinear function (could be extended in
# the future):
s_mat <- apply(x_cont[, 1, drop = FALSE], 2, function(x, b = 2, c = - pi / 4) {
  b * sin(x - c)
})
s_sum <- rowSums(s_mat)
nterms_s <- ncol(s_mat)
eta_gam <- eta_glm + s_sum

# Multiply by 2 because of the baseline linear term as well as the standard
# deviation for the wiggliness around it:
nterms_gam <- nterms_glm + 2L * nterms_s

## GAMMs ------------------------------------------------------------------
## Add nonlinear (smoothed) effects to the GLMMs

eta_gamm <- eta_glmm + s_sum

nterms_gamm <- nterms_glmm + 2L * nterms_s

## Combined dataset -------------------------------------------------------

f_gauss <- gaussian()
f_binom <- f_brnll <- binomial()
f_poiss <- poisson()
dis_tst <- runif(1L, 1, 2)
wobs_tst <- sample(1:4, nobsv, replace = TRUE)
# For the "brnll" `fam_nm`, offsets are simply not added to have some
# scenarios without offsets.
# For GAMs, offsets are not added because of rstanarm issue #546 (see
# also further below).
# For GAMMs, offsets are not added because of rstanarm issue #253 (see
# also further below).
# (The brms "gam" and "gamm" cases are handled in the same way as the rstanarm
# "gam" and "gamm" cases to avoid too many special cases.)
# For the "categ" `fam_nm`, offsets are not added because they are currently not
# supported for it.
offs_expr <- expression(!(fam_nm %in% c("brnll", "categ") ||
                            mod_nm %in% c("gam", "gamm")))
cre_dat <- function(idxs_crr, offs_crr, wobs_crr, dis_crr) {
  nobsv_crr <- length(idxs_crr)
  dat_crr <- lapply(mod_nms, function(mod_nm) {
    lapply(fam_nms, function(fam_nm) {
      pred_link <- get(paste0("eta_", mod_nm))
      pred_link <- pred_link[idxs_crr, , drop = FALSE]
      if (eval(offs_expr)) {
        pred_link <- pred_link + offs_crr
      }
      if (fam_nm %in% fam_nms_ordin) {
        pred_link <- pred_link - icpt
        if (eval(offs_expr)) {
          # The equal-probability thresholds defined above refer to the state
          # before offsets are added, so we need to subtract them here in the
          # data-generating model:
          pred_link <- pred_link - offs_crr
        }
        thres_eta <- sapply(thres, function(thres_k) {
          thres_k - pred_link
        })
      } else if (fam_nm %in% fam_nms_categ) {
        pred_link <- sapply(thres, function(thres_k) {
          thres_k + pred_link
        })
      }
      if (fam_nm == "cumul") {
        pred_resp <- augdat_ilink_cumul(thres_eta, link = link_str)
      } else if (fam_nm %in% fam_nms_ordin) {
        if (fam_nm %in% c("crtio", "adcat")) {
          thres_eta <- -thres_eta
        }
        ilink_crr <- get(paste0("inv_link_", get_fam_long(fam_nm)),
                         asNamespace("brms"), mode = "function",
                         inherits = FALSE)
        pred_resp <- ilink_crr(thres_eta, link = link_str)
      } else if (fam_nm %in% fam_nms_categ) {
        ilink_crr <- get(paste0("inv_link_", get_fam_long(fam_nm)),
                         asNamespace("brms"), mode = "function",
                         inherits = FALSE)
        pred_resp <- ilink_crr(pred_link)
      } else {
        pred_resp <- get(paste0("f_", fam_nm))$linkinv(pred_link)
      }
      if (fam_nm == "gauss") {
        return(rnorm(nobsv_crr, mean = pred_resp, sd = dis_crr))
      } else if (fam_nm == "brnll") {
        return(rbinom(nobsv_crr, 1, pred_resp))
      } else if (fam_nm == "binom") {
        return(rbinom(nobsv_crr, wobs_crr, pred_resp))
      } else if (fam_nm == "poiss") {
        return(rpois(nobsv_crr, pred_resp))
      } else if (fam_nm %in% fam_nms_aug) {
        ryunq <- sapply(seq_len(nobsv_crr), function(i_obs) {
          sample(yunq_num, size = 1L, prob = pred_resp[i_obs, ])
        })
        if (use_fac) {
          ryunq <- factor(ryunq, levels = yunq_num, labels = yunq_chr,
                          ordered = fam_nm %in% fam_nms_ordin)
        }
        return(ryunq)
      } else {
        stop("Unknown `fam_nm`.")
      }
    })
  })
  dat_crr <- unlist(dat_crr, recursive = FALSE)
  names(dat_crr) <- paste("y", gsub("\\.", "_", names(dat_crr)), sep = "_")
  return(dat_crr)
}
dat <- cre_dat(idxs_crr = seq_len(nobsv), offs_crr = offs_tst,
               wobs_crr = wobs_tst, dis_crr = dis_tst)
dat <- data.frame(
  dat,
  xco = x_cont, xca = lapply(x_cate_list, "[[", "x_cate"),
  z = lapply(z_list, "[[", "z"),
  s = s_mat,
  wobs_col = wobs_tst, offs_col = offs_tst,
  check.names = FALSE
)
if (ncol(s_mat) == 1) {
  names(dat)[names(dat) == "s"] <- "s.1"
}

## Modified datasets ------------------------------------------------------

dat_wobs_ones <- within(dat, {
  wobs_col <- NULL
  wobs_col_ones <- rep_len(1, length.out = nobsv)
})
dat_wobs_new <- within(dat, {
  wobs_col <- NULL
  wobs_col_new <- rep_len(2:5, length.out = nobsv)
})

dat_offs_zeros <- within(dat, {
  offs_col <- NULL
  offs_col_zeros <- rep_len(0, length.out = nobsv)
})
dat_offs_new <- within(dat, {
  offs_col <- NULL
  offs_col_new <- seq(-2, 2, length.out = nobsv)
})

nobsv_indep <- tail(nobsv_tst, 1)
dis_indep <- runif(1L, 1, 2)
offs_indep <- rnorm(nobsv_indep)
wobs_indep <- sample(1:4, nobsv_indep, replace = TRUE)
idxs_indep <- sample.int(nobsv, size = nobsv_indep, replace = TRUE)
dat_indep <- cre_dat(idxs_crr = idxs_indep, offs_crr = offs_indep,
                     wobs_crr = wobs_indep, dis_crr = dis_indep)
dat_indep <- cbind(
  as.data.frame(dat_indep),
  dat[idxs_indep,
      grep("^y_", names(dat), value = TRUE, invert = TRUE),
      drop = FALSE]
)
dat_indep$wobs_col <- wobs_indep
dat_indep$offs_col <- offs_indep

# Fits --------------------------------------------------------------------

## Setup ------------------------------------------------------------------

if (!requireNamespace("rstanarm", quietly = TRUE)) {
  warning("Package 'rstanarm' is needed for the rstanarm tests, but could not ",
          "be found. Deactivating the rstanarm tests now. Furthermore, in ",
          "this case, `run_prj`, `run_vs`, and `run_cvvs` are currently set ",
          "to `FALSE`.")
  pkg_nms <- character()
  # TODO: Adapt the tests to avoid the following line, at least if `run_brms` is
  # `TRUE` (better: avoid it in any case, no matter whether `run_brms` is `TRUE`
  # or `FALSE`).
  run_prj <- run_vs <- run_cvvs <- FALSE
} else {
  pkg_nms <- "rstanarm"
}

if (run_brms && !requireNamespace("brms", quietly = TRUE)) {
  warning("Package 'brms' is needed for the brms tests, but could not be ",
          "found. Deactivating the brms tests now.")
  run_brms <- FALSE
}
if (run_brms) {
  pkg_nms <- c(pkg_nms, "brms")
  # For storing "brmsfit"s locally:
  file_pth <- testthat::test_path("bfits")
  if (!dir.exists(file_pth)) dir.create(file_pth)
  # Backend:
  if (identical(Sys.getenv("TESTS_BRMS_BACKEND"), "cmdstanr") &&
      requireNamespace("cmdstanr", quietly = TRUE) &&
      # Relative file paths for cmdstanr's global option
      # `cmdstanr_write_stan_file_dir` didn't work before cmdstanr PR #665.
      # Using the workaround `file.path(getwd(), file_pth)` instead of only
      # `file_pth` also doesn't work in `R CMD check` (it doesn't throw any
      # exceptions, but recompilations take place, causing a huge increase in
      # runtime). At the time of cmdstanr's PR #665, the (development) version
      # number of cmdstanr was 0.5.2.1, so requiring >= 0.5.3 guarantees that
      # the fix is included:
      packageVersion("cmdstanr") >= "0.5.3" &&
      !is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE))) {
    options(brms.backend = "cmdstanr")
    options(cmdstanr_write_stan_file_dir = file_pth)
  }
}
pkg_nms <- setNames(nm = pkg_nms)

chains_tst <- 2L
iter_tst <- 500L
nrefdraws <- chains_tst * iter_tst %/% 2L
seed_fit <- 74345

### Formula ---------------------------------------------------------------

# Notes:
#   * Argument `offset` is not supported by rstanarm::stan_gamm4(). Instead, use
#     offset() in the formula (like for all other models). However, because of
#     rstanarm issue #546 and rstanarm issue #253, omit the offsets in GAMs and
#     GAMMs.
#   * In rstanarm::stan_gamm4(), multilevel terms are specified via argument
#     `random`.

trms_common <- c("xco.1", "xco.2", "xco.3", "xca.1", "xca.2",
                 "offset(offs_col)")
trms_grp <- c("(xco.1 | z.1)")
trms_add <- c("s(s.1)") # , "s(s.2)", "s(s.3)"
if (use_polym) {
  trm_poly <- "polym(xco.1, degree = 2, raw = TRUE)"
} else {
  trm_poly <- "poly(xco.1, 2, raw = TRUE)"
}
trms_common_spcl <- c(trm_poly,
                      "sqrt(abs(xco.3)^2) * I(!as.logical(xco.3 > 0))", "xca.1",
                      "xca.2", "offset(offs_col)")
trms_universe <- unique(c(trms_common, trms_grp, trms_add, trms_common_spcl))
trms_universe_split <- setdiff(trms_universe, "offset(offs_col)")
# Handle interaction terms:
stopifnot(!any(grepl(":", trms_universe_split)))
trms_universe_split_IA <- grep("\\*", trms_universe_split, value = TRUE)
if (length(trms_universe_split_IA)) {
  trms_universe_split_noIA <- setdiff(trms_universe_split,
                                      trms_universe_split_IA)
  # Replace " * " in interaction terms by ":":
  trms_universe_split_IA <- gsub(" \\* ", ":", trms_universe_split_IA)
  # Add main effects from interaction terms:
  trms_universe_split_noIA <- union(
    trms_universe_split_noIA,
    unlist(strsplit(trms_universe_split_IA, split = ":"))
  )
  trms_universe_split <- c(trms_universe_split_noIA, trms_universe_split_IA)
}
# Add lower-order group-level terms:
if ("(xco.1 | z.1)" %in% trms_universe_split) {
  trms_universe_split <- union(trms_universe_split, "(1 | z.1)")
}
# Ensure that for all terms on the left-hand side of the `|` character in
# group-level terms, corresponding population-level terms exist:
if ("(xco.1 | z.1)" %in% trms_universe_split) {
  trms_universe_split <- union(trms_universe_split, "xco.1")
}

# Solution terms for project()-ing from `"refmodel"`s:
solterms_x <- c("xco.1", "xca.1")
solterms_z <- c("(1 | z.1)", "(xco.1 | z.1)") # removing one of them later
solterms_s <- c("s(s.1)") # , "s(s.2)"
solterms_spcl <- c("xca.1", trm_poly,
                   "sqrt(abs(xco.3)^2)", "I(!as.logical(xco.3 > 0))",
                   "sqrt(abs(xco.3)^2):I(!as.logical(xco.3 > 0))")

### Weights (observations) ------------------------------------------------

# Argument `weights` is not needed when using the cbind() syntax (for the
# binomial family with > 1 trials). Furthermore, rstanarm:::kfold.stanreg() does
# not support weights. Thus, we have two possible options for the `weights`
# argument:
wobss_tst <- list(with_wobs = list(weights = wobs_tst),
                  without_wobs = list())

### Offsets ---------------------------------------------------------------

### See the notes above: Due to rstanarm issue #541, the fact that rstanarm
### doesn't support argument `offset` for GAMs and GAMMs, and the fact that
### brms::brm() has no argument `offset`, the easiest way to use offsets is to
### always specify them in the formula (or not at all, see the definition of
### `args_fit` below). Therefore, the following object which could be used to
### use *argument* `offset` is just a dummy (but used nevertheless, namely to
### construct the names for the argument lists):
offss_tst <- list(with_offs = list(), # with_offs = list(offset = offs_tst),
                  without_offs = list())
###

### Argument list ---------------------------------------------------------

# For some arguments, if they are specified via objects,
# rstanarm:::kfold.stanreg() seems to assume these objects to lie in the
# global environment. Since `testthat` uses a new environment for running
# the tests (see `?testthat::test_env`), we need the following code to be
# able to run devtools::test():
for (obj_symb_chr in c(paste0("f_", fam_nms_trad))) {
  if (!exists(obj_symb_chr, envir = .GlobalEnv)) {
    assign(obj_symb_chr, get(obj_symb_chr), envir = .GlobalEnv)
  }
}

args_fit <- lapply(pkg_nms, function(pkg_nm) {
  # Depending on the package used for fitting the reference model, `mod_nms`
  # could be restricted here.
  mod_nms <- setNames(nm = mod_nms)
  lapply(mod_nms, function(mod_nm) {
    if (pkg_nm == "rstanarm") {
      fam_nms <- setdiff(fam_nms, fam_nms_brms_only)
      if (mod_nm != "glm" || !use_fac) {
        # rstanarm::stan_polr() does not support multilevel or additive terms
        # and it also does not support a numeric response:
        fam_nms <- setdiff(fam_nms, "cumul")
      }

      if (mod_nm != "gamm") {
        random_arg <- list()
      } else {
        random_arg <- list(random = as.formula(paste("~", trms_grp)))
      }
    }

    if (mod_nm != "glm") {
      if (mod_nm %in% c("gam", "gamm")) {
        # Additive models are currently not supported by the augmented-data
        # projection:
        fam_nms <- setdiff(fam_nms, fam_nms_aug)
      }
      if (pkg_nm == "brms") {
        # For speed reasons, do not test all families:
        if (mod_nm == "glmm") {
          fam_nms <- intersect(fam_nms, c("brnll", "cumul", "categ"))
        } else {
          fam_nms <- intersect(fam_nms, "binom")
        }
      }
      # Because of issue #207:
      fam_nms <- setdiff(fam_nms, "poiss")
    }

    if (pkg_nm == "rstanarm" && mod_nm == "gamm") {
      # Exclude "binom" from `fam_nms` since there seems to be an issue with
      # get_refmodel.stanreg() in this case (probably issue #148):
      fam_nms <- setdiff(fam_nms, "binom")
      # TODO (GAMMs): Fix this. This exclusion also has the downside that K-fold
      # CV cannot be tested in that case.
    }

    fam_nms <- setNames(nm = fam_nms)
    lapply(fam_nms, function(fam_nm) {
      y_chr <- paste("y", mod_nm, fam_nm, sep = "_")

      if (fam_nm == "gauss" && !(pkg_nm == "rstanarm" && mod_nm == "gamm")) {
        # Here, we also test a special formula (the rstanarm "gamm" case is
        # excluded because of rstanarm issue #545):
        formul_nms <- c("stdformul", "spclformul")
      } else {
        formul_nms <- "stdformul"
      }

      fam_nm_long <- get_fam_long(fam_nm)
      if (pkg_nm == "brms" && !is.na(fam_nm_long)) {
        family_crr <- substitute(
          get(fam_nm_long_subst, envir = asNamespace("brms"))(),
          list(fam_nm_long_subst = fam_nm_long)
        )
      } else {
        family_crr <- as.name(paste0("f_", fam_nm))
      }

      if (eval(offs_expr)) {
        offss_nms <- "with_offs"
      } else {
        offss_nms <- "without_offs"
      }

      formul_nms <- setNames(nm = formul_nms)
      lapply(formul_nms, function(formul_nm) {
        if (formul_nm == "spclformul") {
          trms_common <- trms_common_spcl
          if (fam_nm != "gauss") {
            stop("`y_chr` needs to be adapted for families other than ",
                 "`\"gauss\"`.")
          }
          y_chr <- paste0("log(abs(", y_chr, "))")
        }
        if (fam_nm == "binom") {
          if (pkg_nm == "rstanarm") {
            y_chr <- paste0("cbind(", y_chr, ", wobs_col - ", y_chr, ")")
          } else if (pkg_nm == "brms") {
            y_chr <- paste(y_chr, "| trials(wobs_col)")
          }
        }
        trms <- switch(mod_nm,
                       "glm" = trms_common,
                       "glmm" = c(trms_common, trms_grp),
                       "gam" = c(trms_common, trms_add),
                       "gamm" = switch(
                         pkg_nm,
                         "rstanarm" = c(trms_common, trms_add),
                         "brms" = c(trms_common, trms_add, trms_grp),
                         stop("Unknown `pkg_nm`.")
                       ),
                       stop("Unknown `mod_nm`."))

        if (fam_nm %in% c("brnll", "binom", fam_nms_aug)) {
          # In this case, observation weights are not supported by projpred:
          wobss_nms <- "without_wobs"
        } else {
          wobss_nms <- "with_wobs"
        }

        wobss_nms <- setNames(nm = wobss_nms)
        lapply(wobss_nms, function(wobss_nm) {
          if (pkg_nm == "brms" && wobss_nm == "with_wobs") {
            if (fam_nm == "binom") {
              stop("Because of `\"| trials(wobs_col)\"` above, the code here ",
                   "(for pasting `\"| weights(wobs_col)\"`) needs to be ",
                   "adapted.")
            }
            y_chr <- paste(y_chr, "| weights(wobs_col)")
          }

          offss_nms <- setNames(nm = offss_nms)
          lapply(offss_nms, function(offss_nm) {
            if (offss_nm == "without_offs") {
              trms <- setdiff(trms, "offset(offs_col)")
            }
            formul_crr <- as.formula(paste(
              y_chr, "~", paste(trms, collapse = " + ")
            ))

            if (pkg_nm == "rstanarm") {
              pkg_args <- c(list(QR = TRUE),
                            wobss_tst[[wobss_nm]],
                            offss_tst[[offss_nm]],
                            random_arg)
            } else if (pkg_nm == "brms") {
              pkg_args <- list(file = file_pth,
                               file_refit = "on_change",
                               silent = 2)
              if (identical(getOption("brms.backend", "rstan"), "cmdstanr")) {
                pkg_args <- c(pkg_args, list(diagnostics = NULL))
              }
            }

            return(c(
              nlist(mod_nm, fam_nm, pkg_nm, formula = formul_crr,
                    family = family_crr, data = quote(dat),
                    chains = chains_tst, iter = iter_tst, seed = seed_fit,
                    refresh = 0),
              pkg_args
            ))
          })
        })
      })
    })
  })
})
args_fit <- unlist_cust(args_fit)
stopifnot(length(unique(names(args_fit))) == length(args_fit))
# For "brmsfit"s, set a unique file name (done here because during the creation
# of `args_fit`, these unique names are not easily accessible):
args_fit <- lapply(setNames(nm = names(args_fit)), function(args_fit_nm) {
  if (args_fit[[args_fit_nm]]$pkg_nm == "brms" &&
      "file" %in% names(args_fit[[args_fit_nm]])) {
    args_fit[[args_fit_nm]]$file <- file.path(
      args_fit[[args_fit_nm]]$file,
      paste0("bfit_", args_fit_nm)
    )
  }
  return(args_fit[[args_fit_nm]])
})

if (!run_more && length(pkg_nms)) {
  sel_fits <- c(
    "rstanarm.glm.gauss.stdformul.with_wobs.with_offs",
    "rstanarm.glm.brnll.stdformul.without_wobs.without_offs",
    "rstanarm.glmm.gauss.spclformul.with_wobs.with_offs",
    "rstanarm.gam.gauss.spclformul.with_wobs.without_offs",
    "rstanarm.gamm.brnll.stdformul.without_wobs.without_offs",
    "brms.glm.poiss.stdformul.with_wobs.with_offs",
    "brms.glmm.brnll.stdformul.without_wobs.without_offs",
    # "brms.gam.binom.stdformul.without_wobs.without_offs",
    "brms.gamm.binom.stdformul.without_wobs.without_offs",
    # grep(paste(paste0("\\.", fam_nms_aug, "\\."), collapse = "|"),
    #      names(args_fit), value = TRUE)
    "rstanarm.glm.cumul.stdformul.without_wobs.with_offs",
    "brms.glm.cumul.stdformul.without_wobs.with_offs",
    "brms.glm.srtio.stdformul.without_wobs.with_offs",
    "brms.glm.crtio.stdformul.without_wobs.with_offs",
    "brms.glm.adcat.stdformul.without_wobs.with_offs",
    "brms.glm.categ.stdformul.without_wobs.without_offs",
    "brms.glmm.cumul.stdformul.without_wobs.with_offs",
    "brms.glmm.categ.stdformul.without_wobs.without_offs"
  )
  if (!use_fac) {
    # rstanarm::stan_polr() cannot deal with a numeric response:
    sel_fits <- grep("^rstanarm\\.glm\\.cumul\\.", sel_fits, value = TRUE,
                     invert = TRUE)
  } else {
    # The non-multilevel (and non-additive) brms::cumulative() case is
    # redundant, given the corresponding rstanarm::stan_polr() case and the
    # multilevel (and non-additive) brms::cumulative() case:
    sel_fits <- grep("^brms\\.glm\\.cumul\\.", sel_fits, value = TRUE,
                     invert = TRUE)
  }
  args_fit <- args_fit[names(args_fit) %in% sel_fits]
  if (run_brms && "rstanarm" %in% pkg_nms) {
    stopifnot(setequal(names(args_fit), sel_fits))
  } else if (run_brms && !"rstanarm" %in% pkg_nms) {
    stopifnot(setequal(names(args_fit),
                       grep("^brms\\.", sel_fits, value = TRUE)))
  } else {
    stopifnot(setequal(names(args_fit),
                       grep("^brms\\.", sel_fits, value = TRUE, invert = TRUE)))
  }
}

## Run --------------------------------------------------------------------

fits <- suppressWarnings(lapply(args_fit, function(args_fit_i) {
  fit_fun_nm <- get_fit_fun_nm(args_fit_i)
  if (args_fit_i$pkg_nm == "rstanarm" && args_fit_i$fam_nm == "cumul") {
    fit_fun_nm <- "stan_polr"
    args_fit_i$family <- NULL
    args_fit_i$prior <- quote(rstanarm::R2(location = 0.5, what = "median"))
    args_fit_i$QR <- NULL
  }
  ### Option 1:
  # do.call(fit_fun_nm,
  #         excl_nonargs(args_fit_i),
  #         envir = as.environment(asNamespace(args_fit_i$pkg_nm)))
  ###
  ### Option 2:
  do.call(get(fit_fun_nm, asNamespace(args_fit_i$pkg_nm)),
          excl_nonargs(args_fit_i))
  ###
}))

# projpred ----------------------------------------------------------------

## Setup ------------------------------------------------------------------

seed_tst <- 20411346
seed2_tst <- 866028
seed3_tst <- 1208499

nclusters_tst <- 2L
nclusters_pred_tst <- 3L
if (!run_more) {
  ndr_ncl_pred_tst <- list()
} else {
  ndr_ncl_pred_tst <- list(default_ndr_ncl = list())
}
ndr_ncl_pred_tst <- c(ndr_ncl_pred_tst, list(
  noclust = list(ndraws = nclusters_pred_tst),
  clust = list(nclusters = nclusters_pred_tst),
  clust1 = list(nclusters = 1L)
))
if (any(unlist(lapply(ndr_ncl_pred_tst, "[[", "ndraws")) <= 20)) {
  # Suppress the message concerning small `ndraws` or `ndraws_pred` values:
  options(projpred.mssg_ndraws = FALSE)
}
nresample_clusters_tst <- c(1L, 100L)

meth_tst <- list(
  default_meth = list(),
  L1 = list(method = "L1"),
  forward = list(method = "forward")
)
# Suppress the warning for interaction terms being selected before all involved
# main effects have been selected (only concerns L1 search):
options(projpred.warn_L1_interactions = FALSE)
# Suppress the warning thrown by proj_predict() in case of observation weights
# that are not all equal to `1`:
options(projpred.warn_wobs_ppd = FALSE)
# Suppress the verbose-mode progress bar in project():
options(projpred.verbose_project = FALSE)
# Suppress instability warnings:
options(projpred.warn_instable_projections = FALSE)

search_trms_tst <- list(
  default_search_trms = list(),
  alltrms = list(search_terms = setdiff(trms_common, "offset(offs_col)")),
  fixed = list(search_terms = c("xco.1", "xco.1 + xco.2", "xco.1 + xco.3",
                                "xco.1 + xco.2 + xco.3")),
  excluded = list(search_terms = c("xco.2", "xco.3", "xco.2 + xco.3")),
  empty_size = list(search_terms = c("xco.1 + xco.2", "xco.1 + xco.3",
                                     "xco.2 + xco.3", "xco.1 + xco.2 + xco.3"))
)

K_tst <- 2L
cvmeth_tst <- list(
  default_cvmeth = list(),
  LOO = list(cv_method = "LOO"),
  kfold = list(cv_method = "kfold", K = K_tst)
)

resp_oscale_tst <- list(
  default_r_oscale = list(),
  r_oscale_F = list(resp_oscale = FALSE)
)

vsel_funs <- nlist("summary.vsel", "plot.vsel", "suggest_size.vsel")
# Performance statistics common across all families, when using the traditional
# projection (or the latent projection with `resp_oscale = FALSE` or the latent
# projection with `resp_oscale = TRUE`, but the latter only in combination with
# `<refmodel>$family$cats` being `NULL`):
stats_common <- c("elpd", "mlpd", "mse", "rmse")
# Performance statistics for the binomial() family only, when using the
# traditional projection (or the latent projection with `resp_oscale = TRUE`,
# but the latter only in combination with `<refmodel>$family$cats` being
# `NULL`):
stats_binom <- c(stats_common, "acc", "auc")
# For creating test setups:
stats_tst <- list(
  default_stats = list(),
  common_stats = list(stats = stats_common),
  binom_stats = list(stats = stats_binom),
  augdat_stats = list(stats = c("elpd", "mlpd", "acc"))
)
type_tst <- c("mean", "lower", "upper", "se")

rk_abbv_tst <- list(
  default_abbv = list(),
  abbv3 = list(ranking_abbreviate = TRUE,
               ranking_abbreviate_args = list(minlength = 3))
)

rk_repel_tst <- list(
  default_repel = list(),
  repelText = list(ranking_repel = "text",
                   ranking_repel_args = list(seed = seed3_tst))
)

rk_col_tst <- as.list(setNames(nm = c(FALSE, TRUE)))
names(rk_col_tst) <- paste0("col", names(rk_col_tst))

cumulate_tst <- as.list(setNames(nm = c(FALSE, TRUE)))
names(cumulate_tst) <- paste0("cu", names(cumulate_tst))

angle_tst <- list(
  default_angle = list(),
  angle45 = list(text_angle = 45)
)

### nterms ----------------------------------------------------------------

ntermss <- sapply(mod_nms, function(mod_nm) {
  get(paste("nterms", mod_nm, sep = "_"))
})
# The `nterms_max` setting which will be used throughout the tests, except for
# the special `search_terms` tests:
nterms_max_tst <- min(ntermss)
if (!run_more) {
  nterms_max_tst <- min(nterms_max_tst, 2L)
}

nterms_unavail <- list(
  single = nterms_max_tst + 130L,
  vec = c(nterms_max_tst + 130L, nterms_max_tst + 290L)
)
if (!run_more) {
  nterms_avail <- list()
} else {
  nterms_avail <- list(default_nterms = NULL)
}
nterms_avail <- c(nterms_avail, list(
  empty = 0L,
  single = nterms_max_tst %/% 2L,
  subvec = as.integer(round(seq(0, nterms_max_tst, length.out = 2))),
  full = 0:nterms_max_tst
))

nterms_max_smmry <- list(
  default_nterms_max_smmry = NULL,
  halfway = nterms_max_tst %/% 2L,
  zero = 0L
)

nterms_max_rk <- list(
  default_nterms_max_rk = list(),
  halfway = list(nterms_max = nterms_max_tst %/% 2L),
  zero = list(nterms_max = 0L)
)

rk_max_tst <- list(
  default_rk_max = list(),
  rk_max_NA = list(ranking_nterms_max = NA),
  rk_max_1 = list(ranking_nterms_max = 1L)
)

## Reference model --------------------------------------------------------

args_ref <- lapply(setNames(nm = names(fits)), function(tstsetup_fit) {
  if (args_fit[[tstsetup_fit]]$pkg_nm == "brms" &&
      packageVersion("brms") >= "2.16.4") {
    pkg_args <- list(brms_seed = seed2_tst)
  } else {
    pkg_args <- list()
  }

  if (args_fit[[tstsetup_fit]]$fam_nm == "brnll") {
    # In this case, test the latent projection and the augmented-data projection
    # (the latter only in a slightly reduced set of test setups). For the
    # corresponding traditional projection case (needed for comparison), use the
    # keyword `trad_compare` to be able to find this case more easily later.
    augdat_args <- list(
      trad_compare = list(),
      latent = list(latent = TRUE) # , latent_y_unqs = c("0", "1")
    )
    if (!args_fit[[tstsetup_fit]]$mod_nm %in% c("gam", "gamm")) {
      augdat_args <- c(augdat_args, list(
        augdat = list(augdat_y_unqs = c("0", "1"),
                      augdat_link = quote(augdat_link_binom),
                      augdat_ilink = quote(augdat_ilink_binom))
      ))
    }
  } else if (args_fit[[tstsetup_fit]]$fam_nm %in% fam_nms_aug) {
    augdat_args <- list(augdat = list())
    if (args_fit[[tstsetup_fit]]$fam_nm %in% "cumul") {
      augdat_args <- c(augdat_args, list(
        latent = list(latent = TRUE)
      ))
    }
  } else {
    augdat_args <- list(trad = list())
  }

  lapply(setNames(nm = names(augdat_args)), function(augdat_args_nm) {
    return(c(nlist(tstsetup_fit),
             only_nonargs(args_fit[[tstsetup_fit]]),
             list(prj_nm = augdat_args_nm),
             pkg_args,
             augdat_args[[augdat_args_nm]]))
  })
})
args_ref <- unlist_cust(args_ref)

refmods <- lapply(args_ref, function(args_ref_i) {
  do.call(get_refmodel, c(
    list(object = fits[[args_ref_i$tstsetup_fit]]),
    excl_nonargs(args_ref_i)
  ))
})

## Variable selection -----------------------------------------------------

### varsel() --------------------------------------------------------------

if (run_vs) {
  # Some families are not supported yet, apart from the creation of a `refmodel`
  # object:
  tstsetups_vs_ref <- grep(fam_nms_unsupp_regex, names(refmods), value = TRUE,
                           invert = TRUE)
  if (!run_more) {
    tstsetups_vs_ref <- grep(paste0("\\.glmm", fam_nms_aug_regex),
                             tstsetups_vs_ref, value = TRUE, invert = TRUE)
  }
  tstsetups_vs_ref <- setNames(nm = tstsetups_vs_ref)
  args_vs <- lapply(tstsetups_vs_ref, function(tstsetup_ref) {
    mod_crr <- args_ref[[tstsetup_ref]]$mod_nm
    fam_crr <- args_ref[[tstsetup_ref]]$fam_nm
    prj_crr <- args_ref[[tstsetup_ref]]$prj_nm
    if (prj_crr == "trad" && mod_crr == "glm" && fam_crr == "gauss") {
      # Here, we test the default `method` (which is L1 search here) as well as
      # forward search:
      meth <- meth_tst[setdiff(names(meth_tst), "L1")]
    } else if (prj_crr %in% c("trad_compare", "latent")) {
      # For traditional settings which correspond to an augmented-data setting,
      # choose forward search (needed for comparing the two approaches);
      # correspondingly, we also need forward search for the latent projection
      # (even though in principle, the latent projection can be used with L1
      # search):
      meth <- meth_tst["forward"]
    } else {
      # Here, we only test the default `method`:
      meth <- meth_tst["default_meth"]
    }
    lapply(meth, function(meth_i) {
      if (mod_crr == "glm" && fam_crr == "gauss" &&
          grepl("\\.stdformul\\.", tstsetup_ref) &&
          identical(meth_i$method, "forward")) {
        # Here, we also test non-NULL `search_terms`:
        search_trms <- search_trms_tst
      } else {
        search_trms <- search_trms_tst["default_search_trms"]
      }
      lapply(search_trms, function(search_trms_i) {
        if (length(search_trms_i) &&
            !identical(search_trms_i$search_terms,
                       search_trms_tst$alltrms$search_terms)) {
          nterms_max_tst <- count_terms_chosen(search_trms_i$search_terms) - 1L
        }
        if (mod_crr == "glmm" && fam_crr == "categ") {
          # Quick-and-dirty solution to get some working results (it's probably
          # due to unfortunate test data simulated here that convergence at the
          # default settings is not given):
          extra_args <- list(avoid.increase = TRUE)
        } else {
          extra_args <- list()
        }
        return(c(
          nlist(tstsetup_ref), only_nonargs(args_ref[[tstsetup_ref]]),
          list(
            nclusters = nclusters_tst, nclusters_pred = nclusters_pred_tst,
            nterms_max = nterms_max_tst, verbose = FALSE, seed = seed_tst
          ),
          meth_i, search_trms_i, extra_args
        ))
      })
    })
  })
  args_vs <- unlist_cust(args_vs)
  stopifnot(sum(sapply(args_vs, function(args_vs_i) {
    !is.null(args_vs_i$search_terms)
  })) >= 1)

  vss <- lapply(args_vs, function(args_vs_i) {
    if (args_vs_i$prj_nm == "augdat" && args_vs_i$fam_nm == "cumul") {
      warn_expected <- "non-integer #successes in a binomial glm!"
    } else if (!is.null(args_vs_i$avoid.increase)) {
      warn_expected <- warn_mclogit
    } else {
      warn_expected <- NA
    }
    expect_warning(
      vs_out <- do.call(varsel, c(
        list(object = refmods[[args_vs_i$tstsetup_ref]]),
        excl_nonargs(args_vs_i)
      )),
      warn_expected,
      info = args_vs_i$tstsetup_ref
    )
    return(vs_out)
  })
}

### cv_varsel() -----------------------------------------------------------

if (run_cvvs) {
  tstsetups_cvvs_ref <- tstsetups_vs_ref
  # Even in the `run_more = TRUE` case (which is not run by default), we need to
  # impose some restrictions to have the tests run through in a reasonable
  # amount of time:
  tstsetups_cvvs_ref <- grep(
    paste0("\\.glmm", fam_nms_aug_regex), tstsetups_cvvs_ref, value = TRUE,
    invert = TRUE
  )
  if (!run_more) {
    tstsetups_cvvs_ref <- grep("\\.gam\\.", tstsetups_cvvs_ref, value = TRUE,
                               invert = TRUE)
  }
  # Under the special test settings used here, Bernoulli GAMMs often seem to run
  # into lme4 errors. However, since these Bernoulli GAMMs are basically
  # redundant given the other tested models, we can simply skip them:
  tstsetups_cvvs_ref <- grep("\\.gamm\\.brnll\\.", tstsetups_cvvs_ref,
                             value = TRUE, invert = TRUE)
  tstsetups_cvvs_ref <- setNames(nm = tstsetups_cvvs_ref)
  args_cvvs <- lapply(tstsetups_cvvs_ref, function(tstsetup_ref) {
    pkg_crr <- args_ref[[tstsetup_ref]]$pkg_nm
    mod_crr <- args_ref[[tstsetup_ref]]$mod_nm
    fam_crr <- args_ref[[tstsetup_ref]]$fam_nm
    prj_crr <- args_ref[[tstsetup_ref]]$prj_nm
    if (prj_crr %in% c("trad_compare", "latent")) {
      # For traditional settings which correspond to an augmented-data setting
      # (`trad_compare`), choose forward search (needed for comparing the two
      # approaches; therefore also necessary for the `latent` setting even
      # though in principle, the latent projection can be used with L1 search):
      meth <- meth_tst["forward"]
    } else {
      meth <- meth_tst["default_meth"]
    }
    if (grepl("\\.without_wobs", tstsetup_ref)) {
      # In principle, we want to use K-fold CV here and LOO CV else because
      # rstanarm:::kfold.stanreg() doesn't support observation weights. However,
      # there are some special cases to take care of:
      if (pkg_crr == "brms" && packageVersion("brms") <= "2.16.3") {
        # For brms versions <= 2.16.3, there is a reproducibility issue when
        # using K-fold CV, so use LOO CV:
        cvmeth <- cvmeth_tst["default_cvmeth"]
      } else if (pkg_crr == "brms" && mod_crr == "gamm") {
        # For GAMMs fitted by brms, there is a (random, i.e., only occasional)
        # reproducibility issue when using K-fold CV, so use LOO CV:
        cvmeth <- cvmeth_tst["default_cvmeth"]
      } else if (prj_crr %in% c("latent", "augdat") && fam_crr != "brnll") {
        # We also want to test the latent and the augmented-data projection with
        # LOO CV:
        cvmeth <- cvmeth_tst["default_cvmeth"]
      } else {
        cvmeth <- cvmeth_tst["kfold"]
      }
    } else {
      cvmeth <- cvmeth_tst["default_cvmeth"]
    }
    lapply(meth, function(meth_i) {
      lapply(cvmeth, function(cvmeth_i) {
        if (!run_valsearch_always && !identical(cvmeth_i$cv_method, "kfold") &&
            # Handle augmented-data and corresponding traditional projection:
            (!prj_crr %in% c("latent", "augdat", "trad_compare") ||
             (prj_crr %in% c("latent", "augdat", "trad_compare") &&
              !run_valsearch_aug_lat)) &&
            # Forward search:
            ((length(meth_i) == 0 &&
              (mod_crr != "glm" || prj_crr == "augdat")) ||
             (length(meth_i) > 0 && meth_i$method == "forward"))) {
          # These are cases with forward search, LOO CV, and
          # `!run_valsearch_always` where we want to save time by using
          # `validate_search = FALSE`:
          meth_i <- c(meth_i, list(validate_search = FALSE))
        }
        search_trms <- search_trms_tst["default_search_trms"]
        lapply(search_trms, function(search_trms_i) {
          if (length(search_trms_i) &&
              !identical(search_trms_i$search_terms,
                         search_trms_tst$alltrms$search_terms)) {
            nterms_max_tst <- count_terms_chosen(search_trms_i$search_terms) -
              1L
          }
          return(c(
            nlist(tstsetup_ref), only_nonargs(args_ref[[tstsetup_ref]]),
            list(
              nclusters = nclusters_tst, nclusters_pred = nclusters_pred_tst,
              nterms_max = nterms_max_tst, verbose = FALSE, seed = seed_tst
            ),
            meth_i, cvmeth_i, search_trms_i
          ))
        })
      })
    })
  })
  args_cvvs <- unlist_cust(args_cvvs)

  # Use suppressWarnings() because of occasional warnings concerning Pareto k
  # diagnostics. Additionally to suppressWarnings(), suppressMessages() could be
  # used here (because of the refits in K-fold CV):
  cvvss <- suppressWarnings(lapply(args_cvvs, function(args_cvvs_i) {
    cvvs_expr <- expression(do.call(cv_varsel, c(
      list(object = refmods[[args_cvvs_i$tstsetup_ref]]),
      excl_nonargs(args_cvvs_i)
    )))
    if (args_cvvs_i$mod_nm == "gamm" &&
        !identical(args_cvvs_i$cv_method, "kfold")) {
      # Due to issue #239, we have to wrap the call to cv_varsel() in try():
      return(try(eval(cvvs_expr), silent = TRUE))
    } else {
      return(eval(cvvs_expr))
    }
  }))
  success_cvvs <- !sapply(cvvss, inherits, "try-error")
  err_ok <- sapply(cvvss[!success_cvvs], function(cvvs_err) {
    attr(cvvs_err, "condition")$message ==
      "Not enough (non-NA) data to do anything meaningful"
  })
  expect_true(
    all(err_ok),
    info = paste("Unexpected error for",
                 paste(names(cvvss[!success_cvvs])[!err_ok], collapse = ", "))
  )
  cvvss <- cvvss[success_cvvs]
  args_cvvs <- args_cvvs[success_cvvs]
}

## Projection -------------------------------------------------------------

### From "refmodel" -------------------------------------------------------

if (run_prj) {
  # Some families are not supported yet, apart from the creation of a `refmodel`
  # object:
  tstsetups_prj_ref <- grep(fam_nms_unsupp_regex, names(refmods), value = TRUE,
                            invert = TRUE)
  tstsetups_prj_ref <- setNames(nm = tstsetups_prj_ref)
  args_prj <- lapply(tstsetups_prj_ref, function(tstsetup_ref) {
    pkg_crr <- args_ref[[tstsetup_ref]]$pkg_nm
    mod_crr <- args_ref[[tstsetup_ref]]$mod_nm
    fam_crr <- args_ref[[tstsetup_ref]]$fam_nm
    prj_crr <- args_ref[[tstsetup_ref]]$prj_nm
    if (grepl("\\.spclformul", tstsetup_ref)) {
      solterms_x <- solterms_spcl
    }
    solterms <- nlist(empty = character(), solterms_x)
    if (prj_crr %in% c("augdat", "trad_compare") && fam_crr == "brnll" &&
        mod_crr == "glmm") {
      # We need a single group-level term (which only consists of group-level
      # intercepts) to be able to use `nAGQ` later:
      solterms_z <- setdiff(solterms_z, "(xco.1 | z.1)")
    } else {
      solterms_z <- setdiff(solterms_z, "(1 | z.1)")
    }
    if (mod_crr %in% c("glmm", "gamm")) {
      solterms <- c(solterms,
                    nlist(solterms_z, solterms_xz = c(solterms_x, solterms_z)))
    }
    if (mod_crr %in% c("gam", "gamm")) {
      solterms <- c(solterms,
                    nlist(solterms_s, solterms_xs = c(solterms_x, solterms_s)))
    }
    if (mod_crr == "gamm") {
      solterms <- c(solterms,
                    nlist(solterms_sz = c(solterms_s, solterms_z),
                          solterms_xsz = c(solterms_x, solterms_s, solterms_z)))
    }
    if (!run_more &&
        (fam_crr != "gauss" || grepl("\\.spclformul", tstsetup_ref))) {
      solterms <- tail(solterms, 1)
    }
    lapply(setNames(nm = names(solterms)), function(solterms_nm_i) {
      if (pkg_crr == "rstanarm" && mod_crr == "glm" &&
          fam_crr == "gauss" && solterms_nm_i == "solterms_x") {
        ndr_ncl_pred <- ndr_ncl_pred_tst
      } else if (pkg_crr == "rstanarm" && mod_crr == "glm" &&
                 fam_crr == "gauss" && solterms_nm_i == "empty") {
        ndr_ncl_pred <- ndr_ncl_pred_tst[c("noclust", "clust", "clust1")]
      } else if (
        (run_more && (
          (pkg_crr == "rstanarm" && mod_crr == "glmm" &&
           fam_crr == "brnll" && solterms_nm_i == "solterms_xz") ||
          (pkg_crr == "rstanarm" && mod_crr == "gam" &&
           fam_crr == "binom" && solterms_nm_i == "solterms_xs") ||
          (pkg_crr == "rstanarm" && mod_crr == "gamm" &&
           fam_crr == "brnll" && solterms_nm_i == "solterms_xsz")
        )) ||
        (!run_more && mod_crr %in% c("glmm", "gam", "gamm")) ||
        prj_crr %in% c("latent", "augdat", "trad_compare")
      ) {
        # The `noclust` setting is important for the test "non-clustered
        # projection does not require a seed" in `test_project.R`.
        ndr_ncl_pred <- ndr_ncl_pred_tst[c("noclust", "clust")]
      } else {
        ndr_ncl_pred <- ndr_ncl_pred_tst[c("clust")]
      }
      if (prj_crr %in% c("augdat", "trad_compare") && fam_crr == "brnll" &&
          mod_crr == "glmm" && grepl("z", solterms_nm_i)) {
        # We need an increased accuracy to be able to compare traditional and
        # augmented-data projection:
        divmin_args <- list(nAGQ = 30L)
      } else {
        divmin_args <- list()
      }
      lapply(ndr_ncl_pred, function(ndr_ncl_pred_i) {
        if (mod_crr == "glmm" && fam_crr == "categ") {
          # Quick-and-dirty solution to get some working results (it's probably
          # due to unfortunate test data simulated here that convergence at the
          # default settings is not given):
          divmin_args <- c(divmin_args, list(avoid.increase = TRUE))
        }
        return(c(
          nlist(tstsetup_ref), only_nonargs(args_ref[[tstsetup_ref]]),
          list(solution_terms = solterms[[solterms_nm_i]], seed = seed_tst),
          ndr_ncl_pred_i, divmin_args
        ))
      })
    })
  })
  args_prj <- unlist_cust(args_prj)

  prjs <- lapply(args_prj, function(args_prj_i) {
    if (args_prj_i$prj_nm == "augdat" && args_prj_i$fam_nm == "cumul" &&
        !any(grepl("\\|", args_prj_i$solution_terms))) {
      warn_expected <- "non-integer #successes in a binomial glm!"
    } else if (!is.null(args_prj_i$avoid.increase) &&
               any(grepl("\\|", args_prj_i$solution_terms))) {
      warn_expected <- warn_mclogit
    } else {
      warn_expected <- NA
    }
    expect_warning(
      prj_out <- do.call(project, c(
        list(object = refmods[[args_prj_i$tstsetup_ref]]),
        excl_nonargs(args_prj_i)
      )),
      warn_expected,
      info = args_prj_i$tstsetup_ref
    )
    return(prj_out)
  })
}

### From "vsel" -----------------------------------------------------------

# A helper function to create the argument list for project() for a given
# character vector of test setups (referring to either `vss` or `cvvss`):
cre_args_prj_vsel <- function(tstsetups_prj_vsel) {
  vsel_type <- deparse(substitute(tstsetups_prj_vsel))
  args_obj <- switch(vsel_type,
                     "tstsetups_prj_vs" = args_vs,
                     "tstsetups_prj_cvvs" = args_cvvs,
                     stop("Unexpected `vsel_type`."))
  lapply(tstsetups_prj_vsel, function(tstsetup_vsel) {
    args_out <- c(
      nlist(tstsetup_vsel), only_nonargs(args_obj[[tstsetup_vsel]]),
      list(nclusters = nclusters_pred_tst, seed = seed_tst)
    )
    if (args_obj[[tstsetup_vsel]]$mod_nm != "glm" ||
        !is.null(args_obj[[tstsetup_vsel]]$search_terms) ||
        args_obj[[tstsetup_vsel]]$prj_nm %in% c("latent", "augdat") ||
        grepl("\\.spclformul", tstsetup_vsel)) {
      nterms_avail <- nterms_avail["subvec"]
    }
    if (!is.null(args_obj[[tstsetup_vsel]]$search_terms)) {
      nterms_max_cut <- args_obj[[tstsetup_vsel]]$nterms_max
      if (all(grepl("\\+", args_obj[[tstsetup_vsel]]$search_terms))) {
        # This is the "empty_size" setting, so we have to subtract the skipped
        # model size (see issue #307):
        nterms_max_cut <- nterms_max_cut - 1L
      }
      nterms_avail <- lapply(nterms_avail, function(nterms_avail_i) {
        pmin(nterms_avail_i, nterms_max_cut)
      })
    }
    lapply(nterms_avail, function(nterms_crr) {
      if (!is.null(nterms_crr)) {
        args_out <- c(args_out, list(nterms = nterms_crr))
      }
      return(args_out)
    })
  })
}

#### varsel() -------------------------------------------------------------

if (run_vs) {
  tstsetups_prj_vs <- unlist(lapply(mod_nms, function(mod_nm) {
    if (any(grepl(paste0("\\.", mod_nm, "\\.gauss\\."), names(vss)))) {
      tstsetups_out <- grep(
        paste0("\\.", mod_nm, "\\.gauss\\..*\\.default_meth"), names(vss),
        value = TRUE
      )
    } else {
      tstsetups_out <- grep(
        paste0("\\.", mod_nm, "\\..*\\.default_meth"), names(vss),
        value = TRUE
      )
    }
    if (!run_more) {
      tstsetups_out <- head(tstsetups_out, 1)
    }
    return(tstsetups_out)
  }))
  tstsetups_prj_vs <- union(
    tstsetups_prj_vs,
    grep("\\.default_search_trms", names(vss), value = TRUE, invert = TRUE)
  )
  tstsetups_prj_vs <- union(
    tstsetups_prj_vs,
    grep("\\.(latent|augdat)\\.", names(vss), value = TRUE)
  )
  tstsetups_prj_vs <- setNames(nm = tstsetups_prj_vs)
  stopifnot(length(tstsetups_prj_vs) > 0)
  args_prj_vs <- cre_args_prj_vsel(tstsetups_prj_vs)
  args_prj_vs <- unlist_cust(args_prj_vs)

  prjs_vs <- lapply(args_prj_vs, function(args_prj_vs_i) {
    do.call(project, c(
      list(object = vss[[args_prj_vs_i$tstsetup_vsel]]),
      excl_nonargs(args_prj_vs_i)
    ))
  })
}

#### cv_varsel() ----------------------------------------------------------

if (run_cvvs) {
  tstsetups_prj_cvvs <- unlist(lapply(mod_nms, function(mod_nm) {
    if (any(grepl(paste0("\\.", mod_nm, "\\.gauss\\."), names(cvvss)))) {
      tstsetups_out <- grep(
        paste0("\\.", mod_nm,
               "\\.gauss\\..*\\.default_meth\\.default_cvmeth"),
        names(cvvss), value = TRUE
      )
    } else {
      tstsetups_out <- grep(
        paste0("\\.", mod_nm, "\\..*\\.default_meth\\.default_cvmeth"),
        names(cvvss), value = TRUE
      )
    }
    if (!run_more) {
      tstsetups_out <- head(tstsetups_out, 1)
    }
    return(tstsetups_out)
  }))
  tstsetups_prj_cvvs <- union(
    tstsetups_prj_cvvs,
    grep("\\.default_search_trms", names(cvvss), value = TRUE, invert = TRUE)
  )
  tstsetups_prj_cvvs <- union(
    tstsetups_prj_cvvs,
    grep("\\.(latent|augdat)\\.", names(cvvss), value = TRUE)
  )
  tstsetups_prj_cvvs <- setNames(nm = tstsetups_prj_cvvs)
  stopifnot(length(tstsetups_prj_cvvs) > 0)
  args_prj_cvvs <- cre_args_prj_vsel(tstsetups_prj_cvvs)
  args_prj_cvvs <- unlist_cust(args_prj_cvvs)

  # Use suppressWarnings() because of occasional pwrssUpdate() warnings:
  prjs_cvvs <- suppressWarnings(lapply(
    args_prj_cvvs,
    function(args_prj_cvvs_i) {
      do.call(project, c(
        list(object = cvvss[[args_prj_cvvs_i$tstsetup_vsel]]),
        excl_nonargs(args_prj_cvvs_i)
      ))
    }
  ))
}

## Prediction -------------------------------------------------------------

### From "projection" -----------------------------------------------------

if (run_prj) {
  pls <- lapply(prjs, proj_linpred, .seed = seed2_tst)
  pps <- lapply(prjs, proj_predict, .seed = seed2_tst)
}

### From "proj_list" ------------------------------------------------------

#### varsel() -------------------------------------------------------------

if (run_vs) {
  pls_vs <- lapply(prjs_vs, proj_linpred, .seed = seed2_tst)
  pps_vs <- lapply(prjs_vs, proj_predict, .seed = seed2_tst)
}

#### cv_varsel() ----------------------------------------------------------

if (run_cvvs) {
  pls_cvvs <- lapply(prjs_cvvs, proj_linpred, .seed = seed2_tst)
  pps_cvvs <- lapply(prjs_cvvs, proj_predict, .seed = seed2_tst)
}

## summary.vsel() ---------------------------------------------------------

cre_args_smmry_vsel <- function(args_obj) {
  tstsetups <- names(args_obj)
  # Choose all test setups which are for special `search_terms` settings:
  tstsetups_smmry_vsel <- tstsetups[sapply(tstsetups, function(tstsetup_vsel) {
    !is.null(args_obj[[tstsetup_vsel]]$search_terms)
  })]
  # Choose all test setups which are for the latent or the augmented-data
  # projection or which are corresponding to such:
  tstsetups_smmry_vsel <- union(
    tstsetups_smmry_vsel,
    tstsetups[sapply(tstsetups, function(tstsetup_vsel) {
      args_obj[[tstsetup_vsel]]$prj_nm %in% c("latent", "augdat",
                                              "trad_compare")
    })]
  )

  # Ensure that from each model type (`mod_nm`) and each family (`fam_nm`), we
  # have at least one test setup:
  mods_fams_existing <- sapply(tstsetups_smmry_vsel, function(tstsetup_vsel) {
    paste0(args_obj[[tstsetup_vsel]]$mod_nm, ".",
           args_obj[[tstsetup_vsel]]$fam_nm)
  })
  mods_fams_possible <- apply(expand.grid(mod_nms, fam_nms), 1, paste,
                              collapse = ".")
  mods_fams_missing <- setdiff(mods_fams_possible, mods_fams_existing)
  tstsetups_smmry_vsel <- union(
    tstsetups_smmry_vsel,
    unlist(lapply(mods_fams_missing, function(mod_fam) {
      head(
        grep(paste0(".", mod_fam, "."), tstsetups, value = TRUE, fixed = TRUE),
        1
      )
    }))
  )

  tstsetups_smmry_vsel <- setNames(nm = tstsetups_smmry_vsel)
  stopifnot(length(tstsetups_smmry_vsel) > 0)
  lapply(tstsetups_smmry_vsel, function(tstsetup_vsel) {
    mod_crr <- args_obj[[tstsetup_vsel]]$mod_nm
    fam_crr <- args_obj[[tstsetup_vsel]]$fam_nm
    prj_crr <- args_obj[[tstsetup_vsel]]$prj_nm
    add_stats <- switch(mod_crr,
                        "glm" = switch(prj_crr,
                                       "augdat" = "augdat_stats",
                                       "latent" = "augdat_stats",
                                       switch(fam_crr,
                                              "brnll" = "binom_stats",
                                              "binom" = "binom_stats",
                                              "common_stats")),
                        character())
    if (!run_more && !is.null(args_obj[[tstsetup_vsel]]$search_terms)) {
      add_stats <- character()
    }
    stats_tst <- stats_tst[c("default_stats", add_stats)]
    lapply(stats_tst, function(stats_crr) {
      if (!run_more) {
        if (!is.null(args_obj[[tstsetup_vsel]]$search_terms)) {
          nterms_tst <- nterms_max_smmry["default_nterms_max_smmry"]
        } else {
          nterms_tst <- nterms_max_smmry["halfway"]
        }
      } else {
        if (mod_crr == "glm" && fam_crr == "gauss" &&
            is.null(args_obj[[tstsetup_vsel]]$search_terms) &&
            length(stats_crr) == 0) {
          nterms_tst <- nterms_max_smmry[c("default_nterms_max_smmry",
                                           "halfway", "zero")]
        } else {
          nterms_tst <- nterms_max_smmry["default_nterms_max_smmry"]
        }
      }
      if (length(stats_crr)) {
        cumulate_tst <- cumulate_tst["cuFALSE"]
      }
      lapply(nterms_tst, function(nterms_crr) {
        if (prj_crr != "latent") {
          resp_oscale_tst <- resp_oscale_tst["default_r_oscale"]
        }
        lapply(resp_oscale_tst, function(resp_oscale_crr) {
          if (isFALSE(resp_oscale_crr$resp_oscale) &&
              any(setdiff(stats_binom, stats_common) %in% stats_crr$stats)) {
            return(dummy_glob)
          }
          lapply(cumulate_tst, function(cumulate_crr) {
            return(c(
              nlist(tstsetup_vsel), only_nonargs(args_obj[[tstsetup_vsel]]),
              list(type = type_tst, nterms_max = nterms_crr),
              stats_crr, resp_oscale_crr, list(cumulate = cumulate_crr)
            ))
          })
        })
      })
    })
  })
}

### varsel() --------------------------------------------------------------

if (run_vs) {
  args_smmry_vs <- cre_args_smmry_vsel(args_vs)
  args_smmry_vs <- unlist_cust(args_smmry_vs)
  args_smmry_vs <- rm_dummies(args_smmry_vs)
  has_zero_vs <- any(sapply(lapply(args_smmry_vs, "[[", "nterms_max"),
                            identical, 0L))

  smmrys_vs <- lapply(args_smmry_vs, function(args_smmry_vs_i) {
    if (any(c("rmse", "auc") %in% args_smmry_vs_i$stats)) {
      smmry_seed <- list(seed = seed3_tst)
    } else {
      smmry_seed <- list()
    }
    do.call(summary, c(
      list(object = vss[[args_smmry_vs_i$tstsetup_vsel]]),
      excl_nonargs(args_smmry_vs_i),
      smmry_seed
    ))
  })
}

### cv_varsel() -----------------------------------------------------------

if (run_cvvs) {
  args_smmry_cvvs <- cre_args_smmry_vsel(args_cvvs)
  args_smmry_cvvs <- unlist_cust(args_smmry_cvvs)
  args_smmry_cvvs <- rm_dummies(args_smmry_cvvs)
  has_zero_cvvs <- any(sapply(lapply(args_smmry_cvvs, "[[", "nterms_max"),
                              identical, 0L))

  smmrys_cvvs <- lapply(args_smmry_cvvs, function(args_smmry_cvvs_i) {
    if (any(c("rmse", "auc") %in% args_smmry_cvvs_i$stats)) {
      smmry_seed <- list(seed = seed3_tst)
    } else {
      smmry_seed <- list()
    }
    do.call(summary, c(
      list(object = cvvss[[args_smmry_cvvs_i$tstsetup_vsel]]),
      excl_nonargs(args_smmry_cvvs_i),
      smmry_seed
    ))
  })
}

if (run_more) {
  has_zero_combined <- logical()
  if (run_vs) {
    has_zero_combined <- c(has_zero_combined, has_zero_vs)
  }
  if (run_cvvs) {
    has_zero_combined <- c(has_zero_combined, has_zero_cvvs)
  }
  if (length(has_zero_combined)) {
    stopifnot(any(has_zero_combined))
  }
}

## plot.vsel() ------------------------------------------------------------

cre_args_plot_vsel <- function(args_obj) {
  tstsetups <- grep("\\.brnll\\..*\\.trad", names(args_obj), value = TRUE)
  tstsetups <- union(
    tstsetups,
    grep("\\.default_search_trms|\\.alltrms", names(args_obj), value = TRUE,
         invert = TRUE)
  )
  tstsetups <- union(
    tstsetups,
    head(grep("\\.spclformul", names(args_obj), value = TRUE), 1)
  )
  lapply(
    setNames(nm = tstsetups),
    function(tstsetup_vsel) {
      nterms_max_plot <- nterms_max_smmry[c("default_nterms_max_smmry",
                                            "halfway")]
      lapply(nterms_max_plot, function(nterms_crr) {
        lapply(rk_max_tst, function(rk_max_crr) {
          lapply(rk_abbv_tst, function(rk_abbv_crr) {
            lapply(rk_repel_tst, function(rk_repel_crr) {
              lapply(rk_col_tst, function(rk_col_crr) {
                lapply(cumulate_tst, function(cumulate_crr) {
                  lapply(angle_tst, function(angle_crr) {
                    return(c(
                      nlist(tstsetup_vsel),
                      only_nonargs(args_obj[[tstsetup_vsel]]),
                      list(nterms_max = nterms_crr),
                      rk_max_crr, rk_abbv_crr, rk_repel_crr,
                      list(ranking_colored = rk_col_crr,
                           cumulate = cumulate_crr),
                      angle_crr
                    ))
                  })
                })
              })
            })
          })
        })
      })
    })
}

### varsel() --------------------------------------------------------------

if (run_vs) {
  args_plot_vs <- cre_args_plot_vsel(args_vs)
  args_plot_vs <- unlist_cust(args_plot_vs)

  plots_vs <- lapply(args_plot_vs, function(args_plot_vs_i) {
    do.call(plot, c(
      list(x = vss[[args_plot_vs_i$tstsetup_vsel]]),
      excl_nonargs(args_plot_vs_i)
    ))
  })
}

### cv_varsel() -----------------------------------------------------------

if (run_cvvs) {
  args_plot_cvvs <- cre_args_plot_vsel(args_cvvs)
  args_plot_cvvs <- unlist_cust(args_plot_cvvs)

  plots_cvvs <- lapply(args_plot_cvvs, function(args_plot_cvvs_i) {
    do.call(plot, c(
      list(x = cvvss[[args_plot_cvvs_i$tstsetup_vsel]]),
      excl_nonargs(args_plot_cvvs_i)
    ))
  })
}

## ranking() --------------------------------------------------------------

### varsel() --------------------------------------------------------------

if (run_vs) {
  args_rk_vs <- lapply(setNames(nm = names(vss)), function(tstsetup_vsel) {
    lapply(nterms_max_rk, function(nterms_crr) {
      return(c(
        nlist(tstsetup_vsel), only_nonargs(args_vs[[tstsetup_vsel]]),
        nterms_crr
      ))
    })
  })
  args_rk_vs <- unlist_cust(args_rk_vs)

  rks_vs <- lapply(args_rk_vs, function(args_rk_vs_i) {
    do.call(ranking, c(
      list(object = vss[[args_rk_vs_i$tstsetup_vsel]]),
      excl_nonargs(args_rk_vs_i)
    ))
  })
}

### cv_varsel() -----------------------------------------------------------

if (run_cvvs) {
  args_rk_cvvs <- lapply(setNames(nm = names(cvvss)), function(tstsetup_vsel) {
    lapply(nterms_max_rk, function(nterms_crr) {
      return(c(
        nlist(tstsetup_vsel), only_nonargs(args_cvvs[[tstsetup_vsel]]),
        nterms_crr
      ))
    })
  })
  args_rk_cvvs <- unlist_cust(args_rk_cvvs)

  rks_cvvs <- lapply(args_rk_cvvs, function(args_rk_cvvs_i) {
    do.call(ranking, c(
      list(object = cvvss[[args_rk_cvvs_i$tstsetup_vsel]]),
      excl_nonargs(args_rk_cvvs_i)
    ))
  })
}

## cv_proportions() -------------------------------------------------------

err_no_foldwise_rk <- "Could not find fold-wise predictor rankings"

### varsel() --------------------------------------------------------------

if (run_vs) {
  args_pr_vs <- lapply(setNames(nm = names(rks_vs)), function(tstsetup_rk) {
    return(c(
      nlist(tstsetup_rk), only_nonargs(args_rk_vs[[tstsetup_rk]])
    ))
  })
  args_pr_vs <- unlist_cust(args_pr_vs)

  prs_vs <- lapply(args_pr_vs, function(args_pr_vs_i) {
    err_expected <- err_no_foldwise_rk
    expect_error(
      do.call(cv_proportions, c(
        list(object = rks_vs[[args_pr_vs_i$tstsetup_rk]]),
        excl_nonargs(args_pr_vs_i)
      )),
      err_expected,
      info = args_pr_vs_i$tstsetup_rk
    )
    return(dummy_glob)
  })
  keep_prs_vs <- rm_dummies(prs_vs, return_logical = TRUE)
  prs_vs <- prs_vs[keep_prs_vs]
  args_pr_vs <- args_pr_vs[keep_prs_vs]
}

### cv_varsel() -----------------------------------------------------------

if (run_cvvs) {
  args_pr_cvvs <- lapply(setNames(nm = names(rks_cvvs)), function(tstsetup_rk) {
    lapply(cumulate_tst, function(cumulate_crr) {
      return(c(
        nlist(tstsetup_rk), only_nonargs(args_rk_cvvs[[tstsetup_rk]]),
        list(cumulate = cumulate_crr)
      ))
    })
  })
  args_pr_cvvs <- unlist_cust(args_pr_cvvs)

  prs_cvvs <- lapply(args_pr_cvvs, function(args_pr_cvvs_i) {
    if (isFALSE(args_cvvs[[args_pr_cvvs_i$tstsetup_vsel]]$validate_search)) {
      err_expected <- err_no_foldwise_rk
    } else if (isTRUE(
      args_rk_cvvs[[args_pr_cvvs_i$tstsetup_rk]]$nterms_max == 0
    )) {
      err_expected <- "Needing `nterms_max >= 1`"
    } else {
      err_expected <- NA
    }
    expect_error(
      pr_out <- do.call(cv_proportions, c(
        list(object = rks_cvvs[[args_pr_cvvs_i$tstsetup_rk]]),
        excl_nonargs(args_pr_cvvs_i)
      )),
      err_expected,
      info = args_pr_cvvs_i$tstsetup_rk
    )
    if (is.na(err_expected)) {
      return(pr_out)
    } else {
      return(dummy_glob)
    }
  })
  keep_prs_cvvs <- rm_dummies(prs_cvvs, return_logical = TRUE)
  prs_cvvs <- prs_cvvs[keep_prs_cvvs]
  args_pr_cvvs <- args_pr_cvvs[keep_prs_cvvs]
}

## plot.cv_proportions() --------------------------------------------------

if (run_cvvs) {
  args_plotpr <- lapply(setNames(nm = names(prs_cvvs)), function(tstsetup_pr) {
    return(c(
      nlist(tstsetup_pr), only_nonargs(args_pr_cvvs[[tstsetup_pr]])
    ))
  })
  args_plotpr <- unlist_cust(args_plotpr)

  plotprs <- lapply(args_plotpr, function(args_plotpr_i) {
    do.call(plot, c(
      list(x = prs_cvvs[[args_plotpr_i$tstsetup_pr]]),
      excl_nonargs(args_plotpr_i)
    ))
  })
}

## Output names -----------------------------------------------------------

vsel_nms <- c(
  "refmodel", "nobs_train", "search_path", "solution_terms",
  "solution_terms_cv", "ce", "type_test", "y_wobs_test", "nobs_test",
  "summaries", "nterms_all", "nterms_max", "method", "cv_method", "K",
  "validate_search", "clust_used_search", "clust_used_eval", "nprjdraws_search",
  "nprjdraws_eval", "projpred_version"
)
# Related to prediction (in contrast to selection):
vsel_nms_pred <- c("summaries", "solution_terms", "ce")
vsel_nms_pred_opt <- c("solution_terms")
# Related to `nloo`:
vsel_nms_nloo <- c("summaries", "solution_terms_cv")
vsel_nms_nloo_opt <- c("solution_terms_cv")
# Related to `validate_search`:
vsel_nms_valsearch <- c("validate_search", "summaries", "ce",
                        "solution_terms_cv")
vsel_nms_valsearch_opt <- character()
# Related to `cvfits`:
vsel_nms_cvfits <- c("refmodel", "summaries", "solution_terms_cv")
vsel_nms_cvfits_opt <- c("solution_terms_cv")
vsel_smmrs_sub_nms <- vsel_smmrs_ref_nms <- c("mu", "lppd")

## Defaults ---------------------------------------------------------------

nclusters_default <- 20L
ndraws_pred_default <- 400L
nresample_clusters_default <- 1000L
regul_default <- 1e-4

# Seed --------------------------------------------------------------------

if (run_randRNG) {
  set.seed(NULL)
}
