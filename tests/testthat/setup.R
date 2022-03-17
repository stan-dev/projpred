#_______________________________________________________________________________
# Setup for the unit tests
#_______________________________________________________________________________

# General setup -----------------------------------------------------------

options(warn = 1)

# These switches may be set to `FALSE` to save time (e.g., when debugging
# interactively):
# Run more tests, at the downside of increased runtime?:
run_more <- FALSE
# Run varsel()?:
run_vs <- identical(Sys.getenv("NOT_CRAN"), "true")
# Run cv_varsel()?:
run_cvvs <- run_vs
# Run cv_varsel() with `validate_search = TRUE` always (`TRUE`) or just for L1
# search (`FALSE`)?:
run_valsearch_always <- FALSE
# Run the `cvfits` test for all possible test setups (`TRUE`) or just for the
# first one among the GLMMs (`FALSE`; note that if there is no GLMM available in
# that test, the first test setup among those for K-fold CV is used)?:
run_cvfits_all <- FALSE
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
if (run_snaps) {
  testthat_ed_max2 <- edition_get() <= 2
}
# Run parallel tests?:
# Notes:
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
      stop("Package \"doParallel\" is needed for these tests. Please ",
           "install it.",
           call. = FALSE)
    }
    dopar_backend <- "doParallel"
  } else {
    # This case (which should not be possible by default) is only included
    # here to demonstrate how parallelization should be used on Windows (but
    # currently, this makes no sense, as explained above).
    if (!requireNamespace("doFuture", quietly = TRUE)) {
      stop("Package \"doFuture\" is needed for these tests. Please ",
           "install it.",
           call. = FALSE)
    }
    dopar_backend <- "doFuture"
    if (identical(.Platform$OS.type, "windows")) {
      ### Not used in this case because the 'future.callr' package provides a
      ### faster alternative on Windows (which is still slower than a sequential
      ### run, though):
      # future_plan <- "multisession"
      ###
      if (!requireNamespace("future.callr", quietly = TRUE)) {
        stop("Package \"future.callr\" is needed for these tests. Please ",
             "install it.",
             call. = FALSE)
      }
      future_plan <- "callr"
    } else {
      # This case (which should not be possible by default) is only included
      # here to demonstrate how other systems should be used with the 'doFuture'
      # package.
      future_plan <- "multicore"
    }
  }
}
# Run all test scripts (following this setup script) in a completely random RNG
# state? (The tests should still pass then, because in all situations where RNG
# is used, a specific seed is supposed to be set.):
run_randRNG <- identical(Sys.getenv("NOT_CRAN"), "true")
# Run tests for additive models (GAMs and GAMMs)?:
run_additive <- TRUE

source(testthat::test_path("helpers", "unlist_cust.R"), local = TRUE)
source(testthat::test_path("helpers", "testers.R"), local = TRUE)
source(testthat::test_path("helpers", "args.R"), local = TRUE)
source(testthat::test_path("helpers", "getters.R"), local = TRUE)
source(testthat::test_path("helpers", "formul_handlers.R"), local = TRUE)

mod_nms <- setNames(nm = c("glm", "glmm", "gam", "gamm"))
if (run_additive) {
  # Suppress the warning for additive models (GAMs and GAMMs) stating that their
  # implementation is currently only experimental:
  options(projpred.warn_additive_experimental = FALSE)
} else {
  mod_nms <- setNames(nm = setdiff(mod_nms, c("gam", "gamm")))
}

fam_nms <- setNames(nm = c("gauss", "brnll", "binom", "poiss"))

# Data --------------------------------------------------------------------

## Setup ------------------------------------------------------------------

# Number of observations:
nobsv <- 41L
# Values for testing:
nobsv_tst <- c(1L, nobsv %/% 3L)

# Seed:
seed_dat <- 8541351
set.seed(seed_dat)

## GLMs --------------------------------------------------------------------
## Add nonpooled ("fixed") effects to the intercept-(and-offset-)only model

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
## Add partially pooled ("random") effects to the GLMs

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
  r_icpts <- rnorm(nlvl_ran_i, sd = 2.8)
  r_xco1 <- rnorm(nlvl_ran_i, sd = 2.8)
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
dat <- lapply(mod_nms, function(mod_nm) {
  lapply(fam_nms, function(fam_nm) {
    pred_link <- get(paste0("eta_", mod_nm))
    if (fam_nm != "brnll" && !mod_nm %in% c("gam", "gamm")) {
      # For the "brnll" `fam_nm`, offsets are simply not added to have some
      # scenarios without offsets.
      # For GAMs, offsets are not added because of rstanarm issue #546 (see
      # also further below).
      # For GAMMs, offsets are not added because of rstanarm issue #253 (see
      # also further below).
      pred_link <- pred_link + offs_tst
    }
    pred_resp <- get(paste0("f_", fam_nm))$linkinv(pred_link)
    if (fam_nm == "gauss") {
      return(rnorm(nobsv, mean = pred_resp, sd = dis_tst))
    } else if (fam_nm == "brnll") {
      return(rbinom(nobsv, 1, pred_resp))
    } else if (fam_nm == "binom") {
      return(rbinom(nobsv, wobs_tst, pred_resp))
    } else if (fam_nm == "poiss") {
      return(rpois(nobsv, pred_resp))
    } else {
      stop("Unknown `fam_nm`.")
    }
  })
})
dat <- unlist(dat, recursive = FALSE)
names(dat) <- paste("y", gsub("\\.", "_", names(dat)), sep = "_")
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

## nterms -----------------------------------------------------------------

ntermss <- sapply(mod_nms, function(mod_nm) {
  get(paste("nterms", mod_nm, sep = "_"))
})
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
  subvec = as.integer(round(seq(0, nterms_max_tst, length.out = 3))),
  full = 0:nterms_max_tst
))

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

# Fits --------------------------------------------------------------------

## Setup ------------------------------------------------------------------

if (!requireNamespace("rstanarm", quietly = TRUE)) {
  stop("Package \"rstanarm\" is needed for these tests. Please install it.",
       call. = FALSE)
}
pkg_nms <- "rstanarm"

if (run_brms) {
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop("Package \"brms\" is needed for these tests. Please install it.",
         call. = FALSE)
  }
  pkg_nms <- c(pkg_nms, "brms")
  # For storing "brmsfit"s locally:
  file_pth <- testthat::test_path("bfits")
  if (!dir.exists(file_pth)) dir.create(file_pth)
}
pkg_nms <- setNames(nm = pkg_nms)

chains_tst <- 2L
iter_tst <- 500L
nrefdraws <- chains_tst * iter_tst %/% 2L
seed_fit <- 74345

### Formula ---------------------------------------------------------------

# Notes:
#   * Argument `offset` has an issue for rstanarm::stan_glmer() (see rstanarm
#     issue #541). Instead, use offset() in the formula.
#   * Argument `offset` is not supported by rstanarm::stan_gamm4(). Instead, use
#     offset() in the formula. However, because of rstanarm issue #546 and
#     rstanarm issue #253, omit the offsets in GAMs and GAMMs.
#   * In rstanarm::stan_gamm4(), multilevel terms are specified via argument
#     `random`.

trms_common <- c("xco.1", "xco.2", "xco.3", "xca.1", "xca.2",
                 "offset(offs_col)")
trms_grp <- c("(xco.1 | z.1)")
trms_add <- c("s(s.1)") # , "s(s.2)", "s(s.3)"
trms_common_spcl <- c("xco.1", "I(xco.1^2)",
                      "exp(xco.2) * I(as.numeric(xco.3 > 0))", "xca.1", "xca.2",
                      "offset(offs_col)")

# Solution terms for project()-ing from `"refmodel"`s:
solterms_x <- c("xco.2", "xca.1")
solterms_z <- c("(1 | z.1)", "(xco.1 | z.1)")
solterms_s <- c("s(s.1)") # , "s(s.2)"
solterms_spcl <- c("xca.1", "xco.1", "I(xco.1^2)", "exp(xco.2)",
                   "I(as.numeric(xco.3 > 0))",
                   "exp(xco.2):I(as.numeric(xco.3 > 0))")

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
for (obj_symb_chr in c(paste0("f_", fam_nms))) {
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
      if (mod_nm != "gamm") {
        random_arg <- list()
      } else {
        random_arg <- list(random = as.formula(paste("~", trms_grp)))
      }
    }

    if (mod_nm != "glm") {
      if (pkg_nm == "brms") {
        # For speed reasons, do not test all families:
        if (mod_nm == "glmm") {
          fam_nms <- intersect(fam_nms, "brnll")
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

      if (pkg_nm == "brms" && fam_nm == "brnll") {
        family_crr <- quote(get("bernoulli", envir = asNamespace("brms"))())
      } else {
        family_crr <- as.name(paste0("f_", fam_nm))
      }

      if (fam_nm == "brnll" ||
          (pkg_nm == "rstanarm" && mod_nm %in% c("gam", "gamm"))) {
        # For the "brnll" `fam_nm`, the offsets are simply omitted to have some
        # scenarios without offsets.
        # In the rstanarm "gam" and "gamm" case, the offsets are omitted because
        # of rstanarm issue #546 and rstanarm issue #253.
        offss_nms <- "without_offs"
      } else {
        offss_nms <- "with_offs"
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

        if (fam_nm %in% c("brnll", "binom")) {
          # In this case, observation weights are not supported by projpred (and
          # for rstanarm, the number of trials is specified in the formula via
          # the cbind() syntax (indirectly, because the number of trials is the
          # sum of the two columns)):
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
                               file_refit = "on_change") # , silent = 2
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

if (!run_more) {
  sel_fits <- c(
    "rstanarm.glm.gauss.stdformul.with_wobs.with_offs",
    "rstanarm.glm.brnll.stdformul.without_wobs.without_offs",
    "rstanarm.glmm.gauss.spclformul.with_wobs.with_offs",
    "rstanarm.gam.gauss.spclformul.with_wobs.without_offs",
    "rstanarm.gamm.brnll.stdformul.without_wobs.without_offs",
    "brms.glm.poiss.stdformul.with_wobs.with_offs",
    "brms.glmm.brnll.stdformul.without_wobs.without_offs",
    # "brms.gam.binom.stdformul.without_wobs.with_offs",
    "brms.gamm.binom.stdformul.without_wobs.with_offs"
  )
  args_fit <- args_fit[names(args_fit) %in% sel_fits]
}

## Run --------------------------------------------------------------------

fits <- suppressWarnings(lapply(args_fit, function(args_fit_i) {
  fit_fun_nm <- switch(args_fit_i$pkg_nm,
                       "rstanarm" = switch(args_fit_i$mod_nm,
                                           "glm" = "stan_glm",
                                           "glmm" = "stan_glmer",
                                           "stan_gamm4"),
                       "brms" = "brm",
                       stop("Unknown `pkg_nm`."))
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

K_tst <- 2L
cvmeth_tst <- list(
  default_cvmeth = list(),
  LOO = list(cv_method = "LOO"),
  kfold = list(cv_method = "kfold", K = K_tst)
)

vsel_funs <- nlist("summary.vsel", "plot.vsel", "suggest_size.vsel")
stats_common <- c("elpd", "mlpd", "mse", "rmse")
stats_tst <- list(
  default_stats = list(),
  common_stats = list(stats = stats_common),
  binom_stats = list(stats = c(stats_common, "acc", "auc"))
)
type_tst <- c("mean", "lower", "upper", "se")

seed_tst <- 74341
seed2_tst <- 866028
seed3_tst <- 1208499

## Reference model --------------------------------------------------------

args_ref <- lapply(setNames(nm = names(fits)), function(tstsetup_fit) {
  if (args_fit[[tstsetup_fit]]$pkg_nm == "brms" &&
      packageVersion("brms") >= "2.16.4") {
    pkg_args <- list(brms_seed = seed2_tst)
  } else {
    pkg_args <- list()
  }
  return(c(nlist(tstsetup_fit),
           only_nonargs(args_fit[[tstsetup_fit]]),
           pkg_args))
})

refmods <- lapply(args_ref, function(args_ref_i) {
  do.call(get_refmodel, c(
    list(object = fits[[args_ref_i$tstsetup_fit]]),
    excl_nonargs(args_ref_i)
  ))
})

## Variable selection -----------------------------------------------------

### varsel() --------------------------------------------------------------

if (run_vs) {
  tstsetups_vs_ref <- setNames(nm = names(refmods))
  args_vs <- lapply(tstsetups_vs_ref, function(tstsetup_ref) {
    mod_crr <- args_ref[[tstsetup_ref]]$mod_nm
    fam_crr <- args_ref[[tstsetup_ref]]$fam_nm
    if (mod_crr == "glm" && fam_crr == "gauss") {
      # Here, we test the default `method` (which is L1 search here) as well as
      # forward search:
      meth <- meth_tst[setdiff(names(meth_tst), "L1")]
    } else {
      # Here, we only test the default `method`:
      meth <- meth_tst["default_meth"]
    }
    lapply(meth, function(meth_i) {
      return(c(
        nlist(tstsetup_ref), only_nonargs(args_ref[[tstsetup_ref]]),
        list(
          nclusters = nclusters_tst, nclusters_pred = nclusters_pred_tst,
          nterms_max = nterms_max_tst, verbose = FALSE, seed = seed_tst
        ),
        meth_i
      ))
    })
  })
  args_vs <- unlist_cust(args_vs)

  vss <- lapply(args_vs, function(args_vs_i) {
    do.call(varsel, c(
      list(object = refmods[[args_vs_i$tstsetup_ref]]),
      excl_nonargs(args_vs_i)
    ))
  })
}

### cv_varsel() -----------------------------------------------------------

if (run_cvvs) {
  tstsetups_cvvs_ref <- names(refmods)
  if (!run_more) {
    tstsetups_cvvs_ref <- grep("\\.gam\\.", tstsetups_cvvs_ref, value = TRUE,
                               invert = TRUE)
  }
  tstsetups_cvvs_ref <- setNames(nm = tstsetups_cvvs_ref)
  args_cvvs <- lapply(tstsetups_cvvs_ref, function(tstsetup_ref) {
    pkg_crr <- args_ref[[tstsetup_ref]]$pkg_nm
    mod_crr <- args_ref[[tstsetup_ref]]$mod_nm
    fam_crr <- args_ref[[tstsetup_ref]]$fam_nm
    meth <- meth_tst["default_meth"]
    if (grepl("\\.without_wobs", tstsetup_ref)) {
      # In principle, we want to use K-fold CV here and LOO CV else because
      # rstanarm:::kfold.stanreg() doesn't support observation weights. However,
      # there are some special cases to take care of:
      if (mod_crr == "gamm" && fam_crr == "brnll") {
        # In this case, K-fold CV leads to an error in pwrssUpdate()
        # ("(maxstephalfit) PIRLS step-halvings failed to reduce deviance in
        # pwrssUpdate"). Therefore, use LOO CV:
        cvmeth <- cvmeth_tst["default_cvmeth"]
        # TODO (GAMMs): Fix this.
      } else if (pkg_crr == "brms" && packageVersion("brms") <= "2.16.3") {
        # For brms versions <= 2.16.3, there is a reproducibility issue when
        # using K-fold CV, so use LOO CV:
        cvmeth <- cvmeth_tst["default_cvmeth"]
      } else if (pkg_crr == "brms" && mod_crr == "gamm") {
        # For GAMMs fitted by brms, there is a (random, i.e., only occasional)
        # reproducibility issue when using K-fold CV, so use LOO CV:
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
            ((length(meth_i) == 0 && mod_crr != "glm") ||
             (length(meth_i) > 0 && meth_i$method == "forward"))) {
          # In this case, we have forward search (and LOO CV) and
          # `!run_valsearch_always` indicates that we want to save time by using
          # `validate_search = FALSE`:
          meth_i <- c(meth_i, list(validate_search = FALSE))
        }
        return(c(
          nlist(tstsetup_ref), only_nonargs(args_ref[[tstsetup_ref]]),
          list(
            nclusters = nclusters_tst, nclusters_pred = nclusters_pred_tst,
            nterms_max = nterms_max_tst, verbose = FALSE, seed = seed_tst
          ),
          meth_i, cvmeth_i
        ))
      })
    })
  })
  args_cvvs <- unlist_cust(args_cvvs)

  # Use suppressWarnings() because of occasional warnings concerning Pareto k
  # diagnostics. Additionally to suppressWarnings(), suppressMessages() could be
  # used here (because of the refits in K-fold CV):
  cvvss <- suppressWarnings(lapply(args_cvvs, function(args_cvvs_i) {
    do.call(cv_varsel, c(
      list(object = refmods[[args_cvvs_i$tstsetup_ref]]),
      excl_nonargs(args_cvvs_i)
    ))
  }))
}

## Projection -------------------------------------------------------------

### From "refmodel" -------------------------------------------------------

tstsetups_prj_ref <- setNames(nm = names(refmods))
args_prj <- lapply(tstsetups_prj_ref, function(tstsetup_ref) {
  pkg_crr <- args_ref[[tstsetup_ref]]$pkg_nm
  mod_crr <- args_ref[[tstsetup_ref]]$mod_nm
  fam_crr <- args_ref[[tstsetup_ref]]$fam_nm
  if (grepl("\\.spclformul", tstsetup_ref)) {
    solterms_x <- solterms_spcl
  }
  solterms <- nlist(empty = character(), solterms_x)
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
      (!run_more && mod_crr %in% c("glmm", "gam", "gamm"))
    ) {
      # The `noclust` setting is important for the test "non-clustered
      # projection does not require a seed" in `test_project.R`.
      ndr_ncl_pred <- ndr_ncl_pred_tst[c("noclust", "clust")]
    } else {
      ndr_ncl_pred <- ndr_ncl_pred_tst[c("clust")]
    }
    lapply(ndr_ncl_pred, function(ndr_ncl_pred_i) {
      return(c(
        nlist(tstsetup_ref), only_nonargs(args_ref[[tstsetup_ref]]),
        list(solution_terms = solterms[[solterms_nm_i]], seed = seed_tst),
        ndr_ncl_pred_i
      ))
    })
  })
})
args_prj <- unlist_cust(args_prj)

prjs <- lapply(args_prj, function(args_prj_i) {
  do.call(project, c(
    list(object = refmods[[args_prj_i$tstsetup_ref]]),
    excl_nonargs(args_prj_i)
  ))
})

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
        grepl("\\.spclformul", tstsetup_vsel)) {
      nterms_avail <- nterms_avail["subvec"]
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
  tstsetups_prj_vs <- setNames(
    nm = unlist(lapply(mod_nms, function(mod_nm) {
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
  )
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
  tstsetups_prj_cvvs <- setNames(
    nm = unlist(lapply(mod_nms, function(mod_nm) {
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
  )
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

pls <- lapply(prjs, proj_linpred)
pps <- lapply(prjs, proj_predict, .seed = seed2_tst)

### From "proj_list" ------------------------------------------------------

#### varsel() -------------------------------------------------------------

if (run_vs) {
  pls_vs <- lapply(prjs_vs, proj_linpred)
  pps_vs <- lapply(prjs_vs, proj_predict, .seed = seed2_tst)
}

#### cv_varsel() ----------------------------------------------------------

if (run_cvvs) {
  pls_cvvs <- lapply(prjs_cvvs, proj_linpred)
  pps_cvvs <- lapply(prjs_cvvs, proj_predict, .seed = seed2_tst)
}

## summary.vsel() ---------------------------------------------------------

cre_args_smmry_vsel <- function(tstsetups_smmry_vsel) {
  vsel_type <- deparse(substitute(tstsetups_smmry_vsel))
  args_obj <- switch(vsel_type,
                     "tstsetups_smmry_vs" = args_vs,
                     "tstsetups_smmry_cvvs" = args_cvvs,
                     stop("Unexpected `vsel_type`."))
  lapply(tstsetups_smmry_vsel, function(tstsetup_vsel) {
    mod_crr <- args_obj[[tstsetup_vsel]]$mod_nm
    fam_crr <- args_obj[[tstsetup_vsel]]$fam_nm
    add_stats <- switch(mod_crr,
                        "glm" = switch(fam_crr,
                                       "brnll" = "binom_stats",
                                       "binom" = "binom_stats",
                                       "common_stats"),
                        character())
    stats_tst <- stats_tst[c("default_stats", add_stats)]
    lapply(stats_tst, function(stats_crr) {
      if (!run_more) {
        nterms_tst <- nterms_avail["single"]
      } else {
        if (mod_crr == "glm" && fam_crr == "gauss" && length(stats_crr) == 0) {
          nterms_tst <- nterms_avail[c("default_nterms", "single")]
        } else {
          nterms_tst <- nterms_avail["default_nterms"]
        }
      }
      lapply(nterms_tst, function(nterms_crr) {
        return(c(
          nlist(tstsetup_vsel), only_nonargs(args_obj[[tstsetup_vsel]]),
          list(type = type_tst, nterms_max = nterms_crr),
          stats_crr
        ))
      })
    })
  })
}

### varsel() --------------------------------------------------------------

if (run_vs) {
  tstsetups_smmry_vs <- setNames(nm = unlist(lapply(mod_nms, function(mod_nm) {
    unlist(lapply(fam_nms, function(fam_nm) {
      head(grep(paste0("\\.", mod_nm, "\\.", fam_nm), names(vss), value = TRUE),
           1)
    }))
  })))
  stopifnot(length(tstsetups_smmry_vs) > 0)
  args_smmry_vs <- cre_args_smmry_vsel(tstsetups_smmry_vs)
  args_smmry_vs <- unlist_cust(args_smmry_vs)

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
  tstsetups_smmry_cvvs <- setNames(
    nm = unlist(lapply(mod_nms, function(mod_nm) {
      unlist(lapply(fam_nms, function(fam_nm) {
        head(grep(paste0("\\.", mod_nm, "\\.", fam_nm), names(cvvss),
                  value = TRUE),
             1)
      }))
    }))
  )
  stopifnot(length(tstsetups_smmry_cvvs) > 0)
  args_smmry_cvvs <- cre_args_smmry_vsel(tstsetups_smmry_cvvs)
  args_smmry_cvvs <- unlist_cust(args_smmry_cvvs)

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

## Output names -----------------------------------------------------------

vsel_nms <- c(
  "refmodel", "search_path", "d_test", "summaries", "solution_terms", "kl",
  "nterms_max", "nterms_all", "method", "cv_method", "validate_search",
  "clust_used_search", "clust_used_eval", "nprjdraws_search", "nprjdraws_eval",
  "suggested_size", "summary"
)
vsel_nms_cv <- c(
  "refmodel", "search_path", "d_test", "summaries", "kl", "solution_terms",
  "pct_solution_terms_cv", "nterms_all", "nterms_max", "method", "cv_method",
  "validate_search", "clust_used_search", "clust_used_eval", "nprjdraws_search",
  "nprjdraws_eval", "suggested_size", "summary"
)
# Related to prediction (in contrast to selection):
vsel_nms_pred <- c("summaries", "solution_terms", "kl", "suggested_size",
                   "summary")
vsel_nms_pred_opt <- c("solution_terms", "suggested_size")
# Related to `d_test`:
vsel_nms_dtest <- c("d_test", setdiff(vsel_nms_pred, c("solution_terms", "kl")))
# Related to `nloo`:
vsel_nms_cv_nloo <- c("summaries", "pct_solution_terms_cv", "suggested_size",
                      "summary")
vsel_nms_cv_nloo_opt <- c("pct_solution_terms_cv", "suggested_size")
# Related to `validate_search`:
vsel_nms_cv_valsearch <- c("validate_search", "summaries",
                           "pct_solution_terms_cv", "suggested_size",
                           "summary")
vsel_nms_cv_valsearch_opt <- c("suggested_size")
# Related to `cvfits`:
vsel_nms_cv_cvfits <- c("refmodel", "d_test", "summaries",
                        "pct_solution_terms_cv", "summary", "suggested_size")
vsel_nms_cv_cvfits_opt <- c("pct_solution_terms_cv", "suggested_size")
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
