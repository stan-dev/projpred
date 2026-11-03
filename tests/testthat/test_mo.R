context("monotonic terms mo()")

# Helpers -----------------------------------------------------------------

# Imitates brms::mo(), which returns a list that model.frame() cannot use. The
# reference model formulas below have this function in their environment (as
# if brms was attached), so the tests also check that projpred uses its own
# mo().
mo_brms_like <- function(x, id = NA) {
  structure(list(term = deparse(substitute(x))),
            class = c("mo_term", "sp_term"))
}

# Simulated data and draws of a reference model with the brms parametrization
# b * D * sum(zeta[seq_len(j)]) of monotonic effects (`j` is the position of
# the value in the levels, starting at 0): `x` is an integer with levels
# x_min, ..., x_min + D, `f` an ordered factor with 3 levels, `z` and `w` are
# continuous (`w` has no effect), and `g` is a grouping variable.
mo_sim <- function(n = 120, S = 40, D = 4, x_min = 0L, seed = 5602) {
  set.seed(seed)
  dat <- data.frame(
    x = sample(rep_len(x_min + 0:D, n)),
    f = factor(sample(rep_len(c("lo", "mid", "hi"), n)),
               levels = c("lo", "mid", "hi"), ordered = TRUE),
    z = rnorm(n),
    w = rnorm(n),
    g = factor(sample(rep_len(paste0("gr", 1:6), n)))
  )
  zeta_x <- rev(seq_len(D)) / sum(seq_len(D))
  zeta_f <- c(0.7, 0.3)
  b <- cbind(alpha = rnorm(S, 0.2, 0.05), bx = rnorm(S, 0.8, 0.05),
             bf = rnorm(S, -0.6, 0.05), bz = rnorm(S, 0.5, 0.05))
  # Linear predictor (N x S) for data `nd`:
  eta_fun <- function(nd) {
    mo_x <- D * c(0, cumsum(zeta_x))[nd$x - x_min + 1]
    mo_f <- 2 * c(0, cumsum(zeta_f))[as.integer(nd$f)]
    cbind(1, mo_x, mo_f, nd$z) %*% t(b)
  }
  eta_mean <- rowMeans(eta_fun(dat))
  dat$y_gauss <- eta_mean + rnorm(n, sd = 0.5)
  dat$y_binom <- rbinom(n, 1, plogis(eta_mean))
  dat$y_poiss <- rpois(n, exp(eta_mean))
  # The coefficients of the cumulative indicators (S x number of columns):
  coef_true <- cbind(
    "(Intercept)" = b[, "alpha"],
    b[, "bx"] * D * matrix(zeta_x, S, D, byrow = TRUE,
                           dimnames = list(NULL,
                                           paste0("mo(x)", x_min + seq_len(D)))),
    b[, "bf"] * 2 * matrix(zeta_f, S, 2, byrow = TRUE,
                           dimnames = list(NULL, c("mo(f)mid", "mo(f)hi"))),
    z = b[, "bz"]
  )
  return(nlist(dat, eta_fun, coef_true, S, D))
}

# A reference model for `sim` (see mo_sim()). `ref_predfun_arr` turns the N x S
# linear predictor into the output of `ref_predfun` (for the augmented-data
# projection).
mo_refm <- function(sim, fml, family = gaussian(), data = sim$dat,
                    ref_predfun_arr = identity, ...) {
  environment(fml) <- list2env(list(mo = mo_brms_like),
                               parent = environment(fml))
  ref_predfun <- function(fit, newdata = NULL) {
    nd <- if (is.null(newdata)) fit$data else newdata
    ref_predfun_arr(fit$sim$eta_fun(nd))
  }
  dis <- if (family$family == "gaussian") rep(0.5, sim$S) else NULL
  cvfun <- function(folds) {
    lapply(seq_len(max(folds)), function(k) {
      list(sim = sim, data = data[folds != k, , drop = FALSE],
           omitted = which(folds == k))
    })
  }
  dot_args <- list(...)
  cvrefbuilder <- function(cvfit) {
    do.call(init_refmodel, c(
      list(object = cvfit, data = cvfit$data, formula = fml, family = family,
           ref_predfun = ref_predfun, dis = dis, cvfun = cvfun,
           called_from_cvrefbuilder = TRUE),
      dot_args
    ))
  }
  init_refmodel(
    list(sim = sim, data = data), data = data, formula = fml, family = family,
    ref_predfun = ref_predfun, dis = dis, cvfun = cvfun,
    cvrefbuilder = cvrefbuilder, ...
  )
}

# Evaluate `expr` without the warnings whose message contains `msg`.
muffle_warning <- function(expr, msg) {
  withCallingHandlers(expr, warning = function(w) {
    if (grepl(msg, conditionMessage(w), fixed = TRUE)) {
      invokeRestart("muffleWarning")
    }
  })
}

# MASS::polr() warns about non-integer data. projpred filters this warning for
# the draw-wise fits, but the filter misses calls with `:` (here
# `MASS::polr()`), so the warning appears once per draw.
muffle_polr_warning <- function(expr) {
  muffle_warning(expr, "non-integer #successes in a binomial glm!")
}

# The training rows of `nd` (a subset of the training data):
rows_of <- function(nd) as.integer(rownames(nd))

# `a` without attributes other than the dimensions:
strip_attr <- function(a) {
  if (is.null(dim(a))) {
    return(as.vector(a))
  }
  return(array(as.vector(a), dim = dim(a)))
}

# New data that lack some levels of `x` and `f`:
mo_newdata <- function(sim) {
  dat <- sim$dat
  dat[dat$x <= 2 & dat$f != "hi", , drop = FALSE][1:8, ]
}

# mo() --------------------------------------------------------------------

test_that("mo() gives the cumulative indicators of an integer predictor", {
  x <- c(0, 2, 1, 3, 2)
  b <- mo(x)
  expect_identical(dim(b), c(5L, 3L))
  expect_identical(colnames(b), c("1", "2", "3"))
  expect_identical(attr(b, "mo_levels"), 0:3)
  expect_identical(strip_attr(b), outer(x, 1:3, ">=") * 1)
  # With the levels of the training data:
  b_new <- mo(c(1, 0), mo_levels = 0:3)
  expect_identical(strip_attr(b_new), rbind(c(1, 0, 0), c(0, 0, 0)))
  # As in brms, the levels start at the smallest value (also if negative):
  b_shift <- mo(c(3, 1, 2))
  expect_identical(attr(b_shift, "mo_levels"), 1:3)
  expect_identical(strip_attr(b_shift), rbind(c(1, 1), c(0, 0), c(1, 0)))
  b_neg <- mo(c(0, -1, 1))
  expect_identical(colnames(b_neg), c("0", "1"))
  expect_identical(strip_attr(b_neg), rbind(c(1, 0), c(0, 0), c(1, 1)))
  # Missing values give missing rows:
  expect_true(all(is.na(mo(c(1, NA, 2))[2, ])))
})

test_that("mo() gives the cumulative indicators of an ordered factor", {
  f <- factor(c("b", "a", "c", "b"), levels = c("a", "b", "c"), ordered = TRUE)
  b <- mo(f)
  expect_identical(colnames(b), c("b", "c"))
  expect_identical(attr(b, "mo_levels"), c("a", "b", "c"))
  expect_identical(strip_attr(b), cbind(c(1, 0, 1, 1), c(0, 0, 1, 0)))
  # A factor with dropped levels uses the levels of the training data:
  b_new <- mo(droplevels(f[2]), mo_levels = c("a", "b", "c"))
  expect_identical(dim(b_new), c(1L, 2L))
})

test_that("mo() throws errors for unsupported predictors", {
  expect_error(mo(factor(c("a", "b"))), "ordered factor")
  expect_error(mo(c(0, 1.5)), "integer vector")
  expect_error(mo(c("a", "b")), "integer vector")
  expect_error(mo(c(0, 0)), "at least two levels")
  expect_error(mo(c(0, 4), mo_levels = 0:3), "outside the range")
  expect_error(mo(c(-1, 2), mo_levels = 0:3), "outside the range")
  f_new <- factor("d", levels = c("a", "d"), ordered = TRUE)
  expect_error(mo(f_new, mo_levels = c("a", "b", "c")), "do not exist")
  expect_error(mo(c(0, Inf)), "integer vector")
  expect_error(mo(c(0, Inf), mo_levels = 0:3), "integer vector")
})

test_that("mo() matches new data by the type of the training levels", {
  lvls <- c("lo", "mid", "hi")
  b_ref <- strip_attr(mo(factor(c("mid", "hi"), levels = lvls,
                                ordered = TRUE)))
  # Character vectors and unordered factors are matched by their labels:
  expect_identical(strip_attr(mo(c("mid", "hi"), mo_levels = lvls)), b_ref)
  expect_identical(strip_attr(mo(factor(c("mid", "hi")), mo_levels = lvls)),
                   b_ref)
  # Integer codes are not accepted for an ordered factor, and factors are not
  # accepted for an integer:
  expect_error(mo(c(2L, 3L), mo_levels = lvls), "factor or a character")
  expect_error(mo(factor(1), mo_levels = 0:3), "must be an integer in new")
})

test_that("mo() terms in unsupported positions throw an error", {
  sim <- mo_sim(n = 60, S = 4)
  expect_error(mo_refm(sim, y_gauss ~ z + (mo(x) | g)),
               "not supported in group-level terms")
  expect_error(mo_refm(sim, y_gauss ~ mo(x) + s(z) + (1 | g)),
               "not supported in additive multilevel models")
})

test_that(paste(
  "collapse_ranked_predictors() does not change mo() terms by the",
  "replacements for factor terms"
), {
  dat <- data.frame(y = rnorm(6), gbx = c(0:2, 0:2),
                    g = factor(rep(c("a", "b"), 3)))
  fml <- y ~ mo(gbx) + g
  environment(fml) <- list2env(list(mo = mo), parent = globalenv())
  fr <- model.frame(fml, data = dat)
  expect_identical(
    collapse_ranked_predictors(c("mo(gbx)1", "gb", "mo(gbx)2"),
                               formula = fml, data = fr),
    c("mo(gbx)", "g")
  )
})

test_that("the terms of mo() store the training levels in \"predvars\"", {
  dat <- data.frame(x = c(0, 1, 3, 2), y = 1:4)
  mf <- model.frame(y ~ mo(x), data = dat)
  pv <- attr(attr(mf, "terms"), "predvars")
  expect_identical(pv[[3]][[1]], quote(projpred:::mo))
  expect_identical(pv[[3]]$mo_levels, 0:3)
  # The stored call works from the global environment (as for mgcv::gam(),
  # which sets the environment of the terms to the global environment):
  expect_equal(dim(eval(pv[[3]], list(x = 1), globalenv())), c(1L, 3L))
})

test_that("mo_formula_env() makes `mo` refer to projpred's mo()", {
  env_user <- new.env()
  assign("mo", mo_brms_like, envir = env_user)
  assign("user_obj", 1, envir = env_user)
  fml <- y ~ mo(x) + z
  environment(fml) <- env_user
  fml_new <- mo_formula_env(fml)
  expect_identical(get("mo", envir = environment(fml_new)), mo)
  expect_identical(parent.env(environment(fml_new)), env_user)
  expect_identical(get("user_obj", envir = environment(fml_new)), 1)
  # A second call does not add another environment:
  expect_identical(environment(mo_formula_env(fml_new)),
                   environment(fml_new))
  # Formulas without mo() calls (also with a variable called `mo`) stay as
  # they are:
  fml_no <- y ~ mo + z
  environment(fml_no) <- env_user
  expect_identical(environment(mo_formula_env(fml_no)), env_user)
})

# Sign groups of the monotonicity constraint ------------------------------

# The groups of mo_sign_groups() as column names:
sign_groups_of <- function(fml, dat) {
  fml <- mo_formula_env(fml)
  fr <- model.frame(fml, data = dat)
  mm <- model.matrix(fml, data = fr)
  lapply(mo_sign_groups(fr, mm), function(g) colnames(mm)[g])
}

test_that("mo_sign_groups() gives one group per mo() column set", {
  dat <- data.frame(
    y = rnorm(12), x = rep(1:4, 3), x2 = rep(c(5, 6, 7), 4),
    bin = rep(0:1, 6), z = seq(-1, 1, length.out = 12),
    f = factor(rep(c("a", "b", "c"), each = 4)),
    o = factor(rep(c("lo", "mid", "mid:x", "hi"), 3),
               levels = c("lo", "mid", "mid:x", "hi"), ordered = TRUE)
  )
  # Main effects; a binary predictor has no group:
  expect_identical(
    sign_groups_of(y ~ mo(x) + mo(bin) + z, dat),
    list("mo(x)" = c("mo(x)2", "mo(x)3", "mo(x)4"))
  )
  # Interaction with a numeric variable:
  expect_identical(
    sign_groups_of(y ~ z + mo(x2):z, dat),
    list("z:mo(x2)" = c("z:mo(x2)6", "z:mo(x2)7"))
  )
  # Interaction with a factor: one group per column of the factor:
  expect_identical(
    sign_groups_of(y ~ mo(x) * f, dat),
    list("mo(x)" = c("mo(x)2", "mo(x)3", "mo(x)4"),
         "mo(x):fb" = c("mo(x)2:fb", "mo(x)3:fb", "mo(x)4:fb"),
         "mo(x):fc" = c("mo(x)2:fc", "mo(x)3:fc", "mo(x)4:fc"))
  )
  expect_identical(
    sign_groups_of(y ~ f + mo(x):f:z, dat),
    list("fa:mo(x):z" = c("fa:mo(x)2:z", "fa:mo(x)3:z", "fa:mo(x)4:z"),
         "fb:mo(x):z" = c("fb:mo(x)2:z", "fb:mo(x)3:z", "fb:mo(x)4:z"),
         "fc:mo(x):z" = c("fc:mo(x)2:z", "fc:mo(x)3:z", "fc:mo(x)4:z"))
  )
  # Without the main effect of the factor (indicators of all levels), and with
  # a term whose variables are not in the order of the formula:
  expect_identical(
    sign_groups_of(y ~ z + f:mo(x) + z:mo(x2), dat),
    list("fa:mo(x)" = c("fa:mo(x)2", "fa:mo(x)3", "fa:mo(x)4"),
         "fb:mo(x)" = c("fb:mo(x)2", "fb:mo(x)3", "fb:mo(x)4"),
         "fc:mo(x)" = c("fc:mo(x)2", "fc:mo(x)3", "fc:mo(x)4"),
         "z:mo(x2)" = c("z:mo(x2)6", "z:mo(x2)7"))
  )
  # A binary predictor in an interaction with a factor has no group:
  expect_identical(sign_groups_of(y ~ mo(bin) * f, dat), list())
  # An ordered factor with `:` in a level (also with a variable `x`, so that
  # `mid:x` is both a level and a level-variable pair):
  expect_identical(
    sign_groups_of(y ~ mo(o) + z + mo(o):z, dat),
    list("mo(o)" = c("mo(o)mid", "mo(o)mid:x", "mo(o)hi"),
         "mo(o):z" = c("mo(o)mid:z", "mo(o)mid:x:z", "mo(o)hi:z"))
  )
  expect_identical(
    sign_groups_of(y ~ mo(o) + mo(o):x, dat),
    list("mo(o)" = c("mo(o)mid", "mo(o)mid:x", "mo(o)hi"),
         "mo(o):x" = c("mo(o)mid:x", "mo(o)mid:x:x", "mo(o)hi:x"))
  )
  # The indices refer to the model matrix (here with the intercept column):
  fml <- mo_formula_env(y ~ z + mo(x))
  fr <- model.frame(fml, data = dat)
  mm <- model.matrix(fml, data = fr)
  grps <- mo_sign_groups(fr, mm)
  expect_identical(attr(grps, "keys"), "mo(x)#")
  attr(grps, "keys") <- NULL
  expect_identical(grps, list("mo(x)" = 3:5))
  # The keys distinguish the codings of a factor in an interaction (the
  # difference to the first level, or the effect within a level), but not of
  # numeric variables:
  keys_of <- function(fml) {
    fml <- mo_formula_env(fml)
    fr <- model.frame(fml, data = dat)
    attr(mo_sign_groups(fr, model.matrix(fml, data = fr)), "keys")
  }
  expect_identical(keys_of(y ~ f * mo(x)),
                   c("mo(x)#", "fb:mo(x)#1", "fc:mo(x)#1"))
  expect_identical(keys_of(y ~ f + f:mo(x)),
                   c("fa:mo(x)#2", "fb:mo(x)#2", "fc:mo(x)#2"))
  expect_identical(keys_of(y ~ mo(x) * z), c("mo(x)#", "mo(x):z#"))
  expect_identical(keys_of(y ~ z + mo(x):z), "z:mo(x)#")
  # No mo() terms:
  expect_identical(sign_groups_of(y ~ z + f, dat), list())
})

test_that("mo_sign_groups() throws an error for unsupported mo() terms", {
  dat <- data.frame(y = rnorm(12), x = rep(1:4, 3), x2 = rep(c(5, 6, 7), 4))
  expect_error(sign_groups_of(y ~ mo(x):mo(x2), dat),
               "more than one mo\\(\\) call")
  expect_error(sign_groups_of(y ~ I(2 * mo(x)), dat),
               "inside other calls")
})

# Sign-constrained glm_ridge() --------------------------------------------

# Data with non-monotonic effects of two mo() predictors (columns 1-4 and 5-7
# of `x`), so that the constraints bind:
cons_data <- function(family, n = 150, seed = 3611) {
  set.seed(seed)
  x1 <- sample(rep_len(0:4, n))
  x2 <- sample(rep_len(1:4, n))
  z <- rnorm(n)
  eta <- 0.3 + c(0, 0.8, -0.4, 0.6, 0.2)[x1 + 1] +
    c(0, -0.5, -0.2, -0.9)[x2] + 0.4 * z
  y <- switch(family$family,
              gaussian = eta + rnorm(n, sd = 0.7),
              binomial = rbinom(n, 1, family$linkinv(eta)),
              poisson = rpois(n, family$linkinv(eta)))
  x <- cbind(strip_attr(mo(x1)), strip_attr(mo(x2)), z, rnorm(n))
  list(x = x, y = y, weights = sample(1:3, n, replace = TRUE),
       offset = rnorm(n, sd = 0.2))
}

# Gradient of the penalized loss of glm_ridge() on the scale of the
# standardized `x` (with the intercept first, if any):
ridge_grad <- function(fit, x, y, family, weights, offset, lambda,
                       intercept = TRUE, normalize = TRUE) {
  transf <- standardization(x, center = intercept, scale = normalize,
                            weights = weights)
  x_std <- t((t(x) - transf$shift) / transf$scale)
  b_std <- drop(fit$beta) * transf$scale
  pen <- rep(1, ncol(x))
  if (intercept) {
    x_std <- cbind(1, x_std)
    b_std <- c(fit$beta0 + sum(transf$shift * fit$beta), b_std)
    pen <- c(0, pen)
  }
  pd <- pseudo_data(drop(x_std %*% b_std), y, family, offset = offset,
                    weights = weights)
  drop(crossprod(x_std, pd$grad)) + lambda * pen * b_std
}

test_that(paste(
  "glm_ridge() with sign constraints gives the constrained optimum (KKT",
  "conditions, equality with the fit without the bound columns)"
), {
  lambda <- 1e-4
  fams <- list(extend_family(gaussian()), extend_family(binomial()),
               extend_family(poisson()))
  # Constrained columns and their signs (the first four settings: both mo()
  # predictors with all sign patterns; the last: columns not in order):
  cons_sets <- list(
    list(cols = 1:7, s = rep(c(1, 1), c(4, 3))),
    list(cols = 1:7, s = rep(c(-1, 1), c(4, 3))),
    list(cols = 1:7, s = rep(c(1, -1), c(4, 3))),
    list(cols = 1:7, s = rep(c(-1, -1), c(4, 3))),
    list(cols = c(5L, 2L, 3L, 7L), s = c(-1, 1, 1, -1))
  )
  settings <- list(c(intercept = TRUE, normalize = TRUE),
                   c(intercept = FALSE, normalize = TRUE),
                   c(intercept = TRUE, normalize = FALSE))
  n_bound_total <- 0
  for (fam in fams) {
    d <- cons_data(fam)
    for (cs in cons_sets) {
      for (st in settings) {
        cons <- cs$cols
        s <- cs$s
        fit <- glm_ridge(d$x, d$y, family = fam, lambda = lambda,
                         weights = d$weights, offset = d$offset,
                         intercept = st[["intercept"]],
                         normalize = st[["normalize"]], thresh = 1e-12,
                         cons_cols = cons, cons_signs = s)
        b <- drop(fit$beta)
        expect_true(all(s * b[cons] >= 0))
        bound <- cons[b[cons] == 0]
        n_bound_total <- n_bound_total + length(bound)
        # The KKT conditions:
        grad <- ridge_grad(fit, d$x, d$y, fam, d$weights, d$offset, lambda,
                           intercept = st[["intercept"]],
                           normalize = st[["normalize"]])
        i0 <- as.integer(st[["intercept"]])
        scale_g <- max(abs(grad), 1)
        free <- setdiff(seq_len(ncol(d$x)), bound)
        expect_lt(max(abs(grad[c(seq_len(i0), free + i0)])), 1e-6 * scale_g)
        if (length(bound) > 0) {
          expect_true(all(s[match(bound, cons)] * grad[bound + i0] >
                            -1e-6 * scale_g))
        }
        # The fit without the bound columns (without constraints) has the
        # same coefficients:
        fit_free <- glm_ridge(d$x[, free, drop = FALSE], d$y, family = fam,
                              lambda = lambda, weights = d$weights,
                              offset = d$offset,
                              intercept = st[["intercept"]],
                              normalize = st[["normalize"]], thresh = 1e-12)
        expect_equal(c(fit$beta0, b[free]),
                     c(fit_free$beta0, drop(fit_free$beta)),
                     tolerance = 1e-6)
      }
    }
  }
  # The constraints bind in most of these fits:
  expect_gt(n_bound_total, 60)
})

test_that(paste(
  "glm_ridge() with sign constraints that do not bind gives the",
  "unconstrained fit, also from an infeasible start"
), {
  lambda <- 1e-4
  for (fam in list(extend_family(gaussian()), extend_family(binomial()))) {
    d <- cons_data(fam)
    fit_unc <- glm_ridge(d$x, d$y, family = fam, lambda = lambda,
                         weights = d$weights, offset = d$offset,
                         thresh = 1e-12)
    b_unc <- drop(fit_unc$beta)
    # One constraint per column, with the sign of the unconstrained solution:
    cons <- seq_len(ncol(d$x))
    s <- ifelse(b_unc >= 0, 1, -1)
    fit <- glm_ridge(d$x, d$y, family = fam, lambda = lambda,
                     weights = d$weights, offset = d$offset, thresh = 1e-12,
                     cons_cols = cons, cons_signs = s)
    expect_equal(c(fit$beta0, drop(fit$beta)), c(fit_unc$beta0, b_unc),
                 tolerance = 1e-6)
    # An infeasible start (on the standardized scale) is moved to the bounds:
    fit_start <- glm_ridge(d$x, d$y, family = fam, lambda = lambda,
                           weights = d$weights, offset = d$offset,
                           thresh = 1e-12, cons_cols = cons, cons_signs = s,
                           beta_init = -s)
    expect_equal(c(fit_start$beta0, drop(fit_start$beta)),
                 c(fit_unc$beta0, b_unc), tolerance = 1e-6)
  }
})

test_that("glm_ridge() checks the sign constraints", {
  d <- cons_data(extend_family(gaussian()), n = 20)
  expect_error(glm_ridge(d$x, d$y, cons_cols = 1:2, cons_signs = 1),
               "invalid sign constraints")
  expect_error(glm_ridge(d$x, d$y, cons_cols = 10L, cons_signs = 1),
               "invalid sign constraints")
  expect_error(glm_ridge(d$x, d$y, cons_cols = 1L, cons_signs = 0),
               "invalid sign constraints")
  expect_error(glm_ridge(d$x, d$y, cons_cols = c(2L, 2L), cons_signs = c(1, 1)),
               "invalid sign constraints")
  expect_error(glm_ridge(NULL, d$y, cons_cols = 1L, cons_signs = 1),
               "invalid sign constraints")
})

# Traditional projection --------------------------------------------------

test_that(paste(
  "the projection onto the full submodel reproduces the monotonic effects of",
  "the reference model draw by draw"
), {
  sim <- mo_sim(x_min = -2L)
  trms <- c("mo(x)", "mo(f)", "z", "w")
  for (fam in list(gaussian(), binomial(), poisson())) {
    y_nm <- switch(fam$family, gaussian = "y_gauss", binomial = "y_binom",
                   poisson = "y_poiss")
    fml <- as.formula(paste(y_nm, "~ mo(x) + mo(f) + z + w"))
    refm <- mo_refm(sim, fml, family = fam)
    prj <- project(refm, predictor_terms = trms, ndraws = sim$S,
                   verbose = FALSE)
    coefs <- do.call(rbind, lapply(prj$outdmin, function(o) {
      drop(coef(o))
    }))
    # The ridge penalty of the submodel fitter (`regul`) has a larger effect
    # for the non-Gaussian families:
    tol <- if (fam$family == "gaussian") 1e-4 else 1e-3
    expect_equal(strip_attr(coefs[, colnames(sim$coef_true)]),
                 strip_attr(sim$coef_true), tolerance = tol,
                 info = fam$family)
    expect_equal(unname(coefs[, "w"]), rep(0, sim$S), tolerance = tol,
                 info = fam$family)
    expect_equal(strip_attr(proj_linpred(prj, transform = FALSE)$pred),
                 t(sim$eta_fun(sim$dat)), tolerance = tol, info = fam$family)
  }
})

test_that(paste(
  "projected predictions for new data with mo() terms equal the predictions",
  "for the same training observations"
), {
  sim <- mo_sim()
  nd <- mo_newdata(sim)
  refm <- mo_refm(sim, y_gauss ~ mo(x) + mo(f) + z + w)
  prj <- project(refm, predictor_terms = c("mo(x)", "mo(f)", "z"),
                 ndraws = 25, seed = 3906, verbose = FALSE)
  pred_train <- proj_linpred(prj)$pred
  pred_new <- proj_linpred(prj, newdata = nd)$pred
  expect_equal(unname(pred_new), unname(pred_train[, rows_of(nd)]),
               tolerance = 1e-10)
  # New data with a value of `x` above the training levels:
  nd_out <- nd[1, ]
  nd_out$x <- sim$D + 1L
  expect_error(proj_linpred(prj, newdata = nd_out), "outside the range")
})

test_that("forward search ranks mo() terms as single terms", {
  sim <- mo_sim()
  refm <- mo_refm(sim, y_gauss ~ mo(x) + mo(f) + z + w)
  vs <- varsel(refm, nclusters = 4, nclusters_pred = 8, seed = 3906,
               verbose = FALSE)
  rk <- ranking(vs)$fulldata
  expect_setequal(rk, c("mo(x)", "mo(f)", "z", "w"))
  expect_identical(utils::tail(rk, 1), "w")
})

test_that(paste(
  "L1 search ranks mo() terms as single terms, also with more than 10",
  "levels, and its submodels predict new data"
), {
  sim <- mo_sim(D = 11)
  nd <- sim$dat[sim$dat$x <= 5, , drop = FALSE][1:8, ]
  refm <- mo_refm(sim, y_gauss ~ mo(x) + z + w)
  # Monotonic projections do not support L1 search:
  expect_error(varsel(refm, method = "L1", verbose = FALSE),
               "do not support L1 search")
  opts_old <- options(projpred.mo_monotone = FALSE)
  on.exit(options(opts_old), add = TRUE)
  vs <- varsel(refm, method = "L1", refit_prj = FALSE, seed = 3906,
               verbose = FALSE)
  rk <- ranking(vs)$fulldata
  expect_setequal(rk, c("mo(x)", "z", "w"))
  prj <- project(vs, nterms = 3, refit_prj = FALSE, verbose = FALSE)
  pred_train <- proj_linpred(prj)$pred
  pred_new <- proj_linpred(prj, newdata = nd)$pred
  expect_equal(unname(drop(pred_new)), unname(drop(pred_train)[rows_of(nd)]),
               tolerance = 1e-10)
})

test_that("interactions of mo() terms work in forward and L1 search", {
  sim <- mo_sim()
  refm <- mo_refm(sim, y_gauss ~ mo(x) + z + mo(x):z)
  vs_fw <- varsel(refm, nclusters = 4, nclusters_pred = 8, seed = 3906,
                  verbose = FALSE)
  rk_fw <- ranking(vs_fw)$fulldata
  expect_setequal(rk_fw, c("mo(x)", "z", "mo(x):z"))
  expect_identical(rk_fw[3], "mo(x):z")
  # (L1 search only with unconstrained projections.)
  opts_old <- options(projpred.mo_monotone = FALSE)
  on.exit(options(opts_old), add = TRUE)
  vs_l1 <- varsel(refm, method = "L1", seed = 3906, verbose = FALSE)
  expect_setequal(ranking(vs_l1)$fulldata, c("mo(x)", "z", "mo(x):z"))
})

test_that("multilevel and additive submodels with mo() terms work", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("mgcv")
  sim <- mo_sim()
  nd <- mo_newdata(sim)
  # The binomial multilevel fits to non-integer pseudo-data throw lme4
  # convergence warnings; this test only compares the predictions (of
  # unconstrained projections, because monotonic projections do not support
  # these submodels):
  opts_old <- options(projpred.warn_additive_experimental = FALSE,
                      projpred.warn_instable_projections = FALSE,
                      projpred.warn_proj_drawwise = FALSE,
                      projpred.check_convergence = FALSE,
                      projpred.mo_monotone = FALSE)
  on.exit(options(opts_old), add = TRUE)
  cases <- list(
    glmm = list(fml = y_gauss ~ mo(x) + z + (1 | g),
                trms = c("mo(x)", "z", "(1 | g)"), family = gaussian()),
    glmm_binom = list(fml = y_binom ~ mo(x) + z + (1 | g),
                      trms = c("mo(x)", "z", "(1 | g)"), family = binomial()),
    gam = list(fml = y_gauss ~ mo(x) + s(z), trms = c("mo(x)", "s(z)"),
               family = gaussian())
  )
  for (nm in names(cases)) {
    refm <- mo_refm(sim, cases[[nm]]$fml, family = cases[[nm]]$family)
    prj <- muffle_warning(muffle_warning(
      project(refm, predictor_terms = cases[[nm]]$trms, ndraws = 25,
              seed = 3906, verbose = FALSE),
      "Model failed to converge"
    ), "Model is nearly unidentifiable")
    pred_train <- proj_linpred(prj)$pred
    pred_new <- proj_linpred(prj, newdata = nd)$pred
    expect_equal(unname(pred_new), unname(pred_train[, rows_of(nd)]),
                 tolerance = 1e-8, info = nm)
  }
})

# Monotonic projections ---------------------------------------------------

# A reference model with an effect `eff_x` of the integer `x` (levels 0, ...,
# 4; by default non-monotonic), an effect of `z`, an effect `c_w` of `w` (which
# decreases with `x`), and no effect of the integer `x2` (levels 1, 2, 3), of
# the factor `f` and of the grouping variable `g`. Each draw adds `dip`
# (normal with mean `dip_mean` and SD `dip_sd`) to the effect of the levels 2,
# 3 and 4 of `x`.
mo_nonmono_refm <- function(family = gaussian(), fml = y ~ mo(x) + z + mo(x2),
                            S = 20, n = 100, seed = 5129,
                            eff_x = c(0, 0.8, -0.4, 0.6, 0.2), c_w = 0,
                            dip_mean = 0, dip_sd = 0, sd_g = 0, c_h = 0,
                            ...) {
  set.seed(seed)
  dat <- data.frame(x = sample(rep_len(0:4, n)), z = rnorm(n),
                    x2 = sample(rep_len(1:3, n)),
                    f = factor(sample(rep_len(c("a", "b"), n))),
                    g = factor(sample(rep_len(paste0("gr", 1:5), n))))
  dat$w <- -0.5 * dat$x + rnorm(n, sd = 0.5)
  b <- cbind(alpha = rnorm(S, 0.1, 0.05), bx = rnorm(S, 1, 0.1),
             bz = rnorm(S, 0.4, 0.05), dip = rnorm(S, dip_mean, dip_sd))
  # Group effects (the same in all draws), only if `sd_g > 0`:
  eff_g <- if (sd_g > 0) rnorm(nlevels(dat$g), sd = sd_g) else 0
  # A covariate `h` which is not in the formulas, so that no submodel fits the
  # draws exactly (only if `c_h > 0`):
  dat$h <- if (c_h > 0) rnorm(n) else 0
  # Linear predictor (N x S) for data `nd`:
  eta_fun <- function(nd) {
    outer(eff_x[nd$x + 1], b[, "bx"]) + outer(nd$z, b[, "bz"]) +
      outer(nd$x >= 2, b[, "dip"]) + c_w * nd$w + c_h * nd$h +
      matrix(b[, "alpha"], nrow(nd), S, byrow = TRUE) +
      if (sd_g > 0) eff_g[as.integer(nd$g)] else 0
  }
  eta_mean <- rowMeans(eta_fun(dat))
  dat$y <- switch(family$family,
                  gaussian = eta_mean + rnorm(n, sd = 0.5),
                  binomial = rbinom(n, 1, plogis(eta_mean)))
  environment(fml) <- list2env(list(mo = mo_brms_like),
                               parent = environment(fml))
  init_refmodel(
    list(data = dat), data = dat, formula = fml, family = family,
    ref_predfun = function(fit, newdata = NULL) {
      eta_fun(if (is.null(newdata)) fit$data else newdata)
    },
    dis = if (family$family == "gaussian") rep(0.5, S) else NULL,
    cvfun = function(folds) NULL, cvrefbuilder = function(cvfit) NULL, ...
  )
}

# The coefficients (S x columns, without intercept) of the projected draws:
prj_coefs <- function(prj) {
  out <- t(vapply(prj$outdmin, function(o) as.vector(o$beta),
                  numeric(nrow(prj$outdmin[[1]]$beta))))
  colnames(out) <- rownames(prj$outdmin[[1]]$beta)
  out
}

# Are the rows of `b` monotonic (all elements >= 0 or all <= 0, with the
# tolerance of the constraint check)?
is_mono <- function(b, tol = 1e-8) {
  apply(b, 1, function(v) all(v >= -tol) || all(v <= tol))
}

# The number of calls of the function `fun` of projpred while `expr` runs
# (with `non_null = TRUE`: of the calls that return a value other than
# `NULL`).
count_calls <- function(expr, fun, non_null = FALSE) {
  n <- 0
  fun_orig <- get(fun, envir = asNamespace("projpred"))
  wrapper <- function(...) {
    out <- fun_orig(...)
    if (!non_null || !is.null(out)) {
      n <<- n + 1
    }
    out
  }
  out <- do.call(with_mocked_bindings,
                 c(list(code = substitute(expr)), stats::setNames(list(wrapper),
                                                                  fun)),
                 envir = parent.frame())
  list(value = out, n = n)
}

test_that(paste(
  "monotonic projections give monotonic mo() effects with the signs from",
  "the projection onto the full model"
), {
  cols_x <- paste0("mo(x)", 1:4)
  opts_old <- options(projpred.mo_monotone = TRUE)
  on.exit(options(opts_old), add = TRUE)
  for (fam in list(gaussian(), binomial())) {
    # `w` decreases with `x`, so the effect of `x` in the submodel without `w`
    # has the opposite sign of its effect in the full model:
    refm <- mo_nonmono_refm(fam, fml = y ~ mo(x) + z + w, c_w = 1)
    S <- ncol(refm$mu)
    # (With a small convergence threshold for the comparison with glm_ridge()
    # below.)
    options(projpred.mo_monotone = TRUE)
    prj <- project(refm, predictor_terms = c("mo(x)", "z"), ndraws = S,
                   thresh_conv = 1e-12, verbose = FALSE)
    b_x <- prj_coefs(prj)[, cols_x]
    expect_true(all(is_mono(b_x)), info = fam$family)
    # Unconstrained projections (onto the submodel and onto the full model):
    options(projpred.mo_monotone = FALSE)
    prj_unc <- project(refm, predictor_terms = c("mo(x)", "z"), ndraws = S,
                       thresh_conv = 1e-12, verbose = FALSE)
    prj_full <- project(refm, predictor_terms = c("mo(x)", "z", "w"),
                        ndraws = S, verbose = FALSE)
    b_unc <- prj_coefs(prj_unc)
    binding <- !is_mono(b_unc[, cols_x])
    expect_gt(sum(binding), S / 2)
    s_full <- ifelse(rowSums(prj_coefs(prj_full)[, cols_x]) >= 0, 1, -1)
    s_sum <- ifelse(rowSums(b_unc[, cols_x]) >= 0, 1, -1)
    # The signs of the full model differ from those of the submodel:
    expect_true(all(s_full[binding] != s_sum[binding]), info = fam$family)
    expect_true(all(b_x[binding, ] * s_full[binding] >= 0), info = fam$family)
    # Draws without a binding constraint keep the unconstrained fit:
    expect_equal(prj_coefs(prj)[!binding, , drop = FALSE],
                 b_unc[!binding, , drop = FALSE], tolerance = 1e-12,
                 info = fam$family)
    # The constrained draws equal sign-constrained glm_ridge() fits to the
    # reference model's draws (`ndraws = S` projects all draws). This checks
    # the signs, columns and arguments that the projection passes, not the
    # bounded solver (see the KKT tests above):
    dat <- refm$fetch_data()
    x <- cbind(strip_attr(mo(dat$x)), z = dat$z)
    for (s in which(binding)[1:3]) {
      fit_s <- glm_ridge(x, refm$mu[, s], family = refm$family,
                         weights = refm$wobs, lambda = 1e-4, thresh = 1e-12,
                         cons_cols = 1:4, cons_signs = rep(s_full[s], 4))
      expect_equal(unname(c(prj$outdmin[[s]]$alpha, prj_coefs(prj)[s, ])),
                   unname(c(fit_s$beta0, drop(fit_s$beta))),
                   tolerance = 1e-6, info = fam$family)
    }
  }
})

test_that(paste(
  "mo_cons_signs() takes the signs of the full model and otherwise the sign",
  "of the sum of the coefficients"
), {
  groups <- structure(list(1:2, 3:4, 5:6), keys = c("a#", "b#", "c#1"))
  beta <- c(1, -2, 3, -1, -0.5, 0.2)
  expect_identical(mo_cons_signs(beta, groups, NULL), c(-1, 1, -1))
  expect_identical(mo_cons_signs(beta, groups, c("a#" = 1, "c#" = 1)),
                   c(1, 1, -1))
})

test_that(paste(
  "a failure of the projection onto the full model gives a warning and no",
  "full-model signs (once per set of projected draws)"
), {
  refm <- mo_nonmono_refm()
  refm$div_minimizer <- function(...) stop("full fit failed")
  p_ref <- get_refdist(refm, ndraws = 4)
  expect_warning(
    signs <- mo_full_signs(p_ref, refm, args_divmin = list()),
    "failed with the error \"full fit failed\""
  )
  expect_identical(dim(signs), c(4L, 0L))
  # The second call uses the cache:
  expect_silent(signs2 <- mo_full_signs(p_ref, refm, args_divmin = list()))
  expect_identical(signs2, signs)
})

test_that(paste(
  "monotonic projections equal the unconstrained ones if no constraint binds"
), {
  sim <- mo_sim()
  refm <- mo_refm(sim, y_gauss ~ mo(x) + mo(f) + z + w)
  prj <- project(refm, predictor_terms = c("mo(x)", "mo(f)", "z"),
                 ndraws = 10, verbose = FALSE)
  opts_old <- options(projpred.mo_monotone = FALSE)
  on.exit(options(opts_old), add = TRUE)
  prj_unc <- project(refm, predictor_terms = c("mo(x)", "mo(f)", "z"),
                     ndraws = 10, verbose = FALSE)
  expect_identical(prj_coefs(prj), prj_coefs(prj_unc))
})

test_that(paste(
  "the joint Gaussian solve with refits of the binding draws gives the",
  "draw-wise monotonic projections"
), {
  # About half of the draws have a non-monotonic effect of `x`:
  refm <- mo_nonmono_refm(eff_x = c(0, 0.4, 0.8, 1.2, 1.6), dip_mean = -0.4,
                          dip_sd = 0.4)
  prj_trms <- c("mo(x)", "z")
  cnt_joint <- count_calls(
    project(refm, predictor_terms = prj_trms, ndraws = 20, verbose = FALSE),
    "glm_ridge_gauss_multi", non_null = TRUE
  )
  # The joint solve ran (for the submodel and for the full model):
  expect_identical(cnt_joint$n, 2)
  cnt_refit <- count_calls(
    project(refm, predictor_terms = prj_trms, ndraws = 20, verbose = FALSE),
    "mo_cons_args"
  )
  # Only the binding draws were fitted again:
  expect_gt(cnt_refit$n, 0)
  expect_lt(cnt_refit$n, 20)
  prj_dw <- with_mocked_bindings(
    project(refm, predictor_terms = prj_trms, ndraws = 20, verbose = FALSE),
    glm_ridge_gauss_multi = function(...) NULL
  )
  expect_equal(prj_coefs(cnt_joint$value), prj_coefs(prj_dw),
               tolerance = 1e-10)
})

test_that(paste(
  "the parallel case of divmin() gives the same monotonic projections"
), {
  skip_if_not_installed("foreach")
  skip_if_not_installed("iterators")
  # Use a registered backend, otherwise the sequential one:
  if (!foreach::getDoParRegistered()) {
    foreach::registerDoSEQ()
  }
  opts_old <- options(projpred.parallel_proj_trigger = Inf)
  on.exit(options(opts_old), add = TRUE)
  for (fam in list(gaussian(), binomial())) {
    refm <- mo_nonmono_refm(fam, eff_x = c(0, 0.4, 0.8, 1.2, 1.6),
                            dip_mean = -0.4, dip_sd = 0.4)
    options(projpred.parallel_proj_trigger = Inf)
    prj_seq <- project(refm, predictor_terms = c("mo(x)", "z"), ndraws = 20,
                       verbose = FALSE)
    options(projpred.parallel_proj_trigger = 1)
    cnt <- new.env()
    cnt$n <- 0
    trace(foreach::foreach,
          bquote(assign("n", get("n", envir = .(cnt)) + 1, envir = .(cnt))),
          print = FALSE, where = asNamespace("foreach"))
    prj_par <- tryCatch(
      project(refm, predictor_terms = c("mo(x)", "z"), ndraws = 20,
              verbose = FALSE),
      finally = untrace(foreach::foreach, where = asNamespace("foreach"))
    )
    expect_gt(cnt$n, 0)
    expect_equal(prj_coefs(prj_par), prj_coefs(prj_seq), tolerance = 1e-12,
                 info = fam$family)
  }
})

test_that(paste(
  "the signs from the full model are computed once per set of projected",
  "draws"
), {
  refm <- mo_nonmono_refm()
  cnt <- count_calls(
    varsel(refm, nclusters = 4, nclusters_pred = 6, seed = 3906,
           verbose = FALSE),
    "mo_fit_coefs"
  )
  # One full-model projection for the search (4 clusters) and one for the
  # performance evaluation (6 clusters):
  expect_identical(cnt$n, 4 + 6)
  expect_setequal(ranking(cnt$value)$fulldata, c("mo(x)", "z", "mo(x2)"))
})

test_that(paste(
  "monotonic projections work in LOO CV with validated search, in the",
  "latent projection and with a multilevel reference model"
), {
  opts_old <- options(projpred.warn_psis = FALSE)
  on.exit(options(opts_old), add = TRUE)
  cols_x <- paste0("mo(x)", 1:4)
  refm <- mo_nonmono_refm()
  cv_loo <- muffle_warning(muffle_warning(cv_varsel(
    refm, cv_method = "LOO", validate_search = TRUE, nloo = 3, nclusters = 3,
    nclusters_pred = 6, seed = 3906, verbose = FALSE
  ), "Some Pareto k diagnostic values are too high"), "Not enough tail samples")
  expect_true(all(is.finite(summary(cv_loo)$perf_sub$elpd)))
  prj <- project(cv_loo, nterms = 3, ndraws = 20, verbose = FALSE)
  expect_true(all(is_mono(prj_coefs(prj)[, cols_x])))
  # Latent projection:
  refm_lat <- suppressMessages(mo_nonmono_refm(binomial(), latent = TRUE))
  prj_lat <- project(refm_lat, predictor_terms = c("mo(x)", "z"), ndraws = 20,
                     verbose = FALSE)
  b_lat <- prj_coefs(prj_lat)[, cols_x]
  expect_true(all(is_mono(b_lat)))
  expect_gt(sum(rowSums(abs(b_lat)) > 0), 0)
  # Multilevel reference model (the full model is projected with lme4), GLM
  # submodel. `g` has no effect, so lme4 warns about the full-model fits:
  skip_if_not_installed("lme4")
  refm_ml <- mo_nonmono_refm(fml = y ~ mo(x) + z + (1 | g), sd_g = 0.7,
                            c_h = 0.5)
  prj_ml <- muffle_warning(muffle_warning(muffle_warning(
    project(refm_ml, predictor_terms = c("mo(x)", "z"), ndraws = 20,
            verbose = FALSE),
    "Model failed to converge"
  ), "Model is nearly unidentifiable"), "might not have converged")
  expect_true(all(is_mono(prj_coefs(prj_ml)[, cols_x])))
})

# KKT conditions of the generalized least-squares step of the constrained
# multilevel fit `fit` (an lmerMod object) with the constraints
# signs[k] * beta[cols[k]] >= 0: the largest gradient of the free coefficients
# (relative to the largest gradient) and the smallest signed gradient of the
# coefficients at the bound 0 (must be nonnegative).
lmm_kkt <- function(fit, cols, signs) {
  X <- lme4::getME(fit, "X")
  sw <- sqrt(stats::weights(fit))
  LtZt <- lme4::getME(fit, "Lambdat") %*% lme4::getME(fit, "Zt") %*%
    Matrix::Diagonal(x = sw)
  L <- Matrix::Cholesky(Matrix::tcrossprod(LtZt), LDL = FALSE, Imult = 1)
  vinv <- function(a) {
    as.matrix(a) - as.matrix(Matrix::crossprod(
      LtZt, Matrix::solve(L, LtZt %*% as.matrix(a), system = "A")
    ))
  }
  Xw <- X * sw
  beta <- lme4::fixef(fit)
  g <- drop(crossprod(Xw, vinv(Xw)) %*% beta -
              crossprod(Xw, vinv(lme4::getME(fit, "y") * sw)))
  bound <- seq_along(beta) %in% cols & beta == 0
  list(free = max(abs(g[!bound])) / max(abs(g)),
       bound = min(c(Inf, (g * signs[match(seq_along(beta), cols)])[bound])))
}

test_that(paste(
  "monotonic projections onto Gaussian multilevel submodels satisfy the KKT",
  "conditions and predict new data"
), {
  skip_if_not_installed("lme4")
  # (The submodels include the full model, for which projpred warns.)
  opts_old <- options(projpred.warn_instable_projections = FALSE)
  on.exit(options(opts_old), add = TRUE)
  cols_x <- paste0("mo(x)", 1:4)
  trms <- c("mo(x)", "z", "(1 | g)")
  refm <- mo_nonmono_refm(fml = y ~ mo(x) + z + (1 | g), sd_g = 0.7,
                          c_h = 0.5, eff_x = c(0, 0.4, 0.8, 1.2, 1.6),
                          dip_mean = -0.4,
                          dip_sd = 0.4)
  prj <- project(refm, predictor_terms = trms, ndraws = 10, verbose = FALSE)
  fixefs <- t(vapply(prj$outdmin, lme4::fixef, numeric(6)))
  expect_true(all(is_mono(fixefs[, cols_x])))
  opts_old2 <- options(projpred.mo_monotone = FALSE)
  on.exit(options(opts_old2), add = TRUE)
  prj_unc <- project(refm, predictor_terms = trms, ndraws = 10,
                     verbose = FALSE)
  fixefs_unc <- t(vapply(prj_unc$outdmin, lme4::fixef, numeric(6)))
  binding <- !is_mono(fixefs_unc[, cols_x])
  expect_gt(sum(binding), 0)
  expect_lt(sum(binding), 10)
  # Draws without a binding constraint keep the lmer() fit:
  expect_equal(fixefs[!binding, ], fixefs_unc[!binding, ], tolerance = 1e-12)
  # The submodel is the full model, so the signs are those of the
  # unconstrained fits:
  s_full <- ifelse(rowSums(fixefs_unc[, cols_x]) >= 0, 1, -1)
  for (s in which(binding)) {
    fit <- prj$outdmin[[s]]
    kkt <- lmm_kkt(fit, cols = 2:5, signs = rep(s_full[s], 4))
    expect_lt(kkt$free, 1e-8)
    expect_gt(kkt$bound, -1e-8)
    # The lmerMod object: fitted values, conditional modes and sigma for the
    # constrained coefficients:
    X <- lme4::getME(fit, "X")
    Zt <- as.matrix(lme4::getME(fit, "Zt"))
    Lt <- as.matrix(lme4::getME(fit, "Lambdat"))
    w <- stats::weights(fit)
    y <- lme4::getME(fit, "y")
    beta <- lme4::fixef(fit)
    u <- lme4::getME(fit, "u")
    expect_equal(unname(fitted(fit)),
                 unname(drop(X %*% beta +
                               t(Zt) %*% as.matrix(lme4::getME(fit, "b")))),
                 tolerance = 1e-10)
    lhs <- Lt %*% Zt %*% (w * t(Zt)) %*% t(Lt) + diag(nrow(Lt))
    rhs <- Lt %*% Zt %*% (w * (y - drop(X %*% beta)))
    expect_equal(drop(lhs %*% u), drop(rhs), tolerance = 1e-8)
    pwrss <- sum(w * (y - fitted(fit))^2) + sum(u^2)
    expect_equal(sigma(fit)^2, pwrss / (nrow(X) - ncol(X)), tolerance = 1e-10)
  }
  # Predictions for new data equal the predictions for the same training
  # observations:
  nd <- refm$fetch_data()[c(3, 8, 20, 41), ]
  pred_train <- proj_linpred(prj)$pred
  pred_new <- proj_linpred(prj, newdata = nd)$pred
  expect_equal(unname(pred_new), unname(pred_train[, rows_of(nd)]),
               tolerance = 1e-10)
})

test_that(paste(
  "the constrained ML fit of a multilevel submodel equals the lmer() ML fit",
  "without the columns at the bound (with observation weights)"
), {
  skip_if_not_installed("lme4")
  refm <- mo_nonmono_refm(fml = y ~ mo(x) + z + (z | g), sd_g = 0.7,
                          c_h = 0.5)
  dat <- refm$fetch_data()
  dat$y1 <- refm$mu[, 1]
  set.seed(4117)
  wts <- runif(nrow(dat), 0.5, 2)
  for (re in c("(1 | g)", "(z | g)")) {
    fml <- mo_formula_env(as.formula(paste("y1 ~ mo(x) + z +", re)))
    args <- list(data = dat, weights = wts, REML = FALSE,
                 control = lme4::lmerControl())
    fit_u <- do.call(lme4::lmer, c(list(formula = fml), args))
    b_u <- lme4::fixef(fit_u)[2:5]
    expect_false(all(b_u >= 0) || all(b_u <= 0), info = re)
    # Positive signs (the reference model's effect of `x`):
    fit_c <- mo_lmer_cons(fit_u, fml, args, c("mo(x)#" = 1))
    beta <- lme4::fixef(fit_c)
    bound <- which(seq_along(beta) %in% 2:5 & beta == 0)
    expect_gt(length(bound), 0)
    expect_true(all(beta[2:5] >= 0))
    # The ML fit without the columns at the bound:
    dat$xr <- I(lme4::getME(fit_u, "X")[, -bound, drop = FALSE])
    fit_r <- lme4::lmer(as.formula(paste("y1 ~ 0 + xr +", re)), data = dat,
                        weights = wts, REML = FALSE)
    expect_equal(unname(beta[-bound]), unname(lme4::fixef(fit_r)),
                 tolerance = 1e-5, info = re)
    expect_equal(unname(lme4::getME(fit_c, "theta")),
                 unname(lme4::getME(fit_r, "theta")), tolerance = 1e-4,
                 info = re)
    expect_equal(deviance(fit_c), deviance(fit_r), tolerance = 1e-8,
                 info = re)
  }
})

test_that(paste(
  "without signs from the full model, monotonic projections use the signs",
  "of the sums of the unconstrained coefficients (GLM and LMM submodels)"
), {
  skip_if_not_installed("lme4")
  opts_old <- options(projpred.warn_instable_projections = FALSE)
  on.exit(options(opts_old), add = TRUE)
  cols_x <- paste0("mo(x)", 1:4)
  # `w` decreases with `x`, so the signs from the full model and the signs of
  # the sums differ:
  cases <- list(
    glm = list(fml = y ~ mo(x) + z + w, trms = c("mo(x)", "z")),
    lmm = list(fml = y ~ mo(x) + z + w + (1 | g),
               trms = c("mo(x)", "z", "(1 | g)"))
  )
  for (nm in names(cases)) {
    refm <- mo_nonmono_refm(fml = cases[[nm]]$fml, c_w = 1, sd_g = 0.7,
                            c_h = 0.5)
    prj <- with_mocked_bindings(
      project(refm, predictor_terms = cases[[nm]]$trms, ndraws = 10,
              verbose = FALSE),
      mo_full_signs = function(p_ref, ...) {
        matrix(numeric(), nrow = NCOL(p_ref$mu), ncol = 0)
      }
    )
    opts_old2 <- options(projpred.mo_monotone = FALSE)
    prj_unc <- project(refm, predictor_terms = cases[[nm]]$trms, ndraws = 10,
                       verbose = FALSE)
    options(opts_old2)
    coefs_of <- function(p) {
      t(vapply(p$outdmin, function(f) mo_fit_coefs(f)[cols_x], numeric(4)))
    }
    b_x <- coefs_of(prj)
    b_unc <- coefs_of(prj_unc)
    binding <- !is_mono(b_unc)
    expect_gt(sum(binding), 0)
    s_sum <- ifelse(rowSums(b_unc) >= 0, 1, -1)
    expect_true(all(b_x[binding, ] * s_sum[binding] >= 0), info = nm)
    expect_true(all(is_mono(b_x)), info = nm)
  }
})

test_that(paste(
  "monotonic projections onto multilevel submodels work in the latent",
  "projection and throw an error for non-Gaussian families"
), {
  skip_if_not_installed("lme4")
  # (The submodels include the full model, for which projpred warns.)
  opts_old <- options(projpred.warn_instable_projections = FALSE)
  on.exit(options(opts_old), add = TRUE)
  cols_x <- paste0("mo(x)", 1:4)
  trms <- c("mo(x)", "z", "(1 | g)")
  refm_lat <- suppressMessages(mo_nonmono_refm(
    binomial(), fml = y ~ mo(x) + z + (1 | g), sd_g = 0.7, c_h = 0.5,
    latent = TRUE
  ))
  prj_lat <- project(refm_lat, predictor_terms = trms, ndraws = 10,
                     verbose = FALSE)
  fixefs <- t(vapply(prj_lat$outdmin, lme4::fixef, numeric(6)))
  expect_true(all(is_mono(fixefs[, cols_x])))
  opts_old2 <- options(projpred.mo_monotone = FALSE)
  prj_lat_unc <- project(refm_lat, predictor_terms = trms, ndraws = 10,
                         verbose = FALSE)
  options(opts_old2)
  fixefs_unc <- t(vapply(prj_lat_unc$outdmin, lme4::fixef, numeric(6)))
  expect_gt(sum(!is_mono(fixefs_unc[, cols_x])), 0)
  refm_bin <- mo_nonmono_refm(binomial(), fml = y ~ mo(x) + z + (1 | g),
                              sd_g = 0.7, c_h = 0.5)
  expect_error(
    project(refm_bin, predictor_terms = trms, ndraws = 4, verbose = FALSE),
    "do not support multilevel submodels with a non-Gaussian family"
  )
})

test_that(paste(
  "monotonic projections throw errors for unsupported cases before the",
  "search or the projection starts"
), {
  opts_old <- options(projpred.warn_additive_experimental = FALSE,
                      projpred.warn_instable_projections = FALSE,
                      projpred.mo_monotone = TRUE)
  on.exit(options(opts_old), add = TRUE)
  # get_refdist() runs at the start of each search and projection:
  expect_error_early <- function(expr, msg) {
    with_mocked_bindings(
      expect_error(expr, msg),
      get_refdist = function(...) stop("The search started.")
    )
  }
  refm_s <- mo_nonmono_refm(fml = y ~ mo(x) + s(z))
  expect_error_early(varsel(refm_s, verbose = FALSE),
                     "do not support reference models with smooth terms")
  expect_error_early(project(refm_s, predictor_terms = "mo(x)",
                             verbose = FALSE),
                     "do not support reference models with smooth terms")
  refm <- mo_nonmono_refm()
  expect_error_early(varsel(refm, method = "L1", verbose = FALSE),
                     "do not support L1 search")
  refm_c <- refm
  refm_c$div_minimizer <- function(...) divmin(...)
  expect_error_early(varsel(refm_c, verbose = FALSE),
                     "do not support custom divergence minimizers")
  opts_old2 <- options(projpred.glm_fitter = "fit_glm_callback")
  on.exit(options(opts_old2), add = TRUE)
  expect_error_early(varsel(refm, verbose = FALSE),
                     "do not support other submodel fitters")
  # Multilevel submodels in project() do not use this fitter:
  expect_silent(check_mo_monotone(refm, c("mo(x)", "(1 | g)")))
  expect_error(check_mo_monotone(refm, c("mo(x)", "(1 | g)"),
                                 method = "forward"),
               "do not support other submodel fitters")
  options(opts_old2)
  expect_error(check_mo_monotone(refm, c("z", "I(2 * mo(x))")),
               "do not support mo\\(\\) calls inside other calls")
  # Terms that contain `+` (also without mo() terms):
  expect_silent(check_mo_monotone(refm, c("(1 + z | g)", "I(x2 + z)")))
  expect_silent(check_mo_monotone(refm, c("mo(x) + I(x2 + z)", "1")))
  refm_ia <- mo_nonmono_refm(fml = y ~ mo(x) + mo(x2) + mo(x):mo(x2))
  expect_error_early(varsel(refm_ia, verbose = FALSE),
                     "do not support terms with more than one mo\\(\\) call")
  refm_bin <- mo_nonmono_refm(binomial(), fml = y ~ mo(x) + z + (1 | g),
                              sd_g = 0.7, c_h = 0.5)
  expect_error_early(
    varsel(refm_bin, verbose = FALSE),
    "do not support multilevel submodels with a non-Gaussian family"
  )
  # Without such terms in `search_terms`, the search runs:
  vs_ia <- varsel(refm_ia, search_terms = c("mo(x)", "mo(x2)"),
                  nclusters = 2, nclusters_pred = 3, seed = 3906,
                  verbose = FALSE)
  expect_setequal(ranking(vs_ia)$fulldata, c("mo(x)", "mo(x2)"))
  # A reference model without mo() terms, with a term that contains `+`:
  refm_i <- mo_nonmono_refm(fml = y ~ I(x2 + z) + z)
  vs_i <- varsel(refm_i, nclusters = 2, nclusters_pred = 3, seed = 3906,
                 verbose = FALSE)
  expect_setequal(ranking(vs_i)$fulldata, c("I(x2 + z)", "z"))
  # With unconstrained projections, the checks do not apply:
  options(projpred.mo_monotone = FALSE)
  expect_s3_class(project(refm_s, predictor_terms = "mo(x)", ndraws = 4,
                          verbose = FALSE), "projection")
})

test_that(paste(
  "the checks in proj_to_submodl() stop unsupported monotonic projections",
  "also without the early checks"
), {
  opts_old <- options(projpred.warn_instable_projections = FALSE)
  on.exit(options(opts_old), add = TRUE)
  refm_bin <- mo_nonmono_refm(binomial(), fml = y ~ mo(x) + z + (1 | g),
                              sd_g = 0.7, c_h = 0.5)
  expect_error(
    with_mocked_bindings(
      project(refm_bin, predictor_terms = c("mo(x)", "z", "(1 | g)"),
              ndraws = 4, verbose = FALSE),
      check_mo_monotone = function(...) invisible(TRUE)
    ),
    "multilevel submodels with a non-Gaussian family"
  )
  refm_c <- mo_nonmono_refm()
  refm_c$div_minimizer <- function(...) divmin(...)
  expect_error(
    with_mocked_bindings(
      project(refm_c, predictor_terms = c("mo(x)", "z"), ndraws = 4,
              verbose = FALSE),
      check_mo_monotone = function(...) invisible(TRUE)
    ),
    "custom divergence minimizers"
  )
})

test_that(paste(
  "monotonic projections throw an error for terms with two mo() calls, but",
  "not for submodels without such terms"
), {
  refm <- mo_nonmono_refm(fml = y ~ mo(x) + mo(x2) + mo(x):mo(x2))
  trms_ia <- c("mo(x)", "mo(x2)", "mo(x):mo(x2)")
  expect_error(project(refm, predictor_terms = trms_ia, ndraws = 4,
                       verbose = FALSE),
               "more than one mo\\(\\) call")
  prj <- project(refm, predictor_terms = c("mo(x)", "mo(x2)"), ndraws = 4,
                 verbose = FALSE)
  expect_true(all(is_mono(prj_coefs(prj)[, paste0("mo(x)", 1:4)])))
  opts_old <- options(projpred.mo_monotone = FALSE)
  on.exit(options(opts_old), add = TRUE)
  expect_s3_class(project(refm, predictor_terms = trms_ia, ndraws = 4,
                          verbose = FALSE), "projection")
})

# Latent and augmented-data projection ------------------------------------

test_that("the latent projection reproduces the monotonic effects", {
  sim <- mo_sim()
  refm <- suppressMessages(mo_refm(sim, y_binom ~ mo(x) + mo(f) + z,
                                   family = binomial(), latent = TRUE))
  prj <- project(refm, predictor_terms = c("mo(x)", "mo(f)", "z"),
                 ndraws = sim$S, verbose = FALSE)
  expect_equal(strip_attr(proj_linpred(prj, transform = FALSE)$pred),
               t(sim$eta_fun(sim$dat)), tolerance = 1e-4)
})

test_that(paste(
  "the augmented-data projection (cumulative family) reproduces the",
  "monotonic effects and predicts new data"
), {
  skip_if_not_installed("MASS")
  # Levels 1, ..., 5 (as in a Likert scale):
  sim <- mo_sim(x_min = 1L)
  thres <- c(-0.5, 0.8)
  dat <- sim$dat
  eta_mean <- rowMeans(sim$eta_fun(dat))
  dat$y_cumul <- factor(
    1L + (runif(nrow(dat)) > plogis(thres[1] - eta_mean)) +
      (runif(nrow(dat)) > plogis(thres[2] - eta_mean)),
    levels = 1:3
  )
  # The output of `ref_predfun`: an S x N x C_lat array of thres_c - eta:
  thres_arr <- function(eta) {
    arr <- vapply(thres, function(th) t(th - eta), t(eta))
    return(arr)
  }
  refm <- mo_refm(
    sim, y_cumul ~ mo(x) + mo(f) + z,
    family = structure(list(family = "cumulative_rstanarm", link = "logit"),
                       class = "family"),
    data = dat, ref_predfun_arr = thres_arr,
    augdat_y_unqs = as.character(1:3),
    augdat_link = augdat_link_cumul, augdat_ilink = augdat_ilink_cumul,
    augdat_args_link = list(link = "logit"),
    augdat_args_ilink = list(link = "logit")
  )
  # Monotonic projections do not support the augmented-data projection (the
  # error comes before the projection starts, that is, before get_refdist()):
  with_mocked_bindings(
    expect_error(
      project(refm, predictor_terms = c("mo(x)", "mo(f)", "z"), ndraws = 4,
              verbose = FALSE),
      "do not support the augmented-data projection"
    ),
    get_refdist = function(...) stop("The projection started.")
  )
  opts_old <- options(projpred.mo_monotone = FALSE)
  on.exit(options(opts_old), add = TRUE)
  prj <- muffle_polr_warning(project(
    refm, predictor_terms = c("mo(x)", "mo(f)", "z"), ndraws = sim$S,
    verbose = FALSE
  ))
  pred_train <- strip_attr(proj_linpred(prj, transform = FALSE)$pred)
  expect_equal(pred_train, thres_arr(sim$eta_fun(dat)), tolerance = 1e-3)
  nd <- mo_newdata(sim)
  pred_new <- proj_linpred(prj, newdata = dat[rows_of(nd), ],
                           transform = FALSE)$pred
  expect_equal(strip_attr(pred_new),
               pred_train[, rows_of(nd), , drop = FALSE], tolerance = 1e-10)
})

# Cross-validation --------------------------------------------------------

test_that("cv_varsel() works with mo() terms (LOO and K-fold CV)", {
  # The Pareto k warnings come from the small number of draws:
  opts_old <- options(projpred.warn_psis = FALSE)
  on.exit(options(opts_old), add = TRUE)
  sim <- mo_sim()
  refm <- mo_refm(sim, y_gauss ~ mo(x) + mo(f) + z + w)
  cv_loo <- muffle_warning(cv_varsel(
    refm, cv_method = "LOO", validate_search = TRUE, nloo = 5, nclusters = 3,
    nclusters_pred = 6, seed = 3906, verbose = FALSE
  ), "Some Pareto k diagnostic values are too high")
  expect_setequal(ranking(cv_loo)$fulldata, c("mo(x)", "mo(f)", "z", "w"))
  cv_kf <- cv_varsel(
    refm, cv_method = "kfold", K = 2, nclusters = 3, nclusters_pred = 6,
    seed = 3906, verbose = FALSE
  )
  expect_length(ranking(cv_kf)$foldwise, 2 * 4)
  expect_true(all(is.finite(summary(cv_kf)$perf_sub$elpd)))
})

test_that(paste(
  "K-fold CV throws an error if the training data of a fold lack the largest",
  "value of an integer mo() predictor"
), {
  sim <- mo_sim()
  # Only one observation with the largest value:
  is_max <- sim$dat$x == sim$D
  sim$dat$x[is_max] <- sim$D - 1L
  sim$dat$x[which(is_max)[1]] <- sim$D
  refm <- mo_refm(sim, y_gauss ~ mo(x) + z)
  expect_error(
    cv_varsel(refm, cv_method = "kfold", K = 2, nclusters = 3,
              nclusters_pred = 6, seed = 3906, verbose = FALSE),
    "outside the range of the training data"
  )
})
