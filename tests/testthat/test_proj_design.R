# The cache entries of proj_to_submodl() from the column blocks of one model
# matrix of all main-effect terms (proj_prepare_design()) must equal the
# entries of the general path (proj_prepare_general()).

# Data with numeric, factor (with sum contrasts), ordered, character and
# logical columns:
dat_design <- function(n = 120, seed = 6173) {
  set.seed(seed)
  dat <- as.data.frame(matrix(rnorm(n * 5), n, 5,
                              dimnames = list(NULL, paste0("x", 1:5))))
  dat$x5 <- abs(dat$x5) + 0.1
  dat$f1 <- factor(sample(c("a", "b", "c"), n, replace = TRUE))
  dat$f2 <- factor(sample(c("p", "q", "r", "s"), n, replace = TRUE))
  dat$o1 <- factor(sample(c("lo", "mid", "hi"), n, replace = TRUE),
                   levels = c("lo", "mid", "hi"), ordered = TRUE)
  dat$c1 <- sample(c("u", "v"), n, replace = TRUE)
  dat$l1 <- sample(c(TRUE, FALSE), n, replace = TRUE)
  dat$y <- dat$x1 - 0.5 * dat$x2 + c(a = 0, b = 0.7, c = -0.4)[dat$f1] +
    rnorm(n)
  contrasts(dat$f2) <- contr.sum(4)
  return(dat)
}

fml_design <- y ~ x1 + x2 + poly(x3, 2) + log(x5) + I(x4^2) + f1 + f2 + o1 +
  c1 + l1 + x1:x2

# A list with the parts of a reference model that proj_prepare_design() and
# proj_prepare_general() use:
refm_parts <- function(dat, fml) {
  ctrs <- contrasts_from_data(dat, vars = all.vars(fml[[3]]))
  dat_store <- strip_contrasts(dat)
  list(formula = fml, fetch_data = function(...) dat_store, contrasts = ctrs)
}

# A cache store in which the next miss builds the design:
new_store <- function() {
  store <- new.env(parent = emptyenv())
  store$design_misses <- 2L
  return(store)
}

# The cache entry of the general path:
entry_general <- function(refm, terms_, y) {
  return(proj_prepare_general(refm, terms_, y)[c("fml_divmin", "prep")])
}

# Checks that the cache entry `ent` equals the entry `ent_ref`. Both paths call
# poly(), which computes a QR decomposition. Some BLAS libraries (for example
# MKL) can round two calls with the same data differently, so the model matrix
# and the coefficients of poly() in the `predvars` of the terms can differ in
# the last bits. All other parts must be identical:
expect_entry_equal <- function(ent, ent_ref, info) {
  blas_parts <- function(e) {
    pv <- as.list(attr(e$prep$shared_x$terms, "predvars"))[-1]
    pv <- lapply(pv, function(el) if (is.call(el)) as.list(el) else el)
    list(x = e$prep$shared_x$x, predvars = pv)
  }
  other_parts <- function(e) {
    e$prep$shared_x$x <- NULL
    attr(e$prep$shared_x$terms, "predvars") <- NULL
    e
  }
  expect_identical(other_parts(ent), other_parts(ent_ref), info = info)
  expect_equal(blas_parts(ent), blas_parts(ent_ref), tolerance = 1e-12,
               info = info)
}

# Checks that proj_prepare_design() gives the entry of the general path for
# each element of `sets` that it covers, and returns the number of these sets.
expect_design_entries <- function(refm, sets, y, store = new_store()) {
  n_design <- 0L
  for (s in sets) {
    ent <- proj_prepare_design(store, refm, s, y)
    if (is.null(ent)) {
      next
    }
    n_design <- n_design + 1L
    expect_entry_equal(ent, entry_general(refm, s, y),
                       info = paste(c(s, ncol(y)), collapse = " "))
  }
  return(n_design)
}

test_that(paste(
  "proj_prepare_design() gives the cache entries of the general path, and",
  "NULL for interactions"
), {
  dat <- dat_design()
  refm <- refm_parts(dat, fml_design)
  cand <- labels(terms(fml_design))
  set.seed(9031)
  sets <- c(list(character(), cand[cand != "x1:x2"]),
            lapply(1:60, function(i) sample(cand, sample.int(length(cand), 1))))
  sets_main <- Filter(function(s) !"x1:x2" %in% s, sets)
  y_k <- matrix(runif(nrow(dat) * 5), ncol = 5)
  y_rn <- y_k
  rownames(y_rn) <- paste0("obs", seq_len(nrow(dat)))
  for (y in list(y_k, y_rn, y_k[, 1, drop = FALSE])) {
    expect_identical(expect_design_entries(refm, sets, y), length(sets_main))
    store <- new_store()
    expect_null(proj_prepare_design(store, refm, c("x1", "x1:x2"), y))
    expect_true(store$design_ok)
  }
  # One store, responses with other row names in turn:
  store <- new_store()
  for (y in list(y_rn, y_k, y_rn, y_k[, 1, drop = FALSE])) {
    expect_identical(expect_design_entries(refm, sets_main[1:10], y,
                                           store = store), 10L)
  }
})

test_that("proj_prepare_design() builds the design at the third miss", {
  dat <- dat_design()
  refm <- refm_parts(dat, fml_design)
  y <- matrix(runif(nrow(dat) * 3), ncol = 3)
  store <- new.env(parent = emptyenv())
  expect_null(proj_prepare_design(store, refm, "x1", y))
  expect_null(proj_prepare_design(store, refm, "x2", y))
  expect_null(store$design_ok)
  expect_identical(proj_prepare_design(store, refm, "x1", y),
                   entry_general(refm, "x1", y))
})

test_that("proj_prepare_design() covers other columns and formulas", {
  dat <- dat_design()
  y <- matrix(runif(nrow(dat) * 3), ncol = 3)
  # Missing values only in a column that is not a candidate:
  dat_na <- dat
  dat_na$z <- dat$x1
  dat_na$z[3] <- NA
  refm_na <- refm_parts(dat_na, y ~ x1 + x2 + f1)
  store <- new_store()
  expect_identical(
    expect_design_entries(refm_na, list("x1", c("f1", "x2")), y, store), 2L
  )
  expect_true(store$design_ok)
  # A data-dependent transformation (scale() keeps the center and the scale of
  # the training data in "predvars"), a matrix column, a factor with a level
  # that no row has (and a contrast matrix):
  dat_sp <- dat
  dat_sp$m <- cbind(dat$x1, dat$x2)
  dat_sp$f3 <- factor(sample(c("a", "b"), nrow(dat), replace = TRUE),
                      levels = c("a", "b", "c"))
  contrasts(dat_sp$f3) <- contr.sum(3)
  refm_sp <- refm_parts(dat_sp, y ~ scale(x3) + m + f3 + x1)
  expect_identical(
    expect_design_entries(refm_sp, list("scale(x3)", "m", c("f3", "x1"),
                                        c("m", "scale(x3)")),
                          y),
    4L
  )
  # A reference model formula with a group-level term: the design covers the
  # submodels without it:
  dat_g <- dat
  dat_g$g <- factor(rep_len(paste0("g", 1:6), nrow(dat)))
  refm_g <- refm_parts(dat_g, y ~ x1 + x2 + (1 | g))
  store <- new_store()
  expect_identical(expect_design_entries(refm_g, list("x1", c("x1", "x2")), y,
                                         store),
                   2L)
  expect_null(proj_prepare_design(store, refm_g, c("x1", "(1 | g)"), y))
  # Another `options("contrasts")` gives another design:
  refm <- refm_parts(dat, fml_design)
  store <- new_store()
  expect_identical(expect_design_entries(refm, list(c("o1", "f1")), y, store),
                   1L)
  opts_old <- options(contrasts = c("contr.sum", "contr.treatment"))
  on.exit(options(opts_old), add = TRUE)
  expect_identical(expect_design_entries(refm, list(c("o1", "f1")), y, store),
                   1L)
  options(opts_old)
})

test_that("proj_prepare_design() takes the general path outside its gate", {
  dat <- dat_design()
  y <- matrix(runif(nrow(dat) * 3), ncol = 3)
  # A missing value in a candidate column:
  dat_na <- dat
  dat_na$x2[3] <- NA
  store <- new_store()
  expect_null(proj_prepare_design(store, refm_parts(dat_na, fml_design),
                                  c("x1", "f1"), y))
  expect_false(store$design_ok)
  # A factor with one level (model.matrix() fails for it, but not for the
  # submodels without it):
  dat_one <- dat
  dat_one$f0 <- factor(rep("a", nrow(dat)))
  refm_one <- refm_parts(dat_one, y ~ x1 + x2 + f0)
  store <- new_store()
  expect_null(proj_prepare_design(store, refm_one, "x1", y))
  expect_false(store$design_ok)
  expect_false(is.null(entry_general(refm_one, "x1", y)$prep))
  # Column names that data.frame() changes:
  dat_nm <- dat
  names(dat_nm)[names(dat_nm) == "x1"] <- "x 1"
  fml_nm <- y ~ `x 1` + x2
  expect_null(proj_prepare_design(new_store(), refm_parts(dat_nm, fml_nm),
                                  "x2", y))
  # A term that the reference model formula does not have:
  expect_null(proj_prepare_design(new_store(), refm_parts(dat, y ~ x1 + x2),
                                  "x3", y))
})

test_that(paste(
  "varsel() and cv_varsel() give the same results with and without the",
  "design path"
), {
  skip_on_cran()
  # A Bernoulli reference model (a Gaussian one would use the search of
  # search_forward_gauss(), which does not call proj_to_submodl()):
  dat <- dat_design()
  set.seed(2207)
  dat$y <- rbinom(nrow(dat), 1, plogis(dat$y))
  fml_ref <- y ~ x1 + x2 + poly(x3, 2) + log(x5) + I(x4^2) + f1 + f2 + o1 +
    c1 + l1
  dat_dflt <- strip_contrasts(dat)
  tt <- delete.response(terms(fml_ref))
  mm <- model.matrix(tt, dat_dflt)
  ndraws <- 50
  beta_hat <- glm.fit(mm, dat$y, family = binomial())$coefficients
  expect_false(anyNA(beta_hat))
  beta <- beta_hat + matrix(rnorm(ncol(mm) * ndraws, sd = 0.05), ncol(mm))
  xlev <- .getXlevels(tt, model.frame(tt, dat_dflt))
  fit <- list(beta = beta, data = dat_dflt, tt = tt, xlev = xlev)
  ref_predfun <- function(fit, newdata = NULL) {
    nd <- strip_contrasts(if (is.null(newdata)) fit$data else newdata)
    model.matrix(fit$tt, model.frame(fit$tt, nd, xlev = fit$xlev)) %*%
      fit$beta
  }
  refm <- suppressWarnings(init_refmodel(
    fit, data = dat, formula = fml_ref, family = binomial(),
    ref_predfun = ref_predfun
  ))
  run_vs <- function() {
    list(
      vs = varsel(refm, method = "forward", nclusters = 3,
                  nclusters_pred = 5, verbose = FALSE, seed = 5502),
      cvvs = suppressWarnings(cv_varsel(
        refm, method = "forward", cv_method = "LOO", nloo = 8,
        validate_search = TRUE, nclusters = 3, nclusters_pred = 5,
        verbose = FALSE, seed = 5502
      ))
    )
  }
  # Count the calls of proj_prepare_design() that gave an entry:
  calls <- new.env()
  calls$n <- 0L
  ns <- asNamespace("projpred")
  suppressMessages(trace(
    "proj_prepare_design", where = ns, print = FALSE,
    exit = bquote(if (!is.null(returnValue(NULL))) {
      assign("n", get("n", envir = .(calls)) + 1L, envir = .(calls))
    })
  ))
  on.exit(suppressMessages(untrace("proj_prepare_design", where = ns)),
          add = TRUE)
  res_design <- run_vs()
  expect_gt(calls$n, 0L)
  res_general <- with_mocked_bindings(
    run_vs(),
    proj_prepare_design = function(...) NULL
  )
  for (nm in names(res_design)) {
    expect_identical(ranking(res_design[[nm]])$fulldata,
                     ranking(res_general[[nm]])$fulldata, info = nm)
    expect_identical(summary(res_design[[nm]])$perf_sub,
                     summary(res_general[[nm]])$perf_sub, info = nm)
    expect_identical(res_design[[nm]]$summaries,
                     res_general[[nm]]$summaries, info = nm)
  }
})
