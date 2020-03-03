context("datafit")
suppressWarnings(RNGversion("3.5.0"))


# tests for data based estimates (no actual reference model)

if (!requireNamespace("glmnet", quietly = TRUE)) {
  stop("glmnet needed for this function to work. Please install it.",
    call. = FALSE
  )
}

set.seed(1235)
n <- 40
nv <- 5
x <- matrix(rnorm(n * nv, 0, 1), n, nv)
b <- runif(nv) - 0.5
dis <- runif(1, 1, 2)
weights <- sample(1:4, n, replace = T)
offset <- rnorm(n)
chains <- 2
seed <- 1235
iter <- 500
source(file.path("helpers", "SW.R"))

f_gauss <- gaussian()
df_gauss <- data.frame(y = rnorm(n, f_gauss$linkinv(x %*% b), dis), x = x)
f_binom <- binomial()
df_binom <- data.frame(y = rbinom(n, weights, f_binom$linkinv(x %*% b)), x = x)
f_poiss <- poisson()
df_poiss <- data.frame(y = rpois(n, f_poiss$linkinv(x %*% b)), x = x)

formula <- y ~ x.1 + x.2 + x.3 + x.4 + x.5

dref_gauss <- init_refmodel(fit = NULL, df_gauss, df_gauss$y, formula, gaussian(),
                                offset = offset, weights = weights)
dref_binom <- init_refmodel(fit = NULL, df_binom, df_binom$y / weights, formula,
                                binomial(), offset = offset, weights = weights)
dref_poiss <- init_refmodel(fit = NULL, df_poiss, df_poiss$y, formula, poisson(),
                                offset = offset, weights = weights)

dref_list <- list(gauss = dref_gauss, binom = dref_binom, poiss = dref_poiss)

# varsel
vsd_list <- lapply(dref_list, varsel, nv_max = nv + 1, verbose = FALSE)

# cv_varsel
cvvsd_list <- lapply(dref_list, cv_varsel, nv_max = nv + 1, verbose = FALSE)

#
predd_list <- lapply(vsd_list, proj_linpred, xnew = data.frame(x=x, weights=weights),
                     offsetnew = offset, weightsnew = weights, nv = 3,
                     seed = seed)

test_that("predict fails for 'datafit' objects", {
  expect_error(
    predict(dref_gauss, df_gauss),
    "Cannot make predictions with data reference only"
  )
})

test_that("output of varsel is sensible with only data provided as reference model", {
  for (i in seq_along(vsd_list)) {
    # vind seems legit
    expect_equal(length(vsd_list[[i]]$vind), nv)

    # kl seems legit
    expect_equal(length(vsd_list[[i]]$kl), nv + 1)

    # kl decreasing
    expect_equal(vsd_list[[i]]$kl, cummin(vsd_list[[i]]$kl), tolerance=1e-2)

    # summaries seems legit
    expect_named(vsd_list[[i]]$summaries, c("sub", "ref"))
    expect_equal(length(vsd_list[[i]]$summaries$sub), nv + 1)
    expect_named(vsd_list[[i]]$summaries$sub[[1]], c("mu", "lppd"))
    expect_named(vsd_list[[i]]$summaries$ref, c("mu", "lppd"))
  }
})

test_that("output of cv_varsel is sensible with only data provided as reference model", {
  for (i in seq_along(cvvsd_list)) {
    # vind seems legit
    expect_equal(length(cvvsd_list[[i]]$vind), nv)

    # kl seems legit
    expect_equal(length(cvvsd_list[[i]]$kl), nv + 1)

    # kl decreasing
    expect_equal(cvvsd_list[[i]]$kl, cummin(cvvsd_list[[i]]$kl), tolerance=1e-2)

    # summaries seems legit
    expect_named(cvvsd_list[[i]]$summaries, c("sub", "ref"))
    expect_equal(length(cvvsd_list[[i]]$summaries$sub), nv + 1)
    expect_named(cvvsd_list[[i]]$summaries$sub[[1]], c("mu", "lppd"))
    expect_named(cvvsd_list[[i]]$summaries$ref, c("mu", "lppd"))
  }
})

test_that("varsel_stats stops if baseline = 'ref' and deltas = TRUE", {
  expect_error(
    varsel_stats(vsd_list[[1]], baseline = "ref", deltas = TRUE),
    "Cannot use deltas = TRUE and baseline = 'ref' when there is no reference model"
  )
})

test_that("output of project is sensible with only data provided as reference model", {
  for (i in 1:length(vsd_list)) {

    # length of output of project is legit
    p <- project(vsd_list[[i]], nv = 0:nv)
    expect_equal(length(p), nv + 1)

    for (j in 1:length(p)) {
      expect_named(p[[j]], c("kl", "weights", "dis", "vind",
                             "sub_fit", "p_type", "family"),
      ignore.order = TRUE
      )
      # number of draws should equal to the number of draw weights
      ns <- length(p[[j]]$weights)
      expect_equal(length(p[[j]]$sub_fit$alpha), ns)
      expect_equal(length(p[[j]]$dis), ns)
      if (j > 1)
        expect_equal(ncol(p[[j]]$sub_fit$beta), ns)
      # j:th element should have j-1 variables
      expect_equal(length(which(p[[j]]$sub_fit$beta != 0)), j - 1)
      expect_equal(length(p[[j]]$vind), j - 1)
      # family kl
      expect_equal(p[[j]]$family, vsd_list[[i]]$family)
    }
    # kl should be non-increasing on training data
    klseq <- sapply(p, function(e) e$kl)
    expect_equal(klseq, cummin(klseq), tolerance=1e-2)

    # all submodels should use the same clustering/subsampling
    expect_equal(p[[1]]$weights, p[[nv]]$weights)
  }
})


test_that("output of proj_linpred is sensible with only data provided as reference model", {
  for (i in 1:length(vsd_list)) {

    # length of output of project is legit
    pred <- proj_linpred(vsd_list[[i]],
      xnew = data.frame(x=x), seed = seed,
      offsetnew = offset, weightsnew = weights, nv = 3
    )
    expect_equal(length(pred$pred), nrow(x))

    pred <- proj_linpred(vsd_list[[i]],
      xnew = data.frame(x=x), ynew = dref_list[[i]]$y, seed = seed,
      offsetnew = offset, weightsnew = weights, nv = 3
    )

    expect_equal(length(pred$pred), nrow(x))
    expect_equal(length(pred$lpd), nrow(x))
  }
})


# below are some tests that check Lasso solution computed with varsel is the same
# as that of glmnet. (notice that glm_ridge and glm_elnet are already tested separately, so
# these would only check that the results do not change due to varsel/cv_varsel etc.)


set.seed(1235)
n <- 100
nv <- 10
x <- matrix(rnorm(n * nv, 0, 1), n, nv)
b <- seq(0, 1, length.out = nv)
dis <- runif(1, 0.3, 0.5)
weights <- sample(1:4, n, replace = T) #
offset <- 0.1 * rnorm(n)
seed <- 1235
source(file.path("helpers", "SW.R"))



fams <- list(gaussian(), binomial(), poisson())
x_list <- lapply(fams, function(fam) x)
y_list <- lapply(fams, function(fam) {
  if (fam$family == "gaussian") {
    y <- rnorm(n, x %*% b, 0.5)
    weights <- NULL
    y_glmnet <- y
  } else if (fam$family == "binomial") {
    y <- rbinom(n, weights, fam$linkinv(x %*% b))
    y <- y / weights
    y_glmnet <- cbind(1 - y, y) # different way of specifying binomial y for glmnet
    weights <- weights
  } else if (fam$family == "poisson") {
    y <- rpois(n, fam$linkinv(x %*% b))
    y_glmnet <- y
    weights <- NULL
  }
  nlist(y, y_glmnet, weights)
})


test_that("L1-projection with data reference gives the same results as Lasso from glmnet.", {
  for (i in seq_along(fams)) {
    x <- x_list[[i]]
    y <- y_list[[i]]$y
    y_glmnet <- y_list[[i]]$y_glmnet
    fam <- fams[[i]]
    weights <- y_list[[i]]$weights
    if (is.null(weights))
      weights <- rep(1, NROW(y))

    lambda_min_ratio <- 1e-7
    nlambda <- 1500

    df <- data.frame(y=y, x=x)
    formula <- y ~ x.1 + x.2 + x.3 + x.4 + x.5 + x.6 + x.7 + x.8 + x.9 + x.10
    # Lasso solution with projpred
    ref <- init_refmodel(NULL, df, y, formula, family = fam, weights = weights, offset = offset)
    vs <- varsel(ref, method = "l1", lambda_min_ratio = lambda_min_ratio, nlambda = nlambda, thresh = 1e-12)
    pred1 <- proj_linpred(vs, xnew = data.frame(x=x), nv = 0:nv,
                              transform = FALSE, offsetnew = offset)

    # compute the results for the Lasso
    lasso <- glmnet::glmnet(x, y_glmnet,
      family = fam$family, weights = weights, offset = offset,
      lambda.min.ratio = lambda_min_ratio, nlambda = nlambda, thresh = 1e-12)
    vind <- predict(lasso, type = "nonzero", s = lasso$lambda)
    nselected <- sapply(vind, function(e) length(e))
    lambdainds <- sapply(unique(nselected), function(nv) max(which(nselected == nv)))
    lambdaval <- lasso$lambda[lambdainds]
    pred2 <- predict(lasso, newx = x, type = "link", s = lambdaval, newoffset = offset)

    # check that the predictions agree (up to nv-2 only, because glmnet terminates the coefficient
    # path computation too early for some reason...)
    for (j in 1:(nv - 2)) {
      expect_true(median(abs(pred1[[j]]$pred - pred2[, j])) < 3e-1)
    }

    # check that the coefficients are similar
    ind <- match(vs$vind, setdiff(split_formula(formula), "1"))
    betas <- sapply(vs$spath$sub_fits, function(x) x$beta %||% 0)
    delta <- sapply(seq_len(nv), function(i)
                    abs(t(betas[[i + 1]]) - lasso$beta[ind[1:i], lambdainds[i + 1]]))
    expect_true(median(unlist(delta)) < 6e-2)
    expect_true(median(abs(sapply(vs$spath$sub_fits, function(x) x$alpha) -
                           lasso$a0[lambdainds])) < 1.5e-1)
  }
})
