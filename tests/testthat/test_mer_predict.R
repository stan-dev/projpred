context("mer_predict")

if (exists(".Random.seed", envir = .GlobalEnv)) {
  rng_old <- get(".Random.seed", envir = .GlobalEnv)
}
set.seed(13579)
n_mp <- 150
dat_mp <- data.frame(
  x1 = rnorm(n_mp), x2 = rnorm(n_mp),
  f = factor(sample(c("a", "b", "c"), n_mp, replace = TRUE)),
  g = factor(paste0("g", sample(8, n_mp, replace = TRUE))),
  h = factor(paste0("h", sample(4, n_mp, replace = TRUE)))
)
eta_mp <- with(dat_mp, 0.5 * x1 - 0.3 * x2 + rnorm(8)[as.integer(g)] +
                 0.5 * rnorm(4)[as.integer(h)])
# Draw-wise fits of one submodel (as from a projection):
fits_mp <- function(fml, family = NULL, n_draws = 4) {
  lapply(seq_len(n_draws), function(s) {
    dat_s <- dat_mp
    dat_s$y <- eta_mp + rnorm(n_mp, 0, 0.3)
    if (is.null(family)) {
      suppressMessages(lme4::lmer(fml, data = dat_s))
    } else {
      dat_s$y <- plogis(dat_s$y)
      suppressWarnings(lme4::glmer(fml, data = dat_s, family = family))
    }
  })
}
# The per-draw predictions of subprd() on master:
prd_per_draw <- function(fits, newdata) {
  do.call(cbind, lapply(fits, function(fit) {
    predict(fit, newdata = newdata, re.form = NULL, allow.new.levels = TRUE) +
      repair_re(fit, newdata = newdata)
  }))
}

test_that("subprd_mer_direct() gives the predictions of the draw-wise path", {
  cases <- list(
    lmer_slope = fits_mp(y ~ x1 + x2 + f + (1 + x1 | g) + (1 | h)),
    lmer_icpt = fits_mp(y ~ x1 + (1 | g)),
    lmer_nofix = fits_mp(y ~ 1 + (1 | g) + (1 | h)),
    glmer = fits_mp(y ~ x1 + (1 | g) + (1 | h), family = binomial()),
    # Special formula (functions in the fixed-effects terms):
    lmer_spcl = fits_mp(y ~ poly(x1, 2, raw = TRUE) +
                          sqrt(abs(x2)^2) * I(!as.logical(x2 > 0)) + f +
                          (x1 | g))
  )
  newdat <- dat_mp[sample(n_mp, 30), ]
  for (nm in names(cases)) {
    for (nd in list(dat_mp, newdat)) {
      set.seed(1)
      prd_expected <- prd_per_draw(cases[[nm]], nd)
      seed_expected <- .Random.seed
      set.seed(1)
      prd <- subprd_mer_direct(cases[[nm]], nd)
      expect_false(is.null(prd), info = nm)
      expect_equal(prd, prd_expected, tolerance = 1e-12, info = nm)
      # No random numbers are drawn:
      expect_identical(.Random.seed, seed_expected, info = nm)
      expect_equal(subprd(cases[[nm]], nd), prd_expected, tolerance = 1e-12,
                   info = nm)
    }
  }
})

test_that("subprd_mer_direct() returns `NULL` for other cases", {
  fits <- fits_mp(y ~ x1 + (1 | g))
  # New group level:
  nd_new <- dat_mp[1:10, ]
  nd_new$g <- factor(c("gnew", as.character(nd_new$g[-1])))
  expect_null(subprd_mer_direct(fits, nd_new))
  # All group levels are treated as new:
  opts_old <- options(projpred.mlvl_pred_new = TRUE)
  expect_null(subprd_mer_direct(fits, dat_mp))
  options(opts_old)
  # Other fits, no `newdata`, fits of different submodels:
  expect_null(subprd_mer_direct(list(lm(x1 ~ x2, data = dat_mp)), dat_mp))
  expect_null(subprd_mer_direct(fits, NULL))
  fits_2 <- c(fits[1], fits_mp(y ~ x2 + (1 | g), n_draws = 1))
  expect_null(subprd_mer_direct(fits_2, dat_mp))
})

if (exists("rng_old")) assign(".Random.seed", rng_old, envir = .GlobalEnv)
