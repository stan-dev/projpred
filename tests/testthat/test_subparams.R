# These tests run on CRAN on purpose. They fit small lme4 and ordinal models
# directly (no reference model) and check the group-level variance components
# that get_subparams() returns for both naming schemes.

# The variance components of get_subparams() output:
vc_subparams <- function(subparams) {
  subparams[grep("^Sigma\\[|^sd_|^cor_", names(subparams))]
}

# lme4::sleepstudy with a quadratic term and a second grouping factor `grp2`
# (crossed with `Subject`) that has a group-level standard deviation of 20.
# The state of the random number generator is restored at the end.
dat_sleep <- function() {
  withr::local_seed(5)
  d <- lme4::sleepstudy
  d$Days2 <- (d$Days - 4.5)^2 / 10
  d$grp2 <- factor(paste0("h", rep(1:5, length.out = nrow(d))))
  d$Reaction <- d$Reaction + rnorm(nlevels(d$grp2), sd = 20)[d$grp2]
  d
}

test_that(paste(
  "get_subparams() returns the variances and covariances in the rstanarm",
  "order for `nm_scheme = \"rstanarm\"` (lmer() fit)"
), {
  fit <- lme4::lmer(Reaction ~ Days + (Days + Days2 | Subject) + (1 | grp2),
                    data = dat_sleep())
  vc <- lme4::VarCorr(fit)
  vc_s <- vc$Subject
  vc_g <- vc$grp2
  expect_gt(vc_g[1, 1], 1)
  vc_expected <- c(
    "Sigma[Subject:(Intercept),(Intercept)]" = vc_s[1, 1],
    "Sigma[Subject:Days,(Intercept)]" = vc_s[2, 1],
    "Sigma[Subject:Days2,(Intercept)]" = vc_s[3, 1],
    "Sigma[Subject:Days,Days]" = vc_s[2, 2],
    "Sigma[Subject:Days2,Days]" = vc_s[3, 2],
    "Sigma[Subject:Days2,Days2]" = vc_s[3, 3],
    "Sigma[grp2:(Intercept),(Intercept)]" = vc_g[1, 1]
  )
  # The order of the groups follows lme4::VarCorr():
  vc_expected <- vc_expected[order(
    match(sub("^Sigma\\[(.*):.*$", "\\1", names(vc_expected)), names(vc))
  )]
  expect_equal(vc_subparams(get_subparams(fit, nm_scheme = "rstanarm")),
               vc_expected, tolerance = 1e-12)
  # The variances are the squared standard deviations:
  expect_equal(unname(vc_expected[c(1, 4, 6)]),
               unname(attr(vc_s, "stddev")^2), tolerance = 1e-12)
})

test_that(paste(
  "get_subparams() returns the standard deviations and correlations for",
  "`nm_scheme = \"brms\"` (lmer() fit)"
), {
  fit <- lme4::lmer(Reaction ~ Days + (Days | Subject), data = dat_sleep())
  vc_s <- lme4::VarCorr(fit)$Subject
  sds <- attr(vc_s, "stddev")
  vc_expected <- c(
    "sd_Subject__Intercept" = sds[[1]],
    "sd_Subject__Days" = sds[[2]],
    "cor_Subject__Intercept__Days" = attr(vc_s, "correlation")[2, 1]
  )
  expect_equal(vc_subparams(get_subparams(fit, nm_scheme = "brms")),
               vc_expected, tolerance = 1e-12)
})

test_that(paste(
  "get_subparams() returns the group-level variance for",
  "`nm_scheme = \"rstanarm\"` (glmer() fit)"
), {
  fit <- lme4::glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
                     data = lme4::cbpp, family = binomial())
  sd_herd <- attr(lme4::VarCorr(fit)$herd, "stddev")[[1]]
  expect_gt(sd_herd, 0.1)
  expect_equal(vc_subparams(get_subparams(fit, nm_scheme = "rstanarm")),
               c("Sigma[herd:(Intercept),(Intercept)]" = sd_herd^2),
               tolerance = 1e-12)
  expect_equal(vc_subparams(get_subparams(fit, nm_scheme = "brms")),
               c("sd_herd__Intercept" = sd_herd),
               tolerance = 1e-12)
})

test_that(paste(
  "get_subparams() returns the group-level variance for",
  "`nm_scheme = \"rstanarm\"` (clmm() fit)"
), {
  fit <- ordinal::clmm(rating ~ temp + (1 | judge), data = ordinal::wine)
  sd_judge <- attr(ordinal::VarCorr(fit)$judge, "stddev")[[1]]
  expect_gt(sd_judge, 0.1)
  expect_equal(vc_subparams(get_subparams(fit, nm_scheme = "rstanarm")),
               c("Sigma[judge:(Intercept),(Intercept)]" = sd_judge^2),
               tolerance = 1e-12)
  expect_equal(vc_subparams(get_subparams(fit, nm_scheme = "brms")),
               c("sd_judge__Intercept" = sd_judge),
               tolerance = 1e-12)
})
