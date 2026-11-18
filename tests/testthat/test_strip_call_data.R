context("strip_call_data")

# Data with the response columns of several draws, as the submodel fitters get
# them:
n_scd <- 150L
idx_scd <- seq_len(n_scd)
dat_scd <- data.frame(x1 = sin(idx_scd), x2 = cos(1.3 * idx_scd),
                      g = factor(paste0("g", idx_scd %% 6L)))
u_scd <- c(0.8, -0.5, 0.3, -0.9, 0.6, -0.2)
eta_scd <- dat_scd$x1 + 0.5 * dat_scd$x2 + u_scd[as.integer(dat_scd$g)]
n_draws_scd <- 30L
for (s in seq_len(n_draws_scd)) {
  dat_scd[[paste0(".y.", s)]] <- eta_scd +
    0.3 * sin(2.3 * (s + 0.37) * idx_scd + 0.5)
}
dat_scd$yb <- as.integer(eta_scd + 1.2 * sin(2.9 * idx_scd + 0.2) > 0)
dat_scd$yo <- cut(eta_scd + 1.2 * sin(2.9 * idx_scd + 0.2),
                  c(-Inf, -0.5, 0.5, Inf),
                  labels = c("lo", "mid", "hi"), ordered_result = TRUE)
dat_scd$yc <- factor(dat_scd$yo, ordered = FALSE)
newdat_scd <- dat_scd[c(3, 50, 77, 120), c("x1", "x2", "g")]

# Whether a call of `x` (or of its list elements, one level deep) holds a
# data frame as element `data`:
has_call_data <- function(x, depth = 0L) {
  cl <- if (isS4(x)) {
    if (methods::.hasSlot(x, "call")) x@call
  } else if (is.list(x) && !is.data.frame(x)) {
    x[["call"]]
  }
  if (is.call(cl) && is.data.frame(cl$data)) {
    return(TRUE)
  }
  if (depth < 1L && is.list(x) && !isS4(x) && !is.data.frame(x)) {
    return(any(vapply(x, has_call_data, logical(1), depth = depth + 1L)))
  }
  return(FALSE)
}

test_that(paste(
  "strip_call_data() removes the data frame from the calls of the submodel",
  "fits and keeps their predictions for new data"
), {
  # Each case: the fit (through do.call(), as in the fitters of projpred) and
  # the predictions for new data (`NULL`: no predict() method).
  cases <- list(
    lmer = list(pkg = "lme4", fit = function() {
      do.call(lme4::lmer, list(formula = .y.1 ~ x1 + x2 + (1 | g),
                               data = dat_scd))
    }, pred = function(f) predict(f, newdata = newdat_scd)),
    glmer = list(pkg = "lme4", fit = function() {
      do.call(lme4::glmer, list(formula = yb ~ x1 + (1 | g), data = dat_scd,
                                family = binomial()))
    }, pred = function(f) predict(f, newdata = newdat_scd)),
    gam = list(pkg = "mgcv", fit = function() {
      do.call(mgcv::gam, list(formula = .y.1 ~ s(x1) + x2, data = dat_scd))
    }, pred = function(f) predict(f, newdata = newdat_scd)),
    # gamm4() fits keep no data frame in their calls:
    gamm4 = list(pkg = "gamm4", fit = function() {
      do.call(gamm4::gamm4, list(formula = .y.1 ~ s(x1) + x2,
                                 random = ~ (1 | g), data = dat_scd))
    }, pred = function(f) predict(f$gam, newdata = newdat_scd),
    no_data = TRUE),
    glmmPQL = list(pkg = "MASS", fit = function() {
      do.call(MASS::glmmPQL, list(fixed = yb ~ x1, random = ~ 1 | g,
                                  family = binomial(), data = dat_scd,
                                  verbose = FALSE))
    }, pred = function(f) predict(f, newdata = newdat_scd, level = 0)),
    polr = list(pkg = "MASS", fit = function() {
      do.call(MASS::polr, list(formula = yo ~ x1 + x2, data = dat_scd))
    }, pred = function(f) predict(f, newdata = newdat_scd, type = "probs")),
    clmm = list(pkg = "ordinal", fit = function() {
      do.call(ordinal::clmm, list(formula = yo ~ x1 + (1 | g),
                                  data = dat_scd))
    }, pred = NULL),
    multinom = list(pkg = "nnet", fit = function() {
      do.call(nnet::multinom, list(formula = yc ~ x1 + x2, data = dat_scd,
                                   trace = FALSE))
    }, pred = function(f) predict(f, newdata = newdat_scd, type = "probs"))
  )
  for (nm in names(cases)) {
    cs <- cases[[nm]]
    if (!requireNamespace(cs$pkg, quietly = TRUE)) next
    fit <- suppressMessages(suppressWarnings(cs$fit()))
    if (isTRUE(cs$no_data)) {
      expect_false(has_call_data(fit), info = nm)
      expect_identical(strip_call_data(fit), fit, info = nm)
      next
    }
    expect_true(has_call_data(fit), info = nm)
    fit_s <- strip_call_data(fit)
    expect_false(has_call_data(fit_s), info = nm)
    expect_lt(length(serialize(fit_s, NULL)), length(serialize(fit, NULL)),
              label = nm)
    if (!is.null(cs$pred)) {
      expect_identical(cs$pred(fit_s), cs$pred(fit), info = nm)
    }
  }
})

test_that("strip_call_data() leaves other objects as they are", {
  expect_identical(strip_call_data(NULL), NULL)
  expect_identical(strip_call_data(dat_scd), dat_scd)
  x <- list(a = 1, call = quote(f(data = d)), b = list(call = quote(g())))
  expect_identical(strip_call_data(x), x)
})

test_that("divmin() gives multilevel and additive fits without the data", {
  fml_resp <- paste0("cbind(", paste0(".y.", seq_len(n_draws_scd),
                                      collapse = ", "), ")")
  divmin_scd <- function(rhs) {
    suppressWarnings(divmin(
      formula = as.formula(paste(fml_resp, "~", rhs)), data = dat_scd,
      family = extend_family(gaussian()), weights = rep(1, n_scd),
      projpred_var = matrix(NA, n_scd, n_draws_scd),
      projpred_ws_aug = matrix(NA, n_scd, n_draws_scd), verbose_divmin = FALSE
    ))
  }
  for (rhs in c("x1 + x2 + (1 | g)", "s(x1) + x2")) {
    fits <- divmin_scd(rhs)
    expect_length(fits, n_draws_scd)
    expect_false(any(vapply(fits, has_call_data, logical(1))), info = rhs)
  }
})
