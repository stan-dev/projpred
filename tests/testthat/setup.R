seed <- 1235
set.seed(seed)
n <- 40L
nterms <- 5L
x <- matrix(rnorm(n * nterms, 0, 1), n, nterms)
b <- runif(nterms) - 0.5
dis <- runif(1L, 1, 2)
weights <- sample(1:4, n, replace = TRUE)
offset <- rnorm(n)
chains <- 2L
iter <- 500L
source(testthat::test_path("helpers", "SW.R"))

f_gauss <- gaussian()
df_gauss <- data.frame(y = rnorm(n, f_gauss$linkinv(x %*% b), dis), x = x)
f_binom <- binomial()
df_binom <- data.frame(y = rbinom(n, weights, f_binom$linkinv(x %*% b)),
                       x = x, weights = weights, offset = offset)
fam_nms <- setNames(nm = c("gauss", "binom"))
ys <- lapply(fam_nms, function(fam_nm) {
  get(paste0("df_", fam_nm))$y
})
SW({
  fit_gauss <- rstanarm::stan_glm(
    y ~ x.1 + x.2 + x.3 + x.4 + x.5,
    family = f_gauss, data = df_gauss,
    weights = weights, offset = offset,
    chains = chains, seed = seed, iter = iter, QR = TRUE
  )
  fit_binom <- rstanarm::stan_glm(
    cbind(y, weights - y) ~ x.1 + x.2 + x.3 + x.4 + x.5,
    family = f_binom, data = df_binom,
    offset = offset, # `weights` is not needed when using the cbind() syntax
    chains = chains, seed = seed, iter = iter
  )
})
fit_list <- lapply(fam_nms, function(fam_nm) {
  get(paste0("fit_", fam_nm))
})
