# tests for varsel and cv_varsel:

set.seed(1235)
n <- 40
nv <- 5
x <- matrix(rnorm(n*nv, 0, 1), n, nv)
b <- runif(nv)-0.5
dis <- runif(1, 1, 2)
weights <- sample(1:4, n, replace = T)
chains <- 1
cores <- 1
seed <- 1235
iter <- 500
# change this to something else once offsets work
offset <- rep(0, n)

f_gauss <- gaussian()
df_gauss <- data.frame(y = rnorm(n, f_gauss$linkinv(x%*%b), dis), x = x)
f_binom <- binomial()
df_binom <- data.frame(y = rbinom(n, weights, f_binom$linkinv(x%*%b)), x = x)
f_poiss <- poisson()
df_poiss <- data.frame(y = rpois(n, f_poiss$linkinv(x%*%b)), x = x)

fit_gauss <- stan_glm(y ~ x, family = f_gauss, data = df_gauss, QR = T,
                      weights = weights, offset = offset,
                      chains = chains, cores = cores, seed = seed, iter = iter)
fit_binom <- stan_glm(cbind(y, weights-y) ~ x, family = f_binom, QR = T,
                      data = df_binom, weights = weights, offset = offset,
                      chains = chains, cores = cores, seed = seed, iter = iter)
fit_poiss <- stan_glm(y ~ x, family = f_poiss, data = df_poiss, QR = T,
                      weights = weights, offset = offset,
                      chains = chains, cores = cores, seed = seed, iter = iter)
suppressWarnings(
  fit_glmer <- stan_glmer(mpg ~ wt + (1|cyl), data = mtcars, chains = chains,
                          cores = cores, seed = seed, iter = iter)
)
fit_list <- list(fit_gauss, fit_binom, fit_poiss)

vsf <- function(x, m) varsel(x, method = m, nv_max = nv, verbose = FALSE)
vs_list <- list(l1 = lapply(fit_list, vsf, 'L1'),
                fs = lapply(fit_list, vsf, 'forward'))

# kfold not tested currently
cvsf <- function(x, m) cv_varsel(x, method = m, nv_max = nv, verbose = FALSE)
suppressWarnings(
  cvs_list <- list(l1 = lapply(fit_list, cvsf, 'L1'),
                   fs = lapply(fit_list, cvsf, 'forward'))
)


context("varsel")
test_that("varsel returns an object with a field named 'varsel'", {
  for(i in length(vs_list))
    for(j in length(vs_list[[i]]))
      expect_true('varsel' %in% names(vs_list[[i]][[j]]))
})

test_that("object retruned by varsel contains relevant fields with L1-search", {

  for(i in length(vs_list)) {
    for(j in length(vs_list[[i]])) {
      # chosen seems legit
      expect_length(vs_list[[i]][[j]]$varsel$chosen, nv)
      # chosen_names seems legit
      expect_equal(names(coef(fit_gauss)[-1])[vs_list[[i]][[j]]$varsel$chosen],
                   vs_list[[i]][[j]]$varsel$chosen_names)
      # kl seems legit
      expect_length(vs_list[[i]][[j]]$varsel$kl, nv + 1)
      # d_test seems legit
      expect_length(vs_list[[i]][[j]]$varsel$d_test$y, n)
      expect_length(vs_list[[i]][[j]]$varsel$d_test$weights, n)
      expect_length(vs_list[[i]][[j]]$varsel$d_test$weights, n)
      expect_type(vs_list[[i]][[j]]$varsel$d_test$type, 'character')
      expect_equal(cvs_list[[i]][[j]]$varsel$d_test$type, 'loo')
      # summaries seems legit
      expect_named(vs_list[[i]][[j]]$varsel$summaries, c('sub', 'full'))
      expect_length(vs_list[[i]][[j]]$varsel$summaries$sub, nv + 1)
      expect_named(vs_list[[i]][[j]]$varsel$summaries$sub[[1]], c('mu', 'lppd'))
      expect_named(vs_list[[i]][[j]]$varsel$summaries$full, c('mu', 'lppd'))
      # family_kl seems legit
      expect_equal(vs_list[[i]][[j]]$family$family,
                   vs_list[[i]][[j]]$varsel$family_kl$family)
      expect_equal(vs_list[[i]][[j]]$family$link,
                   vs_list[[i]][[j]]$varsel$family_kl$link)
      expect_gt(length(vs_list[[i]][[j]]$varsel$family_kl),
                length(vs_list[[i]][[j]]$family$family))
    }
  }
})

test_that("nv_max has an effect on varsel for gaussian models", {
  vs1 <- varsel(fit_gauss, method = 'forward', nv_max = 3, verbose = FALSE)
  expect_length(vs1$varsel$chosen, 3)
})

test_that("nv_max has an effect on varsel for non-gaussian models", {
  vs1 <- varsel(fit_binom, method = 'forward', nv_max = 3, verbose = FALSE)
  expect_length(vs1$varsel$chosen, 3)
})

test_that("Having something else than stan_glm as the fit throws an error", {
  expect_error(varsel(fit_glmer, verbose = FALSE), regexp = 'supported')
  expect_error(varsel(1, verbose = FALSE), regexp = 'not a stanreg')
})


context('cv_varsel')
test_that("cv_varsel returns an object with a field named 'varsel'", {
  for(i in length(cvs_list))
    for(j in length(cvs_list[[i]]))
      expect_true('varsel' %in% names(cvs_list[[i]][[j]]))
})

test_that("object retruned by cv_varsel contains relevant fields", {
  for(i in length(cvs_list)) {
    for(j in length(cvs_list[[i]])) {
      # chosen seems legit
      expect_length(cvs_list[[i]][[j]]$varsel$chosen, nv)
      # chosen_names seems legit
      expect_equal(names(coef(fit_gauss)[-1])[cvs_list[[i]][[j]]$varsel$chosen],
                   cvs_list[[i]][[j]]$varsel$chosen_names)
      # kl seems legit
      expect_length(cvs_list[[i]][[j]]$varsel$kl, nv + 1)
      # d_test seems legit
      expect_length(cvs_list[[i]][[j]]$varsel$d_test$y, n)
      expect_length(cvs_list[[i]][[j]]$varsel$d_test$weights, n)
      expect_length(cvs_list[[i]][[j]]$varsel$d_test$weights, n)
      expect_type(cvs_list[[i]][[j]]$varsel$d_test$type, 'character')
      expect_equal(cvs_list[[i]][[j]]$varsel$d_test$type, 'loo')
      # summaries seems legit
      expect_named(cvs_list[[i]][[j]]$varsel$summaries, c('sub', 'full'))
      expect_length(cvs_list[[i]][[j]]$varsel$summaries$sub, nv + 1)
      expect_named(cvs_list[[i]][[j]]$varsel$summaries$sub[[1]], c('mu', 'lppd'),
                   ignore.order = TRUE)
      expect_named(cvs_list[[i]][[j]]$varsel$summaries$full, c('mu', 'lppd'),
                   ignore.order = TRUE)
      # family_kl seems legit
      expect_equal(cvs_list[[i]][[j]]$family$family,
                   cvs_list[[i]][[j]]$varsel$family_kl$family)
      expect_equal(cvs_list[[i]][[j]]$family$link,
                   cvs_list[[i]][[j]]$varsel$family_kl$link)
      expect_gt(length(cvs_list[[i]][[j]]$varsel$family_kl),
                length(cvs_list[[i]][[j]]$family$family))
      # pctch seems legit
      expect_equal(dim(cvs_list[[i]][[j]]$varsel$pctch), c(nv, nv + 1))
      expect_true(all(cvs_list[[i]][[j]]$varsel$pctch[,-1] <= 1 &&
                        cvs_list[[i]][[j]]$varsel$pctch[,-1] >= 0))
      expect_equal(cvs_list[[i]][[j]]$varsel$pctch[,1], 1:nv)
      expect_equal(colnames(cvs_list[[i]][[j]]$varsel$pctch),
                   c('size', cvs_list[[i]][[j]]$varsel$chosen_names))
      # ssize seems legit
      expect_true(cvs_list[[i]][[j]]$varsel$ssize>=0)
    }
  }
})


test_that("nv_max has an effect on cv_varsel for gaussian models", {
  suppressWarnings(
    vs1 <- cv_varsel(fit_gauss, method = 'forward', nv_max = 3, verbose = FALSE)
  )
  expect_length(vs1$varsel$chosen, 3)
})

test_that("nv_max has an effect on cv_varsel for non-gaussian models", {
  suppressWarnings(
    vs1 <- varsel(fit_binom, method = 'forward', nv_max = 3, verbose = FALSE)
  )
  expect_length(vs1$varsel$chosen, 3)
})

test_that("Having something else than stan_glm as the fit throws an error", {
  expect_error(cv_varsel(fit_glmer, verbose = FALSE), regexp = 'supported')
  expect_error(varsel(1, verbose = FALSE), regexp = 'not a stanreg')
})

test_that("kfold cv throws an error", {
  expect_error(cv_varsel(fit_gauss, cv_method = 'kfold', verbose = FALSE), regexp = 'unavailable')
})
