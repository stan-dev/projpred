# tests for varsel and cv_varsel:

set.seed(1235)
n <- 50
nv <- 5
x <- matrix(rnorm(n*nv, 0, 1), n, nv)
b <- runif(nv)-0.5
dis <- runif(1, 1, 2)
weights <- sample(1:4, n, replace = T)
chains <- 2
seed <- 1235
iter <- 500
offset <- rnorm(n)
source(file.path('helpers', 'SW.R'))

f_gauss <- gaussian()
df_gauss <- data.frame(y = rnorm(n, f_gauss$linkinv(x%*%b), dis), x = x)
f_binom <- binomial()
df_binom <- data.frame(y = rbinom(n, weights, f_binom$linkinv(x%*%b)), x = x)
f_poiss <- poisson()
df_poiss <- data.frame(y = rpois(n, f_poiss$linkinv(x%*%b)), x = x)

SW(
  fit_gauss <- stan_glm(y ~ x, family = f_gauss, data = df_gauss, QR = T,
                        weights = weights, offset = offset,
                        chains = chains, seed = seed, iter = iter)
)
SW(
  fit_binom <- stan_glm(cbind(y, weights-y) ~ x, family = f_binom, QR = T,
                        data = df_binom, weights = weights, offset = offset,
                        chains = chains, seed = seed, iter = iter)
)
SW(
  fit_poiss <- stan_glm(y ~ x, family = f_poiss, data = df_poiss, QR = T,
                        weights = weights, offset = offset,
                        chains = chains, seed = seed, iter = iter)
)
SW(
  fit_lm <- stan_lm(y ~ x, data = df_gauss, weights = weights, offset = offset,
                    prior = R2(0.3),
                    chains = chains, seed = seed, iter = iter)
)
SW(
  fit_glmer <- stan_glmer(mpg ~ wt + (1|cyl), data = mtcars,
                          chains = chains, seed = seed, iter = iter)
)
fit_list <- list(gauss = fit_gauss, binom = fit_binom, poiss = fit_poiss,
                 lm = fit_lm)

vsf <- function(x, m) varsel(x, method = m, nv_max = nv, verbose = FALSE)
vs_list <- list(l1 = lapply(fit_list, vsf, 'L1'),
                fs = lapply(fit_list, vsf, 'forward'))

# kfold not tested currently
cvsf <- function(x, m)
  cv_varsel(x, method = m, cv_method = 'LOO', nv_max = nv, verbose = FALSE)

SW(
  cvs_list <- list(l1 = lapply(fit_list, cvsf, 'L1'),
                   fs = lapply(fit_list, cvsf, 'forward'))
)


context('varsel')
test_that('varsel returns an object with a field named "varsel"', {
  for(i in length(vs_list)) {
    i_inf <- names(vs_list)[i]
    for(j in length(vs_list[[i]])) {
      j_inf <- names(vs_list[[i]])[j]
      expect_true('varsel' %in% names(vs_list[[i]][[j]]))
    }
  }
})

test_that('object retruned by varsel contains the relevant fields', {
  for(i in length(vs_list)) {
    i_inf <- names(vs_list)[i]
    for(j in length(vs_list[[i]])) {
      j_inf <- names(vs_list[[i]])[j]
      # chosen seems legit
      expect_equal(length(vs_list[[i]][[j]]$varsel$chosen), nv,
                   info = paste(i_inf, j_inf))
      # chosen_names seems legit
      expect_equal(names(coef(fit_gauss)[-1])[vs_list[[i]][[j]]$varsel$chosen],
                   vs_list[[i]][[j]]$varsel$chosen_names,
                   info = paste(i_inf, j_inf))
      # kl seems legit
      expect_equal(length(vs_list[[i]][[j]]$varsel$kl), nv + 1,
                   info = paste(i_inf, j_inf))
      # decreasing
      expect_equal(vs_list[[i]][[j]]$varsel$kl,
                   cummin(vs_list[[i]][[j]]$varsel$kl),
                   info = paste(i_inf, j_inf))
      # d_test seems legit
      expect_equal(length(vs_list[[i]][[j]]$varsel$d_test$y), n,
                   info = paste(i_inf, j_inf))
      expect_equal(length(vs_list[[i]][[j]]$varsel$d_test$weights), n,
                   info = paste(i_inf, j_inf))
      expect_equal(length(vs_list[[i]][[j]]$varsel$d_test$weights), n,
                   info = paste(i_inf, j_inf))
      expect_equal(typeof(vs_list[[i]][[j]]$varsel$d_test$type), 'character',
                   info = paste(i_inf, j_inf))
      expect_equal(cvs_list[[i]][[j]]$varsel$d_test$type, 'loo',
                   info = paste(i_inf, j_inf))
      # summaries seems legit
      expect_named(vs_list[[i]][[j]]$varsel$summaries, c('sub', 'full'),
                   info = paste(i_inf, j_inf))
      expect_equal(length(vs_list[[i]][[j]]$varsel$summaries$sub), nv + 1,
                   info = paste(i_inf, j_inf))
      expect_named(vs_list[[i]][[j]]$varsel$summaries$sub[[1]], c('mu', 'lppd'),
                   info = paste(i_inf, j_inf))
      expect_named(vs_list[[i]][[j]]$varsel$summaries$full, c('mu', 'lppd'),
                   info = paste(i_inf, j_inf))
      # family_kl seems legit
      expect_equal(vs_list[[i]][[j]]$family$family,
                   vs_list[[i]][[j]]$varsel$family_kl$family,
                   info = paste(i_inf, j_inf))
      expect_equal(vs_list[[i]][[j]]$family$link,
                   vs_list[[i]][[j]]$varsel$family_kl$link,
                   info = paste(i_inf, j_inf))
      expect_true(length(vs_list[[i]][[j]]$varsel$family_kl) >=
                    length(vs_list[[i]][[j]]$family$family),
                  info = paste(i_inf, j_inf))
    }
  }
})

test_that('nv_max has an effect on varsel for gaussian models', {
  vs1 <- varsel(fit_gauss, method = 'forward', nv_max = 3, verbose = FALSE)
  expect_equal(length(vs1$varsel$chosen), 3)
})

test_that('nv_max has an effect on varsel for non-gaussian models', {
  vs1 <- varsel(fit_binom, method = 'forward', nv_max = 3, verbose = FALSE)
  expect_equal(length(vs1$varsel$chosen), 3)
})

test_that('Having something else than stan_glm as the fit throws an error', {
  expect_error(varsel(fit_glmer, verbose = FALSE), regexp = 'supported')
  expect_error(varsel(1, verbose = FALSE), regexp = 'not recognized')
})


context('cv_varsel')
test_that('cv_varsel returns an object with a field named "varsel"', {
  for(i in length(cvs_list)){
    i_inf <- names(cvs_list)[i]
    for(j in length(cvs_list[[i]])) {
      j_inf <- names(cvs_list[[i]])[j]
      expect_true('varsel' %in% names(cvs_list[[i]][[j]]),
                  info = paste(i_inf, j_inf))
    }
  }

})

test_that('object retruned by cv_varsel contains the relevant fields', {
  for(i in length(cvs_list)) {
    i_inf <- names(cvs_list)[i]
    for(j in length(cvs_list[[i]])) {
      j_inf <- names(cvs_list[[i]])[j]
      # chosen seems legit
      expect_equal(length(cvs_list[[i]][[j]]$varsel$chosen), nv,
                   info = paste(i_inf, j_inf))
      # chosen_names seems legit
      expect_equal(names(coef(fit_gauss)[-1])[cvs_list[[i]][[j]]$varsel$chosen],
                   cvs_list[[i]][[j]]$varsel$chosen_names,
                   info = paste(i_inf, j_inf))
      # kl seems legit
      expect_equal(length(cvs_list[[i]][[j]]$varsel$kl), nv + 1,
                   info = paste(i_inf, j_inf))
      # decreasing
      expect_equal(cvs_list[[i]][[j]]$varsel$kl,
                   cummin(cvs_list[[i]][[j]]$varsel$kl),
                   info = paste(i_inf, j_inf))
      # d_test seems legit
      expect_equal(length(cvs_list[[i]][[j]]$varsel$d_test$y), n,
                   info = paste(i_inf, j_inf))
      expect_equal(length(cvs_list[[i]][[j]]$varsel$d_test$weights), n,
                   info = paste(i_inf, j_inf))
      expect_equal(length(cvs_list[[i]][[j]]$varsel$d_test$weights), n,
                   info = paste(i_inf, j_inf))
      expect_equal(typeof(cvs_list[[i]][[j]]$varsel$d_test$type), 'character',
                   info = paste(i_inf, j_inf))
      expect_equal(cvs_list[[i]][[j]]$varsel$d_test$type, 'loo',
                   info = paste(i_inf, j_inf))
      # summaries seems legit
      expect_named(cvs_list[[i]][[j]]$varsel$summaries, c('sub', 'full'),
                   info = paste(i_inf, j_inf))
      expect_equal(length(cvs_list[[i]][[j]]$varsel$summaries$sub), nv + 1,
                   info = paste(i_inf, j_inf))
      expect_named(cvs_list[[i]][[j]]$varsel$summaries$sub[[1]], c('mu', 'lppd'),
                   ignore.order = TRUE, info = paste(i_inf, j_inf))
      expect_named(cvs_list[[i]][[j]]$varsel$summaries$full, c('mu', 'lppd'),
                   ignore.order = TRUE, info = paste(i_inf, j_inf))
      # family_kl seems legit
      expect_equal(cvs_list[[i]][[j]]$family$family,
                   cvs_list[[i]][[j]]$varsel$family_kl$family,
                   info = paste(i_inf, j_inf))
      expect_equal(cvs_list[[i]][[j]]$family$link,
                   cvs_list[[i]][[j]]$varsel$family_kl$link,
                   info = paste(i_inf, j_inf))
      expect_true(length(cvs_list[[i]][[j]]$varsel$family_kl) >=
                    length(cvs_list[[i]][[j]]$family$family),
                  info = paste(i_inf, j_inf))
      # pctch seems legit
      expect_equal(dim(cvs_list[[i]][[j]]$varsel$pctch), c(nv, nv + 1),
                   info = paste(i_inf, j_inf))
      expect_true(all(cvs_list[[i]][[j]]$varsel$pctch[,-1] <= 1 &&
                        cvs_list[[i]][[j]]$varsel$pctch[,-1] >= 0),
                  info = paste(i_inf, j_inf))
      expect_equal(cvs_list[[i]][[j]]$varsel$pctch[,1], 1:nv,
                   info = paste(i_inf, j_inf))
      expect_equal(colnames(cvs_list[[i]][[j]]$varsel$pctch),
                   c('size', cvs_list[[i]][[j]]$varsel$chosen_names),
                   info = paste(i_inf, j_inf))
      # ssize seems legit
      expect_true(cvs_list[[i]][[j]]$varsel$ssize>=0,
                  info = paste(i_inf, j_inf))
    }
  }
})


test_that('nv_max has an effect on cv_varsel for gaussian models', {
  suppressWarnings(
    vs1 <- cv_varsel(fit_gauss, method = 'forward', nv_max = 3, verbose = FALSE)
  )
  expect_equal(length(vs1$varsel$chosen), 3)
})

test_that('nv_max has an effect on cv_varsel for non-gaussian models', {
  suppressWarnings(
    vs1 <- varsel(fit_binom, method = 'forward', nv_max = 3, verbose = FALSE)
  )
  expect_equal(length(vs1$varsel$chosen), 3)
})

test_that('Having something else than stan_glm as the fit throws an error', {
  expect_error(cv_varsel(fit_glmer, verbose = FALSE), regexp = 'supported')
  expect_error(varsel(1, verbose = FALSE), regexp = 'not recognized')
})

test_that('kfold cv throws an error', {
  expect_error(cv_varsel(fit_gauss, cv_method = 'kfold', verbose = FALSE), regexp = 'unavailable')
})
