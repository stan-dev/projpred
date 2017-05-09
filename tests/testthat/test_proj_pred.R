# tests for proj_linpred and proj_predict

set.seed(1235)
n <- 40
nv <- 5
x <- matrix(rnorm(n*nv, 0, 1), n, nv)
b <- runif(nv)-0.5
dis <- runif(1, 1, 2)
weights <- sample(1:4, n, replace = T)
offset <- rnorm(n)
chains <- 2
seed <- 1235
iter <- 500
source(file.path('helpers', 'SW.R'))


f_gauss <- gaussian()
df_gauss <- data.frame(y = rnorm(n, f_gauss$linkinv(x%*%b), dis), x = x)
f_binom <- binomial()
df_binom <- data.frame(y = rbinom(n, weights, f_binom$linkinv(x%*%b)), x = x)
f_poiss <- poisson()
df_poiss <- data.frame(y = rpois(n, f_poiss$linkinv(x%*%b)), x = x)
ys <- list()
ys[[1]] <- df_gauss$y
ys[[2]] <- df_binom$y/weights
ys[[3]] <- df_poiss$y

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

fit_list <- list(gauss = fit_gauss, binom = fit_binom, poiss = fit_poiss)
vs_list <- lapply(fit_list, varsel, nv_max = nv, verbose = FALSE)
proj_vind_list <- lapply(vs_list, project, vind = c(2,3), seed = seed)
proj_all_list <- lapply(vs_list, project, intercept = FALSE, seed = seed)

context('proj_linpred')
test_that("output of proj_linpred is sensible with fit-object as input", {
  for(i in 1:length(vs_list)) {
    i_inf <- names(vs_list)[i]
    pl <- proj_linpred(vs_list[[i]], xnew = x)
    expect_equal(length(pl), nv + 1, info = i_inf)
    for(j in 1:length(pl))
      expect_equal(ncol(pl[[j]]), n, info = i_inf)
  }
})

test_that("output of proj_linpred is sensible with project-object as input", {
  for(i in 1:length(proj_vind_list)) {
    i_inf <- names(proj_vind_list)[i]
    pl <- proj_linpred(proj_vind_list[[i]], xnew = x)
    expect_equal(ncol(pl), n, info = i_inf)
  }
  for(i in 1:length(proj_all_list)) {
    i_inf <- names(proj_all_list)[i]
    pl <- proj_linpred(proj_all_list[[i]], xnew = x)
    expect_equal(length(pl), nv + 1, info = i_inf)
    for(j in 1:length(pl))
      expect_equal(ncol(pl[[j]]), n, info = i_inf)
  }
})

test_that("proj_linpred: error when varsel has not been performed for the object", {
  expect_error(proj_linpred(1))
  expect_error(proj_linpred(fit_gauss))
})

test_that("proj_linpred: specifying ynew has an expected effect", {
  for(i in 1:length(vs_list)) {
    i_inf <- names(vs_list)[i]
    pl <- proj_linpred(vs_list[[i]], xnew = x, ynew = ys[[i]], weightsnew=weights)
    pl2 <- proj_linpred(vs_list[[i]], xnew = x, weightsnew=weights)
    for(j in 1:length(pl)) {
      expect_equal(names(pl[[j]]), c('pred', 'lpd'))
      expect_equal(ncol(pl[[j]]$pred), n, info = i_inf)
      expect_equal(ncol(pl[[j]]$lpd), n, info = i_inf)
      expect_equal(ncol(pl2[[j]]), n, info = i_inf)
    }
  }
})

test_that("proj_linpred: specifying weights has an expected effect", {
  for(i in 1:length(proj_vind_list)) {
    # for binomial models weights have to be specified
    if (proj_vind_list[[i]]$family_kl$family != 'binomial') {
      i_inf <- names(proj_vind_list)[i]
      plw <- proj_linpred(proj_vind_list[[i]], xnew = x, ynew = ys[[i]],
                          weightsnew = weights)
      pl <- proj_linpred(proj_vind_list[[i]], xnew = x, ynew = ys[[i]])
      expect_equal(names(plw), c('pred', 'lpd'), info = i_inf)
      expect_equal(ncol(plw$pred), n, info = i_inf)
      expect_equal(ncol(plw$lpd), n, info = i_inf)
      expect_true(sum(plw$lpd != pl$lpd) > 0, info = i_inf)
    }
  }
})

test_that("proj_linpred: specifying offset has an expected effect", {
  for(i in 1:length(proj_vind_list)) {
    i_inf <- names(proj_vind_list)[i]
    plo <- proj_linpred(proj_vind_list[[i]], xnew = x, ynew = ys[[i]], weightsnew=weights,
                        offsetnew = offset)
    pl <- proj_linpred(proj_vind_list[[i]], xnew = x, ynew = ys[[i]], weightsnew=weights)
    expect_equal(names(plo), c('pred', 'lpd'), info = i_inf)
    expect_equal(ncol(plo$pred), n, info = i_inf)
    expect_equal(ncol(plo$lpd), n, info = i_inf)
    expect_true(sum(plo$lpd != pl$lpd) > 0, info = i_inf)
  }
})

test_that("proj_linpred: specifying transform has an expected effect", {
  for(i in 1:length(proj_vind_list)) {
    i_inf <- names(proj_vind_list)[i]
    plt <- proj_linpred(proj_vind_list[[i]], xnew = x, transform = TRUE)
    plf <- proj_linpred(proj_vind_list[[i]], xnew = x, transform = FALSE)
    expect_equal(proj_vind_list[[i]]$family_kl$linkinv(plf), plt, info = i_inf)
  }
})

test_that("proj_linpred: specifying integrated has an expected effect", {
  for(i in 1:length(proj_vind_list)) {
    i_inf <- names(proj_vind_list)[i]
    plt <- proj_linpred(proj_vind_list[[i]], xnew = x, integrated = TRUE)
    plf <- proj_linpred(proj_vind_list[[i]], xnew = x, integrated = FALSE)
    expect_equal(drop(proj_vind_list[[i]]$weights%*%plf), plt, info = i_inf)
  }
})

test_that("proj_linpred: adding more regularization has an expected effect", {
  regul <- c(1e-6, 1e-1, 1e2)
  for(i in 1:length(vs_list)) {
    i_inf <- names(vs_list)[i]
    norms <- rep(0, length(regul))
    for (j in 1:length(regul)) {
      pred <- proj_linpred(vs_list[[i]], xnew = x, nv = 2, transform = FALSE,
                           integrated = TRUE, regul=regul[j])
      norms[j] <- sum(pred^2)
    }
    for (j in 1:(length(regul)-1))
      expect_true(all(norms[j] >= norms[j+1]), info = i_inf)
  }
})


test_that("proj_linpred: arguments passed to project work accordingly", {
  for(i in 1:length(vs_list)) {
    i_inf <- names(vs_list)[i]
    pr <- project(vs_list[[i]], nv = c(2, 4), nc = 2, ns = 20,
                  intercept = FALSE, regul = 1e-8, seed = 12)
    prl1 <- proj_linpred(pr, xnew = x)
    prl2 <- proj_linpred(vs_list[[i]], xnew = x, nv = c(2, 4), nc = 2, ns = 20,
                         intercept = FALSE, regul = 1e-8, seed = 12)
    expect_equal(prl1, prl2, info = i_inf)
  }
})


context('proj_predict')
test_that("output of proj_predict is sensible with fit-object as input", {
  for(i in 1:length(vs_list)) {
    i_inf <- names(vs_list)[i]
    pl <- proj_predict(vs_list[[i]], xnew = x)
    expect_equal(length(pl), nv + 1, info = i_inf)
    for(j in 1:length(pl))
      expect_equal(ncol(pl[[j]]), n, info = i_inf)
  }
})

test_that("output of proj_predict is sensible with project-object as input", {
  for(i in 1:length(proj_vind_list)) {
    i_inf <- names(proj_vind_list)[i]
    pl <- proj_predict(proj_vind_list[[i]], xnew = x)
    expect_equal(ncol(pl), n, info = i_inf)
  }
  for(i in 1:length(proj_all_list)) {
    i_inf <- names(proj_all_list)[i]
    pl <- proj_predict(proj_all_list[[i]], xnew = x)
    expect_equal(length(pl), nv + 1, info = i_inf)
    for(j in 1:length(pl))
      expect_equal(ncol(pl[[j]]), n, info = i_inf)
  }
})

test_that("proj_predict: error when varsel has not been performed for the object", {
  expect_error(proj_predict(1))
  expect_error(proj_predict(fit_gauss))
})

test_that("proj_predict: specifying weightsnew has an expected effect", {
  pl <- proj_predict(proj_vind_list[['binom']], xnew = x, seed = seed)
  plw <- proj_predict(proj_vind_list[['binom']], xnew = x, seed = seed,
                      weightsnew = weights)
  expect_true(sum(pl != plw)>0)
})

test_that("proj_predict: specifying offsetnew has an expected effect", {
  for(i in 1:length(proj_vind_list)) {
    i_inf <- names(proj_vind_list)[i]
    pl <- proj_predict(proj_vind_list[[i]], xnew = x, draws = iter,
                       seed = seed)
    plo <- proj_predict(proj_vind_list[[i]], xnew = x, draws = iter,
                        seed = seed, offsetnew = offset)
    expect_true(sum(pl != plo) > 0, info = i_inf)
  }
})

test_that("proj_predict: specifying draws has an expected effect", {
  for(i in 1:length(proj_vind_list)) {
    i_inf <- names(proj_vind_list)[i]
    pl <- proj_predict(proj_vind_list[[i]], xnew = x, draws = iter)
    expect_equal(dim(pl), c(iter, n))
  }
})

test_that("proj_predict: specifying seed_sam has an expected effect", {
  for(i in 1:length(proj_vind_list)) {
    i_inf <- names(proj_vind_list)[i]
    pl1 <- proj_predict(proj_vind_list[[i]], xnew = x, seed_samp = seed)
    pl2 <- proj_predict(proj_vind_list[[i]], xnew = x, seed_samp = seed)
    expect_equal(pl1, pl2, info = i_inf)
  }
})

test_that("proj_predict: arguments passed to project work accordingly", {
  for(i in 1:length(vs_list)) {
    i_inf <- names(vs_list)[i]
    pr1 <- project(vs_list[[i]], nv = c(2, 4), nc = 2, ns = 20,
                  intercept = FALSE, regul = 1e-8, seed = 12)
    prp1 <- proj_predict(pr1, xnew = x, draws = 100, seed_samp = 11)
    prp2 <- proj_predict(vs_list[[i]], xnew = x, draws = 100, seed_samp = 11,
                         nv = c(2, 4), nc = 2, ns = 20, intercept = FALSE,
                         regul = 1e-8, seed = 12)
    expect_equal(prp1, prp2, info = i_inf)
  }
})

