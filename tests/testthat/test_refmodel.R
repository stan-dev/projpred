context('refmodel')


# tests for generic reference model

if (require(rstanarm) && require(brms)) {
  

  seed <- 1235
  set.seed(seed)
  n <- 50
  nv <- 5
  x <- matrix(rnorm(n*nv, 0, 1), n, nv)
  b <- runif(nv)-0.5
  dis <- runif(1, 1, 2)
  weights <- sample(1:4, n, replace = T)
  chains <- 2
  iter <- 500
  offset <- rnorm(n)
  source(file.path('helpers', 'SW.R'))
  
  f_gauss <- gaussian()
  df_gauss <- data.frame(y = rnorm(n, f_gauss$linkinv(x%*%b), dis), x = x)
  f_binom <- binomial()
  df_binom <- data.frame(y = rbinom(n, weights, f_binom$linkinv(x%*%b)), x = x, weights=weights)
  
  SW({
    fit_gauss <- stan_glm(y ~ x, family = f_gauss, data = df_gauss,
                          chains = chains, seed = seed, iter = iter)
    fit_binom <- brm(y | trials(weights) ~ x, family = f_binom,
                     data = df_binom chains = chains, seed = seed, iter = iter)
    ref_gauss <- get_refmodel(fit_gauss)
    ref_binom <- get_refmodel(fit_binom)
  })
  
  test_that('get_refmodel produces sensible results', {
    expect_s3_class(ref_gauss, "refmodel")
    expect_s3_class(ref_binom, "refmodel")
  })
  
  test_that('get_refmode checks for the absence of data', {
    SW({
    fit_nodata <- stan_glm(df_gauss$y ~ x, family = f_gauss, QR = T,
                           weights = weights, offset = offset,
                           chains = chains, seed = seed, iter = iter)
    })
    expect_error(get_refmodel(fit_nodata),
                 'Model was fitted without a \'data\' argument')
  })
  
  test_that('predict checks the \'type\' argument', {
    expect_error(predict(ref_gauss, df_gauss, type = 'zzz'),
                 '\'arg\' should be one of')
  })
  
  test_that('predict produces sensible results for gaussian models', {
    out.resp <- predict(ref_gauss, df_gauss, type = 'response')
    expect_vector(out.resp)
    expect_length(out.resp, nrow(df_gauss))
  
    out.link <- predict(ref_gauss, df_gauss, type = 'link')
    expect_equal(out.resp, out.link)
  })
  
  test_that('predict produces sensible results for binomial models', {
    out.resp <- predict(ref_binom, df_binom, type = 'response')
    expect_vector(out.resp)
    expect_length(out.resp, nrow(df_binom))
    expect_true(all(out.resp >= 0 & out.resp <= 1))
  
    out.link <- predict(ref_binom, df_binom, type = 'link')
    expect_length(out.resp, nrow(df_binom))
  })
  
  test_that('predict produces sensible results when specifying ynew', {
    out <- predict(ref_gauss, df_gauss, ynew = df_gauss$y)
    expect_vector(out)
    expect_length(out, length(df_gauss$y))
  
    expect_error(predict(ref_gauss, df_gauss, ynew = df_gauss),
                 'must be a numerical vector')
  })

}
