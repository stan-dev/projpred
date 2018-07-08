# tests for generic reference model

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


dref_gauss <- init_refmodel(x,df_gauss$y,gaussian(),offset=offset,wobs=weights)
dref_binom <- init_refmodel(x,df_binom$y/weights,binomial(),offset=offset,wobs=weights)
dref_poiss <- init_refmodel(x,df_poiss$y,poisson(),offset=offset,wobs=weights)
dref_list <- list(gauss = dref_gauss, binom = dref_binom, poiss = dref_poiss)

# varsel
vsd_list <- lapply(dref_list, varsel, nv_max = nv, verbose = FALSE)

# cv_varsel
cvvsd_list <- lapply(dref_list, cv_varsel, nv_max = nv, verbose = FALSE)

# 
predd_list <- lapply(vsd_list, proj_linpred, xnew=x, seed = seed,
										 offsetnew=offset, weightsnew=weights, nv=3)


context('init_refmodel')

test_that('output of varsel is sensible with only data provided as reference model', {
  for(i in seq_along(vsd_list)) {
    # vind seems legit
    expect_equal(length(vsd_list[[i]]$vind), nv)
    
    # kl seems legit
    expect_equal(length(vsd_list[[i]]$kl), nv + 1)
    
    # kl decreasing
    expect_equal(vsd_list[[i]]$kl, cummin(vsd_list[[i]]$kl))
    
    # summaries seems legit
    expect_named(vsd_list[[i]]$summaries, c('sub', 'full'))
    expect_equal(length(vsd_list[[i]]$summaries$sub), nv + 1)
    expect_named(vsd_list[[i]]$summaries$sub[[1]], c('mu', 'lppd'))
    expect_named(vsd_list[[i]]$summaries$full, c('mu', 'lppd'))
  }
})

test_that("output of cv_varsel is sensible with only data provided as reference model", {
  for(i in seq_along(cvvsd_list)) {
    # vind seems legit
    expect_equal(length(cvvsd_list[[i]]$vind), nv)
    
    # kl seems legit
    expect_equal(length(cvvsd_list[[i]]$kl), nv + 1)
    
    # kl decreasing
    expect_equal(cvvsd_list[[i]]$kl, cummin(cvvsd_list[[i]]$kl))
    
    # summaries seems legit
    expect_named(cvvsd_list[[i]]$summaries, c('sub', 'full'))
    expect_equal(length(cvvsd_list[[i]]$summaries$sub), nv + 1)
    expect_named(cvvsd_list[[i]]$summaries$sub[[1]], c('mu', 'lppd'))
    expect_named(cvvsd_list[[i]]$summaries$full, c('mu', 'lppd'))
  }
})

test_that("output of project is sensible with only data provided as reference model", {
  for(i in 1:length(vsd_list)) {
    
    # length of output of project is legit
    p <- project(vsd_list[[i]], nv=0:nv)
    expect_equal(length(p), nv + 1)
    
    for(j in 1:length(p)) {
      expect_named(p[[j]], c('kl', 'weights', 'dis', 'alpha', 'beta', 'vind',
                             'p_type', 'intercept', 'family_kl'),
                   ignore.order = T)
      # number of draws should equal to the number of draw weights
      ns <- length(p[[j]]$weights)
      expect_equal(length(p[[j]]$alpha), ns)
      expect_equal(length(p[[j]]$dis), ns)
      expect_equal(ncol(p[[j]]$beta), ns)
      # j:th element should have j-1 variables
      expect_equal(nrow(p[[j]]$beta), j-1)
      expect_equal(length(p[[j]]$vind), j-1)
      # family kl
      expect_equal(p[[j]]$family_kl, vsd_list[[i]]$family_kl)
    }
    # kl should be non-increasing on training data
    klseq <- sapply(p, function(e) e$kl)
    expect_equal(klseq, cummin(klseq))
    
    # all submodels should use the same clustering/subsampling
    expect_equal(p[[1]]$weights, p[[nv]]$weights)
  }
})


test_that("output of proj_linpred is sensible with only data provided as reference model", {
  for(i in 1:length(vsd_list)) {
    
    # length of output of project is legit
    # p <- project(vsd_list[[i]], nv=0:nv)
    pred <- proj_linpred(vsd_list[[i]], xnew=x, seed = seed,
                        offsetnew=offset, weightsnew=weights, nv=3)
    expect_equal(length(pred), nrow(x))
    
    pred <- proj_linpred(vsd_list[[i]], xnew=x, ynew=dref_list[[i]]$y, seed = seed,
                         offsetnew=offset, weightsnew=weights, nv=3)
    
    expect_equal(length(pred$pred), nrow(x))
    expect_equal(length(pred$lpd), nrow(x))
  }
})
