context('project')


# tests for project

if (require(rstanarm)) {
  

  seed <- 1235
  set.seed(seed)
  n <- 40
  nv <- 5
  x <- matrix(rnorm(n*nv, 0, 1), n, nv)
  b <- runif(nv)-0.5
  dis <- runif(1, 1, 2)
  weights <- sample(1:4, n, replace = T)
  offset <- rnorm(n)
  chains <- 2
  iter <- 500
  source(file.path('helpers', 'SW.R'))
  
  f_gauss <- gaussian()
  df_gauss <- data.frame(y = rnorm(n, f_gauss$linkinv(x%*%b), dis), x = x)
  f_binom <- binomial()
  df_binom <- data.frame(y = rbinom(n, weights, f_binom$linkinv(x%*%b)), x = x)
  f_poiss <- poisson()
  df_poiss <- data.frame(y = rpois(n, f_poiss$linkinv(x%*%b)), x = x)
  
  SW({
  fit_gauss <- stan_glm(y ~ x, family = f_gauss, data = df_gauss, QR = T,
                        weights = weights, offset = offset,
                        chains = chains, seed = seed, iter = iter)
  fit_binom <- stan_glm(cbind(y, weights-y) ~ x, family = f_binom, QR = T,
                        data = df_binom, weights = weights, offset = offset,
                        chains = chains, seed = seed, iter = iter)
  fit_poiss <- stan_glm(y ~ x, family = f_poiss, data = df_poiss, QR = T,
                        weights = weights, offset = offset,
                        chains = chains, seed = seed, iter = iter)
  })
  fit_list <- list(fit_gauss, fit_binom, fit_poiss)
  vs_list <- lapply(fit_list, varsel, nv_max = nv, verbose = FALSE)
  
  
  test_that("project: relaxing has the expected effect", {
    
    vs_list <- lapply(fit_list, varsel, nv_max = nv, verbose = FALSE, method='l1')
    for (i in seq_along(vs_list)) {
      
      p0 <- project(vs_list[[i]], relax=F, nv=1:nv)
      p1 <- project(vs_list[[i]], relax=T, nv=1:nv, nc=1, regul=1e-9)
      
      for (j in seq_along(p1)) {
        # L1-penalised coefficients should have smaller L1-norm
        expect_true( sum(abs(p0[[j]]$beta)) < sum(abs(p1[[j]]$beta)) )
      }
    }
  })
  
  test_that("object returned by project contains the relevant fields", {
    for(i in 1:length(vs_list)) {
      i_inf <- names(vs_list)[i]
      p <- project(vs_list[[i]], nv=0:nv)
      expect_type(p, "list")
      expect_length(p, nv + 1)
  
      for(j in 1:length(p)) {
        expect_s3_class(p[[j]], "projection")
        expect_named(p[[j]], c('kl', 'weights', 'dis', 'alpha', 'beta', 'vind',
                               'p_type', 'intercept', 'family_kl'),
                     ignore.order = T, info = i_inf)
        # number of draws should equal to the number of draw weights
        ns <- length(p[[j]]$weights)
        expect_length(p[[j]]$alpha, ns)
        expect_length(p[[j]]$dis, ns)
        expect_equal(ncol(p[[j]]$beta), ns, info = i_inf)
        # j:th element should have j-1 variables
        expect_equal(nrow(p[[j]]$beta), j - 1, info = i_inf)
        expect_length(p[[j]]$vind, j - 1)
        # family kl
        expect_equal(p[[j]]$family_kl, vs_list[[i]]$family_kl,
                     info = i_inf)
      }
      # kl should be non-increasing on training data
      klseq <- sapply(p, function(x) x$kl)
      expect_equal(klseq, cummin(klseq), info = i_inf)
  
      # all submodels should use the same clustering
      expect_equal(p[[1]]$weights, p[[nv]]$weights, info = i_inf)
    }
  })
  
  test_that("project: error when varsel has not been performed for the object", {
    expect_error(project(1, xnew = x),
                 'is not a variable selection object')
    expect_error(project(fit_gauss, xnew = x),
                 'is not a variable selection object')
  })
  
  test_that("project: nv is checked", {
    expect_error(project(vs_list[[1]], nv = 1000),
                 'Cannot perform the projection with 1000 variables')
    expect_error(project(vs_list[[1]], nv = -1),
                 'must contain non-negative values')
    expect_error(project(vs_list[[1]], nv = 'a'),
                 'must contain non-negative values')
    expect_error(project(vs_list[[1]], nv = df_gauss),
                 'must contain non-negative values')
  })
  
  test_that("project: setting nv = NULL has the expected effect", {
    for(i in 1:length(vs_list)) {
      i_inf <- names(vs_list)[i]
      p <- project(vs_list[[i]], nv = NULL)
      # if only one model size is projected, do not return a list of length one
      expect_true(length(p) >= 1, info = i_inf)
      # beta has the correct number of rows
      expect_equal(nrow(p$beta), vs_list[[i]]$ssize, info = i_inf)
      expect_length(p$vind, vs_list[[i]]$ssize)
    }
  })
  
  test_that("project: setting nv = 0 has an expected effect", {
    for(i in 1:length(vs_list)) {
      i_inf <- names(vs_list)[i]
      nv <- 0
      p <- project(vs_list[[i]], nv = nv)
      # if only one model size is projected, do not return a list of length one
      expect_true(length(p) >= 1, info = i_inf)
      # beta has the correct number of rows
      expect_equal(nrow(p$beta), nv, info = i_inf)
      expect_length(p$vind, nv)
    }
  })
  
  test_that("project: setting nv = 3 has an expected effect", {
    for(i in 1:length(vs_list)) {
      i_inf <- names(vs_list)[i]
      nv <- 3
      p <- project(vs_list[[i]], nv = nv)
      # if only one model is projected, do not return a list of length one
      expect_true(length(p) >= 1, info = i_inf)
      # beta has the correct number of rows
      expect_equal(nrow(p$beta), nv, info = i_inf)
      expect_length(p$vind, nv)
    }
  })
  
  test_that("project: setting vind to 4 has an expected effect", {
    for(i in 1:length(vs_list)) {
      vind <- 4
      p <- project(vs_list[[i]], vind = vind)
      expect_equivalent(p$vind, vind)
      expect_equal(nrow(p$beta), 1)
      exp_ind <- which(vs_list[[i]]$vind == vind)
      expect_named(p$vind, names(vs_list[[i]]$vind)[exp_ind])
    }
  })
  
  test_that("project: setting vind to 1:2 has an expected effect", {
    for(i in 1:length(vs_list)) {
      # i_inf <- names(vs_list)[i]
      vind <- 1:2
      # names(vind) <- names(coef(vs_list[[i]]))[vind+1]
      p <- project(vs_list[[i]], vind = vind)
      expect_equivalent(p$vind, vind)
      expect_equal(nrow(p$beta), length(vind), info = i_inf)
      exp_ind <- sapply(vind, function(x) which(vs_list[[i]]$vind == x))
      expect_named(p$vind, names(vs_list[[i]]$vind)[exp_ind])
    }
  })
  
  test_that("project: setting vind to something nonsensical returns an error", {
    # variable selection objects
    expect_error(project(vs_list[[1]], vind = 1:10),
                 'vind contains an index larger than')
    expect_error(project(vs_list[[1]], vind = 17),
                 'vind contains an index larger than')
  
    # fit objects
    expect_error(project(fit_list[[1]], vind = 1:10),
                 'vind contains an index larger than')
    expect_error(project(fit_list[[1]], vind = 17),
                 'vind contains an index larger than')
  })
  
  test_that("project: setting ns to 1 has an expected effect", {
    for(i in 1:length(vs_list)) {
      i_inf <- names(vs_list)[i]
      ns <- 1
      p <- project(vs_list[[i]], ns = ns, nv = nv)
      # expected number of draws
      expect_length(p$weights, ns)
      expect_length(p$alpha, ns)
      expect_length(p$dis, ns)
      expect_equal(ncol(p$beta), ns, info = i_inf)
      expect_equal(p$weights, 1, info = i_inf)
    }
  })
  
  test_that("project: setting ns to 40 has an expected effect", {
    for(i in 1:length(vs_list)) {
      i_inf <- names(vs_list)[i]
      ns <- 40
      p <- project(vs_list[[i]], ns = ns, nv = nv)
      # expected number of draws
      expect_length(p$weights, ns)
      expect_length(p$alpha, ns)
      expect_length(p$dis, ns)
      expect_equal(ncol(p$beta), ns, info = i_inf)
  
      # no clustering, so draw weights should be identical
      expect_true(do.call(all.equal, as.list(p$weights)), info = i_inf)
    }
  })
  
  test_that("project: setting nc to 1 has an expected effect", {
    for(i in 1:length(vs_list)) {
      i_inf <- names(vs_list)[i]
      nc <- 1
      p <- project(vs_list[[i]], nc = nc, nv = nv)
      # expected number of draws
      expect_length(p$weights, nc)
      expect_length(p$alpha, nc)
      expect_length(p$dis, nc)
      expect_equal(ncol(p$beta), nc, info = i_inf)
    }
  })
  
  test_that("project: setting nc to 20 has an expected effect", {
    for(i in 1:length(vs_list)) {
      i_inf <- names(vs_list)[i]
      nc <- 20
      p <- project(vs_list[[i]], nc = nc, nv = nv)
      # expected number of draws
      expect_length(p$weights, nc)
      expect_length(p$alpha, nc)
      expect_length(p$dis, nc)
      expect_equal(ncol(p$beta), nc, info = i_inf)
    }
  })
  
  test_that("project: setting ns or nc to too big throws an error", {
    expect_error(project(vs_list[[1]], ns = 400000, nv = nv),
                 'exceed the number of columns')
    expect_error(project(vs_list[[1]], nc = 400000, nv = nv),
                 'exceed the number of columns')
    expect_error(project(fit_list[[1]], vind = 1:nv, ns = 400000),
                 'exceed the number of columns')
    expect_error(project(fit_list[[1]], vind = 1:nv, nc = 400000),
                 'exceed the number of columns')
  })
  
  test_that("project: specifying intercept has an expected effect", {
    for(i in 1:length(vs_list)) {
      i_inf <- names(vs_list)[i]
      p <- project(vs_list[[i]], nv = nv, intercept = TRUE)
      expect_true(p$intercept, info = i_inf)
      p <- project(vs_list[[i]], nv = nv, intercept = FALSE)
      expect_true(!p$intercept, info = i_inf)
      expect_true(all(p$alpha==0), info = i_inf)
    }
  })
  
  test_that("project: specifying the seed does not cause errors", {
    for(i in 1:length(vs_list)) {
      i_inf <- names(vs_list)[i]
      p <- project(vs_list[[i]], nv = nv, seed = seed)
      expect_named(p, c('kl', 'weights', 'dis', 'alpha', 'beta', 'vind',
                        'p_type', 'intercept', 'family_kl'),
                   ignore.order = T, info = i_inf)
    }
  })
  
  test_that("project: adding more regularization has an expected effect", {
      regul <- c(1e-6, 1e-3, 1e-1, 1e1, 1e4)
      for(i in 1:length(vs_list)) {
          #i_inf <- names(vs_list)[i]
          norms <- rep(0, length(regul))
          for (j in 1:length(regul))
              norms[j] <- sum(project(vs_list[[i]], nv = 3, seed = seed, nc=1, regul=regul[j])$beta^2)
          for (j in 1:(length(regul)-1))
              expect_gt(norms[j],norms[j+1])
      }
  })
  
  test_that("project: projecting full model onto itself does not change results", {
  	
  	tol <- 1e-3
    
    for (i in 1:length(fit_list)) {
      fit <- fit_list[[i]]
      draws <- as.data.frame(fit)
      alpha_ref <- draws$`(Intercept)`
      beta_ref <- draws[,1+(1:nv),drop=F]
      S <- nrow(draws)
      proj <- project(fit, vind = 1:nv, seed = seed, ns=S, regul=0)
      
      # test alpha and beta
      dalpha <- max(abs(proj$alpha - alpha_ref))
      dbeta <- max(abs(proj$beta - t(beta_ref)))
      expect_lt(dalpha, tol)
      expect_lt(dbeta, tol)
      
      if (ncol(draws) > nv+1) {
        # test dispersion
        dis_ref <- draws[,ncol(draws)]
        ddis <- max(abs(proj$dis - draws$sigma))
        expect_lt(ddis, tol)
      }
    }
  })
  
  test_that("project: works as expected from a cvsel object", {
    SW({
    cvs <- cv_varsel(fit_gauss, nv_max = 3, verbose = FALSE)
    p <- project(cvs, nv=3)
    })
    expect_equal(nrow(p$beta), 3)
    expect_length(p$vind, 3)
  })

}