context("elnet")

# tests for glm_elnet

if (!requireNamespace("glmnet", quietly = TRUE)) {
  stop("glmnet needed for this function to work. Please install it.",
    call. = FALSE
  )
}

set.seed(1235)
n <- 40
nterms <- 10
nterms_fit <- nterms - 5
x <- matrix(rnorm(n * nterms, 0, 1), n, nterms)
b <- c(seq(0, 1, length.out = nterms_fit), rep(0, nterms - nterms_fit)) # runif(nterms)-0.5
dis <- runif(1, 0.3, 0.5) # runif(1, 1, 2)
x_tr <- x[, 1:nterms_fit]
b_tr <- b[1:nterms_fit]
weights <- sample(1:4, n, replace = TRUE)
weights_norm <- weights / sum(weights) * n
offset <- 0.1 * rnorm(n) # rnorm(n)
penalty <- runif(ncol(x_tr)) + 0.5
## must scale the penalties to be comparable to glmnet
penalty <- penalty / sum(penalty) * ncol(x_tr)

tol <- 1e-04
extra_thresh <- 1e-10

test_that(paste("glm_elnet: various families and setups, glm_elnet and glmnet",
                "should give same result"), {
  fams <- list(gaussian(), binomial(), poisson())
  x_list <- lapply(fams, function(fam) x_tr)
  y_list <- lapply(fams, function(fam) {
    if (fam$family == "gaussian") {
      y <- rnorm(n, x_tr %*% b_tr, 0.1)
      y_glmnet <- y
    } else if (fam$family == "binomial") {
      y <- rbinom(n, weights, fam$linkinv(x_tr %*% b_tr))
      y <- y / weights
      ## different way of specifying binomial y for glmnet
      y_glmnet <- cbind(1 - y, y)
    } else if (fam$family == "poisson") {
      y <- rpois(n, fam$linkinv(x_tr %*% b_tr))
      y_glmnet <- y
    }
    list(y = y, y_glmnet = y_glmnet)
  })

  for (intercept in c(FALSE, TRUE)) {
    ## cannot test 0 < alpha < 1 because it seems glmnet uses the 'unnaive'
    ## elastic net
    for (alpha in c(0, 1)) {
      ## it seems glmnet does the normalization differently so we can't test
      ## normalize TRUE
      for (normalize in c(FALSE)) {
        for (use_offset in c(FALSE, TRUE)) {
          for (use_weights in c(FALSE, TRUE)) {
            for (i in seq_along(fams)) {
              x <- x_list[[i]]
              y <- y_list[[i]]$y
              y_glmnet <- y_list[[i]]$y_glmnet
              fam <- fams[[i]]
              nlam <- 500
              lambda_min_ratio <- 1e-3

              if (use_offset) {
                os <- offset
              } else {
                os <- rep(0, n)
              }
              if (use_weights) {
                w <- weights
              } else {
                w <- rep(1, n)
              }

              # compute the whole solution paths
              fit1 <- glm_elnet(x, y, fam,
                alpha = alpha,
                lambda_min_ratio = 0.1 * lambda_min_ratio, nlambda = nlam,
                weights = w, offset = os,
                normalize = normalize, thresh = 1e-12, intercept = intercept
              )
              fit2 <- glmnet::glmnet(x, y_glmnet,
                family = fam$family, alpha = alpha,
                lambda.min.ratio = lambda_min_ratio, nlambda = nlam,
                weights = w, offset = os, standardize = normalize,
                thresh = 1e-12, intercept = intercept
              )
              ## check that with a given L1-norm, the coefficient values are the
              ## same (need to check it this way since the lambda values are not
              ## comparable between glm_elnet and glmnet)
              b1 <- rbind(fit1$beta0, fit1$beta)
              enorm1 <- colSums(abs(b1))
              b2 <- as.matrix(rbind(fit2$a0, fit2$beta))
              enorm2 <- colSums(abs(b2))
              for (j in 1:nrow(b1)) {
                # loop through each coefficient (including intercept)
                infostr <- paste0(
                  "alpha = ", alpha, ", normalization = ", normalize,
                  ", weights = ", use_weights, ", offset = ", use_offset,
                  ", family = ", fam$family
                )
                b1j_interp <- approxfun(enorm1, b1[j, ])
                magn <- max(abs(b2[j, ]) + 1e-9)
                max_rel_diff <- max(abs(b1j_interp(enorm2) - b2[j, ]) / magn,
                                    na.rm = TRUE)
                expect_true(max_rel_diff < 2 * 1e-2, info = infostr)
              }

              ## plot coefficient path of some variable (useful when debugging
              ## so we leave this here commented)
              ## ds <- 0.3
              ## j <- 1
              ## ggplot() +
              ##   geom_point(aes(x=enorm1, y=b1[j,]), color='black', size=ds) +
              ##   geom_line(aes(x=enorm1, y=b1[j,]), color='black') +
              ##   geom_point(aes(x=enorm2, y=b2[j,]), color='red', size=ds) +
              ##   geom_line(aes(x=enorm2, y=b2[j,]), color='red')
              ##    + xlim(0.5,1) + ylim(0,0.1)
            }
          }
        }
      }
    }
  }
})

test_that(paste("glm_elnet: poisson, log-link, normalization should not affect",
                "the maximum likelihood solution"), {
  fam <- extend_family(poisson(link = "log"))
  y <- rpois(n, fam$linkinv(x %*% b))

  nlam <- 100
  elnetfit1 <- glm_elnet(x_tr, y,
    family = fam, nlambda = nlam, lambda_min_ratio = 1e-7,
    offset = offset, weights = weights_norm,
    intercept = TRUE, normalize = FALSE
  )
  elnetfit2 <- glm_elnet(x_tr, y,
    family = fam, nlambda = nlam, lambda_min_ratio = 1e-7,
    offset = offset, weights = weights_norm,
    intercept = TRUE, normalize = TRUE
  )

  expect_equal(c(elnetfit1$beta0[nlam], elnetfit1$beta[, nlam]),
    c(elnetfit2$beta0[nlam], elnetfit2$beta[, nlam]),
    tolerance = tol
  )
})

test_that("glm_elnet with alpha=0 and glm_ridge give the same result.", {
  for (famstr in c("gaussian", "binomial", "poisson")) {
    if (famstr == "gaussian") {
      fam <- gaussian(link = "identity")
      y <- rnorm(n, x_tr %*% b_tr, 0.5)
    } else if (famstr == "binomial") {
      fam <- binomial(link = "probit")
      y <- rbinom(n, weights, fam$linkinv(x_tr %*% b_tr)) / weights
    } else if (famstr == "poisson") {
      fam <- poisson(link = "log")
      y <- rpois(n, fam$linkinv(x_tr %*% b_tr))
    }

    for (intercept in c(TRUE, FALSE)) {
      for (normalize in c(TRUE, FALSE)) {

        # compute the L2-path with glm_elnet
        elnetfit <- glm_elnet(x_tr, y,
          family = fam, nlambda = 50, alpha = 0,
          offset = offset, weights = weights, penalty = penalty,
          intercept = intercept, normalize = normalize, thresh = 1e-15
        )
        b1 <- rbind(elnetfit$beta0, elnetfit$beta)

        # compute the solutions using glm_ridge in the same lambda grid
        b2 <- array(dim = dim(b1))
        for (j in seq_along(elnetfit$lambda)) {
          lam <- elnetfit$lambda[j]
          ridgefit <- glm_ridge(x_tr, y,
            family = fam, lambda = lam,
            offset = offset, weights = weights, penalty = penalty,
            intercept = intercept, normalize = normalize, thresh = 1e-15
          )
          b2[1, j] <- ridgefit$beta0
          b2[2:nrow(b2), j] <- ridgefit$beta
        }

        infostr <- paste0(
          "intercept = ", intercept,
          ", normalize = ", normalize, ", family = ", fam$family
        )
        expect_true(max(abs(b1 - b2)) < 1e-6, info = infostr)
      }
    }
  }
})
