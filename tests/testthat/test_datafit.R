# tests for data based estimates 

set.seed(1235)
n <- 40
nv <- 5
x <- matrix(rnorm(n*nv, 0, 1), n, nv)
b <- seq(0,1,length.out = nv)
dis <- runif(1, 0.3,0.5)
weights <-  sample(1:4, n, replace = T)#
offset <- 0.1*rnorm(n)
seed <- 1235
# source(file.path('helpers', 'SW.R'))



fams <- list(gaussian(), binomial(), poisson())
x_list <- lapply(fams, function(fam) x)
y_list <- lapply(fams, function(fam) {
  if (fam$family == 'gaussian') {
    y <- rnorm(n, x%*%b, 0.5)
    y_glmnet <- y
  } else if (fam$family == 'binomial') {
    y <- rbinom(n, weights, fam$linkinv(x%*%b))
    y <- y/weights
    y_glmnet <- cbind(1-y,y) # different way of specifying binomial y for glmnet
  } else if (fam$family == 'poisson') {
    y <- rpois(n, fam$linkinv(x%*%b))
    y_glmnet <- y
  }
  list(y=y, y_glmnet=y_glmnet)
})





context('datafit')


# test_that("L1-projection with data reference gives the same results as Lasso from glmnet.", {
#   
#   for (i in seq_along(fams)) {
#     
#     x <- x_list[[i]]
#     y <- y_list[[i]]$y
#     y_glmnet <- y_list[[i]]$y_glmnet
#     fam <- fams[[i]]
#     
#     lambda_min_ratio <- 1e-4
#     nlambda <- 200
#     
#     # Lasso solution with projpred
#     ref <- init_refmodel(x,y,family = fam, wobs = weights, offset = offset)
#     vs <- varsel(ref, method='l1', lambda_min_ratio = lambda_min_ratio, nlambda = nlambda, thresh = 1e-12)
#     pred1 <- proj_linpred(vs, xnew = x, nv=0:nv, transform = F)
#     
#     # compute the results for the Lasso
#     lasso <- glmnet(x,y_glmnet,family=fam$family, weights = weights, offset = offset,
#                     lambda.min.ratio = lambda_min_ratio, nlambda = nlambda, thresh = 1e-12, dfmax = nv+5)
#     vind <- predict(lasso, type='nonzero', s=lasso$lambda)
#     nselected <- sapply(vind, function(e) length(e))
#     lambdainds <- sapply(unique(nselected), function(nv) max(which(nselected==nv)))
#     lambdaval <- lasso$lambda[lambdainds]
#     pred2 <- predict(lasso, newx=x, type='link', s=lambdaval, newoffset=rep(0,n))
#     
#     
#     k <- 3
#     qplot(pred1[[k]], pred2[,k]) + geom_abline(slope=1)
#     
#     # plot(lasso, xlim=c(0,0.3), ylim=c(0,0.2))
#     plot(lasso, xvar='lambda')
#     plot(lasso)
#     b <- vs$varsel$spath$beta
#     l1norm <- colSums(abs(b))
#     for (j in 1:nrow(b))
#       lines(l1norm, b[j,], type = 'p')
#     
#   }
# })
# 
# 
# 
# lam <- 2
# out <- projpred:::glm_elnet(x,y,projpred:::kl_helpers(fam), lambda=lam*n, weights = weights, offset=offset)
# out2 <- glmnet(x,y,family = 'poisson', lambda=lam, weights = weights, offset=offset)
# 
# 
# 
# 
# 
# # out <- projpred:::glm_elnet(x,y,projpred:::kl_helpers(fam), nlambda=100,  offset=offset)
# # out2 <- glmnet(x,y,family = 'poisson', nlambda=100,  offset=offset)
# 
# nlam <- 300
# tic(); out <- projpred:::glm_elnet(x,y,projpred:::kl_helpers(fam), nlambda=nlam, 
#                             weights = weights/sum(weights), offset=offset, normalize = F); toc()
# tic(); out2 <- glmnet(x,y,family = 'poisson', nlambda=nlam, weights = 10*weights, offset=offset, standardize = F); toc()
# 
# b1 <- out$beta
# l1norm1 <- colSums(abs(b1))
# b2 <- out2$beta
# l1norm2 <- colSums(abs(b2))
# 
# ds <- 0.3
# j <- 1
# ggplot() +
#   geom_point(aes(x=l1norm1, y=b1[j,]), color='black', size=ds) +
#   geom_line(aes(x=l1norm1, y=b1[j,]), color='black') +
#   geom_point(aes(x=l1norm2, y=b2[j,]), color='red', size=ds) +
#   geom_line(aes(x=l1norm2, y=b2[j,]), color='red') 

