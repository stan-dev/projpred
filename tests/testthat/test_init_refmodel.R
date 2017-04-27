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

SW(
    fit_gauss <- stan_glm(y ~ x, family = f_gauss, data = df_gauss, QR = T,
                          weights = weights, offset = offset,
                          chains = chains, seed = seed, iter = iter)
)
e <- extract(fit_gauss$stanfit)
mu <- t(posterior_linpred(fit_gauss, newdata = df_gauss, transform = T, offset=offset))
dis <- e$aux
ref_gauss <- init_refmodel(x,df_gauss$y,gaussian(),mu=mu,dis=dis,offset=offset,wobs=weights,loglik=log_lik(fit_gauss))

SW(
    fit_binom <- stan_glm(cbind(y, weights-y) ~ x, family = f_binom, QR = T,
                          data = df_binom, weights = weights, offset = offset,
                          chains = chains, seed = seed, iter = iter)
)
e <- extract(fit_binom$stanfit)
mu <- t(posterior_linpred(fit_binom, newdata = df_binom, transform = T, offset=offset))
ref_binom <- init_refmodel(x,df_binom$y/weights,binomial(),mu=mu,offset=offset,wobs=weights,loglik=log_lik(fit_binom))

SW(
    fit_poiss <- stan_glm(y ~ x, family = f_poiss, data = df_poiss, QR = T,
                          weights = weights, offset = offset,
                          chains = chains, seed = seed, iter = iter)
)
e <- extract(fit_poiss$stanfit)
mu <- t(posterior_linpred(fit_poiss, newdata = df_poiss, transform = T, offset=offset))
ref_poiss <- init_refmodel(x,df_poiss$y,poisson(),mu=mu,offset=offset,wobs=weights,loglik=log_lik(fit_poiss))

fit_list <- list(gauss = fit_gauss, binom = fit_binom, poiss = fit_poiss)
ref_list <- list(gauss = ref_gauss, binom = ref_binom, poiss = ref_poiss)
vs_list <- lapply(fit_list, varsel, nv_max = nv, verbose = FALSE)
vs2_list <- lapply(ref_list, varsel, nv_max = nv, verbose = FALSE)
SW({
cvvs_list <- lapply(fit_list, cv_varsel, nv_max = nv, verbose = FALSE)
cvvs2_list <- lapply(ref_list, cv_varsel, nv_max = nv, verbose = FALSE)
})
proj_vind_list <- lapply(vs_list, project, vind = c(2,3))
proj2_vind_list <- lapply(vs2_list, project, vind = c(2,3))
pred_list <- lapply(vs_list, proj_linpred, xnew=x, offsetnew=offset, weightsnew=weights, nv=3)
pred2_list <- lapply(vs2_list, proj_linpred, xnew=x, offsetnew=offset, weightsnew=weights, nv=3)


context('init_refmodel')
test_that("output of varsel is the same when using stanfit and init_refmodel", {
    for(i in 1:length(vs_list)) {
        expect_equal(vs_list[[i]]$varsel, vs2_list[[i]]$varsel)
    }
})

test_that("output of cv_varsel is the same when using stanfit and init_refmodel", {
    for(i in 1:length(vs_list)) {
        expect_equal(cvvs_list[[i]]$varsel, cvvs2_list[[i]]$varsel)
    }
})

test_that("output of project is the same when using stanfit and init_refmodel", {
    for(i in 1:length(vs_list)) {
        expect_equal(proj_vind_list[[i]]$varsel, proj2_vind_list[[i]]$varsel)
    }
})

test_that("output of proj_linpred is the same when using stanfit and init_refmodel", {
    for(i in 1:length(vs_list)) {
        expect_equal(pred_list[[i]]$varsel, pred2_list[[i]]$varsel)
    }
})

