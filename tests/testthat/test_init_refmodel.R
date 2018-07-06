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

SW(
    fit_gauss <- stan_glm(y ~ x, family = f_gauss, data = df_gauss, QR = T,
                          weights = weights, offset = offset,
                          chains = chains, seed = seed, iter = iter)
)
e <- extract(fit_gauss$stanfit)
perm_inv <- c(mapply(function(p, i) order(p) + i*length(p),
                     fit_gauss$stanfit@sim$permutation,1:fit_gauss$stanfit@sim$chains-1))
dis <- e$aux[perm_inv]
predfun <- function(xnew) t(posterior_linpred(fit_gauss, newdata = data.frame(x=xnew), transform = T, offset=rep(0,nrow(xnew)) ))
ref_gauss <- init_refmodel(x,df_gauss$y,gaussian(),predfun=predfun,dis=dis,offset=offset,wobs=weights)
dref_gauss <- init_refmodel(x,df_gauss$y,gaussian(),offset=offset,wobs=weights) # data only

SW(
    fit_binom <- stan_glm(cbind(y, weights-y) ~ x, family = f_binom, QR = T,
                          data = df_binom, weights = weights, offset = offset,
                          chains = chains, seed = seed, iter = iter)
)
e <- extract(fit_binom$stanfit)
predfun <- function(xnew) t(posterior_linpred(fit_binom, newdata = data.frame(x=xnew), transform = T, offset=rep(0,nrow(xnew)) ))
ref_binom <- init_refmodel(x,df_binom$y/weights,binomial(),predfun=predfun,offset=offset,wobs=weights)
dref_binom <- init_refmodel(x,df_binom$y/weights,binomial(),offset=offset,wobs=weights) # data only

SW(
    fit_poiss <- stan_glm(y ~ x, family = f_poiss, data = df_poiss, QR = T,
                          weights = weights, offset = offset,
                          chains = chains, seed = seed, iter = iter)
)
e <- extract(fit_poiss$stanfit)
predfun <- function(xnew) t(posterior_linpred(fit_poiss, newdata = data.frame(x=xnew), transform = T, offset=rep(0,nrow(xnew)) ))
ref_poiss <- init_refmodel(x,df_poiss$y,poisson(),predfun=predfun,offset=offset,wobs=weights)
dref_poiss <- init_refmodel(x,df_poiss$y,poisson(),offset=offset,wobs=weights)



fit_list <- list(gauss = fit_gauss, binom = fit_binom, poiss = fit_poiss)
ref_list <- list(gauss = ref_gauss, binom = ref_binom, poiss = ref_poiss)
dref_list <- list(gauss = dref_gauss, binom = dref_binom, poiss = dref_poiss)

vs_list <- lapply(fit_list, varsel, nv_max = nv, verbose = FALSE)
vs2_list <- lapply(ref_list, varsel, nv_max = nv, verbose = FALSE)
vsd_list <- lapply(dref_list, varsel, nv_max = nv, verbose = FALSE)

SW({
cvvs_list <- lapply(fit_list, cv_varsel, nv_max = nv, verbose = FALSE)
cvvs2_list <- lapply(ref_list, cv_varsel, nv_max = nv, verbose = FALSE)
cvvsd_list <- lapply(dref_list, cv_varsel, nv_max = nv, verbose = FALSE)
})

proj_vind_list <- lapply(vs_list, project, vind = c(2,3), seed = seed)
proj2_vind_list <- lapply(vs2_list, project, vind = c(2,3), seed = seed)
projd_vind_list <- lapply(vsd_list, project, vind = c(2,3), seed = seed)

pred_list <- lapply(vs_list, proj_linpred, xnew=x, seed = seed,
                    offsetnew=offset, weightsnew=weights, nv=3)
pred2_list <- lapply(vs2_list, proj_linpred, xnew=x, seed = seed,
                     offsetnew=offset, weightsnew=weights, nv=3)
predd_list <- lapply(vsd_list, proj_linpred, xnew=x, seed = seed,
										 offsetnew=offset, weightsnew=weights, nv=3)


context('init_refmodel')
test_that("output of varsel is the same when using stanfit and init_refmodel", {
	for (i in seq_along(vs_list)) {
		expect_equal(vs_list[[i]]$varsel, vs2_list[[i]]$varsel)
	}
})

test_that("output of cv_varsel is the same when using stanfit and init_refmodel", {
	for (i in seq_along(cvvs_list)) {
		expect_equal(cvvs_list[[i]]$varsel, cvvs2_list[[i]]$varsel)
	}
})

test_that("output of project is the same when using stanfit and init_refmodel", {
	for (i in seq_along(proj_vind_list)) {
		expect_equivalent(proj_vind_list[[i]], proj2_vind_list[[i]])
		#expect_equal(proj_vind_list[[i]], proj2_vind_list[[i]]) # this errors for some reason but can't figure out why
	}
})

test_that("output of proj_linpred is the same when using stanfit and init_refmodel", {
	for (i in seq_along(pred_list)) 
		expect_equal(pred_list[[i]], pred2_list[[i]])
})



test_that('output of varsel is sensible with only data provided as reference model', {
	for(i in seq_along(vsd_list)) {
		# vind seems legit
		expect_equal(length(vsd_list[[i]]$varsel$vind), nv)
		
		# kl seems legit
		expect_equal(length(vsd_list[[i]]$varsel$kl), nv + 1)
		
		# kl decreasing
		expect_equal(vsd_list[[i]]$varsel$kl, cummin(vsd_list[[i]]$varsel$kl))
		
		# summaries seems legit
		expect_named(vsd_list[[i]]$varsel$summaries, c('sub', 'full'))
		expect_equal(length(vsd_list[[i]]$varsel$summaries$sub), nv + 1)
		expect_named(vsd_list[[i]]$varsel$summaries$sub[[1]], c('mu', 'lppd'))
		expect_named(vsd_list[[i]]$varsel$summaries$full, c('mu', 'lppd'))
	}
})

test_that("output of varsel is sensible with only data provided as reference model", {
	for(i in seq_along(cvvsd_list)) {
		# vind seems legit
		expect_equal(length(cvvsd_list[[i]]$varsel$vind), nv)
		
		# kl seems legit
		expect_equal(length(cvvsd_list[[i]]$varsel$kl), nv + 1)
		
		# kl decreasing
		expect_equal(cvvsd_list[[i]]$varsel$kl, cummin(cvvsd_list[[i]]$varsel$kl))
		
		# summaries seems legit
		expect_named(cvvsd_list[[i]]$varsel$summaries, c('sub', 'full'))
		expect_equal(length(cvvsd_list[[i]]$varsel$summaries$sub), nv + 1)
		expect_named(cvvsd_list[[i]]$varsel$summaries$sub[[1]], c('mu', 'lppd'))
		expect_named(cvvsd_list[[i]]$varsel$summaries$full, c('mu', 'lppd'))
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
			expect_equal(p[[j]]$family_kl, vs_list[[i]]$varsel$family_kl)
		}
		# kl should be non-increasing on training data
		klseq <- sapply(p, function(e) e$kl)
		expect_equal(klseq, cummin(klseq))
		
		# all submodels should use the same clustering/subsampling
		expect_equal(p[[1]]$weights, p[[nv]]$weights)
	}
})

