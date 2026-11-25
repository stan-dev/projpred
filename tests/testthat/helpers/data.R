# Simulated data for the tests. setup.R sources this file with the seed
# `seed_dat`, after it has set the switches and the model and family
# names.

## Setup ------------------------------------------------------------------

# Number of observations:
nobsv <- 39L
# Values for testing:
nobsv_tst <- c(1L, nobsv %/% 3L)

# For ordinal models (but also used for categorical models):
nthres <- 2L
ncat <- nthres + 1L
yunq_num <- seq_len(ncat)
yunq_chr <- paste0("y", yunq_num)
# The intercepts at centered predictors, also known as thresholds:
thres <- qlogis(seq_len(nthres) / ncat)
link_str <- "logit"

## GLMs --------------------------------------------------------------------
## Add population-level effects to the intercept-(and-offset-)only model

nterms_cont <- 3L
x_cont <- matrix(rnorm(nobsv * nterms_cont), nobsv, nterms_cont)
b_cont <- runif(nterms_cont, min = -0.5, max = 0.5)

nlvl_fix <- c(3L, 2L)
nlvl_fix <- setNames(nlvl_fix, seq_along(nlvl_fix))
if (length(nlvl_fix) <= 1) {
  names(nlvl_fix) <- paste0("z.", names(nlvl_fix))
}
nterms_cate <- length(nlvl_fix)
x_cate_list <- lapply(nlvl_fix, function(nlvl_fix_i) {
  x_cate <- gl(n = nlvl_fix_i, k = floor(nobsv / nlvl_fix_i), length = nobsv,
               labels = paste0("lvl", seq_len(nlvl_fix_i)))
  b_cate <- runif(nlvl_fix_i, min = -0.5, max = 0.5)
  ### Using a model.matrix() approach:
  # x_cate_mat <- model.matrix(~ 0 + x_cate)
  # eta_cate <- x_cate_mat %*% b_cate
  ###
  ### Using an indexing approach:
  eta_cate <- b_cate[x_cate]
  ###
  return(nlist(x_cate, eta_cate, b_cate))
})

# Intercept and offsets (offsets are added later):
icpt <- -0.42
offs_tst <- rnorm(nobsv)

eta_glm <- icpt +
  x_cont %*% b_cont +
  do.call("+", lapply(x_cate_list, "[[", "eta_cate"))

nterms_glm <- nterms_cont + nterms_cate

## GLMMs ------------------------------------------------------------------
## Add group-level effects to the GLMs, yielding GLMMs

nlvl_ran <- c(6L)
nlvl_ran <- setNames(nlvl_ran, seq_along(nlvl_ran))
if (length(nlvl_ran) <= 1) {
  names(nlvl_ran) <- paste0("z.", names(nlvl_ran))
}
# Multiply by 2 because of the random intercept as well as the random slope (the
# latter for `x_cont[, 1]`):
nterms_z <- length(nlvl_ran) * 2L
z_list <- lapply(nlvl_ran, function(nlvl_ran_i) {
  z <- gl(n = nlvl_ran_i, k = floor(nobsv / nlvl_ran_i), length = nobsv,
          labels = paste0("lvl", seq_len(nlvl_ran_i)))
  r_icpts <- rnorm(nlvl_ran_i, sd = 0.4)
  r_xco1 <- rnorm(nlvl_ran_i, sd = 0.4)
  eta_z <- r_icpts[z] + r_xco1[z] * x_cont[, 1]
  return(nlist(z, eta_z, r_icpts, r_xco1))
})
eta_glmm <- eta_glm + do.call("+", lapply(z_list, "[[", "eta_z"))

nterms_glmm <- nterms_glm + nterms_z

## GAMs -------------------------------------------------------------------
## Add nonlinear (smoothed) effects to the GLMs

# For simplicity, always use the same nonlinear function (could be extended in
# the future):
s_mat <- apply(x_cont[, 1, drop = FALSE], 2, function(x, b = 2, c = - pi / 4) {
  b * sin(x - c)
})
s_sum <- rowSums(s_mat)
nterms_s <- ncol(s_mat)
eta_gam <- eta_glm + s_sum

# Multiply by 2 because of the baseline linear term as well as the standard
# deviation for the wiggliness around it:
nterms_gam <- nterms_glm + 2L * nterms_s

## GAMMs ------------------------------------------------------------------
## Add nonlinear (smoothed) effects to the GLMMs

eta_gamm <- eta_glmm + s_sum

nterms_gamm <- nterms_glmm + 2L * nterms_s

## Combined dataset -------------------------------------------------------

f_gauss <- gaussian()
f_binom <- f_brnll <- binomial()
f_poiss <- poisson()
dis_tst <- runif(1L, 1, 2)
wobs_tst <- sample(1:4, nobsv, replace = TRUE)
# For the "brnll" `fam_nm`, offsets are simply not added to have some
# scenarios without offsets.
# For GAMs, offsets are not added because of rstanarm issue #546 (see
# also further below).
# For GAMMs, offsets are not added because of rstanarm issue #253 (see
# also further below).
# (The brms "gam" and "gamm" cases are handled in the same way as the rstanarm
# "gam" and "gamm" cases to avoid too many special cases.)
# For the "categ" `fam_nm`, offsets are not added because they are currently not
# supported for it.
offs_expr <- expression(!(fam_nm %in% c("brnll", "categ") ||
                            mod_nm %in% c("gam", "gamm")))
cre_dat <- function(idxs_crr, offs_crr, wobs_crr, dis_crr) {
  nobsv_crr <- length(idxs_crr)
  dat_crr <- lapply(mod_nms, function(mod_nm) {
    lapply(fam_nms, function(fam_nm) {
      pred_link <- get(paste0("eta_", mod_nm))
      pred_link <- pred_link[idxs_crr, , drop = FALSE]
      if (eval(offs_expr)) {
        pred_link <- pred_link + offs_crr
      }
      if (fam_nm %in% fam_nms_ordin) {
        pred_link <- pred_link - icpt
        if (eval(offs_expr)) {
          # The equal-probability thresholds defined above refer to the state
          # before offsets are added, so we need to subtract them here in the
          # data-generating model:
          pred_link <- pred_link - offs_crr
        }
        thres_eta <- sapply(thres, function(thres_k) {
          thres_k - pred_link
        })
      } else if (fam_nm %in% fam_nms_categ) {
        pred_link <- sapply(thres, function(thres_k) {
          thres_k + pred_link
        })
      }
      if (fam_nm == "cumul") {
        pred_resp <- augdat_ilink_cumul(thres_eta, link = link_str)
      } else if (fam_nm %in% fam_nms_ordin) {
        if (fam_nm %in% c("crtio", "adcat")) {
          thres_eta <- -thres_eta
        }
        ilink_crr <- get(paste0("inv_link_", get_fam_long(fam_nm)),
                         asNamespace("brms"), mode = "function",
                         inherits = FALSE)
        pred_resp <- ilink_crr(thres_eta, link = link_str)
      } else if (fam_nm %in% fam_nms_categ) {
        ilink_crr <- get(paste0("inv_link_", get_fam_long(fam_nm)),
                         asNamespace("brms"), mode = "function",
                         inherits = FALSE)
        pred_resp <- ilink_crr(pred_link)
      } else {
        pred_resp <- get(paste0("f_", fam_nm))$linkinv(pred_link)
      }
      if (fam_nm == "gauss") {
        return(rnorm(nobsv_crr, mean = pred_resp, sd = dis_crr))
      } else if (fam_nm == "brnll") {
        return(rbinom(nobsv_crr, 1, pred_resp))
      } else if (fam_nm == "binom") {
        return(rbinom(nobsv_crr, wobs_crr, pred_resp))
      } else if (fam_nm == "poiss") {
        return(rpois(nobsv_crr, pred_resp))
      } else if (fam_nm %in% fam_nms_aug) {
        ryunq <- sapply(seq_len(nobsv_crr), function(i_obs) {
          sample(yunq_num, size = 1L, prob = pred_resp[i_obs, ])
        })
        if (use_fac) {
          ryunq <- factor(ryunq, levels = yunq_num, labels = yunq_chr,
                          ordered = fam_nm %in% fam_nms_ordin)
        }
        return(ryunq)
      } else {
        stop("Unknown `fam_nm`.")
      }
    })
  })
  dat_crr <- unlist(dat_crr, recursive = FALSE)
  names(dat_crr) <- paste("y", gsub("\\.", "_", names(dat_crr)), sep = "_")
  return(dat_crr)
}
dat <- cre_dat(idxs_crr = seq_len(nobsv), offs_crr = offs_tst,
               wobs_crr = wobs_tst, dis_crr = dis_tst)
dat <- data.frame(
  dat,
  xco = x_cont, xca = lapply(x_cate_list, "[[", "x_cate"),
  z = lapply(z_list, "[[", "z"),
  s = s_mat,
  wobs_col = wobs_tst, offs_col = offs_tst,
  check.names = FALSE
)
if (ncol(s_mat) == 1) {
  names(dat)[names(dat) == "s"] <- "s.1"
}

## Modified datasets ------------------------------------------------------

dat_wobs_ones <- within(dat, {
  wobs_col <- NULL
  wobs_col_ones <- rep_len(1, length.out = nobsv)
})
dat_wobs_new <- within(dat, {
  wobs_col <- NULL
  wobs_col_new <- rep_len(2:5, length.out = nobsv)
})

dat_offs_zeros <- within(dat, {
  offs_col <- NULL
  offs_col_zeros <- rep_len(0, length.out = nobsv)
})
dat_offs_new <- within(dat, {
  offs_col <- NULL
  offs_col_new <- seq(-2, 2, length.out = nobsv)
})

nobsv_indep <- tail(nobsv_tst, 1)
dis_indep <- runif(1L, 1, 2)
offs_indep <- rnorm(nobsv_indep)
wobs_indep <- sample(1:4, nobsv_indep, replace = TRUE)
idxs_indep <- sample.int(nobsv, size = nobsv_indep, replace = TRUE)
dat_indep <- cre_dat(idxs_crr = idxs_indep, offs_crr = offs_indep,
                     wobs_crr = wobs_indep, dis_crr = dis_indep)
dat_indep <- cbind(
  as.data.frame(dat_indep),
  dat[idxs_indep,
      grep("^y_", names(dat), value = TRUE, invert = TRUE),
      drop = FALSE]
)
dat_indep$wobs_col <- wobs_indep
dat_indep$offs_col <- offs_indep
