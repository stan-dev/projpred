# Extends an abbreviated family name to its long form (only for categorical and
# ordinal families):
get_fam_long <- function(fam_nm) {
  switch(fam_nm,
         brnll = "bernoulli",
         categ = "categorical",
         cumul = "cumulative",
         srtio = "sratio",
         crtio = "cratio",
         adcat = "acat",
         NA_character_)
}

# As get_fam_long(), but including the families from the traditional projection:
get_fam_long_full <- function(fam_nm) {
  switch(fam_nm,
         "gauss" = "gaussian",
         "binom" = "binomial",
         "poiss" = "poisson",
         get_fam_long(fam_nm))
}

get_f_cumul <- function(link_nm = link_str) {
  structure(list(family = "cumulative_rstanarm",
                 link = link_nm),
            class = "family")
}

# Standardize the left-hand side of a formula, i.e., get a slightly modified
# formula as well as the response variable names before and after any
# evaluations of special expressions:
stdize_lhs <- function(formul_crr) {
  formul_crr <- rm_cbind(formul_crr)
  formul_crr <- rm_addresp(formul_crr)
  y_nm_orig <- as.character(formul_crr)[2]
  y_nm <- gsub("\\(|\\)", "", y_nm_orig)
  return(nlist(fml = formul_crr, y_nm_orig, y_nm))
}

# A function to retrieve the formula from a fit (`fit_obj`):
get_formul_from_fit <- function(fit_obj) {
  formul_out <- formula(fit_obj)
  if (inherits(fit_obj, "brmsfit")) {
    formul_out <- formula(formul_out)
  }
  return(formul_out)
}

# A function to adapt a given dataset (`dat`) appropriately to a given formula
# (`formul_crr`):
get_dat_formul <- function(formul_crr, needs_adj, dat_crr = dat,
                           add_offs_dummy = FALSE) {
  if (needs_adj) {
    stdized_lhs <- stdize_lhs(formul_crr)
    dat_crr[[stdized_lhs$y_nm]] <- eval(str2lang(stdized_lhs$y_nm_orig),
                                        dat_crr)
  }
  if (add_offs_dummy) {
    # Needed in some latent projection cases where posterior_linpred.stanreg()
    # complains about `offs_col` not found, even though its argument `offset` is
    # used. The specific value doesn't matter:
    dat_crr$offs_col <- 42
  }
  return(dat_crr)
}

# A function to adapt a given dataset (`dat`) appropriately to a given "test
# setup" (`tstsetup`):
get_dat <- function(tstsetup, dat_crr = dat, offs_ylat = 0, ...) {
  dat_crr <- get_dat_formul(
    args_fit[[args_prj[[tstsetup]]$tstsetup_fit]]$formula,
    needs_adj = grepl("\\.spclformul", tstsetup), dat_crr = dat_crr, ...
  )
  if (args_prj[[tstsetup]]$prj_nm == "latent") {
    if (args_prj[[tstsetup]]$pkg_nm == "rstanarm" &&
        grepl("\\.with_offs\\.", tstsetup)) {
      dat_crr$projpred_internal_offs_stanreg <- offs_ylat
    }
    y_nm <- stdize_lhs(prjs[[tstsetup]]$refmodel$formula)$y_nm
    # Use `ref_predfun_usr` here (instead of `ref_predfun`) to include
    # offsets:
    refprd_with_offs <- get(
      "ref_predfun_usr",
      envir = environment(prjs[[tstsetup]]$refmodel$ref_predfun)
    )
    dat_crr[[y_nm]] <- rowMeans(unname(
      refprd_with_offs(fit = prjs[[tstsetup]]$refmodel$fit, newdata = dat_crr)
    ))
  }
  return(dat_crr)
}

# A function to get the elements which may be supplied to argument `penalty`:
get_penal_possbl <- function(formul_crr) {
  return(setdiff(colnames(model.matrix(formul_crr, data = dat)), "(Intercept)"))
}
