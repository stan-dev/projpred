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
get_dat_formul <- function(formul_crr, needs_adj, dat_crr = dat) {
  if (needs_adj) {
    stdized_lhs <- stdize_lhs(formul_crr)
    dat_crr[[stdized_lhs$y_nm]] <- eval(str2lang(stdized_lhs$y_nm_orig),
                                        dat_crr)
  }
  return(dat_crr)
}

# A function to adapt a given dataset (`dat`) appropriately to a given "test
# setup" (`tstsetup`):
get_dat <- function(tstsetup, dat_crr = dat) {
  get_dat_formul(args_fit[[args_prj[[tstsetup]]$tstsetup_fit]]$formula,
                 needs_adj = grepl("\\.spclformul", tstsetup),
                 dat_crr = dat_crr)
}

# A function to get the elements which may be supplied to argument `penalty`:
get_penal_possbl <- function(formul_crr) {
  return(setdiff(colnames(model.matrix(formul_crr, data = dat)), "(Intercept)"))
}
