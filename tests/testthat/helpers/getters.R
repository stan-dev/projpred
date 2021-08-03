get_dat_formul <- function(formul_crr, needs_adj, dat_crr = dat) {
  if (needs_adj) {
    y_nm_orig <- as.character(formul_crr)[2]
    y_nm <- gsub("\\(|\\)", "", y_nm_orig)
    dat_crr[[y_nm]] <- eval(str2lang(y_nm_orig), dat_crr)
  }
  return(dat_crr)
}

get_dat <- function(tstsetup, dat_crr = dat) {
  get_dat_formul(args_fit[[args_prj[[tstsetup]]$tstsetup_fit]]$formula,
                 needs_adj = grepl("\\.spclformul", tstsetup),
                 dat_crr = dat_crr)
}

# A function to get the possible length of the vector supplied to argument
# `penalty`:
get_penal_possbl <- function(formul_crr) {
  trms_crr <- labels(terms(formul_crr))
  if (length(trms_crr) > 0 && any(grepl("xca\\.", trms_crr))) {
    contr_crr <- lapply(
      setNames(nm = grep("xca\\.", trms_crr, value = TRUE)),
      function(x_nm) {
        # Note: As mentioned in issue #149, the reference level of a
        # categorical predictor actually should not have its own coefficient.
        # But to be able to run the tests for now, we also use `contrasts =
        # FALSE` here:
        contrasts(get(x_nm, envir = as.environment(dat)),
                  contrasts = FALSE)
      }
    )
  } else {
    contr_crr <- NULL
  }
  mm_crr <- model.matrix(update(formul_crr, . ~ . + 0),
                         data = dat,
                         contrasts.arg = contr_crr)
  return(colnames(mm_crr))
}
