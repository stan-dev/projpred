get_dat <- function(tstsetup, dat_crr = dat) {
  if (grepl("\\.spclformul", tstsetup)) {
    y_nm_orig <- as.character(
      args_fit[[args_prj[[tstsetup]]$tstsetup_fit]]$formula
    )[2]
    y_nm <- gsub("\\(|\\)", "", y_nm_orig)
    dat_crr[[y_nm]] <- eval(str2lang(y_nm_orig), dat_crr)
  }
  return(dat_crr)
}
