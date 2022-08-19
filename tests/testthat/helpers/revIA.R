# A function for reversing the order of the individual terms in ":" interaction
# terms:
revIA <- function(trms) {
  trms_split <- strsplit(grep(":", trms, value = TRUE), ":")
  return(unlist(lapply(trms_split, function(trm_split) {
    paste(rev(trm_split), collapse = ":")
  })))
}
