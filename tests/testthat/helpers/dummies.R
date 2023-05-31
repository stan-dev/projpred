dummy_glob <- "REMOVE THIS DUMMY ENTRY"

# Remove undesired dummy entries from a list.
#
# @param x A list.
# @param dummy_txt The text that identifies a dummy entry.
#
# @return `x` with dummy elements removed.
rm_dummies <- function(x, dummy_txt = dummy_glob, return_logical = FALSE) {
  keep_this <- !sapply(x, identical, dummy_txt)
  if (return_logical) {
    return(keep_this)
  }
  return(x[keep_this])
}
