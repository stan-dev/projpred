# A custom unlist() wrapper
#
# @param x A list, at a deeper level containing an unnamed sublist or a sublist
#   with element "mod_nm".
#
# @return The recursively "unlisted" `x` (recursively up to the point where the
#   next sublist is unnamed or contains an element called "mod_nm").
#
unlist_cust <- function(x) {
  stopifnot(is.list(x) && length(x) > 0)
  if (is.null(names(x[[1]])) || "mod_nm" %in% names(x[[1]])) {
    return(x)
  } else {
    return(unlist_cust(unlist(x, recursive = FALSE)))
  }
}
