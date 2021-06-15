# A custom unlist() wrapper
#
# @param x A list, at a deeper level containing an unnamed sublist or a named
#   sublist containing the name given in `nm_stop`.
# @param nm_stop A character of length 1, giving the name which signals to stop
#   the recursion.
#
# @return The recursively "unlisted" `x` (recursively up to the point where the
#   next sublist is unnamed or contains an element with name given in
#   `nm_stop`).
#
unlist_cust <- function(x, nm_stop = "mod_nm") {
  stopifnot(is.list(x) && length(x) > 0)
  if (is.null(names(x[[1]])) || nm_stop %in% names(x[[1]])) {
    return(x)
  } else {
    return(unlist_cust(unlist(x, recursive = FALSE)))
  }
}
