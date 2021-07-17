# A helper function for retrieving details about the actually used `ndraws` or
# `nclusters` argument as well as about associated objects
#
# @param args_i A list of arguments supplied to a function with arguments
#   `ndraws` and `nclusters` (project(), varsel(), etc.).
#
# @return A list with elements:
#   * `ndr_ncl_nm`: The name of the actually used argument (`ndraws` or
#     `nclusters`).
#   * `nprjdraws`: The value of the actually used argument.
#   * `p_type`: The corresponding expected element `p_type` of a `"projection"`
#     object.
ndr_ncl_dtls <- function(args_i) {
  ndr_ncl_nm <- intersect(names(args_i), c("ndraws", "nclusters"))
  if (length(ndr_ncl_nm) == 0) {
    ndr_ncl_nm <- "ndraws"
    nprjdraws <- ndraws_pred_default
  } else {
    stopifnot(length(ndr_ncl_nm) == 1)
    nprjdraws <- args_i[[ndr_ncl_nm]]
  }
  return(nlist(
    ndr_ncl_nm,
    nprjdraws,
    p_type = (ndr_ncl_nm == "nclusters" || nprjdraws <= 20)
  ))
}
