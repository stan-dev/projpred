# A helper function for picking only those elements from an argument list which
# are not arguments in the narrower sense.
#
# @param args_i A list of arguments.
#
# @return Those elements of `args_i` named as follows (if existing):
#   * `"mod_nm"`,
#   * `"fam_nm"`,
#   * everything starting with "tstsetup_" (regexp: "^tstsetup_").
only_nonargs <- function(args_i) {
  nms_only <- c("mod_nm", "fam_nm", "pkg_nm", "prj_nm",
                grep("^tstsetup_", names(args_i), value = TRUE))
  return(args_i[intersect(names(args_i), nms_only)])
}

# A helper function for excluding elements from an argument list which are not
# arguments in the narrower sense.
#
# @param args_i A list of arguments.
# @param nms_excl_add A character vector of element names which should also be
#   excluded from `args_i`.
#
# @return `args_i` with the following elements excluded:
#   * `"mod_nm"`,
#   * `"fam_nm"`,
#   * everything starting with "tstsetup_" (regexp: "^tstsetup_").
excl_nonargs <- function(args_i, nms_excl_add = character()) {
  nms_excl <- c("mod_nm", "fam_nm", "pkg_nm", "prj_nm",
                grep("^tstsetup_", names(args_i), value = TRUE),
                nms_excl_add)
  return(args_i[setdiff(names(args_i), nms_excl)])
}

# A helper function for creating expectations for the actually used `ndraws` or
# `nclusters` argument of project() as well as for argument `nresample_clusters`
# of proj_predict().
#
# @param args_i A list of arguments supplied to project().
# @param nresample_clusters_crr The value of proj_predict()'s argument
#   `nresample_clusters` which should be used for calculating `nprjdraws_out`.
#
# @return A list with elements:
#   * `ndr_ncl_nm`: The name of the actually used argument (`ndraws` or
#     `nclusters`).
#   * `nprjdraws`: The value of the actually used argument.
#   * `clust_used`: A single logical value indicating whether the projected
#     draws will have nonconstant weights (`TRUE`) or not (`FALSE`).
#   * `nprjdraws_out`: The number of projected draws in the output. In contrast
#     to `nprjdraws`, this also takes proj_predict()'s argument
#     `nresample_clusters` into account.
ndr_ncl_dtls <- function(args_i,
                         nresample_clusters_crr = nresample_clusters_default) {
  ndr_ncl_nm <- intersect(names(args_i), c("ndraws", "nclusters"))
  if (length(ndr_ncl_nm) == 0) {
    ndr_ncl_nm <- "ndraws"
    nprjdraws <- ndraws_pred_default
  } else {
    stopifnot(length(ndr_ncl_nm) == 1)
    nprjdraws <- args_i[[ndr_ncl_nm]]
  }
  clust_used <- ndr_ncl_nm == "nclusters" && nprjdraws > 1
  if (clust_used) {
    nprjdraws_out <- nresample_clusters_crr
  } else {
    nprjdraws_out <- nprjdraws
  }
  return(nlist(
    ndr_ncl_nm,
    nprjdraws,
    clust_used,
    nprjdraws_out
  ))
}
