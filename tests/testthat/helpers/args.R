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
#
# @return A list with elements:
#   * `nprjdraws`: The value of the actually used argument, i.e., the number of
#     (resulting) projected draws.
#   * `clust_used`: A single logical value indicating whether the projected
#     draws will have nonconstant weights (`TRUE`) or not (`FALSE`).
#   * `clust_used_gt1`: A single logical value indicating whether `clust_used`
#     is `TRUE` *and* `nprjdraws` is greater than `1`.
ndr_ncl_dtls <- function(args_i) {
  ndr_ncl_nm <- intersect(names(args_i), c("ndraws", "nclusters"))
  if (length(ndr_ncl_nm) == 0) {
    ndr_ncl_nm <- "ndraws"
    nprjdraws <- ndraws_pred_default
  } else {
    stopifnot(length(ndr_ncl_nm) == 1)
    nprjdraws <- args_i[[ndr_ncl_nm]]
  }
  clust_used <- ndr_ncl_nm == "nclusters" && nprjdraws < nrefdraws
  clust_used_gt1 <- clust_used && nprjdraws > 1
  return(nlist(nprjdraws, clust_used, clust_used_gt1))
}

# A helper function for retrieving project()'s output element
# `const_wdraws_prj`, but taking into account that project() output may be a
# `proj_list`.
#
# @param prj_out Output from project().
#
# @return The same as project()'s output element `const_wdraws_prj`.
has_const_wdr_prj <- function(prj_out) {
  if (!is_proj_list(prj_out)) prj_out <- list(prj_out)
  return(prj_out[[1]]$const_wdraws_prj)
}

# A helper function for creating expectations for the number of projected draws
# in the output of proj_predict().
#
# @param args_i Passed to ndr_ncl_dtls().
# @param prj_out Passed to has_const_wdr_prj().
# @param nresample_clusters_crr The value which was used for proj_predict()'s
#   argument `nresample_clusters`.
#
# @return The number of projected draws in the output of proj_predict().
ndr_pp_out <- function(args_i, prj_out,
                       nresample_clusters_crr = nresample_clusters_default) {
  if (has_const_wdr_prj(prj_out)) {
    ndr_ncl <- ndr_ncl_dtls(args_i)
    nprjdraws_out <- ndr_ncl$nprjdraws
  } else {
    nprjdraws_out <- nresample_clusters_crr
  }
  return(nprjdraws_out)
}
