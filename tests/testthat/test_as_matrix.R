context("as.matrix.projection")

# Skip GAMs:
args_gam <- grep("^gam\\.", names(args_prj))
if (length(args_gam)) {
  warning("Skipping GAMs because of issue #151. Note that for GAMs, the ",
          "current expectations in `test_as_matrix.R` refer to a mixture of ",
          "brms's and rstanarm's naming scheme; as soon as issue #152 ",
          "is solved, these expectations need to be adopted.")
  args_prj <- args_prj[-args_gam]
}
# Skip GAMMs:
args_gamm <- grep("^gamm\\.", names(args_prj))
if (length(args_gamm)) {
  warning("Skipping GAMMs because of issue #131.")
  args_prj <- args_prj[-args_gamm]
}

test_that("as.matrix.projection() works", {
  for (tstsetup in names(prjs)) {
    tstsetup_ref <- args_prj[[tstsetup]]$tstsetup_ref
    mod_crr <- args_prj[[tstsetup]]$mod_nm
    fam_crr <- args_prj[[tstsetup]]$fam_nm
    solterms <- args_prj[[tstsetup]]$solution_terms
    ndr_ncl <- ndr_ncl_dtls(args_prj[[tstsetup]])

    # Expected warning (more precisely: regexp which is matched against the
    # warning; NA means no warning) for as.matrix.projection():
    if (ndr_ncl$clust_used) {
      # Clustered projection, so we expect a warning:
      warn_prjmat_expect <- "the clusters might have different weights"
    } else {
      warn_prjmat_expect <- NA
    }
    expect_warning(m <- as.matrix(prjs[[tstsetup]]),
                   warn_prjmat_expect, info = tstsetup)

    if (fam_crr == "gauss") {
      npars_fam <- "sigma"
    } else if (fam_crr == "binom") {
      npars_fam <- character()
    }

    colnms_prjmat_expect <- c(
      "Intercept",
      grep("\\|", grep("x(co|ca)\\.[[:digit:]]", solterms, value = TRUE),
           value = TRUE, invert = TRUE)
    )
    xca_idxs <- as.integer(
      sub("^xca\\.", "",
          grep("^xca\\.", colnms_prjmat_expect, value = TRUE))
    )
    for (xca_idx in xca_idxs) {
      colnms_prjmat_expect <- grep(paste0("^xca\\.", xca_idx, "$"),
                                   colnms_prjmat_expect,
                                   value = TRUE, invert = TRUE)
      colnms_prjmat_expect <- c(
        colnms_prjmat_expect,
        paste0("xca.", xca_idx, "lvl", seq_len(nlvl_fix[xca_idx])[-1])
      )
    }
    poly_trms <- grep("poly\\(.*\\)", colnms_prjmat_expect, value = TRUE)
    if (length(poly_trms) > 0) {
      poly_degree <- sub(".*(poly\\(.*\\)).*", "\\1", poly_trms)
      if (length(unique(poly_degree)) != 1) {
        stop("This test needs to be adopted. Info: ", tstsetup)
      }
      poly_degree <- unique(poly_degree)
      poly_degree <- sub(".*,[[:blank:]]+(.*)\\)", "\\1", poly_degree)
      poly_degree <- as.integer(poly_degree)
      colnms_prjmat_expect <- c(
        setdiff(colnms_prjmat_expect, poly_trms),
        unlist(lapply(poly_trms, function(poly_trms_i) {
          paste0(poly_trms_i, seq_len(poly_degree))
        }))
      )
    }
    colnms_prjmat_expect <- paste0("b_", colnms_prjmat_expect)
    if ("(1 | z.1)" %in% solterms) {
      colnms_prjmat_expect <- c(colnms_prjmat_expect, "sd_z.1__Intercept")
      colnms_prjmat_expect <- c(
        colnms_prjmat_expect,
        paste0("r_z.1[lvl", seq_len(nlvl_ran[1]), ",Intercept]")
      )
    }
    if ("(xco.1 | z.1)" %in% solterms) {
      colnms_prjmat_expect <- c(colnms_prjmat_expect, "sd_z.1__xco.1")
      colnms_prjmat_expect <- c(
        colnms_prjmat_expect,
        paste0("r_z.1[lvl", seq_len(nlvl_ran[1]), ",xco.1]")
      )
    }
    if (all(c("(1 | z.1)", "(xco.1 | z.1)") %in% solterms)) {
      colnms_prjmat_expect <- c(colnms_prjmat_expect,
                                "cor_z.1__Intercept__xco.1")
    }
    s_nms <- sub("\\)$", "",
                 sub("^s\\(", "",
                     grep("^s\\(.*\\)$", solterms, value = TRUE)))
    if (length(s_nms) > 0) {
      stopifnot(inherits(refmods[[tstsetup_ref]]$fit, "stanreg"))
      # Get the number of basis coefficients:
      s_info <- refmods[[tstsetup_ref]]$fit$jam$smooth
      s_terms <- sapply(s_info, "[[", "term")
      s_dfs <- setNames(sapply(s_info, "[[", "df"), s_terms)
      ### Alternative:
      # par_nms_orig <- colnames(
      #   as.matrix(refmods[[tstsetup_ref]]$fit)
      # )
      # s_dfs <- sapply(s_nms, function(s_nm) {
      #   sum(grepl(paste0("^s\\(", s_nm, "\\)"), par_nms_orig))
      # })
      ###
      # Construct the expected column names for the basis coefficients:
      for (s_nm in s_nms) {
        colnms_prjmat_expect <- c(
          colnms_prjmat_expect,
          paste0("b_s(", s_nm, ").", seq_len(s_dfs[s_nm]))
        )
      }
      # Needed for the names of the `smooth_sd` parameters:
      s_nsds <- setNames(
        lapply(lapply(s_info, "[[", "sp"), names),
        s_terms
      )
      # Construct the expected column names for the SDs of the smoothing
      # terms:
      for (s_nm in s_nms) {
        colnms_prjmat_expect <- c(
          colnms_prjmat_expect,
          paste0("smooth_sd[", s_nsds[[s_nm]], "]")
        )
      }
    }
    colnms_prjmat_expect <- c(colnms_prjmat_expect, npars_fam)

    expect_identical(dim(m), c(ndr_ncl$nprjdraws, length(colnms_prjmat_expect)),
                     info = tstsetup)
    ### expect_setequal() does not have argument `info`:
    # expect_setequal(colnames(m), colnms_prjmat_expect)
    expect_true(setequal(colnames(m), colnms_prjmat_expect),
                info = tstsetup)
    ###
  }
})
