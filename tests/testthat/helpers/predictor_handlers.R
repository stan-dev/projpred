# A function for reversing the order of the individual terms in ":" interaction
# terms:
revIA <- function(trms) {
  trms_split <- strsplit(grep(":", trms, value = TRUE), ":")
  return(unlist(lapply(trms_split, function(trm_split) {
    paste(rev(trm_split), collapse = ":")
  })))
}

# Expand poly() terms, e.g., `poly(x, 2, raw = TRUE)` to
# `poly(x, 2, raw = TRUE)1` and `poly(x, 2, raw = TRUE)2`:
expand_poly <- function(trms, info_str) {
  poly_trms <- grep("poly\\(.*\\)", trms, value = TRUE)
  if (length(poly_trms)) {
    poly_degree <- sub(
      "poly\\(.*,[[:blank:]]*([[:digit:]]+)[[:blank:]]*,.*\\)", "\\1",
      poly_trms
    )
    poly_degree <- unique(poly_degree)
    if (length(poly_degree) != 1) {
      stop("This test needs to be adapted. Info: ", info_str)
    }
    poly_degree <- as.integer(poly_degree)
    trms <- c(
      setdiff(trms, poly_trms),
      unlist(lapply(poly_trms, function(poly_trms_i) {
        paste0(poly_trms_i, seq_len(poly_degree))
      }))
    )
  }
  return(trms)
}
