# Skips -----------------------------------------------------------------------

# Skip a test that only loops over the test setups in `x` (for example
# `refmods`) if `x` is empty. On CRAN, the setup fits no reference models, so
# such a test would run no expectation; the skip reports the reason.
skip_if_no_setups <- function(x) {
  skip_if(length(x) == 0,
          paste0("`", deparse(substitute(x)), "` is empty (no test setups)"))
}
