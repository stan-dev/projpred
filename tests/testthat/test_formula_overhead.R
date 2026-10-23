context("formula and model-matrix shortcuts of the projections")

test_that(paste(
  "validate_response_formula() gives the same formulas as one update() per",
  "response"
), {
  rhs_list <- c(
    "x1 + x2", "x1 * f", "x1 + offset(o)", "x1 + (1 | g)",
    "x1 + (x2 | g) + (1 | h)", "s(x1) + x2", "1"
  )
  resps <- paste0(".y_draw.", 1:4)
  for (rhs in rhs_list) {
    fml <- as.formula(paste0("cbind(", paste(resps, collapse = ", "), ") ~ ",
                             rhs))
    fmls_expected <- lapply(resps, function(r) update(fml, paste0(r, " ~ .")))
    expect_identical(validate_response_formula(fml), fmls_expected,
                     info = rhs)
  }
})
