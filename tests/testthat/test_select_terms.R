context("select_possible_terms_size() with precomputed term information")

test_that(paste(
  "select_possible_terms_size() gives the same candidates with and without",
  "`term_info`"
), {
  term_sets <- list(
    union("1", split_formula(
      y ~ x1 + x2 + x3 + f + x1:x2 + poly(x4, 2) + I(x5^2) + x3:x6,
      data = NULL
    )),
    c("1", "x1", "x2", "x1:x2", "x3", "f", "poly(x4,2)"),
    c("1", paste0("x", 1:15))
  )
  for (st in term_sets) {
    info <- search_term_info(st)
    expect_false(is.null(info))
    for (rep_i in 1:10) {
      set.seed(rep_i)
      chosen <- character()
      for (size in 1:8) {
        cands <- select_possible_terms_size(chosen, st, size)
        cands_info <- select_possible_terms_size(chosen, st, size,
                                                 term_info = info)
        expect_identical(cands_info, cands)
        if (is.null(cands)) next
        chosen <- c(chosen, cands[sample.int(length(cands), 1)])
      }
    }
  }
})

test_that("search_term_info() returns `NULL` for smooth or group-level terms", {
  expect_null(search_term_info(c("1", "x1", "s(x2)")))
  expect_null(search_term_info(c("1", "x1", "(1 | g)")))
})
