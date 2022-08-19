context("formula")

test_that(paste(
  "check that we recover the correct terms for a simple linear",
  "model without interactions or group terms"
), {
  formula <- y ~ x + z
  tt <- extract_terms_response(formula)
  expect_length(tt$individual_terms, 2)
  expect_length(tt$interaction_terms, 0)
  expect_length(tt$group_terms, 0)
  expect_length(tt$response, 1)

  formula <- y ~ 1
  tt <- extract_terms_response(formula)
  expect_length(tt$individual_terms, 0)
  expect_length(tt$interaction_terms, 0)
  expect_length(tt$group_terms, 0)
  expect_length(tt$response, 1)
})

test_that(paste(
  "check that we recover the correct terms for a simple ",
  "multilevel model with interactions"
), {
  formula <- y ~ x + z + x:z + (x:z | g)
  tt <- extract_terms_response(formula)
  expect_length(tt$individual_terms, 2)
  expect_length(tt$interaction_terms, 1)
  expect_length(tt$group_terms, 1)
  expect_length(tt$response, 1)
})

test_that(paste(
  "check that we recover the correct terms for a simple ",
  "additive model without interactions"
), {
  formula <- y ~ x + s(z)
  tt <- extract_terms_response(formula)
  expect_length(tt$individual_terms, 1)
  expect_length(tt$interaction_terms, 0)
  expect_length(tt$additive_terms, 1)
  expect_length(tt$group_terms, 0)
  expect_length(tt$response, 1)
})

test_that(paste(
  "check that we recover the correct terms for a simple ",
  "additive model with multidimensional interactions"
), {
  formula <- y ~ t2(x, z)
  tt <- extract_terms_response(formula)
  expect_length(tt$individual_terms, 0)
  expect_length(tt$interaction_terms, 0)
  expect_length(tt$additive_terms, 1)
  expect_length(tt$group_terms, 0)
  expect_length(tt$response, 1)
})

test_that(paste(
  "check that we return a list of length one containing the same formula for a",
  "single response"
), {
  formula <- y ~ x + z
  valrespformul <- validate_response_formula(formula)
  expect_type(valrespformul, "list")
  expect_length(valrespformul, 1)
  expect_equal(formula, valrespformul[[1]])
})

test_that("check that we return a list of formulas for multiple responses", {
  formula <- cbind(y.1, y.2) ~ x + z
  valrespformul <- validate_response_formula(formula)
  expect_type(valrespformul, "list")
  expect_length(valrespformul, 2)
  expect_equal(y.1 ~ x + z, valrespformul[[1]])
  expect_equal(y.2 ~ x + z, valrespformul[[2]])
})

test_that("check that we properly flatten a formula with duplicated terms", {
  formula <- (y ~ x + z + x:z + (1 | g) + (x | g) + (z | g) + (x + z | g) +
                (x + z + x:z | g))
  flat <- projpred:::flatten_formula(formula)
  # don't check 'flat' directly as sorting of terms is OS specific
  terms <- attr(terms(flat), "term.labels")
  expect_length(terms, 4)
  expect_true(all(terms %in% c("x", "z", "x + z + x:z | g", "x:z")))
})

test_that("check that we properly split a formula", {
  formula <- y ~ x + z + x:z + (x + z + x:z | g)
  sp <- split_formula(formula)
  expect_length(sp, 11)
  expect_length(setdiff(
    c(
      "1", "x", "z", "x + z + x:z", "(1 | g)", "x + (x | g)", "z + (z | g)",
      "x + (1 | g)", "z + (1 | g)", "x + z + x:z + (1 | g)",
      "x + z + x:z + (x + z + x:z | g)"
    ),
    sp
  ), 0)

  formula <- y ~ 0 + (x + z + x:z | g)
  sp <- split_formula(formula)
  expect_length(sp, 7)
  expect_length(setdiff(
    c(
      "(1 | g) + 0", "x + (x | g) + 0", "z + (z | g) + 0",
      "x + (1 | g) + 0", "z + (1 | g) + 0", "x + z + x:z + (1 | g) + 0",
      "x + z + x:z + (x + z + x:z | g) + 0"
    ),
    sp
  ), 0)

  formula <- y ~ (x + z + x:z | g)
  sp <- split_formula(formula)
  expect_length(sp, 8)
  expect_length(setdiff(
    c(
      "(1 | g)", "x + (x | g)", "z + (z | g)",
      "x + (1 | g)", "z + (1 | g)", "x + z + x:z + (1 | g)",
      "x + z + x:z + (x + z + x:z | g)"
    ),
    sp
  ), 0)

  formula <- y ~ (0 + x + z + x:z | g)
  sp <- split_formula(formula)
  expect_length(sp, 4)
  expect_length(setdiff(
    c(
      "1", "x + (0 + x | g)", "z + (0 + z | g)",
      "x + z + x:z + (0 + x + z + x:z | g)"
    ),
    sp
  ), 0)

  formula <- y ~ x + z + x:z
  sp <- split_formula(formula)
  expect_length(sp, 4)
  expect_length(setdiff(
    c(
      "1", "x", "z", "x + z + x:z"
    ),
    sp
  ), 0)

  formula <- y ~ 0 + x + z + x:z
  sp <- split_formula(formula)
  expect_length(sp, 3)
  expect_length(setdiff(
    c(
      "x + 0", "z + 0", "x + z + x:z + 0"
    ),
    sp
  ), 0)

  formula <- y ~ s(x) + s(z)
  sp <- split_formula(formula)
  expect_length(sp, 5)
  expect_length(setdiff(
    c(
      "1", "s(x)", "s(z)", "x", "z"
    ),
    sp
  ), 0)

  formula <- y ~ t2(x, z)
  sp <- split_formula(formula)
  expect_length(sp, 4)
  expect_length(setdiff(
    c(
      "1", "t2(x, z)", "x", "z"
    ),
    sp
  ), 0)
})

test_that("check that we can identify formulas with group terms", {
  formula <- y ~ x + z + x:z + (x + z + x:z | g)
  expect_true(formula_contains_group_terms(formula))

  formula <- y ~ x + z
  expect_false(formula_contains_group_terms(formula))
})

test_that(paste(
  "check that we can subset a formula and update the data",
  "columns properly"
), {
  data <- data.frame(y = rnorm(20), x = matrix(rnorm(40), 20, 4))
  fake_y <- matrix(rnorm(20), 20, 1)
  formula <- y ~ x.1 + x.2 + x.3 + x.4
  s <- subset_formula_and_data(formula, c("x.1", "x.3"), data, y = fake_y)

  cols <- colnames(s$data)
  expect_equal(cols[1], ".y")
  expect_equal(s$data[, 1], as.vector(fake_y))
  expect_equal(s$formula, .y ~ x.1 + x.3)

  formula <- y ~ x.1 + x.2 + x.3 + x.4
  fake_y <- matrix(rnorm(40), 20, 2)
  s <- subset_formula_and_data(formula, c("x.1", "x.3"), data, y = fake_y)
  cols <- colnames(s$data)
  expect_equal(cols[1:2], c(".y.1", ".y.2"))
  expect_true(all(s$data[, 1:2] == fake_y))
  expect_equal(s$formula, cbind(.y.1, .y.2) ~ x.1 + x.3)
})

test_that("check that we count terms correctly", {
  formula <- y ~ x + z
  expect_equal(count_terms_in_formula(formula), 3)

  formula <- y ~ x + z + x:z
  expect_equal(count_terms_in_formula(formula), 4)

  formula <- y ~ x + z + x:z + (1 | g)
  expect_equal(count_terms_in_formula(formula), 5)

  formula <- y ~ x + z + x:z + (x | g)
  expect_equal(count_terms_in_formula(formula), 6)

  expect_equal(count_terms_chosen(c("x", "z")), 3)
  expect_equal(count_terms_chosen(c("x", "z", "x:z")), 4)
  expect_equal(count_terms_chosen(c("x", "z", "x:z", "x + (x | g)")), 6)
})

test_that("check that we correctly sort models by size", {
  submodels <- c(
    y ~ x + z,
    y ~ x + z + x:z,
    y ~ x + z + x:z + (1 | g),
    y ~ x + z + x:z + (x | g)
  )
  s <- sort_submodels_by_size(submodels)
  expect_null(s[[1]])
  expect_equal(unlist(s), as.character(submodels))
})

test_that("select_possible_terms_size() avoids redundant models", {
  chosen <- "x1 * x2"
  size_chosen <- count_terms_chosen(chosen)
  expect_equal(size_chosen, 4)
  allterms <- c("x1 * x2", "x1", "x2")
  expect_null(
    select_possible_terms_size(chosen, allterms, size = size_chosen + 1)
  )

  chosen <- "x + (x | g)"
  size_chosen <- count_terms_chosen(chosen)
  expect_equal(size_chosen, 4)
  allterms <- c("x + (x | g)", "x", "(1 | g)", "(x | g)")
  expect_null(
    select_possible_terms_size(chosen, allterms, size = size_chosen + 1)
  )

  chosen <- "s(x1)"
  size_chosen <- count_terms_chosen(chosen)
  expect_equal(size_chosen, 2)
  allterms <- c("s(x1)", "x1")
  expect_null(
    select_possible_terms_size(chosen, allterms, size = size_chosen + 1)
  )
})

test_that("select_possible_terms_size() works for non-redundant models", {
  chosen <- "x1"
  size_chosen <- count_terms_chosen(chosen)
  expect_equal(size_chosen, 2)
  allterms <- c("x1 * x2", "x1", "x2")
  expect_identical(
    select_possible_terms_size(chosen, allterms, size = size_chosen + 1),
    "x2"
  )
  chosen <- c(chosen, "x2")
  size_chosen <- count_terms_chosen(chosen)
  expect_equal(size_chosen, 3)
  expect_identical(
    select_possible_terms_size(chosen, allterms, size = size_chosen + 1),
    "x1:x2"
  )

  chosen <- "x"
  size_chosen <- count_terms_chosen(chosen)
  expect_equal(size_chosen, 2)
  allterms <- c("x", "(1 | g)", "(x | g)")
  expect_identical(
    select_possible_terms_size(chosen, allterms, size = size_chosen + 1),
    "(1 | g)"
  )

  chosen <- "x1"
  size_chosen <- count_terms_chosen(chosen)
  expect_equal(size_chosen, 2)
  allterms <- c("s(x1)", "x1")
  expect_identical(
    select_possible_terms_size(chosen, allterms, size = size_chosen + 1),
    "s(x1)"
  )
})
