test_that("check that we recover the correct terms for a simple linear model without interactions or group terms", {
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

test_that("check that we recover the correct terms for a simple multilevel model with interactions", {
  formula <- y ~ x + z + x:z + (x:z | g)
  tt <- extract_terms_response(formula)
  expect_length(tt$individual_terms, 2)
  expect_length(tt$interaction_terms, 1)
  expect_length(tt$group_terms, 1)
  expect_length(tt$response, 1)
})

test_that("check that we return the same formula for a single response", {
  formula <- y ~ x + z
  expect_equal(formula, validate_response_formula(formula))
})

test_that("check that we return a list of formulas for multiple responses", {
  formula <- cbind(y.1, y.2) ~ x + z
  expect_true(length(validate_response_formula(formula)) == 2)
  expect_equal(y.1 ~ x + z, validate_response_formula(formula)[[1]])
})

test_that("check that we properly flatten a formula with duplicated terms", {
  formula <- (y ~ x + z  + x:z + (1 | g) + (x | g) + (z | g) + (x + z | g) +
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
      "x + z + x:z + (x:z | g)"
    ),
    sp), 0)

  formula <- y ~ 0 + (x + z + x:z | g)
  sp <- split_formula(formula)
  expect_length(sp, 7)
  expect_length(setdiff(
    c(
      "(1 | g) + 0", "x + (x | g) + 0", "z + (z | g) + 0",
      "x + (1 | g) + 0", "z + (1 | g) + 0", "x + z + x:z + (1 | g) + 0",
      "x + z + x:z + (x:z | g) + 0"
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
      "x + z + x:z + (x:z | g)"
    ),
    sp), 0)

  formula <- y ~ (0 + x + z + x:z | g)
  sp <- split_formula(formula)
  expect_length(sp, 4)
  expect_length(setdiff(
    c(
      "1", "x + (0 + x | g)", "z + (0 + z | g)", "x + z + x:z + (0 + x:z | g)"
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
    sp), 0)

  formula <- y ~ 0 + x + z + x:z
  sp <- split_formula(formula)
  expect_length(sp, 3)
  expect_length(setdiff(
    c(
      "x + 0", "z + 0", "x + z + x:z + 0"
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

test_that("check that we can subset a formula and update the data columns properly", {
  data <- data.frame(y=rnorm(20), x=matrix(rnorm(40), 20, 4))
  fake_y <- matrix(rnorm(20), 20, 1)
  formula <- y ~ x.1 + x.2 + x.3 + x.4
  s <- subset_formula_and_data(formula, c("x.1", "x.3"), data, y=fake_y)

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
  expect_equal(count_terms_in_subformula(formula), 3)

  formula <- y ~ x + z + x:z
  expect_equal(count_terms_in_subformula(formula), 4)

  formula <- y ~ x + z + x:z + (1 | g)
  expect_equal(count_terms_in_subformula(formula), 6)

  formula <- y ~ x + z + x:z + (x | g)
  expect_equal(count_terms_in_subformula(formula), 7)

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

test_that("check redundant models", {
  current <- c("x + (x | g)")
  expect_true(is_next_submodel_redundant(current, "x"))
  expect_true(is_next_submodel_redundant(current, "(1 | g)"))
  expect_true(is_next_submodel_redundant(current, "(x | g)"))
  expect_false(is_next_submodel_redundant(current, "(z | g)"))
  expect_false(is_next_submodel_redundant(current, "z"))
})

test_that("check reduce models", {
  chosen <- c("x + (x | g)", "x", "(1 | g)", "(x | g)")
  r <- reduce_models(chosen)
  expect_equal(r, chosen[1])

  chosen <- c("x", "(1 | g)", "(x | g)")
  r <- reduce_models(chosen)
  expect_equal(r, chosen)
})
