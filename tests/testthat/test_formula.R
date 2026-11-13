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
  flat <- flatten_formula(formula)
  # don't check 'flat' directly as sorting of terms is OS specific
  terms <- attr(terms(flat), "term.labels")
  expect_length(terms, 4)
  expect_true(all(terms %in% c("x", "z", "x + z + x:z | g", "x:z")))
})

test_that("split_interaction_term() works for 2-way interaction terms", {
  ia_split <- split_interaction_term("X1:X2")
  expect_identical(ia_split, "X1 + X2 + X1:X2")
  ia_noop <- split_interaction_term("X1:X2", add_lower_terms = FALSE)
  expect_identical(ia_noop, "X1:X2")
})

test_that("split_interaction_term() works for 3-way interaction terms", {
  ia_split <- split_interaction_term("X1:X2:X3")
  expect_identical(ia_split, "X1 + X2 + X3 + X1:X2 + X1:X3 + X2:X3 + X1:X2:X3")
  ia_noop <- split_interaction_term("X1:X2:X3", add_lower_terms = FALSE)
  expect_identical(ia_noop, "X1:X2:X3")
})

test_that("split_interaction_term() works for 4-way interaction terms", {
  ia_split <- split_interaction_term("X1:X2:X3:X4")
  expect_identical(ia_split, paste0(
    "X1 + X2 + X3 + X4 + ",
    "X1:X2 + X1:X3 + X1:X4 + X2:X3 + X2:X4 + X3:X4 + ",
    "X1:X2:X3 + X1:X2:X4 + X1:X3:X4 + X2:X3:X4 + ",
    "X1:X2:X3:X4"
  ))
  ia_noop <- split_interaction_term("X1:X2:X3:X4", add_lower_terms = FALSE)
  expect_identical(ia_noop, "X1:X2:X3:X4")
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

test_that("split_plus() splits term strings only at the `+` between terms", {
  expect_identical(split_plus("x1"), "x1")
  expect_identical(split_plus(c("x1 + x2", "x3")), c("x1", "x2", "x3"))
  expect_identical(
    split_plus(c("log(p + 1)", "log(q + 1) + d")),
    c("log(p + 1)", "log(q + 1)", "d")
  )
  expect_identical(
    split_plus("I(a+b) + x1:x2 + (1 | g) + s(x, k = 3) + offset(o)"),
    c("I(a + b)", "x1:x2", "(1 | g)", "s(x, k = 3)", "offset(o)")
  )
  expect_identical(split_plus("1 + x"), c("1", "x"))
  expect_identical(split_plus("poly(x3,2)"), "poly(x3, 2)")
  expect_identical(split_plus("`a b` + x2 + I(`a b`^2)"),
                   c("`a b`", "x2", "I(`a b`^2)"))
  expect_identical(split_plus(character(0)), character(0))
  expect_identical(split_plus("a + (b + c)"), c("a", "(b + c)"))
  expect_identical(split_plus("+x + y"), c("+x", "y"))
  expect_identical(split_plus("a + NULL"), c("a", "NULL"))
  # The order of the parts with `+` calls as right operands:
  expect_identical(split_plus("a + `+`(b, `+`(c, d)) + e"),
                   c("a", "b", "c", "d", "e"))
  # A long chain of `+` (a left-deep call tree of 4999 levels):
  long_terms <- paste0("x", seq_len(5000))
  expect_identical(split_plus(paste(long_terms, collapse = " + ")), long_terms)
})

test_that(paste(
  "collapse_ranked_predictors() gives the term of each model-matrix column,",
  "also if a column name is the start of another column name"
), {
  withr::local_seed(7714)
  n <- 60
  dat <- data.frame(
    x1 = rnorm(n), x3 = rnorm(n),
    f = factor(sample(c("c", "a", "ab"), n, replace = TRUE),
               levels = c("c", "a", "ab")),
    g = factor(sample(c("u", "v", "w"), n, replace = TRUE))
  )
  dat$y <- rnorm(n)
  # The column "fa" is the start of "fab", and "poly(x3, 10)1" is the start of
  # "poly(x3, 10)10". The last four formulas have no such column names:
  fmls <- list(y ~ x1 + f, y ~ x1 * f, y ~ g * f, y ~ f - 1,
               y ~ x1 + poly(x3, 10), y ~ poly(x3, 10) * x1,
               y ~ x1 + g + poly(x3, 2), y ~ poly(x3, 2) * g,
               y ~ poly(x3, 2):g, y ~ x1:g)
  for (fml in fmls) {
    info <- paste(deparse(fml), collapse = "")
    fr <- model.frame(fml, dat)
    cols <- setdiff(colnames(model.matrix(fml, data = fr)), "(Intercept)")
    labs <- attr(terms(fml), "term.labels")
    # The columns of a term are next to each other, so the terms come in the
    # order of the columns:
    expect_identical(collapse_ranked_predictors(cols, fml, fr), labs,
                     info = info)
    expect_identical(collapse_ranked_predictors(rev(cols), fml, fr), rev(labs),
                     info = info)
  }
  # An element that is not a column stays, and each term comes once:
  fml <- y ~ x1 + f
  fr <- model.frame(fml, dat)
  expect_identical(
    collapse_ranked_predictors(c("fab", NA, "x1", "fa"), fml, fr),
    c("f", NA, "x1")
  )
})

test_that(paste(
  "collapse_ranked_predictors() gives the term for a term with several",
  "columns and without contrasts"
), {
  skip_if_not_installed("splines")
  withr::local_seed(7714)
  dat <- data.frame(x1 = rnorm(60), x3 = rnorm(60), y = rnorm(60))
  fml <- y ~ x1 + splines::ns(x3, df = 3)
  fr <- model.frame(fml, dat)
  cols <- setdiff(colnames(model.matrix(fml, data = fr)), "(Intercept)")
  expect_identical(collapse_ranked_predictors(rev(cols), fml, fr),
                   c("splines::ns(x3, df = 3)", "x1"))
})

test_that(paste(
  "select_possible_terms_size() returns a long search term without line",
  "breaks"
), {
  # as.character() of a formula breaks a right-hand side of this length (1289
  # characters) into lines:
  long_term <- paste0("v", 1:200, collapse = " + ")
  allterms <- c("z", long_term)
  for (use_term_info in c(FALSE, TRUE)) {
    expect_identical(
      select_possible_terms_size("z", allterms, size = 202,
                                 use_term_info = use_term_info),
      long_term,
      info = as.character(use_term_info)
    )
  }
  # A candidate that overlaps with the chosen terms:
  expect_identical(
    select_possible_terms_size("v1", allterms, size = 201),
    paste0("v", 2:200, collapse = " + ")
  )
  expect_identical(to_character_rhs(as.formula(paste(". ~", long_term))),
                   long_term)
  # One term above 500 characters: split_plus() and the terms() label give the
  # same text.
  long_i <- paste0("I(", long_term, ")")
  expect_identical(split_plus(paste("z +", long_i)), c("z", long_i))
  lbls <- attr(terms(as.formula(paste("~ z +", long_i))), "term.labels")
  expect_identical(norm_term_labels(lbls), c("z", long_i))
})

test_that(paste(
  "smooth_covered_terms() finds the linear terms next to a smooth term of the",
  "same variable, and flatten_formula() drops them with `duplicates = FALSE`"
), {
  expect_identical(smooth_covered_terms(y ~ x + s(x) + z), "x")
  expect_identical(smooth_covered_terms(y ~ x + z + t2(x, z)), c("x", "z"))
  expect_identical(smooth_covered_terms(y ~ x + s(log(z), x)), "x")
  expect_identical(smooth_covered_terms(y ~ x + s(z)), character(0))
  expect_identical(smooth_covered_terms(y ~ x + z), character(0))
  expect_identical(
    labels(terms(flatten_formula(y ~ x + s(x) + z + (1 | g),
                                 duplicates = FALSE))),
    c("z", "s(x)", "1 | g")
  )
})

test_that(paste(
  "parse_additive_terms() fails for a named argument of s() or t2() at any",
  "position"
), {
  for (st in c("s(x, by = g)", "s(by = g, x)", "s(k = 5, x)", "s(x = x1)",
               "t2(x, z, k = 4)", "t2(k = 4, x, z)", "log(s(x, k = 5))")) {
    expect_error(parse_additive_terms(c("z", st)),
                 "arguments other than predictors are not allowed",
                 info = st)
  }
  expect_identical(parse_additive_terms(c("z", "s(x)", "t2(x, z)")),
                   c("s(x)", "t2(x, z)"))
  expect_identical(parse_additive_terms(c("s(log(z), x)", "s(X[, 1])")),
                   c("s(log(z), x)", "s(X[, 1])"))
  expect_error(extract_terms_response(y ~ x + g + s(by = g, x)),
               "arguments other than predictors are not allowed")
})

test_that(paste(
  "formula.gamm4() for an rstanarm::stan_gamm4() fit does not add the",
  "variables of the smooth terms as linear terms"
), {
  # The elements of an rstanarm::stan_gamm4() fit that formula.gamm4() uses:
  fit <- structure(list(formula = y ~ x + s(z),
                        glmod = list(formula = y ~ x + s(z) + (1 + x | g))),
                   class = "gamm4")
  expect_setequal(labels(terms(formula(fit))), c("x", "s(z)", "x | g"))
  # Without group-level terms, the formula does not change:
  fit$glmod <- NULL
  expect_identical(formula(fit), y ~ x + s(z))
})
