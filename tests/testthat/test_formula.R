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

test_that(paste(
  "select_possible_terms_size() adds a smooth term only after the linear",
  "terms of its variables that are search terms"
), {
  allterms <- c("1", "s(x1)", "x2", "x1")
  for (use_ti in c(FALSE, TRUE)) {
    expect_identical(
      select_possible_terms_size("1", allterms, size = 2,
                                 use_term_info = use_ti),
      c("x2", "x1"),
      info = paste("use_term_info:", use_ti)
    )
  }
  expect_identical(select_possible_terms_size(c("1", "x2"), allterms, size = 3),
                   "x1")
  expect_identical(select_possible_terms_size(c("1", "x1"), allterms, size = 3),
                   c("s(x1)", "x2"))
  # The submodel with all search terms has the same size on all search paths:
  expect_identical(nterms_search_all(allterms), 3L)
  dat <- data.frame(y = 1:4, x1 = 1:4, x2 = 1:4, g = gl(2, 2))
  fml <- y ~ x2 + s(x1) + (1 | g)
  expect_identical(
    nterms_search_all(union("1", split_formula(fml, data = dat))), 4L
  )
  # Without the linear term in the search terms:
  allterms <- c("1", "s(x1)", "x2")
  expect_identical(select_possible_terms_size("1", allterms, size = 2),
                   c("s(x1)", "x2"))
  expect_identical(nterms_search_all(allterms), 2L)
  # A smooth term of two variables:
  allterms <- c("1", "t2(x1, x3)", "x1", "x3")
  expect_identical(select_possible_terms_size(c("1", "x1"), allterms, size = 3),
                   "x3")
  expect_identical(
    select_possible_terms_size(c("1", "x1", "x3"), allterms, size = 4),
    "t2(x1, x3)"
  )
  allterms <- c("1", "t2(x1, x3)", "x1")
  expect_identical(select_possible_terms_size(c("1", "x1"), allterms, size = 3),
                   "t2(x1, x3)")
  # A search term with the linear and the smooth term:
  expect_identical(
    select_possible_terms_size("1", c("1", "x1 + s(x1)", "x2"), size = 3),
    "x1 + s(x1)"
  )
  # Linear terms in search terms with a smooth term do not count (otherwise,
  # these two search terms would wait for each other):
  expect_identical(
    select_possible_terms_size("1", c("1", "x2 + s(x1)", "x1 + s(x2)"),
                               size = 3),
    c("x2 + s(x1)", "x1 + s(x2)")
  )
  # Non-syntactic names:
  allterms <- c("1", "s(`my x`)", "`my x`")
  expect_identical(select_possible_terms_size("1", allterms, size = 2),
                   "`my x`")
  expect_identical(
    select_possible_terms_size(c("1", "`my x`"), allterms, size = 3),
    "s(`my x`)"
  )
})

test_that(paste(
  "the warning about the projection onto the full multilevel model uses the",
  "submodel with all search terms"
), {
  withr::local_options(projpred.warn_instable_projections = TRUE)
  dat <- data.frame(y = 1:12, x1 = sin(1:12), x2 = cos(1:12), g = gl(3, 4))
  refmod_fake <- function(fml) {
    list(formula = fml, fetch_data = function(...) dat,
         family = list(family = "gaussian", for_latent = FALSE,
                       for_augdat = FALSE))
  }
  parse_args <- function(fml, search_terms) {
    parse_args_varsel(refmod_fake(fml), method = "forward", refit_prj = TRUE,
                      nterms_max = NULL, nclusters = 1,
                      search_terms = search_terms, verbose = 0)
  }
  warn_full <- "the projection onto the full model can be instable"
  # varsel() and cv_varsel():
  expect_warning(args <- parse_args(y ~ x2 + s(x1) + (1 | g), NULL),
                 warn_full)
  expect_identical(args$nterms_all, 4L)
  expect_equal(args$nterms_max, 4)
  expect_warning(args <- parse_args(y ~ x2 + s(x1) + (1 | g),
                                    c("s(x1)", "x2", "(1 | g)")),
                 warn_full)
  expect_identical(args$nterms_all, 3L)
  # Search terms without group-level terms:
  expect_no_warning(args <- parse_args(y ~ x1 + x2 + (1 | g), c("x1", "x2")))
  expect_identical(args$nterms_all, 2L)
  # project(): a submodel with a smooth term and without the linear term of its
  # variable is the full model:
  refmod <- refmod_fake(y ~ s(x1) + (1 | g))
  expect_true(is_full_projection(c("s(x1)", "(1 | g)"), refmod, refmod))
  expect_true(is_full_projection(c("x1", "s(x1)", "(1 | g)"), refmod, refmod))
  expect_false(is_full_projection(c("x1", "(1 | g)"), refmod, refmod))
  vs <- structure(list(args_search = list(search_terms = NULL)),
                  class = "vsel")
  expect_true(is_full_projection(c("x1", "s(x1)", "(1 | g)"), vs, refmod))
  vs$args_search$search_terms <- c("x1", "(1 | g)")
  expect_true(is_full_projection(c("x1", "(1 | g)"), vs, refmod))
  # Group-level slopes:
  refmod <- refmod_fake(y ~ x1 + (x1 | g))
  expect_true(is_full_projection(c("x1", "(x1 | g)"), refmod, refmod))
  expect_false(is_full_projection(c("x1", "(1 | g)"), refmod, refmod))
  vs$args_search$search_terms <- NULL
  expect_true(is_full_projection(c("(1 | g)", "x1 + (x1 | g)"), vs, refmod))
  expect_false(is_full_projection(c("(1 | g)", "x1 + (1 | g)"), vs, refmod))
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
