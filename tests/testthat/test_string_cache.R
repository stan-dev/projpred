context("cache of results that depend only on strings")

test_that("cached_by_strings() returns the cached value and empties the cache", {
  skip_if(getRversion() < "4.2.0", "the cache needs R >= 4.2.0")
  n_calls <- 0
  compute <- function() {
    n_calls <<- n_calls + 1
    NULL
  }
  key <- c("a:b", "test_string_cache", as.character(runif(1)))
  expect_null(cached_by_strings("test", key, compute))
  expect_null(cached_by_strings("test", key, compute))
  # A `NULL` value is cached, too:
  expect_identical(n_calls, 1)
  # Another type or key gives another entry:
  expect_null(cached_by_strings("test2", key, compute))
  expect_identical(n_calls, 2)
  # `NA` and "NA" are different keys:
  expect_identical(cached_by_strings("test", NA_character_, function() 1), 1)
  expect_identical(cached_by_strings("test", "NA", function() 2), 2)
  # A key longer than the limit of R variable names (10000 bytes):
  long_key <- strrep("x1 + ", 3000)
  expect_identical(cached_by_strings("test", long_key, function() 3), 3)
  expect_identical(cached_by_strings("test", long_key, function() 4), 3)
  # With `max_n` entries, the cache is emptied before the new entry:
  n_now <- utils::numhash(string_cache$tab)
  expect_identical(cached_by_strings("test", "x", function() 5, max_n = n_now),
                   5)
  expect_identical(utils::numhash(string_cache$tab), 1L)
})

test_that("the cached functions give the results of the uncached ones", {
  if (!is.null(string_cache$tab)) {
    utils::clrhash(string_cache$tab)
  }
  lbls <- c("x2:x1", "x1", " x3 : x2 ")
  ref_lbls <- vapply(strsplit(lbls, ":", fixed = TRUE), function(parts) {
    paste(sort(trimws(parts)), collapse = ":")
  }, character(1))
  for (k in 1:2) {
    expect_identical(norm_term_labels(lbls), ref_lbls)
    expect_identical(norm_term_labels(setNames(lbls, c("a", "b", "c"))),
                     setNames(ref_lbls, c("a", "b", "c")))
  }
  expect_identical(norm_term_labels(NA_character_), "")
  expect_identical(norm_term_labels("NA"), "NA")
  fmls <- list(y ~ x1 + x2:x1 + s(x3) + (1 | g), y ~ x1 + offset(o),
               ~ x1 + x2 - 1)
  for (fml in fmls) {
    for (k in 1:2) {
      expect_identical(extract_terms_response(fml),
                       extract_terms_response_uncached(fml))
    }
  }
  # The environment of a formula does not matter:
  fml_env <- local(y ~ x1 + x2)
  expect_identical(extract_terms_response(fml_env),
                   extract_terms_response_uncached(y ~ x1 + x2))
  terms <- c("x1", "x2", "x3", "x1:x2", "x1 + x2")
  info <- search_term_info(terms)
  expect_identical(info, search_term_info_uncached(terms))
  for (size in 1:4) {
    for (chosen in list(character(0), "x1", c("x1", "x2"))) {
      for (uti in c(FALSE, TRUE)) {
        for (k in 1:2) {
          expect_identical(
            select_possible_terms_size(chosen, terms, size,
                                       use_term_info = uti),
            select_possible_terms_size_uncached(chosen, terms, size,
                                                use_term_info = uti)
          )
        }
      }
    }
  }
})

test_that("the cached functions store their results under their keys", {
  skip_if(getRversion() < "4.2.0", "the cache needs R >= 4.2.0")
  cached_by_strings("test", "init", function() NULL)
  tab <- string_cache$tab
  utils::clrhash(tab)
  lbls <- c("b:a", "c")
  fml <- y ~ x1 + x2
  terms <- c("x1", "x2")
  calls <- list(
    list(fun = function() norm_term_labels(lbls),
         key = list("norm_term_labels", lbls)),
    list(fun = function() extract_terms_response(fml),
         key = list("extract_terms_response",
                    c(class(fml), deparse(fml, width.cutoff = 500L)))),
    list(fun = function() search_term_info(terms),
         key = list("search_term_info", terms)),
    list(fun = function() select_possible_terms_size(character(0), terms, 2),
         key = list("select_possible_terms_size",
                    list(2, FALSE, character(0), terms))),
    list(fun = function() select_possible_terms_size(character(0), terms, 2,
                                                     use_term_info = TRUE),
         key = list("select_possible_terms_size",
                    list(2, TRUE, character(0), terms)))
  )
  for (k in seq_along(calls)) {
    value <- calls[[k]]$fun()
    hit <- utils::gethash(tab, calls[[k]]$key)
    expect_false(is.null(hit), info = as.character(k))
    expect_identical(hit[[1]], value, info = as.character(k))
    # The second call adds no entry:
    n_entries <- utils::numhash(tab)
    expect_identical(calls[[k]]$fun(), value, info = as.character(k))
    expect_identical(utils::numhash(tab), n_entries, info = as.character(k))
  }
})

test_that("without the cache (R < 4.2.0), the values are computed each time", {
  tab_old <- string_cache$tab
  available_old <- string_cache$available
  on.exit({
    string_cache$tab <- tab_old
    string_cache$available <- available_old
  }, add = TRUE)
  string_cache$tab <- NULL
  string_cache$available <- FALSE
  n_calls <- 0
  compute <- function() {
    n_calls <<- n_calls + 1
    "value"
  }
  for (k in 1:2) {
    expect_identical(cached_by_strings("test", "key", compute), "value")
  }
  expect_identical(n_calls, 2)
  expect_null(string_cache$tab)
  expect_identical(norm_term_labels("x2:x1"), "x1:x2")
  expect_identical(extract_terms_response(y ~ x1 + x2),
                   extract_terms_response_uncached(y ~ x1 + x2))
})
