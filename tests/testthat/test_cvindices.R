context('cv-indices')

# tests for cvfolds and cvind

test_that('k is valid', {
  cvfuns <- c(cvfolds, cvind)
  for (cvfun in cvfuns) {
    expect_error(cvfun(n = 10, k = 1000),
                 'cannot exceed n')
    expect_error(cvfun(n = 10, k = 1),
                 'must be at least 2')
    expect_error(cvfun(n = 10, k = c(4, 9)),
                 'a single integer value')
    expect_error(cvfun(n = 10, k = 'a'),
                 'a single integer value')
  }
})
