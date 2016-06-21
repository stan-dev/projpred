library(glmproj)

context("Model-specific helper functions")
test_that("gaussian family returns the correct helper functions",{
  family <- gaussian()
  funs <- kl_helpers(family)
  expect_type(funs$kl,'closure')
  expect_type(funs$dkl,'closure')
  expect_type(funs$dme,'closure')
  expect_type(funs$dis,'closure')
})

test_that("binomial family with logit link returns the correct helper functions",{
  family <- binomial(link = 'logit')
  funs <- kl_helpers(family)
  expect_type(funs$kl,'closure')
  expect_type(funs$dkl,'closure')
  expect_type(funs$dme,'closure')
  expect_type(funs$dis,'closure')
})

test_that("binomial family with probit link returns the correct helper functions",{
  family <- binomial(link = 'probit')
  funs <- kl_helpers(family)
  expect_type(funs$kl,'closure')
  expect_type(funs$dkl,'closure')
  expect_type(funs$dme,'closure')
  expect_type(funs$dis,'closure')
})

test_that("poisson family returns the correct helper functions",{
  family <- poisson()
  funs <- kl_helpers(family)
  expect_type(funs$kl,'closure')
  expect_type(funs$dkl,'closure')
  expect_type(funs$dme,'closure')
  expect_type(funs$dis,'closure')
})
