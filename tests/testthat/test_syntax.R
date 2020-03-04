context('syntax')

# test for simple fits but with varying syntax (e.g. varying formula etc.)


if (require(brms)) {


  # load the mesquite data
  data('mesquite', package = 'projpred')
  
  
  # fit the model with some transformations on the target variable and the original inputs
  SW( 
    fit <- brm(LeafWt ~ log(Diam1) + log(Diam2) + log(CanHt) + log(TotHt) + log(Dens) +
                  log(Diam1)*log(Diam2) + Group, 
                   data = mesquite, refresh=0, chain=2)
  )
  
  # selection
  SW(cvs <- cv_varsel(fit, verbose=F))
  suggested_size <- suggest_size(cvs)
  
  # project onto some model size
  proj <- project(cvs, nv = 3)
  
  
  
  test_that('varsel/cv_varsel/project return objects with correct types', {
    expect_true('vsel' %in% class(vs))
    expect_true('cvsel' %in% class(cvs))
    expect_true('projection' %in% class(proj))
  })
  
  test_that('suggested model size is ok', {
    expect_true(!is.na(suggested_size))
  })
  
}  
