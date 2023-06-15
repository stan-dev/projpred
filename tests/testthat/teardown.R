#_______________________________________________________________________________
# Teardown for the unit tests
#_______________________________________________________________________________

if (run_prll) {
  if (dopar_backend == "doParallel") {
    doParallel::stopImplicitCluster()
  } else if (dopar_backend == "doFuture") {
    future::plan(future::sequential)
    options(export_default)
    rm(export_default)
  } else {
    stop("Unrecognized `dopar_backend`.")
  }
}
