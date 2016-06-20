.onAttach <- function(...) {
  ver <- utils::packageVersion("glmproj")
  packageStartupMessage("This is glmproj version ", ver)
}
