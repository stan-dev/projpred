# A function to remove the "cbind" part from the response in a formula (in fact,
# only relevant for formulas from "stanreg" fits with a binomial family
# with > 1 trials, but doesn't harm to apply it on other formulas as well):
rm_cbind <- function(formul) {
  formul_chr <- as.character(formul)
  stopifnot(length(formul_chr) == 3)
  formul_chr[2] <- sub("^cbind\\(", "", formul_chr[2])
  formul_chr[2] <- sub(",.*\\)$", "", formul_chr[2])
  formul <- update(formul, paste(formul_chr[c(2, 1, 3)], collapse = " "))
  return(formul)
}

# A function to remove additional response information from a formula (in fact,
# only relevant for formulas from "brmsfit"s, but doesn't harm to apply it on
# other formulas as well):
rm_addresp <- function(formul) {
  formul_chr <- as.character(formul)
  stopifnot(length(formul_chr) == 3)
  formul_chr[2] <- sub("[[:blank:]]*\\|[[:blank:]]*weights\\(wobs_col\\)$",
                       "", formul_chr[2])
  formul <- update(formul, paste(formul_chr[c(2, 1, 3)], collapse = " "))
  return(formul)
}
