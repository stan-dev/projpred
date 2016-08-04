#' @export cv_varsel 
cv_varsel <- function(fit, fits = NULL, ...) {
  UseMethod('cv_varsel')
}

cv_varsel.stanreg <- function(fit, fits = NULL, ...) {

  if(is.null(fits)) fits <- cv_fit(fit)
  k <- nrow(fits)
  verbose <- ifelse(is.null(list(...)$verbose), F, list(...)$verbose)

  params <- extract_params(fit)
  if(ncol(params$x) < 2)
    stop('Data must have at least 2 features.')
  if(!(family(fit)$family %in% c('gaussian','binomial','poisson')))
    stop(paste0(family(fit)$family, 'family not yet supported.'))


  msgs <- paste('Forward selection for the', c('full model.', paste0('fold number ', 1:k,'/',k,'.')))
  allfits <- rbind(list(fit=fit,d_test=NA),fits)

  sel <- mapply(function(fit, d_test, msg, verbose) {
    if(verbose) print(msg)
    varsel(fit, d_test, ...)
  }, allfits[,'fit'], allfits[,'d_test'], msgs, MoreArgs = list(verbose = verbose))

  sel_sub <- simplify2array(sel['submodel',-1])
  varnames <- setdiff(names(sel[['submodel', 2]]), c('chosen','kl'))

  sub_avg <- as.data.frame(sapply(varnames,
    function(varname, sel_sub) rowMeans(simplify2array(sel_sub[varname,])), sel_sub))

  pctch <- sapply(1:nrow(sub_avg), function(ind, sel_sub, chosen_full, k) {
    char <- as.character(chosen_full[ind])
    count <- unname(table(simplify2array(sel_sub['chosen',])[1:ind,])[char])
    ifelse(is.na(count), 0, count/k)
  }, sel_sub, sel[['submodel',1]]$chosen, k)

  res <- list(submodel = cbind(as.data.frame(sel[['submodel', 1]]), sub_avg, pctchosen = pctch))

  res$full <- apply(simplify2array(sel['full',-1]), 1, function(x) mean(unlist(x)))

  res$family <- family(fit)

  structure(res, class = 'varsel')

}
