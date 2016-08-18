project.varsel <- function(fit, obj, size, ...) {
  params <- extract_params(fit)
  ns_total <- ncol(params$b)
  args <- list(...)
  args$ns <- min(ifelse(is.null(args$ns), 400, args$ns), ns_total)
  family_kl <- kl_helpers(family(fit))
  v_inds <- obj$submodel$chosen[1:size]

  d <- list(x = params$x[,varinds], w = params$w, offset = params$offset)

  # Sample indices to be used with forward selection and final projection
  s_ind <- round(seq(1, ns_total, length.out  = args$ns))
  b0 <- matrix(rowMeans(params$b[v_inds, , drop = F]), ncol = 1)
  p <- list(mu = family_kl$linkinv(d$x%*%params$b[v_inds, s_ind]),
            dis = params$dis[s_ind])

  if(family_kl$family == 'gaussian') {
    q <- list(b = solve(crossprod(d$x), crossprod(d$x, p$mu)))
    q$dis <- sqrt(p$dis^2 + colMeans((p$mu - d$x%*%q$b)^2))
  } else {
    res <- sapply(1:args$ns, function(s_ind, p, d, b0, family_kl) {
      NR(list(mu = p$mu[, s_ind, drop = F], dis = p$dis[s_ind]), d, b0, family_kl)
      }, p, d, b0, family_kl)
    q <- list(b = do.call(cbind, res['b',]))
    if('dis' %in% rownames(res)) q$dis <- unlist(res['dis',])
  }
  q$chosen <- v_inds
  q$family <- family(fit)
  structure(q, 'glmproj')
}
