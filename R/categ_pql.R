# PQL fit of multilevel submodels for the brms::categorical() family
#
# The algorithm of mclogit::mblogit() with `method = "PQL"`, `estimator = "ML"`,
# `catCov = "free"` and the default inner optimizer nlminb() with an analytic
# gradient, but with base-R dense matrices.  Two differences from mclogit (up
# to at least version 0.9.17.1):
# * The gradient of the inner objective sums the diagonal blocks of the inverse
#   of `Z'WZ + Sigma^-1` correctly.  mclogit:::sum_blockDiag() uses a wrong
#   offset for blocks larger than 1 x 1 (more than two response categories or
#   group-level slopes), which stops nlminb() early.
# * The data have one row per observation and a response matrix of the
#   augmented-data weights (see categ_mat_resp()).  The approximate deviance
#   of the augmented data is this deviance plus the constant
#   -2 sum(W * log(W / rowSums(W))) >= 0.
# The returned object has the elements of an `mmblogit` object that projpred
# and mclogit:::predict.mmblogit() use.

# Convert augmented data (N observations in C blocks, one block per response
# category, see `?augdat-internals`) and the augmented-length weights into N
# rows with an N x C response matrix of the weights.  Returns `NULL` if `data`
# does not have this layout.
categ_mat_resp <- function(formula, data, weights) {
  y_nm <- extract_terms_response(formula)$response
  y <- data[[y_nm]]
  if (!is.factor(y) || length(weights) != nrow(data)) return(NULL)
  ncat <- nlevels(y)
  nobs <- nrow(data) %/% ncat
  if (nobs == 0 || nobs * ncat != nrow(data) ||
      !identical(as.integer(y), rep(seq_len(ncat), each = nobs))) {
    return(NULL)
  }
  idx1 <- seq_len(nobs)
  nms_pred <- setdiff(names(data), y_nm)
  for (k in seq_len(ncat)[-1]) {
    idxk <- idx1 + (k - 1L) * nobs
    for (nm in nms_pred) {
      if (!identical(data[[nm]][idx1], data[[nm]][idxk])) return(NULL)
    }
  }
  resp_nm <- "projpred_internal_resp_mat"
  if (resp_nm %in% names(data)) {
    stop("Need to write to column `", resp_nm, "` of `data`, but that column ",
         "already exists. Please rename this column in `data` and try again.")
  }
  data_mat <- data[idx1, nms_pred, drop = FALSE]
  data_mat[[resp_nm]] <- matrix(weights, nrow = nobs, ncol = ncat,
                                dimnames = list(NULL, levels(y)))
  formula[[2]] <- as.name(resp_nm)
  return(list(formula = formula, data = data_mat))
}

# Split a one-sided formula `~ <terms> | <group>` as mclogit does.
categ_pql_random <- function(formula) {
  fo <- delete.response(terms(formula))
  attributes(fo) <- NULL
  if (length(fo[[2]]) < 2 || as.character(fo[[2]][1]) != "|") {
    stop("Missing '|' operator in a group-level term.")
  }
  groups <- fo
  fo[2] <- fo[[2]][2]
  groups[2] <- groups[[2]][3]
  return(list(formula = structure(fo, class = "formula"),
              groups = all.vars(groups)))
}

# The PQL fit.  `formula` has the response matrix on the left-hand side;
# `random` is a list of formulas `~ <terms> | <group>`.  Returns `NULL` for
# cases that this function does not handle (nested group terms, a group
# variable in more than one term, a rank-deficient model matrix, more than
# `max_cells` cells in the group-level design matrix, a failure in the first
# iteration); the caller then uses mclogit::mblogit().
categ_pql_fit <- function(formula, random, data, epsilon = 1e-8, maxit = 25,
                          maxit.inner = 100, max_cells = 1e7) {
  if (!is.numeric(epsilon) || length(epsilon) != 1 || !is.finite(epsilon) ||
      epsilon <= 0) {
    stop("`epsilon` must be a single number > 0.")
  }
  if (!is.numeric(maxit) || length(maxit) != 1 || !is.finite(maxit) ||
      maxit < 1) {
    stop("`maxit` must be a single number >= 1.")
  }
  # predict.mmblogit() is registered when the namespace is loaded:
  loadNamespace("mclogit")
  random_raw <- random
  random <- lapply(random, categ_pql_random)
  if (any(lengths(lapply(random, "[[", "groups")) != 1)) return(NULL)
  grp_nms <- vapply(random, "[[", character(1), "groups")
  if (anyDuplicated(grp_nms)) return(NULL)

  # Model frame and design matrices, as in mclogit::mblogit():
  mt <- attr(model.frame(formula, data = data, drop.unused.levels = TRUE),
             "terms")
  rf <- as.formula(paste(c(". ~ .", unique(unlist(lapply(random_raw,
                                                          all.vars)))),
                         collapse = " + "))
  mf <- model.frame(update(formula, rf), data = data,
                    drop.unused.levels = TRUE)
  Y <- model.response(mf, "any")
  X <- model.matrix(mt, mf)
  if (qr(X)$rank < ncol(X)) return(NULL)
  contrasts <- attr(X, "contrasts")
  xlevels <- .getXlevels(mt, mf)
  ncat <- ncol(Y)
  nobs <- nrow(Y)
  D <- diag(ncat)[, -1, drop = FALSE]
  dimnames(D) <- list(colnames(Y), colnames(Y)[-1])
  w_obs <- rowSums(Y)
  y_rel <- Y / w_obs
  y_rel[w_obs == 0, ] <- 0
  y <- as.vector(t(y_rel))
  w <- rep(w_obs, each = ncat)
  s <- rep(seq_len(nobs), each = ncat)
  XD <- X %x% D
  colnames(XD) <- paste0(rep(colnames(D), ncol(X)), "~",
                         rep(colnames(X), each = ncol(D)))
  nlong <- nrow(XD)

  # Group-level design matrices (columns: group levels, and within a level the
  # d coefficients) with the column names of mclogit:::mkZ():
  d <- integer(length(random))
  m <- integer(length(random))
  groups <- list()
  vc_nms <- list()
  Z_list <- list()
  for (k in seq_along(random)) {
    Zk <- suppressWarnings(model.matrix(terms(random[[k]]$formula), mf,
                                        contrasts.arg = contrasts))
    ZDk <- Zk %x% D
    colnames(ZDk) <- gsub("(Intercept)", "1", fixed = TRUE, paste0(
      rep(colnames(D), ncol(Zk)), "~", rep(colnames(Zk), each = ncol(D))
    ))
    g <- rep(as.factor(mf[[grp_nms[k]]]), each = ncat)
    d[k] <- ncol(ZDk)
    m[k] <- nlevels(g)
    q_k <- sum(d[seq_len(k)] * m[seq_len(k)])
    if (as.numeric(nlong) * q_k > max_cells || as.numeric(q_k)^2 > max_cells) {
      return(NULL)
    }
    Zfull <- matrix(0, nrow = nlong, ncol = m[k] * d[k])
    jk <- rep((as.integer(g) - 1L) * d[k], d[k]) +
      rep(seq_len(d[k]), each = nlong)
    Zfull[cbind(rep(seq_len(nlong), d[k]), jk)] <- ZDk
    colnames(Zfull) <- if (d[k] > 1) {
      as.vector(outer(colnames(ZDk), levels(g), paste, sep = "|"))
    } else {
      levels(g)
    }
    Z_list[[k]] <- Zfull
    groups[[grp_nms[k]]] <- g
    vc_nms[[k]] <- colnames(ZDk)
  }
  Z <- do.call(cbind, Z_list)
  z_idx <- split(seq_len(ncol(Z)), rep(seq_along(d), d * m))

  # Start values as in mclogit:
  eta <- local({
    n_ <- w * y + 0.5
    f <- n_ / rowsum(n_, s, reorder = FALSE)[s]
    log(f) - stats::ave(log(f), s)
  })
  prb <- categ_pql_probs(eta, s)

  parms <- list(alpha = NULL, Phi = NULL)
  converged <- FALSE
  step_back <- FALSE
  fit <- NULL
  last_parms <- NULL
  last_deviance <- NULL
  prev_last_deviance <- NULL
  last_eta <- eta
  for (iter in seq_len(maxit)) {
    y_star <- eta + (y - prb) / prb
    prev_last_parms <- last_parms
    last_parms <- parms
    parms <- tryCatch(categ_pql_inner(
      parms, y_star = y_star, X = XD, Z = Z, z_idx = z_idx, d = d, m = m,
      w = w, prb = prb, s = s, maxit_inner = maxit.inner
    ), categ_pql_numeric_error = function(e) e)
    step_back <- FALSE
    if (inherits(parms, "categ_pql_numeric_error")) {
      if (is.null(fit)) return(NULL)
      msg <- conditionMessage(parms)
      if (length(prev_last_deviance) && last_deviance > prev_last_deviance &&
          length(prev_last_parms$alpha)) {
        warning("Numeric problems in inner iteration (", msg, ") and ",
                "previous step increased deviance, stepping back twice",
                call. = FALSE)
        parms <- prev_last_parms
      } else {
        warning("Numeric problems in inner iteration (", msg, "), stepping ",
                "back", call. = FALSE)
        parms <- last_parms
      }
      step_back <- TRUE
      fit <- categ_pql_eval(parms, XD, Z, z_idx, d, y, w, s)
    } else {
      last_fit <- fit
      fit <- categ_pql_eval(parms, XD, Z, z_idx, d, y, w, s)
      if (!is.finite(fit$deviance)) {
        if (is.null(last_fit)) return(NULL)
        warning("Non-finite deviance, backing up", call. = FALSE)
        step_back <- TRUE
        parms <- last_parms
        fit <- last_fit
      }
    }
    eta <- fit$eta
    prb <- fit$prb
    if (step_back) break
    if (iter > 1) last_eta <- last_fit$eta
    prev_last_deviance <- last_deviance
    last_deviance <- fit$deviance
    if (sum((eta - last_eta)^2) / sum(eta^2) <= epsilon) {
      converged <- TRUE
      break
    }
  }
  if (!converged && !step_back) {
    warning("Algorithm did not converge", call. = FALSE)
  }
  if (step_back) {
    warning("Algorithm stopped without convergence", call. = FALSE)
  }
  eps <- 10 * .Machine$double.eps
  if (any(prb < eps) || any(1 - prb < eps)) {
    warning("Fitted probabilities numerically 0 or 1 occurred", call. = FALSE)
  }
  # The check of mclogit, on the deviance scale of the augmented rows:
  dev_const <- -2 * sum(Y[Y > 0] * log(y_rel[Y > 0]))
  if (fit$deviance + dev_const < 0) {
    warning("Approximate deviance is negative.\nYou might be overfitting ",
            "your data or the group size is too small.", call. = FALSE)
  }

  coefs <- setNames(drop(parms$alpha), colnames(XD))
  VarCov <- setNames(lapply(seq_along(d), function(k) {
    matrix(parms$Phi[[k]], d[k], d[k],
           dimnames = list(vc_nms[[k]], vc_nms[[k]]))
  }), grp_nms)
  random_effects <- lapply(z_idx, function(idx) {
    matrix(parms$b[idx], ncol = 1, dimnames = list(colnames(Z)[idx], NULL))
  })
  names(random_effects) <- NULL
  out <- list(
    coefficients = coefs,
    coefmat = matrix(coefs, nrow = ncol(D),
                     dimnames = list("Logit eqn." = colnames(D),
                                     "Predictors" = colnames(X))),
    random.effects = random_effects,
    VarCov = VarCov,
    lambda = parms$lambda,
    linear.predictors = eta,
    deviance = fit$deviance,
    iter = iter,
    converged = converged,
    method = "PQL",
    estimator = "ML",
    formula = formula,
    terms = mt,
    random = random,
    groups = groups,
    contrasts = contrasts,
    xlevels = xlevels,
    D = D,
    N = sum(w_obs),
    response.type = "matrix",
    catCov = "free"
  )
  class(out) <- c("mmblogit", "mblogit", "mmclogit", "mclogit", "lm")
  return(out)
}

# Category probabilities from the linear predictors (`s`: the observation of
# each element).
categ_pql_probs <- function(eta, s) {
  expeta <- exp(eta - stats::ave(eta, s, FUN = max))
  expeta / rowsum(expeta, s, reorder = FALSE)[s]
}

# One inner fit: the variance parameters by nlminb() for the working response
# `y_star` with the working weights of `prb`, then the coefficients.
categ_pql_inner <- function(parms, y_star, X, Z, z_idx, d, m, w, prb, s,
                            maxit_inner) {
  p <- ncol(X)
  q <- ncol(Z)
  A <- cbind(X, Z, y_star)
  # W A, with W = diag(w prb) - (sqrt(w) prb)(sqrt(w) prb)' within each
  # observation:
  swp <- sqrt(w) * prb
  WA <- w * prb * A -
    swp * rowsum(swp * A, s, reorder = FALSE)[s, , drop = FALSE]
  G <- crossprod(A, WA)
  ix <- seq_len(p)
  iz <- p + seq_len(q)
  iy <- p + q + 1L
  XWX <- G[ix, ix, drop = FALSE]
  XWy <- G[ix, iy]
  ZWX <- G[iz, ix, drop = FALSE]
  ZWZ <- G[iz, iz, drop = FALSE]
  ZWy <- G[iz, iy]
  yWy <- G[iy, iy]

  alpha_start <- parms$alpha
  if (!length(alpha_start)) {
    alpha_start <- categ_pql_numeric(solve(XWX, XWy))
  }
  Phi_start <- parms$Phi
  if (!length(Phi_start)) {
    r <- drop(y_star - X %*% alpha_start)
    Phi_start <- categ_pql_numeric(lapply(seq_along(d), function(k) {
      Zk <- Z[, z_idx[[k]], drop = FALSE]
      b_k <- solve(crossprod(Zk) + diag(ncol(Zk)), crossprod(Zk, r))
      S_k <- tcrossprod(matrix(b_k, nrow = d[k]))
      if (qr(S_k)$rank < d[k]) S_k <- diag(diag(S_k), nrow = d[k])
      S_k / (m[k] - 1)
    }))
  }
  lambda_start <- categ_pql_numeric(unlist(lapply(Phi_start, function(Phi_k) {
    L <- chol(solve(Phi_k))
    L[upper.tri(L, diag = TRUE)]
  })))
  pll <- categ_pql_objective(XWX = XWX, XWy = XWy, ZWX = ZWX, ZWZ = ZWZ,
                             ZWy = ZWy, yWy = yWy, z_idx = z_idx, d = d, m = m)
  devfunc <- function(lambda) -2 * pll(lambda)$logLik
  gradfunc <- function(lambda) -2 * pll(lambda, gradient = TRUE)$gradient
  opt <- categ_pql_numeric(stats::nlminb(
    lambda_start, objective = devfunc, gradient = gradfunc,
    control = list(iter.max = maxit_inner)
  ))
  if (opt$convergence != 0) {
    warning(sprintf("Inner iterations did not converge - nlminb message: %s",
                    opt$message), call. = FALSE)
  }
  res <- categ_pql_numeric(pll(opt$par))
  if (!is.finite(res$logLik)) {
    categ_pql_numeric(stop("the inner PQL objective is not finite"))
  }
  return(list(alpha = res$alpha, b = res$b, Psi = res$Psi,
              Phi = categ_pql_numeric(lapply(res$Psi, solve)),
              lambda = opt$par, log_det_iSigma = res$log_det_iSigma,
              log_det_H = res$log_det_H))
}

# Evaluate `expr` and turn an error into an error of class
# `categ_pql_numeric_error`, which the outer PQL loop handles (numerical
# failures such as a singular matrix).
categ_pql_numeric <- function(expr) {
  tryCatch(expr, error = function(e) {
    stop(structure(
      class = c("categ_pql_numeric_error", "error", "condition"),
      list(message = conditionMessage(e), call = NULL)
    ))
  })
}

# The pseudo log-likelihood of the inner fit as a function of `lambda` (the
# upper triangles of the Cholesky factors of the group-level precision
# matrices), with its gradient.
categ_pql_objective <- function(XWX, XWy, ZWX, ZWZ, ZWy, yWy, z_idx, d, m) {
  n_lambda <- d * (d + 1) / 2
  lambda_idx <- split(seq_len(sum(n_lambda)), rep(seq_along(d), n_lambda))
  lambda_to_mats <- function(lambda) {
    lapply(seq_along(d), function(k) {
      L <- matrix(0, d[k], d[k])
      L[upper.tri(L, diag = TRUE)] <- lambda[lambda_idx[[k]]]
      L
    })
  }
  # nlminb() asks for the objective and the gradient at the same point, so
  # keep the last result:
  cache <- new.env(parent = emptyenv())
  pll <- function(lambda, gradient = FALSE) {
    if (identical(lambda, cache$lambda) &&
        (!gradient || !is.null(cache$res$gradient))) {
      return(cache$res)
    }
    Lambda <- lambda_to_mats(lambda)
    Psi <- lapply(Lambda, crossprod)
    # H = Z'WZ + Sigma^-1, with Sigma^-1 block diagonal (Psi_k for each level
    # of group k):
    H <- ZWZ
    for (k in seq_along(d)) {
      for (j in seq_len(m[k])) {
        ii <- z_idx[[k]][(j - 1L) * d[k] + seq_len(d[k])]
        H[ii, ii] <- H[ii, ii] + Psi[[k]]
      }
    }
    R <- tryCatch(chol(H), error = function(e) NULL)
    if (is.null(R)) {
      res <- list(logLik = -Inf, gradient = rep(NA_real_, length(lambda)))
    } else {
      RA <- backsolve(R, ZWX, transpose = TRUE)
      Ry <- backsolve(R, ZWy, transpose = TRUE)
      XiVX <- XWX - crossprod(RA)
      XiVX <- (XiVX + t(XiVX)) / 2
      XiVy <- XWy - drop(crossprod(RA, Ry))
      alpha <- solve(XiVX, XiVy)
      b <- drop(backsolve(R, Ry - RA %*% alpha))
      quad <- yWy - sum(XWy * alpha) - sum(ZWy * b)
      log_det_iSigma <- sum(vapply(seq_along(d), function(k) {
        2 * m[k] * sum(log(abs(diag(Lambda[[k]]))))
      }, numeric(1)))
      log_det_H <- 2 * sum(log(diag(R)))
      res <- list(logLik = (log_det_iSigma - log_det_H - quad) / 2,
                  alpha = alpha, b = b, Psi = Psi,
                  log_det_iSigma = log_det_iSigma, log_det_H = log_det_H)
      if (gradient) {
        K <- chol2inv(R)
        res$gradient <- unlist(lapply(seq_along(d), function(k) {
          dk <- d[k]
          mk <- m[k]
          # Sum of the mk diagonal dk x dk blocks of K = H^-1:
          K_kk <- K[z_idx[[k]], z_idx[[k]], drop = FALSE]
          offs <- rep((seq_len(mk) - 1L) * dk, each = dk * dk)
          i <- rep(rep(seq_len(dk), dk), mk) + offs
          j <- rep(rep(seq_len(dk), each = dk), mk) + offs
          SumK <- matrix(rowSums(matrix(K_kk[cbind(i, j)], nrow = dk * dk)),
                         dk, dk)
          S <- tcrossprod(matrix(b[z_idx[[k]]], nrow = dk))
          Gr <- Lambda[[k]] %*% (mk * solve(Psi[[k]]) - SumK - S)
          Gr[upper.tri(Gr, diag = TRUE)]
        }))
      }
    }
    assign("lambda", lambda, envir = cache)
    assign("res", res, envir = cache)
    return(res)
  }
  return(pll)
}

# Linear predictors, probabilities and the approximate deviance of mclogit.
categ_pql_eval <- function(parms, X, Z, z_idx, d, y, w, s) {
  eta <- drop(X %*% parms$alpha + Z %*% parms$b)
  rand_ssq <- sum(vapply(seq_along(d), function(k) {
    B <- matrix(parms$b[z_idx[[k]]], nrow = d[k])
    sum(B * (parms$Psi[[k]] %*% B))
  }, numeric(1)))
  prb <- categ_pql_probs(eta, s)
  dev_resids <- ifelse(y > 0, 2 * w * y * (log(y) - log(prb)), 0)
  deviance <- -parms$log_det_iSigma + parms$log_det_H + sum(dev_resids) +
    rand_ssq
  return(list(eta = eta, prb = prb, deviance = deviance))
}
