fetch_data <- function(data, obs = NULL, newdata = NULL) {
  if (is.null(obs)) {
    if (is.null(newdata)) {
      return(data)
    } else {
      return(newdata)
    }
  } else if (is.null(newdata)) {
    return(data[obs, , drop = FALSE])
  } else {
    return(newdata[obs, , drop = FALSE])
  }
}

linear_mle <- function(formula, data, weights = NULL, regul = NULL) {
  formula <- validate_response_formula(formula)
  fit_lm_ridge_callback <- function(f) {
    if (count_terms_in_subformula(f) == 1) {
      lm(f, data = data, weights = weights)
    } else {
      fit <- lm.ridge(f, data = data, weights = weights, lambda = regul)
      fit$data <- data
      fit$formula <- f
      fit$weights <- weights
      fit
    }
  }
  if (inherits(formula, "formula")) {
    return(fit_lm_ridge_callback(formula))
  } else if (inherits(formula, "list")) {
    return(lapply(formula, function(f) (fit_lm_ridge_callback(f))))
  } else {
    stop("The provided formula is neither a formula object nor a list")
  }
}

#' Use lmer to fit the projection to the posterior draws for multilevel models.
#' Note that we don't use glmer because the target is a pseudo-Gaussian
#' transformation.
linear_multilevel_mle <- function(formula, data, weights = NULL, regul = NULL) {
  formula <- validate_response_formula(formula)
  fit_lmer_callback <- function(f) {
    tryCatch(lme4::lmer(f, data = data, weights = weights),
      error = function(e) {
        if (grepl("No random effects", as.character(e))) {
          lm.ridge(f, data = data, weights = weights, lambda = regul)
        } else if (grepl("not positive definite", as.character(e))) {
          lme4::lmer(f,
            data = data, weights = weights,
            control = lmerControl(
              optimizer = "optimx",
              optCtrl = list(method = "nlminb")
            )
          )
        } else {
          e
        }
      }
    )
  }
  if (inherits(formula, "formula")) {
    return(fit_lmer_callback(formula))
  } else if (inherits(formula, "list")) {
    return(lapply(formula, function(f) (fit_lmer_callback(f))))
  } else {
    stop("The provided formula is neither a formula object nor a list")
  }
}

linear_multilevel_proj_predfun <- function(fit, newdata = NULL, weights = NULL) {
  if (is.null(weights)) {
    weights <- 1
  }
  if (inherits(fit, "list")) {
    return(do.call(cbind, lapply(fit, function(fit) {
      predict(fit, newdata = newdata, allow.new.levels = TRUE, weights = weights)
    })))
  } else {
    return(predict(fit,
      newdata = newdata, allow.new.levels = TRUE,
      weights = weights
    ))
  }
}

linear_proj_predfun <- function(fit, newdata = NULL, weights = NULL) {
  if (is.null(weights)) {
    weights <- 1
  }
  if (inherits(fit, "list")) {
    if (!is.null(newdata)) {
      return(do.call(cbind, lapply(fit, function(fit) {
        predict(fit, newdata = newdata, weights = weights)
      })))
    } else {
      return(do.call(cbind, lapply(fit, function(fit) {
        predict(fit)
      })))
    }
  }
  else {
    if (!is.null(newdata)) {
      return(predict(fit, newdata = newdata, weights = weights))
    } else {
      return(predict(fit))
    }
  }
}

# importing MASS solves a weird error caused by coef(fit) sometimes
# returning NULL for some reason that seems platform dependent
#' @import MASS
predict.ridgelm <- function(fit, newdata = NULL, weights = NULL) {
  b <- coef(fit)
  center <- fit$xm
  scales <- fit$scales
  if (!is.null(newdata)) {
    if (is.null(weights)) {
      weights <- 1
    }
    x <- model.matrix(delete.response(terms(fit$formula)), newdata)
  } else {
    if (is.null(fit$weights)) {
      weights <- 1
    } else {
      weights <- fit$weights
    }
    x <- model.matrix(delete.response(terms(fit$formula)), fit$data)
  }
  if (NCOL(x) > 1) {
    x <- cbind(x[, 1, drop = FALSE],
               scale(x[, -1, drop = FALSE], center = center, scale = scales))
  }
  x <- weights * x
  return(x %*% b)
}

## taken from MASS, but fix a small bug when removing the intercept
lm.ridge <- function(formula, data, subset, na.action, lambda = 0, model = FALSE,
                     x = FALSE, y = FALSE, contrasts = NULL, ...) {
  m <- match.call(expand.dots = FALSE)
  m$model <- m$x <- m$y <- m$contrasts <- m$... <- m$lambda <- NULL
  m[[1L]] <- quote(stats::model.frame)
  m <- eval.parent(m)
  Terms <- attr(m, "terms")
  Y <- model.response(m)
  X <- model.matrix(Terms, m, contrasts)
  n <- nrow(X)
  p <- ncol(X)
  offset <- model.offset(m)
  if (!is.null(offset)) {
    Y <- Y - offset
  }
  if (Inter <- attr(Terms, "intercept")) {
    Xm <- colMeans(X[, -Inter, drop = FALSE])
    Ym <- mean(Y)
    p <- p - 1
    X <- X[, -Inter, drop = FALSE] - rep(Xm, rep(n, p))
    Y <- Y - Ym
  }
  else {
    Ym <- Xm <- NA
  }
  Xscale <- drop(rep(1 / n, n) %*% X^2)^0.5
  X <- X / rep(Xscale, rep(n, p))
  Xs <- svd(X)
  rhs <- t(Xs$u) %*% Y
  d <- Xs$d
  lscoef <- Xs$v %*% (rhs / d)
  lsfit <- X %*% lscoef
  resid <- Y - lsfit
  s2 <- sum(resid^2) / (n - p - Inter)
  HKB <- (p - 2) * s2 / sum(lscoef^2)
  LW <- (p - 2) * s2 * n / sum(lsfit^2)
  k <- length(lambda)
  dx <- length(d)
  div <- d^2 + rep(lambda, rep(dx, k))
  a <- drop(d * rhs) / div
  dim(a) <- c(dx, k)
  coef <- Xs$v %*% a
  dimnames(coef) <- list(names(Xscale), format(lambda))
  GCV <- colSums((Y - X %*% coef)^2) / (n - colSums(matrix(
    d^2 / div,
    dx
  )))^2
  res <- list(
    coef = drop(coef), scales = Xscale, Inter = Inter,
    lambda = lambda, ym = Ym, xm = Xm, GCV = GCV, kHKB = HKB,
    kLW = LW
  )
  class(res) <- "ridgelm"
  res
}
