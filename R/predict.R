predict.varsel <- function(fit, obj, size, newdata, ...) {
  predict(project(fit, obj, size), newdata, ...)
}

predict.glmproj <- function(q, newdata) {
  if(!is.list(newdata)) newdata <- list(x = newdata)
  if(!is.null(newdata$w)) newdata$x <- newdata$x*newdata$w
  eta <- newdata$x[, q$chosen, drop = F]%*%q$b
  if(!is.null(newdata$offset)) eta <- eta + newdata$offset

  drop(q$family$linkinv(eta))
}
