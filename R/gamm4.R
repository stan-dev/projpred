## taken from gam4
gamm4.setup <- function(formula, pterms, data = NULL, knots = NULL) {
  ## first simply call `gam.setup'....
  G <- mgcv:::gam.setup(formula, pterms,
      data = data, knots = knots, sp =
          NULL, min.sp = NULL, H = NULL, absorb.cons = TRUE, sparse.cons = 0,
      gamm.call = TRUE
  )
  if (!is.null(G$L)) {
    stop(
      "gamm can not handle linked smoothing parameters",
      "(probably from use of `id' or adaptive smooths)"
    )
  }

  first.f.para <- G$nsdf + 1

  random <- list()

  if (G$nsdf > 0) ind <- 1:G$nsdf else ind <- rep(0, 0)
  X <- G$X[, ind, drop = FALSE] # accumulate fixed effects into here

  xlab <- rep("", 0)
  ## sparse version of full matrix, treating smooths as fixed
  G$Xf <- as(X, "dgCMatrix")
  first.para <- G$nsdf + 1
  used.names <- names(data) ## keep track of all variable names already used

  if (G$m) {
    for (i in 1:G$m) { ## work through the smooths
      sm <- G$smooth[[i]]
      sm$X <- G$X[, sm$first.para:sm$last.para, drop = FALSE]
      ## convert smooth to random effect and fixed effects
      rasm <- mgcv::smooth2random(sm, used.names, type = 2)
      used.names <- c(used.names, names(rasm$rand))

      sm$fixed <- rasm$fixed

      ## deal with creation of sparse full model matrix
      if (!is.null(sm$fac)) {
        flev <- levels(sm$fac) ## grouping factor for smooth
        n.lev <- length(flev)
        for (k in 1:n.lev) {
          G$Xf <- cbind2(G$Xf, as(
              sm$X * as.numeric(sm$fac == flev[k]),
              "dgCMatrix"
          ))
        }
      } else {
        n.lev <- 1
        G$Xf <- cbind2(G$Xf, as(sm$X, "dgCMatrix"))
      }

      ## now append random effects to main list
      n.para <- 0 ## count random coefficients
      if (!sm$fixed) {
        for (k in 1:length(rasm$rand)) n.para <- n.para + ncol(rasm$rand[[k]])
        sm$lmer.name <- names(rasm$rand)
        random <- c(random, rasm$rand)
        sm$trans.D <- rasm$trans.D
        sm$trans.U <- rasm$trans.U ## matrix mapping fit coefs back to original
      }

      ## ensure stored first and last para relate to G$Xf in expanded version
      sm$last.para <- first.para + ncol(rasm$Xf) + n.para - 1
      sm$first.para <- first.para
      first.para <- sm$last.para + 1

      if (ncol(rasm$Xf)) {
        Xfnames <- rep("", ncol(rasm$Xf))
        k <- length(xlab) + 1
        for (j in 1:ncol(rasm$Xf)) {
          xlab[k] <- Xfnames[j] <-
            mgcv::new.name(paste(sm$label, "Fx", j, sep = ""), xlab)
          k <- k + 1
        }
        colnames(rasm$Xf) <- Xfnames
      }

      X <- cbind(X, rasm$Xf) # add fixed model matrix to overall fixed X

      sm$first.f.para <- first.f.para
      first.f.para <- first.f.para + ncol(rasm$Xf)
      ## note less than sm$first.f.para => no fixed
      sm$last.f.para <- first.f.para - 1

      ## store indices of random parameters in smooth specific array
      sm$rind <- rasm$rind
      sm$rinc <- rasm$rinc
      ## pen.ind==i TRUE for coef penalized by ith penalty
      sm$pen.ind <- rasm$pen.ind

      sm$n.para <- n.para

      sm$X <- NULL ## delete model matrix

      G$smooth[[i]] <- sm ## replace smooth object with extended version
    }
  }

  G$random <- random ## named list of random effect matrices
  G$X <- X ## fixed effects model matrix

  G
}

## refactored from gamm4 to return the model matrix for generating predictions
## with fit$mer and newdata
model.matrix.gamm4 <- function(formula, random = NULL, data = NULL,
                               family = gaussian()) {
  if (!is.null(random)) {
    if (!inherits(random, "formula")) {
      stop("gamm4 requires `random' to be a formula")
    }
    random.vars <- all.vars(random)
  } else {
    random.vars <- NULL
  }
  # create model frame.....
  gp <- mgcv::interpret.gam(formula) # interpret the formula
  mf <- match.call(expand.dots = FALSE)
  mf$formula <- gp$fake.formula
  mf$REML <- mf$verbose <- mf$control <- mf$start <- mf$family <- mf$scale <-
    mf$knots <- mf$random <- mf$... <- NULL ## mf$weights?
  mf[[1]] <- as.name("model.frame")
  pmf <- mf
  gmf <- eval(mf, parent.frame())
  gam.terms <- attr(gmf, "terms")

  if (length(random.vars)) {
    mf$formula <- as.formula(paste(paste(deparse(gp$fake.formula,
      backtick = TRUE
    ), collapse = ""), "+", paste(random.vars,
      collapse = "+"
    )))
    mf <- eval(mf, parent.frame())
  } else {
    mf <- gmf
  }
  rm(gmf)

  if (nrow(mf) < 2) {
    stop("Not enough (non-NA) data to do anything meaningful")
  }

  ## summarize the *raw* input variables
  ## note can't use get_all_vars here -- buggy with matrices
  vars <- all.vars(gp$fake.formula[-2]) ## drop response here
  inp <- parse(text = paste("list(", paste(vars, collapse = ","), ")"))
  dl <- eval(inp, data, parent.frame())
  names(dl) <- vars ## list of all variables needed
  ## summarize the input data
  var.summary <- mgcv:::variable.summary(gp$pf, dl, nrow(mf))

  ## lmer offset handling work around...
  ## variables not in mf raw -- can cause lmer problem
  mvars <- vars[!vars %in% names(mf)]
  if (length(mvars) > 0) {
    for (i in 1:length(mvars)) {
      mf[[mvars[i]]] <- dl[[mvars[i]]]
    }
  } ## append raw versions to mf
  pmf$formula <- gp$pf
  pmf <- eval(pmf, parent.frame()) # pmf contains all data for non-smooth part
  pTerms <- attr(pmf, "terms")

  G <- gamm4.setup(gp, pterms = pTerms, data = mf)
  G$var.summary <- var.summary

  ## number of random smooths (i.e. s(...,fx=FALSE,...) terms)
  n.sr <- length(G$random)

  if (is.null(random) && n.sr == 0) {
    return(mf)
  }

  yname <- mgcv::new.name("y", names(mf))
  eval(parse(text = paste("mf$", yname, "<-G$y", sep = "")))
  Xname <- mgcv::new.name("X", names(mf))
  eval(parse(text = paste("mf$", Xname, "<-G$X", sep = "")))

  offset.name <- attr(mf,"names")[attr(attr(mf,"terms"),"offset")]
  lme4.formula <- paste(yname, "~", Xname, "-1")
  if (length(offset.name)) {
    lme4.formula <- paste(lme4.formula, "+", offset.name)
  }

  ## Add the random effect dummy variables for the smooth
  r.name <- names(G$random)
  if (n.sr) {
    ## adding the constructed variables to the model frame avoiding name
    ## duplication
    for (i in 1:n.sr) {
      mf[[r.name[i]]] <- factor(rep(1:ncol(G$random[[i]]),
        length = nrow(G$random[[i]])
      ))
      lme4.formula <- paste(lme4.formula, "+ (1|", r.name[i], ")")
    }
  }
  if (!is.null(random)) { ## append the regular random effects
    lme4.formula <- paste(
      lme4.formula, "+",
      substring(paste(deparse(random, backtick = TRUE), collapse = ""),
        first = 2
      )
    )
  }
  lme4.formula <- as.formula(lme4.formula)
  if (family$family == "gaussian" && family$link == "identity")
    linear <- TRUE else linear <- FALSE

  control <- if (linear) lme4::lmerControl()
             else lme4::glmerControl()

  ## NOTE: further arguments should be passed here...
  b <- if (linear) {
    lme4::lFormula(lme4.formula,
      data = mf, weights = G$w, REML = TRUE,
      control = control,
    )
  } else {
    lme4::glFormula(lme4.formula,
      data = mf, family = family, weights = G$w,
      control = control,
    )
  }
  if (n.sr) {
    tn <- names(b$reTrms$cnms)
    ind <- 1:length(tn)
    sn <- names(G$random) ## names of smooth random components
    for (i in 1:n.sr) { ## loop through random effect smooths
      k <- ind[sn[i] == tn] ## which term should contain G$random[[i]]
      ii <- (b$reTrms$Gp[k] + 1):b$reTrms$Gp[k + 1]
      b$reTrms$Zt[ii, ] <- as(t(G$random[[i]]), "dgCMatrix")
      b$reTrms$cnms[[k]] <- attr(G$random[[i]], "s.label")
    }
  }

  return(nlist(mf, b))
}
