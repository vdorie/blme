optWrapWrap <- function(
  optimizer,
  fn,
  par,
  lower=-Inf,
  upper=Inf,
  control=list(),
  adj=FALSE,
  calc.derivs=TRUE,
  force.calc.derivs=FALSE,
  use.last.params=FALSE,
  verbose=0L
) {
  ## for backward compatibility, only call with valid arguments

  optwrap <- get("optwrap", getNamespace("lme4"))
  args <- list(
    optimizer,
    fn,
    par,
    lower=lower,
    control=control,
    adj=adj,
    verbose=verbose
  )

  optionalArgs <- names(formals(optwrap))
  optionalArgs <- optionalArgs[sapply(formals(optwrap), function(x) x != '')]

  missingArgs <- optionalArgs[optionalArgs %not_in% names(args)]
  for (missingArg in missingArgs)
    args[[missingArg]] <- get(missingArg)
  do.call(optwrap, args)
}

## derived from lme4/R/modular.R

optimizeLmer <- function(devfun,
                         optimizer    = formals(lmerControl)$optimizer,
                         restart_edge = formals(lmerControl)$restart_edge,
                         boundary.tol = formals(lmerControl)$boundary.tol,
                         start   = NULL,
                         verbose = 0L,
                         control = list(),
                         ...) {
  verbose <- as.integer(verbose)
  rho <- environment(devfun)

  lme4Version <- packageVersion("lme4")
  lme4Namespace <- getNamespace("lme4")
  
  parInfo <- rho$parInfo
  startingValues <- getStartingValues(start, rho, parInfo)
  lowerBounds <- getBounds(parInfo, "lower")
  upperBounds <- if (lme4Version >= "2.0.0") getBounds(parInfo, "upper") else rep(Inf, length(lowerBounds))
  ## b/c bounds are pulled from rho to check convergence
  rho$lower <- lowerBounds
  rho$upper <- upperBounds
  thetaLowerBounds <- lowerBounds[seq_along(rho$pp$theta)]
  thetaUpperBounds <- upperBounds[seq_along(rho$pp$theta)]

  opt <- optWrapWrap(
    optimizer,
    devfun,
    startingValues,
    lower=lowerBounds,
    upper=upperBounds,
    control=control,
    adj=FALSE,
    verbose=verbose,
    ...
  )
  
  if (restart_edge) {
    ## FIXME: should we be looking at rho$pp$theta or opt$par
    ##  at this point???  in koller example (for getData(13)) we have
    ##   rho$pp$theta=0, opt$par=0.08
    theta0 <- if (lme4Version >= "2.0.0") rho$mkPar(rho$pp$theta) else rho$pp$theta
    par0 <- opt$par
    par0[seq_along(theta0)] <- theta0

    if (any(
      length(wl <- which(par0 == lowerBounds)) > 0L |
      length(wu <- which(par0 == upperBounds)) > 0L
    ))
    {
      d0 <- devfun(par0)
      btol <- 1e-5  ## FIXME: make user-settable?
      bgrad <- mapply(
        function(i, bval, btol) {
          par <- par0
          par[i] <- bval + btol
          (devfun(par) - d0)/btol
        },
        i = c(wl, wu),
        bval = c(lowerBounds[wl], upperBounds[wu]),
        btol = rep(c(btol, -btol),
        c(length(wl), length(wu)))
      )
      ## what do I need to do to reset rho$pp$theta to original value???
      devfun(par0) ## reset rho$pp$theta after tests
      ## FIXME: allow user to specify ALWAYS restart if on boundary?
      if (any(is.na(bgrad))) {
        warning("some gradient components are NA near boundaries, skipping boundary check")
        return(opt)
      } else {
        if (any(bgrad < 0)) {
          if (verbose) message("some theta parameters on the boundary, restarting")
          opt <- optWrapWrap(
            optimizer,
            devfun,
            opt$par,
            lower=lowerBounds,
            upper=upperBounds,
            control=control,
            adj=FALSE,
            verbose=verbose,
            ...
          )
        }
      } ## bgrad not NA
    }
  }
  if (!is.null(boundary.tol) && boundary.tol > 0) {
    if (exists("check.boundary", lme4Namespace))
      opt <- get("check.boundary", lme4Namespace)(rho, opt, devfun, boundary.tol)
  }
  
  opt
}

