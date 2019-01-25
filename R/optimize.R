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

  lme4Env <- asNamespace("lme4")
  
  parInfo <- rho$parInfo
  startingValues <- getStartingValues(start, rho, parInfo)
  lowerBounds <- getLowerBounds(parInfo)
  rho$lower <- lowerBounds ## b/c lower bounds are pulled from devfunenv to check convergence
  thetaLowerBounds <- lowerBounds[seq_along(rho$pp$theta)]
  
  optwrap <- get("optwrap", lme4Env)
  lme4IsOld <- is.null(formals(optwrap)[["calc.derivs"]])
  opt <-
    if (!lme4IsOld)
      optwrap(optimizer, devfun, startingValues, lower = lowerBounds, control = control, adj = FALSE, verbose = verbose, ...)
    else
      optwrap(optimizer, devfun, startingValues, lower = lowerBounds, control = control, adj = FALSE, verbose = verbose)
  
  if (restart_edge) {
    ## FIXME: should we be looking at rho$pp$theta or opt$par
    ##  at this point???  in koller example (for getData(13)) we have
    ##   rho$pp$theta=0, opt$par=0.08
    if (length(bvals <- which(rho$pp$theta == thetaLowerBounds)) > 0L) {
      par <- opt$par
      ## *don't* use numDeriv -- cruder but fewer dependencies, no worries
      ##  about keeping to the interior of the allowed space
      theta0 <- new("numeric", rho$pp$theta) ## 'deep' copy ...
      d0 <- devfun(par)
      btol <- 1e-5  ## FIXME: make user-settable?
      bgrad <- sapply(bvals,
                      function(i) {
                        bndval <- rho$lower[i]
                        par[seq_along(theta0)] <- theta0
                        par[i] <- bndval + btol
                        (devfun(par) - d0) / btol
                      })
      ## what do I need to do to reset rho$pp$theta to original value???
      par[seq_along(theta0)] <- theta0
      devfun(par) ## reset rho$pp$theta after tests
      ## FIXME: allow user to specify ALWAYS restart if on boundary?
      if (any(bgrad < 0)) {
        if (verbose) message("some theta parameters on the boundary, restarting")
        opt <- if (!lme4IsOld)
          optwrap(optimizer, devfun, opt$par, lower = lowerBounds, control = control, adj = FALSE, verbose = verbose, ...)
        else
          optwrap(optimizer, devfun, opt$par, lower = lowerBounds, control = control, adj = FALSE, verbose = verbose)
      }
    }
  }
  if (!is.null(boundary.tol) && boundary.tol > 0) {
    if (exists("check.boundary", lme4Env))
      opt <- get("check.boundary", lme4Env)(rho, opt, devfun, boundary.tol)
  }
  
  opt
}
