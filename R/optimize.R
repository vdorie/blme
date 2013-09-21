foldPars <- function(pars) {
  if (is.numeric(pars)) return(pars);
  as.numeric(unlist(pars));
}

foldLowers <- function(pars) {
  if (!is.list(pars) && is.numeric(pars)) return(attr(pars, "lower"));
  unlist(sapply(pars, function(par) attr(par, "lower")));
}

expandPars <- function(parsVector, parsList) {
  if (!is.list(parsList)) return(invisible(NULL));
  
  parentEnv <- parent.frame();
  parNames <- names(parsList);

  index <- 0;
  for (i in 1:length(parsList)) {
    parLength <- length(parsList[[i]]);
    parName <- parNames[[i]];

    parentEnv[[parName]] <- parsVector[index + 1:parLength];
    index <- index + parLength;
  }
  return(invisible(NULL));
}

padStartList <- function(start, rho)
{
  if (is(rho$resp, "lmerResp") && rho$blmerControl$sigmaOptimizationType == SIGMA_OPTIM_NUMERIC) {
    if (is.null(start)) return(list(theta = numeric(), sigma = numeric()));
    if (is.numeric(start)) return(list(theta = start, sigma = numeric()));
    if (is.list(start)) {
      if (is.null(start$sigma)) {
        start$sigma <- numeric();
        return(start);
      }
      return(start);
    }
    stop("start argument must be numeric vector or list");
  }
  return(start);
}

setLowerBounds <- function(pars, rho)
{
  if (!is.list(pars)) {
    attr(pars, "lower") <- rho$lower;
    return(pars);
  }
  numTheta <- length(pars$theta);
  attr(pars$theta, "lower") <- rho$lower[1:numTheta];
  if (!is.null(pars$fixef)) {
    numFixef <- length(pars$fixef);
    attr(pars$fixef, "lower") <- rho$lower[numTheta + 1:numFixef];
  }
  if (is(rho$resp, "lmerResp") && rho$blmerControl$sigmaOptimizationType == SIGMA_OPTIM_NUMERIC) {
    attr(pars$sigma, "lower") <- 0;
  }
  return(pars);
}

getStartingValues <- function(start, rho) {
  start <- padStartList(start, rho);
  start <- getStart(start, rho);
  start <- setLowerBounds(start, rho);
  rho$pars <- start;
  
  return(start);
}

getStart <- function(pars, rho) {
  pp <- rho$pp; resp <- rho$resp;
  
  getDefault <- function(par, parName) {
    if (!is.numeric(par)) stop("start for ", parName, " must be a numeric vector");
    
    defaultFunction <- attr(par, "default");
    if (!is.null(defaultFunction)) {
      default <- defaultFunction();
      attributes(default) <- attributes(par);
      if (length(par) == 0) {
        par <- default;
      } else if (length(par) != length(default)) {
        stop("incorrect number of ", parName, " components (!=", length(default), ")");
      }
    }

    if (length(par) == 0) stop("could not find starting values for ", parName);
    return(par);
  }

  if (is.null(pars)) pars <- numeric(0);
  
  if (!is.list(pars)) {
    defaultFunction <- attr(pars, "default");
    if (is.null(defaultFunction)) {
      defaultFunction <- function() pp$theta;
      environment(defaultFunction)$pp <- pp;
      attr(pars, "default") <- defaultFunction;
    }

    return(getDefault(pars, "theta"));
  }

  if (any(is.null(names(pars)))) stop("par list requires elements to be named");
  if (length(unique(names(pars))) != length(pars)) stop("par list contains duplicate elements");
  
  for (i in 1:length(pars)) {
    par <- pars[[i]];
    parName <- names(pars)[[i]];
    
    defaultFunction <- attr(par, "default");
    if (is.null(defaultFunction)) {
      if (parName == "theta") {
        defaultFunction <- function() pp$theta;
        environment(defaultFunction)$pp <- pp;
        attr(par, "default") <- defaultFunction;
      } else if (parName == "beta") {
        defaultFunction <- function() pp$delb;
        environment(defaultFunction)$pp <- pp;
        attr(par, "default") <- defaultFunction;
      } else if (parName == "sigma") {
        defaultFunction <- function() sd(resp$y);
        environment(defaultFunction)$resp <- resp;
        attr(par, "default") <- defaultFunction;
      }
    }
    
    pars[[i]] <- getDefault(par, parName);
  }

  return (pars);
}

setLower <- function(pars, lowers) {
  if (!is.list(lowers)) {
    if (!is.numeric(lowers)) stop("supplied lower bounds not numeric or list");
    
    if (!is.list(pars)) {
      if (!is.numeric(pars)) stop("expected numeric vector for theta");
      if (length(pars) != length(lowers)) stop("lower bounds for theta not of length equal to theta");
      
      attr(pars, "lower") <- lowers;
      return(pars);
    } else {
      if (is.null(pars$theta)) stop("expected theta to be part of pars");
      if (length(pars$theta) != length(lowers)) stop("lower bounds for theta not of length equal to theta");
      
      attr(pars$theta, "lower") <- lowers;
      return(pars);
    }
  } else {
    if (!is.list(pars)) stop("multiple lower bound vectors supplied, but only one parameter vector (theta)");
    
    if (length(lowers) == length(pars)) {
      nullNames <- is.null(names(lowers));
      if (any(nullNames))
        names(lowers)[nullNames] <- names(pars)[nullNames];
    }
    lowerNames <- names(lowers);
    parNames <- names(pars);
    for (i in 1:length(lowers)) {
      lower <- lowers[[i]];
      lowerName <- lowerNames[[i]];
      
      if (!is.numeric(lower)) stop("supplied lower bound for ", lowerName, " not numeric or list");

      matchRows <- lowerName == parNames;
      if (!any(matchRows)) stop("supplied lower bound for ", lowerName, " does not match any named parameter vector");

      matchRow <- min(which(matchRows));

      if (length(pars[[matchRow]]) != length(lower)) stop("lower bounds for theta not of length equal to theta");
      
      attr(pars[[matchRow]], "lower") <- lower;
    }
  }

  return(pars);
}

##' @rdname modular
##' @inheritParams lmer
##' @inheritParams lmerControl
##' @param devfun a deviance function, as generated by \code{\link{mkLmerDevfun}}
##' @return \bold{optimizeLmer}: Results of an optimization.
##' \cr
##' \cr
##' @export
optimizeLmer <- function(devfun,
                         optimizer="Nelder_Mead",
                         restart_edge=FALSE,
                         start = NULL,
                         verbose = 0L,
                         control = list()) {
    verbose <- as.integer(verbose)
    rho <- environment(devfun)

    if (is.null(start)) {
      start <- getStart(start, rho);
      start <- setLowerBounds(start, rho);
    }

    startingValues <- foldPars(start);
    lowerBounds <- foldLowers(start);
    opt <- get("optwrap", asNamespace("lme4"))(optimizer,
                          devfun,
                          startingValues,
                          lower = lowerBounds,
                          control=control,
                          adj=FALSE, verbose=verbose)

    if (restart_edge) {
        ## FIXME: should we be looking at rho$pp$theta or opt$par
        ##  at this point???  in koller example (for getData(13)) we have
        ##   rho$pp$theta=0, opt$par=0.08
        if (length(bvals <- which(rho$pp$theta==rho$lower))>0) {
            par <- opt$par;
            ## *don't* use numDeriv -- cruder but fewer dependencies, no worries
            ##  about keeping to the interior of the allowed space
            theta0 <- new("numeric",rho$pp$theta) ## 'deep' copy ...
            d0 <- devfun(par)
            btol <- 1e-5  ## FIXME: make user-settable?
            bgrad <- sapply(bvals,
                            function(i) {
                                bndval <- rho$lower[i]
                                par[1:length(theta0)] <- theta0;
                                par[i] <- bndval+btol
                                (devfun(par)-d0)/btol
                            })
            ## what do I need to do to reset rho$pp$theta to original value???
            par[1:length(theta0)] <- theta0;
            devfun(par) ## reset rho$pp$theta after tests
            ## FIXME: allow user to specify ALWAYS restart if on boundary?
            if (any(bgrad<0)) {
                if (verbose) message("some theta parameters on the boundary, restarting")
                opt <- get("optwrap", asNamespace("lme4"))(optimizer,
                                      devfun,
                                      opt$par,
                                      lower=lowerBounds, control=control,
                                      adj=FALSE, verbose=verbose)
            }
        }
    }
    return(opt)
}
