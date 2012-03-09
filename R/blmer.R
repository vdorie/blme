blmer <-
  function(formula, data, family = NULL, REML = TRUE,
           control = list(), start = NULL, verbose = FALSE, doFit = TRUE,
           subset, weights, na.action, offset, contrasts = NULL,
           model = TRUE, x = TRUE, cov.prior = "wishart",
           fixef.prior = "normal", var.prior = "inverse.gamma",
           ...)
{
  mc <- match.call()
  if (!is.null(family)) {             # call bglmer
    mc[[1]] <- as.name("bglmer")
    return(eval.parent(mc))
  }
  stopifnot(length(formula <- as.formula(formula)) == 3)
  
  fr <- lmerFrames(mc, formula, contrasts) # model frame, X, etc.
  FL <- lmerFactorList(formula, fr, 0L, 0L) # flist, Zt, dims
  largs <- list(...)
  if (!is.null(method <- list(...)$method)) {
    warning(paste("Argument", sQuote("method"),
                  "is deprecated.  Use", sQuote("REML"),
                  "instead"))
    REML <- match.arg(method, c("REML", "ML")) == "REML"
    largs <- largs[names(largs) != "method"]
  }
  if (length(largs) > 0) {
    warning("the following '...' arguments have  *not* been used: ",
            sub("^list", "", deparse(largs, control=NULL)))
  }
  ### FIXME: issue a warning if the control argument has an msVerbose component
  cv <- do.call(lmerControl, control)
  if (missing(verbose)) verbose <- cv$msVerbose
  FL$dims["mxit"] <- cv$maxIter
  FL$dims["mxfn"] <- cv$maxFN
  
  ans <- list(fr = fr, FL = FL, start = start, REML = REML,
              verbose = verbose,
              covariancePrior = cov.prior,
              unmodeledCoefficientPrior = fixef.prior,
              commonScalePrior = var.prior,
              callingEnvironment = parent.frame());
  
  if (doFit) {
    ans <- do.call(blmer_finalize, ans)
    ans@call <- mc
  }
  return(ans);
}

bglmer <-
function(formula, data, family = gaussian, start = NULL,
         verbose = FALSE, nAGQ = 1, doFit = TRUE, subset, weights,
         na.action, offset, contrasts = NULL, model = TRUE,
         control = list(), cov.prior = "wishart",
         fixef.prior = "normal",
         ...)
### Fit a generalized linear mixed model
{
    mc <- match.call()
    ## Evaluate and check the family [[hmm.. have  famType() for that ...]]
    if(is.character(family))
        family <- get(family, mode = "function", envir = parent.frame(2))
    if(is.function(family)) family <- family()
    if(!is.list(family) || is.null(family$family))
      stop(gettextf("family '%s' not recognized", deparse(substitute(family)),
                    domain = "R-lme4"))
    if(family$family == "gaussian" && family$link == "identity") {
        mc[[1]] <- as.name("blmer")      # use blmer not bglmer
        mc$family <- NULL
        return(eval.parent(mc))
    }
    stopifnot(length(formula <- as.formula(formula)) == 3)

    ## Check for method argument which is no longer used
    if (!is.null(method <- list(...)$method)) {
        msg <- paste("Argument", sQuote("method"),
                     "is deprecated.\nUse", sQuote("nAGQ"),
                     "to choose AGQ.  PQL is not available.")
        if (match.arg(method, c("Laplace", "AGQ")) == "Laplace") {
            warning(msg)
        } else stop(msg)
    }

    fr <- lmerFrames(mc, formula, contrasts) # model frame, X, etc.
    offset <- wts <- NULL
    if (length(fr$wts)) wts <- fr$wts
    if (length(fr$off)) offset <- fr$off
    glmFit <- glm.fit(fr$X, fr$Y, weights = wts, # glm on fixed effects
                      offset = offset, family = family,
                      intercept = attr(attr(fr$mf, "terms"), "intercept") > 0)
    FL <- lmerFactorList(formula, fr, 0L, 0L) # flist, Zt, dims
    
### FIXME: issue a warning if the control argument has an msVerbose component
    cv <- do.call(lmerControl, control)
    if (missing(verbose)) verbose <- cv$msVerbose
### FIXME: issue a warning if the model argument is FALSE.  It is ignored.
    FL$dims["mxit"] <- cv$maxIter
    FL$dims["mxfn"] <- cv$maxFN

    ans <- list(fr = fr, FL = FL, glmFit = glmFit, start = start,
                nAGQ = nAGQ, verbose = verbose,
                covariancePrior = cov.prior,
                unmodeledCoefficientPrior = fixef.prior,
                callingEnvironment = parent.frame());
    if (doFit) {
        ans <- do.call(bglmer_finalize, ans)
        ans@call <- mc
    }
    ans
}

blmer_finalize <- function(fr, FL, start, REML, verbose,
                           covariancePrior, unmodeledCoefficientPrior,
                           commonScalePrior, callingEnvironment)
{
  Y <- as.double(fr$Y)
  if (is.list(start) && all(sort(names(start)) == sort(names(FL))))
    start <- list(ST = start)
  if (is.numeric(start)) start <- list(STpars = start)
  dm <- mkZt(FL, start[["ST"]])
  
  dm$dd["REML"] <- as.logical(REML)
  dm$dd["verb"] <- as.integer(verbose)
  swts <- sqrt(unname(fr$wts))
  p <- dm$dd["p"]
  n <- length(Y)
  
  ans <- new(Class = "bmer",
             env = new.env(),
             nlmodel = (~I(x))[[2]],
             frame = fr$mf,
             call = call("foo"),      # later overwritten
             flist = dm$flist,
             X = fr$X,
             Zt = dm$Zt,
             pWt = unname(fr$wts),
             offset = unname(fr$off),
### FIXME: Should y retain its names? As it stands any row names in the
### frame are dropped.  Really?  Are they part of the frame slot (if not
### reduced to 0 rows)?
             y = unname(Y),
             Gp = unname(dm$Gp),
             dims = dm$dd,
             ST = dm$ST,
             A = dm$A,
             Cm = dm$Cm,
             Cx = if (length(swts) > 0) (dm$A)@x else numeric(0),
             L = dm$L,
             deviance = dm$dev,
             fixef = fr$fixef,
             ranef = numeric(dm$dd[["q"]]),
             u = numeric(dm$dd[["q"]]),
             eta = numeric(n),
             mu = numeric(n),
             resid = numeric(n),
             sqrtrWt = swts,
             sqrtXWt = as.matrix(swts),
             RZX = matrix(0, dm$dd[["q"]], p),
             RX = matrix(0, p, p),
             cov.prior = list(),
             fixef.prior = createFlatPriorObject(),
             var.prior = createFlatPriorObject())
  if (!is.null(stp <- start$STpars) && is.numeric(stp)) {
    STp <- .Call(mer_ST_getPars, ans)
    if (length(STp) == length(stp))
      .Call(mer_ST_setPars, ans, stp)
  }

  ans <- setPrior(ans, covariancePrior, unmodeledCoefficientPrior,
                  commonScalePrior, callingEnvironment);

  ### This checks that the number of levels in a grouping factor < n
  ### Only need to check the first factor because it is the one with
  ### the most levels.
  #
  # blme edit:
  # if the common scale is fixed, can work with num groups = n
  numGroups <- length(levels(dm$flist[[1]]));
  numObserv <- length(Y);
  if (numGroups >= numObserv) {
    if (numGroups > numObserv ||
        (ans@var.prior@type == getEnumOrder(typeEnumeration, DIRECT_TYPE_NAME) &&
         ans@var.prior@families[1] != getEnumOrder(familyEnumeration, POINT_FAMILY_NAME)))
      stop(paste("Number of levels of a grouping factor for the random effects",
                 "must be less than the number of observations", sep = "\n"))
  }
  
  return(mer_finalize(ans));
}

bglmer_finalize <- function(fr, FL, glmFit, start, nAGQ, verbose,
                            covariancePrior, unmodeledCoefficientPrior,
                            callingEnvironment)
{
  if (is.list(start) && all(sort(names(start)) == sort(names(FL))))
    start <- list(ST = start)
  if (is.numeric(start)) start <- list(STpars = start)
  dm <- mkZt(FL, start[["ST"]])
  ft <- famType(glmFit$family)
  dm$dd[names(ft)] <- ft
  useSc <- as.integer(!(famNms[dm$dd[["fTyp"]]] %in%
                        c("binomial", "poisson")))
  dm$dd["useSc"] <- useSc
  ## Only need to check the first factor because it is the one with
  ## the most levels.
  M1 <- length(levels(dm$flist[[1]]))
  n <- ncol(dm$Zt)
  if (M1 >= n) {
    msg1 <- "Number of levels of a grouping factor for the random effects\n"
    msg3 <- "n, the number of observations"
    if (useSc)
      stop(msg1, "must be less than ", msg3)
    else if (M1 == n)
      message(msg1, "is *equal* to ", msg3)
  }
  if ((nAGQ <- as.integer(nAGQ)) < 1) nAGQ <- 1L
  if (nAGQ %% 2 == 0) nAGQ <- nAGQ + 1L # reset nAGQ to be an odd number
  dm$dd[["nAGQ"]] <- as.integer(nAGQ)
  AGQlist <- .Call(lme4_ghq, nAGQ)
  y <- unname(as.double(glmFit$y))
  ##    dimnames(fr$X) <- NULL
  p <- dm$dd[["p"]]
  dm$dd["verb"] <- as.integer(verbose)
  fixef <- fr$fixef
  fixef[] <- coef(glmFit)
  if (!is.null(ff <- start$fixef) && is.numeric(ff) &&
      length(ff) == length(fixef)) fixef <- ff

  ans <- new(Class = "bmer",
             env = new.env(),
             nlmodel = (~I(x))[[2]],
             frame = fr$mf,
             call = call("foo"),      # later overwritten
             flist = dm$flist,
             Zt = dm$Zt, X = fr$X, y = y,
             pWt = unname(glmFit$prior.weights),
             offset = unname(fr$off),
             Gp = unname(dm$Gp),
             dims = dm$dd, ST = dm$ST, A = dm$A,
             Cm = dm$Cm, Cx = (dm$A)@x, L = dm$L,
             deviance = dm$dev,
             fixef = fixef,
             ranef = numeric(dm$dd[["q"]]),
             u = numeric(dm$dd[["q"]]),
             eta = unname(glmFit$linear.predictors),
             mu = unname(glmFit$fitted.values),
             muEta = numeric(dm$dd[["n"]]),
             var = numeric(dm$dd[["n"]]),
             resid = unname(glmFit$residuals),
             sqrtXWt = as.matrix(numeric(dm$dd[["n"]])),
             sqrtrWt = numeric(dm$dd[["n"]]),
             RZX = matrix(0, dm$dd[["q"]], p),
             RX = matrix(0, p, p),
             ghx = AGQlist[[1]],
             ghw = AGQlist[[2]],
             cov.prior = list(),
             fixef.prior = createFlatPriorObject(),
             var.prior = createFlatPriorObject());
  if (!is.null(stp <- start$STpars) && is.numeric(stp)) {
    STp <- .Call(mer_ST_getPars, ans)
    if (length(STp) == length(stp))
      .Call(mer_ST_setPars, ans, stp)
  }

  ans <- setPrior(ans,
                  covariancePrior,
                  unmodeledCoefficientPrior,
                  NULL,
                  callingEnvironment);
  
  return (mer_finalize(ans));
}

validateRegressionArgument <- function(regression, regressionName) {
  if (missing(regression)) stop("'regression' missing.");
  
  # check for existence and null-ness
  if (is.null(regression)) stop("object '", regressionName, "' is null.");
  if (!inherits(regression, "bmer")) stop("object '", regressionName, "' does not inherit from S4 class 'bmer'.");
}

setPrior <- function(regression, cov.prior = NULL,
                     fixef.prior = NULL, var.prior = NULL, env = parent.frame())
{
  validateRegressionArgument(regression, match.call()$regression);
  
  if (!is.null(cov.prior)) {
    regression@cov.prior <- parseCovariancePriorSpecification(regression, cov.prior, env);
  }
  if (!is.null(fixef.prior)) {
    regression@fixef.prior <- parseUnmodeledCoefficientPriorSpecification(regression, fixef.prior, env);
    if (is.na(regression@deviance[["sigmaREML"]])) regression@deviance[["sigmaREML"]] <- 1.0;
    if (is.na(regression@deviance[["sigmaML"]])) regression@deviance[["sigmaML"]] <- 1.0;
  }
  if (!is.null(var.prior)) {
    regression@var.prior <- parseCommonScalePriorSpecification(regression, var.prior, env);
  }
  return (regression);
}

parsePrior <- function(regression, cov.prior = NULL,
                       fixef.prior = NULL, var.prior = NULL, env = parent.frame())
{
  validateRegressionArgument(regression, match.call()$regression);
  
  if (!is.null(cov.prior)) {
    return (parseCovariancePriorSpecification(regression, cov.prior, env));
  }
  if (!is.null(fixef.prior)) {
    return (parseUnmodeledCoefficientPriorSpecification(regression, fixef.prior, env));
  }
  if (!is.null(var.prior)) {
    return (parseCommonScalePriorSpecification(regression, var.prior, env));
  }
}

runOptimizer <- function(regression, verbose=FALSE)
{
  validateRegressionArgument(regression, match.call()$regression);
  
  if (verbose) {
    regression@dims[["verb"]] <- as.integer(1)
  } else {
    regression@dims[["verb"]] <- as.integer(0)
  }
  return (mer_finalize(regression));
}

runOptimizerWithPrior <- function(regression, cov.prior = NULL,
                                  fixef.prior = NULL, var.prior = NULL,
                                  verbose = FALSE, env = parent.frame())
{
  validateRegressionArgument(regression, match.call()$regression);
  
  regression <- setPrior(regression, cov.prior, fixef.prior, var.prior, env);
  
  return(runOptimizer(regression, verbose));
}

getDeviance <- function(regression)
{
  validateRegressionArgument(regression, match.call()$regression);
  
  return (.Call(bmer_get_dev, regression));
}

approximateDeviance <- function(regression)
{
  validateRegressionArgument(regression, match.call()$regression);
  
  return (.Call(bmer_approximate_dev, regression));
}

getDerivatives <- function(regression)
{
  validateRegressionArgument(regression, match.call()$regression);
  
  return (.Call(bmer_getDerivatives, regression));
}

getPriorPenalty <- function(regression) {
  validateRegressionArgument(regression, match.call()$regression);

  return(.Call(bmer_calculatePriorPenalty, regression));
}
