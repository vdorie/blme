## for R-check
wishart <- "ignored";
blmer <- function(formula, data = NULL, REML = TRUE,
                  control = lmerControl(), start = NULL,
                  verbose = 0L, subset, weights, na.action, offset,
                  contrasts = NULL, devFunOnly = FALSE,
                  cov.prior = wishart, fixef.prior = NULL,
                  resid.prior = NULL,
                  ...)
{
  mc <- mcout <- match.call()
  missCtrl <- missing(control)
  missCovPrior <- missing(cov.prior)
  ## see functions in modular.R for the body ...
  if (!missCtrl && !inherits(control, "lmerControl")) {
    if(!is.list(control)) stop("'control' is not a list; use lmerControl()")
    ## back-compatibility kluge
    warning("passing control as list is deprecated: please use lmerControl() instead",
            immediate.=TRUE)
    control <- do.call(lmerControl, control)
  }
  if (!is.null(list(...)[["family"]])) {
    warning("calling lmer with 'family' is deprecated; please use glmer() instead")
    mc[[1]] <- quote(lme4::glmer)
    if(missCtrl) mc$control <- glmerControl()
    return(eval(mc, parent.frame(1L)))
  }

  fixef.prior <- mc$fixef.prior; ## for delayed evaluation, get quoted versions
  cov.prior <- if (!missCovPrior) mc$cov.prior else formals(blmer)$cov.prior;
  resid.prior <- mc$resid.prior;
  if (!is.null(mc$var.prior)) resid.prior <- parse(text = mc$var.prior)[[1]];
  mc$fixef.prior <- NULL;
  mc$cov.prior <- NULL;
  mc$resid.prior <- NULL;
  mc$var.prior <- NULL;
  
  sigmaIsFixed <-
    !is.null(resid.prior) && (grepl("^\\W*point", resid.prior) || (is.call(resid.prior) && resid.prior[[1]] == "point"));
  
  if (sigmaIsFixed) {
    control$checkControl$check.nobs.vs.nlev <- "ignore";
    control$checkControl$check.nobs.vs.rankZ <- "ignore";
    control$checkControl$check.nobs.vs.nRE <- "ignore";
  }
  mc$control <- control ## update for  back-compatibility kluge
  
  ## https://github.com/lme4/lme4/issues/50
  ## parse data and formula
  mc[[1]] <- quote(lme4::lFormula)
  lmod <- eval(mc, parent.frame(1L));
  mcout$formula <- lmod$formula
  lmod$formula <- NULL

  ## peel off the starting values lmer stuff expects to see
  lmerStart <- NULL;
  if (!is.null(start) && is.list(start) && length(start) > 1)
    lmerStart <- start$theta;
  
  devfun <- do.call(mkBlmerDevfun,
                    c(lmod, lmod$X, lmod$reTrms,
                      list(priors = list(covPriors = cov.prior, fixefPrior = fixef.prior, residPrior = resid.prior),
                           start = lmerStart, verbose = verbose, control = control)));
  
  if (devFunOnly) return(devfun);
  
  devFunEnv <- environment(devfun);
  ## start <- getStartingValues(start, devFunEnv);

  opt <- optimizeLmer(devfun,
                      optimizer=control$optimizer,
                      restart_edge=control$restart_edge,
                      control=control$optCtrl,
                      verbose=verbose,
                      start=start,
                      boundary.tol=control$boundary.tol,
                      calc.derivs=control$calc.derivs,
                      use.last.params=control$use.last.params)
  cc <- NULL;
  lme4Env <- asNamespace("lme4");
  if (exists("checkConv", lme4Env)) {
    cc <- get("checkConv", lme4Env)(attr(opt,"derivs"),opt$par,
                                    ctrl = control$checkConv,
                                    lbound=environment(devfun)$lower)
  }

  args <- list(rho = devFunEnv, opt = opt, reTrms = lmod$reTrms, fr = lmod$fr, mc = mcout);
  if (!is.null(formals(mkMerMod)$lme4conv)) args$lme4conv <- cc;
  result <- do.call("mkMerMod", args, TRUE, sys.frame(0));
  result <- repackageMerMod(result, opt, devFunEnv);
  
  return(result);
}

bglmer <- function(formula, data=NULL, family = gaussian,
                   control = glmerControl(), start = NULL, verbose = 0L, nAGQ = 1L,
                   subset, weights, na.action, offset,
                   contrasts = NULL, mustart, etastart, devFunOnly = FALSE,
                   cov.prior = wishart, fixef.prior = NULL,
                   ...)
{
  covPriorMissing <- missing(cov.prior);
  
  if (!inherits(control, "glmerControl")) {
    if(!is.list(control)) stop("'control' is not a list; use glmerControl()")
    ## back-compatibility kluge
    msg <- "Use control=glmerControl(..) instead of passing a list"
    if(length(cl <- class(control))) msg <- paste(msg, "of class", dQuote(cl[1]))
    warning(msg, immediate.=TRUE)
    control <- do.call(glmerControl, control)
  }
  mc <- mcout <- match.call()
  
  fixef.prior <- mc$fixef.prior; ## for delayed evaluation, store as quoted
  cov.prior <- if (!covPriorMissing) mc$cov.prior else formals(bglmer)$cov.prior;
  mc$fixef.prior <- NULL;
  mc$cov.prior <- NULL;
  
  ## family-checking code duplicated here and in glFormula (for now) since
  ## we really need to redirect at this point; eventually deprecate formally
  ## and clean up
  if (is.character(family))
    family <- get(family, mode = "function", envir = parent.frame(2))
  if( is.function(family)) family <- family()
  if (isTRUE(all.equal(family, gaussian()))) {
    ## redirect to lmer (with warning)
    warning("calling glmer() with family=gaussian (identity link) as a shortcut to lmer() is deprecated;",
            " please call lmer() directly")
    mc[[1]] <- quote(lme4::lmer)
    mc["family"] <- NULL            # to avoid an infinite loop
    return(eval(mc, parent.frame()))
  }
  
  ## see https://github.com/lme4/lme4/issues/50
  ## parse the formula and data
  mc[[1]] <- quote(lme4::glFormula)
  glmod <- eval(mc, parent.frame(1L))
  mcout$formula <- glmod$formula
  glmod$formula <- NULL
  
  ## create deviance function for covariance parameters (theta)
  
  devfun <- do.call(mkBglmerDevfun, c(glmod, glmod$X, glmod$reTrms,
                                      list(priors = list(covPriors = cov.prior, fixefPrior = fixef.prior),
                                           verbose = verbose, control = control, nAGQ = 0)))
  if (nAGQ==0 && devFunOnly) return(devfun)
  ## optimize deviance function over covariance parameters
  
  if (is.list(start) && !is.null(start$fixef))
    if (nAGQ==0) stop("should not specify both start$fixef and nAGQ==0")

  args <- list(devfun = devfun, optimizer = control$optimizer[[1]],
               restart_edge = if (nAGQ == 0) control$restart_edge else FALSE,
               start = start, verbose = verbose, control = control$optCtrl, nAGQ = 0);
  if (!is.null(formals(optimizeGlmer)$boundary.tol)) args$boundary.tol <- if (nAGQ == 0) control$boundary.tol else 0;
  if (!is.null(formals(optimizeGlmer)[["..."]])) {
    args$calc.derivs <- FALSE;
  }
  opt <- do.call("optimizeGlmer", args, TRUE, sys.frame(0));
  
  if(nAGQ > 0L) {
    
    start <- get("updateStart", asNamespace("lme4"))(start,theta=opt$par)
    
    ## update deviance function to include fixed effects as inputs
    devfun <- updateBglmerDevfun(devfun, glmod$reTrms, nAGQ = nAGQ)
    
    if (devFunOnly) return(devfun)
    
    args$devfun <- devfun;
    args$optimizer <- control$optimizer[[2]];
    args$restart_edge <- control$restart_edge;
    args$nAGQ <- nAGQ;
    args$stage <- 2;
    if (!is.null(formals(optimizeGlmer)$boundary.tol)) args$boundary.tol <- control$boundary.tol;
    if (!is.null(formals(optimizeGlmer)[["..."]])) {
      args$calc.derivs <- FALSE;
      args$use.last.params <- control$use.last.params;
    }
    ## reoptimize deviance function over covariance parameters and fixed effects
    opt <- do.call("optimizeGlmer", args, TRUE, sys.frame(0));
  }
  lme4Env <- asNamespace("lme4");
  cc <- if (!is.null(control$calc.derivs) && !control$calc.derivs) NULL else {
    if (verbose>10) cat("checking convergence\n")
    if (exists("checkConv", lme4Env))
      get("checkConv", lme4Env)(attr(opt,"derivs"), opt$par,
                ctrl = control$checkConv,
                lbound=environment(devfun)$lower)
    else NULL
  }
  
  ## prepare output
  args <- list(rho = environment(devfun), opt = opt, reTrms = glmod$reTrms, fr = glmod$fr, mc = mcout);
  if (!is.null(formals(mkMerMod)$lme4conv)) args$lme4conv <- cc;
  result <- do.call("mkMerMod", args, TRUE, sys.frame(0));
  result <- repackageMerMod(result, opt, environment(devfun));
  
  return(result);
}

lmmObjective <- function(pp, resp, sigma, exponentialTerms, polynomialTerm, blmerControl) {
  sigma.sq <- sigma^2;
  result <- resp$objective(pp$ldL2(), pp$ldRX2(), pp$sqrL(1.0), sigma.sq);

  exponentialTerm <- 0;
  for (i in 1:length(exponentialTerms)) {
    power <- as.numeric(names(exponentialTerms)[[i]]);
    value <- exponentialTerms[[i]];
    if (!is.finite(value)) return(value);
    
    exponentialTerm <- exponentialTerm + value * sigma^power;
  }

  priorPenalty <- exponentialTerm + polynomialTerm + blmerControl$constant + blmerControl$df * log(sigma.sq);
  
  result <- result + priorPenalty;

  return(result);
}

repackageMerMod <- function(merMod, opt, devFunEnv) {
  isLMM <- is(merMod, "lmerMod");

  blmerControl <- devFunEnv$blmerControl;
  priors <- devFunEnv$priors;

  if (isLMM) {
    expandParsInCurrentFrame(opt$par, devFunEnv$parInfo);
    if (blmerControl$fixefOptimizationType != FIXEF_OPTIM_NUMERIC) beta <- merMod@pp$beta(1.0)
    else merMod@beta <- beta;
  } else {
    beta <- opt$par[-devFunEnv$dpars];
  }
  
  Lambda.ts <- getCovBlocks(merMod@pp$Lambdat, devFunEnv$ranefStructure);
  exponentialTerms <- calculatePriorExponentialTerms(priors, beta, Lambda.ts, sigma);

  if (isLMM) {
    if (blmerControl$fixefOptimizationType == FIXEF_OPTIM_NUMERIC) {
      fixefExponentialTerm <- calculateFixefExponentialTerm(beta, merMod@pp$beta(1.0), merMod@pp$RX());
      if (is.null(exponentialTerms[["-2"]])) {
        exponentialTerms[["-2"]] <- fixefExponentialTerm;
      } else {
        exponentialTerms[["-2"]] <- exponentialTerms[["-2"]] + fixefExponentialTerm;
      }
    }

    if (!is.null(exponentialTerms[["-2"]]))
      merMod@devcomp$cmp[["pwrss"]] <- merMod@devcomp$cmp[["pwrss"]] + as.numeric(exponentialTerms[["-2"]]);
  
    ## recover sigma
    sigmaOptimizationType <- blmerControl$sigmaOptimizationType;
    if (sigmaOptimizationType == SIGMA_OPTIM_POINT) {
      sigma <- priors$residPrior@value;
    } else if (sigmaOptimizationType != SIGMA_OPTIM_NUMERIC) {
      profileSigma <- getSigmaProfiler(priors, blmerControl);
      sigma <- profileSigma(merMod@pp, merMod@resp, exponentialTerms, blmerControl);
    }
    ## set sigma in final object
    numObs   <- merMod@devcomp$dims[["n"]];
    numFixef <- merMod@devcomp$dims[["p"]];
    if (merMod@devcomp$dims[["REML"]] > 0L) {
      merMod@devcomp$cmp[["sigmaREML"]] <- sigma;
      merMod@devcomp$cmp[["sigmaML"]] <- sigma * sqrt((numObs - numFixef) / numObs);
    } else {
      merMod@devcomp$cmp[["sigmaML"]] <- sigma;
      merMod@devcomp$cmp[["sigmaREML"]] <- sigma * sqrt(numObs / (numObs - numFixef));
    }

    
    objectiveValue <- merMod@resp$objective(merMod@pp$ldL2(), merMod@pp$ldRX2(), merMod@pp$sqrL(1.0), sigma^2);
    if (blmerControl$fixefOptimizationType == FIXEF_OPTIM_NUMERIC)
      objectiveValue <- objectiveValue + fixefExponentialTerm / sigma^2;
    
    if (merMod@devcomp$dims[["REML"]] > 0L) {
      priorPenalty <- merMod@devcomp$cmp[["REML"]] - objectiveValue;
      merMod@devcomp$cmp[["REML"]] <- objectiveValue;
    } else {
      priorPenalty <- merMod@devcomp$cmp[["dev"]] - objectiveValue;
      merMod@devcomp$cmp[["dev"]] <- objectiveValue;
    }
    merMod@devcomp$cmp[["penalty"]] <- priorPenalty;

    return(new("blmerMod",
               resp    = merMod@resp,
               Gp      = merMod@Gp,
               call    = merMod@call,
               frame   = merMod@frame,
               flist   = merMod@flist,
               cnms    = merMod@cnms,
               lower   = merMod@lower,
               theta   = merMod@theta,
               beta    = beta,
               u       = merMod@u,
               devcomp = merMod@devcomp,
               pp      = merMod@pp,
               optinfo = merMod@optinfo,
               priors  = priors));
  } else {
    if (length(exponentialTerms) > 0)
      priorPenalty <- exponentialTerms[[1]] + calculatePriorPolynomialTerm(priors$covPriors, Lambda.ts) + blmerControl$constant
    else
      priorPenalty <- 0;
    merMod@devcomp$cmp[["dev"]] <- merMod@devcomp$cmp[["dev"]] - priorPenalty;
    merMod@devcomp$cmp[["penalty"]] <- priorPenalty;

    return(new("bglmerMod",
               resp    = merMod@resp,
               Gp      = merMod@Gp,
               call    = merMod@call,
               frame   = merMod@frame,
               flist   = merMod@flist,
               cnms    = merMod@cnms,
               lower   = merMod@lower,
               theta   = merMod@theta,
               beta    = merMod@beta,
               u       = merMod@u,
               devcomp = merMod@devcomp,
               pp      = merMod@pp,
               optinfo = merMod@optinfo,
               priors  = priors));
  }
}

validateRegressionArgument <- function(regression, regressionName) {
  if (missing(regression)) stop("'regression' missing.");
  
  # check for existence and null-ness
  if (is.null(regression)) stop("object '", regressionName, "' is null.");
  if (!is(regression, "bmerMod")) stop("object '", regressionName, "' does not inherit from S4 class 'bmerMod'.");
}

setPrior <- function(regression, cov.prior = NULL,
                     fixef.prior = NULL, resid.prior = NULL, envir = parent.frame(1L), ...)
{
  matchedCall <- match.call();

  covMissing   <- missing(cov.prior);
  fixefMissing <- missing(fixef.prior);
  residMissing <- missing(resid.prior);
  
  validateRegressionArgument(regression, matchedCall$regression);
  
  if (residMissing && !is.null(matchedCall$var.prior)) {
    matchedCall$resid.prior <- matchedCall$var.prior;
    residMissing <- FALSE;
  }

  priors <- evaluatePriorArguments(matchedCall$cov.prior, matchedCall$fixef.prior, matchedCall$resid.prior,
                                   regression@devcomp$dim, regression@cnms,
                                   as.integer(diff(regression@Gp) / sapply(regression@cnms, length)),
                                   envir);

  if (!covMissing) regression@covPriors <- priors$covPriors;
  if (!fixefMissing) regression@fixefPrior <- priors$fixefPrior;
  if (!residMissing) regression@residPrior <- priors$residPrior;
  
  return (regression);
}

parsePrior <- function(regression, cov.prior = NULL,
                       fixef.prior = NULL, resid.prior = NULL, envir = parent.frame(), ...)
{
  matchedCall <- match.call();

  covMissing   <- missing(cov.prior);
  fixefMissing <- missing(fixef.prior);
  residMissing <- missing(resid.prior);
  
  validateRegressionArgument(regression, matchedCall$regression);
  
  if (residMissing && !is.null(matchedCall$var.prior)) {
    matchedCall$resid.prior <- matchedCall$var.prior;
    residMissing <- FALSE;
  }

  priors <- evaluatePriorArguments(matchedCall$cov.prior, matchedCall$fixef.prior, matchedCall$resid.prior,
                                   regression@devcomp$dim, regression@cnms,
                                   as.integer(diff(regression@Gp) / sapply(regression@cnms, length)),
                                   envir);

  result <- list();
  if (!covMissing) result$covPriors <- priors$covPriors;
  if (!fixefMissing) result$fixefPrior <- priors$fixefPrior;
  if (!residMissing) result$residPrior <- priors$residPrior;

  if (length(result) == 1) return(result[[1]]);
  return(result);
}

if (FALSE) {
runOptimizer <- function(regression, verbose = FALSE)
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
                                  verbose = FALSE, envir = parent.frame())
{
  validateRegressionArgument(regression, match.call()$regression);
  
  regression <- setPrior(regression, cov.prior, fixef.prior, var.prior, envir);
  
  return(runOptimizer(regression, verbose));
}
}
