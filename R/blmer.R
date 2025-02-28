## copyright note:
##   a lot of of this was copy/pasted from the lme4 package (http://cran.r-project.org/web/packages/lme4/index.html, GPL-2)
##   a lot of it was not
## ideally, blme wouldn't have to recreate the functions with minor tweaks but we're just not there yet

## for R-check
wishart <- "ignored"

blmer <- function(formula, data = NULL, REML = TRUE,
                  control = lmerControl(), start = NULL,
                  verbose = 0L, subset, weights, na.action, offset,
                  contrasts = NULL, devFunOnly = FALSE,
                  cov.prior = wishart, fixef.prior = NULL,
                  resid.prior = NULL,
                  ...)
{
  mc <- mcout <- match.call()
  callingEnv <- parent.frame(1L)
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
  
  fixef.prior <- mc$fixef.prior ## for delayed evaluation, get quoted versions
  cov.prior <- if (!missCovPrior) mc$cov.prior else formals(blmer)$cov.prior
  resid.prior <- mc$resid.prior
  if (!is.null(mc$var.prior)) resid.prior <- parse(text = mc$var.prior)[[1]]
  mc$fixef.prior <- NULL
  mc$cov.prior <- NULL
  mc$resid.prior <- NULL
  mc$var.prior <- NULL
  
  sigmaIsFixed <-
    !is.null(resid.prior) && ((is.character(resid.prior) && grepl("^\\W*point", resid.prior)) ||
                              (is.call(resid.prior) && resid.prior[[1]] == "point"))
  
  if (sigmaIsFixed) {
    control$checkControl$check.nobs.vs.nlev  <- "ignore"
    control$checkControl$check.nobs.vs.rankZ <- "ignore"
    control$checkControl$check.nobs.vs.nRE   <- "ignore"
  }

  hasPseudoData <-
    !is.null(fixef.prior) && ((is.character(fixef.prior) && grepl("^\\W*normal", fixef.prior)) ||
                              (is.call(fixef.prior) && fixef.prior[[1]] == "normal"))

  if (hasPseudoData) {
    control$checkControl$check.rankX <- "ignore"
  }
  
  mc$control <- control ## update for  back-compatibility kluge
  
  ## https://github.com/lme4/lme4/issues/50
  ## parse data and formula
  mc[[1]] <- quote(lme4::lFormula)
  lmod <- eval(mc, parent.frame(1L))
  mcout$formula <- lmod$formula
  lmod$formula <- NULL

  ## peel off the starting values lmer stuff expects to see
  lmerStart <- NULL
  if (!is.null(start) && is.list(start) && length(start) > 1)
    lmerStart <- start$theta
  devfun <- do.call(mkBlmerDevfun,
                    c(lmod, lmod$X, lmod$reTrms,
                      list(priors = list(covPriors = cov.prior, fixefPrior = fixef.prior, residPrior = resid.prior),
                           start = lmerStart, verbose = verbose, control = control, env = callingEnv)))
  
  if (devFunOnly) return(devfun)
  
  devFunEnv <- environment(devfun)
  opt <- if (control$optimizer=="none")
    list(par=NA,fval=NA,conv=1000,message="no optimization")
  else {
        optimizeLmer(devfun, optimizer = control$optimizer,
                     restart_edge = control$restart_edge,
                     boundary.tol = control$boundary.tol,
                     control = control$optCtrl,
                     verbose=verbose,
                     start=start,
                     calc.derivs=control$calc.derivs,
                     use.last.params=control$use.last.params)
  }
  
  ## dirty hacks to give some backwards lme4 compatibility
  cc <- NULL
  lme4Namespace <- getNamespace("lme4")
  if (exists("checkConv", lme4Namespace)) {
    checkConv <- get("checkConv", lme4Namespace)
    cc <- checkConv(attr(opt, "derivs"), opt$par,
                    ctrl = control$checkConv,
                    lbound = environment(devfun)$lower)
  }

  args <- list(rho = devFunEnv, opt = opt, reTrms = lmod$reTrms, fr = lmod$fr, mc = mcout)
  if ("lme4conv" %in% names(formals(mkMerMod))) args$lme4conv <- cc
  result <- do.call(mkMerMod, args, TRUE)
  result <- repackageMerMod(result, opt, devFunEnv)
  
  return(result)
}

bglmer <- function(formula, data = NULL, family = gaussian,
                   control = glmerControl(), start = NULL, verbose = 0L, nAGQ = 1L,
                   subset, weights, na.action, offset,
                   contrasts = NULL, mustart, etastart, devFunOnly = FALSE,
                   cov.prior = wishart, fixef.prior = NULL,
                   ...)
{
  covPriorMissing <- missing(cov.prior)
  callingEnv <- parent.frame(1L)
  
  if (!inherits(control, "glmerControl")) {
    if(!is.list(control)) stop("'control' is not a list; use glmerControl()")
    ## back-compatibility kluge
    msg <- "Use control=glmerControl(..) instead of passing a list"
    if(length(cl <- class(control))) msg <- paste(msg, "of class", dQuote(cl[1]))
    warning(msg, immediate.=TRUE)
    control <- do.call(glmerControl, control)
  }
  mc <- mcout <- match.call()
  
  fixef.prior <- mc$fixef.prior ## for delayed evaluation, store as quoted
  cov.prior <- if (!covPriorMissing) mc$cov.prior else formals(bglmer)$cov.prior
  mc$fixef.prior <- NULL
  mc$cov.prior <- NULL
  
  ## family-checking code duplicated here and in glFormula (for now) since
  ## we really need to redirect at this point; eventually deprecate formally
  ## and clean up
 if (is.character(family))
   family <- get(family, mode = "function", envir = parent.frame(2))
  if( is.function(family)) family <- family()
  if (isTRUE(all.equal(family, gaussian()))) {
    ## redirect to lmer (with warning)
    warning("calling bglmer() with family=gaussian (identity link) as a shortcut to blmer() is deprecated;",
            " please call blmer() directly")
    mc[[1]] <- quote(blme::blmer)
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
  nAGQinit <- if(control$nAGQ0initStep) 0L else 1L
  devfun <- do.call(mkBglmerDevfun, c(glmod,
                                      list(priors = list(covPriors = cov.prior, fixefPrior = fixef.prior),
                                           verbose = verbose,
                                           control = control,
                                           nAGQ = nAGQinit,
                                           env = callingEnv)))
  if (nAGQ==0 && devFunOnly) return(devfun)
  ## optimize deviance function over covariance parameters
  
  ## FIXME: perhaps should be in glFormula instead??
  if (is.list(start)) {
    start.bad <- setdiff(names(start),c("theta","fixef"))
    if (length(start.bad)>0) {
      stop(sprintf("bad name(s) for start vector (%s); should be %s and/or %s",
                   paste(start.bad,collapse=", "),
                   shQuote("theta"),
                   shQuote("fixef")),call.=FALSE)
    }
    if (!is.null(start$fixef) && nAGQ==0)
      stop("should not specify both start$fixef and nAGQ==0")
  }
  
  if (packageVersion("lme4") <= "1.1-7" || identical(control$nAGQ0initStep, TRUE)) {
    args <- list(devfun = devfun,
                 optimizer = control$optimizer[[1]],
                 ## DON'T try fancy edge tricks unless nAGQ=0 explicitly set
                 restart_edge = if (nAGQ == 0) control$restart_edge else FALSE,
                 control = control$optCtrl,
                 start = start,
                 nAGQ = 0,
                 verbose = verbose)
    if (!is.null(formals(optimizeGlmer)$boundary.tol)) args$boundary.tol <- if (nAGQ == 0) control$boundary.tol else 0
    if (!is.null(formals(optimizeGlmer)[["..."]])) args$calc.derivs <- FALSE
    
    opt <- do.call(optimizeGlmer, args, TRUE)
  }
  
  if(nAGQ > 0L) {
    start <- get("updateStart", getNamespace("lme4"))(start,theta=opt$par)
    
    ## update deviance function to include fixed effects as inputs
    devfun <- updateBglmerDevfun(devfun, glmod$reTrms, nAGQ = nAGQ)
    
    if (devFunOnly) return(devfun)
    ## reoptimize deviance function over covariance parameters and fixed effects 
    args <- list(devfun = devfun,
                 optimizer = control$optimizer[[2]],
                 restart_edge = control$restart_edge,
		 control = control$optCtrl,
                 start = start,
                 nAGQ = nAGQ,
                 verbose = verbose,
                 stage = 2)
    if (!is.null(formals(optimizeGlmer)$boundary.tol)) args$boundary.tol <- control$boundary.tol
    if (!is.null(formals(optimizeGlmer)[["..."]])) {
      args$calc.derivs <- control$calc.derivs
      args$use.last.params <- control$use.last.params
    }
    ## reoptimize deviance function over covariance parameters and fixed effects
    opt <- do.call(optimizeGlmer, args, TRUE)
  }

  lme4Namespace <- getNamespace("lme4")
  cc <- if (!is.null(control$calc.derivs) && !control$calc.derivs) NULL else {
    if (exists("checkConv", lme4Namespace)) {
      if (verbose > 10) cat("checking convergence\n")
      checkConv <- get("checkConv", lme4Namespace)
      checkConv(attr(opt,"derivs"), opt$par,
                ctrl = control$checkConv,
                lbound = environment(devfun)$lower)
    }
    else NULL
  }
  
  ## prepare output
  args <- list(rho = environment(devfun), opt = opt, reTrms = glmod$reTrms, fr = glmod$fr, mc = mcout)
  if ("lme4conv" %in% names(formals(mkMerMod))) args$lme4conv <- cc
  result <- do.call(mkMerMod, args, TRUE)
  result <- repackageMerMod(result, opt, environment(devfun))
  
  return(result)
}

lmmObjective <- function(pp, resp, sigma, exponentialTerms, polynomialTerm, blmerControl) {
  sigma.sq <- sigma^2
  result <- resp$objective(pp$ldL2(), pp$ldRX2(), pp$sqrL(1.0), sigma.sq)

  exponentialTerm <- 0
  for (i in seq_along(exponentialTerms)) {
    power <- as.numeric(names(exponentialTerms)[[i]])
    value <- exponentialTerms[[i]]
    if (!is.finite(value)) return(value)
    
    exponentialTerm <- exponentialTerm + value * sigma^power
  }

  priorPenalty <- exponentialTerm + polynomialTerm + blmerControl$constant + blmerControl$df * log(sigma.sq)
  
  result <- result + priorPenalty

  return(result)
}

repackageMerMod <- function(merMod, opt, devFunEnv) {
  isLMM <- is(merMod, "lmerMod")

  blmerControl <- devFunEnv$blmerControl
  priors <- devFunEnv$priors
  
  sigma <- NULL
  if (isLMM) {
    expandParsInCurrentFrame(opt$par, devFunEnv$parInfo)
    if (blmerControl$fixefOptimizationType != FIXEF_OPTIM_NUMERIC) beta <- merMod@pp$beta(1.0)
    else merMod@beta <- beta
    
    if (blmerControl$sigmaOptimizationType == SIGMA_OPTIM_POINT) sigma <- priors$residPrior@value
  } else {
    beta <- opt$par[-devFunEnv$dpars]
  }
  
  if (!is.null(merMod@optinfo)) {
    parLength <- devFunEnv$parInfo$theta$length + if (!isLMM) devFunEnv$parInfo$beta$length else 0
    if (parLength != length(merMod@optinfo$val)) {
      merMod@optinfo$val_full    <- merMod@optinfo$val
      merMod@optinfo$derivs_full <- merMod@optinfo$derivs
      
      merMod@optinfo$val <- merMod@optinfo$val[parLength]
      merMod@optinfo$derivs$gradient <- merMod@optinfo$derivs$gradient[parLength]
      merMod@optinfo$derivs$Hessian <- merMod@optinfo$derivs$Hessian[parLength, parLength, drop = FALSE]
    }
  }
  
  Lambda.ts <- getCovBlocks(merMod@pp$Lambdat, devFunEnv$ranefStructure)
  exponentialTerms <- calculatePriorExponentialTerms(priors, beta, Lambda.ts, sigma)

  if (isLMM) {
    if (blmerControl$fixefOptimizationType == FIXEF_OPTIM_NUMERIC) {
      ## when beta can't be profiled, a term is added that is a weighted distance of the numeric optimum
      ## and a least-squares solution (what the merMod calculates)
      fixefExponentialTerm <- calculateFixefExponentialTerm(beta, merMod@pp$beta(1.0), merMod@pp$RX())
      if (is.null(exponentialTerms[["-2"]])) {
        exponentialTerms[["-2"]] <- fixefExponentialTerm
      } else {
        exponentialTerms[["-2"]] <- exponentialTerms[["-2"]] + fixefExponentialTerm
      }
    }

    if (!is.null(exponentialTerms[["-2"]]))
      merMod@devcomp$cmp[["pwrss"]] <- merMod@devcomp$cmp[["pwrss"]] + as.numeric(exponentialTerms[["-2"]])
  
    ## recover sigma
    sigmaOptimizationType <- blmerControl$sigmaOptimizationType
    if (sigmaOptimizationType %not_in% c(SIGMA_OPTIM_NA, SIGMA_OPTIM_POINT, SIGMA_OPTIM_NUMERIC)) {
      profileSigma <- getSigmaProfiler(priors, blmerControl)
      sigma <- profileSigma(merMod@pp, merMod@resp, exponentialTerms, blmerControl)
    }
    ## set sigma in final object
    numObs   <- merMod@devcomp$dims[["n"]]
    numFixef <- merMod@devcomp$dims[["p"]]
    if (merMod@devcomp$dims[["REML"]] > 0L) {
      merMod@devcomp$cmp[["sigmaREML"]] <- sigma
      merMod@devcomp$cmp[["sigmaML"]] <- sigma * sqrt((numObs - numFixef) / numObs)
    } else {
      merMod@devcomp$cmp[["sigmaML"]] <- sigma
      merMod@devcomp$cmp[["sigmaREML"]] <- sigma * sqrt(numObs / (numObs - numFixef))
    }

    
    objectiveValue <- merMod@resp$objective(merMod@pp$ldL2(), merMod@pp$ldRX2(), merMod@pp$sqrL(1.0), sigma^2)
    if (blmerControl$fixefOptimizationType == FIXEF_OPTIM_NUMERIC)
      objectiveValue <- objectiveValue + fixefExponentialTerm / sigma^2
    
    if (merMod@devcomp$dims[["REML"]] > 0L) {
      priorPenalty <- merMod@devcomp$cmp[["REML"]] - objectiveValue
      merMod@devcomp$cmp[["REML"]] <- objectiveValue
    } else {
      priorPenalty <- merMod@devcomp$cmp[["dev"]] - objectiveValue
      merMod@devcomp$cmp[["dev"]] <- objectiveValue
    }
    merMod@devcomp$cmp[["penalty"]] <- priorPenalty

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
               priors  = priors))
  } else {
    if (length(exponentialTerms) > 0)
      priorPenalty <- exponentialTerms[[1]] + calculatePriorPolynomialTerm(priors$covPriors, Lambda.ts) + blmerControl$constant
    else
      priorPenalty <- 0
    merMod@devcomp$cmp[["dev"]] <- merMod@devcomp$cmp[["dev"]] - priorPenalty
    merMod@devcomp$cmp[["penalty"]] <- priorPenalty

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
               priors  = priors))
  }
}

validateRegressionArgument <- function(regression, regressionName) {
  if (missing(regression)) stop("'regression' missing.")
  
  # check for existence and null-ness
  if (is.null(regression)) stop("object '", regressionName, "' is null.")
  if (!is(regression, "bmerMod")) stop("object '", regressionName, "' does not inherit from S4 class 'bmerMod'.")
}

setPrior <- function(regression, cov.prior = NULL,
                     fixef.prior = NULL, resid.prior = NULL, envir = parent.frame(1L), ...)
{
  matchedCall <- match.call()

  covMissing   <- missing(cov.prior)
  fixefMissing <- missing(fixef.prior)
  residMissing <- missing(resid.prior)
  
  validateRegressionArgument(regression, matchedCall$regression)
  
  if (residMissing && !is.null(matchedCall$var.prior)) {
    matchedCall$resid.prior <- matchedCall$var.prior
    residMissing <- FALSE
  }
  
  priors <- evaluatePriorArguments(matchedCall$cov.prior, matchedCall$fixef.prior, matchedCall$resid.prior,
                                   regression@devcomp$dim, colnames(regression@pp$X), regression@cnms,
                                   as.integer(diff(regression@Gp) / sapply(regression@cnms, length)),
                                   envir)

  if (!covMissing) regression@covPriors <- priors$covPriors
  if (!fixefMissing) regression@fixefPrior <- priors$fixefPrior
  if (!residMissing) regression@residPrior <- priors$residPrior
  
  return (regression)
}

parsePrior <- function(regression, cov.prior = NULL,
                       fixef.prior = NULL, resid.prior = NULL, envir = parent.frame(), ...)
{
  matchedCall <- match.call()

  covMissing   <- missing(cov.prior)
  fixefMissing <- missing(fixef.prior)
  residMissing <- missing(resid.prior)
  
  validateRegressionArgument(regression, matchedCall$regression)
  
  if (residMissing && !is.null(matchedCall$var.prior)) {
    matchedCall$resid.prior <- matchedCall$var.prior
    residMissing <- FALSE
  }

  priors <- evaluatePriorArguments(matchedCall$cov.prior, matchedCall$fixef.prior, matchedCall$resid.prior,
                                   regression@devcomp$dim, colnames(regression@pp$X), regression@cnms,
                                   as.integer(diff(regression@Gp) / sapply(regression@cnms, length)),
                                   envir)

  result <- list()
  if (!covMissing) result$covPriors <- priors$covPriors
  if (!fixefMissing) result$fixefPrior <- priors$fixefPrior
  if (!residMissing) result$residPrior <- priors$residPrior

  if (length(result) == 1) return(result[[1]])
  return(result)
}

if (FALSE) {
runOptimizer <- function(regression, verbose = FALSE)
{
  validateRegressionArgument(regression, match.call()$regression)
  
  if (verbose) {
    regression@dims[["verb"]] <- as.integer(1)
  } else {
    regression@dims[["verb"]] <- as.integer(0)
  }
  return (mer_finalize(regression))
}

runOptimizerWithPrior <- function(regression, cov.prior = NULL,
                                  fixef.prior = NULL, var.prior = NULL,
                                  verbose = FALSE, envir = parent.frame())
{
  validateRegressionArgument(regression, match.call()$regression)
  
  regression <- setPrior(regression, cov.prior, fixef.prior, var.prior, envir)
  
  return(runOptimizer(regression, verbose))
}
}

getRefitControl <- function(x, controlArg) {
  if (!is.null(controlArg)) {
    if (length(controlArg$optCtrl) == 0) { ## use object's version:
      obj.control <- x@optinfo$control
      ignore.pars <- c("xst", "xt")
      if (any(ign <- names(controlArg) %in% ignore.pars))
        obj.control <- obj.control[!ign]
      controlArg$optCtrl <- obj.control
    }
    controlArg
  } else if (isGLMM(x)) {
    glmerControl()
  } else {
    lmerControl()
  }
}

refit.bmerMod <- function(
  object,
  newresp = NULL,
  newweights = NULL,
  rename.response = FALSE,
  maxit = 100L,
  ...
)
{
  lme4Namespace <- getNamespace("lme4")
  lme4Version   <- packageVersion("lme4")

  l... <- list(...)
  
  ctrl.arg <- NULL
  if ("control" %in% names(l...)) ctrl.arg <- l...$control
  
  if (!all(names(l...) %in% c("control", "verbose")))
    warning("additional arguments to refit.bmerMod ignored")
  
  ## TODO: not clear whether we should reset the names
  ##       to the new response variable.  Maybe not.
  
  ## retrieve name before it gets mangled by operations on newresp
  newrespSub <- substitute(newresp)
  
  ## for backward compatibility/functioning of refit(fit,simulate(fit))
  if (is.list(newresp)) {
    if (length(newresp) == 1) {
      na.action <- attr(newresp,"na.action")
      newresp <- newresp[[1]]
      attr(newresp, "na.action") <- na.action
    } else {
      stop("refit not implemented for 'newresp' lists with length > 1: ",
           "consider ", sQuote("lapply(object, refit)"))
    }
  }
    
  ## oldresp <- object@resp$y # need to set this before deep copy,
  ##                          # otherwise it gets reset with the call
  ##                          # to setResp below
  
  
  ## somewhat repeated from profile.merMod, but sufficiently
  ##  different that refactoring is slightly non-trivial
  ## "three minutes' thought would suffice ..."
  control <- getRefitControl(object, ctrl.arg) # NOTE: blme change

  if (object@optinfo$optimizer == "optimx") {
    control$optCtrl <- object@optinfo$control
  }
  ## we need this stuff defined before we call .glmerLaplace below ...
  pp        <- object@pp$copy()
  dc        <- object@devcomp
  nAGQ      <- unname(dc$dims["nAGQ"]) # possibly NA # NOTE: blme change
  nth       <- dc$dims[["nth"]]
  verbose <- l...$verbose; if (is.null(verbose)) verbose <- 0L
  if (!is.null(newresp)) {
    ## update call and model frame with new response
    rcol <- attr(attr(model.frame(object), "terms"), "response")
    if (rename.response) {
      attr(object@frame,"formula")[[2L]] <- object@call$formula[[2L]] <- newrespSub
      names(object@frame)[rcol] <- deparse(newrespSub)
    }
    if (!is.null(na.act <- attr(object@frame,"na.action")) &&
        is.null(attr(newresp, "na.action"))) {
      ## will only get here if na.action is 'na.omit' or 'na.exclude'
      ## *and* newresp does not have an 'na.action' attribute
      ## indicating that NAs have already been filtered
      newresp <- if (is.matrix(newresp))
        newresp[-na.act, ]
      else newresp[-na.act]
    }
    object@frame[,rcol] <- newresp
  }

  if (!is.null(newweights)) {
    ## DRY ...
    if (!is.null(na.act <- attr(object@frame,"na.action")) &&
      is.null(attr(newweights, "na.action"))) {
      newweights <- newweights[-na.act]
    }
    object@frame[["(weights)"]] <- newweights
    oc <- attr(attr(object@frame, "terms"), "dataClasses")
    attr(attr(object@frame, "terms"), "dataClasses") <- c(oc, `(weights)` = "numeric")
    
    object@call$weights <- substitute(newweights)

    ## try to make sure new weights are findable later
    assign(deparse(substitute(newweights)),
           newweights,
           environment(formula(object)))
  }


  rr <- if (isLMM(object))
    mkRespMod(model.frame(object), REML = isREML(object))
  else if (isGLMM(object)) {
    modelFrame <- model.frame(object) ## NOTE: blme change
    if (lme4Version <= "1.1-6") modelFrame$mustart <- object@resp$mu
    mkRespMod(modelFrame, family = family(object))
  } else
    stop("refit.bmerMod not working for nonlinear mixed models.")
  
  if (!is.null(newresp)) {
    if (family(object)$family == "binomial") {
      ## re-do conversion of two-column matrix and factor
      ##  responses to proportion/weights format
      if (is.matrix(newresp) && ncol(newresp) == 2) {
        ntot <- rowSums(newresp)
        ## FIXME: test what happens for (0,0) rows
        newresp <- newresp[,1] / ntot
        rr$setWeights(ntot)
      }
      if (is.factor(newresp)) {
        ## FIXME: would be better to do this consistently with
        ## whatever machinery is used in glm/glm.fit/glmer ... ??
        newresp <- as.numeric(newresp) - 1
      }
    }
    
    ## if (isGLMM(object) && rr$family$family=="binomial") {
    ## }
    stopifnot(length(newresp <- as.numeric(as.vector(newresp))) ==
              length(rr$y))
    
  }
  
  glmerPwrssUpdate <- get("glmerPwrssUpdate", lme4Namespace)
  if (isGLMM(object)) {
    GQmat <- GHrule(nAGQ)
    if (nAGQ <= 1) {
      if (lme4Version <= "1.1-7")
        glmerPwrssUpdate(pp, rr, control$tolPwrss, GQmat)
      else
        glmerPwrssUpdate(pp, rr, control$tolPwrss, GQmat, maxit = maxit)
    } else {
      if (lme4Version <= "1.1-7")
        glmerPwrssUpdate(pp, rr, control$tolPwrss, GQmat, grpFac = object@flist[[1]])
      else
        glmerPwrssUpdate(pp, rr, control$tolPwrss, GQmat, maxit = maxit, grpFac = object@flist[[1]])
    }
  }
  ## .Call(glmerLaplace, pp$ptr(), rr$ptr(), nAGQ,
  ## control$tolPwrss, as.integer(30), verbose)
  ##              nAGQ,
  ##              control$tolPwrss, as.integer(30), # maxit = 30
  ##              verbose)
  ##        lp0         <- pp$linPred(1) # each pwrss opt begins at this eta

  devlist <- if (isGLMM(object)) {
    baseOffset <- get("forceCopy", lme4Namespace)(object@resp$offset)

    list(tolPwrss    = dc$cmp [["tolPwrss"]],
	 compDev     = dc$dims[["compDev"]],
	 nAGQ        = unname(nAGQ),
	 lp0         = pp$linPred(1), ## object@resp$eta - baseOffset,
	 baseOffset  = baseOffset,
	 pwrssUpdate = glmerPwrssUpdate,
	 ## save GQmat in the object and use that instead of nAGQ
	 GQmat       = GHrule(nAGQ),
	 fac         = object@flist[[1]],
	 pp          = pp,
         resp        = rr,
         u0          = pp$u0,
         verbose     = verbose,
         dpars       = seq_len(nth))
  } else
    list(pp      = pp,
         resp    = rr,
         u0      = pp$u0,
         verbose = verbose,
         dpars   = seq_len(nth))
  
  ## NOTE: blme changes start
  ff <- makeRefitDevFun(list2env(devlist), nAGQ = nAGQ, verbose = verbose, maxit = maxit, object = object)
  reTrms <- list(flist = object@flist, cnms = object@cnms, Gp = object@Gp, lower = object@lower)
  if (isGLMM(object))
    ff <- updateBglmerDevfun(ff, reTrms, nAGQ)
  
  
  ## commenting out xst (not used) and x0, which we grab elsewhere
  ## xst       <- rep.int(0.1, nth)
  ## x0        <- pp$theta
  ## lower     <- object@lower
  lower <- environment(ff)$lower
  ## if (!is.na(nAGQ) && nAGQ > 0L) {
  ##    xst   <- c(xst, sqrt(diag(pp$unsc())))
  ##    x0    <- c(x0, unname(fixef(object)))
  ##     lower <- c(lower, rep(-Inf, length(x0) - length(lower)))
  ##}
  
  ## NOTE: blme changes end
  
  ## control <- c(control, list(xst = 0.2 * xst, xt = xst * 0.0001))
  ## FIX ME: allow use.last.params to be passed through
  calc.derivs <- !is.null(object@optinfo$derivs)
  ## if (isGLMM(object)) {
  ##   rho$resp$updateWts()
  ##   rho$pp$updateDecomp()
  ##   rho$lp0 <- rho$pp$linPred(1)
  ## }

  optimizer <- object@optinfo$optimizer
  if (!is.null(newopt <- ctrl.arg$optimizer)) {
    ## we might end up with a length-2 optimizer vector ...
    ##  use the *last* element
    optimizer <- newopt[length(newopt)]
  }
  
  
  ## NOTE: blme changes start
  opt <-
    if (isLMM(object)) {
      optimizeLmer(ff,
                   optimizer = object@optinfo$optimizer,
                   control = control$optCtrl,
                   verbose = verbose,
                   start = extractParameterListFromFit(object, environment(ff)$blmerControl),
                   calc.derivs = calc.derivs,
                   use.last.params = if (!is.null(control$use.last.params)) control$use.last.params else FALSE)
    } else {
      args <- list(devfun = ff,
                   optimizer = object@optinfo$optimizer,
                   restart_edge = control$restart_edge,
                   control = control$optCtrl,
                   start = extractParameterListFromFit(object, environment(ff)$blmerControl),
                   nAGQ = nAGQ,
                   verbose = verbose,
                   stage = 2)
      if (!is.null(formals(optimizeGlmer)$boundary.tol)) args$boundary.tol <- control$boundary.tol
      if (!is.null(formals(optimizeGlmer)[["..."]])) {
        args$calc.derivs <- control$calc.derivs
        args$use.last.params <- if (!is.null(control$use.last.params)) control$use.last.params else FALSE
      }
      do.call(optimizeGlmer, args, TRUE)
    }
  # NOTE: blmer changes meaningfully end
  cc <- NULL
  if (exists("checkConv", lme4Namespace)) {
    cc <- get("checkConv", lme4Namespace)(attr(opt,"derivs"), opt$par,
                                          ctrl = control$checkConv,
                                          lbound = lower)
  }
  
  if (isGLMM(object)) rr$setOffset(baseOffset)
  
  args <- list(rho = environment(ff), opt = opt,
               reTrms = reTrms,
               fr = object@frame, mc = getCall(object))
  if ("lme4conv" %in% names(formals(mkMerMod))) args$lme4conv <- cc
  result <- do.call(mkMerMod, args, TRUE, sys.frame(0))
  repackageMerMod(result, opt, environment(ff)) 
}

refitML.bmerMod <- function (x, optimizer="bobyqa", ...) {
  # NOTE: blme changes start
  l... <- list(...)

  if (!all(names(l...) %in% c("control", "verbose", "maxit")))
    warning("additional arguments to refitML.bmerMod ignored")

  ctrl.arg <- NULL
  if ("control" %in% names(l...)) ctrl.arg <- l...$control
  control <- getRefitControl(x, ctrl.arg)
  verbose <- l...$verbose; if (is.null(verbose)) verbose <- 0L
  maxit <- l...$maxit; if (is.null(maxit)) maxit <- 100L

  lme4Namespace <- getNamespace("lme4")
  # NOTE: blme changes end


  ## FIXME: optimizer is set to 'bobyqa' for back-compatibility, but that's not
  ##  consistent with lmer (default NM).  Should be based on internally stored 'optimizer' value
  if (!isREML(x)) return(x)
  stopifnot(is(rr <- x@resp, "lmerResp"))
  rho <- new.env(parent=parent.env(environment()))
  rho$resp <- new(class(rr), y=rr$y, offset=rr$offset, weights=rr$weights, REML=0L)
  xpp <- x@pp$copy()
  rho$pp <- new(class(xpp), X=xpp$X, Zt=xpp$Zt, Lambdat=xpp$Lambdat,
                Lind=xpp$Lind, theta=xpp$theta, n=nrow(xpp$X))
  # NOTE: blme changes start
  devfun <- makeRefitDevFun(rho, verbose = verbose, maxit = maxit, control = control, object = x)
  ## NOTE: blme changes end
  
  optwrap <- get("optwrap", lme4Namespace)
  opt <- ## "smart" calc.derivs rules
      if(optimizer == "bobyqa" && !any("calc.derivs" == ...names()))
          optwrap(optimizer, devfun, x@theta, lower=x@lower, calc.derivs=TRUE, ...)
      else
          optwrap(optimizer, devfun, x@theta, lower=x@lower, ...)
  ## FIXME: Should be able to call mkMerMod() here, and be done
  n <- length(rr$y)
  pp <- rho$pp
  p <- ncol(pp$X)
  dims <- c(N=n, n=n, p=p, nmp=n-p, q=nrow(pp$Zt), nth=length(pp$theta),
            useSc=1L, reTrms=length(x@cnms),
            spFe=0L, REML=0L, GLMM=0L, NLMM=0L)#, nAGQ=NA_integer_)
  wrss <- rho$resp$wrss()
  ussq <- pp$sqrL(1)
  pwrss <- wrss + ussq
  cmp <- c(ldL2=pp$ldL2(), ldRX2=pp$ldRX2(), wrss=wrss, ussq=ussq,
           pwrss=pwrss, drsum=NA, dev=opt$fval, REML=NA,
           sigmaML=sqrt(pwrss/n), sigmaREML=sqrt(pwrss/(n-p)))
  ## modify the call  to have REML=FALSE. (without evaluating the call!)
  cl <- x@call
  cl[["REML"]] <- FALSE
  result <- new("lmerMod", call = cl, frame=x@frame, flist=x@flist,
      cnms=x@cnms, theta=pp$theta, beta=pp$delb, u=pp$delu,
      optinfo = get(".optinfo", lme4Namespace)(opt),
      lower=x@lower, devcomp=list(cmp=cmp, dims=dims), pp=pp, resp=rho$resp,
      Gp=x@Gp)
  # NOTE: blme change
  repackageMerMod(result, opt, rho)
}
