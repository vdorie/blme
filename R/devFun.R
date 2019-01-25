mkBlmerDevfun <- function(fr, X, reTrms, REML = TRUE, start = NULL,
                          verbose = 0L, control = lmerControl(), priors = NULL, ...) {
  devfun <- mkLmerDevfun(fr, X, reTrms, REML, start, verbose, control, ...)
  devFunEnv <- environment(devfun)
  pred <- devFunEnv$pp
  resp <- devFunEnv$resp
  
  devFunEnv$ranefStructure <- getRanefStructure(pred, resp, reTrms)
  
  if (is.null(priors)) priors <- list()
  devFunEnv$priors <-
    evaluatePriorArguments(priors$covPriors, priors$fixefPrior, priors$residPrior,
                           c(n = nrow(X), p = ncol(X), GLMM = 0L, REML = if (REML) 1L else 0L),
                           colnames(X),
                           reTrms$cnms, devFunEnv$ranefStructure$numGroupsPerFactor,
                           parent.frame(2L))
  
  devFunEnv$blmerControl <- createBlmerControl(pred, resp, devFunEnv$priors)
  devFunEnv$parInfo <- getParInfo(pred, resp, devFunEnv$ranefStructure, devFunEnv$blmerControl)
  devFunBody <- getBlmerDevianceFunctionBody(devFunEnv)

  if (!is.null(devFunBody)) body(devfun) <- parse(text = devFunBody)

  devfun
}

mkBglmerDevfun <- function(fr, X, reTrms, family, nAGQ = 1L, verbose = 0L,
                           maxit = 100L,
                           control=glmerControl(),
                           priors = NULL, ...) {
  devfun <-
    if (packageVersion("lme4") <= "1.1.7") {
     mkGlmerDevfun(fr, X, reTrms, family, nAGQ, verbose, control, ...)
    } else {
     mkGlmerDevfun(fr, X, reTrms, family, nAGQ, verbose, maxit, control, ...)
    }
  
  devFunEnv <- environment(devfun)
  pred <- devFunEnv$pp
  resp <- devFunEnv$resp

  devFunEnv$ranefStructure <- getRanefStructure(pred, resp, reTrms)
  if (is.null(priors)) priors <- list()
  devFunEnv$priors <-
    evaluatePriorArguments(priors$covPriors, priors$fixefPrior, NULL,
                           c(n = nrow(X), p = ncol(X), GLMM = 1L),
                           colnames(X),
                           reTrms$cnms, devFunEnv$ranefStructure,
                           parent.frame(2L))

  
  devFunEnv$blmerControl <- createBlmerControl(pred, resp, devFunEnv$priors)
  devFunEnv$parInfo <- getParInfo(pred, resp, devFunEnv$ranefStructure, devFunEnv$blmerControl)

  devFunBody <- getBglmerDevianceFunctionBody(devFunEnv, nAGQ != 0L)

  if (!is.null(devFunBody)) body(devfun) <- parse(text = devFunBody)
  
  devfun
}

makeRefitDevFun <- function(env, nAGQ = 1L, verbose = 0, maxit=100L,
                            control = list(),
                            object) {
  lme4Namespace <- asNamespace("lme4")
  devfun <-
    if (packageVersion("lme4") <= "1.1.7") {
      get("mkdevfun", lme4Namespace)(env, nAGQ, verbose, control)
    } else
      get("mkdevfun", lme4Namespace)(env, nAGQ, maxit, verbose, control)
  
  pred <- env$pp
  resp <- env$resp
  
  env$maxit <- as.integer(maxit)
  env$lower <- object@lower
  env$priors <- object@priors
  env$ranefStructure <- getRanefStructure(pred, resp, list(cnms = object@cnms, Gp = object@Gp))
  env$blmerControl <- createBlmerControl(pred, resp, env$priors)
  env$parInfo <- getParInfo(pred, resp, env$ranefStructure, env$blmerControl)

  devFunBody <-
    if (is(resp, "lmerResp"))
      getBlmerDevianceFunctionBody(env)
    else
      getBglmerDevianceFunctionBody(env, nAGQ != 0L)
  
  if (!is.null(devFunBody)) body(devfun) <- parse(text = devFunBody)
  
  devfun
}

## environment already populated at this point
updateBglmerDevfun <- function(devfun, reTrms, nAGQ = 1L) {
  devfun <- updateGlmerDevfun(devfun, reTrms, nAGQ = nAGQ)
  devFunEnv <- environment(devfun)
  devFunBody <- getBglmerDevianceFunctionBody(devFunEnv, nAGQ != 0L)

  if (!is.null(devFunBody)) body(devfun) <- parse(text = devFunBody)

  devfun
}


getBlmerDevianceFunctionBody <- function(devFunEnv)
{
  priors <- devFunEnv$priors
  
  if (!anyPriorsApplied(priors)) return(NULL)

  blmerControl <- devFunEnv$blmerControl
  
  sigmaOptimizationType <- blmerControl$sigmaOptimizationType
  fixefOptimizationType <- blmerControl$fixefOptimizationType

  fixefPrior <- priors$fixefPrior

  devFunBody <- NULL
  stringConnection <- textConnection("devFunBody", "w", local=TRUE)
  sink(stringConnection)

  cat("{\n")
  cat("  expandParsInCurrentFrame(theta, parInfo)\n",
      "  pp$setTheta(as.double(theta))\n\n", sep = "")
  devFunEnv$expandParsInCurrentFrame <- expandParsInCurrentFrame
  
  if (sigmaOptimizationType == SIGMA_OPTIM_POINT)
    cat("  sigma <- priors$residPrior@value\n")

  if (is(fixefPrior, "bmerNormalDist")) {
    if (fixefPrior@commonScale == FALSE) {
      cat("  pp$updateDecomp(sigma * priors$fixefPrior@R.cov.inv)\n")
    } else {
      cat("  pp$updateDecomp(priors$fixefPrior@R.cov.inv)\n")
    }
  } else {
    cat("  pp$updateDecomp()\n")
  }

  cat("\n")

  cat("  resp$updateMu(pp$linPred(0.0))\n",
      "  pp$updateRes(resp$wtres)\n",
      "  pp$solve()\n",
      "  resp$updateMu(pp$linPred(1.0))\n\n", sep = "")

  if (fixefOptimizationType != FIXEF_OPTIM_NUMERIC) {
    cat("  beta <- pp$beta(1.0)\n")
  }
  cat("  Lambda.ts <- getCovBlocks(pp$Lambdat, ranefStructure)\n")
  if (sigmaOptimizationType == SIGMA_OPTIM_NUMERIC ||
      sigmaOptimizationType == SIGMA_OPTIM_POINT) {
    cat("  exponentialTerms <- calculatePriorExponentialTerms(priors, beta, Lambda.ts, sigma)\n")
  } else {
    cat("  exponentialTerms <- calculatePriorExponentialTerms(priors, beta, Lambda.ts)\n")
  }
  cat("  polynomialTerm <- calculatePriorPolynomialTerm(priors$covPriors, Lambda.ts)\n\n")
  devFunEnv$calculatePriorExponentialTerms <- calculatePriorExponentialTerms
  devFunEnv$calculatePriorPolynomialTerm <- calculatePriorPolynomialTerm
  devFunEnv$getCovBlocks <- getCovBlocks

  if (fixefOptimizationType == FIXEF_OPTIM_NUMERIC) {
    cat("  exponentialTerms <- calculateFixefExponentialTerm(beta, pp$beta(1.0), pp$RX(), exponentialTerms)\n")
    devFunEnv$calculateFixefExponentialTerm <- calculateFixefExponentialTerm
  }
  
  if (sigmaOptimizationType != SIGMA_OPTIM_NUMERIC &&
      sigmaOptimizationType != SIGMA_OPTIM_POINT) {
    cat("  sigma <- profileSigma(pp, resp, exponentialTerms, blmerControl)\n\n", sep = "")
    devFunEnv$profileSigma <- getSigmaProfiler(priors, blmerControl)
  }

  cat("  lmmObjective(pp, resp, sigma, exponentialTerms, polynomialTerm, blmerControl)\n")
  devFunEnv$lmmObjective <- lmmObjective
  
  cat("}\n")
  
  sink()
  close(stringConnection)
  
  devFunBody
}

anyPriorsApplied <- function(priors) {
  !is.null(priors$fixefPrior) || any(sapply(priors$covPriors, function(cov.prior.i) !is.null(cov.prior.i))) ||
  !is.null(priors$residPrior)
}

getSigmaProfiler <- function(priors, blmerControl) {
  sigmaOptimizationType <- blmerControl$sigmaOptimizationType
  if (sigmaOptimizationType == SIGMA_OPTIM_SQ_LINEAR) {
    return (function(pp, resp, exponentialTerms, blmerControl) {
      pwrss <- resp$wrss() + pp$sqrL(1.0)
      if (!is.null(exponentialTerms[["-2"]])) pwrss <- pwrss + exponentialTerms[["-2"]]

      df <- nrow(pp$X) - resp$REML + blmerControl$df

      sqrt(pwrss / df)
    })
  } else if (sigmaOptimizationType == SIGMA_OPTIM_SQ_QUADRATIC) {
    return (function(pp, resp, exponentialTerms, blmerControl) {
      pwrss <- resp$wrss() + pp$sqrL(1.0)
      if (!is.null(exponentialTerms[["-2"]])) pwrss <- pwrss + exponentialTerms[["-2"]]
      a <- exponentialTerms[["2"]]
      
      df <- nrow(pp$X) - resp$REML + blmerControl$df

      disc <- sqrt(df^2 + 4 * pwrss * a)
      
      sqrt((disc - df) / (2 * a))
    })
  } else if (sigmaOptimizationType == SIGMA_OPTIM_QUADRATIC) {
    return (function(pp, resp, exponentialTerms, blmerControl) {
      pwrss <- resp$wrss() + pp$sqrL(1.0)
      if (!is.null(exponentialTerms[["-2"]])) pwrss <- pwrss + exponentialTerms[["-2"]]
      a <- exponentialTerms[["-1"]]
      
      df <- nrow(pp$X) - resp$REML + blmerControl$df

      disc <- sqrt(a^2 + 16 * df * pwrss)
      
      0.25 * (disc + a) / df
    })
  } else stop("illegal sigma optimization type")
}

calculatePriorExponentialTerms <- function(priors, beta, Lambda.ts, sigma = NULL)
{
  result <- list()
  fixefPrior <- priors$fixefPrior
  covPriors  <- priors$covPriors
  residPrior <- priors$residPrior

  if (!is.null(fixefPrior)) {
    if ((is(fixefPrior, "bmerTDist") || is(fixefPrior, "bmerHorseshoeDist")) && fixefPrior@commonScale == TRUE) {
      term <- getExponentialTerm(fixefPrior, beta, sigma)
    } else {
      term <- getExponentialTerm(fixefPrior, beta)
    }
    result[[toString(term[1L])]] <- term[2L]
  }

  for (i in seq_along(covPriors)) {
    if (is.null(covPriors[[i]])) next
    covPrior.i <- covPriors[[i]]

    if (is(covPrior.i, "bmerCustomDist") && covPrior.i@commonScale == FALSE) {
      term <- getExponentialTerm(covPrior.i, Lambda.ts[[i]] * sigma)
    } else {
      term <- getExponentialTerm(covPrior.i, Lambda.ts[[i]])
    }
    power <- toString(term[1L])
    exponential <- term[2L]
    if (is.null(result[[power]])) result[[power]] <- exponential
    else result[[power]] <- result[[power]] + exponential
  }

  if (is.null(residPrior)) return(result)
  
  term <- getExponentialTerm(residPrior)
  power <- toString(term[1L])
  exponential <- term[2L]
  if (is.null(result[[power]])) result[[power]] <- exponential
  else result[[power]] <- result[[power]] + exponential

  result
}

calculatePriorPolynomialTerm <- function(covPriors, Lambda.ts)
{
  sum(sapply(seq_along(covPriors), function(i)
      if (!is.null(covPriors[[i]])) getPolynomialTerm(covPriors[[i]], Lambda.ts[[i]]) else 0))
}

calculateFixefExponentialTerm <- function(beta, beta.tilde, RX, exponentialTerms = NULL)
{
  exponential <- crossprod(RX %*% (beta - beta.tilde))[1L]
  if (is.null(exponentialTerms)) return(exponential)
  
  if (is.null(exponentialTerms[["-2"]])) {
    exponentialTerms[["-2"]] <- exponential
  } else {
    exponentialTerms[["-2"]] <- exponentialTerms[["-2"]] + exponential
  }
  exponentialTerms
}

testGetBglmerDevianceFunctionBody <- function(devFun)
{
  devFunEnv <- environment(devFun)
  priors <- devFunEnv$priors
  
  if (!anyPriorsApplied(priors)) return(NULL)

  fixefPrior <- priors$fixefPrior
  
  devFunBody <- NULL
  stringConnection <- textConnection("devFunBody", "w", local = TRUE)
  sink(stringConnection)

  cat("{\n",
      "  Lambda.ts        <- getCovBlocks(pp$Lambdat, ranefStructure)\n",
      "  exponentialTerms <- calculatePriorExponentialTerms(priors, spars, Lambda.ts)\n",
      "  polynomialTerm   <- calculatePriorPolynomialTerm(priors$covPriors, Lambda.ts)\n\n",
      "  ",
      sep = "")
  
  oldDevFunBody <- deparse(body(devFun))
  cat(oldDevFunBody[-length(oldDevFunBody)], sep = "\n") ## cut trailing "}"

  cat("  } + exponentialTerms[[1]] + polynomialTerm + blmerControl$constant\n",
      "}", sep = "")

  sink()
  close(stringConnection)
  
  devFunBody
}

getBglmerDevianceFunctionBody <- function(devFunEnv, fixefAreParams)
{
  priors <- devFunEnv$priors
  
  if (!anyPriorsApplied(priors)) return(NULL)

  fixefPrior <- priors$fixefPrior

  devFunBody <- NULL
  stringConnection <- textConnection("devFunBody", "w", local=TRUE)
  sink(stringConnection)

  cat("{\n")

  if (fixefAreParams)
    cat("  resp$setOffset(baseOffset)\n")

  cat("  resp$updateMu(lp0)\n")

  if (!fixefAreParams) {
    cat("  spars <- rep(0, ncol(pp$X))\n",
        "  pp$setTheta(as.double(theta))\n", sep = "")
    if (packageVersion("lme4") <= "1.1.7") {
      cat("  p <- pwrssUpdate(pp, resp, tolPwrss, GHrule(0L), compDev, verbose=verbose)\n")
    } else {
      cat("  p <- pwrssUpdate(pp, resp, tolPwrss, GHrule(0L), compDev, maxit=maxit, verbose=verbose)\n")
    }
  } else {
    cat("  pp$setTheta(as.double(pars[dpars]))\n",
        "  spars <- as.numeric(pars[-dpars])\n",
        "  offset <- if (length(spars) == 0) baseOffset else baseOffset + pp$X %*% spars\n",
        "  resp$setOffset(offset)\n\n", sep = "")
    if (packageVersion("lme4") <= "1.1.7") {
      cat("  p <- pwrssUpdate(pp, resp, tolPwrss, GQmat, compDev, fac, verbose=verbose)\n")
    } else {
      cat("  p <- pwrssUpdate(pp, resp, tolPwrss, GQmat, compDev, fac, maxit=maxit, verbose=verbose)\n")
    }
  }
  
  cat("  resp$updateWts()\n\n",
      
      "  Lambda.ts <- getCovBlocks(pp$Lambdat, ranefStructure)\n",
      "  exponentialTerms <- calculatePriorExponentialTerms(priors, spars, Lambda.ts)\n",
      "  polynomialTerm <- calculatePriorPolynomialTerm(priors$covPriors, Lambda.ts)\n\n",
      
      "  p + exponentialTerms[[1]] + polynomialTerm + blmerControl$constant\n",
      "}\n",
      sep = "")
  devFunEnv$getCovBlocks <- getCovBlocks
  devFunEnv$calculatePriorExponentialTerms <- calculatePriorExponentialTerms
  devFunEnv$calculatePriorPolynomialTerm <- calculatePriorPolynomialTerm
  
  sink()
  close(stringConnection)
  
  devFunBody
}
