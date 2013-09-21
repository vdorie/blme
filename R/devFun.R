mkBlmerDevfun <- function(fr, X, reTrms, REML = TRUE, start = NULL,
                          verbose = 0L, control = lmerControl(), priors = NULL, ...) {
  devfun <- mkLmerDevfun(fr, X, reTrms, REML, start, verbose, control, ...);

  if (is.null(priors)) priors <- list();
  priors <- evaluatePriorArguments(priors$covPriors, priors$fixefPrior, priors$residPrior,
                                   c(n = nrow(X), p = ncol(X), GLMM = 0L),
                                   reTrms$cnms, parent.frame(2L));
  
  devFunEnv <- environment(devfun);
  devFunEnv$priors <- priors;
  devFunEnv$blmerControl <- createBlmerControl(reTrms, priors);
  devFunBody <- getBlmerDevianceFunctionBody(priors, devFunEnv);

  if (!is.null(devFunBody)) body(devfun) <- parse(text = devFunBody);

  return(devfun);
}

mkBglmerDevfun <- function(fr, X, reTrms, family, nAGQ = 1L, verbose = 0L,
                          control=glmerControl(), priors = NULL, ...) {
  devfun <- mkGlmerDevfun(fr, X, reTrms, family, nAGQ, verbose, control, ...);

  priors <- evaluatePriorArguments(priors$covPriors, priors$fixefPrior, NULL,
                                   c(n = nrow(X), p = ncol(X), GLMM = 1L),
                                   reTrms$cnms, parent.frame(2L));
  devFunEnv <- environment(devfun);
  devFunEnv$priors <- priors;
  devFunEnv$blmerControl <- createBlmerControl(reTrms, priors);
  devFunBody <- getBglmerDevianceFunctionBody(priors, devFunEnv, nAGQ != 0L);

  if (!is.null(devFunBody)) body(devfun) <- parse(text = devFunBody);
  
  return(devfun);
}


## environment already populated at this point
updateBglmerDevfun <- function(devfun, reTrms, nAGQ = 1L) {
  devfun <- updateGlmerDevfun(devfun, reTrms, nAGQ = nAGQ);
  devFunEnv <- environment(devfun);
  devFunBody <- getBglmerDevianceFunctionBody(devFunEnv$priors, devFunEnv, nAGQ != 0L);

  if (!is.null(devFunBody)) body(devfun) <- parse(text = devFunBody);

  return (devfun);
}


getBlmerDevianceFunctionBody <- function(priors, devFunEnv)
{
  if (!requiresPiecewiseOptimization(priors)) return(NULL);

  blmerControl <- devFunEnv$blmerControl;
  
  sigmaOptimizationType <- blmerControl$sigmaOptimizationType;

  fixefPrior <- priors$fixefPrior;

  devFunBody <- NULL;
  stringConnection <- textConnection("devFunBody", "w", local=TRUE);
  sink(stringConnection);

  cat("{\n");

  cat ("  expandPars(theta, pars);\n",
       "  pp$setTheta(as.double(theta));\n\n", sep = "");
  devFunEnv$expandPars <- expandPars;
  
  if (sigmaOptimizationType == SIGMA_OPTIM_POINT)
    cat("  sigma <- priors$residPrior@value;\n");
  
  if (is(fixefPrior, "bmerNormalDist")) {
    if (fixefPrior@commonScale == FALSE) {
      cat("  pp$updateDecomp(sigma * priors$fixefPrior@R.cov.inv);\n");
    } else {
      cat("  pp$updateDecomp(priors$fixefPrior@R.cov.inv);\n");
    }
  } else {
    cat("  pp$updateDecomp();\n");
  }

  cat("\n");

  cat("  resp$updateMu(pp$linPred(0.0));\n",
      "  pp$updateRes(resp$wtres);\n",
      "  pp$solve();\n",
      "  resp$updateMu(pp$linPred(1.0));\n\n", sep = "");

  cat("  beta <- pp$beta(1.0);\n",
      "  Lambda.ts <- getCovBlocks(pp$Lambdat, blmerControl$ranefStructure);\n",
      "  exponentialTerms <- calculatePriorExponentialTerms(priors, beta, Lambda.ts);\n",
      "  polynomialTerm <- calculatePriorPolynomialTerm(priors$covPriors, Lambda.ts);\n\n", sep = "");
  devFunEnv$calculatePriorExponentialTerms <- calculatePriorExponentialTerms;
  devFunEnv$calculatePriorPolynomialTerm <- calculatePriorPolynomialTerm;
  devFunEnv$getCovBlocks <- getCovBlocks;
  
  if (sigmaOptimizationType != SIGMA_OPTIM_NUMERIC &&
      sigmaOptimizationType != SIGMA_OPTIM_POINT) {
    cat("  sigma <- profileSigma(pp, resp, exponentialTerms, blmerControl);\n\n", sep = "");
    devFunEnv$profileSigma <- getSigmaProfiler(priors, blmerControl);
  }

  cat("  lmmObjective(pp, resp, sigma, exponentialTerms, polynomialTerm, blmerControl);\n");
  devFunEnv$lmmObjective <- lmmObjective;
  
  cat("}\n");
  
  sink();
  close(stringConnection);
  
  return(devFunBody);
}

requiresPiecewiseOptimization <- function(priors) {
  !is.null(priors$fixefPrior) || any(sapply(priors$covPriors, function(cov.prior.i) !is.null(cov.prior.i))) ||
  !is.null(priors$residPrior);
}

getSigmaProfiler <- function(priors, blmerControl) {
  sigmaOptimizationType <- blmerControl$sigmaOptimizationType;
  if (sigmaOptimizationType == SIGMA_OPTIM_SQ_LINEAR) {
    return (function(pp, resp, exponentialTerms, blmerControl) {
      pwrss <- resp$wrss() + pp$sqrL(1.0);
      if (!is.null(exponentialTerms[["-2"]])) pwrss <- pwrss + exponentialTerms[["-2"]];

      df <- nrow(pp$X) - resp$REML + blmerControl$df;

      return (sqrt(pwrss / df));
    });
  } else if (sigmaOptimizationType == SIGMA_OPTIM_SQ_QUADRATIC) {
    return (function(pp, resp, exponentialTerms, blmerControl) {
      pwrss <- resp$wrss() + pp$sqrL(1.0);
      if (!is.null(exponentialTerms[["-2"]])) pwrss <- pwrss + exponentialTerms[["-2"]];
      a <- exponentialTerms[["2"]];
      
      df <- nrow(pp$X) - resp$REML + blmerControl$df;

      disc <- sqrt(df^2 + 4 * pwrss * a);
      
      return (sqrt((disc - df) / (2 * a)));
    });
  } else if (sigmaOptimizationType == SIGMA_OPTIM_QUADRATIC) {
    return (function(pp, resp, exponentialTerms, blmerControl) {
      pwrss <- resp$wrss() + pp$sqrL(1.0);
      if (!is.null(exponentialTerms[["-2"]])) pwrss <- pwrss + exponentialTerms[["-2"]];
      a <- exponentialTerms[["-1"]];
      
      df <- nrow(pp$X) - resp$REML + blmerControl$df;

      disc <- sqrt(a^2 + 16 * df * pwrss);
      return (0.25 * (disc + a) / df);
    });
  } else stop("wtf flow control");
}

calculatePriorExponentialTerms <- function(priors, beta, Lambda.ts)
{
  result <- list();
  fixefPrior <- priors$fixefPrior;
  covPriors  <- priors$covPriors;
  residPrior <- priors$residPrior;

  if (!is.null(fixefPrior)) {
    term <- getExponentialTerm(fixefPrior, beta);
    result[[toString(term[1])]] <- term[2];
  }

  for (i in 1:length(covPriors)) {
    if (is.null(covPriors[[i]])) next;
    term <- getExponentialTerm(covPriors[[i]], Lambda.ts[[i]]);
    power <- toString(term[1]);
    exponential <- term[2];
    if (is.null(result[[power]])) result[[power]] <- exponential
    else result[[power]] <- result[[power]] + exponential;
  }

  if (is.null(residPrior)) return(result);
  
  term <- getExponentialTerm(residPrior);
  power <- toString(term[1]);
  exponential <- term[2];
  if (is.null(result[[power]])) result[[power]] <- exponential
  else result[[power]] <- result[[power]] + exponential;

  return (result);
}

calculatePriorPolynomialTerm <- function(covPriors, Lambda.ts)
{
  sum(sapply(1:length(covPriors), function(i)
      if (!is.null(covPriors[[i]])) getPolynomialTerm(covPriors[[i]], Lambda.ts[[i]]) else 0));
}

getBglmerDevianceFunctionBody <- function(priors, devFunEnv, fixefAreParams)
{
  if (!requiresPiecewiseOptimization(priors)) return(NULL);

  fixefPrior <- priors$fixefPrior;

  devFunBody <- NULL;
  stringConnection <- textConnection("devFunBody", "w", local=TRUE);
  sink(stringConnection);

  cat("{\n");

  if (fixefAreParams)
    cat("  resp$setOffset(baseOffset);\n");

  cat("  resp$updateMu(lp0);\n");

  if (fixefAreParams) {
    cat("  pp$setTheta(as.double(pars[dpars]));\n",
        "  spars <- as.numeric(pars[-dpars]);\n",
        "  offset <- if (length(spars) == 0) baseOffset else baseOffset + pp$X %*% spars;\n",
        "  resp$setOffset(offset);\n\n",
        "  p <- pwrssUpdate(pp, resp, tolPwrss, GQmat, compDev, fac, verbose);\n",
        sep = "");
  } else {
    cat("  spars <- rep(0, ncol(pp$X));\n",
        "  pp$setTheta(as.double(theta));\n",
        "  p <- pwrssUpdate(pp, resp, tolPwrss, GHrule(0L), compDev, verbose);\n",
        sep = "");
  }
  
  cat("  resp$updateWts();\n\n",
      
      "  Lambda.ts <- getCovBlocks(pp$Lambdat, blmerControl$ranefStructure);\n",
      "  exponentialTerms <- calculatePriorExponentialTerms(priors, spars, Lambda.ts);\n",
      "  polynomialTerm <- calculatePriorPolynomialTerm(priors$covPriors, Lambda.ts);\n\n",
      
      "  p + exponentialTerms[[1]] + polynomialTerm + blmerControl$constant",
      "}\n",
      sep = "");
  devFunEnv$getCovBlocks <- getCovBlocks;
  devFunEnv$calculatePriorExponentialTerms <- calculatePriorExponentialTerms;
  devFunEnv$calculatePriorPolynomialTerm <- calculatePriorPolynomialTerm;
  
  sink();
  close(stringConnection);
  
  return(devFunBody);
}
