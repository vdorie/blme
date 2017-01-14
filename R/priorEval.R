evaluateFixefPrior <- function(fixefPrior, defnEnv, evalEnv) {
  if (is.character(fixefPrior)) fixefPrior <- parse(text = fixefPrior)[[1L]]
  
  if (is.symbol(fixefPrior) && exists(toString(fixefPrior), envir = evalEnv) &&
      !(as.character(fixefPrior) %in% fixefDistributions)) {
    fixefPrior <- get(toString(fixefPrior), envir = evalEnv)
    if (is.character(fixefPrior)) fixefPrior <- parse(text = fixefPrior)[[1L]]
  }
  
  if (!is.null(fixefPrior)) {
    if (is.symbol(fixefPrior)) fixefPrior <- call(as.character(fixefPrior))
    fixefDistributionName <- as.character(fixefPrior[[1L]])
    if (!(fixefDistributionName %in% fixefDistributions)) stop("unrecognized fixef distribution: '", fixefDistributionName, "'")

    return(eval(fixefPrior, envir = evalEnv))
  }
  
  NULL
}

evaluateCovPriors <- function(covPriors, factorColumnNames, numGroupsPerFactor, defnEnv, evalEnv) {
  numFactors <- length(factorColumnNames)
  factorNames <- names(factorColumnNames)
  result <- vector("list", numFactors)
  defaultCovPrior <- NULL

  if (is.null(covPriors)) return(result)
  
  if (is.character(covPriors)) {
    covPriors <- gsub("inverse.wishart", "invwishart", covPriors)
    covPriors <- gsub("inverse.gamma", "invgamma", covPriors)
    covPriors <- parse(text = covPriors)[[1L]]
  }
  
  if (is.call(covPriors) && covPriors[[1L]] == "list") covPriors[[1L]] <- NULL
  
  if (!is.list(covPriors)) covPriors <- list(covPriors)

  for (i in seq_along(covPriors)) {
    covPrior.i <- covPriors[[i]]
    ## can't just let 'em re-define "wishart", or use the built-in gamma
    if (is.symbol(covPrior.i) && exists(toString(covPrior.i), envir = evalEnv) &&
        !(as.character(covPrior.i) %in% covDistributions)) {
      covPrior.i <- get(toString(covPrior.i), envir = evalEnv)
      if (is.character(covPrior.i)) covPrior.i <- parse(text = covPrior.i)[[1]]
      covPriors[[i]] <- covPrior.i
    }
  }
  for (i in seq_along(covPriors)) {
    covPrior.i <- covPriors[[i]]
    if (is.character(covPrior.i)) {
      covPrior.i <- gsub("inverse.wishart", "invwishart", covPrior.i)
      covPrior.i <- gsub("inverse.gamma", "invgamma", covPrior.i)
      covPrior.i <- parse(text = covPrior.i)[[1]]
    }

    ## turn 'wishart' into 'wishart()'
    if (is.symbol(covPrior.i)) covPrior.i <- call(as.character(covPrior.i))
    
    if (is.formula(covPrior.i)) {
      factorName <- as.character(covPrior.i[[2L]])

      if (!(factorName %in% factorNames))
        stop("grouping factor '", factorName, "' for covariance prior not in model formula")
      
      ## turn 'group ~ wishart' into 'group ~ wishart()'
      if (is.symbol(covPrior.i[[3L]])) covPrior.i[[3L]] <- call(as.character(covPrior.i[[3L]]))

      ## for each grouping factor with the given name, store function call for later
      matchingFactors <- which(factorName == factorNames)
      for (j in seq_along(matchingFactors)) result[[matchingFactors[j]]] <- covPrior.i[[3L]]
      
    } else {
      ## default
      if (!is.null(defaultCovPrior)) warning("more than one default covariance prior specified, only using the last one")
      defaultCovPrior <- covPrior.i
    }
  }
  
  for (i in seq_len(numFactors)) {
    if (is.null(result[[i]]) && is.null(defaultCovPrior)) next

    result.i <- result[[i]]
    if (is.null(result[[i]]) && !is.null(defaultCovPrior)) result.i <- defaultCovPrior

    covDistributionName <- as.character(result.i[[1L]])
    if (!(covDistributionName %in% covDistributions)) stop("unrecognized ranef covariance distribution: '", covDistributionName, "'")

    defnEnv$q.k <- defnEnv$level.dim <- length(factorColumnNames[[i]])
    defnEnv$j.k <- defnEnv$n.grps <- numGroupsPerFactor[i]
    
    result.i <- eval(result.i, envir = evalEnv)
    
    if (!is.null(result.i)) result[[i]] <- result.i
  }
  
  result
}

evaluateResidualPrior <- function(residPrior, defnEnv, evalEnv) {
  if (is.character(residPrior)) {
    residPrior <- gsub("inverse.gamma", "invgamma", residPrior)
    residPrior <- parse(text = residPrior)[[1L]]
  }

  if (is.symbol(residPrior) && exists(toString(residPrior), envir = evalEnv)) {
    fixefPrior <- get(toString(residPrior), envir = evalEnv)
    if (is.character(residPrior)) residPrior <- parse(text = residPrior)[[1]]
  }
  
  if (!is.null(residPrior)) {
    if (is.symbol(residPrior)) residPrior <- call(as.character(residPrior))
    residDistributionName <- as.character(residPrior[[1L]])
    if (!(residDistributionName %in% residDistributions)) stop("unrecognized residual variance distribution: '", residDistributionName, "'")

    return(eval(residPrior, envir = evalEnv))
  }
  
  NULL
}
  
evaluatePriorArguments <- function(covPriors, fixefPrior, residPrior,
                                   dims, factorColumnNames, numGroupsPerFactor, parentEnv) {
  result <- list()
  evalEnv <- new.env(parent = parentEnv)
  defnEnv <- new.env()
  
  defnEnv$p <- defnEnv$n.fixef <- dims[["p"]]
  defnEnv$n <- defnEnv$n.obs   <- dims[["n"]]

  isLMM <- dims[["GLMM"]] == 0L
  
  ## add the names of dist functs to the evaluating env
  for (distributionName in names(lmmDistributions)) {
    distributionFunction <- lmmDistributions[[distributionName]]

    environment(distributionFunction) <- defnEnv
    if (!isLMM) {
      ## need both copies to have their envs tweaked, but only one called
      distributionFunction <- glmmDistributions[[distributionName]]
      if (!is.null(distributionFunction)) environment(distributionFunction) <- defnEnv
    }
    if (!is.null(distributionFunction)) assign(distributionName, distributionFunction, envir = evalEnv)
  }
  
  result$fixefPrior <- evaluateFixefPrior(fixefPrior, defnEnv, evalEnv)
  if ((is(result$fixefPrior, "bmerTDist") || is(result$fixefPrior, "bmerHorseshoeDist")) && isLMM && dims[["REML"]] > 0L)
    stop("t/horseshoe distribution for fixed effects only supported when REML = FALSE")
  result$covPriors  <- evaluateCovPriors(covPriors, factorColumnNames, numGroupsPerFactor, defnEnv, evalEnv)

  if (isLMM) {
    environment(residualVarianceGammaPrior) <- defnEnv
    environment(residualVarianceInvGammaPrior) <- defnEnv
    assign("gamma", residualVarianceGammaPrior, envir = evalEnv)
    assign("invgamma", residualVarianceInvGammaPrior, envir = evalEnv)
    
    result$residPrior <- evaluateResidualPrior(residPrior, defnEnv, evalEnv)
  }
  
  result
}
