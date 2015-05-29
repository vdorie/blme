## "parInfo" is a named list, each member also being a list containing
##   length
##   default - function of rho, aka the environment of the dev fun
##   lower   - lower bounds for par

## the order of elements in it is the order that they are passed to the optimizer
## not that that should matter to a downstream user. if they care, construct it by
## index

## for models with a single parameter vector, the name "theta" is sufficient
## so we can return without further analysis
expandParsInCurrentFrame <- function(parVector, parInfo) {
  if (length(parInfo) == 1) return(invisible(NULL))
  
  parentEnv <- parent.frame()
  parNames <- names(parInfo)
  
  index <- 0
  for (i in 1:length(parInfo)) {
    parLength <- parInfo[[i]]$length
    parName <- parNames[[i]]

    parentEnv[[parName]] <- parVector[index + 1:parLength]
    index <- index + parLength
  }
  invisible(NULL)
}

getStartingValues <- function(userStart, devFunEnv, parInfo) {
  if (is.null(userStart)) userStart <- list()
  if (is.numeric(userStart)) userStart <- list(theta = userStart)
  if (is.list(userStart) && length(userStart) == 1 &&
      is.null(names(userStart))) names(userStart) <- "theta"
  
  invalidStartingValues <- !(names(userStart) %in% names(parInfo))
  if (any(invalidStartingValues))
    warning("starting values for parameter(s) '", toString(names(userStart)[invalidStartingValues]),
            "' not part of model and will be ignored")

  start <- numeric(sum(sapply(parInfo, function(par.i) par.i$length)))
  offset <- 0L
  for (i in 1:length(parInfo)) {
    parName <- names(parInfo)[[i]]
    parLength <- parInfo[[i]]$length

    userValue <- userStart[[parName]]
    useDefault <- TRUE
    if (!is.null(userValue)) {
      if (length(userValue) != parLength) {
        warning("parameter '", parName, "' is of length ", parLength, ", yet supplied vector is of length ",
                length(userValue), ". start will be ignored")
      } else {
        start[offset + 1:parLength] <- userValue
        useDefault <- FALSE
      }
    }
    if (useDefault)
      start[offset + 1:parLength] <- parInfo[[i]]$default(devFunEnv)
    
    offset <- offset + parLength
  }
  start
}

extractParameterListFromFit <- function(fit, blmerControl) {
  result <- list(theta = fit@theta)
  if (blmerControl$fixefOptimizationType == FIXEF_OPTIM_NUMERIC) {
    if (fit@devcomp$dims[["GLMM"]] != 0L)
      result$fixef <- fit@beta
    else
      result$beta  <- fit@beta
  }
  if (fit@devcomp$dims[["GLMM"]] == 0L && blmerControl$fixefOptimizationType == SIGMA_OPTIM_NUMERIC) {
    result$sigma <- if (fit@devcomp$dims[["REML"]] == 0L) fit@devcomp$cmp[["sigmaML"]] else fit@devcomp$cmp[["sigmaREML"]]
  }
  result
}

getLowerBounds <- function(parInfo) {
  result <- numeric(sum(sapply(parInfo, function(par.i) par.i$length)))
  offset <- 0L
  for (i in 1:length(parInfo)) {
    parName <- names(parInfo)[[i]]
    parLength <- parInfo[[i]]$length
    parLower <- parInfo[[i]]$lower
    if (parLength != length(parLower)) stop("length of lower bounds for parameter '", parName, "' does not equal length of vector")

    result[offset + 1:parLength] <- parLower
    offset <- offset + parLength
  }
  
  result
}

getParInfo <- function(pred, resp, ranefStructure, blmerControl) {
  numPars <- 1
  result <- list(theta = list(length = ranefStructure$numCovParameters,
                   lower = ranefStructure$lower,
                   default = function(devFunEnv) pred$theta))
  
  if (blmerControl$fixefOptimizationType == FIXEF_OPTIM_NUMERIC) {
    numPars <- numPars + 1
    numFixef <- if (length(pred$X) > 0) ncol(pred$X) else 0
    result[[numPars]] <-
      list(length = numFixef,
             lower = rep(-Inf, numFixef),
             default = function(devFunEnv) pred$beta0 + pred$delb)
    names(result)[[numPars]] <- "beta"
  }
  if (blmerControl$sigmaOptimizationType == SIGMA_OPTIM_NUMERIC) {
    numPars <- numPars + 1
    result[[numPars]] <-
      list(length = 1L,
             lower = 0,
             default = function(devFunEnv) sd(resp$y))
    names(result)[[numPars]] <- "sigma"
  }
  
  result
}
