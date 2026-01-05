## "parInfo" is a named list, each member also being a list containing
##   length
##   default - function of rho, aka the environment of the dev fun
##   lower   - lower bounds for par

## the order of elements in it is the order that they are passed to the optimizer
## not that that should matter to a downstream user. if they care, construct it by
## index

expandParsInCurrentFrame <- function(parVector, parInfo) {
  parentEnv <- parent.frame()
  parNames <- names(parInfo)
  
  offset <- 0
  for (i in seq_along(parInfo)) {
    parLength <- parInfo[[i]]$length
    parName <- parNames[[i]]

    parentEnv[[parName]] <- parVector[offset + seq_len(parLength)]
    offset <- offset + parLength
  }
  invisible(NULL)
}

getStartingValues <- function(userStart, devFunEnv, parInfo) {
  if (is.null(userStart)) userStart <- list()
  if (is.numeric(userStart)) userStart <- list(theta = userStart)
  if (is.list(userStart) && length(userStart) == 1 &&
      is.null(names(userStart))) names(userStart) <- "par"

  if (length(parPos <- which(names(userStart) == "par")) >= 1)
  {
    names(userStart)[parPos] <- "theta"
  }
  if (length(fixefPos <- which(names(userStart) == "fixef")) >= 1)
    names(userStart)[fixefPos] <- "beta"
  
  invalidStartingValues <- !(names(userStart) %in% names(parInfo))
  if (any(invalidStartingValues))
    warning("starting values for parameter(s) '", toString(names(userStart)[invalidStartingValues]),
            "' not part of model and will be ignored")

  start <- numeric(sum(sapply(parInfo, function(par.i) par.i$length)))
  offset <- 0L
  for (i in seq_along(parInfo)) {
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
  lme4Version <- packageVersion("lme4")
  result <- if (lme4Version >= "2.0.0") list(par = fit@theta) else list(theta = fit@theta)
  if (blmerControl$fixefOptimizationType == FIXEF_OPTIM_NUMERIC) {
    if (fit@devcomp$dims[["GLMM"]] != 0L)
      result$beta <- fit@beta
    else
      result$beta <- fit@beta
  }
  if (fit@devcomp$dims[["GLMM"]] == 0L && blmerControl$fixefOptimizationType == SIGMA_OPTIM_NUMERIC) {
    result$sigma <- if (fit@devcomp$dims[["REML"]] == 0L) fit@devcomp$cmp[["sigmaML"]] else fit@devcomp$cmp[["sigmaREML"]]
  }
  result
}

getBounds <- function(parInfo, direction) {
  result <- numeric(sum(sapply(parInfo, function(par.i) par.i$length)))
  offset <- 0L
  for (i in seq_along(parInfo)) {
    parName <- names(parInfo)[[i]]
    parLength <- parInfo[[i]]$length
    bound <- parInfo[[i]][[direction]]
    if (parLength != length(bound)) {
      stop(
        "length of ", direction, " bounds for parameter '", parName,
        "' does not equal length of vector"
      )
    }

    result[offset + 1:parLength] <- bound
    offset <- offset + parLength
  }
  
  result
}

getParInfo <- function(pred, resp, ranefStructure, blmerControl) {
  numPars <- 1L
  lme4Version <- packageVersion("lme4")
  result <- list(
    theta=list(
      length=ranefStructure$numCovParameters,
      lower=ranefStructure$lower,
      default=function(devFunEnv) pred$theta
    )
  )

  if (lme4Version >= "2.0.0")
    result$theta$upper <- ranefStructure$upper

  if (blmerControl$fixefOptimizationType == FIXEF_OPTIM_NUMERIC) {
    numPars <- numPars + 1L
    numFixef <- if (length(pred$X) > 0L) ncol(pred$X) else 0L
    result[[numPars]] <- list(
      length=numFixef,
      lower=rep(-Inf, numFixef),
      default=function(devFunEnv) pred$beta0 + pred$delb
    )
    if (lme4Version >= "2.0.0")
      result[[numPars]]$upper <- rep(Inf, numFixef)
    names(result)[[numPars]] <- "beta"
  }
  if (blmerControl$sigmaOptimizationType == SIGMA_OPTIM_NUMERIC) {
    numPars <- numPars + 1L
    result[[numPars]] <- list(
      length=1L,
      lower=0,
      default=function(devFunEnv) sd(resp$y)
    )
    if (lme4Version >= "2.0.0")
      result[[numPars]]$upper <- Inf
    names(result)[[numPars]] <- "sigma"
  }

  result
}

