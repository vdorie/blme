covariancePriorsToString <- function(covPriors, numGroupsPerFactor, digits)
{
  result <- character(0)
  resultIndex <- 1
  
  numFactors <- length(numGroupsPerFactor)
  factorNames <- names(numGroupsPerFactor)
  for (i in 1:numFactors) {
    prior.i <- covPriors[[i]]

    if (is.null(prior.i)) next

    result[resultIndex] <- paste(factorNames[i], " ~ ", toString(prior.i, digits), sep = "")
    resultIndex <- resultIndex + 1
  }

  result
}

printPriors <- function(priors, numGroupsPerFactor, digits) {
  covariancePriorOutput <- covariancePriorsToString(priors$covPriors, numGroupsPerFactor, digits)
  if (length(covariancePriorOutput) > 0) {
    cat("Cov prior  : ", covariancePriorOutput[1], "\n", sep="")
    if (length(covariancePriorOutput) > 1) {
      for (i in 2:length(covariancePriorOutput))
        cat("           : ", covariancePriorOutput[i], "\n", sep="")
    }
  }
  if (!is.null(priors$fixefPrior))
    cat("Fixef prior: ", toString(priors$fixefPrior, digits), "\n", sep="")

  if (!is.null(priors$residPrior))
    cat("Resid prior: ", toString(priors$residPrior, digits, FALSE), "\n", sep="")
}
