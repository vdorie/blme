covariancePriorsToString <- function(covPriors, numGroupsPerFactor, digits)
{
  result <- character(0L)
  resultIndex <- 1L
  
  numFactors <- length(numGroupsPerFactor)
  factorNames <- names(numGroupsPerFactor)
  for (i in seq_len(numFactors)) {
    prior.i <- covPriors[[i]]

    if (is.null(prior.i)) next

    result[resultIndex] <- paste(factorNames[i], " ~ ", toString(prior.i, digits), sep = "")
    resultIndex <- resultIndex + 1L
  }

  result
}

printPriors <- function(priors, numGroupsPerFactor, digits) {
  covariancePriorOutput <- covariancePriorsToString(priors$covPriors, numGroupsPerFactor, digits)
  if (length(covariancePriorOutput) > 0L) {
    cat("Cov prior  : ", covariancePriorOutput[1L], "\n", sep="")
    if (length(covariancePriorOutput) > 1L) {
      for (i in seq.int(2L, length(covariancePriorOutput)))
        cat("           : ", covariancePriorOutput[i], "\n", sep="")
    }
  }
  if (!is.null(priors$fixefPrior))
    cat("Fixef prior: ", toString(priors$fixefPrior, digits), "\n", sep="")

  if (!is.null(priors$residPrior))
    cat("Resid prior: ", toString(priors$residPrior, digits, FALSE), "\n", sep="")
}
