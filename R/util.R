getCovBlocks <- function(cov, ranefStructure) {
  index <- 0
  result <- list()
  for (i in 1:ranefStructure$numFactors) {
    result[[i]] <- as.matrix(cov[index + 1:ranefStructure$numCoefPerFactor[i],
                                 index + 1:ranefStructure$numCoefPerFactor[i]])
    index <- index + ranefStructure$numRanefPerFactor[i]
  }

  return(result)
}

## from lme4
namedList <- function(...) {
  result <- list(...)
  substituteNames <- sapply(substitute(list(...)), deparse)[-1]
  if (is.null(resultNames <- names(result))) resultNames <- substituteNames
  if (any(noNames <- resultNames == "")) resultNames[noNames] <- substituteNames[noNames]
  setNames(result, resultNames)
}

isFormula <- function(term) is.language(term) && term[[1]] == '~'

#quoteInNamespace <- function(name, character.only = FALSE) {
#  result <- quote(a + b)
#  result[[1]] <- as.symbol(":::")
#  result[[2]] <- as.symbol("blme")
#  
#  result[[3]] <- if (character.only) name else match.call()[[2]]
#  result
#}
