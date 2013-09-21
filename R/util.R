getCovBlocks <- function(cov, ranefStructure) {
  index <- 0;
  result <- list();
  for (i in 1:ranefStructure$numFactors) {
    result[[i]] <- as.matrix(cov[index + 1:ranefStructure$numCoefPerFactor[i],
                                 index + 1:ranefStructure$numCoefPerFactor[i]]);
    index <- index + ranefStructure$numRanefPerFactor[i];
  }

  return(result);
}
