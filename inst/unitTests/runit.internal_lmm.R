cat("\n\nRUnit test cases for blme::bmer_lmmTest function\n\n");

test.blme.lmm.internals <- function()
{
  .Call("bmer_lmmTest");
}

# this block will calculate the rotated sparse model matrix for
# a specific mer object
if (FALSE) {
getBlockCovariance <- function(Sigmas, totalNumModeledParameters, numGroups) {
  numLevels <- length(Sigmas);

  numNonZeroEntries <- sum(sapply(1:numLevels, function(k) { return(dim(Sigmas[[k]])[1]^2 * numGroups[k]) }));
  rowIndices <- rep(0, numNonZeroEntries);
  colIndices <- rep(0, numNonZeroEntries);
  values <- rep(0, numNonZeroEntries);

  sparseIndex <- 1;
  upperLeftIndex <- 1;
  for (k in 1:numLevels) {
    numModeledParameters <- dim(Sigmas[[k]])[1];
    numValues <- numModeledParameters^2;

    for (j in 1:numGroups[k]) {
      rowIndices[sparseIndex:(sparseIndex + numValues - 1)] <-
        rep(upperLeftIndex:(upperLeftIndex + numModeledParameters - 1),
            numModeledParameters);
      colIndices[sparseIndex:(sparseIndex + numValues - 1)] <-
        as.vector(t(matrix(rep(upperLeftIndex:(upperLeftIndex + numModeledParameters - 1)),
                           numModeledParameters, numModeledParameters)));
      values[sparseIndex:(sparseIndex + numValues - 1)] <-
        as.vector(Sigmas[[k]]);
      
      sparseIndex <- sparseIndex + numValues;
      upperLeftIndex <- upperLeftIndex + numModeledParameters;
    }
  }
  
  return(sparseMatrix(rowIndices, colIndices, x=values));
}

getPermutationFromLmerToBlock <- function(numFactors, numGroupsPerFactor, numModeledParametersPerFactor)
{
  indices <- rep(0, sum(numGroupsPerFactor * numModeledParametersPerFactor));

  index <- 1;
  offset <- 0;
  for (k in 1:numFactors) {
    for (l in 1:numModeledParametersPerFactor[k]) {
      for (j in 1:numGroupsPerFactor[k]) {
        indices[index] <- offset + l + (j - 1) * numModeledParametersPerFactor[k];
        index <- index + 1;
      }
    }
    offset <- offset + numModeledParametersPerFactor[k] * numGroupsPerFactor[k];
  }
  return(indices);
}

rotateSparseDesignMatrix <- function(model) {
  factorDimensions <- sapply(model@ST, nrow);
  numGroupsPerFactor <- (model@Gp[-1] - model@Gp[-length(model@Gp)]) / factorDimensions;
  
  S <- lapply(model@ST, function(matrix) diag(diag(matrix), nrow(matrix)));
  T <- lapply(model@ST, function(matrix) { diag(matrix) <- rep(1, nrow(matrix)); matrix });

  Lambdas <- lapply(1:length(S), function(i) T[[i]] %*% S[[i]])
  Lambda <- getBlockCovariance(Lambdas, sum(factorDimensions * numGroupsPerFactor), numGroupsPerFactor)
  perm <- as(getPermutationFromLmerToBlock(length(numGroupsPerFactor), numGroupsPerFactor,
                                           factorDimensions), "pMatrix")

  return(perm %*% t(Lambda) %*% t(perm) %*% model@Zt);
}

computeAugmentedDesignFactorizations <- function(model) {
  W <- Diagonal(nrow(model@X), model@sqrtXWt);
  X.w <- W %*% model@X
  C <- model@A %*% W
  P <- as(model@L@perm + 1, "pMatrix")

  L <- Cholesky(tcrossprod(P %*% C), Imult=1, LDL=FALSE, perm=FALSE);
  L@perm <- model@L@perm;
  L@type[1] <- as.integer(2);

  RZX <- as(solve(L, P %*% C %*% X.w, "L"), "matrix");
  RZ <- as(chol(crossprod(X.w) - crossprod(RZX)), "matrix");

  return(list(L = L, RZX = RZX, RZ = RZ));
}

calculateJointMode <- function(model) {
  C <- model@A;
  
  if (length(model@sqrtXWt) > 0) {
    W <- Diagonal(nrow(model@X), model@sqrtXWt);
    C@x <- model@Cx;
  } else {
    W <- Diagonal(nrow(model@X));
  }

  Y.w <- W %*% model@y;
  X.w <- W %*% model@X;
  
  P <- as(model@L@perm + 1, "pMatrix");
  
  theta.tilde <- solve(model@L, P %*% C %*% Y.w, "L");
  beta.tilde <- solve(t(model@RX), crossprod(X.w, Y.w) - crossprod(model@RZX, theta.tilde));

  beta.hat <- as(solve(model@RX, beta.tilde), "numeric");
  theta.hat <- as(solve(model@L, theta.tilde - model@RZX %*% beta.hat, "Lt"), "numeric");

  pwrss <- (crossprod(Y.w) - (crossprod(theta.tilde) +
                              crossprod(beta.tilde)))[1];
  
  return(list(beta.hat = beta.hat, theta.hat = theta.hat,
              beta.tilde = beta.tilde, theta.tilde = theta.tilde,
              pwrss = pwrss));
}

calculateCommonScale <- function(model) {
  C <- model@A;
  
  if (length(model@sqrtXWt) > 0) {
    W <- Diagonal(nrow(model@X), model@sqrtXWt);
    C@x <- model@Cx;
  } else {
    W <- Diagonal(nrow(model@X));
  }
    
  Y.w <- W %*% model@y;
  X.w <- W %*% model@X;

  P <- as(model@L@perm + 1, "pMatrix");

  theta.hat <- model@u;
  beta.hat <- model@fixef;
  
  weightedResidualSumOfSquares <- crossprod(Y.w - crossprod(P %*% C, theta.hat) - X.w %*% beta.hat)[1];
  penaltySumOfSquares <- crossprod(model@u)[1];

  numObservations <- nrow(X.w);
  deviance <- numObservations * (1 + log(weightedResidualSumOfSquares + penaltySumOfSquares) +
                                 log(2 * pi / numObservations)) + 2 * determinant(model@L, logarithm = TRUE)$modulus;
  sigma.hat <- sqrt((weightedResidualSumOfSquares + penaltySumOfSquares) / numObservations);

  return (list(weightedResidualSumOfSquares = weightedResidualSumOfSquares,
          penaltySumOfSquares = penaltySumOfSquares,
          deviance = deviance,
          sigma.hat = sigma.hat));
}
 
}
