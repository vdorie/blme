blme_getBlockCov <- function(blockList, numRepsPerBlock) {
  numFactors <- length(blockList);
  factorDims <- sapply(blockList, nrow);

  # suppose blocks can be decomps, but ignore that for now
  numNonZeroes <- sum(factorDims^2 * numRepsPerBlock);
  rowIndices <- rep(0, numNonZeroes);
  colIndices <- rep(0, numNonZeroes);
  values <- rep(0, numNonZeroes);

  sparseIndex <- 1;
  upperLeftIndex <- 1;
  for (i in 1:numFactors) {
    numParams <- factorDims[i];
    numValues <- numParams^2;

    for (j in 1:numRepsPerBlock[i]) {
      sparseRange <- sparseIndex + 1:numValues - 1;
      upperLeftRange <- upperLeftIndex + 1:numParams - 1;

      rowIndices[sparseRange] <- rep(upperLeftRange, numParams);
      colIndices[sparseRange] <- rep(upperLeftRange, rep(numParams, numParams));
      values[sparseRange] <- as.vector(blockList[[i]]);

      sparseIndex <- sparseIndex + numValues;
      upperLeftIndex <- upperLeftIndex + numParams;
    }
  }

  return(sparseMatrix(rowIndices, colIndices, x = values));
}

blme_stToCov <- function(ST) {
  dimension <- nrow(ST);
  T <- ST;
  diag(T) <- rep(1, dimension);
  S <- diag(diag(ST), dimension);
  return(tcrossprod(T %*% S));
}

blme_stToCovInv <- function(ST) {
  dimension <- nrow(ST);
  T <- ST;
  diag(T) <- rep(1, dimension);
  S.inv <- diag(1 / diag(ST), dimension);
  return(crossprod(S.inv %*% solve(T)));
}

blme_stToChol <- function(ST) {
  dimension <- nrow(ST);
  T <- ST;
  diag(T) <- rep(1, dimension);
  S <- diag(diag(ST), dimension);
  return(T %*% S);
}

blme_covToSdCor <- function(cov) {
  sds <- diag(sqrt(diag(cov)));
  sds.inv <- diag(1 / sqrt(diag(cov)));
  cor <- sds.inv %*% cov %*% sds.inv;
  return(list(sds = sds, cor = cor));
}

blme_covToSt <- function(cov) {
  TS <- t(chol(cov));
  S.inv <- diag(1 / diag(TS));
  ST <- TS %*% S.inv;
  diag(ST) <- diag(TS);
  return(ST);
}


blme_stMatricesToVector <- function(ST) {
  result <- rep(0, sum(sapply(ST, function(ST.i) { n <- nrow(ST.i); n * (n + 1) / 2 })));
  offset <- 0;
  for (i in 1:length(ST)) {
    factorDimension <- nrow(ST[[i]]);
    
    result[offset + 1:factorDimension] <- diag(ST[[i]]);
    offset <- offset + factorDimension;
    
    if (factorDimension == 1) next;
    
    lowerTriangleLength <- factorDimension * (factorDimension - 1) / 2;
    result[offset + 1:lowerTriangleLength] <- ST[[i]][lower.tri(ST[[i]])];
    offset <- offset + lowerTriangleLength;
  }
  
  return(result);
}

blme_stVectorToMatrix <- function(stAsVector, factorDimension) {
  result <- diag(stAsVector[1:factorDimension], factorDimension);
  if (factorDimension > 1) {
    lowerTriangleLength <- factorDimension * (factorDimension - 1) / 2;
    result[lower.tri(result)] <- stAsVector[factorDimension + 1:lowerTriangleLength];
  }
  return(result);
}
  

blme_stVectorToMatrices <- function(parameters, numFactors, factorDimensions) {
  ST <- list();
  offset <- 0;
  for (i in 1:numFactors) {
    factorDimension <- factorDimensions[i];
    stVectorLength <- factorDimension * (factorDimension + 1) / 2;
    
    ST[[i]] <- blme_stVectorToMatrix(parameters[offset + 1:stVectorLength], factorDimension);
    offset <- offset + stVectorLength;
  }
  
  return(ST);
}

blme_ranefVcov <- function(model) {
  Sigmas <- lapply(model@ST, blme_stToCov);

  factorDims <- sapply(model@ST, nrow);
  numRepsPerBlock <- (model@Gp[-1] - model@Gp[-length(model@Gp)]) / factorDims;

  return(blme_getBlockCov(Sigmas, numRepsPerBlock));
}

blme_ranefVcovInv <- function(model) {
  Lambdas <- lapply(model@ST, blme_stToChol);
  Lambdas.inv <- lapply(Lambdas, solve);
  Sigmas.inv <- lapply(Lambdas.inv, crossprod);

  factorDims <- sapply(model@ST, nrow);
  numRepsPerBlock <- (model@Gp[-1] - model@Gp[-length(model@Gp)]) / factorDims;

  return(blme_getBlockCov(Sigmas.inv, numRepsPerBlock));
}

blme_ranefChol <- function(model) {
  Lambdas <- lapply(model@ST, blme_stToChol);
  factorDims <- sapply(model@ST, nrow);
  
  numRepsPerBlock <- (model@Gp[-1] - model@Gp[-length(model@Gp)]) / factorDims;

  return(blme_getBlockCov(Lambdas, numRepsPerBlock));
}

blme_logDetVcov <- function(model) {
  numFactors <- model@dims[["nt"]];
  
  Lambdas <- lapply(model@ST, blme_stToChol);
  
  factorDims <- sapply(model@ST, nrow);
  numRepsPerBlock <- (model@Gp[-1] - model@Gp[-length(model@Gp)]) / factorDims;

  result <- 0;
  for (i in 1:numFactors) {
    factorDet <- sum(log(diag(Lambdas[[i]])));
    result <- result + factorDet * numRepsPerBlock[i];
  }
  return(2.0 * result);
}

# the permutation matrix used internally, produced by cholmod.
# fill-reduces
blme_cholmodPerm <- function(model) {
  as(model@L@perm + 1, "pMatrix");
}

# permutes the block diagonal form into what lmer uses
blme_ranefPerm <- function(model)
{
  numFactors <- model@dims[["nt"]];
  factorDims <- sapply(model@ST, nrow);
  numRepsPerBlock <- (model@Gp[-1] - model@Gp[-length(model@Gp)]) / factorDims;
  
  indices <- rep(0, sum(factorDims * numRepsPerBlock));

  index <- 1;
  offset <- 0;
  for (k in 1:numFactors) {
    for (l in 1:factorDims[k]) {
      for (j in 1:numRepsPerBlock[k]) {
        indices[index] <- offset + l + (j - 1) * factorDims[k];
        index <- index + 1;
      }
    }
    offset <- offset + numRepsPerBlock[k] * factorDims[k];
  }
  return(as(indices, "pMatrix"));
}
