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

blme_getInverseInformation <- function(model) {
  # Xt.sp = sparse part of design matrix, transposed
  # X.dn  = dense part
  # transformations will be *not* applied, aside from weighting
  # observations
  # Xt.sp = L' Z' dmu/deta W^0.5
  # X.dn  = W^0.5 dmu/deta X
  #
  # dmu/deta is 1 if linear, is what it is for glm (depends on link)
  # W is 1 if linear, 1 / var(y) if glm

  Xt.sp <- model@A;
  X.dn  <- model@X;
  if (length(model@sqrtXWt) > 0) {
    if (all(dim(model@Cm)) > 0) {
      Xt.sp <- model@Cm;
    } else {
      Xt.sp@x <- model@Cx;
    }
    X.dn <- diag(as.vector(model@sqrtXWt)) %*% X.dn;
  }

  # at this point, we have the components of the augmented design matrix that
  # we need to take the crossproduct of, and then invert.

  # block-wise matrix inversion;
  # http://en.wikipedia.org/wiki/Invertible_matrix#Blockwise_inversion
  # A, B, C, and D correspond to there, although C = B' so we don't need
  # to create it
  #
  # To figure out where this comes from, you can check the docs for
  # blmer, there should be some calculations there.
  A <- tcrossprod(Xt.sp) + Diagonal(nrow(Xt.sp));
  B <- Xt.sp %*% X.dn;
  D <- crossprod(X.dn);

  A.inv <- solve(A);

  # X.i here are blocks of the inverse
  temp <- crossprod(B, A.inv); # CA^-1, comes up a lot
  D.i <- solve(D - temp %*% B);
  C.i <- -1 * D.i %*% temp;
  A.i <- A.inv - crossprod(temp, C.i);

  # since everything above is "spherical", we have to bring it to the
  # appropriate scale/rotation
  P.ranef <- blme_ranefPerm(model);
  Lambda <- P.ranef %*% blme_ranefChol(model) %*% t(P.ranef);
  
  A.i <- Lambda %*% tcrossprod(A.i, Lambda);
  C.i <- tcrossprod(C.i, Lambda);

  numRanef <- model@dims[["q"]];
  numFixef <- model@dims[["p"]];
  result <- matrix(NA, numRanef + numFixef, numRanef + numFixef);

  ranefRange <- 1:numRanef;
  fixefRange <- numRanef + 1:numFixef;
  
  result[ranefRange, ranefRange] <- as.matrix(A.i);
  result[fixefRange, ranefRange] <- as.matrix(C.i);
  result[ranefRange, fixefRange] <- t(result[fixefRange, ranefRange]);
  result[fixefRange, fixefRange] <- as.matrix(D.i);

  return(result);
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

blme_stVectorToMatrices <- function(parameters, numFactors, factorDimensions) {
  ST <- list();
  offset <- 0;
  for (i in 1:numFactors) {
    factorDimension <- factorDimensions[i];
    
    ST.i <- diag(parameters[offset + 1:factorDimension], factorDimension);
    offset <- offset + factorDimension;
    
    if (factorDimension > 1) {
      lowerTriangleLength <- factorDimension * (factorDimension - 1) / 2;
        ST.i[lower.tri(ST.i)] <- parameters[offset + 1:lowerTriangleLength];
      offset <- offset + lowerTriangleLength;
    }
    
    ST[[i]] <- ST.i;
  }
  
  return(ST);
}

blme_rotateSparseDesign <- function(model) {
  # except for rounding errors, the following are equivalent
  .Call("mer_ST_setPars", model, blme_stMatricesToVector(model@ST));
  return(model@A);
#  P.ranef <- blme_ranefPerm(model);
#  Lambda <- P.ranef %*% blme_ranefChol(model) %*% t(P.ranef);
  
#  return(t(Lambda) %*% model@Zt);
}

blme_getWeightMatrix <- function(model) {
  if (length(model@sqrtXWt) > 0) {
    W <- Diagonal(nrow(model@X), model@sqrtXWt);
  } else {
    W <- Diagonal(nrow(model@X));
  }
  return(W);
}

blme_computeDesignFactors <- function(model) {
  W <- blme_getWeightMatrix(model);
  
  X.w <- W %*% model@X;
  C <- model@A %*% W;
  P <- blme_cholmodPerm(model);

  # except for rounding errors, the following are equivalent
  modelCopy <- model;
  modelCopy@deviance[["sigmaML"]] <- 0.0; # just something to get it to do a deep copy
  .Call("mer_update_L", modelCopy);
  L <- modelCopy@L;
#  L <- Cholesky(tcrossprod(P %*% C), Imult=1, LDL=FALSE, perm=FALSE);
#  L@perm <- model@L@perm;
#  L@type[1] <- as.integer(2);
  
  RZX <- as(solve(L, P %*% C %*% X.w, "L"), "matrix");
  
  if (model@fixef.prior@type == 0) {
    RX <- as(chol(crossprod(X.w) - crossprod(RZX)), "matrix");
  } else {
    RXPartial <- crossprod(X.w) - crossprod(RZX);
    
    sigma.sq <- ifelse(model@dims[["REML"]],
                       model@deviance[["sigmaREML"]], model@deviance[["sigmaML"]])^2;
    
    # first param is log det - not needed here
    hyperparameters <- model@fixef.prior@hyperparameters[-1];
    if (length(hyperparameters) == 1) {
      Sigma.beta.inv <- diag(hyperparameters^2, ncol(model@X))
    } else if (length(hyperparameters) == ncol(model@X)) {
      Sigma.beta.inv <- diag(hyperparameters^2);
    } else {
      matrixLength <- ncol(model@X)^2;
      Sigma.beta.inv <- matrix(hyperparameters[1:matrixLength + matrixLength],
                               ncol(model@X), ncol(model@X));
    }
    if (model@fixef.prior@scales[1] == 2)
      Sigma.beta.inv <- sigma.sq * Sigma.beta.inv;
    RX <- as(chol(RXPartial + Sigma.beta.inv), "matrix");
  }
  
  return(list(L = L, RZX = RZX, RX = RX));
}

blme_calculateJointMode <- function(model) {
  C <- model@A;

  W <- blme_getWeightMatrix(model);
  if (length(model@sqrtXWt) > 0) C@x <- model@Cx;
  
  Y.w <- W %*% model@y;
  X.w <- W %*% model@X;
    
  P <- blme_cholmodPerm(model);

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

blme_calculateDeviances <- function(model) {
  C <- model@A;

  W <- blme_getWeightMatrix(model);
  if (length(model@sqrtXWt) > 0) C@x <- model@Cx;
    
  Y.w <- W %*% model@y;
  X.w <- W %*% model@X;
    
  P <- blme_cholmodPerm(model);
    
  theta.hat <- model@u;
  beta.hat <- model@fixef;

  deviance <- model@deviance;
  pwrss <- deviance[["pwrss"]];
  deviance[["usqr"]] <- crossprod(model@u)[1];
  deviance[["wrss"]] <- deviance[["disc"]] <- pwrss - deviance[["usqr"]];
    
  degreesOfFreedomML   <- nrow(X.w);
  degreesOfFreedomREML <- nrow(X.w) - ncol(X.w);

  
  if (model@fixef.prior@type == 0) {
    deviance[["sigmaML"]]   <- sqrt(pwrss / degreesOfFreedomML);
    deviance[["sigmaREML"]] <- sigmaREML <- sqrt(pwrss / degreesOfFreedomREML);
    
    deviance[["ML"]]   <- deviance[["ldL2"]] +
      degreesOfFreedomML   * (1 + log(pwrss) + log(2 * pi / degreesOfFreedomML));
    deviance[["REML"]] <- deviance[["ldL2"]] + deviance[["ldRX2"]] +
      degreesOfFreedomREML * (1 + log(pwrss) + log(2 * pi / degreesOfFreedomREML));
  } else {
    sigmaML.sq   <- deviance[["sigmaML"]]^2;
    sigmaREML.sq <- deviance[["sigmaREML"]]^2;
    
    deviance[["ML"]]   <- deviance[["ldL2"]] +
      degreesOfFreedomML   * log(2 * pi *   sigmaML.sq) + pwrss / sigmaML.sq;
    deviance[["REML"]] <- deviance[["ldL2"]] + deviance[["ldRX2"]] +
      degreesOfFreedomREML * log(2 * pi * sigmaREML.sq) + pwrss / sigmaREML.sq;
  }

  return(deviance);
}

blme_deviance <- function(parameters, model) {
  if (model@fixef.prior@type[1] != 0) {
    stParameters <- parameters[-length(parameters)];
    if (model@dims[["REML"]]) {
      model@deviance[["sigmaREML"]] <- parameters[length(parameters)];
    } else {
      model@deviance[["sigmaML"]]   <- parameters[length(parameters)];
    }
  } else {
    stParameters <- parameters;
  }
  
  model@ST <- blme_stVectorToMatrices(stParameters, model@dims[["nt"]], sapply(model@ST, nrow));
  model@A <- blme_rotateSparseDesign(model);
  factorizations <- blme_computeDesignFactors(model);
  model@L <- factorizations$L;
  model@RZX <- factorizations$RZX;
  model@RX <- factorizations$RX;
  
  model@deviance[["ldL2"]]  <- 2.0 * determinant(model@L,  logarithm=TRUE)$modulus;
  model@deviance[["ldRX2"]] <- 2.0 * determinant(model@RX, logarithm=TRUE)$modulus;
  
  modes <- blme_calculateJointMode(model);
  
  model@fixef <- modes$beta.hat;
  model@u     <- modes$theta.hat;
  model@deviance[["pwrss"]] <- modes$pwrss;
  
  model@deviance <- blme_calculateDeviances(model);
  
  if (model@dims[["REML"]]) return(model@deviance[["REML"]]);
  return(model@deviance[["ML"]]);
}
