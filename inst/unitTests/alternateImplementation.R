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

blme_rotateSparseDesign <- function(model) {
  # except for rounding errors, the following are equivalent
  .Call("mer_ST_setPars", model, blme_stMatricesToVector(model@ST));
  A <- model@A;
  
#  P.ranef <- blme_ranefPerm(model);
#  Lambda <- P.ranef %*% blme_ranefChol(model) %*% t(P.ranef);
#  A <- t(Lambda) %*% model@Zt;

  if (length(model@sqrtXWt) == 0) {
    Cx <- A@x;
  } else {
    temp <- A %*% blme_getWeightMatrix(model);
    Cx <- temp@x;
  }
  
  return(list(A = A, Cx = Cx));
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
  modelCopy@deviance[["NULLdev"]] <- 0.0; # just something to get it to do a deep copy
  .Call("mer_update_L", modelCopy);
  L <- modelCopy@L;
#  L <- Cholesky(tcrossprod(P %*% C), Imult=1, LDL=FALSE, perm=FALSE);
#  L@perm <- model@L@perm;
#  L@type[1] <- as.integer(2);
  
  RZX <- as(solve(L, P %*% C %*% X.w, "L"), "matrix");

  priorTypeNone     <- blme:::getEnumOrder(blme:::typeEnum, blme:::NONE_TYPE_NAME);
  priorTypeDirect   <- blme:::getEnumOrder(blme:::typeEnum, blme:::DIRECT_TYPE_NAME);
  priorFamilyNormal <- blme:::getEnumOrder(blme:::familyEnum, blme:::NORMAL_FAMILY_NAME);

  if (model@fixef.prior@type == priorTypeNone) {
    temp <- crossprod(X.w) - crossprod(RZX);
    if (any(is.nan(temp@x))) browser();
    RX <- as(chol(temp), "matrix");
  } else if (model@fixef.prior@type == priorTypeDirect &&
             model@fixef.prior@families[1] == priorFamilyNormal) {
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

    onCommonScale <- blme:::getScaleInt(blme:::getEnumOrder(blme:::posteriorScaleEnum, blme:::defaultUnmodeledCoefficientPosteriorScale),
                                        blme:::getEnumOrder(blme:::commonScaleEnum, blme:::COMMON_SCALE_TRUE_NAME));
    if (model@fixef.prior@scales[1] != onCommonScale)
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
  beta.tilde  <- solve(t(model@RX), crossprod(X.w, Y.w) - crossprod(model@RZX, theta.tilde));
    
  beta.hat  <- as(solve(model@RX, beta.tilde), "numeric");
  theta.hat <- as(solve(model@L, theta.tilde - model@RZX %*% beta.hat, "Lt"), "numeric");
  
  return(list(beta.hat = beta.hat, theta.hat = theta.hat,
              beta.tilde = beta.tilde, theta.tilde = theta.tilde));
}

blme_calculateDeviances <- function(model, modes, stParameters) {
  C <- model@A;

  W <- blme_getWeightMatrix(model);
  if (length(model@sqrtXWt) > 0) C@x <- model@Cx;
    
  Y.w <- W %*% model@y;
  X.w <- W %*% model@X;
    
  P <- blme_cholmodPerm(model);
    
  theta.hat <- model@u;
  beta.hat <- model@fixef;

  totalSumOfSquares <- (crossprod(Y.w) - crossprod(modes$theta.tilde) - crossprod(modes$beta.tilde))[1];

  deviance <- model@deviance;
  deviance[["usqr"]] <- crossprod(model@u)[1];

  priorDF <- blme_getPriorDegreesOfFreedom(model);
  
  degreesOfFreedomML   <- model@dims[["n"]];
  degreesOfFreedomREML <- model@dims[["n"]] -  model@dims[["p"]];
  
  if (blme_canProfileCommonScale(model)) {
    priorExpPart <- blme_getCommonScaleExponentialPart(model, stParameters);

    deviance[["sigmaML"]]   <- sigmaML   <- sqrt((totalSumOfSquares + 2.0 * priorExpPart$mTwo) / (degreesOfFreedomML + priorDF));
    deviance[["sigmaREML"]] <- sigmaREML <- sqrt((totalSumOfSquares + 2.0 * priorExpPart$mTwo) / (degreesOfFreedomREML + priorDF));
  } else {
    sigmaML   <- deviance[["sigmaML"]];
    sigmaREML <- deviance[["sigmaREML"]];
  }
    
  sigmaML.sq   <- sigmaML^2;
  sigmaREML.sq <- sigmaREML^2;

  deviance[["wrss"]] <- crossprod(Y.w - X.w %*% modes$beta.hat - crossprod(P %*% C, modes$theta.hat))[1];
  pwrss <- deviance[["wrss"]] + deviance[["usqr"]];
  deviance[["ML"]]   <- deviance[["ldL2"]] +
    degreesOfFreedomML   * log(2 * pi *   sigmaML.sq) + pwrss / sigmaML.sq;
  deviance[["REML"]] <- deviance[["ldL2"]] + deviance[["ldRX2"]] +
    degreesOfFreedomREML * log(2 * pi * sigmaREML.sq) + pwrss / sigmaREML.sq;

  return(deviance);
}

blme_getObjectiveFunction <- function(model) {
  parameters <- blme_stMatricesToVector(model@ST);
  return(blme_getObjectiveFunctionForParameters(parameters, model));
}

blme_getObjectiveFunctionForParameters <- function(parameters, model) {
  if (blme_parametersIncludeCommonScale(model)) {
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
  if (any(sapply(model@ST, function(ST) any(diag(ST) < 0)))) return(.Machine$double.xmax * .Machine$double.eps);
  sparseParts <- blme_rotateSparseDesign(model);
  model@A <- sparseParts$A;
  model@Cx <- sparseParts$Cx;
  
  factorizations <- blme_computeDesignFactors(model);
  model@L <- factorizations$L;
  model@RZX <- factorizations$RZX;
  model@RX <- factorizations$RX;
  
  model@deviance[["ldL2"]]  <- 2.0 * determinant(model@L,  logarithm=TRUE)$modulus;
  model@deviance[["ldRX2"]] <- 2.0 * determinant(model@RX, logarithm=TRUE)$modulus;

  modes <- blme_calculateJointMode(model);
  
  model@fixef <- modes$beta.hat;
  model@u     <- modes$theta.hat;

  model@deviance <- blme_calculateDeviances(model, modes, stParameters);

  result <- ifelse(model@dims[["REML"]], model@deviance[["REML"]], model@deviance[["ML"]]);
  result <- result + blme_getPriorPenalty(model);

  if (is.nan(result) || !is.finite(result)) return(.Machine$double.xmax * .Machine$double.eps);
  return(result);
}
