getGLMMData <- function() {
  set.seed(3, "Mersenne-Twister", "Inversion");
  J <- 4;
  n <- 8;
  N <- J * n;
  x.1 <- rnorm(N);
  x.2 <- rnorm(N);
  theta <- rnorm(J, 0, 2);
  theta.g <- rep(theta, rep(n, J));
  eta <- 3 * x.1 + 2 * x.2 + theta.g;
  mu <- exp(eta) / (1 + exp(eta));
  y <- rbinom(N, 1, mu);

  g <- gl(J, n);

  return(data.frame(y = y, x.1 = x.1, x.2 = x.2, g = g));
}
testData <- getGLMMData();
rm(getGLMMData);
