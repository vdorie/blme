getLMMData <- function() {
  set.seed(0, "Mersenne-Twister", "Inversion");
  
  N <- 100;
  J.1 <- 5;
  J.2 <- 5;
  beta <- c(5, 2, 4);
  theta.1 <- matrix(rnorm(J.1 * 2), J.1, 2);
  theta.2 <- matrix(rnorm(J.2 * 3), J.2, 3);
  
  x.1 <- rnorm(N);
  x.2 <- rnorm(N);
  g.1 <- rmultinom(N, 1, runif(J.1));
  g.2 <- rmultinom(N, 1, runif(J.2));
  g.1 <- sapply(1:N, function(i) which(g.1[,i] == 1));
  g.2 <- sapply(1:N, function(i) which(g.2[,i] == 1));
  
  y <- 1   * (beta[1] + theta.1[g.1,1] + theta.2[g.2,1]) +
       x.1 * (beta[2] + theta.1[g.1,2] + theta.2[g.2,2]) +
       x.2 * (beta[3] +                  theta.2[g.2,3]) +
       rnorm(N);
  
  weights <- runif(N);
  weights <- weights / sum(weights);
  
  return(data.frame(y = y, x.1 = x.1, x.2 = x.2, g.1 = g.1, g.2 = g.2, w = weights));
}
testData <- getLMMData();
rm(getLMMData);
