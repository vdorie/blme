cat("\n\nRUnit test cases for bmer::blmer function with priors on the common scale\n\n");

test.bmer.blmer.varPrior <- function()
{
  y <- c(28, 8, -3, 7, -1, 1, 18, 12);
  sigma <- c(15, 10, 16, 11, 9, 11, 10, 18);

  y.z <- (y - mean(y)) / sigma;

  g <- 1:8
  
  model1 <- blmer(y.z ~ 1 + (1 | g), var.prior = "point(1)", cov.prior = NULL, fixef.prior = NULL);
}
