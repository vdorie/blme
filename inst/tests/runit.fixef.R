cat("\n\nRUnit test cases for blme:::parsePrior function, fixef.prior argument\n\n");


test.blme.parsePrior.fixef.prior <- function()
{
  source(system.file("common", "lmmData.R", package = "blme"));
  source(system.file("common", "checkWarning.R", package = "blme"));

  model1 <- blmer(y ~ x.1 + (1 | g.1), testData,
                  cov.prior = NULL, fixef.prior = NULL, var.prior = NULL);
  
  # exceptions and warnings
  RUnitOptions <- getOption("RUnit");
  RUnitOptions$silent <- TRUE;
  options("RUnit" = RUnitOptions);
  
  checkException(blme:::parsePrior(model1, fixef.prior = "normal(common.scale = 'crazy')"));
  checkException(blme:::parsePrior(model1, fixef.prior = "normal(cov = diag(3))"));
  negDefiniteMatrix <- matrix(c(1, 0, 0, -0.1), 2, 2);
  checkException(blme:::parsePrior(model1, fixef.prior = "normal(cov = negDefiniteMatrix)"));
  asymmetricMatrix <- matrix(c(1, 0.5, 0.3, 0.7), 2, 2);
  checkException(blme:::parsePrior(model1, fixef.prior = "normal(cov = asymmetricMatrix)"));
}
