cat("\n\nRUnit test cases for blme::blmer function\n\n");

# simple test to see if blmer can match lmer in output
test.blme.blmer <- function()
{
  generateData <- FALSE;
  testRoot <- file.path(.path.package(package="blme"), "unitTests");  
  source(file.path(testRoot, "lmmData.R"), TRUE);

  lmerFit  <-  lmer(y ~ x.1 + x.2 + (1 + x.1 | g.1) + (1 + x.1 + x.2 | g.2));
  blmerFit <- blmer(y ~ x.1 + x.2 + (1 + x.1 | g.1) + (1 + x.1 + x.2 | g.2),
                    cov.prior = NULL, fixef.prior = NULL, var.prior = NULL);

  # as.numeric to strip names from blmer ST
  # off-by-one bug fix in blmer means that the two are no longer exactly equal
  checkEquals(as.numeric(unlist(lmerFit@ST)), as.numeric(unlist(blmerFit@ST)), tolerance = 1.5e-5);
  checkEquals(lmerFit@ranef, blmerFit@ranef, tolerance=2.0e-5);
  checkEquals(lmerFit@fixef, blmerFit@fixef, tolerance=1.5e-6);
}
