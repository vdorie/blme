cat("\n\nRUnit test cases for bmer:::blmer function.\n\n");

test.blme.blmer <- function() {
  source(system.file("common", "lmmData.R", package = "blme"));

  control <- lmerControl(optimizer = "bobyqa");
  
  lmerFit  <-  lmer(y ~ x.1 + x.2 + (1 + x.1 | g.1) + (1 + x.1 + x.2 | g.2), testData, control = control);
  blmerFit <- blmer(y ~ x.1 + x.2 + (1 + x.1 | g.1) + (1 + x.1 + x.2 | g.2), testData, control = control,
                    cov.prior = NULL, fixef.prior = NULL, resid.prior = NULL);

  checkEquals(lmerFit@pp$theta, blmerFit@pp$theta);
  checkEquals(lmerFit@pp$u(1.0), blmerFit@pp$u(1.0));
  checkEquals(lmerFit@pp$beta(1.0), blmerFit@pp$beta(1.0));
}
