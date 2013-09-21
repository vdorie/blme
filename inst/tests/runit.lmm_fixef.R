cat("\n\nRUnit test cases for bmer::blmer function with priors on the fixed effects\n\n");

test.blme.blmer.fixefPrior <- function()
{
  source(system.file("common", "lmmData.R", package = "blme"));
  
  fixef.prior <- "normal(sd = 1, common.scale = TRUE)";
  ## the default optimizers seem to be sucking for this model
  ## under the assumption that they'll improve in the future,
  ## start it at the mode as computed by the old version
  startingValues <- c(0.648196454969222, -0.251155871182205, 1.41982711952805, 5.30818893266426, 1.64487709178204, 2.42515000251243, 0.934438088183127, 0.218132399295426, 0);
  fit <- blmer(y ~ x.1 + x.2 + (1 + x.1 | g.1) + (1 + x.1 + x.2 | g.2), testData,
                    cov.prior = NULL, fixef.prior = fixef.prior, start = startingValues);
  checkEquals(fit@pp$theta, c(0.648195706527889, -0.251155765769507, 1.41982725882128, 5.30817842042119, 1.6448714780241, 2.42514919244429, 0.93443679163898, 0.218131815825292, 8.25408519583847e-07));


  fixef.prior <- "normal(sd = 1, common.scale = FALSE)";
  startingValues <- list(theta = c(0.696498273543799, -0.173433928372341, 1.49823772191009, 5.25514420225188,
                                   1.54058807242222, 3.37965560881278, 1.04257820385898, -0.0394352712532918,
                                   0.854358283868515),
                         sigma = 0.970090843027481);
  fit <- blmer(y ~ x.1 + x.2 + (1 + x.1 | g.1) + (1 + x.1 + x.2 | g.2), testData,
                    cov.prior = NULL, fixef.prior = fixef.prior, start = startingValues);
  checkEquals(c(fit@pp$theta, fit@devcomp$cmp[["sigmaREML"]]),
              c(0.696497912161098, -0.173434135295791, 1.49823739668995, 5.25513473016839, 1.54058911169323, 3.37964801540207, 1.04257692653163, -0.0394332012063832, 0.854358710427337, 0.970090843027481));
  
  
  # check even that these run without error
  options(warn = 2);
  ignored <- blmer(Reaction ~ Days + (1 + Days|Subject), sleepstudy,
                   cov.prior = NULL, resid.prior = NULL,
                   fixef.prior = "normal");
  ignored <- blmer(Reaction ~ Days + (1 + Days|Subject), sleepstudy,
                   cov.prior = NULL, resid.prior = NULL,
                   fixef.prior = "normal(cov = diag(0.5, 2))");
  options(warn = 0);
}
