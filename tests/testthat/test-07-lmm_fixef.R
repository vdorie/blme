context("blmer numerical results with fixef prior")

source(system.file("common", "lmmData.R", package = "blme"))
lme4Version <- packageVersion("lme4")
control <- lmerControl(optimizer = "bobyqa")

test_that("blmer fits test data with normal(7, TRUE) prior, matching previous version", {
  fixef.prior <- "normal(sd = 7, common.scale = TRUE)"

  startingValues <- c(0.714336877636958, -0.242234853872256, 1.56142829865131, 0.931702840718855, 0.456177995916484, -0.174861679569041, 1.0585277913399, 0.121071648252222, 0.215801873693294)
  result <- if (lme4Version < "1.1-4") c(0.714336904883696, -0.242233333549434, 1.56142849039447, 0.931702729108028, 0.456177204451304, -0.174861811614276, 1.05852821195682, 0.121071547240353, 0.215801842870277) else c(0.714336904883696, -0.242233333549434, 1.56142849039447, 0.931702729108028, 0.456177204451304, -0.174861811614276, 1.05852821195682, 0.121071547240353, 0.215801842870277)
  
  fit <- blmer(y ~ x.1 + x.2 + (1 + x.1 | g.1) + (1 + x.1 + x.2 | g.2), testData, control = control,
                    cov.prior = NULL, fixef.prior = fixef.prior, start = startingValues)
  expect_equal(fit@theta, result, tolerance = 5.0e-5)
})

test_that("blmer fits test data with normal(10, FALSE) prior, matching previous version", {
  fixef.prior <- "normal(sd = 10, common.scale = FALSE)"
  
  startingValues <- list(theta = c(0.705301445472825, -0.236130064856711, 1.54070576284237, 0.919298480793096, 0.444958591085821, -0.162201425613492, 1.04498858978601, 0.121905334663798, 0.204897688209115),
                         sigma = 0.969103097682058)
  result <- if (lme4Version < "1.1-4") c(0.705369855182081, -0.236759905121764, 1.54063251814471, 0.919250008248663, 0.444836570608055, -0.162132239807962, 1.04497528986881, 0.121858574203024, 0.204725931113902) else c(0.705369855182081, -0.236759905121764, 1.54063251814471, 0.919250008248663, 0.444836570608055, -0.162132239807962, 1.04497528986881, 0.121858574203024, 0.204725931113902)
  
  fit <- blmer(y ~ x.1 + x.2 + (1 + x.1 | g.1) + (1 + x.1 + x.2 | g.2), testData, control = control,
                    cov.prior = NULL, fixef.prior = fixef.prior, start = startingValues)
  expect_equal(fit@theta, result, tolerance = 5.0e-5)
  expect_equal(fit@devcomp$cmp[["sigmaREML"]], if (lme4Version < "1.1-4") 0.969074276597577 else 0.969074276597577, tolerance = 1.0e-6)
})
  
test_that("blmer fits test data with t prior, matching previous version", {
  fixef.prior <- "t(3, c(10^2, 2.5^2), common.scale = FALSE)"
  
  startingValues <- list(theta = c(0.645289664330177, -0.151604332140352, 1.39404761930357, 0.788435718441722, 0.312013729923666, -0.0155461916762167, 0.949082870229164, 0.117100582888698, 0),
                         beta = c(5.32508665168687, 1.16859904165051, 4.0443701271478))
  result <- if (lme4Version < "1.1-4") c(0.645289146996319, -0.151634501090343, 1.39403793373549, 0.788432069261316, 0.312010137757441, -0.0155458970707687, 0.949081665570772, 0.117100684805151, 3.13476220325792e-07) else c(0.645289146996319, -0.151634501090343, 1.39403793373549, 0.788432069261316, 0.312010137757441, -0.0155458970707687, 0.949081665570772, 0.117100684805151, 0)
  fixefResult <- if (lme4Version < "1.1-4") c(5.32507818836626, 1.16860398465568, 4.04437041491386) else c(5.32507818836626, 1.16860398465568, 4.04437041491386)
  
  fit <- blmer(y ~ x.1 + x.2 + (1 + x.1 | g.1) + (1 + x.1 + x.2 | g.2), testData, REML = FALSE, control = control,
                    cov.prior = NULL, fixef.prior = fixef.prior, start = startingValues)
  expect_equal(fit@theta, result, tolerance = 5.0e-5)
  expect_equal(fit@beta, fixefResult, tolerance = 5.0e-5)
})

test_that("blmer fits sleep study example in documentation", {
  oldWarnings <- options()$warn
  options(warn = 2)
  
  data("sleepstudy", package = "lme4")
  
  fit <- blmer(Reaction ~ Days + (1 + Days|Subject), sleepstudy,
               cov.prior = NULL, resid.prior = NULL,
               fixef.prior = "normal")
  expect_is(fit, "blmerMod")
  
  fit <- blmer(Reaction ~ Days + (1 + Days|Subject), sleepstudy,
               cov.prior = NULL, resid.prior = NULL,
               fixef.prior = "normal(cov = diag(0.5, 2))")
  expect_is(fit, "blmerMod")
  
  options(warn = oldWarnings)
})
