# defaults
loadCovarianceDefaults <- function(env) {
  env$defaultSpectralPosteriorScale <- VAR_SCALE_NAME;
  env$defaultSpectralDataScale <- FREE_SCALE_NAME;
  env$defaultSpectralFamily <- GAMMA_FAMILY_NAME;
  env$defaultSpectralGammaShape <-
    expression(ifelse(posteriorScale == SD_SCALE_NAME, 2.0, 1.5));
  env$defaultSpectralGammaRate <-
    expression(ifelse(shape > 1,
                      (shape - 1) / ifelse(posteriorScale == SD_SCALE_NAME,
                                           10^2, 10^4),
                      shape / ifelse(posteriorScale == SD_SCALE_NAME,
                                     10^2, 10^4)));
  env$defaultSpectralInverseGammaShape <- 0.5;
  env$defaultSpectralInverseGammaScale <- 
    expression(ifelse(posteriorScale == SD_SCALE_NAME, 10^2, 10^4) *
               (shape + 1));

  # correlation priors
  # note that since the correlation is S * R * S, if the 's' components peak around
  # the sqrt(sd(y)) and the 'r' components peak around sd(y), all multiplied together
  # gives a rough peak near the data variance
  env$defaultCorrelationPosteriorScale <- SD_SCALE_NAME;
  env$defaultCorrelationDataScale <- FREE_SCALE_NAME;
  env$defaultCorrelationCoordinateFamily <- GAMMA_FAMILY_NAME;
  env$defaultCorrelationMatrixFamily <- WISHART_FAMILY_NAME;
  env$defaultCorrelationGammaShape <-
    expression(ifelse(posteriorScale == SD_SCALE_NAME, 2.0, 1.5));
  env$defaultCorrelationGammaRate <-
    expression(ifelse(shape > 1,
                      (shape - 1) / ifelse(posteriorScale == SD_SCALE_NAME,
                                           10^1, 10^2),
                      shape / ifelse(posteriorScale == SD_SCALE_NAME,
                                     10^1, 10^2)));

  env$defaultCorrelationInverseGammaShape <- 0.5;
  env$defaultCorrelationInverseGammaScale <-
    expression(ifelse(posteriorScale == SD_SCALE_NAME, 10^1, 10^2) *
               (shape + 1));
  env$defaultCorrelationWishartDegreesOfFreedom <- expression(factorDimension + 1);
  env$defaultCorrelationWishartScale <-
    expression(diag(10^2 / ifelse(degreesOfFreedom > factorDimension + 1,
                                     degreesOfFreedom - factorDimension - 1,
                                     degreesOfFreedom),
                    factorDimension));

  env$defaultCorrelationInverseWishartDegreesOfFreedom <- expression(factorDimension - 0.5);
  env$defaultCorrelationInverseWishartInverseScale <-
    expression(diag(10^-2 / (degreesOfFreedom + factorDimension + 1), factorDimension));



  # direct priors
  env$defaultDirectPosteriorScale <- SD_SCALE_NAME;
  env$defaultDirectDataScale <- FREE_SCALE_NAME;
  
  env$defaultDirectGammaShape <-
    expression(ifelse(posteriorScale == SD_SCALE_NAME, 2.0, 1.5));
  env$defaultDirectGammaRate <-
    expression(ifelse(shape > 1,
                      (shape - 1) / ifelse(posteriorScale == SD_SCALE_NAME,
                                           10^2, 10^4),
                      shape / ifelse(posteriorScale == SD_SCALE_NAME,
                                     10^2, 10^4)));
     
  env$defaultDirectInverseGammaShape <- 0.5;
  env$defaultDirectInverseGammaScale <-
    expression(ifelse(posteriorScale == SD_SCALE_NAME, 10^2, 10^4) *
        (shape + 1));

  env$defaultDirectWishartDegreesOfFreedom <- expression(factorDimension + 1);
  env$defaultDirectWishartScale <-
    expression(diag(10^4 / ifelse(degreesOfFreedom > factorDimension + 1,
                                       degreesOfFreedom - factorDimension - 1,
                                       degreesOfFreedom),
                    factorDimension));

  env$defaultDirectInverseWishartDegreesOfFreedom <- expression(factorDimension - 0.5);
  env$defaultDirectInverseWishartInverseScale <-
    expression(diag(10^4 / (degreesOfFreedom + factorDimension + 1), factorDimension));
}
