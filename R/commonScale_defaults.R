loadCommonScaleDefaults <- function(env) {
  env$defaultCommonScalePointPosteriorScale <- SD_SCALE_NAME;
  env$defaultCommonScalePointPriorValue <- 1.0;

  env$defaultCommonScaleInverseGammaShape <- 0;
  env$defaultCommonScaleInverseGammaScale <- 0;
  env$defaultCommonScaleInverseGammaPosteriorScale <- VAR_SCALE_NAME;
}
