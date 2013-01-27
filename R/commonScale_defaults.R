loadCommonScaleDefaults <- function(env) {
  env$defaultCommonScaleCommonScale <- COMMON_SCALE_TRUE_NAME ; # not used, just need something valid and consistent
  
  env$defaultCommonScalePointPosteriorScale <- SD_SCALE_NAME;
  env$defaultCommonScalePointPriorValue <- 1.0;

  env$defaultCommonScaleInverseGammaShape <- 0;
  env$defaultCommonScaleInverseGammaScale <- 0;
  env$defaultCommonScaleInverseGammaPosteriorScale <- VAR_SCALE_NAME;

  env$defaultCommonScaleGammaShape <- 0;
  env$defaultCommonScaleGammaRate <- 0;
  env$defaultCommonScaleGammaPosteriorScale <- VAR_SCALE_NAME;
}
