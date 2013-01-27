loadUnmodeledCoefficientDefaults <- function(env) {
  env$defaultUnmodeledCoefficientPosteriorScale <- VAR_SCALE_NAME; # not used, just need something valid and consistent
  
  env$defaultUnmodeledCoefficientNormalCommonScale <- COMMON_SCALE_TRUE_NAME;
  env$defaultUnmodeledCoefficientNormalSD <- c(10, 2.5);
  env$defaultUnmodeledCoefficientFamily <- NORMAL_FAMILY_NAME;
  env$defaultUnmodeledCoefficientDataScale <- FREE_SCALE_NAME;
}
