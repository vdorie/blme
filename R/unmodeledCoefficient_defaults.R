loadUnmodeledCoefficientDefaults <- function(env) {
  env$defaultUnmodeledCoefficientNormalCovarianceScale <- ABSOLUTE_SCALE_NAME;
  env$defaultUnmodeledCoefficientNormalSD <- c(10, 2.5);
  env$defaultUnmodeledCoefficientFamily <- NORMAL_FAMILY_NAME;
  env$defaultUnmodeledCoefficientDataScale <- FREE_SCALE_NAME;
}
