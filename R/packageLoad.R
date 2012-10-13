# Most of the internal names are defined in the C portion of the package.
# Rather than enter them twice and risk them getting out of sync,
# we use the package load hook to pull them into the namespace.
.onLoad <- function(libname, pkgname)
{
  namespace <- getNamespace(pkgname);
  if (is.null(namespace)) stop("cannot find namespace environment, cannot load constants from C");
  
  loadConstants(namespace);
  loadCovarianceDefaults(namespace);
  loadUnmodeledCoefficientDefaults(namespace);
  loadCommonScaleDefaults(namespace);
}
