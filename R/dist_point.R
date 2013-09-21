setClass("bmerPointDist", contains = "bmerDist",
         slots = c(value = "numeric"));

toString.bmerPointDist <- function(x, digits = getOption("digits"), ...)
  paste("point(value = ", round(x@value, digits), ")", sep = "");

