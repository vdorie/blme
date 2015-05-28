setClass("bmerPointDist",
         representation(value = "numeric"),
         contains = "bmerDist")

toString.bmerPointDist <- function(x, digits = getOption("digits"), ...)
  paste("point(value = ", round(x@value, digits), ")", sep = "")
