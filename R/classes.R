# "AllClass" name is to bump up this file in load order
setClass("bmerPrior",
         representation = list(type            = "integer",
                               families        = "integer",
                               scales          = "integer",
                               hyperparameters = "numeric"),
         validity = function(object) TRUE);

setClass("bmer",
	 representation(cov.prior   = "list", # list of bmerPriors, one per grouping factor
                        fixef.prior = "bmerPrior",
                        var.prior   = "bmerPrior"),
         contains = "mer",
         validity = function(object) .Call(mer_validate, object))

setClass("summary.bmer",
         contains = c("bmer", "summary.mer"))

# Constructors, more or less. Not necessary, but should help with
# cleanliness.
createFlatPriorObject <- function() {
  return(new(Class           = "bmerPrior",
             type            = getEnumOrder(typeEnumeration, NONE_TYPE_NAME),
             families        = integer(0),
             scales          = integer(0),
             hyperparameters = double(0)));
}

createPriorObject <- function(type, fields) {
  return(new(Class           = "bmerPrior",
             type            = getEnumOrder(typeEnumeration, type),
             families        = fields$families,
             scales          = fields$scales,
             hyperparameters = fields$hyperparameters));
}
