setClass("bmerMod",
         representation(priors = "list"), contains = "merMod");
setClass("blmerMod", contains=c("lmerMod", "bmerMod"));
setClass("bglmerMod", contains=c("glmerMod", "bmerMod"));
