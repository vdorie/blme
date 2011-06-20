
trimWhitespaceRegularExpression <- "^\\s*((?:\\S+\\s*)*\\S+)\\s*$";
# removes whitespace from either end of a string
trim <- function(string) {
  if (!typeof(string) == "character") return (string);

  if (length(string) > 1) {
    return(sapply(1:length(string),
                  function(i) sub(trimWhitespaceRegularExpression, "\\1",
                                  string[i], perl=TRUE)));
  }
  
  return(sub(trimWhitespaceRegularExpression, "\\1", string, perl=TRUE));
}


# in the highly, highly unlikely circumstance that the user feeds us a \1, this will bork.
# a nicer implementation has a regex function that returns a vector of the matches
subSplit <- function(pattern, replacement, x, ignore.case = FALSE, perl = FALSE, 
    fixed = FALSE, useBytes = FALSE, splitChar = '\1')
{
  substitutedString <- sub(pattern, replacement, x, perl=TRUE);
  
  return (strsplit(substitutedString, splitChar)[[1]]);
}

startsWith <- function(string, prefix) {
  return (grepl(paste("^", prefix, sep=""), string, perl=TRUE));
}

endsWith <- function(string, suffix) {
  return (grepl(paste(suffix, "$", sep=""), string, perl=TRUE));
}

contains <- function(string, substring) {
  return (grepl(substring, string, perl=TRUE));
}

flattenStrings <- function(stringVector, sep = " ") {
  if (typeof(stringVector) != "character" || length(stringVector)  < 2) return(stringVector);

  result <- sprintf("%s", stringVector[1]);
  for (i in 2:length(stringVector)) {
    result <- sprintf("%s%s%s", result, sep, stringVector[i]);
  }
  return (result);
}

# splits on commas, but respects parentheses
splitListOnComma <- function(string)
{
  commaIndex <- .Call(bmer_findCommaInList, string);
  if (!is.na(commaIndex)) {
    head <- substr(string, 1, commaIndex - 1);
    tail <- substr(string, commaIndex + 1, nchar(string));
    return (c(head, tail));
  } else {
    return (string);
  }
}
  
splitListOnCommas <- function(string)
{
  splitResult <- splitListOnComma(string);
  if (length(splitResult) == 1) return (list(string));

  return (append(list(splitResult[1]), splitListOnCommas(splitResult[2])));
}
