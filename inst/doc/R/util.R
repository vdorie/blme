functionToString <- function(name, fn) {
  maxColWidth <- 75;

  indent <- nchar(name) + 1;
  indentString <- "";
  if (indent > 0) {
    for (i in 1:indent) indentString <- sprintf("%s%s", indentString, " ")
  }

  result <- paste(name, "(", sep = "");
  
  values <- formals(fn)
  names <- names(values)

  if (is.symbol(values[[1]])) {
    result <- paste(result, names[[1]], sep = "");
  } else {
    value <- values[[1]];
    if (is.character(value)) value <- paste("\"", value, "\"", sep = "")
    if (is.language(value)) value <- deparse(value)
    result <- paste(result, names[[1]], " = ", value, sep = "")
  }
  
  if (length(values) > 1) for (i in 2:length(values)) {
    if (is.symbol(values[[i]])) {
      result <- paste(result, ", ", names[[i]], sep = "")
    } else {
      value <- values[[i]];
      if (is.character(value)) value <- paste('"', value, '"', sep = "")
      if (is.language(value)) value <- deparse(value)
      result <- paste(result, ", ", names[[i]], " = ", value, sep = "")
    }
  }

  result <- commaColumnWrap(result, maxColWidth);

  if (length(result) == 2) {
    cat("\\verb?", result[1], "?\\\\\n", sep = "");
    cat("\\verb?", indentString, result[2], ")?\n", sep = "");
        ##result <- paste(head, ",\\\\\\\\\n", indentString, tail, sep = "");
  } else {
    cat("\\verb?", result, ")?\n", sep = "");
  }
}


commaColumnWrap <- function(string, maxColWidth) {
  if (nchar(string) <= maxColWidth) return(string);

  ## replaces last ", " by a question mark and then splits on that
  ## OK, because question marks aren't allowed aside from help queries (I think)
  split <- strsplit(sub("(.*), (.*)", "\\1?\\2", string, perl = TRUE), "\\?")
  head <- split[[1]][1];
  tail <- split[[1]][2];

  while (nchar(head) > maxColWidth) {
    split <- strsplit(sub("(.*), (.*)", "\\1?\\2", head, perl = TRUE), "\\?")
    head <- split[[1]][1];
    tail <- paste(split[[1]][2], ", ", tail, sep = "");
  }
  
  return(c(paste(head, ",", sep = ""), tail));
}

makeOrCheckDir <- function(dirName)
{
  if (!file.exists(dirName)) dir.create(dirName)
  info <- file.info(dirName)
  if (!info$isdir) stop("'", dirName, "' is not a directory")
}
