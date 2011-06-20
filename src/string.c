#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <Rdefines.h>
#include <R_ext/RS.h>  /* for Calloc/Free */
#include <wchar.h>
// #include <wctype.h>    /* for wctrans_t */

#include <R_ext/libextern.h>
LibExtern Rboolean mbcslocale;

// since in R it is bloody hard to get at a string character by character, and even if we do so it is
// really slow, here is an alternative
SEXP bmer_findCommaInList(SEXP input)
{
  if (!isString(input)) error("non-character argument");
  
  int inputLength = LENGTH(input);
  
  if (inputLength == 0) return(R_NilValue);
  
  Rboolean useUTF8 = FALSE;
  const char *inputCharacters = NULL;
  
  SEXP resultExpression;
  PROTECT(resultExpression = allocVector(INTSXP, inputLength));
  int *result = INTEGER(resultExpression);
  
  int i;
  for (int inputStringNum = 0; inputStringNum < inputLength; ++inputStringNum) {
    if (getCharCE(STRING_ELT(input, inputStringNum)) == CE_UTF8) {
      useUTF8 = TRUE;
      break;
    } else {
      useUTF8 = FALSE;
    }

    if (useUTF8) {
      inputCharacters = translateCharUTF8(STRING_ELT(input, inputStringNum));
    } else {
      inputCharacters = translateChar(STRING_ELT(input, inputStringNum));
    }
    
    if (inputCharacters == NULL) {
      result[inputStringNum] = NA_INTEGER;
      continue;
    }
    
    int numParentheses = 0;
    for (i = 0; inputCharacters[i] != '\0'; ++i) {
      switch (inputCharacters[i]) {
        case '(':
          ++numParentheses;
          break;
        case ')':
          --numParentheses;
          break;
        case ',':
          if (numParentheses == 0) goto BMER_SPLIT_ON_COMMAS_LOOP_DONE;
          break;
        default:
          break;
      }
    }
BMER_SPLIT_ON_COMMAS_LOOP_DONE:
    if (inputCharacters != NULL && inputCharacters[i] != ',') {
      result[inputStringNum] = NA_INTEGER;
      continue;
    }
    
    int commaIndex = i;
    
    result[inputStringNum] = commaIndex + 1;
  }
  
  UNPROTECT(1);
  
  return(resultExpression);
}
