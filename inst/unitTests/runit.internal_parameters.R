cat("\n\nRUnit test cases for blme:::bmer_parametersTest function\n\n");

test.blme.parameters.internals <- function()
{
  .Call("bmer_parametersTest");
}
