cat("\n\nRUnit test cases for blme:::bmer_matrixTest function\n\n");

test.blme.matrix.internals <- function()
{
  .Call("bmer_matrixTest");
}
