cat("\n\nRUnit test cases for blme::bmer_lmmTest function\n\n");

test.blme.lmm.internals <- function()
{
  .Call("bmer_lmmTest");
}
