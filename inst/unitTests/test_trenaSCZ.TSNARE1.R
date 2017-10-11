library(RUnit)
library(trenaSCZ.TSNARE1)
#------------------------------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#------------------------------------------------------------------------------------------------------------------------
#PORT.RANGE <- 8000:8020
#if(!exists("tv"))
#   tv <- trenaViz(PORT.RANGE, quiet=TRUE);
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
  test_constructor();
  test_dataSummary();

} # runTests
#------------------------------------------------------------------------------------------------------------------------
test_constructor <- function()
{
   printf("--- test_constructor")
   tscz <- trenaSCZ(gwasLocusNumber=5)
   checkTrue("trenaSCZ" %in% is(tscz))

} # test_constructor
#------------------------------------------------------------------------------------------------------------------------
test_dataSummary <- function()
{
   printf("--- test_dataSummary")
   tscz <- trenaSCZ(gwasLocusNumber=5)
   dataSummary(tscz)

} # test_dataSummary
#------------------------------------------------------------------------------------------------------------------------
