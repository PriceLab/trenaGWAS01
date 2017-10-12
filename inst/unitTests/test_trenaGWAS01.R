library(RUnit)
library(trenaGWAS01)
#------------------------------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#------------------------------------------------------------------------------------------------------------------------
if(!exists("tgwas"))
   tgwas <- trenaGWAS01(gwasLocusNumber=5, targetGene="TSNARE1", targetGene.tss=142354831, quiet=TRUE)

#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
  test_constructor();
  test_dataSummary();
  test_setGetRegion()
  test_loadFootprintsAndSnps_thenIntersection()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
test_constructor <- function()
{
   printf("--- test_constructor")
   checkTrue("trenaGWAS01" %in% is(tgwas))

} # test_constructor
#------------------------------------------------------------------------------------------------------------------------
test_dataSummary <- function()
{
   printf("--- test_dataSummary")
   dataSummary(tgwas)

} # test_dataSummary
#------------------------------------------------------------------------------------------------------------------------
test_setGetRegion <- function()
{
   printf("--- test_dataSummary")
   setRegion(tgwas, "chr1")
   Sys.sleep(1)
   region.at.start <- getRegion(tgwas)
   checkEquals(region.at.start, list(chrom="chr1", start=1, end=248956422))

   setRegion(tgwas, "TSNARE1")
   Sys.sleep(1)
   region.tsnare1 <- getRegion(tgwas)
   checkEquals(region.tsnare1, list(chrom="chr8", start=142211080, end=142404249))

} # test_setGetRegion
#------------------------------------------------------------------------------------------------------------------------
test_loadFootprintsAndSnps_thenIntersection <- function()
{
   printf("--- test_loadFootprintsAndSnps_thenIntersection")

   setRegion(tgwas, "chr8:142,230,914-142,234,913") # rs49776977 is at 142,232,793

   system.time(
      tbl.fp <- loadFootprints(tgwas)
         )   # about 32 seconds locally to whovian, 192 to bddsrds

   checkEquals(ncol(tbl.fp), 10)
   checkTrue(nrow(tbl.fp) > 50)
   expected.colnames <-  c("chrom",  "motifStart",  "motifEnd",  "motifName", "strand",  "score",  "length",
                           "distance.from.tss","id", "fpSource")
   checkTrue(all(expected.colnames %in% colnames(tbl.fp)))

   system.time(
      tbl.snps <- loadSNPs(tgwas)
      )   # 85 seconds with SNPlocs package, 0.3 seconds with extdata/tbl.imputedSnps.chr8.RData

   tbl.snpsInFp <- findSNPsInFootprints(tgwas, tbl.snps, tbl.fp)
   checkEquals(nrow(tbl.snpsInFp), 2)
   checkEquals(sort(tbl.snpsInFp$rsid), c("rs12676511", "rs4976977"))

} # test_loadFootprintsAndSnps_thenIntersection
#------------------------------------------------------------------------------------------------------------------------
test_displayMotif <- function()
{
   printf("--- test_motif")
   setRegion(tgwas, "chr8:142,230,914-142,234,913") # rs49776977 is at 142,232,793


} # test_loadSNPs
#------------------------------------------------------------------------------------------------------------------------
# chr8 142232793      * |   rs4976977  W (A -> T)
demo_NFIX_snp <- function()
{

} # demo_NFIX_snp
#------------------------------------------------------------------------------------------------------------------------
