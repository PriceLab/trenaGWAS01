.trenaSCZ <- setClass("trenaSCZ",
                      representation=representation(
                         dataSummary="character",
                         gwasLocusNumber="numeric",
                         trenaViz="trenaViz"
                         )
                      )
#------------------------------------------------------------------------------------------------------------------------
setGeneric("dataSummary",    signature="obj", function(obj) standardGeneric("dataSummary"))
setGeneric("setRegion",      signature="obj", function(obj, chrom, start, end) standardGeneric("setRegion"))
setGeneric("getRegion",      signature="obj", function(obj) standardGeneric("getRegion"))
#setGeneric("loadFootprints", signature="obj", function(obj,
#------------------------------------------------------------------------------------------------------------------------
trenaViz.PORT.RANGE <- 8000:8020
#------------------------------------------------------------------------------------------------------------------------
trenaSCZ = function(gwasLocusNumber)
{
   tv <- trenaViz(PORT.RANGE, quiet=TRUE)
   setGenome(tv, "hg38")

   obj <- .trenaSCZ(gwasLocusNumber=gwasLocusNumber, trenaViz=tv)

} # trenaSCZ
#------------------------------------------------------------------------------------------------------------------------
setMethod("dataSummary", "trenaSCZ",

     function(obj){
        s <- "---- basic data used in this instance of trenaSCZ";
        s <- paste(s, "", sep="\n")
        s <- paste(s, "tbl.loci: ", sep="\n")
        s <- paste(s, "tbl.credSnp: ", sep="\n")
        s <- paste(s, "tbl.wg.trena: ", sep="\n")
        s <- paste(s, "tbl.geneModel.mtx.cer.TSNARE1.250kb: ", sep="\n")
        cat(s)
        cat("\n")
       })

#------------------------------------------------------------------------------------------------------------------------
setMethod("setRegion", "trenaSCZ",

     function(obj, chrom, start, end){
        chromString <- sprintf("%s:%d-%d", chrom, start, end)
        })

#------------------------------------------------------------------------------------------------------------------------
setMethod("getRegion", "trenaSCZ",

     function(obj){
        parseChromLocString(getGenomicRegion(obj@trenaViz))
        })

#------------------------------------------------------------------------------------------------------------------------

