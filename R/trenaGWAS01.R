library(SNPlocs.Hsapiens.dbSNP150.GRCh38)

.trenaGWAS01 <- setClass("trenaGWAS01",
                      representation=representation(
                         dataSummary="character",
                         gwasLocusNumber="numeric",
                         targetGene="character",
                         targetGene.tss="numeric",
                         trena="Trena",
                         trenaViz="trenaViz",
                         quiet="logical"
                         )
                      )
#------------------------------------------------------------------------------------------------------------------------
setGeneric("dataSummary",     signature="obj", function(obj) standardGeneric("dataSummary"))
setGeneric("setRegion",       signature="obj", function(obj, targetString) standardGeneric("setRegion"))
setGeneric("getRegion",       signature="obj", function(obj) standardGeneric("getRegion"))
setGeneric("loadFootprints",  signature="obj", function(obj) standardGeneric("loadFootprints"))

setGeneric("loadSNPs",        signature="obj", function(obj) standardGeneric("loadSNPs"))
setGeneric("findSNPsInFootprints", signature="obj", function(obj, tbl.snps, tbl.fp) standardGeneric("findSNPsInFootprints"))
setGeneric("displayMotif",    signature="obj", function(obj, motifName) standardGeneric("displayMotif"))

setGeneric("assessSnpsInContextOfGeneModel", signature="obj", function(obj, tbl.snps, geneModel, targetGene,
                                                                       matchThreshold, shoulder)
              standardGeneric("assessSnpsInContextOfGeneModel"))
#------------------------------------------------------------------------------------------------------------------------
trenaViz.PORT.RANGE <- 8000:8020
database.host <- "bddsrds.globusgenomics.org"
#database.host <- "whovian"
#------------------------------------------------------------------------------------------------------------------------
trenaGWAS01 = function(gwasLocusNumber, targetGene, targetGene.tss, quiet=TRUE)
{
   local.trena <- Trena("hg38", quiet)

   local.tv <- trenaViz(trenaViz.PORT.RANGE, quiet=quiet)
   setGenome(local.tv, "hg38")

   obj <- .trenaGWAS01(gwasLocusNumber=gwasLocusNumber,
                       targetGene=targetGene,
                       targetGene.tss=targetGene.tss,
                       trena=local.trena,
                       trenaViz=local.tv,
                       quiet=quiet)

      # user will probably want direct access to these objects also
   printf("  also creating user-level (global scope) 'trena' and 'tv' (trenaViz) objects....")
   trena <<- local.trena
   tv <<- local.tv

   obj

} # trenaGWAS01
#------------------------------------------------------------------------------------------------------------------------
setMethod("dataSummary", "trenaGWAS01",

     function(obj){
        s <- "---- basic data used in this instance of trenaGWAS01";
        s <- paste(s, "", sep="\n")
        s <- paste(s, "tbl.loci: ", sep="\n")
        s <- paste(s, "tbl.credSnp: ", sep="\n")
        s <- paste(s, "tbl.wg.trena: ", sep="\n")
        s <- paste(s, "tbl.geneModel.mtx.cer.TSNARE1.250kb: ", sep="\n")
        cat(s)
        cat("\n")
       })

#------------------------------------------------------------------------------------------------------------------------
setMethod("setRegion", "trenaGWAS01",

     function(obj, targetString){
        showGenomicRegion(obj@trenaViz, targetString)
        })

#------------------------------------------------------------------------------------------------------------------------
setMethod("getRegion", "trenaGWAS01",

     function(obj){
        parseChromLocString(getGenomicRegion(obj@trenaViz))
        })

#------------------------------------------------------------------------------------------------------------------------
# query trenaViz (in the browser) for the currently displayed genomic region
# then use the trena convenience method "getRegulatoryChromosomalRegions" to get and
# display brain DHS footprints reported in that region
#
# targetGene and targetGene.tss are used so that the resulting footprints can be
# annotated with respect to the (presumed) gene of interest
#
setMethod("loadFootprints", "trenaGWAS01",

    function(obj){
       dbNames <- c("brain_hint_16", "brain_hint_20", "brain_wellington_20", "brain_wellington_16")
       brain.hint16.db.uri <- sprintf("postgres://%s/%s", database.host, dbNames[1])
       brain.hint20.db.uri <- sprintf("postgres://%s/%s", database.host, dbNames[2])
       brain.wellington16.db.uri <- sprintf("postgres://%s/%s", database.host, dbNames[3])
       brain.wellington20.db.uri <- sprintf("postgres://%s/%s", database.host, dbNames[4])
       sources <- list(brain.hint16.db.uri, brain.hint20.db.uri,
                       brain.wellington16.db.uri,brain.wellington20.db.uri)

       names(sources) <- dbNames
       current.region <- parseChromLocString(getGenomicRegion(obj@trenaViz))
       printf("span: %d", 1 + current.region$end - current.region$start)
       tbls <- getRegulatoryChromosomalRegions(obj@trena,
                                               current.region$chrom, current.region$start, current.region$end,
                                               sources, obj@targetGene, obj@targetGene.tss)
       names(tbls) <- dbNames
       track.count <- length(tbls)

       #colors <- rainbow_hcl(track.count)
       #colors <- terrain_hcl(track.count)
       colors <- heat_hcl(track.count)
       #colors <- sequential_hcl(track.count)

       next.color <- 1

       for(i in seq_len(length(tbls))){
          track.name <- names(sources)[i]
          # disabled until this gets fixed in igv.js (pshannon (12 oct 2017))
          #if(track.name %in% getTrackNames(obj@trenaViz))
          #   removeTracksByName(obj@trenaViz, track.name)
          columns.of.interest <- c("chrom", "motifStart", "motifEnd", "motifName")

          tbl <- tbls[[i]]
          tbls[[i]]$fpSource <- dbNames[i]   # in preparation for combining them

          stopifnot(all(columns.of.interest %in% colnames(tbl)))
          addBedTrackFromDataFrame(obj@trenaViz, track.name, tbl[, columns.of.interest], color=colors[next.color]);
          next.color <- next.color + 1
          } # for i

       tbl.single <- do.call(rbind, tbls)
       rownames(tbl.single) <- NULL
       tbl.single
       }) # addFootprintsInCurrentRegion

#----------------------------------------------------------------------------------------------------
# both the index snp, and its accompanying imputed ("credible") snps
setMethod("loadSNPs", "trenaGWAS01",

     function(obj){

        load(system.file(package="trenaGWAS01", "extdata", "tbl.loci.RData"))
        load(system.file(package="trenaGWAS01", "extdata", "tbl.credSnp.RData"))
        load(system.file(package="trenaGWAS01", "extdata", "tbl.imputedSnps.chr8.RData"))

           # get the index snps, and all the credible snps, in the current locus
           # note that this may be different from the currently display region in trenaViz
        locus.chrom <- tbl.loci$chrom[obj@gwasLocusNumber]
        locus.start <- tbl.loci$start38[obj@gwasLocusNumber]
        locus.end   <- tbl.loci$end38[obj@gwasLocusNumber]
        locusLoc <- sprintf("%s:%d-%d", locus.chrom, locus.start, locus.end)

        tbl.lociSnps <- subset(tbl.credSnp, chrom==locus.chrom & pos38 >= locus.start & pos38 <= locus.end)

        with(unique(tbl.lociSnps[, c("id", "chrom", "pos38")]),
             addBedTrackFromDataFrame(obj@trenaViz, "indexSNP",
                                      data.frame(chrom=chrom, start=pos38, end=pos38, stringsAsFactors=FALSE), color="red"))

             # now get their hg38 coordinates.  for now, ignore all imputed snps without rsid
        imputed.snps <- grep("^rs", tbl.lociSnps$credibleSNP, value=TRUE)
           # SNPlocs is way too slow.  so for this case study, chr8 snps of interest have been extracted
           # from that bioc data packages, saved in extdata, loaded above
           # printf("retrieveing snp locations from Bioconductor package 'SNPlocs.Hsapiens.dbSNP150.GRch38'")
           # printf("takes a few minutes on the first call, and ~30 seconds on subsequent calls...")
           # tbl.imputed.snps <- as.data.frame(snpsById(SNPlocs.Hsapiens.dbSNP150.GRCh38, imputed.snps))

        tbl.imputed.snps <- subset(tbl.imputedSnps.chr8, RefSNP_id %in% tbl.credSnp$credibleSNP)
              # transform tbl.imputed.snps a bit for rendering in trenaViz: rearrange & rename columns
        tbl.imputedSnp.bed <- tbl.imputed.snps[, c("seqnames", "pos", "pos", "RefSNP_id")]
        colnames(tbl.imputedSnp.bed) <- c("chrom", "start", "end", "rsid")
        tbl.imputedSnp.bed$chrom <- as.character(tbl.imputedSnp.bed$chrom)
        tbl.imputedSnp.bed$chrom <- sprintf("chr%s", tbl.imputedSnp.bed$chrom)
        current.region <- parseChromLocString(getGenomicRegion(obj@trenaViz))
        tbl.imputedSnp.bed <- subset(tbl.imputedSnp.bed, start >= current.region$start &
                                                         end   <= current.region$end &
                                                         chrom == current.region$chrom)
        printf("found %d imputed snps in region", nrow(tbl.imputedSnp.bed))
        addBedTrackFromDataFrame(obj@trenaViz, "imputed", tbl.imputedSnp.bed, color="purple")
        return(tbl.imputedSnp.bed)
        }) # loadSNPs

#----------------------------------------------------------------------------------------------------
setMethod("findSNPsInFootprints", "trenaGWAS01",

     function(obj, tbl.snps, tbl.fp){
        tbl.fp.simple <- unique(tbl.fp[, 1:4])
        colnames(tbl.fp.simple) <- c("chrom", "start", "end", "motifName")
        gr.fp <- GRanges(tbl.fp.simple)
        gr.snp <- GRanges(tbl.snps)
        tbl.ov <- as.data.frame(findOverlaps(gr.snp, gr.fp))
        if(nrow(tbl.ov) == 0){
           printf("no snps in footprints")
           return(data.frame())
           }
        colnames(tbl.ov) <- c("snp", "fp")
        tbl.snpFp <- tbl.snps[unique(tbl.ov$snp),]
        rownames(tbl.snpFp) <- NULL
        addBedTrackFromDataFrame(obj@trenaViz, "imputedInFP", tbl.snpFp, color="magenta")
        tbl.snpFp
        })

#----------------------------------------------------------------------------------------------------
setMethod("assessSnpsInContextOfGeneModel", "trenaGWAS01",

     function(obj, tbl.snps, geneModel, targetGene, matchThreshold, shoulder){

        stopifnot(all(c("chrom", "start", "end", "rsid") %in% colnames(tbl.snps)))

        jaspar2016.pfms <- query(MotifDb, "jaspar2016")
        human.pfms <- query(jaspar2016.pfms, "sapiens")
        mouse.pfms <- query(jaspar2016.pfms, "mus")
        rat.pfms   <- query(jaspar2016.pfms, "rnorvegicus")
        pfms <- as.list(c(human.pfms, mouse.pfms, rat.pfms))   # 619

        motifMatcher <- MotifMatcher("hg38", pfms)

        for(r in seq_len(nrow(tbl.snps))){
          rsid <- tbl.snps$rsid[r]
          printf("assessing %s", rsid)
          tbl.region <- with(tbl.snps[r,],
                             data.frame(chrom=chrom, start=start-shoulder, end=end+shoulder, stringsAsFactors=FALSE))

          tbl.wt  <- findMatchesByChromosomalRegion(motifMatcher, tbl.region, matchThreshold)
          tbl.mut <- findMatchesByChromosomalRegion(motifMatcher, tbl.region, matchThreshold, variants=rsid)

          tfs.wt <- c()
          tfs.mut <- c()

          if(nrow(tbl.wt) > 0){
             tbl.wt$shortMotif <- unlist(lapply(tbl.wt$motifName,
                                                function(s) {tokens <- strsplit(s, "-")[[1]]; return(tokens[length(tokens)])}))

             tfs.wt <- c(tfs.wt, associateTranscriptionFactors(MotifDb, tbl.wt, source="TFClass", expand.rows=TRUE)$geneSymbol)
             tfs.wt <- unique(c(tfs.wt, associateTranscriptionFactors(MotifDb, tbl.wt, source="MotifDb", expand.rows=TRUE)$geneSymbol))
             if(any(is.na(tfs.wt)))
                tfs.wt <- tfs.wt[-which(is.na(tfs.wt))]
             } # nrow(tbl.wt) > 0

          if(nrow(tbl.mut) > 0){
             tbl.mut$shortMotif <- unlist(lapply(tbl.mut$motifName,
                                                 function(s) {tokens <- strsplit(s, "-")[[1]]; return(tokens[length(tokens)])}))
             tfs.mut <- c(tfs.mut, associateTranscriptionFactors(MotifDb, tbl.mut, source="TFClass", expand.rows=TRUE)$geneSymbol)
             tfs.mut <- unique(c(tfs.mut, associateTranscriptionFactors(MotifDb, tbl.mut, source="MotifDb", expand.rows=TRUE)$geneSymbol))
             if(any(is.na(tfs.mut)))
                tfs.mut <- tfs.mut[-which(is.na(tfs.mut))]
             } # nrow(tbl.mut) > 0


          model.tfs <- geneModel$gene

          tfs.wt <- unique(tfs.wt)
          tfs.wt.inModel <- intersect(tfs.wt, model.tfs)

          tfs.mut <- unique(tfs.mut)
          tfs.mut.inModel <- intersect(tfs.mut, model.tfs)

          tfs.lost <- setdiff(tfs.wt, tfs.mut)
          tfs.gained <- setdiff(tfs.mut, tfs.wt)

          tfs.inModel.lost   <- intersect(tfs.lost,   model.tfs)
          tfs.inModel.gained <- intersect(tfs.gained, model.tfs)

          printf("%s %s, tfs.wt,     %2d, in model: %2d", targetGene, rsid, length(tfs.wt),     length(tfs.wt.inModel))
          printf("%s %s, tfs.mut,    %2d, in model: %2d", targetGene, rsid, length(tfs.mut),    length(tfs.mut.inModel))
          printf("%s %s, tfs.lost,   %2d, in model: %2d", targetGene, rsid, length(tfs.lost),   length(tfs.inModel.lost))
          printf("%s %s, tfs.gained, %2d, in model: %2d", targetGene, rsid, length(tfs.gained), length(tfs.inModel.gained))

          if(length(tfs.inModel.lost) > 0){
             tf.model.indices <- match(tfs.inModel.lost, model.tfs)
             if(any(is.na(tf.model.indices)))
                tf.model.indices <- tf.model.indices[-which(is.na(tf.model.indices))]
             for(index in tf.model.indices)
                printf("  tf %s rank %2d lost in model", model.tfs[index], index)
             } # if tfs.inModel.lost

          if(length(tfs.inModel.gained) > 0){
             tf.model.indices <- match(tfs.inModel.gained, model.tfs)
             if(any(is.na(tf.model.indices)))
                tf.model.indices <- tf.model.indices[-which(is.na(tf.model.indices))]
             for(index in tf.model.indices)
                printf("  tf %s rank %2d gained in wg model", model.tfs[index], index)
            } # if tfs.inModel.gained
          } # for r
     }) # assessSnpsInContextOfGeneModel

#----------------------------------------------------------------------------------------------------
setMethod("displayMotif", "trenaGWAS01",

     function(obj, motifName){
        browser()
        xyz <- 99
        })

#----------------------------------------------------------------------------------------------------
# collapse all footprints from multiple sources into one non-redundant, non-overlapping track
# this may be a useful visualization, as well as the basis for an exhaustive search for binding
# motifs, using trena's MotifMatcher class
#
# the "expandFootprintsBy" parameter will lead to the unification of more individual footprints
# this will sometimes be useful:  calling footprints is an inexact art; allowing for extra
# width may help in generating hypotheses - implicating snps - which would otherwise be missed.
createCollapsedFootprintTrack <- function(expandFootprintsBy=10)
{
   dbNames <- c("brain_hint_16", "brain_hint_20", "brain_wellington_20", "brain_wellington_16")
   brain.hint16.db.uri <- sprintf("postgres://%s/%s", database.host, dbNames[1])
   brain.hint20.db.uri <- sprintf("postgres://%s/%s", database.host, dbNames[2])
   brain.wellington16.db.uri <- sprintf("postgres://%s/%s", database.host, dbNames[3])
   brain.wellington20.db.uri <- sprintf("postgres://%s/%s", database.host, dbNames[4])

   sources <- list(brain.hint16.db.uri, brain.hint20.db.uri, brain.wellington16.db.uri, brain.wellington20.db.uri)
   names(sources) <- dbNames

   current.region <- parseChromLocString(getGenomicRegion(tv))

   printf("span: %d", 1 + current.region$end - current.region$start)
   tbls <- getRegulatoryChromosomalRegions(obj@trena,
                                           current.region$chrom, current.region$start, current.region$end,
                                           sources, targetGene, targetGene.tss)
     # make a GRangesList out of the 4 GRanges, converting each data.frame to a GRange object
   grl <- GRangesList(lapply(tbls, GRanges))
   gr.all <- Reduce(c, grl)
   strand(gr.all) <- "*"  # ignore strand differences
   start(gr.all) <- start(gr.all) - expandFootprintsBy
   end(gr.all) <- end(gr.all) + expandFootprintsBy

   tbl.fpCollapsed <- as.data.frame(reduce(gr.all))  # metadata discarded
   trackName <- "fpCollapsed"
   if(trackName %in% getTrackNames(tv))
      removeTracksByName(tv, trackName)
   addBedTrackFromDataFrame(tv, trackName, tbl.fpCollapsed, color="green")

     # find the imputed snp, collapsed fp overlaps
   #colnames(tmp.bed)[1] <- "seqnames"
   #tmp.bed$seqnames <- sprintf("chr%s", tmp.bed$seqnames)

   #tmp.bed$seqnames <- "chr8"
   #tbl.overlaps <- as.data.frame(findOverlaps(GRanges(tmp.bed), GRanges(tbl.fpCollapsed)))
   #colnames(tbl.overlaps) <- c("snp", "fp")
   #trackName <- "snpInFp"
   #if(trackName %in% getTrackNames(tv))
   #   removeTracksByName(tv, trackName)

   #tbl.snpsInFP <<- tmp.bed[tbl.overlaps$snp,]
   #addBedTrackFromDataFrame(tv, trackName, tmp.bed[tbl.overlaps$snp,], color="magenta")

   invisible(tbl.fpCollapsed)

} # createCollapsedFootprintTrack
#----------------------------------------------------------------------------------------------------

