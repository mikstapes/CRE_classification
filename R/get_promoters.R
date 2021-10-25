#! /home/phan/miniconda3/envs/r_env/bin/Rscript
library(here)
source(here("R", "functions.R"))

# get proximal proms based on TSS, chrom profile and RNA read counts


args <- commandArgs(trailingOnly=T)
if (length(args) != 6) stop('Usage: ./get_promoters.R <ref genome> <prom_state.bed> <atac_peaks.bed> <RNA bam dir> <outdir> <sample_name>')

bioc_pkg <- c("GenomicFeatures", "rtracklayer", "GenomicAlignments", "Rsamtools")
tidy_pkg <- c("ggplot2", "tidyr", "dplyr", "readr")
genomes <- c("mm10", "hg38")

suppressMessages(lapply(bioc_pkg, library, character.only = TRUE)) 
suppressMessages(lapply(tidy_pkg, library, character.only = TRUE))

##---Parsing inputs

ref_genome <- args[1]

prom_state <- import.bed(args[2])

atac <- import.bed(args[3])

rna_bam_path <- args[4]

out_path <- args[5]

sample_id <- args[6]


##---- Load txdb, input data & filter only standard chroms (autosomes + chrX or chrZ(gg))

if (ref_genome %in% genomes) {
  if (ref_genome == "mm10") pkg <- "TxDb.Mmusculus.UCSC.mm10.knownGene"
  if (ref_genome == "hg38") pkg <- "TxDb.Hsapiens.UCSC.hg38.knownGene"
  suppressMessages(require(pkg, character.only = TRUE))
  assign("txdb", eval(parse(text = pkg)))
  atac <- getStandardChrom(atac, ref = ref_genome)
} else if (ref_genome=="galGal6") {
  txdb <- loadDb("/project/MDL_Ibrahim/MP_all/annotations/gg6/TxDb.galGal6.ncbiRefseq.sqlite")
  atac <- getStandardChrom(atac, ref = ref_genome)
} else {
  stop('genome not supported')
}

######--- Get Tx with active promoters ---######

#- Filter 1: unique TSS

Tx <- GenomicFeatures::transcripts(txdb)
  Tx <- trim(Tx[!duplicated(start(Tx))])
  Tx <- getStandardChrom(Tx, ref = ref_genome)

#- Filter 2: take overlap of TSS and prom-like state

Tx_with_active_prom <- Tx[findOverlaps(resize(Tx, width = 1, fix = "start"), 
                                       prom_state, maxgap = 0)@from]
  

#- Filter 3: Get unique & (longest) isoform of each Tx for counting

#get length in kb
Tx_with_active_prom$length <- width(Tx_with_active_prom)/1e3

#filter by length and overlap with other isoform

Tx_with_active_prom_filt <- getUniqueTx(Tx_with_active_prom)


######--- Load RNAseq BAM file from both reps and get counts ---######

bam_path <- list.files(path = rna_bam_path, 
                       pattern = ".bam$", full.names = TRUE) %>% 
  setNames(gsub(".bam$", "", basename(.)))

if (ref_genome == "galGal6") {
  
  bam_files <- BamFileList(bam_path, asMates = TRUE, yieldSize = 5e6)
  tss_se <- summarizeOverlaps(features = Tx_with_active_prom_filt, 
                              reads = bam_files, 
                              singleEnd=FALSE, fragments=TRUE, mode = "IntersectionNotEmpty",
                              BPPARAM = BiocParallel::MulticoreParam(20))
} else if (ref_genome != "galGal6") {
  
  bam_files <- BamFileList(bam_path, asMates = FALSE, yieldSize = 2e6)
  tss_se <- summarizeOverlaps(features = Tx_with_active_prom_filt, 
                              reads = bam_files, 
                              singleEnd=TRUE, fragments=FALSE, mode = "IntersectionNotEmpty",
                              BPPARAM = BiocParallel::MulticoreParam(10))
} else { stop('genome not supported')}
  
  
# add Tx names to counts 
  
#rownames(se_bam_counts)  <- rowRanges(se_bam_counts)$tx_name

# add counts & TPM to filtered Tx GRanges 

Tx_with_active_prom_filt$count.1  <- assay(tss_se[,1])[,1]; Tx_with_active_prom_filt$count.1.TPM  <- calculateTPM(tss_se[,1])[,1]
Tx_with_active_prom_filt$count.2  <- assay(tss_se[,2])[,1]; Tx_with_active_prom_filt$count.2.TPM  <- calculateTPM(tss_se[,2])[,1]

# add average TPM and save as RData

Tx_with_active_prom_filt$avg.TPM <- (Tx_with_active_prom_filt$count.1.TPM + Tx_with_active_prom_filt$count.2.TPM)/2
    
# Filter 4: Keep TSS where Tx has >= 2 TPM

prom_regions <- promoters(Tx_with_active_prom_filt[which(Tx_with_active_prom_filt$avg.TPM >= 2)], 
                       upstream = 1000, downstream = 200)

# save as RData for references
save(prom_regions, file = file.path(out_path, paste0(sample_id, "_promoters_from_annotated_TSS.RData")))

# removes all counts mcols for final prom object
mcols(prom_regions) <- mcols(prom_regions)["tx_name"]
# rename to mcol for consistency
names(mcols(prom_regions)) <- "name"
# add annotation status
prom_regions$status <- "annotated"
# add ranges name
names(prom_regions) <- prom_regions$name

######--- Get unannoted promoters ---######

#- get promoters without filtered annotated TSS

prom_noRef <- prom_state[!overlapsAny(prom_state,
                                      start(Tx_with_active_prom_filt),
                                      maxgap = 0)]

#- get ATAC peaks as possible un-annotated promoters

atac_prom <- atac[overlapsAny(atac, prom_noRef, maxgap = 0)]

#- Get BAM counts from rep

if (ref_genome == "galGal6") {
  
  bam_files <- BamFileList(bam_path, asMates = TRUE, yieldSize = 5e6)
  atac_se <- summarizeOverlaps(features = atac_prom, 
                               reads = bam_files, 
                               singleEnd=FALSE, fragments=TRUE, mode = "Union",
                               BPPARAM = BiocParallel::MulticoreParam(20))
  
  # get lib size (PE)
  bam.param  <- ScanBamParam(flag = scanBamFlag(isPaired  = T, isDuplicate = F))
  bam_counts  <- countBam(bam_files, param = bam.param)
  bam_counts  <- bam_counts$records/2
  
} else if (ref_genome != "galGal6") {
  
  bam_files <- BamFileList(bam_path, asMates = FALSE, yieldSize = 2e6)
  atac_se <- summarizeOverlaps(features = atac_prom, 
                               reads = bam_files, 
                               singleEnd=TRUE, fragments=FALSE, mode = "Union",
                               BPPARAM = BiocParallel::MulticoreParam(10))
  # get lib size (SE)
  bam.param  <- ScanBamParam(flag = scanBamFlag(isDuplicate = F))
  bam_counts  <- countBam(bam_files, param = bam.param)
  bam_counts  <- bam_counts$records
  
} else {stop('genome not supported')}
  

# add counts + CPM to orginal gr
  
atac_prom$count.1  <- assay(atac_se[,1])[,1]
  atac_prom$count.1.CPM  <- getCPM(atac_se[,1], lib.size = bam_counts[1])[,1]
atac_prom$count.2  <- assay(atac_se[,2])[,1] 
  atac_prom$count.2.CPM  <- getCPM(atac_se[,2], lib.size = bam_counts[2])[,1]

# remove peaks with no reads in either reps & get avg CPM

atac_prom  <- atac_prom[which(atac_prom$count.1 > 0 & atac_prom$count.2 > 0)]

atac_prom$avg.CPM <- (atac_prom$count.1.CPM + atac_prom$count.2.CPM)/2

# Final filter: keep only regions where average CPM >= 2

atac_prom <- atac_prom[which(round(atac_prom$avg.CPM) >= 2)]
 
# sort and merge final atac regions: resize to the same width as annotated prom regions, and merge/remove those within the set prom width(i.e. 1.2kb)

atac_prom  <- resize(atac_prom, fix = "center", width = 1200)
atac_prom <- sort_ranges(atac_prom_sorted, sort_by=atac_prom_sorted$avg.CPM, gap_size = 1200) 

  save(atac_prom, file = file.path(out_path, paste0(sample_id, "_promoters_unannotated.RData")))

# removes all counts mcols + add names + annotation status
mcols(atac_prom) <- mcols(atac_prom)["name"]
names(atac_prom) <- atac_prom$name
atac_prom$status <- "unannotated"

# combine both and sort by seq levels
prom_regions <- c(prom_regions, atac_prom)
prom_regions <- sortSeqlevels(prom_regions)
prom_regions <- sort(prom_regions)


##--- Export all promoters

  export.bed(prom_regions, file.path(out_path, paste0(sample_id, "_PROM_final.bed")))






