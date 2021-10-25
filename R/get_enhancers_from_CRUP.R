#! /home/phan/miniconda3/envs/r_env/bin/Rscript
library(here)
source(here("scripts", "functions.R"))

# get 1100bp predicted enh candidate after incorporating ATAC pks 
args <- commandArgs(trailingOnly=T)
if (length(args) != 5) stop('Usage: ./get_enh_set.R <atac_peaks.bed> <crup_predictions.rds> <genome_build> <outdir/> <sample_name>')

bioc_pkg <- c("GenomicFeatures", "rtracklayer")
tidy_pkg <- c("ggplot2", "tidyr", "dplyr", "readr")
genomes <- c("mm10", "hg38")

suppressMessages(lapply(bioc_pkg, library, character.only = TRUE)) 
suppressMessages(lapply(tidy_pkg, library, character.only = TRUE))

##---Parsing inputs

atac <- import.bed(args[1])

crup <- readRDS(args[2])

ref_genome <- args[3]


##---- Load txdb, input data & filter only standard chroms (autosomes + chrX or chrZ(gg))

if (ref_genome %in% genomes) {
  if (ref_genome == "mm10") pkg <- "TxDb.Mmusculus.UCSC.mm10.knownGene"
  if (ref_genome == "hg38") pkg <- "TxDb.Hsapiens.UCSC.hg38.knownGene"
  suppressMessages(require(pkg, character.only = TRUE))
  assign("txdb", eval(parse(text = pkg)))
  txdb <- getStandardChrom(txdb, ref = ref_genome)
  atac <- getStandardChrom(atac, ref = ref_genome)
} else if (ref_genome=="galGal6") {
  txdb <- loadDb("/project/MDL_Ibrahim/MP_all/annotations/gg6/TxDb.galGal6.ncbiRefseq.sqlite")
  txdb <- getStandardChrom(txdb, ref = ref_genome)
  atac <- getStandardChrom(atac, ref = ref_genome)
} else {
  stop('genome not supported')
}

##--- Define TSS/promoters from TxDb

tss <- GenomicFeatures::promoters(txdb, upstream = 1000, downstream = 200)
tss <- trim(tss)  #trim 6 out of bound tss/promoters
tss <- granges(tss[!duplicated.GenomicRanges(tss)])

##--- Get 'distal' eatac peaks

atac_distal <- atac[!overlapsAny(atac, tss)] 

##--- Load crup bins over cutoff prob as 'peaks'

cutoff <- 0.5

pks_filt_all <- crup[which(crup$prob > cutoff)]

# get all 100bp bins overlapping atac pks outside proximal region

pks_filt_distal <- pks_filt_all[findOverlaps(pks_filt_all, atac_distal, maxgap = 250)@from]

# create 1100 pks regions (windows N = 5)

extend_count <- 5

start(pks_filt_distal) <- start(pks_filt_distal) - extend_count*100
width(pks_filt_distal) <- extend_count*100*2 + 100

## final peak set (sort and merge like CRUP)

# split bins into chroms for sorting and merging + parallelization 

out <- mclapply(split(pks_filt_distal, seqnames(pks_filt_distal)),
                sort_peaks,
                mc.cores = 8)

# merged all chroms into a single gr
out <- do.call("c", unname(out))

# sort ENH by seqlevels

out <- sortSeqlevels(out)
out <- sort(out)

# give ENH 'names'

if (ref_genome == "mm10") names(out) <- paste0("mENH_", seq_len(length(out)))
if (ref_genome == "hg38") names(out) <- paste0("hENH_", seq_len(length(out)))
if (ref_genome == "galGal6") names(out) <- paste0("gENH_", seq_len(length(out)))
            
##--- Export filtered enhancers

export.bed(out, paste0(args[4], args[5],"_ENH_final.bed"))

##--- Export stats 

# Table stats
gr_count <- c(length(tss), length(atac), length(atac_distal),
              length(pks_filt_all), length(out))
gr_label <- c("Promoters", "ATAC peaks", "ATAC peaks(distal)", 
              "CRUP peaks", "ENH")

tibble(gr_label, gr_count) %>% write_tsv(paste0(args[4], args[5], "_ENH_stats.tsv"))








