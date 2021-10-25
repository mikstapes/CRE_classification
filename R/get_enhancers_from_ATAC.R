#! /home/phan/miniconda3/envs/r_env/bin/Rscript
library(here)
source(here("scripts", "functions.R"))


args <- commandArgs(trailingOnly=T)
if (length(args) != 4) stop('Usage: ./get_enh_set.R <atac_peaks.bed> <crup_enhancer.bedGraph> <genome_build> <outfile.bed>')

bioc_pkg <- c("GenomicFeatures", "rtracklayer")
tidy_pkg <- c("ggplot2", "tidyr", "dplyr", "readr")
genomes <- c("mm10", "hg38")

suppressMessages(lapply(bioc_pkg, library, character.only = TRUE)) 
suppressMessages(lapply(tidy_pkg, library, character.only = TRUE))

##---Parsing inputs

atac <- import.bed(args[1])
crup <- import.bedGraph(args[2])
  names(crup) <- paste0("crup.enh_", seq_len(length(crup)))
  
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

##--- Get 'distal' enhancers

atac_distal <- atac[!overlapsAny(atac, tss)] 
enhancers <- atac_distal[overlapsAny(atac_distal, crup, maxgap = 250)] 

##--- Export filtered enhancers

export.bed(enhancers, args[4])







