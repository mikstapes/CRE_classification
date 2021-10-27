library(here)
source(here("scripts", "functions.R"))

args <- commandArgs(trailingOnly=T)
if (length(args) != 5) stop('Usage: ./get_enh_set.R <genome_build> <atac_peaks.narrowPeak> <crup_predictions.rds> <outdir/> <sample_name>')

bioc_pkg <- c("GenomicFeatures", "rtracklayer")
tidy_pkg <- c("ggplot2", "tidyr", "dplyr", "readr")

suppressMessages(lapply(bioc_pkg, library, character.only = TRUE)) 
suppressMessages(lapply(tidy_pkg, library, character.only = TRUE))

genomes <- c("mm10", "hg38")

suppressMessages(lapply(bioc_pkg, library, character.only = TRUE)) 
suppressMessages(lapply(tidy_pkg, library, character.only = TRUE))

##---Parsing inputs

ref_genome <- args[1]

atac <- import(args[2],
               format = "narrowPeak",
               genome = ref_genome)

crup <- readRDS(args[3])


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

##--- Get 'distal' atac peaks

atac_distal <- atac[!overlapsAny(atac, tss)] 

##--- Load crup bins over cutoff prob as 'peaks'

cutoff <- 0.5

pks_all <- crup[which(crup$prob > cutoff)]

# get bins outside 

pks_distal <- pks_all[findOverlaps(pks_all, atac_distal, maxgap = 250)@from]

#compare ENH sizes, export stats as tsv

# pre-set different sizes for ENH
enh.width <- c(500, 700, 900, 1100)

# empty vector and tibble for loops to fill in 
count <- numeric(length(enh.width))
dups_count <- numeric(length(enh.width))
dups_df <- tibble(count, dups_count)

# for each set enh size/width, calculate # called ENH, and # of dups

for (i in seq_along(enh.width)) {
  
  # resize 
  pks_out <- resize(pks_distal, fix = "center", width = enh.width[i])
  
  # parallel sort by descending score and merge by overlap
  pks_out <- mclapply(split(pks_out, seqnames(pks_out)),
                         sort_peaks,
                         mc.cores = 8)
  
  # re-combine split gr into a single gr
  pks_out <- do.call("c", unname(pks_out))
    pks_out <- sort(pks_out)
  
  # get dups stats df
  dups_df[i,] <-  get_dups_count(peaks = pks_out, 
                                 size = enh.width[i], 
                                 atac.gr = atac_distal,
                                 gap = TRUE)
}

## writes the number of ENH: total, dups, final (no dups)
dups_df <- mutate(dups_df, final_count = count - dups_count/2) 
    
    # write tsv
    dups_df %>% write_tsv(paste0(args[4], args[5], "_ENH_size_overlap_stats.tsv"))


