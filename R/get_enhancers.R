#! /home/phan/miniconda3/envs/r_env/bin/Rscript

####################

## Get enhancers AFTER calling active promoter regions
# take ATAC peaks as active enhancer regions, centered by summit


args <- commandArgs(trailingOnly=T)
if (length(args) != 6) stop('Usage: ./get_enhancers.R <genome_build> <atac_peaks.narrowPeak> 
                            <crup.bed/bedGraph> <promoters.bed> <outdir/> <sample_name>')

bioc_pkg <- c("GenomicFeatures", "rtracklayer")
cran_pkg <- c("here", "tidyverse")

invisible(suppressMessages(lapply(c(bioc_pkg,cran_pkg), library, character.only = TRUE)))

source(here("scripts", "functions.R"))

##---Parsing inputs

ref_genome <- args[1]

atac <- rtracklayer::import(args[2], genome = ref_genome)
  atac <- getStandardChrom(atac, ref = ref_genome)

crup <- rtracklayer::import(args[3])

called_promoters <- import.bed(args[4])

out_path <- args[5]

sample_id <- args[6]


##--- Resize atac pk to center by summit (+/- 500bp)

atac <- getPkSummit(atac, win.size = 1000)

##--- Get 'distal' atac peaks

atac_distal <- atac[!overlapsAny(atac, reduce(called_promoters), maxgap = 500)] 

##-- Take overlaps with crup single-enh

out <- atac_distal[findOverlaps(atac_distal, crup, maxgap = 0)@from]

# give ENH 'names'

if (ref_genome == "mm10") names(out) <- paste0("mENH_", seq_len(length(out)))
if (ref_genome == "hg38") names(out) <- paste0("hsENH_", seq_len(length(out)))
if (ref_genome == "galGal6") names(out) <- paste0("gENH_", seq_len(length(out)))

##--- Export filtered enhancers

export.bed(out, file.path(out_path, paste0(sample_id, "_ENH_filtered.bed")))

##--- Export stats 

# Table stats

gr_count <- c(length(called_promoters), 
              length(atac), length(atac_distal), length(crup),
              length(out))

gr_label <- c("Promoters", "ATAC peaks", "ATAC peaks(distal)", 
              "CRUP peaks", "filtered ENH")

tibble(gr_label, gr_count) %>% 
  write_tsv(file.path(out_path, paste0(sample_id, "_ENH_filt_from_atac_stats.tsv")))




