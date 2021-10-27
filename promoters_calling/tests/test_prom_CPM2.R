# get counts for all promoter regions without annotation to determine cut-offs


library(here)
source(here("scripts", "functions.R"))


#args <- commandArgs(trailingOnly=T)
#if (length(args) != 6) stop('Usage: ./test_prom_CPM.R <genome_build> <prom_state.bed> <atac.bed> <RNA_bam_path> <sample_name> <outdir/>')

bioc_pkg <- c("GenomicFeatures", "rtracklayer", "GenomicAlignments", "Rsamtools")
tidy_pkg <- c("ggplot2", "tidyr", "dplyr", "readr")
genomes <- c("mm10", "hg38")

suppressMessages(lapply(bioc_pkg, library, character.only = TRUE)) 
suppressMessages(lapply(tidy_pkg, library, character.only = TRUE))

##---Parsing inputs

ref_genome <- "galGal6"

prom_state <- import.bed(paste0(here(), "/chromseg/HH22h/bed/HH22h_4-states_st_1.bed"))

atac <- import.bed(paste0(here(), "/data/peaks/atac-HH22h_peaks.bed"))

rna_bam_path <- "/project/RNAseq_Mdl_data/Mikie/STAR/chicken/bam-files"


sample_id <- "HH22h"

out_path <- paste0(here("prom_calling"), "/")




# load TxDb

genomes <- c("mm10", "hg38")

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


# Get transcripts with unique TSS
Tx <- GenomicFeatures::transcripts(txdb)
  Tx <- trim(Tx[!duplicated(start(Tx))])
  
# Load RNAseq BAM

bam_path <- list.files(rna_bam_path, 
		       pattern = ".bam$", full.names = TRUE)


prom_noRef <- prom_state[!overlapsAny(prom_state,
                                      resize(Tx, width = 1, fix = "start"), maxgap = 0)]

#- get ATAC peaks as possible annotated promoters

atac_prom <- atac[overlapsAny(atac, prom_noRef, maxgap = 0)]

# Get counts at peaks

if (ref_genome == "galGal6") {
  
  bam_files <- BamFileList(bam_path, asMates = TRUE)
  
  atac_se <- summarizeOverlaps(features = atac_prom, 
                               reads = bam_files, 
                               singleEnd=FALSE, fragments=TRUE, mode = "Union",
                               BPPARAM = BiocParallel::MulticoreParam(16))
 # get lib size 
  
  bam.param  <- ScanBamParam(flag = scanBamFlag(isPaired  = T, isDuplicate = F))
  bam_counts  <- countBam(bam_files, param = bam.param)
  bam_counts  <- bam_counts$records/2

} else {
  
  bam_files <- BamFileList(bam_path, asMates = FALSE)
  atac_se <- summarizeOverlaps(features = atac_prom, 
                               reads = bam_files, 
                               singleEnd=TRUE, fragments=FALSE, mode = "Union",
                               BPPARAM = BiocParallel::MulticoreParam(8))
  bam.param  <- ScanBamParam(flag = scanBamFlag(isDuplicate = F))
  bam_counts  <- countBam(bam_files, param = bam.param)
  bam_counts  <- bam_counts$records


}


# add counts to original gr

atac_prom$count.1  <- assay(atac_se[,1]); atac_prom$count.2  <- assay(atac_se[,2])

# filter out peaks with no read counts in either reps
atac_prom  <- atac_prom[which(atac_prom$count.1 > 0 & atac_prom$count.2 > 0)]

# add CPM counts
atac_prom$count.1.CPM  <- getCPM(atac_se[,1], lib.size=bam_counts[1])
atac_prom$count.2.CPM  <- getCPM(atac_se[,2], lib.size=bam_counts[2])

  atac_prom$avg.CPM <- (atac_prom$count.1.CPM + atac_prom$count.2.CPM)/2


# plot counts distribution
  
count1.raw <- unname(atac_prom$count.1[,1])
count2.raw <- unname(atac_prom$count.2[,1])
count1.CPM <- unname(atac_prom$count.1.CPM[,1])
count2.CPM <- unname(atac_prom$count.2.CPM[,1])

counts_all <- c(count1.raw, count1.CPM, count2.raw, count2.CPM)
rep <- c(rep("rep1", length(atac_prom)*2), rep("rep2", length(atac_prom)*2))

grp <- c(rep("raw_count", length(atac_prom)), rep("CPM_count", length(atac_prom)),
	 rep("raw_count", length(atac_prom)), rep("CPM_count", length(atac_prom)))

counts_df <- tibble(rep, grp, counts_all)
counts_df$grp <- factor(counts_df$grp)
counts_df$rep <- factor(counts_df$rep)


pdf(paste0(out_path, sample_id, "_counts_dist.pdf"))

counts_df %>% ggplot(aes(x=rep, y=counts_all, fill = rep)) + 
  geom_violin() +
  geom_boxplot(width = 0.2) +
  facet_wrap(vars(grp))

dev.off()

#plot only CPM and scale the counts where 97.5th quantile = 100
cpm_df  <- counts_df %>% filter(grp == "CPM_count") %>% select(-grp)
cpm_df  <- cpm_df %>% 
	mutate(counts_scaled=case_when(counts_all > quantile(cpm_df$counts_all, 0.975) ~ quantile(cpm_df$counts_all, 0.975),
				       TRUE ~ counts_all))


pdf(paste0(out_path, sample_id, "_CPM_scaled_dist.pdf"))

cpm_df %>% ggplot(aes(x=rep, y=counts_scaled, fill = rep)) +
	geom_violin() +
	geom_boxplot(width = 0.3)
dev.off()

# save RData
save(atac_prom, 
     file = paste0(out_path, sample_id, "_prom_no_ref_counts.RData"))


# set cutoff avg CPM >= 2

unannotated_prom  <- atac_prom[which(ceiling(atac_prom$avg.CPM) >= 2)]

export.bed(unannotated_prom, here("prom_calling", "HH22h_unannoted_prom.bed"))
