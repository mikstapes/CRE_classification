#! /home/phan/miniconda3/envs/r_env/bin/Rscript
library(here)
source(here("R", "functions.R"))


bioc_pkg <- c("GenomicFeatures", "rtracklayer", "GenomicAlignments", "Rsamtools", 
              "TxDb.Mmusculus.UCSC.mm10.knownGene", "org.Mm.eg.db")
tidy_pkg <- c("ggplot2", "tidyr", "dplyr", "readr")

suppressMessages(lapply(tidy_pkg, library, character.only = TRUE))
suppressMessages(lapply(bioc_pkg, library, character.only = TRUE)) 


##---Parsing inputs

prom_state <- import.bed(here("segmentation", "e105h", "e105h_4-states_st_1.bed"))

rna_bam_path <- "/project/RNAseq_Mdl_data/Mikie/STAR/chicken/bam-files"

out_path <- here("prom_calling","HH22h")

sample_id <- "HH22h"


# promoters from https://epd.epfl.ch/wwwtmp/mouse_epdnew_p6rZV.bed
# only loaded promoters most representative for each gene

epd_prom <- import.bed(here("promoters_calling", "mm10_epd3.bed"))

# tidy gene name
epd_prom$name <- gsub("_[0-9]*$", "", epd_prom$name)

symb2entrez <- AnnotationDbi::select(org.Mm.eg.db, keys = epd_prom$name, columns = "ENTREZID", keytype = "SYMBOL")

epd_prom$gene_id <- symb2entrez$ENTREZID

# resize prom, keep those with prom-like chromatin
epd_prom_gr <- promoters(epd_prom, upstream = 1000, downstream = 200)
  export.bed(epd_prom_gr, here("promoters_calling", "mm10_epd3_resized.bed"))
  
epd_prom_gr_filt <- epd_prom_gr[overlapsAny(epd_prom_gr, prom_state)]
  export.bed(epd_prom_gr_filt, here("promoters_calling", "mm10_epd3_resized_filt_by_state.bed"))

#- Get TSS with unique TSS

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

Tx <- GenomicFeatures::transcripts(txdb, columns = c("tx_name", "gene_id"))
  Tx <- trim(Tx[!duplicated(start(Tx))]) #131932

###### keep those which overlaps with epd
tx_by_epd <- Tx[findOverlaps(resize(Tx, width = 1000, fix = "start"), epd_prom_gr_filt)@from] #41113

# get width
tx_by_epd$length <- width(tx_by_epd)/1e3

# split by geneid & keep longest

tx_by_epd$gene_id <- as.character(tx_by_epd$gene_id)
split_gr <- splitAsList(tx_by_epd, tx_by_epd$gene_id)
split_max <-   lapply(split_gr, function(x) x[which.max(width(x))])

longest_tx_by_epd <- do.call("c", unname(split_max))
  

######--- Load RNAseq BAM file from both reps and get counts ---######

bam_path <- list.files(path = rna_bam_path, 
                       pattern = ".bam$", full.names = TRUE) %>% 
  setNames(gsub(".bam$", "", basename(.)))

bam_files <- BamFileList(bam_path, asMates = FALSE, yieldSize = 2e6)
tss_se <- summarizeOverlaps(features = longest_tx_by_epd, 
                            reads = bam_files, 
                            singleEnd=TRUE, fragments=FALSE, mode = "Union",
                            BPPARAM = BiocParallel::MulticoreParam(16))


# add counts & TPM to filtered Tx GRanges 

longest_tx_by_epd$count.1  <- assay(tss_se[,1])[,1]; longest_tx_by_epd$count.1.TPM  <- calculateTPM(tss_se[,1])[,1]
longest_tx_by_epd$count.2  <- assay(tss_se[,2])[,1]; longest_tx_by_epd$count.2.TPM  <- calculateTPM(tss_se[,2])[,1]

# add average TPM and save as RData

longest_tx_by_epd$avg.TPM <- (longest_tx_by_epd$count.1.TPM + longest_tx_by_epd$count.2.TPM)/2
    
# Filter 4: Keep TSS where Tx has >= 2 TPM and get promoter region

prom_regions <- promoters(longest_tx_by_epd[which(round(longest_tx_by_epd$avg.TPM) >= 2)],
                          upstream = 1000, downstream = 200)

# save as RData for references
save(prom_regions, file = paste0(out_path,"/active_promoters_from_EPD_TSS.RData")) 

# removes all counts mcols for final prom object
mcols(prom_regions) <- mcols(prom_regions)["tx_name"]
# rename to mcol for consistency
names(mcols(prom_regions)) <- "name"

# add ranges name
names(prom_regions) <- prom_regions$name


##--- Export all promoters

export.bed(prom_regions, paste0(out_path,"/", sample_id, "_PROM_final.bed"))




