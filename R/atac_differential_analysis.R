## Differential analysis for tissue-specificity 'enhancers'
## Input are genrich peaks from replicated and filt.bam file after mapping (dedups, mapq10)

########## Load packages #######

library(bamsignals)
library(Rsamtools)
library(GenomicRanges)
library(rtracklayer)
library(GenomicAlignments)
library(EnrichedHeatmap)

library(circlize)

library(csaw)
library(edgeR)
library(limma)

library(here)
library(dplyr)
library(ggplot2)


##----------Global var and functions-----------

# standard chroms only
chr_gg <- paste0("chr", c(1:33, "Z"))
chr_mm <- paste0("chr", c(1:19, "X"))

#read parameters for PE atac
param_gg <- readParam(pe = 'both', restrict = chr_gg)
param_mm <- readParam(pe = 'both', restrict = chr_mm)

##----------Load Data-----------
#pks_heart <- import.bed(here("data", "atac-HH24h_peaks.bed"))
#pks_limb <-  import.bed(here("data", "atac-HH24fl_peaks.bed"))

#pks_all <- union(pks_heart, pks_limb)

bam_files <- list.files(here("data", "bam", "mm"), 
                        pattern = "bam$", full.names=TRUE) %>% 
setNames(gsub(".bam$", "", basename(.)))
#name atac-sample-r1/2

###--------------------------

#                          param = param_mm, 
#                          BPPARAM = BiocParallel::MulticoreParam(4))
#
##count BAM in 200bp windows
atac_win.counts <- csaw::windowCounts(bam_files, 
                                param = param_mm,
                                width = 200,
                                BPPARAM = BiocParallel::MulticoreParam(6))
#
##count BAM within a 2kb neighborhood for local enrichment
#atac_local.bkgr <- csaw::regionCounts(bam_files, 
#                                      regions = suppressWarnings(resize(rowRanges(atac_win.counts),
#                                                                        width = 2000, 
#                                                                        fix = "center")),
#                                      param = param_mm,
#                                      BPPARAM = BiocParallel::MulticoreParam(6))

atac_global.bkgr<- csaw::windowCounts(bam_files,
                                        param = param_mm, 
                                        width = 2000, bin = TRUE,
                                        BPPARAM = BiocParallel::MulticoreParam(6))

# Get filter stats based on read counts of flanking 2kb regs around window counts, 
# stat = win.counts - bkground.counts at each window

#      local.stat <- filterWindowsLocal(atac_win.counts, atac_local.bkgr)
#      keep.local <- local.stat$filter > log2(3)
#      
      global.stat <- filterWindowsGlobal(atac_win.counts, atac_global.bkgr)
      keep.global <- global.stat$filter > log2(3)
#
#    # Plot logFC from local control
#    pdf(here("plots", "hist_local_FC.pdf"))
#    hist(local.stat$filter, xlab = "LogFC from 2kb local background",
#         breaks=100, main="", 
#         col="blue", xlim=c(0,5))
#    abline(v=log2(3), col="red", lwd=2)
#    dev.off()
    
    # Plot logFC from global control
#    pdf(here("plots", "hist_global_FC.pdf"))
#    hist(global.stat$filter, xlab = "LogFC from global background",
#         breaks=100, main="", 
#         col="blue", xlim=c(0,5))
#    abline(v=log2(3), col="red", lwd=2)
#    dev.off()




# Keep only regions within at least a 3 fold change above global background
keep  <- keep.global
atac_filt.counts <- atac_win.counts[keep,]

##----------Normalize Counts-----------      

# Loess normalization of peaks
#pks.count_norm <- csaw::normOffsets(pks.count, se.out = TRUE)

# Normalization of 'local' enriched bam, i.e. de novo peaks
atac_filt.counts_norm <- csaw::normOffsets(atac_filt.counts, se.out = TRUE)

##--------------DA--------------  

#pks <- asDGEList(pks.count_norm)
aln <- asDGEList(atac_filt.counts_norm)

##---(1)Design matrix from aln counts
aln$samples$group <- c("E105-heart", "E105-heart", "E105-FL", "E105-FL") 

aln_design <- model.matrix(~0+group, data=aln$samples)
aln_design


## Plot MDS from limma ##
    # Plot samples on 2-min scatter plot based on log2 FC 
    # input is normalized count from aln, with CPM scaling
aln_norm.cpm <- calculateCPM(atac_filt.counts_norm, use.offsets = TRUE)

    pdf(here("plots", "MDS_mouse-atac-CPMnorm.pdf"))
    plotMDS(aln_norm.cpm, labels=names(bam_files),
            col=c("#1A936F", "#C6DABF", "#C1666B", "#BF211E"))
    dev.off()

# Estimate binomial dispersion of all samples
# input is DGElist of norm counts

aln <-  edgeR::estimateDisp(aln, aln_design)


    ## Plot BCV from edgeR ##
 #   pdf(here("plots", "BCV_24H-atac.pdf"))
 #   plotBCV(aln)
 #   dev.off()

## Fit negative binomial lin model over counts
aln.fit <- glmQLFit(aln, aln_design, robust=TRUE)
    
  ## Plot the fit
#pdf(here("plots", "Disp_24H-atac-fit.pdf"))
#plotQLDisp(aln.fit)
#dev.off()


## here's the actual meat of DA analysis

colnames(aln_design) <- c("E105.FL", "E105.Heart")
grp_names <-  colnames(aln_design)
contrast_pairs <- c(paste0(grp_names[1], "-", grp_names[2]),
                   paste0(grp_names[2], "-", grp_names[1]))


test_contrast <- function(x, mat, fit, counts) {
  con <- makeContrasts(contrasts = x, levels = colnames(mat))
  message(con)
  results <-  glmQLFTest(fit, contrast = con)
  merged <- csaw::mergeResults(counts, results$table, tol=250, 
                               merge.args = list(max.width = 1000))
  return(merged)
}

DA.list <- mapply(test_contrast, x = contrast_pairs, 
                  MoreArgs = list(mat = aln_design,
                                  fit = aln.fit, 
                                  counts = atac_filt.counts_norm))


names(DA.list) <- contrast_pairs

# extract edgr diff counts as gr with stats as metadata
edgr_to_gr <- function(x, counts) {
  gr <- x$regions #region within the results DF is a gr obj
  mcols(gr) <- DataFrame(x$combined,
                     best.logFC = x$best$logFC,
                     logCPM = x$best$logCPM) 
  return(gr)
}

DA.gr <- lapply(DA.list, edgr_to_gr, counts = atac_filt.counts_norm) #list of pks + stats
## Write diff and mixed regions as bed file
#n: names(gr) to vectorize over

lapply(names(DA.gr), function(x) {
  for (d in c("up", "down")) {
    gr <- DA.gr[[x]]
    gr <-  gr[gr$FDR < 0.05]
    gr <- gr[gr$direction == d]
    gr <- gr[order(gr$best.logFC)]
    diff_name <- paste0("DA_atac-", x, "_", d,".bed")
    export.bed(gr, con = file.path(here(),"/diff_analysis", diff_name))
}})

lapply(names(DA.gr), function(x) {
    gr <- DA.gr[[x]]
    gr <-  gr[gr$FDR > 0.05]
    gr <- gr[order(gr$best.logFC)]
    ns_name <- paste0("DA_atac-", x, "_","ns.bed")
    export.bed(gr, con = file.path(here(),"/diff_analysis", ns_name))
})

########### Plotting #######

### (1) MA plot

#heart_mpks <- DA.gr[["HH24.Heart-HH24.FL"]]
#heart_mpks <- heart_mpks[order(heart_mpks$best.logFC)]
#heart_mpks$sig <- "n.s"
#heart_mpks$sig[heart_mpks$FDR < 0.05] <- "significant"
#
#pdf(here("plots", "MA_heart_limb.pdf"))
#data.frame(heart_mpks) %>% 
#  ggplot(aes(x = logCPM, y = best.logFC, col = factor(sig, levels = c("n.s", "significant")))) +
#  geom_point() + scale_color_manual(values = c("#99ADAC", "#E0BE15")) + 
#  geom_smooth(inherit.aes=F, aes(x = logCPM, y = best.logFC), method = "loess", se = FALSE) + 
#  geom_hline(yintercept = 0) + labs(col = NULL)
#dev.off()
