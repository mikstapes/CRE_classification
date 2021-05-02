#! /home/phan/miniconda3/envs/r_env/bin/Rscript

args <- commandArgs(trailingOnly=T)
if (length(args) != 4) stop('Usage: ./plot_heatmaps.R <histone.bam or (normalized) hitones.bw> <window size(bp)> <outfile.pdf> <target.bed or genome build for TSS>')

######## Load Libraries #######
bioc_pkg <- c("GenomicFeatures", 
              "rtracklayer", 
              "EnrichedHeatmap", 
              "GenomicAlignments", 
              "Rsamtools",
              "bamsignals")

cran_pkg <- c("circlize")

suppressMessages(lapply(bioc_pkg, require, character.only = TRUE, quietly = TRUE)) 
suppressMessages(lapply(cran_pkg, require, character.only = TRUE, quietly = TRUE))

######## Functions & vars #######
genomes <- c("mm10", "hg38") 
col_fun <-  colorRamp2(quantile(mat2, c(0, 0.99)), c("#424874", "#83AFAA"))

# get_standard_chroms <- function(gr) {
#   si <- seqinfo(gr)
#   seqlevels(gr, pruning.mode = "coarse") <- seqlevels(si)[grep("^chr[0-9]{,2}$|chrX$|chrZ$",
#                                                                           seqlevels(si))]
#   return(seqlevels(gr))
# }

get_tss <- function(ref_genome) {
  #--Load txdb for appropriate genome build
  if (ref_genome %in% genomes) {
    if (ref_genome == "mm10") pkg <- "TxDb.Mmusculus.UCSC.mm10.knownGene"
    if (ref_genome == "hg38") pkg <- "TxDb.Hsapiens.UCSC.hg38.knownGene"
    suppressMessages(require(pkg, character.only = TRUE))
    assign("txdb", eval(parse(text = pkg)))
  } else (ref_genome=="galGal6") {
    txdb <- loadDb("/project/MDL_Ibrahim/MP_all/annotations/gg6/TxDb.galGal6.ncbiRefseq.sqlite")
  }
  #--Get TSS
  tss <- resize(transcripts(txdb), width = 1, fix = "start")
  seqlevels(tss, pruning.mode = "coarse") <- seqlevels(seqinfo(tss))[grep("^chr[0-9]{,2}$|chrX$|chrZ$",
                                                                          seqlevels(seqinfo(tss)))]
  return(tss)
}
  

######## Load target and signals #######
signals.path <- args[1]
extend_size <- args[2]

target <- args[4]

#--(1) Get the target regions, import either as bed or default to TSS
if (endsWith(target, '.bed')) {
    target_gr <- rtracklayer::import.bed(target)
    target_gr <- resize(target_gr, width = 1, fix = "center")
    seqlevels(target_gr, pruning.mode = "coarse") <- seqlevels(seqinfo(target_gr))[grep("^chr[0-9]{,2}$|chrX$|chrZ$",
                                                                            seqlevels(seqinfo(target_gr)))]
} else {
    print('No bed file of target regions provided, heatmaps will center around TSS \
            from provided genome build')
    target_gr <- get_tss(target)
}

  # Get Target extended regions for quick loading of bw
  target_gr.extended <- resize(target_gr, 
                               width = extend_size*2,
                               fix = "center")

#--(2) Get signals as either norm. counts, i.e. bigwigs or filtered bam alignments
if (!endsWith(signals.path, '.bam') | !endsWith(signals.path, '.bw')) {
  stop('Wrong input file for signals. Requires either a bigwig .bw or .bam file')
} else if (endsWith(signals.path, '.bam')) {
  reads <- readGAlignments(signals.path)
  reads.bf <- BamFile(signals.path)
  reads.si <- seqinfo(reads.bf)  #Seqinfo 
  #get coverage after extending reads to 200bp output is a RleList obj
  reads_cov <- coverage(resize(granges(reads), 
                               width = 200, fix = "start"))
  #covert to GRanges
  signals_gr <- as(reads_cov, "GRanges")
  #keep only 'standard chromosomes'
  seqlevels(signals, pruning.mode = "coarse") <- seqlevels(reads.si)[grep("^chr[0-9]{,2}$|chrX$|chrZ$",
                                                                          seqlevels(reads.si))]
} else {
  signals_gr <- rtracklayer::import.bw(signals.path,
                                    selection = BigWigSelection(target_gr.extended))
}
  
#--(3) Norm to matrix and plot heatmap around target
mat.norm <- normalizeToMatrix(signals_gr, 
                              target_gr, 
                              value_column = "score", 
                              background = 0,
                              mean_mode = "w0", 
                              extend = extend_size,
                              smooth = TRUE)
ehmap <- EnrichedHeatmap(mat.norm,
                       col = col_fun,
                       pos_line = FALSE,
                       border = FALSE,
                       use_raster = TRUE, raster_quality = 10, raster_device = "png",
                       rect_gp = gpar(col = "transparent"),
                       heatmap_legend_param = list(
                         legend_direction = "horizontal",
                         title = "normalized counts",
                         labels_gp = gpar(fontsize = 15)),
                       top_annotation = HeatmapAnnotation(
                         enriched = anno_enriched(
                           gp = gpar(col = "dark blue", lty = 1, lwd=3),
                           col="dark blue", 
                           axis_param = list(gp = gpar(fontsize=15))
                         )))

##---(4) Save plot to disk 

pdf(args[3])

draw(ehmap,
     heatmap_legend_side = "bottom",
     annotation_legend_side = "bottom",
     padding = unit(c(4, 4, 4, 4), "mm"))

dev.off()

  





