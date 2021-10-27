### Getting bam counts around projection$

cran_pkg <- c("here", "dplyr", "readr", "tidyr", "ggplot2", "purrr", "forcats", "scales")
bioc_pkg <- c("rtracklayer", "bamsignals")

suppressMessages(lapply(bioc_pkg, library, character.only = TRUE)) 
suppressMessages(lapply(cran_pkg, library, character.only = TRUE))

source(here("scripts", "functions.R"))


#####- FUNCTIONS -#####




#####- Load IPP GRanges -#####

load(here("IPP", "proj_dmax_1_5kb.Rdata"))

#####- Read BAM -#####

h22.bamf  <- list.files(here("data", "bam", "HH22h"), pattern = ".bam$", full.names = TRUE) %>%
	setNames(gsub(".bam$", "", basename(.)))

h10.bamf  <- list.files(here("data", "bam", "e105h", "rep1"), pattern = ".bam$", full.names = TRUE) %>%
	setNames(gsub(".bam$", "", basename(.)))

ggr  <- gENH_proj.GRlist[['1000']]
ggr  <- resize(ggr, fix = "center", width = 1000)

k27ac.sigs  <- bamsignals::bamProfile(bampath = h10.bamf[['e105h.H3K27AC']], 
		       ggr, shift = 100 , filteredFlag = 1024, binsize = 1000, verbose = F)
k27ac.mat  <- bamsignals::alignSignals(k27ac.sigs)

ggr$k27ac.sigs  <- bamsignals::alignSignals(k27ac.sigs)

ggr1kb.df  <- tibble(ref = ggr$ref, gr = ggr$cons.gr, sigs = ggr$k27ac.sigs)

pdf(here("IPP", "plots", "test_gENH_d1000_k27ac_rescaled_75th.pdf"))


ggr1kb.df %>% filter(sigs <= quantile(ggr1kb.df$sigs, probs = 0.75)) %>%
	ggplot(aes(x=fct_rev(gr), y = scales::rescale(sigs, to = c(0,1)), fill = fct_rev(gr))) +
	geom_boxplot()

dev.off()


getIPPbam <- function(bam, gr, target_genome) {
  
  bf  <- names(which(bam.list == bam)); print(bf)

  # filter chroms
  gr <- getStandardChrom(gr, ref = target_genome)
  
  # get count as vector
  if (target_genome == "mm10") { #SE BAM
    sigs <- sapply(bamProfile(bampath = bam,
                              gr = gr,
                              binsize = width(gr[1]), 
                              shift = 100,
                              mapqual = 10, filteredFlag = 1024, ss = F, 
			      verbose =F)@signals, 
                   function(x) x) 
    } else { #PE BAM
    sigs <- sapply(bamProfile(bampath = bam,
                              gr = gr,
                              binsize = win.size, 
                              paired.end = "midpoint",
                              mapqual = 10, filteredFlag = 1024, ss = F, 
			      verbose = F)@signals, 
                   function(x) x)  
    }
    return(sigs)
  }

gr.list  <- lapply(gENH_proj.GRlist, function(x) {
  x <- x[which(x$cons.gr == "PC")]
  x <- resize(x, fix="center", width = 1000)})


bam.list  <- h10.bamf
count.list <- lapply(h10.bamf, function(x){
  r <- mclapply(gr.list, 
                   getIPPbam,
                   bam = x, 
		   target_genome = "mm10", mc.cores = 4)
})


names(count.list)  <- c("H3K27AC", "H3K27ME3", "H3K4ME1", "H3K4ME3")


save(count.list, file =  here("data", "gENH_proj_bamcounts.RData"))
