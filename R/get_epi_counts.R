# given a genomic regions, plot epigenetic BAM profile over a +/- 1kb window

cran_pkg <- c("here", "tidyverse")
bioc_pkg <- c("rtracklayer", "bamsignals")

suppressMessages(lapply(bioc_pkg, library, character.only = TRUE))
suppressMessages(lapply(cran_pkg, library, character.only = TRUE))

source(here("scripts", "functions.R"))

### Load BAM files 

mBAM  <- list.files(here("data", "bam", "e105h", "rep1"), pattern = ".bam$", full.names = TRUE) %>%
        setNames(gsub(".bam$", "", basename(.)))

gBAM <- list.files(here("data", "bam", "HH22h"), pattern = ".bam$", full.names = TRUE) %>%
        setNames(gsub(".bam$", "", basename(.)))

gBAM  <- gBAM[-6]


### Load enh/prom as GRanges

prom.dir <- here("prom_calling")
enh.dir <- here("enhancer_sets")

mPROM.path <- file.path(prom.dir, "e105h", "e105h_PROM.bed")
gPROM.path <- file.path(prom.dir, "HH22h", "HH22h_PROM.bed")

mENH.path <- file.path(enh.dir, "e105h", "e105h_ENH_filt_from_atac.bed")
gENH.path <- file.path(enh.dir, "HH22h", "HH22h_ENH_filt_from_atac.bed")

mGR  <- lapply(c(mPROM.path, mENH.path), function(x) {
    bed  <- import.bed(x)
    bed  <- getStandardChrom(bed, ref = "mm10")
    return(bed)
})


loadFixedRanges <- function(path, ref.genome, win.size) {
    bed  <- import.bed(path)
    bed  <- getStandardChrom(bed, ref = ref.genome)
    bed <- resize(bed, width = win.size*2, fix = "center")
    return(bed)
}

mGR  <- lapply(c(mPROM.path, mENH.path), 
               loadFixedRanges,
               ref.genome = "mm10", 
               win.size = 1000)
names(mGR)  <- c("prom", "enh")
gGR  <- lapply(c(gPROM.path, gENH.path), 
               loadFixedRanges,
               ref.genome = "galGal6", 
               win.size = 1000)
names(gGR)  <- c("prom", "enh")


#####------ Count Function ------#####


get_bamProfile <- function(bam, bam.list, gr, target_genome) {

  bf  <- names(which(bam.list == bam)); print(bf)

  # filter chroms
  gr <- getStandardChrom(gr, ref = target_genome)

  # get count as vector
  if (target_genome == "mm10") { #SE BAM
    sigs <- sapply(bamProfile(bampath = bam,
    gr = gr,
    shift = 100, binsize = 10,
    mapqual = 10, filteredFlag = 1024, ss = F, verbose =F)@signals,
    function(x) x)
    } else { #PE BAM
    sigs <- sapply(bamProfile(bampath = bam,
    gr = gr, 
    paired.end = "midpoint", binsize = 10,
    mapqual = 10, filteredFlag = 1024, ss = F, verbose = F)@signals,
    function(x) x)
    }
    return(sigs)
  }

#####------ Get Counts as nested Lists ------#####

mSigs  <- lapply(mBAM, function(x){
    mclapply(mGR,
             get_bamProfile,
             bam = x, bam.list = mBAM, target_genome = "mm10",
             mc.cores = 6)
    })

names(mSigs)  <- c("ATAC", "H3K27AC", "H3K27ME3", "H3K4ME1", "H3K4ME3") 

gSigs  <- lapply(gBAM, function(x){
    mclapply(gGR,
             get_bamProfile,
             bam = x, bam.list = gBAM, target_genome = "galGal6",
             mc.cores = 6)
    })

names(gSigs)  <- c("ATAC", "H3K27AC", "H3K27ME3", "H3K4ME1", "H3K4ME3")

save(mSigs, gSigs, file = here("data", "RData", "EP_bamProfile_10bin.RData"))


load(here("data", "RData", "EP_bamProfile_10bin.RData"))

mSigs.wide <- unlist(mSigs, recursive = F)

mSigs.df <- map_dfr(names(mSigs.wide), function(x){
  tibble(sig = rowMeans(mSigs.wide[[x]]), 
         col=rep(x, nrow(mSigs.wide[[x]])),
         bp = c(seq(-1000, -10, by = 10), seq(10, 1000, by = 10))) %>% 
    separate(col, into = c("bam", "cre")) %>% 
    group_by(bam) %>% 
    mutate(sigx = round(scales::rescale(sig, to = c(0,1)), digits = 3)) %>% 
    ungroup()
})

mSigs.df  <- mutate(mSigs.df, sp = "mouse")

gSigs.wide <- unlist(gSigs, recursive = F)

gSigs.df <- map_dfr(names(gSigs.wide), function(x){
  tibble(sig = rowMeans(gSigs.wide[[x]]), 
         col=rep(x, nrow(gSigs.wide[[x]])),
         bp = c(seq(-1000, -10, by = 10), seq(10, 1000, by = 10))) %>% 
    separate(col, into = c("bam", "cre")) %>% 
    group_by(bam) %>% 
    mutate(sigx = round(scales::rescale(sig, to = c(0,1)), digits = 3)) %>% 
    ungroup()
})

gSigs.df  <- mutate(gSigs.df, sp = "chicken")

sigs.df <- rbind(mSigs.df, gSigs.df)

save(sigs.df, file = here("data", "RData","EP_bamProfile_10bin_df.RData"))
