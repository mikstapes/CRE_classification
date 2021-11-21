##
####################### data wrangling functions ####################
## 


makeDFfromGranges <-  function(gr) {
  df <-  tibble(chrom=as.character(seqnames(gr)),
                start=start(gr)-1, # gr is 1-based and BED is 0-based
                end=end(gr)-1, 
                names=names(gr))
  return(df)
}

loadBED4asDF <- function(file_path) {
  df <- read_tsv(file_path, col_names = F) %>% select(1:4)
  colnames(df) <- c("chr", "start", "end", "name")
  return(df)
}

# getRgbBED <- function(df)


##
####################### IPP related functions ####################
## 

#- Helper functions #1: import and tidy raw IPP output 
loadTidyIPP <- function(file_path) {
  x <- read_tsv(file_path, col_names = F)[-1:-2,]
  x <- x %>% 
    dplyr::select(name=X1, ref=X2, coords_direct=X3,coords_djikstra=X4,
           score_direct=X5, score_djikstra=X6,bridges=X15) %>% 
    pivot_longer(cols = 3:6, 
                 names_to = c(".value", "proj_method"),
                 names_pattern = "(.+)_(.+)") %>% 
    separate(ref, into = c("ref_chr", "ref_coords"), sep = ":") %>% 
    separate(coords, into = c("proj_chr", "proj_coords"), sep = ":")
  x$score <- as.numeric(x$score)
  
  return(x)
  
}

#- get projected coords by threshold and projection method, removes unprojected 
# as DF

getProjCoordAsDF <- function(df, direct_thr, dj_thr) {
  
  SC_proj <- df %>% 
    filter(score >= direct_thr & proj_method == "direct" ) %>% 
    mutate(gr = "SC") 
  
  PC_proj <- df %>% 
    filter(score>= dj_thr) %>% 
    filter(proj_method == "djikstra") %>% 
    anti_join(SC_proj, by = "name") %>% mutate(gr = "PC") # takes only those not already called SC, i.e score >= dir_thr
  
  NC_proj <- df %>% 
    filter(proj_method == "djikstra" & score < dj_thr & !is.na(proj_coords)) %>% 
    mutate(gr = "NC")
  
  bind_proj <- rbind(SC_proj, PC_proj, NC_proj) %>% arrange(proj_chr)
  
  return(bind_proj)
}

# Take only djikstra projection of all "unlifted" regions

getProjCoordAsDF_new <- function(df, threshold) {

  PC_proj <- df %>% 
  filter(score>= threshold & proj_method == "djikstra") %>% 
  mutate(gr = "PC")   
  
  NC_proj <- df %>% 
  filter(proj_method == "djikstra" & score < threshold & !is.na(proj_coords)) %>% 
  mutate(gr = "NC")
  
  bind_proj <- rbind(PC_proj, NC_proj) %>% arrange(proj_chr)
  
  return(bind_proj)
}



#- Calculate threshold with max distance input

getIPPThreshold <- function(max_distance, dh) {
  
  thr <- exp(max_distance*log(0.5)/dh)
  thr <- round(thr, 3)
  return(thr)
  
}

#- Convert proj df into GR, proj method, score, and original ref stored as mcols. 
# Output gr contains only autosomes and chrX or Z

IPPtoGRanges <- function(proj_df, target_genome) {
  
  gr <- GRanges(seqnames = proj_df$proj_chr,
                ranges = IRanges(start = as.numeric(proj_df$proj_coords), 
                                 width = 1),
                strand = rep_len("*", length(proj_df$proj_coords)),
                name = proj_df$name,
                score = proj_df$score,
                proj = proj_df$proj_method,
                cons.gr = proj_df$gr) 
  
  gr <- getStandardChrom(gr, target_genome)
  
  return(gr)
  
}
  





#- Helper functions #2: extend projected point +/- 500bp 
extend_proj <- function(x, window_size) {
  
  x$proj_center <- as.numeric(x$proj_center)
  
  x <- x %>% 
    mutate(offset=window_size,  
           proj_start=proj_center-offset, 
           proj_end=proj_center+offset) %>%
    select(-proj_center)
  
  return(x)
  
}

#-  Main function #1: extract projections by threshold and proj method
get_thresholded_proj_new <- function(x, threshold, proj) {
  
  #- tidy IPP, extract only relevant cols 
  x <- tidy_IPP(x)
  
  # -  get coords of thresholded projs based on input
  if (proj=="djikstra") {
    x <- x %>% 
      filter(proj_method=="djikstra" & distance_score >=threshold) %>%
      select(name, coords_djikstra, proj_method) %>% 
      tidyr::separate(coords_djikstra, into = c("proj_chr", "proj_center"), sep = ":")
  }
  
  else if (proj=="direct") {
    x <- x %>% 
      filter(proj_method=="direct" & distance_score >=threshold) %>%
      select(name, coords_direct, proj_method) %>% 
      tidyr::separate(coords_direct, into = c("proj_chr", "proj_center"), sep = ":")
  }
  
  else if (proj=="both") {
    x <- x %>% 
      filter(distance_score >=threshold) %>%
      select(name, coords_direct, coords_djikstra) %>% 
      pivot_longer(cols = starts_with("coords"), ## this will duplicates rows with both projs
                   names_to = "proj_method", 
                   names_prefix = "coords_", 
                   values_to = "proj_coords") %>% 
      tidyr::separate(proj_coords, into = c("proj_chr", "proj_center"), sep = ":")
  }
  ## Extend projected point +/- 250bp
  
  x <- extend_proj(x, window_size = 250) %>% 
    select(proj_chr, proj_start, proj_end, name, proj_method) %>% 
    group_by(name) %>% # removes un-projected points and keep only the highest score at each proj/method
    filter(!is.na(proj_chr) & distance_score == max(distance_score))  %>% ungroup()
  
  return(distinct(x))
  
}

#-  Main function #2: Get all projected points as tidy df
get_all_proj <- function(x){
  
  x <- x %>% select(name, coords_direct, coords_djikstra, distance_score) %>% 
    pivot_longer(cols = starts_with("coords"), ## this will duplicates rows with both projs
                 names_to = "proj_method", 
                 names_prefix = "coords_", 
                 values_to = "proj_coords") %>% 
    tidyr::separate(proj_coords, into = c("proj_chr", "proj_center"), sep = ":") 
  
  x <- extend_proj(x, window_size = 250) %>% 
    select(proj_chr, proj_start, proj_end, name, proj_method, distance_score) %>% 
    group_by(name) %>% # removes un-projected points and keep only the highest score at each proj/method
    filter(!is.na(proj_chr) & distance_score == max(distance_score)) %>% ungroup()
  
  return(distinct(x))
  
}



##
####################### GRanges functions ####################
## 

# Functions: get autosomes + chrom X/Z

getStandardChrom  <- function(x, ref) {
  chrom_list  <- seqlevels(seqinfo(x))
  
  if (ref=="mm10" | ref=="hg38") {
    seqlevels(x, pruning.mode="coarse")  <- chrom_list[grep("^chr[0-9]{,2}$|chrX$",
                                                            chrom_list)]
  } else if (ref == "galGal6") {
    seqlevels(x, pruning.mode="coarse")  <- chrom_list[grep("^chr[0-9]{,2}$|chrZ$",
                                                            chrom_list)]
  } else {
    stop('ref not included')
  }
  
  return(x)
}

# Function: load txdb for standard chroms only

loadtxdb <- function(ref_genome) {
  
  if (ref_genome %in% c("mm10", "hg38")) {
    if (ref_genome == "mm10") pkg <- "TxDb.Mmusculus.UCSC.mm10.knownGene"
    if (ref_genome == "hg38") pkg <- "TxDb.Hsapiens.UCSC.hg38.knownGene"
    suppressMessages(require(pkg, character.only = TRUE))
    assign("txdb", eval(parse(text = pkg)))
    txdb <- getStandardChrom(txdb, ref = ref_genome)
  } else if (ref_genome=="galGal6") {
    txdb <- loadDb(here("txdb", "TxDb.galGal6.ncbiRefseq.sqlite"))
    txdb <- getStandardChrom(txdb, ref = ref_genome)
  } else {
    stop('genome not supported')
  }
  
return(txdb)
    
}


# Function: get the center/midpoint of gr

gr.center <- function(gr1, gr2) {
  
  start(gr) <- rowMeans(cbind(start(gr), end(gr)))
  
  width(gr) <- 1
  
  return(gr)
  
}


# Function: get peak summit and resize as needed

getPkSummit <- function(gr, win.size) {
  
  gr_summit <- GRanges(seqnames = seqnames(gr),
                            ranges = IRanges(start = start(gr) + gr$peak-1,
                                             width = win.size),
                            strand = strand(gr), 
                            name = gr$name,
                            pVal = gr$pValue) 
  
  return(gr_summit)
  
}

# Function: sort gr by a specified mcol values, keep those with highest and merge overlapping/adjacent ranges with lower value
# adapted from CRUP sort_peaks()

sort_ranges <- function(gr, sort_by, gap_size){
  
  # sort gr according to specified mcol (should add check point for int./numeric values)
  gr <- gr[sort(sort_by, decreasing = T, index.return = T)$ix]
  
  count <- 0
  while (length(gr) > (count + 1)) {
    
    count <- count + 1
    overlap.to <- findOverlaps(query = gr[count], subject = gr, maxgap = gap_size)@to
    
    if (length(overlap.to) == 1) next
    
    delete.index <- sort(overlap.to, decreasing = F)[-1]
    gr <- gr[-delete.index]
  }
  #resort by seq levels
  gr <- sort(gr)
  return(gr)
}


# Function CRUP::sort_peaks()

sort_crup_bin <- function(peaks){
  
  # sort peaks according to specified mcol (should add check point for int./numeric values)
  peaks <- peaks[sort(peaks$prob, decreasing = T, index.return = T)$ix]
  
  count <- 0
  while (length(peaks) > (count + 1)) {
    
    count <- count + 1
    overlap.to <- findOverlaps(query = peaks[count], subject = peaks)@to
    
    if (length(overlap.to) == 1) next
    
    delete.index <- sort(overlap.to, decreasing = F)[-1]
    peaks <- peaks[-delete.index]
  }
  
  return(peaks)
}

merge_crup_bin <- function(peaks){
  
  count <- 0
  while (length(peaks) > (count + 1)) {
    
    count <- count + 1
    overlap.to <- findOverlaps(query = peaks[count], subject = peaks)@to
    
    if (length(overlap.to) == 1) next
    
    delete.index <- sort(overlap.to, decreasing = F)[-1]
    peaks <- peaks[-delete.index]
  }
  
  return(peaks)
}

# Function: load in genomic regions as BED and resize to fix by center to fixed with

loadFixedRanges <- function(path, ref.genome, win.size) {
  bed  <- import.bed(path)
  bed  <- getStandardChrom(bed, ref = ref.genome)
  bed <- resize(bed, width = win.size*2, fix = "center")
  return(bed)
}



##
####################### Promoter functions ####################
## 

# Function: calculate TPM of RNAseq counts

calculateTPM <- function(se_count) {
  tx.c  <- assay(se_count)
  tx.len  <- rowRanges(se_count)$length
  
  mat  <- tx.c/tx.len
  tpm.mat  <- t(t(mat)*1e6/colSums(mat))
  
  return(tpm.mat)
}

# Function: calculate CPM of read counts

getCPM <- function(se_count, lib.size) {
  
  count_mat <- assay(se_count)
  cpm.count <- t(t(count_mat)*1e6)/lib.size
  
  retunr(cpm.count)
}

# Function: get unique and longest isoforms for each Tx

getUniqueTx <- function(tx_gr) {
  
  #counter
  i <- 0
  
  #loops through each tx, to find overlapping TSS 
  while (length(tx_gr) > i + 1){
    
      i <- i + 1
      overlap.to.idx <- findOverlaps(resize(tx_gr[i], width = 100, fix = "start"),
                                     tx_gr)@to
      
      if (length(overlap.to.idx) == 1) next
      
      # keep only those with ~same length
      keep.length.idx <- which(ceiling(tx_gr[overlap.to.idx]$length) == ceiling(tx_gr[i]$length))
      overlap.to.idx <- overlap.to.idx[keep.length.idx]
      
      # find the longest of all 'isoforms' and delete the rest if there is more than 1
      delete.idx <- sort(tx_gr[overlap.to.idx]$length, 
                         decreasing = T, index.return  = T)$ix
      
      if (length(delete.idx) == 1)
          delete.idx <- delete.idx
      else 
          delete.idx <- delete.idx[-1]
      
      delete.idx <- overlap.to.idx[delete.idx]
      
      tx_gr <- tx_gr[-delete.idx]
  } 
  
  return(tx_gr)
  
}







##
####################### Other functions ####################
## 

###- Function: get_dups_count()

## Get duplicated/overlapping ENH based on size and overlapping ATAC pks
# input is crup pks, the set width of the final called ENH, and ATAC pks gr
# output is a df of pks count, all vs those overlapping the same ATAC pk

get_dups_count <- function(peaks, size, atac.gr, gap){
  
  #peaks <- resize(peaks, fix = "center", width = size)
  
  # get dups idx
  
  if (gap == TRUE) atac.hits <- findOverlaps(peaks, atac.gr, maxgap = 100)@to
  if (gap == FALSE) atac.hits <- findOverlaps(peaks, atac.gr)@to
    
  atac_dups_idx <- atac.hits[which(duplicated(atac.hits))]
  
  # count ENH overlapping
  dups_count <- sum(countOverlaps(atac.gr[atac_dups_idx], peaks, maxgap = 100))
  
  dups_stats <- tibble(length(peaks), dups_count)
  
  return(dups_stats)
  
}
