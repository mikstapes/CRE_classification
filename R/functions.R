##Granges is 1-based and BED is 0-based
makeDFfromGranges <-  function(gr) {
  df <-  tibble(chrom=as.character(seqnames(gr)),
                start=start(gr)-1,
                end=end(gr)-1, 
                names=names(gr))
  return(df)
}


tidy_IPP <- function(x) {
  x <- x[-c(1:3),]
  x <- x %>% 
    select(name=X1, ref=X2, coords_direct=X3,coords_djikstra=X4,
           score_direct=X5, score_djikstra=X6,bridging_species=X15) %>% 
    pivot_longer(cols = starts_with("score"), 
                 names_to = "proj_method", 
                 names_prefix = "score_", 
                 values_to = "distance_score")
  return(x)
  
}


threshold_IPP <- function(x) {
  x <- x %>% 
    filter(proj_method=="djikstra" & distance_score >=0.975) %>%
    select(name, coords_djikstra) %>% 
    tidyr::separate(coords_djikstra, into = c("chrom_proj", "center_proj"), sep = ":")
  
  x$center_proj <- as.numeric(x$center_proj) 
  
  x <- x %>% 
    mutate(offset=250,  
           start_proj=center_proj-offset, 
           end_proj=center_proj+offset) %>%
    select(chrom_proj, start_proj, end_proj, name)
  
  return(x)
  
}