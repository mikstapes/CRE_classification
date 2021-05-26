##Granges is 1-based and BED is 0-based
makeDFfromGranges <-  function(gr) {
  df <-  tibble(chrom=as.character(seqnames(gr)),
                start=start(gr)-1,
                end=end(gr)-1, 
                names=names(gr))
  return(df)
}

#- Helper functions #1: tidy raw IPP output 
tidy_IPP <- function(x) {
  x <- x %>% 
    select(name=X1, ref=X2, coords_direct=X3,coords_djikstra=X4,
           score_direct=X5, score_djikstra=X6,bridges=X15) %>% 
    pivot_longer(cols = starts_with("score"), 
                 names_to = "proj_method", 
                 names_prefix = "score_", 
                 values_to = "distance_score")
  return(x)
  
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
  
  x <- tidy_IPP(x) %>% 
    select(name, coords_direct, coords_djikstra, distance_score) %>% 
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