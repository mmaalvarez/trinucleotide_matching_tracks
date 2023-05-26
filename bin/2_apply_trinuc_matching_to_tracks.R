library(tidyverse)
library(data.table)
library(dtplyr)
library(GenomicRanges)
library(rtracklayer)
library(valr)
library("BSgenome.Hsapiens.UCSC.hg19")
library(spgs)
library(conflicted)
conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("lag", "dplyr")
conflict_prefer("between", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("slice", "dplyr")
conflict_prefer("map", "purrr")
conflict_prefer("extract", "magrittr")
conflict_prefer("reduce", "IRanges")
conflict_prefer("desc", "dplyr")
conflict_prefer("reverseComplement", "spgs")

args = commandArgs(trailingOnly=TRUE)

## load files from 2nd process

full_tracks_trinuc32_freq = ifelse(interactive(),
                                   yes = "full_tracks_trinuc32_freq.tsv",
                                   no = args[1]) %>% 
  read_tsv()

matched_tracks = ifelse(interactive(),
                        yes = "matched_tracks.tsv",
                        no = args[2]) %>% 
  read_tsv()

euclidean_score = ifelse(interactive(),
                         yes = "euclidean_score.tsv",
                         no = args[3]) %>% 
  read_tsv() %>% pull(euclidean_score) %>% as.numeric()

sequences = ifelse(interactive(),
                   yes = "sequences.tsv",
                   no = args[4]) %>% 
  read_tsv()


## within each coordinate, randomly remove as many trinucleotides from each type as have been removed with trinuc_matching, and keep track of their positions

rm_n_trinucs_at_random_indices = function(removed_trinucs, trinuc32, sequences){
  
  # init stuff for while loop
  matches = c(1,2)
  niters = 0
  # remove sample_n_matches in case it exists from before, otherwise if not created in the while loop it will pass again as if it were new
  if(exists("sample_n_matches")){
    rm(sample_n_matches)
  }
  
  # "min(diff(matches))<3" should ensure that the removed trinucs do not overlap (maybe not so important anyway)
  while(removed_trinucs>0 & min(diff(matches))<3 & niters<1e6){
    
    matches = as.numeric(gregexpr(trinuc32, sequences)[[1]])
    
    # if there are many "NNNNNN.." in the sequence, it could be that less trinuc match than the "removed_trinucs"
    n_trinucs_to_remove = min(removed_trinucs, length(matches))
    
    sample_n_matches = sort(sample(matches,
                                   n_trinucs_to_remove,
                                   replace = F))
    niters = niters + 1
  }
  
  if(exists("sample_n_matches") && !is.null(sample_n_matches) && niters<1e6){
    return(sample_n_matches)
  } else {
    return(NA)
  }
}

matched_tracks_granges = full_tracks_trinuc32_freq %>%
  column_to_rownames("bin") %>% 
  # subtract orig from remaining counts
  map2_dfc(matched_tracks, ~ .x - .y) %>% 
  # bin names back, and sequences
  rownames_to_column("bin") %>% 
  separate(bin, into = c("seqnames", "start", "end", "name"), extra = "merge") %>% 
  bind_cols(sequences) %>% 
  ### WARNING: updating the 'end' since the actual length of some sequences (obtained with BSGenome::getSeq()) do not match the start and end difference
  mutate(end = as.numeric(start) + nchar(`sequences`)) %>% 
  relocate(sequences) %>% relocate(name) %>% relocate(end) %>% relocate(start) %>% relocate(seqnames) %>% 
  pivot_longer(cols = !matches("seqnames") & !matches("start") & !matches("end") & !matches("name") & !matches("sequences"),
               names_to = 'trinuc32',
               values_to = 'removed_trinucs') %>% 
  # include reverse complement for each trinuc32, so it's actually trinuc64
  mutate(trinuc32 = paste0(trinuc32, "|", spgs::reverseComplement(trinuc32, case="upper"))) %>% 
  ## do the removing thing
  rowwise() %>% 
  mutate(new_end = list(rm_n_trinucs_at_random_indices(removed_trinucs, trinuc32, sequences))) %>% 
  group_by(seqnames, start, end, name) %>% 
  # new_end is now unique positions 1bp before each trinucleotide removed
  summarise(new_end = list((unique(unlist(new_end)) + as.numeric(start)) - 1)) %>% 
  ungroup %>% 
  unnest(new_end) %>% 
  distinct() %>% 
  # in rows just before/at a region in which matching was not performed (typically "target" regions, as they are shorter), or at ending of chr, the new_end == end
  mutate(new_end = ifelse(is.na(new_end),
                          as.numeric(end),
                          new_end)) %>% 
  arrange(seqnames, start, new_end) %>%
  group_by(seqnames, start, end, name) %>% 
  # new_start is 4 nts after the new_end (i.e. skipping the removed trinucleotide) OF THE PREVIOUS ROW
  mutate(new_start = lag(new_end) + 4,
         # in rows at start of chr, the new_start == start
         new_start = ifelse(is.na(new_start),
                            as.numeric(start),
                            new_start),
         width = new_end - new_start + 1,
         strand = "*") %>% 
  ungroup %>% 
  # remove negative or 0 widths (due to overlapping of removed trinucleotides)
  filter(width >= 1) %>% 
  select(seqnames, new_start, new_end, width, strand, name) %>%
  rename("start" = "new_start",
         "end" = "new_end") %>%
  arrange(seqnames, start, end) %>%
  makeGRangesFromDataFrame(keep.extra.columns = T)

export.bed(matched_tracks_granges,
           paste0(filename, "__3ntMatched_euclidean-", round(euclidean_score, 4), ".bed.gz"))
