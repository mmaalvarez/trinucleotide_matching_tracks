library(tidyverse)
library(data.table)
library(dtplyr)
library(GenomicRanges)
library(rtracklayer)
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

interactive_work_path = "../work/3c/923348ba7ab03ef48c3e7a6255c9b7/"

full_tracks_trinuc32_freq = ifelse(interactive(),
                                   yes = Sys.glob(paste0(interactive_work_path, "*_full_tracks_trinuc32_freq.tsv")),
                                   no = args[1]) %>% 
  read_tsv()
bin_names = select(full_tracks_trinuc32_freq, bin)

matched_tracks = ifelse(interactive(),
                        yes = Sys.glob(paste0(interactive_work_path, "*_matched_tracks.tsv")),
                        no = args[2]) %>% 
  read_tsv() %>%
  column_to_rownames("bin")

euclidean_score = ifelse(interactive(),
                         yes = Sys.glob(paste0(interactive_work_path, "*_euclidean_score.tsv")),
                         no = args[3]) %>% 
  read_tsv() %>% pull(euclidean_score) %>% as.numeric()

sequences = ifelse(interactive(),
                   yes = Sys.glob(paste0(interactive_work_path, "*_sequences.tsv")),
                   no = args[4]) %>% 
  read_tsv()

filename = ifelse(interactive(),
                   yes = Sys.glob(paste0(interactive_work_path, "*_filename.tsv")),
                   no = args[5]) %>% 
  read_tsv() %>% pull(filename)

# source of rm_n_trinucs_at_random_indices() function
trinuc_matching_source = ifelse(interactive(),
                                yes = "utils.R",
                                no = args[6]) %>% 
  source()


## within each coordinate, randomly remove as many trinucleotides from each type as have been removed with trinuc_matching, and keep track of their positions

matched_tracks_granges = full_tracks_trinuc32_freq %>%
  column_to_rownames("bin") %>% 
  # subtract orig from remaining counts
  map2_dfc(matched_tracks, ~ .x - .y) %>% 
  # bin names back, and sequences
  bind_cols(bin_names) %>% relocate(bin) %>%  
  separate(bin, into = c("seqnames", "start", "end", "name"), extra = "merge") %>% 
  bind_cols(sequences) %>% 
  ### WARNING: updating the 'end' since the actual length of some sequences (obtained with BSGenome::getSeq()) do not match the start and end difference
  mutate(end = as.numeric(start) + nchar(`sequences`)) %>% 
  relocate(sequences) %>% relocate(name) %>% relocate(end) %>% relocate(start) %>% relocate(seqnames) %>% 
  pivot_longer(cols = !matches("seqnames") & !matches("start") & !matches("end") & !matches("name") & !matches("sequences"),
               names_to = 'trinuc32',
               values_to = 'removed_trinucs') %>% 
  # include reverse complement for each trinuc32, so it's actually trinuc64
  mutate(trinuc32 = paste0(trinuc32, "|", reverseComplement(trinuc32, case="upper"))) %>% 
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
