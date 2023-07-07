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


### NEW fraction of the input files to parse in this iteration
start_int = ifelse(interactive(),
                   yes = "1",
                   no = args[1]) %>% 
  as.numeric

end_int = ifelse(interactive(),
                 yes = "140000",
                 no = args[2]) %>% 
  as.numeric


## load files from 2nd process

interactive_work_path = "../../work/54/2c44d8fe769335dc4d4a73ad46a0e3/"

full_tracks_trinuc32_freq = ifelse(interactive(),
                                   yes = Sys.glob(paste0(interactive_work_path, "*_full_tracks_trinuc32_freq.tsv")),
                                   no = args[3]) %>% 
  read_tsv() %>%
  ### NEW fraction of the input files to parse in this iteration
  slice(start_int:end_int)
  
bin_names = select(full_tracks_trinuc32_freq, bin)

matched_tracks = ifelse(interactive(),
                        yes = Sys.glob(paste0(interactive_work_path, "*_matched_tracks.tsv")),
                        no = args[4]) %>% 
  read_tsv() %>%
  ### NEW fraction of the input files to parse in this iteration
  slice(start_int:end_int) %>% 
  column_to_rownames("bin")

euclidean_score = as.numeric(args[5])

sequences = ifelse(interactive(),
                   yes = Sys.glob(paste0(interactive_work_path, "*_sequences.tsv")),
                   no = args[6]) %>% 
  read_tsv() %>%
  ### NEW fraction of the input files to parse in this iteration
  slice(start_int:end_int)

filename = args[7]

# source of rm_n_trinucs_at_random_indices() function
trinuc_matching_source = ifelse(interactive(),
                                yes = "../../bin/utils.R",
                                no = args[8]) %>% 
  source()
gc()



## within each coordinate, randomly remove as many trinucleotides from each type as have been removed with trinuc_matching, and keep track of their positions

# subtract orig from remaining counts
matched_tracks_granges = full_tracks_trinuc32_freq %>%
  column_to_rownames("bin") %>% 
  map2_dfc(matched_tracks, ~ .x - .y)
## free up memory
rm(full_tracks_trinuc32_freq) ; rm(matched_tracks) ; gc()
## continue
matched_tracks_granges = matched_tracks_granges %>% 
  # bin names back, and sequences
  bind_cols(bin_names) %>% relocate(bin) %>%  
  separate(bin, into = c("seqnames", "start", "end", "name"), extra = "merge") %>% 
  bind_cols(sequences)
## free up memory
rm(bin_names) ; rm(sequences) ; gc()
## continue
matched_tracks_granges = matched_tracks_granges %>% 
  lazy_dt() %>% 
  ### WARNING: updating the 'end' since the actual length of some sequences (obtained with BSGenome::getSeq()) do not match the start and end difference
  mutate(end = as.numeric(start) + nchar(`sequences`)) %>% 
  relocate(sequences) %>% relocate(name) %>% relocate(end) %>% relocate(start) %>% relocate(seqnames) %>% 
  pivot_longer(cols = !matches("seqnames") & !matches("start") & !matches("end") & !matches("name") & !matches("sequences"),
               names_to = 'trinuc32',
               values_to = 'removed_trinucs')
## free up memory
gc()
## continue
matched_tracks_granges = matched_tracks_granges %>%
  # include reverse complement for each trinuc32, so it's actually trinuc64
  mutate(trinuc32 = paste0(trinuc32, "|", reverseComplement(trinuc32, case="upper"))) %>% 
  ## do the removing thing
  as_tibble %>% 
  rowwise() %>% 
  mutate(new_end = list(rm_n_trinucs_at_random_indices(`removed_trinucs`, `trinuc32`, `sequences`)))
## free up memory
gc()
## continue
matched_tracks_granges = matched_tracks_granges %>%
  lazy_dt %>% 
  group_by(seqnames, start, end, name) %>% 
  # new_end is now unique positions 1bp before each trinucleotide removed
  summarise(new_end = list(unique(unlist(new_end)) + as.numeric(start) - 1)) %>% 
  ungroup
## free up memory
gc()
## continue
matched_tracks_granges = matched_tracks_granges %>%
  as_tibble %>% 
  unnest(new_end) %>% 
  lazy_dt %>% 
  distinct()
## free up memory
gc()
## continue
matched_tracks_granges = matched_tracks_granges %>%
  # in rows just before/at a region in which matching was not performed (typically "target" regions, as they are shorter), or at ending of chr, the new_end == end
  mutate(new_end = ifelse(is.na(new_end),
                          as.numeric(end),
                          new_end))
## free up memory
gc()
## continue
matched_tracks_granges = matched_tracks_granges %>%
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
  ungroup
## free up memory
gc()
## continue
matched_tracks_granges = matched_tracks_granges %>%
  # remove negative or 0 widths (due to overlapping of removed trinucleotides)
  filter(width >= 1) 
## free up memory
gc()
## continue
matched_tracks_granges = matched_tracks_granges %>%
  select(seqnames, new_start, new_end, width, strand, name)
## free up memory
gc()
## continue
matched_tracks_granges = matched_tracks_granges %>%
  rename("start" = "new_start",
         "end" = "new_end") %>%
  arrange(seqnames, start, end)
## free up memory
gc()
## continue
matched_tracks_granges = matched_tracks_granges %>%
  as_tibble() %>% 
  makeGRangesFromDataFrame(keep.extra.columns = T)
gc()

export.bed(matched_tracks_granges,
           paste0(filename, "__3ntMatched_euclidean-", round(euclidean_score, 4), "__int_", start_int, "-", end_int, ".part.bed"))
