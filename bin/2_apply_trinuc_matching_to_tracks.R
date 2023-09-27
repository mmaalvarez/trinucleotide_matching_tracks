library(tidyverse)
library(data.table)
library(dtplyr)
library(GenomicRanges)
library(rtracklayer)
library(spgs)
library(Biostrings)
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


## load files from 1st process

interactive_work_path = "../work/0d/8150eb171672e746c785dfb41d2af8/"

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

euclidean_score = ifelse(interactive(),
                         yes = Sys.glob(paste0(interactive_work_path, "*_euclidean_score.tsv")),
                         no = args[5]) %>% 
  read_tsv() %>% pull(euclidean_score) %>% as.numeric()

sequences = ifelse(interactive(),
                   yes = Sys.glob(paste0(interactive_work_path, "*_sequences.tsv")),
                   no = args[6]) %>% 
  read_tsv() %>%
  ### NEW fraction of the input files to parse in this iteration
  slice(start_int:end_int)

filename = ifelse(interactive(),
                   yes = Sys.glob(paste0(interactive_work_path, "*_filename.tsv")),
                   no = args[7]) %>% 
  read_tsv() %>% pull(filename)

# source of rm_n_trinucs_at_random_indices() function
trinuc_matching_source = ifelse(interactive(),
                                yes = "utils.R",
                                no = args[8]) %>% 
  source()

gc()


## within each coordinate, randomly remove as many trinucleotides from each type as have been removed with trinuc_matching, and keep track of their positions

# subtract orig from remaining counts
matched_tracks_granges = full_tracks_trinuc32_freq %>%
  column_to_rownames("bin") %>% 
  map2_dfc(matched_tracks, ~ .x - .y)
## free up memory
rm(full_tracks_trinuc32_freq) ; gc()
## continue
matched_tracks_granges = matched_tracks_granges %>% 
  # bin names back, and sequences
  bind_cols(bin_names) %>% relocate(bin) %>%  
  separate(bin, into = c("seqnames", "start", "end", "name"), extra = "merge") %>% 
  bind_cols(sequences)
## free up memory
rm(bin_names) ; rm(sequences) ; gc()

### include reverse complement for each trinuc32, so it's actually trinuc64
cols_to_transform = match(colnames(matched_tracks), colnames(matched_tracks_granges))
# Apply the reverseComplement to those columns
transformed_cols = sapply(colnames(matched_tracks_granges)[cols_to_transform], 
                          function(name) {
                            paste0(name, "|", as.character(reverseComplement(DNAString(name))))
                            })
# Replace those column names in the dataframe
colnames(matched_tracks_granges)[cols_to_transform] <- transformed_cols

### WARNING: updating the 'end' since the actual length of some sequences (obtained with BSGenome::getSeq()) do not match the start and end difference
matched_tracks_granges = matched_tracks_granges %>% 
  mutate(end = as.numeric(start) + nchar(`sequences`)) %>% 
  relocate(sequences) %>% relocate(name) %>% relocate(end) %>% relocate(start) %>% relocate(seqnames) %>% 
  pivot_longer(cols = !matches("seqnames") & !matches("start") & !matches("end") & !matches("name") & !matches("sequences"),
               names_to = 'trinuc64',
               values_to = 'removed_trinucs')
gc()


# keep separate for now the rows without removed trinucs
matched_tracks_granges_no_rm_trinucs = matched_tracks_granges %>% 
  filter(removed_trinucs == 0) %>% 
  select(seqnames, start, end, name) %>% 
  mutate(new_end = end)
gc()

## do the removing thing
matched_tracks_granges = matched_tracks_granges %>%
  # only do it in the ones that actually need to have trinucs removed
  filter(removed_trinucs != 0) %>% #slice_sample(n=2) %>% 
  rowwise() %>% 
  mutate(new_end = list(rm_n_trinucs_at_random_indices(`removed_trinucs`, `trinuc64`, `sequences`)))
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
  distinct()  %>% 
  ## add back the rows without removed trinucs
  bind_rows(matched_tracks_granges_no_rm_trinucs) %>% 
  mutate(seqnames = factor(seqnames, levels = paste0("chr", c(seq(1:22), "X", "Y")))) %>% 
  arrange(seqnames, start, new_end) %>%
  group_by(seqnames, start, end, name) %>% 
  lazy_dt 
## free up memory
gc()
## continue
matched_tracks_granges = matched_tracks_granges %>%
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
