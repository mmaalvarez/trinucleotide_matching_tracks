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


## load tracks

file_path = ifelse(interactive(),
                  yes = "../input/example.bed",
                  no = args[1])
filename = gsub(".*\\/", "", file_path) %>% gsub("\\..*", "", .)

# try different extensions
full_tracks = tryCatch(import.bw(file_path), # could also be .bigWig
                       error = function(e) tryCatch(import.bedGraph(file_path),
                                                    error = function(e) tryCatch(import.bed(file_path), # also accepts .bed.gz
                                                                                 error = function(e) makeGRangesFromDataFrame(read_tsv(file_path),
                                                                                                                              keep.extra.columns = T))))
colnames(elementMetadata(full_tracks)) = "name"
gc()


## parameter values for trinuc_matching()

stoppingCriterion = ifelse(interactive(),
                           yes = "0.01",
                           no = args[2]) %>% 
  as.numeric()

maxTime = ifelse(interactive(),
                 yes = "8",
                 no = args[3]) %>% 
  as.numeric()

max_fraction_removed_trinucs = ifelse(interactive(),
                                      yes = "0.5",
                                      no = args[4]) %>% 
  as.numeric()

acceleration_score = ifelse(interactive(),
                            yes = "1",
                            no = args[5]) %>% 
  as.numeric()

euclidean_change_ratio = ifelse(interactive(),
                                yes = "0.1,1.1",
                                no = args[6]) %>% 
  strsplit(., split=",", fixed = T) %>% 
  magrittr::extract2(1) %>% 
  as.numeric()

## source of trinuc_matching() function
trinuc_matching_source = ifelse(interactive(),
                                yes = "utils.R",
                                no = args[7])

## keep SNVs in good mappability regions
good_mappability_regions = ifelse(interactive(),
                                  yes = "../crg75/CRG75_nochr.bed",
                                  no = args[8]) %>%
  import.bed() %>% data.frame %>%
  rename("chrom" = "seqnames") %>% 
  mutate(chrom = gsub("^", "chr", chrom))
gc()

full_tracks_good_mappability = full_tracks %>% 
  data.frame %>% rename("chrom" = "seqnames") %>% 
  bed_intersect(good_mappability_regions, suffix = c("_full_tracks", "_crg75")) %>% 
  ## parse the start and end of the intersects
  group_by(chrom) %>% 
  # important to sort because we will use lag() and lead()
  arrange(start_crg75) %>% 
  mutate_at(vars(start_crg75), 
            # to deal with good_mappability_regions that overlap >=2 adjacent target-background-.. (or background-target-..) coordinates
            funs(dist_with_next = . - lead(.),
                 dist_with_previous = . - lag(.))) %>%
                           # if 1st region of the overlap, its 'start' is either the region's start or the good_mappability's start, whichever is the largest
  mutate(start = case_when(dist_with_next == 0 & (dist_with_previous != 0 | is.na(dist_with_previous)) ~ pmax(start_full_tracks, start_crg75),
                           # elif middle/last region of the overlap, its 'start' is the region's start
                           dist_with_previous == 0 & !is.na(dist_with_previous) ~ start_full_tracks,
                           # elif no overlap between regions, BUT the good_mappability's start is BEFORE the region's start, its 'start' is the region's
                           start_crg75 < start_full_tracks ~ start_full_tracks,
                           # else, keep the good_mappability's start
                           TRUE ~ start_crg75),
                         # if 1st region of the overlap, its 'end' is the region's end
         end = case_when(dist_with_next == 0 & (dist_with_previous != 0 | is.na(dist_with_previous)) ~ end_full_tracks,
                         # elif middle/last region of the overlap, its 'end' is either the region's end or the good_mappability's end, whichever is the smallest
                         dist_with_previous == 0 & !is.na(dist_with_previous) ~ pmin(end_full_tracks, end_crg75),
                         # elif no overlap between regions, BUT the good_mappability's end is AFTER the region's end, its 'end' is the region's
                         end_crg75 > end_full_tracks ~ end_full_tracks,
                         # else, keep the good_mappability's end
                         TRUE ~ end_crg75)) %>% 
  ungroup %>% 
  mutate(seqnames = factor(chrom, ordered = T, levels = paste0("chr", c(seq(1,22), "X", "Y")))) %>% 
  arrange(seqnames, start) %>% 
  select(seqnames, start, end, name_full_tracks) %>% 
  rename_all(~str_replace_all(., "_crg75|_full_tracks", "")) %>%
  makeGRangesFromDataFrame(keep.extra.columns = T)
gc()


## get trinuc frequencies from each coordinate in tracks

# extract the sequences
sequences = getSeq(BSgenome.Hsapiens.UCSC.hg19, ## WARNING: Assuming that coordinates are in hg19 -- if hg38, change to 'BSgenome.Hsapiens.UCSC.hg38'
                   names = full_tracks_good_mappability)
gc()

# get frequency of each trinucleotide per sequence, moving 1 nucleotide downstream each time
trinuc32_freq = trinucleotideFrequency(sequences) %>%
  as_tibble %>%
  rownames_to_column("id") %>% 
  lazy_dt %>%
  pivot_longer(cols = -id,
               names_to = 'trinuc32',
               values_to = "freq") %>%
  group_by(id) %>%
  ## collapse into 32 trinucleotides
  # trinucs that do not have a C or T in center, convert to reverse complement
  mutate(trinuc32 = ifelse(substr(trinuc32, start = 2, stop = 2) %in% c('A', 'G'),
                           reverseComplement(trinuc32, case="upper"),
                           trinuc32)) %>% 
  # sum up frequencies of each N(C|T)N & reverse complement pair, within id
  group_by(id, trinuc32) %>%
  summarise(freq = sum(freq)) %>%
  ungroup %>% 
  # back to original format (each row ('id') maps to the same row in map_features)
  pivot_wider(names_from = 'trinuc32',
              values_from = 'freq') %>% 
  as_tibble %>% 
  arrange(as.numeric(id)) %>% 
  select(-id)
gc()

# bind trinuc32 freqs to original tracks (good mappability)
full_tracks_trinuc32_freq = full_tracks_good_mappability %>%
  data.frame %>%
  unite("bin", seqnames, start, end, name) %>% 
  select(bin) %>% 
  bind_cols(trinuc32_freq) %>%
  column_to_rownames("bin")
gc()


## run matching
source(trinuc_matching_source)

matched_tracks = trinuc_matching(full_tracks_trinuc32_freq, 
                                 stoppingCriterion = stoppingCriterion, # desired Euclidean score (max. overall distance between any bin's trinuc frequencies and all-bin-average trinuc frequencies)
                                 maxIter = 20000*length(full_tracks_trinuc32_freq), # to prevent endless loops
                                 maxTime = maxTime, # max hours (default 8)
                                 max_fraction_removed_trinucs = max_fraction_removed_trinucs, # don't allow to remove more total trinucleotide counts than this fraction of the total original trinucleotide counts (default 0.5)
                                 acceleration_score = acceleration_score, # multiplied to the n of counts to be removed at each iteration (default 1)
                                 euclidean_change_ratio = euclidean_change_ratio) # if the ratio of the current euclidean score compared to the previous iteration's is between the lower term of the range and 1 (i.e. too slow), increase acceleration_score proportionally; decrease accel.._score accordingly (min. 1) if the ratio is larger than the range (i.e. euc. score changes too erratically)                                 
gc()

euclidean_score = matched_tracks[2][[1]]

matched_tracks = matched_tracks[1][[1]]


## pass files to 2nd process
rownames_to_column(full_tracks_trinuc32_freq, "bin") %>% 
  write_tsv(paste0(filename,"_full_tracks_trinuc32_freq.tsv"))

write_tsv(matched_tracks, paste0(filename,"_matched_tracks.tsv"))

data.frame(euclidean_score) %>%
  write_tsv(paste0(filename,"_euclidean_score.tsv"))

data.frame(sequences) %>% 
  write_tsv(paste0(filename,"_sequences.tsv"))

data.frame(filename) %>%
  write_tsv(paste0(filename,"_filename.tsv"))
