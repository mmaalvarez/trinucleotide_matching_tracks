library(tidyverse)
library(data.table)
library(dtplyr)
library(GenomicRanges)
library(rtracklayer)
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

filename = ifelse(interactive(),
                  yes = "/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/TCGA_PCAWG_Hartwig_CPTAC_POG_MMRFCOMMPASS/SHM/CTCF_cohesin/1_CTCF_cohesin_peaks_coords/CTCF_cohesin_peaks_coords",
                  no = args[1])
# try different formats
full_tracks = tryCatch(import.bw(paste0(filename, ".bw")), # could also be .bigWig
                       error = function(e) tryCatch(import.bedGraph(paste0(filename, ".bedGraph")),
                                                    error = function(e) tryCatch(import.bed(paste0(filename, ".bed")), # also accepts .bed.gz
                                                                                 error = function(e) makeGRangesFromDataFrame(read_tsv(paste0(filename, ".tsv")),
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
                                 yes = "/g/strcombio/fsupek_data/users/malvarez/projects/RepDefSig/bin/utils.R",
                                 no = args[7])

## keep SNVs in good mappability regions
good_mappability_regions = ifelse(interactive(),
                                 yes = "/g/strcombio/fsupek_home/mmunteanu/reference/CRG75_nochr.bed",
                                 no = args[8]) %>%
  import.bed() %>% data.frame %>%
  mutate(seqnames = gsub("^", "chr", seqnames))



## get trinuc frequencies from each coordinate in tracks

# extract the sequences
sequences = getSeq(BSgenome.Hsapiens.UCSC.hg19, ## WARNING: Assuming that coordinates are in hg19 -- if hg38, change to 'BSgenome.Hsapiens.UCSC.hg38'
                   names = full_tracks)

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

# bind trinuc32 freqs to original tracks
full_tracks_trinuc32_freq = full_tracks %>%
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
  write_tsv("full_tracks_trinuc32_freq.tsv")

write_tsv(matched_tracks, "matched_tracks.tsv")

data.frame(euclidean_score) %>%
  write_tsv("euclidean_score.tsv")

data.frame(sequences) %>% 
  write_tsv("sequences.tsv")
