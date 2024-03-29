library(tidyverse)
library(data.table)
library(dtplyr)
library(GenomicRanges)
library(rtracklayer)
library(valr)
library("BSgenome.Hsapiens.UCSC.hg19")
library(spgs)
library(rlang)
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
conflict_prefer("strsplit", "base")


args = commandArgs(trailingOnly=TRUE)


## tracks file path and name

file_path = ifelse(interactive(),
                   yes = '/g/strcombio/fsupek_cancer3/malvarez/chromatin_info/DNA_repair__protein_binding/DNA_repair/BER/OGG1/chipseq/bed_hg19/gaps_score0/hg19_GSM_CP-Sample_Flag-OGG1-GOX-60-1_gaps-to-score-0.bed', #"example.bed",
                   no = args[1])
filename = gsub(".*\\/", "", file_path) %>% gsub("\\..*", "", .)

## parameter values for trinuc_matching()

stoppingCriterion = ifelse(interactive(),
                           yes = "0.01",
                           no = args[2]) %>% 
  as.numeric()

maxTime = ifelse(interactive(),
                 yes = "8",
                 no = args[3]) %>% 
  as.numeric()

unitsTime = ifelse(interactive(),
                   yes = "hours",
                   no = args[4])

max_fraction_removed_trinucs = ifelse(interactive(),
                                      yes = "0.25",
                                      no = args[5]) %>% 
  as.numeric()

acceleration_score = ifelse(interactive(),
                            yes = "1",
                            no = args[6]) %>% 
  as.numeric()

euclidean_change_ratio = ifelse(interactive(),
                                yes = "0.1,1.1",
                                no = args[7]) %>% 
  strsplit(., split=",", fixed = T) %>% 
  magrittr::extract2(1) %>% 
  as.numeric()

fast_progress_lfc_cutoff = ifelse(interactive(),
                                  yes = "-0.00001",
                                  no = args[8]) %>%
  as.numeric()

progress_its = ifelse(interactive(),
                      yes = "1000",
                      no = args[9]) %>% 
  as.numeric()

## source of trinuc_matching() function
trinuc_matching_source = ifelse(interactive(),
                                yes = "utils.R",
                                no = args[10]) %>% 
  source()

## OPTIONAL: keep SNVs in good mappability regions
good_mappability_regions = ifelse(interactive(),
                                  yes = "",
                                  no = args[11])


##########################

### load tracks file

full_tracks = load_tracks(file_path)
gc()


## binarize continuous (i.e. numeric) scores

scores = elementMetadata(full_tracks)[[1]]

if(is.numeric(scores)){

  # calculate median score for the feature, for binarization
  median_score = median(scores)

  # check if "score" only contains numbers, even though it is labelled as a "character" column and these numbers are formatted as strings
  } else if(all(grepl("^[-]?[0-9]+(\\.[0-9]+)?$",
                      scores))){
    median_score = median(as.numeric(scores))

  }else{
  cat (sprintf("Tracks file metadata is not continuous (%s), keeping it unchanged\n", paste(unique(mcols(full_tracks))[["name"]], collapse = ", ")))
}

if(exists("median_score")){
  ## binarize weighted average feature value by being lower or larger than the across-genome median
  full_tracks = full_tracks %>%
    data.frame %>% 
    lazy_dt %>% 
    #### WARNING first do the average score at duplicated (start end) ranges, this is due to the (in some features) hg38-->hg19 lift dividing some ranges into 2 alternative ranges with the same score
    group_by(seqnames, start, end) %>% 
    summarise(name = mean(name)) %>% 
    ungroup %>% 
    as_tibble %>% 
    rowwise %>% 
    lazy_dt %>% 
    mutate(name = ifelse(name <= median_score,
                         "low",
                         "high")) %>% 
    as_tibble %>% 
    makeGRangesFromDataFrame(keep.extra.columns = T)
  gc()
  
  ## collapse contiguous ranges if they have same metadata levels
  full_tracks = unlist(reduce(split(full_tracks, ~name)))
  mcols(full_tracks) = names(full_tracks)
  full_tracks = full_tracks %>% 
    as_tibble %>% 
    filter(seqnames %in% paste0("chr", c(seq(1,22),"X","Y"))) %>% 
    mutate(seqnames = factor(seqnames, levels = paste0("chr", c(seq(1,22),"X","Y")))) %>% 
    arrange(seqnames, start) %>% 
    rename("name" = "X") %>% 
    makeGRangesFromDataFrame(keep.extra.columns = T)
  gc()
}


## OPTIONAL: keep SNVs in good mappability regions
if(!good_mappability_regions %in% c("", "None", "none", "NONE", "NULL", "Na", "NA")){
  
  # import good_mappability_regions
  import.bed() %>% data.frame %>%
  rename("chrom" = "seqnames") %>% 
  mutate(chrom = gsub("^", "chr", chrom))

  # filter out bad mappability regions
  full_tracks = full_tracks %>% 
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
}
gc()


## get trinuc frequencies from each coordinate in tracks

# extract the sequences
sequences = getSeq(BSgenome.Hsapiens.UCSC.hg19, ## WARNING: Assuming that coordinates are in hg19 -- if hg38, change to 'BSgenome.Hsapiens.UCSC.hg38'
                   names = full_tracks)
gc()

## collapse into 32 trinucleotides
# get frequency of each trinucleotide per sequence, moving 1 nucleotide downstream each time
trinuc32_freq = trinucleotideFrequency(sequences) %>%
  as_tibble

# Identify columns where the 2nd character is A or G
# trinucs that do not have a C or T in center, convert to reverse complement
cols_to_transform = grep("^.([AG]).$", colnames(trinuc32_freq))

# Apply the reverseComplement to those columns
transformed_cols = sapply(colnames(trinuc32_freq)[cols_to_transform], function(name) {
  as.character(Biostrings::reverseComplement(DNAString(name)))
})

# Replace those column names in the dataframe
colnames(trinuc32_freq)[cols_to_transform] <- transformed_cols

trinuc32_freq = data.frame(trinuc32_freq) %>% 
  rownames_to_column("id") %>% 
  #lazy_dt %>%
  pivot_longer(cols = -id,
               names_to = 'trinuc32',
               values_to = "freq") %>%
  # sum up frequencies of each N(C|T)N & reverse complement pair, within id
  mutate(trinuc32 = gsub("\\..", "", trinuc32)) %>% 
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
matched_tracks = trinuc_matching(full_tracks_trinuc32_freq, 
                                 stoppingCriterion = stoppingCriterion, # desired Euclidean score (max. overall distance between any bin's trinuc frequencies and all-bin-average trinuc frequencies)
                                 maxIter = 20000*length(full_tracks_trinuc32_freq), # to prevent endless loops
                                 maxTime = maxTime, # max time (default 8)
                                 unitsTime = unitsTime, # time units, default "hours", if >24h should specify "days"                                 max_fraction_removed_trinucs = max_fraction_removed_trinucs, # don't allow to remove more total trinucleotide counts than this fraction of the total original trinucleotide counts (default 0.5)
                                 acceleration_score = acceleration_score, # multiplied to the n of counts to be removed at each iteration (default 1)
                                 euclidean_change_ratio = euclidean_change_ratio, # if the ratio of the current euclidean score compared to the previous iteration's is between the lower term of the range and 1 (i.e. too slow), increase acceleration_score proportionally; decrease accel.._score accordingly (min. 1) if the ratio is larger than the range (i.e. euc. score changes too erratically)
                                 fast_progress_lfc_cutoff = fast_progress_lfc_cutoff, # minimum degree of Euclidean score decrease (LFC; e.g. log2(0.0615/0.062)) allowed for the mean of progress_its
                                 progress_its = progress_its) # n last iterations (that reached a mineuclidean_score) used to calculate the progress
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
