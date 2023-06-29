library(tidyverse)
library(data.table)
library(dtplyr)
library(conflicted)
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("between", "dplyr")


###########################################################################################################################
### match trinucleotide(32) proportions between genomic coordinates
# takes trinucleotide(32) count dataframe with N rows (1 per genomic coordinate, e.g. chr1_1_237572_background, chr1_237573_237942_target...) × 32 trinuc type columns
# removes trinucs(32) counts, so that trinucs(32) proportions are "matched" across all bins
# modified from data/users/malvarez/projects/RepDefSig/resources/trinuc_matching/marina/sampler_fran.R

# normalize each row to sum to 1
rowNorm = function(m){
  m / rowSums(m)
}

# euclidean distance
euclidean = function(a, b){
  sqrt(sum((a-b)^2))
} 
 
trinuc_matching = function(full_tracks_trinuc32_freq, 
                           stoppingCriterion = 0.01, # desired Euclidean score (max. overall distance between any bin's trinuc frequencies and all-bin-average trinuc frequencies)
                           maxIter = 20000*length(full_tracks_trinuc32_freq), # to prevent endless loops
                           maxTime = 8, # max time (default 8)
                           unitsTime = "hours", # time units, default "hours", if >24h should specify "days"
                           max_fraction_removed_trinucs = 0.5, # don't allow to remove more total trinucleotide counts than this fraction of the total original trinucleotide counts (default 1)
                           acceleration_score = 1, # multiplied to the n of counts to be removed at each iteration (default 1)
                           euclidean_change_ratio_range = c(0.1, 1.1), # if the ratio of the current euclidean score compared to the previous iteration's is between the lower term of the range and 1 (i.e. too slow), increase acceleration_score proportionally; decrease accel.._score accordingly (min. 1) if the ratio is larger than the range (i.e. euc. score changes too erratically)
                           fast_progress_lfc_cutoff = -0.00001, # minimum degree of Euclidean score decrease (LFC; e.g. log2(0.0615/0.062)) allowed for the mean of progress_its
                           progress_its = 1000){ # n last iterations (that reached a mineuclidean_score) used to calculate the progress
  ## initialize constants/variables
  euclidean_score = sqrt(2) # we want to make this decrease until 'stoppingCriterion'
  mineuclidean_score = euclidean_score
  progress = c() # to keep track of how the euclidean_score decrease is progressing (and early-stop if it shows signs of stagnation)
  bin_names = data.frame("bin" = rownames(full_tracks_trinuc32_freq)) # store for re-adding it in the end, since data.table format loses it
  counts = as.data.table(full_tracks_trinuc32_freq)
  min_counts = min(counts)
  counts_mineuclidean = counts
  n_bins = length(rownames(counts))
  removed_trinucs = 0
  rowSums_counts = rowSums(counts)
  total_orig_trinucs = sum(rowSums_counts)
  max_removed_trinucs = total_orig_trinucs * max_fraction_removed_trinucs
  iter = 0
  
  ## detect which bins have too few trinucs, to ignore them for offender search (i.e. less than the total current sum of trinucs (across rows AND columns) divided by the number of bins)
  too_few_trinuc_bins = data.frame(rowSums_counts) %>% 
    rownames_to_column("bin") %>% 
    filter(`rowSums_counts` <= (total_orig_trinucs / n_bins / 2)) %>% # WARNING: dividing by 2 so that not too many bins are excluded
    pull(bin)
  ## ... and check whether every row has been declared as unusable (would result in an emtpy 'offender')
  if(n_bins - length(too_few_trinuc_bins) <= 0){
    stop(sprintf("Cannot start optimization - 'too_few_trinuc_bins' comprises all possible bins, so no offender bin can be defined -- Exiting..."))
  }
  
  start_time = Sys.time()

  while ( TRUE ) {
    
    iter = iter + 1
    
    ## check whether we have reached max. num. iterations
    if (iter == maxIter) {
      counts = counts_mineuclidean
      cat( sprintf("Stopping optimization - maximum number of iterations reached - Returning min Euclidean score results (%f)\n", mineuclidean_score) );
      break
    }  
    
    ## check whether we have already removed too many trinucs
    if(removed_trinucs > max_removed_trinucs){
      counts = counts_mineuclidean
      cat(sprintf("Stopping optimization at iter %i/%i - %.02f%% of the %f original trinucs have already been removed -- Exiting and returning min Euclidean score results (%f)\nAnalysis terminated after %s\n", iter, maxIter, removed_trinucs/total_orig_trinucs*100, total_orig_trinucs, mineuclidean_score, paste(round(Sys.time() - start_time, 2), units(Sys.time() - start_time))))
      break
    }

    ###################################################
    ## get OFFENDER ROW WITH CORRECTABLE COLUMN
    
    # a data.frame version of the counts data.table, with a "bin" column containing the rownames
    counts_df = as.data.frame(counts) %>% rownames_to_column("bin")
    
    ## all rows are used for the "baseline mean frequencies"
    meanFreqs = colMeans(rowNorm(counts), na.rm = T)
    # ...but the 'too_few_trinuc_bins' are excluded from being candidates for offender...
    freqs_not_too_few_trinuc_bins = counts_df %>% 
      ## ...so filter out bins (rows) with too few trinucs
      filter(! bin %in% too_few_trinuc_bins) %>% 
      column_to_rownames("bin") %>% 
      rowNorm()
    
    ## LOOP HERE UNTIL WE HAVE AN OFFENDER ROW WITH CORRECTABLE COLUMN
    
    got_offender_row_corr_col = F
    
    while (got_offender_row_corr_col == F){
      
      # find the 'offender' row, which is the row most different from the meanFreqs vector
      offender_name = which.max(apply(freqs_not_too_few_trinuc_bins, 1, function(x){ euclidean(x, meanFreqs) })) %>% # note that the meanFreqs DOES use ALL rows
        names()
      offender_row = counts_df %>% 
        filter(bin == offender_name) %>% 
        column_to_rownames("bin")
      offender_row_freqs = rowNorm(offender_row)

      diffs = data.frame(meanFreqs - offender_row_freqs)
      
      if(sum(offender_row) <= 0){
        counts = counts_mineuclidean
        stop(sprintf("Cannot continue optimization at iter %i/%i - no more nts at the offender bin '%s' -- Exiting and returning min Euclidean score results (%f)\nAnalysis terminated after %s\n", iter, maxIter, offender_name, mineuclidean_score, paste(round(Sys.time() - start_time, 2), units(Sys.time() - start_time))))
      }
  
      ## in that row, find the column which is most responsible for the difference; however importantly we care ONLY about the negative differences in this vector! i.e. those are the cases where the offending row has HIGHER freqs (meaning we can correct that by removing sites... we can't add sites!!)
      
      is.correctableCol.zero = T
      
      while(is.correctableCol.zero == T){
        
        correctableColIndex = which.min(diffs)
        correctableCol = names(correctableColIndex)
        correctableColCounts = select(offender_row, all_of(correctableCol))
  
        if(correctableColCounts <= 0) {
          cat( sprintf("WARNING: At iter %i/%i, counts exhausted at bin '%s's 'correctableCol' (trinuc %s), trying next `which.min(diffs)` trinuc - Euclidean score: %f\n%s have passed\n", iter, maxIter, offender_name, correctableCol, as.numeric(euclidean_score), paste(round(Sys.time() - start_time, 2), units(Sys.time() - start_time))) )
          
          diffs = mutate_at(diffs, 
                            vars(all_of(correctableCol)), 
                            ~gsub(".*", NA, .))
          
        } else {
          is.correctableCol.zero = F
          
          # check if we have removed all trinuc from diffs (so we should move on to next row to be used as offender)
          if(length(diffs) == 0){
  
            # don't use this bin anymore as offender
            too_few_trinuc_bins = c(too_few_trinuc_bins, offender_name)
            
            cat(sprintf("WARNING: At iter %i/%i, removed all trinuc from 'diffs'; Bin '%s' is not used anymore as 'offender' bin - Euclidean score: %f\n%s have passed\n", iter, maxIter, offender_name, as.numeric(euclidean_score), paste(round(Sys.time() - start_time, 2), units(Sys.time() - start_time))))
          
          } else {
            # end loop
            got_offender_row_corr_col = T
          }
        }
      }
  
      if(n_bins - length(too_few_trinuc_bins) <= 0){
        counts = counts_mineuclidean
        stop(sprintf("Cannot continue optimization at iter %i/%i - 'too_few_trinuc_bins' comprises all possible bins, so no offender bin can be defined -- Exiting and returning min Euclidean score results (%f)\nAnalysis terminated after %s\n", iter, maxIter, mineuclidean_score, paste(round(Sys.time() - start_time, 2), units(Sys.time() - start_time))))
      }
    }
    
    # store euclidean_score from previous round
    previous_euclidean_score = euclidean_score
    
    ## calculate euclidean score (how diff. are the offender row freqs. from the mean freqs. across the table)
    euclidean_score = euclidean(offender_row_freqs, meanFreqs)
    
    # calculate ratio of change in this iteration
    euclidean_change_ratio = euclidean_score / previous_euclidean_score
    
    # if euclidean_score is new minimum:
    if(euclidean_score < mineuclidean_score){
      
      # 1- store counts table
      mineuclidean_score = euclidean_score
      counts_mineuclidean = counts
      
      # 2- calculate amount of change between previous and current iteration (since it reached a mineuclidean_score)...
      LFC_euc_score = log2(euclidean_change_ratio)
      
      # ...and append it to the 'amount_change's of the last 'progress_its' iterations (that reached a mineuclidean_score)
      progress = c(tail(progress, progress_its-1),
                   LFC_euc_score)
      mean_progress_lfc = mean(progress)
    }
    
    # did we reduce the euclidean_score enough?
    if (euclidean_score <= stoppingCriterion) {
      cat(sprintf("Successfully completed optimization: Euclidean score (%f) lower than %f - Returning current results\n", euclidean_score, stoppingCriterion) )
      break
    } else if (mean_progress_lfc>fast_progress_lfc_cutoff  &  length(progress)==progress_its) {
      # if not, is the euclidean_score stagnated? (i.e. very small/slow progress)
      cat(sprintf("Cannot continue optimization at iter %i/%i - progress (mean LFC) of last %i iterations (that decreased the min. Euclidean score) is too slow: %f > %f (fast_progress_lfc_cutoff) -- Exiting and returning min Euclidean score results (%f)\nAnalysis terminated after %s\n", iter, maxIter, progress_its, mean_progress_lfc, fast_progress_lfc_cutoff, mineuclidean_score, paste(round(Sys.time() - start_time, 2), units(Sys.time() - start_time))))
      break
    }
    
    # if the ratio of the current euclidean score compared to the previous iteration's is between the lower term of the range and 1 (i.e. too slow), increase acceleration_score proportionally; decrease accel.._score accordingly (min. 1) if the ratio is larger than the range (i.e. euc. score changes too erratically)
    if (between(euclidean_change_ratio, 
                min(euclidean_change_ratio_range), 1)){

      # too slow: increase acceleration_score proportionally
      acceleration_score = acceleration_score + (1 - euclidean_change_ratio)
      
    } else if (euclidean_change_ratio > max(euclidean_change_ratio_range)){
      
      # "excess" of euclidean_change_ratio
      excess_euclidean_change_ratio = euclidean_change_ratio / max(euclidean_change_ratio_range)
        
      # too fast/erratic: decrease acceleration_score proportionally to excess (min 1)
      acceleration_score = acceleration_score / excess_euclidean_change_ratio
      
      # ensure that acceleration_score is not lower than 1
      if (acceleration_score < 1){
        acceleration_score = 1
      }
    }

    ###############################
    ### subtract counts from the responsible trinuc (column) of the offender bin (row) so that the column frequencies of the latter get closer to the mean column frequencies (across the whole table)
    
    ### note that this adjustment (subtraction) is too conservative, but by iterating it should converge to the right value
    # the acceleration_score may speed it up a bit
    subtractThis = round(as.numeric(diffs[correctableColIndex] * sum(offender_row)) * acceleration_score)
    
    # total counts currently at the offender_name row × correctableCol trinuc intersection
    counts_at_correctableCol_offender_name = counts[[correctableCol]][as.integer(offender_name)]
    
    # total counts that will remain at the offender_name row × correctableCol trinuc intersection after the removal of 'subtractThis' counts
    remaining_counts_at_correctableCol_offender_name = counts_at_correctableCol_offender_name + subtractThis
    
    # make sure that, possibly due to the acceleration score, there were not more counts removed than the total available
    if(remaining_counts_at_correctableCol_offender_name < 0){
      
      # if so, let there be min_counts (probably 0) left at that row×col intersection
      remaining_counts_at_correctableCol_offender_name = min_counts
        
      # update subtractThis just to keep track of the actual # of counts removed
      subtractThis = min_counts - counts_at_correctableCol_offender_name
    }
    
    # keep track of # of counts removed
    removed_trinucs = removed_trinucs + abs(subtractThis)
    
    ### apply this change to the actual counts table
    set(counts, 
        i = as.integer(offender_name), 
        j = correctableCol,
        value = remaining_counts_at_correctableCol_offender_name)
    
    
    ## output log every 100th iter
    if(iter %% 100 == 0){

      ## first check whether we have reached max. time
      time_passed = str_split(paste(round(Sys.time() - start_time), units(Sys.time() - start_time)), " ")[[1]]
      if (as.numeric(time_passed[1]) >= maxTime & time_passed[2] == unitsTime){
        counts = counts_mineuclidean
        cat( sprintf("Stopping optimization - %s %s have passed - Returning min Euclidean score results (%f)\n", time_passed[1], time_passed[2], mineuclidean_score) )
        break
      }  
      
      ## also update which bins have too few trinucs, to ignore them for offender search (i.e. less than the total current sum of trinucs (across rows AND columns) divided by the number of bins)
      rowSums_counts = rowSums(counts)
      too_few_trinuc_bins = data.frame(rowSums_counts) %>% 
        rownames_to_column("bin") %>% 
        filter(`rowSums_counts` <= (sum(rowSums_counts) / n_bins / 2)) %>% # WARNING: dividing by 2 so that not too many bins are excluded
        pull(bin)
      ## ... and check whether every row has been declared as unusable (would result in an emtpy 'offender')
      if(n_bins - length(too_few_trinuc_bins) <= 0){
        counts = counts_mineuclidean
        stop(sprintf("Cannot continue optimization at iter %i/%i - 'too_few_trinuc_bins' comprises all possible bins, so no offender bin can be defined -- Exiting and returning min Euclidean score results (%f)\nAnalysis terminated after %s\n", iter, maxIter, mineuclidean_score, paste(round(Sys.time() - start_time, 2), units(Sys.time() - start_time))))
      }
      
      # print log message
      cat( sprintf("Iteration %i/%i:\n\tSubtracted %i '%s's at bin #'%s'\n\t%.02f%% of the original trinucleotides have been removed\n\tEuclidean score: %f\n\tAcceleration score: %f\n\tMean progress LFC: %f for the last %i iterations that updated the min. Euclidean score\n\t%s %s have passed\n\n", iter, maxIter, abs(as.numeric(subtractThis)), correctableCol, offender_name, removed_trinucs/total_orig_trinucs*100, euclidean_score, acceleration_score, mean_progress_lfc, progress_its, time_passed[1], time_passed[2]))
    }
    
  } ## keep iterating...
  
  ## return final counts_mineuclidean tables, and its mineuclidean_score
  return(list(as_tibble(counts) %>%
                # put back the bin names
                bind_cols(bin_names) %>% 
                relocate("bin"),
              mineuclidean_score))
}



###########################################################################################################################
### loaded by 2_apply_trinuc_matching_to_tracks, to randomly remove as many trinucleotides from each type as have been removed with trinuc_matching, and keep track of their positions

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



#######################################################################################
## load tracks file, keep metadata without NAs
load_tracks <- function(path_file) {
  
  # load it
  tracks_file = tryCatch(import.bw(path_file),
                         error = function(e) tryCatch(import.bedGraph(path_file),
                                                      error = function(e) tryCatch(import.bed(path_file),
                                                                                   error = function(e) tryCatch(makeGRangesFromDataFrame(read_tsv(path_file), keep.extra.columns = T),
                                                                                                                error = function(e) import.wig(path_file)))))
  # keep only metadata without NAs
  for(metadata_col in names(elementMetadata(tracks_file))){
    if(unique(is.na(unique(elementMetadata(tracks_file)[[metadata_col]])))){
      mcols(tracks_file)[[metadata_col]] <- NULL
    }
  }
  # call the metadata column "name"
  names(mcols(tracks_file)) = "name"
  
  # sort tracks file by chromosome
  tracks_file = sortSeqlevels(tracks_file)
  gc()
  
  return(tracks_file)
}