#!/usr/bin/env nextflow

Channel
    .fromPath(params.file_paths)
    .splitCsv(header:false)
    .set{ file_paths }


process run_trinuc_matching {

    queue = 'normal_prio_long'
    time = { params.time_process1.hour }
    memory = { (params.memory_process1 + 5*(task.attempt-1)).GB }

    input:
    set file_path from file_paths
    val stoppingCriterion from params.stoppingCriterion
    val maxTime from params.maxTime
    val unitsTime from params.unitsTime
    val max_fraction_removed_trinucs from params.max_fraction_removed_trinucs
    val acceleration_score from params.acceleration_score
    val euclidean_change_ratio from params.euclidean_change_ratio
    val fast_progress_lfc_cutoff from params.fast_progress_lfc_cutoff
    val progress_its from params.progress_its
    path utils from params.utils
    val good_mappability_regions from params.good_mappability_regions

    output:
    file '*_full_tracks_trinuc32_freq.tsv' into full_tracks_trinuc32_freq
    file '*_matched_tracks.tsv' into matched_tracks
    file '*_euclidean_score.tsv' into euclidean_score
    file '*_sequences.tsv' into sequences
    file '*_filename.tsv' into filename

    """
    #!/usr/bin/env bash

    if command -v conda &> /dev/null
    then
        if conda env list | grep "^R " >/dev/null 2>/dev/null
        then
            # there is a conda environment named "R"
            conda activate R
            Rscript $PWD/bin/1_run_trinuc_matching.R ${file_path} ${stoppingCriterion} ${maxTime} ${unitsTime} ${max_fraction_removed_trinucs} ${acceleration_score} ${euclidean_change_ratio} ${fast_progress_lfc_cutoff} ${progress_its} ${utils} ${good_mappability_regions}
        else
            # no conda environment named "R"
            Rscript $PWD/bin/1_run_trinuc_matching.R ${file_path} ${stoppingCriterion} ${maxTime} ${unitsTime} ${max_fraction_removed_trinucs} ${acceleration_score} ${euclidean_change_ratio} ${fast_progress_lfc_cutoff} ${progress_its} ${utils} ${good_mappability_regions}
        fi
    else
        # no conda installed
        Rscript $PWD/bin/1_run_trinuc_matching.R ${file_path} ${stoppingCriterion} ${maxTime} ${unitsTime} ${max_fraction_removed_trinucs} ${acceleration_score} ${euclidean_change_ratio} ${fast_progress_lfc_cutoff} ${progress_its} ${utils} ${good_mappability_regions}
    fi
    """
}



// create a channel that emits a sequence of intervals from 1 to the number of rows (ranges) in the input file, with each interval being of length a 100th part of the original (i.e. rows, ranges)

def createRange(start, end, step) {
    def range = []
    for(def i = start; i <= end; i+=step) {
        range << i
    }
    return range
}

Channel
    .from(createRange(1, params.nranges, params.fraction_length))
    .map { it -> [it, (it+(params.fraction_length-1) > params.nranges) ? params.nranges : it+(params.fraction_length-1)] }
    .set{ input_fraction }

process apply_trinuc_matching_to_tracks {

    time = { params.time_process2.hour }
    memory = { (params.memory_process2 + 5*(task.attempt-1)).GB }

    input:
    set start_int,end_int from input_fraction
    file full_tracks_trinuc32_freq from full_tracks_trinuc32_freq
    file matched_tracks from matched_tracks
    val euclidean_score from euclidean_score
    file sequences from sequences
    file filename from filename
    path utils from params.utils

    output:
    file '*__3ntMatched_euclidean-*.bed.gz' into res

    """
    #!/usr/bin/env bash

    if command -v conda &> /dev/null
    then
        if conda env list | grep "^R " >/dev/null 2>/dev/null
        then
            # there is a conda environment named "R"
            conda activate R
            Rscript $PWD/bin/2_apply_trinuc_matching_to_tracks.R ${start_int} ${end_int} ${full_tracks_trinuc32_freq} ${matched_tracks} ${euclidean_score} ${sequences} ${filename} ${utils}
        else
            # no conda environment named "R"
            Rscript $PWD/bin/2_apply_trinuc_matching_to_tracks.R ${start_int} ${end_int} ${full_tracks_trinuc32_freq} ${matched_tracks} ${euclidean_score} ${sequences} ${filename} ${utils}
        fi
    else
        # no conda installed
        Rscript $PWD/bin/2_apply_trinuc_matching_to_tracks.R ${start_int} ${end_int} ${full_tracks_trinuc32_freq} ${matched_tracks} ${euclidean_score} ${sequences} ${filename} ${utils}
    fi
    """
}


// CORRECT THIS so that ${filename} and ${euclidean_score} are the actual parameter values, not literary '${filename}' and '${euclidean_score}'
filename = params.filename
euclidean_score = params.euclidean_score

// ALSO look up if there is a way to sort collected files after collection, based on 1st (chr) and 2nd (start) columns, otherwise the 100 subtables are concatenated at random
res
    .collectFile(name: 'res/${filename}__3ntMatched_euclidean-${euclidean_score}.bed')
    .println { "Finished, results saved in res/" }
