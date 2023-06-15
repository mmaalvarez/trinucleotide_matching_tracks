#!/usr/bin/env nextflow

Channel
    .fromPath(params.file_paths)
    .splitCsv(header:false)
    .set{ file_paths }

process run_trinuc_matching {

    queue = 'normal_prio_long'
    time = { (params.maxTime + 3).hour }
    memory = { (params.memory_process1 + 5*(task.attempt-1)).GB }

    input:
    set file_path from file_paths
    val stoppingCriterion from params.stoppingCriterion
    val maxTime from params.maxTime
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
            Rscript $PWD/bin/1_run_trinuc_matching.R ${file_path} ${stoppingCriterion} ${maxTime} ${max_fraction_removed_trinucs} ${acceleration_score} ${euclidean_change_ratio} ${fast_progress_lfc_cutoff} ${progress_its} ${utils} ${good_mappability_regions}
        else
            # no conda environment named "R"
            Rscript $PWD/bin/1_run_trinuc_matching.R ${file_path} ${stoppingCriterion} ${maxTime} ${max_fraction_removed_trinucs} ${acceleration_score} ${euclidean_change_ratio} ${fast_progress_lfc_cutoff} ${progress_its} ${utils} ${good_mappability_regions}
        fi
    else
        # no conda installed
        Rscript $PWD/bin/1_run_trinuc_matching.R ${file_path} ${stoppingCriterion} ${maxTime} ${max_fraction_removed_trinucs} ${acceleration_score} ${euclidean_change_ratio} ${fast_progress_lfc_cutoff} ${progress_its} ${utils} ${good_mappability_regions}
    fi
    """
}

process apply_trinuc_matching_to_tracks {

    publishDir "$PWD/res/", mode: 'copy'

    queue = 'normal_prio_long'
    time = { 48.hour }
    memory = { (params.memory_process2 + 5*(task.attempt-1)).GB }

    input:
    file full_tracks_trinuc32_freq from full_tracks_trinuc32_freq
    file matched_tracks from matched_tracks
    val euclidean_score from euclidean_score
    file sequences from sequences
    file filename from filename
    path utils from params.utils

    output:
    file '*__3ntMatched_euclidean-*.bed.gz'

    """
    #!/usr/bin/env bash

    if command -v conda &> /dev/null
    then
        if conda env list | grep "^R " >/dev/null 2>/dev/null
        then
            # there is a conda environment named "R"
            conda activate R
            Rscript $PWD/bin/2_apply_trinuc_matching_to_tracks.R ${full_tracks_trinuc32_freq} ${matched_tracks} ${euclidean_score} ${sequences} ${filename} ${utils}
        else
            # no conda environment named "R"
            Rscript $PWD/bin/2_apply_trinuc_matching_to_tracks.R ${full_tracks_trinuc32_freq} ${matched_tracks} ${euclidean_score} ${sequences} ${filename} ${utils}
        fi
    else
        # no conda installed
        Rscript $PWD/bin/2_apply_trinuc_matching_to_tracks.R ${full_tracks_trinuc32_freq} ${matched_tracks} ${euclidean_score} ${sequences} ${filename} ${utils}
    fi
    """
}
