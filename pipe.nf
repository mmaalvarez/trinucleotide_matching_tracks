#!/usr/bin/env nextflow

process run_trinuc_matching {

    time = { (params.maxTime + 0.25).hour }
    memory = { (params.memory_process1 + 5*(task.attempt-1)).GB }

    input:
    path filename from params.filename
    val stoppingCriterion from params.stoppingCriterion
    val maxTime from params.maxTime
    val max_fraction_removed_trinucs from params.max_fraction_removed_trinucs
    val acceleration_score from params.acceleration_score
    val euclidean_change_ratio from params.euclidean_change_ratio
    path utils from params.utils
    path good_mappability_regions from params.good_mappability_regions

    output:
    file 'full_tracks_trinuc32_freq.tsv' into full_tracks_trinuc32_freq
    file 'matched_tracks.tsv' into matched_tracks
    file 'euclidean_score.tsv' into euclidean_score
    file 'sequences.tsv' into sequences

    """
    #!/usr/bin/env bash

    if command -v conda &> /dev/null
    then
        if conda env list | grep "^R " >/dev/null 2>/dev/null
        then
            # there is a conda environment named "R"
            conda activate R
            Rscript $PWD/bin/1_run_trinuc_matching.R ${filename} ${stoppingCriterion} ${maxTime} ${max_fraction_removed_trinucs} ${acceleration_score} ${euclidean_change_ratio} ${utils} ${good_mappability_regions}
        else
            # no conda environment named "R"
            Rscript $PWD/bin/1_run_trinuc_matching.R ${filename} ${stoppingCriterion} ${maxTime} ${max_fraction_removed_trinucs} ${acceleration_score} ${euclidean_change_ratio} ${utils} ${good_mappability_regions}
        fi
    else
        # no conda installed
        Rscript $PWD/bin/1_run_trinuc_matching.R ${filename} ${stoppingCriterion} ${maxTime} ${max_fraction_removed_trinucs} ${acceleration_score} ${euclidean_change_ratio} ${utils} ${good_mappability_regions}
    fi
    """
}

process apply_trinuc_matching_to_tracks {

    time = { 2.hour }
    memory = { (params.memory_process2 + 5*(task.attempt-1)).GB }

    input:
    file full_tracks_trinuc32_freq from full_tracks_trinuc32_freq
    file matched_tracks from matched_tracks
    val euclidean_score from euclidean_score
    file sequences from sequences

    output:
    file '*__3ntMatched_euclidean-*.bed.gz' into matched_tracks_granges

    """
    #!/usr/bin/env bash

    if command -v conda &> /dev/null
    then
        if conda env list | grep "^R " >/dev/null 2>/dev/null
        then
            # there is a conda environment named "R"
            conda activate R
            Rscript $PWD/bin/2_apply_trinuc_matching_to_tracks.R ${full_tracks_trinuc32_freq} ${matched_tracks} ${euclidean_score} ${sequences}
        else
            # no conda environment named "R"
            Rscript $PWD/bin/2_apply_trinuc_matching_to_tracks.R ${full_tracks_trinuc32_freq} ${matched_tracks} ${euclidean_score} ${sequences}
        fi
    else
        # no conda installed
        Rscript $PWD/bin/2_apply_trinuc_matching_to_tracks.R ${full_tracks_trinuc32_freq} ${matched_tracks} ${euclidean_score} ${sequences}
    fi
    """
}

matched_tracks_granges
    .println { "Finished! Trinucleotide-matched tracks saved in res/" }
