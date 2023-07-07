#!/usr/bin/env nextflow


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
    path full_tracks_trinuc32_freq from params.full_tracks_trinuc32_freq
    path matched_tracks from params.matched_tracks
    val euclidean_score from params.euclidean_score
    path sequences from params.sequences
    val filename from params.filename
    path utils from params.utils

    output:
    file '*__3ntMatched_euclidean-*.part.bed' into res

    """
    #!/usr/bin/env bash

    conda activate R
    Rscript $PWD/apply_trinuc_matching_to_tracks.R ${start_int} ${end_int} ${full_tracks_trinuc32_freq} ${matched_tracks} ${euclidean_score} ${sequences} ${filename} ${utils}
    """
}

filename = params.filename
euclidean_score = params.euclidean_score

res
    .collectFile(name: 'res/${filename}__3ntMatched_euclidean-${euclidean_score}.bed')
    .println { "Finished, results saved in res/" }
