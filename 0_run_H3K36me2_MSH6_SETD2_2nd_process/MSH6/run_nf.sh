#!/bin/bash

conda activate nextflow

mkdir -p log/


filename="SRR13314132_sorted_FE"
input_relative_path="../../work/5c/1cc7e0fd218e12d3bb617caadc3555/"$filename"_"

nranges=`wc -l "$input_relative_path"full_tracks_trinuc32_freq.tsv | cut -d" " -f1`
# split input tables into n parts
nparts=100
fraction_length=`echo $(($nranges / $nparts))`

euclidean_score="0.4732095732911133"



nextflow -log $PWD/log/nextflow.log run pipe.nf --full_tracks_trinuc32_freq=$PWD/"$input_relative_path"full_tracks_trinuc32_freq.tsv \
												--matched_tracks=$PWD/"$input_relative_path"matched_tracks.tsv \
												--sequences=$PWD/"$input_relative_path"sequences.tsv \
												--nranges=$nranges \
												--fraction_length=$fraction_length \
												--filename=$filename \
												--euclidean_score=$euclidean_score \
												--time_process2=8 \
												--memory_process2=30 \
												--utils=$PWD/../../bin/utils.R \
												--good_mappability_regions=none \
												-resume
