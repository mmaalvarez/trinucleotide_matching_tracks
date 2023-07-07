#!/bin/bash

conda activate nextflow

mkdir -p log/


filename="SRR13314134_sorted_FE"
input_relative_path="../../work/ee/823f9be3ff4101082ace1624d48875/"$filename"_"

nranges=`wc -l "$input_relative_path"full_tracks_trinuc32_freq.tsv | cut -d" " -f1`
# split input tables into n parts
nparts=100
fraction_length=`echo $(($nranges / $nparts))`

euclidean_score="0.5767566448037298"



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
