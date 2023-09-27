#!/bin/bash

conda activate nextflow

#export TOWER_ACCESS_TOKEN=""

mkdir -p log/


## THE nranges and fraction_length HAVE TO BE OBTAINED WITHIN pipe.nf from each input files in file_paths

# filename="H3K36me2_pooled_GSE118954-149670-175750"
# input_relative_path="../../work/54/2c44d8fe769335dc4d4a73ad46a0e3/"$filename"_"
# nranges=`wc -l "$input_relative_path"full_tracks_trinuc32_freq.tsv | cut -d" " -f1`
## split input tables into n parts
# nparts=100
# fraction_length=`echo $(($nranges / $nparts))`


nextflow -log $PWD/log/nextflow.log run pipe.nf --file_paths=$PWD/input/file_paths.csv \
												--stoppingCriterion=0.01 \
												--maxTime=2 \
												--unitsTime="days" \
												--max_fraction_removed_trinucs=0.25 \
												--acceleration_score=100 \
												--euclidean_change_ratio=0.1,1.1 \
												--fast_progress_lfc_cutoff=-0.00001 \
												--progress_its=1000 \
												--time_process1=60 \
												--time_process2=1 \
												--memory_process1=100 \
												--memory_process2=20 \
												--utils=$PWD/bin/utils.R \
												--good_mappability_regions=none \
												-resume #\ -with-tower
												#--nranges=$nranges \
												#--fraction_length=$fraction_length \