#!/bin/bash

conda activate nextflow

#export TOWER_ACCESS_TOKEN=""

mkdir -p log/

nextflow -log $PWD/log/nextflow.log run pipe.nf --file_paths=$PWD/input/file_paths.csv \
												--stoppingCriterion=0.2 \
												--maxTime=2 \
												--unitsTime="days" \
												--max_fraction_removed_trinucs=0.25 \
												--acceleration_score=100 \
												--euclidean_change_ratio=0.1,1.1 \
												--fast_progress_lfc_cutoff=-0.00001 \
												--progress_its=1000 \
												--time_process1=60 \
												--time_process2=8 \
												--memory_process1=100 \
												--memory_process2=350 \
												--utils=$PWD/bin/utils.R \
												--good_mappability_regions=none \
												-resume #\ -with-tower
