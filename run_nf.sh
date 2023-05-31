#!/bin/bash

conda activate nextflow

#export TOWER_ACCESS_TOKEN=""

mkdir -p log/

nextflow -log $PWD/log/nextflow.log run pipe.nf --file_paths=$PWD/input/file_paths.csv \
												--stoppingCriterion=0.01 \
												--maxTime=8 \
												--max_fraction_removed_trinucs=0.5 \
												--acceleration_score=1 \
												--euclidean_change_ratio=0.1,1.1 \
												--fast_progress_lfc_cutoff=-0.00001 \
												--progress_its=1000 \
												--memory_process1=15 \
												--memory_process2=40 \
												--utils=$PWD/bin/utils.R \
												--good_mappability_regions=none \
												-resume #\ -with-tower
