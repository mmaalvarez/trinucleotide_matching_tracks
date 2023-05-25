#!/bin/bash

conda activate nextflow

#export TOWER_ACCESS_TOKEN="eyJ0aWQiOiA0ODk1fS4yNDdkNzEyYmY5NGUxMmFlOTQ1OGYwYWJlYmI2MjY0YmU2Y2E4Yzdl"

mkdir -p log/

nextflow -log $PWD/log/nextflow.log run pipe.nf --filename "/g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/TCGA_PCAWG_Hartwig_CPTAC_POG_MMRFCOMMPASS/SHM/CTCF_cohesin/1_CTCF_cohesin_peaks_coords/CTCF_cohesin_peaks_coords" \
												--stoppingCriterion 0.01 \
												--maxTime 8 \
												--max_fraction_removed_trinucs 0.5 \
												--acceleration_score 1 \
												--euclidean_change_ratio 0.1,1.1 \
												--memory_process1 9 \
												--memory_process2 40 \
												--utils /g/strcombio/fsupek_data/users/malvarez/projects/RepDefSig/bin/utils.R \
												--good_mappability_regions /g/strcombio/fsupek_home/mmunteanu/reference/CRG75_nochr.bed \
												-resume #\ -with-tower
