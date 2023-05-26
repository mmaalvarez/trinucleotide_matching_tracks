#!/usr/bin/env bash

if command -v conda &> /dev/null
then
    if conda env list | grep "^R " >/dev/null 2>/dev/null
    then
        # there is a conda environment named "R"
        conda activate R
        Rscript /g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/TCGA_PCAWG_Hartwig_CPTAC_POG_MMRFCOMMPASS/SHM/CTCF_cohesin/1_CTCF_cohesin_peaks_coords/trinucleotide_matching_tracks/bin/1_run_trinuc_matching.R CTCF_cohesin_peaks_coords 0.01 8 0.5 1 0.1,1.1 utils.R CRG75_nochr.bed
    else
        # no conda environment named "R"
        Rscript /g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/TCGA_PCAWG_Hartwig_CPTAC_POG_MMRFCOMMPASS/SHM/CTCF_cohesin/1_CTCF_cohesin_peaks_coords/trinucleotide_matching_tracks/bin/1_run_trinuc_matching.R CTCF_cohesin_peaks_coords 0.01 8 0.5 1 0.1,1.1 utils.R CRG75_nochr.bed
    fi
else
    # no conda installed
    Rscript /g/strcombio/fsupek_cancer3/malvarez/WGS_tumors/somatic_variation/TCGA_PCAWG_Hartwig_CPTAC_POG_MMRFCOMMPASS/SHM/CTCF_cohesin/1_CTCF_cohesin_peaks_coords/trinucleotide_matching_tracks/bin/1_run_trinuc_matching.R CTCF_cohesin_peaks_coords 0.01 8 0.5 1 0.1,1.1 utils.R CRG75_nochr.bed
fi
