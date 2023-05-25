## In a track file (with .bed-like format), remove trinucleotide counts from genomic coordinates to approximate them to the mean proportions of the remaining coordinates

## It first converts the track file into a data.frame with N rows (1 per genomic coordinate, e.g. chr1_1_237572, chr1_237573_237942) and 32 columns (the 64 possible trinucleotides are collapsed into 32 by complementarity)

## Not downsampling coordinates with too few counts

## Iterate until:

- Trinucleotide count frequencies match across coordinates (Euclidean distance)
- Too many counts removed
- Max. number of iterations reached
- Max. time reached


# How to use

## In run_nf.sh, at --filename "XXX" you have to type the path+name of the bed-formatted tracks file WITHOUT its extension (.bed, .bw, .tsv...), and prepend the relative path if it's in other folder, e.g. if it is in the parent folder:
`--filename "/g/strcombio/rest/of/path/my_bed_file_name" \`

## To run the pipeline, execute run_nf.sh:
`bash run_nf.sh`


# Other

## It uses a container at "/g/strcombio/fsupek_data/users/malvarez/projects/RepDefSig/containers/regressions/container.sif"
- Its .def file is container_def/container.def
- It actually includes libraries not used in this pipeline

## The 'trinuc_matching' function is sourced from /g/strcombio/fsupek_data/users/malvarez/projects/RepDefSig/bin/utils.R

## good_mappability_regions are from /g/strcombio/fsupek_home/mmunteanu/reference/CRG75_nochr.bed
