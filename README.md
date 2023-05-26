In a track file (with .bed-like format), remove trinucleotide counts from genomic coordinates to approximate them to the mean proportions of the remaining coordinates

It first converts the track file into a data.frame with N rows (1 per genomic coordinate, e.g. chr1_1_237572, chr1_237573_237942) and 32 columns (the 64 possible trinucleotides are collapsed into 32 by complementarity)

Not downsampling coordinates with too few counts

You can process several track files in parallel

### Iterate until:

- Trinucleotide count frequencies match across coordinates (Euclidean distance)
- Too many counts removed
- Max. number of iterations reached
- Max. time reached

## How to use

You have to add to $PWD/input/file_paths.tsv the absolute path + file name of each of the bed-formatted track files you want to process, one per row (see $PWD/input/file_paths_example.txt)

To run the pipeline, execute run_nf.sh:
`bash run_nf.sh`

## Other

"good_mappability_regions" (called $PWD/crg75/CRG75_nochr.bed, see example of how the first lines should look) were obtained from https://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign75mer.bigWig

To build the container, run in $PWD/container: 
`sudo singularity build container.sif container.def`
