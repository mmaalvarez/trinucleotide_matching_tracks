In a track file (with .bed-like format), remove trinucleotide counts from genomic coordinates to approximate them to the mean proportions of the remaining coordinates

It first converts the track file into a data.frame with N rows (1 per genomic coordinate, e.g. chr1_1_237572, chr1_237573_237942) and 32 columns (the 64 possible trinucleotides are collapsed into 32 by complementarity)

Not downsampling coordinates with too few counts

You can process several track files in parallel

### Iterate until:

- Trinucleotide count frequencies match across coordinates (Euclidean distance)
- Progress is too slow
- Too many counts removed
- Max. number of iterations reached
- Max. time reached

## How to use

You have to add to $PWD/input/file_paths.csv the absolute path + file name of each of the bed-formatted track files you want to process, one per row (see $PWD/input/file_paths_example.txt)

If using the container, you also need to edit the '/g/' in the nextflow.config file 'runOptions' to match the name of the parent directory where the file paths are: `runOptions = '-B /<HERE>/'`


To run the pipeline, execute run_nf.sh:
`bash run_nf.sh`

## Other

If you don't want to use a container, comment out the following block in the nextflow.config file: `container = "$PWD/container/container.sif"`

To build the container, run in $PWD/container: 
`sudo singularity build container.sif container.def`

Optional: filter out genomic regions not included in the "good_mappability_regions" CRG75 file (called $PWD/crg75/CRG75_nochr.bed, see example of how the first lines should look like) were obtained from https://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign75mer.bigWig
