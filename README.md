# RCMC Analysis Code
This repository contains source code for the article [Region Capture Micro-C reveals coalescence of enhancers and promoters into nested microcompartments](https://www.nature.com/articles/s41588-023-01391-1), used in the analysis of RCMC data. 

Code is provided either in the form of Python/R scripts or as Jupyter notebooks to be run in conda environments containing the required packages. Additionally, genomic positions of microcompartments identified in the paper are included in bedpe format.

## Code summary
### Micro-C alignment (microcbowtie2.py)
Required packages:
-	bowtie2
-	samtools
-	sambamba
-	pairtools
-	cooler
-	pairix

Python script used to align reads in .fastq format from paired-end sequencing of Micro-C experiments and produces as output .pairs, .cool and .mcool files compatible with downstream applications such as HiGlass.

Example usage:

```
python /path/to/script/microcbowtie2.py --file_1 pair1.fastq --file_2 pair2.fasq -g mm39 -t 36 -o exampleoutput
```

### ChIP-seq alignment (spikeinChIP_PE_alignment.py)
Required packages:
-	bowtie2
-	samtools
-	sambamba

Python script used to align reads in .fastq format from paired-end sequencing of ChIP-seq experiments and produces aligned .bam files. .fastq information is input as a .tsv with each line containing the path to the first pair .fastq, the path to the second pair .fastq, and the desired output name of the aligned file.

Example usage:

```
python /path/to/script/spikeinChIP_PE_alignment.py -f list_of_fastqs.tsv -g mm39 -t 36 -o alignmentcountsout
```

### Finding chromatin features overlapping microcompartment anchors (loopFeatureOverlap.R)
Required packages:
-	plyr
-	dplyr
-	reshape2
-	purrr
-	grid
-	IRanges
-	GenomicRanges
-	arrangements
-	foreach

R script used to classify microcompartment interactions by finding overlap between identified microcompartments (.bedpe) and chromatin features (.bed) such as promoters, enhancers, CTCF binding sites, etc. It outputs individual .bedpe files of interactions according to combinatorial classification of chromatin features (e.g for enhancer (E) and promoter (P): P-P, E-E, E-P, E-null, P-null, null-null), including interactions which have no overlap (null category). Classification can be mutually exclusive (E-P cannot also be P-P) or inclusive.

Example usage:

```
Rscript /path/to/script/loopFeatureOverlap.R -l interactions.bedpe -b promoter.bed,enhancer.bed -i P,E -o outputdirectory/
```

### Calculating strength of individual interactions (LoopStrengthRCMC.ipynb) or aggregate pileup analysis (PileupsRCMC.ipynb)
Required packages:
-	seaborn
-	coolpuppy
-	cooltools
-	cooler

Jupyter notebooks used to calculate strengths of individual microcompartments (LoopStrengthRCMC.ipynb) or generate aggregate pileup analysis figures (PileupsRCMC.ipynb). Each one takes .mcool files of contacts from RCMC, a list of interactions to calculate for (.bedpe format), expected files generated by cooltools for each .mcool, and the captured region, and calculates background corrected observed/expected interaction strengths, either for each interaction individually (output as a .bedpe with additional columns for strengths for each .mcool) or as a pileup (output as a .pdf of the aggregate interaction annotated with the calculated strength).

### Visualization of contact maps and genomic tracks (ContactMapVisualizationExampleNotebook.ipynb)
Required packages:
-	cooltools
-	cooler
-	coolbox
-	matplotlib

Jupyter notebook used to generate visualizations of contact maps and genomics tracks for figures. Contact map visualization is accomplished using cooltools and requires a .mcool file of contacts from RCMC or a comparable method. Genomic track visualization is accomplished using coolbox and requires a .mcool file of contacts, gene annotations (.gtf format or similar), and ChIP-seq, RNA-seq, and ATAC-seq datasets (.bw format).

### Calculation of read-containing bin fraction by contact distance (CalculatingFilledBinFractionByDistance.ipynb)
Required packages:
-	cooltools
-	cooler
-	matplotlib

Jupyter notebook used to calculate the fraction of bins in .mcool-derived contact maps which contain at least one read pair at a given resolution and contact distance from the diagonal. The notebook takes an unbalanced .mcool of contacts from RCMC or a comparable method, tabulates the occupied contact bin fraction at specified contact distances, and generates a plot of occupied bin fraction by contact distance.

### Calculation of row sums in ICE-balanced .mcools (BalancedRowsumsCalculation.ipynb)
Required packages:
-	cooltools
-	cooler
-	matplotlib

Jupyter notebook used to confirm successful ICE balancing of .mcool files by calculating and plotting row sums. Two variations of the calculation are provided in the notebook – one for calculating row sums across an entire region, and one for calculating row sums only for bins containing contact anchor sites. Both variations take ICE-balanced .mcool files of contacts from RCMC or a comparable method, and the latter additionally takes a list of contact anchor sites (.bed format). The distributions of calculated row sums in either variation are plotted as histograms.

### List of manually-annotated microcompartment loops (MicrocompartmentLoops_PlusMin1kb.bedpe)
BEDPE format file listing all 1091 manually-annotated microcompartment loops across the Ppm1g (chr5) and Klf1 (chr8) regions used in the microcompartment analysis scripts above. Coordinates are provided for the mm39 reference genome, and loop anchors are listed as plus-and-minus 1kb from each anchor’s point coordinate. Columns in the file are as follows: the first is the chromosome of the left loop anchor, the second is the coordinate of the left loop anchor minus 1 kb, the third is the coordinate of the left loop anchor plus 1 kb, and the remaining three columns are the same for the right loop anchor.

### Lists of probes used for capturing regions of interest (captureprobes_mm10.bed, captureprobes_mm39.bed)
BED format file listing the genomic locations of all probes used for capturing the Sox2 (chr3), Ppm1g (chr5), Nanog (chr6), Klf1 (chr8), and Fbn2 (chr18) regions used in capture. Coordinates are provided for both the mm10 and mm39 reference genomes, and loop anchors are listed as plus-and-minus 1kb from each anchor’s point coordinate. Columns in the file are as follows: the first is the chromosome the region is located on, the second is the start coordinate of the probe, and the third is the end coordinate of the probe.

## How to cite
This work is shared under an MIT license. If you make use of analysis scripts or data from this work, please cite as follows:

Goel, V.Y., Huseyin, M.K. & Hansen, A.S. Region Capture Micro-C reveals coalescence of enhancers and promoters into nested microcompartments. *Nat Genet* (2023). https://doi.org/10.1038/s41588-023-01391-1

Also refer to our deposited and citable code on Zenodo:

Goel, Viraat Y, Huseyin, Miles K, & Hansen, Anders S. (2023). Code supporting Region Capture Micro-C reveals coalescence of enhancers and promoters into nested microcompartments (1.0). Zenodo. https://doi.org/10.5281/zenodo.7641852
