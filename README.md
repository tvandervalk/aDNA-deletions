# aDNA-deletions-scripts
Scripts used to detect genomic deletions in woolly mammoth genomes (van der Valk et. al 2022, in submision)

These scripts were specifically written to analyse ancient mammoth genomes >10X coverage and are not tested on other datasets. Use with caution on your own genomic data.

**Running requirements:**
  - Python3 (https://www.python.org/downloads/)
  - samtools (https://github.com/samtools/)
  - bwa (https://github.com/lh3/bwa)
  - bedtools2 (https://github.com/arq5x/bedtools2)

**How to run:**

The script first builds a reference mappability map. Next, regions where a set of outgroup genomes (e.g. asian elephants) have normal read coverage but a set of target genomes (e.g. woolly mammoths) do not show read support are identified. The output is a bedfile containing the deleted regions fixed in the target set, and plots of the average coverage around those regions. To run, a set of target and outgroup genomes has to be provided in a tab-separated file, with the 1st column depicting outgroup/target and the 2nd column the path to the bam-file (see filepath_example.txt for an example file):

```
aDNA-deletions/run_deletions.sh
    --contig        name of the contig to be analysed (e.g. Chr1)
    --reference     path to fasta reference used for mapping (e.g. /home/reference.fa)
    --samtools      path to samtools command (e.g. /home/samtools/samtools)
    --bwa           path to bwa command (e.g. /home/bwa/bwa)
    --bedtools      path to bedtools command (e.g. /home/bedtools2/bedtools)
    --scripts       path to folder containing the python scripts (e.g. /home/aDNA-deletions/scripts)
    --filepaths     path to txt file containing sample information
  ```
  
  The script analyses one chromosome/contig at a time. To speed up the analysis, one job per contig can be run simultaneously
