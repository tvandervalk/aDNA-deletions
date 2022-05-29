# aDNA-deletions-scripts
Scripts used to detect genomic deletions in woolly mammoth genomes (van der Valk et. al 2022, iScience))

Note that these scripts are specifically written to analyse (ancient) mammoth genomes > 10X coverage and are not tested on other datasets. Use with caution on your own genomic data.

Running requirements:
  - Python3
  - samtools
  - bwa
  - BEDTools

How to run:
The scripts analyse one chromosome/contig at a time. To speed up the analysis one job per contig can be submitted simultaneously
