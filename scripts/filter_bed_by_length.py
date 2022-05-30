#!/usr/bin/python
import sys
import gzip
import sys

#filter deletions below 500 basepairs and add deletion ID to BED file
counter = 0
fasta_read = []

for line in sys.stdin:
    splitted = line.strip().split("\t")
    start,end = int(splitted[1]),int(splitted[2])
    if (end-start) > 500:
        counter += 1
        print(line.strip() + "\t" + str(counter))
