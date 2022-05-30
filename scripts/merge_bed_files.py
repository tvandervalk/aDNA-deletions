#!/usr/bin/python
from optparse import OptionParser
import gzip
from Bio import SeqIO

"""
Author: Tom van der Valk
Contact: tvdvalk1989@gmail.com
Date: May 18,2022
Usage: merge_bed_files.py --target contiq_name
"""

usage = "python merge_bed_files.py --target Chromosome/contig name"
parser = OptionParser(usage=usage)
parser.add_option("--target", action="store", type="string", dest="target",help="Chromosome/contig name")
(options, args) = parser.parse_args()

outputfile_reference_outgroup = open(options.target + ".intermediate.reference_outgroup.bed","w")
outputfile_target = open(options.target + ".intermediate.target.bed","w")

#Find last position in BED files
bed_dict = {}
final_end = 0
with open(options.target + ".merged.reference.bed", 'r') as f1:
    for line in f1:
        splitted = line.strip().split("\t")
        start,end = int(splitted[1]),int(splitted[2])
        if end > final_end:
            final_end = end + 1000

with open(options.target + ".merged.outgroup.bed", 'r') as f1:
    for line in f1:
        splitted = line.strip().split("\t")
        start,end = int(splitted[1]),int(splitted[2])
        if end > final_end:
            final_end = end + 1000

with open(options.target + ".merged.target.bed", 'r') as f1:
    for line in f1:
        splitted = line.strip().split("\t")
        start,end = int(splitted[1]),int(splitted[2])
        if end > final_end:
            final_end = end + 1000

#Merge BED files into dictionary
for i in range(final_end):
    bed_dict[i] = [0,0,0]

with open(options.target + ".merged.reference.bed", 'r') as f1:
    for line in f1:
        splitted = line.strip().split("\t")
        start,end = int(splitted[1]),int(splitted[2])
        for i in range(start,end):
            bed_dict[i][0] = 1

with open(options.target + ".merged.outgroup.bed", 'r') as f1:
    for line in f1:
        splitted = line.strip().split("\t")
        start,end = int(splitted[1]),int(splitted[2])
        for i in range(start,end):
            bed_dict[i][1] = 1

with open(options.target + ".merged.target.bed", 'r') as f1:
    for line in f1:
        splitted = line.strip().split("\t")
        chromname,start,end = splitted[0],int(splitted[1]),int(splitted[2])
        for i in range(start,end):
            bed_dict[i][2] = 1

#Write overlapping BED regions into BED file
for j in range(1,final_end):
    coverage = bed_dict[j]
    reference = coverage[0]
    outgroup = coverage[1]
    target = coverage[2]
    if reference == 1 or outgroup == 1:
        outputfile_reference_outgroup.write(chromname + "\t" + str(j) + "\t" + str(j+1) + "\n")
    if reference == 0 and outgroup == 0 and target == 1:
        outputfile_target.write(chromname + "\t" + str(j) + "\t" + str(j+1) + "\n")

outputfile_reference_outgroup.close()
outputfile_target.close()
