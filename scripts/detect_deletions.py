#!/usr/bin/python
from optparse import OptionParser
import gzip
import sys
from Bio import SeqIO
import os

"""
Author: Tom van der Valk
Contact: tvdvalk1989@gmail.com
Date: May 18,2022
Usage: python detect_deletions.py --depthfile depth.txt --samples samples.txt
"""

usage = "python detect_deletions.py --depthfile depth.txt --samples samples.txt"
parser = OptionParser(usage=usage)
parser.add_option("--depthfile", action="store", type="string", dest="depthfile",help="output from samtools depth")
parser.add_option("--samples", action="store", type="string", dest="samples",help="tab seperated file containing sample names (column 1) and sample groups (column 2)")
(options, args) = parser.parse_args()

#store group names into list
contig_name = (options.depthfile).split(".")[0]
group_list = []
with open(options.samples) as f1:
    for line in f1:
        group_name = line.strip().split("\t")[0]
        group_list += [group_name]

#obtain list indices for each sample group
group_set = set(group_list)
if len(group_set) > 3:
    print("too many groups list in the sample-file, maximum of 3 groups allowed")
elif len(group_set) < 3:
    print("too few groups listed in the sample-file, 3 groups are neccesary")
else:
    reference_list = []
    outgroup_list = []
    target_list = []
    for i in range(len(group_list)):
        if group_list[i].lower() == "reference":
            reference_list += [i]
        elif group_list[i].lower() == "outgroup":
            outgroup_list += [i]
        elif group_list[i].lower() == "target":
            target_list += [i]
        else:
            print("unrecognised grouping:" + i + " in sample-file, only group names: reference, outgroup, or targetgroup are allowed")

#Write the average coverage for each group to file
reference_coverage_dict = {}
for i in range(1000):
    reference_coverage_dict[i] = 0
outgroup_sum_coverage = [0,0]
target_sum_coverage = [0,0]
with open(options.depthfile) as f1:
    for line in f1:
        splitted = line.strip().split("\t")
        depth = [int(i) for i in splitted[2:]]
        reference_coverage = sum(list(map(depth.__getitem__, reference_list)))
        outgroup_coverage = sum(list(map(depth.__getitem__, outgroup_list)))
        target_coverage = sum(list(map(depth.__getitem__, target_list)))
        if reference_coverage > 0:
            reference_coverage_dict[reference_coverage] += 1
        if outgroup_coverage > 0:
            outgroup_sum_coverage[0] += outgroup_coverage
            outgroup_sum_coverage[1] += 1
        if target_coverage > 0:
            target_sum_coverage[0] += target_coverage
            target_sum_coverage[1] += 1

reference_average_coverage = max(reference_coverage_dict, key=reference_coverage_dict.get)
outgroup_average_coverage = round(outgroup_sum_coverage[0]/outgroup_sum_coverage[1],2)
target_average_coverage = round(target_sum_coverage[0]/target_sum_coverage[1],2)
print("average coverage (sum of all samples per group):" + "\n" + \
    "reference:" + " " + str(reference_average_coverage) + "X" + "\n" + \
        "outgroup:" + " " + str(outgroup_average_coverage) + "X" + "\n" + \
            "target:" + " " + str(target_average_coverage) + "X")
print(" ")
outputfile_coverage = open(contig_name + ".average-coverage","w")
outputfile_coverage.write("reference" + "\t" + str(reference_average_coverage) + "\n" + \
                            "outgroup" + "\t" + str(outgroup_average_coverage) + "\n" + \
                            "target" + "\t" + str(target_average_coverage) + "\n")
outputfile_coverage.close()

#Find deletions as regions with ultra low coverage and write to BED file
outputfile_reference = open(contig_name + ".reference.bed", "w")
outputfile_outgroup = open(contig_name + ".outgroup.bed", "w")
outputfile_target = open(contig_name + ".target.bed", "w")

    #set up matrix to obtain coverage in 100bp windows
n = 100
sample_size = len(group_list)
depth_matrix = [0] * n
for i in range(n):
    depth_matrix[i] = [0] * sample_size

outputfile_window_depth = open("windowed." + options.depthfile,"w")
    #read depth file and find regions with low coverage
site_counter = 0
with open(options.depthfile) as f1:
    for line in f1:
        splitted = line.strip().split("\t")
        chromosome,position,depth = splitted[0],int(splitted[1]),[int(i) for i in splitted[2:]]
        site_counter += 1
        if site_counter % 100000 == 0:
            print("processed " + str(site_counter) + " sites")
        depth_matrix += [depth]
        del depth_matrix[0]
        average_depth_list = [0] * sample_size
        for i in range(100):
            for j in range(sample_size):
                average_depth_list[j] += (depth_matrix[i][j]/100)
        reference_depth = sum(list(map(average_depth_list.__getitem__, reference_list)))
        outgroup_depth = sum(list(map(average_depth_list.__getitem__, outgroup_list)))
        target_depth = sum(list(map(average_depth_list.__getitem__, target_list)))
        outputfile_window_depth.write(chromosome + "\t" + str(position) + "\t" + str(round(reference_depth,2)) + "\t" + str(round(outgroup_depth,2)) + "\t" + str(round(target_depth,2)) + "\n")
        if reference_depth < (0.9 * reference_average_coverage) or reference_depth > (1.1 * reference_average_coverage):
            outputfile_reference.write(chromosome + "\t" + str(position-1) + "\t" + str(position) + "\n")
        if outgroup_depth < (0.10 * outgroup_average_coverage):
            outputfile_outgroup.write(chromosome + "\t" + str(position-1) + "\t" + str(position) + "\n")
        if target_depth < (0.10 * target_average_coverage):
            outputfile_target.write(chromosome + "\t" + str(position-1) + "\t" + str(position) + "\n")

outputfile_reference.close()
outputfile_outgroup.close()
outputfile_target.close()
