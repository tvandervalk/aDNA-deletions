#!/usr/bin/python
from __future__ import division
import matplotlib.pyplot as plt
import gzip
import time
import numpy as np
from string import ascii_lowercase
import itertools
from optparse import OptionParser
plt.rcParams['figure.figsize'] = [25, 5]
from sys import argv

"""
Author: Tom van der Valk
Contact: tvdvalk1989@gmail.com
Date: May 18,2022
Usage: python plot_deletions.py --depth_file windowed_depth_file --average_coverage coverage_file --outgroup_bed outgroup_bed --target_bed target_bed_file
"""

#set options
usage = "usage: %prog python plot_deletions.py --depth_file windowed_depth_file --average_coverage coverage_file --outgroup_bed outgroup_bed --target_bed target_bed_file"
parser = OptionParser(usage=usage)
parser.add_option("--depth_file", action="store", type="string", dest="depth_file",help="windowed depth outputfile from detect_deletions.py script")
parser.add_option("--average_coverage", action="store", type="string", dest="average_coverage",help="average coverage outputfile from detect_deletions.py script")
parser.add_option("--outgroup_bed", action="store", type="string", dest="outgroup_bed",help="Reference + outgroup output bed file from detect_deletions.py script")
parser.add_option("--target_bed", action="store", type="string", dest="target_bed",help="Target output bed file from detect_deletions.py script")
(options, args) = parser.parse_args()
chromname = options.depth_file.split(".")[1]

#safe average coverage per group in list
coverages = []
with open(options.average_coverage) as f1:
    for line in f1:
        splitted = line.strip().split("\t")
        coverages += [float(splitted[1])]

#create list of characters to name output plots
letters = list(ascii_lowercase)
letter_combi_list = list(itertools.combinations_with_replacement(letters, 3))
letter_indexer = 0

#save depth per group into dictionary
coverage_dict = {}
start = 10 ** 100
with open(options.depth_file) as f1:
    for line in f1:
        splitted = line.strip().split("\t")
        parsed_splitted = []
        pos = int(splitted[1])
        if pos < start:
            start = pos
        chromlength = pos
        splitted = [float(i) for i in splitted[2:]]
        cov_reference,cov_outgroup,cov_target = splitted[0],splitted[1],splitted[2]
        if cov_reference < 1:
            parsed_splitted += [0]
        else:
            parsed_splitted += [cov_reference/coverages[0]]
        if cov_outgroup < 0.10 * coverages[1]:
            parsed_splitted += [0]
        else:
            parsed_splitted += [cov_outgroup/coverages[1]]
        if cov_target < 0.10 * coverages[2]:
            parsed_splitted += [0]
        else:
            parsed_splitted += [cov_target/coverages[2]]
        if (cov_reference+cov_outgroup+cov_target) > 0:
            coverage_dict[pos] = parsed_splitted


print("loaded coverage")

#check if there is a deletion in the focal window
windowsize = 100000
for i in range(start,chromlength,windowsize):
    start_window = i
    end_window = i + windowsize
    counter = 0
    plotter = False
    x_axis_scatter = []
    y_axis_scatter = []
    with open(options.target_bed) as f1:
        for line in f1:
            splitted = line.strip().split("\t")
            start,end = int(splitted[1]),int(splitted[2])
            if start >= start_window and end <= end_window:
                plotter = True
            if start < start_window and end < end_window and end > start_window:
                plotter = True
            if start < start_window and end > end_window:
                plotter = True
            if start > start_window and start < end_window and end > end_window:
                plotter = True

    #If deletion is detected, make a plot
    if plotter == True:
        x_axis = []
        y_axis_reference = []
        y_axis_outgroup = []
        y_axis_target = []
        for focal_position in range(start_window,end_window):
            x_axis += [focal_position]
            if focal_position in coverage_dict:
                parsed_splitted = coverage_dict[focal_position]
            else:
                parsed_splitted = [0,0,0]
            y_axis_reference += [parsed_splitted[0]]
            y_axis_outgroup += [parsed_splitted[1]]
            y_axis_target += [parsed_splitted[2]]

        bed_list_outgroup = []
        with open(options.outgroup_bed) as f1:
            for line in f1:
                splitted = line.strip().split("\t")
                start,end = int(splitted[1]),int(splitted[2])
                if start >= start_window and end <= end_window:
                    bed_list_outgroup += [(int(splitted[1]),int(splitted[2]))]
                elif start < start_window and end <= end_window and end > start_window:
                    bed_list_outgroup += [(start_window,int(splitted[2]))]
                elif start >= start_window and start < end_window and end > end_window:
                    bed_list_outgroup += [((int(splitted[1]),end_window))]
                elif start < start_window and end > end_window:
                    bed_list_outgroup += [(start_window,end_window)]


        bed_list_target = []
        with open(options.target_bed) as f1:
            for line in f1:
                splitted = line.strip().split("\t")
                start,end = int(splitted[1]),int(splitted[2])
                start,end,bed_number = int(splitted[1]),int(splitted[2]),str(splitted[3])
                if start >= start_window and end <= end_window:
                    bed_list_target += [(int(splitted[1]),int(splitted[2]),bed_number)]
                elif start > start_window and start < end_window and end > end_window:
                    bed_list_target += [(int(splitted[1]),end_window,bed_number)]
                elif start < start_window and end < end_window and end > start_window:
                    bed_list_target += [(start_window,int(splitted[2]),bed_number)]
                elif start < start_window and end > end_window:
                    bed_list_target += [(start_window,end_window,bed_number)]

            fig, (ax1, ax2, ax3) = plt.subplots(3, 1,sharex=True)
            ax1.plot(x_axis,y_axis_reference,color="black")
            ax2.plot(x_axis,y_axis_outgroup,color="black")
            ax3.plot(x_axis,y_axis_target,color="black")
            ax3.xaxis.set_ticks(np.arange(min(x_axis), max(x_axis)+1, 10000))
            ax3.ticklabel_format(axis="x", style="sci", scilimits=(0,0),useOffset=False)
            ax1.set_ylim(0,3)
            ax2.set_ylim(0,3)
            ax3.set_ylim(0,3)
            ax3.set_xlabel("genomic position",fontsize=14)
            ax2.set_ylabel("relative coverage",fontsize=14)
            ax1.set_title("reference mappability", y=1.0, pad=-14)
            ax2.set_title("Outgroup coverage", y=1.0, pad=-14)
            ax3.set_title("Target coverage", y=0.975, pad=-14)
            for j in bed_list_outgroup:
                if j:
                    start,end = j[0],j[1]
                    ax1.axvspan(start,end, facecolor='red', alpha=0.5)
            for j in bed_list_target:
                if j:
                    start,end = j[0],j[1]
                    start,end,plot_number = j[0],j[1],j[2]
                    ax3.axvspan(start,end, facecolor='green', alpha=0.5)
                    ax3.text(start,2,str(plot_number),fontsize=10)
            file_letter = "".join(letter_combi_list[letter_indexer])
            plt.savefig(file_letter + "_" + chromname + ".window-" + str(start_window) + ".pdf",format="pdf")
            plt.close()
            letter_indexer += 1
