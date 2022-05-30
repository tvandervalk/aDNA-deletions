#!/usr/bin/python
from optparse import OptionParser
import gzip
from Bio import SeqIO

"""
Author: Tom van der Valk
Contact: tvdvalk1989@gmail.com
Date: May 18,2022
Usage: python fasta_to_reads.py --fasta reference.fasta --readlength NN --target contiq_name --out outputfile_name
"""

#set options
usage = "usage: %prog --fasta reference.fasta --readlength NN --target contiq_name --out outputfile_name"
parser = OptionParser(usage=usage)
parser.add_option("--fasta", action="store", type="string", dest="inputfile",help="reference genome in fasta format")
parser.add_option("--readlength", action="store", type="string", dest="readlength",help="the read length used to build a reference mappability map")
parser.add_option("--target_contiq", action="store", type="string", dest="target",help="Chr/contiq name to create the mappability map for")
parser.add_option("--out", action="store", type="string", dest="output",help="Write output to file")
(options, args) = parser.parse_args()

#load fasta into memory
print("loading fasta reference into RAM....")
if options.inputfile.endswith(".gz"):
    fasta_dict = SeqIO.to_dict(SeqIO.parse(gzip.open(options.inputfile,"rt"), "fasta"))
else:
    fasta_dict = SeqIO.to_dict(SeqIO.parse(options.inputfile, "fasta"))

print("parsing fasta to reads...")

#write fasta reads to output file
readlength = int(options.readlength)
if options.target in fasta_dict:
    outputfile = gzip.open(options.output + ".fastq.gz","w")
    contig_value = fasta_dict[options.target]
    contig_sequence = str(contig_value.seq)
    readcounter = 0
    chromname = options.target
    for i in range(0,len(contig_sequence)-readlength,1):
        readcounter += 1
        read_sequence = contig_sequence[i:i+readlength].upper()
        if "N" not in read_sequence:
            outputfile.write(("@" + chromname + "_" + str(readcounter) + "\n").encode())
            outputfile.write((str(read_sequence) + "\n").encode())
    outputfile.close()
else:
    print("contig name:" + options.target + " not found in fasta reference")
