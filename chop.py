#!/usr/bin/python3

'''Useage: python chop.py -f input.fasta -r <number of reads to generate> -l <mean fragment length> -sd <standard deviation of fragment sizes>'''

import sys
import numpy as np
import argparse
import random
from Bio import SeqIO

#set up script parameters
parser = argparse.ArgumentParser()
parser.add_argument('-f', '--fasta', help='input fasta file')
parser.add_argument('-r', '--reads', help='number of reads to simulate')
parser.add_argument('-sd', '--stddev', help='standard deviation of fragment length distribution')
parser.add_argument('-l', '--meanlength', help='mean fragment length')

args = parser.parse_args()

####TO DO: Need to institute paired end, insert size (as in the paired end data??)

#dictionary of fasta headers that are in the file (in case of multi fasta input)
names = [rec.name for rec in SeqIO.parse(args.fasta, "fasta")]
lengths = [len(rec.seq) for rec in SeqIO.parse(args.fasta, "fasta")]
fastaDict = dict(zip(names, lengths))

#use loop to randomly pick one of the reads, get the name and length to get starting coordinate
for i in range(1, int(args.reads)+1):
        read, length = random.choice(list(fastaDict.items()))
        start = np.random.randint(1, int(length)+1)
        fragLen = int(np.random.normal(int(args.meanlength), int(args.stddev), 1))

        for record in SeqIO.parse(open(args.fasta, "r"), "fasta"):
                if record.name == read:
                        if int(start)+int(fragLen) >= length:
                                newSeq = record.seq[int(start):int(start)-int(fragLen)]
                                revComp = newSeq.reverse_complement()
                                print(">%s_%i|%i\n%s" % (read, i, len(newSeq), newSeq))
                        elif int(start)+int(fragLen) <= length:
                                newSeq = record.seq[int(start):int(start)+int(fragLen)]
                                print(">%s_%i|%i\n%s" % (i, read, len(newSeq), newSeq))


