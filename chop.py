#!/usr/bin/python3

'''Useage: python chop.py -f input.fasta -r <number of reads to generate> -i <mean insert length> -t <sequencing type (paired or single)> [optional: -sd <standard deviation of insert size normal distribution> -l <maximum sequencing length> -a <adapter sequence>]'''

import sys
import numpy as np
import argparse
import random
import os
from Bio import SeqIO


###TO DO: SET CIRCULAR GENOME OPTION (END OF FASTA NOT END OF CONTIG), GENERATE QUALITY SCORE SYSTEM -- OUTPUT FASTQ FILES


#set up script parameters
parser = argparse.ArgumentParser()
parser.add_argument('-f', '--fasta', help='input fasta file', required=True)
parser.add_argument('-r', '--reads', help='number of reads to simulate', required=True)
parser.add_argument('-sd', '--stddev', help='standard deviation of fragment length distribution', default='50')
parser.add_argument('-i', '--insertSize', help='mean fragment length (insert size)', required=True)
parser.add_argument('-l', '--maxlen', help='maximum length of read based on sequencing chemistry', default='200')
parser.add_argument('-t', '--seqtype', help='paired or single', required=True)

args = parser.parse_args()

#check if file exists  
assert os.path.exists(args.fasta), bcolors.WARNING + 'Yeah, either %s does not exist or you need to include the path' % args.fasta

#color text
class bcolors:
	WARNING = '\033[91m'
	COMPLETE = '\033[92m'

#dictionary of fasta headers that are in the file (in case of multi fasta input)
names = [rec.name for rec in SeqIO.parse(args.fasta, "fasta")]
lengths = [len(rec.seq) for rec in SeqIO.parse(args.fasta, "fasta")]
fastaDict = dict(zip(names, lengths))

def singleEnd(fastaDict):
	for i in range(int(args.reads)):
		read, length = random.choice(list(fastaDict.items()))
		start = np.random.randint(1, int(length))
		fragLen = int(np.random.normal(int(args.insertSize), int(args.stddev), 1))

		for record in SeqIO.parse(open(args.fasta, "r"), "fasta"):
			if record.name == read:
				if int(start)+int(fragLen) >= length: #make sure end position is not past the total fasta entry length
					newSeq = record.seq[int(start)-int(fragLen):int(start)]
					print(">%s_%i|%s|%i|%i:%i\n%s" % (read, i, args.fasta, len(newSeq), int(start)-int(fragLen), int(start), newSeq[0:int(args.maxlen)]))
				elif int(start)+int(fragLen) <= length:
					newSeq = record.seq[int(start):int(start)+int(fragLen)]
					print(">%s_%i|%s|%i|%i:%i\n%s" % (read, i, args.fasta, len(newSeq), int(start), int(start)+int(fragLen), newSeq[0:int(args.maxlen)]))


def pairedEnd(fastaDict):
        #use loop to randomly pick one of the reads, get the name and length to get starting coordinate
	for i in range(int(args.reads)):
		read, length = random.choice(list(fastaDict.items()))
		start = np.random.randint(1, int(length))
		fragLen = int(np.random.normal(int(args.insertSize), int(args.stddev), 1))

		for record in SeqIO.parse(open(args.fasta, "r"), "fasta"):
			if record.name == read:
				if int(start)+int(fragLen) >= length:
						readone = record.seq[int(start)-int(fragLen):int(start)]
						readtwo = record.seq[int(start):int(start)+int(fragLen)]
						print(">%s_%i|%s|%i|%i:%i\t1\n%s" % (read, i, args.fasta, len(readone), int(start), int(start)+int(fragLen), readone[0:int(args.maxlen)]))
						print(">%s_%i|%s|%i|%i:%i\t2\n%s" % (read, i, args.fasta, len(readtwo), int(start), int(start)+int(fragLen), readtwo[0:int(args.maxlen)]))
				elif int(start)+int(fragLen) <= length:
						readone = record.seq[int(start):int(start)+int(fragLen)]
						readtwo = record.seq[int(start)-int(fragLen):int(start)]
						print(">%s_%i|%s|%i|%i:%i\t1\n%s" % (read, i, args.fasta, len(readone), int(start), int(start)+int(fragLen), readone[0:int(args.maxlen)]))
						print(">%s_%i|%s|%i|%i:%i\t2\n%s" % (read, i, args.fasta, len(readtwo), int(start), int(start)+int(fragLen), readtwo[0:int(args.maxlen)]))

#initialize program
if args.seqtype == 'paired':
        pairedEnd(fastaDict)
elif args.seqtype == 'single':
        singleEnd(fastaDict)
else:
	print("SeqType must be either paired or single")

