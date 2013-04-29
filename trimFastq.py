#!/usr/bin/env python

"""
Trim length of fastq reads
"""

import sys, argparse

__author__ = "Colby Chiang (chiang@chgr.mgh.harvard.edu)"
__version__ = "$Revision: 0.0.1 $"
__date__ = "$Date: 2010/11/10 13:25 $"

parser = argparse.ArgumentParser(description='Trim length of fastq reads')
parser.add_argument('-i', '--input', type=file, help='Input fastq file')
parser.add_argument('-s', '--softclip', action='store_true', help='Soft clip trailing bases instead of removing them.')
parser.add_argument('-l', '--length', type=int, help='Bases to trim to')

args = parser.parse_args()

fqIn = args.input
length = args.length
softClip = args.softclip

if fqIn == None:
	fqIn = sys.stdin

lineNumber = 0

for line in fqIn:
	line = line.rstrip()
	lineNumber += 1
	
	if not softClip:
		if lineNumber % 2 == 0: 
			print line[0:length]
		else: print line
	elif softClip:
		if lineNumber % 4 == 0: 
			print line[0:length] + "#"*(len(line)-length)
		else: print line
