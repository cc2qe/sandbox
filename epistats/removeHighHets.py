#!/usr/bin/env python

import argparse
import sys

parser = argparse.ArgumentParser(description='remove entries with heterozygosity rate over specified percent')
parser.add_argument('-i', '--input', default=sys.stdin, type=argparse.FileType('r'), help='input file (default stdin)')
parser.add_argument('-m', '--max_het_rate', type=float, required=True, help='max heterozygosity rate of the population [0-1]')

args = parser.parse_args()

f = args.input
max_het_rate = args.max_het_rate

for l in f:
    l = l.rstrip()
    v = l.split('\t')
    numHets = v[7:].count('1')
    numSamples = len(v[7:])
    
    if (float(numHets) / numSamples) < max_het_rate:
        print l
    else:
        continue
