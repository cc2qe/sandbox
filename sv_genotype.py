#!/usr/bin/env python

import argparse, sys
import pysam
from collections import Counter
from argparse import RawTextHelpFormatter

__author__ = "Author (email@site.com)"
__version__ = "$Revision: 0.0.1 $"
__date__ = "$Date: 2014-03-19 14:31 $"

# --------------------------------------
# define functions

def get_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="\
sv_genotype.py\n\
author: " + __author__ + "\n\
version: " + __version__ + "\n\
description: Basic python script template")
    parser.add_argument('-a', '--regionA', required=True, help='breakpoint region A (chrom:start-end)')
    parser.add_argument('-b', '--regionB', required=True, help='breakpoint region B (chrom:start-end)')
    parser.add_argument('-i', '--id', required=False, default=1, help='id for the SV')
    #parser.add_argument('-c', '--landing', required=False, default=None,  help='landing area (chrom:start-end)')
    parser.add_argument('-f', '--flank', type=int, required=False, default=20, help='min number of query bases flanking breakpoint on either side [20]')
    parser.add_argument('-l', '--readlength', type=int, required=False, default=101, help='max read length in bam file')
    parser.add_argument('-B', '--bam', type=pysam.Samfile, required=True, default=None, help='full bam file for sample')
    parser.add_argument('-S', '--splitters', type=pysam.Samfile, required=False, default=None, help='split-read bam file for sample')
    parser.add_argument('-D', '--discordants', type=pysam.Samfile, required=False, default=None, help='discordant-read bam file for sample')

    # parse the arguments
    args = parser.parse_args()

    # send back the user input
    return args

# primary function
def sv_genotype(sv_id, regionA, regionB, flank, readlength, bam, splitters, discordants):
    (chromA, pos_stringA) = regionA.split(':')
    (startA, endA) = map(int, pos_stringA.split('-'))

    (chromB, pos_stringB) = regionB.split(':')
    (startB, endB) = map(int, pos_stringB.split('-'))

    split_counter = Counter() # counts the number of splitters over each breakpoint
    for split_read in splitters.fetch(chromA, startA - 1, endA):
        if split_read.cigar[0][0] == 0: # and split_read.pnext > startB - readlength:
            split_counter[split_read.pos + split_read.cigar[0][1]] += 1

        #print split_read
        #i = split_read.pos + split_read.inferred_length
        #print split_read.pos, i, split_read.cigar


    #print split_counter[2911653]

    for posA in range(startA, endA + 1):
        ref_counter = 0
        # make sure these are actually ref
        for ref_read in bam.fetch(chromA, posA - 1, posA):
            #print ref_read
            #print ref_read.pos, ref_read.pos + ref_read.inferred_length, ref_read.aend
            if not ref_read.is_duplicate and posA - ref_read.pos >= flank and ref_read.aend - posA >= flank:
                ref_counter += 1
        print '\t'.join(map(str, (chromA, posA, ref_counter, split_counter[posA], 'A', sv_id)))

    
    # now do the same over regionB
    split_counter = Counter() # counts the number of splitters over each breakpoint
    for split_read in splitters.fetch(chromB, startB - 1, endB):
        #print split_read
        #print split_read.positions[-1]
        if split_read.cigar[-1][0] == 0: # and split_read.pnext > startB - readlength:
            split_counter[split_read.positions[-1] - split_read.cigar[-1][1] + 1] += 1


    # maybe use AlignedRead.positions to get the overlap.
    # or AlignedRead.overlap
    for posB in range(startB, endB + 1):
        ref_counter = 0
        # make sure these are actually ref
        for ref_read in bam.fetch(chromB, posB - 1, posB):
            if not ref_read.is_duplicate and posB - ref_read.pos >= flank and ref_read.aend - posB >= flank:
                ref_counter += 1
        print '\t'.join(map(str, (chromB, posB, ref_counter, split_counter[posB], 'B', sv_id)))

    #for split_read in splitters.fetch(break_chrom, break_pos - 1, break_pos):
    #    print split_read.pos + split_read.inferred_length

def calc_somthing(blocklist):
    print blocklist
    for read in blocklist:
       print read
       print

# --------------------------------------
# main function

def main():
    # parse the command line args
    args = get_args()

    # call primary function
    sv_genotype(args.id, args.regionA, args.regionB, args.flank, args.readlength, args.bam, args.splitters, args.discordants)

    # close the bam files
    args.bam.close()
    if args.splitters:
        args.splitters.close()
    if args.discordants:
        args.discordants.close()

    

# initialize the script
if __name__ == '__main__':
    sys.exit(main())
