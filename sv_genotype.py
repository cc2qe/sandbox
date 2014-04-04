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
    parser.add_argument('-z', '--z', type=float, required=False, default=2, help='z-score of inner-span to be considered discordant (default: 2)')
    parser.add_argument('-B', '--bam', type=pysam.Samfile, required=True, default=None, help='full bam file for sample')
    parser.add_argument('-S', '--splitters', type=pysam.Samfile, required=False, default=None, help='split-read bam file for sample')
    parser.add_argument('-D', '--discordants', type=pysam.Samfile, required=False, default=None, help='discordant-read bam file for sample')

    # parse the arguments
    args = parser.parse_args()

    # send back the user input
    return args

# primary function
def sv_genotype(sv_id, regionA, regionB, flank, readlength, z, bam, splitters, discordants):
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
        print '\t'.join(map(str, (chromA, posA, ref_counter, split_counter[posA], 'spl_A', sv_id)))

    
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
        print '\t'.join(map(str, (chromB, posB, ref_counter, split_counter[posB], 'spl_B', sv_id)))

    #for split_read in splitters.fetch(break_chrom, break_pos - 1, break_pos):
    #    print split_read.pos + split_read.inferred_length

    # now do the spanning coverage over region A
    mean_ispan = 118
    sd_ispan = 80
    mean_ospan = mean_ispan + 2 * readlength
    sd_ospan = 80

    ref_span_counter = Counter()
    disc_span_counter = Counter()
    for read in bam.fetch(chromA, startA - 1 - (mean_ospan + sd_ospan * z), endA + (mean_ospan + sd_ospan * z)):
        if not read.is_reverse and read.mate_is_reverse and not read.is_secondary:
            if ispan(read, readlength) > 0 and ispan(read, readlength) < mean_ispan + sd_ispan * z:
                #print read
                #print ispan(read, readlength)
                for posA in range(read.positions[-1] + 1, read.pnext + 1):
                    ref_span_counter[posA] += 1

                    # check mate on same chrom                    
            elif ispan(read, readlength) > 0 and read.pnext + 1 + readlength > startB and read.pnext + 1 < endB and read.cigar[-1][0] == 0:
                #print read
                #print read.cigar[-1][0]
                #print read.alen, read.qlen, read.rlen
                for posA in range(read.positions[-1] + 1, read.pnext + 1):
                    disc_span_counter[posA] += 1

    '''
    disc_span_counter = Counter()
    for disc_read in discordants.fetch(chromA, startA - 1 - (mean_ospan + sd_ospan * z), endA + 1 + (mean_ospan + sd_ospan * z)):
        if not disc_read.is_duplicate and not disc_read.is_reverse and not disc_read.is_secondary and disc_read.pnext >= startB and disc_read.pnext < endB + 1 + (mean_ospan + sd_ospan * z):
            for posA in range(disc_read.positions[-1], disc_read.pnext):
                disc_span_counter[posA] += 1

    ref_span_counter = Counter()
    for ref_read in bam.fetch(chromA, startA - 1 - (mean_ospan + sd_ospan * z), endA + (mean_ospan + sd_ospan * z) + 1):
        if not ref_read.is_duplicate and ref_read.is_proper_pair and not ref_read.is_secondary and not ref_read.is_reverse and ref_read.pnext - (ref_read.positions[-1] + 1) > 0:
            print ref_read
            #print ref_read.positions[-1], ref_read.pnext, ref_read.pnext - (ref_read.positions[-1] + 1), ref_read.tlen
            for posA in range(ref_read.positions[-1], ref_read.pnext):
                ref_span_counter[posA] +=1
    '''
    for posA in range(startA, endA + 1):
        print '\t'.join(map(str, (chromA, posA, ref_span_counter[posA], disc_span_counter[posA], 'disc_A', sv_id)))
    
    # do the same over region B
    ref_span_counter = Counter()
    disc_span_counter = Counter()
    for read in bam.fetch(chromB, startB - 1 - (mean_ospan + sd_ospan * z), endB + (mean_ospan + sd_ospan * z)):
        if read.is_reverse and not read.mate_is_reverse and not read.is_secondary:
            if ispan(read, readlength) > 0 and ispan(read, readlength) < mean_ispan + sd_ispan * z:
                for posB in range(read.pnext + readlength + 1, read.positions[0] + 1):
                    ref_span_counter[posB] += 1
            
            # check mate on same chrom                    
            elif ispan(read, readlength) > 0 and read.pnext + 1 + readlength > startA and read.pnext + 1 < endA and read.cigar[0][0] == 0:
                for posB in range(read.pnext + readlength + 1, read.positions[0] + 1):
                    disc_span_counter[posB] += 1

    '''
    disc_span_counter = Counter()
    for disc_read in discordants.fetch(chromB, startB - 1 - (mean_ospan + sd_ospan * z), endB + (mean_ospan + sd_ospan * z)):
        # NOTE: below only works for deletions
        if not disc_read.is_duplicate and disc_read.is_reverse and not disc_read.is_secondary and disc_read.pnext + 100 < endA + 1:
            #print disc_read
            for posB in range(disc_read.pnext + readlength, disc_read.positions[0]):
                disc_span_counter[posB] += 1

    # print disc_span_counter

    ref_span_counter = Counter()
    for ref_read in bam.fetch(chromB, startB - 1 - (mean_ospan + sd_ospan * z), endB + (mean_ospan + sd_ospan * z) + 1):
        # NOTE: BELOW ONLY WORKS FOR DELETIONS, NOT WEIRDLY ORIENTED STUFF
        if not ref_read.is_duplicate and ref_read.is_proper_pair and not ref_read.is_secondary and ref_read.is_reverse and ref_read.pnext + readlength - (ref_read.positions[0] + 1) < 0:
            #print ref_read
            #print ref_read.pnext + readlength, ref_read.positions[0]
            for posB in range(ref_read.pnext + readlength, ref_read.positions[0]):
                #print posB
                ref_span_counter[posB] +=1
    '''

    for posB in range(startB, endB + 1):
        print '\t'.join(map(str, (chromB, posB, ref_span_counter[posB], disc_span_counter[posB], 'disc_B', sv_id)))


def ispan(read, readlength):
    if not read.is_reverse:
        # these are (a,b] intervals
        left = read.positions[-1] + 1
        right = read.pnext
    else:
        left = read.pnext + readlength
        right = read.positions[0]
    return right - left


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
    sv_genotype(args.id, args.regionA, args.regionB, args.flank, args.readlength, args.z, args.bam, args.splitters, args.discordants)

    # close the bam files
    args.bam.close()
    if args.splitters:
        args.splitters.close()
    if args.discordants:
        args.discordants.close()

    

# initialize the script
if __name__ == '__main__':
    sys.exit(main())
