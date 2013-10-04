#!/usr/bin/env python

import argparse, sys
import numpy as np
from argparse import RawTextHelpFormatter

__author__ = "Author (email@site.com)"
__version__ = "$Revision: 0.0.1 $"
__date__ = "$Date: 2013-05-09 14:31 $"

# --------------------------------------
# define functions

def get_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="\
deconvolute.py\n\
author: " + __author__ + "\n\
version: " + __version__ + "\n\
description: separate heterogenous tumor sample into subpopulations")
    parser.add_argument('-m', '--method', nargs=1, type=str, choices=['allele', 'copy', 'hybrid'], default='allele', help='method to run')
    parser.add_argument('segfile', nargs='?', type=argparse.FileType('r'), default=None, help='segmentation file. If \'-\' or absent then defaults to stdin.')

    # parse the arguments
    args = parser.parse_args()

    # if no input, check if part of pipe and if so, read stdin.
    if args.segfile == None:
        if sys.stdin.isatty():
            parser.print_help()
            exit(1)
        else:
            args.input = sys.stdin

    # send back the user input
    return args

# primary function
def decon(method, seg_file):
    # print the output header
    print '\t'.join(('# chrom', 'start', 'end', 'seg_id', 'copy_count', 'r', 's', 'pop1', 'pop2', 'copy_ratio', 'y_copy_ratio', 'min_frac', 'y_min_frac'))

    # parse the file (tab delimited)
    # chrom, start, end, id, num_probes, log2copyratio, maj_allele_reads, min_allele_reads, min_allel_fraction
    # 1 3218610 42108999 1 20226 -0.2905 34476 53299 0.432432
    for l in seg_file:
        v = l.rstrip().split('\t')
        chrom = v[0]
        (start, end, seg_id, num_probes, log_copy, maj_reads, min_reads, min_frac) =  v[1:]
        (start, end, seg_id, num_probes, maj_reads, min_reads) = map(int, (start, end, seg_id, num_probes,maj_reads, min_reads))
        (log_copy, min_frac) = map(float, (log_copy, min_frac))
        copy_ratio = 2 ** log_copy

        # let R be the list of aberrant copy number ratios. These
        # fractions are rational numbers based on the assumption of
        # integer number of chromosomes
        R = [0.5, 3/2.0, 4/2.0]
        
        # let S be the list of aberrant minor allele fractions. These
        # fractions are rational numbers based on the assumption of
        # integer number of chromosomes
        S = [0, 1/3.0, 1/4.0]
    
        # make blank matrices for the residuals
        x_copy_resids = [0]*len(R)
        y_copy_resids = [0]*len(R)
        x_min_frac_resids = [0]*len(R)
        y_min_frac_resids = [0]*len(R)
        x_sols = [None]*len(R)
        y_sols = [None]*len(R)

        # iterate over each aberrant possibility
        for i in range(len(R)):
            r = R[i]
            s = S[i]

            # Solve based on copy number first (alpha is coeff for copy number, beta is
            # coeff for minor allele fraction
            alpha = 1
            beta = 0

            # A * x = B, solve for x.
            A = np.matrix([[alpha*1 + beta*0.5 ,alpha*r + beta*s],[1,1]])
            B = np.matrix([[alpha*copy_ratio + beta*min_frac],[1]])
            x = np.linalg.solve(A,B)
            
            # if any of the subpops are > 1 or < 0, then peg them to 1 and 0
            for m in range(x.shape[0]):
                for n in range(x.shape[1]):
                    if x[m,n] < 0:
                        x[m,n] = 0
                    elif x[m,n] > 1:
                        x[m,n] = 1

            # Solve based on minor allele fraction
            alpha = 0
            beta = 1

            # A * y = B, solve for y.
            A = np.matrix([[alpha*1 + beta*0.5 ,alpha*r + beta*s],[1,1]])
            B = np.matrix([[alpha*copy_ratio + beta*min_frac],[1]])
            y = np.linalg.solve(A,B)

            # if any of the subpops are > 1 or < 0, then peg them to 1 and 0
            for m in range(y.shape[0]):
                for n in range(y.shape[1]):
                    if y[m,n] < 0:
                        y[m,n] = 0
                    elif y[m,n] > 1:
                        y[m,n] = 1

            # Now we have copy number solutions (x) and minor allele fraction solutions (y)
            # Calculate the residuals off the other metric for each
            x_copy_ratio = x[0,0]*1 + x[1,0]*r
            x_min_frac = x[0,0]*0.5 + x[1,0]*s
            x_copy_resids[i] = abs(x_copy_ratio - copy_ratio)    # this should always be zero
            x_min_frac_resids[i] = abs(x_min_frac - min_frac)
            x_sols[i] = x

            # do the same for the minor allele fraction solutions
            y_copy_ratio = y[0,0]*1 + y[1,0]*r
            y_min_frac = y[0,0]*0.5 + y[1,0]*s
            y_copy_resids[i] = abs(y_copy_ratio - copy_ratio)
            y_min_frac_resids[i] = abs(y_min_frac - min_frac)    # this should always be zero
            y_sols[i] = y
            
        # now we have the copy number and allele solutions. determine the best one
        x_best = x_min_frac_resids.index(min(x_min_frac_resids))
        y_best = y_copy_resids.index(min(y_copy_resids))


        y = y_sols[y_best]
        s = S[y_best]
        r = R[y_best]
        y_copy_ratio = y[0,0]*1 + y[1,0]*r
        y_min_frac = y[0,0]*0.5 + y[1,0]*s

        # just for the output string, to convert y_best to number of copies
        copy_count = [1,3,4]

        print '\t'.join(map(str, (chrom, start, end, seg_id, copy_count[y_best], r, "%.3f" % s, "%.3f" % y[0,0], "%.3f" % y[1,0], "%.3f" % copy_ratio, "%.3f" % y_copy_ratio, "%.3f" % min_frac, "%.3f" % y_min_frac)))
    
    return

# --------------------------------------
# main function

def main():
    # parse the command line args
    args = get_args()

    # store into global values
    seg_file = args.segfile
    
    # call primary function
    decon(args.method, seg_file)

    # close the input file
    seg_file.close()

# initialize the script
if __name__ == '__main__':
    sys.exit(main())
