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
    parser.add_argument('-a', '--alpha', type=float, default=1, help='weight for copy number statistic (only applies to hybrid method, default: 1)')
    parser.add_argument('-b', '--beta', type=float, default=1, help='weight for allele statistic (only applies to hybrid method, default: 1)')
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
def decon(method, seg_file, argAlpha, argBeta):
    # print the output header
    print '\t'.join(('# chrom', 'start', 'end', 'seg_id', 'copy_count', 'r', 's', 'pop1', 'pop2', 'copy_ratio', 'y_copy_ratio', 'min_frac', 'y_min_frac'))
    
    het = 0.45
    
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

        if method == 'copy':
            # let R be the list of aberrant copy number ratios. These
            # fractions are rational numbers based on the assumption of
            # integer number of chromosomes
            R = [0.5, 3/2.0, 4/2.0]

            # let S be the list of aberrant minor allele fractions. These
            # fractions are rational numbers based on the assumption of
            # integer number of chromosomes
            S = [0, 1/3.0, 1/4.0]

            # just for the output string, to convert y_best to number of copies
            copy_count = [1,3,4]

            # make blank matrices for the residuals
            x_copy_resids = [0]*len(R)
            x_min_frac_resids = [0]*len(R)
            x_sols = [None]*len(R)
            
            # iterate over each aberrant possibility
            for i in range(len(R)):
                r = R[i]
                s = S[i]

                # Solve based on copy number first (alpha is coeff for copy number, beta is
                # coeff for minor allele fraction
                alpha = 1
                beta = 0

                # A * x = B, solve for x.
                A = np.matrix([[alpha*1 + beta*het ,alpha*r + beta*s],[1,1]])
                B = np.matrix([[alpha*copy_ratio + beta*min_frac],[1]])
                x = np.linalg.solve(A,B)

                # if any of the subpops are > 1 or < 0, then peg them to 1 and 0
                for m in range(x.shape[0]):
                    for n in range(x.shape[1]):
                        if x[m,n] < 0:
                            x[m,n] = 0
                        elif x[m,n] > 1:
                            x[m,n] = 1
        
                # Now we have copy number solutions (x) and minor allele fraction solutions (y)
                # Calculate the residuals off the other metric for each
                x_copy_ratio = x[0,0]*1 + x[1,0]*r
                x_min_frac = x[0,0]*het + x[1,0]*s
                x_copy_resids[i] = abs(x_copy_ratio - copy_ratio)    # this should always be zero
                x_min_frac_resids[i] = abs(x_min_frac - min_frac)
                x_sols[i] = x
        
            # now we have the copy number and allele solutions. determine the best one
            x_best = x_min_frac_resids.index(min(x_min_frac_resids))

            # the best calculations for the copy number method
            x = x_sols[x_best]
            s = S[x_best]
            r = R[x_best]
            x_copy_ratio = x[0,0]*1 + x[1,0]*r
            x_min_frac = x[0,0]*het + x[1,0]*s
            print '\t'.join(map(str, (chrom, start, end, seg_id, copy_count[x_best], r, "%.3f" % s, "%.3f" % x[0,0], "%.3f" % x[1,0], "%.3f" % copy_ratio, "%.3f" % x_copy_ratio, "%.3f" % min_frac, "%.3f" % x_min_frac)))
            
        elif method == 'allele':
            # let R be the list of aberrant copy number ratios. These
            # fractions are rational numbers based on the assumption of
            # integer number of chromosomes
            R = [0.5, 3/2.0, 4/2.0]

            # let S be the list of aberrant minor allele fractions. These
            # fractions are rational numbers based on the assumption of
            # integer number of chromosomes
            S = [0, 1/3.0, 1/4.0]

            # just for the output string, to convert y_best to number of copies
            copy_count = [1,3,4]

            # make blank matrices for the residuals
            y_copy_resids = [0]*len(R)
            y_min_frac_resids = [0]*len(R)
            y_sols = [None]*len(R)
            
            # iterate over each aberrant possibility
            for i in range(len(R)):
                r = R[i]
                s = S[i]

                # Solve based on minor allele fraction
                alpha = 0
                beta = 1

                # A * y = B, solve for y.
                A = np.matrix([[alpha*1 + beta*het ,alpha*r + beta*s],[1,1]])
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
                y_copy_ratio = y[0,0]*1 + y[1,0]*r
                y_min_frac = y[0,0]*het + y[1,0]*s
                y_copy_resids[i] = abs(y_copy_ratio - copy_ratio)
                y_min_frac_resids[i] = abs(y_min_frac - min_frac)    # this should always be zero
                y_sols[i] = y
            
            # now we have the copy number and allele solutions. determine the best one
            y_best = y_copy_resids.index(min(y_copy_resids))
            
            # the best calculations for the allele method
            y = y_sols[y_best]
            s = S[y_best]
            r = R[y_best]
            y_copy_ratio = y[0,0]*1 + y[1,0]*r
            y_min_frac = y[0,0]*het + y[1,0]*s
            print '\t'.join(map(str, (chrom, start, end, seg_id, copy_count[y_best], r, "%.3f" % s, "%.3f" % y[0,0], "%.3f" % y[1,0], "%.3f" % copy_ratio, "%.3f" % y_copy_ratio, "%.3f" % min_frac, "%.3f" % y_min_frac)))
        
        elif method == 'hybrid':
            # let R be the list of aberrant copy number ratios. These
            # fractions are rational numbers based on the assumption of
            # integer number of chromosomes
            R = [1/2.0, 2/2.0, 3/2.0, 4/2.0, 4/2.0, 5/2.0, 5/2.0] # max chroms: 5
            # R = [1/2.0, 2/2.0, 3/2.0, 4/2.0, 4/2.0, 5/2.0, 5/2.0, 6/2.0, 6/2.0, 6/2.0]
            # R = [1/2.0, 2/2.0, 3/2.0, 3/2.0, 4/2.0, 4/2.0, 4/2.0, 5/2.0, 5/2.0, 5/2.0, 6/2.0, 6/2.0, 6/2.0, 6/2.0]

            # let S be the list of aberrant minor allele fractions. These
            # fractions are rational numbers based on the assumption of
            # integer number of chromosomes
            S = [0/1.0, 0/2.0, 1/3.0, 1/4.0, 2/4.0, 1/5.0, 2/5.0] # max chroms: 5
            # S = [0/1.0, 0/2.0, 1/3.0, 1/4.0, 2/4.0, 1/5.0, 2/5.0, 1/6.0, 2/6.0, 3/6.0]
            # S = [0/1.0, 0/2.0, 0.03/3.0, 1/3.0, 0/4.0, 1/4.0, 2/4.0, 0/5.0, 1/5.0, 2/5.0, 0/6.0, 1/6.0, 2/6.0, 3/6.0]

            # just for the output string, to convert y_best to number of copies
            copy_count = [1,2,3,4,4,5,5] # max chroms: 5
            # copy_count = [1,2,3,4,4,5,5,6,6,6]
            # copy_count = [1,2,3,3,4,4,4,5,5,5,6,6,6,6]

            # make blank matrices for the residuals
            z_copy_resids = [0]*len(R)
            z_min_frac_resids = [0]*len(R)
            z_combo_resids = [0]*len(R)
            z_sols = [None]*len(R)

            # iterate over each aberrant possibility
            for i in range(len(R)):
                r = R[i]
                s = S[i]
                # Solve based on combo
                alpha = argAlpha
                beta = argBeta
                
                # A * y = B, solve for z.
                A = np.matrix([[alpha*1 + beta*het ,alpha*r + beta*s],[1,1]])
                B = np.matrix([[alpha*copy_ratio + beta*min_frac],[1]])
                z = np.linalg.solve(A,B)
                
                # if any of the subpops are > 1 or < 0, then peg them to 1 and 0
                for m in range(z.shape[0]):
                    for n in range(z.shape[1]):
                        if z[m,n] < 0:
                            z[m,n] = 0
                        elif z[m,n] > 1:
                            z[m,n] = 1

                # Now we have copy number solutions and minor allele fraction solutions
                # Calculate the residuals off the other metric for each
                z_copy_ratio = z[0,0]*1 + z[1,0]*r
                z_min_frac = z[0,0]*het + z[1,0]*s
                z_copy_resids[i] = abs(z_copy_ratio - copy_ratio)
                z_min_frac_resids[i] = abs(z_min_frac - min_frac)
                z_combo_resids[i] = z_copy_resids[i] + z_min_frac_resids[i] # do I need to multiply these by alpha and beta????
                z_sols[i] = z

            # now we have the copy number and allele solutions. determine the best one
            z_best = z_combo_resids.index(min(z_combo_resids))

            # the best calculations for the allele method
            z = z_sols[z_best]
            s = S[z_best]
            r = R[z_best]
            z_copy_ratio = z[0,0]*1 + z[1,0]*r
            z_min_frac = z[0,0]*het + z[1,0]*s
            print '\t'.join(map(str, (chrom, start, end, seg_id, copy_count[z_best], r, "%.3f" % s, "%.3f" % z[0,0], "%.3f" % z[1,0], "%.3f" % copy_ratio, "%.3f" % z_copy_ratio, "%.3f" % min_frac, "%.3f" % z_min_frac)))      
    
    return

# --------------------------------------
# main function

def main():
    # parse the command line args
    args = get_args()

    # store into global values
    seg_file = args.segfile
    
    # call primary function
    decon(args.method[0], seg_file, args.alpha, args.beta)

    # close the input file
    seg_file.close()

# initialize the script
if __name__ == '__main__':
    sys.exit(main())
