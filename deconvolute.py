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
pythonTemplate.py\n\
author: " + __author__ + "\n\
version: " + __version__ + "\n\
description: Basic python script template")
    parser.add_argument('copyratio', type=float, help='log 2 copy ratio of tumor to normal')
    parser.add_argument('mafrac', type=float, help='minor allele fraction')
    # parser.add_argument('-c', '--flagC', required=False, action='store_true', help='sets flagC to true')
    # parser.add_argument('input', nargs='?', type=argparse.FileType('r'), default=None, help='file to read. If \'-\' or absent then defaults to stdin.')

    # parse the arguments
    args = parser.parse_args()

    # if no input, check if part of pipe and if so, read stdin.
    # if args.input == None:
    #    if sys.stdin.isatty():
    #        parser.print_help()
    #        exit(1)
    #    else:
    #        args.input = sys.stdin

    # send back the user input
    return args

# primary function
def decon(log_ratio, ma_frac):
    raw_ratio = 2 ** log_ratio

    # let R be the list of possible copy number ratios. These
    # fractions are rational numbers based on the assumption of
    # integer number of chromosomes
    R = [0.5, 3/2.0, 4/2.0]

    # let S be the list of possible minor allele fractions. These
    # fractions are rational numbers based on the assumption of
    # integer number of chromosomes
    S = [0, 1/3.0, 1/4.0]

    # now we'll assume that the tumor population is comprised of a subpopulation of normal
    # and a subpopulation of mutant cells. So we'll solve a system of equations where
    # subpop1 + subpop2 = 1 and raw_ratio = 0.5 * subpop1 + t * subpop2 where t is in s.

    for i in range(len(R)):
        r = R[i]
        s = S[i]
        
        A = np.matrix([[1, r],[0.5, s]])
        b = np.matrix([[raw_ratio],[ma_frac]])

        x = np.linalg.solve(A,b)

        for m in range(x.shape[0]):
            for n in range(x.shape[1]):
                if x[m,n] < 0:
                    x[m,n] = 0
                elif x[m,n] > 1:
                    x[m,n] = 1

        x[0,0]*0.5 + x[1,0]*s

        if abs(1-sum(x)) < 0.1:
            print 's = %s' % s
            print x[0,0], x[1,0]
            print 'sum is %s' % sum(x)

    return

# --------------------------------------
# main function

def main():
    # parse the command line args
    args = get_args()

    # store into global values
    log_ratio = args.copyratio
    ma_frac = args.mafrac
    
    # call primary function
    decon(log_ratio, ma_frac)

    # close the input file

# initialize the script
if __name__ == '__main__':
    sys.exit(main())
