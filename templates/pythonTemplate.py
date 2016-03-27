#!/usr/bin/env python

import argparse, sys
# import math, time, re
import gzip
# import numpy as np
# from scipy import stats
# from collections import Counter
from argparse import RawTextHelpFormatter

__author__ = "Colby Chiang (cchiang@genome.wustl.edu)"
__version__ = "$Revision: 0.0.1 $"
__date__ = "$Date: 2016-03-27 09:43 $"

# --------------------------------------
# define functions

def get_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="\
fastqtl_qvalue.py\n\
author: " + __author__ + "\n\
version: " + __version__ + "\n\
description: compute q-values from FastQTL nominal p-values (stdin) \n\
             based on the beta distribution of permutations")
    parser.add_argument('-i', '--input',
                        metavar='FILE', dest='input_path',
                        type=str, default=None,
                        help='FastQTL file of nominal p-values [stdin]')
    parser.add_argument('-p', '--permutation',
                        metavar='FILE', dest='permutation_path',
                        required=True,
                        type=str, default=None,
                        help='FastQTL file of permutation p-values [stdin]')
    parser.add_argument('-a', '--argA',
                        metavar='FLOAT', dest='argA',
                        type=float, required=False,
                        help='description of argument')
    parser.add_argument('-c', '--flagC',
                        required=False, action='store_true',
                        help='sets flagC to true')
    parser.add_argument('input', nargs='?', type=argparse.FileType('r'),
                        default=None,
                        help='file to read. If \'-\' or absent then defaults to stdin.')


    # parse the arguments
    args = parser.parse_args()

    # if no input file, check if part of pipe and if so, read stdin.
    if args.input_path == None:
        if sys.stdin.isatty():
            parser.print_help()
            exit(1)

    # send back the user input
    return args

# open file (either plaintext or zip)
def get_file(filename):
    if filename.endswith('.gz'):
        data = gzip.open(filename, 'rb')
    else:
        data = open(filename, 'r')
    return data    

# primary function
def calc_q(nominal, permutation):

    return

# --------------------------------------
# main function

def main():
    # parse the command line args
    args = get_args()

    # if no input file, check if part of pipe and if so, read stdin.
    if args.input_path == None:
        input_file = sys.stdin
    else:
        input_file = get_file(args.input_path)

    # get permutation data
    permutation_file = get_file(args.permutation_path)

    # call primary function
    calc_q(
        input_file,
        permutation_file
        )

    # close the files
    input_file.close()
    


# initialize the script
if __name__ == '__main__':
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise 
