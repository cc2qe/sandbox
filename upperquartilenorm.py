#!/usr/bin/env python

import argparse, sys
from argparse import RawTextHelpFormatter
import math
import numpy as np

__author__ = "Colby Chiang (cc2qe@virginia.edu)"
__version__ = "0.0.2"
__date__ = "$Date: 2015-04-21 14:31 $"

# --------------------------------------
# define functions

def get_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="\
upperquartilenorm.py " + __version__ + "\n\
author: " + __author__ + "\n\
description: upper quartile normalize expression data")
    # parser.add_argument('-a', '--abs', action='store_true', help='take absolute values of input before calculating stats')
    parser.add_argument('-q', '--quantile', type=str, default=0.75, help='quantile for normalization [0.75]')
    parser.add_argument('-r', '--skip_rows', type=int, default=0, help='number of rows to skip [0]')
    parser.add_argument('-c', '--skip_columns', type=int, default=0, help='number of columns to skip [0]')
    parser.add_argument('data', nargs='?', type=argparse.FileType('r'), default=None, help='input data [stdin]')

    # parse the arguments
    args = parser.parse_args()

    # if no input, check if part of pipe and if so, read stdin.
    if args.data == None:
        if sys.stdin.isatty():
            parser.print_help()
            exit(1)
        else:
            args.data = sys.stdin
    return args

# primary function
def upperquartilenorm(data,
                      quantile,
                      skip_rows,
                      skip_columns):
    # store the row and column headers so we can
    # output them at the end
    d_str = np.genfromtxt(data, dtype=str, comments=None)
    row_head = d_str[:skip_rows,:]
    col_head = d_str[:,:skip_columns]

    d = d_str[skip_rows:,skip_columns:].astype(float)
    num_rows = d.shape[0]
    num_cols = d.shape[1]

    # get the (0.75) quantile for each column for
    # the non-zero rows
    quantile_val = [None] * num_cols
    non_zero_rows = d[np.any(d, axis=1)]
    for i in xrange(num_cols):
        quantile_val[i] = np.percentile(non_zero_rows[:,i], 100 * quantile)

    # the mean of the quantile values
    mean_quantile_val = sum(quantile_val) / num_cols

    # print the header
    for i in xrange(row_head.shape[0]):
        print '\t'.join(row_head[i,:])

    # print the name columns plus the actual data
    for i in xrange(num_rows):
        print '\t'.join(col_head[i + skip_rows,:].tolist() +
                        ["%0.6f" % x for x in d[i,:] / quantile_val]
                        )
    return

# --------------------------------------
# main function

def main():
    # parse the command line args
    args = get_args()

    # call primary function
    upperquartilenorm(args.data,
                      args.quantile,
                      args.skip_rows,
                      args.skip_columns)
    
    # close the file
    args.data.close()

# initialize the script
if __name__ == '__main__':
    try:
        main()
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise 
