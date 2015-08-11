#!/usr/bin/env python

import argparse, sys
from argparse import RawTextHelpFormatter

__author__ = "Colby Chiang (cchiang@genome.wustl.edu)"
__version__ = "$Revision: 0.0.1 $"
__date__ = "$Date: 2015-08-11 10:59 $"

# --------------------------------------
# define functions

def get_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="\
column_select.py\n\
author: " + __author__ + "\n\
version: " + __version__ + "\n\
description: select columns from a file by header names")
    parser.add_argument('-c', '--col', metavar='FILE', required=True, type=argparse.FileType('r'), help='list of column headers to extract')
    parser.add_argument('-l', '--leading', metavar='INT', required=False, type=int, default=0, help='number of leading columns to print [0]')
    parser.add_argument('-p', '--pass', metavar='STR', dest='pass_prefix', required=False, default=None, help='prefix for comment lines in INPUT to pass unfiltered')
    parser.add_argument('-m', '--missing', metavar='STR', dest='missing_fill', type=str, required=False, default=None, help="fill missing columns with string (e.g.: NA)")
    parser.add_argument('input', nargs='?', type=argparse.FileType('r'), default=None, help='phenotype file')

    # parse the arguments
    args = parser.parse_args()

    # if no input, check if part of pipe and if so, read stdin.
    if args.input == None:
        if sys.stdin.isatty():
            parser.print_help()
            exit(1)
        else:
            args.input = sys.stdin

    # send back the user input
    return args

# primary function
def extract_cols(col, lead_cols, pass_prefix, missing_fill, source):
    # get_columns = range(lead_cols)
    select = []
    for line in col:
        select.append(line.rstrip())
    
    in_header = True
    for line in source:
        if pass_prefix is not None and  line.startswith(pass_prefix):
            print line.rstrip()
            continue
        v = line.rstrip().split('\t')
        if in_header:
            column_map = {c: v.index(c) for c in v}
            if missing_fill is None:
                get_columns = range(lead_cols) + [column_map[x] for x in select if x in column_map]
                header_v = [v[x] for x in range(lead_cols) + [column_map[x] for x in select if x in column_map]]
            else:
                get_columns = range(lead_cols) + [column_map[x] if x in column_map else None for x in select ]
                header_v = [v[x] for x in range(lead_cols)] + select

            print '\t'.join(header_v)
            in_header = False
        else:
            print '\t'.join(v[x] if x is not None else missing_fill for x in get_columns)

    source.close()
    col.close()    
    return

# --------------------------------------
# main function

def main():
    # parse the command line args
    args = get_args()

    # call primary function
    extract_cols(args.col, args.leading, args.pass_prefix, args.missing_fill, args.input)

# initialize the script
if __name__ == '__main__':
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise
