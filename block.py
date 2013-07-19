#!/usr/bin/env python

import argparse, sys
from argparse import RawTextHelpFormatter
import numpy as np

__author__ = "Author (email@site.com)"
__version__ = "$Revision: 0.0.1 $"
__date__ = "$Date: 2013-05-09 14:31 $"

# --------------------------------------
# define functions

def get_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="\
block.py\n\
author: " + __author__ + "\n\
version: " + __version__ + "\n\
description: Basic python script template")
    #    parser.add_argument('-a', '--argA', metavar='argA', type=str, required=True, help='description of argument')
    #parser.add_argument('-b', '--argB', metavar='argB', required=False, help='description of argument B')
    #parser.add_argument('-c', '--flagC', required=False, action='store_true', help='sets flagC to true')
    parser.add_argument('input', nargs='?', type=argparse.FileType('r'), default=None, help='file to read. If \'-\' or absent then defaults to stdin.')

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
def block(file):
    prevGt = -1
    diffs = list()
    for line in file:
        v = line.rstrip().split('\t')
        gt = v[14]


        if gt != prevGt and prevGt != -1:
            print "dump data %s" % prevGt
            print "mean counts %s" % np.mean(diffs)
            print "stdev %s" % np.std(diffs)
            for i in range(len(diffs)):
                if abs(diffs[i] - np.mean(diffs)):
                    diffs[0:i] + diffs[i+1:len(diffs)]
                # omitArr = diffs[0:i] + diffs[i+1:len(diffs)]
                # print(np.std(omitArr))
            diffs = list()
        
        print line.rstrip()
        counts = v[15].split('|')
        # difference between observed and expected
        diffs.append(float(counts[2]))
        
        prevGt = gt
    
    return

# --------------------------------------
# main function

def main():
    # parse the command line args
    args = get_args()

    # store into global values
    file = args.input
    
    # call primary function
    block(file)

    # close the input file
    file.close()

# initialize the script
if __name__ == '__main__':
    sys.exit(main())
