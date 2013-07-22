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
    obs = list()
    exp = list()
    diffs = list()

    ref1 = list()
    alt1 = list()
    ref2 = list()
    alt2 = list()

    # print header
    print '\t'.join(("# gt", "mean", "stdev", "numPairs", "meanObs", "meanExp", "meanRef1", "meanRef2"))
    for line in file:
        v = line.rstrip().split('\t')
        gt = v[14]

        if gt != prevGt and prevGt != -1:
            for d in diffs:
                if abs(d - np.mean(diffs) > 10):
                    index = diffs.index(d)
                    diffs.remove(d)
                    obs.pop(index)
                    exp.pop(index)
                    ref1.pop(index)
                    alt1.pop(index)
                    ref2.pop(index)
                    alt2.pop(index)
            refMaj1 = list()
            refMaj2 = list()
            for i in range(len(ref1)):
                refMaj1.append( ref1[i] - alt1[i] )
                refMaj2.append( ref2[i] - alt2[i] )
                    
                # omitArr = diffs[0:i] + diffs[i+1:len(diffs)]
                # print(np.std(omitArr))
            print '%s\t%.2f\t%.2f\t%s\t%.2f\t%.2f\t%.2f\t%.2f' % (prevGt,np.mean(diffs),np.std(diffs),len(diffs),np.mean(obs),np.mean(exp),np.mean(refMaj1),np.mean(refMaj2))

            diffs = list()
            obs = list()
            exp = list()
            ref1 = list()
            alt1 = list()
            ref2 = list()
            alt2 = list()
        
        # print line.rstrip()
        counts = v[15].split('|')
        # difference between observed and expected
        diffs.append(float(counts[2]))
        obs.append(float(counts[0]))
        exp.append(float(counts[1]))

        ref1.append(float(v[4]))
        alt1.append(float(v[6]))
        ref2.append(float(v[11]))
        alt2.append(float(v[13]))
        
        prevGt = gt
    
    for d in diffs:
        if abs(d - np.mean(diffs) > 10):
            index = diffs.index(d)
            diffs.remove(d)
            obs.pop(index)
            exp.pop(index)
            ref1.pop(index)
            alt1.pop(index)
            ref2.pop(index)
            alt2.pop(index)
            
        for i in range(len(ref1)):
            refMaj1.append( ref1[i] - alt1[i] )
            refMaj2.append( ref2[i] - alt2[i] )

            # omitArr = diffs[0:i] + diffs[i+1:len(diffs)]
            # print(np.std(omitArr))
    print '%s\t%.2f\t%.2f\t%s\t%.2f\t%.2f\t%.2f\t%.2f' % (prevGt,np.mean(diffs),np.std(diffs),len(diffs),np.mean(obs),np.mean(exp),np.mean(refMaj1),np.mean(refMaj2))
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
