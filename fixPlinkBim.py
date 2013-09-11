#!/usr/bin/env python

import argparse, sys
from argparse import RawTextHelpFormatter

__author__ = "Author (email@site.com)"
__version__ = "$Revision: 0.0.1 $"
__date__ = "$Date: 2013-05-09 14:31 $"

# --------------------------------------
# define functions

def get_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="\
fixPlinkBim.py\n\
author: " + __author__ + "\n\
version: " + __version__ + "\n\
description: Plink bim (extended map) files use columns 5 and 6 to\n\
delineate alleleA and alleleB at a locus. However, they are in\n\
arbitrary order with respect to the reference genome, which causes\n\
downstream issues with converting it into VCF or other standard formats.\n\
This script takes a bim file and a file with the ref base at each locus\n\
as input, and if outputs a bim with the ref allele in column 5 and the\n\
alternate allele in column 6.")
    parser.add_argument('-b', '--bimFile', type=argparse.FileType('r'), required=True, help='plink bim file to fix')
    parser.add_argument('-r', '--refFile', type=argparse.FileType('r'), required=True, help="""Line delimited reference base at each
position in the bim file. Must be the
same number of lines as the bim file with
lines in corresponding order""")

    # parse the arguments
    args = parser.parse_args()

    return args

# primary function
def fixBim(bimFile, refFile):
    while 1:
        bimLine = bimFile.readline()
        refLine = refFile.readline()
        if not bimLine or not refLine:
            break
        bimV = bimLine.rstrip().split('\t')
        refV = refLine.rstrip().split('\t')

        alleleA = bimV[4]
        alleleB = bimV[5]

        refBase = refV[0]

        # if neither of the alleles match the reference
        # genome there is problem. Print it anyway so we don't
        # mess up the plink-seq import, but kick an error and 
        # remember to remove it later.
        if refBase != alleleA and refBase != alleleB:
            sys.stderr.write("Ref base doesn't match allele (ref: %s)\n%s" % (refV[0], bimLine))
            print bimLine.rstrip()

        elif alleleA == refBase:
            print bimLine.rstrip()
        elif alleleB == refBase:
            print '\t'.join(map(str,(bimV[0], bimV[1], bimV[2], bimV[3], bimV[5], bimV[4])))

    return

# --------------------------------------
# main function

def main():
    # parse the command line args
    args = get_args()

    # call primary function
    fixBim(args.bimFile, args.refFile)

    # close the input files
    args.bimFile.close()
    args.refFile.close()

# initialize the script
if __name__ == '__main__':
    sys.exit(main())
