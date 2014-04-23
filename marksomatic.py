#!/usr/bin/env python

import argparse, sys
from argparse import RawTextHelpFormatter
import vcf as pyvcf

__author__ = "Colby Chiang (cc2qe@virginia.edu)"
__version__ = "$Revision: 0.0.1 $"
__date__ = "$Date: 2014-04-23 07:57 $"

# --------------------------------------
# define functions

def get_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="\
marksomatic.py\n\
author: " + __author__ + "\n\
version: " + __version__ + "\n\
description: Basic python script template")
    # parser.add_argument('-a', '--argA', metavar='argA', type=str, required=True, help='description of argument')
    # parser.add_argument('-b', '--argB', metavar='argB', required=False, help='description of argument B')
    # parser.add_argument('-c', '--flagC', required=False, action='store_true', help='sets flagC to true')
    parser.add_argument('-s', '--strict', required=False, action='store_true', help='Require 0 non-reference reads in normal')
    parser.add_argument('-t', '--tumor', required=True, type=str, help='Name of tumor sample(s) (comma separated list)')
    parser.add_argument('-n', '--normal', required=True, type=str, help='Name of normal sample')
    parser.add_argument('-v', '--vcf', type=argparse.FileType('rb'), default=None, help='VCF file to read. If \'-\' or absent then defaults to stdin.')

    # parse the arguments
    args = parser.parse_args()

    # if no input, check if part of pipe and if so, read stdin.
    if args.vcf == None:
        if sys.stdin.isatty():
            parser.print_help()
            exit(1)
        else:
            args.vcf = sys.stdin

    # send back the user input
    return args

# primary function
def get_somatic(strict, normal, tumor_list, vcf, outvcf):
    for rec in vcf:
        if rec.genotype(normal).gt_type == 0 and (not strict or sum(map(int, rec.genotype(normal).data.AO)) == 0):
            for tumor in tumor_list:
                if rec.genotype(tumor).is_variant:
                    rec.add_info('SOMATIC')
                    break
        outvcf.write_record(rec)
    return

# --------------------------------------
# main function

def main():
    # parse the command line args
    args = get_args()

    # format the arguments
    vcf = pyvcf.Reader(args.vcf)
    tumor_list = args.tumor.split(',')

    outvcf = pyvcf.Writer(open('/dev/stdout', 'w'), vcf)
    #outvcf = pyvcf.Writer(open('test.vcf', 'w'), vcf)
    # call primary function
    get_somatic(args.strict, args.normal, tumor_list, vcf, outvcf)

    # close the files
    args.vcf.close()
    outvcf.close()

# initialize the script
if __name__ == '__main__':
    sys.exit(main())
