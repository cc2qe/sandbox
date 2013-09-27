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
    R = [0.5, 3/2.0, 4/2.0, 4/2.0]

    # let S be the list of possible minor allele fractions. These
    # fractions are rational numbers based on the assumption of
    # integer number of chromosomes
    S = [0, 1/3.0, 1/4.0, 2/4.0]

    # now we'll assume that the tumor population is comprised of a subpopulation of normal
    # and a subpopulation of mutant cells. So we'll solve a system of equations where
    # subpop1 + subpop2 = 1 and raw_ratio = 0.5 * subpop1 + t * subpop2 where t is in s.

    # # First try. This matches for CN and ma_frac but doesn't require the populations to
    # # sum to 1.
    # for i in range(len(R)):
    #     r = R[i]
    #     s = S[i]
    #     
    #     A = np.matrix([[1, r],[0.5, s]])
    #     b = np.matrix([[raw_ratio],[ma_frac]])
    # 
    #     x = np.linalg.solve(A,b)
    # 
    #     for m in range(x.shape[0]):
    #         for n in range(x.shape[1]):
    #             if x[m,n] < 0:
    #                 x[m,n] = 0
    #             elif x[m,n] > 1:
    #                 x[m,n] = 1
    # 
    #     x[0,0]*0.5 + x[1,0]*s
    # 
    #     if abs(1-sum(x)) < 0.1:
    #         print 's = %s' % s
    #         print x[0,0], x[1,0]
    #         print 'sum is %s' % sum(x)


    print "\nnew line with %s %s" % (raw_ratio, ma_frac)
    copy_resids = [0]*len(R)
    ma_frac_resids = [0]*len(R)
    sols = [None]*len(R)
    for i in range(len(R)):
        r = R[i]
        s = S[i]

        # first calculate based on the copy number ratio
        A = np.matrix([[1,r],[1,1]])
        b = np.matrix([[raw_ratio],[1]])

        x = np.linalg.solve(A,b)

        for m in range(x.shape[0]):
            for n in range(x.shape[1]):
                if x[m,n] < 0:
                    x[m,n] = 0
                elif x[m,n] > 1:
                    x[m,n] = 1

        # then calculate based on the minor allele fraction
        A = np.matrix([[0.5,s],[1,1]])
        b = np.matrix([[ma_frac],[1]])

        y = np.linalg.solve(A,b)
        for m in range(y.shape[0]):
            for n in range(y.shape[1]):
                if y[m,n] < 0:
                    y[m,n] = 0
                elif y[m,n] > 1:
                    y[m,n] = 1

        print "s is %s, p1=%s, p2=%s" % (s, x[0,0], x[1,0])
        m_ratio = x[0,0]*1 + x[1,0]*r
        print m_ratio, raw_ratio, m_ratio - raw_ratio

        m_ma_frac = x[0,0]*0.5 + x[1,0]*s
        print m_ma_frac, ma_frac, m_ma_frac - ma_frac

        copy_resids[i] = abs(m_ratio - raw_ratio)
        ma_frac_resids[i] = abs(m_ma_frac - ma_frac)
        sols[i] = x

    # b = copy_resids.index(min(copy_resids))
    # print copy_resids
    b = ma_frac_resids.index(min(ma_frac_resids))
    print ma_frac_resids
    
    x = sols[b]
    s = S[b]
    r = R[b]
    m_ratio = x[0,0]*1 + x[1,0]*r
    m_ma_frac = x[0,0]*0.5 + x[1,0]*s
    
    print "b ", b
    print "raw_ratio %s ; ma_frac %s " % (raw_ratio, ma_frac)
    # print "r=%s, s=%s, p1=%s, p2=%s" % (r, s, x[0,0], x[1,0])
    print "r=%s, s=%s" % (r, s)
    print "p1\t%s" % x[0,0]
    print "p2\t%s" % x[1,0]
    print m_ratio, raw_ratio, m_ratio - raw_ratio
    print m_ma_frac, ma_frac, m_ma_frac - ma_frac
    
    print 

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
