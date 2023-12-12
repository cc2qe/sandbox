#!/usr/bin/env python

import argparse, sys
import tempfile
# import math, time, re
# import gzip
# import numpy as np
# from scipy import stats
# from collections import Counter
from argparse import RawTextHelpFormatter
from subprocess import Popen, PIPE, STDOUT

__author__ = "Colby Chiang (colby.chiang@childrens.harvard.edu)"
__version__ = "$Revision: 0.0.1 $"
__date__ = "$Date: 2023-12-12 11:33 $"

# --------------------------------------
# define functions

def get_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="\
o2run.py\n\
author: " + __author__ + "\n\
version: " + __version__ + "\n\
description: wrapper for submitted jobs on HMS O2 cluster. \n\
             See https://harvardmed.atlassian.net/wiki/spaces/O2/overview \n\
             for details.")
    parser.add_argument('--cmd',
                        metavar='STRING', dest='cmd_string',
                        type=str, default=None, required=True,
                        help='Command to run (required)')
    parser.add_argument('-t', '--time',
                        metavar="STRING", dest='time',
                        type=str, default="0-04:00",
                        help="Runtime in D-HH:MM format [0-04:00]")
    parser.add_argument('-m', '--mem',
                        metavar="STRING", dest='mem',
                        type=str, default="8G",
                        help="Memory to reserve. Preferred to included a size unit \n\
(K, M, G, T for kibibyte, mebibyte, etc.) \n\
If you don't specify a unit, requests mebibytes by default. [8G]")
    parser.add_argument('-q', '--queue',
                        metavar="STRING", dest='queue',
                        type=str, default="short",
                        help="Submission queue (short: <12h, medium: <5d, long <30d, \n\
priority, transfer, etc.) [short]")
    parser.add_argument('-n', '--tasks',
                        metavar="INT", dest='ntasks',
                        type=int, default=1,
                        help="Number of tasks [1]")
    parser.add_argument('-c', '--cores',
                        metavar="INT", dest='ncores',
                        type=int, default=1,
                        help="Number of cores per task [1]")
    # parser.add_argument('-o', '--outfile',
    #                     metavar="STRING", dest="outfile",
    #                     type=str, default=None,
    #                     help="Send screen output (STDOUT) to file outfile.")
    # parser.add_argument('-e', '--errfile',
    #                     metavar="STRING", dest="errfile",
    #                     type=str, default=None,
    #                     help="Send errors (STDERR) to file errfile")

    # parser.add_argument('-p', '--permutation',
    #                     metavar='FILE', dest='permutation_path',
    #                     required=True,
    #                     type=str, default=None,
    #                     help='FastQTL file of permutation p-values [stdin]')
    # parser.add_argument('-a', '--argA',
    #                     metavar='FLOAT', dest='argA',
    #                     type=float, required=False,
    #                     help='description of argument')
    # parser.add_argument('-c', '--flagC',
    #                     required=False, action='store_true',
    #                     help='sets flagC to true')
    # parser.add_argument('input', nargs='?', type=argparse.FileType('r'),
    #                     default=None,
    #                     help='file to read. If \'-\' or absent then defaults to stdin.')


    # parse the arguments
    args = parser.parse_args()

    # # if no input file, check if part of pipe and if so, read stdin.
    # if args.input_path == None:
    #     if sys.stdin.isatty():
    #         parser.print_help()
    #         exit(1)

    # send back the user input
    return args

# open file (either plaintext or zip)
def get_file(filename):
    if filename.endswith('.gz'):
        data = gzip.open(filename, 'rb')
    else:
        data = open(filename, 'r')
    return data    

def join_arg(flag, value):
    return (' '.join([str(flag), str(value)]))

# primary function
def generate_slurm_cmd(args, temp_file):
    q = ["sbatch"]
    q.extend(['-c', str(args.ncores)])
    q.extend(['-n', str(args.ntasks)])
    q.extend(['--mem=' + args.mem])
    q.extend(['-p', args.queue])
    # arglist.append(join_arg('-o', args.outfile))
    # arglist.append(join_arg('-e', args.errfile))
    q.extend(['-t', args.time])

    q.append(temp_file.name)


    print(q)

    # cmd = 'sbatch ' + ' '.join(arglist) + ' ' + '/home/cc514/this.sh' #tmp.name
    # return(cmd)
    return(q)

# --------------------------------------
# main function

def main():
    # parse the command line args
    args = get_args()

    # # if no input file, check if part of pipe and if so, read stdin.
    # if args.input_path == None:
    #     input_file = sys.stdin
    # else:
    #     input_file = get_file(args.input_path)

    # # get permutation data
    # permutation_file = get_file(args.permutation_path)

 
    # Creating a temporary file using NamedTemporaryFile()
    temp_file = tempfile.NamedTemporaryFile(
        mode='w+t',
        dir='.',
        prefix='o2run-',
        suffix='.sh',
        delete=False
    )
 
    print('Named file name:', temp_file.name)

    # Writing text data into the file
    temp_file.write('#!/bin/bash\n')
    temp_file.write(args.cmd_string)

    cmd = generate_slurm_cmd(args, temp_file)
    p = Popen(cmd)


    # Closing the file
    temp_file.close()

    # create temporary file for bash script
    # with tempfile.NamedTemporaryFile(dir='.', delete=False) as tmp:
    #     # call primary function
    #     cmd = generate_slurm_cmd(args, tmp)
    #     print(cmd)

    #     p = Popen(cmd, stdout=PIPE, stderr=STDOUT)

    # # close the files
    # input_file.close()
    
# initialize the script
if __name__ == '__main__':
    try:
        sys.exit(main())
    except IOError as e:
        if e.errno != 32:  # ignore SIGPIPE
            raise 
