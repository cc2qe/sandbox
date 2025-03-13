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
    parser.add_argument('cmd',
                        metavar='STRING',
                        type=str, default=None,
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
    parser.add_argument('-J', '--jobname',
                        metavar="STRING", dest='jobname',
                        type=str, default=None,
                        help="Jobname [random string]")
    parser.add_argument('-k', '--keep',
                        action='store_true', default=False,
                        dest='keep',
                        help="Keep temporary bash script file (deleted by default)")
    # parser.add_argument('-o', '--outfile',
    #                     metavar="STRING", dest="outfile",
    #                     type=str, default=None,
    #                     help="Send screen output (STDOUT) to file outfile.")
    # parser.add_argument('-e', '--errfile',
    #                     metavar="STRING", dest="errfile",
    #                     type=str, default=None,
    #                     help="Send errors (STDERR) to file errfile")

    # parse the arguments
    args = parser.parse_args()

    return args

# parse command line arguments to Slurm command
def generate_slurm_cmd(args, temp):
    q = ["sbatch"]
    q.extend(['-c', str(args.ncores)])
    q.extend(['-n', str(args.ntasks)])
    q.extend(['--mem=' + args.mem])
    q.extend(['-p', args.queue])
    # arglist.append(join_arg('-o', args.outfile))
    # arglist.append(join_arg('-e', args.errfile))
    q.extend(['-t', args.time])

    if args.jobname is not None:
        q.extend(['-J', args.jobname])

    # append the bash script to the end of slurm cmd
    q.append(temp.name)

    return(q)

# --------------------------------------
# main function

def main():
    # parse the command line args
    args = get_args()

    # Creating a temporary file using NamedTemporaryFile()
    with tempfile.NamedTemporaryFile(
        mode='w+t',
        dir='.',
        prefix='o2-',
        suffix='.sh',
        delete=(not args.keep)
    ) as temp:
        print('Temporary bash script file:', temp.name)

        # Writing text data into the file
        temp.write('#!/bin/bash\n')
        temp.write(args.cmd_string)
        temp.flush() # close and delete the file when out of scope

        # generate Slurm command string from arguments
        cmd = generate_slurm_cmd(args, temp)

        # print the sbatch command
        print(' '.join(cmd))

        # run the command
        p = Popen(cmd)
        p.wait() # wait until subprocess is complete
    
# initialize the script
if __name__ == '__main__':
    try:
        sys.exit(main())
    except IOError as e:
        if e.errno != 32:  # ignore SIGPIPE
            raise 
