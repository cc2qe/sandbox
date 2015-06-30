#!/usr/bin/env python

import argparse, sys
from argparse import RawTextHelpFormatter

__author__ = "Colby Chiang (cc2qe@virginia.edu)"
__version__ = "$Revision: 0.0.1 $"
__date__ = "$Date: 2015-02-16 10:46 $"

# --------------------------------------
# define functions

def get_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="\
clusterBedpe.py\n\
author: " + __author__ + "\n\
version: " + __version__ + "\n\
description: cluster a BEDPE file by position. input must be sorted")
    parser.add_argument('-d', '--distance',
                        metavar='INT', type=int,
                        dest='max_distance',
                        required=False,
                        default=100000,
                        help='max separation distance (bp) of adjacent loci in cluster [100000]')
    parser.add_argument('-m', '--min_cluster_size',
                        metavar='INT', type=int,
                        dest='min_cluster_size',
                        required=False,
                        default=2,
                        help='minimum cluster size to report feature [2]')
    parser.add_argument('-r', '--min_reciprocal',
                        metavar='INT', type=int,
                        dest='min_reciprocal',
                        required=False,
                        default=0,
                        help='min unique loci on each partner to report feature [0]')
    parser.add_argument('-u', '--union',
                        required=False,
                        action='store_true',
                        help='merge coordinates into their union')
    parser.add_argument('-c', '--cols',
                        metavar='STR', type=str,
                        required=False,
                        help='specify columns from the input file to operate upon (see -o option, below).\nMultiple columns can be specified in a comma-delimited list.')
    parser.add_argument('-o', '--oper',
                        metavar='STR', type=str,
                        required=False,
                        help='specify the operation that should be applied to -c.\nValid operations:\n  collapse (i.e., print a delimited list (duplicates allowed))')

    parser.add_argument('input', nargs='?',
                        type=argparse.FileType('r'),
                        default=None,
                        help='BEDPE file to read. If \'-\' or absent then defaults to stdin.')

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

class Cluster(object):
    def __init__(self):
        self.id = None
        self.elements = []

        self.chrom_a = None
        self.min_a = float("inf")
        self.max_a = 0

        self.chrom_b = None
        self.min_b = float("inf")
        self.max_b = 0

        self.size = 0
        self.uniq_a = set()
        self.uniq_b = set()

    # check whether a bedpe object is addable to this
    # cluster given the max_distance
    def can_add(self, bedpe, max_distance):
        if self.size == 0:
            return True

        if (self.chrom_a != bedpe.chrom_a
            or self.min_a - max_distance > bedpe.end_a
            or self.max_a + max_distance < bedpe.start_a):
            # print "a", self.chrom_a, self.min_a, self.max_a
            return False

        if (self.chrom_b != bedpe.chrom_b
            or self.min_b - max_distance > bedpe.end_b
            or self.max_b + max_distance < bedpe.start_b):
        #     # print "b"
            return False

        else:
            return True

    def add(self, bedpe):
        self.chrom_a = bedpe.chrom_a
        self.min_a = min(self.min_a, bedpe.start_a)
        self.max_a = max(self.max_a, bedpe.end_a)

        self.chrom_b = bedpe.chrom_b
        self.min_b = min(self.min_b, bedpe.start_b)
        self.max_b = max(self.max_b, bedpe.end_b)

        self.elements.append(bedpe)
        self.uniq_a.add((bedpe.chrom_a,
                         bedpe.start_a,
                         bedpe.end_a))
        self.uniq_b.add((bedpe.chrom_b,
                         bedpe.start_b,
                         bedpe.end_b))
        self.size += 1

    def get_cluster_string(self, cols, oper, union):
        summary = []

        if union:
            for c in cols:
                summary.append(','.join(map(str, [s.vec[c - 1] for s in self.elements])))

            c_string = "\t".join(map(str,
                                     [self.chrom_a,
                                      self.min_a,
                                      self.max_a,
                                      self.chrom_b,
                                      self.min_b,
                                      self.max_b,
                                      self.id,
                                      '\t'.join(summary)]
                                     ))
        else:
            c_string = '\n'.join('\t'.join(s.vec[:6] + [str(self.id)] + [s.vec[c - 1] for c in cols]) for s in self.elements)
        return c_string

class Bedpe(object):
    def __init__(self, line):
        v = line.rstrip().split("\t")
        self.chrom_a = v[0]
        self.start_a = int(v[1])
        self.end_a = int(v[2])
        self.chrom_b = v[3]
        self.start_b = int(v[4])
        self.end_b = int(v[5])
        self.id = v[6]
        self.vec = v

# prints and removes clusters from cluster_list that are beyond
# distance window
def prune(cluster_list, bedpe, max_distance, min_cluster_size, min_reciprocal, cols, oper, union):
    new_cluster_list = []
    global cluster_counter

    for cluster in cluster_list:
        # cluster is beyond updatable window:
        if (bedpe is None
            or cluster.chrom_a != bedpe.chrom_a
            or cluster.min_a - max_distance > bedpe.end_a
            or cluster.max_a + max_distance < bedpe.start_a):

            # print the cluster if eligible
            if (cluster.size >= min_cluster_size
                and len(cluster.uniq_a) > min_reciprocal
                and len(cluster.uniq_b) > min_reciprocal):
                cluster_counter += 1
                cluster.id = cluster_counter
                print cluster.get_cluster_string(cols, oper, union)

        # cluster is still within updatable window,
        # leave it in the cluster list
        else:
            new_cluster_list.append(cluster)

    return new_cluster_list

# primary function
def cluster_bedpe(in_file, max_distance, min_cluster_size, min_reciprocal, cols, oper, union):
    # line number
    line_counter = 0
    # the number of clusters that have been output
    global cluster_counter
    cluster_counter = 0

    # a list of clusters in the buffer
    cluster_list = []
    for line in in_file:
        if line.startswith('#'):
            continue
        line_counter += 1
        bedpe = Bedpe(line)

        matched_to_cluster = False
        for cluster in cluster_list:
            if cluster.can_add(bedpe, max_distance):
                cluster.add(bedpe)
                matched_to_cluster = True
                break

        if not matched_to_cluster:
            new_cluster = Cluster()
            new_cluster.add(bedpe)
            cluster_list.append(new_cluster)

        # prune and print eligible clusters
        if line_counter % 1000 == 0:
            cluster_list = prune(cluster_list,
                                 bedpe,
                                 max_distance,
                                 min_cluster_size,
                                 min_reciprocal,
                                 cols,
                                 oper,
                                 union)

    # prune the cluster and print eligible
    # features
    cluster_list = prune(cluster_list,
                         None,
                         max_distance,
                         min_cluster_size,
                         min_reciprocal,
                         cols,
                         oper,
                         union)

    # close the input file
    in_file.close()
    return

# --------------------------------------
# main function

def main():
    # parse the command line args
    args = get_args()

    if args.cols is not None:
        cols = map(int, args.cols.split(','))
    else:
        cols = []

    if args.oper is not None:
        oper = args.oper.split(',')
    else:
        oper = ['collapse'] * len(cols)

    # call primary function
    cluster_bedpe(args.input, args.max_distance, args.min_cluster_size, args.min_reciprocal, cols, oper, args.union)

# initialize the script
if __name__ == '__main__':
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise
