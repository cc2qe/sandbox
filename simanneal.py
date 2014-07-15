#!/usr/bin/env python

import argparse, sys, random
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
    # parser.add_argument('-a', '--argA', metavar='argA', type=str, required=True, help='description of argument')
    # parser.add_argument('-b', '--argB', metavar='argB', required=False, help='description of argument B')
    # parser.add_argument('-c', '--flagC', required=False, action='store_true', help='sets flagC to true')
    # parser.add_argument('input', nargs='?', type=argparse.FileType('r'), default=None, help='file to read. If \'-\' or absent then defaults to stdin.')

    # parse the arguments
    args = parser.parse_args()

    # if no input, check if part of pipe and if so, read stdin.
    # if args.input == None:
    #     if sys.stdin.isatty():
    #         parser.print_help()
    #         exit(1)
    #     else:
    #         args.input = sys.stdin

    # send back the user input
    return args

# return the euclidian distance between point a and point b
def dist(a, b):
    # if not the same vector length then error
    if len(a) != len(b):
        return -1

    sq_sum = 0
    for i in xrange(len(a)):
        sq_sum += (a[i] - b[i])**2
    return sq_sum**(0.5)

def total_dist(points, route):
    d = 0
    a = points[route[0]]
    for r in xrange(1,len(route)):
        b = points[route[r]]
        d += dist(a,b)
        a = b
    return d
    

# primary function
def myFunction():

    size = 20
    flips_per_epoch = 100

    points = []
    for i in xrange(size):
        coord = [random.random(), random.random()]
        points.append(coord)

    # print points[5], points[9]
    # print dist(points[5], points[9])

    min_dist = float('inf')
    best_route = None

    for i in xrange(10000):
        route = range(size)
        random.shuffle(route)
        route.append(route[0])

        d = total_dist(points, route)
        if d < min_dist:
            best_route = route
            min_dist = d


    print best_route
    print min_dist


    f = open('route.txt', 'w')
    for r in best_route:
        f.write('\t'.join(map(str, points[r])) + '\n')

    f.close()
        

    # print route


        

    
    return

# --------------------------------------
# main function

def main():
    # parse the command line args
    args = get_args()



    # call primary function
    myFunction()

        

    # close the input file
    # args.input.close()

# initialize the script
if __name__ == '__main__':
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise
