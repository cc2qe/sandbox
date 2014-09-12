#!/usr/bin/env python

# http://www.theprojectspot.com/tutorial-post/simulated-annealing-algorithm-for-beginners/6
# http://www.psychicorigami.com/2007/06/28/tackling-the-travelling-salesman-problem-simmulated-annealing/


import argparse, sys, random, math
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
    parser.add_argument('-s', '--seed', required=False, help='random seed')
    parser.add_argument('-n', '--num_cities', type=int, required=True, help='number of cities')
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
    
def binary_flip(route, a, b):
    new_route = route[:]
    new_route[a], new_route[b] = route[b], route[a]
    return new_route

# decide if to switch
def switch(dist, new_dist, temp):
    if new_dist < dist:
        p = 1.0
    else:
        p = math.exp(-abs(new_dist - dist)/temp)

    print temp, dist, new_dist, p
    if random.random() < p:
            return True
    else:
        return False

def kirkpatrick_cooling(start_temp, alpha):
    T = start_temp
    while True:
        yield T
        T = alpha * T

# primary function
def myFunction(size):
    start_temp = 200
    alpha = 0.9999
    max_evals = 100000

    points = []
    for i in xrange(size):
        coord = [random.random(), random.random()]
        points.append(coord)

    # print points
    # print points[5], points[9]
    # print dist(points[5], points[9])

    route = range(size)
    random.shuffle(route)
    route.append(route[0])
    print route
    d = total_dist(points, route)

    print 'initial dist', d

    best_route = route
    min_dist = d

    i = 0

    cooling_schedule = kirkpatrick_cooling(start_temp, alpha)
    path_exhaust = 0 # number of temps that have iterated because of exhausted paths
    for temp in cooling_schedule:
        done = False
        tried = set()
        while not done:
            while 1:
                a,b = random.sample(xrange(1, len(route)-1),2)
                if (a,b) not in tried: break
            flip_route = binary_flip(route, a, b)
            flip_d = total_dist(points, flip_route)
            tried.add((a,b))

            print i
            print len(tried)
            if switch(d, flip_d, temp):
                route = flip_route
                d = flip_d

                # only iterate temp if accept
                i += 1
                done = True

            if d < min_dist:
                min_dist = d
                best_route = route[:]

            if len(tried) >= (size - 2) * (size - 1):
                i += 1
                path_exhaust += 1
                done = True

        if i % max_evals/50 == 0:
            f = open('routes/route_%s.txt' % i, 'w')
            for r in route:
                f.write('\t'.join(map(str, points[r])) + '\n')
            f.close()
        if path_exhaust > 5 or i >= max_evals:
            break

    f = open('routes/route_best.txt', 'w')
    for r in best_route:
        f.write('\t'.join(map(str, points[r])) + '\n')
    f.close()

    print best_route
    print min_dist
        




        
    
    return

# --------------------------------------
# main function

def main():

    # parse the command line args
    args = get_args()

    random.seed(args.seed)

    # call primary function
    myFunction(args.num_cities)

        

    # close the input file
    # args.input.close()

# initialize the script
if __name__ == '__main__':
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise
