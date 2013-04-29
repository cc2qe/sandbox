#!/usr/bin/env python

import sys, random, argparse


parser = argparse.ArgumentParser(description='Sample a pseudo-random set of a desired number of lines from an input file and return them in order.')
parser.add_argument('-l', '--lines', type=int, required=True, help='Number of lines to sample from the file.')
parser.add_argument('-n', '--input_file_length', type=int, required=False, help='Number of lines in the input file. Improves speed and allows stdin (optional).')
parser.add_argument('-s', '--seed', type=int, required=False, help='Seed for random number generator. Useful for debugging')
parser.add_argument('-i', '--input', required=False, type=argparse.FileType('r'), help='Input file. (\'-\' for stdin, but -n is required for stdin)')

args = parser.parse_args()

f = args.input
num_lines_wanted = args.lines
num_lines = args.input_file_length
seed = args.seed

# if input length not user specified, count lines in file
if num_lines == None:
    num_lines = 0
    
    for line in f:
        num_lines += 1

    # Go back to the top of the file
    f.seek(0)

# Seed randomness
if seed != None:
    random.seed(seed)
else:
    random.seed()

lines_to_get = random.sample(range(num_lines), num_lines_wanted)
lines_to_get.sort()

for current_line_num in range(num_lines):
    if current_line_num in lines_to_get:
        print f.readline().rstrip()
    else:
        f.readline()


    
