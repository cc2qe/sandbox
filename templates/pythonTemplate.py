#!/usr/bin/env python

import argparse, sys

# --------------------------------------
# define functions

def myFunction(argA, argB, flagC, file):
    for line in file:
        print line.rstrip()
    
    return

# --------------------------------------
# argument parsing

def main():
    parser = argparse.ArgumentParser(description="Basic python script template")
    parser.add_argument('-a', '--argA', metavar='argA', type=str, required=True, help='description of argument')
    parser.add_argument('-b', '--argB', metavar='argB', required=False, help='description of argument B')
    parser.add_argument('-c', '--flagC', required=False, action='store_true', help='sets flagC to true')
    parser.add_argument('file', nargs='?', type=argparse.FileType('r'), default=sys.stdin, help='file to read. If \'-\' or absent then defaults to stdin.')
    
    # parse the arguments
    args = parser.parse_args()

    # store into global values
    argA = args.argA
    argB = args.argB
    flagC = args.flagC
    file = args.file
    
    myFunction(argA, argB, flagC, file)
    file.close()

# initialize the script
if __name__ == '__main__':
    sys.exit(main())
