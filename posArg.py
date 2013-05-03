#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser(description='find out argparse formatting')

parser.add_argument('arg1', help='this is a positional arg')

args = parser.parse_args()

