#!/usr/bin/python

import sys, argparse

f = sys.stdin

parser = argparse.ArgumentParser(description='Crawl a trio vcf for haplotype data. Trio must be sorted by chrom and position, with a sample identifier in column 1.')


args = parser.parse_args()
