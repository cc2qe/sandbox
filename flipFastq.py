#!/usr/bin/env python

import sys
from string import *

# gets the reverse complement of a DNA strand
def revcomp(dna):
    """ reverse complement of a DNA sequence """
    
    comp = dna.translate(maketrans("AGCTagct", "TCGAtcga"))
    lcomp = list(comp)
    lcomp.reverse()

    return ''.join(lcomp)


f = sys.stdin

i = 1

for l in f:
	l = l.rstrip()
	if i % 4 == 2:
		print revcomp(l)
	elif i % 4 == 0:
		print l[::-1]
	else:
		print l
	
	i += 1
