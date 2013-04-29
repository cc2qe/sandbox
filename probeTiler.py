#!/usr/bin/env python

"""
designProbes

This script generates probes of a given length and
tiling density for use in capture experiments.

"""

import shlex, subprocess, sys, os, argparse, re, math, uuid
from string import *

__author__ = "Colby Chiang (chiang@chgr.mgh.harvard.edu)"
__version__ = "$Revision: 0.0.3 $"
__date__ = "$Date: 2011/06/22 13:22 $"

############# FUNCTIONS ###################

#runs a command in the bash shell
def bash(cmd,cwd=None):
	# print(cmd)	# diagnostic
	retVal = subprocess.Popen(cmd, shell=True, \
		stdout=subprocess.PIPE, cwd=cwd).stdout.read().strip('\n').split('\n')
	if retVal==['']:
		return(0)
	else:
		return(retVal)

# gets the reverse complement of a DNA strand
def revcomp(dna):
    """ reverse complement of a DNA sequence """
    
    comp = dna.translate(maketrans("AGCTagct", "TCGAtcga"))
    lcomp = list(comp)
    lcomp.reverse()

    return ''.join(lcomp)

############# END FUNCTIONS ################

parser = argparse.ArgumentParser(description='Design capture probes over multiple genomic regions.')
parser.add_argument('-g', '--genome', required=True, help='Genome 2bit file')
parser.add_argument('-l', '--length', required=True, help='Length of each probe (bp)')
parser.add_argument('-s', '--separation', help='Separation distance between probes. (Either this or probeNumber is required)')
parser.add_argument('-n', '--probeNumber', help='total number of probes to tile to.  (Either this or separation is required)')
parser.add_argument('-str', '--strand', choices=['plus', 'minus', 'both'], default='plus', help='Tile over the plus strand, the minus strand or both strands')
parser.add_argument( '-c', '--coordinates', help='Chromosomal position for region of interest. Coordinates are zero-based. Must be of the form chr6:5302012-5703254')
parser.add_argument('-bed', '--bedFile', help='BED file of regions')
#parser.add_argument('-o', '--output', required=True, help='Output prefix for probe files')


args = parser.parse_args()

genome = args.genome
coords = args.coordinates
probeLength = int(args.length)
strand = args.strand
refGenome = args.genome	# 2bit fasta file (downloaded from here http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/)

myID = uuid.uuid1().hex		# session id
coordList = []

if coords != None and args.bedFile != None:
	parser.error('input must be either BED file or coordinates, but not both')

elif coords != None:
	# Input control for bad values
	if re.search(r'.*:.*-.*', coords) == None:
		parser.error('coordinates must be of the form \'chr#:start-end\'')
		
	coords = re.sub(r'[,\s]', '', coords)		# strip whitespace and commas
	coords = re.split('[:-]', coords)			# splits chr5:2039-23800 into array ['chr5', '2039', '23800']
	
	chr = coords[0].lower()						# escape characters, convert to lowercase
	
	if chr[:3] != 'chr':
		parser.error('could not understand chromosome \'%s\'' % chr)
	
	try:
		n = int(chr[3:])
		if n > 22 or n < 1:
			parser.error('could not understand chromosome \'%s\'\nchromosomes must be of the form \'chr1-22\', \'chrX\', \'chrY\', or \'chrM\'' % chr)
	except ValueError:
		if re.search(chr[3:], '[xym]'):
			chr = chr[:3] + chr[3:].upper()
		else:
			parser.error('could not understand chromosome \'%s\'\nchromosomes must be of the form \'chr1-22\', \'chrX\', \'chrY\', or \'chrM\'' % chr)
	
	start = coords[1]
	try:
		start = int(start)
		if start < 0: parser.error('start position	 cannot be less than 0')
	except ValueError:
		parser.error('start position must be integer greater than 0')
	
	end = coords[2]
	try:
		end = int(end)
		if end < 0: parser.error('end position cannot be less than 0')
	except ValueError:
		parser.error('end position must be integer greater than 0')
	if start > end:
		parser.error('end position must be greater than start position')
	
	coordList.append([chr, start, end])

elif args.bedFile != None:
	bedFile = open(args.bedFile, 'r')
	
	for l in bedFile:
		bedData = l.rstrip().split('\t')
		bedData[1] = int(bedData[1])
		bedData[2] = int(bedData[2])
		coordList.append(bedData)
	
	bedFile.close()

fullLength = 0
for r in coordList:
	fullLength = fullLength + r[2] - r[1]

if args.separation != None:
	stepSize = int(args.separation)
elif args.probeNumber != None:
	probeNumber = int(args.probeNumber)
	stepSize = (float(fullLength) - (probeLength)*len(coordList)) / probeNumber
else:
	parser.error('either --separation or --probeNumber is required')

carryOver = 0
counter = 1
for region in coordList:
	chr = region[0]
	start = region[1]
	end = region[2]

	cmd = "twoBitToFa %s %s_region.fa -seq=%s -start=%s -end=%s" % (refGenome, myID, chr, start,end)
	bash(cmd)
	
	cmd = "tail -n +2 %s_region.fa | tr -d '\n' > %s_region2.fa; mv %s_region2.fa %s_region.fa" % (myID, myID, myID, myID)
	bash(cmd)
	
	fasta = open(myID + "_region.fa", 'r')
	
	seq = fasta.read()
	seq = seq.upper()
	
	i = carryOver
	j = 0
	while round(i,4) < (len(seq) - probeLength):	# have to round here to protect against a tiny bit of rounding error
		roundCoord = int(round(i))
		
		if strand == 'plus':
			probeSeq = seq[roundCoord:roundCoord+probeLength]
		elif strand == 'minus':
			probeSeq = revcomp(seq[roundCoord:roundCoord+probeLength])
		elif strand == 'both':
			if j % 2 == 1:
				probeSeq = revcomp(seq[roundCoord:roundCoord+probeLength])
			else:
				probeSeq = seq[roundCoord:roundCoord+probeLength]
		else:
			parser.error('something went wrong with the strands')
	
		print '%s\t%s' % (counter, probeSeq)
		counter += 1
		i += stepSize
		j += 1
	
	carryOver = i - (len(seq)-probeLength)
#	print 'carryOver ', carryOver, stepSize
	fasta.close()
	os.remove(myID + "_region.fa")

